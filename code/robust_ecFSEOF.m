function [mutantStrain,filtered] = robust_ecFSEOF(model,rxnTarget,expYield,CS_MW,resultsFolder)
mkdir(resultsFolder)
current      = pwd;
tol          = 1E-12;
OE           = 2;
thresholds   = [0.05 5];
% clone GECKO
git ('clone https://github.com/SysBioChalmers/GECKO')
cd GECKO
git checkout feat/add_FSEOF_utilities
%Get model parameters
cd geckomat
parameters = getModelParameters;
c_source   = parameters.c_source;
bioRXN     = parameters.bioRxn;
% Run FSEOF to find gene candidates
Nsteps    = 16;
alphaLims = [0.5*expYield 2*expYield];
cd utilities/ecFSEOF
mkdir('results')
file1   = 'results/genesResults_ecFSEOF.txt';
file2   = 'results/rxnsResults_ecFSEOF.txt';
results    = run_ecFSEOF(model,rxnTarget,c_source,alphaLims,Nsteps,file1,file2);
genes      = results.genes;
disp(['There are ' num2str(length(genes)) ' targets'])
geneShorts = results.geneNames;
actions    = results.k_genes;
actions(actions<0.5) = 0;
actions(actions>1)   = 1;
MWeigths             = [];
%Identify candidate genes in model enzymes
[~,iB]    = ismember(genes,model.enzGenes);
candidates = {};
for i=1:numel(iB)
    if iB(i)>0
        candidates = [candidates; model.enzymes(iB(i))];
        MWeigths   = [MWeigths; model.MWs(iB(i))];
    else
        candidates = [candidates; {''}];
        MWeigths   = [MWeigths; nan];
    end
end
%Get results files structures
candidates = table(genes,candidates,geneShorts,MWeigths,actions,results.k_genes,'VariableNames',{'genes' 'enzymes' 'shortNames' 'MWs' 'actions' 'k_scores'});
candidates = candidates(((candidates.actions==1)|(candidates.actions==0)),:); 
cd (current)
disp('First filter OE or DR?')
disp(['There are ' num2str(height(candidates)) ' targets'])
writetable(candidates,[resultsFolder '/candidates_ecFSEOF.txt'],'Delimiter','\t','QuoteStrings',false);
% Get constraints values
tempModel   = model;
%Get relevant rxn indexes
targetIndx  = find(strcmpi(tempModel.rxns,rxnTarget));
CUR_indx    = find(strcmpi(tempModel.rxnNames,c_source));
growth_indx = find(strcmpi(tempModel.rxns,bioRXN));
%Fix suboptimal experimental biomass yield conditions
Yield = expYield;
V_bio = Yield*CS_MW;
tempModel.lb(growth_indx) = V_bio;
%Fix unit C source uptake
tempModel.lb(CUR_indx)    = (1-tol)*1;
tempModel.ub(CUR_indx)    = (1+tol)*1;
%Get and fix optimal production rate
tempModel = setParam(tempModel, 'obj', targetIndx, +1);
sol       = solveLP(tempModel,1);
WT_prod   = sol.x(targetIndx);
WT_CUR    = sol.x(CUR_indx);
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = (1+tol)*WT_prod;
%Calculate WT yields
WT_prod_yield = WT_prod/WT_CUR;
% Run FVA for all enzyme usages subject to fixed CUR and Grates
disp(' ')
disp('Running enzyme usage variability analysis')
FVAtable = enzymeUsage_FVA(tempModel,candidates.enzymes);
candidateUsages = FVAtable.pU;
minUsages       = FVAtable.minU;
maxUsages       = FVAtable.maxU;
% Identify enzymes with "room" for direct overexpression
for i=1:length(candidates.enzymes)
    if maxUsages(i)~=0 && candidates.actions(i)>0
        candidates.actions(i) = 1;
        %For those enzymes which flexib
        if maxUsages(i)< OE*candidateUsages(i)
            candidates.actions(i) = 2;
        end 
    end
end
candidates.OE(candidates.actions>0)  = OE;
candidates.OE(candidates.actions==0) = 0;
candidates.minUsage = minUsages;
candidates.maxUsage = maxUsages;
candidates.pUsage   = candidateUsages;
candidates.pUsage   = candidateUsages;
%Generate table with FVA results
t = table(candidates.enzymes,minUsages,maxUsages,FVAtable.ranges,candidateUsages,'VariableNames',{'enzNames' 'minUsages' 'maxUsages' 'ranges' 'pUsages'});
writetable(candidates,[resultsFolder '/candidates_enzUsageFVA.txt'],'Delimiter','\t','QuoteStrings',false);
%Discard enzymes whose usage LB = UB = 0
tempMat  = table2array(t(:,2:3));
unused   = find(sum(tempMat,2)==0);
toRemove = intersect(unused,find(candidates.actions>0));
candidates(toRemove,:) = [];
disp(['There are ' num2str(height(candidates)) ' candidates'])
% Mechanistic validations of FSEOF results
disp(' ')
disp('Mechanistic validation of results')
%Relevant rxn indexes
relIndexes = [CUR_indx, targetIndx];
%relax target rxn bounds
tempModel.lb(targetIndx)  = 0;
tempModel.ub(targetIndx)  = 1000;
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',targetIndx,+1);
[~,~,FCs,validated]  = testAllmutants(candidates,tempModel,relIndexes,WT_prod_yield,1E-6);
%Discard genes with a negative impact on production yield
candidates.foldChange = FCs; 
candidates            = candidates(validated,:);
disp(['There are ' num2str(height(candidates)) ' candidates'])
writetable(candidates,[resultsFolder '/candidates_mech_validated.txt'],'Delimiter','\t','QuoteStrings',false);
% Assess linear dependencies
%Get rxn mets network
[GeneMetMatrix,~,Gconect] = getGeneMetMatrix(tempModel,candidates.genes);
%Get linearly independent genes from GeneMetMatrix
[LI_Genes,G2Gmatrix,~] = getGeneDepMatrix(GeneMetMatrix);
%Append algebraic results to candidates table
candidates.LI          = LI_Genes;
candidates.conectivity = Gconect.mets_number;
% Keep top results
toKeep            = find((candidates.k_scores>=thresholds(2)|candidates.k_scores<=thresholds(1)));
candidates        = candidates(toKeep,:);
GeneMetMatrix     = GeneMetMatrix(:,toKeep);
G2Gmatrix         = G2Gmatrix(toKeep,toKeep);
[~,groups]        = getGenesGroups(G2Gmatrix,candidates.genes);
candidates.groups = groups;
disp(['There are ' num2str(height(candidates)) ' candidates'])
% Rank candidates by priority
%%% 1st. LI=1 OEs with both min and pUsage>0 & Deletions with pUsage=0
priority = zeros(height(candidates),1);
%LI Enzymes that are necesarily used
cond1    = (candidates.actions>0 & candidates.pUsage>0 & candidates.minUsage>0);
%LI enzymes that are not used
cond2    = (candidates.actions==0 & candidates.pUsage==0 & candidates.maxUsage>0);
indexes  = (candidates.LI==1 & (cond2 | cond1));
priority(indexes) = 1; 
%%% 2nd. LI=0, for OEs pick the enzyme with the lowest MW for each group 
for i=1:max(candidates.groups)
    %Find group genes
    groupIndxs = find(candidates.groups==i);
    groupTable = candidates(groupIndxs,:); 
    if sum(candidates.actions(groupIndxs))==0
        nonZeros = (candidates.pUsage(groupIndxs)>0);
        if sum(nonZeros)==0
            priority(groupIndxs) = 2;
        else
            priority(groupIndxs(~nonZeros)) = 3;
        end
    else
        if all(candidates.actions(groupIndxs)>0)
            %Select those with pUsage>0  and minUsage>0 and the lowest MW
            groupTable = groupTable((groupTable.pUsage>0),:);
            [~,minMW]  = min(groupTable{:,'MWs'});
            groupGenes = groupTable.genes(minMW);
            if ~isempty(groupGenes)                
                prtyIndx   = (strcmpi(candidates.genes,groupGenes));
                priority(prtyIndx) = 1;
            end
        end
    end
end
candidates.priority = priority;
%Keep priority genes and sort them accordingly
candidates = candidates(priority>0,:);
disp(['There are ' num2str(height(candidates)) ' candidates'])
candidates = sortrows(candidates,'priority','ascend');
writetable(candidates,[resultsFolder '/candidates_priority.txt'],'Delimiter','\t','QuoteStrings',false);
% get optimal strain according to priority candidates
[~,filtered,~,iB] = getOptimalStrain(tempModel,candidates,[targetIndx CUR_indx],CS_MW);
[mutantStrain,filtered,] = getOptimalStrain(tempModel,filtered,[targetIndx CUR_indx],CS_MW);
cd (current)
actions = cell(height(filtered),1);
actions(filtered.actions==0)= {'deletion'};
actions(filtered.actions>0) = {'OE'};
filtered.actions = actions;
writetable(filtered,[resultsFolder '/compatible_genes_results.txt'],'Delimiter','\t','QuoteStrings',false);
origin = 'GECKO/geckomat/utilities/ecFSEOF/results/*';
copyfile(origin,resultsFolder)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxVal,TOPgene,FoldChanges,positive] = testAllmutants(candidates,tempModel,indexes,WTval,tol)
FoldChanges = [];
CUR_indx    = indexes(1);
targetIndx  = indexes(2);
medianUsage = (candidates.maxUsage-candidates.minUsage)/2; 
%Index to minimize (bi-level optimization)
minIndex    = find(contains(tempModel.rxnNames,'prot_pool'));
for i=1:height(candidates)
    gene     = candidates.genes{i};
    short    = candidates.shortNames{i};
    action   = candidates.actions(i);
    OEf      = candidates.OE(i);
    modifications   = {gene action OEf};
    if action == 0
        pUsage = medianUsage(i);
    else
        pUsage = [];
    end
    mutantModel     = getMutant(tempModel,modifications,pUsage);
    [mutSolution,~] = solveECmodel(mutantModel,mutantModel,'pFBA',minIndex);
    if ~isempty(mutSolution)
        yield = mutSolution(targetIndx)/mutSolution(CUR_indx);
        FC    = yield/WTval;
    else
        FC = 0;
    end
    FoldChanges = [FoldChanges; FC];
    %disp(['Ready with genetic modification #' num2str(i) '[' short ': ' num2str(action) '] FC: ' num2str(FC)])
end
positive   = FoldChanges>(1-tol);
[maxVal,I] = max(FoldChanges);
if ~(maxVal>1)
    maxVal = [];
    gene   = [];
else 
    TOPgene = candidates.genes{I};
    FC   = FoldChanges(I);
    %disp(['candidate gene: ' short ' FC: ' num2str(FC)])
end
end

