%function optimize_HemeProd
current      = pwd;
hemeModel    = getHeme_ecYeastGEM; 
tol          = 1E-12;
OE           = 2;
fID          = fopen('../data/previous_results.txt');
prev_results = textscan(fID,'%s %s %f','Delimiter','\t','HeaderLines',1,'TreatAsEmpty',{'NA','na','NaN'});
expValidated = prev_results{1}(prev_results{3}==1);
modelTargets = prev_results{1};
%% clone GECKO
%git ('clone https://github.com/SysBioChalmers/GECKO')
cd GECKO
git checkout feat/add_FSEOF_utilities
%Set minimal media
cd geckomat/kcat_sensitivity_analysis
hemeModel = changeMedia_batch(hemeModel,'D-glucose exchange (reversible)');
%% Run FSEOF to find gene candidates
cd ../utilities/ecFSEOF
newFolder = '../../../../../results';
mkdir(newFolder)
file1     = [newFolder '/genesResults.txt'];
file2     = [newFolder '/rxnsResults.txt'];
rxnTarget = 'r_s3714_Ex';
cSource   = 'D-glucose exchange (reversible)';
alphaLims = [0.01098 0.04392];
Nsteps    = 16;

results              = run_ecFSEOF(hemeModel,rxnTarget,cSource,alphaLims,Nsteps,file1,file2);
genes                = results.genes;
geneShorts           = results.geneNames;
actions              = results.k_genes;
actions(actions<0.5) = 0;
actions(actions>1)   = 1;
actions              = logical(actions);
MWeigths             = [];
%Identify candidate genes in model enzymes
[iA,iB]    = ismember(genes,hemeModel.enzGenes);
candidates = {};
for i=1:numel(iB)
    if iB(i)>0
        candidates = [candidates; hemeModel.enzymes(iB(i))];
        MWeigths   = [MWeigths; hemeModel.MWs(iB(i))];
    else
        candidates = [candidates; {''}];
        MWeigths   = [MWeigths; nan];
    end
end
candidates = table(genes,candidates,geneShorts,MWeigths,actions,results.k_genes,'VariableNames',{'genes' 'enzymes' 'shortNames' 'MWs' 'actions' 'k_scores'});
%% Get constraints values
cd (current)
tempModel   = hemeModel;
heme_yield  = [];
%Get relevant rxn indexes
targetIndx  = find(strcmpi(tempModel.rxnNames,'heme exchange'));
GUR_indx    = find(strcmpi(tempModel.rxnNames,'D-glucose exchange (reversible)'));
growth_indx = find(strcmpi(tempModel.rxnNames,'growth'));
prot_indxs  = find(contains(tempModel.rxnNames,'prot_'));
pool_indxs  = prot_indxs(end);
prot_indxs  = prot_indxs(1:end-1);
%Fix suboptimal experimental biomass yield conditions
Yield                     = 0.122;
V_bio                     = Yield*0.18;
tempModel.lb(growth_indx) = V_bio;
tempModel.lb(GUR_indx)    = (1-tol)*1;
tempModel.ub(GUR_indx)    = (1+tol)*1;
%Max heme production
tempModel = setParam(tempModel, 'obj', targetIndx, +1);
sol       = solveLP(tempModel,1);
WT_prod   = sol.x(targetIndx);
pUsages   = sol.x(prot_indxs);
WT_GUR    = sol.x(GUR_indx);
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = (1+tol)*WT_prod;
%Calculate WT yields
WT_bio_yield  = V_bio/(0.18*WT_GUR);
WT_prod_yield = WT_prod/WT_GUR;
%% Run FVA for all enzyme usages subject to fixed GUR and Grates
ranges    = [];
minUsages = []; 
maxUsages = []; 
FVAprots  = candidates.enzymes;
tempFVA   = tempModel;
tempFVA.lb(growth_indx) = 0;
tempFVA.lb(GUR_indx)    = 0;
tempFVA.ub(GUR_indx)    = 1000;
candidateUsages = [];
for i=1:length(FVAprots)
    if ~isempty(FVAprots{i})
        rxnIndx = find(contains(tempFVA.rxnNames,FVAprots{i}));
        enzIndx = find(strcmpi(tempFVA.enzymes,FVAprots{i}));
        tempFVA = setParam(tempFVA, 'obj', rxnIndx, -1);
        sol     = solveLP(tempFVA);
        if ~isempty(sol.f)
            minFlux = sol.x(rxnIndx); 
            tempFVA = setParam(tempFVA, 'obj', rxnIndx, +1);
            sol     = solveLP(tempFVA);
            if ~isempty(sol.f)
               disp(['Ready with rxn #' num2str(i)])
               maxFlux = sol.x(rxnIndx); 
            else
               maxFlux = nan; 
            end
        else
            minFlux = nan;
            maxFlux = nan; 
        end
        ranges          = [ranges; (maxFlux-minFlux)];
        minUsages       = [minUsages; minFlux]; 
        maxUsages       = [maxUsages; maxFlux];
        candidateUsages = [candidateUsages;pUsages(enzIndx)];
    else
        ranges    = [ranges; 0];
        minUsages = [minUsages; 0]; 
        maxUsages = [maxUsages; 0];
        candidateUsages = [candidateUsages;0];
    end
end
%Identify enzymes with "room" for direct overexpression
for i=1:length(FVAprots)
    if maxUsages(i)~=0 & actions(i)==1
        overExp_enzymes(i) = 1;
        %For those enzymes which flexib
        if maxUsages(i)< OE*candidateUsages(i)
            actions(i) = 2;
        end
    else
        overExp_enzymes(i) = 0; 
    end
end
%% write FVA results on enzyme usage rxns on a file and filter out unused enzymes
overExp_enzymes     = overExp_enzymes';
overExp_enzymes     = overExp_enzymes*OE;
candidates.OE       = overExp_enzymes;
candidates.minUsage = minUsages;
candidates.maxUsage = maxUsages;
candidates.pUsage   = candidateUsages;
%Generate table with FVA results
t = table(candidates.enzymes,minUsages,maxUsages,ranges,candidateUsages,'VariableNames',{'enzNames' 'minUsages' 'maxUsages' 'ranges' 'pUsages'});
writetable(t,'../results/enzUsageRanges_hemeGenes.txt','Delimiter','\t','QuoteStrings',false);
%Discard unusable enzymes for OE
tempMat  = table2array(t(:,2:end));
unused   = find(sum(tempMat,2)==0);
toRemove = intersect(unused,find(candidates.actions==1));
candidates(unused,:) = [];
%Keep those enzymes 
%% Mechanistic validations of FSEOF results
%Relevant rxn indexes
relIndexes = [GUR_indx, targetIndx];
tempModel  = hemeModel;
%Fix suboptimal biomass yield
tempModel.lb(growth_indx) = V_bio;
tempModel.lb(GUR_indx)    = 0;
tempModel.ub(GUR_indx)    = (1+tol)*1;
%set Max heme production as objective function
tempModel = setParam(tempModel,'obj',targetIndx,+1);
%validate candidates in a mechanistic way
[~,~,~,validated] = testAllmutants(candidates,tempModel,relIndexes,WT_prod_yield,tol);
%Discard production phenotype affecting genes
candidates = candidates(validated,:);
%% Assess linear dependencies
%Get rxn mets network
[GeneMetMatrix,Mconect,Gconect] = getGeneMetMatrix(tempModel,candidates.genes);
%Get linearly independent genes from GeneMetMatrix
[LI_Genes,EQ_Gmatrix,~] = getGeneDepMatrix(GeneMetMatrix);
%Append algebraic results to candidates table
candidates.LI          = LI_Genes;
candidates.conectivity = Gconect.mets_number;
%% Keep top results
candidates.k_scores(candidates.k_scores<0.01) = 0;
candidates2 = candidates;
%toKeep        = find((candidates.k_scores==1000|candidates.k_scores==0));
toKeep        = find((candidates.k_scores==1000));
candidates    = candidates(toKeep,:);
GeneMetMatrix = GeneMetMatrix(:,toKeep);
EQ_Gmatrix    = EQ_Gmatrix(toKeep,toKeep);
[~,groups]    = getGenesGroups(EQ_Gmatrix,candidates.genes);
candidates.groups = groups;
candidates.expVal = ismember(candidates.genes,expValidated);
%% Rank candidates by priority
%%% 1st. LI=1 OE's with both min and pUsage>0 & Del's with pUsage=0
priority = zeros(height(candidates),1);
%LI Enzymes that are necesarily used
cond1    = (candidates.actions==1 & candidates.pUsage>0 & candidates.minUsage>0);
%LI enzymes that are not used
cond2    = (candidates.actions==0 & candidates.pUsage==0 & candidates.maxUsage>0);
indexes  = find(candidates.LI==1 & (cond2 | cond1));
priority(indexes) = 1; 

%%% 2nd. LI=0, for OE's pick the enzyme with the lowest MW for each group 
for i=1:max(candidates.groups)
    %Find group genes
    groupIndxs = find(candidates.groups==i);
    groupGenes = candidates.genes(groupIndxs);
    groupTable = candidates(groupIndxs,:); 
    if sum(candidates.actions(groupIndxs))==0
        nonZeros = (candidates.pUsage(groupIndxs)>0);
        if sum(nonZeros)==0
            priority(groupIndxs) = 2;
        else
            priority(groupIndxs(~nonZeros)) = 3;
        end
    else
        if sum(candidates.actions(groupIndxs)==1) == length(groupIndxs)
            %Select those with pUsage>0  and minUsage>0 and the lowest MW
            cond1      = candidates.pUsage(groupIndxs)>0;
            cond2      = candidates.minUsage(groupIndxs)>0;
            groupTable = groupTable(find(cond1 & cond2),:);
            [~,minMW]  = min(groupTable.MWs);
            groupGenes = groupGenes(minMW);
            if ~isempty(groupGenes)
                prtyIndx   = find(strcmpi(candidates.genes,groupGenes));
                priority(prtyIndx) = 1;
            end
        end
    end
end
candidates.priority = priority;
%% Get met-met digraph
MetsIndxs = find(~contains(tempModel.metNames,'prot_'));
nodeMets  = tempModel.mets(MetsIndxs);
sumCols   = sum(GeneMetMatrix,2);
toKeep    = find(sumCols);
nodeMets  = nodeMets(toKeep);
tempGMmatrix = GeneMetMatrix(toKeep,:);

[metGeneGraph,metNames] = getMGgraph(tempGMmatrix,nodeMets,tempModel,candidates.shortNames,candidates.MWs,'force');
%%
function [maxVal,gene,FC,positive] = testAllmutants(candidates,tempModel,indexes,WTval,tol)
FoldChanges = [];
GUR_indx    = indexes(1);
targetIndx  = indexes(2);
%Index to minimize (bi-level optimization)
minIndex    = GUR_indx;
for i=1:height(candidates)
    gene     = candidates.genes{i};
    short    = candidates.shortNames{i};
    enzyme   = candidates.enzymes{i};
    enzUsage = candidates.pUsage(i);
    action   = candidates.actions(i);
    OEf      = candidates.OE(i);
    modifications   = {gene action OEf};
    mutantModel     = getMutant(tempModel,modifications);
    [mutSolution,~] = solveECmodel(mutantModel,mutantModel,'pFBA',minIndex);
    if ~isempty(mutSolution)
        yield = mutSolution(targetIndx)/mutSolution(GUR_indx);
        FC    = yield/WTval;
    else
        FC = 0;
    end
    FoldChanges = [FoldChanges; FC];
    disp(['Ready with genetic modification #' num2str(i) ' (' num2str(action) ') FC: ' num2str(FC)])
end
positive  = FoldChanges>(1+tol);
[maxVal,I] = max(FoldChanges);
if ~(maxVal>1)
    maxVal = [];
    gene   = [];
    FC     = [];
else 
    gene = candidates.genes{I};
    FC   = FoldChanges(I);
    disp(['candidate gene: ' gene ' FC: ' num2str(FC)])
end
end
%%%
function [mutantModel,I] = getGeneMutant(tempModel,candidates,gene)
    I             = find(strcmpi(candidates.genes,gene));
    modifications = {candidates.genes{I} candidates.actions(I) candidates.OE(I)};
    mutantModel   = getMutant(tempModel,modifications);
end