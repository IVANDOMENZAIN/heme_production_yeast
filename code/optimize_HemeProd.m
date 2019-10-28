%function optimize_HemeProd
current       = pwd;
hemeModel     = getHeme_ecYeastGEM; 
tol           = 1E-12;
OE            = 2;
growthDeffect = 0.5;
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
%Identify candidate genes in model enzymes
[iA,iB]    = ismember(genes,hemeModel.enzGenes);
candidates = {};
for i=1:numel(iB)
    if iB(i)>0
        candidates = [candidates; hemeModel.enzymes(iB(i))];
    else
        candidates = [candidates; {''}];
    end
end
candidates = table(genes,candidates,geneShorts,actions,results.k_genes,'VariableNames',{'genes' 'enzymes' 'shortNames' 'actions' 'k_scores'});
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

%Get parsimonious proteome allocation solution
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = (1+tol)*WT_prod;
tempModel = setParam(tempModel, 'obj', pool_indxs, -1);
sol       = solveLP(tempModel,1);
pUsages   = sol.x(prot_indxs);
% WT_GUR    = sol.x(GUR_indx);

% %Get minimal GUR
% tempModel = setParam(tempModel, 'obj', GUR_indx, -1);
% sol       = solveLP(tempModel,1);
% WT_GUR    = sol.x(GUR_indx);
% pUsages   = sol.x(prot_indxs);

%Fix maximum production rate
tempFVA                = tempModel;
tempFVA.lb(targetIndx) = (1-tol)*WT_prod;
tempFVA.ub(targetIndx) = (1+tol)*WT_prod;
WT_bio_yield           = V_bio/(0.18*WT_GUR);
WT_prod_yield          = WT_prod/WT_GUR;
%% Run FVA for all enzyme usages subject to fixed GUR and Grates
ranges    = [];
minUsages = []; 
maxUsages = []; 
FVAprots  = candidates.enzymes;
candidateUsages = [];
for i=1:length(FVAprots)
    if ~isempty(FVAprots{i})
        rxnIndx = find(contains(tempFVA.rxnNames,FVAprots{i}));
        enzIndx = find(strcmpi(tempFVA.enzymes,FVAprots{i}));
        %Fix parsimonious usages
        %tempFVA = setParam(tempModel, 'lb', prot_indxs, (1-tol)*pUsages);
        %tempFVA = setParam(tempFVA, 'ub', prot_indxs, (1+tol)*pUsages);
        %Set i-th enzyme free
        %tempFVA = setParam(tempFVA, 'lb', rxnIndx, 0);
        %tempFVA = setParam(tempFVA, 'ub', rxnIndx, Inf);
        tempFVA = setParam(tempFVA, 'obj', rxnIndx, -1);
        sol     = solveLP(tempFVA);
        if ~isempty(sol.f)
            minFlux   = sol.x(rxnIndx); 
            tempFVA   = setParam(tempFVA, 'obj', rxnIndx, +1);
            sol       = solveLP(tempFVA);
            if ~isempty(sol.f)
               disp(['Ready with rxn #' num2str(i)])
               maxFlux = sol.x(rxnIndx); 
            else
               maxFlux = nan; 
            end
        else
            disp(['Nooot Ready with rxn #' num2str(i)])
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
%%
overExp_enzymes = overExp_enzymes';
overExp_enzymes  = overExp_enzymes*OE;
%overExp_enzymes(overExp_enzymes>0) = maxUsages(overExp_enzymes>0)./candidateUsages(overExp_enzymes>0);
candidates.OE     = overExp_enzymes;
candidates.usages = candidateUsages;
%Generate table with FVA results
t = table(candidates.enzymes,minUsages,maxUsages,ranges,candidateUsages,'VariableNames',{'enzNames' 'minUsages' 'maxUsages' 'ranges' 'pUsages'});
writetable(t,'../results/enzUsageRanges_hemeGenes.txt','Delimiter','\t','QuoteStrings',false);
%%
%Relevant rxn indexes
relIndexes = [GUR_indx, targetIndx];
tempModel  = hemeModel;
%Fix suboptimal biomass yield
tempModel.lb(growth_indx) = V_bio;
tempModel.lb(GUR_indx)    = 0;
tempModel.ub(GUR_indx)    = (1+tol)*1;
%Max heme production
tempModel    = setParam(tempModel, 'obj', targetIndx, +1);
%Get first genetic modification
[maxVal,gene,FC,validated] = testAllmutants(candidates,tempModel,relIndexes,WT_prod_yield);
%Discard production phenotype affecting genes
candidates = candidates(validated,:);
remaining  = candidates;
%Get rxn mets network
[GeneMetMatrix,Mconect,Gconect] = getGeneMetMatrix(tempModel,candidates.genes);
%Get linearly independent genes from GeneMetMatrix
[LI_Genes,EQ_Gmatrix,IndGenes] = getGeneDepMatrix(GeneMetMatrix);
%Append algebraic results to candidates table
candidates.LI          = LI_Genes;
candidates.Independent = IndGenes;
candidates.conectivity = Gconect.mets_number;
%GeneDepMat    = getGeneDepMatrix();
%Get top gene associated information
[mutantModel,I] = getGeneMutant(tempModel,candidates,gene);
remaining(I,:)  = [];
genesOpt        = gene;
FoldChanges     = FC;

%%
while ~isempty(remaining.genes)
    [maxVal,gene,FC]   = testAllmutants(remaining,mutantModel,relIndexes,WT_prod_yield);
    if ~isempty(gene)
        disp(['Optimal target found: ' gene ' FC: ' num2str(FC)])
        [mutantModel,I] = getGeneMutant(mutantModel,remaining,gene);
        remaining(I,:)  = [];
        genesOpt        = [genesOpt;gene];
        FoldChanges     = [FoldChanges; FC];
    else
        remaining.genes = [];
    end
end


%%
function [maxVal,gene,FC,positive] = testAllmutants(candidates,tempModel,indexes,WTval)
FoldChanges = [];
GUR_indx    = indexes(1);
targetIndx  = indexes(2);
%Index to minimize (bi-level optimization)
minIndex    = GUR_indx;
for i=1:height(candidates)
    gene     = candidates.genes{i};
    short    = candidates.shortNames{i};
    enzyme   = candidates.enzymes{i};
    enzUsage = candidates.usages(i);
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
positive  = FoldChanges>1;
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