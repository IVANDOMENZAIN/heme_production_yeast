%Clone GECKO and go to the right branch
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO 
git checkout feat/add_FSEOF_utilities
cd ..
clc
%Load model
load('../models/ecYeastGEM_batch.mat')
model = getHeme_ecYeastGEM(ecModel_batch);
targetIndex = find(strcmpi(model.rxnNames,'heme exchange'));
%Set media conditions
cd GECKO/geckomat/kcat_sensitivity_analysis
c_source = 'D-glucose exchange (reversible)';
tempModel = changeMedia_batch(model,c_source);
clc
CS_MW     = 0.18015;
CS_index  = find(strcmpi(tempModel.rxnNames,c_source));
growthPos = find(strcmpi(tempModel.rxnNames,'growth'));
cd ../../../
%Unconstrain growth
tempModel = setParam(tempModel,'lb',growthPos,0);
tempModel = setParam(tempModel,'ub',growthPos,1000);
%Get biomass yield for a unit glucose uptake rate 
tempModel = setParam(tempModel,'obj',growthPos,1);
tempModel = setParam(tempModel,'ub',CS_index,1);
solution  = solveLP(tempModel,1);
WT_yield  = solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
disp(['The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);
resultsFolder = '../results/ecModel_targets';
mkdir(resultsFolder)
expYield = 0.122;
[optStrain,candidates] = robust_ecFSEOF(tempModel,tempModel.rxns(targetIndex),0.122,CS_MW,resultsFolder);
