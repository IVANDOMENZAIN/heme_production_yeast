%function getHeme_ecYeastGEM(ecModel_batch)
%Block alternative phosphofruktokinase:
load('../models/ecYeastGEM.mat');
ecModel_batch.ub(strcmp(ecModel_batch.rxns,'r_0887')) = 0; %ATP + sedohept-7P -> ADP + H+ + sedohept-1,7biP
%Add heme reaction:
posH                   = strcmp(ecModel_batch.mets,'s_3714'); %heme a [cytoplasm]
rxnsToAdd.mets         = {'s_3714'};                    %Already existing in yeastGEM
rxnsToAdd.rxns         = {'r_s3714_Ex'};
rxnsToAdd.rxnNames     = {'heme exchange'};
rxnsToAdd.stoichCoeffs = {-1};
rxnsToAdd.lb           = 0;
rxnsToAdd.ub           = 1000;
rxnsToAdd.c            = 0;
rxnsToAdd.subSystems   = {''};
rxnsToAdd.grRules      = {''};
newModel               = addRxns(ecModel_batch,rxnsToAdd,1,'c',false);
rxnTarget               = 'r_s3714_Ex';
cSource   = 'D-glucose exchange (reversible)';
alphaLims = [0.2 0.8];
Nsteps    = 16;
git ('clone https://github.com/SysBioChalmers/GECKO')
cd GECKO
git checkout feat/add_FSEOF_utilities
cd geckomat/utilities/ecFSEOF
newFolder = '../../../../../results';
mkdir(newFolder)
file1 = [newFolder '/genesResults.txt'];
file2 = [newFolder '/rxnsResults.txt'];
results = run_ecFSEOF(newModel,rxnTarget,cSource,alphaLims,Nsteps,file1,file2);
%end