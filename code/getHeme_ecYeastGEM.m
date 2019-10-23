function newModel = getHeme_ecYeastGEM
load('../models/ecYeastGEM.mat');
%Block alternative phosphofruktokinase:
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
end