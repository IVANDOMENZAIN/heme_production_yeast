function [optStrain,optGenes,FChanges,iB] = getOptimalStrain(model,candidates,rxnIndxs,CS_MW)
%Get WT yield
sol = solveLP(model,1);
if ~isempty(sol.x)
    targetIndex = rxnIndxs(1);
    GURindex    = rxnIndxs(2);
    WTyield     = sol.x(targetIndex)/(sol.x(GURindex)*CS_MW);
end
medianUsage = (candidates.maxUsage-candidates.minUsage)/2; 
%Create mutants iteratively
optStrain  = model;
FChanges   = [];
genesFC    = [];
counter    = 0;
previousFC = 1;
for i=[1 2 3]
    levelCandidates = candidates(candidates.priority==i,:);
    levelCandidates = sortrows(levelCandidates,'foldChange','descend');
    for j=1:length(levelCandidates.genes)
        gene   = levelCandidates.genes{j};
        short  = levelCandidates.shortNames(j);
        action = levelCandidates.actions(j);
        OEf    = levelCandidates.OE(j);
        modifications = {gene action OEf};
        if action == 0
            enzUsage = medianUsage(i);
        else
            enzUsage = [];
        end
        tempMutant    = getMutant(optStrain,modifications,enzUsage);
        [mutSol,~]    = solveECmodel(tempMutant,model,'pFBA',GURindex);
        if ~isempty(mutSol)
            yield = mutSol(targetIndex)/(mutSol(GURindex)*CS_MW);
            FC    = yield/WTyield;
            %Just keep those genes that don't affect the production phenotype
            if (FC-previousFC)>-1E-6
                FChanges   = [FChanges; FC];
                genesFC    = [genesFC;gene];
                optStrain  = tempMutant;
                previousFC = FC;
                counter = counter+1;
                disp(['Ready with gene #' num2str(counter) ' (' short{1} ')' '  FC:' num2str(FC)])
            end
        else
            FC = 0;
        end
    end
end
[~,iB]   = ismember(genesFC,candidates.genes);
iB       = sort(iB,'ascend');
FChanges = table(genesFC,FChanges,'VariableNames',{'genes' 'FC'});
optGenes = candidates(iB,:);
end