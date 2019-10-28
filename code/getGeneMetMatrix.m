function [GeneMetMatrix,metsConectivity,genesConectivity] = getGeneMetMatrix(model,genes)
%getGeneMetMatrix
%   
%   Function that obtains a binary matrix in which rows represent
%   metabolites and columns genes. Each non-zero coeffiecient represents
%   a relationship between a gene and a metabolite (i.e. the metabolite is
%   either consumed/produced by a reaction encoded by a given gene).
%
%   model            (struct) A MATLAB GEM structure
%   genes            List of gene IDs or gene indexes to take from the model 
%
%   GeneMetMatrix    (sparse) Boolean matrix representing the relationships 
%                    metabolites and genes
%   metsConectivity  (vector) Vector that indicates the amount of
%                    metabolites related to each gene
%   genesConectivity (vector) Indicates the amount of genes related to each
%                    metabolite
%
%   Usage:  [GeneMetMatrix, metsConectivity,genesConectivity] = getGeneMetMatrix(model,genes)
% 
%   Last modified.  Iván Domenzain 2019-10-28 
%

%Manage exceptions
if nargin<2
    genes = model.genes;
end

if ~isempty(genes)
    if ~isnumeric(genes)
        [iA,iB] = ismember(genes,model.genes);
        if numel(iA)~=numel(genes)
            disp('Not all provided genes were found in model')
        end
        genes = iB;
        genes = sort(genes);
    end
else
    genes = 1:numel(model.genes);
end
%Standardize rxnGeneMat
[~,rxnGeneMat]= standardizeGrRules(model,true);
%Get relevant variables size
M = sum((~contains(model.metNames,'prot_') & ~contains(model.metNames,'pmet_')));
Smat          = model.S;
G             = length(model.genes);
RGmat         = rxnGeneMat;
GeneMetMatrix = sparse(M,G);

for i=1:M
    %Get indexes of associated rxns for i-th metabolite
    rxnIndxs = find(Smat(i,:));
    %Get associated genes for each associated rxn
    for j=rxnIndxs
        if sum(RGmat(j,:))~=0
            geneIndxs                  = find(RGmat(j,:));
            GeneMetMatrix(i,geneIndxs) = 1;
        end
    end
end
GeneMetMatrix    = GeneMetMatrix(:,genes);
metsConectivity  = table(model.metNames(1:M),model.mets(1:M),model.metComps(1:M),sum(GeneMetMatrix,2),'VariableNames',{'metNames' 'mets' 'metComps' 'genes_number'});
genesConectivity = table(model.genes(iB),model.geneShortNames(iB),sum(GeneMetMatrix,1)','VariableNames',{'genes' 'short_names' 'mets_number'});
end

