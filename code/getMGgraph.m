function [metGeneGraph,metNames] = getMGgraph(GeneMetMatrix,mets,model,geneShorts,MWs,algorithm)
exclude   = {'ATP' 'ADP' 'AMP' 'GTP' 'GDP' 'GMP' 'CTP' 'CMP' 'CDP' 'oxygen' 'acetyl-CoA' 'acetaldehyde' 'ethanol' ...
             'dGDP' 'dCDP' 'dATP' 'dGTP' 'dTTP' 'dTMP' 'dGMP' 'dAMP' 'dADP' 'dUMP' 'dUTP' 'diphosphate' 'phosphate' 'H+' 'H2O' 'NADPH' 'NADP(+)' 'NADH' 'NAD' 'pmet_' 'carbon dioxide' 'coenzyme A'};
tempGMmat = GeneMetMatrix;
[~,iB]    = ismember(mets,model.mets);
metNames  = model.metNames(iB);

if nargin<5
    algorithm = [];
end
%Exclude highly connected metabolites from GeneMetMatrix
for compound=exclude
    if strcmpi(compound,'pmet_')
        metIndxs = find(contains(metNames,compound));
    else
        metIndxs = find(strcmpi(metNames,compound));
    end
    metNames(metIndxs)    = [];
    mets(metIndxs)        = [];
    tempGMmat(metIndxs,:) = [];
    iB(metIndxs) = [];
end

for i=1:length(iB)
    index = iB(i);
    comp = model.metComps(index);
    str  = ['_[' num2str(comp) ']'];
    metNames{i} = [metNames{i} str];
end

[M,G]        = size(tempGMmat);
metGeneGraph = zeros(M,M);
edgeNames    = {};
for i=1:G
    geneRow  = tempGMmat(:,i);
    metIndxs = find(geneRow);
    Mweigth  = MWs(i);
    geneS    = geneShorts(i);
    for j=1:length(metIndxs)
        index        = metIndxs(j);
        tempIndxs    = metIndxs;
        tempIndxs(j) = [];
        %if metGeneGraph(metIndxs,index)==0
            edgeNames = [edgeNames,geneS];
        %end
        metGeneGraph(metIndxs,index) = i;%metGeneGraph(metIndxs,index)+Mweigth;
    end
end

metGeneGraph = graph(metGeneGraph,mets,'OmitSelfLoops');
if ~isempty(algorithm)
    figure                 % Creates a figure
    set(gca,'FontSize',26) % Creates an axes and sets its FontSize to 18
    p = plot(metGeneGraph);
    labelnode(p,1:height(metGeneGraph.Nodes),metNames)
    layout(p,'force')
    length(edgeNames)
    numedges(metGeneGraph)
    labeledge(p,1:numedges(metGeneGraph),(metGeneGraph.Edges.Weight))
end
end