function PlotMuaMesh(CtrlVar,MUA,ElementList)

%%
% PlotMuaMesh(CtrlVar,MUA,ElementList)
%
% Examples:
%
% figure ; PlotMuaMesh(CtrlVar,MUA,ElementList)
%
% figure ; PlotMuaMesh([],MUA);  % CtrlVar is an optional input
%
% CtrlVar.NodeColor='r';
% figure ; PlotMuaMesh(CtrlVar,MUA,1:100);  % Plot nodes in red
%
% CtrlVar.PlotLabels=0;
% CtrlVar.PlotXYscale=1000;
% figure ; PlotMuaMesh(CtrlVar,MUA,1:100);  % Show only elemetns 1 to 100
%%

if isempty(CtrlVar)
    CtrlVar.PlotLabels=0;
    CtrlVar.MeshColor='k';
    CtrlVar.NodeColor='k';
    CtrlVar.PlotXYscale=1;
    CtrlVar.PlotNodesSymbolSize=3;
    CtrlVar.PlotNodesSymbol='o';
    CtrlVar.PlotNodes=1;
    CtrlVar.time=NaN;
    CtrlVar.FEmeshPlotTitle=[];
    CtrlVar.PlotFEmeshAndSaveMesh=0;
    CtrlVar.PlotsXaxisLabel='x';
    CtrlVar.PlotsYaxisLabel='y';
else




if ~isfield(CtrlVar,'WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo')
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
end

if nargin<3
    
    ElementList=1:MUA.Nele;
end

PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar,ElementList)



if CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo
    hold on
    PlotGmshGeometryDefinition(CtrlVar);
end

end