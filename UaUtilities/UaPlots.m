function [cbar,xGL,yGL,xCF,yCF]=UaPlots(CtrlVar,MUA,F,Variable,options)

%%
%
% cbar=UaPlots(CtrlVar,MUA,F,Variable,options)
%
% Simple plot utility to plot variables and calving fronts and grounding lines as well.
%
% Returns grounding lines (xGL,yGL) and calving fronts (xCF,yCF).
%
% Calving fronts
%
%
% Examples:
%
%   UaPlots(CtrlVar,MUA,F,F.h)
%
%   UaPlots(CtrlVar,MUA,F,"-speed-")
%
%   UaPlots(CtrlVar,MUA,F,"-ubvb-")
%
%   cbar=UaPlots(CtrlVar,MUA,F,F.h,CalvingFrontColor="b",GroundingLineColor="k",GetRidOfValuesDownStreamOfCalvingFronts=false,ColorMap=jet);
%   title(cbar,"h (m)")
%
%   UaPlots(CtrlVar,MUA,F,"-log10speed-",CalvingFrontColor="b",GroundingLineColor="k",GetRidOfValuesDownStreamOfCalvingFronts=false,ColorMap=othercolor("YlGnBu8",100),PlotUnderMesh=true) ;
%
%
%   figetaInt=FindOrCreateFigure("eta Int") ; clf(figetaInt) ;
%   UaPlots(CtrlVar,MUA,F,"eta int",GetRidOfValuesDownStreamOfCalvingFronts=false) ;
%
%%

arguments
    CtrlVar struct
    MUA     struct
    F       {mustBeA(F,{'struct','UaFields','numeric'})}
    Variable   {string, double}
    options.PlotGroundingLines  logical = true
    options.PlotCalvingFronts  logical = true
    options.CalvingFrontColor char = "k"
    options.GroundingLineColor char = "r"
    options.GetRidOfValuesDownStreamOfCalvingFronts=true;
    options.PlotOverMesh=false;
    options.PlotUnderMesh=false;
    options.PlotMuaBoundary=true;
    options.FigureTitle string="UaPlots";  % this is the figure title, not the plot title 


    % options.ColorMap double=othercolor('YlGnBu6',1028)
    % options.ColorMap double=othercolor("Mlightterrain",1028)
    % options.ColorMap double=othercolor("Mdarkterrain",1028)
    % options.ColorMap double=othercolor("Mtemperaturemap",1028)
    options.ColorMap double=othercolor("YlGnBu8",1028)  % See othercolor.m for more options
end

% if fig title has not been set, use by default the variable name
if options.FigureTitle=="UaPlots"
    if isstring(Variable)
       options.FigureTitle=Variable;
    elseif ~isempty(inputname(4))
        options.FigureTitle=inputname(4) ; 
    end

end


fFig=FindOrCreateFigure(options.FigureTitle)  ; clf(fFig)  ; 

if islogical(Variable)
    Variable=double(Variable) ;
end

if isnumeric(Variable)
    Variable=Variable(:);
end

if isempty(F)
    F=UaFields;
end

if isscalar(Variable)

    if contains(Variable,"int")
        options.GetRidOfValuesDownStreamOfCalvingFronts=false;
    end

end



if options.GetRidOfValuesDownStreamOfCalvingFronts  && ~isempty(F.LSF)

    if isempty(F.LSFMask)
        F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    end

    F.ub(~F.LSFMask.NodesIn)=NaN;
    F.vb(~F.LSFMask.NodesIn)=NaN;
    if isnumeric(Variable)
        if numel(Variable)==MUA.Nnodes
            Variable(~F.LSFMask.NodesIn)=NaN;
        end
    end

end


xGL=nan ; yGL=nan ; xCF=nan ; yCF=nan ;

isModifyColormap=true;

if isModifyColormap
    colormap(options.ColorMap);
end

if options.PlotOverMesh
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    PlotMuaMesh(CtrlVar,MUA) ;
    hold on

end


if isnumeric(Variable)

    [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Variable);
    title(cbar,inputname(4)) ;

else


    switch lower(Variable)

        case {"speed","-speed-"}

            speed=sqrt(F.ub.*F.ub+F.vb.*F.vb) ;

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\| \mathbf{v} \|$",Interpreter="latex")
            title(cbar,"$(\mathrm{m\,yr^{-1}})$",interpreter="latex")

        case {"log10speed","-log10speed-"}

            speed=log10(sqrt(F.ub.*F.ub+F.vb.*F.vb)) ;
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,speed);
            title("$\log_{10}(\| \mathbf{v} \|)$",Interpreter="latex")
            title(cbar,"$\log_{10}(m/a)$",Interpreter="latex")


        case {"ubvb","-ubvb-","uv"}

            CtrlVar.VelColorMap=jet(100) ;
            cbar=QuiverColorGHG(F.x,F.y,F.ub,F.vb,CtrlVar) ;
            title(cbar,"(m/a)",Interpreter="latex")
            title(sprintf("velocities at t=%f",CtrlVar.time),Interpreter="latex")


        case "dhdt"


            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.dhdt);
            title(cbar,"(m/a)",Interpreter="latex")
            title(sprintf("$dh/dt$ at t=%f",CtrlVar.time),Interpreter="latex")
            title(cbar,"$(\mathrm{m\,yr^{-1}})$",interpreter="latex")

        case "basal drag"




            [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,[],MUA,F) ;

            CtrlVar.VelColorMap=jet(100) ;
            cbar=QuiverColorGHG(F.x,F.y,txzb,tyzb,CtrlVar) ;
            title(cbar,"(m/a)",Interpreter="latex")
            title(sprintf("basal drag vectors at t=%f",CtrlVar.time),Interpreter="latex")


        case "e"  % effective strain rate


            [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,[],MUA,F) ;

            % e(e<0)=eps ; % the projection onto nodes does not preserve positivy
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,e);
            title(cbar,"(1/a)",Interpreter="latex")
            title(sprintf("effective strain rates at t=%f",CtrlVar.time),Interpreter="latex")


        case "e int"  % effective strain rate at integration points



            [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); % returns integration point values

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,e);
            title(cbar,"(1/a)",Interpreter="latex")
            title(sprintf("effective strain rates at integration points at t=%f",CtrlVar.time),Interpreter="latex")



        case "eta"  % effective strain rate


            [txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e,eta]=CalcNodalStrainRatesAndStresses(CtrlVar,[],MUA,F) ;

            % e(e<0)=eps ; % the projection onto nodes does not preserve positivy
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,eta);
            title(cbar,"(kPa yr)",Interpreter="latex")
            title(sprintf("effective viscosity eta at t=%f",CtrlVar.time),Interpreter="latex")


        case "eta int"  % effective strain rate at integration points



            [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); % returns integration point values


            fFigHist=FindOrCreateFigure(options.FigureTitle+"Hist")  ; clf(fFigHist)  ; 
            histogram((log10(etaInt(:))),Normalization="probability") ; 
            hold on ; 
            xline(log10(CtrlVar.etaZero),'r',LineWidth=2)


            fFig=FindOrCreateFigure(options.FigureTitle)  ; clf(fFig)  ; 
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,log10(etaInt));
            title(cbar,"(kPa yr)",Interpreter="latex")
            title(sprintf("log10 of effective viscosity at integration points at t=%f",CtrlVar.time),Interpreter="latex")



    

        case "surface slope"  % effective strain rate at integration points

            [dfdx,dfdy,xint,yint]=calcFEderivativesMUA(F.s,MUA,CtrlVar) ;
            slope=sqrt(dfdx.*dfdx+dfdy.*dfdy) ;
            slope=ProjectFintOntoNodes(MUA,slope) ;
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.rho.*F.h.*slope);
            title(cbar,"()",Interpreter="latex")
            title(sprintf("surface slope at t=%f",CtrlVar.time),Interpreter="latex")


        otherwise

            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.(Variable));
            title(cbar,"$("+Variable+")$",Interpreter="latex")


    end
end

hold on ;

if options.PlotUnderMesh
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    PlotMuaMesh(CtrlVar,MUA,[],"w") ;
    hold on

end



if options.PlotGroundingLines
    [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],color=options.GroundingLineColor);
end

if options.PlotCalvingFronts
    [xCF,yCF]=PlotCalvingFronts(CtrlVar,MUA,F,color=options.CalvingFrontColor);
end


if options.PlotMuaBoundary
    PlotMuaBoundary(CtrlVar,MUA,"k");
end

% Just guessing that this might be the most common case, the user can easily change afterwards anyhow.

if CtrlVar.PlotXYscale==1000
    xlabel("xps (km)",Interpreter="latex")
    ylabel("yps (km)",Interpreter="latex")
else
    xlabel("x (m)",Interpreter="latex")
    ylabel("y (m)",Interpreter="latex")
end



axis tight

end