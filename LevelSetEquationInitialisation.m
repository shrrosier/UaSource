function [UserVar,RunInfo,LSF,l,LSFqx,LSFqy,BCs]=LevelSetEquationInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)

nargoutchk(7,7)

LSFqx=[]; LSFqy=[] ;

fprintf("LevelSetEquationInitialisation: Initialising the level set. \n ")

switch lower(CtrlVar.LevelSetInitialisationMethod)
    
    
    
    case {"-geometric-","geometric","-geo-","geo"}
        
        options.subdivide = true;
        options.lineup  = true;
        options.plot  = false;
        Value=0 ;  [Xc,Yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Value,options) ;
        CFPD=sqrt(2*min(MUA.EleAreas))/CtrlVar.LevelSetGeometricInitialisationDistanceFactor;
        options.method="InputPoints";
        options.ResampleCalvingFront=true;
        options.CalvingFrontPointDistance=CFPD;
        options.GetRidOfCalvingFrontOutsideComputationalDomain = false;
        options.CalvingFrontPointDistance = 1e3; % this is the (default) distance between the points defining the calving front
        [~,~,LSF]=CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,F0.LSF,options);
        
        
        %%
        
    otherwise
        
        % CtrlVar.LevelSetReinitializePDist=false ;
        [UserVar,RunInfo,LSF,l,LSFqx,LSFqy]=LevelSetEquationAnalyticalInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l);
        F0.LSF=LSF ; F1.LSF=LSF ;
end






end





