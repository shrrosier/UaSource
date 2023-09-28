
function [RunInfo,varargout]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin)

nVar=length(varargin) ;
varargout=cell(nVar,1);



switch CtrlVar.MapOldToNew.method
    
    
    case "scatteredInterpolant"
        
        
        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingScatteredInterpolant(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:}) ;
        
    case "FE form functions"
        
        x2 = MUAnew.coordinates(:,1);
        y2 = MUAnew.coordinates(:,2);
        [varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingFEShapeFunctions(CtrlVar,MUAold,x2,y2,varargin{:});

    case "ShapeAndScattered"  

        [RunInfo,varargout{:}]=MapNodalVariablesFromMesh1ToMesh2UsingShapeAndScattered(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,varargin{:});
        

    otherwise

        error(" case not found")
end




end
