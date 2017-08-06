function [UserVar,RunInfo,r,ruv,rh,rl]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0)



narginchk(15,15)

F1.ub=F1.ub+gamma*dub;
F1.vb=F1.vb+gamma*dvb;
F1.h=F1.h+gamma*dh;
luvh=luvh+gamma*dl;

[UserVar,RunInfo,R]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);




if ~isempty(L)
    
    frhs=-R-L'*luvh;
    grhs=cuvh-L*[F1.ub;F1.vb;F1.h];
    
    %frhs=-R-L'*(luv+gamma*dl);
    %grhs=cuvh-L*[ub+gamma*dub;vb+gamma*dvb;h+gamma*dh];
else
    frhs=-R;
    grhs=[];
end

[r,rl,ruv,rh]=ResidualCostFunction(frhs,grhs,Fext0,MUA.Nnodes);



if isequal(CtrlVar.Residual.uvh,'uv')
    r=ruv; % debug August 2017
end

end

