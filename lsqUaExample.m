
%  R=(x1,x2)


problemtype="[x1,x2]" ;
problemtype="[x1+x2,x2]";
problemtype="[x1^2+x2,x2]";
problemtype="[x1^2,x2]";
% problemtype="[x1^2+x2,x2^2+x1]";
% problemtype="[x1^2-x2,x2^2-x1]";

isConstraint=true; 

CtrlVar.lsqUa.ItMax=10 ;
CtrlVar.lsqUa.tol=1e-30;
CtrlVar.lsqUa.isLSQ=false ;
CtrlVar.lsqUa.LevenbergMarquardt="auto" ; % "fixed"
CtrlVar.lsqUa.LMlambda0=0.01 ;
CtrlVar.lsqUa.Normalize=false;
CtrlVar.lsqUa.SaveIterate=true;


ProgressOverplot=true;



fun = @(x) fRK(x,problemtype)  ;
x0=[-5 ; -5] ; 

lambda= []  ;

if isConstraint
    L=[1 1 ];  c= 5  ;
else
    L=[]; c=[];
end

[xSol,lambda,resnorm,residual,g,h,output] = lsqUa(CtrlVar,fun,x0,lambda,L,c) ;




x1Vector=linspace(-10,10);
x2Vector=linspace(-10,10);
r2=nan(numel(x1Vector),numel(x2Vector));

for I=1:numel(x1Vector)
    for J=1:numel(x2Vector)
        x=[x1Vector(I) x2Vector(J)];
        R=fRK(x,problemtype);
        r2(I,J)=R'*R ;
    end
end


flsqUa=FindOrCreateFigure("lsqUa test") ; clf(flsqUa) ;
contourf(x1Vector,x2Vector,r2',10) ; axis equal tight; colorbar ; axis([-10 10 -10 10])
hold on  ;
if isConstraint
    plot(x1Vector,c-x1Vector,'r')
end

plot(xSol(1),xSol(2),'o',MarkerFaceColor='r',MarkerEdgeColor="w",MarkerSize=12)

for I=1:output.nIt

    % plot(output.xVector(1,I),output.xVector(2,I),'+r')
    text(output.xVector(1,I),output.xVector(2,I),num2str(I-1),color="r")
    
   % txt = input("RET to continue\n") ; 
    
end


flsqUaProg=FindOrCreateFigure("lsqUa progress") ; 

if ProgressOverplot
    hold on
end

npoints=numel(output.rVector);
itVector=0:npoints-1; 
yyaxis left
semilogy(itVector, output.gVector,'o-')
ylabel("$\|g\|^2$",Interpreter="latex")
yyaxis right
semilogy(itVector, output.rVector,'o-')
ylabel("$\|R\|^2$",Interpreter="latex")
xlabel("iteration",Interpreter="latex")


