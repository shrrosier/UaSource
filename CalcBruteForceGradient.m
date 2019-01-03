function dJ = CalcBruteForceGradient(func,p0,CtrlVar,iRange,deltaStep)




fprintf(' Calculating gradients using brute-force method. \n')


J0=func(p0);

% Testing gradient using brute force method

% deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize*norm(p0);

dJ=p0*0+NaN;
dJtemp=dJ;

switch  lower(CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType)
    
    case {'forward','first-order'}
        
        
        parfor k=1:numel(iRange)
            I=iRange(k);
            p1=p0;
            p1(I)=p1(I)+deltaStep;
            J1=func(p1);
            dJtemp(k)=(J1-J0)/deltaStep;
            
        end
        
        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end
        
    case {'central','second-order'}
        
        
        parfor k=1:numel(iRange)
            I=iRange(k);
            p1=p0;
            pm1=p0;
            p1(I)=p1(I)+deltaStep;
            J1=func(p1);
            
            pm1(I)=pm1(I)-deltaStep;
            Jm1=func(pm1);
            
            dJtemp(k)=(J1-Jm1)/deltaStep/2;
        end
        
        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end
        
        
    case 'fourth-order'
        
        parfor k=1:numel(iRange)
            I=iRange(k);
            
            p1=p0;
            pm1=p0;
            p2=p0;
            pm2=p0;
            p1(I)=p1(I)+deltaStep;
            p2(I)=p2(I)+2*deltaStep;
            pm1(I)=pm1(I)-deltaStep;
            pm2(I)=pm2(I)-2*deltaStep;
            
            J1=func(p1);
            J2=func(p2);
            Jm1=func(pm1);
            Jm2=func(pm2);
            
            dJtemp(k)=(Jm2/12-2*Jm1/3+2*J1/3-J2/12)/deltaStep;
        end
        
        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end
        
        
    case 'complex step differentiation'
        
        fprintf('Error: complex step differentiation not possible because now using a linear solver\n')
        fprintf('that does not work with complex matrices.\n')
        error('CalcBruteForceGradient:CSD','complex step differentiation not possible')
        
        
        for k=1:numel(iRange)
            
            I=iRange(k);
            p1=p0;
            p1(I)=p1(I)+1i*deltaStep;
            J1=func(p1);
            dJtemp(k)=imag(J1)/deltaStep;
            
        end
        
        for k=1:numel(iRange)
            dJ(iRange(k))=dJtemp(k);
        end
        
    otherwise
        
        fprintf(' CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType has invalid value. \n')
        error('which case')
end

% %%
% calculates d kv/ d b using the complex number method
% using the fact that df/dx=Im(f(x+i dx))/dx
% delta=1e-6*norm(p0);
% dJ=p0*0;
%
% for I=iRange
%
%     p1=p0;
%     p1(I)=p1(I)+1i*delta;
%     J1=func(p1);
%     dJ(I)=imag(J1)/delta;
%
% end

%%
%DkvDb=imag(kv)/db;
%DrhDb=imag(rh)/db;

end


