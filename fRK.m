function [R,K]=fRK(x,problemtype)

x=x(:) ;

% nargout

switch problemtype

    case "[x1,x2]"

        R=x  ;
        K=eye(numel(x)) ;


    case "[x1^2,x2^2]"

        R=x.^2;
        K=2*[x(1)  0     ;
            0   x(2) ] ;


    case "[x1+x2,x2]"

        R(1)=x(1)+x(2);
        R(2)=x(2);
        R=R(:) ;
        K=[1  1 ; 0 1] ;

    case "[x1^2,x2]"

        R(1)=x(1)^2;
        R(2)=x(2);
        R=R(:) ;
        K=[2*x(1)  0 ; 0 1] ;


    case "[x1^2+x2,x2]"

        R(1)=x(1)^2+x(2);  
        R(2)=x(2);
        R=R(:) ;

        K=[2*x(1)  1 ; ...   % \nabla R1^T
            0     1] ;       % \nabla R2^T

    case "[x1^2+x2,x2^2+x1]"

        R(1)=x(1)^2+x(2);  
        R(2)=x(2)^2+x(1);
        R=R(:) ;

        K=[2*x(1)   1   ; ...   % \nabla R1^T
            1   2*x(2)] ;       % \nabla R2^T

   case "[x1^3-100 x2,-x2^2+10 x1]"

        R(1)=x(1)^3-100*x(2);  
        R(2)=-x(2)^2+10*x(1);
        R=R(:) ;

        K=[3*x(1)^2   -100   ; ...   % \nabla R1^T
            +10   -2*x(2)] ;         % \nabla R2^T

    case "Rosenbrock"

        [f,g,H] = RosenbrockFunction(x) ; 

         R=g ; K = H ; 
    
    otherwise

        error("case not found")


end