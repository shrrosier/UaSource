function cint=DefineCalvingAtIntegrationPoints(UserVar,CtrlVar,cint,nx,ny,uint,vint)


if  contains(CtrlVar.CalvingLaw.Evaluation,"-int-")

    factor=1;
    cint=factor*(uint.*nx+vint.*ny);
    
end

end