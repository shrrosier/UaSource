function PlotRunInfo(RunInfo)
    
    %%
    
    
    FindOrCreateFigure("RunInfo: time step and iterations")
    yyaxis left
    semilogy(RunInfo.Forward.time,RunInfo.Forward.dt,'o-') ; 
    ylabel('time step')
    

    
    yyaxis right 
    stairs(RunInfo.Forward.time,RunInfo.Forward.uvhIterations) ; 
    ylabel('uvh iterations')

    
    xlabel('time') ; 
    legend("time step","#uvh iterations")
    ylim([0 inf])
    
     FindOrCreateFigure("RunInfo: time step histogram and iterations")
     histogram(RunInfo.Forward.dt) ; xlabel('dt')
     title('dt Histogram')
    
     
     
    
end