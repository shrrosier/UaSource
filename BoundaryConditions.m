classdef BoundaryConditions
    %%
    % 
    %
    %
    %%
    
    properties
        
        % basal velocities
        ubFixedNode=[];
        ubFixedValue=[];
        vbFixedNode=[];
        vbFixedValue=[];
        
        ubTiedNodeA=[];
        ubTiedNodeB=[];
        vbTiedNodeA=[];
        vbTiedNodeB=[];
        
        ubvbFixedNormalNode=[];
        ubvbFixedNormalValue=[];
        
        % deformational velocity, only used when solving SIA or hybrid
        udFixedNode=[];
        udFixedValue=[];
        vdFixedNode=[];
        vdFixedValue=[];
        
        udTiedNodeA=[];
        udTiedNodeB=[];
        vdTiedNodeA=[];
        vdTiedNodeB=[];
        
        udvdFixedNormalNode=[];
        udvdFixedNormalValue=[];
        
        % thickness  - only used if solving for thickness, for example in a uvh step
        hFixedNode=[];
        hFixedValue=[];
        hTiedNodeA=[];
        hTiedNodeB=[];

        % positive thickness constraints - these are introduced automatically, do not prescribe direclty
        hPosNode=[];
        hPosValue=[];
         
        
        % rate of thickness change - only used when using calculating dh/dt in combination with the ajoint methods
        %                            for example when using measurements of dh/dt in an inversion 
        dhdtFixedNode=[];
        dhdtFixedValue=[];
        dhdtTiedNodeA=[];
        dhdtTiedNodeB=[];
         
        
        % Boundary conditions for the Level Set Field (LSF)
        
        LSFFixedNode=[];
        LSFFixedValue=[];
        LSFTideNodeA=[];
        LSFTideNodeB=[];
        
    end
 
end