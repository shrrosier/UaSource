function Ua2D(UserRunParameters)
    

    if nargin==0
        UserRunParameters=[];
    end
            

    SetUaPath() %% set path


    warning('off','MATLAB:TriRep:PtsNotInTriWarnId')
    
    
    %% initialize some variables
    Info=UaRunInfo; 
    l=UaLagrangeVariables; % these are the Lagrange variables assosiated with boundary conditions.
                           % not to be confused with the Lagrange variables assosiated with solving the 
                           % adjoint problem
    
    RunInfo=[];
    Lubvb=[];
    GLdescriptors=[] ;
    da0dt=[];
    dsdt=NaN; dbdt=NaN; dhdt=NaN; Ruv=[];
    dGFdt=[];  % get rid of this at a later stage
    
    
    %% Define default values
    CtrlVar=Ua2D_DefaultParameters();
    CtrlVar.UserParameters=UserRunParameters;
    
    
    %% Get user-defined values for Ctrl Variable and FE mesh outline
    %  CtrlVar,UsrVar,Info,UaOuts
    [Experiment,CtrlVar,time,dt,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(CtrlVar);
    
    CtrlVar.Experiment=Experiment;
    CtrlVar.time=time;
    CtrlVar.dt=dt;
    CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates;
    clear Experiment time dt MeshBoundaryCoordinates;
    
    %CtrlVar.Experiment=Experiment;
    CtrlVar.CurrentRunStepNumber=0 ; 
    %CtrlVar.MeshBoundaryCoordinates=MeshBoundaryCoordinates; 
    
    if ~isfield(CtrlVar,'fidlog')
        CtrlVar.fidlog=1;
    end
    
    
    CtrlVar=CtrlVarValidityCheck(CtrlVar);
    
    
    CtrlVar.Logfile=[CtrlVar.Experiment,'.log'];
    [pathstr,name,ext]=fileparts(CtrlVar.GmeshFile); CtrlVar.GmeshFile=[pathstr,name]; % get rid of eventual file extension
   
    if CtrlVar.Restart || CtrlVar.AdjointRestart
        CtrlVar.InfoFile = fopen([CtrlVar.Experiment,'-RunInfo.txt'],'a');
    else
        CtrlVar.InfoFile = fopen([CtrlVar.Experiment,'-RunInfo.txt'],'w');
    end

    
    CtrlVar.MeshChanged=0;  % true if mesh changed in last adapt-meshing stage
    
    [status,message]=copyfile(CtrlVar.Logfile,[CtrlVar.Logfile,'~']);  %  copy potential previous logfile
    
  
        
    %% write out some basic information about the type of run selected
    
    PrintRunInfo(CtrlVar);
    
    %%
    %
    % if CtrlVar.Parallel
    %     if matlabpool('size') >0 ; matlabpool close force local ;end
    %     str=computer;
    %
    %     if  strcmp(str,'GLNXA64')
    %         matlabpool open local 8
    %     else
    %         matlabpool open local 2
    %     end
    % end
    
    %%
    tTime=tic;
    
  
    %% Get input data 
    if ~CtrlVar.InverseRun %  forward run
        
        if CtrlVar.Restart  % Forward restart run
            
            [CtrlVar,MUA,BCs,time,dt,s,b,S,B,ub,vb,ud,vd,l,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
                dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
                GF,GLdescriptors,CurrentRunStepNumber]=GetInputsForForwardRestartRun(CtrlVar);
            CtrlVar.time=time;
            CtrlVar.dt=dt;
            CtrlVar.CurrentRunStepNumber=CurrentRunStepNumber;
            clear time dt CurrentRunStepNumber
            
            
        else % New forward run (ie not a restart)
            [CtrlVar,MUA,BCs,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
                dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,...
                GF]=GetInputsForForwardRun(CtrlVar);
            
            
            if CtrlVar.OnlyMeshDomainAndThenStop
                return
            end
        end
        
    else % inverse run
        
        if CtrlVar.AdjointRestart==1  % inverse restarted run

            [CtrlVar.Experiment,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info]=...
                GetInputsForInverseRestartRun(CtrlVar);
            
        else % New inverse run
            
            % first get variables defining the forward run
            [CtrlVar,MUA,BCs,s,b,S,B,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,alpha,as,ab,...
                dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1,GF]=GetInputsForForwardRun(CtrlVar);

            % now get the additional variables specific to an inverse run
            [InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(CtrlVar.Experiment,CtrlVar,MUA,BCs,CtrlVar.time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF);
            
        end
    end
    
    if CtrlVar.TestUserInputs==1
        CtrlVar.TestUserInputs=0;
    end
    
    
    
    %%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%    now all data specific to particular runs should have been defined %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %  For convenience I assume that the user defines S, B, s and b.  The
    %  program then calculates h=s-b and then s and b from h, B and S given the
    %  ice and ocean specific density.  The thickness is preserved, and s and b
    %  are consistent with the floating condition for a given ice tickness h, rho and rhow.
    
    h=s-b;
    [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
    GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
  
    %[nip,niph]=NrOfIntegrationPoints(CtrlVar);
    
    % pointers to the elements of Boundary.Edges where u and v are fixed
    %Boundary.uFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,ufixednode)')));
    %Boundary.vFixedEdgesAllPtrs=logical(prod(double(ismember(Boundary.Edges,vfixednode)')));
    % xint, yint :  matrices of coordinates of integration points. Nele  x nip
    % Xint, Yint :  vectors of unique coordinates of integration points
    %[DTxy,TRIxy,DTint,TRIint,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(MUA);
    
    
    %%
    if CtrlVar.doInverseStep==1;   % -inverse
        
        %         [db,dc] = Deblurr2D(NaN,...
        %             s,u,v,b,B,...
        %             sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,xMeas,yMeas,...
        %             Experiment,...
        %             coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,...
        %             Luv,Luvrhs,lambdauv,n,m,alpha,rho,rhow,g,nStep);
        %
        %x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
        %figure(21) ; trisurf(TRI,x/CtrlVar.PlotXYscale,y/CtrlVar.PlotXYscale,h) ;  title(' h')
       
        [InvFinalValues,ub,vb,ud,vd,l.ubvb,l.udvd,xAdjoint,yAdjoint,Info]=...
            InvertForModelParameters(CtrlVar.Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l.ubvb,l.udvd,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info);
        
        C=InvFinalValues.C          ; fprintf(CtrlVar.fidlog,' C set equal to InvFinalValues.C \n ');
        AGlen=InvFinalValues.AGlen  ; fprintf(CtrlVar.fidlog,' AGlen set equal InvFinalValues.AGlen \n ');
        m=InvFinalValues.m ; n=InvFinalValues.n ; 
        
        
        % this calculation not really needed as AdjointNR2D should return converged ub,vb,ud,vd values for Cest and AGlenEst
        [ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,L]= uv(CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);
        

        if CtrlVar.doplots
            AdjointResultsPlots(CtrlVar.Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,...
                InvStartValues,Priors,Meas,BCsAdjoint,Info,InvFinalValues,xAdjoint,yAdjoint);
        end

        if CtrlVar.AdjointWriteRestartFile==1
            
            fprintf(CtrlVar.fidlog,' %s saving adjoint restart file: %s \n ',CtrlVar.Experiment,CtrlVar.NameOfAdjointRestartFiletoWrite);
            save(CtrlVar.NameOfAdjointRestartFiletoWrite,...
                'CtrlVar','MUA','BCs','s','b','h','S','B','ub','vb','ud','vd','l','alpha','rho','rhow','g','GF',...
                'InvStartValues','Priors','Meas','BCsAdjoint','Info','InvFinalValues','xAdjoint','yAdjoint');
            
            xEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
            yEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
            fprintf(CtrlVar.fidlog,' saving C and m  in file %s \n ',CtrlVar.NameOfFileForSavingSlipperinessEstimate)        ;
            save(CtrlVar.NameOfFileForSavingSlipperinessEstimate,'C','m','xEle','yEle','MUA')
            fprintf(CtrlVar.fidlog,' saving AGlen and m in file %s \n ',CtrlVar.NameOfFileForSavingAGlenEstimate) ;
            save(CtrlVar.NameOfFileForSavingAGlenEstimate,'AGlen','n','xEle','yEle','MUA')
            
        end
        
        
        SayGoodbye(CtrlVar);
        return  % This is the end of the run
        
    end
    %%
    
    CtrlVar.FE2dTransientPlotsCounter=0;
    if ReminderFraction(CtrlVar.time,CtrlVar.TransientPlotDt)<1e-5 || CtrlVar.TransientPlotDt==0
        
        CtrlVar.FE2dTransientPlotsCounter=CtrlVar.FE2dTransientPlotsCounter+1;
        CtrlVar.FE2dTransientPlotsInfostring='First transient plot';
        [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
        [wSurf,wSurfInt,wBedInt,wBed]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,ub,vb,as,ab,exx,eyy,xint,yint,MUA.coordinates,MUA.connectivity,MUA.nip,CtrlVar);
        [DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);
        FE2dTransientPlots(CtrlVar,DTxy,TRIxy,CtrlVar.MeshBoundaryCoordinates,GF,dGFdt,MUA.coordinates,MUA.connectivity,...
            b,B,S,s,h,ub,vb,wSurf,dhdt,dsdt,dbdt,C,AGlen,m,n,xint,yint,wSurfInt,etaInt,exx,eyy,exy,e,CtrlVar.time,...
            rho,rhow,as+ab,as,ab,MUA.Boundary,MUA.nip);
        
        if CtrlVar.OnlyDoFirstTransientPlotAndThenStop
            return
        end
        
    end
    
    % UaOutputs
    CtrlVar.UaOutputsCounter=0;
    if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
        CtrlVar.UaOutputsInfostring='First call ';
        CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
          
        UaOutputs(CtrlVar,MUA,CtrlVar.time,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,GF,BCs,l);
        
        if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls;
            fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
            return
        end
    end
    %

    %nStep0=CurrentRunStepNumber;  
    CtrlVar.CurrentRunStepNumber0=CtrlVar.CurrentRunStepNumber;
 %   if CtrlVar.PlotWaitBar ;     multiWaitbar('CloseAll'); end
    
    %%  time loop
    while 1 
        
        if CtrlVar.CurrentRunStepNumber >=( CtrlVar.nTimeSteps+CtrlVar.CurrentRunStepNumber0)
            fprintf('Exiting time loop because total number of steps reached. \n')
            break
        end

        if (CtrlVar.TotalTime - CtrlVar.time) <= CtrlVar.dtmin
            fprintf('Exiting time loop because total time reached. \n')
            break
        end
        
        if CtrlVar.TimeDependentRun && CtrlVar.dt <= CtrlVar.dtmin % I limit dt some small value for numerical reasons
            fprintf('Exiting time loop because time step too small (%g<%g)\n',CtrlVar.dt,CtrlVar.dtmin)
            TempFile=[CtrlVar.Experiment,'-UaDumpTimeStepTooSmall.mat']; fprintf(CtrlVar.fidlog,' saving variables in %s \n ',TempFile) ; save(TempFile)
            break
        end
        
        %CurrentRunStepNumber=CurrentRunStepNumber+1;  CtrlVar.CurrentRunStepNumber=CurrentRunStepNumber;
        CtrlVar.CurrentRunStepNumber=CtrlVar.CurrentRunStepNumber+1;  
        
        if CtrlVar.PlotWaitBar ;
            multiWaitbar('Time Steps','Value',(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.nTimeSteps);
            multiWaitbar('Model Run Time','Value',CtrlVar.time/CtrlVar.TotalTime);
        end
        
        MUA=UpdateMUA(CtrlVar,MUA);
        
        % -adapt time step
        if CtrlVar.doPrognostic
            dtOld=CtrlVar.time;
            CtrlVar.dt=AdaptiveTimeStepping(CtrlVar,CtrlVar.time,CtrlVar.dt,RunInfo,dubdt,dvbdt,dhdt);
            dtRatio=CtrlVar.dt/dtOld;
           
        else
            CtrlVar.dt=0;
        end
        
        
        if CtrlVar.doDiagnostic  && CtrlVar.InDiagnosticRunsDefineIceGeometryAtEveryRunStep
            % in a diagnostic run I always always use the user-defined geometry
            % for each and every step of the calculation
            [s,b,S,B,alpha]=GetGeometry(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,'sbSB');
            h=s-b;
        elseif CtrlVar.DefineOceanSurfaceAtEachTimeStep
            [~,~,S,~,~]=GetGeometry(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,'S');
        end
        
        [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
        GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);  
        [C,m]=GetSlipperyDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
        [AGlen,n]=GetAGlenDistribution(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
        if CtrlVar.UpdateBoundaryConditionsAtEachTimeStep
            BCs=GetBoundaryConditions(CtrlVar.Experiment,CtrlVar,MUA,BCs,CtrlVar.time,s,b,h,S,B,ub,vb,ud,vd,GF);
        end

        
        %% adaptive meshing, adapt mesh, adapt-mesh
        if CtrlVar.AdaptMesh || CtrlVar.TimeGeometries.Flag ||  CtrlVar.FEmeshAdvanceRetreat
            
            % looks like far too many parameters here, but I need the whole model definition 
            % and on return everything will be defined on a new mesh
            [CtrlVar,MUA,BCs,CtrlVar.MeshBoundaryCoordinates,GF,GLdescriptors,...
                s,b,h,S,B,ub,vb,ud,vd,l,rho,rhow,g,AGlen,n,C,m,ab,as,dhdt,dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1]=...
                AdaptMesh(CtrlVar,CtrlVar.Experiment,CtrlVar.MeshBoundaryCoordinates,MUA,BCs,CtrlVar.time,CtrlVar.CurrentRunStepNumber,...
                GF,GLdescriptors,alpha,...
                s,b,h,S,B,ub,vb,ud,vd,Ruv,Lubvb,l,rho,rhow,g,AGlen,n,C,m,ab,as,dhdt,dhdtm1,dubdt,dvbdt,dubdtm1,dvbdtm1,duddt,dvddt,duddtm1,dvddtm1);
          
            
            if MUA.Nele==0 ;
                fprintf('FE mesh is empty \n ')
                break ;
            end
            
            if CtrlVar.AdaptMeshAndThenStop
                save(CtrlVar.SaveInitialMeshFileName,'MUA') ; 
                fprintf(CtrlVar.fidlog,' MUA was saved in %s .\n',CtrlVar.SaveInitialMeshFileName);
                fprintf('Exiting \n ')
                return
            end
            
        end
        
        %%
        
        ub0=ub ; vb0=vb; ud0=ud ; vd0=vd ; h0=h; s0=s ; b0=b;
        [as0,ab0]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
        a0=as0+ab0;
        
        
        if CtrlVar.doDiagnostic==1 || ( CtrlVar.doPrognostic==1 && CtrlVar.Implicituvh==0)
            %% Diagnostic calculation
            tdiagnostic=tic;                  % -uv
            [ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,Lubvb]= uv(CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);
           
            
            tdiagnostic=toc(tdiagnostic);
            
            dubdtm1=dubdt ; dvbdtm1=dvbdt; duddtm1=duddt ; dvddtm1=dvddt;
            
            if CtrlVar.dt==0 ;
                dubdt=zeros(MUA.Nnodes,1) ;  dvbdt=zeros(MUA.Nnodes,1);
                duddt=zeros(MUA.Nnodes,1) ;  dvddt=zeros(MUA.Nnodes,1);
            else
                dubdt=(ub-ub0)/CtrlVar.dt ; dvbdt=(vb-vb0)/CtrlVar.dt;
                duddt=(ud-ud0)/CtrlVar.dt ; dvddt=(vd-vd0)/CtrlVar.dt;
            end
            
            
        end
        
        
        
        if CtrlVar.doPrognostic==1 && CtrlVar.Implicituvh==1
            
            fprintf(CtrlVar.fidlog,...
                '\n ===== Implicit uvh going from t=%-.15g to t=%-.15g with dt=%-g. Done %-g %% of total time, and  %-g %% of steps \n ',...
                CtrlVar.time,CtrlVar.time+CtrlVar.dt,CtrlVar.dt,100*CtrlVar.time/CtrlVar.TotalTime,100*(CtrlVar.CurrentRunStepNumber-1-CtrlVar.CurrentRunStepNumber0)/CtrlVar.nTimeSteps);
            
            %% possibly start with one diagnostic step
            if CtrlVar.InitialDiagnosticStep==1
                
                RunInfo.converged=0; dIt=0;
                while RunInfo.converged==0
                    dIt=dIt+1;
                    
                    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-.15g \n ',CtrlVar.time);
                    
                    tdiagnostic=tic;    
                    
                    [ub,vb,ud,vd,l,Kuv,Ruv,RunInfo,Lubvb]= uv(CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF);
                                                           
                
                    tdiagnostic=toc(tdiagnostic);
                    if  RunInfo.converged==1 || dIt==2
                        break
                    end
                    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-g did not converge. Resetting (u,v) to zero and trying again \n ',CtrlVar.time);
                    TempFile=[CtrlVar.Experiment,'-UaDumpDiagnosticStepNotConverged.mat']; fprintf(CtrlVar.fidlog,' saving variables in %s \n ',TempFile) ; save(TempFile)
                    ub=ub*0 ; vb=vb*0 ; ud=ud*0 ; vd=vd*0;  l.ubvb=l.ubvb*0; % if did not converge try to reset
                end
                
                if dIt==2 && RunInfo.converged==0
                    fprintf(CtrlVar.fidlog,' initial diagnostic step at t=%-g did not converge despite resetting (u,v) to zero \n ',CtrlVar.time);
                    fprintf(CtrlVar.fidlog,' saving variables in UaDump2 \n ') ; save UaDump2
                end
                
                ub0=ub ; ud0=ud ; vb0=vb ; vd0=vd;
                CtrlVar.InitialDiagnosticStep=0;
            end
            
            if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
                CtrlVar.UaOutputsInfostring='Diagnostic step';
                CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
                UaOutputs(CtrlVar,MUA,CtrlVar.time,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,GF,BCs,l);
                if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls;
                    fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                        CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
                    return
                end
            end
            
            
            
            %% get an explicit estimate for u, v and h at the end of the time step
            
            [ub1,vb1,ud1,vd1,h1]=ExplicitEstimation(CtrlVar.dt,dtRatio,CtrlVar.CurrentRunStepNumber,ub,dubdt,dubdtm1,vb,dvbdt,dvbdtm1,ud,duddt,duddtm1,vd,dvddt,dvddtm1,h,dhdt,dhdtm1);
            
            %% advance the solution by dt using a fully implicit method with respect to u,v and h
            uvhStep=1;
            while uvhStep==1  && CtrlVar.dt > CtrlVar.dtmin % if uvh step does not converge, it is repeated with a smaller dt value
                
                [as,ab]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time+CtrlVar.dt,s,b,h,S,B,rho,rhow,GF);
                as1=as ; ab1=ab; 
                %        0  : values at t
                %        1  : explicit guess for values at t+dt
                % on output : converged values at t+dt
                [ub,vb,ud,vd,h,l.ubvb,l.h,RunInfo,CtrlVar,BCs,CtrlVar.dt]=...
                    FIuvh2D(CtrlVar,MUA,BCs,CtrlVar.dt,S,B,ub0,vb0,ud0,vd0,h0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,...
                    dubdt,dvbdt,duddt,dvddt,...
                    l.ubvb,l.h,...
                    AGlen,C,n,m,alpha,rho,rhow,g);
                uvhStep=0;  % assuming it converged
                
                if RunInfo.converged==0
                    
                    fprintf(CtrlVar.fidlog,' Warning : Reducing time step from %-g to %-g \n',CtrlVar.dt,CtrlVar.dt/10);
                    CtrlVar.dt=CtrlVar.dt/10; 
                    uvhStep=1;  % continue within while loop
                    
                    fprintf(CtrlVar.fidlog,'Also resetting u1, v1, h1 to ub0, vb0 and h0, and setting estimates for Lagrange parameters to zero. \n');
                    [b,s,h0]=Calc_bs_From_hBS(h0,S,B,rho,rhow,CtrlVar,MUA.coordinates);
                    ub1=ub0*0 ; vb1=vb0*0 ;  ud1=ud0*0 ; vd1=vd0*0 ; h1=h0;
                    l.ubvb=l.ubvb*0; l.h=l.h*0;
                end
            end
            
            CtrlVar.time=CtrlVar.time+CtrlVar.dt; 
            [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
            
            dhdtm1=dhdt ; dubdtm1=dubdt ; dvbdtm1=dvbdt;
            if CtrlVar.dt==0 ;
                dhdt=zeros(MUA.Nnodes,1) ;  dubdt=zeros(MUA.Nnodes,1); dvbdt=zeros(MUA.Nnodes,1); dsdt=zeros(MUA.Nnodes,1) ; dbdt=zeros(MUA.Nnodes,1);
            else
                dhdt=(h-h0)/CtrlVar.dt; dubdt=(ub-ub0)/CtrlVar.dt ; dvbdt=(vb-vb0)/CtrlVar.dt; duddt=(ud-ud0)/CtrlVar.dt ; dvddt=(vd-vd0)/CtrlVar.dt; dsdt=(s-s0)/CtrlVar.dt; dbdt=(b-b0)/CtrlVar.dt;
            end
            
        end
        
        
        if CtrlVar.doPrognostic==1 && CtrlVar.Implicituvh==0
            
            if CtrlVar.InfoLevel>0 ; fprintf(CtrlVar.fidlog,'Prognostic, semi-implicit, advancing time from t=%-g to t=%-g \n',CtrlVar.time,CtrlVar.time+dt);end
            
            
            tprognostic=tic;
            
            [ub1,vb1]=ExplicitEstimation(CtrlVar.dt,dtRatio,CtrlVar.CurrentRunStepNumber,ub,dubdt,dubdtm1,vb,dvbdt,dvbdtm1);

            dub1dt=dubdt; dvb1dt=dvbdt ;  dub0dt=dubdt; dvb0dt=dvbdt ; % could possibly be done a bit better
            [as,ab]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time+dt,s,b,h,S,B,rho,rhow,GF);
            a1=as+ab; da0dt=(a1-a0)/CtrlVar.dt ; da1dt=da0dt; 
            if dt==0 ; da0dt=zeros(MUA.Nnodes,1); da1dt=zeros(MUA.Nnodes,1) ; end
            [h,l.h]=SSS2dPrognostic(CtrlVar.dt,h0,ub0,vb0,dub0dt,dvb0dt,a0,da0dt,ub1,vb1,a1,da1dt,dub1dt,dvb1dt,MUA.coordinates,MUA.connectivity,MUA.Boundary,MUA.niph,Lh,Lhrhs,l.h,CtrlVar.CurrentRunStepNumber,CtrlVar);
            
            CtrlVar.time=CtrlVar.time+CtrlVar.dt; %CtrlVar.time=time;
            
            
            [b,s,h]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar,MUA.coordinates);
            
            dhdtm1=dhdt ;
            if CtrlVar.dt==0 ;
                dhdt=zeros(MUA.Nnodes,1) ; dsdt=zeros(MUA.Nnodes,1) ; dbdt=zeros(MUA.Nnodes,1);
            else
                dhdt=(h-h0)/CtrlVar.dt; dsdt=(s-s0)/CtrlVar.dt; dbdt=(b-b0)/CtrlVar.dt;
            end
            
            
            
            %CtrlVar.hChange=1; CtrlVar.s=1; % h and s just changed
            tprognostic=toc(tprognostic);
            if CtrlVar.InfoLevel >= 1 && fprintf(CtrlVar.fidlog,'SSTREAM semi-implicit step in %-g sec, \t prognostic in %-g and diagnostic in %-g sec \n ',tprognostic+tdiagnostic,tprognostic,tdiagnostic) ; end
        end
        
        %% calculations for this step are now done, only some plotting/writing issues do deal with
        
        % calculating derived quantities
      
        GF = GL2d(B,S,h,rhow,rho,MUA.connectivity,CtrlVar);
                
        
        if CtrlVar.doPrognostic==1 && CtrlVar.InfoLevel>=10 
            dhdtMean=(dhdt+dhdtm1)/2;
            [maxdhdt,imaxdhdt]=max(dhdtMean);
            [mindhdt,imindhdt]=min(dhdtMean);
            fprintf(CtrlVar.fidlog,'max(h) %-g \t min(h) %-g \t max(dhdt) %-g \t min(dhdt) %-g \t mean(dhdt) %-g \t median(dhdt) %-g \t rms(h) %-g \t h(max(dhdt)) %-g h(min(dhdt)) %-g\n ',...
                max(h),min(h),max(dhdtMean),min(dhdtMean),mean(dhdt),median(dhdt),sqrt(norm(dhdt)/numel(h)),h(imaxdhdt),h(imindhdt));
            
        end
        %% plotting results
        if ReminderFraction(CtrlVar.time,CtrlVar.TransientPlotDt)<1e-5 || CtrlVar.TransientPlotDt==0 
            
            
            
            CtrlVar.FE2dTransientPlotsCounter=CtrlVar.FE2dTransientPlotsCounter+1;
            CtrlVar.FE2dTransientPlotsInfostring='Time-loop transient plot';
            [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
            [wSurf,wSurfInt,wBedInt,wBed]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,ub,vb,as,ab,exx,eyy,xint,yint,MUA.coordinates,MUA.connectivity,MUA.nip,CtrlVar);
            [DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);
            FE2dTransientPlots(CtrlVar,DTxy,TRIxy,CtrlVar.MeshBoundaryCoordinates,GF,dGFdt,MUA.coordinates,MUA.connectivity,...
                b,B,S,s,h,ub,vb,wSurf,dhdt,dsdt,dbdt,C,AGlen,m,n,xint,yint,wSurfInt,etaInt,exx,eyy,exy,e,CtrlVar.time,...
                rho,rhow,as+ab,as,ab,MUA.Boundary,MUA.nip);
        end
        
           % UaOutputs
    
    if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 ) 
        CtrlVar.UaOutputsInfostring='inside transient loop';
        CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;    
        
        if CtrlVar.MassBalanceGeometryFeedback>0
            [as,ab]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time+dt,s,b,h,S,B,rho,rhow,GF);
        end
        
        UaOutputs(CtrlVar,MUA,CtrlVar.time,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,GF,BCs,l);
        
        if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls;
            fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
            return
        end
    end
    %
        
%%  Better done outside Ua as a part of post-analysis     
%         if CtrlVar.CompareWithAnalyticalSolutions
%             
%             switch Experiment
%                 
%                 case 'TestGaussPeak'
%                     CompareRestultsWithAnalyticalTransferFunctions(Experiment,MUA.coordinates,MUA.connectivity,ub,vb,s,b,S,B,time,CtrlVar.dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar);
%                 case 'Test1dIceStream'
%                     CompareRestultsWith1dIceStreamSolutions(Experiment,MUA.coordinates,MUA.connectivity,ub,vb,s,b,S,B,time,CtrlVar.dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar);
%                 case 'Test1dIceShelf'
%                     CompareRestultsWith1dIceShelfSolutions(Experiment,MUA.coordinates,MUA.connectivity,ub,vb,s,b,S,B,time,CtrlVar.dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar);
%                 otherwise
%                     fprintf(CtrlVar.fidlog,' case not found for analytical comparison \n');
%             end
%             
%         end
%%
% 
%         if CtrlVar.CompareResultsWithPreviouslyObtainedResults
%             CompareResultsWithPreviouslyObtainedResults(Experiment,MUA.coordinates,MUA.connectivity,ub,vb,s,b,S,B,time,CtrlVar.dt,AGlen,C,n,m,rho,rhow,g,alpha,nip,DTxy,TRIxy,DTint,TRIint,CtrlVar);
%         end
%         
        
        if CtrlVar.WriteRestartFile==1 && mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)==0
            WriteRestartFile()
        end
        
        %% write a dump file
        %if CtrlVar.WriteDumpFile && (mod(nStep,CtrlVar.WriteDumpFileStepInterval)==0 || mod(time,CtrlVar.WriteDumpFileTimeInterval)==0 )
        if CtrlVar.WriteDumpFile && (ReminderFraction(CtrlVar.time,CtrlVar.WriteDumpFileTimeInterval)<1e-5 || ReminderFraction(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteDumpFileStepInterval)==0)
            
            
            if ReminderFraction(CtrlVar.time,CtrlVar.WriteDumpFileTimeInterval)<1e-5
                dumpfile=sprintf('%s-DumpFile%-gT-',CtrlVar.Experiment,CtrlVar.time);  dumpfile=regexprep(dumpfile,'\.','k');
            else
                dumpfile=sprintf('%sDumpFile%i',CtrlVar.Experiment,CtrlVar.CurrentRunStepNumber);
            end
            
            fprintf(CtrlVar.fidlog,' Saving everything in file %s  at t=%-g \n',dumpfile,CtrlVar.time);
            
            try
                save(dumpfile','s','h','b','B','S','ub','vb','wSurf','wBed','as','ab','time','MUA','AGlen','C',...
                    'Luv','Luvrhs','lambdauv','Lh','Lhrhs','lambdah','n','m','alpha','rhow','rho','g',...
                    'dubdt','dvbdt','a0','da0dt','dhdt','dhdtm1','dubdtm1','dvbdtm1','DTxy','TRIxy','CtrlVar')
            catch exception
                fprintf(CtrlVar.fidlog,' Could not save file %s \n ',dumpfile);
                fprintf(CtrlVar.fidlog,'%s \n',exception.message);
            end
        end
                
    end
    
    if CtrlVar.PlotWaitBar 
        multiWaitbar('Time Steps','Value',(CtrlVar.CurrentRunStepNumber-CtrlVar.CurrentRunStepNumber0)/CtrlVar.nTimeSteps);
        multiWaitbar('Model Run Time','Value',CtrlVar.time/CtrlVar.TotalTime);
    end
    
        
    %% plotting results

    [etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
    [wSurf,wSurfInt,wBedInt,wBed]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,ub,vb,as,ab,exx,eyy,xint,yint,MUA.coordinates,MUA.connectivity,MUA.nip,CtrlVar);
    
    if ReminderFraction(CtrlVar.time,CtrlVar.TransientPlotDt)<1e-5 || CtrlVar.TransientPlotDt==0 
        
        
        CtrlVar.FE2dTransientPlotsCounter=CtrlVar.FE2dTransientPlotsCounter+1;
        CtrlVar.FE2dTransientPlotsInfostring='Last transient plot';
        [DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);
        FE2dTransientPlots(CtrlVar,DTxy,TRIxy,CtrlVar.MeshBoundaryCoordinates,GF,dGFdt,MUA.coordinates,MUA.connectivity,...
            b,B,S,s,h,ub,vb,wSurf,dhdt,dsdt,dbdt,C,AGlen,m,n,xint,yint,wSurfInt,etaInt,exx,eyy,exy,e,CtrlVar.time,...
            rho,rhow,as+ab,as,ab,MUA.Boundary,MUA.nip);
    end
    
    if (ReminderFraction(CtrlVar.time,CtrlVar.UaOutputsDt)<1e-5 || CtrlVar.UaOutputsDt==0 )
        CtrlVar.UaOutputsInfostring='Last call';
        CtrlVar.UaOutputsCounter=CtrlVar.UaOutputsCounter+1;
        if CtrlVar.MassBalanceGeometryFeedback>0
            [as,ab]=GetMassBalance(CtrlVar.Experiment,CtrlVar,MUA,CtrlVar.time+dt,s,b,h,S,B,rho,rhow,GF);
        end
        
        UaOutputs(CtrlVar,MUA,CtrlVar.time,s,b,S,B,h,ub,vb,ud,vd,dhdt,dsdt,dbdt,C,AGlen,m,n,rho,rhow,g,as,ab,GF,BCs,l);
        
        if CtrlVar.UaOutputsCounter>=CtrlVar.UaOutputsMaxNrOfCalls;
            fprintf(' Exiting because number of calls to UaOutputs (%i) >= CtrlVar.UaOutputsMaxNrOfCalls (%i) /n',...
                CtrlVar.UaOutputsCounter,CtrlVar.UaOutputsMaxNrOfCalls)
            return
        end
    end
    
    if CtrlVar.doplots && CtrlVar.FE2dPlots
        [DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);
        FE2dPlots(CtrlVar,DTxy,TRIxy,CtrlVar.MeshBoundaryCoordinates,GF,dGFdt,MUA.coordinates,MUA.connectivity,...
            b,B,S,s,h,ub,vb,wSurf,dhdt,C,AGlen,m,n,xint,yint,wSurfInt,etaInt,exx,eyy,exy,e,CtrlVar.time,...
            rho,rhow,as+ab,as,ab,txx,tyy,txy);
    end
    
    
    
    
    %% saving outputs
    if CtrlVar.SaveSurfaceData==1
        fprintf(CtrlVar.fidlog,' Saving surface data in file %s \n ',CtrlVar.SurfaceDataFile);
        sMeas=s ; uMeas=ub ; vMeas=vb ; wMeas=wSurf ; bMeas=b ; BMeas=B ; xMeas=coordinates(:,1) ; yMeas=coordinates(:,2);
        save(CtrlVar.SurfaceDataFile,'sMeas','uMeas','vMeas','wMeas','bMeas','BMeas','xMeas','yMeas')
    end
    
    if CtrlVar.WriteRestartFile==1 &&  mod(CtrlVar.CurrentRunStepNumber,CtrlVar.WriteRestartFileInterval)~=0
        
        WriteRestartFile()

    end
        
    if CtrlVar.PlotWaitBar ;     multiWaitbar('CloseAll'); end
    tTime=toc(tTime); fprintf(CtrlVar.fidlog,' Total time : %-g sec \n',tTime) ;
  
    
    if CtrlVar.fidlog~= 1 ; fclose(CtrlVar.fidlog); end
    fclose(CtrlVar.InfoFile);
    
    SayGoodbye(CtrlVar)
    
    %% nested functions
    
    
    function WriteRestartFile
        RestartFile=CtrlVar.NameOfRestartFiletoWrite;
        fprintf(CtrlVar.fidlog,' \n ################## %s %s ################### \n Writing restart file %s  at t=%-g \n %s \n ',CtrlVar.Experiment,datestr(now),RestartFile,CtrlVar.time);
        %[DTxy,TRIxy]=TriangulationNodesIntegrationPoints(MUA);
        CtrlVarInRestartFile=CtrlVar;
        nStep=CtrlVar.CurrentRunStepNumber;  % later get rid of nStep from all restart files
        Itime=CtrlVar.CurrentRunStepNumber;  % later get rid of Itime from all restart files
        time=CtrlVar.time;
        dt=CtrlVar.dt;
        try
            save(RestartFile,'CtrlVarInRestartFile','MUA','BCs','time','dt','s','b','S','B','h','ub','vb','ud','vd','dhdt','dsdt','dbdt','C','AGlen','m','n','rho','rhow','as','ab','GF',...
                'nStep','Itime','a0','da0dt','dhdtm1','dubdt','dvbdt','dubdtm1','dvbdtm1','duddt','dvddt','duddtm1','dvddtm1',...
                'GLdescriptors','l','alpha','g');
            
        catch exception
            fprintf(CtrlVar.fidlog,' Could not save restart file %s \n ',RestartFile);
            fprintf(CtrlVar.fidlog,'%s \n',exception.message);
        end
        
    end

    
    
end
