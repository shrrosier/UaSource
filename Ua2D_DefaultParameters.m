function CtrlVar=Ua2D_DefaultParameters


%%
% CtrlVar=Ua2D_DefaultParameters
%
%  sets the fields of the CtrlVar to their default values
%
%



%%
%  
%  Most likely when running �a, only a fairly limited number of the parameters listed below need to be set/changed. 
%  Changing the parameter values from their default values should be done by the user in `Ua2D_InitialUserInput.m'. 
%  That user m-file should be located in a separate run-directory, together with all the other user m-files


%%

CtrlVar.Experiment='UaDefaultRun';
CtrlVar.time=NaN;             % In a transient run this variable is the (model) time. Then set this to some 
                              % reasonable initial value, for example CtralVar.time=0;
%% Types of run
% 
CtrlVar.TimeDependentRun=0 ;  % either [0|1].  
                              % If true (i.e. set to 1) then the run is a forward transient one, if not
                              % then velocities based on the current geometry are calculated. 
CtrlVar.InverseRun=0;         % if true then a surface-to-bed inversion is to be performed.
                              % (in an inverse run the value of CtrlVar.TimeDependentRun is irrelevant)
                              
CtrlVar.Restart=0;            % If true then the run is a restart run. Note that this variable applies to both forward and inverse runs.
                              % For example setting: 
                              %       CtrlVar.InverseRun=1; 
                              %       CtrlVar.Restart=1;
                              % will give a restarted inverse run. (make sure a corresponding restart file does exist, see below.)
                              %

                              
CtrlVar.TotalNumberOfForwardRunSteps=1;   % maximum number of forward run steps.  In a transient run this will be the maximum number of time steps.
                                          % In a non-transient (stationary) run, this will be the maximum number of diagnostic calculations.
                                          % (Typically, the number of forward run steps in a non-transient run will be 1, and the user must make sure to set 
                                          % the value accordingly, i.e.  CtrlVar.TotalNumberOfForwardRunSteps=1;)
                                          % In a restart run, TotalNumberOfForwardRunSteps is the total number of run steps done within that restart run, i.e.
                                          % not the total accumulated number of forward run steps.
                              
%% Ice flow approximation
CtrlVar.FlowApproximation='SSTREAM' ;  % any of ['SSTREAM'|'SSHEET'|'Hybrid']  
                                       % Note, both SSTREAM and SSHEET are implemented.
                                       % But Hybrid is still in development and should not be used for the time being.

%% Boundary conditions
CtrlVar.UpdateBoundaryConditionsAtEachTimeStep=0;  % if true, `DefineBCs' is called at the beginning of each time step and boundary conditions are updated
                                                   % otherwise boundary conditions are only updated at the beginning of the run (also at the beginning or a restart run).
CtrlVar.BCsWeights=1;  % testing parameter, do not change

%
%% Manually updating geometry in the course of a run.
% By default DefineGeometry is only called at the beginning of a run, and after
% a any mesh modifications.
%
% However, it is possible to force additional `manual' updates to geometry
% during a diagnostic run. This is often usefull, for example, for doing
% sensitivity tests with respect to geometry.
%
% Also, one can redefine the ocean surface elevation at each run step in both
% transient and non-transiten (diagnostic) runs.
%
CtrlVar.DefineOceanSurfaceAtEachTimeStep=0;   % if true,  `DefineGeometry.m' is called at each time step, returning S, and only S.
                                              % if false, `DefineGeometry.m' is only called at the beginning of a run
                                              %            and when the FE-mesh changes
CtrlVar.InDiagnosticRunsDefineIceGeometryAtEveryRunStep=1;  % If true, the user-defined geometry (sbSB) is always used at each run step in a 
                                                            % diagnostic calculation.
                                                        
%%
CtrlVar.TestUserInputs=1;  % By default user inputs will be tested at the start of the run
                           % to suppress set TestUserInputs=0
                           % if user inputs are always to be tested throughout the run, set TestUserInputs=2 
CtrlVar.TestForRealValues=1;
%% Element type
%
% The options are: linear, quadratic, or cubic Lagrangian triangle elements
CtrlVar.TriNodes=6 ;  % Possible values are 3, 6, 10 node (linear/quadradic/cubic)

%% Control on transient runs
% Once either the number of time steps or total time modelled reaches prescribed values
% the run stops.

CtrlVar.TotalTime=1e10;          % maximum model time
CtrlVar.dt=1;                    % time step (usually overwritten by user by defining dt in the Ua2D_InitialUserInputFile
CtrlVar.dtmin=1e-12;             % for numerical reasons the time step should always be larger than some very small value

CtrlVar.InitialDiagnosticStep=0; % Start a transient run with an initial diagnostic step, even if the step is a restart step.
                                 % Irrespective of the value of this variable, an initial diagnostic step is always performed at the beginning of a transient run if it is not a restart run.
                                 % An initial diagnostic step is therefore done at the beginning of a transient run if:
                                 % 1) so asked by the user, i.e. if the user sets CtrlVar.InitialDiagnosticStep=1, and
                                 % 2) at the start of an implicit uvh transient run.
                                 % Unless asked by the user, no initial diagnostic step is done at the beginning of a transient restart run.

CtrlVar.InitialDiagnosticStepAfterRemeshing=1 ; % after each remeshing, do an initial diagnostic step before continuing with further prognostic steps. (Always a good idea.) 

%% Restart option
CtrlVar.Restart=0;                       % either 0/false or 1/true.  Set to 1 for a restart run. (This also work for inverse runs. See below.)
CtrlVar.WriteRestartFile=1;              % if true, a restart file is written
CtrlVar.WriteRestartFileInterval=100;    % restart file written at this time-step interval  (note, these are run steps, not model time)
CtrlVar.ResetTime=0 ;                    % set to 1 to reset (model) time at start of restart run
CtrlVar.RestartTime=NaN;                 % if ResetTime is true, then this is the model time at the start of the restart run
CtrlVar.ResetTimeStep=0;                 % 1 if time step should be reset to dt given in the Ua2D_InitialUserInputFile
CtrlVar.NameOfRestartFiletoRead='Ua2D_Restartfile.mat';
CtrlVar.NameOfRestartFiletoWrite='Ua2D_Restartfile.mat';

%%
CtrlVar.SaveAdaptMeshFileName=[];          % file name for saving adapt mesh. If left empty, no file is written

%% Plotting
%
% Most plotting is typically done by the user using his own version of the `UaOutputs.m',
% or in a separate post-processing step
% However, some basic plots can be generated directly from within Ua.
%

CtrlVar.doplots=1;          % if true then plotting during runs by Ua are allowed, set to 0 to suppress all plots
CtrlVar.PlotWaitBar=1;      % a waitbar is plotted
CtrlVar.doAdaptMeshPlots=1; % if true and if CtrlVar.doplots true also, then do some extra plotting related to adapt meshing
CtrlVar.PlotOceanLakeNodes=0;        % Shows which nodes are considered a part of the `ocean' and which are within `lakes' that have no connection the ocean
CtrlVar.PlotMeltNodes=0;
CtrlVar.PlotXYscale=1;     % used to scale x and y axis of some of the figures, only used for plotting purposes
                           % (if spatial units are in meters, setting this to 1000 produces xy axis with the units km)
CtrlVar.PlotsXaxisLabel='x' ; CtrlVar.PlotsYaxisLabel='y' ; %
CtrlVar.MinSpeedWhenPlottingVelArrows=0;    % when plotting vel arrows with smaller speed are scaled so that their speed its
                                            % equal to this value  (setting this to a large value makes all arrows
                                            % equally long)

CtrlVar.BoundaryConditionsFixedNodeArrowScale=1;  % Determines the size of arrows indicating boundary conditions when plotting boundary conditions. 
                                                  % The arrows are automatically scales with respect to mesh size, but if they are
                                                  % too small or too large this parameter can be used to affect their size. 

%% Plotting mesh
% The mesh can be plotted within Ua by setting CtrlVar.PlotMesh=1, or by calling 
% either PlotFEmesh or PlotMuaMesh (see help PlotFEmesh)
CtrlVar.PlotMesh=0;        % If true then FE mesh is shown every time a new mesh is generated
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=1; 
CtrlVar.FEmeshPlotTitle=[]; % Title for FE mesh plot, if left empty then something sensible is used instead
CtrlVar.PlotFEmeshAndSaveMesh=0 ; % when plotting mesh also save mesh to a file
CtrlVar.PlotBCs=0;         % If true then boundary conditions are shown at the beginning of the run
CtrlVar.PlotNodes=0;       % If true then nodes are plotted when FE mesh is shown
CtrlVar.PlotLabels=0 ;     % If true elements and nodes are labelled with their respective numbers
CtrlVar.LabelNodes=0;      % Nodal labels are plotted
CtrlVar.LabelElements=0;   % Element labels are plotted
CtrlVar.PlotNodesSymbol='o';
CtrlVar.PlotNodesSymbolSize=3;
CtrlVar.MeshColor='k'; CtrlVar.NodeColor='k';


%% Numerical variables related to transient runs
% In general there should be no need to ever change these values except for testing purposes
%
% Transient runs can be done either (fully) implicitly, or semi-implicitly
% In a (fully) implicit approach, the time-integration is done implicitly with respect to both velocities and thickness.
% In a semi-implict approach, the time-integration is done implicitly with respect to thickness, and explicitly with respect to velocities.
%
% There are currently two fully-implicit time-stepping methods implemented: The 'theta' and the 'supg' methods.
%
% The 'theta' method uses a weighted sum of the values at the beginning and the end of a time step.
% The weighting is controlled by CtrlVar.theta and depending on the value of theta different types of
% approximations are obtained: 0,1/2,1 gives forward Euler, Lax-Wendroff and backwards Euler, respectively.
% The 'supg' method is a Streamline-Upwind Petrov-Galerkin method. The supg-method uses the same
% weighting as the 'theta' method, but the test function for the mass-conservation equation is different.
%
% The default time-stepping method is: Fully implicit Streamline-Upwind Petrov-Galerkin with theta=0.5 (Lax Wendroff).
%
%
CtrlVar.Implicituvh=1;           % 0: prognostic run is semi-implicit (implicit with respect to h only)
                                 % 1: prognostic run is fully-implicit (implicit with respect to uvh)

CtrlVar.uvhTimeSteppingMethod='supg'; % 'theta'|'supg'

CtrlVar.SUPG.beta0=0.5 ; CtrlVar.SUPG.beta1=0 ; % parameters related to the SUPG method.
CtrlVar.theta=0.5;    % theta=0 is forward Euler, theta=1 is backward Euler, theta=1/2 is Lax-Wendroff and is most accurate

% Note: An additional time-stepping method is the Third-Order Taylor-Galerkin (TG3) method.
% It has not been fully tested but seems to work very well for fully implicit transient calculation.
% This option that can be obtained by setting:
% CtrlVar.TG3=1 ;  CtrlVar.Test1=1;  CtrlVar.Test0=0;   CtrlVar.theta=0.5;  
% and using the fully-implicit time-stepping option (CtrlVar.Implicituvh=1)); 
CtrlVar.TG3=0 ; % if true, the prognostic steps uses a third-order Taylor-Galerkin method
                % currently only implemented for periodic boundary conditions                         
                % Note, only theta=0.5 is strictly consistent with TG3=1, so
                % for backward Euler set theta=1 and TG3=0                 
CtrlVar.IncludeTG3uvhBoundaryTerm=0;                     % keep zero (only used for testing)
CtrlVar.IncludeDirichletBoundaryIntegralDiagnostic=0;    % keep zero (only used for testing)
  



%% Numerical Regularisation Parameters  (note: these are not related to inverse modelling regularisation)
CtrlVar.SpeedZero=1e-4;     % needs to be larger than 0 but should also be much smaller than any velocities of interest.
CtrlVar.EpsZero=1e-10;      % needs to be larger than 0 but should also be much smaller than any effective strain rates of interest.
CtrlVar.Czero=1e-10;        % 
CtrlVar.CAdjointZero=CtrlVar.Czero; % used as a regularisation parameter when calculating dIdCq.
CtrlVar.dbdxZero=1;   % when calculating basal shear stresses in the hybrid approximation, a very large bed slope causes errors.
CtrlVar.dbdyZero=1;   % a crude solution is to limit bed slopes to 45 degrees. 
CtrlVar.AGlenAdjointZero=100*eps; 
CtrlVar.AdjointEpsZero=CtrlVar.EpsZero;
%% Constraints on viscosity and slipperiness
% These constraints are always enforced, but only really of any importance when inverting for A and/or C.
% (Using SIA or the hybrid approximation Cmin MUST be set to 0, or at least to a value much less than Czero!)
%
switch lower(CtrlVar.FlowApproximation)
    case 'sstream'
        CtrlVar.Cmin=1e-6;          % a reasonable lower estimate of C is u=C tau^m with min u=1 m/a and max tau=100 dPa => C=u/tau^m=1e-6
    otherwise
        CtrlVar.Cmin=0;          % a reasonable lower estimate of C is u=C tau^m with min u=1 m/a and max tau=100 dPa => C=u/tau^m=1e-6
end
CtrlVar.Cmax=1e10;
CtrlVar.AGlenmin=100*eps;
CtrlVar.AGlenmax=1e10;

%% Non-linear iteration-loop parameters
% The non-linear system is considered solved once the residuals are smaller than NLtol,
% and the normalised chances in u,h and \lambda smaller than du, dh and dl.
%
% The most (arguably even the only) important number is NLtol.
% NLtol is a tolerance on the norm of the solution residuals, i.e. the resulting residuals once the solution is
% plugged back into the equation. So NLtol should ALWAYS be set to a small value (for example <1e-10)
%
% The CtrlVar.du/dh/dl are tolerances for the chance in u,h, and \lambda, respectively, between subsequent non-linear iteration steps.
% Although one would expect these to go to zero with increasing iteration number, these are not very reliable
% estimates of the actual error.  Generally set du and dh to not too large value, but do not focus too much on those numbers
% (The error in solving the boundary conditions is always checked internally.)
CtrlVar.NLtol=1e-15; % tolerance for the square of the norm of the residual error
CtrlVar.du=1e-2;     % tolerance for change in (normalised) speed
CtrlVar.dh=1e-2;     % tolerance for change in (normalised) thickness
CtrlVar.dl=100;      % tolerance for change in (normalised) lambda variables used to enforced BCs
%    Note: there is no need to put any constrains on the Lagrange variables
%    used to enforce the BCs because 1) the BCs are currently always linear,
%    and 2) it is always checked internally that the BCs have been solved correctly.
%    In fact, it can be a bad idea to enforce a limit on this change because
%    sometimes the change in lambda between non-linear iteration steps is just a
%    direct response to how the primary variables (u,v,h) change.  The norm
%    of these changes can then be large despite the BCs being exactly fulfilled.)

CtrlVar.NR=1;             % 1 gives Newton-Raphson (use Newton-Raphson whenever possible)
% Modified Newton-Raphson only evaluates the left-hand side (the stiffness matrix) if certain
% criteria are fullfilled. This will reduced time spend with matrix assembly but
% also reduced the rate of convergence. Depending on the problem using the
% modified NR method may, or may not, lead to an overall reduction in
% computational time.
%
% There are two criteria that determine if the left-hand side is updated or not:
% 1) interval and 2) (residual) reduction criteria. The interval criteria
% determines the number of iterations between updates. (The matris is always
% updated at the beginning of the non-linear iteration.) The reduction criteria
% forces re-assembly if the reduciton in last iteration was not greater than
% a given fraction.
%

CtrlVar.ModifiedNRuvIntervalCriterion=1;  % interval between matrix updates, always a positive integer number.
CtrlVar.ModifiedNRuvReductionCriterion=1; % fractional reduction forcing an update
CtrlVar.ModifiedNRuvhIntervalCriterion=1;  
CtrlVar.ModifiedNRuvhReductionCriterion=1;
% Settingn for example:
% CtrlVar.ModifiedNRuvIntervalCriterion=10;     
% CtrlVar.ModifiedNRuvReductionCriterion=0.95;
% will cause the matrix only to be updated every 10-th non-linear iteration, unless
% the fractional reduction r/r0 over previous iteration was less than 0.95.
%

CtrlVar.Piccard=0;        % 1 gives Piccard iteration, otherwise NR iteration (always use NR whenever possible).
CtrlVar.NRviscosity=1;    % if 1 derivatives with respect to viscosity are included in the NR method
CtrlVar.NRbeta2=1;        % if 1 derivatives with respect to slipperiness are included in the NR method
                          % Note: if Piccard=0 then the NRviscosity and NRbeta2 values are overwritten and set to 0. 
CtrlVar.NRitmax=50;       % maximum number of NR iteration
CtrlVar.Piccarditmax=30;  % maximum number of Piccard iterations
CtrlVar.iarmmax=10;       % maximum number of backtracking steps in NR and Piccard iteration
CtrlVar.NRitmin=1;        % minimum number of NR iteration
CtrlVar.NewtonAcceptRatio=0.5;  % accepted reduction in NR without going into back-stepping
CtrlVar.NewtonBacktrackingBeta=1e-4;  %  affects the Amarijo exit criteria in the back-stepping
CtrlVar.LineSeachAllowedToUseExtrapolation=1; % If true, backtracking algorithm may start with an extrapolation step.
CtrlVar.BacktrackingGammaMin=1e-10;  % smallest step-size in Newton/Piccard backtracking as a fraction of the full Newton/Picard step.
CtrlVar.BacktrackingGammaMinAdjoint=1e-20; % smallest step-size allowed while backtracking in adjoint step. (This is an absolut step size, i.e. not a fraction of initial step size.)



%% Backtracking parameters  -line search 
% Parameters affecting the backtracking algorithm
CtrlVar.BackTrackBeta=0.1 ;               % beta in the Armijo�Goldstein exit condition
CtrlVar.BackTrackMaxIterations=50 ;       % this is plenty
CtrlVar.BackTrackMaxExtrapolations=50  ;  % if set to zero no extrapolation is done (i.e. pure backtracking)
CtrlVar.BackTrackExtrapolationRatio=2.5 ; % ratio between new and old step size in each extrapolation step
CtrlVar.BackTrackMinXfrac=1e-10 ;         % exit backtracking if pos. of minimum is changing by less than this fraction of initial step 
CtrlVar.BackTrackMaxFuncSame=3 ;          % exit backtracking if this many evaluations of cost function resulted in no further decrease of cost function
    


%% Lin equation solver parameters


% Linear symmetrical solver is either Matlab \ operator, or Uzawa (outer) iteration
% CtrlVar.SymmSolver can be one of {'Backslash','Uzawa','AugmentedLagrangian'}
CtrlVar.SymmSolver='AugmentedLagrangian';  %

% Linear asymmetrical solver is either Matlab \ operator or Augmented Lagrangian Solver (ALS)
% CtrlVar.AsymmSolver='Backslash';
CtrlVar.AsymmSolver='AugmentedLagrangian';  %
% For asymmetrical indefinite block-structured systems
% the ALS method is almost always better than the default Matlab backslash operator. 
% ALS is an iterative method with an inner and outer iteration. Convergence depends on
% ALSpower. If ALS does not converge then tying a smaller ALSpower
% usually does the trick.
CtrlVar.ALSIterationMin=3;     CtrlVar.ALSIterationMax=25;   CtrlVar.ALSpower=5;  % ALS parameters
CtrlVar.UzawaIterationMin=3;   CtrlVar.UzawaIterationMax=25; CtrlVar.UzawaPower=5;  % Uzawa parameters


CtrlVar.LinSolveTol=1e-10;  % Residual when solving linear system.
                            % If the standard Matlab backslash algorithm is used, default Matlab values apply and this number is not used
                            % For indefinite block-structured systems of the type [A B' ; B 0] [x;y]=[f;g]
                            % the relative residual is defined in standard way as: 
                            % Residual=norm([A B' ; B sparse(m,m)]*[x;y]-[f ; g])/norm([f;g]);   
                            % A value of 1e-10 is arguably an overly small number, in many cases 1e-6 would be considered acceptable

%% Internal variables related to matrix assembly
% These variables are only for testing purposes. Do not change from default
% values.
CtrlVar.CalvingFrontFullyFloating=0;  % if true then the natural BC is only covers a freely floating calving front (do not change, only for testing)
CtrlVar.GroupRepresentation=0;
%% Number of integration points
% if left empty, the number of integration points is set automatically

CtrlVar.niph=[] ;  % number of integration points for uvh in implicit runs, and for the h-solver in semi-implicit runs
CtrlVar.nip=[] ;   % number of integration points for the uv solver
                   % Possible Nr of integration points: 1,3,4,6,7,12,16
                   % The default values are: 
                   % nip=3 and niph=3 for linear elemetns (three node elements)
                   % nip=7 and niph=7 for quadric elements (six node elements)
                   % nip=12 and niph=12 for cubic elements (ten node elements)
                   % The defaul values are usually fine, but sometimes increasing the number of
                   % intergration points improves convergence of the Newton-Raphson iteration.
%% Level of information given during a run
% A number of variables affect the information given during a run.
% Generally the higher the number, the more information is given.
% 
% Depending on info levels, figures might be plotted as well. However, this is only done
% if corresponding plotting logicals such as CtrlVar.doplots, CtrlVar.doAdaptMeshPlot, etc, are also true.
%
CtrlVar.InfoLevel=1;        % Overall level of information (forward runs)  

CtrlVar.InfoLevelAdjoint=1; % Overall level of information (inverse runs). Note: generally good to combine with CtrlVar.InfoLevelNonLinIt=0; CtrlVar.InfoLevel=0;

CtrlVar.InfoLevelNonLinIt=1; % Info level for non-line solve. Generally:
%   0   : no information on non-linear step printed.
%   1  : prints basic convergence information at end of non-linear step.
%  >1  : detailed info on residuals given at the end of non-linear step.
% >=2  : info on backtracking step as well.
% >=10 : calculates/plots additional info on residuals as a function of step size within line search, and rate of convergence
% >=100 : plots residual vectors

CtrlVar.InfoLevelAdaptiveMeshing=1;  

CtrlVar.InfoLevelLinSolve=0;  % If the linear solver does not converge (it sometimes uses a inner and outer loop to deal with indefinite systems)
                              % then increasing this number will give further information. G

CtrlVar.ThicknessConstraintsInfoLevel=1 ;
                              
CtrlVar.Report_if_b_less_than_B=0; % 

CtrlVar.SymmSolverInfoLevel=0 ;
CtrlVar.InfoLevelBackTrack=1;

CtrlVar.InfoLevelCPU=1;  % if 1 then some info on CPU time usage is given

CtrlVar.StandartOutToLogfile=false ; % if true standard output is directed to a logfile
% name of logfile is  $Experiment.log



%% Inversion 
% 
% Inversion can currently be done for A and C. 
% 
% One can invert for either A or C individually, or both.
%
% The default option is to invert for log(A) and log(C) simultaneously.
%
% The objective function J (i.e. the function to be minimized) has the form
%
%  J=  I + R
%
% where I is a misift term, and R a regularisation term.
%
%
% The misfit term is:
%
%  I= (1/Area)   \int  (((u-uMeas)/uErrors)^2 + ((v-vMeas)/vErrors)^2) ) dx dy
%
% and the regularisation term can be either (Bayesian)
%
%  R= (C-Cprior) inv(KC) (C-Cprior)  +  (A-Aprior) inv(KA) (A-Aprior)  
%
% where KC and KA are covariance matrices, or (Tikhonov)
%
%  R= (1/Area)  \int (  gs^2 (grad (p-prior))^2  + ga^2 (p-prior)^2) dx dy
%
% where p is A or log(A), C or log(C)
%
% There are number of different minimisation methods implemented. Although the
% methodology behind the inversion is rigorous, in practice when working with
% real data the inversions sometimes get stuck in some local minima. The
% different optimisations methods implemented use slighlty different search
% directions, and switching methods may help getting out of a local minima as
% seen by one particular method. (When using synthetic data this is hardly ever
% an issue).
%
%
% The inversion for C and A can be done with C and A defined on nodes or
% elements. See: CtrlVar.AGlenisElementBased and CtrlVar.CisElementBased. In the
% past only inversion for element-based variables was possible, but now (as of
% Jan 2017) one can invert for either nodal or element values. By default, the
% inversion is done on nodal values.
%
%
% Hint: Often starting inverting for C using the fix-point method (see
% "FixPointEstimationOfSlipperiness" below) drives the misfit initially quite
% significantly down. Once that method stagnates (which it almost always will
% because the gradient used in that method is just a rough estimate and
% generally not exact), switch to another minimisation approach, for example the
% UaOptimisation using the adjoint gradients.
%
% Ua has some inbuilt optimisation methods and these are used by default.
% However, if the matlab optimisation toolbox is installed, the matlab routines
% can be used instead.
%
% Note #1: Some parameter combinations can be inconsistent. For example inverting
% for A only and applying regularisation on A and C, i.e.
%
%   CtrlVar.Inverse.InvertFor='logAGlen' ;
%   CtrlVar.Inverse.Regularize.Field='logAGlenlogC'
%
% is considered inconsistent (although in principle possible.) Also using the
% `FixPointC' gradient calculation, which only works for C inversion, and
% inverting for both A and C, i.e. 
%
%   CtrlVar.Inverse.DataMisfit.GradientCalculation='Adjoint' ; % {'Adjoint','FixPointC'}
%   CtrlVar.Inverse.InvertFor='logAGlenlogC' ; % {'C','logC','AGlen','logAGlen','logAGlenlogC'}
%
% is inconsistent. Ua tries to spot these input parameter mistakes and correct
% for them, but it is better to try to keep all inputs consistent.
%
% Note #2: It is possible to invert for any combinatin of log(A) or A  and log(C)
% or C. So for example one can invert for log(A) and C by setting 
%
%   CtrlVar.Inverse.InvertFor='logAGlenC' ;
%
% Also one can invert for log(C) and log(A) and regularise A and C by setting
%
%
%   CtrlVar.Inverse.InvertFor='logAGlenlogC' ;
%   CtrlVar.Inverse.Regularize.Field='AGlenC'
%


CtrlVar.Inverse.MinimisationMethod='UaOptimization'; % {'MatlabOptimization','UaOptimization'}


CtrlVar.Inverse.Iterations=1; % Number of inverse iterations

CtrlVar.Inverse.WriteRestartFile=1;  % always a good idea to write a restart file. 
CtrlVar.Inverse.NameOfRestartOutputFile='InverseRestart.mat';
CtrlVar.Inverse.NameOfRestartInputFile=CtrlVar.Inverse.NameOfRestartOutputFile;
CtrlVar.NameOfFileForSavingSlipperinessEstimate='C-Estimate.mat';
CtrlVar.NameOfFileForSavingAGlenEstimate='AGlen-Estimate.mat';
    
% It is usually better to invert for log(A) and log(C) rather than A and C.
% The default is to invert for log(A) and log(C) simultaneously.
CtrlVar.Inverse.InvertFor='logAGlenlogC' ; % {'C','logC','AGlen','logAGlen','logAGlenlogC'}

% The gradient of the objective function is calculated using the adjoint method.
% When inverting for C only, one can also use a gradient based on a `FixPoint'
% iteration, which is often a very good initial approach. 
CtrlVar.Inverse.DataMisfit.GradientCalculation='Adjoint' ; % {'Adjoint','FixPointC'}

% The gradient of the objective function can be premultiplied with the inverse
% of the mass matrix. This creates a `mesh independent' gradient. This has both
% advantages and disadvantages. 
CtrlVar.Inverse.AdjointGradientPreMultiplier='I'; % {'I','M'}


% Regularisation can be applied on A and C or log(A) and log(C). Also possible
% to use a covariance matrix for A and C. 
%
% Select Bayesian motivated regularisation by setting 
% CtrlVar.Inverse.Regularize.Field='cov' and Tikhonov regularisation
% by setting CtrlVar.Inverse.Regularize.Field to either 'C','logC','AGlen','logAGlen',or 'logAGlenlogC'
%
% Default is Tikhonov regularisation on log(A) and log(C)
CtrlVar.Inverse.Regularize.Field='logAGlenlogC' ; % {'cov','C','logC','AGlen','logAGlen','logAGlenlogC'}


% [ -- Parameters specific to Tikhonov regularisation See the above definition
% of the regularisation term R in the case of Tikhonov regularisation. The
% values of these parameters can be expected to be highly problem dependent. By
% default regularisation is switched on, but can the switched off by setting the
% gs and the ga parameters to zero.
CtrlVar.Inverse.Regularize.C.gs=1; 
CtrlVar.Inverse.Regularize.C.ga=1;
CtrlVar.Inverse.Regularize.logC.ga=1;
CtrlVar.Inverse.Regularize.logC.gs=1 ; 

CtrlVar.Inverse.Regularize.AGlen.gs=1;
CtrlVar.Inverse.Regularize.AGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.gs=1 ;
%  -]

% I and R are multiplied by these followign DataMisit and Regularisation
% multipliers. This is a convening shortcut of getting rid of either the misfit
% (I) or the regularization term (R) in the objective function (J) altogether.
CtrlVar.Inverse.DataMisfit.Multiplier=1;
CtrlVar.Inverse.Regularize.Multiplier=1;


% [----------  The following parameters are only relevant if using the
% UaOptimization i.e. only if
% CtrlVar.Inverse.MinimisationMethod='UaOptimization'; 
%
% The Ua optimisation is a simple non-linear conjugate-gradient method with
% automated resets, combined with a (one-sided) line search. The reset is done
% if the angle between subsequent steepest decent directions is to far from 90
% degrees, or if the update parameter becomes negative (only relevant for
% Polak-Ribiere and Hestens-Stiefel).
CtrlVar.Inverse.GradientUpgradeMethod='ConjGrad' ; %{'SteepestDecent','ConjGrad'}
CtrlVar.Inverse.InitialLineSearchStepSize=[];
CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize=1e-20; % minimum step size in backtracking
CtrlVar.Inverse.MinimumRelativelLineSearchStepSize=1e-5; % minimum fractional step size relative to initial step size
CtrlVar.Inverse.MaximumNumberOfLineSeachSteps=50;
CtrlVar.ConjugatedGradientsRestartThreshold=40 ; % degrees!
CtrlVar.ConjugatedGradientsUpdate='PR'; % (FR|PR|HS|DY)
                                        % FR ;Fletcher-Reeves
                                        % PR :Polak-Ribi\`ere
                                        % HR: Hestenes-Stiefel
                                        % DY :Dai-Yan
% end, UaOptimization parameters
% ------------]

% [------  The following parameters are only relevant if using the MatlabOptimisation option 
% i.e. only if CtrlVar.Inverse.MinimisationMethod='MatlabOptimization'
%
% Refer to the matlab documentation for further information. 
%
% The optimisation used is the matlab routine fminunc.
% You will need to have the matlab optimisation toolbox to be able to do this.
%  
% The Matlab optimisation toolbox has various algorithms to choose from, each of
% which has large number of parameters.
%
% Define the algorithm and set the options defining:  
%
%   CtrlVar.Inverse.MatlabOptimisationParameters
%
% Below are three examples that you might want to start with.
% The last one (which is the default one) uses fmincon with lbfgs update.
% 

CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'FunctionTolerance',1e-10,...
    'StepTolerance',1e-20,...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'SpecifyObjectiveGradient',true,...
    'HessianFcn','objective');

CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',1e-20,...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'SpecifyObjectiveGradient',true);


CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fmincon',...
    'Algorithm','interior-point',...
    'CheckGradients',false,...
    'ConstraintTolerance',1,...
    'Diagnostics','on',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','iter-detailed',...
    'FunValCheck','off',...
    'MaxFunctionEvaluations',100000,...
    'MaxIterations',CtrlVar.Inverse.Iterations,...,...
    'OptimalityTolerance',1e-2,...
    'OutputFcn',@fminuncOutfun,...
    'PlotFcn',[],... % 'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'StepTolerance',1e-10,...
    'FunctionTolerance',1,...
    'UseParallel',true,...
    'HessianApproximation',{'lbfgs',30},...
    'HessianFcn',[],...
    'HessianMultiplyFcn',[],...
    'InitBarrierParam',1,...
    'ScaleProblem','none',...
    'InitTrustRegionRadius',1,...
    'SpecifyConstraintGradient',false,...
    'SpecifyObjectiveGradient',true,...
    'SubproblemAlgorithm','cg');


% end, MatlabOptimisation parameters.   
% ------------]

% Some less-often used parameters related to inversion: 
CtrlVar.Inverse.InfoLevel=1;  % Set to 1 to get some basic information on J, R and I for each iteration, 
                              % >=2 for additional info on backtracking,
                              % >=100 for further info and plots
% In an inversion it it generally better to set other infolevels to a low value. So
% consider setting:
% CtrlVar.InfoLevelNonLinIt=0; CtrlVar.InfoLevel=0;

% [ ------------- Testing the adjoint gradients
% The derivatives obtained with the adjoint method can be
% compared with those obtained from brute force finite difference calculations.
% Only do this for small problems!
CtrlVar.Inverse.TestAdjoint.isTrue=0; % If true then perform a brute force calculation 
                                      % of the directinal derivative of the objective function.  
CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType='second-order' ; % {'first-order','second-order','fourth-order'}
                                         % The brute-force gradient can be calculated using first-order foward
                                         % differences, second-order central differences, or fourth-order central
                                         % differences. 
CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize=1e-8 ;
CtrlVar.Inverse.TestAdjoint.iRange=[] ;  % range of nodes/elements over which brute force gradient is to be calculated.
                                         % if left empty, values are calulated for every node/element within the mesh. 
                                         % If set to for example [1,10,45] values are calculated for these three
                                         % nodes/elements.
% end, testing adjoint parameters. 
% -------------------]          

CtrlVar.Inverse.DataMisfit.HessianEstimate='0'; % {'0','I','MassMatrix'} Do not use, just for testing. 
CtrlVar.Inverse.CalcGradI=true;   % do not change, just for testing
CtrlVar.Inverse.DataMisfit.FunctionEvaluation='integral';   % do not change, just for testing
CtrlVar.Inverse.DataGradient.FunctionEvaluation='integral'; % do not change, just for testing
CtrlVar.Inverse.StoreSolutionAtEachIteration=0; % if true then inverse solution at each iteration is saved in the RunInfo variable.



%% Numbering of nodes and elements
CtrlVar.sweep=1;              % renumber nodes using a `sweep' plane
CtrlVar.SweepAngle=0.01;      % angle of sweep plane with respect to x axis
CtrlVar.CuthillMcKee=0;       % renumber nodes using sparse reverse Cuthill-McKee ordering

%% Creation of a dumpfile
% Mainly used for testing purposes, but can in principle also be used to generate output data files.
CtrlVar.WriteDumpFile=0;                      % a dumpfile is created containing all variables
CtrlVar.WriteDumpFileStepInterval=1000;       % number of time steps between writing a dump file
CtrlVar.WriteDumpFileTimeInterval=0;          % time interval between writing a dump file

%%  Outputs
%
% For outputs Ua calls a routine called 'UaOutputs.m'
% Write your own version of this routine to fit your own output/plotting needs and keep 
% the routine into you local run directory, i.e. the directory from which you run Ua
% Start by copying the example UaOutput.m routine from the Ua source installation folder
% to you local run directory.
%
%
%
%


CtrlVar.UaOutputsDt=0; % model time interval between calling UaOutputs.m
                       % if set to zero UaOutputs is called at every time/run step
                       % if set to a negative number, or NaN, UaOutputs is never called
CtrlVar.UaOutputsMaxNrOfCalls=NaN;  % maximum nr of calls to UaOutputs
                                    % Once this limit is reached, the run stops. (Setting this to 1 or some low number
                                    % can sometimes be useful for testing/control purposes)
                                    % NaN implies no limit to the number of calls 
                                    
          
%% Optaining information about the run, during the run.
%
% A simply way of getting information about the run from within the user m-files
% is by inspecting the fields of the CtrlVar.  The CtrlVar is given an in input
% to all such m-files.
%
% For example, the counter CtrlVar.CurrentRunStepNumber gives the current
% run-step number. 
%
CtrlVar.CurrentRunStepNumber=0 ;  % This is a counter that is increased by one at each run step.
                                   
%% General Meshing Options
% There are various ways of meshing the computational domain.
%
% In almost all cases the simplest option tends to be to define the outlines of
% the computational domain in Ua2D_InitialUserInput.
%
% In that case �a will call an external mesh generator. The external mesh
% generator used by Ua is "gmsh" which is a well known and a well supported open
% source mesh generator (http://geuz.org/gmsh/) The outlines of the mesh are
% defined by the variable 'MeshBoundaryCoordinates' set in
% Ua2D_InitialUserInput.m. This approach is quite flexible and allows for
% complicated computational domains containing holes and/or separated domains.
%
% *For examples of how to generate different type of meshes look at* *ExamplesOfMeshGeneration.m*
%
% Both when done from within �a or externally, generating a FE mesh with the mesh generator `gmsh' typically involves:
%
% *             a) create an input file for gmsh (.geo)
% *             b) call gmsh for that input file (.geo). gmsh in turn generates an output file (.msh)
% *             c) read into �a the resulting gmsh output file (.msh) with the mesh
% All, or some of these three steps can be done withing �a.
%
% More specifically the options are:
%
% *    i)  Directly read existing gmsh output file (.msh)
% *   ii)  First run gmsh with an existing gmsh input file (.geo) and then read the resulting gmsh output file (.msh)
% *   iii) First generate gmsh input file (geo), then run gmsh for that input file, and finally read the resulting gmsh output file (.msh)
%
% Option iii is the default option, in which case �a generates the gmsh input file (.geo), calls gmsh, and then reads the resulting gmsh output file with the mesh.
%
% To select between i, ii and iii set CtrlVar.GmshMeshingMode={'load .msh','mesh domain and load .msh file','create new gmsh .geo input file and mesh domain and load .msh file'}
%
% CtrlVar.GmshMeshingMode='load .msh'                                                               % option i
% CtrlVar.GmshMeshingMode='mesh domain and load .msh file'                                          % option ii
 CtrlVar.GmshMeshingMode='create new gmsh .geo input file and mesh domain and load .msh file';     % option iii, which is the default option
% 
% After having generated a FE mesh, that FE mesh can then be read in as an initial mesh at the start of other runs.
% 
CtrlVar.ReadInitialMesh=0;    % if true then read FE mesh (i.e the MUA variable) directly from a .mat file 
                              % unless the adaptive meshing option is used, no further meshing is done.
CtrlVar.ReadInitialMeshFileName='ExistingMeshFile.mat';
CtrlVar.SaveInitialMeshFileName='NewMeshFile.mat';
% By default, the mesh is always saved into a file, and that file can later be re-read.
% But to generate a new mesh file from, for example a result file or a restart file, is easy. Just load the restart/result file and save MUA to a file 
% So for example:  load Restartfile ; save MyNewMeshFile MUA
% Now `MyNewMeshFile.mat' is a file that can be used as an initial mesh file by setting CtrlVar.ReadInitialMesh=0; CtrlVar.ReadInitialMeshFileName='MyNewMeshFile.mat';

CtrlVar.MeshGenerator='gmsh';  % possible values: {mesh2d|gmsh}
CtrlVar.GmshFile='GmshFile';  % name of gmsh input/output files (no file extensions)

CtrlVar.GmshMeshingAlgorithm=1;     % see gmsh manual
                                    % 1=MeshAdapt
                                    % 2=Automatic
                                    % 5=Delaunay
                                    % 6=Frontal
                                    % 7=bamg
                                    % 8=DelQuad (experimental)
                                    
CtrlVar.GmshVerbosityLevel=1;    % see gmsh manual, higher values give more information
CtrlVar.GmshPause=0;      % very occasionally gmsh returns an error when run from within matlab
                          % but runs OK if run outside of matlab for exactly the same problem (!).
                          % The reasons for this are not clear, possibly related to delayed writing of
                          % files and some syncronisation issues. Possibly remedy is to introduced a short
                          % pause before calling gmsh. GmshPause>0 creates such a pause.
                          % The duration of the pause is measured in seconds.
                          
                          
CtrlVar.GmshInputFormat=1; % When using �a to call Gmsh, the input to Gmsh as defined in Ua2D_InitialUserInput 
                           % can be given in two different ways, i.e. GmshInputFormat=1 or 2. 
                           % Format 1 is simpler
                           % Format 2 is closer to the actual input format of Gmsh (.geo) and is more
                           % flexible. See ExamplesOfMeshGeneration.m for further description and examples.
CtrlVar.GmshBoundaryType='lines';   % (spline|lines)
CtrlVar.GmshCharacteristicLengthExtendFromBoundary=0;
CtrlVar.GmshCharacteristicLengthFromCurvature = 0 ;
CtrlVar.GmshGeoFileAdditionalInputLines{1}='   ';  % these lines are added to the gmsh .geo input file each time such a file is created

CtrlVar.OnlyMeshDomainAndThenStop=0; % if true then only meshing is done and no further calculations. Useful for checking if mesh is reasonable
CtrlVar.AdaptMeshAndThenStop=0;      % if true, then mesh will be adapted but no further calculations performed

%% Controlling element sizes
% 
% if no adaptive meshing is used then the element size is given by
CtrlVar.MeshSize=10e3;                       % over-all desired element size (however if gmsh is used without adaptive meshing
                                             % only CtrlVar.MeshSizeMin and CtrlVar.MeshSizeMax are used)
                                             % 
CtrlVar.MeshSizeMin=0.1*CtrlVar.MeshSize;    % min element size
CtrlVar.MeshSizeMax=CtrlVar.MeshSize;        % max element size

CtrlVar.MaxNumberOfElements=100e3;           % max number of elements. If #elements larger then CtrlMeshSize/min/max are changed
CtrlVar.MaxNumberOfElementsUpperLimitFactor=1.3;  % if actual number of elements is larger than CtrlVar.MaxNumberOfElements by this factor
                                                  % the domain is remeshed by modifying MeshSizeMin 
CtrlVar.MaxNumberOfElementsLowerLimitFactor=0.0;
% Note that the `MeshSize' part of the names of these variables is possibly somewhat
% misleading. These variables relate to the size of the elements not the overall
% size of the computational mesh (which is determined by
% MeshBoundaryCoordinatates).

%% Options related to the Ua mesh structure variable MUA
CtrlVar.MUA.MassMatrix=false;       % true if the mass matrix is to be computed and stored as a part of MUA
CtrlVar.MUA.StiffnessMatrix=false;  % true if the stiffness matrices is to be computed and stored as a part of MUA

%% Pos. thickness constraints,          (-active set-)
% A minimum ice thickness can be enforced in different ways using the following methods:
%  1) `reset method' : simply resetting the thickness to min thickness at node where thickness is less than a prescribed value.
%  2) `active-set' method.
%  3) `thickness-barrier' method
%
% The active-set method is the preferred option and is arguably the only correct way of enforcing min ice thickness.
% The active-set method should therefore be used whenever possible.
% However, in some cases the active set method does not converge, in which case
% options 1) or 3), or combinations thereof, must be used.  If the differences between approach 1) and 2) are small, then using 1)
% allows for shortest computation times
%
% The thickness-barrier method introduces a fictitious surface mass balance term.
% The thickness-barrier method can be used on its own, but should primarily be used in combination with the active-set method
% to improve convergence.
%


CtrlVar.ThickMin=1;                      % minimum allowed thickness without (potentially) doing something about it

% reset method, option 1
CtrlVar.ResetThicknessToMinThickness=0;  % set to 1 to reset thickness values less than ThickMin to ThickMin at each time step (Option 1, not recommended)
CtrlVar.ResetThicknessInNonLinLoop=0;    % if true, thickness in the non-linear iteration of the uvh implicit approach
                                         % is set to zero, provided CtrlVar.ResetThicknessToMinThickness is also true (usually not a good idea)


% active-set method, option 2 
CtrlVar.ThicknessConstraints=1;             % set to 1 to use the active-set method (Option 2, the recommended option).
CtrlVar.ThicknessConstraintsItMax=10  ;     % maximum number of active-set iterations.
                                            % if the maximum number of active-set iterations is reached, a warning is give, but
                                            % the calculation is not stopped. (In many cases there is no need to wait for
                                            % full convergence of the active-set method for each time step.)
                                            % if set to 0, then the active set is updated once and then proceed to next time step.
CtrlVar.ThicknessConstraintsLambdaPosThreshold=0;  % if Thickconstraints are larger than this value they are inactivated, should be zero
CtrlVar.NumberOfActiveThicknessConstraints=0;      % The number of active thickness constraints (just for information, always set initially to zero)
CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints=1000 ; %

% thickness barrier, option 3
CtrlVar.ThicknessBarrier=0;                   % set to 1 for using the barrier method  (Option 3)
CtrlVar.ThicknessBarrierThicknessScale=CtrlVar.ThickMin;     %
CtrlVar.ThicknessBarrierDiagonalFraction=1;   % size of barrier term in comparison to mean abs of diagonal elements
CtrlVar.ThicknessBarrierMinThickMultiplier=2; % exp. barrier is 1 at ThickMin * MinThickMuliplier
CtrlVar.ThicknessBarrierAccumulation=0.01;

%% Advance/Retreat mesh and activation/deactivation of elements
% This option allows for deactivation/activation of elements based on ice thickness.
% A `background' FE mesh is required. In most cases this background FE mesh will simply be the initial FE mesh
% used at the start of the calculation.
% For advancing glaciers this option must be combined with the active-set method (set CtrlVar.ThicknessConstraints=1)
%
% Note: When combined with the active-set method
% then CtrlVar.ThickMinDeactivateElements must be >= CtrlVar.ThickMin.
% It is usually good for this value to be slightly larger than CtrlVar.ThickMin
% Nodes are only included in the active set if thickness at a node < CtrlVar.ThickMin.
% If CtrlVar.ThickMinDeactivateElements>CtrlVar.ThickMin elements will be eliminated
% before all nodes of that element have (potentially) been included in the active set.
% This reduces somewhat the number of nodes in the active set from what it would otherwise be
% if CtrlVar.ThickMinDeactivateElements=CtrlVar.ThickMin.
%
% Note: Elements are only inactivated if ALL nodes have thickness <CtrlVar.ThickMinDeactivateElements
% Elements are activated once at least one of the nodes leaves the active set. It is therefore
% possible to have a new element where some of the nodes have thickness larger than CtrlVar.ThickMin.
% but less than CtrlVar.ThickMinDeactivateElements.
CtrlVar.FEmeshAdvanceRetreat=0;     % activates the Advance/Retreating mesh option
CtrlVar.FEmeshAdvanceRetreatDT=0.5; % activation/deactivation done at this time interval
                                    % for CtrlVar.FEmeshAdvanceRetreatDT=0 the activation/deactivation is done
                                    % at every time step (in many cases the best approach)
CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName='BackgroundMeshfile.mat'; % This file is needed for the advance/retreat option
                                                                             % It must contain the variable `MUA_Background'
                                                                             
CtrlVar.ThickMinDeactivateElements=1.01*CtrlVar.ThickMin;% Elements where thickness at all nodes is less than this value are deactivated
CtrlVar.SelectElementsToDeactivateAlgorithm=1; % (1|2)  There are two different methods implemented for selecting
                                                % elements to be deactivated in conjunction with the FEmeshAdvanceRetreat option:
                                                % 1: Eliminate an element if none of the nodes of that element belong to an element
                                                %    where any of the nodal thicknesses are greater than CtrlVar.ThickMinDeactivateElements
                                                % 2: Eliminate elements where all the nodes are less or equal to CtrlVar.ThickMinDeactivateElements
                                                % Method 1 eliminates less (or equal) number of elements than Method 2.
                                                % Method 2 is OK for retreating cases and some advancing cases but can fail                                             
                                                % if the advance is `too quick'.

CtrlVar.MinSurfAccRequiredToReactivateNodes=0;  % If surface accumulation is larger than this, then a node is considered to have positive ice thickness
                                                % and not eliminated.  This is important in cases where, for example, with time the ELA drops down below 
                                                % the top of a mountain peak that is not included in the current FE-mesh.
                                                % This allows for the formation of new isolated glaciated areas.
                                                % Although the default value is zero, it is presumably better to set this to a small positive value.


%% Mesh refinement: Uniform global mesh refinement
% Mesh can be refined at a start of a run or the start of a restart run by subdividing all triangles into four
% can be useful, for example, for an error estimation
CtrlVar.RefineMeshOnRestart=0;
CtrlVar.RefineMeshOnStart=0;
%% Mesh refinement: Global and local adaptive mesh refinement
% There are various adapt meshing options.
% The most general one is global remeshing using explicit error estimate
%
% Global remeshing can be based on one or more of the following
% relative RefineCriteria:
%
% * 'effective strain rates'
% *          '|dhdt|'
% *        '||grad(dhdt)||'
% *          'dhdt curvature'
% *         'thickness gradient'
% *         'thickness curvature'
% *         'flotation'
% *         'f factor'
%
% These (relative) criteria can be combined. When two or more criteria are combined
% RefineCriteria is given as a cell array
%
% The relative importance of different RefineCriteria can be specified by
% defining `RefineCriteriaWeights'. These weights affect how small the smallest
% element will be for a given refinement criteria.
%
% If CtrlVar.RefineCriteriaWeights=1, the whole range CtrlVar.MeshSizeMin to
% CtrlVar.MeshSizeMax is used.
%
% If CtrlVar.RefineCriteriaWeights=0.5 element size will range from
% CtrlVar.MeshSizeMax down to
% CtrlVar.MeshSizeMin+(CtrlVar.MeshSizeMax-CtrlVar.MeshSizeMin)*(1-RefineCriteriaWeight)
%
% If CtrlVar.RefineCriteriaWeights=0 the criterion is effectively ignored.
%
% Examples:
%
%  CtrlVar.RefineCriteria='effective strain rates';  % specifies 'effective strain rates' as the only criterion
%  CtrlVar.RefineCriteriaWeights=[1];                % with a relative weight of unity
%
%  CtrlVar.RefineCriteria={'flotation','||grad(dhdt)||','dhdt curvature','thickness curvature'}; % several criteria used
%  CtrlVar.RefineCriteriaWeights=[1,1,1];  %
%
% In addition the refinement can be limited to an area within a given flotation distance by defining
% CtrlVar.RefineCriteriaFlotationLimit
% For example, for CtrlVar.RefineCriteriaFlotationLimit=[100,NaN]
% the first refinement criteria will only be applied to area where the glacier bed (b) is within 100 vertical distance units from flotation
%
% One can also specify directly the desired element sizes (explicit:global
% option) or the elements to be refined (explicit:local option), using the user
% m-file `DefineDesireEleSizes.m'
%
% Note: Adapt meshing can only be done in a combination with a forward run.
% Adapt meshing can not be done in an inverse run. Usually before starting with
% an inverse run you may find it usefull to do a number of forward (time
% independent) forward runs and use those to optimize the mesh prior to start of
% any inverse runs.
%
% In addition to the above listed realtive refinedment criteria, one can also
% specify one type of an absolut mesh criterion based on distance from grounding
% lines.
%
% The desired sizes of elements within a given distance from any grounding lines
% can be specified through CtrlVar.MeshAdapt.GLrange
%
% CtrlVar.MeshAdapt.GLrange is an n x 2 array. Each value pair in every line
% specifies a distance from the grounding line and the desired element size
% within that distance. For example setting
%
%   CtrlVar.MeshAdapt.GLrange=[5000 1000];
%
% specifies that all elements witin 5000 meters should be 1000 m large (here
% assuming the distance unit is meters) if all elements.
%
% Setting
%
%   CtrlVar.MeshAdapt.GLrange=[5000 2000 ; 1000 500  ; 250 50];
%
% specifies that elements within 5000 meters from the grounding line should be
% at the most 2000 meters large, those within 1000 m at most 500 m larger, and
% those within 250 m, 50 meters in size. 
%
% If no such mesh criterion is to be specified, set to an empty value, i.e 
%
%   CtrlVar.MeshAdapt.GLrange=[];                                                    
%
% Note: This absolut mesh criterion requires the matlab function rangesearch
% which is a part of the machine learning matlab toolbox.
%
CtrlVar.AdaptMesh=0;          % true if adapt meshing is used, no remeshing is done unless this variable is true
CtrlVar.MeshRefinementMethod='explicit:global';    % can have any of these values:
                                                                                             % 'explicit:local'
                                                   % 'implicit:global'  (broken at the moment, do not use)
                                                   % 'implicit:local'   (broken at the moment, do not use)

% `explicit:global' implies a global remeshing of the whole domain. This is a very flexible approach 
%  allowing for both increasing and decreasing mesh resolution.
% 'explicit:local' implies a local adaptive mesh refinement obtained by splitting
% individual triangles up into four sub-triangles. This is often a very elegant
% way of refining the mesh, but does not allow for subsequent mesh coarsening.
%
                                                                                                      
                                                   
CtrlVar.AdaptMeshInitial=1  ; % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshInterval)~=0.
CtrlVar.AdaptMeshInterval=1 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0

CtrlVar.hpower=1;         % used to go from an error estimate to a size estimate for an element
                          % h=1/error^hpower ,  where `h' is the desired element size and `error' a local
                          % error estimate.

CtrlVar.AdaptMeshIterations=1;  % Number of global adapt mesh iterations within each adapt step
                                % This is seldom anything else but 1, except potentially in a either diagnostic calculation (time independent)
                                % or at the start of a transient (prognostic) calculation where the initial mesh is very coarse
                                % Note that when using gmsh the mesh refinement is always based on indicators given at the nodal points of the original mesh
                                % therefore using a few AdaptMeshIterations can sometimes be be a good idea.

                                
CtrlVar.LocalAdaptMeshSmoothingIterations=5;  % Number of Laplace mesh smoothing iterations used in local mesh refinement
CtrlVar.LocalAdaptMeshRatio=0.25;             % The maximum number of elements subdivided during each local mesh refinement step
                                              % as a fraction of the total number of elements.

CtrlVar.MaxRatioOfChangeInEleSizeDuringAdaptMeshing=5;   % put a strict limit on how much ele sizes change during single
CtrlVar.MinRatioOfChangeInEleSizeDuringAdaptMeshing=1/5; % adaptive meshing step to avoid excessive changes.
                                                         % This does not apply to local mesh refinement where in each adapt step
                                                         % the elements are always only refined, and then always by a factor of two.

CtrlVar.RefineCriteria={'flotation','||grad(dhdt)||','dhdt curvature'};
CtrlVar.RefineCriteriaWeights=[0.1,1,1];                %  
CtrlVar.RefineCriteriaFlotationLimit=[NaN,NaN,NaN];     % Refine criteria is only applied to elements which are this close to flotation
                                                        %        useful to restrict refinement to an area in the (vertical) vicinity of the grounding line
                                                        %        Set to NaN if to be ignored and applied to all regions irrespective of how close to flotation

CtrlVar.NumberOfSmoothingErrorIndicatorIterations=1;    % each of the error indicators can be smooth over neighbouring elements
                                                        % this is done by calculating average values for each element based its nodal values
                                                        % and then interpolating from elements to nodes using number of elements that a node is
                                                        % attached to as a weighting factor.  This introduces a smoothing that is related to
                                                        % connectivity as opposed to spatial distance.
                                                        % This kind of smoothing is never done for the 'flotation' and the `f factor' cases
                                                        % as the spread/smoothing can be determined directly by CtrlVar.RefineDiracDeltaWidth

% absolut mesh adapt criterion                                                        
CtrlVar.MeshAdapt.GLrange=[];                                                    
                                                        
CtrlVar.RefineDiracDeltaWidth=100;  % for `flotation' and 'f factor' the zone within this vertical distance from flotation is refined
CtrlVar.RefineDiracDeltaOffset=0;   %
                              
                                

%% `Time geometries' are boundary geometries that change with time.
% This can be used, for example, to simulate a calving event.
CtrlVar.TimeGeometries.Flag=0;             % true if domain geometry is changed during the run (e.g prescribed calving event)

%% Mesh adjustments: Mesh morphing:
%
% (mesh morphing around a moving grounding line is currently broken. This 
% looked like a good idea, but really is only going to work if the grounding line has
% a simple shape and the topology of the grounding lines does not change.
% Basically a too limited option for practical use.)
CtrlVar.MeshMorphing=0;       % true for mesh-morphing where the mesh is morphed onto moving grounding line
                              % this is a very elegant method, but only works if there is just one grounding line
                              % and therefore not really that useful in a general 2HD situation.
CtrlVar.GLmeshing=0;          % GL meshing based on morphing
CtrlVar.GLtension=1;          % tension of spline used in GL morphing, 1: no smoothing; 0: straight line
CtrlVar.GLds=CtrlVar.MeshSizeMin ; % edge length along GL when using GL meshing




%% Parameters affecting the floating mask

CtrlVar.kH=1;   % kH -> infty gives an exact Heaviside and delta functions.
                % kH=1 implies a grounding line "width" of 1 m up and down from floating condition
                % kH=10 implies a grounding line "width" of 1/10 m up and down from floating condition
CtrlVar.Hh0=0;  % offset is Heaviside function when calculating GF field

%% Parameters affecting calculation of grounding line
% The grounding line position does not enter any calculations done by �a. 
% The grounding line is primarily calculated for plotting purposes.
CtrlVar.GLthreshold=0.5;  % used to define position of GL with respect to the values of the Heaviside function (1 fully grounded, 0 fully floating)
CtrlVar.GLsubdivide=0;    % If 0/false the grounding line is determined based on GL.node values at corners only (using GLthreshold). If 1/true
                          % then all nodal values of 6-node and 10-node triangles are also used. This is done by splitting those into 4 and 9 triangles, respectively



%% A and C as element or nodal variables
% AGlen and C can be either nodal or element variables.
%
% By default A and C are nodal variables.
%
%
CtrlVar.AGlenisElementBased=0; 
CtrlVar.CisElementBased=0;
CtrlVar.AutomaticallyMapAGlenBetweenNodesAndEleIfEnteredIncorrectly=1;


%% Adaptive Time Stepping Algorithm (ATSA)   (adapt time step)
% The adaptive-time-stepping algorithm is based on the idea of keeping the number of non-linear iterations
% close to a certain target (CtrlVar.ATSTargetIterations).
% This is a simple but highly effective method.  However, as the ATSA is not based on any kind of error estimates,
% it does not guarantee that errors will not accumulate, etc, and ATSA does not work for linear problems.
%
% The main idea is to aim at a time step that limits the number of non-linear iteration to a relatively small number.
% If so, then most likely the Newton-Raphson iteration is in the quadratic regime.
% Experience has shown that a target number of iterations (CtrlVar.ATSTargetIterations) within 3 to 5 is good for this purpose
%
% Time step is increased if r<1 where
%
%     r=N/M
%
%   where
%   N is the max number of non-linear iteration over last n time steps
%   M is the target number of iterations
%
%   here
%     M=CtrlVar.ATSTargetIterations
%   and
%     n=CtrlVar.ATSintervalUp
%
%   (N does not need to be specified.)
%
%  Currently the time step is only decreased if either:
%        a) the number of non-linear iterations in last time step was larger than 25
%        b) number of iterations over last n times steps were all larger than 10
%  where n=CtrlVar.ATSintervalDown
%
% There are some further modifications possible:
%  -time step is adjusted so that time interval for making transient plots (CtrlVar.TransientPlotDt) is not skipped over
%  -time step is not increased further than the target time step CtrlVar.ATStimeStepTarget
%  -time step is adjusted so that total simulation time does not exceed CtrlVar.TotalTime
%
%
%
CtrlVar.AdaptiveTimeStepping=1 ;    % true if time step should potentially be modified
CtrlVar.ATStimeStepTarget=1000.0;   % maximum time step size allowed
CtrlVar.ATStimeStepFactorUp=2 ;     % when time step is increased, it is increased by this factor
CtrlVar.ATStimeStepFactorDown=10 ;  % when time step is decreased, it is decreased by this factor
CtrlVar.ATSintervalUp=5 ;           %
CtrlVar.ATSintervalDown=3 ;         %
CtrlVar.ATSTargetIterations=4;      % if number of non-lin iterations has been less than ATSTargetIterations for
                                    % each and everyone of the last ATSintervalUp iterations, the time step is
                                    % increased by the factor ATStimeStepFactorUp
                                    
                                    
%% Mass-balance geometry feedback
% If the mass balance is a function of geometry, an additional non-linearity is introduced to transient runs.
% This non-linearity can be solved in a fully consistent way using the Newton-Raphson method provided the user
% supplies the gradient of the mass balance with respect to thickness.
%
CtrlVar.MassBalanceGeometryFeedback=0;  % If the mass balance depends on geometry then
                                        % setting this parameter to either 1, 2 or 3 has the effect of 
                                        % the mass-balance being updated within the non-linear transient-loop.
                                        % In principle this parameter should always be set to 3, but in practice the
                                        % feedback can often be sufficiently well accounted for by simply updating 
                                        % mass balance at each and every time step (i.e. option 0).
                                        %  
                                        %  0 : no mass-balance geometry feedback considered within non-lin iteration loop 
                                        %      (however, as always, mass balance is updated at each time step)
                                        %  1 : mass-balance feedback included at the start of each non-lin iteration, 
                                        %      but not within the backtracking step.
                                        %  2 : Feedback included in non-lin loop, both at the beginning of each NR iteration, 
                                        %      and within backtracking step.  
                                        %  3 : Consistent mass-balance feedback algorithm. As option 2, but with 
                                        %      the gradient of the mass-balance with respect to thickness added to the
                                        %      left-hand side of the NR system. Requires the user to supply this gradient through 
                                        %      `DefineMassBalance.m'. Doing so can lead to a drastic reduction 
                                        %      in the number of NR steps required.
                                        %
                                        %  If mass balance depends on geometry then always try to use:
                                        %  CtrlVar.MassBalanceGeometryFeedback=3
                                        %  and then also give dasdh and dabdh as return arguments in DefineMassBalance.
                                        %
                                        %  However, if there is no dependency of the mass balance on geometry, always set
                                        %  CtrlVar.MassBalanceGeometryFeedback=0
                                        %  as doing so avoids calls to DefineMassBalance.m within the non-line loop.
                                        %
CtrlVar.MassBalanceGeometryFeedbackDamping=0;  % Dampens the update in surface mass balance.
                                               % If not equal to zero, then the actual mass-balance value used at the end of the time step,
                                               % becomes a weighted average of that at the beginning and the (correct) value at the 
                                               % end of the time step.
                                               % The value must be in the range [0,1]
                                               % Only use this if encountering convergence problems.  
                                               % Should always be equal to 0 if possible.
                                               % If not equal to 0, the algorithm converges to a wrong solution (!),
                                               % although the error might be very small if mass-balance geometry feedback is not that strong.
      

%% Sea ice/melange                                               
%
% �a has some (simple) ice-melange/sea-ice physics that allow for ocean and athmospheric
% drag acting over the floating sections.
%
% If used, then the drag parameters are defined in 'DefineSeaIceParameters'
%
CtrlVar.IncludeMelangeModelPhysics=0;


%%
CtrlVar.MeltNodesDefinition='Edge-Wise';

%% Mapping from Mesh1 to Mesh2
% when a new FE is created, the values from the old mesh need to be mapped onto the new mesh
% if the boundary of the mesh has not changed this only involves interpolation
% but if the boundary of Mesh2 is different from that of Mesh1 then extrapolation might be involved
% In this case one must check for such points and give them some sensible outside value.
% 
CtrlVar.Mesh1To2CheckForPointsOutsideMesh1AndInsideConvexHull=1 ; % for non evolving mesh boundaries, can be set to 0/false
CtrlVar.InpolyTol=0.1;       % tolerance when checking inside outpoints using the `inpoly' m-file, should be small compared to size of any element

%%
% Parallel options:
%
% 
% The parallel profile is not modified within �a. Set the properties of the local
% profile through the general Matlab settings. See the matlab manual for further
% information. If needed, the properties of the local profile can be adjusted in
% the Ua2D_InitialInput file
%
% For example, to change the number of local workers to 6, one can do the
% following: myCluster = parcluster('local') myCluster.NumWorkers = 6;
% saveProfile(myCluster)
%
% Consult the matlab manual for further information
%
CtrlVar.ParallelAssembly=1;

%%
CtrlVar.fidlog=1;  % unit number for standard output, no need to change.

%%

 CtrlVar.DevelopmentVersion=1;  % Internal variable, always set to 0 (unless you want to use some untried, untested and unfinished features....) 

end



