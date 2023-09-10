function [LSF,xcEdges,ycEdges,ShapeDifference]=CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,LSF,options)

%
% Takes initial calving front data (Xc,Yc) and a level set function LSF, which on input only must have the correct sign, and
% initializes the level set based on gemetrical distance from the (Xc,Yc) calving front. Then recomputes the calving front
% from the level set.
%
% The level-set provided on input is initialized based on the (Xc,Yc) values, and the LSF provided on input will be
% recalculated.
%
% Only the sign of the LSF function on input is used. This sign must be correct in the sense that LSF<0 for nodes
% downstream of calving fronts, and LSF>0 for nodes upstream of calving fronts.
%
% The returned calving front, (xc,yc), may not be identical to the calving front (Xc,Yc) provided as input. There are a
% number of (subtle) reasons for this. The (xc,yc) calving front is (re) calculated from the nodal values of level set field
% (LSF), calculated (ie initialized) using the (Xc,Yc) values. Despite best atempts, this operation, ie calculating c=(xc,yc)
% from LSF, which we can write c=c(LSF), is not an exact inverse of the function LSF=LSF(c).
% 
% This can be understood as being related to lack degrees-of-freedom when calculating LSF nodal values, as a single LSF nodal
% value will always affect the positions of the calving fronts across all elements to which it belongs. We can therefore have
% a situation where it is impossible to find distribution of LSF nodal values for which the calculated calving front is
% identical to the calving front as given on input.
%
%
% If the calving front is provided on input as (Xc,Yc), it can optionally be subsampled within elements. Generally this is a
% must if one is performing a re-initialisation during a run, as otherwise the geometric re-initialzation can shift the level
% set. But for the first initialisation this is hardly requried as it only impacts the initial location of the calving front,
% and if the calving front is alread known at a high resolution as (Xc,Yc) from some observations, such resampling is
% superfluous.
%
% There is a special case when Xc and Yc are entered as empty arrays, in which case Xc and Yc are simply determined by
% calculating the zero contours of LSF. This can be OK, but the resulting calving front (xc,yc) will, in general, be impacted
% by the mesh resolution with the resulting calving fronts typically cutting through the mid edges of the elements. If
% (Xc,Yc) is known, the return calving front (xc,yc) is much more likely to cut across the elements similar as (Xc,Yc).
%
%
% On input only the sign of LSF must be correct. This is, for nodes inside/upstream of calving fronts LSF must be positive,
% and for nodes outside/downstream of calving fronts, negative.
%
% On return LSF has be initialized geometrically based on the distance from nodes to the (xcEdges,ycEdges) intersections
% along element edges
%
% (xc,yc) are the calving front positions corresponding to the zero level of the level-set function
%
% Ideally, on  return (xc,yc) should be close to (Xc,Yc)
%
% The "EleEdges" method is slower than "InputPoints" method, and possibly in practice the resulting differences are not that
% significant.
%
% When using the "InputPoints" method, it might be good to resample the calving front at equal intervals by setting
%
%
%   ResampleCalvingFront=true
%   CalvingFrontPointDistance=1000 ; % the distance between points along the re-sampled calving front
%
%
% Examples:
%
%
%   [xc,yc,LSF,xcEdges,ycEdges,ShapeDifference]=CalvingFrontLevelSetGeometricalInitialisation(CtrlVar,MUA,Xc,Yc,F0.LSF,plot=true) ;
%
%


ShapeDifference=[];


UserVar=[] ; RunInfo=[];




if nargin > 4 && ~isempty(Xc)
    P1=[Xc(:) Yc(:)] ;
else
    fprintf("CalvingFrontLevelSetGeometricalInitialisation: Xc and Yc empty on input. \n \t \t \t Xc and Yc calculated as the zero contour lines of LSF.\n")
    Value=0;
    options.lineup  = true;
    options.plot  = false;
    options.subdivide  = false;
    [Xc,Yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,LSF,Value,options) ; 
    P1=[Xc(:) Yc(:)] ;
end



con=MUA.connectivity;
coo=MUA.coordinates;

% To speed up the calculation consider to:
% 1) get rid of any sections of the calving front outside of the computational boundary
% 2) subsample points along the calving front.
%
% But presumably all such manipulations with the input data are better done ahead of
% the call.

if options.GetRidOfCalvingFrontOutsideComputationalDomain
    % 1) Get rid of any points of the calving profile outside the mesh domain
    % But be carfulle with this as it does not deal with more than one calving front!
    io=inpoly2(P1, [MUA.Boundary.x MUA.Boundary.y]);
    P1=P1(io,:) ;
end

if options.ResampleCalvingFront
    % and 2) subsample the points along the calving profile
%%
    ds=options.CalvingFrontPointDistance;

    nPoints=size(P1,1);

    % Check if there are several seperate calving fronts, and resample each individually
    iNaN=find(isnan(P1(:,1)));
    ii=[0;iNaN;nPoints+1] ;

    xi=[]; yi=[];
    for k=1:numel(ii)-1
        
        px=P1(ii(k)+1:ii(k+1)-1,1);
        py=P1(ii(k)+1:ii(k+1)-1,2);


        arclength = sum(sqrt(sum(diff([px py],[],1).^2,2)));
        Npoints=round(arclength/ds);

        pt= interparc(Npoints,px,py,'linear') ;
        xi=[xi;nan;pt(:,1)];
        yi=[yi;nan;pt(:,2)];
    end

    xi(end)=[]; yi(end)=[];
%%
    P1=[P1 ; xi(:) yi(:)] ;  % add the new ones, keep the original ones as well
    %P1=[xi(:) yi(:)] ;


end




switch options.method

    case "EleEdges"

        %
        %  Here the intersections between individual ele edges and the calving front profile as defined by (Xc,Yc)
        %  are determined, and the LSF is then initialized as the distance to those intersection points

        switch MUA.nod

            case 3

                EleEdges=con(:,[1 2 3 1]); % These are element edges, this is only for 3-node elements.
                % In general, these should be the nodal numbers of the corner nodes
                % so that this traces a closed loop for every element.

            case 6

                EleEdges=con(:,[1 3 5 1]); % These are element edges

                % This does not make it any more accurate, only slower
                %con=TriFE(MUA.connectivity);
                %EleEdges=con(:,[1 2 3 1]); % These are element edges, this is only for 3-node elements.



            case 10

                EleEdges=con(:,[1 4 7 1]); % These are element edges

                % This does not make it any more accurate, only slower
                %con=TriFE(MUA.connectivity);
                %EleEdges=con(:,[1 2 3 1]); % These are element edges, this is only for 3-node elements.

            otherwise

                error('not implemented')
        end

        nTri=size(con,1);
        P2=NaN(5*nTri,2);
        ii=1;
        for I=1:nTri

            P2(ii:ii+3,:)=coo(EleEdges(I,:),:);
            P2(ii+4,:)=[NaN NaN ];
            ii=ii+5 ;

        end
        P2(end,:)=[] ;

        % Now calculate the intersection points using this matlab function
        [xcEdges,ycEdges]=polyxpoly(P1(:,1),P1(:,2),P2(:,1),P2(:,2));

    case "InputPoints"

        % Here I simply take the input calving profile.
        % But this may above have been resampled already, in which case this
        % is potentially a finer subdivision

        xcEdges=P1(:,1); ycEdges=P1(:,2);

end



LSF=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,LSF,xcEdges,ycEdges);

CtrlVar.PlotGLs=0;  % to suppress the plot


end