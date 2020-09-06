%% Introduction and geometry of the studied system
% This code can be used to perform simulations of zippering under various
% external influences (i.e. increase in tension, external force etc.) and
% various friction/dissipation strength. The code was not intended for
% external use, and therefore does not include extensive help or user input
% verification. The model data and methods are contained in the class
% zipperModel, while the dynamics calculation is performed by an external
% function modelDynamics. The class can be controlled using the GUI by
% running function zipperModelGUI.
% Following configuration of the network is simulated. The point V1 is
% usually held fixed and only the dynamics of the vertex V2 is calculated

%     X2       
%     ||      
%     ||      
%     V1      OX
%    /  \    /  
%   /    \  O    
%  X1     \/      
%         V2
%         ||       
%         ||
%         X3

% X - fixed end, O - possibly pulled point, OX - fixed end

%% Class begins here
% The class covers the data and functions of the model of vertex dynamics;
% It calls 'modelDynamics' function to calculate time development, and
% itself contains auxilliary functions
classdef zipperModel < handle

    %% Variables
    
    % these are static properties for all instances of the class
    % evaluated only once
    properties (Constant = true)
        tdisc = 100;        % time discretisation; steps per second
        epsilon = 0.000001; % precision for two floats to be considered equal
    end
    
    % calculated from the constant values by constructor;
    properties (Access = public)
       name;
       dt;              % inverse of tdisc
       teq = 50;        % time for equilibrium in seconds; constructor converts to iterations
    end
    
    % geometry and constant parameters
    properties (Access = public)
       % initial geometry [in um] as illustrated on top
       X1 = [3,-1];        
       X2 = [10.9,25];
       X3 = [17.3089,-25];
       V1X = [10.3089,20];             % initial position of vertex 1
       V2X = [17.3089,5.5468];         % initial position of vertex 2
       OX = [24.3089,20];
       V2OXfrac = 0.5;                 % default value for the positioning of point O, as fraction of OX-V2 distance
        
       % basic fixations of vertex V1, point O (where axon can be
       % potentially bend) and even remote fixation point OX; In the
       % model output presented, none of these is mobile
       mobileV1 = false;
       mobileOX = false;
       mobileO = false;                
       % i.e. point O moves even when false, but remains on the V2OX line
       % to perform lateral pulling or deformation experiment, this
       % parameter must be set to true!
        
       % biophysical parameters
       adh = 0.2;          % adhesive parameter [nN]
       sigma0 = 1.0;       % basal tension in axons [nN]
       sigmaL = 1.0;       % basal tension in axon A2
       sigmaR = 1.0;       % basal tension in axon A3
    end
    
    % calculation related parameters, default initial values
    properties (Access = public)
       
        outerForce = 0;     % magnitude of outer force (if any)
        stepTension = 0;    % increase in tension (right axon)
        stepTensionL = 0;   % increase in tension (left axon)
        hookeStiffness = 0; % axon's Hookean stiffness
        normalFriction = 0.2;   % segment normal substrate friction coef
        axialFriction = 0.2;    % segment parallel substrate friction coef
        internalFriction = 3000;% internal friction within the axon when stretching
        zipperFriction = 1.5;
        tmax = 100000;          % number of iterations ...
        tsec;                   % ... in seconds
        tfrac = 0.5;
        toff;                   % time to switch off force
        ramp = 0.5;             % ramp of increasing tension
        movieSwitch = false;    % switch to generate animation or not (time demanding)
        moviepath;              % animation output file full name
        moviesetting = struct('frame', false, 'tension', false, 'force', false,...
                              'velocity', false, 'tfactor', 1, 'ffactor', 1,...
                              'vfactor', 1, 'ftensor',false, 'frfactor', 7,...
                              'mtensor', false, 'mfactor', 7);
        % settings of the output movie in the struct order to show: axon shafts;
        % tension vectors; force vec; velocity vec; vector-length-scaling factors for
        % tension, force and velocity; friction tenson (as an ellipse);
        % scale fric. ten.; mobility tensor (inverse of fric.ten.); scale
        % mobility tensor        
    end
    
    % arrays to make record every second
    properties (Access = public)   
        
        V1trajectory;   % trajectory of V1
        V2trajectory;   % ... V2
        Otrajectory;    % ... point O
        FVecs;          % friction eigenvectors
        Fratio;         % friction ratio of eigenvalues
        V1vel;          % velocity of V1
        V2vel;          % ... V2
        Ovel;           % ... point O
        angle;          % zipper angle
        A1tension;      % tension in axon 1 (tension uniform in axon)
        A2tension;      % ... 2
        A3tension;      % ... 3
        convergence;    % distance to final equilibrium
        V2length;       % length of zipper at V2, i.e. V2X3
        A2length;       % length of middle axon
        A3length;       % length of the rightmost axon
        Energyscape;    % energy landscape of conservative forces (elastic)
        energyscapeAnimation = struct('cdata',[],'colormap',[]);
        energyscapeSlowTrail;
        
    end
        
    % counts number of already assigned instances
    methods(Static)
        
        function out = Counter(cmd)
            
            persistent counter;
            
            if nargin==1; 
                switch cmd
                    case 'add'
                        if numel(counter)
                            counter = counter + 1;
                        else
                            counter = 1;
                        end
                    case 'reset'
                        counter = 0;
                end
            end
            
            out = counter;
        end
    end
                
    
    %% Methods
    methods
        % constructor
        function obj = zipperModel(name_)
            if nargin == 0;
                counter = obj.Counter('add');
                obj.name = strcat('experiment',num2str(counter));
            elseif nargin == 1;
                obj.name = name_;
            end
            
            obj.dt = 1/obj.tdisc;            % set discretization
            obj.teq = obj.teq * obj.tdisc;   % default period of equilibration
            obj.setTime(obj.tmax);           % set times according to tmax
            obj.moviepath = strcat(pwd,'/',obj.name,'.avi');    % set path to save movie
            
        end
    
        % updates tensions so that the input coordinates are a static
        % configuration
        function positionStaticVertex(obj, coor_)
            
            obj.V2X = coor_;
            
            vC = obj.V2X - obj.X3;  % common branch
            vL = obj.V1X - obj.V2X; % left branch
            vR = obj.OX  - obj.V2X; % right branch
            
            aL = acos( dot(vL,vC)/norm(vL)/norm(vC) );
            aR = acos( dot(vR,vC)/norm(vR)/norm(vC) );
            
            geoMat = [ sin(aL) , -sin(aR) ; 1-cos(aL) , 1-cos(aR) ];
            RS = [ 0; obj.adh ];
            
            sigma = geoMat\RS;
            
            obj.sigmaL = sigma(1);
            obj.sigmaR = sigma(2);
            
        end
        
        % finds static coordinates for input tensions
        function adjustStaticVertex(obj, tension_)
           
            obj.sigmaL = tension_(1);
            obj.sigmaR = tension_(2);
            
            finalVel = modelDynamics(obj, 'Equilibrate', true );   % callibrates vertex positions
            
            while( norm(finalVel) > 1e-4 )
                disp(strcat('Final vertex velocity, v=',num2str(norm(finalVel)),...
                    ',exceeds the threshold 1e-4. Repeating equilibration.'));
                finalVel = modelDynamics(obj, 'Equilibrate', true );
            end;
            
            disp('Vertex equilibration finished');
            
        end
        
        % resets number of iterations; reinitialize all arrays
        function setTime(obj, time_)
            
            obj.tmax = time_;
            obj.setToff(obj.tfrac);
            obj.tsec = round(obj.tmax*obj.dt);
            
            % reinitialize the containers
            obj.V1trajectory = zeros(obj.tsec,2);
            obj.V2trajectory = zeros(obj.tsec,2);
            obj.Otrajectory = zeros(obj.tsec,2);
            obj.FVecs = zeros(obj.tsec, 4);
            obj.Fratio = zeros(obj.tsec,2);
            obj.V1vel = zeros(obj.tsec,2);
            obj.V2vel = zeros(obj.tsec,2);
            obj.Ovel = zeros(obj.tsec,2);
            obj.angle = zeros(obj.tsec,1);
            obj.A1tension = zeros(obj.tsec,1);
            obj.A2tension = zeros(obj.tsec,1);
            obj.A3tension = zeros(obj.tsec,1);
            obj.convergence = zeros(obj.tsec,1);
            obj.V2length = zeros(obj.tsec,1);
            obj.A2length = zeros(obj.tsec,1);
            obj.A3length = zeros(obj.tsec,1);
        end
        
        % sets time parameter to switch off outer force
        function setToff(obj, tFrac)
            obj.tfrac = tFrac;
            obj.toff = obj.tmax * obj.tfrac;
        end
        
        % runs the model dynamics for the given parameters            
        function runDynamics(obj)
                        
            modelDynamics(obj, 'Equilibrate', true );   % equilibrates vertex positions
                        
            modelDynamics(obj, 'Equilibrate', false, 'StaticForce', ~(obj.outerForce==0) );   % runs simulation with perturbations
            
            obj.getConvergence();   % calculates convergence
            
        end
    
        % calculate the opproach to the final point
        function getConvergence(obj)            
            for i=1:obj.tsec
                obj.convergence(i) = norm(obj.V2trajectory(obj.tsec,:) - obj.V2trajectory(i,:));
            end
        end
        
        %% Support methods
        % calculate the etaP, etaN, etaI and etaZ to get specific eigenvalues
        function setFrictions(obj,ratio)
            
            % initial values
            start = [100,1,1];
            % model handle
            model = @getEvals;
            
            est = fminsearch(model,start);
            disp(strcat('Optimal values are:',num2str(est)));

            % calculates reatio of friction eigenvalues for given param
            function [err] = getEvals(eta)
                
                    etaI = eta(1);
                    etaP = eta(2);
                    etaN = eta(3);
                    etaZ = eta(4);

                    V2OX  = obj.OX - obj.V2X;
                    lV2OX = norm(V2OX);
                    tV2OX = V2OX/lV2OX;
                    nV2OX = [ tV2OX(2), -tV2OX(1)];
                    V1V2  = obj.V1X - obj.V2X;
                    lV1V2 = norm(V1V2);
                    tV1V2 = V1V2/lV1V2; 
                    nV1V2 = [ tV1V2(2), -tV1V2(1) ];
                    V2X3  = obj.X3 - obj.V2X;
                    lV2X3 = norm(V2X3);
                    tV2X3 = V2X3/lV2X3;
                    nV2X3 = [ -tV2X3(2), tV2X3(1) ];        
                    lA3 = lV2OX + lV2X3;        % length of the axon 3, at the moment
                    lA2 = lV2X3 + lV1V2;% + lV1X2;% length of the axon 2, at the moment

                    % calculate the friction tensor

                    QA3 = [ tV2OX + tV2X3; nV2OX; nV2X3; tV2X3 ];       % projector to the axon 3 and its normals
                    EtaA3 = 1/6*[ etaP * ( (lV2X3^3+lV2OX^3)/lA3^2) + etaI * (3/lA3), 0, 0, 0;...
                                    0, etaN * lV2OX, 0, 0; 0, 0 , etaN * lV2X3, 0;...       % friction coefficients for V2OX & V2X3
                                    0, 0, 0, 3*etaZ ];
                    QA2 = [ tV1V2 + tV2X3; nV1V2; nV2X3; tV2X3 ];       % projector to the axon 2 and its normals
                    EtaA2 = 1/6*[ etaP * ( (lV2X3^3+lV1V2^3)/lA2^2) + etaI * (3/lA2), 0, 0, 0;
                                    0, etaN * lV1V2, 0, 0; 0, 0, etaN * lV2X3, 0;...        % friction coefficients for V2X3 & V1V2
                                    0, 0, 0, 3*etaZ ];
                    fV2Mat = QA3' * EtaA3 * QA3 + QA2' * EtaA2 * QA2;
                    [~,eVals] = eigs(fV2Mat);    % calculate eigenvectors and eigenvalues
                    err = abs(eVals(1,1)/eVals(2,2) - ratio);
                    if (eta(1)<0)||(eta(2)<0)||(eta(3)<0)||(eta(4)<0) 
                        err = 10e10;
                    end
            end
        end        
        
        % simple method to visualize energy landscape
        function frameOut = getEnergyLandscape(obj,sR,sL,varargin)
            
           inp = inputParser;
           defaultScale = 10;    % level of plot detail
           defaultCentered = false;
           defaultSingle = false;
           isZip = @(object) isa(object,'zipperModel');
           
           addRequired(inp,'obj',isZip);
           addRequired(inp,'sR',@isnumeric);
           addRequired(inp,'sL',@isnumeric);
           addParameter(inp,'scale', defaultScale, @isnumeric);
           addParameter(inp,'centered', defaultCentered, @islogical);
           addParameter(inp,'single',defaultSingle, @islogical);
           
           parse(inp, obj, sR, sL, varargin{:});
                      
           scale = inp.Results.scale;       % discretization; pieces per micron
           centered = inp.Results.centered; % zooms fields on the central part
           single = inp.Results.single;     % if generating single standalone, or animation
           
           % ===========================================================
           % discretization of the field
           if ~centered
               xrange = (round(min([obj.V1X(1), obj.X3(1), obj.OX(1)])-4))...
                        :1/scale:...
                        (round(max([obj.X3(1), obj.V1X(1), obj.OX(1)])+4));        
               yrange = (round(min([obj.V1X(2), obj.X3(2), obj.OX(2)])-4))...
                        :1/scale:...
                        (round(max([obj.V1X(2), obj.X3(2), obj.OX(2)])+4));
           else
               xrange =(obj.V2X(1) - 0.5):1/scale:(obj.V2X(1) + 2);
               yrange = (obj.V2X(2) - 4):1/scale:(obj.V2X(2) + 1);
           end
                    
           % calculate the energy matrix
           E = zeros(numel(xrange), numel(yrange));
           
           for x = 1:numel(xrange)
               for y = 1:numel(yrange)
                   C = [xrange(x),yrange(y)];
                   E(x,y) =  scale * (...
                                    sR * ( norm(obj.OX-C) + norm(obj.X3-C) ) + ...
                                    sL * ( norm(obj.V1X-C) + norm(obj.X3-C) ) - ...
                                    obj.adh * ( norm(obj.X3 - C) ) + ...
                                    0.5 * obj.hookeStiffness * ( ...
                                    ( norm(obj.OX-C) + norm(obj.X3-C) - ...
                                    norm(obj.OX-obj.V2X) - norm(obj.X3-obj.V2X) + 1 )^2 + ...
                                    ( norm(obj.V1X-C) + norm(obj.X3-C) - ...
                                    norm(obj.V1X-obj.V2X) - norm(obj.X3-obj.V2X) + 1 )^2 ));
               end
           end
                      
           % calculate the new equilibrium
           [minVal, minInd] = min(E);
           [~, newEq(2)] = min(minVal);
           newEq(1) = minInd(newEq(2));
           newEq = [ xrange(newEq(1)), yrange(newEq(2)) ];
           
           % for a single frame, generate also the gradiant descent
           if single;
               % calculate the gradient for the landscape
               [fx,fy] = gradient(-1*E,1,1);

               % initiate steep descent path
               head(1,:) = obj.V2X;
               envir = [ 0,0; 1,0; 0,1; 1,1; -1,0; 0,-1; -1,-1; 1,-1; -1,1 ];  % environment of the trailhead

               % calculate the steepest descend path, olny if single

               while norm(head(end,:) - newEq) > 1.2/scale
                   ind = [ floor(scale*(head(end,1)-xrange(1)))+1, floor(scale*(head(end,2)-yrange(1)))+1];
                   vec = [0,0];
                   weightsum = 0;
                   for i = 1:9
                       index = ind + envir(i,:);
                       coor = [ xrange(index(1)), yrange(index(2)) ];   % coordinate of the indexed node
                       if i~=1; weight = 1/norm( head(end,:) - coor );   % weight inverse to distance
                       else     weight = 2*scale; end;
                       weightsum = weightsum +  weight;
                       vec = vec + [fy(index(1), index(2)), fx(index(1), index(2))] * weight;
                   end
                   vec = vec/weightsum;     % average the gradient
                   head(end+1,:) = head(end,:) + vec;   % new point on gradient path
               end;
               head(end+1,:) = newEq;
           end;
           
           % generate a figure
           if ~single; energyFig = figure('Visible','off');
           else energyFig = figure('Visible','on'); end;
           ax = gca;
           ax.NextPlot = 'replaceChildren';
           hold on;
           
           scatter(newEq(1), newEq(2), 'blue', 'filled');       % plot new equilibrium 
           contour(xrange, yrange, E', 100);                    % plot the energy landscape
           
           if ~centered     % plot fixation point only in non-centered case
                scatter([obj.V1X(1), obj.X3(1), obj.OX(1)], [obj.V1X(2), obj.X3(2), obj.OX(2)],...
                    'k','d','filled');
           end;

           if single;   % only for single image
                scatter(obj.V2X(1), obj.V2X(2),'r','filled');        % plot initial equilibrium
                line(head(:,1), head(:,2));                          % plot the gradient path
           else
                obj.energyscapeSlowTrail(end+1,:) = newEq;
                line(obj.energyscapeSlowTrail(:,1),obj.energyscapeSlowTrail(:,2));
           end;
           
           daspect([1,1,1]);    % 1:1 aspect ratio
           
           % return the frame
           frameOut = getframe(energyFig);
           
           % save as object data to make accessible from base workspace for
           % a single snapshot; not much sense for video
           if single;
               obj.Energyscape.xrange = xrange;
               obj.Energyscape.yrange = yrange;
               obj.Energyscape.energy = E';
               obj.Energyscape.trail  = head;
           end;
        end        
        
        % calculates friction tensor at each point of the canvas
        function [bigquiver,smallquiver,bigcont,smallcont]=getFrictionField(obj)
            
            scale = 2;
            centered = false;
            %lV1X2 = norm(obj.X2 - obj.V1X);
            etaP = obj.axialFriction;
            etaN = obj.normalFriction;
            etaI = obj.internalFriction;
            etaZ = obj.zipperFriction;
            %fvec = [1;-1];
            
%             if ~centered
%                 xrange = (floor(min([obj.V1X(1), obj.X3(1), obj.OX(1)]))-1)...
%                          :1/scale:...
%                          (ceil(max([obj.X3(1), obj.V1X(1), obj.OX(1)]))+1);        
%                 yrange = (round(min([obj.V1X(2), obj.X3(2), obj.OX(2)])))...
%                          :1/scale:...
%                          (round(max([obj.V1X(2), obj.X3(2), obj.OX(2)])));
%             else
%                 xrange =(obj.V2X(1) - 0.5):1/scale:(obj.V2X(1) + 2);
%                 yrange = (obj.V2X(2) - 4):1/scale:(obj.V2X(2) + 1);
%             end
            
            xrange= 10:1/scale:25;
            yrange= -25:1/scale:20;
                    
            % calculate the energy matrix
            fTensor = zeros(numel(xrange), numel(yrange),6);
            
            for x = 1:numel(xrange)
                for y = 1:numel(yrange)
                    C = [xrange(x),yrange(y)];
                    [fTensor(x,y,1:4),fTensor(x,y,5:6)] =  getTensor(C);
                end
            end
            
            figure
            subplot(1,2,1);
            hold on;
            scatter(obj.V2X(1), obj.V2X(2),'r','filled');
            scatter(obj.V1X(1), obj.V1X(2),'k','filled');
            scatter(obj.OX(1), obj.OX(2),'k','filled');
            scatter(obj.X3(1), obj.X3(2),'k','filled');
            quiver(xrange,yrange, fTensor(:,:,1)', fTensor(:,:,2)');
            xind=1:4:numel(xrange);
            yind=1:4:numel(yrange);
            bigquiver=zeros(numel(xind)*numel(yind),5);
            smallquiver=zeros(numel(xind)*numel(yind),5);
            i=1;
            for xc=xind
                for yc=yind
                    coor=[xrange(xc),yrange(yc)];
                    bigdirs=[fTensor(xc,yc,1),fTensor(xc,yc,2)]*sign(fTensor(xc,yc,1));
                    smalldirs=[fTensor(xc,yc,3),fTensor(xc,yc,4)]*sign(fTensor(xc,yc,4));
                    bigquiver(i,:)=[coor,bigdirs,fTensor(xc,yc,5)];
                    smallquiver(i,:)=[coor,smalldirs,fTensor(xc,yc,6)];
                    i=i+1;
                end
            end            
            contour(xrange, yrange, real(fTensor(:,:,5)'),'ShowText','on');   % plot the frinction contourplot
            cdata=contour(xrange,yrange,real(fTensor(:,:,5)'));
            bigcont=cdata';
            daspect([1 1 1]);
            
            subplot(1,2,2)
            hold on;            
            scatter(obj.V2X(1), obj.V2X(2),'r','filled');
            scatter(obj.V1X(1), obj.V1X(2),'k','filled');
            scatter(obj.OX(1), obj.OX(2),'k','filled');
            scatter(obj.X3(1), obj.X3(2),'k','filled');
            quiver(xrange,yrange, fTensor(:,:,3)', fTensor(:,:,4)');
            contour(xrange, yrange, real(fTensor(:,:,6)'),'ShowText','on');
            cdata=contour(xrange,yrange,real(fTensor(:,:,6)'));
            smallcont=cdata';
            daspect([1 1 1]);
            
            % calculates friction tensor for each coordinate point, returns
            % the largest eigenvalue and corresponding eigenvector
            % (normalized)            
            function [vec,vals] = getTensor(V2)
                
                V2OX  = obj.OX - V2;
                lV2OX = norm(V2OX);
                tV2OX = V2OX/lV2OX;
                nV2OX = [ tV2OX(2), -tV2OX(1)];
                V1V2  = obj.V1X - V2;
                lV1V2 = norm(V1V2);
                tV1V2 = V1V2/lV1V2; 
                nV1V2 = [ tV1V2(2), -tV1V2(1) ];
                V2X3  = obj.X3 - V2;
                lV2X3 = norm(V2X3);
                tV2X3 = V2X3/lV2X3;
                nV2X3 = [ -tV2X3(2), tV2X3(1) ];        
                lA3 = lV2OX + lV2X3;            % length of the axon 3, at the moment
                lA2 = lV2X3 + lV1V2;% + lV1X2;  % length of the axon 2, at the moment

                % calculate the friction tensor

                QA3 = [ tV2OX + tV2X3; nV2OX; nV2X3; tV2X3 ];       % projector to the axon 3 and its normals
                EtaA3 = 1/6*[ etaP * ( (lV2X3^3+lV2OX^3)/lA3^2) + etaI * (3/lA3), 0, 0, 0;...
                    0, etaN * lV2OX, 0, 0; 0, 0 , etaN * lV2X3, 0;...       % friction coefficients for V2OX & V2X3
                    0, 0, 0, 3/2*etaZ ];

                QA2 = [ tV1V2 + tV2X3; nV1V2; nV2X3; tV2X3 ];       % projector to the axon 2 and its normals
                EtaA2 = 1/6*[ etaP * ( (lV2X3^3+lV1V2^3)/lA2^2) + etaI * (3/lA2), 0, 0, 0;
                    0, etaN * lV1V2, 0, 0; 0, 0, etaN * lV2X3, 0;...        % friction coefficients for V2X3 & V1V2
                    0, 0, 0, 3/2*etaZ ];

                fV2Mat = QA3' * EtaA3 * QA3 + QA2' * EtaA2 * QA2;
                
%                fV2MatInv = inv(fV2Mat);
                
                [eVecs,eVals] = eigs(fV2Mat);    % calculate eigenvectors and eigenvalues
                vals = [eVals(1,1),eVals(2,2)];
                vec = [eVecs(:,1)',eVecs(:,2)'];
                
%                 vecs(:,1) = fV2Mat \ [1;0];  % velocity for the given force
%                 vecs(:,2) = fV2Mat \ [0;1];
%                 vals = [1/norm(vecs(:,1)),1/norm(vecs(:,2))];
%                 vecs(:,1) = vecs(:,1)/norm(vecs(:,1));
%                 vecs(:,2) = vecs(:,2)/norm(vecs(:,2));
%                 vec = [vecs(:,1)',vecs(:,2)'];% vectors for particular forces
                     
            end            
        end
        
        % estimate velocity of the vertex
        function estimateVelocity(obj)
            
            tX3V2 = (obj.V2X - obj.X3);
            tX3V2 = tX3V2/norm(tX3V2);
            tV2O  = (obj.OX - obj.V2X);
            tV2O  = tV2O/norm(tV2O);
            tV2V1 = (obj.V1X - obj.V2X);
            tV2V1 = tV2V1/norm(tV2V1);
           
            estVel = abs( obj.adh - (obj.sigmaL+obj.stepTensionL)*(1-dot(tX3V2,tV2V1))-...
                        (obj.sigmaR+obj.stepTension)*(1-dot(tX3V2,tV2O)))/...
                        ( obj.normalFriction * norm(obj.OX(1) - obj.X3(1)) * ...
                         sqrt( 1 - dot(tX3V2,tV2O)^2) );
                     
            disp(strcat('Estimated initial velocity (for rapid onset):', num2str(estVel),'um/s'));
                     
        end
        
    end
end


