%% CALCULATES VERTEX EVOLUTION FOR GIVEN PARAMETERS
% is called by from the zipperModel class
% please note that this software was intended for internal use, and is
% therefore not in quality otherwise expected. Many functionalities are
% rather experimental (like mobile point OX or V1) and could be quite
% tricky to correctly use.

function V2v = modelDynamics(obj, varargin)
%% Input parser
% obj: zipperModel class object 
% Equilibrate: boolean, equilibration or comp. experiemnt
% StaticFroce: boolean, is static external force applied
% StaticVel: boolean, prescribe static velocity of point O
% StaticVelMagnitude: magnitude of prescribed static velocity of point O
% StaticVelTime: period for which static vel. is prescribed

inp = inputParser;
defaultEquilibrium = false;
defaultStaticOForce = true;     
isZip = @(object) isa(object,'zipperModel');

addRequired(inp, 'obj', isZip);
addParameter(inp, 'Equilibrate',defaultEquilibrium,@islogical);
addParameter(inp, 'StaticForce',defaultStaticOForce,@islogical);

parse(inp, obj, varargin{:});

Equilibrate = inp.Results.Equilibrate;
StaticForce = inp.Results.StaticForce;

%% ======================================================================
% sets initial mobile points; set static force vector if necessary
O = obj.V2X + (obj.OX-obj.V2X) * obj.V2OXfrac;
V1 = obj.V1X;
V2 = obj.V2X;
if StaticForce;     % if static external force is applied
    OXO = O - obj.OX;
    tOXO = OXO/norm(OXO);
    nOXOstatic = [-tOXO(2), tOXO(1)];
else                % if none applied (explicit to avoid numerical issues)
    nOXOstatic = [0,0];
end
eps = obj.epsilon;      % num. difference thresh to consider equal
refVecs = [1,0;0,1];    % reference friction vectors (i.e. where [1,0] and [0,1] are mapped by the friction vector)

% --- calculate initial lengths for Hookes law ---
lV2Oi  = norm(obj.V2X-O);
lV2V1i = norm(obj.V2X-obj.V1X);
lV2X3i = norm(obj.V2X-obj.X3);
lV1X2i = norm(obj.V1X-obj.X2);
lV1X1i = norm(obj.V1X-obj.X1);
lOXOi  = norm(O-obj.OX);

lA1i = lV1X2i + lV1X1i;             % initial length of axon 1
lA2i = lV2X3i + lV2V1i + lV1X2i;    % initial length of axon 2
lA3i = lOXOi + lV2Oi + lV2X3i;      % initial length of axon 3


% --- if video is on and not equilibrating, set up video writer and display progressbar ---
if ~Equilibrate
    hwb = waitbar(0,'1', 'Name', 'Dynamics running', 'CreateCancelBtn', 'setappdata(gcbf,''cancelling'',1)');
    setappdata(hwb,'cancelling',0);
    if obj.movieSwitch;
        animationWrite = VideoWriter(obj.moviepath);    % uses default 'Motion JPEG AVI' video profile
        animationWrite.Quality = 100;   % set maximal quality
        animationWrite.FrameRate = 10;  % 10 frames per second;
        open(animationWrite);
    end
end;

% --- equilibration configuration ---
if Equilibrate;              
    tmax = obj.teq;
    toff = 0;           % no outer force
    mobileO = true;     % point O can move freely
    sigmaA2 = obj.sigmaL;   % base tension in axons
    sigmaA3 = obj.sigmaR;
    outerForce = 0;         % zeros while searching static equilibrium
    stepTension = 0;        % constant tensions
    stepTensionL = 0;
    hookeStiffness = 0;     % no Hooke
    etaP = 0.001;       % avoid singular matrix for particular...
    etaN = 0.001;       % ... calculation settings, where some values...
    etaI = 10;          % ... of friction are zero ...
    etaZ = 0;           % ... converge faster w/o friction
% --- not equilibration; recorded experiment run ---
else
    sigmaA2 = obj.sigmaL + obj.stepTensionL/(obj.tmax*obj.ramp);
    sigmaA3 = obj.sigmaR + obj.stepTension/(obj.tmax*obj.ramp);     % initial tension in the axon 3
    obj.V1trajectory(1,:) = obj.V1X;    % initial positions
    obj.V2trajectory(1,:) = obj.V2X;
    obj.Otrajectory(1,:)  = O;
    tmax = obj.tmax;
    toff = obj.toff;    % when external force switches off
    mobileO = (obj.mobileO || (obj.outerForce~=0)); % point O mobile if set to be or if under external force
    outerForce = obj.outerForce;    % outer force magnitude
    stepTension = obj.stepTension;  % change in tension in the axon 3 ...
    stepTensionL = obj.stepTensionL;% ...and axon 2
    hookeStiffness = obj.hookeStiffness;    % stiffness (0 in our experiemnts)
    etaP = obj.axialFriction;       % parallel (axial) substrate friction
    etaN = obj.normalFriction;      % normal substrate friction
    etaI = obj.internalFriction;    % interval axial axon friction (due stretchnig)
    etaZ = obj.zipperFriction;      % dissipative processes due zipper advance
end

%% Dynamics calculation
% note that some experimental, unfinished parts, were removed before before
% the submission for sake of code clarity
for t=2:tmax
    
    % --- progress bar update or cancel the task ---
    if ~Equilibrate;
        if getappdata(hwb,'cancelling'); 
            delete(hwb);
            return; 
        end;
        if mod(t,obj.tdisc)==0;
            waitbar((t*obj.dt)/obj.tsec, hwb, strjoin({'Timelapse:', strcat(num2str(t*obj.dt), '/', num2str(obj.tsec))}));
        end;
    end;
    
    % --- recalculate distances and vectors for this timestep ---   
    OXO = O - obj.OX;
    V2O = O - V2;
    V2OX = obj.OX - V2;
    lOXO = norm(OXO);
    lV2O = norm(V2O);
    lV2OX = norm(V2OX);
    tOXO = OXO/lOXO;
    nOXO = [-tOXO(2), tOXO(1)];
    tV2O = V2O/lV2O;           
    nV2O = [tV2O(2),-tV2O(1)];
    %tV2OX = V2OX/lV2OX;
    
    % --- gradual tension modification ---
    if (t <= tmax*obj.ramp); 
        sigmaA2 = obj.sigmaL + stepTensionL * t/(tmax*obj.ramp); 
        sigmaA3 = obj.sigmaR + stepTension * t/(tmax*obj.ramp); 
    end;
    
    % --- switch forces off when toff is reached ---    
    if t < toff && abs(outerForce) > eps;
        Fo = outerForce * nOXOstatic;                      % force acting on point O
        Fox = - sigmaA3 * tOXO + outerForce * nOXOstatic;  % force acting on point OX
    else            
        Fo = [0,0];
        Fox = - sigmaA3 * tOXO;      % only basal tension exerted by the remaining network
    end;
    
    % --- calculate points distances ---
    lX1V1 = norm(obj.X1-V1);
    lV1X2 = norm(obj.X2-V1);
    lV1V2 = norm(V1-V2);
    lV2X3 = norm(V2-obj.X3);
    lA3 = lOXO + lV2O + lV2X3;  % length of the axon 3, at the moment
    lA2 = lV2X3 + lV1V2 + lV1X2;% length of the axon 2, at the moment
    lA1 = lX1V1 + lV1X2;        % length of the axon 1, at the moment

    %% Equations of motion for points V1,V2,O,OX (usually only V2 is mobile)
    % --- calculate velocities from geometry for O - mobile point ---
    % NOTE: point O was never independently manipulated in the presented
    % experiments.
    if mobileO;
        
        % --- Force on the loose end, total, outer force + tension + Hooke
        forceO = transpose( (Fo - (tV2O+tOXO)*hookeStiffness*(lA3-lA3i) - (tV2O + tOXO) * sigmaA3) );

        if norm(forceO) > eps;  % force is not negligible

            if norm(tOXO + tV2O) < eps; % axon is nearly straight; only normal friction

               Ov = 2*dot(forceO,(nOXO+nV2O))/(etaN*(lOXO+lV2O)) * (nOXO+nV2O);

            else

                % --- simplified friction matrix for the mobile point O

                Qv2ox = [ tOXO + tV2O; nV2O; nOXO ];            % projectors
                Etav2ox = 0.25 * [ etaP*(((lV2X3+lV2O)^3+lOXO^3)/lA3^2) + etaI/lA3 , 0, 0;  % friction coefficients
                                    0, etaN * lA3, 0
                                    0, 0, etaN * lA3 ];

                OMat = Qv2ox' * Etav2ox * Qv2ox;

                % --- O velocity
                Ov = transpose( OMat \ forceO );

            end

              %  Ov = Ov - (dot(Ov,tV2O)*tV2O);      % to forbid point sliding towards the vertex

            if t < toff && (norm(Fo) > eps)
                Ov = (Ov*nOXOstatic') * nOXOstatic; % perpendicular movement restriction during applied force or deformation
            end

        else    % outer force < eps
            Ov = [0,0];
        end;
        
    else % mobileO = false
        Ov = [0,0];     % calculate nothing, later update position so it remains on the V2OX line;
                        % not blocking the O movement can for low friction
                        % parameters result in its movement onto the V2 and
                        % result in a chaotic tV2O direction vector
    end

    if (~Equilibrate && mod(t,obj.tdisc) == 0); obj.Ovel(t*obj.dt,:) = Ov; end;  % keeping track of velocity of O in time
    
    
    % --- OX velocity, the end point ---
    % NOTE: never mobile in presented experiments
    if obj.mobileOX;
       
        % --- Total force on the end ---
        forceOX = transpose( (Fox + tOXO*hookeStiffness*(lV2O-lV2Oi+lV2X3-lV2X3i+lOXO-lOXOi)) + tOXO * sigmaA3 );
        
        if norm(forceOX) > eps

            % normal friction depends on point O mobility
            if mobileO; lV2X = lOXO;
            else lV2X = lV2OX;
            end
            
            % --- Friction matrix for the loose end ---
            Qoxo = [ tOXO; nOXO ];            % projector
            Etaoxo = 0.25 * [ etaP*lA3 + etaI/lA3 , 0; 0, etaN * lV2X ];
            OXMat = Qoxo' * Etaoxo * Qoxo;
                        
            % --- OX velocity
            OXv = transpose( OXMat \ forceOX );

        else % mobileOX = false
            OXv = [0,0];
        end;

    end;
                
    
    % --- V2 VERTEX ---
    % equations to calculate velocity of the vertex V2 under substrate
    % friction, internal dissipation and dissipation in the vertex
    V1V2 = V1-V2;
    V2X3 = obj.X3-V2;
    tV1V2 = V1V2/lV1V2;
    tV2X3 = V2X3/lV2X3;
    nV2X3 = [ -tV2X3(2), tV2X3(1) ];
    nV1V2 = [ tV1V2(2), -tV1V2(1) ];
    forceV2 = (tV1V2 + tV2X3) * sigmaA2 + (tV2O + tV2X3) * sigmaA3;  % tensile force vector on V2
    forceV2 = forceV2 + hookeStiffness * (tV1V2*(lA2 - lA2i) + tV2X3*(lA2 - lA2i)+tV2X3*(lA3-lA3i)+ tV2O*(lA3-lA3i));  % Hooke's law tension increase
    forceV2 = forceV2 - obj.adh * tV2X3;               % overall force on V2
    
    % store information about tension and V2 zipper angle
    if (~Equilibrate && mod(t,100) == 0);   % save information
        obj.A3tension(t/100) = sigmaA3 + hookeStiffness * (lV2X3-lV2X3i+lV2O-lV2Oi+lOXO-lOXOi);
        obj.A2tension(t/100) = sigmaA2 + hookeStiffness * (lV1V2-lV2V1i+lV2X3-lV2X3i+lV1X2-lV1X2i);
        obj.angle(t/100) = real( acos(dot(tV2O,-tV2X3))/pi*180 + ...    % right semiangle
                            acos(dot(tV1V2,-tV2X3))/pi*180 );           % left semiangle
    end;
    
    if norm(forceV2) > eps; % non-negligible force
        
        if ~mobileO; 
            lV2OOX = lV2OX; % if O is rigid, use the whole branch to calculate the friction
        else
            lV2OOX = lV2O + lOXO;
        end;    
        
        % --- Complete friction with transformation ---

        QA3 = [ tV2O + tV2X3; nV2O; nV2X3; tV2X3 ];       % projector to the axon 3 and its normals
        EtaA3 = 1/6*[ etaP * ( (lV2X3^3+lV2OOX^3)/lA3^2) + etaI * (3/lA3), 0, 0, 0;...
            0, etaN * lV2O, 0, 0; 0, 0 , etaN * lV2X3, 0;...       % friction coefficients for V2OX & V2X3
            0, 0, 0, 3*etaZ ];
            
        QA2 = [ tV1V2 + tV2X3; nV1V2; nV2X3; tV2X3 ];     % projector to the axon 2 and its normals
        EtaA2 = 1/6*[ etaP * ( (lV2X3^3+lV1V2^3)/lA2^2) + etaI * (3/lA2), 0, 0, 0;
            0, etaN * lV1V2, 0, 0; 0, 0, etaN * lV2X3, 0;...        % friction coefficients for V2X3 & V1V2
            0, 0, 0, 3*etaZ ];
            
        fV2Mat = QA3' * EtaA3 * QA3 + QA2' * EtaA2 * QA2;   % complete friction matrix
             
        % --- Velocity by matrix inversion ---
        V2v = transpose( fV2Mat \ transpose(forceV2) );
        
    else % force on V2 is negligible//this approach avoids numerical instability
        V2v = [0,0];
        fV2Mat = [1,0;0,1];
    end
    
    if (~Equilibrate && mod(t,100) == 0); obj.V2vel(t/100,:) = V2v; end;  % keeping track of velocity in time
    
    % --- V1 VERTEX ---    
    % NOTE: never mobile in presented experiments
    if obj.mobileV1;
    
        V1X2 = obj.X2-V1;
        V1X1 = obj.X1-V1;
        V2V1 = -V1V2;
        tV1X2 = V1X2/lV1X2; 
        tV1X1 = V1X1/lX1V1;
        tV2V1 = V2V1/lV1V2;
        nV1X1 = [ tV1X1(2), -tV1X1(1) ];
        nV1X2 = [ tV1X2(2), -tV1X2(1) ];
        forceV1 = (tV2V1 + tV1X2) * sigmaA2 + (tV1X2 + tV1X1) * obj.sigma0;   % tensile force vector on V1
        forceV1 = forceV1 + hookeStiffness * (tV2V1*(lA2-lA2i) + tV1X2*(lA2-lA2i)+ tV1X2*(lA1-lA1i)+ tV1X1*(lA1-lA1i));
        forceV1 = forceV1 - obj.adh * tV1X2;                % overall force on V1

        if (~Equilibrate && mod(t,obj.tdisc) == 0);
            obj.A1tension(t/obj.tdisc) = obj.sigma0 + hookeStiffness * (lV1X2-lV1X2i+lX1V1-lV1X1i);
        end

        if norm(forceV1) > eps;

            
            QA1 = [ tV1X1 + tV1X2; nV1X1; nV1X2; tV1X2 ];       % projector to the axon 1 and its normals
            EtaA1 = 1/6*[ etaP * ( (lV1X2^3+lX1V1^3)/lA1^2) + etaI * (3/lA1), 0, 0, 0;...
                0, etaN * lX1V1, 0, 0; 0, 0 , etaN * lV1X2, 0;...       % friction coefficients for V1X1 & V1X2
                0, 0, 0, 3*etaZ ];

            QA2b = [ tV2V1 + tV1X2; -nV1V2; nV1X2; tV1X2 ];     % projector to the axon 2 and its normals
            EtaA2b = 1/6*[ etaP * ( (lV1X2^3+lV1V2^3)/lA2^2) + etaI * (3/lA2), 0, 0, 0;
                0, etaN * lV1V2, 0, 0; 0, 0, etaN * lV1X2, 0;...        % friction coefficients for V2X3 & V1V2
                0, 0, 0, 3*etaZ ];

            % --- Friction force matrix
            fV1Mat = QA1' * EtaA1 * QA1 + QA2b' * EtaA2b * QA2b;
            
            V1v = transpose( fV1Mat \ transpose(forceV1) );

        else
            V1v = [0,0];
        end

        if (~Equilibrate && mod(t,100) == 0); obj.V1vel(t/100,:) = V1v; end;  % record V1 velocity

    else
        V1v = [0,0];
    end;
        
    %% Update position of mobile points
    if obj.mobileOX; obj.OX = obj.OX + OXv*obj.dt; end
    if obj.mobileV1; V1 = V1 + V1v*obj.dt; end
    if mobileO; O = O + Ov*obj.dt;              % make sure O doesn't behave improperly for low frictions
    else O = V2 + (obj.OX - V2)*obj.V2OXfrac;   % to suppress the sliding, if the O motion is not necessary, fix it on the V2OX line
    end; 

    V2 = V2 + V2v*obj.dt;
    
    
    %% Record currect system information (positions, forces, etc.)
    if ~Equilibrate && mod(t,100) == 0; 
        [eVecs,eVals] = eig(fV2Mat);    % find eigenvectors and eigenvalues of V2 friction tensor
        % verify the eigenvalues and eigenvectors do not swap
        if (abs(dot(eVecs(:,1),refVecs(:,1))) < abs(dot(eVecs(:,1),refVecs(:,2))) && ...
            abs(dot(eVecs(:,2),refVecs(:,2))) < abs(dot(eVecs(:,2),refVecs(:,1))))
            temp = eVals(1,1);
            eVals(1,1) = eVals(2,2);
            eVals(2,2) = temp;
            refVecs(:,1) = eVecs(:,2);
            refVecs(:,2) = eVecs(:,1);
        else
            refVecs = eVecs;
        end;
        obj.FVecs(t/100,:) = [refVecs(:,1)',refVecs(:,2)'];
        obj.V1trajectory(t/100,:) = V1;
        obj.V2trajectory(t/100,:) = V2;
        obj.Otrajectory(t/100,:) = O;
        obj.Fratio(t/100,:) = [eVals(1,1), eVals(2,2)];    % ratio of V2 friction tensor eigenvalues
        obj.V2length(t/100,:) = lV2X3;
        obj.A2length(t/100,:) = lA2;
        obj.A3length(t/100,:) = lA3;
        if(obj.movieSwitch && ~Equilibrate); 
            writeVideo(animationWrite, vertexAnimator(V1, V2, obj.X3, obj.OX, O, Fo,Fox, V2v, (tV1V2)*obj.A2tension(t/100), (tV2O)*obj.A3tension(t/100), forceV2, fV2Mat, diag(eVals), obj.moviesetting )); 
        end
    end;

end
    
    % update equilibrated vertex positions after equilibrium run
    if Equilibrate; 
        obj.V1X = V1;
        obj.V2X = V2;
    end

    if ~Equilibrate; 
        delete(hwb);    % delete waitbar handle
        if obj.movieSwitch; close(animationWrite); end;   % close video saving object
    end
end
