function [ varargout ] = zipperModelGUI( varargin )
% zipperModelGUI simplifies the work with zipperModel class parameters
% GUI takes one optional argument
% if zipperModel object is passed, GUI connects to the object
% empty call results in creation of new zipperModel class instance (object)
% GUI returns passed object or the newly created object handle

inp = inputParser;
defaultObj = zipperModel( );
isZip = @(object) isa(object,'zipperModel');

addOptional(inp,'obj',defaultObj,isZip);
parse(inp, varargin{:} );

obj = inp.Results.obj;

if numel(inp.UsingDefaults) == 1; varargout{1} = obj; end;

%% ===================================================================
toPlot = 1;
V2x = obj.V2X(1);
V2y = obj.V2X(2);
sL  = obj.sigmaL;
sR  = obj.sigmaR;

obj.normalFriction = 0;
obj.axialFriction = 0;
obj.internalFriction = 0;
obj.zipperFriction = 0;

% GUI figure
hf = figure('Visible','off','Position',[0,0,1200,800],'Resize','off');

% PANELS
hcb = uibuttongroup('Parent', hf, 'Title','Dynamics parameters', 'Position', [0.65, 0.65, 0.3, 0.33]);
hacb = uibuttongroup('Parent', hf, 'Title', 'Animated content options', 'Position', [0.65, 0.1, 0.3, 0.2]);
hmp = uipanel('Parent', hf, 'Title', 'Animated output options', 'Position', [0.65, 0.3, 0.3, 0.15]);
hgp = uipanel('Parent', hf, 'Title', 'Graph options', 'Position', [0.65, 0.45, 0.3, 0.1]);
hbtn = uipanel('Parent', hf, 'Title', 'Dynamics time and controls', 'Position', [0.65, 0.55, 0.3, 0.1]);
hV2 = uibuttongroup('Parent', hf, 'Title', 'Vertex location', 'Position', [0.05, 0.85, 0.5, 0.15],...
                    'SelectionChangedFcn', @vertexselection);

% BUTTONS
hrun = uicontrol(hbtn,'Style','pushbutton', 'String', 'run', 'Position', [10,35,70,25], ...
                 'Callback',{@runbutton_callback});
hplot = uicontrol(hbtn, 'Style','pushbutton', 'String', 'plot', 'Position', [10,5,70,25],...
                   'Callback',{@plotbutton_callback},'Enable','off');

% STATIC TEXT
hdatatext = uicontrol(hbtn,'Style','text','String','No data available','Position',[100,5,120,15]);
hitText   = uicontrol(hbtn, 'Style', 'text', 'String', 'sec', 'Position', [220, 35, 50, 15]);    
hitVtxt = uicontrol(hacb,'Style', 'text', 'String', 'rescale vectors', 'Position', [190,90,120,20]);
hitDtxt = uicontrol(hacb,'Style', 'text', 'String', 'display', 'Position', [10,90,120,20]);
hramplbl = uicontrol(hcb, 'Style', 'text', 'String', 'Slide to select duration of linear tension increase', 'Position', [10, 190, 330, 20]);
hfriclbl = uicontrol(hcb, 'Style', 'text', 'String', 'Check type of friction to use (edit magnitude)', 'Position', [10, 145, 330, 20]);
htenlbl = uicontrol(hcb, 'Style', 'text', 'String', 'Check axon tension to change (edit magnitude)', 'Position', [10, 40, 330, 20]);

% EDITABLE TEXT FIELDS -- INPUTS
hitSL = uicontrol(hcb, 'Style', 'edit', 'String', num2str(obj.stepTensionL), 'Position', [210, 5, 130, 20],...
                'Enable', 'off','Callback',{@it_callback, 'stepTensionL'});
hitS = uicontrol(hcb, 'Style', 'edit', 'String', num2str(obj.stepTension), 'Position', [210, 25, 130, 20],...
                'Enable', 'off','Callback',{@it_callback, 'stepTension'});
hitEn = uicontrol(hcb, 'Style', 'edit', 'String', num2str(obj.normalFriction), 'Position', [210, 65, 130, 20],...
                'Enable', 'off','Callback',{@it_callback, 'normalFriction'});   
hitEp = uicontrol(hcb, 'Style', 'edit', 'String', num2str(obj.axialFriction), 'Position', [210, 85, 130, 20],...
                'Enable', 'off','Callback',{@it_callback, 'axialFriction'});   
hitEi = uicontrol(hcb, 'Style', 'edit', 'String', num2str(obj.internalFriction), 'Position', [210, 105, 130, 20],...
                'Enable', 'off','Callback',{@it_callback, 'internalFriction'});
hitEz = uicontrol(hcb, 'Style', 'edit', 'String', num2str(obj.zipperFriction), 'Position', [210, 125, 130, 20],...
                'Enable', 'off','Callback',{@it_callback, 'zipperFriction'});    
hitT = uicontrol(hbtn, 'Style', 'edit', 'String', num2str(obj.tmax * obj.dt), 'Position', [100,35,120,25],...
                'Callback', {@it_callback, 'tmax'});        
hitV = uicontrol(hacb,'Style', 'edit', 'String', num2str(obj.moviesetting.vfactor), 'Position', [190,70,120,20],...
                'Enable', 'off', 'Callback', {@it_callback, 'vfactor'});
hitFo = uicontrol(hacb,'Style', 'edit', 'String', num2str(obj.moviesetting.ffactor), 'Position', [190,50,120,20],...
                'Enable', 'off','Callback', {@it_callback, 'ffactor'});
hitTe = uicontrol(hacb,'Style', 'edit', 'String', num2str(obj.moviesetting.tfactor), 'Position', [190,30, 120,20],...
                'Enable', 'off','Callback', {@it_callback, 'tfactor'});      
hitFr = uicontrol(hacb,'Style','edit', 'String', '1', 'Position', [190,10,170,20],...
                'Enable', 'off', 'Visible', 'off');
            
% CHECK BOXES
hcbSL = uicontrol(hcb, 'Style', 'checkbox', 'String', 'Left Tension change (T1)', 'Value', (obj.stepTensionL~=0), ...
                 'Position', [10, 5, 190, 20],'Callback',{@cb_callback,hitSL,'stepTensionL'});
hcbS = uicontrol(hcb, 'Style', 'checkbox', 'String', 'Right Tension change (T2)', 'Value', (obj.stepTension~=0), ...
                 'Position', [10, 25, 190, 20],'Callback',{@cb_callback,hitS,'stepTension'});
hcbEn = uicontrol(hcb, 'Style', 'checkbox', 'String', 'Normal Friction (En)', 'Value', (obj.normalFriction~=0), ...
                 'Position', [10, 65, 190, 20],'Callback',{@cb_callback,hitEn,'normalFriction'});
hcbEp = uicontrol(hcb, 'Style', 'checkbox', 'String', 'Axial Friction (Ep)', 'Value', (obj.axialFriction~=0), ...
                 'Position', [10, 85, 190, 20],'Callback',{@cb_callback,hitEp,'axialFriction'});       
hcbEi = uicontrol(hcb, 'Style', 'checkbox', 'String', 'Internal Friction (Ei)', 'Value', (obj.internalFriction~=0), ...
                 'Position', [10, 105, 190, 20],'Callback',{@cb_callback,hitEi,'internalFriction'});    
hcbEz = uicontrol(hcb, 'Style', 'checkbox', 'String', 'Zipper Friction (Ez)', 'Value', (obj.zipperFriction~=0), ...
                 'Position', [10, 125, 190, 20],'Callback',{@cb_callback,hitEz,'zipperFriction'});                 
hacbT = uicontrol(hacb, 'Style', 'checkbox', 'String', 'tension', 'Value', (obj.moviesetting.tension), ...
                 'Position', [10, 30, 170, 20],'Callback',{@acb_callback,hitTe,'aTension'});
hacbF = uicontrol(hacb, 'Style', 'checkbox', 'String', 'axon shafts', 'Value', (obj.moviesetting.frame), ...
                 'Position', [10, 10, 170, 20],'Callback',{@acb_callback,hitFr,'aFrame'});
hacbFo = uicontrol(hacb, 'Style', 'checkbox', 'String', 'force', 'Value', (obj.moviesetting.force), ...
                 'Position', [10, 50, 170, 20],'Callback',{@acb_callback,hitFo,'aForce'});
hacbV = uicontrol(hacb, 'Style', 'checkbox', 'String', 'velocity', 'Value', (obj.moviesetting.velocity), ...
                 'Position', [10, 70, 170, 20],'Callback',{@acb_callback,hitV,'aVelocity'});             

if hcbS.Value; hitS.Enable = 'on';end;
if hcbSL.Value; hitSL.Enable = 'on';end;
if hcbEn.Value; hitEn.Enable = 'on';end;
if hcbEp.Value; hitEp.Enable = 'on';end;
if hcbEi.Value; hitEi.Enable = 'on';end;
if hcbEz.Value; hitEz.Enable = 'on';end;
if hacbT.Value; hitTe.Enable = 'on';end;
if hacbF.Value; hitFr.Enable = 'on';end;
if hacbFo.Value; hitFo.Enable = 'on';end;
if hacbV.Value; hitV.Enable = 'on';end;
             
% SLIDER BARS
hramptxt = uicontrol(hcb, 'Style', 'text', 'String', strcat('tension onset = ',num2str(obj.ramp)), 'Position', [190, 170, 150, 20]);
hramp = uicontrol(hcb, 'Style', 'slider', 'Max', 1, 'Min', 0.01, 'Value', obj.ramp, ...
               'TooltipString', 'Franction of overall experiment duration to linearly change the tension', ...
               'SliderStep', [0.01, 0.1], 'Position', [10, 170, 170, 20], 'Callback', {@slider_callback,hramptxt,'ramp'});

% DROP DOWN MENU
htext = uicontrol(hgp,'Style', 'text', 'String', 'Quantity to plot', 'Position', [10, 40, 120, 25]);
hdrop = uicontrol(hgp,'Style', 'popupmenu', 'String', {'Velocity', 'Convergence', 'Angle', 'Track', 'Zip Length',...
                  'Left Tension', 'Right Tension', 'Friction Anisotropy', 'Left Length', 'Right Length'},...
                  'Position', [10,10,150,25], 'Callback', {@dropdown_callback});
hcurrentplot = uicontrol(hgp,'Style','text','String','Velocity','Position',[170,10,100,25]);         

% MOVIE CONTROLS
hmovcb = uicontrol(hmp, 'Style', 'checkbox', 'String', 'Generate Movie', 'Value', obj.movieSwitch, 'Position', [10, 60, 200,20],...
                 'Callback', {@movie_callback});             
hmovtxt = uicontrol(hmp, 'Style','text','String','Movie file full name (including path):', 'Position', [10, 35, 300, 20]);
hmovpath = uicontrol(hmp, 'Style', 'edit', 'String', obj.moviepath, 'Position', [10,10,330,20], ...
                     'Callback', {@moviepath_callback}, 'Enable','off');

if hmovcb.Value; hmovpath.Enable = 'on'; end;

% RE-SET V2 POSITION
uicontrol(hV2,'Style', 'text', 'String', 'Vertex coordinates [um]', 'Position', [10, 45, 180, 20]);
uicontrol(hV2,'Style', 'text', 'String', 'Basal tensions (T) [nN]', 'Position', [190, 45, 180, 20]);
uicontrol(hV2,'Style', 'text', 'String', 'Adhesion (S) [nN]', 'Position', [360, 45, 120, 20]); 
hrbS = uicontrol(hV2, 'Style', 'radiobutton', 'String', 'tensions', 'Position', [500, 45, 100, 20]);
hrbC = uicontrol(hV2, 'Style', 'radiobutton', 'String', 'coordinates', 'Position', [500, 25, 100,20]);
hcbSym = uicontrol(hV2, 'Style', 'checkbox', 'String', 'symmetric', 'Position', [500, 5, 100, 20], ...
                    'Value', false, 'Callback', {@cbsym_callback});

hitA = uicontrol(hV2, 'Style', 'edit', 'String', num2str(round(obj.adh,2)),'Enable','on',...
                 'Position', [370, 5, 50, 20], 'Callback', {@it_callback, 'adhesion'});
hitsL = uicontrol(hV2, 'Style', 'edit', 'String', num2str(round(obj.sigmaL,2)),'Enable','on',...
                 'Position', [210, 5, 50, 20], 'Callback', {@it_callback,'sigmaL'});
hitsR = uicontrol(hV2, 'Style', 'edit', 'String', num2str(round(obj.sigmaR,2)),'Enable','on', ...
                 'Position', [265, 5, 50, 20], 'Callback', {@it_callback,'sigmaR'});

hV2x = uicontrol(hV2, 'Style', 'edit', 'String', num2str(round(V2x,2)), 'Position', [10,5,50,20],...
                'Callback',{@it_callback,'V2x'},'Enable','off');
hV2y = uicontrol(hV2, 'Style', 'edit', 'String', num2str(round(V2y,2)), 'Position', [65,5,50,20],...
                'Callback',{@it_callback,'V2y'},'Enable','off');

hgetEnergy = uicontrol(hV2, 'Style', 'pushbutton', 'String', 'E-landscape','Position', [495, 65, 100, 20],...
                'Callback', {@getenergy_callback}, 'Enable', 'on');            
hsetCoor = uicontrol(hV2, 'Style','pushbutton', 'String', 'set', 'Position', [10,25,70,20],...
                'Callback',{@setcoor_callback, hitsL, hitsR },'Enable','off');
hsetTen = uicontrol(hV2, 'Style','pushbutton', 'String', 'set', 'Position', [210,25,70,20],...
                'Callback',{@settension_callback, hV2x, hV2y },'Enable','on');           

% AXES OF the GUI
ha = axes('Units','pixels','Position',[50,100,700,550]);

align([hrun,hplot,htext,hdrop,hcb,hmovcb,hmovtxt,hmovpath,hacb,hsetCoor],'Left','none');

% initialize GUI
set([hf,ha, hrun,hplot,htext,hdrop,hcb,hacb,hitT,hmp,hgp,hbtn,hV2],'Units','normalized');

hf.Name = 'Model GUI';
movegui(hf,'center');

hf.Visible = 'on';

%% Callback functions
% dropdown menu callback
    function dropdown_callback(source, ~)
        % determine the source of call
        str = source.String;
        val = source.Value;
        
        switch str{val}
            case 'Velocity'
                toPlot = 1;
            case 'Convergence'
                toPlot = 2;
            case 'Angle'
                toPlot = 3;
            case 'Track'
                toPlot = 4;
            case 'Zip Length'
                toPlot = 5;
            case 'Left Tension'
                toPlot = 6;
            case 'Right Tension'
                toPlot = 7;
            case 'Friction Anisotropy'
                toPlot = 8;
            case 'Left Length'
                toPlot = 9;
            case 'Right Length'
                toPlot = 10;
        end      
        
    end

% buttons callback
    % run button callback - starts dynamics calculation
    function runbutton_callback(~,~)
        zipfix=false;
        if obj.normalFriction==0 && obj.axialFriction==0 && ...
           obj.zipperFriction==0 && obj.internalFriction==0
            msgbox('At least one type of friction must be non-zero; please check the corresponding box',...
                   'No friction set');
            return;
        end
        if (obj.zipperFriction~=0 && obj.normalFriction==0 && ...
           obj.axialFriction==0 &&  obj.internalFriction==0)
            obj.axialFriction=0.001;
            obj.normalFriction=0.001;
            zipfix = true;
            disp('Substrate friction set to 1 Pas to avoid numerical instability');
        end
        obj.runDynamics();
        hdatatext.String = 'Data generated';
        hplot.Enable = 'on';
        if zipfix;
            obj.axialFriction=0;
            obj.normalFriction=0;
        end;
        assignin('base','velocity',obj.V2vel);
        assignin('base','trajectory',obj.V2trajectory);
    end

    % plotter callback - plots information about quantities during numeric
    % experiment. This means the 'run' must have occured before
    function plotbutton_callback(~,~)
        
        switch toPlot
            case 1
                plot((obj.V2vel(:,1).^2 + obj.V2vel(:,2).^2).^0.5, 'LineWidth',1);
                xlabel('time [s]');
                ylabel('velocity [um/s]');
                str = 'Velocity';
            case 2
                plot(obj.convergence,'-','LineWidth',1);
                xlabel('time [s]');
                ylabel('convergence [um]');
                str = 'Convergence';
                
            case 3
                plot(obj.angle,'-','LineWidth',1);
                xlabel('time [s]');
                ylabel('angle [deg]');
                str = 'Angle';
            case 4
                plot(obj.V2trajectory(:,1),obj.V2trajectory(:,2),'-','LineWidth',1);
                hold on;
                alignstring = 'right';
                plot(obj.V2trajectory(1,1),obj.V2trajectory(1,2),'s','MarkerSize',10);    % mark the initial point
                for ts=round(linspace(1,numel(obj.V2trajectory(:,1)),30))       % mark time stamps along the track
                    if strcmp(alignstring,'left'); alignstring = 'right';    
                    else alignstring = 'left'; end;
                    text(obj.V2trajectory(ts,1), obj.V2trajectory(ts,2), ...
                        strcat('t=', num2str(ts),'s'), 'HorizontalAlignment', alignstring );
                end
                    
                xlabel('x-coordinate [um]');
                ylabel('y-coordinate [um]');
                str = 'Trajectory';
                hold off;
            case 5
                plot(obj.V2length,'-','LineWidth',1);
                xlabel('time [s]');
                ylabel('zip length [um]');
                str = 'Zip length';
            case 6
                plot(obj.A2tension,'-','LineWidth',1);
                xlabel('time [s]');
                ylabel('left axon tension [nN]');
                str = 'Left tension';
            case 7
                plot(obj.A3tension,'-','LineWidth',1);
                xlabel('time [s]');
                ylabel('right axon tension [nN]');
                str = 'Right tension';
            case 8
                plot(obj.Fratio(:,1), obj.Fratio(:,2),'-','LineWidth',1);
                hold on;
                alignstring = 'right';
                plot(obj.Fratio(1,1),obj.Fratio(1,2),'s','MarkerSize',10);    % mark the initial point
                for ts=round(linspace(1,numel(obj.Fratio(:,1)),10))     % mark time stamps along the track
                    if strcmp(alignstring,'left'); alignstring = 'right';    
                    else alignstring = 'left'; end;
                    %text(obj.Fratio(ts,1), obj.Fratio(ts,2), ...
                    %    strcat('t=', num2str(ts),'s'), 'HorizontalAlignment', alignstring );
                    X = [obj.Fratio(ts,1), obj.Fratio(ts,1)];
                    Y = [obj.Fratio(ts,2), obj.Fratio(ts,2)];
                    U = -[obj.FVecs(ts,1), obj.FVecs(ts,3)];
                    V = -[obj.FVecs(ts,2), obj.FVecs(ts,4)];
                    quiver(X,Y,U,V,1/10,'-b','LineWidth',1);
                    daspect([1,1,1]);
                end
                xlabel('x-eigenval');
                ylabel('y-eigenval');
                str = 'Fric. anisotropy';
                hold off;
            case 9
                plot(obj.A2length, '-', 'LineWidth',1);
                xlabel('time[s]');
                ylabel('left axon length [um]');
                str = 'Left length';
            case 10
                plot(obj.A3length, '-', 'LineWidth',1);
                xlabel('time[s]');
                ylabel('right axon length [um]');
                str = 'Right length';
                
        end
        
        hcurrentplot.String = [];
        hcurrentplot.String = str;  % set information text about the current plot
        
    end

% check box callback
    function cb_callback(source, ~, input, varName)
        state = source.Value;
        hplot.Enable = 'off';
        
        if state; 
            input.Enable = 'on';
            switch varName
                case 'stepTension'
                    obj.stepTension = 0.5;
                    input.String = '0.5';
                case 'stepTensionL'
                    obj.stepTensionL = 0.5;
                    input.String = '0.5';
                case 'normalFriction'
                    obj.normalFriction = 0.2;
                    input.String = '200';
                case 'axialFriction'
                    obj.axialFriction = 0.2;
                    input.String = '200';
                case 'internalFriction'
                    obj.internalFriction = 3000;
                    input.String = '3000';
                case 'zipperFriction'
                    obj.zipperFriction = 1.5;
                    input.String = '1.5';
                case 'Asymmetric'
                    obj.V2X(1) = 17.3089;
                    V2x = obj.V2X(1);
                    input.String = round(V2x,2);
            end
        else
            input.Enable = 'off';
            input.String = '0';
            switch varName
                case 'outerForce'
                    obj.outerForce = 0;
                case 'stepTension'
                    obj.stepTension = 0;
                case 'stepTensionL'
                    obj.stepTensionL = 0;
                case 'hookeStiffness'
                    obj.hookeStiffness = 0;
                case 'normalFriction'
                    obj.normalFriction = 0;
                case 'axialFriction'
                    obj.axialFriction = 0;
                case 'internalFriction'
                    obj.internalFriction = 0;
                case 'zipperFriction'
                    obj.zipperFriction = 0;
                case 'Asymmetric'
                    obj.V2X(1) = 17.3089;
                    V2x = obj.V2X(1);
                    input.String = round(V2x,2);
            end
        end;
        
    end

% sets coordinate of V2 vertex and calculates appropriate tensions
    function setcoor_callback(~, ~, hsL, hsR)
        
       newcoor = [V2x, V2y];
       obj.positionStaticVertex(newcoor);
       
       hsL.String = num2str(round(obj.sigmaL,2));
       sL = obj.sigmaL;
       hsR.String = num2str(round(obj.sigmaR,2));
       sR = obj.sigmaR;
       
    end

% sets tensions in vertex V2 and calculates appropriate location
    function settension_callback(~, ~, hcX, hcY)
        
        newtens = [sL, sR];
        obj.adjustStaticVertex(newtens);
        
        hcX.String = num2str(round(obj.V2X(1),2));
        V2x = obj.V2X(1);
        hcY.String = num2str(round(obj.V2X(2),2));
        V2y = obj.V2X(2);
        
    end

% generates and displays energy landscape for the chosen 
    function getenergy_callback(~,~)
        
        span = 1;
        
        obj.energyscapeAnimation = struct('cdata',[],'colormap',[]);
        obj.energyscapeAnimation(span) = struct('cdata',[],'colormap',[]);
        obj.energyscapeSlowTrail = [];
        obj.energyscapeSlowTrail(1,:) = obj.V2X;
        
        disp('Starting calculation');
        right = sR + obj.stepTension;
        left  = sL + obj.stepTensionL;
        obj.energyscapeAnimation = obj.getEnergyLandscape(right, left,'single',true,'scale',4);
        disp('Finished');
        
    end

% animation check box callback
    function acb_callback(source, ~, input, variable)
        
       state = source.Value;
       
       if state; input.Enable = 'on';
       else input.Enable = 'off';
            input.String = '1';  
       end;
       
       switch variable
           case 'aTension'
               obj.moviesetting.tension = state;
           case 'aFrame'
               obj.moviesetting.frame = state;
           case 'aForce'
               obj.moviesetting.force = state;
           case 'aVelocity'
               obj.moviesetting.velocity = state;
           case 'aFriction'
               obj.moviesetting.ftensor = state;
           case 'aMotion'
               obj.moviesetting.mtensor = state;
       end
      
       hitFr.Enable = 'off';
        
    end

% symmetric checkbox callback
    function cbsym_callback(source,~)
        
        if source.Value;
            hV2x.Enable = 'off';
            hitsR.Enable = 'off';
            sR = sL;
            V2x = round(17.3089,2);
            hV2x.String = num2str(round(V2x,2));
            hitsR.String = num2str(round(sR,2));
        else
            hV2x.Enable = 'on';
            hV2y.Enable = 'on';
            hsetCoor.Enable = 'on';
            hrbS.Value = false;
            hrbC.Value = true;
            hitsL.Enable = 'off';
            hitsR.Enable = 'off';
            hsetTen.Enable = 'off';            
        end
    end


% radiobutton selection callback

    function vertexselection(~, callbackdata)
       
        if strcmp(callbackdata.NewValue.String, 'tensions')
            hV2x.Enable = 'off';
            hV2y.Enable = 'off';
            hsetCoor.Enable = 'off';
            hitsL.Enable = 'on';
            if (~hcbSym.Value); hitsR.Enable = 'on'; end;
            hsetTen.Enable = 'on';
        elseif strcmp(callbackdata.NewValue.String, 'coordinates')
            hitsL.Enable = 'off';
            hitsR.Enable = 'off';
            hsetTen.Enable = 'off';
            if (~hcbSym.Value); hV2x.Enable = 'on'; end;
            hV2y.Enable = 'on';
            hsetCoor.Enable = 'on';
        end
            
    end
    
% input text callback
    function it_callback(source, ~, varName)
        
        val = str2double(source.String);
        hplot.Enable = 'off';
        
        switch varName
            case 'outerForce'
                obj.outerForce = val;
            case 'stepTension'
                obj.stepTension = val;
            case 'stepTensionL'
                obj.stepTensionL = val;
            case 'hookeStiffness'
                obj.hookeStiffness = val;
            case 'tmax'
                obj.setTime(round(val / obj.dt));
            case 'normalFriction'
                obj.normalFriction = val/1000;  % not in Pas
            case 'axialFriction'
                obj.axialFriction = val/1000;   % not in Pas
            case 'internalFriction'
                obj.internalFriction = val;
            case 'zipperFriction'
                obj.zipperFriction = val;
            case 'tfactor'
                obj.moviesetting.tfactor = val;
            case 'ffactor'
                obj.moviesetting.ffactor = val;
            case 'vfactor'
                obj.moviesetting.vfactor = val;
            case 'V2x'
                V2x = val;
            case 'V2y'
                V2y = val;
            case 'sigmaL'
                sL = val;
                if (hcbSym.Value); 
                    sR = sL; 
                    hitsR.String = num2str(round(sR,2));
                end;
            case 'sigmaR'
                sR = val;
            case 'adhesion'
                obj.adh = val;
        end
        
    end
   
% slider callback
    function slider_callback(source, ~, label, variable )
        
        hplot.Enable = 'off';
       
        switch variable
            case 'ramp'
                obj.ramp = round(source.Value,2);
                label.String = strjoin({'tension onset =', num2str(obj.ramp)});
            case 'toff'
                obj.setToff(round(source.Value,2));
                label.String = strjoin({'toff =', num2str(obj.tfrac)});
        end
        
    end

% film callback
    function movie_callback(source,~)
       
        val = source.Value;
        
        if val; 
            hmovpath.Enable = 'on';
            obj.movieSwitch = true;
        else
            hmovpath.Enable = 'off';
            obj.movieSwitch = false;
        end
    end

% filmpath callback
    function moviepath_callback(source,~)
        
        obj.moviepath = source.String;
        
    end


end

