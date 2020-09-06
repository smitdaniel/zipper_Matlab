%% This function constructs an animation of the vertex movement
% V1 = coordinate of vertex 1 (left vertex)
% V2 = coordinate of vertex 2 (right vertex)
% X3 = bottom fixation point coordinate for vertex V2
% OX = right fixation point coordinate for vertex V2
% O = mobile point on the right arm of the axon (Fo acts here)
% Fo = outer force on point O
% Fox = outer force on point OX
% V2v = velocity of the vertex V2
% Tx = tension in axon x
% V2f = force on vertex 2
% ftensor = friction tensor for V2
% eVec, eVals = eigen vectors and eigenvalues of V2 friction matrix
% settings = booleans as what to plot
%
% frame = generated video frame passed back
% ==================================================

function [frame] = vertexAnimator(V1, V2, X3, OX, O, Fo, Fox, V2v, T2, T3, V2f, ftensor, eVals, setting)

    margin = 5;
    
    fh = figure('Visible','off','Position',[0,0,900,1200]);
    
    X = [ V1(1),V2(1),X3(1),OX(1),O(1)];
    Y = [ V1(2),V2(2),X3(2),OX(2),O(2)];
    xrange = [min(X)-margin, max(X)+margin];
    yrange = [min(Y)-margin, max(Y)+margin];
    
    ax = gca;   % get axis handle
    cla(ax);    % clear axis
    
    ax.XLim = xrange;   % set x axis limits
    ax.YLim = yrange;   % set y axis limits
    
    hold on;
    
    if setting.frame;
        tinyPlot(V1,V2,'line');
        tinyPlot(V2,X3,'line');
        tinyPlot(V2,O,'line');
        tinyPlot(O,OX,'line');
    end;
    
    if setting.tension;
        tinyPlot(V2, T2*setting.tfactor, 'vector','g');
        tinyPlot(V2, T3*setting.tfactor, 'vector','g');
    end;
    
    if setting.force
        tinyPlot(V2, V2f*setting.ffactor, 'vector','b');
        tinyPlot(O, Fo*setting.ffactor, 'vector','b');
        tinyPlot(OX, Fox*setting.ffactor, 'vector','b');
    end;
    
    if setting.velocity
        tinyPlot(V2, V2v*setting.vfactor, 'vector','k');
    end;
    
    if setting.ftensor
        tensorPlot( V2, ftensor, eVals, setting.frfactor );
    end
    
    if setting.mtensor
        tensorPlot( V2, -inv(ftensor), -[1/eVals(1),1/eVals(2)], setting.mfactor );
    end
    
    daspect([1,1,1]);
    frame = getframe(gcf);
    close(fh);
    
end

% call to plot tensor representation
function tensorPlot( V2, tensor, eVals, factor )

        maxEval = max(abs(eVals));
        sampling = 60;
        time = linspace(0,2*pi,sampling);
        coor = zeros(sampling,2);
        i=1;
        for t=time;
            COOR = [cos(t); sin(t)];
            coor(i,:) = V2 - factor/maxEval*transpose(tensor * COOR);
            i=i+1;
        end;
        for q=[1,16]
            plot([V2(1), coor(q,1)], [V2(2), coor(q,2)],'m','LineWidth',1.5);
        end;
        plot(coor(:,1),coor(:,2),'m','LineWidth',1.5);
end

% call to plot line/vector representation
function tinyPlot( A, B, type, color )

    switch type
        case 'line'
            plot([A(1),B(1)], [A(2),B(2)], '-r', 'LineWidth',2.5);
        case 'vector'
            quiver(A(1),A(2),B(1),B(2),color,'LineWidth',1.5,'AutoScale','off');
        otherwise

    end;
    
end
