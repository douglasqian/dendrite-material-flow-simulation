% path example = 'C:\\Users\\douglas.qian\\Documents\\neuron_morphologies\\test\\ganglion\\retina\\chalupa\\cell-124-trace.CNG.swc'

function matflowsim(path)
    % Loading a neuron tree from an .swc file
    sneuron = load_tree(path);
    plot_tree(sneuron, [0,0,0], [0,0,0], [], 8, '-2q');
    hold on
    
    % All ptle are stored in struct array (ptles)
    nptles = 2;
    radius = 4;
    step_size = shortest_segment(sneuron);
    ptles(nptles).target = NaN;
    
    for pt = 1: nptles
        ptles(pt).graphics = plot(sneuron.X(1), sneuron.Y(1), 'ro', 'MarkerSize', radius);
        ptles(pt).target = NaN;
    end
  
    id = transpose(1:length(sneuron.dA));
    endpoints = [sneuron.X(T_tree(sneuron)), sneuron.Y(T_tree(sneuron)), id(T_tree(sneuron))];
    bifpoints = [sneuron.X(B_tree(sneuron)), sneuron.Y(B_tree(sneuron)), id(B_tree(sneuron))];
    contpoints = [sneuron.X(C_tree(sneuron)), sneuron.Y(C_tree(sneuron)), id(C_tree(sneuron))];
    
    % SPLITTING RULE = BRANCH ASYMMETRY (VOLUME)
    asym_data = asym_tree(sneuron, cvol_tree(sneuron));
    
    %% ZOOM HANDLING %%
    orig_ax = axis;
    h = zoom; % handle to zoom utility
    set(h,'ActionPostCallback',@zoomCallBack);
    function zoomCallBack(~, evd)      
        new_ax = axis(evd.Axes);    
        % Change in y-axis lengths
        orig_len = orig_ax(4)-orig_ax(3);
        new_len = new_ax(4)-new_ax(3);
        for i = 1: length(ptles)
            set(ptles(i).graphics, 'MarkerSize', 8*(orig_len/new_len))
        end
    end

    iterator = 0;
    %% SIMULATION LOOP %% (continues until allreachedend is true)
    allreachedend = 0;
    while ~allreachedend
       % move ptles
       for pt = 1: length(ptles)
           [close2B, bid] = is_close_member([ptles(pt).graphics.XData, ptles(pt).graphics.YData], bifpoints, step_size, ptles(pt).target);
           [close2C, cid] = is_close_member([ptles(pt).graphics.XData, ptles(pt).graphics.YData], contpoints, step_size, ptles(pt).target);
           if close2B % reached bifurcation point
               % jump to bif point
               ptles(pt).graphics.XData = sneuron.X(bid);
               ptles(pt).graphics.YData = sneuron.Y(bid);
               % choose branch according to rules               
               nchild = size(asym_data(bid).vals, 1);
               cumprob = asym_data(bid).vals(:,1);
               indices = asym_data(bid).vals(:,2);
               
               result = rand;
               % using asymmetry ratios
               for c = 1: nchild
                   if result <= cumprob(c)
                       ptles(pt).target = indices(c);
                       break
                   end
               end
                              
           elseif close2C % on continuation point
                              
               % jump to cont point
               ptles(pt).graphics.XData = sneuron.X(cid);
               ptles(pt).graphics.YData = sneuron.Y(cid);
               % continue to id's only child (find(sneuron.dA(:,iB(id)))
               ptles(pt).target = find(sneuron.dA(:, cid));
                                         
           end
                      
           % take a step towards the target
           x1 = ptles(pt).graphics.XData;
           y1 = ptles(pt).graphics.YData;
           x2 = sneuron.X(ptles(pt).target);
           y2 = sneuron.Y(ptles(pt).target);
           
           ydiff = y2-y1;
           xdiff = x2-x1;
           
           theta = rad2deg(atan2(ydiff, xdiff)) - 90;
           if ~mod(theta, 180) == 0
              theta = -theta;
           end
           
           ptles(pt).graphics.XData = ptles(pt).graphics.XData + step_size * sin(deg2rad(theta));
           ptles(pt).graphics.YData = ptles(pt).graphics.YData + step_size * cos(deg2rad(theta));
           
       end
       
       drawnow % animation (updating screen)
              
       allreachedend = 1;
       for pt = 1: length(ptles)
           [close2T, tid] = is_close_member([ptles(pt).graphics.XData, ptles(pt).graphics.YData], endpoints, step_size, ptles(pt).target);
           if ~close2T
               allreachedend = 0;
               break
           end
       end
       
       iterator = iterator + 1;
       
    end
end