% MATFLOWSIM -- simulates unidrectional cargo (particle) transport in neuron dendrites
% 
% matflowsim (path, varargin)
% -----------------------------------
% 
% Input
% -----
%
% - Mandatory Parameters
%
%     - path::    string: full path to neuron morphology file, including extension 
%                 (see function 'load_tree' for accepted file types)
%         
% Optional Parameters
%
%     - rule::    string: rule determining probabilities for cargo choosing 
%                 children neurites at branch points (must be a TREES toolbox 
%                 function that returns a column vector) {DEFAULT: 'cvol_tree'}
%
%     - rtime::   scalar: total run time of simulation {DEFAULT: 60}
% 
%     - freq::    scalar: frequency of new cargo generation (sec) {DEFAULT: 5}
% 
%     - amp::     scalar: number of cargo generated each wave {DEFAULT: 3}
% 
%     - mvel::    scalar: mean velocity of cargo {DEFAULT: 1}
%
%     - sprob::   scalar: probability for cargo's direction to change {DEFAULT: 0.005}
%
%     - xspeed::  scalar: simulation speed-up factor {DEFAULT: 1}
%
%     - options:: string: using special simulation features {DEFAULT: '-m'}
%                 '-m': make .mp4 movie out of frame data
%                 '-s': simplify neuron tree morphology for fast simulation
%
% Example
% -------
% matflowsim ('C:\\neurons\\cell-124-trace.CNG.swc')
% 
% Written by Douglas Qian

function matflowsim(path, varargin)
    %% PARSING INPUT ARGUMENTS %%
    p = inputParser;
    % Setting validation functions
    valfcn1 = @(x) validateattributes(x, {'char'},   {'nonempty'});
    valfcn2 = @(x) validateattributes(x, {'double'}, {'scalar', 'positive', 'nonempty'});
    % Parameter names, default values, and corresponding validation functions
    params = {'path'  'rule'       'rtime'  'freq'  'amp'  'mvel'  'sprob'  'xspeed'  'options'};
    dfvals = {NaN     'cvol_tree'   60       5       3      1       0.005    1        '-m'     };
    valfcn = [1        1            2        2       2      2       2        2         1       ];
    % Adding arguments to input parser scheme
    for i = 1: length(params)
        if valfcn(i) == 1
            if i == 1
                addRequired(p, 'path', valfcn1);
            else
                addParameter(p, params{i}, dfvals{i}, valfcn1);
            end
        else
            addParameter(p, params{i}, dfvals{i}, valfcn2);
        end
    end
    % Parsing function inputs
    parse(p, path, varargin{:});
    % Extracting variable values
    path    = p.Results.path;
    rule    = str2func(p.Results.rule);
    rtime   = p.Results.rtime;
    freq    = p.Results.freq;
    amp     = p.Results.amp;
    mvel    = p.Results.mvel;
    sprob   = p.Results.sprob;
    xspeed  = p.Results.xspeed;
    options = p.Results.options;
    % fprintf('path: %s\nrule: %s\nrtime: %d\nfreq: %d\namp: %d\nmvel: %f\nsprob: %f\noptions: %s\n', path, rule, rtime, freq, amp, mvel, sprob, options)
    
    %% SIMULATION SETUP %%
    sneuron = load_tree(path);
    rule_data = asym_tree(sneuron, rule(sneuron));
    if ~isempty(options) && ~isempty(strfind(options, '-s'))
        sneuron = simplify_tree(sneuron);
        rule_data = asym_tree(sneuron, rule(sneuron));
    end
    % 2D plotting of neuron morphology
    f = figure;
    set(f, 'Position', [0 0 1000 1000], 'Color', 'white')
    plot_tree(sneuron, [0 0 0], [0 0 0], [], 8, '-2q');
    axis off, hold on,
    % Cargo (ptle) parameters
    radius = 4;
    step_size = shortest_segment(sneuron);
    ptles = struct([]);
    % mean velocity (~1 micron/s) standardization
    framerate = mvel/step_size*xspeed;
    
    typetree = typeN_tree(sneuron);
    
    %% MOVIE MAKER %%
    if ~isempty(options) && ~isempty(strfind(options, '-m'))
        if ~isempty(options) && ~isempty(strfind(options, '-s'))
            vname = sprintf('%s simplified (%s %0.2fmps)', sneuron.name(1:end-4), func2str(rule), mvel*xspeed);
        else
            vname = sprintf('%s (%s %0.2fmps)', sneuron.name(1:end-4), func2str(rule), mvel*xspeed);
        end
        v = VideoWriter(vname, 'MPEG-4');
        v.FrameRate = framerate;
        open(v);
    end

    %% ZOOM HANDLING %%
    orig_ax = axis;
    h = zoom;
    set(h,'ActionPostCallback',@zoomCallBack);
    function zoomCallBack(~, evd)      
        new_ax = axis(evd.Axes);    
        % Zoom based on changes in y-axis lengths
        orig_len = orig_ax(4)-orig_ax(3);
        new_len = new_ax(4)-new_ax(3);
        for pt = 1: length(ptles)
            set(ptles(pt).graphics, 'MarkerSize', radius*(orig_len/new_len))
        end
    end

    %% SIMULATION LOOP %%
    iterations = ceil(rtime * framerate);
    
    for i = 0: iterations
       fprintf('%d\n', i)
       % Generating new ptles
       if mod(i, round(freq*framerate)) == 0
           for j = 1: amp
               ptles(end+1).graphics = plot(sneuron.X(1), sneuron.Y(1), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', radius);
               ptles(end).target = [NaN 1];
               ptles(end).dir = 1;
               ptles(end).speed = normrnd(step_size, step_size*0.25, 1);
           end
       end
       % Move all ptles
       for pt = 1: length(ptles)
           x1 = ptles(pt).graphics.XData;
           y1 = ptles(pt).graphics.YData;
           dir = ptles(pt).dir;
           target = ptles(pt).target(dir+1);
           x2 = sneuron.X(target);
           y2 = sneuron.Y(target);
%            fprintf('-----------\n')
%            fprintf('X: %f Y: %f\n', x1, y1)
%            fprintf('DIR: %d target: %d\n',dir, target)
               
           if pdist([x1 y1; x2 y2], 'euclidean') < step_size
               % Jump to point first
               ptles(pt).graphics.XData = x2;
               ptles(pt).graphics.YData = y2;
               % TYPEN_TREE VALUES
               % 0: terminating point
               % 1: continuation point
               % 2: branch point
               tartype = typetree(target);
               % FORWARD MOVING
               if dir == 1
                   if tartype == 2
                       result = rand;
                       cumprob = rule_data(target).vals(:,1);
                       indices = rule_data(target).vals(:,2);
                       % Choosing a child
                       for c = 1: length(indices)
                           if result <= cumprob(c)
                               ptles(pt).target = [target indices(c)];
                               break
                           end
                       end
                   elseif tartype == 1
                       ptles(pt).target = [target find(sneuron.dA(:, target))];
                   else
                       delete(ptles(pt).graphics);
                       ptles(pt) = [];
                       break
                   end
               % BACKWARD MOVING
               else
                   if target == 1
                       delete(ptles(pt).graphics);
                       ptles(pt) = [];
                       break
                   end
                   ptles(pt).target = [find(sneuron.dA(target, :)) target];
               end
           end    
           
           % Take a step in the direction towards the target
           x1 = ptles(pt).graphics.XData;
           y1 = ptles(pt).graphics.YData;
           target = ptles(pt).target(dir+1);
           x2 = sneuron.X(target);
           y2 = sneuron.Y(target);
           
           ydiff = y2-y1;
           xdiff = x2-x1;
           
           theta = rad2deg(atan2(ydiff, xdiff)) - 90;
           if ~mod(theta, 180) == 0
              theta = -theta;
           end
           
           ptles(pt).graphics.XData = x1 + ptles(pt).speed * sin(deg2rad(theta));
           ptles(pt).graphics.YData = y1 + ptles(pt).speed * cos(deg2rad(theta));
           
           % Slight probability of changing directions
           change = rand;
           if change < sprob
               ptles(pt).dir = ~ptles(pt).dir;
           end
           
       end
              
       drawnow
       
       if ~isempty(options) && ~isempty(strfind('-m', options))
           writeVideo(v, getframe(gcf));
       end
       
    end
    
    if ~isempty(options) && ~isempty(strfind('-m', options))
        close(v);
    end
end