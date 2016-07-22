% MATFLOWSIM (v2.3) - dendritic bidrectional cargo transport simulator
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
%     - rule::    string: name of branching rule {DEFAULT: 'cvol_tree'}
%
%     - rtime::   scalar: total run time of simulation movie [s] {DEFAULT: 30}
% 
%     - freq::    scalar: frequency of new cargo generation [s] {DEFAULT: 2}
% 
%     - amp::     scalar: number of cargo generated in each wave {DEFAULT: 5}
% 
%     - mvel::    scalar: mean velocity of cargo [micron/s] {DEFAULT: 1}
%
%     - radius::  scalar: radius of cargo {DEFAULT: 4}
%
%     - sprob::   scalar: probability for cargo's direction to change {DEFAULT: 0.01}
%
%     - xspeed::  scalar: simulation speed-up factor {DEFAULT: 8}
%
%     - options:: string: using special simulation features {DEFAULT: '-s'}
%                 '-s': resample neuron tree for faster simulation
%                 '-m': make .mp4 movie out of frame data
%                 '-h': ptle density heat map
%                 '-k': kymograph
%                 '-d': display statistics
%
% Example
% -------
% matflowsim ('C:\\neurons\\cell-124-trace.CNG.swc')
% 
% Written by Douglas Qian

function matflowsim(path, varargin)
    %% PARSING INPUT ARGUMENTS %%
    
    p = inputParser;
    
    % Input validation functions
    valfcn1 = @(x) validateattributes(x, {'char'},   {'nonempty'});
    valfcn2 = @(x) validateattributes(x, {'double'}, {'scalar', 'positive', 'nonempty'});
    
    % Parameter names, default values, and corresponding validation functions
    params = {'path'  'rule'       'rtime'  'freq'  'amp'  'mvel'  'radius'  'sprob'  'xspeed'  'options'};
    dfvals = {NaN     'cvol_tree'   30       2       5      1       4         0.00     8        '-s'     };
    valfcn = [1        1            2        2       2      2       2         2        2         1       ];
    
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
    
    % Parsing inputs
    parse(p, path, varargin{:});
    
    % Extracting variable values
    path    = p.Results.path;
    rule    = str2func(p.Results.rule);
    rtime   = p.Results.rtime;
    freq    = p.Results.freq;
    amp     = p.Results.amp;
    mvel    = p.Results.mvel;
    radius  = p.Results.radius;
    sprob   = p.Results.sprob;
    xspeed  = p.Results.xspeed;
    options = p.Results.options;
    
    %% SIMULATION SETUP %%
    
    % Extracting neuron morphology from file
    sneuron = load_tree(path);
    
    % Probabilities generated based off of rule (Nx1 vector)
    rule_data = asym_tree(sneuron, rule(sneuron));
    
    % Simplifying (using TREEs resample_tree) the morphology
    if ~isempty(strfind(options, '-s'))
        sneuron = resample_tree(sneuron, 10);
        rule_data = asym_tree(sneuron, rule(sneuron));
    end
    
    % Neuron morphology plotting
    f1 = figure;
    f1.Name = 'Simulation'; f1.Position = [0 0 1000 1000]; f1.Color = 'white';
    plot_tree(sneuron);
    axis off, hold on,
    
    % Calculating step size 
    step_size = shortest_segment(sneuron);
    
    % Setting up kymograph parameters
    kymo = 0;
    if ~isempty(strfind(options, '-k'))
        % Kymograph region of interest is longest segment
        [~,kroi] = longest_segment(sneuron);
        % Kymograph only samples middle-third of simulation
        kdur = [round(0.33*rtime) round(0.67*rtime)];
        kymo = 1;
    end
    
    % Options switches
    dstats = 0;
    if ~isempty(strfind(options, '-d'))
        dstats = 1;
    end
    
    hmap = 0;
    if ~isempty(strfind(options, '-h'))
        hmap = 1;
    end
    
    % Unique cargo ID's
    ptleidct = 0;
    ptids    = transpose(1:length(sneuron.dA)); % Nx1 vector storing id's of points
    
    % Pre-existing cargo
    branchpts = [sneuron.X(B_tree(sneuron)) sneuron.Y(B_tree(sneuron)) ptids(B_tree(sneuron))];
    branchpts(1,:) = [];
    ptles(length(branchpts)).target = NaN; % struct array storing information about all cargo particles
    for pt = 1: length(branchpts)
       if rand >= 0.5
           ptles(pt).graphics = plot(branchpts(pt,1), branchpts(pt,2), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', radius);
           ptles(pt).dir = 1;
       else
           ptles(pt).graphics = plot(branchpts(pt,1), branchpts(pt,2), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', radius);
           ptles(pt).dir = 0;
       end
       ptles(pt).target = [find(sneuron.dA(branchpts(pt,3),:)) branchpts(pt,3)];
       ptles(pt).speed  = normrnd(step_size, step_size*0.25, 1);
       ptles(pt).id     = ptleidct;
       ptleidct         = ptleidct + 1;
    end
    
    % Other data structures
    if hmap
        ptcoords  = []; % Zx2 vector of ptle coordinates
    end
    if dstats
        density   = zeros(length(sneuron.dA), 1); % Nx1 vector of ptle counts
    end
    typetree  = typeN_tree(sneuron); % Nx1 vector of point types
    endpts    = [sneuron.X(T_tree(sneuron)) sneuron.Y(T_tree(sneuron)) ptids(T_tree(sneuron))]; % Mx2 vector of endpoint coordinates
    
    %% MOVIE MAKER %%
    
    % Simulation speed standardization
    framerate = round(mvel/step_size*xspeed);
    
    % Creating movie
    if ~isempty(strfind(options, '-m'))
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
    z = zoom;
    set(z,'ActionPreCallBack', @zoomPreCallBack, 'ActionPostCallback', @zoomPostCallBack);
    last_ax = axis;
    function zoomPreCallBack(~, evd)
       last_ax = axis(evd.Axes); 
    end
    function zoomPostCallBack(~, evd)      
        new_ax = axis(evd.Axes);    
        % Zoom based on changes in y-axis lengths
        last_len = last_ax(4)-last_ax(3);
        new_len = new_ax(4)-new_ax(3);
        radius = radius * (last_len/new_len);
        for pt = 1: length(ptles)
            set(ptles(pt).graphics, 'MarkerSize', radius)
        end
    end

    %% SIMULATION LOOP %%
    
    % NOTE: 1 iteration = 1 frame
    iterations = rtime * framerate;

    for i = 0: iterations
       fprintf('%d\n', i)
       if mod(i, freq*framerate) == 0
           % Generate cargo coming from soma (root) and endpts
           tpoints = randperm(size(endpts, 1), amp);
           for j = 1: 2*amp
               if j <= amp
                   % Coming from endpts (endocytosis) --> red, dir = 0
                   tid = endpts(tpoints(j), 3);
                   ptles(end+1).graphics = plot(endpts(tpoints(j), 1), endpts(tpoints(j), 2), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', radius);
                   ptles(end).target = [find(sneuron.dA(tid,:)) tid];
                   ptles(end).dir = 0;
               else
                   % Going towards endpts (exocytosis) --> green, dir = 1
                   ptles(end+1).graphics = plot(sneuron.X(1), sneuron.Y(1), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', radius);
                   ptles(end).target = [NaN 1];
                   ptles(end).dir = 1;
               end
               ptles(end).speed = normrnd(step_size, step_size*0.25, 1);
               ptles(end).id = ptleidct;
               ptleidct = ptleidct + 1;
           end
       end
       
       % Start logging data for kymograph
       if kymo && i == kdur(1) * framerate
           kids = []; % tracking ptle id's
           kdata = zeros(diff(kdur)*framerate, 1); % distance of ptles from start of kroi
           % NOTE: length(kids) should equal size(kdata, 3) 
           
       % Stop logging data for kymograph    
       elseif kymo && i == kdur(2)* framerate
           kymo = 0;
       end
       
       % LOOPING THROUGH ALL PTLES %
       for pt = 1: length(ptles)
           % Accounts for ptles deleted within one iteration
           if pt > length(ptles)
               break
           end
           
           % Cargo position data
           x1     = ptles(pt).graphics.XData;
           y1     = ptles(pt).graphics.YData;
           dir    = ptles(pt).dir;
           target = ptles(pt).target(dir+1);
           x2     = sneuron.X(target);
           y2     = sneuron.Y(target);
           
           % Logging kymograph data
           if kymo && i >= kdur(1) * framerate && ptles(pt).target(2) == kroi
               num = find(kids == ptles(pt).id);
               ptarget = ptles(pt).target(1); % or find(sneuron.dA(kroi,:))
               px = sneuron.X(ptarget); py = sneuron.Y(ptarget);
               dist = pdist([x1 y1; px py], 'euclidean');
               if ~isempty(num)
                   kdata(i-kdur(1)*framerate+1, 1, num) = dist;
               else
                   kids(end+1) = ptles(pt).id;
                   nnum = find(kids == ptles(pt).id);
                   kdata(i-kdur(1)*framerate+1, 1, nnum) = dist;
               end
           end
           
           % Logging heat map and statistics data
           if dstats
               density(ptles(pt).target(2)) = density(ptles(pt).target(2)) + 1;
           end
           if hmap
               ptcoords = [ptcoords; x1 y1];
           end
               
           % Checking if cargo is close enough to target
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
                       continue
                   end
               % BACKWARD MOVING
               else
                   if target == 1
                       delete(ptles(pt).graphics);
                       ptles(pt) = [];
                       continue
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
       
       % Write frame into movie object
       if ~isempty(strfind(options, '-m'))
           writeVideo(v, getframe(gcf));
       end
       
    end
    
    %% HEAT MAP %%
    if hmap
        f2 = figure;
        f2.Name = 'Heat Map';
        [values, centers] = hist3(ptcoords, [200,200]);
        imagesc(centers{:}, values);
        colormap jet,
        colorbar,
    end
    
    %% KYMOGRAPH %%
    if ~isempty(strfind(options, '-k'))
        f3 = figure;
        f3.Name = 'Kymograph';
        timescale = kdur(1): 1/framerate: kdur(2)-1/framerate;
        kdata(kdata == 0) = NaN;
        
        for i = 1: size(kdata,3)
           plot(timescale, kdata(:,:,i));
           hold on;
        end
        title('Kymograph'), ylabel('Distance to soma [\mum]'), xlabel('Time [s]'),
    end    
    
    %% DISPLAY STATISTICS %%
    if dstats
        density = density/iterations;
        f4 = figure;
        f4.Name = 'Statistics';
        a = subplot(2,1,1);
        scatter(a, eucl_tree(sneuron), density, 'r.')
        title(a, 'Sholl Analysis'), ylabel(a, 'Number of particles'), xlabel(a, 'Distance from cell body [\mum]'),
        b = subplot(2,1,2);
        scatter(b, cvol_tree(sneuron), density, 'g.')
        title(b, 'Volume'), ylabel(b, 'Number of particles'), xlabel(b, 'Segment volume [1/\mum]'),
    end

    if ~isempty(strfind(options, '-m'))
        close(v);
    end
    
end