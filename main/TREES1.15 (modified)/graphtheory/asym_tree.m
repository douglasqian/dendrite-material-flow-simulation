% ASYM_TREE   Branch point asymmetry.
% (trees package)
% 
% [asym_data] = asym_tree (intree, v, options)
% -------------------------------------
%
% calculates for each branching point the ratio of the sums of the two
% daughter branches. The summed values are given by v which attributes a
% value to each node. Typically this can be the count of terminals
% (default) or the tree length etc... For v1 is smaller summed value of
% sub-trees and v2 the other one: v1/(v1 + v2). Reports NaN where there is
% no branch point. Tree must be BCT (at least trifurcations are forbidden
% of course), use "repair_tree" if necessary.
%
% Input
% -----
% - intree::integer:index of tree in trees structure or structured tree
% - v::vertical vector: {DEFAULT: count child terminals (== "T_tree")}
% - options::string: {DEFAULT: ''}
%     '-s' : show
%     '-m' : explanatory movie 
% 
% Output
% ------
% - asym_data::structure: contains ratios for all branching points (index:
% B_tree) and indices of corresponding children
%           fields: {vals (array)}
%
% Example
% -------
% asym_tree (sample_tree, [], '-m -s')
%
% See also child_tree sub_tree
% Uses ipar_tree B_tree ver_tree dA
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function [asym_data] = asym_tree (intree, v, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length(trees); % {DEFAULT tree: : last tree in trees cell array} 
end;

ver_tree(intree); % verify that input is a tree structure

% use only directed adjacency for this function
if ~isstruct(intree),
    dA = trees{intree}.dA;
else
    dA = intree.dA;
end

if (nargin <2)||isempty(v),
    v = T_tree(intree); % {DEFAULT vector: count termination points}
end

if (nargin < 3)||isempty(options),
    options = ''; % {DEFAULT: no option}
end

iB = find(B_tree(intree));  % index of branching points
ipar = ipar_tree(intree);   % parent index paths (see "ipar_tree")
asym_data(length(dA)).vals = []; % structure array containing asymmetry values for each BP

% modified to account for branching points with more than 1 children

for ward = 1: length(dA),
   if ~any(ward==iB),
       asym_data(ward).vals = [NaN NaN];
   end
end

for ward = 1:length(iB),
    BB = find(dA(:,iB(ward))); % index of children of branch points
    nchild = length(logical(BB)); % number of children
    sums = zeros(nchild,1);

    for i = 1: nchild
        [sub i2] = ind2sub(size(ipar), find(ipar==BB(i)));
        sums(i) = sum(v(sub)); % summed values for all sub-trees
    end

    total = sum(sums);
    cumratio = 0;
    % asymmetry cumulative ratios
    for i = 1: nchild
       ratio = sums(i)/total + cumratio;
       asym_data(iB(ward)).vals = [asym_data(iB(ward)).vals; ratio, BB(i)];
       cumratio = ratio;
    end
    
    if strfind(options,'-m'), % movie option
        clf; hold on; shine;
        HP = plot_tree(intree); set(HP, 'facealpha', 0.2);
        plot_tree(intree,[1 0 0],[],sub1);
        plot_tree(intree,[0 1 0],[],sub2);
        HT = text(0,0,num2str(asym(ward)));
        set(HT,'fontsize',12,'color',[1 0 0]);
        title ('asymmetry at branch points');
        xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
        view(2); grid on; axis image;
        pause(.4);
    end
end


if strfind(options,'-s'), % show option
    clf; hold on; shine;
    HP = plot_tree(intree,[],[],find(~B_tree(intree))); set(HP, 'facealpha', 0.2);
    iB = find(B_tree(intree)); plot_tree(intree,asym_vals(iB),[],iB);
    title (['asymmetry at branch points, mean: ' num2str(nanmean(asym_vals))]);
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view(2); grid on; axis image;
    set(gca,'clim',[0 0.5]); colorbar;
end
