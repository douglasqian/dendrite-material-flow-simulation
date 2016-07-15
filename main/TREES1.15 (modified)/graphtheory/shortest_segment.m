function lens = shortest_segment(intree)
    inv_dA = intree.dA';
    [row, col] = find(inv_dA);
    sects = [row, col]; % all segments between points
    plen = Pvec_tree(intree); % cumulative path lengths to root
    plen(sects);
    slen = diff(plen(sects), [], 2); % all segment lengths
%     lens = [];
%     for i = 1: 7
%         [min_len, id] = min(slen)
%         lens = [lens; min_len, id];
%         slen(id) = [];
%     end
    lens = min(slen(find(slen>0.0001)));
end