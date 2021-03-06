function [minlen, minid] = shortest_segment(intree)
    inv_dA = intree.dA';
    [row, col] = find(inv_dA);
    sects = [row, col]; % all segments between points
    plen = Pvec_tree(intree); % cumulative path lengths to root
    plen(sects);
    slen = diff(plen(sects), [], 2); % all segment lengths
    minlen = min(slen);
    minid = sects(find(slen == minlen), 2);
end