% SIMPLIFY_TREE  Takes a tree and simplifies it (B,T points only)
% 
% simp_tree = simplify_tree(intree)
% -------------------------------------------------------
% reconstructs tree without the continuation points
% does not preserve length
% 
% Input
% ----
% 
% - intree:: struct: structured tree
% 
% Output
% -----
% 
% - simp_tree:: struct: structured tree

function stree = simplify_tree(intree)
    sections = dissect_tree(intree);
    newdim = length(sections);

    stree.dA = sparse(newdim, newdim);

    for id = 2: newdim
        pid = find(sections(:,2) == sections(id,1));
        stree.dA(id, pid) = 1;
    end
    stree.name = intree.name;
    stree.rnames = intree.rnames;
    stree.X = intree.X(sections(:,2));
    stree.Y = intree.Y(sections(:,2));
    stree.Z = intree.Z(sections(:,2));
    stree.D = intree.D(sections(:,2));
    stree.R = intree.R(sections(:,2));
end