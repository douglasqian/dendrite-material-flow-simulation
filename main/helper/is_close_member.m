function [close, id] = is_close_member(pos, points, sradius, target) % search radius determined by step size
    close = 0;
    for i = 1: length(points)
        x1 = pos(1);
        y1 = pos(2);
        x2 = points(i, 1);
        y2 = points(i, 2);
        id = points(i, 3);
        % only sets close to 1 if the point is also a downstream point (child of target)
        if pdist([x1,y1; x2,y2], 'euclidean') < sradius && (isnan(target) || id == target)
            close = 1;
            break
        end
    end
end
