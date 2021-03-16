function scale = normalise2dpts(pts)

    dist = sqrt(pts(1,:).^2 + pts(2,:).^2);
    meandist = mean(dist(:));   
    scale = sqrt(2)/meandist;