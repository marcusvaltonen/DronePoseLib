function y = radialundistort(x, kappa)
% y = radialundistort(x, kappa) removes radial distortion from the
% (homogeneous or inhomogeneous) coordinates x using the parameter kappa
% (should be negative) with the model y = (1 + kappa*|y|^2)*x. for
% coordinates in [-1,1]^2 kappa=-.01 is a mild distortion and kappa=-.5
% is a pretty heavy distortion.

ishom = (size(x, 1) == 3);
if(ishom)
    x = pflat(x);
    x = x(1:2, :);
end

% compute distorted radius
rd2 = sum(x.^2);
rd = sqrt(rd2); 


% compute undistorted coordinates
y = [x; ones(1, size(x, 2)) + kappa*rd2];
y = pflat(y);

if(~ishom)
    y = y(1:2, :);
end
