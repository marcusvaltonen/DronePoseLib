function R = rotm(axis, angle)

if strcmp(axis, 'x')
    R = [1 0 0;
        0 cos(angle) -sin(angle);
        0 sin(angle) cos(angle)];
elseif strcmp(axis, 'y')
    R = [cos(angle) 0 sin(angle);
        0 1 0;
        -sin(angle) 0 cos(angle)];
elseif strcmp(axis, 'z')
    R = [cos(angle) -sin(angle) 0;
        sin(angle) cos(angle) 0;
        0 0 1];
else
    error('Incorrect axis specified.')
end