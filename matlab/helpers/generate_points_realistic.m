function [R1, R2, Ry, f, F, x1, x2, R, t, x1u, x2u] = generate_points_realistic(N, y_angle, distortion_param)
if nargin < 1
    N = 100;
end
if nargin < 2
    y_angle = 0;
end
if nargin < 3
    distortion_param = 0;
end

reproj = inf;
while norm(reproj) > 1e-12
    
    % Relative translation
    t = randn(3,1);
    t = t / norm(t) * 5.5 * 0.1; % The distance between the two cameras is set to be 10 percent of the average scene dept
    f = 300 + rand * 2700;
    K = diag([f, f, 1]);
    Kinv = diag([1 / f, 1 / f, 1]);
    
    % [R1, ~] = qr(randn(3));
    % [R2,  ~] = qr(randn(3));
    R1 = eul2rotm(deg2rad(30*randn(1, 3)), 'XYZ');
    R2 = eul2rotm(deg2rad(30*randn(1, 3)), 'XYZ');
    Ry = eul2rotm([0 y_angle 0], 'XYZ');
    
    R = R2 * Ry * R1';
    
    P1 = K * [R1 zeros(3,1)];
    P2 = K * [R2 * Ry t];
    
    % Fundamental matrix
    F = Kinv * skew(t) * R * Kinv;
    
    % Generate points
    X = [6*rand(2,N)-3; 5*rand(1,N)+3; ones(1,N)];
    
    % Generate point correspondences (pinhole)
    x1 = pflat(P1 * X);
    x2 = pflat(P2 * X);
    
    % Check GT
    reproj = diag(x2' * F * x1);
    
    % Add radial distortion (if desired)
    x1u = x1;
    x2u = x2;
    if distortion_param < 0
        x1 = radialdistort(x1, distortion_param);
        x2 = radialdistort(x2, distortion_param);       
    end
end