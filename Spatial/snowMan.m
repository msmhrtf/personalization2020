function [outL, outR] = snowMan(in, fs, frac, alfa_min, theta_min, theta_flat, a, b, h, rho, source, microphone)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
c = 343; % speed of sound

% initialise snowman positions
M = [microphone.position.x, microphone.position.y, microphone.position.z]; % center of head point
B = [M(1), M(2), M(3)-(h+b)]; % center of torso point
torso = Position(B(1), B(2), B(3)); % torso center Position
if B(3)-b < 0 % check if snowman is in the space
    fprintf('The snowman is out of the box! \n')
end

% find position of the ears D_l and D_r
ears_position = microphone.getEarsPosition(a); % ears at 90 degrees left and right relative to head "heading" (microphone heading)
ear_left = ears_position{1}; % left ear Position
ear_right = ears_position{2}; % right ear Position
D_l = [ear_left.x, ear_left.y, ear_left.z]; % left ear point
D_r = [ear_right.x, ear_right.y, ear_right.z]; % right ear point

% vectors needed for snowman
d_l = D_l - B; % vectors from torso center to ears
d_r = D_r - B;
d_l_norm = norm(d_l); % vector norms
d_r_norm = norm(d_r);
% s = Position.getUnitDirectionVector(microphone.position, source.position); % unit vector from head center to source
s = Position.getUnitDirectionVector(torso, source.position); % unit vector from torso center to source

%%% note for vector s : it is assumed that the source is infinitly far from
%%% the listener, meaning that the angle of incidence is the same for the
%%% torso and the head. TO DO: manage the case where the source is close to
%%% the listener.

% torso-shadow cone limit from left and right ear
shadow_limit_l = -sqrt(d_l_norm^2 - b^2);
shadow_limit_r = -sqrt(d_r_norm^2 - b^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE MANAGEMENT (in or out of the shadow cone for each ear, making 4 cases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% left ear
ds_dot_l = dot(d_l,s);
if ds_dot_l < shadow_limit_l
    outL = torsoShadowSub(in, M, D_l, s, d_l, a, b, alfa_min, theta_min, theta_flat, fs, frac);
else
    outL = torsoReflexionSub(in, M, D_l, s, d_l, a, b, c, rho, alfa_min, theta_min, fs, frac);
end

% right ear
ds_dot_r = dot(d_r,s);
if ds_dot_r < shadow_limit_r
    outR = torsoShadowSub(in, M, D_r, s, d_r, a, b, alfa_min, theta_min, theta_flat, fs, frac);
else
    outR = torsoReflexionSub(in, M, D_r, s, d_r, a, b, c, rho, alfa_min, theta_min, fs, frac);
end

end

