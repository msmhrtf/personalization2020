
% dry violin sound
in = audioread('trvf-open.wav');
in = changeSamplingRate(in, 96000, 44100);
in = in(1:48000*15);

in = [1; zeros(199,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
frac = 0; % process with fractional delay
fs = 44100;
c = 343; % speed of sound
alfa_min = 0.1; % shadow cst
theta_min = 5*pi/6; % shadow cst (=150 degree)
theta_flat = theta_min*(0.5+1/pi*asin(alfa_min/(2-alfa_min))); % in rad

% anthropometric constants
a = 0.0775; % head radius in [m]
b = 0.169; % torso radius
h = 0.053; % neck length
rho = 0.3; % shoulder reflexion coefficient

% set source and listener position in 3D space (SDN objects)
source = Source();
source.position = Position(2,3,100); % source far away on the right of the listener
microphone = Microphone(); 
microphone.position = Position(2,3,3);
microphone.heading = deg2rad(90); % listener heading to have source on the right side

% initialise snowman positions
S = [source.position.x, source.position.y, source.position.z]; % source point 
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
shadow_limit_r = -sqrt(d_l_norm^2 - b^2);

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

soundsc([outL, outR], fs)





% % left-right ITDs
% T_l = (a+a*theta)/c; % left delay
% T_r = (a-a*sin(theta))/c; % right delay
% if theta < 0
%     temp = T_l; % interchange delay if negative azimuth
%     T_l = T_r;
%     T_r = temp;
% end
% 
% % head-shadow filter
% alfa = (1+alfa_min/2) + (1-alfa_min/2)*cos(theta/theta_min*pi);
% alfa_l = 1-sin(theta); 
% alfa_r = 1+sin(theta); % note: alfa_l + alfa_r = 2
% theta_flat = theta_min*(0.5+1/pi*sin(alfa_min/(2-alfa_min)));
% % shadow filter
% tau = 2*c/a;
% b_shadow_l = [2*alfa_l+tau, tau/fs-2*alfa_l]; % shadow filter coefficients
% b_shadow_r = [2*alfa_r+tau, tau/fs-2*alfa_r];
% a_shadow = [2+tau/fs, tau/fs-2]; % denominator is constant
