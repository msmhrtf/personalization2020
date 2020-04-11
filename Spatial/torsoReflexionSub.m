function out = torsoReflexionSub(in, M, D, s, d, a, b, c, rho, alfa_min, theta_min, fs, frac)
% TORSO REFLEXION SUB-MODEL PROCESSING
% in       : signal input
% M        : position listener
% D        : position ear
% s        : vector center of torso to source
% d        : vector center of torso to ear
% a        : radius head
% b        : radius torso
% c        : speed of sound
% rho      : shoulder reflexion coefficient
% alfa_min : algo constant
% theta_min: algo constant
% fs       : sampling frequency
% frac     : boolean for processing with fractional delay

d_norm = norm(d); % norm of vector d (from center torso to ear)

% observation angle for direct path -> equation 7
e = D - M; % vector from head center to ear
theta_D = acos(dot(s,e)/a);

% direct path processing
direct_out = head_shadow(in, theta_D, a, alfa_min, theta_min, fs);
direct_out = head_itd(direct_out, theta_D, a, fs, frac);

% shoulder reflexion variable initialization (Appendix A)
A = d_norm/b;
alfa0 = (A-1)/(2*A-1)*(pi/2);
alfa_max = acos(1/A);
epsilon = pi/2 - acos(dot(d,s)/d_norm);
if epsilon < 0
    alfa = alfa0 - (1-alfa0/alfa_max)*epsilon;
else
    alfa = alfa0 * (1-epsilon/(pi/2));
end
% NOTE:
%%% I don't know what to do with the second epsilon. Should I reiterate the calculation of alfa and find a convergence? 
% epsilon_2 = pi/2-2*alfa-atan(sin(alfa)/(A-cos(alfa)));
% epsilon_2_r = pi/2-2*alfa_r-atan(sin(alfa_r)/(A_r-cos(alfa_r)));
%%%
f = sqrt(b^2 + d_norm^2 - 2*b*d_norm*cos(alfa));
beta = atan(sin(alfa)/(A-cos(alfa)));
psi = alfa + beta;
delta_T = f/c*(1+cos(2*psi)); % additional delays from the shoulder
% observation angle theta_R
d2 = d_norm^2*s - (dot(d,s)*d);
b = b*cos(alfa)*d/d_norm + b*sin(alfa)*d2/norm(d2);
r = (b-d)/norm(b-d);
theta_R = acos(dot(r,e)/a);

% torso reflexion processing
torso_ref = torsoReflexion( in, delta_T, theta_R, a, alfa_min, theta_min, frac, fs );

% torso reflection sub-model (combination of direct path and torso reflexion)
out = ( direct_out + rho*torso_ref ) * 1/(1+rho);

% out = direct_out;

end

