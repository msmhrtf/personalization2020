function out = torsoShadowSub(in, M, D, s, d, a, b, alfa_min, theta_min, theta_flat, fs, frac)
% TORSO SHADOW SUB-MODEL PROCESSING
% in        : signal input
% M         : position listener
% D         : position ear
% s         : vector center of torso to source
% d         : vector center of torso to ear
% a         : radius head
% b         : radius torso
% alfa_min  : algo constant
% theta_min : algo constant
% theta_flat: algo constant
% fs        : sampling frequency
% frac      : boolean for processing with fractional delay

d_norm = norm(d);

% torso shadow variables
zeta = acos(dot(s,d)/d_norm);
zeta_min = pi/2 + acos(b/d_norm);
theta_T = (pi*(zeta-zeta_min) - theta_flat*(zeta-pi)) / (pi-zeta_min); % interpolation

% ITD angle
e = D - M;
theta_HT = acos(dot(s,e)/a);

% head-shadow angle
w1 = (b/d_norm)^2;
w2 = sqrt( w1*(1-w1) / (d_norm^2-(dot(d,s))^2) );
d2 = d_norm^2*s - (dot(d,s)*d);
b_vec = w1*d + w2*d2;
r = (b_vec-d)/norm(b_vec-d);
theta_HS = acos(dot(r,e)/a);

% torso shadow processing
out = torsoShadow(in, theta_T, theta_HS, theta_HT, a, b, alfa_min, theta_min, frac, fs);

end
