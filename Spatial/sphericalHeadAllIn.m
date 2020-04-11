function [out_L, out_R, delay_L, delay_R] = sphericalHeadAllIn(az, el, a, L, fs, frac, shift_ears)
% in       : signal input
% M        : position listener
% D        : position ear
% s        : vector center of head to source
% a        : radius head
% alfa_min : algo constant
% theta_min: algo constant
% fs       : sampling frequency
% frac     : boolean for processing with fractional delay

in = [1; zeros(L,1)];
alfa_min = 0.1; % shadow cst
theta_min = 5*pi/6; % shadow cst (=150 degree)

if (shift_ears)
    e_b = 0.0094; % ear front shift
    e_d = 0.021; % ear down shift
    [D_l(1), D_l(2), D_l(3)] = sph2cart(atan(a/e_b), atan(a/e_d) - pi/2, a);
    [D_r(1), D_r(2), D_r(3)] = sph2cart(-atan(a/e_b), atan(a/e_d) - pi/2, a);
else
    D_l = [0, a, 0];
    D_r = [0, -a, 0];
end
M = zeros(1, 3);
[x, y, z] = sph2cart(-deg2rad(az), deg2rad(el), 1);
% remove noise below eps
x(abs(x)<eps)=0;
y(abs(y)<eps)=0;
z(abs(z)<eps)=0;
S = [x, y, z];
s = (S - M) / norm(S - M); 

% left
e = D_l - M; % vector from head center to ear
theta_D = acos(dot(s, e)/a);
direct_out = head_shadow(in, theta_D, a, alfa_min, theta_min, fs);
[out_L, delay_L] = head_itd(direct_out, theta_D, a, fs, frac);

% right
e = D_r - M; % vector from head center to ear
theta_D = acos(dot(s, e)/a);
direct_out = head_shadow(in, theta_D, a, alfa_min, theta_min, fs);
[out_R, delay_R] = head_itd(direct_out, theta_D, a, fs, frac);

end