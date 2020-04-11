function [hrir_l, hrir_r] = genericPlusSnowman(source, microphone, a, b, h, rho, fs)

% Generic hrtf + snowman to compensente at low frequencies

% Paper: HRTF PERSONALIZATION USING ANTHROPOMETRIC MEASUREMENTS - Dniitry
% N. Zotkin et al. - 2003

% go to line 89 for algorithm

yes = 0; % plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set position of source and listener and get hrtf from database

M = [microphone.position.x, microphone.position.y, microphone.position.z]; % center of head point

% get azimuth and elevation relative to source and mic positions and heading
azimuth = microphone.getAzimuthInteraural(source);
elevation = microphone.getElevationInteraural(source);

[hrir_left, hrir_right] = getGenericHRTF(source, microphone, fs);

% plots
if yes
    % compute FFTs of HRIR (CIPIC function)
    [hrtf_left, sampfreq_l] = freq_resp(hrir_left, 0, Inf); % [hrir_left; zeros(2^15,1)]
    [hrtf_right, sampfreq_r] = freq_resp(hrir_right, 0, Inf); % [hrir_right; zeros(2^15,1)]
    figure
    subplot(2,1,1); plot(hrir_left); hold on; plot(hrir_right); grid on
    subplot(2,1,2); plot(sampfreq_l, hrtf_left); hold on; plot(sampfreq_r, hrtf_right); grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate snowman

% snowman variables
frac = 1; % process with fractional delay
c = 343; % speed of sound
alfa_min = 0.1; % shadow cst
theta_min = 5*pi/6; % shadow cst (=150 degree)
theta_flat = theta_min*(0.5+1/pi*asin(alfa_min/(2-alfa_min))); % in rad

% snowman torso position
B = [M(1), M(2), M(3)-(h+b)]; % center of torso point
torso = Position(B(1), B(2), B(3)); % torso center Position
if B(3)-b < 0 % check if snowman is in the space
    fprintf('The snowman torso is out of the box! \n')
end

L = 200;
% compute snowman
impulse = [1; zeros(L-1,1)];
[hrirSnowman_left, hrirSnowman_right] = snowMan(impulse, fs, frac, alfa_min, theta_min, theta_flat, a, b, h, rho, source, microphone);

if yes
    % snowman fft
    [hrtfSnowman_left, sampfreq_L] = freq_resp(hrirSnowman_left, 0, Inf);
    [hrtfSnowman_right, sampfreq_R] = freq_resp(hrirSnowman_right, 0, Inf);
    figure
    subplot(2,1,1); plot(hrirSnowman_left); hold on; plot(hrirSnowman_right); grid on
    subplot(2,1,2); plot(sampfreq_L, hrtfSnowman_left); hold on; plot(sampfreq_R, hrtfSnowman_right); grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hybrid between measured HRTF and snowman

Hc_r = fft(hrir_right);
Hh_r = fft(hrirSnowman_right);
Hc_phase_r = angle(Hc_r);
Hh_phase_r = angle(Hh_r);
Ac_r = log10(abs(Hc_r));
Ah_r = log10(abs(Hh_r));

Hc_l = fft(hrir_left);
Hh_l = fft(hrirSnowman_left);
Hc_phase_l = angle(Hc_l);
Hh_phase_l = angle(Hh_l);
Ac_l = log10(abs(Hc_l));
Ah_l = log10(abs(Hh_l));

% L = length(Ah_r);
nfreq = floor(L/2);     % DC is not counted
omega = ((0:nfreq)*(fs/(2*nfreq)))';
omega_low = 500;
omega_high = 3000;

Ac_cut_r = Ac_r(1:nfreq+1);
Ah_cut_r = Ah_r(1:nfreq+1);
Ac_cut_l = Ac_l(1:nfreq+1);
Ah_cut_l = Ah_l(1:nfreq+1);

[m_low, idx_low] = min(abs(omega-omega_low));
[m_high, idx_high] =  min(abs(omega-omega_high));

% construct desired log magnitude
As_r = [Ah_cut_r(1:idx_low); zeros(idx_high-idx_low-1, 1); Ac_cut_r(idx_high:end)];
for i = idx_low+1:idx_high-1
    As_r(i) = Ah_cut_r(i) + (Ac_cut_r(i)-Ah_cut_r(i))/(omega(idx_high)-omega(idx_low));
end

As_l = [Ah_cut_l(1:idx_low); zeros(idx_high-idx_low-1, 1); Ac_cut_l(idx_high:end)];
for i = idx_low+1:idx_high-1
    As_l(i) = Ah_cut_l(i) + (Ac_cut_l(i)-Ah_cut_l(i))/(omega(idx_high)-omega(idx_low));
end
         
As_r = [As_r; flipud(conj(As_r(2:end-1,:)))];
As_l = [As_l; flipud(conj(As_l(2:end-1,:)))];

% minimum phase reconstruction
[m, i1] = max(hrirSnowman_left);
[m, i2] = max(hrirSnowman_right);
delay_s = abs(i1-i2);
delay_t = delay_s / fs;
As_mps_l = mps(As_l);
As_mps_r = mps(As_r);
f = ((0:length(As_mps_l)-1)*(fs/(length(As_mps_l))))';
tau_phi = -2 * pi * delay_t * f;
if azimuth < 0
    As_mps_phase = angle(As_mps_r);
    As_mps_phase = As_mps_phase + tau_phi;
    Hs_r = (As_mps_r).*exp(1i*As_mps_phase);
    Hs_l = As_mps_l;
else
    As_mps_phase = angle(As_mps_l);
    As_mps_phase = As_mps_phase + tau_phi;
    Hs_r = As_mps_r;
    Hs_l = (As_mps_l).*exp(1i*As_mps_phase);
end

% !!!!!!!!!! Need to keep the phase of the snowman to keep the personalised ITD.
% Hs_r = (10.^As_r).*exp(1i*Hh_phase_r); % keeping the phase of the model hrtf
% Hs_l = (10.^As_l).*exp(1i*Hh_phase_l); 

hrir_r = ifft(Hs_r,'symmetric');
hrir_l = ifft(Hs_l,'symmetric');

if yes
    [y_fft, ss] = freq_resp(hrir_r, 0, Inf);
    figure
    plot(ss, y_fft)
    hold on
    plot(sampfreq_l, hrtf_right);
    legend('new hrtf', 'original hrtf')
end


end

