% Generic hrtf + snowman to compensente at low frequencies
close all
% Paper: HRTF PERSONALIZATION USING ANTHROPOMETRIC MEASUREMENTS
% N. Zotkin et al. - 2003

% go to line 89 for algorithm

% dry violin sound
in = audioread('trvf-open.wav');
in = changeSamplingRate(in, 96000, 44100);
in = in(1:44100*3);

yes = 1; % plot
fs = 44100;
impulse = [1; zeros(199,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set position of source and listener and get hrtf from database

% set source and listener position in 3D space (SDN objects)
source = Source();
source.position = Position(200,3,3); % source far away on the right of the listener
microphone = Microphone(); 
microphone.position = Position(2,3,3);
microphone.heading = deg2rad(90); % listener heading to have source on the right side

S = [source.position.x, source.position.y, source.position.z]; % source point 
M = [microphone.position.x, microphone.position.y, microphone.position.z]; % center of head point

% get azimuth and elevation relative to source and mic positions and heading
azimuth = rad2deg(microphone.getAzimuth(source))
elevation = rad2deg(microphone.getElevation(source))

[hrir_left,hrir_right] = getGenericHRTF(source, microphone, fs);

% compute FFTs of HRIR (CIPIC function)
[hrtf_left, sampfreq_l] = freq_resp(hrir_left, 0, Inf); % [hrir_left; zeros(2^15,1)]
[hrtf_right, sampfreq_r] = freq_resp(hrir_right, 0, Inf); % [hrir_right; zeros(2^15,1)]

% Xmag1 = abs(fft(hrtf_left)); 
% w1 = [0:length(hrtf_left)-1].*fs/length(hrtf_left);
% figure
% plot(w1(1:length(w1)/2), 20*log10(Xmag1(1:length(Xmag1)/2))); grid on

% plots
if yes
    figure
    subplot(2,1,1); plot(hrir_left); hold on; plot(hrir_right); grid on
    subplot(2,1,2); plot(sampfreq_l, hrtf_left); hold on; plot(sampfreq_r, hrtf_right); grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate snowman

% snowman variables
frac = 0; % process with fractional delay
c = 343; % speed of sound
alfa_min = 0.1; % shadow cst
theta_min = 5*pi/6; % shadow cst (=150 degree)
theta_flat = theta_min*(0.5+1/pi*asin(alfa_min/(2-alfa_min))); % in rad

% anthropometric constants
a = 0.0875; % head radius in [m] a = 0.0875;
b = 0.169; % torso radius
h = 0.053; % neck length
rho = 0.3; % shoulder reflexion coefficient

% snowman torso position
B = [M(1), M(2), M(3)-(h+b)]; % center of torso point
torso = Position(B(1), B(2), B(3)); % torso center Position
if B(3)-b < 0 % check if snowman is in the space
    fprintf('The snowman torso is out of the box! \n')
end

% compute snowman
[hrirSnowman_left, hrirSnowman_right] = snowMan(impulse, fs, frac, alfa_min, theta_min, theta_flat, a, b, h, rho, source, microphone);

% snowman fft
[hrtfSnowman_left, sampfreq_L] = freq_resp(hrirSnowman_left, 0, Inf);
[hrtfSnowman_right, sampfreq_R] = freq_resp(hrirSnowman_right, 0, Inf);

if yes
    figure
    subplot(2,1,1); plot(hrirSnowman_left); hold on; plot(hrirSnowman_right); grid on
    subplot(2,1,2); plot(sampfreq_L, hrtfSnowman_left); hold on; plot(sampfreq_R, hrtfSnowman_right); grid on
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hybrid between measured HRTF and snowman

Hc = fft(hrir_right);
Hh = fft(hrirSnowman_right);
Hh_phase_right = angle(Hh);
Hc_phase = angle(Hc);
Ac = (abs(Hc));
Ah = (abs(Hh));
L = length(Ah);
nfreq = floor(L/2);      % DC is not counted
omega = ((0:nfreq)*(fs/(2*nfreq)))';
omega_low = 500;
omega_high = 3000;

Ac_cut = Ac(1:nfreq+1);       %%%%%%%          MAYBE NEED TO ADD +1 -> Ac(1:nfreq+1)
Ah_cut = Ah(1:nfreq+1);
% figure; plot(omega, Ac)

[m_low, idx_low] = min(abs(omega-omega_low));
[m_high, idx_high] =  min(abs(omega-omega_high));

% construct desired log magnitude
As = [Ah_cut(1:idx_low); zeros(idx_high-idx_low-1, 1); Ac_cut(idx_high:end)];
for i = idx_low+1:idx_high-1
    As(i) = Ah_cut(i) + (Ac_cut(i)-Ah_cut(i))/(omega(idx_high)-omega(idx_low));
end
         
As = [As; flipud(conj(As(2:end-1,:)))];

% minimum phase reconstruction   [0:L-1].*fs/L
% delay_s = 17;
% delay_t = delay_s / fs;
% As_mps = mps(As);
% As_mps_phase = angle(As_mps);
% f = ((0:length(As_mps)-1)*(fs/(length(As_mps))))';
% tau_phi = -2 * pi * delay_t * f;
% % ...
% Hs = As_mps;

Hs = (As).*exp(1i*Hh_phase_right); % snowman phase

y_R = ifft(Hs,'symmetric');

[y_fft, ss] = freq_resp(y_R, 0, Inf);

if yes
    figure
    plot(ss, y_fft)
    hold on
    plot(sampfreq_l, hrtf_right);
    legend('new hrtf', 'original hrtf')
end

%% left

Hc = fft(hrir_left);
Hh = fft(hrirSnowman_left);
Hh_phase_left = angle(Hh);
Hc_phase = angle(Hc);
Ac = (abs(Hc));
Ah = (abs(Hh));
L = length(Ah);
nfreq = floor(L/2);     % DC is not counted
omega = ((0:nfreq)*(fs/(2*nfreq)))';
omega_low = 500;
omega_high = 3000;

Ac_cut = Ac(1:nfreq+1);       %%%%%%%          MAYBE NEED TO ADD +1 -> Ac(1:nfreq+1)
Ah_cut = Ah(1:nfreq+1);
% figure; plot(omega, Ac)

[m_low, idx_low] = min(abs(omega-omega_low));
[m_high, idx_high] =  min(abs(omega-omega_high));

% construct desired log magnitude
As = [Ah_cut(1:idx_low); zeros(idx_high-idx_low-1, 1); Ac_cut(idx_high:end)];
for i = idx_low+1:idx_high-1
    As(i) = Ah_cut(i) + (Ac_cut(i)-Ah_cut(i))/(omega(idx_high)-omega(idx_low));
end

As = [As; flipud(conj(As(2:end-1,:)))];

% minimum phase reconstruction
[m, i1] = max(hrirSnowman_left);
[m, i2] = max(hrirSnowman_right);

% delay_s = abs(i1-i2)
% delay_t = delay_s / fs;
% As_mps = mps(As);
% As_mps_phase = angle(As_mps);
% f = ((0:length(As_mps)-1)*(fs/(length(As_mps))))';
% tau_phi = -2 * pi * delay_t * f;
% As_mps_phase = As_mps_phase + tau_phi;
% Hs = (As_mps).*exp(1i*As_mps_phase);

Hs = (As).*exp(1i*Hh_phase_left); % keeping the phase of the measurement hrtf, not the one of the snowman as desrcibed in the paper

y_L = ifft(Hs,'symmetric');

[y_fft, ss] = freq_resp(y_L, 0, Inf);

% if yes
%     figure
%     plot(ss, y_fft)
%     hold on
%     plot(sampfreq_l, hrtf_right);
%     legend('new hrtf', 'original hrtf')
% end

figure;
plot(y_L);
hold on;
plot(y_R); 
hold on;
plot(hrir_left);
hold on;
plot(hrir_right);
legend('myL', 'myR', 'gL', 'gR')


soundsc([convolveFFT(in, y_L), convolveFFT(in, y_R)], fs)

