close all

fs = 44100;

hrir = load('hrir_final.mat');

hrir_l = squeeze(hrir.hrir_l(4,5,:));
figure; plot(hrir_l); hold on;
[~,delay_samp] = max(hrir_l);
delay_t = (delay_samp-1) / fs;

hrtf_mps = mps(abs(fft(hrir_l)));
hrir_mps = ifft(hrtf_mps);

%plot(hrir_mps);
%plot([zeros(29,1);hrir_mps])

hrtf_mps_phase = angle(hrtf_mps);
hrtf_mps_mag = abs(hrtf_mps);
f = ((0:length(hrtf_mps)-1)*(fs/(length(hrtf_mps))))';
tau_phi = -2 * pi * delay_t * f;
new_phase = hrtf_mps_phase + tau_phi;

% hrtf_reconstruct = (hrtf_mps_mag).*exp(1i*new_phase);
hrtf_reconstruct = hrtf_mps_mag .* (cos(new_phase) + 1i * sin(new_phase));
hrir_reconstruct = ifft(hrtf_reconstruct);

plot(real(hrir_reconstruct));





