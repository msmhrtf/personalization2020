
az_idx = 14; % azi 0
eleSelect = -45 + 5.626*(0:16);
eleMin = -45;
eleMax = 45;
freqMin = 1000;
freqMax = 15000;
magMin = -30;
magMax = 20;

azimuths = [-90 -80 -65 -55 -45:5:45 55 65 80 90]; % azimuth 0 @ index = 14
elevations = [-90 -45 + 5.626*(0:49) 270]; % elevation 0 @ index = 10 / elevation 180 @ index = 42
fs = 44100;
L = 200;
% KEMAR 021 generic ------------------------------------------------------

load('C:\Users\jtissi16\Documents\MATLAB\TH\CIPIC\CIPIC_hrtf_database\standard_hrir_database\subject_021\hrir_final.mat')

% figure;
for i = 1:length(eleSelect)
    [~, el_idx] = min(abs(elevations - eleSelect(i)));
    x1 = [squeeze(hrir_l(az_idx, el_idx, :)); zeros(2^15,1)];
    Xmag(:,i) = abs(fft(x1));
    w(:,i) = [0:length(x1)-1].*fs/length(x1);
%     plot(w(1:length(w)/2, i), 20*log10(Xmag(1:length(Xmag)/2, i)));
%     hold on
end

figure
subplot(1,2,1)
surf(eleSelect', w(748:11215, :), 20*log10(Xmag(748:11215, :)))
% set(gca, 'YScale', 'log');
yticks([2000 4000 6000 8000 10000 12000 14000])
yticklabels({'2','4','6','8','10','12','14'})
view(2)
caxis([magMin magMax]);
xlim([eleMin eleMax]);
ylim([freqMin freqMax]);
% zlim([magMin magMax]);
shading flat;
colormap('hot');
c = colorbar;
axis tight;
xlabel('Elevation (deg)');
ylabel('Frequency (kHz)');
zlabel('Magnitude (dB)');
c.Label.String = 'Magnitude (dB)';
title('CIPIC 021: KEMAR, left ear, median plane')


% Perso of KEMAR -> subject 060 -------------------------------------------

load('C:\Users\jtissi16\Documents\MATLAB\TH\CIPIC\CIPIC_hrtf_database\standard_hrir_database\subject_060\hrir_final.mat')

% figure
for i = 1:length(eleSelect)
    [~, el_idx] = min(abs(elevations - eleSelect(i)));
    x1 = [squeeze(hrir_l(az_idx, el_idx, :)); zeros(2^15,1)];
    Xmag(:,i) = abs(fft(x1));
    w(:,i) = [0:length(x1)-1].*fs/length(x1);
%     plot(w(1:length(w)/2, i), 20*log10(Xmag(1:length(Xmag)/2, i)));
%     hold on
end

subplot(1,2,2)
surf(eleSelect', w(748:11215, :), 20*log10(Xmag(748:11215, :)))
% set(gca, 'YScale', 'log');
yticks([2000 4000 6000 8000 10000 12000 14000])
yticklabels({'2','4','6','8','10','12','14'})
view(2)
caxis([magMin magMax]);
xlim([eleMin eleMax]);
ylim([freqMin freqMax]);
% zlim([magMin magMax]);
shading flat;
colormap('hot');
c = colorbar;
axis tight;
xlabel('Elevation (deg)');
ylabel('Frequency (Hz)');
zlabel('Magnitude (dB)');
c.Label.String = 'Magnitude (dB)';
title('CIPIC 060: Personalisation for KEMAR')