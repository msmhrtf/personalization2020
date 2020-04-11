
load('Perso_Kemar_060.mat'); % matrix_l_new & matrix_r_new: full personalisation of KEMAR 021 giving CIPIC subject 060
load('KEMAR_021_plus.mat'); % hrir_l & hrir_r: cipic subject 021 KEMAR with added points
load('KEMAR_021_perso.mat') % hrtf_l_perso & hrtf_r_perso: cipic subject 021 KEMAR personlaised ITD and ILD according to duda

fs = 44100;
L = 200;

%% ITD (KEMAR 021 vs KEMAR 021 ITD personalisation)
% ITD study in the horizontal plane

azimuths = [-90 -80 -65 -55 -45:5:45 55 65 80 90]; % azimuth 0 @ index = 14
elevations = [-90 -45 + 5.626*(0:49) 270]; % elevation 0 @ index = 10 / elevation 180 @ index = 42
% azimuths = [-80 -65 -55 -45:5:45 55 65 80]; 
% elevations = [-45 + 5.626*(0:49)];
load('horPlaneAzimuth.mat'); % azimuths2: full circle azimuths from 0 to 180 and back to 0
load('aziElePairs.mat'); % aziEle: cipic full circle azimuths with corresponding elevations

in = [1; zeros(L-1,1)];
onsets_generic = zeros(2, length(aziEle));
onsets_perso = zeros(2, length(aziEle));

for i = 1:length(aziEle)
    [~, az_idx] = min(abs(azimuths - aziEle(1,i)));
    [~, el_idx] = min(abs(elevations - aziEle(2,i)));
    onsets_generic(1,i) = hrirs_onset_by_crosscorrelation( squeeze(hrtf_l(az_idx, el_idx, :)), in );
    onsets_generic(2,i) = hrirs_onset_by_crosscorrelation( squeeze(hrtf_r(az_idx, el_idx, :)), in );
    onsets_perso(1,i) = hrirs_onset_by_crosscorrelation( squeeze(hrtf_l_perso(az_idx, el_idx, :)), in );
    onsets_perso(2,i) = hrirs_onset_by_crosscorrelation( squeeze(hrtf_r_perso(az_idx, el_idx, :)), in );
end

% itd in micro seconds
itd_g = (onsets_generic(1,:) - onsets_generic(2,:)) / fs * 1e6; 
itd_p = (onsets_perso(1,:) - onsets_perso(2,:)) / fs * 1e6;

[azimuths2, dum] = sort(azimuths2);
itd_g = itd_g(dum);
itd_p = itd_p(dum);

figure;
plot(azimuths2, itd_g, 'o-')
hold on
plot(azimuths2, itd_p, 'o-')
title('horizontal plane ITD: generic KEMAR vs modified ITD KEMAR')
xlabel('azimuth [degrees]')
ylabel('ITD [micro-second(s)]')
legend('generic KEMAR', 'modified ITD KEMAR')
xlim([-180 180])
ylim([-850 850])
xticks([-180 -90 0 90 180])
grid on

%% ILD horizontal plane (KEMAR 021 vs KEMAR 021 with modified ITD and ILD)
% ILD study in the horizontal plane 8 azimuth position around the listener

aziSelect = [-135, -90, -45, 0, 45, 90, 135, 180]; % azi vertical polar
aziEleSelect = [-45, -90, -45, 0, 45, 90,  45,   0; ... % azi interaural polar
                180,   0,   0, 0,  0,  0, 180, 180]; % ele interaural polar

% KEMAR 021 original
figure;

for i = 1:length(aziEleSelect)
    [~, az_idx] = min(abs(azimuths - aziEleSelect(1,i)));
    [~, el_idx] = min(abs(elevations - aziEleSelect(2,i)));
    x1 = [squeeze(hrtf_l(az_idx, el_idx, :)); zeros(2^15,1)];
    Xmag = abs(fft(x1));
    % Xmag = Xmag - max(Xmag);
    logMag = 20*log10(Xmag(1:length(Xmag)/2));
    
%     logMag = logMag - max(logMag);
%     logMag = logMag - logMag(1);
    w = ([0:length(x1)-1].*fs/length(x1))';
    plot(subplot(2,1,1), w(1:length(w)/2), logMag);
%     plot(w(1:749), 20*log10(Xmag(1:749))); % freq limit = 1000 Hz
%     plot(subplot(2,1,1), w(1:1123), 20*log10(Xmag(1:1123))); % freq limit = 1500 Hz
    hold on
end
grid on
title('Generic KEMAR, left ear: low-end response in the horizontal plane')
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')
xlim([20 20000])
% xlim([0 1000])
legend({['\theta= -135', char(176)], ['\theta= -90', char(176)], ['\theta= -45', char(176)], ['\theta= 0', char(176)], ['\theta= 45', char(176)], ['\theta = 90', char(176)], ['\theta = 135', char(176)], ['\theta=180', char(176)]})


% KEMAR 021 perso

for i = 1:length(aziEleSelect)
    [~, az_idx] = min(abs(azimuths - aziEleSelect(1,i)));
    [~, el_idx] = min(abs(elevations - aziEleSelect(2,i)));
    x1 = [squeeze(hrtf_l_perso(az_idx, el_idx, :)); zeros(2^15,1)];
    Xmag = abs(fft(x1));
    logMag = 20*log10(Xmag(1:length(Xmag)/2));
%     logMag = logMag - max(logMag);
%     logMag = logMag - logMag(1);
    w = ([0:length(x1)-1].*fs/length(x1))';
    plot(subplot(2,1,2), w(1:length(w)/2), logMag);
%     plot(w(1:749), 20*log10(Xmag(1:749))); % freq limit = 1000 Hz
%     plot(subplot(2,1,2), w(1:1123), 20*log10(Xmag(1:1123))); % freq limit = 1500 Hz
    hold on
end
grid on
title('Personalised KEMAR, left ear: low-end response in the horizontal plane')
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')
xlim([20 20000])
% xlim([0 1000])
legend({['\theta=-135', char(176)], ['\theta=-90', char(176)], ['\theta=-45', char(176)], ['\theta=0', char(176)], ['\theta=45', char(176)], ['\theta=90', char(176)], ['\theta=135', char(176)], ['\theta=180', char(176)]})

%% Median Plane (KEMAR 021 vs Personalisation of KEMAR -> subject 060)
% Median plane study for azimuth 0 and elevation -45-45.
% Frequency 1 kHz to 15 kHz.
% No intermpolation in between.

az_idx = 14; % azi 0
eleSelect = -45 + 5.626*(0:16);
eleMin = -45;
eleMax = 45;
freqMin = 1000;
freqMax = 15000;
magMin = -30;
magMax = 20;


% KEMAR 021 generic ------------------------------------------------------

% figure;
for i = 1:length(eleSelect)
    [~, el_idx] = min(abs(elevations - eleSelect(i)));
    x1 = [squeeze(hrtf_l(az_idx, el_idx, :)); zeros(2^15,1)];
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

% figure
for i = 1:length(eleSelect)
    [~, el_idx] = min(abs(elevations - eleSelect(i)));
    x1 = [squeeze(matrix_l_new(az_idx, el_idx, :)); zeros(2^15,1)];
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
ylabel('Frequency (kHz)');
zlabel('Magnitude (dB)');
c.Label.String = 'Magnitude (dB)';
title('CIPIC 060: Personalisation for KEMAR')

%% Horizontal Plane (KEMAR 021 vs Personalisation of KEMAR -> subject 060)


el_idx = 10; % ele 0
aziMin = -90;
aziMax = 90;
freqMin = 0;
freqMax = 1000;
magMin = -15;
magMax = 2;


% KEMAR 021 generic ------------------------------------------------------

% figure;
for i = 1:length(azimuths)
    x1 = [squeeze(hrtf_l(i, el_idx, :)); zeros(2^15,1)];
    Xmag(:,i) = abs(fft(x1));
    w(:,i) = [0:length(x1)-1].*fs/length(x1);
%     plot(w(1:length(w)/2, i), 20*log10(Xmag(1:length(Xmag)/2, i)));
%     hold on
end

figure
subplot(1,2,1)
surf(azimuths', w(1:749, :), 20*log10(Xmag(1:749, :)))
% set(gca, 'YScale', 'log');
% yticks([2000 4000 6000 8000 10000 12000 14000])
% yticklabels({'2','4','6','8','10','12','14'})
view(2)
caxis([magMin magMax]);
xlim([aziMin aziMax]);
ylim([freqMin freqMax]);
% zlim([magMin magMax]);
shading flat;
colormap('hot');
c = colorbar;
axis tight;
xlabel('Azimuth (deg)');
ylabel('Frequency (Hz)');
zlabel('Magnitude (dB)');
c.Label.String = 'Magnitude (dB)';
title('CIPIC 021: KEMAR, left ear, horizontal plane')


% Perso of KEMAR -> subject 060 -------------------------------------------

% figure;
for i = 1:length(azimuths)
    x1 = [squeeze(matrix_l_new(i, el_idx, :)); zeros(2^15,1)];
    Xmag(:,i) = abs(fft(x1));
    w(:,i) = [0:length(x1)-1].*fs/length(x1);
%     plot(w(1:length(w)/2, i), 20*log10(Xmag(1:length(Xmag)/2, i)));
%     hold on
end

subplot(1,2,2)
surf(azimuths', w(1:749, :), 20*log10(Xmag(1:749, :)))
% set(gca, 'YScale', 'log');
% yticks([2000 4000 6000 8000 10000 12000 14000])
% yticklabels({'2','4','6','8','10','12','14'})
view(2)
caxis([magMin magMax]);
xlim([aziMin aziMax]);
ylim([freqMin freqMax]);
% zlim([magMin magMax]);
shading flat;
colormap('hot');
c = colorbar;
axis tight;
xlabel('Azimuth (deg)');
ylabel('Frequency (Hz)');
zlabel('Magnitude (dB)');
c.Label.String = 'Magnitude (dB)';
title('Personalised CIPIC 060: left ear, horizontal plane')
