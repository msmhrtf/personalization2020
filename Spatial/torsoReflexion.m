function out = torsoReflexion( in, deltaT, thetaR, a, alfa_min, theta_min, frac, fs )
% torso reflexion path (without shoulder reflexion coefficient)
% 3 filters in cascade ---> TO DO: the whole process in a nicer and more optimal way.

% deltaT delay
delay_samp = deltaT*fs;
if frac
    % fractional delay 
    base_delay = floor(delay_samp);
    frac_delay = delay_samp - base_delay;
    h = dfilt.delay(base_delay);
    d = fdesign.fracdelay(frac_delay);
    ThirdOrderFrac = design(d,'lagrange','filterstructure','farrowfd');
    Hcas = dfilt.cascade(h,ThirdOrderFrac);
    out = filter(Hcas, in);
else
    % delay filter with rounding the delay to closest integer
    h = dfilt.delay(round(delay_samp));
    out = filter(h, in);
end

% head-shadow
out = head_shadow(out, thetaR, a, alfa_min, theta_min, fs);

% ITD
out = head_itd(out, thetaR, a, fs, frac);

end

