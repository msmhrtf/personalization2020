function [ onset ] = hrirs_onset_by_crosscorrelation( data, signal )
    %% data = hrirs_onset_by_crosscorrelation( data, signal )
    % Compute the onset using the cross correlation between a raw data and
    % the stimolous signal
    %
	%
    % Input:
    %       data: the raw input data
    %       signal: the stimolous signal
    % Output:
    %       onset: the detected onset
    %
    % @author	Michele Marostica
    % @email	michelemaro@gmail.com
    % @insitute	Sound and Music Computing Group, Department of Information Engineering, University of Padua
    % @date		27th Jun 2013
    
    
    [c lags] = xcorr(data, signal);  %cross-correlation of L-canal to find the onset

    [a, pos] = max(c);

    onset = lags(pos);
            
end

