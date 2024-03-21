function [xSignal, xNoise] = extractSN(x, fs, sigStart, sigStop, noiseDist, clipBufferSize, dFilter, units)
% Isolate signal and associated noise samples from a larger audio time
% series vector, given a pre-determined signal location.
%
% Last updated by Wilfried Beslin
% 2024-03-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - Things I might do:
%   -- combine sigStart and sigStop into a single variable (makes for less
%   documentation and fewer input arguments
%   -- add input parsing with inputParser
%   -- add noise range as an output argument (in samples)

    
    import bandpass.noDelayFilt

    % get sigStart, sigStop, and noiseDist as samples, based on "units"
    switch units
        case 'seconds'
            samplesFromInput = @(a) round(a*fs);
        case 'samples'
            samplesFromInput = @(a) a;
        otherwise
            error('Invalid value for "units": must be either ''seconds'' or ''samples''.')
    end
    sigStartSample = samplesFromInput(sigStart);
    sigStopSample = samplesFromInput(sigStop);
    noiseDistSamples = samplesFromInput(noiseDist);
    clipBufferSamples = samplesFromInput(clipBufferSize);
    
    % get signal and noise size
    nSigSamples = sigStopSample - sigStartSample + 1;
    
    % generate shorter clip (easier to filter)
    %%% if it's not possible to generate the clip (i.e., because the signal
    %%% is too close to the beginning of the sequence), then return empties
    clipStartSample = sigStartSample - noiseDistSamples - nSigSamples - clipBufferSamples;
    clipStopSample = sigStopSample + clipBufferSamples;
    
    if clipStartSample > 0 && clipStopSample <= numel(x)
        xClip = x(clipStartSample:clipStopSample);

        % get relative signal and noise samples
        sigStartSampleClip = clipBufferSamples + nSigSamples + noiseDistSamples + 1;
        sigStopSampleClip = sigStartSampleClip + nSigSamples - 1;
        noiseStartSampleClip = sigStartSampleClip - noiseDistSamples - nSigSamples;
        noiseStopSampleClip = sigStartSampleClip - noiseDistSamples - 1;

        % apply digital filter
        xClipFilt = noDelayFilt(dFilter, xClip);

        % isolate signal
        xSignal = xClipFilt(sigStartSampleClip:sigStopSampleClip);

        % isolate noise
        xNoise = xClipFilt(noiseStartSampleClip:noiseStopSampleClip);

        %** DEBUG PLOT
        %{
        tClip = (((1:numel(xClip))-1)./fs)';
        figure;
        ax = axes();
        ax.NextPlot = 'add';
        plot(ax, tClip, xClip)
        plot(ax, tClip, xClipFilt)
        plot(ax, tClip(sigStartSampleClip:sigStopSampleClip), xSignal);
        plot(ax, tClip(noiseStartSampleClip:noiseStopSampleClip), xNoise);
        legend(ax, {'Unfiltered Clip', 'Filtered Clip', 'Signal', 'Noise'})
        grid(ax, 'on')
        box(ax, 'on')
        keyboard
        %}
    else
        xSignal = [];
        xNoise = [];
    end
end