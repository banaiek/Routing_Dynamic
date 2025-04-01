%%--------------------
% This function calculate the spike-triggered MUA using a nonparametric
% approach for spikes in one channel and MUA in another channels

%%%%  Inputs:
% spike_In : a vector of spike indices (for a given trial)
% WB_In : wideband signal from a simultaneously recorded channel
% Par : Input parameters instead of default values:
% par = struct();
% par.Fs = 3E4; % Wideband signal sampling rate (default)
% par.fband = [800 3000]; %bandPass cutOffs (default)
% par.slide_Step = 1; % sliding time step in ms (default)
% par.window_Length = 100; % time around spikes in ms (default)
% par.window_Time = 2.5; % Window length to compute MUA around each step in
% ms (default)
% par.per_post_Window = 20;  % Size for pre and post in ms (default)

%%%%  Oututs:
% out.Lead_lag_Val : a vector of lead-lag values in (default) 20ms around spikes
% out.timpestamp : time vectore for spike_triggered_MUA
% out.spike_triggered_MUA : spike_triggered_MUA for spikes
    

%%%%%%%%%%%%%%------------------Written By
%%%% Kia Banaie Boroujeni, 2023

% Transition During Oscillatory Bursts and Attentional selection. Neuron(2023) 
% Banaie Boroujeni, K., Womelsdorf, T. Routing States 
% doi: https://doi.org/10.1016/j.neuron.2023.06.012
  

function [out] = spike_triggered_MUA(spike_In,wb_In,par)

if ~exist('Par', 'var')
    par = struct();
end

defaultValues = struct();
defaultValues.Fs = 3E4;
defaultValues.fband = [800 3000];
defaultValues.slide_Step = 1;
defaultValues.window_Length = 100;
defaultValues.window_Time = 2.5;
defaultValues.per_post_Window = 20;
defaultValues.smoothen = 0;

fieldNames = fieldnames(defaultValues);

for i = 1:length(fieldNames)
    fieldName = fieldNames{i};
    if ~isfield(par, fieldName)
        par = setfield(par, fieldName, getfield(defaultValues, fieldName));
    end
end

useGPU = false;
if gpuDeviceCount > 0
    useGPU = true;
    wb_In = gpuArray(wb_In);
    spike_In = gpuArray(spike_In);
end

Fs = par.Fs;
fband = par.fband;
slide_Step = par.slide_Step;
window_Length = par.window_Length;
window_Time = par.window_Time;
per_post_Window = par.per_post_Window;

timestamps = [-window_Length:slide_Step:window_Length];
windowed_Segments = timestamps*Fs/1000;
window_size = round(window_Time*Fs/1000);

pre_spike = timestamps < 0 & timestamps >= -per_post_Window;
post_spike = timestamps > 0 & timestamps <= per_post_Window;

Peaks_Val = nan(size(wb_In), 'like', wb_In);
MUA_Filtered = bandpass_filter(wb_In, fband, Fs);
MUA_Rectified = abs(MUA_Filtered);
MUA_Peaks = find(MUA_Rectified - circshift(MUA_Rectified,1) > 0 & MUA_Rectified - circshift(MUA_Rectified,-1) > 0);
Peaks_Val(MUA_Peaks) = MUA_Rectified(MUA_Peaks);

num_steps = length(windowed_Segments);
num_spikes = length(spike_In);

spk_triggered_MUA = nan(num_spikes, num_steps, 'like', wb_In);

for i = 1:length(spike_In)
    spike_Ind = spike_In(i);
    intervals = spike_Ind + windowed_Segments;
    count = 1;
    for j = 1:num_steps
        window_start = intervals(j) - window_size;
        window_end = intervals(j) + window_size;
        if window_start < 1 || window_end > length(Peaks_Val)
            spk_triggered_MUA(i, count) = NaN;
        else
            spk_triggered_MUA(i, count) = median(Peaks_Val(window_start:window_end), 'omitnan');
        end
        count = count + 1;
    end
end

if par.smoothen > 0
    spk_triggered_MUA = movmean(spk_triggered_MUA, par.smoothen, 2);
end

lead_lag_Val = (mean(spk_triggered_MUA(:, post_spike), 2, 'omitnan') - mean(spk_triggered_MUA(:, pre_spike),2, 'omitnan'));

out.lead_lag_Val = lead_lag_Val;
out.spike_triggered_MUA = spk_triggered_MUA;
out.timestamps = timestamps;

if useGPU
    out.lead_lag_Val = gather(out.lead_lag_Val);
    out.spike_triggered_MUA = gather(out.spike_triggered_MUA);
    out.timestamps = gather(out.timestamps);
end

end

function bp_data = bandpass_filter(signal, passband, fs)
[b, a] = butter(2, passband/(fs/2), 'bandpass');
if isa(signal, 'gpuArray')
    b = gpuArray(b);
    a = gpuArray(a);
end
bp_data = filtfilt(b, a, signal);
end
