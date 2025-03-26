function out = detect_burst(X, freq_band, Num_cycle, separating_cycle, Fs, varargin)
% Inputs:
%   X                  - Input signal (vector)
%   freq_band          - Frequency band [f1 f2] in Hz
%   Num_cycle          - Minimum number of cycles for a burst (e.g., 1.5)
%   separating_cycle   - Minimum number of cycles separating bursts (e.g., 5)
%   Fs                 - Sampling frequency in Hz
%
% Name-Value Pair Arguments:
%   'FilterOrder'        - Order of the Butterworth filter (default: 4)
%   'BurstDurationCycles'- Significant burst duration in cycles (default: 2.5)
%
% Outputs:
%   out - Structure containing:
%     - burst_center_points            : Indices of burst center points
%     - burst_phases                   : Phases at burst center points
%     - burst_start_end_points         : Nx2 array of start and end indices of bursts
%     - burst_start_end_zero_crossings : Nx2 array of adjusted start and end indices at zero crossings
%     - band_passed_signal             : Band-passed filtered signal
%     - burst_signal                   : Signal with non-burst values set to NaN
%     - burst_signal_zc                : Signal with non-burst values set to NaN (adjusted to zero crossings)
%     - burst_density                  : Burst density over time (vector same size as input signal)
%
% This function implements an adaptive burst detection method that includes burst density calculation 
% and boundary adjustments to the nearest zero crossings.

% Parse Optional Inputs
p = inputParser;
addParameter(p, 'FilterOrder', 4, @(x) isnumeric(x) && mod(x,2)==0);
addParameter(p, 'BurstDurationCycles', 2.5, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});

filter_order = p.Results.FilterOrder;
burst_duration_cycles = p.Results.BurstDurationCycles;

% Filter the Signal
[b, a] = butter(filter_order / 2, freq_band / (Fs / 2), 'bandpass');
bp_signal = filtfilt(b, a, X);

% Get the Analytic Signal, Amplitude, and Phase
analytic_signal = hilbert(bp_signal);
a_t = abs(analytic_signal);         % Instantaneous amplitude
theta_t = angle(analytic_signal);     % Instantaneous phase

% Compute the Signal's Energy
IE_t = a_t .^ 2;
IE_mean = mean(IE_t);

% Set an Adaptive Energy Threshold
IE_median = median(IE_t);
IE_lower = IE_t(IE_t < IE_median);
std_half_normal = std(IE_lower) / sqrt(1 - 2 / pi);
threshold = IE_mean + 3.3 * std_half_normal;

% Calculate an Amplitude Threshold
[amp_peaks, ~] = findpeaks(a_t);
amp_peaks_rms = rms(amp_peaks);
amp_threshold = sqrt(2) * amp_peaks_rms;

% Mark Potential Burst Points
candidate_points = (IE_t > threshold) & (a_t > amp_threshold);

% Identify Continuous Burst Regions
bursts_labels = bwlabel(candidate_points);
num_bursts = max(bursts_labels);
burst_center_points = [];
burst_phases = [];
burst_start_end_points = [];
burst_signal = nan(size(X));  % Initialize burst signal with NaNs

% Compute Frequency Information
theta_t_unwrapped = unwrap(theta_t);
omega_t = [NaN; diff(theta_t_unwrapped) * Fs];
delta_omega_t = abs([NaN; diff(omega_t)]);
mean_delta_omega = mean(delta_omega_t,'omitmissing');
std_delta_omega = std(delta_omega_t,'omitmissing');

% Process Each Burst Region
upper_bound_freq = freq_band(2);
samples_per_cycle = Fs / upper_bound_freq;
samples_half_cycle = round(0.5 * samples_per_cycle);
min_burst_duration_samples = ceil(Num_cycle * samples_per_cycle);
min_burst_duration_samples2 = ceil(burst_duration_cycles * samples_per_cycle);
min_separation_samples = ceil(separating_cycle * samples_per_cycle);
burst_start_end_zero_crossings = [];
zero_crossings = find_zero_crossings(bp_signal);

for i = 1:num_bursts
    indices = find(bursts_labels == i);
    if length(indices) >= min_burst_duration_samples
        % Find the peak in this burst segment.
        [~, peak_idx_rel] = max(IE_t(indices));
        peak_idx = indices(1) + peak_idx_rel - 1;
        start_idx = indices(1);
        end_idx = indices(end);
        
        % Look backward from the peak to determine the burst's start.
        for idx = peak_idx:-1:indices(1)
            if IE_t(idx) < IE_mean || delta_omega_t(idx) > mean_delta_omega + 2 * std_delta_omega
                start_idx = idx + 1;
                break;
            end
        end
        
        % Look forward from the peak to determine the burst's end.
        for idx = peak_idx:indices(end)
            if IE_t(idx) < IE_mean || delta_omega_t(idx) > mean_delta_omega + 2 * std_delta_omega
                end_idx = idx - 1;
                break;
            end
        end
        
        % Only keep bursts that are long enough.
        if (end_idx - start_idx + 1) >= min_burst_duration_samples2
            center_point = round((start_idx + end_idx) / 2);
            burst_center_points = [burst_center_points; center_point];
            burst_phases = [burst_phases; theta_t(center_point)];
            burst_start_end_points = [burst_start_end_points; start_idx, end_idx];
            burst_signal(start_idx:end_idx) = bp_signal(start_idx:end_idx);
            
            % Adjust boundaries to align with zero crossings.
            extended_start_idx = max(start_idx - samples_half_cycle, 1);
            potential_starts = zero_crossings(zero_crossings <= extended_start_idx);
            if ~isempty(potential_starts)
                start_idx_zero = potential_starts(end);
            else
                start_idx_zero = 1;
            end
            extended_end_idx = min(end_idx + samples_half_cycle, length(X));
            potential_ends = zero_crossings(zero_crossings >= extended_end_idx);
            if ~isempty(potential_ends)
                end_idx_zero = potential_ends(1);
            else
                end_idx_zero = length(X);
            end
            burst_start_end_zero_crossings = [burst_start_end_zero_crossings; start_idx_zero, end_idx_zero];
        end
    end
end

% Filter Out Bursts That Are Too Close
burst_peak_amplitudes = zeros(size(burst_center_points));
for i = 1:length(burst_center_points)
    idx_range = burst_start_end_points(i, 1):burst_start_end_points(i, 2);
    burst_peak_amplitudes(i) = max(a_t(idx_range));
end
valid_bursts = true(size(burst_center_points));
for i = 1:length(burst_center_points) - 1
    if ~valid_bursts(i)
        continue;
    end
    for j = i + 1:length(burst_center_points)
        if ~valid_bursts(j)
            continue;
        end
        separation = burst_start_end_points(j, 1) - burst_start_end_points(i, 2);
        if separation < min_separation_samples
            % Keep the burst with the larger peak amplitude.
            if burst_peak_amplitudes(i) >= burst_peak_amplitudes(j)
                valid_bursts(j) = false;
            else
                valid_bursts(i) = false;
                break;
            end
        else
            break;
        end
    end
end
burst_center_points = burst_center_points(valid_bursts);
burst_phases = burst_phases(valid_bursts);
burst_start_end_points = burst_start_end_points(valid_bursts, :);
burst_start_end_zero_crossings = burst_start_end_zero_crossings(valid_bursts, :);

% Rebuild the Burst Signal for Valid Bursts
burst_signal = nan(size(X));
burst_signal_zc = nan(size(X));
for i = 1:length(burst_center_points)
    idx_range = burst_start_end_points(i, 1):burst_start_end_points(i, 2);
    burst_signal(idx_range) = bp_signal(idx_range);
    idx_range_zc = burst_start_end_zero_crossings(i, 1):burst_start_end_zero_crossings(i, 2);
    burst_signal_zc(idx_range_zc) = bp_signal(idx_range_zc);
end

% Compute the Burst Density
burst_events = zeros(size(X));
burst_events(burst_center_points) = 1;
window_length = round(0.5 * Fs);  % 500 ms window length
gauss_window = gausswin(window_length, 2.5);
gauss_window = gauss_window / max(gauss_window);
burst_density = conv(burst_events, gauss_window, 'same');

% Prepare the Output Structure
out = struct();
out.burst_center_points = burst_center_points;
out.burst_phases = burst_phases;
out.burst_start_end_points = burst_start_end_points;
out.burst_start_end_zero_crossings = burst_start_end_zero_crossings;
out.band_passed_signal = bp_signal;
out.burst_signal = burst_signal;
out.burst_signal_zc = burst_signal_zc;
out.burst_density = burst_density;
end

% Helper Function: Find Zero Crossings
function zero_crossings = find_zero_crossings(signal)
% Returns indices where the signal crosses zero.
sign_changes = diff(sign(signal));
zero_crossings = find(sign_changes ~= 0);
zero_indices = find(signal == 0);
zero_crossings = sort(unique([zero_crossings; zero_indices]));
end