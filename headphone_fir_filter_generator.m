%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is a Headphone FIR-filter generator - a GNU Octave script, intended to 
%   generate FIR filter to use with headphones.
%   Copyright (C) 2021 Alexey Kovalev
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define global parameters...

directory = './';

fname_direct = 'IRC_1006_R_R0195_T030_P000';
fname_opposite = 'IRC_1006_R_R0195_T330_P000';

channel_idx = 1; % Channel index, used from input files: 1 - left, 2 - right.
fs = 44100;
max_impulse_len = 4096; % Suppose, there is no silence at the start of the impulses, which we read from disk


function impulse = read_impulse(fname, channel_idx, max_impulse_len)
	[impulse, fs] = audioread(fname);
	len = min(rows(impulse), max_impulse_len);
	impulse = impulse(1 : len, channel_idx);
endfunction

function view_results(impulse_direct, impulse_opposite, impulse_tf)
	impulse_test = conv(impulse_direct, impulse_tf);
	impulse_test = impulse_test(1 : rows(impulse_opposite)); % match the lengths
	x = linspace(1, rows(impulse_opposite), rows(impulse_opposite));
	plot(x, 6 + 10 * impulse_test, x, 6 + 10 * impulse_opposite, 
        x, 4 + abs(fft(impulse_test)), x, 4 + abs(fft(impulse_opposite)), 
        x, 2 + real(fft(impulse_test)), x, 2 + real(fft(impulse_opposite)), 
        x, imag(fft(impulse_test)), x, imag(fft(impulse_opposite)));
endfunction

function filteredImpulse = lowPassFilter(impulse, fs)
    pkg load signal;
	% Generating low pass FIR filter using Parks-McClellan algorithm.
	% One may also use fir1() function for window-based filter generation, or any other suitable method.
	lpf = remez(2 * round(fs / 600), [0 20000 / (fs / 2) 21000 / (fs / 2) 1], [1 1 0 0]); % For Matlab environment change "remez" to "firpm". 

    filteredImpulse=conv(impulse, lpf, 'same');
endfunction

% Working...

% Load input data...
impulse_direct = read_impulse(strcat(directory, fname_direct, '.wav'), channel_idx, max_impulse_len);
impulse_opposite = read_impulse(strcat(directory, fname_opposite, '.wav'), channel_idx, max_impulse_len);

if size(impulse_direct) ~= size(impulse_opposite)
    error('Size mismatch');
end

% Deconvolution by least mean square method, Levinson-Durbin matrix inversion
impulse_tf = lmsDeconv(impulse_direct, impulse_opposite, size(impulse_direct));

% Optional: cut ultrasonic noise from the impulse.
impulse_tf = lowPassFilter(impulse_tf, fs);

% Windowing
% Approximate delay, introduced by the filter at most frequences - it will be center of the window
avgDelay = round(median(grpdelay(impulse_tf)));
wnd = window(@blackman, 2 * (size(impulse_tf)(1) - avgDelay));
% Cut the window in size
wnd = flipud(wnd(1:size(impulse_tf)));
% Apply window
impulse_tf = impulse_tf .* wnd;

% Make "left" stereo impulse for the convolver
impulse_tf_stereo_L(:, 2) = impulse_tf;
impulse_tf_stereo_L(1, 1) = 1;

% Make "right" stereo impulse for the convolver (the same, but left-right channels reversed)
impulse_tf_stereo_R(:, 1) = impulse_tf;
impulse_tf_stereo_R(1, 2) = 1;

% Write results to disk...
audiowrite(strcat(directory, fname_direct, '_TF_mono_raw.wav'), impulse_tf, fs);
audiowrite(strcat(directory, fname_direct, '_TF_stereo_L.wav'), impulse_tf_stereo_L, fs);
audiowrite(strcat(directory, fname_direct, '_TF_stereo_R.wav'), impulse_tf_stereo_R, fs);

% Now check visually if impulse_tf is similar to impulse_opposite
view_results(impulse_direct, impulse_opposite, impulse_tf);
