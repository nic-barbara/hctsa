function out = SP_Summaries_reduced(y,psdMeth,~,~,dologabs)
% SP_Summaries  Statistics of the power spectrum of a time series
%
% The estimation can be done using a periodogram, using the periodogram code in
% Matlab's Signal Processing Toolbox, or a fast fourier transform, implemented
% using Matlab's fft code.
%
%---INPUTS:
% y, the input time series
%
% psdMeth, the method of obtaining the spectrum from the signal:
%               (i) 'periodogram': periodogram
%               (ii) 'fft': fast fourier transform
%               (iii) 'welch': Welch's method
%
% windowType, the window to use:
%               (i) 'boxcar'
%               (ii) 'rect'
%               (iii) 'bartlett'
%               (iv) 'hann'
%               (v) 'hamming'
%               (vi) 'none'
%
% nf, the number of frequency components to include, if
%           empty (default), it's approx length(y)
%
% dologabs, if 1, takes log amplitude of the signal before
%           transforming to the frequency domain.
%
% doPower, analyzes the power spectrum rather than amplitudes of a Fourier
%          transform
%
%---OUTPUTS:
% Statistics summarizing various properties of the spectrum,
% including its maximum, minimum, spread, correlation, centroid, area in certain
% (normalized) frequency bands, moments of the spectrum, Shannon spectral
% entropy, a spectral flatness measure, power-law fits, and the number of
% crossings of the spectrum at various amplitude thresholds.
%
% This version of this fcn has been reduced and edited by Nicholas Barbara.
% Email: nbar5346@uni.sydney.edu.au

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check that a Curve-Fitting Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox')

% ------------------------------------------------------------------------------
% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y,2) > size(y,1)
    y = y'; % Time series must be a column vector
end
if nargin < 2 || isempty(psdMeth)
    psdMeth = 'fft'; % fft by default
end
if nargin < 5 || isempty(dologabs)
    dologabs = false;
end

if dologabs % a boolean
    % Analyze the spectrum of logarithmic absolute deviations
    y = log(abs(y));
end

Ny = length(y); % time-series length


% ------------------------------------------------------------------------------
% Compute the Fourier Transform
% ------------------------------------------------------------------------------
switch psdMeth
    case 'fft'
        % Fast Fourier Transform
        NFFT = 2^nextpow2(Ny);
        S = fft(y,NFFT); % Fourier Transform
        S = 2*abs(S(1:NFFT/2+1)).^2/Ny; % single-sided power spectral density
        S = S/(2*pi); % convert to angular frequency space
    otherwise
        error('Unknown spectral estimation method ''%s''',psdMeth);
end

if ~any(isfinite(S)) % no finite values in the power spectrum
    % This time series must be really weird -- return NaN (unsuitable operation)...
    warning('NaN in power spectrum? A weird time series.');
    out = NaN; return
end

% Ensure S is a row vector:
if size(S,1) > size(S,2)
    S = S';
end


%-------------------------------------------------------------------------------
% Peaks:
%-------------------------------------------------------------------------------

% Characterize all peaks using findpeaks function:
% Minimum angular separation of 0.02...?
minDist_w = 0.02;
ptsPerw = length(S)/pi;
minPkDist = ceil(minDist_w*ptsPerw);
[pkHeight,~,pkWidth,pkProm] = findpeaks(S,'SortStr','descend','minPeakDistance',minPkDist);
pkWidth = pkWidth/ptsPerw;

% Power in peaks with prominence of at least 2
out.peakPower_prom2 = sum(pkHeight(pkProm > 2).*pkWidth(pkProm > 2));

end
