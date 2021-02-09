function out = CO_AddNoise_reduced(y,tau,~,extraParam,randomSeed)
% CO_AddNoise  Changes in the automutual information with the addition of noise
%
% Adds Gaussian-distributed noise to the time series with increasing standard
% deviation, eta, across the range eta = 0, 0.1, ..., 2, and measures the
% mutual information at each point
% Can be measured using histograms with extraParam bins (implemented using
% CO_HistogramAMI), or using the Information Dynamics Toolkit.
%
% The output is a set of statistics on the resulting set of automutual
% information estimates, including a fit to an exponential decay, since the
% automutual information decreases with the added white noise.
%
% Can calculate these statistics for time delays 'tau', and for a number 'extraParam'
% bins.
%
% This algorithm is quite different, but was based on the idea of 'noise
% titration' presented in: "Titration of chaos with added noise", Chi-Sang Poon
% and Mauricio Barahona P. Natl. Acad. Sci. USA, 98(13) 7107 (2001)
%
% This version of this fcn has been reduced and edited by Nicholas Barbara.
% Email: nbar5346@uni.sydney.edu.au
%
%---INPUTS:
%
% y, the input time series
%
% tau, the time delay for computing AMI
%
% amiMethod, the method for computing AMI:
%      * one of 'std1','std2','quantiles','even' for histogram-based estimation,
%      * one of 'gaussian','kernel','kraskov1','kraskov2' for estimation using JIDT
%
% extraParam, e.g., the number of bins input to CO_HistogramAMI, or parameter
%             for IN_AutoMutualInfo
%
% randomSeed: settings for resetting the random seed for reproducible results
%               (using BF_ResetSeed)

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
% Preliminary checks
% ------------------------------------------------------------------------------
% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = 0; % plot outputs to figure

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
% Expecting a z-scored input time series:
BF_iszscored(y);

if nargin < 2
    tau = []; % set default in CO_HistogramAMI
end
% Set tau to minimum of autocorrelation function
if ~isempty(tau) && ischar(tau) && ismember(tau,{'ac','tau'})
    tau = CO_FirstZero(y,'ac');
end
if nargin < 3
    amiMethod = 'even'; % using evenly spaced bins in CO_HistogramAMI
end
if nargin < 4
    extraParam = []; % number of bins for CO_HistogramAMI
end
if nargin < 5
    randomSeed = [];
end

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------
noiseRange = linspace(0,3,50); % compare properties across this noise range
BF_ResetSeed(randomSeed); % reset the random seed if specified
numRepeats = length(noiseRange);
amis = zeros(numRepeats,1);
noise = randn(size(y)); % uncorrelated additive noise

% ------------------------------------------------------------------------------
% Compute the automutual information across a range of noise levels
% ------------------------------------------------------------------------------
% The *same* noise vector, noise, is added to the signal, with increasing
% standard deviation (one could imagine repeating the calculation with different
% random seeds)...
amiMethod = 'std1';

% Histogram-based methods using my naive implementation in CO_Histogram.m
for i = 1:numRepeats
    amis(i) = CO_HistogramAMI(y+noiseRange(i)*noise,tau,amiMethod,extraParam);
    if isnan(amis(i))
        error('Error computing AMI: Time series too short (?)');
    end
end


% ------------------------------------------------------------------------------
% Fit linear function to output
% ------------------------------------------------------------------------------
p = polyfit(noiseRange',amis,1);
out.fitlinb = p(2); % intercept

end
