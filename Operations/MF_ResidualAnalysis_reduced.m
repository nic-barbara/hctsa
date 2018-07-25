function out = MF_ResidualAnalysis_reduced(e)
% MF_ResidualAnalysis   Analysis of residuals from a model fit.
%
% Given an input residual time series residuals, e, this function returns a
% structure with fields corresponding to statistical tests on the residuals.
% These are motivated by a general expectation of model residuals to be
% uncorrelated.
%
% This version of this fcn has been reduced and edited by Nicholas Barbara.
% Email: nbar5346@uni.sydney.edu.au
%
%---INPUT:
% e, should be raw residuals as prediction minus data (e = yp - y) as a column
%       vector.

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
%% Preliminaries
% ------------------------------------------------------------------------------

% Check that a System Identification Toolbox license is available (for spa):
BF_CheckToolbox('identification_toolbox')

if size(e,2) > size(e,1)
    e = e'; % make sure residuals are a column vector
end
if all(e>0)
    warning('Very weird that all model residuals are positive...')
elseif all(e<0)
    warning('Very weird that all model residuals are negative...')
end
N = length(e);


% ------------------------------------------------------------------------------
%% Analyze autocorrelation in residuals
% ------------------------------------------------------------------------------
% See if there are any linear correlations in residuals.
% Also see if any of these are abnormally large (i.e., may be remnant
% autocorrelation at some level, or may be a characteristic shape in this
% function...)
% Will output values scaled by sqrt(length), as is normal (within a 
% constant).
maxLag = 25;

acs = CO_AutoCorr(e,1:maxLag,'Fourier'); % autocorrelations
sqrtN = sqrt(N);

% Stdev normalized distance from zero
out.acsnd0 = std(abs(acs))*sqrtN;

end
