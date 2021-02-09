function out = SB_MotifTwo_reduced(y,binarizeHow)
% SB_MotifTwo   Local motifs in a binary symbolization of the time series
%
% Coarse-graining is performed by a given binarization method.
%
%---INPUTS:
% y, the input time series
% binarizeHow, the binary transformation method:
%       (i) 'diff': incremental time-series increases are encoded as 1, and
%                   decreases as 0,
%       (ii) 'mean': time-series values above its mean are given 1, and those
%                    below the mean are 0,
%       (iii) 'median': time-series values above the median are given 1, and
%       those below the median 0.
%
%---OUTPUTS:
% Probabilities of words in the binary alphabet of lengths 1, 2, 3, and 4, and
% their entropies.
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

if nargin < 2 || isempty(binarizeHow)
    % Use changes in the time series as the basis for the transformation
    binarizeHow = 'diff';
end

% Generate a binarized version of the input time series:
yBin = BF_Binarize(y,binarizeHow);

% Define the length of the new, symbolized sequence: N
N = length(yBin);

if N < 5
    warning('Time series too short');
    out = NaN; return
end

% Compute proportion of 1, make sure ranges are valid for looking at the 
% next one
r1 = (yBin == 1);
r1 = r1(1:end-1);

% Record up-up
r11 = r1 & yBin(2:end) == 1;
out.uu = mean(r11);

end
