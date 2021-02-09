function out = CO_Embed2_Shapes_reduced(y,tau,shape,r)
% CO_Embed2_Shapes Shape-based statistics in a 2-d embedding space
%
% Takes a shape and places it on each point in the two-dimensional time-delay
% embedding space sequentially. This function counts the points inside this shape
% as a function of time, and returns statistics on this extracted time series.
%
% This version of this fcn has been reduced and edited by Nicholas Barbara.
% Email: nbar5346@uni.sydney.edu.au
%
%---INPUTS:
% y, the input time-series as a (z-scored) column vector
% tau, the time-delay
% shape, has to be 'circle' for now...
% r, the radius of the circle
%
%---OUTPUTS:
% The constructed time series of the number of nearby points, and
% include the autocorrelation, maximum, median, mode, a Poisson fit to the
% distribution, histogram entropy, and stationarity over fifths of the time
% series.

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
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
	fprintf(1,['Setting tau as the first zero crossing ' ...
    			'of the autocorrelation function.\n']);
    tau = 'tau';
end
if nargin < 3 || isempty(shape)
	shape = 'circle'; % use a circle by default
	fprintf(1,'Using a circle.\n');
end
if nargin < 4 || isempty(r)
	r = 1; % default radius of 1
end

% Can set time lag equal to first zero crossing of the autocorrelation function with the 'tau' input
if strcmp(tau,'tau')
    tau = CO_FirstZero(y,'ac');
    % Cannot set the time delay greater than 10% the length of the time series
    if tau > length(y)/10
        tau = floor(length(y)/10);
    end
end

% Ensure y is a column vector:
if size(y,2) > size(y,1)
	y = y';
end

% ------------------------------------------------------------------------------
%% Create the recurrence space, populated by points m
% ------------------------------------------------------------------------------
m = [y(1:end-tau), y(1+tau:end)];
N = length(m);


%% Start the analysis

counts = zeros(N,1); % stores the counts for points within circle

switch shape
	case 'circle'
		% Puts a circle around each point in the embedding space in turn
		% counts how many points are inside this shape, looks at the time series thus formed

		for i = 1:N % across all points in the time series

			m_c = m - ones(N,1)*m(i,:); % points wrt current point i
			m_c_d = sum(m_c.^2,2); % Euclidean distances from point i

		    counts(i) = sum(m_c_d <= r^2); % number of points enclosed in a circle of radius r
		end

    otherwise
        error('Unknown shape ''%s''', shape)
end
counts = counts - 1; % ignore self-counts

% ------------------------------------------------------------------------------
% No meaningful output if never got a count with any other point!
% ------------------------------------------------------------------------------
% (radius, r, is probably too small)
if all(counts == 0)
    fprintf(1,'No counts detected!\n');
    out = NaN; return
end

% ------------------------------------------------------------------------------
% Return basic statistics on the counts
% ------------------------------------------------------------------------------
out.ac1 = CO_AutoCorr(counts,1,'Fourier');

end
