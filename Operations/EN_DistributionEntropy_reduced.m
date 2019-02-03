function out = EN_DistributionEntropy_reduced(y,histOrKS,numBins,~)
% EN_DistributionEntropy    Distributional entropy.
%
% Estimates of entropy from the distribution of a data vector. The
% distribution is estimated either using a histogram with numBins bins, or as a
% kernel-smoothed distribution, using the ksdensity function from Matlab's
% Statistics Toolbox with width parameter, w (specified as the iunput numBins).
%
% An optional additional parameter can be used to remove a proportion of the
% most extreme positive and negative deviations from the mean as an initial
% pre-processing.
%
%---INPUTS:
%
% y, the input time series
%
% histOrKS: 'hist' for histogram, or 'ks' for ksdensity
%
% numBins: (*) (for 'hist'): an integer, uses a histogram with that many bins
%          (*) (for 'ks'): a positive real number, for the width parameter for
%                       ksdensity (can also be empty for default width
%                                       parameter, optimum for Gaussian)
%
% olremp [opt]: the proportion of outliers at both extremes to remove
%               (e.g., if olremp = 0.01; keeps only the middle 98% of data; 0
%               keeps all data. This parameter ought to be less than 0.5, which
%               keeps none of the data).
%               If olremp is specified, returns the difference in entropy from
%               removing the outliers.
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
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(histOrKS)
    histOrKS = 'hist'; % use histogram by default
end
if nargin < 3 % (can be empty for default width for ksdensity)
    numBins = 'fd'; % Changed default, see help for histcounts
end

% ------------------------------------------------------------------------------
% Form the histogram
% ------------------------------------------------------------------------------
switch histOrKS
case 'hist' 
    
    % Use histogram to calculate pdf
    [px,binEdges] = histcounts(y,'BinMethod',numBins,'Normalization',...
        'probability');
    
    % Compute bin widths:
    binWidths = diff(binEdges);

    otherwise
    % error; must specify 'hist'
    error('Unknown distribution method -- specify ''ks'' or ''hist''')
end


% ------------------------------------------------------------------------------
% Compute the entropy sum and return it as output
% ------------------------------------------------------------------------------
out = -sum(px(px>0).*log(px(px>0)./binWidths(px>0)));

end
