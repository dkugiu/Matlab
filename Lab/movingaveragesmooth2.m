function sV = movingaveragesmooth2(xV,maorder) 
% sV = movingaveragesmooth2(xV,maorder)
% MOVINGAVERAGESMOOTH2 makes a smoothing of the given time series using a 
% moving average filter of a given order. It makes use of the 'filtfilt'
% Matlab function and produces a time series of the same length as the
% given time series with the first and last maorder/2 values being somehow
% estimated by the filter.
% INPUTS 
% - xV      : vector of length 'n' of the time series
% - maorder : the maorder of the moving average filter
% OUTPUTS
% - sV      : vector of length 'n' of the smoothed time series

n = length(xV);
xV = xV(:);
if maorder > 1
    b = ones(1,maorder)/maorder;
    sV = filtfilt(b,1,xV);
else
    sV = NaN*ones(n,1);
end
