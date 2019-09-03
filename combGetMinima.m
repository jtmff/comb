function vectorMinima = combGetMinima(trace, bcl, varargin)
%% A function searching for minima in traces from regularly activated cardiac preparations (or similar data).
% trace - a vector containing signal (e.g., intensity of a pixel in optical
% mapping, or an electrophysiological recording). Assumed to be at 1 KHz
% (this is not a hard requirement - the algorithm can run without any
% changes provided the sampling rate is constant).
%
% IN:
%
% bcl - the basic cycle length (number of ms or, in general case, samples,
% between two activations)
%
% varargin{1} - refinementWidth parameter (in ms/samples - the radius of
% local search around comb teeth)
%
% varargin{2} = a 2-by-1 vector determining boundaries from where to where
% are the minima reported (this may be used e.g., to ignore artefactual
% pseudominima at the start of the recording, which may arise from imaging
% artefacts, etc.).
%
% OUT:
% vectorMinima - the vector of local minima in the signal, which are
% approximately bcl ms (or samples) apart.

%% Default parameter initialization and processing of extra inputs.
refinementWidth = 10;
minStart = 1;
maxEnd = length(trace);

if (numel(varargin) >= 1) % 1st extra argument is the width of local search
    refinementWidth = varargin{1};
end

if (numel(varargin) >= 2)
    limits = varargin{2};
    try
    minStart = limits(1);
    maxEnd = limits(2);
    catch
        error('The last parameter must be a vector of length 2 specifying from where and until where are minima extracted');
    end
end

%% Comb positioning
nFrames = length(trace);

candidateMinimaMean = inf*ones(bcl-1, 1);
for iStart = minStart:(minStart + bcl-1)
    candidateMinima = iStart:bcl:nFrames;
    candidateMinimaMean(iStart-minStart+1) = mean(trace(candidateMinima));
end

% The start leading to minimal average value is used as the
% indicator of minima.
[~, bestMinima] = min(candidateMinimaMean);
vectorMinima = minStart - 1 + bestMinima:bcl:nFrames;

%% Comb refinement 
% The vector of teeth (~near-minima) is refined by looking refinementWidth to the
% left and right, taking the actual minimum there.
vectorMinimaRefined = zeros(size(vectorMinima));

for iMinimum = 1:length(vectorMinima)
    leftBound = max(vectorMinima(iMinimum) - refinementWidth, 1);
    rightBound = min(vectorMinima(iMinimum) + refinementWidth, length(trace));

    [~, actualMinimum] = min(trace(leftBound:rightBound));
    vectorMinimaRefined(iMinimum) = leftBound + actualMinimum - 1;
end
vectorMinima = vectorMinimaRefined;

vectorMinima(vectorMinima > maxEnd) = [];
end