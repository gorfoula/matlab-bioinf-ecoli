function [auc] = AUC(classes, predictions, positive)
% AUC Area under the ROC curve
%   [auc] = AUC(classes, predictions, positive) computes an approximation of
%   the integral of the ROC curve via the trapezoidal method.
%
%   The ROC curve is computed as described in 
%   "ROC Graphs: Notes and Practical Considerations for Data Mining
%   Researchers" by Tom Fawcett using a slightly modified version of the 
%   Algorithm 2 (p.11, section 5)
%   Alternatively the Matlab built-in function 'perfcurve' can be used to
%   compute the ROC Curve points using the same input arguments as in
%   RocPoiunts.
%   The main difference is that the currently used function is a lot faster
%   than the Matlab built-in function 'perfcurve'.
%
%   classes:      A vector of the true classes ( 1xN or Nx1 )
%   predictions:  A vector of the classifier's predictions  ( 1xN or Nx1 )
%   positive:     The positive class label  ( 1x1 )
%
%   Notes:
%   1) Checks about the sizes of the vectors are made
%   2) Checks about the correctness of the data are NOT made
%   3) All NaN values in the classes vector are ignored

% Computing the (X,Y) pairs of the ROC Curve
% Alternatively the built-in Matlab function perfcurve can be used by
% calling it with the following arguments:
% [X Y] = perfcurve(classes,predictions,positive);

% Error Checks
if(length(size(classes)) ~= 2 || (size(classes,1) ~= 1 && size(classes,2) ~= 1))
    error('classes not a vector (1xN or Nx1)');
end

if(length(size(predictions)) ~= 2 || (size(predictions,1) ~= 1 && size(predictions,2) ~= 1))
    error('predictions not a vector (1xN or Nx1)');
end

if(size(positive,1) ~= 1 || size(positive,2) ~= 1)
    error('positive not a 1x1 vector');
end

if(size(classes,2) ~= 1)
    classes = classes';
end

if(size(predictions,2) ~= 1)
    predictions = predictions';
end

if(size(classes,1) ~= size(predictions,1))
    error('Number of classes and predictions do not match');
end

[X Y] = RocPoints(classes,predictions, positive);
% Computing the area under the ROC curve (AUC) 
% by integrating the curve using the Trapezoidal rule
figure(10);plot(X,Y);pause; hold off;
auc = trapz(X,Y);

end

function [ X, Y ] = RocPoints( classes, predictions, positive )
maxvalue = max(predictions);

if(maxvalue ~= positive)
    maxindexes = predictions == maxvalue;
    posindexes = predictions == positive;
    predictions(maxindexes) = positive;
    predictions(posindexes) = maxvalue;
    
    maxindexes = classes == maxvalue;
    posindexes = classes == positive;
    classes(maxindexes) = positive;
    classes(posindexes) = maxvalue;
    
    positive = maxvalue;
end

L = sortrows([predictions classes]);
P = length(find(classes == positive));
N = length(find(classes ~= positive)) - length(find(isnan(classes)));

FP = 0;
TP = 0;
X = zeros(2*size(L,1),1);
Y = zeros(2*size(L,1),1);
Lprev = -Inf;

curpoints = 2*size(L,1);

for i = 1:size(L,1)
    if(isnan(L(i,1)))
        continue;
    end
    
    if( L(i,1) ~= Lprev )
        X(curpoints) = FP/N;
        Y(curpoints) = TP/P;
        curpoints = curpoints-1;
        Lprev = L(i,1);
    end
    
    if( L(i,2) == positive )
        TP = TP + 1;
    else
        FP = FP + 1;
    end
    
    X(curpoints) = FP/N;
    Y(curpoints) = TP/P;
    
    curpoints = curpoints-1;
end

X = 1-X(curpoints+1:end);
Y = 1-Y(curpoints+1:end);

end
