function w = getWeightFactors(svmStruct)

nSize = size(svmStruct.SupportVectors, 2);
if (numel(svmStruct.ScaleData) > 0)
    shift = svmStruct.ScaleData.shift;
    scale = svmStruct.ScaleData.scaleFactor;
    w = svmStruct.Alpha' * svmStruct.SupportVectors .* scale;
    w(nSize+1) = svmStruct.Bias + w(1:nSize) * shift' ;
else
    w = svmStruct.Alpha' * svmStruct.SupportVectors;
    w(nSize+1) = svmStruct.Bias ;
end
w = -w'; 