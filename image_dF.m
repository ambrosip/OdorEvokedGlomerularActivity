baselineStack = single(read_file(imgToAnalyzeFileDir));

lowQuantile = double(quantile(baselineStack(:), 0.005));
highQuantile = double(quantile(baselineStack(:), 0.995));
quantiles = [lowQuantile, highQuantile];

baselineMeanImage = mat2gray(mean(baselineStack, ndims(baselineStack)), quantiles);
clear("baselineStack");



odorStack = single(read_file(imgToAnalyzeFileDir));

lowQuantile = double(quantile(odorStack(:), 0.005));
highQuantile = double(quantile(odorStack(:), 0.995));
quantiles = [lowQuantile, highQuantile];

odorMeanImage = mat2gray(mean(odorStack, ndims(odorStack)), quantiles);
clear("odorStack");

%%

subtractedImage = imsubtract(odorMeanImage,baselineMeanImage);
figure; imshow(subtractedImage)
dividedImage = imdivide(subtractedImage,baselineMeanImage);
figure; imshow(dividedImage)