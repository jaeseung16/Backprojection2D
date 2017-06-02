load('rawData.mat'); % example data

twoRR = abs(rawData.spec); % absolute spectra

size(twoRR)

%%
test = BackProjection2D( twoRR, (-100:100) + 512 ); % initialize

%%
test.showImage(); % show the reconstructed image

%%
frames = test.makeMovie(); % show the progress of the image reconstration