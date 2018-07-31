addpath('./direc_6/')

ch0L = fread(fopen('ch0L.pcm', 'r'), inf, 'int32');
ch0R = fread(fopen('ch0R.pcm', 'r'), inf, 'int32');
ch1L = fread(fopen('ch1L.pcm', 'r'), inf, 'int32');
ch1R = fread(fopen('ch1R.pcm', 'r'), inf, 'int32');
ch2L = fread(fopen('ch2L.pcm', 'r'), inf, 'int32');
ch2R = fread(fopen('ch2R.pcm', 'r'), inf, 'int32');
ch3R = fread(fopen('ch3R.pcm', 'r'), inf, 'int32');
fs = 44100;

size_min = min([...
    size(ch0L,1) size(ch0R,1) ...
    size(ch1L,1) size(ch1R,1) ...
    size(ch2L,1) size(ch2R,1)]);

ch0L = ch0L(1:size_min);
ch0R = ch0R(1:size_min);
ch1L = ch1L(1:size_min);
ch1R = ch1R(1:size_min);
ch2L = ch2L(1:size_min);
ch2R = ch2R(1:size_min);
ch3R = ch3R(1:size_min);

ch_2 = [ch1L ch1R ch0L ch0R ch2L ch2R ch3R];

imagesc(cov(ch_2))
cov(ch_2)