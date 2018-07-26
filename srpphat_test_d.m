clear;
figure(3)
wlen = 512;
fluence = 1;
% [ref, fs] = audioread('./ch1L.wav');
% 
% ch1L = delayseq(ref,cosd(45)*0.06928/340,fs);
% ch0L = delayseq(ref,0/340,fs);
% ch2L = delayseq(ref,cosd(45)*0.06928*2/340, fs);
% ch1R = delayseq(ref, 0/340, fs);
% ch2R = delayseq(ref, cosd(45)*0.06928*2/340, fs);
% ch0R = delayseq(ref, cosd(45)*0.06928/340, fs);
addpath('./direc_4/')

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

%  [ch1L, fs] = audioread('./ch1L.wav');
%  [ch1R, ~] = audioread('./ch1R.wav');
%  [ch0L, ~] = audioread('./ch0L.wav');
%  [ch0R, ~] = audioread('./ch0R.wav');
%  [ch2L, ~] = audioread('./ch2L.wav');
%  [ch2R, ~] = audioread('./ch2R.wav');
%ch_all = [ch1L ch1R ch0L ch0R ch2L ch2R];
%ch_all = ch_all(180000:360000,:);
ch_2 = [ch1L ch1R ch0L ch0R ch2L ch2R ch3R];
ch_2=filter([1,-0.97],1,ch_2);
%ch_2 = ch_2(40000:120000,:);
mic_loc_2 = 0.8*[0.1000    0.0000         0
    0.0500    0.0866         0
    -0.05    0.0866         0
    -0.1000    0.0000         0
    -0.05    -0.0866         0
    0.0500    -0.0866         0
    0 0 0
];
%mic_a_2 = [0 60 120 180 240 300];
num_frame = floor(length(ch_2)/wlen);
ch_slice = zeros(wlen,2);


% mic_loc = 0.81*[0.1000    0.0000         0
%     0.0500    0.0866         0
%     -0.05    0.0866         0
%     -0.1000    0.0000         0
%     -0.05    -0.0866         0
%     0.0500    -0.0866         0];
mic_a = [0 60 120 180 240 300];
num_doa = 60;   % 0~360  
num_doa_high = 3;  %  0~60
r = 0.08;
max_id = zeros(1,num_frame-1);
max_m = zeros(1,num_frame-1);
tau = zeros(1,num_frame-1);
tau2 = zeros(1,num_frame-1);
tau3 = zeros(1,num_frame-1);
srp = zeros(num_doa,num_frame-1);
srp_m = zeros(num_doa,num_doa_high, num_frame-1);
srp_ch = zeros(num_frame-1, num_doa, num_doa_high, 21);
tic
for i=1:num_frame-1
   ch_slice = ch_2(i*wlen+1:(i+1)*wlen,:); 
   tau(i) = gccphat(ch_slice(:,1),ch_slice(:,2));
%   tau2(i) = gccphat(ch_slice(:,2),ch_slice(:,3));
%    tau3(i) = gccphat(ch_slice(:,1),ch_slice(:,3));
%    tau4(i) = gccphat(ch_slice(:,1),ch_slice(:,4));
%    tau5(i) = gccphat(ch_slice(:,2),ch_slice(:,4));
%    tau6(i) = gccphat(ch_slice(:,3),ch_slice(:,4));
   %Sx = stft_multi(ch_slice',wlen);
   Sx = fft(ch_slice, wlen);
   %Sx(1:25,:) = Sx(1:25,:).*0.5;
   [srp(:,i), max_id(i), max_m(i), srp_ch(i,:,:,:)] = srp_phat_d(Sx,mic_loc_2, num_doa,num_doa_high,fs);
end
toc
subplot(511)
plot(ch_2(:,2))
axis([1 length(ch_2)-1 -0.02 0.02] );

subplot(512)

plot(max_id)
hold on;
plot(max_m)
hold off;
title("direction expectation")
axis([1 num_frame-1 1 num_doa]);
subplot(513)
temp = zeros(num_doa,num_frame-16);
max_t = zeros(floor(num_frame-16));
temp_vad = sum(srp(:,1:fluence),2);
[thread_v,~] = max(temp_vad(:,1));
for bin = 1:num_frame-fluence
   temp(:,bin) = sum(srp(:,bin:bin+fluence-1),2);
   [e_v,~] = max(temp(:,bin));
   if e_v > thread_v*1.5
        [~,max_t(bin)] = max(temp(:,bin));
   end
end
% 
imagesc(temp);
hold on;
plot(max_t ,'w');
hold off;
axis([1 num_frame-fluence 1 num_doa]);
% imagesc(srp);
% title("srp direction expectation")
% axis([1 num_frame-1 1 num_doa]);
subplot(514)

max_n = zeros(num_frame-fluence);
for bin =1:num_frame-fluence
    [max_n(bin),~] = max(temp(:,bin));
end
plot(max_n ,'r');
axis([1 num_frame-fluence 0 fluence*2]);
subplot(515)

plot(tau);
title("delay expectation")
hold on;

% plot(tau2);
% plot(tau3);
% plot(tau4);
% plot(tau5);
% plot(tau6);
hold off;
axis([1 num_frame-1 -20 20]);

srp_ch_sum = squeeze(sum(sum(srp_ch),3));

srp_dir = zeros(60,6);
srp_dir(:,1) = sum(srp_ch_sum(:, [7, 18, 6, 3, 19]), 2);
srp_dir(:,2) = sum(srp_ch_sum(:, [2, 17]), 2);
srp_dir(:,3) = sum(srp_ch_sum(:, [1, 15, 21, 14, 16]), 2);
srp_dir(:,4) = sum(srp_ch_sum(:, [13, 10]), 2);
srp_dir(:,5) = sum(srp_ch_sum(:, [12, 11, 20, 9, 5]), 2);
srp_dir(:,6) = sum(srp_ch_sum(:, [8, 4]), 2);
