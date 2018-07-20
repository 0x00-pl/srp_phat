clear;
wlen = 4096;
% [ref, fs] = audioread('./ch1L.wav');
% ch1L = delayseq(ref,0.06928/340,fs);
% ch0L = delayseq(ref,0/340,fs);
% ch2L = delayseq(ref,0.06928*2/340, fs);
% ch1R = delayseq(ref, 0/340, fs);
% ch2R = delayseq(ref, 0.06928*2/340, fs);
% ch0R = delayseq(ref, 0.06928/340, fs);

 [ch1L, fs] = audioread('./ch1L.wav');
 [ch1R, ~] = audioread('./ch1R.wav');
 [ch0L, ~] = audioread('./ch0L.wav');
 [ch0R, ~] = audioread('./ch0R.wav');
 [ch2L, ~] = audioread('./ch2L.wav');
 [ch2R, ~] = audioread('./ch2R.wav');
%ch_all = [ch1L ch1R ch0L ch0R ch2L ch2R];
%ch_all = ch_all(180000:360000,:);
ch_2 = [ch1L ch1R ch0L ch0R ch2L ch2R];
%ch_2 = ch_2(40000:120000,:);
mic_loc_2 = 0.8*[0.1000    0.0000         0
    0.0500    0.0866         0
    -0.05    0.0866         0
 -0.1000    0.0000         0
    -0.05    -0.0866         0
    0.0500    -0.0866         0
];
mic_a_2 = [0 60 120 180 240 300];
num_frame = floor(length(ch_2)/wlen);
ch_slice = zeros(wlen,2);


mic_loc = 0.8*[0.1000    0.0000         0
    0.0500    0.0866         0
    -0.05    0.0866         0
    -0.1000    0.0000         0
    -0.05    -0.0866         0
    0.0500    -0.0866         0];
mic_a = [0 60 120 180 240 300];
num_doa = 60;
r = 0.08;
max_id = zeros(1,num_frame-1);
tau = zeros(1,num_frame-1);
tau2 = zeros(1,num_frame-1);
tau3 = zeros(1,num_frame-1);
srp = zeros(num_doa,num_frame-1);
tic
for i=1:num_frame-1
   ch_slice = ch_2(i*wlen+1:(i+1)*wlen,:); 
   tau(i) = gccphat(ch_slice(:,1),ch_slice(:,2));
   tau2(i) = gccphat(ch_slice(:,2),ch_slice(:,3));
%    tau3(i) = gccphat(ch_slice(:,1),ch_slice(:,3));
%    tau4(i) = gccphat(ch_slice(:,1),ch_slice(:,4));
%    tau5(i) = gccphat(ch_slice(:,2),ch_slice(:,4));
%    tau6(i) = gccphat(ch_slice(:,3),ch_slice(:,4));
   %Sx = stft_multi(ch_slice',wlen);
   Sx = fft(ch_slice, wlen);
   [srp(:,i), max_id(i)] = srp_phat(Sx,mic_loc_2,mic_a_2, num_doa, r,fs);
end
toc
subplot(411)
plot(ch_2(:,2))
axis([1 length(ch_2)-1 -0.02 0.02] );

subplot(412)

plot([max_id' 20*max(srp)'])
title("direction expectation")
axis([1 num_frame-1 1 num_doa]);
subplot(413)

imagesc(srp);
title("srp direction expectation")
axis([1 num_frame-1 1 num_doa]);
subplot(414)

plot([tau' tau3']);
title("delay expectation")
% hold on;
% 
% plot(tau2);
% plot(tau3);
% plot(tau4);
% plot(tau5);
% plot(tau6);
% hold off;
axis([1 num_frame-1 -20 20]);

