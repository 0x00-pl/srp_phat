function [srp, max_id, max_m]=srp_phat_d(Sx, mic_loc,mic_a, num_doa,num_doa_high, r, fs)
    %[nbin, num_frames, num_chs] = size(Sx); 
    [nbin, num_chs] = size(Sx); 
%     delay_mic = [0.04 0.0283 0 -0.0283 -0.04 -0.0283 0 0.0283;
%                  0.04 0.0566 0.04 0 -0.04 -0.0566 -0.04 0;
%                  0 0.0283 0.04 0.0283 0 -0.0283 -0.04 -0.0283;
%                  0 0.0283 0.04 0.0283 0 -0.0283 -0.04  -0.0283;
%                  -0.04 0 0.04 0.0566 0.04 0 -0.04 -0.0566;
%                  -0.04 -0.0283 0 0.0283 0.04 0.0283 0 -0.0283];

    delay_mic = compute_delay_mic_loc(mic_loc,mic_a, num_doa,num_doa_high, r);
    shiftTau = delay_mic/340;
    maxshift = floor(nbin/2);
    shiftdelay = round(-fs * shiftTau+ maxshift + 1);
    % GCC-PHAT
    N = num_chs * (num_chs - 1)/2;

    Z= complex(zeros(nbin ,N));
    p= 0;
    for m1=1:num_chs-1
        for m2=m1+1:num_chs
            p= p+1;
            Z(:,p)= Sx(:, m1).*conj(Sx(:,m2));
            Z(:,p) = Z(:,p)./(abs(Z(:,p)) +eps);           
        end
    end
    %R= zeros(nbin ,N,num_frames);
    R= zeros(nbin ,N);
    for p=1:N
        %for q = 1:num_frames
           % R(:,p,q)= fftshift(real(ifft(Z(:,p,q))));
        R(:,p)= fftshift(real(ifft(Z(:,p))));
        %end
    end
    % SRP search
    srp = zeros(num_doa,1);
    srp_m = zeros(num_doa,num_doa_high);
    for q=1:num_doa
        %for m = 1:num_doa_high
        temp = 0;
        temp_m = zeros(num_doa_high,1);
            for p = 1: N
                temp = temp + max(R(shiftdelay(p,q,:),p));
%                 temp = temp + sum(R(shiftdelay(p,q,:),p))/num_doa_high;
                temp_m = temp_m + R(shiftdelay(p,q,:),p);
            end
        %end
        srp(q) = temp;
        srp_m(q,:) = temp_m;
    end
    %[~,max_id] = max(srp);
    [max_x,max_y] =find(srp_m==max(max(srp_m)));
    max_id = max_x(1);
    max_m = max_y(1);
    % SRP search frame
%     srp = zeros(num_doa,num_frames);
%     for q=1:num_doa
%         temp = zeros(1,num_frames);
%         for p = 1: N
%              temp = temp + reshape(R(shiftdelay(p,q),p,:),[1,num_frames]);
%         end
%         srp(q,:) = temp;
%     end
%     [~,max_id] = max(srp);

%     imagesc(srp);
    return 
end
function [delay_mic] = compute_delay_mic_loc(mic_loc,mic_a, num_doa,num_doa_high, r)
    fbin = linspace(0,360,num_doa+1);
    f = fbin(1:end-1);
    f_h = linspace(0,60,num_doa_high);
    
    [num_ch, ~] = size(mic_loc);
    N = num_ch*(num_ch-1)/2;
    delay_mic = zeros(N, num_doa, num_doa_high);
    n = 0;
    for i=1:num_ch-1    
        for j = i+1:num_ch
            n = n+1;
            for k = 1:num_doa
                for m = 1:num_doa_high
                    delay_mic(n,k,m) = r * (cosd(f(k)-mic_a(i)) - cosd(mic_a(j)-f(k))) * cosd(f_h(m));
                end
            end
        end
    end
end
% function [max_id]=srp_phat(Sx, num_doa, fs)
%     [nbin, num_frames, num_chs] = size(Sx); 
%     delay_mic = [0.04 0.0283 0 -0.0283 -0.04 -0.0283 0 0.0283;
%                  0.04 0.0566 0.04 0 -0.04 -0.0566 -0.04 0;
%                  0 0.0283 0.04 0.0283 0 -0.0283 -0.04 -0.0283;
%                  0 0.0283 0.04 0.0283 0 -0.0283 -0.04  -0.0283;
%                  -0.04 0 0.04 0.0566 0.04 0 -0.04 -0.0566;
%                  -0.04 -0.0283 0 0.0283 0.04 0.0283 0 -0.0283];
%     shiftTau = delay_mic/340;
%     aug = zeros(num_frames, num_doa);
%     for d=1:num_doa
%          j= 1;
%         for p = 1:num_chs-1
%             for q = p+1: num_chs
%                 
%                 tau = shiftTau(j, d);
%                 aug(:,d) = aug(:,d) + gcc_phat(Sx(:, :, p), Sx(:, :, q), tau, fs, nbin)';
%                 j = j + 1;
%             end
%            
%         end
%         
%     end
%     [~, max_id]= max(sum(real(aug), 1));
%     imagesc(sum(real(aug), 1));
%     return 
% end
% 
% function R=gcc_phat(x1, x2, tau, fs, nbin) 
% f = linspace(0,(fs/2),nbin)';
% R = x1.*conj(x2);
% R = R.*exp(-2*1i*pi*tau*f);
% R = sum(R,1);
% %R = zeros(nfram, nbin);
% end