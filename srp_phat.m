function [srp, max_id]=srp_phat(Sx, mic_loc, mic_a, num_doa, r, fs)
    %[nbin, num_frames, num_chs] = size(Sx); 
    [nbin, num_chs] = size(Sx); 
%     delay_mic = [0.04 0.0283 0 -0.0283 -0.04 -0.0283 0 0.0283;
%                  0.04 0.0566 0.04 0 -0.04 -0.0566 -0.04 0;
%                  0 0.0283 0.04 0.0283 0 -0.0283 -0.04 -0.0283;
%                  0 0.0283 0.04 0.0283 0 -0.0283 -0.04  -0.0283;
%                  -0.04 0 0.04 0.0566 0.04 0 -0.04 -0.0566;
%                  -0.04 -0.0283 0 0.0283 0.04 0.0283 0 -0.0283];

    delay_mic = compute_delay_mic_loc(mic_loc,mic_a, num_doa, r);
    shiftTau = delay_mic/340;
    maxshift = floor(nbin/2);
    shiftdelay = round(-fs * shiftTau+ maxshift + 1);
    % GCC-PHAT
    N = num_chs * (num_chs - 1)/2;
    %Z= complex(zeros(nbin ,N,num_frames));
    Z= complex(zeros(nbin ,N));
    p= 0;
    for m1=1:num_chs-1
        for m2=m1+1:num_chs
            p= p+1;
            Z(:,p)= Sx(:, m1).*conj(Sx(:,m2));
            Z(:,p) = Z(:,p)./(abs(Z(:,p) +eps));           
%             Z(:,p,:)= Sx(:, :, m1).*conj(Sx(:, :, m2));
%             Z(:,p,:) = Z(:,p,:)./(abs(Z(:,p,:) +eps));
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
    for q=1:num_doa
        temp = 0;
        for p = 1: N
             %temp = temp + sum(R(shiftdelay(p,q),p,:),3);
             temp = temp + (R(shiftdelay(p,q),p,:));
        end
        srp(q) = temp;
    end
    [~, max_id] = max(srp);
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

    imagesc(srp);
    return 
end
function [delay_mic] = compute_delay_mic_loc(mic_loc,mic_a, num_doa, r)
    fbin = linspace(0,360,num_doa+1);
    f = fbin(1:end-1);
    [num_ch, ~] = size(mic_loc);
    N = num_ch*(num_ch-1)/2;
    delay_mic = zeros(N, num_doa);
    n = 0;
    for i=1:num_ch-1    
        for j = i+1:num_ch
            n = n+1;
            for k = 1:num_doa
                delay_mic(n,k) = r * (cosd(f(k)-mic_a(i)) - cosd(mic_a(j)-f(k)));
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