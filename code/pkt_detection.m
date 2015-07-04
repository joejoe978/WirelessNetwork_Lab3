function [lts_ind payload_ind cf] = pkt_detection(rx_ant, THRESH_LTS_CORR, cf)
%LTS_LEN = 160;
LTS_LEN = 320;

% Long preamble (LTS) for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

%lts = [lts_t(33:64) lts_t lts_t lts_t lts_t lts_t(1:32)]

% Hint 1:
% Complex cross correlation of Rx waveform with time-domain LTS 
lts_corr = [];
delta = length(lts_t) -1

for i = 1:(length(rx_ant) - delta)
   max_corr = max(abs(xcorr(rx_ant(i:i+63), lts_t)));
   lts_corr(i) = max_corr;
end

%{
det = length(lts)-1
for i = 1:(length(rx_ant)-det)
   tmp = max(abs(xcorr(rx_ant(i:i+det), lts)));
   lts_corr(i) = tmp;
end
%}

% Skip early and late samples
lts_corr = lts_corr(32:end-32);

cf = cf + 1
figure(cf);
plot(lts_corr);

% Hint 2:
% Get the lts_ind and payload_ind
% Don't forget that "You have skip early and late samples"

[maxVal, max_index] = max(lts_corr);
maxVal
max_index

%if packet is the last one
if((max_index-480) < 1) 
    lts_corr = lts_corr(max_index+4000:end);
    [maxVal, max_index] = max(lts_corr);
elseif ((max_index+4288) > length(lts_corr)) 
    lts_corr = lts_corr(1:max_index-320);
    [maxVal, max_index] = max(lts_corr);
end

cf = cf + 1
figure(cf);
plot(lts_corr);

%which peak is max
%max_index = max_index -1

tmp = max_index;
while 1
    if(abs(lts_corr(tmp)-lts_corr(tmp-64)) > 0.1)
        max_index = tmp;
        break
    else
        tmp = tmp - 64;
    end
end
        
lts_ind = max_index
payload_ind = max_index + LTS_LEN

