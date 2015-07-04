function [] = decode()
evalin('caller','clear all'); 
close all;

global ANT_CNT LTS_LEN SYM_LEN NUM_SYM FFT_OFFSET LTS_CORR_THRESH

DO_CFO_CORRECTION = 1;	% Enable CFO estimation/correction
DO_PHASE_TRACK = 1; % Enable phase tracking
LTS_LEN = 160;
NUM_LTS = 2;
NUM_SYM = 50;
NUM_AC = 0;
LTS_CORR_THRESH = 0.6;
FFT_OFFSET = 1;		% Number of CP samples to use in FFT (on average)
ANT_CNT = 1;
SEGMENT_START = 10000;

% OFDM params
SC_IND_DATA   = [2:7 9:21 23:27 39:43 45:57 59:64]; % Data subcarrier indices
N_SC = 64;          % Number of subcarriers
CP_LEN = 16;        % Cyclic prefix length
SYM_LEN = N_SC + CP_LEN;

% Read tx samples

load('../trace/src_data_1.mat');
%load('src_data_1.mat');

% Read recv samples 
cf = 1;
figure(cf);

rx = read_complex_binary(['../trace/recv_signal.bin']);
%rx = read_complex_binary(['../trace/src_data_1.bin']);
%rx = read_complex_binary(['../trace/recv_signal_v5.bin']);
%rx = read_complex_binary(['../trace/recv_signal_v5.bin']);
%rx = read_complex_binary(['recv_signal2.bin']);
rx = rx(SEGMENT_START:SEGMENT_START+40000-1);
rx_ant = rx;
save(['../trace/src_data_1gen.mat'], 'rx');

plot(real(rx_ant).^2);
raw_title = sprintf( 'Raw Signals %d', i );
title(raw_title);

% Todo
% Packet detection
% implement pkt_detection as the following format

[lts_ind payload_ind, cf] = pkt_detection(rx_ant, LTS_CORR_THRESH, cf);

%payload_ind = 801 
%lts_ind = 481

lts_ind	%display the lts you find
payload_ind %display the payload_ind

% CFO correction
if(DO_CFO_CORRECTION)
	rx_ant = cfo_correction(rx_ant, lts_ind);
end

% Re-extract LTS for channel estimate
cf = cf + 1;
figure(cf);
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
H_ant_lts = zeros(N_SC, ANT_CNT);

rx_lts = rx_ant(lts_ind:lts_ind + LTS_LEN*NUM_LTS - 1);
SC_OFDM = [LTS_LEN - N_SC + 1:LTS_LEN] - FFT_OFFSET;
rx_lts1 = rx_lts(SC_OFDM - N_SC);
rx_lts2 = rx_lts(SC_OFDM);
rx_lts3 = rx_lts(SC_OFDM + N_SC);
rx_lts4 = rx_lts(SC_OFDM + N_SC*2);

rx_lts1_f = fft(rx_lts1);
rx_lts2_f = fft(rx_lts2);
rx_lts3_f = fft(rx_lts3);
rx_lts4_f = fft(rx_lts4);

H_lts1 = rx_lts1_f./ lts_f.';
H_lts2 = rx_lts1_f./ lts_f.';
H_lts3 = rx_lts1_f./ lts_f.';
H_lts4 = rx_lts1_f./ lts_f.';

H_ant_lts = (rx_lts1_f + rx_lts2_f + rx_lts3_f + rx_lts4_f)/4./lts_f.';
H = H_ant_lts;
hold on;
x = [-32:31];
plot(x, real(fftshift(H_lts1)),'r');
%{
plot(x, real(fftshift(H_lts2)),'g');
plot(x, real(fftshift(H_lts3)),'b');
plot(x, real(fftshift(H_lts4)),'k');
%}
plot(x, real(fftshift(H)),'b');
%plot(x, imag(fftshift(H)),'b');
hold off;
grid on;
axis([min(x)+5 max(x)-5 -1.1*max(abs(H)) 1.1*max(abs(H))])
title('Channel Estimates (I-Q)');
xlabel('Subcarrier Index');



rx_ant = rx_ant(payload_ind:payload_ind + SYM_LEN * NUM_SYM - 1);
SC_OFDM = [SYM_LEN - N_SC + 1:SYM_LEN] - FFT_OFFSET;

% Calculate channel estimate
cf = cf + 1;
figure(cf);
H_ant = zeros(N_SC, ANT_CNT);
SC_OFDM = [SYM_LEN - N_SC + 1:SYM_LEN] - FFT_OFFSET;

rx_t1 = rx_ant(SC_OFDM);
rx_t2 = rx_ant(SC_OFDM + SYM_LEN);
rx_t3 = rx_ant(SC_OFDM + 2*SYM_LEN);
rx_f1 = fft(rx_t1);
rx_f2 = fft(rx_t2);
rx_f3 = fft(rx_t3);

H_ant = (rx_f1 + rx_f2 + rx_f3)./ tx_mod_data / 3;
H = H_ant;
hold on;
x = [-32:31];
plot(x, real(fftshift(H)),'r');
plot(x, imag(fftshift(H)),'b');
hold off;
grid on;
axis([min(x)+5 max(x)-5 -1.1*max(abs(H)) 1.1*max(abs(H))])
title('Channel Estimates (I-Q)');
xlabel('Subcarrier Index');


rx_ant = rx_ant(NUM_AC * SYM_LEN + 1:end);

N_OFDM_SYMS = NUM_SYM - NUM_AC;
% Decoding
cf = cf + 1;
figure(cf); clf;
tx_mod_data = repmat(tx_mod_data, 1, N_OFDM_SYMS);

payload_mat = reshape(rx_ant, SYM_LEN, N_OFDM_SYMS);
payload_mat_noCP = payload_mat(SC_OFDM,:);
syms_f_mat = fft(payload_mat_noCP, N_SC, 1);
syms_eq_mat = syms_f_mat./repmat((H_ant), 1, N_OFDM_SYMS);
%syms_eq_mat = syms_f_mat(:,:,i)./repmat(H_ant_lts(:,i), 1, N_OFDM_SYMS);

% Todo
% Phase track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% implement phaseTrack as following usage
if (DO_PHASE_TRACK)
    for sym_i = 1:N_OFDM_SYMS
        syms_eq_pc_mat( :, sym_i ) = phaseTrack( syms_eq_mat( :, sym_i ), tx_mod_data(:, sym_i), cf );
    end
    payload_syms_mat = syms_eq_pc_mat( SC_IND_DATA,: );

else
    payload_syms_mat = syms_eq_mat(SC_IND_DATA,:);
end

%size( syms_eq_pc_mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print 

%figure(10)
%plot(payload_syms_mat,'rx')
%hold on; 
%plot([-1, 1], [0, 0], 'bx');
%plot(rx_ant , 'bx');
%hold off;
%axis([-1.5 1.5 -1.5 1.5]);
%}        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

signal = payload_syms_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Todo
% Calculate average SNR(subcarrier, symbol) here

size(signal) 

noise_pow = abs(signal-tx_mod_data( SC_IND_DATA,: )).^2;

signal_pow = abs(signal).^2 ; 

% SNR = 10*log10( signal_pow./noise_pow );
SNR = 10*log10( signal_pow./noise_pow );

% subcarrier_SNR = ...;	%record average SNR of each subcarrier
subcarrier_SNR =  mean(SNR,2) ;

% SNR = ...;	%record average SNR of each symbol
SNR = mean(SNR) ;
 	     
cf  = cf+1;
figure(cf);
% Plot average SNR of each subcarrier

plot( [[-24:-1] [1:24]], fftshift( subcarrier_SNR ), '-bx' );
xlabel( 'subcarrier' );
ylabel( 'SNR(dB)' );
title( 'Average SNR of each subcarrier' );

cf  = cf+1;
figure(cf);
% Plot average SNR of each subcarrier
plot( [1:50], SNR, '-bo' );
xlabel( 'symbol' );
ylabel( 'SNR(dB)' );
title( 'Average SNR of each symbol' );

cf  = cf+1;
figure(cf);

% Plot constellation points
hold on;
plot(payload_syms_mat,'r.');
axis square; axis(1.5*[-1 1 -1 1]);

plot(tx_mod_data(SC_IND_DATA,:)+0.0001i, 'bo');
legend('Rx','Tx');
title('Tx and Rx Constellations')
grid on;
hold off;

% =================================================================
% CFO correction
% 	-- Use the first antenna to compute CFO --
% =================================================================
function [rx_ant] = cfo_correction(rx_ant, lts_ind)
global ANT_CNT FFT_OFFSET 
% Extract LTS (not yet CFO corrected)
rx_lts = rx_ant(lts_ind:lts_ind+159, 1);
rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);

% Calculate coarse CFO est
rx_cfo_est_lts = mean(unwrap(angle(rx_lts1 .* conj(rx_lts2))));
rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);

%Apply CFO correction to raw Rx waveform
rx_cfo_corr_t = exp(1i*2*pi*rx_cfo_est_lts*[0:length(rx_ant(:,1))-1]');
rx_ant = rx_ant .* repmat(rx_cfo_corr_t, 1, ANT_CNT);

% =================================================================
% initiate settings
% =================================================================
function [] = init_stat()
