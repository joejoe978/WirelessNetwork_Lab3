%%This file is for single-antenna client

clear;
close all;
%Systems params

ANTENNA_CNT = 1;

%Params:
USE_WARPLAB_TXRX = 0;   %Enable WARPLab-in-the-loop (otherwise sim-only)
WRITE_PNG_FILES = 0;    %Enable writing plots to PNG

%Waveform params
N_OFDM_SYMS     = 50;	%Number of OFDM symbols
N_AC			= 0;	%Number of access codes
MOD_ORDER       = 2;	%Modulation order (2/4/16 = BSPK/QPSK/16-QAM)
TX_SCALE        = 1.0;	%Scale for Tx waveform ([0:1])
INTERP_RATE     = 2;	%Interpolation rate (1 or 2)
SIGNAL_ROUND    = 1;

%OFDM params
SC_IND_PILOTS = [8 22 44 58];   %Pilot subcarrier indices
SC_IND_DATA   = [2:7 9:21 23:27 39:43 45:57 59:64]; %Data subcarrier indices
N_SC = 64;          %Number of subcarriers
CP_LEN = 16;        %Cyclic prefix length
N_DATA_SYMS = N_OFDM_SYMS * length(SC_IND_DATA); %Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
USE_PILOT_TONES = 1;    %Enabel phase error correction

%% Define the preamble
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
a = sts_t;
sts_t = sts_t(1:16);


%LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);
b = lts_t;
%Use 30 copies of the 16-sample STS for extra AGC settling margin
sts = repmat(sts_t, 1, 30);
%lts = [lts_t(33:64) lts_t lts_t];
lts = [lts_t(33:64) lts_t lts_t lts_t lts_t lts_t(1:32)];
%preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];

%return;

%% Generate a payload
%ac = randi(MOD_ORDER, length(SC_IND_DATA), N_AC) - 1;
data = randi(MOD_ORDER, length(SC_IND_DATA), 1) - 1;
%tx_data = [ac repmat(data, 1, N_OFDM_SYMS - N_AC)];
tx_data = repmat(data, 1, N_OFDM_SYMS);



%Functions for data -> complex symbol mapping (avoids comm toolbox requirement for qammod)
modvec_bpsk =  (1/sqrt(2))  .* [-1 1];
modvec_16qam = (1/sqrt(10)) .* [-3 -1 +3 +1];

mod_fcn_bpsk = @(x) complex(modvec_bpsk(1+x),0);
mod_fcn_qpsk = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));

%Map the data values on to complex symbols
switch MOD_ORDER
    case 2 %BPSK
        tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
    case 4 %QPSK
        tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
    case 16 %16-QAM
        tx_syms = arrayfun(mod_fcn_16qam, tx_data);      
    otherwise
        fprintf('Invalid MOD_ORDER (%d)!\n', MOD_ORDER);
        return;
end


%Duplicate multiple data symbols
tx_syms_mat = tx_syms;

%Define the pilot tones
if(USE_PILOT_TONES)
    pilots = [1 1 -1 1].';
else
    pilots = [0 0 0 0].';
end

%Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);
size(pilots_mat)

%% IFFT

%Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);


%Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :) = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;

%Perform the IFFT
tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);

%Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
    tx_payload_mat = [tx_cp; tx_payload_mat];
end

%Reshape to a vector
tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));

for ant_i = 1:ANTENNA_CNT
    ANTENNA_CNT
    from = length(sts) + length(lts) * (ant_i-1) + 1;
    to   = length(sts) + length(lts) * ant_i;
    %Construct the full time-domain OFDM waveform
    preamble = [sts, lts];
    preamble = preamble ./ max(abs(preamble));
	max(abs(tx_payload_vec))
    tx_payload_vec = tx_payload_vec ./ max(abs(tx_payload_vec));
	max(abs(tx_payload_vec))
    tx_vec = [preamble tx_payload_vec];
    tx_vec = 0.1*tx_vec;
    %Write signals to the binary file
    fid = fopen(['src_data_' int2str(ant_i) '.bin'],'w');
    fwrite(fid, [real(tx_vec); imag(tx_vec)], 'float');
    fclose(fid);
	
	%Record samples to the mat file for decoding
	tx_mod_data = ifft_in_mat(:, 1:N_AC+1);
	tx_mod_data = reshape(tx_mod_data, numel(tx_mod_data), 1);
	tx_data = reshape(tx_data, 1, numel(tx_data));
	tx_syms = reshape(tx_syms, 1, numel(tx_syms));
	save(['src_data_' int2str(ant_i) '.mat'], 'tx_mod_data', 'tx_data', 'tx_syms');
    
    %Show Signals
    figure;
    plot(abs(tx_vec));
end

return;
