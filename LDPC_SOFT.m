%
% LDPC performance for any LDPC Parity Matrix and Generator Matrix
% Authors: Biniyam Zewede and Shem Kikamaze

% mex soft_decoder.cpp
clc;
clear all;
close all;
H = dlmread('H');
G = dlmread('G');
[m, n] = size(H);
[k, ~]= size(G);

hMod = comm.PSKModulator(2, 'BitInput',true);
hChan = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (SNR)','SNR',0);
hDemod = comm.PSKDemodulator(2, 'BitOutput',true,...
    'DecisionMethod','Approximate log-likelihood ratio', ...
    'Variance', 1/10^(hChan.SNR/10));

EbNoDb = -2:0.5:10;
EbNo = 10.^(EbNoDb/10);
EcNo = EbNo*(k/n);
EcNoDb = 10*log10(EcNo);
Averaging = 100;
error = zeros(1,length(EbNo));
bit_error_mat = zeros(1,length(EbNo));
no_code_error = error;

iter = 100;
for a = 1:Averaging
    for i = 1:length(EbNo)
        b = round(rand(1,k));
        c = mod((b*G),2);
        b = c(n-k+1:n);
        Modsignal = step(hMod, c');
        set(hChan,'SNR',EcNoDb(i));
        set(hDemod,'Variance',1/10^(hChan.SNR/10));
        receivedsignal = step(hChan, Modsignal);
        llr = (step(hDemod, receivedsignal))';
        hard_decision = (sign(llr)-1)./(-2);
        error_before = sum(abs(hard_decision(n-k+1:n)-c(n-k+1:n)));
        final_llr = soft_decoder(H,llr,n,n-m,iter);
        final_bits_2 = sign(final_llr);
        c_new = (final_bits_2-1)./(-2);
        b_new = c_new(n-k+1:n);
        error_after = sum(abs(b_new-b));
        no_code_error(i) = no_code_error(i) + error_before;
        error(i) = error(i)+ error_after;
    end
end
CER = error/(k*Averaging);
semilogy(EbNoDb,CER,'b-');
ylabel('Probability of bit error');
hold on
grid on
theory = qfunc(sqrt(2*EbNo));
semilogy(EbNoDb,theory, 'ko');
clear functions
