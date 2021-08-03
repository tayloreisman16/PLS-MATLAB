clear
dbstop if error
clc

load('Bob2Alice.mat');
load('ref_sig.mat');
bandwidth = 960000;
bin_spacing = 15000;
num_ant = 2;
num_data_bins = 4;
num_subbands = 2;
total_num_symbols = 12;
bit_codebook = 1;

rx_block = RX_Process(bandwidth, bin_spacing, num_ant, num_data_bins, bit_codebook, num_subbands);

buffer_rx_data = Synchronize(rx_block, tx);
chan_est_bins_sort = ChannelEstimate(rx_block, buffer_rx_data, ref_sig);
% chan_est_sb = bins2subbands(rx_block, chan_est_bins_sort); 
% Start with this
% [lsv,sval,rsv] = SVD(rx_block, chan_est_sb);
% PMI_Estimate(rx_block, rsv);

