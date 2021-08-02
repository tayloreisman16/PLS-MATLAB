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

rx_block = RX_Process();
obj=rx_block.ConstructorFunc(bandwidth, bin_spacing, num_ant, num_data_bins, bit_codebook, num_subbands);
buffer_rx_data = rx_block.Synchronize(obj, tx);
chan_est_bins_sort = rx_block.ChannelEstimate(obj, buffer_rx_data, ref_sig);
chan_est_sb = rx_block.bins2subbands(obj, chan_est_bins);
[lsv,sval,rsv] = rx_block.SVD(obj, chan_est_sb);
rx_block.PMI_Estimate(obj, rsv);

