classdef RX_Process
    
    properties
        bandwidth = 960000; %import
        bin_spacing = 15000; %import
        num_ant = 2; %import
        bit_codebook = 2; %import
        NFFT = 64;
        CP = 16;
        OFDMSymb_len = 80;
        num_data_bins = 4; %import
        used_data_bins = [30 31 33 34]; %import - MATLAB differences?
        codebook = [];
        subband_size = 2; % how to tie to num_ant?
        num_bits_symb = (4 / 2) * 1;
        num_subbands = 2; %import
        num_PMI = 2; % num_subbands
        synch = []; % Generate? In RX? Or in another file?
        max_impulse = 64; % I thought this was the length of the CP.
        total_num_symb = 12; % Import
        %         total_symb_len = total_num_symb * OFDMSymb_len;
        total_symb_len = 12 * 80;
        num_data_symb = 4; % Calculate/Import
        num_synch_symb = 8; % Calculate/Import
        SNRdB = 100; % Need to import into function
        %         SNR_type = SNR_type; % Need to import into function
        num_patterns = 3;
        synch_data_pattern = []
        synch_pattern = []
        data_pattern = []
        symb_pattern0 = []
        symb_pattern = []
        
        % symb_pattern how to tile?
    end
    
    methods   %constructor method
        function obj = RX_Process(varargin)
            
            for i = 1: nargin
                switch i
                    case 1
                        obj.bandwidth = varargin{1};
                    case 2
                        obj.bin_spacing = varargin{2};
                        
                    case 3
                        obj.num_ant = varargin{3};
                    case 4
                        obj.num_data_bins = varargin{4};
                    case 5
                        obj.bit_codebook = varargin{5};
                    case 6
                        obj.num_subbands = varargin{6};
                end
            end
            dbg = 1;
            obj.NFFT = floor(obj.bandwidth/obj.bin_spacing);
            obj.CP = obj.NFFT * 0.25;
            obj.OFDMSymb_len = obj.NFFT + obj.CP;
            obj.used_data_bins = [30 31 33 34]; % we will rewrite this
            obj.subband_size = obj.num_ant;
            obj.num_bits_symb = (obj.num_data_bins / obj.num_ant) * obj.bit_codebook;
            obj.num_PMI = obj.num_subbands;
            obj.max_impulse = obj.NFFT;
            
            
            obj.SNRdB = 100;
            obj.num_patterns = 3; % number of times the synch data pattern repeats
            obj.synch_data_pattern = [2, 1];
            obj.synch_pattern = zeros(1, obj.synch_data_pattern(1));
            obj.data_pattern = ones(1, obj.synch_data_pattern(2));
            obj.symb_pattern0 = [obj.synch_pattern, obj.data_pattern];
            obj.symb_pattern = repmat(obj.symb_pattern0, 1, 3); % Make 3 adjustable
            obj.num_data_symb = sum(obj.symb_pattern); % need to automate
            
            obj.num_synch_symb = length(obj.symb_pattern) - obj.num_data_symb; % need to automate
            obj.total_num_symb = obj.num_data_symb + obj.num_synch_symb;
            obj.total_symb_len = obj.total_num_symb * obj.OFDMSymb_len;
            
            
            
            num_precoders = 2^(obj.bit_codebook);
            obj.codebook = cell(1, num_precoders);
            for p = 0: num_precoders-1
                precoder = zeros(obj.num_ant, obj.num_ant); % assuming Alice is tx Bob is Rx
                for m = 0: obj.num_ant-1 % columns in each precoder matrix
                    for n = 0: obj.num_ant-1 % rows in that column
                        w = exp(1j*2*pi*n/obj.num_ant*(m + p/num_precoders));
                        precoder(n+1, m+1)= (1/sqrt(obj.num_ant))*w;
                    end
                end
                
                
                obj.codebook{p+1} = precoder;
                
            end
            
            
        end
        
        
        
        function buffer_rx_data=Synchronize(obj, rx_time_buffer)
            total_symb_count = 0;
            synch_symb_count = 0;
            data_symb_count = 0;
            
            buffer_rx_data = [];
            
            for symb = obj.symb_pattern
                symb_start = (total_symb_count * obj.OFDMSymb_len) + 1;
                symb_end = symb_start + obj.OFDMSymb_len-1;
                
                if symb == 0
                    synch_symb_count = synch_symb_count + 1;
                else
                    data_start = (data_symb_count * obj.NFFT) + 1;
                    data_end = data_start + obj.NFFT - 1;
                    
                    data_with_CP = rx_time_buffer(:, symb_start: symb_end);
                    data_without_CP = data_with_CP(:, obj.CP + 1:obj.CP+obj.NFFT);
                    buffer_rx_data(:, data_start: data_end) = data_without_CP;
                    data_symb_count = data_symb_count + 1;
                end
                total_symb_count = total_symb_count + 1;
            end
        end
        
        function chan_est_bins_sort=ChannelEstimate(obj, buffer_rx_data, ref_signal)
            
            
            chan_est_bins_sort = zeros(obj.num_ant, obj.num_data_symb, obj.num_data_bins/obj.subband_size, obj.subband_size);
            chan_est_bins = zeros(obj.num_ant, obj.num_data_symb * obj.num_data_bins);
            
            for symb = 1:(obj.num_data_symb)
                symb_start = (symb-1) * obj.NFFT + 1;
                symb_end = symb_start + obj.NFFT - 1;
                
                used_symb_start = (symb - 1) * obj.num_data_bins + 1;
                used_symb_end = used_symb_start + obj.num_data_bins - 1;
                for ant = 1:obj.num_ant
                    time_data = buffer_rx_data(ant, symb_start:symb_end);
                    data_fft = fft(time_data, obj.NFFT);
                    data_in_used_bins = data_fft(obj.used_data_bins);
                    
                    est_channel = data_in_used_bins * (conj(ref_signal(symb, :))/(abs(ref_signal(symb, :))))';
                    chan_est_bins(ant, used_symb_start: used_symb_end) = est_channel;
                    for subband_index = 1: (obj.num_data_bins / obj.subband_size)
                        start = (subband_index - 1) * obj.subband_size + 1;
                        stop = start + obj.subband_size - 1;
                        chan_est_bins_sort(ant, symb, subband_index, :) = est_channel(start: stop);
                    end
                end
            end
        end
        
        function chan_est_sb=bins2subbands(obj, chan_est_bins)
            chan_est_sb = [];
            for symb = 1: (obj.num_data_symb + 1)
                symb_start = (symb - 1) * obj.num_data_bins;
                symb_end = symb_start + obj.num_data_bins;
                chan_est = chan_est_bins(:, symb_start: symb_end);
                for sb = 1: (obj.num_subbands + 1)
                    sb_start = (sb - 1) * obj.subband_size;
                    sb_end = sb_start + obj.subband_size;
                    
                    chan_est_sb(symb, sb) = chan_est(:, sb_start: sb_end);
                end
            end
        end
        %         function [lsv,sval,rsv] = SVD(chan_est_sb)
        %             lsv = [];
        %             sval = [];
        %             rsv = [];
        %
        %             for symb = 1:(obj.num_data_symb + 1)
        %                 for sb = 1:(obj.num_subbands + 1)
        %                     [U,S,VH] = svd(chan_est_sb(symb, sb));
        %                     V = conj(VH)';
        %                     ph_shift_u = diag(exp(-1j * angle(U(1, :))));
        %                     ph_shift_v = diag(exp(-1j * angle(V(1, :))));
        %                     lsv(symb, sb) = dot(U, ph_shift_u);
        %                     sval(symb, sb) = S;
        %                     rsv(symb, sb) = dot(V, ph_shift_v);
        %                 end
        %             end
        %         end
        %         function [PMI_sb_estimate, bits_sb_estimate]=PMI_Estimate(rx_precoder)
        %             PMI_sb_estimate = [];
        %             bits_sb_estimate = [];
        %             for symb = 1:(obj.num_data_symb + 1)
        %                 for sb = 1:(obj.num_subbands + 1)
        %                     dist = 0.0;
        %                     for prec = 1:(obj.codebook_sz + 1)
        %                         diff = rx_precoder(symb, sb) - obj.codebook(prec);
        %                         diff_squared = real(diff * conj(diff));
        %                         dist(prec) = sqrt(diff_squared.sum());
        %                     end
        %                     min_dist = min(dist);
        %                     PMI_estimate = find(dist == min_dist);
        %                     PMI_sb_estimate(symb, sb) = PMI_estimate;
        %                     bits_sb_estimate(symb, sb) = de2bi(PMI_estimate, obj.bit_codebook);
        %                 end
        %             end
        %         end
    end
end