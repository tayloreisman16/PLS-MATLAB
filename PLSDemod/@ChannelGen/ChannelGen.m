classdef ChannelGen < handle
    
    properties
        codebook=[]; %null output vector
        codebook_sz = 0;
    end
    
    methods   %constructor method
        function obj=ChannelGen(codebook_sz)
            obj.codebook_sz = codebook_sz;
            num_precoders = 2^obj.codebook_sz;
            matrix0 = [1+0i 1+0i; 1+0i -1+0i];
            matrix1 = [1+0i 1+0i; 0+1i 0-1i];
            scale = 1 / sqrt(2);
            e = cat(3, matrix0, matrix1);
            obj.codebook = e * scale;
        end
    end   
end
