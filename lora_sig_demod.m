function lora_sig_demod(fname,maximum_payload_len,pwr_threshold_ratio)

%% parameters
d_sf = 11; % spreading factor
d_symbol_preamble = 6; % number of symbols for preamble

d_fft_size = 2^d_sf; % fft size
d_sampling_rate = 2e6; % sampling rate [samples per second]
d_bandwidth = 125e3; % RF bandwidth
d_resample_rate = d_bandwidth/d_sampling_rate;
d_len_payload = maximum_payload_len; % suppose a maximum length of payload in symbol number
d_total_symbol = d_symbol_preamble + 2 + 2.5 + d_len_payload;
d_len_frame = d_total_symbol*d_fft_size; % samples in sampling rate

%% Generate upchirp and downchirp
accumulator = 0;
phase = -pi;
sig_upchirp = zeros(2*d_fft_size,1);
for ii= 1:2*d_fft_size
    accumulator = accumulator + phase;
    sig_upchirp(ii) = exp(j*accumulator);
    phase = phase + 2*pi/d_fft_size;
end
sig_downchirp = zeros(2*d_fft_size,1);
for ii= 1:2*d_fft_size
    accumulator = accumulator + phase;
    sig_downchirp(ii) = exp(-j*accumulator);
    phase = phase + 2*pi/d_fft_size;
end

%% Set data file
fid = fopen(fname, 'rb');
len_read = 2*d_len_frame/d_resample_rate; % 2 times of frame length

% % Manually cut off some data at the beginning of the data file
% % Let the start part of signal is background noise
% [data, count] = fread(fid, 2*10^6, 'int8');

% set a sig buffer
sig_buffer = zeros(d_len_frame/d_resample_rate,1);
buffer_flag = 0; % 0 - no data; 1 - has data
pos_buffer_write = 1;

%% start to read data
while(1)
    [data, count] = fread(fid, len_read, 'int8');
    disp(['==== Read ', num2str(count), ' samples ====']);
    if (count < len_read)
        disp(['==== Reach end of file ====']);
        break; % stop reading
    end
    
    sig = data(1:2:end) + j*data(2:2:end);
    sig_len = length(sig);
    
    %% Select the valid data part
    if (buffer_flag==0)
        % locate the start position of the following frame
        avg_abs = mean(abs(sig(1:100))); % Suppose the beginning 100 samples are noise FIXME
        index = find(abs(sig)>pwr_threshold_ratio*avg_abs); % Signal strength reaches the threshold
        if (length(index)==0)
            continue; % No valid data
        end
        start_pos = index(1)-400/d_resample_rate; % shift a little to guarantee including the header
        start_pos = max(1,start_pos);
        if ((sig_len-start_pos+1)<d_len_frame/d_resample_rate)
            sig_buffer(1:sig_len-start_pos+1) = sig(start_pos:end);
            buffer_flag = 1;
            pos_buffer_write = sig_len-start_pos+1 +1;
            len_read = (d_len_frame/d_resample_rate - (sig_len-start_pos+1))*2;
            continue; % read the rest of this frame
        else
            sig_buffer = sig(start_pos:start_pos+d_len_frame/d_resample_rate-1);
        end
    else
        sig_buffer(pos_buffer_write:end) = sig;
        buffer_flag = 0;
        pos_buffer_write = 1;
        len_read = 2*d_len_frame/d_resample_rate;
    end
    sig_frame = sig_buffer;
    disp(['==== Info: Found signal ====']);
    figure();
    plot(abs(sig_frame));
    
    %% Downsampling
    sig_downsampled = downsample(sig_frame,1/d_resample_rate);
    sig_frame = sig_downsampled;
    
    %% synchronization
    pos_start = 0;
    for ii = 1:1000
        sig_window = sig_frame(ii:ii+d_fft_size-1);
        fft_input = sig_window.*sig_downchirp(1:d_fft_size);
        [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
        if (max_index==1)
            pos_start = ii;
            break;
        end
    end
    if (pos_start==0)
        disp(['==== Error: Did not find preamble ====']);
        continue;
    end
    
    %% Verify the preamble part
    disp(['==== Verify Preamble Part ====']);
    sig_frame = sig_frame(pos_start:end);
    for symbol_index = 1:d_symbol_preamble
        sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size+d_fft_size);
        fft_input = sig_window.*sig_downchirp(1:d_fft_size);
        [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
        disp(['The ', num2str(symbol_index),' th symbol peak index: ',num2str(max_index)]);
    end
    % caculate average abs
    avg_abs = mean(abs(sig_window));
    
    %% Check sync words
    disp(['==== Sync Word Part ====']);
    symbol_index = d_symbol_preamble+1;
    sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size+d_fft_size);
    fft_input = sig_window.*conj(sig_upchirp(1:d_fft_size));
    [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
    sync_word_0 = round((max_index-1)/8); % sometimes it is not an integer
    
    symbol_index = d_symbol_preamble+2;
    sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size+d_fft_size);
    fft_input = sig_window.*sig_downchirp(1:d_fft_size);
    [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
    sync_word_1 = round((max_index-1)/8);
    disp(['The Sync Word is 0x ', dec2hex(sync_word_0),' ',dec2hex(sync_word_1)]);
    
    %% SFD part
    disp(['==== Start of Frame Delimiter Part ====']);
    symbol_index = d_symbol_preamble+3;
    sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size+d_fft_size);
    fft_input = sig_window.*conj(sig_downchirp(1:d_fft_size));
    [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
    disp(['The ', num2str(symbol_index),' th symbol peak index: ',num2str(max_index)]);
    
    symbol_index = d_symbol_preamble+4;
    sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size+d_fft_size);
    fft_input = sig_window.*sig_upchirp(1:d_fft_size);
    [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
    disp(['The ', num2str(symbol_index),' th symbol peak index: ',num2str(max_index)]);
    
    % the 0.25 part
    symbol_index = d_symbol_preamble+5;
    sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size+d_fft_size*0.25);
    fft_input = sig_window.*sig_upchirp(1:d_fft_size*0.25);
    [max_value,max_index] = max(abs(fft(fft_input,d_fft_size*0.25)));
    disp(['The ', num2str(symbol_index),' th symbol peak index: ',num2str(max_index)]);
    
    %% Payload part
    disp(['==== Start of Payload Part ====']);
    for payload_symbol_index = 1:d_len_payload
        symbol_index = d_symbol_preamble +4.25 +payload_symbol_index;
        sig_window = sig_frame((symbol_index-1)*d_fft_size+1:(symbol_index-1)*d_fft_size +d_fft_size); % including 0.25 part
        tmp_abs = mean(abs(sig_window));
        if (tmp_abs < avg_abs*0.5) % set 0.5 as threshold
            disp(['End of payload']);
            break;
        end
        fft_input = sig_window.*sig_downchirp(1+0.25*d_fft_size:d_fft_size+0.25*d_fft_size); % shifted 0.25 symbol
        [max_value,max_index] = max(abs(fft(fft_input,d_fft_size)));
        disp(['The ', num2str(payload_symbol_index),' th payload symbol: ',num2str(max_index), ' in Hex 0x', dec2hex(max_index,4)]);
    end
    
    
end

fclose(fid);
end


