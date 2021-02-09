clear
clc
close all

priority('h');

%% trial parameters
date_info = now;
cur_date = ['Date-' datestr(date_info,'mmm-dd-yyyy')];
cur_time = ['Time-' datestr(date_info,'HH-MM-SS')];
freq_type = ['Freq-' 'Measure'];

trial_str = strjoin({freq_type,cur_date,cur_time},'_');

%control parameters
save_data = false;

if ~save_data
    write_file = false;
else
    write_file = true;
end

% if contains(freq_type,'SCOPE')
%     save_data = false;
%     write_file = false;
% end

%% Configure audio box
redbox = audioDeviceReader;
audio_driver = redbox.getAudioDevices;

num_channel = 2;

Fs = 48000;
rec_minutes = 5/60;
rec_time = rec_minutes*60; %in seconds

data_type = 'single';


data_size = Fs*rec_time;
data_buffer = nan(data_size,num_channel,data_type);

frame_size = 1024;
frame_time = (0:(frame_size-1))./Fs*1000;
frame_iter = num_channel;
frame_dur = (frame_size/Fs);
frame_count_max = ceil(rec_time./frame_dur);
frame_time_store = nan(frame_count_max,1);

redbox.Driver = "ASIO"; %
redbox.Device = "Focusrite USB ASIO"; %audio_driver{contains(audio_driver,'Focusrite')};
redbox.SampleRate = Fs;
redbox.SamplesPerFrame = frame_size;
redbox.BitDepth = '32-bit float';
redbox.ChannelMappingSource = 'Auto';
redbox.NumChannels = num_channel;
redbox.OutputDataType = data_type;

%quarter in mic sensitivity: 1.36 mV/Pa
%SPL zero = 20 uPa
setup(redbox);

%% configure plotting
plot_dur = 50/1000; %in ms

frame_rate = 30;
frame_time_monitor = 1000/frame_rate;

num_frames_plot = ceil(plot_dur/frame_dur);
plot_buffer = nan(num_frames_plot*frame_size,num_channel,data_type);

plot_sample_max = floor(plot_dur * Fs);
plot_vec = (1:plot_sample_max);

fig_handle = figure(1);
clf(fig_handle);
axes_handle = axes(fig_handle);
hold on

title(trial_str)

xplot_time = (0:1:(num_frames_plot*frame_size-1))./Fs.*1000;
plot_handle1 = plot(xplot_time,zeros(size(xplot_time),data_type),'r');
plot_handle2 = plot(xplot_time,zeros(size(xplot_time),data_type),'b');

xlim([0 max(xplot_time)]);
ylim([-1 1]);

xlabel('Time (ms)')
ylabel('Samples (arb)');



data_idx = 1;
trial_idx = data_idx:(data_idx+frame_size-1);

plot_iter = 1;
plot_idx = plot_iter:(plot_iter+frame_size-1);

fprintf(1,'Starting data collection\n');


control_var = true;
start_time = tic;
while control_var
    cur_time = toc(start_time);
    
    if frame_iter > frame_count_max
        control_var = false;
        
    else
        [audioIn, nOverrun] = record(redbox);
        if nOverrun > 0
            fprintf(1,'Audio device reader queue was overrun by %d samples.\n',...
                nOverrun);
        end
        
        trial_idx = data_idx:(data_idx+frame_size-1);
        data_buffer(trial_idx,1:num_channel) = audioIn;
        
        data_idx = data_idx + frame_size;
        
        frame_time_store(frame_iter) = cur_time;
        frame_iter = frame_iter + 1;
        
        plot_idx = plot_iter:(plot_iter+frame_size-1);
        plot_buffer(plot_idx,1:num_channel) = audioIn;
        plot_iter = plot_iter + frame_size;
        
        if rem(frame_iter,num_frames_plot) == 0
            plot_handle1.YData = plot_buffer(:,1)+0.2;
            plot_handle2.YData = plot_buffer(:,2)-0.2;
            drawnow;
            plot_iter = 1;
        end
        
    end
end
fprintf(1,'Finished data collection\n');

release(redbox);

nan_vals1 = find(isnan(data_buffer(:,1)),1,'first');
if num_channel == 2
    nan_vals2 = find(isnan(data_buffer(:,2)),1,'first');
end

if ~isempty(nan_vals1)
    
    if num_channel == 2
        nan_start = min([nan_vals1 nan_vals2]);
    else
        nan_start = nan_vals1;
    end
    
    data_buffer(nan_start:data_size,:) = [];
end

near_min = ceil(cur_time/60);


if write_file
    aud_file_name = fullfile(pwd,[trial_str '.wav']);
    
    audio_writer = dsp.AudioFileWriter(aud_file_name,...
        'SampleRate',Fs,'DataType',data_type);
    audio_writer(data_buffer);
    release(audio_writer);
end

priority('n');

% segment_size =
%
% data_frame = reshape(data_buffer,frame_size,

fft_data = double(data_buffer(:,2));

fft_pts = 4.*frame_size; %2.*2048;

spec_resolution = Fs/fft_pts;

fft_spec = double(0:(Fs/fft_pts):(Fs-1/fft_pts));
fft_mag = abs(fft(fft_data,fft_pts));

fft_spec = fft_spec(1:fft_pts/2);
fft_mag = fft_mag(1:fft_pts/2);

cal_offset = (104.6-30); %need calibrator to calculate this

fft_plot = 20.*log10(fft_mag)+cal_offset;
fft_plot = movmean(fft_plot,3);

freq_window = [6500 7500];

fft_spec_idx = fft_spec>=min(freq_window) & fft_spec<=max(freq_window);
fft_offset_idx = find(fft_spec<=median(freq_window),1,'first');

[peak_output,peak_output_idx] = max(fft_plot(fft_spec_idx));
freq_output_max = median(freq_window)./1000; %fft_spec(peak_output_idx+fft_offset_idx);

% fft_spec = fft_spec./1000;

figure(2)
clf(2)
hold on

plot(fft_spec./1000,fft_plot)
plot(freq_window./1000,peak_output.*ones(size(freq_window)),'r')

text(freq_output_max,peak_output,sprintf('Peak: %3.2f',peak_output))
% set(gca,'xscale','log')

xlabel('Freq (kHz)')
ylabel('Mag. (dB SPL)')


















