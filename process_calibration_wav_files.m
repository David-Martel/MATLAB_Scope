
clear
clc
close all

file_list = getAllFiles(pwd);

wav_idx = contains(file_list,'.wav');
measure_idx = contains(file_list,'ChirpMeasure');
measure_file = file_list(wav_idx & measure_idx);

test_file = measure_file{1};
[measure_data,Fs] = audioread(test_file);

source_idx = contains(file_list,'protocolTime');
source_file = file_list(wav_idx & source_idx);
[source_data,source_Fs] = audioread(source_file{1});

% test_url = 'https://www.youtube.com/watch?v=FM7MFYoylVs&list=RDFM7MFYoylVs&start_radio=1';


% aud_envelope = abs(hilbert(aud_data));

frame_size = 1024;
frame_time = frame_size/Fs;

% sig_mean = mean(aud_envelope);
% sig_std = std(aud_envelope);
% 
% aud_zscore = (aud_envelope-sig_mean)./sig_std;


% [pks,locs] = findpeaks(abs(aud_envelope),'MinPeakHeight',0.4);

group_frames = 1;

measure_frames = reshape(measure_data,[],frame_size*group_frames);
[num_measure_frames, group_frame_size] = size(measure_frames);

source_frames = reshape(source_data,[],group_frame_size);

[num_source_frames, ~] = size(source_frames);


time_data = (0:(group_frame_size-1))./Fs*1000;

figure(2)
clf(2)
hold on

measure_frame_max = max(abs(measure_frames),[],2);
source_frame_max = max(abs(source_frames),[],2);

plot(measure_frame_max,'b')
plot(source_frame_max,'r')


%% plot measure data
% fig_handle = figure(1);
% clf(1)
% 
% axes_handle = axes(fig_handle);
% line_handle = plot(time_data,measure_frames(1,:),'color',[1 0 0]);
% 
% axes_handle.YLim = [-0.5 0.5];
% 
% 
% for frame_iter = 1:num_frame
%     
%     cur_time = sprintf('%3.3f sec',frame_iter * frame_time * group_frames);
%     
%     line_handle.YData = measure_frames(frame_iter,:);
%     
%     axes_handle.Title.String = cur_time;
%     
%     drawnow;
%     java_pause(0.1);
% end

    

