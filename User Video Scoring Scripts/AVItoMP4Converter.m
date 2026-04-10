clear count input_video_file video_file_path video_file_name key_frame videoFReader videoFWrite output_image_path

[ video_file_name,video_file_path ] = uigetfile({'*.avi'});
if(video_file_path==0)
    return;
end
% Output path
output_image_path = fullfile(video_file_path,[video_file_name(1:strfind(video_file_name,'.')-1),'.mp4']);
% mkdir(output_image_path);
input_video_file = [video_file_path,video_file_name];
% Read Video
videoFReader = VideoReader(input_video_file);
% Write Video
videoFWrite = VideoWriter(output_image_path,'MPEG-4');
videoFWrite.FrameRate = videoFReader.FrameRate;
open(videoFWrite);
for count = 1:abs(videoFReader.Duration*videoFReader.FrameRate)    
    fprintf('%d / %d\n',count,videoFReader.Duration*videoFReader.FrameRate)
    key_frame = read(videoFReader,count);
    writeVideo(videoFWrite,key_frame);
end
% Release video object
%close(videoFReader);
close(videoFWrite);
disp('COMPLETED');