% Sanjay R. Kharche.
% I used this in 2019 to make sure I had enough, and that enough was Gaussian.
% See that all your random numbers are centered (peak) around the mu you want, and that they are symmetric Gaussians around the mu.
 % 18 November 2020: This is how you may do it, its not a complete program.
 %
clear all
clear all
close all
close all
%
datat = load('job_file');
for i=1:1:11
col1 = datat(:,i);
h1 = histogram(col1, 400,'EdgeColor','none','FaceColor', [0 0 1]);
clear col1;
pause(2);
end;
