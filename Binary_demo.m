%% obtain the gray image
clc; close all;clear all;
Img = imread('110_20184.bmp'); 
Img = mat2gray(rgb2gray(Img));
[row,col] = size(Img);

%% set parameter
r = 40; 
l = 3;
lmada = 0.5;
alpha = 2;
c_u = 0.7;
c_W = 0.1;


%% start level set evolution
u =  binary_test(Img,c_u,r,l,lmada,c_W,alpha);


%% figure
figure;imshow(u,[]);title('final u');colorbar
figure;imshow(u>0,[]);title('binary result');colorbar


%% save result
save_path = 'test_one/'; gt_path = ['gt\']; di2 = dir(gt_path);
number = 1;
u_temp = 255*double(u>0);
imwrite(u_temp,[save_path,num2str(number),'.bmp']);

%% computer measure
path = cell(1,4);
path{1,1} = [save_path,num2str(number),'.bmp'];
path{1,2} = [gt_path,di2((number)*3).name];
path{1,3} = [gt_path,di2((number)*3+1).name];
path{1,4} = [gt_path,di2((number)*3+2).name];
result = system(['DIBCO_metrics',' ', path{1,2},' ',path{1,1},' ',path{1,4},' ',path{1,3}]);








