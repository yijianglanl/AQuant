% Supplementary Matlab code and example data analysis

% Alignment vector analysis method
% Matlab code was implemented by Byoungkoo Lee 
% Matlab version: R2014a, 8.3.0.532

% Multiphoton microscopy images (.tif format) of non-small cell lung 
% carcinoma cell-induced collagen fibers were taken by Jessica Konen
% Image size: 512 x 512 pixels, 1 pixel = 0.83 micro-meter.
% Z stack interval: 1 micro-meter

% Preprocessing step 1: 
% Using ImageJ (1.47v), the brightness/contrast of images are adjusted to
% more clearly see the collagen fibers, using "auto" adjustment option.

% Preprocessing step 2: 
% Using CT-FIRE (win64 Beta 1.2.1, http://loci.wisc.edu/software/ctfire),
% collagen fiber segments were extracted. 
% CT-FIRE parameters
% thresh_im2:15
% s_xlinkbox:4 
% thresh_ext:15 
% thresh_dang_L:5
% thresh_short_L:5
% s_fiberdir:4
% thresh_linkd:3
% thresh_linka:-165, which is equal to 15.
% thresh_flen:5

% curvelet filter parameters in CT-FIRE for denoising
% Percentile of the remaining curvelet coeffs: 0.2
% Number of selected scales: 2

% Output figure control (optional) in CT-FIRE
% Minimum fiber length: 5
% Maximum fiber number: 2999
% Fiber line width: 0.5
% Histogram bins number: 10

% CT-FIRE parameters were selected by massive tests for varying each value
% separately. We quantized fibers into the smallest pixel length but larger
% than noise-like signals. We selected 5 pixel length for the quantization
% process.

% example: H1299, shLKB1, 21hour after seeding in a collagen gel.
%          Z=1 (the bottom of the gel) image.
% shLKB1_21hr_z01c2.tif: example image after preprocessing step 1
% ctFIREout_shLKB1_21hr_z01c2.mat: Matlab data after preprocessing step 2

% Execute the following two Matlab commands, if Matlab return error messages 
% about Matlab toolboxs and paths.
% "Error using Settings
% Undefined function or variable 'matlab.internal.getSettingsRoot'."
restoredefaultpath;
rehash toolboxcache;

clear all;

% Local Alignment Vectors
% Step1: Find the tumor boundary circle manually
% Display the multiphoton microscopy image with extracted fiber
% segments and then select the circular tumor boundary

% change the folder path if you save in the different folder
imgpath = 'D:\Local_Alignment_Vector\';
filepath = 'D:\Local_Alignment_Vector\';

% I used the first z stack image (Z=1), for the demonstration.
% It should extend whole z stacks.
fname=1; 
imgname = ['shLKB1_21hr_z0' num2str(fname) 'c2.tif'];
filename =['ctFIREout_shLKB1_21hr_z0' num2str(fname) 'c2.mat'];
           
fullname2 =  fullfile(imgpath,imgname);
info = imfinfo(fullname2);
pixw = info(1).Width;  
pixh = info(1).Height;
image=imread(fullname2);
image2=image;
for i=1:length(image2(1,:))
    for j=1:length(image2(:,1))
        yp=abs(j-pixh-1);
        image(yp,i)=image2(j,i);
    end
end
   
% A rough circular tumor boundary was manually selected.    
% For the H1299 shLKB1 21hr, the total Z stack number is 76.
% For this example, we used only the first Z stack (Z=fname=1). 
% Tumor boundary can change through Z axis.
tumorradius=265-fname*(265-205)/76; 
tumorcenterx=-65-fname*(-65+50)/76;
tumorcentery=325-fname*(325-365)/76;

% Display the example multiphoton collagen image with extracted fibers and
% selected the circular tumor boundary.
figure();
imshow(image); % display the example collagen image
axis xy; 
hold on;

% display the selected circular tumor boundary
circleang=[0:0.01:2*pi];
circlex=[];
circley=[];
for i=1:length(circleang)
    boundaryx=cos(circleang(i))*tumorradius+tumorcenterx;
    boundaryy=sin(circleang(i))*tumorradius+tumorcentery;
    if (boundaryx>0 && boundaryx<512 && boundaryy>0 && boundaryy<512)
        circlex=[circlex; boundaryx];
        circley=[circley; boundaryy];
    end
end
size1=length(circlex(:,1));
for bline=2:length(circlex(:,1))
    tempxdelta=abs(circlex(bline-1)-circlex(bline));
    tempydelta=abs(circley(bline-1)-circley(bline));
    if(tempxdelta<10 && tempydelta<10)
        plot(circlex(bline-1:bline, 1), circley(bline-1:bline, 1), '-y','linewidth',1.5); 
    end
end
        
% display extracted fiber segments
fullname1 =  fullfile(filepath,filename);
load(fullname1);
FN = find(data.M.L>0); 
FLout = data.M.L(FN);
LFa = length(FN);         
for LL = 1:LFa
    VFa.LL = data.Fa(1,FN(LL)).v;    % field data.Fa includes the index of the points of each fiber 
    XFa.LL = data.Xa(VFa.LL,:);      % field data.Xa includes the coordinate of all the points
    segnum=length(XFa.LL(:,1));
    plot(XFa.LL(:,1),abs(XFa.LL(:,2)-pixh-1), '-c','linewidth',1); 
end
title('Test example: H1299 shLKB1, 21hr, Z=1','FontName','Arial','FontWeight','Bold','FontSize',12);
saveas(gcf,'Step1_Fig3b.tif');
close(gcf);

% Step2: Quantization of extracted fibers
% For 512 x 512 pixels image, we chose 5 pixel length for the quantization
% of extracted fiber segments, which is the smallest fiber length to maximize
% local fiber feature but larger than noise-like signals in our images. 
% Input: CT-FIRE Data
% Output: angpos, quantized fiber segments

% minimum fiber segment length for quantization = 5 pixels
minfsegleng=5; 

% Quantized fiber position data
% Each extracted fiber segments were quantized by 5 pixel length, fitted
% into a straight line: 
% X position and Y position for one end of the line, 
% X position and Y position for the other end of the line
angpos=[]; % x1, y1, x2, y2, fiber segment length

for LL = 1:LFa
    VFa.LL = data.Fa(1,FN(LL)).v;    % field data.Fa includes the index of the points of each fiber 
    XFa.LL = data.Xa(VFa.LL,:);      % field data.Xa includes the coordinate of all the points
    segnum=length(XFa.LL(:,1));
    
    % exclude the fiber segment if it located inside of the tumor boundary circle.
    exclude=0;
    for i=1:segnum
        temppos=[XFa.LL(i,1) abs(XFa.LL(i,2)-pixh-1)];
        dist=norm(temppos-[tumorcenterx, tumorcentery], 2);
        if(dist<tumorradius) 
            exclude=1;
        end
    end
    if(exclude==0)
        xpos=XFa.LL(:,1);
        ypos=abs(XFa.LL(:,2)-pixh-1);
        fdist=norm([xpos(segnum) ypos(segnum)]-[xpos(1) ypos(1)],2);
        tempfnum=floor(fdist/minfsegleng);
        if (tempfnum > 0)
            [temppos]= Quantize_Fiber_Segments(xpos, ypos, minfsegleng);
            angpos=[angpos;temppos];
        end
    end
end

% Plot the quantized and line-fitted fiber segments with tumor boundary
figure();
axis xy;
axis equal;
axis off;
hold on;
Display_Boundary(1, 512, 1, 512);    
circleang=[0:0.01:2*pi];
circlex=[];
circley=[];
for i=1:length(circleang)
    boundaryx=cos(circleang(i))*tumorradius+tumorcenterx;
    boundaryy=sin(circleang(i))*tumorradius+tumorcentery;
    if (boundaryx>0 && boundaryx<512 && boundaryy>0 && boundaryy<512)
        circlex=[circlex; boundaryx];
        circley=[circley; boundaryy];
    end
end
size1=length(circlex(:,1));
for bline=2:length(circlex(:,1))
    tempxdelta=abs(circlex(bline-1)-circlex(bline));
    tempydelta=abs(circley(bline-1)-circley(bline));
    if(tempxdelta<10 && tempydelta<10)
        plot(circlex(bline-1:bline, 1), circley(bline-1:bline, 1), '-', 'color', [0.5 0.5 0.5]);
    end
end
for i=1:length(angpos(:,1))
    plot([angpos(i,1) angpos(i,3)], [angpos(i,2) angpos(i,4)], '-k');
end
title('Quantized fiber segments','FontName','Arial','FontWeight','Bold','FontSize',12);
saveas(gcf,'Step1_Fig3c.tif');
close(gcf);

% Step3: Calculate alignment vectors         
% To calculate local alignment vector, we calculated the alignment vector for
% every 5 pixel interval in X and Y. 
% The total number of local alignment vector for each Z stack image is
% 103 x 103 = 10609 local points.
xbins=[1:5:512];
ybins=[1:5:512];
local_point_num=length(xbins);

% fiber counts inside of a local circular window bin
fnum=zeros(local_point_num,local_point_num); 

% The radius of each local circle
Sradius=25; 

% |alignment vector| = alignment index, where 1 is for the perfectly
% aligned fibers in the local circle, and 0 is for equally distributed
% fibers in the circular statistics (e.g. 90, 180, 270, 360 degree angle
% for four fibers)
AI=zeros(local_point_num,local_point_num); 

% Alignment vetor angle, from 0 to 180 degree
AV_angle=zeros(local_point_num,local_point_num);

for x=1:local_point_num
    for y=1:local_point_num
        [fcount, temp_AI, temp_ang] = Calculate_Alignment_Vector(angpos, xbins(x), ybins(y), Sradius);
        fnum(x,y)=fcount;
        AI(x,y)=temp_AI; 
        AV_angle(x,y)=temp_ang;
    end
end

% Step4: Plot local alignment vectors
% Plot the vector if the fiber number of a local circle > 10, and the
% alignment index of a local circle > 0
minAI=0;  
fmin=10; 

% Local alignment vectors are color-coded (jet) and varying length (from 0 pixel 
% to Sradius pixel)
cval=colormap('jet');
intmrv=(1-minAI)/64;

figure();
axis xy;
axis equal;
axis off;
hold on;
circleang=[0:0.01:2*pi];
circlex=[];
circley=[];
for i=1:length(circleang)
    boundaryx=cos(circleang(i))*tumorradius+tumorcenterx;
    boundaryy=sin(circleang(i))*tumorradius+tumorcentery;
    if (boundaryx>0 && boundaryx<512 && boundaryy>0 && boundaryy<512)
        circlex=[circlex; boundaryx];
        circley=[circley; boundaryy];
    end
end
size1=length(circlex(:,1));
for bline=2:length(circlex(:,1))
    tempxdelta=abs(circlex(bline-1)-circlex(bline));
    tempydelta=abs(circley(bline-1)-circley(bline));
    if(tempxdelta<10 && tempydelta<10)
        plot(circlex(bline-1:bline, 1), circley(bline-1:bline, 1), '-', 'color', [0.5 0.5 0.5]);
    end
end

AI_distribution=[];
for x=1:local_point_num
    for y=1:local_point_num
        if((xbins(x)-Sradius)>=0 && (xbins(x)+Sradius)<=512 && (ybins(y)-Sradius)>=0 && (ybins(y)+Sradius)<=512)
            testdistance=norm([tumorcenterx, tumorcentery]-[xbins(x), ybins(y)],2);
            if(testdistance>=(Sradius+tumorradius))
                if (fnum(x,y)>fmin && AI(x,y)>minAI)
                    AI_distribution=[AI_distribution; AI(x,y)];
                    tempc=round((AI(x,y)-minAI)/intmrv);
                    if (tempc==0) 
                       tempc=1;
                    end
                    R2=Sradius*AI(x,y);
                    x2=R2*cos(AV_angle(x,y)); 
                    y2=R2*sin(AV_angle(x,y));
                    plot([xbins(x)-x2/2, xbins(x)+x2/2], [ybins(y)-y2/2, ybins(y)+y2/2], '-', 'color', cval(tempc,:));
                end
            end
        end
    end
end
    
caxis([minAI 1]);
colorbar;

% Show the local circle size
SampleLocalCircle_x=60;
SampleLocalCircle_y=430;
circleang=[0:0.01:2*pi];
circlex=[];
circley=[];
for i=1:length(circleang)
    boundaryx=cos(circleang(i))*Sradius+SampleLocalCircle_x;
    boundaryy=sin(circleang(i))*Sradius+SampleLocalCircle_y;
    if (boundaryx>0 && boundaryx<512 && boundaryy>0 && boundaryy<512)
        circlex=[circlex; boundaryx];
        circley=[circley; boundaryy];
    end
end
size1=length(circlex(:,1));
for bline=2:length(circlex(:,1))
    tempxdelta=abs(circlex(bline-1)-circlex(bline));
    tempydelta=abs(circley(bline-1)-circley(bline));
    if(tempxdelta<10 && tempydelta<10)
        plot(circlex(bline-1:bline, 1), circley(bline-1:bline, 1), '-k');
    end
end
plot(SampleLocalCircle_x, SampleLocalCircle_y, '.k', 'MarkerSize', 5);
title('Local alignment vectors, R=25, F>10','FontName','Arial','FontWeight','Bold','FontSize',12);
Display_Boundary(1, 512, 1, 512);        
saveas(gcf,'Step4_Fig3h.tif');
close(gcf);

% Step5: Plot the distribution of local alignment vector
interval=0:0.05:1;
[ne, cen]=hist(AI_distribution,interval);

% The local AI distribution was normalized with one single z stack for the
% example case. For the main figure, the distribution was normalized with
% all z stacks.
nen=ne./sum(ne);

linevalue=2;
figure();
plot(cen, nen, '-m', 'linewidth', linevalue);
title('H1299 shLKB1, 21hr, Z=1','FontName','Arial','FontWeight','Bold','FontSize',12);
xlim([0 1]);
ylim([0 0.13]);
set(gca,'XTick',[0 0.5 1]);
set(gca,'YTick',[0 0.12]);
set(gca, 'FontName','Arial','FontSize',12);
ylabel('Local alignment','FontName','Arial','FontWeight','Bold','FontSize',12);
xlabel('Alignment index','FontName','Arial','FontWeight','Bold','FontSize',12);
saveas(gcf,'Step5.tif');
close(gcf);
         
         