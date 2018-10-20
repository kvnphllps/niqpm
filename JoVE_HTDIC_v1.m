% JoVE_HTDIC_v1.m % Overview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated: April 30, 2013
%
% Copyright 2013. 
% Kevin G. Phillips, Sandra M. Baker-Gronberg, Owen J.T. McCarty
% Department of Biomedical Engineering
% Oregon Health & Science University
%
% **Cite the corresponding JoVE article when using this program in
% published studies.**
%
% DEPENDENCIES:
%   1.) hilbert_transform_dic.m [performs the Hilbert transform of an image]
%
%   2.) sobel_edge_detect.m  [performs Sobel-based edge detection]
%
% INPUT:
%   1.) Your directory of .tiff images from high NA illumination through-focus
%   DIC imagery.
%
%   2.) The number of microns per pixel in your images.
%
% BASIC MATLAB PRELIMINARIES:
%
%   1.) Turn on "code folding" to run this program a "cell at a time."
%     
%   2.) Run an individual cell by clicking anywhere inside it and typing:
%        pc: control+enter
%        mac: command+enter   
%
%   3.) Comment a line of code by placing a "%" at the beginning of the
%   line - the text will then turn green. Un-comment the code by deleting
%   the "%" at the beginning of the line.
%
%   4.) Path update: you will need to install the dependency programs "hilbert_transform_dic.m" and "sobel_edge_detection.m" 
%     in the same directory as one another. Then in Section 0 you will
%     update the path to include this directory location.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 0. Clean the slate, update path to dependencies
close all
clc
clear

% Copy and paste the directory of your matlab programs: hilbert_transform_dic.m &
% sobel_edge_detect.m between the single quotes below:
dependencies_directory = '/Users/kevinGphillips/Documents/work/papers/JOVE/resubmission';

% on a mac
% dependencies_directory = '/Users/your/directory/containing/the/matlab/programs/above';

% on a PC
% dependencies_directory = 'c:\your\directory\containing\the\matlab\programs\above';

path(path,dependencies_directory)

%% Section 1.  Run ONCE, CD to the directory of interest and import image cube

% Copy and paste the directory of your .tif images between the single
% quotes below
images_directory = '/Volumes/LITTLE MAN/kevin work/papers/Jove/tiff/na_0p9/sw620/dic';

% on a mac
% images_directory = '/User/your/directory/containing/tiff/files';

% on a PC
% images_directory = 'c:\your\directory\containing\tiff\files';

cd(images_directory)

% the goal is to run this section once only.
% create DIC stack for use throughout the program

% creat list of .tif files
tifList = dir('*.tif');

% Preallocate the image cube "dicStack"
[N_pxl_z, foo] = size(tifList);
[N_pxl_x N_pxl_y N_chan] = size(imread(tifList(1).name));
dicStack = zeros(N_pxl_x, N_pxl_y, N_pxl_z);

for iz = 1:N_pxl_z
    % Read in the DIC intensity image
    Inew = double(imread(tifList(iz).name));
    % take red channel of .tif
    dicStack(:,:,iz) = Inew(:,:,1);     
end

%% Section 2. Choose a cell of interest & the corresponding rotation angle: update dx and dz

% clear previous parameters

clear I3 Inew NX NY Nx NxROI Ny NyROI Nz Segout clean co currPlane currPlaneHil
clear dummy dx dz filter final_clean focal_ImageNum hilbertImage hilbertStack
clear hilbertTrans inFocus inFocus1 iz k list mask_xy plt_area plt_volume poop roi1
clear roi2 roiWidth roiWidth2 rot_inFocus stepX stepY threshold vol_auto xz_hil

% assign parameters
prompt = {' Enter focal plane number', 'Enter lateral resolution [micron per pixel]', 'Enter axial resolution [micron per pixel]', 'Enter rotation angle', 'Enter the region of interest size [pixel]'};
dlg_title = sprintf('Define HTDIC Parameters');
num_lines = 1;
% Default values of:
%    focus no.,   dx,  dz,    rot. angle  roi width
def = {'135',     '0.1',  '0.1',    '-135',     '400'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

% assign user input to program variables
focal_ImageNum   = str2double(answer{1,1});
dx = str2double(answer{2,1});
dz = str2double(answer{3,1});
hilbertAngle = str2double(answer{4,1});
roiWidth = str2double(answer{5,1});

[Nx Ny Nz] = size(dicStack);

% Select in-focus image
currPlane = reshape(dicStack(:,:,focal_ImageNum),Nx,Ny);

% Select the cell of interest
figure(200);
imagesc(currPlane)

colormap(gray)
roiWidth = 250;
h = imrect(gca, [15 15 roiWidth roiWidth]);
roi1 = wait(h);
close(figure(200))
%
cropCurr   = imcrop(currPlane, roi1);

% Rotate the cell- align bias with +x direction
rotCurr = imrotate(cropCurr,hilbertAngle);
figure(200)
imagesc(rotCurr)
colormap(gray)
axis equal image

% Determine the crop for the rotated cell
roiWidth2 = 150;
h = imrect(gca, [10 10 roiWidth2 roiWidth2]);
roi2 = wait(h);
close(figure(200))
clear inFocus

% Apply the crop
crop_rot_curr   = imcrop(rotCurr, roi2);

% Circle the cell of interest
figure(200)
imagesc(crop_rot_curr)
colormap(gray)
axis equal image

%% Section 3. Define xy mask 

figure(300)

imagesc(crop_rot_curr)
colormap(gray)
axis equal image

% use to enclose the cell in a rectangle
%h = imrect;

% use to circle the cell by freehand
h = imfreehand;

wait(h);
mask_xy = createMask(h);
close(figure(300))

%% Section 4. Perform the Hilbert transform and stack

% set maskON = 1 to apply the mask, set = 0 for no mask. 
maskON = 1;

% set plotON = 1 to plot each xy image of the z-stack, set = 0 for no
% plotting (faster).
plotON = 0;

% Preallocation of the Hilbert-transformed DIC stack
[N_pxl_x N_pxl_y N_pxl_z] = size(dicStack);
[Nx Ny] = size(crop_rot_curr);
clear hilbertStack dicStack2
hilbertStack = zeros(Nx, Ny, N_pxl_z);
dicStack2 = zeros(Nx, Ny, N_pxl_z);

for iz = 1:N_pxl_z
    
    % Read in the DIC intensity image
    I3 = reshape(dicStack(:,:,iz),N_pxl_x,N_pxl_y);
    
    % First crop
    inFocus1   = imcrop(I3, roi1);
    
    % Now rotate
    rot_inFocus = imrotate(inFocus1,hilbertAngle);
    
    % Apply final crop
    inFocus   = (imcrop(rot_inFocus, roi2));
    
    % Perform Hilbert Transform
    [hilbertTrans] = hilbert_transform_dic(inFocus);
    
    % Apply mask if desired
    if maskON
    hilbertTrans = hilbertTrans.*(mask_xy);
    inFocus = inFocus.*(mask_xy);
    end
    
    % Plot each xy hilbert transform image if desired
    if plotON
    figure(400); 
    imagesc(hilbertTrans); 
    axis equal image
    colormap(gray)
    end
    
    hilbertStack(:,:,iz) = hilbertTrans;
    dicStack2(:,:,iz) = inFocus;        
end

%% Section 5. Compare cross sectional images of DIC, HT-DIC, F-HT-DIC & segmentation

threshold = .5; %<-------- OPTIMIZE THIS NUMBER FOR EDGE DETECTION
loopON = 0;
[Nx Ny Nz] = size(hilbertStack);

if loopON
    lower_x_limit = 1;
    upper_x_limit = Nx;
else
    lower_x_limit = floor(Nx/2)+5;
    upper_x_limit = floor(Nx/2)+5;
end

for ix = lower_x_limit: upper_x_limit
    
    curr_dic_plane = (reshape(dicStack2(:, ix, :), Nx, Nz)');
    
    curr_hilbert_plane = (abs((reshape(hilbertStack(:, ix, :), Nx, Nz)')  ) );
    
    % filter the Hilbert Stack: high pass filter.
    filter = ones(Nx,Nz)';
    %diag(filter) = 0;
    filter(:,1) = 0;
    filter(1,:) = 0;
    %       filter(1,:) = 0;
    %       diag(filter) = 0;
    curr_ft_hilbert_plane = abs(ifft2(fft2(curr_hilbert_plane).*filter).^2);
  
    figure(500); clf
    subplot(2,3,1)
    imagesc((1:Nx).*dx, (1:Nz).*dz, curr_dic_plane)
    axis equal image
    colormap(gray)
    title('DIC')
    
    subplot(2,3,2)
    imagesc((1:Nx).*dx, (1:Nz).*dz, curr_hilbert_plane)
    axis equal image
    colormap(gray)
    title('HT-DIC')
    
    subplot(2,3,3)
    imagesc((1:Nx).*dx, (1:Nz).*dz, curr_ft_hilbert_plane, [100 700])
    axis equal image
    % colorbar
    colormap(gray)
    title('F-HT-DIC')
    
    % Outline the DIC based on edge detection in DIC
    [mask_xy1,Segout_foo] = sobel_edge_detect(curr_dic_plane,threshold);
    BWoutline0 = bwperim(mask_xy1);
    Segout0 = curr_dic_plane;
    Segout0(BWoutline0) = 0;
    
    % Outline the DIC based on edge detection in HT-DIC
    [mask_xy1,Segout] = sobel_edge_detect(curr_hilbert_plane,threshold);
    BWoutline1 = bwperim(mask_xy1);
    Segout1 = curr_dic_plane;
    Segout1(BWoutline1) = 0;
    
    % Outline the DIC based on the F-HT-DIC
    [mask_xy2,Segout2] = sobel_edge_detect(curr_ft_hilbert_plane,threshold);
    BWoutline2 = bwperim(mask_xy2,4);
    Segout3 = curr_dic_plane;
    Segout3(BWoutline2) = 0;
    
    figure(500)
    subplot(2,3,4)
    imagesc((1:Nx).*dx, (1:Nz).*dz, Segout0)
    axis equal image
    title('Segment with DIC input')
    
    subplot(2,3,5)
    imagesc((1:Nx).*dx, (1:Nz).*dz, Segout1)
    axis equal image
    title('Segment with HT-DIC input')
    
    subplot(2,3,6)
    imagesc((1:Nx).*dx, (1:Nz).*dz, Segout3)
    axis equal image
    colormap(gray)
    title('Segment with F-HT-DIC input')
    
    pause(0.1)
end

%% Section 6. Volume calculation based on DIC input

[Nx Ny Nz] = size(dicStack2);

clear vol_auto;
%  threshold = .8; % uncomment to adjust

dummy = 1;
dz = 0.1;
dx = 0.1;

for iy = 1 : Ny
    
    % define en face DIC image
    en_face_dic = reshape(dicStack2(:,:,focal_ImageNum),Nx, Ny);
    en_face_dic(:,iy) = 0;
    
    % plot en face DIC image
    figure(600)
    subplot(1,2,1)
    imagesc((1:Ny).*dx, (1:Nx).*dx, en_face_dic)
    colormap(gray)
    axis equal image
    
    % define the cross sectional DIC plane
    curr_dic_plane = (reshape(dicStack2(:, iy, :), Nx, Nz)');
    
    
   
    % Outline the DIC based on edge detection in DIC cross sectional image
    % define cross sectional DIC
    curr_dic_plane = (reshape(dicStack2(:, iy, :), Nx, Nz)');
    % threshold = .58;
    [mask_xy1,Segout] = sobel_edge_detect(curr_dic_plane,threshold); % input sagittal image type
    BWoutline1 = bwperim(mask_xy1);
    %     Segout1 = curr_dic_plane;
    %     Segout1 = curr_dic_plane;
    Segout1 = curr_dic_plane;
    
    % make the outline more visible
    if iy <= floor(Ny/2)
        Segout1(BWoutline1) = 100;
    else
        Segout1(BWoutline1) = 0;
    end
    
    vol_auto(dummy) = sum(mask_xy1(:)).*dz*dx*dx;
    dummy = dummy + 1;
    
    % plot
    subplot(1,2,2)
    imagesc( (1:Nx).*dx, (1:Nz).*dz, Segout1, [50 150])
    axis equal image
    colormap(gray)
    title('Segment with DIC input')
    
end
cell_volume = sum(vol_auto);
subplot(1,2,1)
title(sprintf('Volume = %f [fL]',cell_volume))

%% Section 7. Volume calcuation based on Hilbert Transform input

[Nx Ny Nz] = size(hilbertStack);

clear vol_auto;
% threshold = .58; % uncomment to adjust


dummy = 1;
dz = 0.1;
dx = 0.1;

for iy = 1 : Ny
    
    % define en face DIC image
    en_face_dic = reshape(dicStack2(:,:,focal_ImageNum),Nx, Ny);
    en_face_dic(:,iy) = 0;
    
    % plot en face DIC image
    figure(700)
    subplot(1,2,1)
    imagesc((1:Ny).*dx, (1:Nx).*dx, en_face_dic)
    colormap(gray)
    axis equal image
    
    % define the cross sectional DIC plane
    curr_dic_plane = (reshape(dicStack2(:, iy, :), Nx, Nz)');
    
    % define the current Hilbert plane
    curr_hilbert_plane = (abs((reshape(hilbertStack(:, iy, :), Nx, Nz)')  ) );
    
    % Outline the DIC based on edge detection in HT-DIC
    % define cross sectional DIC
    curr_dic_plane = (reshape(dicStack2(:, iy, :), Nx, Nz)');
    % threshold = .58;
    [mask_xy1,Segout] = sobel_edge_detect(curr_hilbert_plane,threshold); % input sagittal image type
    BWoutline1 = bwperim(mask_xy1);
    Segout1 = curr_dic_plane;
    Segout1(BWoutline1) = 1;
    
    vol_auto(dummy) = sum(mask_xy1(:)).*dz*dx*dx;
    dummy = dummy + 1;
    
    % plot
    subplot(1,2,2)
    imagesc( (1:Nx).*dx, (1:Nz).*dz, Segout1, [50 100])
    axis equal image
    colormap(gray)
    title('Segment with HT-DIC input')
%   pause()  
end
cell_volume = sum(vol_auto);
subplot(1,2,1)
title(sprintf('Volume = %f [fL]',cell_volume))

%% Section 8. Volume calcuation based on F-HT-DIC input

[Nx Ny Nz] = size(hilbertStack);

clear vol_auto;
%  threshold = .05; % uncomment to adjust


dummy = 1;
dz = 0.1;
dx = 0.1;

for iy = 1 : Ny
    
    % define en face DIC image
    en_face_dic = reshape(dicStack2(:,:,focal_ImageNum),Nx, Ny);
    en_face_dic(:,iy) = 0;
    
    % plot en face DIC image
    figure(800)
    subplot(1,2,1)
    imagesc((1:Ny).*dx, (1:Nx).*dx, en_face_dic)
    colormap(gray)
    axis equal image
    
    % define the cross sectional DIC plane
    curr_dic_plane = (reshape(dicStack2(:, iy, :), Nx, Nz)');
    
    % define the current Hilbert plane
    curr_hilbert_plane = (abs(reshape(hilbertStack(:, iy, :), Nx, Nz)'));
    
    
    % define the filtered Hilbert plane
    curr_ft_hilbert_plane = abs(ifft2(fft2(curr_hilbert_plane).*filter).^2);
        
    % Outline the DIC based on edge detection in F-HT-DIC
    % define cross sectional DIC
    curr_dic_plane = (reshape(dicStack2(:, iy, :), Nx, Nz)');
    
    [mask_xy1,Segout] = sobel_edge_detect(curr_ft_hilbert_plane,threshold); % input sagittal image type
    BWoutline1 = bwperim(mask_xy1);
    Segout1 = curr_dic_plane;
    Segout1(BWoutline1) = 300;
    
    vol_auto(dummy) = sum(mask_xy1(:)).*dz*dx*dx;
    dummy = dummy + 1;
    
    % plot
    subplot(1,2,2)
    imagesc( (1:Nx).*dx, (1:Nz).*dz, Segout1, [50 200])
    axis equal image
    colormap(gray)
    title('Segment with HT-DIC input')
    
end
cell_volume = sum(vol_auto);
subplot(1,2,1)
title(sprintf('Volume = %f [fL]',cell_volume))

