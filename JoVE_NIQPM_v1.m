% JoVE_NIQPM_v1.m % Overview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Updated: May 1, 2013
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
%   1.) dstn.m     [discrete sine transform in n-dimensions]
%   2.) isdstn.m    [discrete inverse sine transform in n-dimensions] 
%   3.) dctn.m     [discrete cosine transform in n-dimensions]
%   4.) idctn.m     [discrete inverse cosine transform in n-dimensions]
%   3.) phase_algorithm.m  [implementation of FFT based phase algorithm of Frank.]
%
% INPUT:
%   1.) Your directory of .tiff images from low NA illumination
%   through-focus bright field imagery.
%
% BASIC MATLAB PRELIMINARIES:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 0. Clean the slate, update path to dependencies 
close all
clc
clear

% Examples:
% on a mac
% dependencies_directory = '/Users/your/directory/containing/the/matlab/programs/above';

% on a PC
% dependencies_directory = 'c:\your\directory\containing\the\matlab\programs\above';

% Copy and paste the directory of your matlab programs: hilbert_transform_dic.m &
% sobel_edge_detect.m between the single quotes below:
dependencies_directory = '/Users/kevinGphillips/Documents/work/papers/JOVE/resubmission';

% Copy and paste the directory of your bright field images between the single quotes below:
brightfield_directory = '/Volumes/LITTLE MAN/kevin work/papers/Jove/tiff/na_0p1/sw620/bf';

% Copy and paste the directory of your DIC images between the single quotes below:
dic_directory = '/Volumes/LITTLE MAN/kevin work/papers/Jove/tiff/na_0p9/sw620/dic';

% update path
path(path,dependencies_directory)

% go get corresponding DIC imagery information (pray for co-registration)
cd(dic_directory)
dic_tiff_list = dir('*.tif');
[nDIC foo] = size(dic_tiff_list);

cd(brightfield_directory)

%% Section 1. Run ONCE, CD to directory of interest and import the image cube

% on a mac
% directory = '/Volumes/your/directory/containing/tiff/files';

% on a PC
% directory = 'c:\your\directory\containing\tiff\files';

cd(brightfield_directory)

% creat list of .tif files
tifList = dir('*.tif');

% Preallocate the image cube "bfStack"
[N_pxl_z, foo] = size(tifList);
[N_pxl_x N_pxl_y N_chan] = size(imread(tifList(1).name));
bfStack = zeros(N_pxl_x, N_pxl_y, N_pxl_z);

% Assign en face images to each layer of the cube
for iz = 1:N_pxl_z
    % Read in the BF intensity image
    Inew = double(imread(tifList(iz).name));
    % take green channel of .tif
    bfStack(:,:,iz) = Inew(:,:,2); 
end

%% Section 2. Choose a cell of interest, update dx and dz

%%%% Assign parameters %%%%%%%%%%%%%%%%%%%%
prompt = {'Enter focal plane number', 'Enter lateral resolution [micron per pixel]', 'Enter axial resolution [micron per pixel]' 'Enter the region of interest size [pixel]'};
dlg_title = sprintf('Define NIQPM Parameters');
num_lines = 1;
% Default values of:
%   focus no.,   dx,    dz,     roi width
def = {'155',    '0.1',    '0.1',     '200'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

% assign user input to program variables
focal_ImageNum   = str2double(answer{1,1});
dx = str2double(answer{2,1});
dz = str2double(answer{3,1});
roiWidth = str2double(answer{4,1});

% Plot the full field
figure(200); clf
curr_xy = reshape(bfStack(:,:,focal_ImageNum),N_pxl_x,N_pxl_y);
imagesc(curr_xy)
colormap(gray)
h = imrect(gca, [10 10 roiWidth roiWidth]);
roi = wait(h);
close(figure(200))

% Crop image to ROI
inFocus   = imcrop(curr_xy, roi);
[Nx Ny] = size(inFocus);

% Plot the full field cropped to the ROI
figure(100)
imagesc((1:Ny).*dx, (1:Nx).*dx, inFocus); axis square; colormap(gray)
xlabel('x [\mum]')
ylabel('y [\mum]')
title('En face BF profile')
axis equal image

%% Section 3. Construct the bright field stack cropped to the image Perform the axial averaging

% WARNING: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust "aveOn", "planesToAve", "planesToOmit" as a last resort.
%
% The NIQPM method has been validated only with their default values of:
%   aveOn = 0; planesToAve = 1; planesToOmit = 1;
%
% Changing these parameters will require a validation experiment on spheres 
% to ensure these parameters do not provide artificial information. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nx Ny] = size(inFocus);
clear axial_ave_stack
axial_ave_stack = zeros(Nx, Ny, N_pxl_z);

% Perform either axial averaging or omit planes from the cube.
aveOn = 0;

% Adjust parameters to eliminate noise.
planesToAve = 1;
planesToOmit = 1;

if aveOn
    for iz = 1 : planesToAve: N_pxl_z-planesToAve
        currPlane = sum(bfStack( :, :, iz: iz+planesToAve ) ,3);
        cropCurr = imcrop(currPlane, roi);
        axial_ave_stack(:,:,iz) = cropCurr;
    end
else
    for iz = 1 : planesToOmit: N_pxl_z-planesToOmit
        currPlane = reshape(bfStack( :, :, iz ), N_pxl_x, N_pxl_y);
        cropCurr = imcrop(currPlane, roi);
        axial_ave_stack(:,:,iz) = cropCurr;
    end
end

[Nx Ny Nz] = size(axial_ave_stack); 

%% Section 4. Perform the phase algorithm and plot summary images

% Calculation parameters

%%%% Assign phase algorithm parameters %%%%%%%%%%%%%%%%%%%%
prompt = {'Enter bright field plane no.', 'Enter DIC focal plane no.', 'Enter Fourier filter (0-3) '};
dlg_title = sprintf('Define Phase Algorithm Parameters');
num_lines = 1;
% Default values of:
%   focus no.,      DIC focal plane,     Fourier cutoff (0-3)
def = {'100',          '150',             '0'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

% Focal image number to center phase calculation around.
bf_focus = str2double(answer{1,1});
% DIC optimal focal plane
dic_focus = str2double(answer{2,1});
% Fourier Filtering if desired
cutOff = str2double(answer{3,1});


% Warning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "dz" pixel size and basis for PDE solver.
% Adjust "deFocus_stepSize1" as a last resort.
%
% The NIQPM method presented has been validated only with the default values of:
% deFocus_stepSize1 = 11, sineON = 1;
%
% Advanced users can change to sineON = 0 to employ a cosine basis - this
% effectively changes the boundary condition for the T.I.E.
%
% Changing these parameters will require a validation experiment on spheres 
% to ensure the accuracy of the phase reconstruction.
deFocus_stepSize1 = 11;
sineON = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute the derivative in pixel units! No dz in microns, please.
[tot_phase, inFocus] = phase_algorithm(axial_ave_stack, planesToOmit, cutOff, bf_focus, deFocus_stepSize1, sineON);

% Summary image
figure(400); clf

% plot the bright field image
subplot(1,4,1)
imagesc((1:Ny).*dx, (1:Nx).*dx, (inFocus) )
axis equal image
colormap(gray)
title('Bright field intensity')
xlabel('[\mum]')
ylabel('[\mum]')

% plot the reconstructed phase map.
subplot(1,4,2)
if sineON
    imagesc( (1:Ny).*dx, (1:Nx).*dx,  (tot_phase), [0 2] )
else
    imagesc( (1:Ny).*dx, (1:Nx).*dx,  (tot_phase), [-pi/2 pi/2] )
end

% colorbar
title('phase map')
axis equal image
xlabel('[\mum]')
ylabel('[\mum]')

% Factor converting phase to mass density.
convert_to_mass = .54/(2.*pi)/.18;

% define mass density from phase
massDensityMap = ((tot_phase)).*convert_to_mass; 

% Compute "pseudo"-DIC image from phase map
% [Nx, Ny] = size(tot_phase);
pseudo_dic = zeros(Nx, Ny);
stepdiag = 1;
for ix = 1:Nx
    for iy = 1:Ny
        if ix+stepdiag < Nx && iy+stepdiag < Ny 
        pseudo_dic(ix , iy) = -(tot_phase(ix,iy)) + (tot_phase(ix+stepdiag,iy+stepdiag));
        else
        pseudo_dic(ix , iy) = tot_phase(ix,iy);   
        end
    end
end

% plot pseudo DIC
figure(400)
subplot(1,4,3) 
imagesc(  (1:Ny).*dx, (1:Nx).*dx, pseudo_dic )
axis equal image
title('pseudo DIC')
colormap(gray)  
xlabel('[\mum]')
ylabel('[\mum]')
 
% go get corresponding DIC imagery (pray for co-registration)
cd(dic_directory)
dic_tiff_list = dir('*.tif');
[nDIC foo] = size(dic_tiff_list);
Inew = double(imread(dic_tiff_list(dic_focus).name));
DIC = Inew(:,:,2); %
% crop to ROI
inFocus_dic   = imcrop(DIC, roi);

% plot the NA = 0.9 DIC image
figure(400)
subplot(1,4,4)
imagesc((1:Ny).*dx, (1:Nx).*dx, inFocus_dic ) 
title('Orig. DIC')
axis equal image
colormap(gray)
xlabel('[\mum]')
ylabel('[\mum]')

%% Section 5a. Select outline of cell using edge detection: determine mass and density distribution
clear Segout mask_xy

threshold = .8; %<<<----------- Adjust as needed to get proper edge detection

% oulining
[mask_xy,Segout] = sobel_edge_detect(massDensityMap, threshold);
BWoutline = bwperim(mask_xy);
Segout2 = massDensityMap;
Segout2(BWoutline) = 10;

% determine mass
density_in_mask = massDensityMap(mask_xy);

if isempty(density_in_mask)
    f = warndlg('No mask determined. Update threshold or choose a different ROI. ', 'Error found');
    break
else    
final_density = density_in_mask - min(density_in_mask(:));
cell_mass = sum(final_density).*dx^2;
end

figure(500);clf
subplot(1,2,1)
imagesc( (1:Nx).*dx, (1:Nx).*dx, Segout2-mask_xy.*min(density_in_mask(:)), [0 1])
colormap(jet)
axis square
colorbar
title(sprintf('Cell Mass = %f [pg] ',cell_mass))
xlabel('[\mum]')
ylabel('[\mum]')

% bin into histogram
subplot(1,2,2)
hist(final_density, 50);
axis square
xlabel('Density [pg/\mum^2]')
ylabel('Counts')

%% Section 5b. Select outline of cell by hand: determine mass and density distribution

figure(501)
imagesc( (1:Nx).*dx, (1:Nx).*dx, massDensityMap, [0 1])
axis square
xlabel('[\mum]')
ylabel('[\mum]')
h = imfreehand;
wait(h);
mask_xy = createMask(h);
BWoutline = bwperim(mask_xy);
Segout2 = massDensityMap;
Segout2(BWoutline) = 10;
close(figure(501))

% determine mass
density_in_mask = massDensityMap(mask_xy);
final_density = density_in_mask - min(density_in_mask(:));
cell_mass = sum(final_density).*dx^2;

figure(500);clf
subplot(1,2,1)
imagesc( (1:Nx).*dx, (1:Nx).*dx, Segout2, [0 1])
colormap(jet)
colorbar
axis square
title(sprintf('Cell Mass = %f [pg] ',cell_mass))
xlabel('[\mum]')
ylabel('[\mum]')

% bin into histogram
subplot(1,2,2)
hist(final_density, 50);
axis square
xlabel('Density [pg/\mum^2]')
ylabel('Counts')