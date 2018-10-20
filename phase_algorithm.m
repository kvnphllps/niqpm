function [tot_phase, inFocus] = phase_algorithm(stack, dz, cutOff, focal_ImageNum, deFocus_stepSize1,   sineON)

[Nx Ny Nz] = size(stack);


%% Section 1. Define Eigenvalues of Laplace operator & cutoff freq. for image geometry

eigLap_x = zeros(Nx,Ny);
eigLap_y = zeros(Nx,Ny);

% eigLap
for ix = 1:Nx
    for iy = 1:Ny
        kx = 2*pi*ix/(Nx);
        ky = 2*pi*iy/(Ny);
        eigLap_x(ix,iy) =  kx./(kx^2 + ky^2);
        eigLap_y(ix,iy) =  ky./(kx^2 + ky^2);
    end
end
% cutOff = 1;
if cutOff
    eigLap_x(1:cutOff,1:cutOff) = 0;
    eigLap_y(1:cutOff,1:cutOff) = 0;
end

%% Section 3. Perform TBFI algorithm on axial_ave_stack

% Plane in z-stack to compute index map
% focal_ImageNum = ;

iz = focal_ImageNum;

% dz's used in first and second derivative, respectively, in units of micron
dz1 = dz*deFocus_stepSize1; % first derivative approx.
[Nx Ny Nz] = size(stack);
% Compute second order z-derivative


% compute first z-derivatives

% One unit above focal plane
currData1 = reshape(stack(:,:,iz+deFocus_stepSize1),Nx,Ny);
currData1 = currData1./(mean(currData1(:)));

currData2 = reshape(stack(:,:,iz-deFocus_stepSize1),Nx,Ny);
currData2 = currData2./(mean(currData2(:)));

% focal plane
inFocus = reshape(stack(:,:,iz),Nx,Ny);
inFocus = inFocus./(mean(inFocus(:)));

%firstDeriv = -(currData1-currData2)./dz1/2;
d1 = (currData2-inFocus)./dz1;

% First order
% d1 = (firstDeriv(1:roiWidth,1:roiWidth)) ;


% smooth in focus BF image

waveNum = 2*pi/.540;
% TBFI reconstruction algorithm: dstn vs. dctn
% sineON = 1;
size(eigLap_y)
size(d1)

if sineON
    % Discrete Sine Based Calc.
    
     phi_x = idstn(eigLap_x.*dstn((1./(inFocus)).*idstn(eigLap_x.*dstn(d1))));
     phi_y = idstn(eigLap_y.*dstn((1./(inFocus)).*idstn(eigLap_y.*dstn(d1))));
    
    tot_phase = ( phi_x+phi_y).*waveNum;
else
    % Discrete Cosine Based Calc.
    phi_x = idctn(eigLap_x.*dctn((1./(inFocus)).*idctn(eigLap_x.*dctn(d1))));
    phi_y = idctn(eigLap_y.*dctn((1./(inFocus)).*idctn(eigLap_y.*dctn(d1))));
    tot_phase = (phi_x + phi_y).*waveNum;
end

