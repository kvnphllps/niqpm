function [hilbertTrans] = hilbert_transform_dic(inFocus)

% Hilbert Mask
[Nx Ny] = size(inFocus);
hilbertMask = ones(Nx,Ny);
hilbertMask(:,1:floor(Nx/2)) = -1;
hilbertMask(:,floor(Nx/2)) = 0;

% Hilbert Transform
real_inFocus = real(fft2( inFocus ));
imag_inFocus = imag(fft2( inFocus ));
stuff =  (imag_inFocus.*hilbertMask) + 1i.*(real_inFocus.*hilbertMask);
hilbertTrans =  (fliplr(flipud(real((ifft2(stuff))) ))) ;