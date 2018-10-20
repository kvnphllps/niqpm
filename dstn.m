function [y,w] = dstn(y,DIM,w)

%DSTN N-D discrete sine transform.
%   Y = DSTN(X) returns the discrete sine transform of X. The array Y is
%   the same size as X and contains the discrete sine transform
%   coefficients. This transform can be inverted using IDSTN.
%
%   DSTN(X,DIM) applies the DSTN operation across the dimension DIM.
%
%   Class Support
%   -------------
%   Input array can be numeric or logical. The returned array is of class
%   double.
%
%   Reference
%   ---------
%   Adapted from Narasimha M. et al, On the computation of the discrete
%   cosine transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dstn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DST matrix
%   to zero, then reconstruct the image using the inverse DST.
%
%       J(abs(J)<10) = 0;
%       K = idstn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also IDSTN, DCTN
%
%   -- Damien Garcia -- 2009/11
%   Visit my <a
%   href="matlab:web('http://www.biomecardio.com/matlab/dctn.html')">website</a> for discrete cosine transform

% ----------
%   [Y,W] = DSTN(X,DIM,W) uses and returns the weights which are used by
%   the program. If DSTN is required for several large arrays of same size,
%   the weights can be reused to make the algorithm faster. A typical
%   syntax is the following:
%      w = [];
%      for k = 1:10
%          [y{k},w] = dstn(x{k},[],w);
%      end
%   The weights (w) are calculated during the first call of DSTN then
%   reused in the next calls.
% ----------

error(nargchk(1,3,nargin))

y = double(y);
sizy = size(y);

% Test DIM argument
if ~exist('DIM','var'), DIM = []; end
assert(~isempty(DIM) || ~isscalar(DIM),...
    'DIM must be a scalar or an empty array')
assert(isempty(DIM) || DIM==round(DIM) && DIM>0,...
    'Dimension argument must be a positive integer scalar within indexing range.')

% If DIM is empty, a DCT is performed across each dimension

if isempty(DIM), y = squeeze(y); end % Working across singleton dimensions is useless
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1
        if DIM==1, w = []; return
        elseif DIM==2, DIM=1;
        end
        y = y.';
    elseif DIM==2, w = []; return
    end
end

% Weighing vectors
if ~exist('w','var') || isempty(w)
    w = cell(2,dimy);
    for dim = 1:dimy
        if ~isempty(DIM) && dim~=DIM, continue, end
        n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
        w{1,dim} = exp(1i*(0:n-1)'*pi*2/n);
        w{2,dim} = exp(1i*(1:n)'*pi/2/n);
    end
end

% --- DST algorithm ---
if ~isreal(y)
    y = complex(dstn(real(y),DIM,w),dstn(imag(y),DIM,w));
else
    for dim = 1:dimy
        if ~isempty(DIM) && dim~=DIM
            y = shiftdim(y,1);
            continue
        end        
        siz = size(y);
        n = siz(1);
        y = y([1:2:n 2*floor(n/2):-2:2],:);
        y(ceil(n/2)+1:n,:) = -y(ceil(n/2)+1:n,:);
        y = reshape(y,n,[]);
        y = y*sqrt(2*n);
        y = bsxfun(@times,y,w{1,dim});
        y = ifft(y,[],1);
        y = bsxfun(@times,y,w{2,dim});     
        y = imag(y);
        y(n,:) = y(n,:)/sqrt(2);
        y = reshape(y,siz);
        y = shiftdim(y,1);
    end
end
        
y = reshape(y,sizy);



