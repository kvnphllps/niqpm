
function [BWfinal,Segout] = edge_detect_kp_1march2012(hilbertTrans, fudgeFactor)

%% Section 1. Detect Entire Cell

%I = inFocus; 
% [I] = hilbert_transform_dic(inFocus);
% I = (hilbertTrans).^1;
% clean = I;
%  final_clean = clean.*(clean>0);
%     final_clean = final_clean./max(final_clean(:));
I = hilbertTrans;
[junk threshold] = edge(I, 'roberts');
%fudgeFactor = .2;
BWs = edge(I,'sobel', threshold * fudgeFactor);
%figure, imshow(BWs), title('binary gradient mask');

% Dilate the image
se90 = strel('line', 6, 90);
se0 = strel('line', 6, 0);

BWsdil = imdilate(BWs, [se90 se0]);
%figure, imagesc(BWsdil), title('dilated gradient mask');

BWdfill = imfill(BWsdil, 'holes');

%figure, imagesc(BWdfill);
%title('binary image with filled holes');

BWnobord = imclearborder(BWdfill, 4);
%figure, imagesc(BWnobord), title('cleared border image');


seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
%figure, imagesc(BWfinal), title('segmented image');


 BWoutline = bwperim(BWfinal);
 Segout = I;
 Segout(BWoutline) = 1;
% figure, imagesc(Segout), title('outlined original image');
% colormap(gray)
% colorbar
% axis equal image
% 
% figure, imagesc(I)
% axis equal image
% 
% figure, imagesc(inFocus)
% axis equal image




