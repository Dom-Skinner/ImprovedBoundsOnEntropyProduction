function X = get_traj_MT(k)
files = dir;
file = files(k).name;
I = imread(file);
I = I/max(I(:));

% ========================================================================
% Remove weird rows of pixels usually at top or bottom of image
row_mean = mean(I,2);
mean_row_mean = mean(row_mean);
std_row_mean = std(row_mean);

row_keep = abs(row_mean - mean_row_mean) < 3*std_row_mean ;

I = I(row_keep,:);
I = imadjust(I);
% ========================================================================


% ========================================================================
% cut off left of image so origin is at border
T = graythresh(I);
BW = imbinarize(I,T);
BW2 = bwareafilt(BW,1);
idx = find(mean(BW2,1) > 0.4,1, 'first');
if ~any(idx)
    idx = find(mean(BW2,1) > 0.2,1, 'first');
    disp("worth checking to see if error")
end
if ~any(idx)
    idx = find(mean(BW2,1) > 0.1,1, 'first');
    disp("worth checking to see if error")
end
I = imadjust(I(:,idx:end));
% ========================================================================


% ========================================================================
T = graythresh(I);
BW = imbinarize(I,T);

CC = bwconncomp(BW);
rp = regionprops(BW,'Extrema');
rem = [];
for i = 1:length(rp)
    rem(i) =  (min(rp(i).Extrema(:,1)) > 2) & (max(rp(i).Extrema(:,1)) > size(I,2) - 1);
end
BW3 = ismember(labelmatrix(CC),find(~rem)); 
BW3 = imdilate(BW3,strel('square',3));

% ========================================================================

X = zeros(size(BW3,1),1);
for i = 1:size(BW3,1)
    set = 0;
    if all(BW3(i,1:3) == 0)
        set = 1;
    end
    for j = 2:size(BW3,2)
        if (set == 0) && (BW3(i,j-1) == 1) && all(BW3(i,j:min(j+2,size(BW3,2))) == 0)
            X(i) = j-1;
            set = 1;
        end
    end
    if set == 0
        X(i) = size(BW3,2);
    end
end
end