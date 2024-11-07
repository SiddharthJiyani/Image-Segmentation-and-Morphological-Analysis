img = imread("./img1.jpg");

%Step2 : Convert to grayscale

%getting the dimention of the image
[rows, cols,channels] = size(img);
%fprintf('Image dimensions: Rows = %d, Columns = %d, Channels = %d\n', rows, cols, channels);

    
% Initialize the grayscale image matrix
gray_img = zeros(rows, cols);

% Convert RGB to grayscale
for i = 1:rows
    for j = 1:cols
        % Convert RGB to grayscale using the formula: Gray = 0.29 * R + 0.59 * G + 0.11 * B
        gray_img(i,j) = 0.29 * double(img(i,j,1)) + 0.59 * double(img(i,j,2)) + 0.11 * double(img(i,j,3));
    end
end 

% Convert to uint8 data type
% 'uint8(gray_img)' ensures that the grayscale image matrix gray_img is represented as 8-bit integers, which is the standard format for images in MATLAB 
gray_img = uint8(gray_img);

% Display the grayscale image
%figure;
%imshow(gray_img);
%title('Grayscale Image');


% Step3 : Add 30% noice to image

noisy_img = gray_img;
for i=1:rows
    for j=1:cols
        %generate a random number
        random_num = rand(1);
        if random_num < 0.1
            noisy_img(i,j) = 255;
        end

        if random_num >0.1 && random_num < 0.2
            noisy_img(i,j) = 0 ;
        end
        
    end
end

% Step 4: De-noising image using adaptive median filter

% Define maximum window size for the adaptive median filter
max_winsize = 10;

% Initialize denoised image
deNoise_image = noisy_img;

% Iterate over each pixel
for i = 1:rows
    for j = 1:cols
        % Compute maximum window size
        win_s = min([max_winsize, i-1, j-1, rows-i, cols-j]);
        
        % Iterate over window sizes
        for W = 0:win_s
            % Extract local neighborhood
            S = noisy_img(i-W:i+W, j-W:j+W);
            
            % Calculate median and extreme values
            xmed = median(S, 'all');
            xmin = min(S, [], 'all');
            xmax = max(S, [], 'all');
            
            % Check conditions for noise removal
            if (xmed - xmin > 0 && xmax - xmed > 0)
                break; % Exit loop if conditions are met
            end
        end
        
        % Apply median value to denoised image
        if ~(noisy_img(i, j) - xmin > 0 && noisy_img(i, j) - xmax < 0)
            deNoise_image(i, j) = xmed;
        end
    end
end

% Step 5: Otsu's Thresholding

Histogram = imhist(deNoise_image); %computing histogram

max = 0 ;
P = Histogram/sum(Histogram); % Histogram is normalized
var(1) = 0 ; % between class variance
threshold=0;

for k=1:255    % Moving through each and every pixel
    p1 = sum(P(1:k));        %Probability of class 1 (seperated by threshold)
    p2 = sum(P(k+1:255));    %Probability of class 1 (seperated by threshold)
    
    mean1 = sum((1:k)'.*P(1:k))/p1 ; %class mean 1 = (1/pi)*(sum(i*pi))
    mean2 = sum((k+1:255)'.*P(k+1:255))/p2 ; %class mean 

    var(k) = p1*p2*((mean1-mean2)^2) ; % between class variance using simplified formula
    if(var(k) > max)
        max = var(k);
        threshold = k-1 ;
    end
end

% Apply thresholding to obtain binary image
binary_img = deNoise_image > threshold;

% Step 6.1 : Morphological analysis - Erosion

% Define the structuring element for erosion
SE = ones(5); % 5x5 square structuring element

% Create a zero matrix for the eroded image
erosion_img = zeros(size(binary_img));

% Perform erosion operation
for i = 3:size(binary_img, 1)-2
    for j = 3:size(binary_img, 2)-2
        % Extract the neighborhood around the current pixel
        neighborhood = binary_img(i-2:i+2, j-2:j+2);
        
        % Apply the structuring element (logical AND operation) and find the minimum value
        erosion_img(i, j) = min(neighborhood(SE == 1), [], 'all');
    end
end


% Step 6.2 : Morphological analysis - Dilation

% Create a zero matrix for the dilated image
dilation_img = zeros(size(erosion_img));

% Perform dilation operation
for i = 3:size(erosion_img, 1)-2
    for j = 3:size(erosion_img, 2)-2
        % Extract the neighborhood around the current pixel
        neighborhood = erosion_img(i-2:i+2, j-2:j+2);
        
        % Check if any pixel in the neighborhood is 1
        if any(neighborhood(SE == 1))
            dilation_img(i, j) = 1; % Set the corresponding pixel in the dilation image to 1
        end
    end
end

% Step 7 : Extracting largest connected component
X0=zeros(rows,cols,"double");
X0(100,100)=1;
X_prev=zeros(rows,cols,"double");
Xk=zeros(rows,cols,"double");

while(true)
    %Dilation
    X_prev=zeros(rows,cols,"double");
    SE=ones(5);
    for i = 3:(size(X0, 1)-3)
        for j = 3:(size(X0, 2)-3)
            se_win_row = i-2:i+2;
            se_win_col = j-2:j+2;

            if(X0(i,j)==1)
                X_prev(se_win_row,se_win_col)=SE;
            end
        end
    end
    Xk=X_prev.*dilation_img;
    if(Xk==X0)
        break;
    end
    X0=Xk;
end

% Calculate the area of the largest connected component
area = sum(Xk(:));

% Display the largest connected component
figure;
imshow(Xk);
caption = sprintf('Area = %d', area);
title(sprintf('Largest Connected Component (Area = %d)', area));