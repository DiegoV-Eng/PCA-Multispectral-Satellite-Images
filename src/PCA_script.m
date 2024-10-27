clear all;
close all;
clc;

% Loading images in variables
I1 = imread('Fig1138(a)(WashingtonDC_Band1_564).tif');
I2 = imread('Fig1138(b)(WashingtonDC_Band2_564).tif');
I3 = imread('Fig1138(c)(WashingtonDC_Band3_564).tif');
I4 = imread('Fig1138(d)(WashingtonDC_Band4_564).tif');
I5 = imread('Fig1138(e)(WashingtonDC_Band5_564).tif');
I6 = imread('Fig1138(f)(WashingtonDC_Band6_564).tif');

N = size(I1,1);
M = size(I1,2);

% Calculating mean and covariance matrix
Nd = 6;
X = zeros(N,M,Nd);
mx = zeros(Nd,1);
Cx = zeros(Nd,Nd);
for i = 1:N
    for j = 1:M
        X(i,j,:) = [I1(i,j); I2(i,j); I3(i,j); I4(i,j); I5(i,j); I6(i,j)];
        mx = mx + squeeze(X(i,j,:));
        Cx = Cx + squeeze(X(i,j,:))*squeeze((X(i,j,:)))';
    end
end
mx = mx/(N*M); 
Cx = Cx/(N*M)-mx*mx';

% Calculating eigenvalues and eigenvectors
[Q,D] = eig(Cx);
d_v = sort(diag(D),'descend');
Q_or = fliplr(Q);
D_or = diag(d_v);

% Eigenvector Matrix and Principal components of Matrix A (Ak)
A = Q_or'; 

for d = 2:4
    Ak = A(1:d,:); 

    % Projection of the vector x in the new space of characteristics by using
    % the transformation Matrix A
    Y = zeros(N,M,Nd);
    for k= 1:Nd
        for i = 1:N
            for j = 1:M
                X(i,j,:) = [I1(i,j); I2(i,j); I3(i,j); I4(i,j); I5(i,j); I6(i,j)];
                Y(i,j,k) = A(k,:)*(squeeze(X(i,j,:))-mx);
            end
        end
    end

    % Reconstruct vectors (X_hat) by using some eigenvectors stored in "Ak"
    X_hat = zeros(N,M,Nd);
    Y_hat = [];
    for i = 1:N
        for j = 1:M
            X(i,j,:) = [I1(i,j); I2(i,j); I3(i,j); I4(i,j); I5(i,j); I6(i,j)];
            Y_hat(i,j,:) = Ak*(squeeze(X(i,j,:))-mx);
            X_hat(i,j,:) = Ak'*squeeze(Y_hat(i,j,:))+mx;
        end
    end

    % Deploy images
    figure(1); % Original images
    subplot(2,3,1); imshow(I1); axis off
    subplot(2,3,2); imshow(I2); axis off
    subplot(2,3,3); imshow(I3); axis off
    subplot(2,3,4); imshow(I4); axis off
    subplot(2,3,5); imshow(I5); axis off
    subplot(2,3,6); imshow(I6); axis off
    a = axes; t1 = title('Original Images');
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on');

    figure(2); % PCA Images
    subplot(2,3,1); imshow(Y(:,:,1),[]); axis off
    subplot(2,3,2); imshow(Y(:,:,2),[]); axis off
    subplot(2,3,3); imshow(Y(:,:,3),[]); axis off
    subplot(2,3,4); imshow(Y(:,:,4),[]); axis off
    subplot(2,3,5); imshow(Y(:,:,5),[]); axis off
    subplot(2,3,6); imshow(Y(:,:,6),[]); axis off
    a = axes; t1 = title('PCA Images (per component)');
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on');

    figure(3); % Reconstructed images
    subplot(2,3,1); imshow(X_hat(:,:,1),[]); axis off
    subplot(2,3,2); imshow(X_hat(:,:,2),[]); axis off
    subplot(2,3,3); imshow(X_hat(:,:,3),[]); axis off
    subplot(2,3,4); imshow(X_hat(:,:,4),[]); axis off
    subplot(2,3,5); imshow(X_hat(:,:,5),[]); axis off
    subplot(2,3,6); imshow(X_hat(:,:,6),[]); axis off
    a = axes; t1 = title('Reconstructed Images');
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on');

    figure(4); % Measure error
    subplot(2,3,1); imshow(X_hat(:,:,1)-double(I1)); axis off
    subplot(2,3,2); imshow(X_hat(:,:,2)-double(I2)); axis off
    subplot(2,3,3); imshow(X_hat(:,:,3)-double(I3)); axis off
    subplot(2,3,4); imshow(X_hat(:,:,4)-double(I4)); axis off
    subplot(2,3,5); imshow(X_hat(:,:,5)-double(I5)); axis off
    subplot(2,3,6); imshow(X_hat(:,:,6)-double(I6)); axis off
    a = axes; t1 = title('Difference of the Images');
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on');

    % Measure the error and calculating PSNR
    err1 = immse(X_hat(:,:,1),double(I1)); % the format of the matrices should be the same
    PSNR1 = 10*log10(255^2/err1);

    err2 = immse(X_hat(:,:,2),double(I2));
    PSNR2 = 10*log10(255^2/err2);

    err3 = immse(X_hat(:,:,3),double(I3));
    PSNR3 = 10*log10(255^2/err3);

    err4 = immse(X_hat(:,:,4),double(I4));
    PSNR4 = 10*log10(255^2/err4);

    err5 = immse(X_hat(:,:,5),double(I5));
    PSNR5 = 10*log10(255^2/err5);

    err6 = immse(X_hat(:,:,6),double(I6));
    PSNR6 = 10*log10(255^2/err6);

    % Making table to present results
    Variable1 = {'Error 1'; 'Error 2'; 'Error 3'; 'Error 4'; 'Error 5'; 'Error 6'};
    values1 = [err1; err2; err3; err4; err5; err6];

    Variable2 = {'PSNR 1'; 'PSNR 2'; 'PSNR 3'; 'PSNR 4'; 'PSNR 5'; 'PSNR 6'};
    values2 = [PSNR1; PSNR2; PSNR3; PSNR4; PSNR5; PSNR6];
    T = table(Variable1,values1,Variable2,values2);
    disp(['Results using best', string(d), 'principal components'])
    disp(T);

end
