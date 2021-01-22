data=load('data.mat');
face=data.face;
face_neutral=face(:,:,1:3:end);
face_exp=face(:,:,2:3:end);
face_illum=face(:,:,3:3:end);

data=load('pose.mat');
pose=data.pose;

figure;
colormap gray;
for j=1:13
    subplot(3,5,j);
    imagesc(pose(:,:,j,1));
end

figure;
colormap gray;
for j=1:20
    subplot (4,5,j)
    imagesc(face_neutral(:,:,j))
end 

figure;
colormap gray;
for j=1:20
    subplot (4,5,j)
    imagesc(face_exp(:,:,j))
end 

figure; 
colormap gray;
for j=1:20
    subplot (4,5,j)
    imagesc(face_illum(:,:,j))
end 

%Fisher linear discriminant (not used)
[d1,d2,n] = size(face_neutral);
X1=zeros(n,d1*d2); %data for neutral faces
X2=zeros(n,d1*d2); %data for expressions
for j=1:n
    aux = face_neutral(:,:,j);
    X1(j,:)=aux(:)';
    aux=face_exp(:,:,j);
    X2(j,:)=aux(:)';
end

X=[X1;X2]

%PCA to map to higher dimensional space
[U,Sig,V] = svd(X','econ');
nPCA = 20;
%project to nPCA-dimensional space
Y = X*U(:,1:nPCA);
figure;
hold on;
grid;
plot3(Y(D1,1), Y(D1,2), Y(D1,3),'.', 'Markersize', 20, 'color', 'k');
plot3(Y(D2,1), Y(D2,2), Y(D2,3),'.', 'Markersize', 20, 'color', 'r');
view(3);

%Bayesian decision theory
Y1=Y(D1,:);%projected data of neutral face
Y2=Y(D2,:);%projected data of expression face
mu1=mean(Y1,1);
mu2=mean(Y2,1);
fprintf('norm(mu1-mu2)=%d\n', norm(mu1-mu2));

%estimate covariance matrices for each class
Y1c=Y1-ones(n1,1)*mu1; %center the data for neutral face
Y2c=Y2-ones(n2,1)*mu2; %center the data for expression face
S1=Y1c'*Y1c/n1;
S2=Y2c'*Y2c/n2;
figure;
imagesc(S1);
colorbar;
figure;
imagesc(S2);
colorbar;

%define the descriminant functions
iS1=inv(S1);
iS2=inv(S2);
mu1=mu1';
mu2=mu2';
w0=0.5*(log(det(S2)/det(S1)))-0.5*(mu1'*iS1*mu1-mu2'*iS2*mu2);
g=@(x)-0.5*x'*(iS1-iS2)*x+x'*(iS1*mu1-iS2*mu2)+w0;

%classify illuminated faces
Y3=zeros(n,nPCA);
label=zeros(n,1);
for j=1:n
    aux=face_illum(:,:,j);
    y=(aux(:)'*U(:,1:nPCA))'; 
    Y3(j,:)=y';
    label(j)=sign(g(y));
end
iplus=find(label>0);
iminus=find(label<0);
fprintf('#iplus=%d, #iminus=%d\n', length(iplus), length(iminus));