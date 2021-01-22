function BasicPCA()
global fsz
close all
fsz = 16;
dim = 3;
fname = 'PacmanData.mat';
fname = 'CatData.mat';
dat = load(fname);
N = dat.Ndata; % the number of data points
D = dat.m; % the dimension of the space where data points live
X = dat.data; % N-by-D matrix of data points, data points are rows
d = dat.d; % the distance matrix
ang = dat.a; % rotation anglesfor data points
%% center the data so that their center of mass is zero
colmeans = mean(X,1);
X = X - ones(N,1)*colmeans;
[U,Sigma,V] = svd(X','econ');
esort = diag(Sigma);
figure;
plot(esort,'.','Markersize',20);
grid;
%
[col,X] = getcolors(X);
% make angles form 0 to 2*pi
ind = find(ang < 0);
if ~isempty(ind)
    ang(ind) = ang(ind) + 2*pi;
end
Y1 = X*U(:,1);
Y2 = X*U(:,2);
Y3 = X*U(:,3);
figure;
hold on; grid;
if dim == 2 
    for i = 1 : N
        plot(Y1(i),Y2(i),'.','Markersize',20,'color',col(i,:));
    end   
    colorbar
    caxis([0,pi]);
    set(gca,'Fontsize',20);
    daspect([1,1,1])
else
    for i = 1 : N
        plot3(Y1(i),Y2(i),Y3(i),'.','Markersize',20,'color',col(i,:));
    end   
    colorbar
    caxis([0,pi]);
    view(3)
    set(gca,'Fontsize',20);
    daspect([1,1,1])
end    
end

 %%
 function [c,X] = getcolors(X)
 global fsz
 [n,d] = size(X);
%% determine colors
% compute pairwise distances
d = zeros(n);
e = ones(n,1);
for i = 1 : n
    d(i,:) = sqrt(sum((X - e*X(i,:)).^2,2));
end
% find k nearest neighbors and define weighted directed graph
k = 10; % the number of nearest neighbors for computing distances
% for each point, find k nearest neighbors
ineib = zeros(n,k);
dneib = zeros(n,k);
for i = 1 : n
    [dsort,isort] = sort(d(i,:),'ascend');
    dneib(i,:) = dsort(1:k);
    ineib(i,:) = isort(1:k);
end
% find a point closest to the origin
dor = sum(X.^2,2);
[~,imin] = min(dor);
% compute shortest paths in the graph
D = zeros(n);
ee = ones(1,k);
g = ineib';
g = g(:)';
w = dneib';
w = w(:)';
G = sparse(kron((1:n),ee),g,w);
c = zeros(n,3);
[dist,~,~] = graphshortestpath(G,imin);
[dist,isort] = sort(dist,'ascend');
X = X(isort,:);
N = 1000;
col = parula(N);
for i = 1 : n
    if isfinite(dist(i))
        c(i,:) = col(getcolor(dist(i),max(dist(isfinite(dist))),N),:);
    else
        c(i,:) = [1,0,0];
    end
end

end
 %%
function c = getcolor(u,umax,N)
c = max(1,round(N*(u/umax)));
end


    