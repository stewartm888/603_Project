function Cat()
close all
c = imread('cat_in_the_hat.jpg');
[m,n] = size(c);
% m must be equal to n
if m ~= n
    m = min(m,n);
end
c1 = c(1 : m, 1 : m,1);
figure;
imshow(c1);
xx = linspace(-1,1,m);
[x,y] = meshgrid(xx,xx);
im = ones(m,m)*255;
ind = find(c1 < 50);
im(ind) = 0;




figure;
colormap gray
image(xx,xx,im);
%%
figure;
for k = 1 : 20
    a = 2*pi*rand;
    z = [x(ind)';y(ind)'];
    R = [cos(a),-sin(a);sin(a),cos(a)];
    znew = R*z;
    xnew = max(-1,min(znew(1,:),1));
    ynew = max(-1,min(znew(2,:),1));
    inew = sub2ind([m,m],1+round((1 + ynew)*(m/2 - 1)),1+round((1 + xnew)*(m/2 - 1)));
    inew = unique([inew;inew + 1;inew - 1;inew + m;inew - m]);
    im = ones(m,m)*255;
    im(inew) = 0;
    subplot(4,5,k);
    colormap gray
    image(xx,xx,im);
end
%% make data
Ndata = 200;
m2 = m*m;
data = zeros(Ndata,m2);
a = 2*pi*rand(Ndata,1);
for k = 1 : Ndata
    z = [x(ind)';y(ind)'];
    R = [cos(a(k)),-sin(a(k));sin(a(k)),cos(a(k))];
    znew = R*z;
    xnew = max(-1,min(znew(1,:),1));
    ynew = max(-1,min(znew(2,:),1));
    inew = sub2ind([m,m],1+round((1 + ynew)*(m/2 - 1)),1+round((1 + xnew)*(m/2 - 1)));
    inew = unique([inew;inew + 1;inew - 1;inew + m;inew - m]);
    im = ones(m,m)*255;
    im(inew) = 0;
    aux = im(:);
    data(k,:) = aux';
end

for i = 1 : Ndata
    d(:,i) = sum((data - ones(Ndata,1)*data(i,:)).^2,2);
end
save('CatData.mat','data','Ndata','m','d','a');
end





