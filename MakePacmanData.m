function MakeData()
% Makes Pacman data
close all
m = 65;
phi = pi/6; % opening half-angle
r = 0.8;
r2 = r*r;
xx = linspace(-1,1,m);
[x,y] = meshgrid(xx,xx);
im = ones(m,m)*255;
dd2 = x.^2 + y.^2;
az = angle(x + 1i*y);
izplus = find(az >= 0);
izminus = find(az < 0);
ind = find(dd2 < r2 & abs(az) > phi); 
im(ind) = 0;
figure;
colormap gray
image(xx,xx,im);
daspect([1,1,1])
%%
figure;
for k = 1 : 20
    a = 2*pi*rand - pi;
    a1 = phi + a;
    azz = az;
    if a1 > pi
        azz(izminus) = azz(izminus) + 2*pi;
    end
    a2 = -phi + a;
    if a2 < -pi
        azz(izplus) = azz(izplus) - 2*pi;
    end
    im = ones(m,m)*255;
    inew = find(dd2 < r2 & (azz > a1 | azz < a2 )); 
    im(inew) = 0;
    subplot(4,5,k);
    colormap gray
    image(xx,xx,im);
end
%% make data
Ndata = 200;
m2 = m*m;
data = zeros(Ndata,m2);
a = linspace(-pi,pi,Ndata+1);
a(end)=[];
for k = 1 : Ndata
    azz = az;
    a1 = phi + a(k);
    if a1 >= pi
        azz(izminus) = azz(izminus) + 2*pi;
    end
    a2 = -phi + a(k);
    if a2 < -pi
        azz(izplus) = azz(izplus) - 2*pi;
    end
    im = ones(m,m)*255;
    inew = find(dd2 < r2 & (azz > a1 | azz < a2 )); 
    im(inew) = 0;
    aux = im(:);
    data(k,:) = aux';
end

for i = 1 : Ndata
    d(:,i) = sum((data - ones(Ndata,1)*data(i,:)).^2,2);
end
save('PacmanData.mat','data','Ndata','m','d','a');
end





