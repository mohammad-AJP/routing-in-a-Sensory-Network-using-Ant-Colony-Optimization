clc
clear all
close all
t = 0;
for i=1:10
    for j = 1:10
        x(i,j) = floor(110*rand+t);
        y(i,j) = floor(110*rand+t);
    end
    t = t+110;
end
X = reshape(x',1,100);
Y = reshape(y,1,100);
ini_energy = 10*randi([1 , 3],1,100);
nodes = [X;Y;ini_energy];
plot(X,Y,'LineStyle','none','marker' , 's','MarkerSize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
