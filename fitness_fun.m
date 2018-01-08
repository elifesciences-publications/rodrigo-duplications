%..
f0 = 1;
a = 0.17;
k = 0.40;
b = 0.036;
h = 1.8;
n = 4;
x0 = 0.13;
s = 0.5;
sex = 0;

%..
x = x0*lognrnd(0,1,1,1e4);
noise1 = lognrnd(0,s,1,1e4);
noise2 = lognrnd(0,s,1,1e4);
noise0 = lognrnd(0,sex,1,1e4);
y1 = 1./(1 + (x.*noise1.*noise0/x0).^(-n));
y2 = y1/2 + (1/2)./(1 + (x.*noise2.*noise0/x0).^(-n));

%..
xdet = x0*exp(-3:0.1:3);
ydet = 0:0.02:1.2;
ydetmodel = 1./(1 + (xdet/x0).^(-n));
ydetopt = h - sqrt(b*h*(k+xdet)./(a*xdet));
figure('color','w');
semilogx(xdet/x0,ydetmodel,'k','LineWidth',3);
hold on;
semilogx(xdet/x0,ydetopt,'k--','LineWidth',3);
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 1;
axis(ax,'square');
xlabel('lactose (x/x_0)','FontSize',30);
ylabel('LacZ (y)','FontSize',30);
xlim([0.1 10]);
ylim([0 1.2]);

%..
for i=1:length(xdet)
    for j=1:length(ydet)
        Fdet(i,j) = f0*(1 + a*ydet(j).*xdet(i)./(k+xdet(i)) - b*ydet(j)./(h-ydet(j)));
    end
end
figure('color','w');
surf(xdet/x0,ydet,Fdet');
mymap=[(0.2:0.05:1)' (0.2:0.05:1)' (0.2:0.05:1)'].^0.75;
%mymap='jet';
colormap(mymap);
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 1;
ax.XAxis.Scale = 'log';
axis(ax,'square');
xlabel('lactose (x/x_0)','FontSize',30,'rot',20);
ylabel('LacZ (y)','FontSize',30,'rot',-40);
zlabel('fitness (W)','FontSize',30);
xlim([0.1 10]);
ylim([0 1.2]);
zlim([0.95 1.1]);
ax.YTick = [0 0.6 1.2];


