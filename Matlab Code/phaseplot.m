M = dlmread('Phase Plot/M_15_300_pm1.txt');
 M = fliplr(flipud(M));
 M = M';
 close all
figure
imagesc(kron(M/10,ones(20,20)))
hold on
colormap('gray')

S = 0.01:.001:15;
L = 2*2.*S.*(log(200*S/5));
L = (L/4)*20-10;
S = (15-S)*20+10;

S2 = .01:.001:15;
L2 = .25.*2.*S2.*(log(200*S2/5)).^2;
L2 = (L2/4)*20-10;
S2 = (15-S2)*20+10;


% plot(L,S,'r','LineWidth',3)
% x = [0.7 .65];
% y = [.45 .35];
% str = '$$ L = 4S\log(N+K) $$';
% annotation('textarrow',x,y,'Interpreter','latex','FontSize',12,'String',str)

plot(L2,S2,'b','LineWidth',3)
x = [0.65 .68];
y = [.28 .41];
str = '$$ L = 0.25(S_1+S_2)\log^2(N+K) $$';
annotation('textarrow',x,y,'Interpreter','latex','FontSize',24,'String',str)



yticks = 10:40:300;
yticklabels = 15:-2:1;
 set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
 xticks =10:40:750;
 xticklabels = 4:8:140;
 set(gca, 'XTick', xticks,'FontSize',24);
 set(gca,'xticklabel',xticklabels)%,sprintf('%d |',xticklabels));
xlabel('$L$','Interpreter','Latex','FontSize',24)
ylabel('$S=S_1=S_2$','Interpreter','Latex','FontSize',24)  
hold off
