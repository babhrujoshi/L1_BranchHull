close all
result = load('L1_BH_SNR.mat').M;

alpha_list = result(1,:);
rel_error_mean = result(2,:);
rel_error_std = result(3,:);
max_sup_noise = result(4,:);


% prepare it for the fill function
X_plot  = [alpha_list, fliplr(alpha_list)];
Y_plot  = [rel_error_mean-rel_error_std, fliplr(rel_error_mean+rel_error_std)];


hold on 
plot(alpha_list, rel_error_mean, 'blue', 'LineWidth', 1.2)
% plot(alpha_list, log2(37*sqrt(max_sup_noise)), 'black', 'LineWidth', 1.2)
fill(X_plot, Y_plot , 1,....
        'facecolor','blue', ...
        'edgecolor','none', ...
        'facealpha', 0.3);
xlabel('$\alpha$, noise level','Interpreter','Latex','FontSize',24)
ylabel('Relative error','Interpreter','Latex','FontSize',24)  
hold off 

figure
hold on 
plot(alpha_list, log2(rel_error_mean), 'blue', 'LineWidth', 1.2)
plot(alpha_list, log2(37*sqrt(max_sup_noise)), 'black', 'LineWidth', 1.2)
fill(X_plot, log2(Y_plot) , 1,....
        'facecolor','blue', ...
        'edgecolor','none', ...
        'facealpha', 0.3);
xlabel('$\alpha$, noise level','Interpreter','Latex','FontSize',24)
ylabel('Relative error','Interpreter','Latex','FontSize',24)  

set(gca, 'YTickLabel', str)   
xl = get(gca, 'XLim');
yt = get(gca, 'YTick');
str = cellstr( num2str(yt(:),'2^{%d}') );
legend('Relative error', 'Theoritical bound', 'Interpreter','Latex','FontSize',24)

% hTxt = text(yt, xl(ones(size(yt))), str, ...   %# create text at same locations
%     'Interpreter','tex');    

hold off
