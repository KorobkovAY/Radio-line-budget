clear all
close all
clc


%% Входные данные
R  = 6000;     % m
F  = (1:0.1:150)*10^9;    % Hz
ht = 200;    % m
el = 5;    % deg
for ch = 1:length(F)
    for i = 1:length(el)
        Lgas(i,ch) = tropopl(R,F(ch),ht,el(i),'LatitudeModel','Mid','Season','Winter');
    end
end
%% Построение графика
% plot(el,Lgas(:,1));
% grid on
% hold on
% plot(el,Lgas(:,2));
% plot(el,Lgas(:,3));
% xlabel('Угол места, град');
% ylabel('Ослабление, дБ');
% legend({'F = 100 ГГц', 'F = 200 ГГц', 'F = 300 ГГц'},'Interpreter','tex','Location','northeast','FontSize',12);
% % ylim([1e-3 1e5]);

semilogy(F/10^9,Lgas)
xlabel('Частота, ГГц');
ylabel('Ослабление, дБ');
grid on
xlim([1 160]);