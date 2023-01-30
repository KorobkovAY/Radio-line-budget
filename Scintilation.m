clear all
close all
clc


%% Входные данные
KPD  = 0.9;     % КПД зеркала
D = 12;   % Диаметр зеркала
F  = 1:20;    % Hz
el = [4 8 16];    % deg
t = 10; % град С
P = 1013.3; % гПа
p = 0.05;  % Процент времени, за которое рассматривается ослабление 0.01<p<=50

%% Расчет
for ch = 1:length(el)
    for i = 1:length(F)
        es =  (1+10^(-4)*(7.2+P*(0.00320+5.9*10^(-7)*t^2)))*6.1121*exp((t*(18.678-t/234.5))/(t+257.14));     % Давление насыщенного пара
        e = es*80/100; % Давление водяных паров
        T = 273.15+t; % Абсолютная температура
        Nwat = 3.732e5*e/(T^2); % Влажностна составляющая
        sigmaREF = 3.6*10^(-3)+Nwat*10^(-4);   % Стандартное отклонение амплитуды сигнала
        L = 2000/(sqrt(sind(el(ch))^2+2.35*10^(-4))+sind(el(ch)));  % Эффективная длинна трассы
        Deff = sqrt(KPD)*D;
        x = 1.22*Deff^2*(F(i)/L);
        gx = sqrt(3.86*(x^2+1)^(11/12)*sind(11*atand(1/x)/6)-7.08*x^(5/6));
        sigma = sigmaREF*F(i)^(7/12)*gx/(sind(el(ch)))^1.2;
        ap = -0.061*(log10(p)^3)+0.072*(log10(p))^2-1.71*log10(p)+3;  % К-т процента времени
        As(i,ch) = ap*sigma;
    end
end
semilogy(F,As(:,1));
grid on
hold on
semilogy(F,As(:,2));
semilogy(F,As(:,3));
xlabel('Частота, ГГц');
ylabel('Ослабление, дБ');
legend({'\theta = 4\circ', '\theta = 8\circ', '\theta = 16\circ'},'Interpreter','tex','Location','northeast','FontSize',12);
ylim([0 10]);
