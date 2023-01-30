clear all
close all
clc


%% Входные данные
% Дождь
hstation = 0.2; % Высота НС над уровнем моря, км
teta = 5; % Угол места
%Re = effearthradius; % Эффективный радиус Земли (~8500 км)
altitude = 55.75; % Широта НС
Re = 8500;
R001 = 10; % Интенсивность дождя, мм/ч
h0 = 3; % Среднегодовая высота изотермы 0 град
% Мерцания
KPD  = 0.7;     % КПД зеркала
D = 5;   % Диаметр зеркала
F  = 5;    % GHz
t = 10; % град С
P = 1013.3; % гПа
p = 0.1;  % Процент времени, за которое рассматривается ослабление 0.01<p<=50
% Облака
T = 273.15; % Температура жидкой воды в облаках, К
Lred = 0.5; % Общий столбчатый объем жидкой воды, кг/м3 (для января 1% - 0.5)
% Газы
R  = 100000;     % m
ht = 200;    % m
%% Расчет

for ch = 1:length(F)
    for i = 1:length(teta)
        Lgas(i,ch) = tropopl(R,F(ch)*10^9,ht,teta(i),'LatitudeModel','Mid','Season','Winter');
    end
end
for i = 1:length(F)
    Fita = 300/T;
    E0 = 77.66+103.3*(Fita-1);
    E1 = 0.0671*E0;
    E2 = 3.52;
    fp = 20.20-146*(Fita-1)+316*(Fita-1)^2; % Главная частота релаксации, ГГц
    fs = 39.8*fp; % Вторичная частота релаксации, ГГц
    Ad1 = (E0-E1)/(1+(F(i)/fp)^2);
    Ad2 = (E1-E2)/(1+(F(i)/fs)^2);
    Es = Ad1+Ad2+E2;
    Ad1 = (F(i)*(E0-E1))/((1+(F(i)/fp)^2)*fp);
    Ad2 = (F(i)*(E1-E2))/((1+(F(i)/fs)^2)*fs);
    Ess = Ad1+Ad2;
    Eta = (2+Es)/Ess;
    Kl = 0.819*F(i)/(Ess*(1+Eta^2)); % Коэффициент погонного ослабления, ((дБ/км)/(г/м3))
    A_cloud(i) = Lred*Kl/sind(teta); % Ослабление из-за облачности на наклонных трассах, дБ
end

KAlfa = readmatrix('alfa.xlsx');
for coun2 = 1:length(teta)
    for coun = 1:length(F)
        hr = h0+0.36; % Высота дождя
        if teta(coun2) >= 5 % Высота земной станции над средним уровнем моря, км
            Ls = (hr-hstation)/sind(teta(coun2));
        else
            Ls = 2*(hr-hstation)/((sind(teta(coun2))^2+2*(hr-hstation)/Re)^0.5+sind(teta(coun2)));
        end
        LG = Ls*cosd(teta(coun2)); % Горизонтальная проекция трассы
        for i=1:length(KAlfa(:,1))
            if F(coun) >= KAlfa(i,1)
                KamH = KAlfa(i,2);
                alfaH = KAlfa(i,3);
                KamV = KAlfa(i,4);
                alfaV = KAlfa(i,5);
            elseif F(coun) <1
                KamH = KAlfa(1,2);
                alfaH = KAlfa(1,3);
                KamV = KAlfa(1,4);
                alfaV = KAlfa(1,5);
            end
        end
        Kam = (KamH+KamV)/2;
        AlFa = (alfaH*KamH+alfaV*KamV)/(2*Kam);
        gammaR = Kam*R001^AlFa; % Погонное ослабление
        roo1 = 1/(1+0.78*sqrt(LG*gammaR/F(coun))-0.38*(1-exp(-2*LG))); % К-т ослабления по горизонтали
        Etta = atand((hr-hstation)/(LG*roo1)); % К-т подстройки по вертикали
        if Etta > teta(coun2)
            Lr = (LG*roo1)/cosd(teta(coun2));
        else
            Lr = (hr-hstation)/sind(teta(coun2));
        end
        if altitude <36
            Hi = 36-abs(altitude);
        else
            Hi = 0;
        end
        mult = (1-exp(-teta(coun2)/(1+Hi)))*sqrt(Lr*gammaR)/F(coun)^2;
        Nu = 1/(1+sqrt(sind(teta(coun2)))*(31*mult-0.45));
        LE = Lr*Nu; % Эффективная длина трассы
        A001(coun,coun2) = gammaR*LE; % Прогнозируемое значение ослабления
    end
        
end

for ch = 1:length(teta)
    for i = 1:length(F)
        es =  (1+10^(-4)*(7.2+P*(0.00320+5.9*10^(-7)*t^2)))*6.1121*exp((t*(18.678-t/234.5))/(t+257.14));     % Давление насыщенного пара
        e = es*80/100; % Давление водяных паров
        T = 273.15+t; % Абсолютная температура
        Nwat = 3.732e5*e/(T^2); % Влажностна составляющая
        sigmaREF = 3.6*10^(-3)+Nwat*10^(-4);   % Стандартное отклонение амплитуды сигнала
        L = 2000/(sqrt(sind(teta(ch))^2+2.35*10^(-4))+sind(teta(ch)));  % Эффективная длинна трассы
        Deff = sqrt(KPD)*D;
        x = 1.22*Deff^2*(F(i)/L);
        gx = sqrt(3.86*(x^2+1)^(11/12)*sind(11*atand(1/x)/6)-7.08*x^(5/6));
        sigma = sigmaREF*F(i)^(7/12)*gx/(sind(teta(ch)))^1.2;
        ap = -0.061*(log10(p)^3)+0.072*(log10(p))^2-1.71*log10(p)+3;  % К-т процента времени
        As(i,ch) = ap*sigma;
    end
end

Lsv = 92.45+20*log10(F*R);

bar(1, Lsv)
grid on
hold on
bar(2 ,As)
bar(3 ,A001)
bar(4, A_cloud)
bar(5, Lgas)
text(1.4,2.75,'\leftarrow 122 дБ', 'FontSize',11)
xlabel('');
ylabel('Ослабление, дБ');
legend({'Ослабление в свободном пространстве', 'Ослабление в следствии мерцаний', 'Ослабление в следствии дождя', 'Затухание в облаках', 'Затухание в атмосферных газах'},'Interpreter','tex','Location','northeast','FontSize',12);
ylim([0 5])
xlim([0 10])

