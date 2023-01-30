close all
clear all
clc

%% Входные данные
hstation = 0.2; % Высота НС над уровнем моря, км
teta = 5:3:20; % Угол места
F = 0.1:1:15.1; %Частота, ГГц
%Re = effearthradius; % Эффективный радиус Земли (~8500 км)
altitude = 55.75; % Широта НС
Re = 8500;
R001 = 10; % Интенсивность дождя, мм/ч
h0 = 3; % Среднегодовая высота изотермы 0 град

%% Расчет
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

figure(1)
semilogy(F, A001(:,1))
grid on
hold on
semilogy(F, A001(:,2))
semilogy(F, A001(:,3))
semilogy(F, A001(:,4))
semilogy(F, A001(:,5))
semilogy(F, A001(:,6))
xlabel('Частота, ГГц');
ylabel('Ослабление, дБ');
legend({'\theta = 5\circ','\theta = 8\circ','\theta = 11\circ','\theta = 14\circ','\theta = 17\circ','\theta = 20\circ'},'Interpreter','tex','Location','southeast','FontSize',12);


