clear all
close all
clc

tic
el = 5;
r0 = 6371; % радиус Земли, км
hstation = 0;
A0 = 95.571899;
A1 = -4.011801;
A2 = 6.424731*10^(-2);
A3 = -4.789660*10^(-4);
A4 = 1.340543*10^(-6);
Po0 = 7.5; %  плотность водяных паров на уровне Земли, г/м3
for layer = 1:922 % Характеристики слоев
    sigma(layer) = 0.0001*exp((layer-1)/100); % км
    if layer == 1
        heighR(layer) = r0;
        heighL(layer) = sigma(layer);
    else
        heighR(layer) = heighR(layer-1)+sigma(layer-1);
        heighL(layer) = heighL(layer-1)+sigma(layer);
    end
    heighL_pot(layer) = (heighL(layer)*6356.766)/(6356.766+heighL(layer));
    
    if      0<= heighL_pot(layer) &  heighL_pot(layer)<=11
        T(layer) = 288.15-6.5*heighL_pot(layer);
        P(layer) = 1013.25*(288.15/(288.15-6.5*heighL_pot(layer)))^(-34.1632/6.5);
    elseif  11< heighL_pot(layer) &  heighL_pot(layer)<=20
        T(layer) = 216.65;
        P(layer) = 226.3226*exp(-34.1632*(heighL_pot(layer)-11)/216.65);
    elseif  20< heighL_pot(layer) &  heighL_pot(layer)<=32
        T(layer) = 216.65+(heighL_pot(layer)-20);
        P(layer) = 54.7498*(216.65/(216.65+(heighL_pot(layer)-20)))^(34.1632);
    elseif  32< heighL_pot(layer) &  heighL_pot(layer)<=47
        T(layer) = 228.65+2.8*(heighL_pot(layer)-32);
        P(layer) = 8.680422*(228.65/(228.65+2.8*(heighL_pot(layer)-32)))^(34.1632/2.8);
    elseif  47< heighL_pot(layer) &  heighL_pot(layer)<=51
        T(layer) = 270.65;
        P(layer) = 1.109106*exp(-34.1632*(heighL_pot(layer)-47)/270.65);
    elseif  51< heighL_pot(layer) &  heighL_pot(layer)<=71
        T(layer) = 270.65-2.8*(heighL_pot(layer)-51);
        P(layer) = 0.6694167*(270.65/(270.65-2.8*(heighL_pot(layer)-51)))^(-34.1632/2.8);
    elseif  71< heighL_pot(layer) &  heighL_pot(layer)<=84.852
        T(layer) = 214.65-2.0*(heighL_pot(layer)-71);
        P(layer) = 0.03956649*(214.65/(214.65-2.0*(heighL_pot(layer)-71)))^(-34.1632/2.0);
    elseif  86<= heighL(layer) &  heighL(layer)<=91
        T(layer) = 186.8673;
    elseif  91< heighL(layer) &  heighL(layer)<=101
        T(layer) = 263.1905-76.3232*sqrt(1-((heighL(layer)-91)/19.9429)^2);
    end
    if 86<= heighL(layer) &  heighL(layer)<=101
        P(layer) = exp(A0+A1*heighL(layer)+A2*heighL(layer)^2+A3*heighL(layer)^3+A4*heighL(layer)^4);
    end
    Po(layer) = Po0*exp(-heighL(layer)/2);
    e(layer)=Po(layer)*T(layer)/216.7;
    if e(layer)/P(layer) <= 2*10^-6
        Po(layer)=Po(layer-1);
    end
    Nref(layer) = 77.6*(P(layer)+4810*e(layer)/T(layer))/T(layer);
    nref(layer) = 1+Nref(layer)*10^-6;
end
figure(1)
subplot(1,2,1);
plot(T,heighL);
grid on
xlabel('Температура, К')
ylabel('Высота, км')
subplot(1,2,2);
semilogx(P,heighL);
grid on
xlabel('Давление, гПа')
ylabel('Высота, км')
for layer = 1:922 % Геометрия слоев
    if layer == 1
        betaL(layer) = 90-el;
        New_betaL(layer) = betaL(layer);
    end
    Apas(layer) = -heighR(layer)*cosd(betaL(layer))+0.5*sqrt(4*heighR(layer)^2*cosd(betaL(layer))^2+8*heighR(layer)*sigma(layer)+4*sigma(layer)^2);
    alfaL(layer) = 180-acosd((-(Apas(layer))^2-2*heighR(layer)*sigma(layer)-(sigma(layer))^2)/(2*Apas(layer)*heighR(layer)+2*Apas(layer)*sigma(layer)));
    if layer <=921
        betaL(layer+1) = asind(nref(layer)*sind(alfaL(layer))/nref(layer+1));
    end
    if layer ==1
        gamma(layer) = betaL(layer)-alfaL(layer);
    else
        gamma(layer) = gamma(layer-1)+(betaL(layer)-alfaL(layer));
    end
end

for count = 1:923
    if count == 1
        goriz(count) = 0;
    elseif count == 2
        goriz(count) = tand(betaL(count-1))*heighL(count-1);
    else
        goriz(count) = tand(gamma(count-2)+betaL(count-1))*heighL(count-1);        
    end
    if count == 1
        vert(count) = 0;
    else
        vert(count) = heighL(count-1);
    end
end
tfs1 = 1.728 + 0.5411*el + 0.03723*el^2;
tfs2 = 0.1815 + 0.06272*el + 0.01380*el^2;
tfs3 =  0.01727 + 0.008288*el;
tfs = 1/(tfs1+hstation*tfs2+hstation^2*tfs3);
Te = el-tfs;
figure(2)
goriz1 =  heighL(922)/0.0874886;
Xx = [0, goriz1];
Yy = [0, heighL(922)];
goriz2 =  heighL(922)/tand(Te);
Xx1 = [0, goriz2];
plot(goriz,vert,'-','LineWidth',1.5)
hold on
plot(Xx,Yy,'--')
plot(Xx1,Yy,'-.','LineWidth',1.2)
grid on
grid minor
ylabel('Высота, км')
xlabel('Дальность пути по горизонтали, км')
legend('МСЭ-R','Прямолинейное распространение','Видимый угол места')
toc


