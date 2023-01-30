clear all
close all
clc

R0 = 6371;
Ri = R0 + 20;

c = @(h) 2*h;
fun = @(h) 1./(h.^2)+c(h);
q = integral(fun,R0,Ri);


