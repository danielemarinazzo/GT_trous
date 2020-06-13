clear;clc;
load('dati_scimmia_grezzi.mat')
J=4;m=5;
x=data(2000:3000,:);
type='poly';par=1;
[cb, co]=causality_scale(x,J,type,par,m);

