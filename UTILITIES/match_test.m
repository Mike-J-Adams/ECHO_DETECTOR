%match_test.m

clear 
close all

x = ones(10,1);
b = x(end:-1:1);
y = filter(b,1,x);