clear all;
clear all;

dd = load('sond.out');
t  = dd(:,1);
r  = dd(:,2);
f  = dd(:,3);

figure(1)
loglog(t,r,'b')
grid

figure(2)
loglog(t,f,'ro')
grid