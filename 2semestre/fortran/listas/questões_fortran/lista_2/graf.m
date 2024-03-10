fi_N=load(laplace.dat)
fi_A=load(laplace_sol_analitica.dat)

figure
 surf(Fi_N)

figure
 surf(fi_A)
 zlim([0 1])
   
 x=1:50;
 L=25;
 
figure
plot(x,fi_A(L,:),x,fi_N(L,:),'o')  
