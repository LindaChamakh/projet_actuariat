//clear
stacksize max

//Input:
//Taux d'interet risque neutre par calibration
k_r = 0.2080239;
r_inf = 0.0402051;
r_0 = 0.005;
sigma_r  = 0.0018370;
//Processus x_t
k_x = 0.3;
x_inf = -0.01;
x_0 = 0.02;
sigma_x = 0.008;
//Log return of the asset
mu_a = 0.04;
sigma_a = 0.06;
//Correlation
rho_xr = 0.;
rho_as = 0.95;
rho_ar = 0.25;
//Maturité
T = 1;



//nb de pas dans l'EDP
N = 60;
M_1 = 15;
M_2 = 15;

//Résolution de l'EDP pour approcher phi
[phi,dx,dr,x,r] = edp(k_r,r_inf,r_0,sigma_r,k_x,x_inf,x_0,sigma_x,rho_xr,T,N,M_1,M_2);
//phi = Phi(:,N+1);



M=100000;//taille échantillon Monte Carlo 
A_0 = 101479200;//ACTIF: obligations, actions, immobilier à t=0
E_0 = 57238200;//PASSIF: dettes vis-à-vis des actionnaires
//L0 = ;//PASSIF: dettes vis-à-vis des assurés //INUTILE
P = 36000000;


//Deux gaussiennes indépendantes
Gr = rand(M,1,"normal");
Gx = rand(M,1,"normal");
Ga = rand(M,1,"normal");

//Modélisation des dynamiques r_1 et x_1
r_1 = r_0*exp(-k_r) + r_inf*(1-exp(-k_r)) + sigma_r*sqrt((1-exp(-2*k_r))/(2*k_r))*Gr;
x_1 = x_0*exp(-k_x) + x_inf*(1-exp(-k_x)) + sigma_x*sqrt((1-exp(-2*k_x))/(2*k_x))*((rho_as/sqrt(1-rho_ar*rho_ar))*Ga + sqrt((1-rho_ar**2- rho_as**2)/(1-rho_ar*rho_ar))*Gx);
R_1 = mu_a + rho_ar*sigma_a*Gr + sqrt(1-rho_ar**2)*sigma_a*Ga;
A_1 = A_0.*(1 + R_1);

//Calcul de BE
BE = zeros(M,1);
g_0 = g(find(x<=x_0 & x>(x_0 -dx)));
for i=1:M
    ind_r = find(r<r_1(i) & r>(r_1(i) -dr));
    ind_x = find(x<x_1(i) & x>(x_1(i) -dx));
    phi_1 = phi(ind_r*(M_1-1) + ind_x);
    BE(i) = P*exp(f(x_0,r_0) - g(x_0))*phi_1;
end

//Affichage de la moyenne et de l'écart type de BE
mprintf("Moyenne de BE = %f \n Ecart type de BE = %f \n",mean(BE),stdev(BE));

//Calcul de E_1
E_1 = A_1 - BE - P*g(x_0); 
mprintf("Moyenne de E1 = %f \n Ecart type de E1 = %f \n",mean(E_1),stdev(E_1));

mprintf("Moyenne de F1 = %f \n Ecart type de F1 = %f \n",mean(P*g(x_0)),stdev(P*g(x_0)));

//Calcul de L
L = exp(-r_0).*E_1;

//Statistique d'ordre de L (trier L)
L = gsort(L,'g','i');
mprintf("VaR = %f \n",L(ceil(M*0.005)));

//Calcul de SCR_0
SCR_0 = E_0 - L(ceil(M*0.005));
mprintf("SCR_0 = %f \n",SCR_0);

//Intervalle de confiance pour la VAR
N = linspace(1,M,M);
y = cdfbin("PQ",N,M*ones(1,M),(1-0.005)*ones(1,M),0.005*ones(1,M));

mprintf("borne inf et sup à environ 95: inf = %f,\t sup = %f\n",L(M-find(y>=0.9,1)),L(M-find(y>=0.05,1)));

plot(E_0 - L);
plot(SCR_0,'ro')
xlabel('N',"fontsize",6)
ylabel('$ E_0 - e^{-r_0}E_1 $',"fontsize",6)
