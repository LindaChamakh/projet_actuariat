clear all
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
//Log return of the asser
mu_a = 0.04;
sigma_a = 0.06;
//Correlation
rho_xr = 0.;
rho_as = 0.95;
rho_ar = 0.25;
//Parametres taux de rachat
alpha = -0.05;
beta_ = -0.01;
gamma_ = 0.01;
delta = 0.03;
mu_min = -0.05;
mu_max = 0.3;
//Rachat structurel
mu_i = 0.05;
eta = 2;
//Taux minimum garanti
TMG = 0.015;
//Maturité
T = 10 - 1;

//Nb de pas de temps
N = 15;
dt = T/N;
//Nb de pas en x
L_1 = 45*sigma_x*sqrt((1-exp(-2*k_x))/2*k_x);
a_x = exp(-k_x)*(x_0 - x_inf) + x_inf - L_1/2;
b_x = exp(-k_x)*(x_0 - x_inf) + x_inf + L_1/2;
M_1 = 20;
dx = L_1/M_1;
//Nb de pas en r
L_2 = 45*sigma_r*sqrt((1-exp(-2*k_r))/2*k_r);
a_r = exp(-k_r)*(r_0 - r_inf) + r_inf - L_2/2;
b_r = exp(-k_r)*(r_0 - r_inf) + r_inf + L_2/2;
M_2 = 20;
dr = L_2/M_2;

//discrétisation de l'espace
x = linspace(a_x + dx,b_x - dx,M_1-1);
r = linspace(a_r + dr,b_r - dr,M_2-1);

//foncton f et g
deff('[x]=f(x,r,TMG)','x = max(x+r,TMG)');
g = mu_i + eta*max((-x)',0);

//vecteur contenant la solution au temps actuel et au temps précédent
phi_prec = ones((M_1-1)*(M_2-1),1);//cond initiale
phi = zeros((M_1-1)*(M_2-1),1);

//Remplissage des matrices du système linéaire
A = zeros((M_1-1)*(M_2-1),(M_1-1)*(M_2-1));
B = zeros((M_1-1)*(M_2-1),(M_1-1)*(M_2-1));
i = 1;
j = 1;
for l=1:(M_1-1)*(M_2-1)
    //disp(i)
    //disp(j)
    A(l,l) = 1/dt + k_r*(r_inf-r(j))*0.5/dr + k_x*(x_inf-x(i))*0.5/dx + sigma_r*sigma_r*0.5/(dr*dr) + sigma_x*sigma_x*0.5/(dx*dx) + rho_xr*sigma_x*sigma_r/dx/dr - 0.5*(f(x(i),r(j),TMG)-r(j)-g(i));
    B(l,l) = 1/dt - k_r*(r_inf-r(j))*0.5/dr - k_x*(x_inf-x(i))*0.5/dx - sigma_r*sigma_r*0.5/(dr*dr) - sigma_x*sigma_x*0.5/(dx*dx) - rho_xr*sigma_x*sigma_r/dx/dr + 0.5*(f(x(i),r(j),TMG)-r(j)-g(i));
    if j > 1
        //disp(l)
        A(l,l-M_1+1) = - sigma_r*sigma_r*0.25/dr/dr;
        B(l,l-M_1+1) = sigma_r*sigma_r*0.25/dr/dr;
    end
    if i > 1
        A(l,l-1) = - sigma_x*sigma_x*0.25/dx/dx;
        B(l,l-1) = sigma_x*sigma_x*0.25/dx/dx;
    end
    if  i < M_1-1
        A(l,l+1) = - sigma_x*sigma_x*0.25/dx/dx -k_x*(x_inf-x(i))*0.5/dx - rho_xr*sigma_x*sigma_r/dx/dr;
        B(l,l+1) = sigma_x*sigma_x*0.25/dx/dx + k_x*(x_inf-x(i))*0.5/dx + rho_xr*sigma_x*sigma_r/dx/dr;
    end
    if j < M_2-1
        A(l,l+M_1-1) = - sigma_r*sigma_r*0.25/dr/dr -k_r*(r_inf-r(j))*0.5/dr - rho_xr*sigma_x*sigma_r/dx/dr;
        B(l,l+M_1-1) =  sigma_r*sigma_r*0.25/dr/dr + k_r*(r_inf-r(j))*0.5/dr + rho_xr*sigma_x*sigma_r/dx/dr;
    end
    if j < M_2-1 & i < M_1-1
        A(l,l+M_1) = rho_xr*sigma_x*sigma_r/dx/dr;
        B(l,l+M_1) = -rho_xr*sigma_x*sigma_r/dx/dr;
    end
    i=i+1;
    if l == j*(M_1-1)
        j = j+1;
        i = 1;
    end
end

//Approximation de l'EDP par différences finie en utilisant le schéma de crank-nicolson
// A*phi^{n+1} = B*phi^n + G
G = g;
for i=1:M_2-2
    G = [G;g];
end

for t=1:N
    b = B*phi_prec + G;
    phi = A\b;
    phi_prec = phi;
end



M=10000;//taille échantillon Monte Carlo 
A_0 = 101479200;//ACTIF: obligations, actions, immobilier à t=0
E_0 = 57238200;//PASSIF: dettes vis-à-vis des actionnaires
//L0 = ;//PASSIF: dettes vis-à-vis des assurés //INUTILE
P = 100000000;


//Deux gaussiennes indépendantes
Gr = rand(M,1,"normal");
Gx = rand(M,1,"normal");
Ga = rand(M,1,"normal");

//Modélisation des dynamiques r_1 et x_1
r_1 = r_0*exp(-k_r) + r_inf*(1-exp(-k_r)) + sigma_r*sqrt((1-exp(-2*k_r))/(2*k_r))*Gr;
x_1 = x_0*exp(-k_x) + x_inf*(1-exp(-k_x)) + sigma_x*sqrt((1-exp(-2*k_x))/(2*k_x))*((rho_as/sqrt(1-rho_ar*rho_ar))*Ga + sqrt((1-rho_ar**2- rho_as**2)/(1-rho_ar*rho_ar))*Gx);
R_1 = mu_a + rho_ar*sigma_a*Gr + sqrt(1-rho_ar**2)*sigma_a*Ga;
A_1 = A_0*(1 + R_1);

//Calcul de BE
BE = zeros(M,1);
g_0 = g(find(x<=x_0 & x>(x_0 -dx)));
for i=1:M
    ind_r = find(r<=r_1(i) & r>(r_1(i) -dr));
    ind_x = find(x<=x_1(i) & x>(x_1(i) -dx));
    phi_1 = phi(ind_r*(M_1-1) + ind_x);
    BE(i) = P*exp(f(x_0,r_0) - g_0)*phi_1;
    if (BE(i) == 0)
        disp(ind_r);
        disp(ind_x);
        disp(phi_1);
    end
end

//Affichage de la moyenne et de l'écart type de BE
mprintf("Moyenne = %f \n Ecart type = %f \n",mean(BE),stdev(BE));

//Calcul de E_1
E_1 = A_1 - BE - P*g_0; 

//Calcul de L
L = exp(-r_0).*E_1;

//Statistique d'ordre de L (trier L)
L = gsort(L,'g','i');
mprintf("VaR = %f \n",L(ceil(N*0.95)));

//Calcul de SCR_0
SCR_0 = E_0 - L(ceil(N*0.95));
mprintf("SCR_0 = %f \n",SCR_0);

 
