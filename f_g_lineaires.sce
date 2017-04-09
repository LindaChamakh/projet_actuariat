function [z]=zeta(k,t)
    z=(1-exp(-k*t))/k
endfunction

function [m, v, dm, dv]=Z_esp_var(x_0, x_inf, k, sigma, t)
    m=(x_0-x_inf)*zeta(k,t)+x_inf*t
    v=sigma**2*(t-zeta(k,t))/k**2 -sigma**2*zeta(k,t)**2/(2*k)
    dm = (x_0-x_inf)*exp(-k*t)+x_inf
    dv = sigma**2*(1-exp(-k*t))/k**2 - sigma**2*zeta(k,t)*exp(-k*t)/(2*k)
endfunction

function [theta1, theta2]=thetas(x_0, x_inf, k, sigma, t, eta)
    [m, v, dm, dv]=Z_esp_var(x_0, x_inf, k, sigma, t)
    theta1 = exp((1+eta)*m + (1+eta)**2 *v/2)
    theta2 = -eta*(dm+(1+eta)*dv/2)*theta1
endfunction

//discretisation: calcul de phi(t,T)
function [phi_t_T] = discrete_phi_t_T(t, T, N_pas, x_0, x_inf, k, sigma, eta, mu_i)
    phi_t_T = 0
    delta_t = int((T-t)/N_pas)
    t_i = t
    for i=0:(T-t)
        [theta1, theta2]=thetas(x_0, x_inf, k, sigma, t_i - t, eta)
        phi_t_T = phi_t_T + delta_t* exp(-mu_i*(t-t_i))*(mu_i*theta1+theta2)
        t_i = t_i + delta_t
    end
    
    [theta1, theta2]=thetas(x_0, x_inf, k, sigma, T - t, eta)

    phi_t_T = phi_t_T + exp(-mu_i*(T-t))*theta1
endfunction

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
T = 10 ;
//Nb de pas de temps
N = 15;


function [f]=f(x,r)
    f = x+r
endfunction

function [g]=g(x)
    g = mu_i - eta*x
endfunction

M=10000;//taille échantillon Monte Carlo 
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

A_1 = A_0*(1 + R_1);

PM_1 = P*exp(f(x_0,r_0) - g(x_0))

t = 1
BE_1 = PM_1*discrete_phi_t_T(t, T, N, x_0, x_inf, k_x, sigma_x, eta, mu_i)
test = discrete_phi_t_T(t, T, N, x_0, x_inf, k_x, sigma_x, eta, mu_i)

BE = BE_1*ones(M,1)
//Calcul de E_1
F_1 = P*g(x_0)
E_1 = A_1 - BE - F_1; 

//Calcul de L
L = exp(-r_0).*E_1;

//Statistique d'ordre de L (trier L)
L = gsort(L,'g','i');
mprintf("VaR = %f \n",L(ceil(M*0.005)));

//Calcul de SCR_0
SCR_0 = E_0 - L(ceil(M*0.005));
mprintf("SCR_0 = %f \n",SCR_0);

Nb_simul = linspace(1,M,M)'
plot(Nb_simul, E_0-L)
//hold on
plot(M*0.005, SCR_0, 'r*')
xlabel('N',"fontsize",6)
ylabel('$ E_0 - e^{-r_0}E_1 $',"fontsize",6)
