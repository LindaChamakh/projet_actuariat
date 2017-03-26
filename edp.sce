clear all
//Input:
//Taux d'interet risque neutre
k_r = 0.18;
r_inf = 0.03;
r_0 = 0.005;
sigma_r = 0.02;
//Processus x_t
k_x = 0.3;
x_inf = -0.01;
x_0 = 0.02;
sigma_x = 0.008;
//Correlation
rho = 0.5;
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
T = 10;

//Nb de pas de temps
N = 15;
dt = T/N;
//Nb de pas en x
L_1 = 6*sigma_x/sqrt(2*k_x);
M_1 = 20;
dx = L_1/M_1;
//Nb de pas en r
L_2 = 6*sigma_r/sqrt(2*k_r);
M_2 = 20;
dr = L_2/M_2;

x = linspace(x_inf - L_1/2,x_inf + L_1/2,M_1);
r = linspace(r_inf - L_2/2,r_inf + L_2/2,M_2);

//foncton f et g
deff('[x]=f(x,r,TMG)','x = max(x+r,TMG)');
g = mu_i + eta*x';

//vecteur contenant la solution au temps actuel et au temps précédent
phi_prec = ones(M_1*M_2,1);//cond initiale
phi = zeros(M_1*M_2,1);

//Remplissage des matrices du système linéaire
A = zeros(M_1*M_2,M_1*M_2);
B = zeros(M_1*M_2,M_1*M_2);
i = 1;
j = 1;
for l=1:M_1*M_2
    //disp(i)
    //disp(j)
    A(l,l) = 1/dt + k_r*(r_inf-r(j))*0.5/dr + k_x*(x_inf-x(i))*0.5/dx + sigma_r*sigma_r*0.5/(dr*dr) + sigma_x*sigma_x*0.5/(dx*dx) + rho*sigma_x*sigma_r/dx/dr - 0.5*(f(x(i),r(j),TMG)-r(j)-g(i));
    B(l,l) = 1/dt - k_r*(r_inf-r(j))*0.5/dr - k_x*(x_inf-x(i))*0.5/dx - sigma_r*sigma_r*0.5/(dr*dr) - sigma_x*sigma_x*0.5/(dx*dx) - rho*sigma_x*sigma_r/dx/dr + 0.5*(f(x(i),r(j),TMG)-r(j)-g(i));
    if j > 1
        //disp(l)
        A(l,l-M_1) = - sigma_r*sigma_r*0.25/dr/dr;
        B(l,l-M_1) = sigma_r*sigma_r*0.25/dr/dr;
    end
    if i > 1
        A(l,l-1) = - sigma_x*sigma_x*0.25/dx/dx;
        B(l,l-1) = sigma_x*sigma_x*0.25/dx/dx;
    end
    if  i < M_1
        A(l,l+1) = - sigma_x*sigma_x*0.25/dx/dx -k_x*(x_inf-x(i))*0.5/dx - rho*sigma_x*sigma_r/dx/dr;
        B(l,l+1) = sigma_x*sigma_x*0.25/dx/dx + k_x*(x_inf-x(i))*0.5/dx + rho*sigma_x*sigma_r/dx/dr;
    end
    if j < M_2
        A(l,l+M_1) = - sigma_r*sigma_r*0.25/dr/dr -k_r*(r_inf-r(j))*0.5/dr - rho*sigma_x*sigma_r/dx/dr;
        B(l,l+M_1) =  sigma_r*sigma_r*0.25/dr/dr + k_r*(r_inf-r(j))*0.5/dr + rho*sigma_x*sigma_r/dx/dr;
    end
    if j < M_2 & i < M_1
        A(l,l+M_1+1) = rho*sigma_x*sigma_r/dx/dr;
        B(l,l+M_1+1) = -rho*sigma_x*sigma_r/dx/dr;
    end
    i=i+1;
    if l == j*(M_1)
        j = j+1;
        i = 1;
    end
end

//Approximation de l'EDP par différences finie en utilisant le schéma de crank-nicolson
// A*phi^{n+1} = B*phi^n + G
G = g;
for i=1:M_2-1
    G = [G;g];
end

for t=1:N
    b = B*phi_prec + G;
    phi = A\b;
    phi_prec = phi;
end


