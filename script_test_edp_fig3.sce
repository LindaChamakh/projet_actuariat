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
T = 10;

//nb de pas dans l'EDP
N = 180;
M_1 = 15;
M_2 = 15;

//Résolution de l'EDP pour approcher phi
N_x = 15;//nb de points pour k_x
N_sigma = 15//nb de points pour sigma_x
abscisse = zeros(N_sigma,1);
ordonnee = zeros(N_x,1);
resultat = zeros(N_sigma,N_x);
dk = 0.3/(N_x-1);
dsigma = 0.008/(N_sigma-1);
tic();
for i=1:N_x
    k_x = 0.01 + (i-1)*dk;
    ordonnee(i) = k_x;
    for j=1:N_sigma
        sigma_x = 0.0045 + (j-1)*dsigma;
        abscisse(j) = sigma_x;
        [phi,dx,dr,x,r] = edp(k_r,r_inf,r_0,sigma_r,k_x,x_inf,x_0,sigma_x,rho_xr,T,N,M_1,M_2);
        //phi = Phi(:,N+1);
        ind_r = find((r<=r_0) & r>(r_0 -dr));
        ind_x = find(x<=x_0 & x>(x_0 -dx));
        if ind_x == M_1 - 1
            resultat(j,i) = (phi(ind_r*(M_1-1) + ind_x) - phi(ind_r*(M_1-1) + ind_x - 1))/(-dx);
        else
            resultat(j,i) = (phi(ind_r*(M_1-1) + ind_x +1) - phi(ind_r*(M_1-1) + ind_x))/(dx);
        end
    end
end
printf('\ntemps nécessaire  : %.3f\n',toc());

f=gcf();
f.color_map = hotcolormap(128);
surf(abscisse,ordonnee,resultat)
a = gca();
a.rotation_angles=[60,310];
xtitle('$\partial_{x} \psi\ à\ t=0\ en\ (x_0,r_0)$')
xlabel("$\sigma_x$","fontsize",6);
ylabel("$\kappa_x$","fontsize",6);
zlabel("$\partial_{x} \psi\$","fontsize",6);
