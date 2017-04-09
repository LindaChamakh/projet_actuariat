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

//sensibilité pour la différence finie
epsilon = 0.00001

//nb de pas dans l'EDP
N = 50;
M_1 = 30;
M_2 = 30;

//Résolution de l'EDP pour approcher phi
N_x_0 = 20;//nb de points pour x_0
N_r_0 = 20//nb de points pour r_0
abscisse = zeros(N_x_0,1);
ordonnee = zeros(N_r_0,1);
resultat = zeros(N_x_0,N_r_0);
dx_0 = 0.02/(N_x_0-1);
dr_0 = 0.007/(N_r_0-1);
for i=1:N_r_0
    r_0 = 0.001 + (i-1)*dr_0;
    ordonnee(i) = r_0;
    for j=1:N_x_0
        x_0 = -0.01 + (j-1)*dx_0;
        abscisse(j) = x_0;
        [Phi,dx,dr,x,r] = edp(k_r,r_inf,r_0,sigma_r,k_x,x_inf,x_0 + epsilon,sigma_x,rho_xr,T,N,M_1,M_2);
        phi = Phi(:,N+1);
        ind_r = find((r<=r_0) & r>(r_0 -dr));
        ind_x = find(x<=x_0 & x>(x_0 -dx));
        resultat(j,i) = (phi(ind_r*(M_1-1) + ind_x +1) - phi(ind_r*(M_1-1) + ind_x))/(x(ind_x+1)-x(ind_x));
    end
end


f=gcf();
f.color_map = hotcolormap(128);
surf(abscisse,ordonnee,resultat)
a = gca();
a.rotation_angles=[60,310];
xtitle('$\partial_{x} \psi\ à\ t=0\ en\ (x_0,r_0)$')
xlabel("$x_0$","fontsize",3);
ylabel("$r_0$","fontsize",3);
zlabel("$\partial_{x} \psi\$","fontsize",3);
