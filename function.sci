clear
function [y] = f(x,r)
    //Taux minimum garanti
    TMG = 0.025;
    y = max(x+r,TMG);
//    y = x+r;
endfunction

function [y] = g(x)
    //Rachat structurel
    mu_i = 0.05;
    eta = 2;
    //y = min(1,max(0,mu_c(x')+mu_i));
    //y = mu_i - eta*x';
    y = mu_i + min(x',0);
endfunction

function [y] = mu_c(x)
    //Parametres taux de rachat
    alpha = -0.05;
    beta_ = -0.01;
    gamma_ = 0.01;
    delta = 0.03;
    mu_min = -0.05;
    mu_max = 0.3;
    y(find(x<=alpha)) = mu_max;
    y(find(x>=alpha & x <= beta_)) = mu_max*((x(find(x>=alpha & x <= beta_))-beta_)/(alpha-beta_));
    y(find(x>=beta_ & x <= gamma_)) = 0;
    y(find(x>=gamma_ & x <= delta)) = mu_min*((x(find(x>=delta & x <= gamma_))-gamma_)/(delta-gamma_));
    y(find(x>=delta)) = mu_min;
endfunction

function [phi,dx,dr,x,r] = edp(k_r,r_inf,r_0,sigma_r,k_x,x_inf,x_0,sigma_x,rho_xr,T,N,M_1,M_2)
    //Nb de pas en x
    L_1 = 45*sigma_x*sqrt((1-exp(-2*k_x))/2*k_x);
    a_x = exp(-k_x)*(x_0 - x_inf) + x_inf - L_1/2;
    b_x = exp(-k_x)*(x_0 - x_inf) + x_inf + L_1/2;
//    L_1 = 10*sigma_x/sqrt(2*k_x);
//    a_x = x_inf - L_1/2;
//    b_x = x_inf + L_1/2;
    dx = L_1/M_1;
    //Nb de pas en r
    L_2 = 45*sigma_r*sqrt((1-exp(-2*k_r))/2*k_r);
    a_r = exp(-k_r)*(r_0 - r_inf) + r_inf - L_2/2;
    b_r = exp(-k_r)*(r_0 - r_inf) + r_inf + L_2/2;
//    L_2 = 10*sigma_r*sqrt(2*k_r);
//    a_r = r_inf - L_2/2;
//    b_r = r_inf + L_2/2;
    dr = L_2/M_2;

    //discrétisation de l'espace
    x = linspace(a_x + dx,b_x - dx,M_1-1);
    r = linspace(a_r + dr,b_r - dr,M_2-1);
    
    //Nb de pas de temps avec condition CFL
    //dt = T/N;
//    printf("ici : %f\n",- k_r*(r_inf-max(r))*0.5/dr - k_x*(x_inf-max(x))*0.5/dx - sigma_r*sigma_r*0.5/(dr*dr) - sigma_x*sigma_x*0.5/(dx*dx) - 0.5*rho_xr*sigma_x*sigma_r/dx/dr + 0.5*(max(f(x,r))-max(r)-max(g(x))));
    dt = 1/(- k_r*(r_inf-max(r))*0.5/dr - k_x*(x_inf-max(x))*0.5/dx - sigma_r*sigma_r*0.5/(dr*dr) - sigma_x*sigma_x*0.5/(dx*dx) - 0.5*rho_xr*sigma_x*sigma_r/dx/dr + 0.5*(max(f(x,r))-max(r)-max(g(x))));
    if dt <= 0 
        printf("dt négatif : %f\n",dt);
    end
    dt = abs(dt);/// - 0.000001;
    //vecteur contenant la solution à tous les temps
    phi_prec = ones((M_1-1)*(M_2-1),1);//cond initiale
    phi = ones((M_1-1)*(M_2-1),1);
    //phi = ones((M_1-1)*(M_2-1),N+1);

    //Remplissage des matrices du système linéaire
    //A = zeros((M_1-1)*(M_2-1),(M_1-1)*(M_2-1));
    //B = zeros((M_1-1)*(M_2-1),(M_1-1)*(M_2-1));
    A = sparse ([], [], [(M_1-1)*(M_2-1),(M_1-1)*(M_2-1)]);
    B = sparse ([], [], [(M_1-1)*(M_2-1),(M_1-1)*(M_2-1)]);
    i = 1;
    j = 1;
    for l=1:(M_1-1)*(M_2-1)
        if j > 1
            A(l,l-M_1+1) = - sigma_r*sigma_r*0.25/dr/dr;
            B(l,l-M_1+1) = sigma_r*sigma_r*0.25/dr/dr;
        end
        if i > 1
            A(l,l-1) = - sigma_x*sigma_x*0.25/dx/dx;
            B(l,l-1) = sigma_x*sigma_x*0.25/dx/dx;
        end
        A(l,l) = 1/dt + k_r*(r_inf-r(j))*0.5/dr + k_x*(x_inf-x(i))*0.5/dx + sigma_r*sigma_r*0.5/(dr*dr) + sigma_x*sigma_x*0.5/(dx*dx) + 0.5*rho_xr*sigma_x*sigma_r/dx/dr - 0.5*(f(x(i),r(j))-r(j)-g(x(i)));
        B(l,l) = 1/dt - k_r*(r_inf-r(j))*0.5/dr - k_x*(x_inf-x(i))*0.5/dx - sigma_r*sigma_r*0.5/(dr*dr) - sigma_x*sigma_x*0.5/(dx*dx) - 0.5*rho_xr*sigma_x*sigma_r/dx/dr + 0.5*(f(x(i),r(j))-r(j)-g(x(i)));
        if  i < M_1-1
            A(l,l+1) = - sigma_x*sigma_x*0.25/dx/dx -k_x*(x_inf-x(i))*0.5/dx - rho_xr*sigma_x*sigma_r/dx/dr;
            B(l,l+1) = sigma_x*sigma_x*0.25/dx/dx + k_x*(x_inf-x(i))*0.5/dx + rho_xr*sigma_x*sigma_r/dx/dr;
        end
        if j < M_2-1
            A(l,l+M_1-1) = - sigma_r*sigma_r*0.25/dr/dr -k_r*(r_inf-r(j))*0.5/dr - rho_xr*sigma_x*sigma_r/dx/dr;
            B(l,l+M_1-1) =  sigma_r*sigma_r*0.25/dr/dr + k_r*(r_inf-r(j))*0.5/dr + rho_xr*sigma_x*sigma_r/dx/dr;
        end
        if j < M_2-1 & i < M_1-1
            A(l,l+M_1) = 0.5*rho_xr*sigma_x*sigma_r/dx/dr;
            B(l,l+M_1) = -0.5*rho_xr*sigma_x*sigma_r/dx/dr;
        end
        i=i+1;
        if l == j*(M_1-1)
            j = j+1;
            i = 1;
        end
    end
    //A=sparse(A);
    //B=sparse(B);

    //Approximation de l'EDP par différences finie en utilisant le schéma de crank-nicolson
    // A*phi^{n+1} = B*phi^n + G
    G = g(x);
    for i=1:M_2-2
        G = [G;g(x)];
    end
    t = dt;
    while t < T
        //b = B*phi(:,t) + G;
        //tic();
        //phi(:,t+1) = umfpack(A,'\',b);
        //printf('\ntemps nécessaire à la résolution du système : %.3f\n tn = %.3f',toc(),t*T/N);
        b = B*phi_prec + G;
        phi = umfpack(A,'\',b);
        phi_prec = phi;
        t =t +dt;
    end
endfunction
