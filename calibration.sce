[fd,SST,Sheetnames,Sheetpos] = xls_open('STTI2013_base.xls')

[Value,TextInd] = xls_read(fd,Sheetpos(1))

//acces a la i eme colonne: Value(i,:)
//acces a la i eme ligne: Value(,:i)

//cration matrice taux zeros coupons R(t,t+T)
R = Value'

mclose(fd)

//suppression des premières ligne et colone
R(1,:) = []
R(:,1) = []

nb_de_termes  = length(R(1,:))
nb_de_periodes = length((R(:,1)))

function [phit]=phi(t)
    phit = (1-exp(-t))/t
endfunction

function [psit]=psi(t)
    psit = phi(t) - exp(-t)
endfunction

function [erreur] = erreur_quadratique(a,b,R,tau)
    erreur = 0
    for t=1:nb_de_periodes
        for T=1:nb_de_termes
            R_NS_t_T = b(t)+a(t,1)*phi(T/tau)+a(t,2)*psi(T/tau)
            erreur = erreur + (R_NS_t_T - R(t,T))**2
        end
    end
endfunction

tau_init = 0.1
tau_final = 10
step = 0.5

function [tau_min, taux_courts] = recherche_tau_opt(R, tau_init, tau_final, step)
    
    nb_de_termes  = length(R(1,:))
    //create for each t, parameters mu_1, mu_2, mu_3
    tau_tab = []
    erreur_quadra_tab = []
    N_tau = int((tau_final - tau_init)/step)
    
    erreur_min = 0
        
    for i=0:N_tau
        tau = tau_init + i*step
        tau_tab = [tau_tab; tau]
        phi_tab = []
        psi_tab = []
        
        for T=1:nb_de_termes
            psiT = psi(T/tau)
            phiT = phi(T/tau)
            psi_tab = [psi_tab; psiT]
            phi_tab = [phi_tab; phiT]
        end
        //regression linéaire de R(t,T) T=1...nb_of_terms sur 1, psi_tab, phi_tab
        X = [ phi_tab' ; psi_tab']
        [a,b,sig] = reglin(X, R)
        erreur_quadra = erreur_quadratique(a,b,R,tau)
        erreur_quadra_tab = [erreur_quadra_tab; erreur_quadra] 
        
        if tau == tau_init then
            erreur_min = erreur_quadra
            tau_min = tau
        end
        
        if erreur_quadra < erreur_min then
            erreur_min = erreur_quadra
            tau_min = tau
            taux_courts = a(:,1) + b
        end
        
 
    end    
    //plot2d(tau_tab, erreur_quadra_tab)
    //xtitle('Calibration de tau_1', 'tau_1', 'Erreur globale')
    
endfunction


function plotNSvsReel(R, tau, t)
    phi_tab = []
    psi_tab = []
    for T=1:nb_de_termes
        psiT = psi(T/tau)
        phiT = phi(T/tau)
        psi_tab = [psi_tab; psiT]
        phi_tab = [phi_tab; phiT]
    end
    //regression linéaire de R(t,T) T=1...nb_of_terms sur 1, psi_tab, phi_tab
    X = [ phi_tab' ; psi_tab']
    [a,b,sig] = reglin(X, R)
    
    //taux zero-coupons de Nelson Siegel
    NSrates = a*X+b*ones(1,nb_de_termes)

    temps = linspace(1,nb_de_periodes, nb_de_periodes)'
    plot(temps,R(:,t), 'b', temps,  NSrates(:,t), 'g') 
    legends(['Taux réels' 'Taux Nelson Siegel'], [2,3],4)
    xtitle('Comparaison en t =' + string(t), 'temps', 'taux')
endfunction
//remarque: taux négatifs

//calibration du vasicek
function [r_infini, kappa, sigma_r]=calibration_vasicek(taux_courts)
    nb_de_periodes = length(taux_courts)
    
    [a,b,sigma] = reglin(taux_courts(1:nb_de_periodes-1)',taux_courts(2:nb_de_periodes)')
    
    if a>0 then
        kappa = -log(a)
        r_infini = b/(1-a)
        sigma_r = sigma*sqrt(-2*log(a)/(1-a**2))    
    else
        r_infini = 0
        kappa = 0
        sigma_r = 0
    end
    
    
endfunction

