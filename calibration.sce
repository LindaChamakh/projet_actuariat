[fd,SST,Sheetnames,Sheetpos] = xls_open('/Users/lindachamakh/Documents/3A_MAF/actuariat/documents/STTI2013_base.xls')
[Value,TextInd] = xls_read(fd,Sheetpos(1))

//acces a la i eme colonne: Value(i,:)
//acces a la i eme ligne: Value(,:i)

//cration matrice taux zeros coupons R(t,t+T)
R = Value'

mclose(fd)

//delete first row and first column
R(1,:) = []
R(:,1) = []

nb_of_terms  = length(R(1,:))
nb_of_periods = length((R(:,1)))

function [phit]=phi(t)
    phit = (1-exp(-t))/t
endfunction

function [psit]=psi(t)
    psit = phi(t) - exp(-t)
endfunction

function [erreur] = erreur_quadratique(a,b,R,tau)
    erreur = 0
    for t=1:nb_of_periods
        for T=1:nb_of_terms
            R_NS_t_T = b(t)+a(t,1)*phi(T/tau)+a(t,2)*psi(T/tau)
            erreur = erreur + (R_NS_t_T - R(t,T))**2
        end
    end
endfunction

//create for each t, parameters mu_1, mu_2, mu_3
tau_tab = []
erreur_quadra_tab = []
tau_init = 0.1
tau_final = 10
step = 0.5
N_tau = int((tau_final - tau_init)/step)

for i=0:N_tau
    tau = tau_init + i*step
    tau_tab = [tau_tab; tau]
    phi_tab = []
    psi_tab = []
    
    for T=1:nb_of_terms
        [psiT] = psi(T/tau)
        [phiT] = phi(T/tau)
        psi_tab = [psi_tab; psiT]
        phi_tab = [phi_tab; phiT]
    end
    //regression linéaire de R(t,T) T=1...nb_of_terms sur 1, psi_tab, phi_tab
    X = [ phi_tab' ; psi_tab']
    [a,b,sig] = reglin(X, R)
    erreur_quadra = erreur_quadratique(a,b,R,tau)
    erreur_quadra_tab = [erreur_quadra_tab; erreur_quadra] 
end

plot(tau_tab, erreur_quadra_tab)



