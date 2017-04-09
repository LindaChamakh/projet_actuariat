# projet_actuariat

Contrat d'épargne en euro avec participation aux bénéfices
On cherche à mesurer l'impact sur les besoins en fonds propres de la possibilité de rachat du contrat.

Rapport:

https://www.overleaf.com/8779860zwsckxkztfyh#/31341214/

Fonctionnement du code:

- calibration.sce: calibration du modèle de Vasicek à partir du fichier excel des zéro-coupons STTI2013_base.xls

Résolution avec taux linéaires et formules fermées:
- f_g_lineaires.sce: contient les fonctions intermédiaires (formules fermées) et d'exécution (calcul des éléments du SCR) dans le cas de taux linéaires

Résolution par EDP
- functions.sci: contient les fonctions des taux f et g et la résolution  de l'EDP par  différences finies.
- script.sce: script utilisant  les fonctions du fichier functions.sci pour calculer le SCR à t = 0 avec intervalle de confiance pour la VaR.

Code annexe:
- script-test... : création de graphes de sensibilité par différences finies.

