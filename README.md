# TMDSpatial
 
L'idée est ici de pouvoir simuler des jeux de données proches de ceux utilisés en écologie des commuanutés, et éprouver les différentes méthodes disponibles face ces données.  
On réutilise le package {mcomsimr} (Thompson et al., 2020, EcoLet., voir https://github.com/plthompson/mcomsimr) mais légérement modifié. Les modifications actuelles sont : 
+ __Autocorrelation temporelle__ (au niveau des patches): Celle-ci est initialement inclue dans le package {ecospat}, i.e., il existe une variabilité temporelle, plus ou moins autocorrélée, dans les conditions environnementales des patches. Cette variabilité ajoute un niveau de complexité supplémentaire. Ainsi, celle-ci n'est pas considérée dans cette version du package. Voir l'argument __temporal_autocorr__ de la fonction simulate_MC() pour l'activer. A noter; le même niveau d'autocorrélation est considéré à la fois pour le temps et pour l'espace.
+ __Extinction de patche__: Voir les paramètres __extrip_prob__ et __extirp_all_pop__ de la fonction simulate_MC(). Le premier argument est un argument de la fonction de base, il donne la probabilité qu'à chaque pas de temps, une population s'éteigne dans un patch. Le second argument est lui en revanche nouveau. Il permet de spécifier si l'extinction concerne __toutes__ les espèces présentes dans un patche ou si ce tirage est réalisé pour chacune des espèces. A noter; si le probabilité est de 0.01 et que le patche contient 3 espèces, dans le cas où __extirp_all_pop = TRUE__, les trois espèces peuvent s'éteindre simultanément avec une prob. de 0.01 mais les extinctions individuelles ne sont pas possible via ce paramètre d'extirpation (mais possible par la dynamique locale). En revanche, si __extirp_all_pop = FALSE__, alors la probabilité pour les trois espèces de s'éteindre simultanément est ici : $0.01^3=10^{-6}$.


Il y a également quelques fonctions en plus pour faciliter l'utilisation des sorties de la fonction simulate_MC().
+ sim_to_matrix()
+ abund_to_occ()
+ plots_occupancies()



## Simulations

Les dynamiques locales sont gérées par : $equation$.

Pour les simulations, on peut jouer avec les paramètres suivants: 

### Environnement et niches (fund.)

L'environnement est constitué uniquement d'un axe. 

Niche : auto ou manuel (optima_spacing, env_niche_breadth, min_env, max_env, env_optima, env1Scale)

+ Local dynamique : max_r

### Landscape and dispersal

landscape (data.Frame), torus, dispersal, disp_mat, kernet_exp

### Competition



## Trucs testables

+ 

### Références

Thompson, Patrick L., Laura Melissa Guzman, Luc De Meester, Zsófia Horváth, Robert Ptacnik, Bram Vanschoenwinkel, Duarte S. Viana, et Jonathan M. Chase. «A Process‐based Metacommunity Framework Linking Local and Regional Scale Community Ecology ». Ecology Letters 23, nᵒ 9 (septembre 2020): 1314‑29. https://doi.org/10.1111/ele.13568.

