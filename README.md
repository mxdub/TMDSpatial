# TMDSpatial
 
L'idée est ici de pouvoir simuler des jeux de données proches de ceux utilisés en écologie des commuanutés, et éprouver les différentes méthodes disponibles face à ces données.  
On réutilise le package {mcomsimr} (Thompson et al., 2020, EcoLet., obtenu depuis https://github.com/plthompson/mcomsimr) mais légérement modifié. Les modifications actuelles sont : 
+ __Autocorrelation temporelle__ (au niveau des patchs): Celle-ci est initialement inclue dans le package {ecospat}, i.e., il existe une variabilité temporelle, plus ou moins autocorrélée, dans les conditions environnementales des patchs. Cette variabilité ajoute un niveau de complexité supplémentaire. Ainsi, celle-ci n'est pas considérée dans cette version du package. Voir l'argument __temporal_autocorr__ de la fonction simulate_MC() pour l'activer. A noter; le même niveau d'autocorrélation est considéré à la fois pour le temps et pour l'espace.
+ __Extinction de patche__: Voir les paramètres __extrip_prob__ et __extirp_all_pop__ de la fonction simulate_MC(). Le premier argument est un argument de la fonction de base, il donne la probabilité qu'à chaque pas de temps, une population s'éteigne dans un patch. Le second argument est lui en revanche nouveau. Il permet de spécifier si l'extinction concerne __toutes__ les espèces présentes dans un patch ou si ce tirage est réalisé pour chacune des espèces. A noter; si la probabilité est de 0.01 et que le patch contient trois espèces, dans le cas où __extirp_all_pop = TRUE__, les trois espèces peuvent s'éteindre simultanément avec une prob. de 0.01 mais les extinctions individuelles ne sont pas possibles via ce paramètre d'extirpation (mais possible par la dynamique locale). En revanche, si __extirp_all_pop = FALSE__, alors la probabilité pour les trois espèces de s'éteindre simultanément est ici : $0.01^3=10^{-6}$.


Il y a également quelques fonctions en plus pour faciliter l'utilisation des sorties de la fonction simulate_MC().
+ sim_to_matrix() : produit une matrice d'abondances à trois dimensions (espèces x sites x temps) depuis la sortie de simulate_MC()
+ abund_to_occ() : produit une matrice d'occurrences (3D) à partir d'une matrice d'abondances
+ plots_occupancies() : produit un graphique du taux d'occupation des patchs à partir d'une matrice d'occurrences. 
+ plots_envt() : produit un graphique avec la répartion des patchs dans l'espace et leur valeur environnementale

## Simulations

Les dynamiques locales suivent l'équation (1) de Thompson et al. (2020). A noter qu'il existe de la stochasticité démographique (i.e., abondances sont des tirages dans une loi de Poisson).

Pour les simulations, les paramètres suivants sont disponibles : 

### Environnement et niches (fund.)

Un unique axe de variation environnementale existe (valeurs de 0 à 1). Comme précisé plus haut, la "valeur environnementale" de chaque site est fixe dans le temps par défaut. L'autocorrelation spatiale dépend du paramètre __env1Scale__ (entre 0 et 1000). 

La niche fondamentale des espèces peut être définie automatiquement (si __env_optima__ n'est pas spécifié) ou manuellement. Si les niches sont définies automatiquement, l'argument __optima_spacing__ définit la manière dont les optimums sont tirés (soit aléatoirement, soit de manière uniforme - *random* vs. *even*). Sinon, un vecteur de la taille du nombre d'espèces donnent les valeurs des optimums pour chacune des espèces. Trois arguments supplémentaires peuvent être spécifiés : __max_env__ et __min_env__ donnent les valeurs minimales et maximales pour les optimums (lorsque tirés aléatoirement), et __env_niche_breadth__ donne l'étendue de la tolérance environnementale (peut être soit une valeur pour toutes les espèces, ou un vecteur de la taille du nombre d'espèces). Le taux de croissance effectif au sein d'un patch est une fonction Gaussienne du gradient envionnemental (et indépendant de la densité locale). Le taux de croissance maximal étant défini par l'argument __max_r__. (voir Eq. 2 dans Thompson et al. 2020).

### Landscape and dispersal

Par défaut, les positions des patchs sont tirées aléatoirement (positions x et y entre 1 et 100), alternativement, il est possible de fournir un dataframe avec les coordonnées x/y des patches en colonnes via l'argument __landscape__. Par defaut, le paysage n'est pas torique cela peut néanmoins être modifié via l'argument __torus__ (passé à *TRUE*).

La dispersion est gérée via trois arguments. L'argument __dispersal__ donne la probabilité que chaque individu disperse à chaque pas de temps. Les probabilités de migration d'un patch à un autre sont gérées soit via (i) l'argument __kernel_exp__, dans ce cas, une matrice de probabilité de dispersion est générée entre chaque paire de patchs (i.e. avec décroissance exponentielle, avec __kernel_exp__ la vitesse de décroissance, en fonction de la distance), soit (ii) en fournissant directement la matrice de probabilité de dispersion interpatchs (__disp_mat__).

### Competition

Les interactions interspecifiques peuvent être ici encore soit spécifiées manuellement (via l'argument __int_mat__) ou automatiquement. Dans ce dernier cas, l'argument __intra__ donne les valeurs pour les coefficients de compétition intraspécifique (soit une valeur pour toutes les espèces, ou un vecteur de la taille du nombre d'espèces), __min_inter__ et __max_inter__ donnent les valeurs minimales et maximales pour les coefficients interspécifiques. Finalement, l'argument __comp_scaler__ permet de changer l'échelle pour tous les coefficients (i.e., augmente ou réduit les effets d'interactions p/r aux effets environnementaux, cf. Eq. 2). A noter; l'effet des interactions n'est pas dépendant des conditions environnementales.

## Trucs testables

+ 

### Références

Thompson, Patrick L., Laura Melissa Guzman, Luc De Meester, Zsófia Horváth, Robert Ptacnik, Bram Vanschoenwinkel, Duarte S. Viana, et Jonathan M. Chase. «A Process‐based Metacommunity Framework Linking Local and Regional Scale Community Ecology ». Ecology Letters 23, nᵒ 9 (septembre 2020): 1314‑29. https://doi.org/10.1111/ele.13568.

