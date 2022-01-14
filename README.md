# TMDSpatial
 
L'idée est ici de pouvoir simuler des jeux de données typiques à ceux utilisés en écologie des commuanutés, et éprouver les différentes méthodes disponibles face ces données.  
On réutilise le package {ecospat} mais légérement modifié. Les modifications notables sont actuellement: 
+ __Autocorrelation temporelle__ (au niveau des patches): Celle-ci est initialement inclue dans le package {ecospat}, i.e., il existe une variabilité temporelle, plus ou moins autocorrélée, dans les conditions environnementales des patches. Cette variabilité ajoute un niveau de complexité supplémentaire. Ainsi, celle-ci n'est pas considérée dans cette version du package. Voir l'argument __temporal_autocorr__ de la fonction simulate_MC() pour l'activer. A noter; le même niveau d'autocorrélation est considéré à la fois pour le temps et pour l'espace.
+ __Extinction de patche__: Voir les paramètres __extrip_prob__ et __extirp_all_pop__ de la fonction simulate_MC(). Le premier argument est un argument de la fonction de base, il donne la probabilité qu'à chaque pas de temps, une population s'éteigne dans un patch. Le second argument est lui en revanche nouveau. Il permet de spécifier si l'extinction concerne __toutes__ les espèces présentes dans un patche ou si ce tirage est réalisé pour chacune des espèces. A noter; si le probabilité est de 0.01 et que le patche contient 3 espèces, dans le cas où __extirp_all_pop = FALSE__, les trois espèces peuvent s'éteindre simultanément avec une prob. de 0.01 mais les extinctions individuelles ne sont pas possible via ce paramètre d'extirpation (mais possible par la dynamique locale). En revanche, si __extirp_all_pop = TRUE__, alors la probabilité pour les trois espèces de s'éteindre simultanément est ici : $0.01^3=10^{-6}$.


Il y a également quelques fonctions en plus pour faciliter l'utilisation des sorties de la fonction simulate_MC().
+ sim_to_matrix()
+ abund_to_occ()
+ plots_occupancies
