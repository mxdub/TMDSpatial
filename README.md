# TMDSpatial
 
L'idée est ici de pouvoir simuler des jeux de données typiques à ceux utilisés en écologie des commuanutés, et éprouver les différentes méthodes disponibles face ces données.  
On réutilise le package {ecospat} mais un peu modifié. Les modifications notables sont actuellement: 
+ __Autocorrelation temporelle__ (au niveau des patches): Celle-ci est initialement inclue dans le package {ecospat}, i.e., il existe une variabilité temporelle, plus ou moins autocorrélée, dans les conditions environnementales des patches. Cette variabilité ajoute un niveau de complexité supplémentaire. Ainsi, celle-ci n'est pas considérer dans cette version du package. Voir l'argument __temporal_autocorr__ de la fonction simulate_MC() pour l'activer. A noter, que le même niveau d'autocorrélation est considéré à la fois pour le temps et pour l'espace.
+ __Extinction de patche__:
