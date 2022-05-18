library(tidyverse)
library(ecospat)
library(TMDSpatial)
library(vegan)
library(scatterpie)

#### Simulate data ####

N_patch = 100
N_species=10

# Case : NO envt control, limited dispersal & NO envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 1,
                env_niche_breadth = 1,
                dispersal = 0.01, kernel_exp = 0.2,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100)

# Case : envt control, limited dispersal & NO envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 1,
                env_niche_breadth = 0.1,
                dispersal = 0.01, kernel_exp = 0.2,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100)

# Case : NO envt control, limited dispersal &  envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 999,
                env_niche_breadth = 1,
                dispersal = 0.01, kernel_exp = 0.2,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100)

# Case : envt control, limited dispersal & envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 999,
                env_niche_breadth = 0.1,
                dispersal = 0.01, kernel_exp = 0.2,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100, initialization = 0, local_start = F)

# Case : NO envt control, NO limited dispersal & NO envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 1,
                env_niche_breadth = 1,
                dispersal = 0.01, kernel_exp = 0.01,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100, initialization = 0, local_start = F)

# Case : envt control, NO limited dispersal & NO envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 1,
                env_niche_breadth = 0.1,
                dispersal = 0.01, kernel_exp = 0.01,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100, initialization = 0, local_start = F)

# Case : NO envt control, NO limited dispersal & envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 999,
                env_niche_breadth = 1,
                dispersal = 0.01, kernel_exp = 0.01,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100, initialization = 0, local_start = F)

# Case : envt control, NO limited dispersal & envt autocorrelation
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 0.8,
                temporal_autocorr = F, env1Scale = 999,
                env_niche_breadth = 0.1,
                dispersal = 0.01, kernel_exp = 0.01,
                extirp_prob = c(0.05),
                timesteps = 200, burn_in = 100, initialization = 0, local_start = F)

### Check spatial autocorrelation
plots_envt(t)

### Gather environmental & geo info
data_geo = get_geoposition(t)
data_envt = get_envt(t)

#### Transform output ####
abundances = sim_to_matrix(t)
occupancies = abund_to_occ(abundances)

## Plots occupancies
plots_occupancies(occupancies)

#### plot in time ######
# ## Plots in time
# n_steps = length(unique(t$dynamics.df$time))
# for(i in seq(1,n_steps,n_steps/10)){
#   snapshot = t(abundances[,,i])
# 
#   tp_ = tibble(x = t$landscape$x, y = t$landscape$y)
#   for(j in 1:N_species)
#     tp_ = tp_ %>% add_column(!!paste0('S', j) := snapshot[,j])
# 
#   ggplot(data = tp_)+
#     geom_point(aes(x=x,y=y), size = 9, shape = 1)+
#     geom_scatterpie(aes(x = x, y = y), data = tp_, cols = paste0("S", 1:N_species))+
#     theme_bw()
# 
#   ggsave(filename = paste0('p_', i, '.jpeg'))
# }
###########################


## Get Snapshot
position_  = 200
snapshot_occ = t(occupancies[,,position_]) %>% as_tibble()
snapshot_abund = t(abundances[,,position_]) %>% as_tibble()

####### SARs #########
## SAR
# compute_sa = function(comm, pos, i_min = 2, i_max = 10, max_plot = 10){
#   splits = i_min:i_max
#   
#   results = tibble(area = numeric(), S = numeric())
#   
#   for(i in splits){
#     
#     bands = 100 / i
#     start = bands / 2
#     end = 100 - start
#     centers = seq(start, end, bands)
#     
#     grid = expand.grid(x = centers, y = centers) %>% as_tibble()
#     grid = grid %>% slice( sample(1:dim(grid)[1], min(dim(grid)[1], max_plot)) )
#     
#     r = (bands / 2) * 0.95
#     
#     S = grid %>% mutate(S = map2_dbl(x, y,
#                                      .f = function(x, y) {  s = sqrt( (x-pos$x)^2+(y-pos$y)^2 ) <= r; sum(apply(comm[s,], 2, any))})) %>%
#       pull(S)
#     
#     results = results %>% add_row(area = 3.14*r^2, S = S )
#   }
#   
#   results
# }
# 
# sar = compute_sa(snapshot_occ, pos = data_geo)
# 
# ggplot(sar, aes(x=log(area), y=S))+geom_point()
# 
# fit <- sar_average(data = sar, grid_start = "none")
# summary(fit)
# plot(fit)
####### SARs #########

# Bray-Curtis : si les diff?rences absolue d'abondances sont importantes
# Hellinger (transfo) : si diff. relatives d'abondances sont diff?rentes et esp?ce communes plus importantes
# Chi2 : diff. relatives + esp?ces rares


#### Var. part ####

# Space description 
pcnms = scores(pcnm(dist(data_geo))) 

# Select pcnms
mod0 = rda(decostand(snapshot_abund, "hel") ~ 1, data = as_tibble(pcnms))
mod1 = rda(decostand(snapshot_abund, "hel") ~ ., data = as_tibble(pcnms))
os = ordistep(mod0, scope = formula(mod1))
# Which pcnm :
pcnm_tokeep = names(os$terminfo$ordered)

pcnms_reduced = pcnms[,pcnm_tokeep]
pcnms_reduced = pcnms

mod = varpart(snapshot_abund,
              ~env1+I(env1^2),
              pcnms_reduced,
              data = data_envt, transfo = 'hel')
mod

plot(mod, bg = c("hotpink","skyblue"))

# Environmental table 
frac_envt = rda(decostand(snapshot_abund, "hel")~env1+I(env1^2),
                data = data_envt)
frac_envt
RsquareAdj(frac_envt)
anova(frac_envt)

# Geo. table
frac_geo  = rda(decostand(snapshot_abund, "hel")~pcnms_reduced,
                data = data_envt)
frac_geo
anova(frac_geo)


# Envt indp from space
frac_envt_idp = rda(decostand(snapshot_abund, "hel")~env1+I(env1^2)+Condition(pcnms_reduced),
                data = data_envt)
frac_envt_idp
anova(frac_envt_idp)

# Space indp from envt
frac_geo_idp  = rda(decostand(snapshot_abund, "hel")~Condition(env1+I(env1^2))+pcnms_reduced,
            data = data_envt)
frac_geo_idp
anova(frac_geo_idp)


# BrayCurtis

mod0 = capscale(snapshot_abund ~ 1 ,data = as_tibble(pcnms), dist = 'bray')
mod1 = capscale(snapshot_abund ~ . ,data = as_tibble(pcnms), dist = 'bray')
os = ordistep(mod0, scope = formula(mod1))
# Which pcnm :
pcnm_tokeep = names(os$terminfo$ordered)
pcnms_reduced = pcnms[,pcnm_tokeep]

# For Bray-Curtis & Chi2 - no empty patch allows (zero division)
to_remove = which(apply(snapshot_abund, 1, sum) == 0)
snapshot_abund = snapshot_abund[-to_remove,] 
data_envt = data_envt[-to_remove,]
pcnms_reduced = pcnms_reduced[-to_remove,]

mod = varpart(vegdist(snapshot_abund, method = "bray"),
              ~env1+I(env1^2),
              pcnms_reduced,
              data = data_envt)
mod

plot(mod, bg = c("hotpink","skyblue"))


# Environmental table 
frac_envt = dbrda(snapshot_abund ~ env1+I(env1^2),data = data_envt, dist = 'bray')
frac_envt
RsquareAdj(frac_envt)
anova(frac_envt)

# Geo. table
frac_geo  = dbrda(snapshot_abund ~ pcnms_reduced, dist = 'bray')
frac_geo
RsquareAdj(frac_geo)
anova(frac_geo)

# Envt indp from space
frac_envt_idp = dbrda(snapshot_abund ~ env1+I(env1^2)+Condition(pcnms_reduced),data = data_envt, dist = 'bray')
frac_envt_idp
anova(frac_envt_idp)

# Space indp from envt
frac_geo_idp  = dbrda(snapshot_abund ~ Condition(env1+I(env1^2))+pcnms_reduced,data = data_envt, dist = 'bray')
frac_geo_idp
anova(frac_geo_idp)


# ChiSquare

mod0 = cca(snapshot_abund ~ 1 ,data = as_tibble(pcnms))
mod1 = cca(snapshot_abund ~ . ,data = as_tibble(pcnms))
os = ordistep(mod0, scope = formula(mod1))
# Which pcnm :
pcnm_tokeep = names(os$terminfo$ordered)
pcnms_reduced = pcnms[,pcnm_tokeep]


mod = varpart(snapshot_abund,
              ~env1+I(env1^2),
              pcnms_reduced,
              data = data_envt, chisquare = T)
mod

plot(mod, bg = c("hotpink","skyblue"))

# Environmental table 
frac_envt = cca(snapshot_abund ~ env1+I(env1^2),data = data_envt)
frac_envt
anova(frac_envt)

# Geo. table
frac_geo  = cca(snapshot_abund ~ pcnms_reduced)
frac_geo
anova(frac_geo)

# Envt indp from space
frac_envt_idp = cca(snapshot_abund ~ env1+I(env1^2)+Condition(pcnms_reduced),data = data_envt)
frac_envt_idp
anova(frac_envt_idp)

# Space indp from envt
frac_geo_idp  = cca(snapshot_abund ~ Condition(env1+I(env1^2))+pcnms_reduced,data = data_envt)
frac_geo_idp
anova(frac_geo_idp)


#### C-score #### (on occupancies)

# Need at least one co-occ...
# If observed C-score > expected C-score : competition, otherwise, "something else".
ecospat.Cscore(snapshot, nperm = 1000, outpath = "./outputs/", verbose = T)

# Constrained C-scores -  by envt.

# Fit glm (with quadratic effect - remember, envt effect on fitness is Gaussian)
inv.logit = function(x) 1/(1+exp(-x))
preds = tibble(dummy = rep(NA, N_patch))
for(i in 1:dim(snapshot)[2])
  preds = preds %>% add_column(!!colnames(snapshot)[i] := inv.logit(predict(glm(snapshot %>% pull(i)~ env1 + I(env1^2), family = 'binomial', data = data_envt))))
preds = preds %>% select(-dummy)

# Constrainted cscore
ecospat.cons_Cscore(snapshot,
                    pred = preds, nperm = 1000, outpath = "./outputs/", verbose = T)

# Utiliser cooc_null_model() de EcoSimR plutôt !
m = EcoSimR::cooc_null_model(snapshot)
summary(m)

#### JSDM (HMSC) ####

# Hmsc
library(Hmsc)

# Créer un tableau qui donne les id des sites pour l'effet random (de taille nrow() de ta matrice de presence/absence)
studyDesign = data.frame( samples = as.factor(rownames(snapshot))  )

# Créer un liste d'effet random, ici, juste un, l'effet 'site'
ranEff = list()
ranEff[['samples']] = HmscRandomLevel(units = unique(studyDesign$samples))

# Var env.
XData = tibble(env1 = t$env.df %>% filter(time_run == min(time_run)) %>% pull(env1))

# Model
m=Hmsc(Y = abund_to_occ(as.matrix(snapshot)),
       XData = as.data.frame(XData), XFormula = ~ 1+env1+I(env1^2),
       studyDesign = studyDesign, ranLevels = ranEff,
       distr = "probit")

m = sampleMcmc(m, samples = 1000, nChains = 3, nParallel = 3)

preds = computePredictedValues(m)
evaluateModelFit(hM = m, predY = preds)

computeAssociations(m)[[1]]$support
Hmsc::getPostEstimate(m, "Beta")

# With Poisson

m3=Hmsc(Y = as.matrix(snapshot),
        XData = as.data.frame(XData), XFormula = ~ 1+env1+I(env1^2),
        studyDesign = studyDesign, ranLevels = ranEff,
        distr = "poisson")

m3 = sampleMcmc(m3, samples = 1000, nChains = 3, nParallel = 3)

preds = computePredictedValues(m3)
evaluateModelFit(hM = m3, predY = preds)

computeAssociations(m3)
Hmsc::getPostEstimate(m3, "Beta")



#### JSDM (PLNnet, see others? ) ####

library(PLNmodels)
XData = as.data.frame(XData)
rownames(XData) = paste0("P", rownames(XData))

occ_cscore = occ_cscore[apply(occ_cscore, 1, sum)!=0,]

dt = prepare_data(occ_cscore, as.data.frame(XData)%>%mutate("env1.2" = env1^2))

O = log(rowSums(occ_cscore) %o% rep(1, ncol(occ_cscore)))

nm=PLNnetwork(Abundance ~ 1 + env1 + env1.2 + offset(O), data = dt)
plot(nm)
coefficient_path(nm, corr = TRUE) %>%
  ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) +
  geom_line(show.legend = FALSE) +  coord_trans(x="log10") + theme_bw()
bm = getBestModel(nm, "BIC")
plot(bm, type = "support", output = "corrplot")
bm$model_par$Omega

stability_selection(nm)
plot(nm, "stability")
bm = getBestModel(nm, "StARS")
plot(bm, type = "support", output = "corrplot")
bm$model_par$Omega
bm$model_par$Sigma


nm=PLN(Abundance ~ . + offset(O), data = dt, control = list(covariance = "diagonal"))
nm2=PLN(Abundance ~ . + offset(O), data = dt)
nm$model_par$Sigma
nm2$model_par$Sigma

nm$BIC
nm2$BIC

### Occupancy model

library(R2jags)

serie_length = 20
start_point = 100
species = 1

obs = occupancies[,,start_point:(start_point+serie_length-1)]
obs = obs[species,,]

plot(apply(obs[,], 2, mean), ylim = c(0,1))
plot(apply(occupancies[species,,], 2, mean), ylim = c(0,1))

n_site = dim(obs)[1]
n_year = dim(obs)[2]

datax <- list(obs = obs,
              n_site = dim(obs)[1],
              n_year = dim(obs)[2])

x = obs
x[is.na(x)] = sample(c(0,1), sum(is.na(x)), replace = T)
inits = list(list(x = x))

parameters <- c('e', 'c', 'psi', 'detec')

jm = jags(model.file = "./MacKenzie_mp.txt",
          data = datax,
          inits = inits,
          n.chains = 1,
          n.iter = 1000,
          n.burnin = 500,
          n.thin = 5,
          parameters.to.save = parameters)
# traceplot(jm)

e = jm$BUGSoutput$sims.list$e
c = jm$BUGSoutput$sims.list$c


n_rep = 50
rep = matrix(0, ncol = serie_length, nrow = n_rep)
for(r in 1:n_rep){
  z = obs
  for(i in 2:serie_length)
    z[,i] = sapply(z[,i-1]*(1-sample(e, 1))+(1-z[,i-1])*sample(c, 1), FUN = function(x) rbinom(1,1,x) )
  rep[r,] = apply(z, 2, mean)
}

rep = as_tibble(rep) %>% 
  add_rownames() %>%
  pivot_longer(-rowname) %>% 
  mutate(name = map_dbl(.x = name,
                        .f = function(x) as.numeric(str_remove(x, 'V'))))

obs_ = tibble(year = 1:serie_length,
              x = apply(obs, 2, mean))

ggplot(rep)+
  geom_boxplot(aes(x = name, y = value, group = name))+
  geom_point(data = obs_, aes(x=year, y = x), color = "red", size = 4)+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()


## Levins like

parameters <- c('e', 'gamma', 'psi', 'detec')
jm = jags(model.file = "./Levins_like.txt",
          data = datax,
            inits = inits,
          n.chains = 1,
          n.iter = 1000,
          n.burnin = 500,
          n.thin = 5,
          parameters.to.save = parameters)
# traceplot(jm)

e = jm$BUGSoutput$sims.list$e
gamma = jm$BUGSoutput$sims.list$gamma

time_out_fit = 100

rep = matrix(0, ncol = serie_length+time_out_fit, nrow = n_rep)
z = matrix(0, ncol = serie_length+time_out_fit, nrow = n_site)
for(r in 1:n_rep){
  z[,1] = obs[,1]
  for(i in 2:(serie_length+time_out_fit))
    z[,i] = sapply(z[,i-1]*(1-sample(e, 1))+(1-z[,i-1])*(1-exp(-sample(gamma, 1)*mean(z[,i-1]))), FUN = function(x) rbinom(1,1,x) )
  rep[r,] = apply(z, 2, mean)
}

rep = as_tibble(rep) %>% 
  add_rownames() %>%
  pivot_longer(-rowname) %>% 
  mutate(name = map_dbl(.x = name,
                        .f = function(x) as.numeric(str_remove(x, 'V'))))

obs_ = occupancies[species,,start_point:(start_point+serie_length-1+time_out_fit)]
obs_ = tibble(year = 1:(serie_length+time_out_fit),
              x = apply(obs_, 2, mean))

ggplot(rep)+
  geom_boxplot(aes(x = name, y = value, group = name))+
  geom_point(data = obs_, aes(x=year, y = x), color = "red", size = 4)+
  scale_y_continuous(limits = c(0,1))+
  theme_bw()
