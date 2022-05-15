library(tidyverse)
library(ecospat)
library(TMDSpatial)
library(vegan)
library(scatterpie)

#### Simulate data ####

N_patch = 400
N_species=3

t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 2,
                temporal_autocorr = F, env1Scale = 1,
                env_niche_breadth = 0.1, env_optima = c(0.3,0.5,0.7),
                int_mat = matrix(c(0.5,0.0,0.0,
                                   0.0,0.5,0.5,
                                   0.0,0.5,0.5), byrow = T, nrow = 3),
                dispersal = 0.01, kernel_exp = 0.001,
                extirp_prob = c(0),
                timesteps = 500, burn_in = 100, initialization = 0)

N_patch = 250
N_species=25
t = simulate_MC(N_patch, species = N_species,
                timesteps = 500, burn_in = 100, initialization = 0, env_niche_breadth = 0.2,
                env1Scale = 1, kernel_exp = 0.001, dispersal=0.001)

### Check spatial autocorrelation
plots_envt(t)

### Gather environmental & geo info
data_geo = t$landscape %>% as_tibble()
data_envt = tibble( env1 = t$env.df$env1[1:N_patch] ) # We assume no temporal variability for envt., so only keep data for first year

#### Transform output ####
abundances = sim_to_matrix(t)
occupancies = abund_to_occ(abundances)

## Plots occupancies
plots_occupancies(occupancies)

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




## Get Snapshot
position_  = 350
snapshot_occ = t(occupancies[,,position_]) %>% as_tibble()
snapshot_abund = t(abundances[,,position_]) %>% as_tibble()


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

# Bray-Curtis : si les différences absolue d'abondances sont importantes
# Hellinger (transfo) : si diff. relatives d'abondances sont différentes et espèce communes plus importantes
# Chi2 : diff. relatives + espèces rares

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

# Caution, empty site are not allow in CCA (zero division)
to_delete = which(apply(snapshot_abund, 1, sum) == 0)
snapshot_abund=snapshot_abund[-to_delete,]
pcnms_reduced=pcnms_reduced[-to_delete,]
data_envt=data_envt[-to_delete,]

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

# Utiliser cooc_null_model() de EcoSimR plutÃ´t !
m = EcoSimR::cooc_null_model(snapshot)
summary(m)

#### JSDM (HMSC) ####

# Hmsc
library(Hmsc)

# CrÃ©er un tableau qui donne les id des sites pour l'effet random (de taille nrow() de ta matrice de presence/absence)
studyDesign = data.frame( samples = as.factor(rownames(snapshot))  )

# CrÃ©er un liste d'effet random, ici, juste un, l'effet 'site'
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


#### MP mod. (unmarked) ####

library(unmarked)
N <- 1000
nspecies <- 3
J <- 5

occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('occ_cov',1:10,sep='')

det_covs <- list()
for (i in 1:nspecies){
  det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
}
names(det_covs) <- paste('det_cov',1:nspecies,sep='')

#True vals
beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
f1 <- beta[1] + beta[2]*occ_covs$occ_cov1
f2 <- beta[3] + beta[4]*occ_covs$occ_cov2
f3 <- beta[5] + beta[6]*occ_covs$occ_cov3
f4 <- beta[7]
f5 <- beta[8]
f6 <- beta[9]
f7 <- beta[10]
f <- cbind(f1,f2,f3,f4,f5,f6,f7)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
colnames(z) <- paste('sp',1:nspecies,sep='')
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

psi <- exp(f %*% t(dm))
psi <- psi/rowSums(psi)

#True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
  ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
}

p_true <- c(0.6,0.7,0.5)

# fake y data
y <- list()

for (i in 1:nspecies){
  y[[i]] <- matrix(NA,N,J)
  for (j in 1:N){
    for (k in 1:J){
      y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
    }
  }
}
names(y) <- c('coyote','tiger','bear')

#Create the unmarked data object
data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

#Summary of data object
summary(data)
plot(data)

#### MP mod. (hand-made) ####

library(R2jags)

occupancies = occ
occupancies[occ>0]=1
occupancies = occupancies[,,seq(1,200,length.out = 20)]
occupancies = structure(occupancies, .Dim = c(1,N_patch,20))

n_sp = dim(occupancies)[1]
n_site = dim(occupancies)[2]
n_year = dim(occupancies)[3]

datax <- list(obs = occupancies,
              p = occupancies,
              X= cbind(XData$env1, XData$env1^2),
              mask = 1-diag(rep(1,n_sp)),
              n_site = n_site,
              n_year = n_year,
              n_sp = n_sp)

inits <- list(list(mu_g = rnorm(n_sp), mu_e = rnorm(n_sp),
                   beta_g = matrix(rnorm(n_sp*2), nrow = n_sp),
                   beta_e = matrix(rnorm(n_sp*2), nrow = n_sp),
                   delta_g = matrix(rnorm(n_sp*n_sp), nrow = n_sp),
                   delta_e = matrix(rnorm(n_sp*n_sp), nrow = n_sp)))

parameters <- c("mu_g","mu_e",
                "beta_g","beta_e",
                "delta_g", "delta_e")


jm = jags.parallel(model.file = "../../models/mod1",
          data = datax,
          n.chains = 3,
          inits = inits,
          n.iter = 5000,
          n.burnin = 1000,
          n.thin = 5,
          parameters.to.save = parameters)
jm
plot(jm)
jm = as.mcmc(jm)

jm[[1]][,c(8,1,2)]
jm[[1]][,c(9,3,4)]
