library(tidyverse)
library(ecospat)
library(TMDSpatial)

#### Simulate data ####

N_patch = 1000
N_species=30

t = simulate_MC(N_patch, N_species, temporal_autocorr = F, timesteps = 500, burn_in = 100)

t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 2, env1Scale = 50,
                temporal_autocorr = F,
                env_niche_breadth = 0.1, env_optima = c(0.7,0.2,0.2),
                int_mat = matrix(c(0.5,0.0,0.0,
                                   0.0,0.5,0.6,
                                   0.0,0.6,0.5), byrow = T, nrow = 3),
                dispersal = 0.01, kernel_exp = 0.3,
                extirp_prob = c(0),
                timesteps = 100, burn_in = 100, initialization = 0)

N_species=1
t = simulate_MC(N_patch, species = N_species,
                min_inter = 0, max_inter = 2,
                env_niche_breadth = 10, env_optima = c(0.5),
                dispersal = 0.01, kernel_exp = 0.00001,
                extirp_prob = c(0.05),
                timesteps = 100, burn_in = 100, initialization = 0)

occ = t$dynamics.df %>% select(-env_niche_breadth, -max_r, -optima, -env)

#### Transform output ####

to_array = function(x){
  x_ = x %>%spread(time, N)
  tmp = array(0, dim = c(length(unique(x$species)),
                         length(unique(x$patch)),
                         length(unique(x$time))))
  for(s in 1:length(unique(x$species))){
    tmp[s,,] = x_ %>% filter(species == s)%>%arrange(patch)%>%select(-patch, -species)%>%as.matrix()
  }
  tmp
}

occ = to_array(occ)

occupancies = occ
occupancies[occ>0]=1

## Proportion de site occupés
occupancies = apply(occupancies, c(1,3), sum) / N_patch
occupancies = occupancies %>% as_tibble() %>% rowid_to_column("species") %>%
  pivot_longer(-species) %>%
  mutate(name = str_replace(name, "V", "")) %>%
  mutate_at(.vars = c("name"), as.numeric)

ggplot(occupancies%>%filter(name>0)%>%filter(name%%2==0), aes(x=name, y=value, color = as.factor(species)))+
  geom_line()+
  scale_y_continuous(limits=c(0,1))

#### C-score ####

## Computing C-score
position_  = 125
occupancies = occ
occupancies[occ>0]=1

occ_cscore = t(occupancies[,,position_])
colnames(occ_cscore) = paste0("S", 1:N_species)
rownames(occ_cscore) = paste0("P", 1:N_patch)
ecospat.Cscore(occ_cscore, nperm = 1000, outpath = "./outputs/", verbose = T)


#### Var. part ####

library(vegan)

mod = varpart(occ_cscore, ~., pcnm(dist(t$landscape))$vectors, data = data.frame(env1 = t$env.df$env1[1:100]), transfo = 'hel')
mod

showvarparts(2, bg = c("hotpink","skyblue"))
plot(mod, bg = c("hotpink","skyblue"))

afrac = rda(decostand(occ_cscore, "hel"), model.matrix(~., data.frame(env1 = t$env.df$env1[1:100]))[,-1], pcnm(dist(t$landscape))$vectors)
anova(afrac, step = 200, perm.max = 200)

afrac = rda(decostand(occ_cscore, "hel"), pcnm(dist(t$landscape))$vectors, model.matrix(~., data.frame(env1 = t$env.df$env1[1:100]))[,-1])
anova(afrac, step = 200, perm.max = 200)


# BrayCurtis

mod = varpart(vegdist(occ_cscore), ~., pcnm(dist(t$landscape))$vectors, data = data.frame(env1 = t$env.df$env1[1:100]))
mod

showvarparts(2, bg = c("hotpink","skyblue"))
plot(mod, bg = c("hotpink","skyblue"))

afrac = rda(decostand(occ_cscore, "hel"), model.matrix(~., data.frame(env1 = t$env.df$env1[1:100]))[,-1], pcnm(dist(t$landscape))$vectors)
anova(afrac, step = 200, perm.max = 200)

afrac = rda(decostand(occ_cscore, "hel"), pcnm(dist(t$landscape))$vectors, model.matrix(~., data.frame(env1 = t$env.df$env1[1:100]))[,-1])
anova(afrac, step = 200, perm.max = 200)


#### JSDM (HMSC) ####

# Hmsc
library(Hmsc)

# Créer un tableau qui donne les id des sites pour l'effet random (de taille nrow() de ta matrice de presence/absence)
studyDesign = data.frame( samples = as.factor(rownames(occ_cscore))  )

# Créer un liste d'effet random, ici, juste un, l'effet 'site'
ranEff = list()
ranEff[['samples']] = HmscRandomLevel(units = unique(studyDesign$samples))

# Var env.
XData = tibble(env1 = t$env.df %>% filter(time_run == min(time_run)) %>% pull(env1))

# Et finalement, pour définir le modèle (du coup, a adapter avec tes variables à toi, c, x et t - je n'utilisais pas de traits pour ma part)
m=Hmsc(Y = as.matrix(occ_cscore),
       XData = as.data.frame(XData), XFormula = ~ 1+env1+I(env1^2),
       studyDesign = studyDesign, ranLevels = ranEff,
       distr = "probit")

m = sampleMcmc(m, samples = 1000, nChains = 3, nParallel = 3)

preds = computePredictedValues(m)
evaluateModelFit(hM = m, predY = preds)

computeAssociations(m)
Hmsc::getPostEstimate(m, "Beta")

# With Poisson
occupancies = occ

occ_cscore = t(occupancies[,,position_])
colnames(occ_cscore) = paste0("S", 1:N_species)
rownames(occ_cscore) = paste0("P", 1:N_patch)

m3=Hmsc(Y = as.matrix(occ_cscore),
        XData = as.data.frame(XData), XFormula = ~ 1+env1+I(env1^2),
        studyDesign = studyDesign, ranLevels = ranEff,
        distr = "poisson")

m3 = sampleMcmc(m3, samples = 1000, nChains = 3, nParallel = 3)
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


#### Testing exec time
library(profvis)
profvis({
  simulate_MC <- function(patches, species, dispersal = 0.01,
                          plot = TRUE,
                          torus = FALSE, kernel_exp = 0.1,
                          env1Scale = 500, temporal_autocorr = TRUE, timesteps = 1200, burn_in = 800, initialization = 200,
                          max_r = 5, min_env = 0, max_env = 1, env_niche_breadth = 0.5, optima_spacing = "random",
                          intra = 1, min_inter = 0, max_inter = 1, comp_scaler = 0.05,
                          extirp_prob = 0, extirp_all_pop = F,
                          landscape, disp_mat, env.df, env_optima, int_mat){
    if (missing(landscape)){
      landscape <- landscape_generate(patches = patches, plot = plot)
    } else {
      landscape <- landscape_generate(patches = patches, xy = landscape, plot = plot)
    }

    if (missing(disp_mat)){
      disp_mat <- dispersal_matrix(landscape = landscape,torus = torus, kernel_exp = kernel_exp, plot = plot)
    } else {
      disp_mat <- dispersal_matrix(landscape = landscape, disp_mat = disp_mat, torus = torus, kernel_exp = kernel_exp, plot = plot)
    }

    if (missing(env.df)){
      env.df <- env_generate_wrp(landscape = landscape, env1Scale = env1Scale, temporal_autocorr = temporal_autocorr, timesteps = timesteps+burn_in, plot = plot)
    } else {
      env.df <- env_generate_wrp(landscape = landscape, env.df = env.df, env1Scale = env1Scale, temporal_autocorr = temporal_autocorr, timesteps = timesteps+burn_in, plot = plot)
    }

    if (missing(env_optima)){
      env_traits.df <- env_traits(species = species, max_r = max_r, min_env = min_env, max_env = max_env, env_niche_breadth = env_niche_breadth, optima_spacing = optima_spacing, plot = plot)
    } else {
      env_traits.df <- env_traits(species = species, max_r = max_r, min_env = min_env, max_env = max_env, env_niche_breadth = env_niche_breadth, optima_spacing = optima_spacing, optima = env_optima, plot = plot)
    }

    if (missing(int_mat)){
      int_mat <- species_int_mat(species = species, intra = intra, min_inter = min_inter, max_inter = max_inter, comp_scaler = comp_scaler, plot = TRUE)
    } else {
      int_mat <- species_int_mat(species = species, int_mat = int_mat, intra = intra, min_inter = min_inter, max_inter = max_inter, comp_scaler = comp_scaler, plot = TRUE)
    }

    dynamics.df <- list()
    N <- matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
    pb <- txtProgressBar(min = 0, max = initialization + burn_in + timesteps, style = 3)
    for(i in 1:(initialization + burn_in + timesteps)){
      if(i <= initialization){
        if(i %in% seq(10,100, by = 10)){
          N <- N + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
        }
        env <- env.df$env1[env.df$time == 1]
      } else {
        env <- env.df$env1[env.df$time == (i-initialization)]
      }

      r <- max_r*exp(-(t((env_traits.df$optima - matrix(rep(env, each = species), nrow = species, ncol = patches))/(2*env_traits.df$env_niche_breadth)))^2)
      N_hat <- N*r/(1+N%*%int_mat)
      N_hat[N_hat < 0] <- 0
      N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches)

      E <- matrix(rbinom(n = patches * species, size = N_hat, prob = rep(dispersal, each = patches)), nrow = patches, ncol = species)
      dispSP <- colSums(E)
      I_hat_raw <- disp_mat%*%E
      I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
      I_hat[is.nan(I_hat)] <- 1
      I <- sapply(1:species, function(x) {
        if(dispSP[x]>0){
          table(factor(sample(x = patches, size = dispSP[x], replace = TRUE, prob = I_hat[,x]), levels = 1:patches))
        } else {rep(0, patches)}
      })

      N <- N_hat - E + I

      if(extirp_all_pop){
        N[rbinom(n = patches, size = 1, prob = extirp_prob)>0,] <- 0
      }else{
        N[rbinom(n = species * patches, size = 1, prob = extirp_prob)>0] <- 0
      }

      dynamics.df[[i]] = data.frame(N = c(N), patch = 1:patches, species = rep(1:species, each = patches), env = env, time = i-initialization-burn_in)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    dynamics.df = do.call(rbind, dynamics.df) %>% left_join(env_traits.df)
    env.df$time_run <- env.df$time - burn_in

    env.df_init <- data.frame(env1 = env.df$env1[env.df$time == 1], patch = 1:patches, time = NA, time_run = rep(seq(-(burn_in + initialization), -burn_in, by = 1), each = patches))
    env.df <- rbind(env.df_init,env.df)

    if(plot == TRUE){
      sample_patches <- sample(1:patches, size = min(c(patches,6)), replace = FALSE)
      g <- dynamics.df %>%
        filter(time %in% seq(min(dynamics.df$time),max(dynamics.df$time), by =10)) %>%
        filter(patch %in% sample_patches) %>%
        ggplot(aes(x = time, y = N, group = species, color = optima))+
        geom_line()+
        facet_wrap(~patch)+
        scale_color_viridis_c()+
        geom_path(data = filter(env.df, patch %in% sample_patches), aes(y = -5, x = time_run, color = env1, group = NULL), size = 3)

      print(g)
    }

    return(list(dynamics.df = dynamics.df, landscape = landscape, env.df = env.df, env_traits.df = env_traits.df, disp_mat = disp_mat, int_mat = int_mat))
  }
  simulate_MC(1000,30,temporal_autocorr = F)
})

