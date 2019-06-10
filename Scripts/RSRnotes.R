tct.order <- geog.tct[order(geog.tct$GEOID),] %>% 
  mutate(phiocc = phi.occ,
         estphi = phi.occ.bias.summary$mn,
         diff = phiocc - estphi)

# Generate RSR matrices ---------------------------------------------------
Q.occ <- err.prec.matx.tct
A.occ <- diag(diag(Q.occ)) - Q.occ
A.occ[A.occ != 0] <- 1 #necessary because rho != 1
n.occ <- NCOL(Q.occ)
P.occ <- diag(n.occ) - X.tct %*% solve(crossprod(X.tct), t(X.tct))
Op.occ <- (nrow(A.occ)/sum(A.occ)) * (P.occ %*% (A.occ %*% P.occ))
numSpatre.occ <- round(nrow(geog.tct) * 0.10)
e <- rARPACK::eigs(Op.occ, numSpatre.occ) # e<-eigen(Op, symmetric = TRUE)

K.occ <- e$vectors[,1:numSpatre.occ] #we take the first q eigenvalues
Minv.occ <- crossprod(K.occ,Q.occ)%*%K.occ #Using the entire spatial lattice

Q.p <- err.prec.matx.bg
A.p <- diag(diag(Q.p)) - Q.p
A.p[A.p != 0] <- 1 #necessary because rho != 1
n.p <- NCOL(Q.p)
P.p <- diag(n.p) - X.bg %*% solve(crossprod(X.bg), t(X.bg))
Op.p <- (nrow(A.p)/sum(A.p)) * (P.p %*% (A.p %*% P.p))
numSpatre.det <- round(nrow(geog.bg) * 0.10)
e <- rARPACK::eigs(Op.p, numSpatre.det) # e<-eigen(Op, symmetric = TRUE)

K.p <- e$vectors[,1:numSpatre.det] #we take the first q eigenvalues
Minv.p <- crossprod(K.p,Q.p)%*%K.p #Using the entire spatial lattice

stan_data <- list(
  n_site = nrow(geog.tct),
  m_tct = m.tct,
  X_tct = X.tct,
  total_surveys = nrow(geog.bg),
  m_bg = m.bg,
  X_bg = X.bg, 
  site = survey.df$site,
  y = survey.df$y,
  start_idx = start_idx,
  end_idx = end_idx,
  any_seen = any_seen,
  n_survey = n_survey,
  numSpatreOcc = numSpatre.occ,
  numSpatreDet = numSpatre.det,
  K_occ = K.occ,
  K_det = K.p,
  Minv_occ = Minv.occ,
  Minv_det = Minv.p
)