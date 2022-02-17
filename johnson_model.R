# Johnson model

# Inputs
## PAR: PAR, umol PPFD m-2 s-1
## Temp: Leaf temperature, C
## CO2: Mesophyll CO2, ubar
## O2: Atmospheric O2, mbar
## Abs: Total leaf absorptance to PAR, mol PPFD absorbed mol-1 PPFD incident
## beta: PSII fraction of total leaf absorptance, mol PPFD absorbed by PSII mol-1 PPFD absorbed
## CB6F: Cyt b6f density, mol sites m-2
## RUB: Rubisco density, mol sites m-2
## Rds: Scalar for dark respiration, dimensionless
## Ku2: Rate constant for exciton sharing at PSII, s-1
## theta1: Curvature parameter for Aj/Ac transition, dimensionless
## eps1: PSI transfer function, mol PSI F to detector mol-1 PSI F emitted
## eps2: PSII transfer function, mol PSII F to detector mol-1 PSII F emitted
## alpha_opt: option for static or dymanic absorption cross-sections of PSI and PSII

model <- function(PAR = 800, Temp = 25, CO2 = 200, O2 = 209, Abs = 0.85, beta = 0.52, 
                  CB6F = ((350 / 300) / 1e6), RUB = ((100 / 3.6) / 1e6), Rds = 0.01,
                  Ku2 = 0e09, theta1 = 1, eps1 = 0, eps2 = 1, 
                  alpha_opt = "static"){
  
  # Environmental variables
  Q <- PAR / 1e6      # PPFD, mol PAR m-2 s-1
  Tc <- Temp / 1e6    # Leaf temperature, degrees C
  C <- CO2 / 1e6      # Partial pressure of CO2, bar
  O <- O2 / 1e3       # Partial pressure of O2, bar
  
  # Photochemical constants
  Kf <- 0.05e09       # Rate constant for fluoresence at PSII and PSI, s-1
  Kd <- 0.55e09       # Rate constant for constitutive heat loss at PSII and PSI, s-1
  Kp1 <- 14.5e09      # Rate constant for photochemistry at PSI, s-1
  Kn1 <- 14.5e09      # Rate constant for regulated heat loss at PSI, s-1
  Kp2 <- 4.5e09       # Rate constant for photochemistry at PSII, s-1
  
  # Biochemical constants
  kq <- 300           # Cyt b6f kcat for PQH2, mol e-1 mol sites-1 s-1
  nl <- 0.75          # ATP per e- in linear flow, ATP/e-
  nc <- 1             # ATP per e- in cyclic flow, ATP/e-
  kc <- 3.6           # Rubisco kcat for CO2, mol CO2 mol sites-1 s-1
  ko <- 3.6 * 0.27    # Rubisco kcat for O2, mol O2 mol sites-1 s-1
  Kc <- 260 / 1e6     # Rubisco Km for CO2, bar
  Ko <- 179000 / 1e6  # Rubisco Km for O2, bar
  
  # Calculate derived variables
  Vqmax <- CB6F * kq # maximum Cyt b6f activity, mol e-1 m-2 s-1
  Vcmax <- RUB * kc # maximum Rubisco activity, mol CO2 m-2 s-1
  Rd <- Vcmax * Rds # mitochondrial respiration, mol CO2 m-2 s-1
  S <- (kc / Kc) * (Ko / ko) # Rubisco specificity for CO2/O2, dimensionless
  gammas <- O / (2 * S) # CO2 compensation point in the absence of Rd, bar
  eta <- (1 - (nl / nc) + (3 + 7 * gammas / C) / ((4 + 8 * gammas / C) * nc)) # PSI/II ETR
  phi1P_max <- Kp1 / (Kp1 + Kd + Kf) # maximum photochemical yield PSI
  
  # Establish PSII and PSI cross-sections, mol PPFD abs PSII/PSI mol-1 PPFD
  if(alpha_opt == "static"){
    a2 <- Abs * beta
    a1 <- Abs - a2
  }else if(alpha_opt == "dynamic"){
    a2 <- solve_xcs(Abs, CB6F, Kd, Kf, Kp2, Ku2, Q, eta, kq, phi1P_max)
    a1 <- Abs - a2
  } else {
    stop("argument alpha_opt not idetified")
  }
  
  # Calculate limiting rates for gas-exchange and electron transport
  ## Expressions for potential Cytochrome b6f-limited rates (_j) N.B., see Eqns. 30-31
  JP700_j <- (Q * Vqmax) / (Q + Vqmax / (a1 * phi1P_max)) # rate of electron transport through PSI
  JP680_j <- JP700_j / eta # rate of electron transport through PSII
  Vc_j <- JP680_j / (4 * (1 + 2 * gammas / C))
  Vo_j <- Vc_j * 2 * gammas / C
  Ag_j <- Vc_j - Vo_j / 2 # potential rate of net CO2 assimilation under Cyt b6f limitation
  
  ## Expressions for potential Rubisco-limited rates (_c) N.B., see Eqns. 32-33
  Vc_c <- C * Vcmax / (C + Kc * (1 + O / Ko))
  Vo_c <- Vc_c * 2 * gammas / C
  Ag_c <- Vc_c - Vo_c / 2 # potential rate of net CO2 assimilation under Rubisco limitation
  JP680_c <- Ag_c * 4 * (1 + 2 * gammas / C) / (1 - gammas / C) # total rate of LEF through PSII
  JP700_c <- JP680_c * eta # total rate of LEF and CEF1 through PSI
  
  # Select min of Rubisco-limited and Cyt b6f-limited rates
  
  ## Define anonymous function for quadratic to smooth transitions N.B., 
  ## this returns an array with the two roots of the quadratic
  tr <- function(l1, l2, th){
    
    data.frame("q1" = ((l1 + l2) + sqrt((l1 + l2)^2 - 4 * th * l1 * l2)) / (2 * th),
               "q2" = ((l1 + l2) - sqrt((l1 + l2)^2 - 4 * th * l1 * l2)) / (2 * th))
  }
  
  ## N.B., min(X,[],2) returns column vector with minimum value of each row; 
  ## these minimum rates are interpreted as the actual rates (_a)
  
  # Select minimum PSI ETR 
  JP700_a <- apply(tr(JP700_j, JP700_c, theta1), 1, min)
  
  # Select minimum PSII ETR 
  JP680_a <- apply(tr(JP680_j, JP680_c, theta1), 1, min)
  
  # Select minimum Ag_a 
  ## gross rate of CO2 assimilation
  Ag_a  <- apply(tr(Ag_j, Ag_c, theta1), 1, min) * as.integer(C > gammas) +
    apply(tr(Ag_j, Ag_c, theta1), 1, max) * as.integer(C <= gammas)
  ## net rate of CO2 assimilation
  An_a <- Ag_a - Rd 
  
  # Derive a2/a1 at light saturation point and update N.B.,
  # this represents dynamic optimization of a2/a1 under limiting light,
  # and then under saturating light holds a2/a1 at the values attained at the light saturation point.
  
  if(alpha_opt == "dynamic"){
    
    JP700_j_c <- data.frame("JP700_j" = JP700_j, "JP700_c" = JP700_c)
    I <- apply(JP700_j_c, 1, which.min)
    I_diff <- diff(I, 1, 1)
    position <- which(I_diff != 0)
    a2_new <- rep(a2[position], length(a2))
    a2_old <- a2
    a2_update <- data.frame("a2_old" = a2_old, "a2_new" = a2_new)
    a2 <- apply(a2_update, 1, max)
    a1 <- Abs - a2
    
  }
  
  # Derive fluorescence parameters from gas-exchange and electron transport
  
  ## Primary fluorescence parameters
  CB6F_a <- JP700_j / kq        # Concentration of Cyt b6f with the Qp site occupied by PQH (mol m−2), Eqns. 21, 30a, 34
  phi1P_a <- JP700_a / (Q * a1) # First-order turnover constant for open PSII centers (s−1), Eqn. 20
  q1_a <- phi1P_a / phi1P_max   # Ratio of the concentration of open PSI reaction centers to the total concentration of PSI reaction centers, Eqn. 19a
  phi2P_a <- JP680_a / (Q * a2) # First-order turnover constant for open PSII centers (s−1), Eqn. 26
  q2_a <- 1 - CB6F_a / CB6F     # Ratio of the concentration of closed PSII reaction centers to the total concentration of PSII reaction centers, Eqns. 28 and 34
  
  ## N.B., rearrange Eqn. 25a to solve for Kn2_a, 
  Kn2_a <- ((Kp2^2 * phi2P_a^2 - 2 * Kp2^2 * phi2P_a * q2_a + Kp2^2 * q2_a^2 -
               4 * Kp2 * Ku2 * phi2P_a^2 * q2_a + 2 * Kp2 * Ku2 * phi2P_a^2 +
               2 * Kp2 * Ku2 * phi2P_a * q2_a + Ku2^2 * phi2P_a^2)^(1 / 2) -
              Kp2 * phi2P_a + Ku2 * phi2P_a + Kp2 * q2_a) / (2 * phi2P_a) - Kf - Ku2 - Kd
  
  # Derived fluorescence parameters -- 'True values'
  
  ## Photosystem II (Eqns. 23a-23e and 25a-25d)
  ### internal yields of the whole bed of PSII units for photochemistry, 
  ### regulated heat loss in the antennae, constitutive heat loss in the antennae, 
  ### fluorescence, and inter-unit exciton sharing (mol energy dissipated mol−1 energy absorbed)
  phi2p_a <- (q2_a) * Kp2 / (Kp2 + Kn2_a + Kd + Kf + Ku2)
  phi2n_a <- (q2_a) * Kn2_a / (Kp2 + Kn2_a + Kd + Kf + Ku2) + (1 - q2_a) * Kn2_a / (Kn2_a + Kd + Kf + Ku2) 
  phi2d_a <- (q2_a) * Kd / (Kp2 + Kn2_a + Kd + Kf + Ku2) + (1 - q2_a) * Kd / (Kn2_a + Kd + Kf + Ku2)
  phi2f_a <- (q2_a) * Kf / (Kp2 + Kn2_a + Kd + Kf + Ku2) + (1 - q2_a) * Kf / (Kn2_a + Kd + Kf + Ku2)
  phi2u_a <- (q2_a) * Ku2 / (Kp2 + Kn2_a + Kd + Kf + Ku2) + (1 - q2_a) * Ku2 / (Kn2_a + Kd + Kf + Ku2)
  ### overall yields of the whole bed of PSII units for photochemistry, 
  ### regulated heat loss in the antennae, constitutive heat loss in the antennae,
  ### and fluorescence (mol energy dissipated mol−1 energy absorbed).
  phi2P_a <- phi2p_a / (1 - phi2u_a)
  phi2N_a <- phi2n_a / (1 - phi2u_a)
  phi2D_a <- phi2d_a / (1 - phi2u_a)
  phi2F_a <- phi2f_a / (1 - phi2u_a)
  
  ## For Photosystem I (Eqns. 19a-19d)
  ### overall yields of the whole bed of PSI units for photochemistry,
  ### constitutive heat loss from closed reaction centers, 
  ### constitutive heat loss from the antennae, and fluorescence
  ### (mol energy dissipated mol−1 energy absorbed)
  phi1P_a <- (q1_a) * Kp1 / (Kp1 + Kd + Kf)
  phi1N_a <- (1 - q1_a) * Kn1 / (Kn1 + Kd + Kf)
  phi1D_a <- (q1_a) * Kd / (Kp1 + Kd + Kf) + (1 - q1_a) * Kd / (Kn1 + Kd + Kf)
  phi1F_a <- (q1_a) * Kf / (Kp1 + Kd + Kf) + (1 - q1_a) * Kf / (Kn1 + Kd + Kf)
  
  # Derived fluorescence parameters -- 'Observed values'
  
  ## PAM measured fluorescence levels (Eqns. 38-42)
  ## N.B., hardcoding of a2(1) for dark-adapted value
  Fm_a <- a2[1] * Kf / (Kd + Kf) * eps2 + a1[1] * Kf / (Kn1 + Kd + Kf) * eps1 # maximum fluorescence level in a dark-adapted leaf
  Fo_a <- a2[1] * Kf / (Kp2 + Kd + Kf) * eps2 + a1[1] * Kf / (Kp1 + Kd + Kf) * eps1 # minimum fluorescence level in a dark-adapted leaf
  Fmp_a <- a2 * Kf / (Kn2_a + Kd + Kf) * eps2 + a1 * Kf / (Kn1 + Kd + Kf) * eps1 # maximum fluorescence level in a light-adapted leaf
  Fop_a <- a2 * Kf / (Kp2 + Kn2_a + Kd + Kf) * eps2 + a1 * Kf / (Kp1 + Kd + Kf) * eps1 # minimum fluorescence level in a light-adapted leaf
  Fs_a <- a2 * phi2F_a * eps2 + a1 * phi1F_a * eps1 # steady-state fluorescence
  
  ## PAM indices used in plotter_forward_fun.m
  PAM1_a <- 1 - Fs_a / Fmp_a # PhiP
  PAM2_a <- Fs_a * (1 / Fmp_a - 1 / Fm_a) # PhiN
  PAM3_a <- Fs_a / Fm_a # PhiD + PhiF
  
  # Other PAM indices used in paper
  ## PAM4_a <- Q * 0.85./2 * PAM1_a; % ETR
  ## PAM5_a <- (Fmp_a - Fs_a)./(Fmp_a - Fop_a); % qP
  ## PAM6_a <- (Fmp_a - Fs_a) * Fop_a./((Fmp_a - Fop_a) * Fs_a); % qL
  ## PAM7_a <- PAM4_a./(1-PAM5_a); % kPuddle
  ## PAM8_a <- PAM4_a./(1-PAM6_a); % kLake
  ## PAM9_a <- Fm_a./Fmp_a - 1; % NPQ
  
  results <- data.frame("PAR"       = PAR,        # PAR, umol PPFD m-2 s-1
                        "Temp"      = Temp,       # Leaf temperature, C
                        "CO2"       = CO2,        # Mesophyll CO2, ubar
                        "O2"        = O2,         # Atmospheric O2, mbar
                        "Q"         = Q,          # PPFD, mol PAR m-2 s-1
                        "C"         = C,          # Partial pressure of CO2, bar
                        "O"         = O,          # Partial pressure of O2, bar
                        "Abs"       = Abs,        # Total leaf absorptance to PAR, mol PPFD absorbed mol-1 PPFD incident
                        "beta"      = beta,       # PSII fraction of total leaf absorptance, mol PPFD absorbed by PSII mol-1 PPFD absorbed
                        "CB6F"      = CB6F,       # Cyt b6f density, mol sites m-2
                        "RUB"       = RUB,        # Rubisco density, mol sites m-2
                        "Rds"       = Rds,        # Scalar for dark respiration, dimensionless
                        "Ku2"       = Ku2,        # Rate constant for exciton sharing at PSII, s-1
                        "theta1"    = theta1,     # Curvature parameter for Aj/Ac transition, dimensionless
                        "eps1"      = eps1,       # PSI transfer function, mol PSI F to detector mol-1 PSI F emitted
                        "eps2"      = eps2,       # PSII transfer function, mol PSII F to detector mol-1 PSII F emitted
                        "Ag_a"      = Ag_a,       # Gross rate of CO2 assimilation
                        "Ag_c"      = Ag_c,       # Potential rate of net CO2 assimilation under Rubisco limitation
                        "Ag_j"      = Ag_j,       # Potential rate of net CO2 assimilation under Cyt b6f limitation
                        "An_a"      = An_a,       # Net rate of CO2 assimilation
                        "CB6F_a"    = CB6F_a,     # Concentration of Cyt b6f with the Qp site occupied by PQH (mol m−2), 
                        "Fm_a"      = Fm_a,       # Maximum fluorescence level in the dark
                        "Fmp_a"     = Fmp_a,      # Maximum fluorescence level in the light
                        "Fo_a"      = Fo_a,       # Minimum fluorescence level in the dark
                        "Fop_a"     = Fop_a,      # Minimum fluorescence level in the light
                        "Fs_a"      = Fs_a,       # Steady-state fluorescence
                        "JP680_a"   = JP680_a,    # Minimum PSII ETR
                        "JP680_c"   = JP680_c,    # Total rate of LEF through PSII, mol e− m−2 s−1
                        "JP680_j"   = JP680_j,    # Rate of electron transport through PSII
                        "JP700_a"   = JP700_a,    # Minimum PSI ETR
                        "JP700_c"   = JP700_c,    # Total rate of LEF and CEF1 through PSI, mol e− m−2 s−1
                        "JP700_j"   = JP700_j,    # Rate of electron transport through PSI
                        "Kn2_a"     = Kn2_a,      # Inter-unit exciton sharing for regulated heat loss in the antennae (s-1)
                        "PAM1_a"    = PAM1_a,     # Photochemistry yield
                        "PAM2_a"    = PAM2_a,     # Regulated NPQ yield
                        "PAM3_a"    = PAM3_a,     # Residual yield
                        "Rd"        = Rd,         # Mitochondrial respiration, mol CO2 m-2 s-1
                        "S"         = S,          # Rubisco specificity for CO2/O2, dimensionless
                        "Vcmax"     = Vcmax,      # Maximum Rubisco activity, mol CO2 m-2 s-1
                        "Vqmax"     = Vqmax,      # Maximum Cyt b6f activity, mol e-1 m-2 s-1
                        "a1"        = a1,         # PSI cross-section, mol PPFD abs PSII/PSI mol-1 PPFD
                        "a2"        = a2,         # PSII cross-section, mol PPFD abs PSII/PSI mol-1 PPFD
                        "eta"       = eta,        # PSI/II ETR
                        "gammas"    = gammas,     # CO2 compensation point in the absence of Rd, bar
                        "phi1P_max" = phi1P_max,  # Maximum photochemical yield PSI
                        "phi1D_a"   = phi1D_a,    # Overall yield of constitutive heat loss from the antennae, mol energy dissipated mol−1 energy absorbed
                        "phi1F_a"   = phi1F_a,    # Overall yield of fluorescence, mol energy dissipated mol−1 energy absorbed
                        "phi1N_a"   = phi1N_a,    # Overall yield of constitutive heat loss from closed reaction centers, mol energy dissipated mol−1 energy absorbed
                        "phi1P_a"   = phi1P_a,    # Overall yield of photochemistry, mol energy dissipated mol−1 energy absorbed
                        "phi2D_a"   = phi2D_a,    # Overall yield of constitutive heat loss in the antennae, mol energy dissipated mol−1 energy absorbed
                        "phi2F_a"   = phi2F_a,    # Overall yield of fluorescence, mol energy dissipated mol−1 energy absorbed
                        "phi2N_a"   = phi2N_a,    # Overall yield of regulated heat loss in the antennae, mol energy dissipated mol−1 energy absorbed
                        "phi2P_a"   = phi2P_a,    # Overall yield of photochemistry, mol energy dissipated mol−1 energy absorbed
                        "phi2d_a"   = phi2d_a,    # Internal yield of constitutive heat loss in the antennae, mol energy dissipated mol−1 energy absorbed
                        "phi2f_a"   = phi2f_a,    # Internal yield of fluorescence, mol energy dissipated mol−1 energy absorbed
                        "phi2n_a"   = phi2n_a,    # Internal yield of regulated heat loss in the antennae, mol energy dissipated mol−1 energy absorbed
                        "phi2p_a"   = phi2p_a,    # Internal yield of photochemistry, mol energy dissipated mol−1 energy absorbed
                        "phi2u_a"   = phi2u_a,    # Internal yield of inter-unit exciton sharing, mol energy dissipated mol−1 energy absorbed
                        "q1_a"      = q1_a,       # Ratio of the concentration of open PSI reaction centers to the total concentration of PSI reaction centers
                        "q2_a"      = q2_a)       # Ratio of the concentration of closed PSII reaction centers to the total concentration of PSII reaction centers
  
  return(results)
  
}

# Symbolic solution for optimal absorption cross-sections
solve_xcs <- function(Abs, CB6F, Kd, Kf, Kp2, Ku2, Q, eta, kq, phi1P_max){
  
  xcs <- (-sqrt((Kd + Kf + Ku2) * (CB6F^2 * Kd^3 * kq^2 * phi1P_max^2 + CB6F^2 * Kf^3 * kq^2 * phi1P_max^2 + 
                                     CB6F^2 * Kd * Kp2^2 * eta^2 * kq^2 + CB6F^2 * Kf * Kp2^2 * eta^2 * kq^2 + 
                                     CB6F^2 * Kp2^2 * Ku2 * eta^2 * kq^2 + CB6F^2 * Kd * Kf^2 * kq^2 * 
                                     phi1P_max^2 * 3.0 + CB6F^2 * Kd^2 * Kf * kq^2 * phi1P_max^2 * 3.0 + 
                                     CB6F^2 * Kd * Kp2^2 * kq^2 * phi1P_max^2 + CB6F^2 * Kd^2 * Kp2 * kq^2 * 
                                     phi1P_max^2 * 2.0 + CB6F^2 * Kf * Kp2^2 * kq^2 * phi1P_max^2 + CB6F^2 * 
                                     Kf^2 * Kp2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kd^2 * Ku2 * kq^2 * 
                                     phi1P_max^2 + CB6F^2 * Kf^2 * Ku2 * kq^2 * phi1P_max^2 + CB6F^2 * Kp2^2 * 
                                     Ku2 * kq^2 * phi1P_max^2 + Abs^2 * Kd * Kp2^2 * Q^2 * eta^2 * phi1P_max^2 + 
                                     Abs^2 * Kf * Kp2^2 * Q^2 * eta^2 * phi1P_max^2 + Abs^2 * Kp2^2 * Ku2 * 
                                     Q^2 * eta^2 * phi1P_max^2 + CB6F^2 * Kd * Kf * Kp2 * kq^2 * phi1P_max^2 * 
                                     4.0 + CB6F^2 * Kd * Kf * Ku2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kd * 
                                     Kp2 * Ku2 * kq^2 * phi1P_max^2 * 2.0 + CB6F^2 * Kf * Kp2 * Ku2 * kq^2 * 
                                     phi1P_max^2 * 2.0 + CB6F^2 * Kd * Kp2^2 * eta * kq^2 * phi1P_max * 2.0 + 
                                     CB6F^2 * Kd^2 * Kp2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kf * 
                                     Kp2^2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kf^2 * Kp2 * eta * kq^2 * 
                                     phi1P_max * 2.0 + CB6F^2 * Kp2^2 * Ku2 * eta * kq^2 * phi1P_max * 2.0 + 
                                     CB6F^2 * Kd * Kf * Kp2 * eta * kq^2 * phi1P_max * 4.0 + CB6F^2 * Kd * 
                                     Kp2 * Ku2 * eta * kq^2 * phi1P_max * 2.0 + CB6F^2 * Kf * Kp2 * Ku2 * eta * 
                                     kq^2 * phi1P_max * 2.0 + Abs * CB6F * Kd * Kp2^2 * Q * eta * kq * 
                                     phi1P_max^2 * 2.0 + Abs * CB6F * Kd * Kp2^2 * Q * eta^2 * kq * phi1P_max * 
                                     2.0 + Abs * CB6F * Kd^2 * Kp2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * 
                                     CB6F * Kf * Kp2^2 * Q * eta * kq * phi1P_max^2 * 2.0 + Abs * CB6F * Kf * 
                                     Kp2^2 * Q * eta^2 * kq * phi1P_max * 2.0 + Abs * CB6F * Kf^2 * Kp2 * Q * 
                                     eta * kq * phi1P_max^2 * 2.0-Abs * CB6F * Kp2^2 * Ku2 * Q * eta * kq * 
                                     phi1P_max^2 * 2.0 + Abs * CB6F * Kp2^2 * Ku2 * Q * eta^2 * kq * 
                                     phi1P_max * 2.0 + Abs * CB6F * Kd * Kf * Kp2 * Q * eta * kq * 
                                     phi1P_max^2 * 4.0 + Abs * CB6F * Kd * Kp2 * Ku2 * Q * eta * kq * 
                                     phi1P_max^2 * 2.0 + Abs * CB6F * Kf * Kp2 * Ku2 * Q * eta * kq * 
                                     phi1P_max^2 * 2.0)) + 
            (CB6F * Kd^2 * kq * phi1P_max + CB6F * Kf^2 * kq * phi1P_max + Abs * Kd^2 * Q * phi1P_max^2 * 2.0 + 
               Abs * Kf^2 * Q * phi1P_max^2 * 2.0 + CB6F * Kd * Kp2 * eta * kq + CB6F * Kf * Kp2 * eta * kq + 
               CB6F * Kp2 * Ku2 * eta * kq + CB6F * Kd * Kf * kq * phi1P_max * 2.0 + CB6F * Kd * Kp2 * kq * 
               phi1P_max + CB6F * Kf * Kp2 * kq * phi1P_max + CB6F * Kd * Ku2 * kq * phi1P_max + CB6F * Kf * 
               Ku2 * kq * phi1P_max + CB6F * Kp2 * Ku2 * kq * phi1P_max + Abs * Kd * Kf * Q * phi1P_max^2 * 
               4.0 + Abs * Kd * Kp2 * Q * phi1P_max^2 * 2.0 + Abs * Kf * Kp2 * Q * phi1P_max^2 * 2.0 + Abs * 
               Kd * Ku2 * Q * phi1P_max^2 * 2.0 + Abs * Kf * Ku2 * Q * phi1P_max^2 * 2.0 + Abs * Kd * Kp2 * Q * 
               eta * phi1P_max + Abs * Kf * Kp2 * Q * eta * phi1P_max + Abs * Kp2 * Ku2 * Q * eta * phi1P_max)) / 
    (Q * phi1P_max * (Kd^2 * phi1P_max + Kf^2 * phi1P_max + Kd * Kp2 * eta + Kf * Kp2 * eta + Kp2 * 
                        Ku2 * eta + Kd * Kf * phi1P_max * 2.0 + Kd * Kp2 * phi1P_max + Kf * Kp2 * 
                        phi1P_max + Kd * Ku2 * phi1P_max + Kf * Ku2 * phi1P_max) * 2.0)
  
  return(xcs)
}

