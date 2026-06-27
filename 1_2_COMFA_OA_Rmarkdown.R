# [Set up function]
comfa_oa_compute <- function(input_df, output_path = NULL) {
  
  ###################Set up input parameters################
  # --- helpers ---
  # Vr - Xiaoyu Li's COMFA-O method (Mathematical Expectation for 0-360 wind)
  Vr_function <- function(Vac, Ws) {
    mapply(function(Vac_i, Ws_i) {
      integrand <- function(theta) sqrt(Vac_i^2 + Ws_i^2 - 2 * Vac_i * Ws_i * cos(theta))
      (1 / (2 * pi)) * integrate(integrand, lower = 0, upper = 2 * pi)$value
    }, Vac, Ws)
  }
  # [Purpose] tsf (clothes surface temperature) ~ (METc,Ta,Pv) Ominvar method (only used for basic sweat)
  # [目的] 皮肤表面温度 (仅用于估算基础出汗)
  Tsf_Ominivar_function <- function(METc, Pv, Ta) {
    35.7 - 0.028 * (METc) -
      0.155 * (METc - 3.05e-3 * (5733 - 6.99 * (METc) - 1000 * Pv) -
                 0.42 * ((METc) - 58.15) - 1.7e-5 * METc * (5867 - 1000 * Pv) -
                 0.0014 * METc * (34 - Ta))
  }
  ###################Set up input parameters finished#######
  
  # ----- Internal core calculator so we can reuse it for CET -----
  compute_core <- function(df) {
    # Input:
    Ta          <- df$Ta            # Celsius
    Rh          <- df$Rh            # %
    clo         <- df$clo           # clo
    Tmrt        <- df$Tmrt          # Celsius
    METs        <- df$METs          # in met
    Age         <- df$Age
    sex_raw     <- df$sex           # "male" = 1, "female" = 2
    Ws          <- df$Ws            # m/s
    Height      <- df$Height        # meter
    Weight      <- df$Weight        # kg
    bodyPosture <- df$bodyPosture   # "sitting" or "standing" or "crouching", default value = "standing"
    
    # Input data correction part
    # ---- enforce minimum clothes insulation ----
    clo_min <- 0.01   # example floor; adjust (e.g., 0.30 for sports, 0.50 for indoor summer tests)
    df$clo_raw  <- df$clo
    df$clo[is.na(df$clo)] <- clo_min  # Replace NAs with the floor, then clamp all values to at least CLO_MIN
    df$clo <- pmax(df$clo, clo_min)
    clo         <- df$clo
    # ---- enforce minimum wind speed (Ws) >= 0.1 m/s ----
    ws_min <- 0.1
    df$Ws_raw <- df$Ws
    df$Ws[is.na(df$Ws)] <- ws_min
    df$Ws <- pmax(df$Ws, ws_min)
    Ws <- df$Ws
    
    #####################COMFA_OA Calculation##############################
    # Stefan-Boltzmann constant
    kB <- 5.67e-08
    
    # Radiation efficiency [辐射效率]
    bp <- tolower(as.character(bodyPosture))
    rad_eff <- ifelse(is.na(bp), 0.895,
                      ifelse(bp == 'sitting', 0.691,
                             ifelse(bp == 'standing', 0.895,
                                    ifelse(bp == 'crouching', 0.6, 0.895))))
    
    # Barometric pressure [大气压力]
    pBaroHPA_setting <- 1013.25
    pBaroHPA <- 1013.25
    
    # Vapour pressure
    Pv <- ((0.6108 * (exp((17.269 * Ta) / (Ta + 237.3)))) * (Rh / 100)) # [Unit] kPa, Tetens formula
    
    # Revised Harris-benedict
    sex_char <- tolower(as.character(sex_raw))
    sex_is_male <- (sex_char %in% c("male", "m")) | (sex_raw %in% c(1, "1"))
    sex_std <- ifelse(sex_is_male, 1, 2)
    
    rmr <- ifelse(sex_is_male, 9.65 * Weight + 573 * Height - 5.8 * Age + 260,
                  7.38 * Weight + 607 * Height - 2.31 * Age + 43)
    Q_met_min <- rmr / 24 * 1.163 #Q_met_min: minimum meta heat by RMR
    
    ### Optimized Kozey formula
    VO2_original  <- 3.5
    VO2_corrected <- rmr / 1440 / 5 / Weight * 1000  # Harris-Benedict predicted RMR, Kcal/day to ml·kg⁻¹·min⁻¹
    METc <- 1 + (METs - 1) * (VO2_original / VO2_corrected)
    
    # BSA: 3D scan 2010
    BSA <- ifelse(sex_is_male,
                  79.8106 * (Height * 100)^0.7271 * Weight^0.3980 / 10000,
                  84.4673 * (Height * 100)^0.6997 * Weight^0.4176 / 10000)
    
    # Relative speed function
    # Effect of temperature and clothing thermal resistance on human sweat at low activity levels
    Vac <- ifelse(METs >= 1.5, 0.0052 * (METs * 58.15 - 58.15 * 1.5), 0) # activity speed
    
    # Vr - Xiaoyu Li's COMFA-O method (Mathematical Expectation for 0-360 wind)
    Vr <- Vr_function(Vac, Ws)
    
    ###  get the data from clo
    Rco  <- 0.155 * clo / 0.082 * 100
    Rcvo <- Rco * 2.74524064171123
    
    Trans <- 100 # Transmissivity (%)
    Rcv <- Rcvo * (-0.8 * (1 - exp(-Vr / 1.095)) + 1) # Clothing Vapour Resistance, from human clothing table
    
    Rc <- Rco * (-0.37 * (1 - exp(-Vr / 0.72)) + 1) # Clothing resistance, from energy table, Row AH
    
    # M: Metabolic energy(W/m2)
    f <- (0.15 - (0.0173 * Pv)) - (0.0014 * Ta) # Heat loss consumed through breathing (%)
    
    ### Qmet, avoid underestimation of Q_met by applying the RMR as minimum
    Q_met_raw <- (1 - f) * (rmr / 24 * 1.163) * (METc * 58.15) / 80
    Q_met <- ifelse(Q_met_raw <= Q_met_min, Q_met_min, Q_met_raw)
    
    ##################COMFA-O Tmrt Method#######################
    # COMFA - kabs
    kabs_labs_totl <- (Tmrt + 273.15)^4 * (5.67e-8 * 0.98)  # 0.98 is the skin emissivity
    
    # COMFA - Rabs
    Rabs <- kabs_labs_totl * rad_eff
    
    # Tc: core temperature
    Tc <- ifelse(Age < 20,
                 36.5 + (0.0043 * Q_met),
                 ifelse(sex_is_male,
                        (36.5 + (0.0043 * Q_met)) - 0.0060 * (Age - 20),    ###0.0060
                        (36.5 + (0.0043 * Q_met)) - 0.0066 * (Age - 20)))   ###0.0066
    
    # [Purpose] Tsk (skin temperature) ~ (Ta,Tmrt,Pv,Vr,METc,Tc)
    # New Tsk formula
    Tsk_high <- 12.2 + 0.020 * Ta + 0.044 * Tmrt + 0.194 * Pv - 0.253 * Vr + 0.00297 * METc * 58.15 + 0.513 * Tc
    Tsk_low  <- 7.2  + 0.064 * Ta + 0.061 * Tmrt + 0.198 * Pv - 0.348 * Vr + 0.616 * Tc
    Tsk_mid  <- (7.2 + 0.064 * Ta + 0.061 * Tmrt + 0.198 * Pv - 0.348 * Vr + 0.616 * Tc) +
      2.5 * ((12.2 + 0.020 * Ta + 0.044 * Tmrt + 0.194 * Pv - 0.253 * Vr + 0.513 * Tc) -
               (7.2 + 0.064 * Ta + 0.061 * Tmrt + 0.198 * Pv - 0.348 * Vr + 0.616 * Tc)) * (clo - 0.2)
    Tsk <- ifelse(clo >= 0.6, Tsk_high,
                  ifelse(clo > 0.2, Tsk_mid, Tsk_low))
    
    # TRemitted: The emitted terrestrial radiation (W/m2)
    Re <- 11333 * Vr
    # Ra is
    Ra <- ifelse(Re > 40000, (0.17 / ((0.71^0.33) * 0.000022 * 0.0266 * (Re^0.805))),
                 ifelse(Re < 4000, (0.17 / ((0.71^0.33) * 0.000022 * 0.683 * (Re^0.466))),
                        (0.17 / ((0.71^0.33) * 0.000022 * 0.193 * (Re^0.618)))))
    
    #############################Sweat heat loss##################
    f_clo_area <- (173.51 * clo - 2.36 - 100.76 * clo * clo + 19.28 * clo^3.0) / 100.0 # facl, in decimal %, percentage of body covered by clothes
    f_clo_area <- ifelse (f_clo_area >= 1, 1, f_clo_area)
    f_clo_area <- pmax(f_clo_area, 0) # To address <0 value for extreme low clo condition
    
    fCloAreaExpansion <- 1 + (0.31 * clo) # fcl, Spasic’s expansion formula, According to ISO 110793) and ISO 79331) (based on McCullough et al.9))
    
    areaClo  <- BSA * f_clo_area + BSA * (fCloAreaExpansion - 1.0)
    areaNude <- BSA * (1.0 - f_clo_area)
    
    # clothes surface temperature [Purpose: COMFA convective heat transfer model]
    # 服装表面温度 [目的] 用于COMFA 对流换热模型
    Tsf <- (((Tsk - Ta) / (Rc + Ra)) * Ra) + Ta # Edit 2025 01 24,
    
    # convective heat transfer coefficients for free convection modes
    # 对流换热系数 （自由对流）
    h_c_free <- 2.38 * (Tsf - Ta)^0.25
    h_c_free[is.nan(h_c_free)] <- -999 # Set the NAN value to -999 for the h_c_free
    
    # Convective heat transfer coefficient (forced convection)
    # 对流换热系数（强制对流）
    h_c_forced <- 2.67 + 6.5 * Vr^0.67 # convective heat transfer coefficients for forced convection
    h_c_forced <- h_c_forced * (pBaroHPA / pBaroHPA_setting)^0.55
    
    # Final hConv
    hConv <- pmax(h_c_free, h_c_forced) # compare free vs. forced and take the bigger value as the convective heat transfer coefficient    
    
    # Convection
    Q_conv_bare <- hConv * (Ta - Tsk) * areaNude
    Q_conv_clo  <- hConv * (Ta - Tsf) * areaClo
    Conv <- -(Q_conv_bare + Q_conv_clo) / BSA
    
    # Radeff method (SELECTED)
    TRemitted <- rad_eff * ((0.95 * (kB) * ((Tsk + 273.15)^4)))
    
    # [Purpose] basal sweat estimation (non-activity triggered sweat)
    # [目的] 基础出汗估算 (非运动刺激产生的出汗)
    # [Method] the sweat estimation by E_req and E_max, Gonzalez et al. adjusted by Omidvar
    Tsf_Ominivar <- Tsf_Ominivar_function(METc, Pv, Ta) # Clothes surface temperature, Ominivar method
    T_cl <- 29.187 - clo * (17.447 + 0.0011 * Pv * 1000 + 0.0504 * Ta) # Clothes surface temperature when met = 4
    h_cl <- 6.45 / clo # Heat conduction coefficient
    
    # Clothed vs. nude ratio
    f_cl <- ifelse(clo < 0.5, 1 + 0.2 * clo, 1.05 + 0.1 * clo) # Clothed vs. nude ratio
    
    # convective heat transfer coefficients for forced convection
    E_max <- (16.7 + 0.371 * clo * h_cl^2) * (4.13 - 0.001 * Pv * 1000) # Maximum evaporation limited by environmental parameters when met =4
    
    E_req <- 232.6 - f_cl * hConv * (T_cl - Ta) - 3.96 * 10^(-8) * f_cl * ((T_cl + 273.15)^4 - (Tmrt + 273.15)^4) # Required evaporation by energy balance when met =4
    
    # Es_b: basal sweat rate before adjustment
    # 基础出汗速率(调整前)
    Es_b <- 36.75 + 0.3818 * E_req - 0.2175 * E_max
    Es_b <- ifelse(E_req <= 0, 0, Es_b)
    
    # Beta coefficient to correct the basal sweat rate
    # Beta系数用于调整基础出汗率
    beta_sweat <- ifelse(Es_b >= 200, 1, ifelse(Es_b > 150, 0.5, ifelse(Es_b >= 0, 0.333, 0)))
    
    # Final: Total sweat heat loss calculation
    ## Basal sweat heat loss (environmental triggered)
    ## 由环境因素导致的基础出汗量
    q_sweat_basal <- ifelse(
      Age >= 20,
      (0.675 * beta_sweat / ((1 + 0.155 * clo * (4.6 + hConv + 0.046 * Tmrt))) * Es_b) * (1 - (0.005 * (Age - 20))),
      (0.675 * beta_sweat / ((1 + 0.155 * clo * (4.6 + hConv + 0.046 * Tmrt))) * Es_b)
    )
    
    #############################################################################################################################
    ## Define basal sweat initialization point
    q_sweat_basal <- ifelse((Q_met + Rabs - Conv - TRemitted) < 20, 0, q_sweat_basal) # When the net budget is less than 20 w/m2, we assume no basal sweat here
    
    q_sweat_act <- ifelse(Age >= 20,
                          0.42 * (Q_met - 58.15) * (1 - (0.005 * (Age - 20))),
                          0.42 * (Q_met - 58.15))
    q_sweat_act <- ifelse(q_sweat_act <= 0, 0, q_sweat_act) # basal sweat always >= 0
    
    ## Evaporative heat loss through sweat
    ## 出汗导致的总蒸发散热
    Es <- q_sweat_act + q_sweat_basal
    Es <- ifelse(E_req <= 0, 0, Es)
    
    qs <- 0.6108 * (exp((17.269 * Tsk) / (Tsk + 237.3))) # Specific humidity at skin
    ############Evap: The evaporative heat loss (W/m2)##################
    # [Purpose] Evaporative heat loss limited by environmental maximum evap ability
    Ei <- (1.16 * 2442) * ((qs - Pv) / (7700 + Rcv + (0.92 * Ra))) # Evaporative heat loss through diffusion (invisible)
    Em <- (1.16 * 2442) * ((qs - Pv) / (Rcv + (0.92 * Ra)))        # Maximal evaporative heat loss
    
    Evap <- Es + Ei
    Evap <- ifelse(Evap < Em, Evap, Em)
    ##############################################################################################################
    ##############################################################################################################
    
    BMI <- Weight / (Height)^2 # height in m, Weight in kg
    
    BF <- ifelse(sex_is_male,  (1.39 * BMI) + (0.16 * Age) - 19.34,
                 (1.39 * BMI) + (0.16 * Age) - 9)
    
    Tc_diff  <- ifelse(Tc <= 37, 37 - Tc, 0)
    Tsk_diff <- ifelse(Tsk <= 33, 33 - Tsk, 0)
    # Set the threshold of shivering (When COMFA < -20)
    cold_signal <- ifelse ((Q_met + Rabs - Conv - Evap - TRemitted) < (-15), 1, 0) # when the net COMFA-OA < -15, start to shivering
    
    Qshiv <- cold_signal * (155.5 * Tc_diff + 47 * Tsk_diff - 1.57 * Tsk_diff^2) / BF^0.5
    
    # --- compute COMFA-OA (as a valid R object name) ---
    comfa_oa <- Q_met + Rabs - Conv - Evap - TRemitted + Qshiv
    
    data.frame(COMFA_OA = comfa_oa, check.names = FALSE)
  }
  
  # ---------- Compute COMFA_OA for the provided data ----------
  df_main <- input_df
  out_core <- compute_core(df_main)
  comfa_vec <- out_core$COMFA_OA
  
  # ---------- CET: solve for Ta that reproduces COMFA_OA under the fixed indoor reference ----------
  # CET definition:
  # ---- robust root finder: safe bracketing without NA in while/if conditions ----
  # ---------- CET: solve for Ta that reproduces COMFA_OA under UTCI-style reference ----------
  # Reference condition: Tmrt = Ta, Ws = 0.3 m/s, Rh = 50%, clo = 0.9, sex = 1, Height = 1.75, Weight = 75,
  # Age = 65, METs = 2.3, bodyPosture = "standing"
  
  ref_constants <- list(
    Rh = 50, clo = 0.9, sex = 1, Height = 1.75, Weight = 75,
    Age = 65, METs = 2.3, bodyPosture = "standing",
    Ws = 0.3  # <-- use ~0.3 m/s at body height, 0.5m for 10m ws
  )
  
  # Build a template DF for the reference condition (one row; we'll vary Ta; Tmrt = Ta)
  make_ref_df <- function(Ta_guess, n = 1L) {
    data.frame(
      Ta = rep(Ta_guess, n),
      Rh = rep(ref_constants$Rh, n),
      clo = rep(ref_constants$clo, n),
      Tmrt = rep(Ta_guess, n),                # <-- UTCI: Tmrt equals Ta
      METs = rep(ref_constants$METs, n),
      Age = rep(ref_constants$Age, n),
      sex = rep(ref_constants$sex, n),
      Ws = rep(ref_constants$Ws, n),          # <-- fixed 0.3 m/s
      Height = rep(ref_constants$Height, n),
      Weight = rep(ref_constants$Weight, n),
      bodyPosture = rep(ref_constants$bodyPosture, n),
      check.names = FALSE
    )
  }
  
  # ---- robust root finder: safe bracketing without NA in while/if conditions ----
  solve_cet_for_target <- function(target_comfa, ta_center = 21, init_width = 20, max_expand = 5) {
    f <- function(ta_scalar) {
      ref_df <- make_ref_df(ta_scalar, n = 1L)
      compute_core(ref_df)$COMFA_OA[1] - target_comfa
    }
    
    # start with a bracket around ta_center
    low  <- ta_center - init_width
    high <- ta_center + init_width
    
    eval_f <- function(x) {
      val <- suppressWarnings(f(x))
      if (!is.finite(val)) NA_real_ else val
    }
    
    f_low  <- eval_f(low)
    f_high <- eval_f(high)
    
    # expand bracket a few times if needed (avoid NA in logical checks)
    expanded <- 0L
    sign_change <- is.finite(f_low) && is.finite(f_high) && (sign(f_low) != sign(f_high))
    
    while (!sign_change && expanded < max_expand) {
      low  <- low  - init_width
      high <- high + init_width
      f_low  <- eval_f(low)
      f_high <- eval_f(high)
      sign_change <- is.finite(f_low) && is.finite(f_high) && (sign(f_low) != sign(f_high))
      expanded <- expanded + 1L
    }
    
    if (!sign_change) {
      # final fallback: coarse scan to try to locate a sign change
      grid <- seq(low, high, by = 1)
      vals <- vapply(grid, eval_f, numeric(1))
      ok <- which(is.finite(vals))
      if (length(ok) >= 2) {
        idx <- which(sign(vals[ok][-1]) != sign(vals[ok][-length(ok)]))
        if (length(idx) >= 1) {
          i1 <- ok[idx[1]]
          i2 <- ok[idx[1] + 1]
          low  <- grid[i1];  high <- grid[i2]
          f_low  <- vals[i1]; f_high <- vals[i2]
          sign_change <- TRUE
        }
      }
    }
    
    if (!sign_change) return(NA_real_)
    
    # now safe to call uniroot
    out <- try(uniroot(f, lower = low, upper = high, tol = 1e-04, maxiter = 200), silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else out$root
  }
  
  # vectorize CET over rows; use observed Ta as center to help convergence
  ta_center_vec <- df_main$Ta
  cet_vals <- mapply(solve_cet_for_target, target_comfa = comfa_vec, ta_center = ta_center_vec)
  
  # --- append to input and (optionally) write CSV ---
  out <- df_main
  out[["COMFA_OA"]] <- comfa_vec
  out[["CET"]] <- as.numeric(cet_vals)
  
  if (!is.null(output_path)) {
    write.csv(out, output_path, row.names = FALSE)
  }
  return(out)
  
}

# ---- example (uncomment to run) ----
input <- read.csv("./2_COMFA_OA_example_input.csv")
output <- comfa_oa_compute(input, output_path = "./2_COMFA_OA_output.csv")
