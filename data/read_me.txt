OD_Benchmark_foreignBloc
    {"A":     Eq1_coeff_A,         # <-- matrix; 
     "mu":    Eq1_SigmaMu,         # <-- matrix;
     "R2":    Eq1_r2               # <-- array.
     }


OD_Benchmark_Seq
    "uncorrected_est":         Eq3bm_seq_impact_byCountry_med
    "uncorrected_mad":         Eq3bm_seq_impact_byCountry_mad
    "ssBias_est":              Tab2_seq_biasEst_byCountry_med
    "ssBias_mcN":              Tmontecarlo
    "corrected_est":           Tab2_seq_corrImpact_byCountry_med
    "corrected_mad":           Tab2_seq_corrImpact_byCountry_mad


OD_Benchmark_Jnt
    uncorrected_est         Eq3bm_joint_impact_byCountry_med
    uncorrected_mad         Eq3bm_joint_impact_byCountry_mad
    ssBias_est              Tab2_joint_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab2_joint_corrImpact_byCountry_med
    corrected_mad           Tab2_joint_corrImpact_byCountry_mad
    

OD_IntRate_foreignBloc
    A        Pg11_coeff_A            
    mu       Pg11_SigmaMu
    R2       Pg11_r2


OD_wIntRate_Seq
    uncorrected_est         Eq3tr_seq_impact_byCountry_med  @  Uncorrected Median Variance Share of Sequential Est
    uncorrected_mad         Eq3tr_seq_impact_byCountry_mad    
    ssBias_est              Tab3_seq_biasEst_byCountry_med  @  Sequential (w/ interest) Estimate of Bias
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab3_seq_corrImpact_byCountry_med   @  Sequential (w/ interest) Corrected Estimates: 
    corrected_mad           Tab3_seq_corrImpact_byCountry_mad   @  Sequential (w/ interest) MAD of Corrected Est.:


OD_wIntRate_Jnt
    uncorrected_est         Eq3tr_joint_impact_byCountry_med
    uncorrected_mad         Eq3tr_joint_impact_byCountry_mad    @@ round 2 removed
    ssBias_est              Tab3_joint_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab3_joint_corrImpact_byCountry_med
    corrected_mad           Tab3_joint_corrImpact_byCountry_mad


******* 换莫奈特卡洛数量 *******

OD_wGlobalY_Seq
    uncorrected_est         Tab4C_seq_impact_byCountry_med
    uncorrected_mad         Tab4C_seq_impact_byCountry_mad
    ssBias_est              Tab4C_seq_biasEst_byCountry
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab4C_seq_corrImpact_byCountry_med
    corrected_mad           Tab4C_seq_corrImpact_byCountry_mad


OD_1GlobalY_Seq
    uncorrected_est         Tab4D_seq_impact_byCountry_med
    uncorrected_mad         Tab4D_seq_impact_byCountry_mad
    ssBias_est              Tab4D_seq_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab4D_seq_corrImpact_byCountry_med
    corrected_mad           Tab4D_seq_corrImpact_byCountry_mad



OD_1Agri_Seq
    uncorrected_est         Tab5ii_seq_impact_byCountry_med
    uncorrected_mad         Tab5ii_seq_impact_byCountry_mad
    ssBias_est              Tab5ii_seq_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab5ii_seq_corrImpact_byCountry_med
    corrected_mad           Tab5ii_seq_corrImpact_byCountry_mad


OD_1Metal_Seq
    uncorrected_est         Tab5iii_seq_impact_byCountry_med
    uncorrected_mad         Tab5iii_seq_impact_byCountry_mad
    ssBias_est              Tab5iii_seq_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab5iii_seq_corrImpact_byCountry_med
    corrected_mad           Tab5iii_seq_corrImpact_byCountry_mad


OD_1Fuel_Seq
    uncorrected_est         Tab5iv_seq_impact_byCountry_med
    uncorrected_mad         Tab5iv_seq_impact_byCountry_mad
    ssBias_est              Tab5iv_seq_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab5iv_seq_corrImpact_byCountry_med
    corrected_mad           Tab5iv_seq_corrImpact_byCountry_mad


OD_1IntRate_Seq
    uncorrected_est         Tab5v_seq_impact_byCountry_med
    uncorrected_mad         Tab5v_seq_impact_byCountry_mad
    ssBias_est              Tab5v_seq_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab5v_seq_corrImpact_byCountry_med
    corrected_mad           Tab5v_seq_corrImpact_byCountry_mad


OD_1ToT_Seq
    uncorrected_est         Tab5iix_seq_impact_byCountry_med
    uncorrected_mad         Tab5iix_seq_impact_byCountry_mad
    ssBias_est              Tab5iix_seq_biasEst_byCountry_med
    ssBias_mcN              Tmontecarlo
    corrected_est           Tab5iix_seq_corrImpact_byCountry_med
    corrected_mad           Tab5iix_seq_corrImpact_byCountry_mad


                 

#####################################
###  @@@@@ @@@@@ @    @ @@@@@@@   ###
###  @     @   @ @   @  @         ###
###  @@@@@ @@@@@ @  @   @@@@@     ###
###      @ @   @ @ @    @         ###
###  @@@@@ @   @ @@     @@@@@@@   ###
###    Benchmark Seq. Results     ###
#####################################




## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
##                                                           ##
##        @@@@@@@@@   @@@@@@@@   @@@@@@@@   @@@@@@@@@        ## 
##           @        @          @             @             ##
##           @        @          @             @             ##
##           @        @@@@@@     @@@@@@@@      @             ##
##           @        @                 @      @             ##
##           @        @                 @      @             ##
##           @        @@@@@@@@@  @@@@@@@@      @             ##
##                                                           ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###




## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
##                                                           ##
##        @@@@@@@@      @@      @      @@@@@@@               ## 
##        @             @ @     @      @      @              ##
##        @             @  @    @      @       @             ##
##        @@@@@@        @   @   @      @        @            ##
##        @             @    @  @      @        @            ##
##        @             @     @ @      @       @             ##
##        @@@@@@@@      @      @@      @@@@@@@@              ##
##                                                           ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ###