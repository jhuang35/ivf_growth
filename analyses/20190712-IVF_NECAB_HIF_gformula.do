use "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\DATA\ANALYTIC\20190125-IVF_Gen_MI_Drop.dta", clear
keep if _mi_m == 1
gen MOM_AGE = mother_age_delivery
gen MOM_HT = m_height_pw26
gen FAT_HT = f_height_m24
gen FAT_WT = f_weight_m24
gen FAST = ogtt_fasting_pw26
gen OGTT = ogtt_2hour_pw26

gen sbp_avg = (sbp1 + sbp2 + sbp3) / 3
	replace sbp_avg = (sbp2 + sbp3) / 2 if sbp_avg == .	
	replace sbp_avg = (sbp1 + sbp3) / 2 if sbp_avg == .
	replace sbp_avg = (sbp1 + sbp2) / 2 if sbp_avg == .
gen dbp_avg = (dbp1 + dbp2 + dbp3) / 3
	replace dbp_avg = (dbp2 + dbp3) / 2 if dbp_avg == .	
	replace dbp_avg = (dbp1 + dbp3) / 2 if dbp_avg == .
	replace dbp_avg = (dbp1 + dbp2) / 2 if dbp_avg == .

gen ASBP = sbp_avg
gen ADBP = dbp_avg

gen HIF = cg27146050
gen NEC = cg03904042
save "C:\Users\JHUANGYH\Desktop\merged_MI.dta", replace

/**************************************
RUN G-COMPUTATIONS
***************************************/
matrix out = J(96,9,.)
// 6 measures x 16 visits
set more off
scalar n = 1
tokenize 0 1 3 6 9 12 15 18 24 36 48 54 60 66 72 78
quietly forvalues i = 1/16 {
	noisily di _newline _newline "Month ``i''"
	use "C:\Users\JHUANGYH\Desktop\merged_MI.dta", clear
	keep if visit == ``i''
	set more off
	// CHOOSE WHICH CPG MEDIATOR
	local CPG HIF
	*local CPG NEC
	foreach var of varlist zlen height zwei log_weight log_bmi zbmi {
		noisily di _newline _newline "`CPG': `var'" 
		gformula IVF `var' `CPG' ///
			MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT FAST ASBP GA, mediation ///
			outcome(`var') exposure(IVF) mediator(`CPG')					///
			post_confs(FAST ASBP GA)								///
			base_confs(MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT)	/// baseline(IVF:0) 																///
			control(`CPG':0)																	/// 
			commands(`CPG':regress, FAST:regress, ASBP:regress, GA:regress, `var':regress) ///
			equations(	FAST: IVF MOM_AGE i.MOM_EDU i.ETHNIC i.HH_INC MOM_HT ppBMI MALE SCORE parity i.smk_home FAT_HT FAT_WT, 	///
						ASBP: IVF MOM_AGE i.MOM_EDU i.ETHNIC i.HH_INC MOM_HT ppBMI MALE SCORE parity i.smk_home FAT_HT FAT_WT, 	///
						GA: FAST ASBP IVF MOM_AGE i.MOM_EDU i.ETHNIC i.HH_INC MOM_HT ppBMI MALE SCORE parity i.smk_home FAT_HT FAT_WT, 	///
						`CPG': FAST ASBP GA IVF MOM_AGE i.MOM_EDU i.ETHNIC i.HH_INC MOM_HT ppBMI MALE SCORE parity i.smk_home FAT_HT FAT_WT, ///
						`var': `CPG' FAST ASBP IVF MOM_AGE i.MOM_EDU i.ETHNIC i.HH_INC MOM_HT ppBMI MALE SCORE parity i.smk_home FAT_HT FAT_WT) ///
			seed(42782) samples(100) obe
			matrix out[n,1] = r(tce)
			matrix out[n,2] = r(se_tce)
			matrix out[n,3] = r(nde)
			matrix out[n,4] = r(se_nde)
			matrix out[n,5] = r(nie)
			matrix out[n,6] = r(se_nie)
			matrix out[n,7] = r(cde)
			matrix out[n,8] = r(se_cde)
			matrix out[n,9] = ``i''
			scalar n = n + 1
			}
		}
svmat out, names(gcomp)
export excel gcomp1-gcomp9 using "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\IVF and Growth\Output/$S_DATE-gformula_`CPG'.xls", firstrow(variables) replace			
			
			
/*
// OTHER STUFF			
capture drop no_cpg
capture drop cpg
logistic IVF MOM_AGE ppBMI MOM_HT i.MOM_EDU SCORE FAT_WT FAT_HT i.ETHNIC
predict no_cpg, xb
logistic IVF NEC MOM_AGE ppBMI MOM_HT i.MOM_EDU SCORE FAT_WT FAT_HT i.ETHNIC
predict cpg, xb
roccomp IVF no_cpg cpg, graph summary title("ROC for IVF status with or without NECAB3") ///
	subtitle("(N = 83 IVF; 1035 spontaneous)") rlopts(lpattern(dash)) ///
	scheme(s1color) plot1opts(msize(vsmall)) plot2opts(msize(vsmall)) ///
	note("Predictors: maternal age, ethnicity, BMI, height, edu, PRS; paternal height, weight; fetal cg03904042")			
		
//rocfit IVF p, continuous(10)
			
			
gformula IVF zlen NEC ///
	MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT ///
	GA FAST OGTT ASBP ADBP, mediation ///
	outcome(zlen) exposure(IVF) mediator(NEC)					///
	post_confs(FAST ASBP GA)								///
	base_confs(MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT)	///
	baseline(IVF:0) 																///
	control(NEC:0)																	/// 
	commands(NEC:regress, GA:regress, FAST:regress, ASBP:regress, zlen:regress) ///
	equations(	FAST: IVF MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT, 	///
				ASBP: IVF MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT, 	///
				NEC: IVF GA FAST ASBP MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT, ///
				GA: IVF FAST ASBP MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT, 	///
				zlen: NEC IVF GA FAST ASBP MOM_AGE MOM_EDU ETHNIC HH_INC MOM_HT ppBMI MALE SCORE parity smk_home FAT_HT FAT_WT) ///
			seed(42782) samples(50) 			
			
capture drop len_78			
predict len_78

xi: regress HIF NEC IVF NECxIVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
			m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity ///
			i.smk_home GA days if visit == 0, robust
			
lincom NEC + NECxIVF
*/
