use "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\DATA\ANALYTIC\20190723-merge_CpG_dad.dta", clear

//sum father_*
//codebook father_*

gen GDM = (gdm_who_1999 == "Yes") / (1 - (gdm_who_1999 == ""))
gen MALE = (sex == "Male") / (sex != "")
gen ETHNIC = ((mother_ethnicity == "indian") + 2*(mother_ethnicity == "malay")) / (1 - (mother_ethnicity == "" | mother_ethnicity == "others"))
gen MOM_EDU = (	(mother_highest_education == "gce") + 				///
				(mother_highest_education == "secondary")*2 + 		///
				(mother_highest_education == "ite_ntc")*3 + 		///
				(mother_highest_education == "primary")*4 + 		///
				(mother_highest_education == "no_education")*5) / 	///
				(mother_highest_education != "")
gen HH_INC = (	(household_income == "4000_5999") + 	///
				(household_income == "2000_3999")*2 + 	///
				(household_income == "1000_1999")*3 + 	///
				(household_income == "0_999")*4) / 		///
				(household_income != "")
gen log_bmi = ln(bmi)
gen log_weight = ln(weight)
gen log_tri = ln(triceps)
gen log_bi = ln(biceps)
gen log_sub = ln(subscapular)
gen log_supra = ln(suprailiac)
gen log_fat_5 = ln(fat_5_qmr)
gen log_fat_6 = ln(fat_6_qmr)
				
// IMPUTATION HERE *******************************

// Begin imputation 
drop if MALE == .
replace MOM_EDU = 4 if MOM_EDU == 5 // no education to primary, thus category = "primary or less"
replace father_HBP = 0 if father_HBP == 1 & father_HBP_age > father_age_delivery & father_HBP_age != . & father_age_delivery != . // N = 2 IVF, 14 SC
replace father_diabetes = 0 if father_diabetes == 1 & father_diabetes_age > father_age_delivery & father_diabetes_age != . & father_age_delivery != . // N = 2 IVF, 5 SC

mi set flong

// IMPUTED VARS
// impute by regression
mi register imputed mother_age_delivery m_height_pw26 f_height_m24 f_weight_m24 ppBMI /// 
					ogtt_fasting_pw26 ogtt_2hour_pw26 GA SCORE `cpgs' sbp1 sbp2 sbp3 dbp1 dbp2 dbp3
mi register imputed father_age_delivery // NEW
mi register imputed zlen height log_weight zwei log_bmi zbmi ///
					EFW_19 EFW_26 EFW_32 SBP DBP glucose_1 ///
					log_tri log_bi log_sub log_supra // outcome measures					
// impute by PMM
mi register imputed	HH_INC ETHNIC MOM_EDU parity smk_home hi_bp 
mi register imputed father_HBP father_diabetes // NEW
				
// COMPLETE OR CALCULATED VARS
mi register regular MALE IVF 

gen zlen_miss = missing(zlen)
gen height_miss = missing(height)
gen weight_miss = missing(log_weight)
gen zwei_miss = missing(zwei)
gen bmi_miss = missing(log_bmi)
gen zbmi_miss = missing(zbmi)
gen EFW_19_miss = missing(EFW_19)
gen EFW_26_miss = missing(EFW_26)
gen EFW_32_miss = missing(EFW_32)
gen SBP_miss = missing(SBP)
/*
gen DBP_miss = missing(DBP)
gen ZSBP_miss = missing(ZSBP)
gen ZDBP_miss = missing(ZDBP)
*/
gen glucose_miss = missing(glucose_1)
gen tri_miss = missing(log_tri)
gen bi_miss = missing(log_bi)
gen sub_miss = missing(log_sub)
gen supra_miss = missing(log_supra)
				
set more off
mi impute chained 	///
	(regress) 	mother_age_delivery m_height_pw26 f_height_m24 f_weight_m24 ppBMI /// regress
				ogtt_fasting_pw26 ogtt_2hour_pw26 GA SCORE `cpgs' sbp1 sbp2 sbp3 dbp1 dbp2 dbp3 ///
				father_age_delivery ///
				zlen height log_weight zwei log_bmi zbmi ///
				EFW_19 EFW_26 EFW_32 SBP DBP glucose_1 ///
				log_tri log_bi log_sub log_supra ///
	(pmm, knn(3)) 	HH_INC ETHNIC MOM_EDU parity smk_home hi_bp /// pmm
					father_HBP father_diabetes ///
	= 		MALE IVF, ///
	replace add(10) rseed(42782) // augment

save "D:\SICS\Archive\20190723-IVF_Dad_MI.dta", replace
use "D:\SICS\Data and Instruments\Archive\20190723-IVF_Dad_MI.dta", clear

// DROP MISSING DEPENDENT VARIABLES
replace zlen = . if zlen_miss //; sleep 500
replace height = . if height_miss //; sleep 500
replace log_weight = . if weight_miss //; sleep 500
replace zwei = . if zwei_miss //; sleep 500
replace log_bmi = . if bmi_miss //; sleep 500
replace zbmi = . if zbmi_miss //; sleep 500
/*
mi xeq: replace log_fat_5 = . if fat_5_miss
mi xeq: replace lean_5_qmr = . if lean_5_miss
mi xeq: replace log_fat_6 = . if fat_6_miss
mi xeq: replace lean_6_qmr = . if lean_6_miss
*/
replace EFW_19 = . if EFW_19_miss //; sleep 500
replace EFW_26 = . if EFW_26_miss //; sleep 500
replace EFW_32 = . if EFW_32_miss //; sleep 500

replace SBP = . if SBP_miss //; sleep 500
replace DBP = . if SBP_miss //; sleep 500
replace ZSBP = . if SBP_miss //; sleep 500
replace ZDBP = . if SBP_miss //; sleep 500

replace glucose_1 = . if glucose_miss //; sleep 500

replace log_tri = . if tri_miss //; sleep 500
replace log_bi = . if bi_miss //; sleep 500
replace log_sub = . if sub_miss //; sleep 500
replace log_supra = . if supra_miss //; sleep 500

save "D:\SICS\Archive\20190723-IVF_Dad_MI_Drop.dta", replace
use "D:\SICS\Archive\20190723-IVF_Dad_MI_Drop.dta", clear
	
mi passive: gen sbp_avg = (sbp1 + sbp2 + sbp3)/3
mi passive: gen dbp_avg = (dbp1 + dbp2 + dbp3)/3

save "D:\SICS\Archive\20190723-IVF_Dad_MI_Drop.dta", replace
use "D:\SICS\Archive\20190723-IVF_Dad_MI_Drop.dta", clear
		
		
		
		
//************************************************
// ESTIMATE HERE
// test
/*
xi: mi estimate: regress zlen IVF mother_age_delivery father_age_delivery ///
					i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE parity ///
					father_HBP father_diabetes ///
					i.smk_home ogtt_fasting_pw26 ogtt_2hour_pw26 GA if visit == 72
*/
/*
foreach var of varlist SBP log_tri log_bi log_sub log_supra {
	di "`var:'"
	codebook visit if `var' != . & _mi_m == 0, tab(20)
	}
*/
// SBP visits: 			36 48 60 72
// tri, sub visits: 	0, 18-78
// bi visits:			18-78
// supra visits: 		48-78

global est_list "" // Store parameter estimates

tokenize 0 1 3 6 9 12 15 18 24 36 48 54 60 66 72 78
foreach y of varlist zlen height zwei log_weight log_bmi zbmi {
	di _newline _newline 
	di "***************************"
	di "`y':"
		quietly forvalues t = 1/16 {
		noisily di "``t'': " 
		*crude models
		tempfile `y'_crude
		quietly: xi: mi estimate: regress `y' IVF ///
					if visit == ``t''
		parmest, flist(est_list) saving(``y'_crude', replace) stars(0.05 0.01) idstr("regress `y' IVF crude") escal(N)
		*min models - SES related + anthro
		tempfile `y'_min
		quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
					father_age_delivery father_HBP father_diabetes ///
					if visit == ``t''
		parmest, flist(est_list) saving(``y'_min', replace) stars(0.05 0.01) idstr("regress `y' IVF min") escal(N)
		*max models - SES and prenatal health / related to fertility
		tempfile `y'_adj
		quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
					father_age_delivery father_HBP father_diabetes ///
					ogtt_fasting_pw26 ogtt_2hour_pw26 sbp_avg dbp_avg GA days if visit == ``t''
		parmest, flist(est_list) saving(``y'_adj', replace) stars(0.05 0.01) idstr("regress `y' IVF adj") escal(N)
		*count if IVF == 1 & e(sample) == 1
		*lincom IVF
		}
	}
foreach y of varlist EFW_19 EFW_26 EFW_32 glucose_1 {
	di _newline _newline 
	di "***************************"
	di "`y':"
	*crude models
		tempfile `y'_crude
		quietly: xi: mi estimate: regress `y' IVF ///
					if visit == 0
		parmest, flist(est_list) saving(``y'_crude', replace) stars(0.05 0.01) idstr("regress `y' IVF crude") escal(N)
	*min models - SES related + anthro
		tempfile `y'_min
		quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
					father_age_delivery father_HBP father_diabetes ///
					if visit == 0
		parmest, flist(est_list) saving(``y'_min', replace) stars(0.05 0.01) idstr("regress `y' IVF min") escal(N)
	*max models - SES and prenatal health / related to fertility
		tempfile `y'_adj
		quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
					father_age_delivery father_HBP father_diabetes ///
					ogtt_fasting_pw26 ogtt_2hour_pw26 sbp_avg dbp_avg GA days if visit == 0
		parmest, flist(est_list) saving(``y'_adj', replace) stars(0.05 0.01) idstr("regress `y' IVF adj") escal(N)
	}
	
tokenize 36 48 60 72
foreach y of varlist SBP ZSBP DBP ZDBP {
	di _newline _newline 
	di "***************************"
	di "`y':"
		quietly forvalues t = 1/4 {
		noisily di "``t'': " 
		*crude models
		tempfile `y'_crude
		quietly: xi: mi estimate: regress `y' IVF ///
					if visit == ``t''
		parmest, flist(est_list) saving(``y'_crude', replace) stars(0.05 0.01) idstr("regress `y' IVF crude") escal(N)
		*min models - SES related + anthro
		tempfile `y'_min
		quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
					father_age_delivery father_HBP father_diabetes ///
					if visit == ``t''
		parmest, flist(est_list) saving(``y'_min', replace) stars(0.05 0.01) idstr("regress `y' IVF min") escal(N)
		*max models - SES and prenatal health / related to fertility
		tempfile `y'_adj
		quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
					m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
					father_age_delivery father_HBP father_diabetes ///
					ogtt_fasting_pw26 ogtt_2hour_pw26 sbp_avg dbp_avg GA days if visit == ``t''
		parmest, flist(est_list) saving(``y'_adj', replace) stars(0.05 0.01) idstr("regress `y' IVF adj") escal(N)
		*count if IVF == 1 & e(sample) == 1
		*lincom IVF
		}
	}

tokenize 0 18 24 36 48 54 60 66 72 78
foreach y of varlist log_tri log_bi log_sub log_supra {
	di _newline _newline 
	di "***************************"
	di "`y':"
		quietly forvalues t = 1/10 {
		noisily di "``t'': " 
		*crude models
		if ("`y'" == "log_tri") | ("`y'" == "log_sub") | ///
		(("`y'" == "log_bi") & ``t'' != 0) | ///
		(("`y'" == "log_supra") & ``t'' >= 48) {
			tempfile `y'_crude
			quietly: xi: mi estimate: regress `y' IVF ///
						if visit == ``t''
			parmest, flist(est_list) saving(``y'_crude', replace) stars(0.05 0.01) idstr("regress `y' IVF crude") escal(N)
			*min models - SES related + anthro
			tempfile `y'_min
			quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
						father_age_delivery father_HBP father_diabetes ///
						if visit == ``t''
			parmest, flist(est_list) saving(``y'_min', replace) stars(0.05 0.01) idstr("regress `y' IVF min") escal(N)
			*max models - SES and prenatal health / related to fertility
			tempfile `y'_adj
			quietly: xi: mi estimate: regress `y' IVF mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home ///
						father_age_delivery father_HBP father_diabetes ///
						ogtt_fasting_pw26 ogtt_2hour_pw26 sbp_avg dbp_avg GA days if visit == ``t''
			parmest, flist(est_list) saving(``y'_adj', replace) stars(0.05 0.01) idstr("regress `y' IVF adj") escal(N)
			}	
		}
	}
	
tempfile alldata
save `alldata'
/********************************************/
* dsconcat grabs all the datasets with names in the global macro regr
	dsconcat ${est_list}
* parse out the information in the id string
	gen model 		= word(idstr,1)
	gen outcome		= word(idstr,2)
	gen exposure	= word(idstr,3)
	gen adjust		= word(idstr,4)
	rename es_1 N
* create concatenated CI for output. I just like to do this because it makes it easier to copy into tables
	gen CI			= string(round(estimate,0.01),"%3.2f")+ " (" + string(round(min95,0.01),"%3.2f") + ///
					", " + string(round(max95,0.01),"%3.2f")+ ")" + stars
* create a main coefficient flag 
	gen main_coef = 0 if exposure == parm // reverse 1/0 order to allow ascending sort
	replace main_coef = 1 if main_coef == .
* label variables
	label var N "N"
	label var model "Regression Type"
	label var exposure "Exposure"
	label var main_coef "Primary exposure?"
	label var outcome "Outcome"
	label var CI "Estimate and 95% CI"
	label var adjust "Adjustment Model (Crude, Adj)"
* get data ready for output
	order model adjust outcome exposure N CI parm
* export output
	cd "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\IVF and Growth\"
	export excel using "Output/$S_DATE-IVF-dad_MI.xls" if main_coef==0 & adjust=="crude", sheet("crude") firstrow(variables) sheetreplace
	export excel using "Output/$S_DATE-IVF-dad_MI.xls" if main_coef==0 & adjust=="min", sheet("min") firstrow(variables) sheetreplace
	export excel using "Output/$S_DATE-IVF-dad_MI.xls" if main_coef==0 & adjust=="adj", sheet("adjusted") firstrow(variables) sheetreplace



// RESTORE ORIGINAL DATA
use `alldata', clear
