clear 
clear matrix

set memory 20m
set maxvar 10000

adopath + "C:\ado\CDC igrowup_stata"
adopath + "C:\ado\who2007_stata" 

use "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\DATA\ANALYTIC\20180904-Anthro.dta", clear
cd "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\DATA\ANALYTIC\"


gen str80 reflib = "C:\ado\CDC igrowup_stata"
*gen str80 reflib = "C:\ado\who2007_stata"
lab var reflib "Directory of reference tables"

gen str100 datalib = "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\DATA\ANALYTIC\Zscores"
lab var datalib "Directory for datafiles"

gen str30 datalab = "WHOout"
lab var datalab "Working file"

gen ageunit = "days"
gen str1 measure = " "
gen str1 oedema = "n"
gen sw = 1

replace sex = substr(sex,1,1)
capture gen c_age_birth = 0

save "20180904-Anthro.dta", replace			

* WHO 2006 standards for Age 0-5
tokenize birth day1 wk3 m3 m6 m9 m12 m15 
forvalues i = 1/8 {
	use "20180904-Anthro.dta", clear
	set more off
	adopath + "C:\ado\CDC igrowup_stata"
	replace datalab = "WHOout_``i''"
	replace measure = "l"
	igrowup_restricted reflib datalib datalab sex c_age_``i'' ageunit c_weight_``i'' c_length_``i'' measure oedema sw 
	}
	
tokenize m18 m24 m36 yr4 yr4_5 yr5
forvalues i = 1/6 {
	use "20180904-Anthro.dta", clear
	set more off
	adopath + "C:\ado\CDC igrowup_stata"
	replace datalab = "WHOout_``i''"
	replace measure = " "
	if "``i''" == "m18" | "``i''" == "m24" {
		replace c_length_``i'' = c_height_``i'' + 0.7 if c_length_``i'' == . & c_height_``i'' != .
		replace c_height_``i'' = c_length_``i''
		replace measure = "l"
		}
	igrowup_restricted reflib datalib datalab sex c_age_``i'' ageunit c_weight_``i'' c_height_``i'' measure oedema sw 
	}

* WHO 2007 Standards for AGE 5.5 & 6
use "20180904-Anthro.dta", clear
replace reflib = "C:\ado\who2007_stata"
replace datalab = "WHOout_yr5_5"
who2007 reflib datalib datalab sex c_age_yr5_5 ageunit c_weight_yr5_5 c_height_yr5_5 oedema sw
replace datalab = "WHOout_yr6"
who2007 reflib datalib datalab sex c_age_yr6 ageunit c_weight_yr6 c_height_yr6 oedema sw

* TABULATE birth parameters
use RNA_sub_prep_12.dta, clear
gen c_bmi_00 = c_weight_00 / (c_length_00/100)^2
tabstat c_length_00 c_weight_00 c_bmi_00 c_head_00, col(stat) stat(N mean sd p25 med p75 min max)

* TABULATE 0.03 - 15 mo, 18 - 60 mo
tokenize 0_03 01 03 06 09 12 15 // m18 m24 m36 yr4 yr4_5 yr5
forvalues i = 1/7 {
	di "---------------------------------------------------------------------------------"
	di " ``i'': "
	use "Output\WHOout_``i''_z_st.dta", clear
	/*quietly {
		gen ID = substr(PSCID,5,5)
		destring ID, replace
		gen SUBSET = inlist(ID, ///
				18006, 04188, 04314, 50118, 66138, 21409, 21494, 21575, 04026, ///
				21646, 21650, 21921, 21990, 50135, 20486, 20537, 20640, 20951, ///
				20993, 21308, 20515, 20556, 21153, 21353, 21478, 21614, 21640, ///
				21784, 08044, 74131, 20839, 20899, 21578, 21931, 04214, 21956, ///
				20501, 20694, 20861, 21288, 21471, 21796, 21951, 22106, 66123, ///
				74140, 20638)
		}*/
		tabstat _cbmi _zbmi _zwfl _zlen _zac _zts _zss, col(stat) stat(N mean sd p25 med p75 min max)
		di "WFL:"
		count if _zwfl < -1
		count if _zwfl < -1 & SUBSET == 1
		count if _zwfl < -2
		count if _zwfl < -2 & SUBSET == 1
		di "LFA:"
		count if _zlen < -1
		count if _zlen < -1 & SUBSET == 1
		count if _zlen < -2
		count if _zlen < -2 & SUBSET == 1
		di "MUAC:"
		count if _zac < -1
		count if _zac < -1 & SUBSET == 1
		count if _zac < -2
		count if _zac < -2 & SUBSET == 1
		di "BMI: "
		count if _zbmi > 1 & _zbmi != .
		count if _zbmi > 1 & SUBSET == 1 & _zbmi != .
		count if _zbmi > 2 & _zbmi != .
		count if _zbmi > 2 & SUBSET == 1 & _zbmi != .
	di _newline _newline
	}

tokenize yr5_5 yr6
forvalues i = 1/2 {
	di "---------------------------------------------------------------------------------"
	di " ``i'': "
	use "Output\WHOout_``i''_z.dta", clear
	quietly {
		gen ID = substr(PSCID,5,5)
		destring ID, replace
		gen SUBSET = inlist(ID, ///
				18006, 04188, 04314, 50118, 66138, 21409, 21494, 21575, 04026, ///
				21646, 21650, 21921, 21990, 50135, 20486, 20537, 20640, 20951, ///
				20993, 21308, 20515, 20556, 21153, 21353, 21478, 21614, 21640, ///
				21784, 08044, 74131, 20839, 20899, 21578, 21931, 04214, 21956, ///
				20501, 20694, 20861, 21288, 21471, 21796, 21951, 22106, 66123, ///
				74140, 20638)
		}
	tabstat _cbmi _zbfa _zhfa c_midarm_``i'', col(stat) stat(N mean sd p25 med p75 min max)
	di "LOW BMI:"
	count if _zbfa < -1
	count if _zbfa < -1 & SUBSET == 1
	count if _zbfa < -2
	count if _zbfa < -2 & SUBSET == 1
	di "HFA:"
	count if _zhfa < -1
	count if _zhfa < -1 & SUBSET == 1
	count if _zhfa < -2
	count if _zhfa < -2 & SUBSET == 1
	di "BMI: "
	count if _zbfa > 1 & _zbfa != .
	count if _zbfa > 1 & SUBSET == 1 & _zbfa != .
	count if _zbfa > 2 & _zbfa != .
	count if _zbfa > 2 & SUBSET == 1 & _zbfa != .
	di _newline _newline
	}

* Hypertension classification
* classify height
centile c_height_72 if sex == "Male", centile(5 10 25 50 75 90 95)
forvalues i = 1/7 {
	scalar m_`i' = r(c_`i')
	}
centile c_height_72 if sex == "Female", centile(5 10 25 50 75 90 95)
forvalues i = 1/7 {
	scalar f_`i' = r(c_`i')
	}
* enter BP thresholds
scalar m_s_PRE_1 = 105 
scalar m_s_PRE_2 = 106 
scalar m_s_PRE_3 = 108 
scalar m_s_PRE_4 = 110 
scalar m_s_PRE_5 = 111 
scalar m_s_PRE_6 = 113 
scalar m_s_PRE_7 = 113 

scalar m_s_HTN_1 = 109 
scalar m_s_HTN_2 = 110 
scalar m_s_HTN_3 = 112 
scalar m_s_HTN_4 = 114 
scalar m_s_HTN_5 = 115 
scalar m_s_HTN_6 = 117 
scalar m_s_HTN_7 = 117
	
scalar m_d_PRE_1 = 68
scalar m_d_PRE_2 = 68 
scalar m_d_PRE_3 = 69 
scalar m_d_PRE_4 = 70
scalar m_d_PRE_5 = 71
scalar m_d_PRE_6 = 72
scalar m_d_PRE_7 = 72 

scalar m_d_HTN_1 = 72
scalar m_d_HTN_2 = 72
scalar m_d_HTN_3 = 73
scalar m_d_HTN_4 = 74
scalar m_d_HTN_5 = 75
scalar m_d_HTN_6 = 76
scalar m_d_HTN_7 = 76 	
	
scalar f_s_PRE_1 = 104
scalar f_s_PRE_2 = 105
scalar f_s_PRE_3 = 106
scalar f_s_PRE_4 = 108
scalar f_s_PRE_5 = 109
scalar f_s_PRE_6 = 110
scalar f_s_PRE_7 = 111 

scalar f_s_HTN_1 = 108
scalar f_s_HTN_2 = 109
scalar f_s_HTN_3 = 110
scalar f_s_HTN_4 = 111
scalar f_s_HTN_5 = 113
scalar f_s_HTN_6 = 114
scalar f_s_HTN_7 = 115 

scalar f_d_PRE_1 = 68
scalar f_d_PRE_2 = 68
scalar f_d_PRE_3 = 69
scalar f_d_PRE_4 = 70
scalar f_d_PRE_5 = 70
scalar f_d_PRE_6 = 71
scalar f_d_PRE_7 = 72

scalar f_d_HTN_1 = 72
scalar f_d_HTN_2 = 72
scalar f_d_HTN_3 = 73
scalar f_d_HTN_4 = 74
scalar f_d_HTN_5 = 74
scalar f_d_HTN_6 = 75
scalar f_d_HTN_7 = 76 
	
gen ht_cent = .
forvalues i = 1/7 {
	replace ht_cent = `i' if (sex == "Male" & c_height_72 > m_`i') | (sex == "Female" & c_height_72 > f_`i')
	}
replace ht_cent = . if c_height_72 == .

gen PreHTN = 0
gen HTN = 0

forvalues i = 1/7 {
	replace PreHTN = 1 if sex == "Male" & ht_cent == `i' & ((q1_3_cardio_sbp_y6 > (m_s_PRE_`i' + 5)) | (q1_4_cardio_dbp_y6 > (m_d_PRE_`i' + 5)))
	replace HTN = 1 if sex == "Male" & ht_cent == `i' & ((q1_3_cardio_sbp_y6 > (m_s_HTN_`i' + 5)) | (q1_4_cardio_dbp_y6 > m_d_HTN_`i' + 5))
	replace PreHTN = 1 if sex == "Female" & ht_cent == `i' & ((q1_3_cardio_sbp_y6 > (f_s_PRE_`i' + 5)) | (q1_4_cardio_dbp_y6 > (f_d_PRE_`i' + 5)))
	replace HTN = 1 if sex == "Female" & ht_cent == `i' & ((q1_3_cardio_sbp_y6 > (f_s_HTN_`i' + 5)) | (q1_4_cardio_dbp_y6 > (f_d_HTN_`i' + 5)))
	}
replace PreHTN = . if (q1_3_cardio_sbp_y6 == . & q1_4_cardio_dbp_y6 == .)
replace HTN = . if (q1_3_cardio_sbp_y6 == . & q1_4_cardio_dbp_y6 == .)

gen OVER = 0
replace OVER = 1 if (q1_3_cardio_sbp_y6 > 120 | q1_4_cardio_dbp_y6 > 80)
replace OVER = . if (q1_3_cardio_sbp_y6 == . & q1_4_cardio_dbp_y6 == .)

sum PreHTN HTN OVER
sum PreHTN HTN OVER if SUBSET == 1
foreach var of varlist PreHTN HTN OVER {
	di "`var':" _
	count if `var' == 1
	di "`var' (SUBSET):" _
	count if `var' == 1 & SUBSET == 1
	}

tabstat q1_3_cardio_sbp_y6 q1_4_cardio_dbp_y6 q1_5_cardio_map_y6 q1_6_cardio_pp_y6 q1_7_cardio_hr_y6, ///
	col(stat) stat(N mean sd p25 med p75 min max)
tabstat q1_3_cardio_sbp_y6 q1_4_cardio_dbp_y6 q1_5_cardio_map_y6 q1_6_cardio_pp_y6 q1_7_cardio_hr_y6 if SUBSET == 1, ///
	col(stat) stat(N mean sd p25 med p75 min max)
