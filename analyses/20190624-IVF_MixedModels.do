use "C:\Users\JHUANGYH\Dropbox\00 - Singapore\GUSTO\01 - PROJECTS\DATA\ANALYTIC\20190125-IVF_Gen_MI_Drop.dta", clear

/****************************************/
// ZLEN
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed zlen IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins zlen
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed zlen IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted HFAZ") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted mean height-for-age Z-Scores (HFAZ), by IVF status.")

/****************************************/
// HEIGHT
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed height IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins height
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed height IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted height (cm)") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted mean height (cm), by IVF status.")

/****************************************/
// ZWEI
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed zwei IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins zwei
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed zwei IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted WFAZ") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted weight-for-age Z-Score (WFAZ), by IVF status.")

/****************************************/
// WEIGHT
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed weight IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins weight
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed weight IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted weight (kg)") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted weight (kg), by IVF status.")	
	
/****************************************/
// ZBMI
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed zbmi IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins zbmi
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed zbmi IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted ZBMI") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted BMI-for-age Z-score (ZBMI), by IVF status.")	
	
/****************************************/
// BMI
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed bmi IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins bmi
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed bmi IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted BMI") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted BMI (kg/m{superscript:2}), by IVF status.")		
	
	
/****************************************/
// ZSBP
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed ZSBP IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins ZSBP
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed ZSBP IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted ZSBP") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted standardized systolic blood pressure (ZSBP), by IVF status.")	


/****************************************/
// ZDBP
/****************************************/
set matsize 5000
// Estimate on MI set
mi estimate: xtmixed ZDBP IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days, cov(un) var reml
mi estimate, saving(out, replace): emargins ZDBP
// Save marginal estimates
mat b = e(b_mi)
mat V = e(V_mi)
// re-rerun "dummy" margins in non-imputed set
	quietly: xtmixed ZDBP IVF##visit mother_age_delivery i.MOM_EDU i.ETHNIC i.HH_INC ///
						m_height_pw26 ppBMI f_height_m24 f_weight_m24 MALE SCORE parity i.smk_home || SubjectID: days if _mi_m == 0, cov(un) var reml
	quietly: margins IVF#visit
// re-post margins estimates from MI set
	myret
// plot as if non-MI
mata: st_global("e(cmd)", "margins")

marginsplot, x(visit) recast(line) recastci(rarea) ytitle("predicted ZDBP") xtitle("visit month") ///
	scheme(plotplain) ///
	ci1opts(color(blue%25) lwidth(none)) plot1opts(lcolor(blue)) ///
	ci2opts(color(red%25) lwidth(none)) plot2opts(lpattern(dash) lcolor(red)) ///
	plot(, labels("SC" "IVF")) legend(pos(6) rows(1)) ///
	title("Predicted standardized diastolic blood pressure (ZDBP), by IVF status.")	
	
	
/****************************************/
/****************************************/	
// PROGRAMS
/****************************************/
capture program drop emargins
program emargins `1', eclass properties(mi)
	xtmixed `1' IVF##visit || SubjectID: days, cov(un) var reml
	margins IVF#visit, post
end

capture program drop myret
program myret, rclass
	return add
	return matrix b = b
	return matrix V = V
end
