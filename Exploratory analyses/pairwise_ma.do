clear all
version 16
set more off

*Set working directory and load data
cd "R:\SCFreeman\scf20\Survival\Data"
import delimited melanoma_ipd.csv, clear
drop v1

* Set the colour scheme
set scheme s1color

*Declare survival data
stset time, fail(event=1)

* Pairiwse meta-analysis of ipilimumab versus nivolumab plus ipilimumab
preserve
	keep if studycode==4 | studycode==5
	keep if txcode==8 | txcode==10
	
	gen trt=0 if txcode==8
	replace trt=1 if txcode==10
	
	ipdmetan, study(study) eform effect("HR") forestplot(xlabel(0.4(0.2)1.4)): stcox i.trt
restore	