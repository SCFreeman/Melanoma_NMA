clear all
version 16
set more off

* Set the working directory and import data
cd "R:\SCFreeman\scf20\Melanoma_NMA\Data"
import delimited melanoma_ipd.csv, clear
drop v1

* set the colour scheme
set scheme s1color

* Declare survival data
stset time, fail(event=1)

*log-log plots
stphplot if studycode==1, by(txcode) legend(label(1 "Dacarbazine") label(2 "Dabrafenib")) title("BREAK-3")
stphplot if studycode==2, by(txcode) legend(label(1 "Dacarbazine") label(2 "Vemurafenib")) title("BRIM-3")
stphplot if studycode==3, by(txcode) legend(label(1 "Dacarbazine") label(2 "Nivolumab")) title("CheckMate 066")
stphplot if studycode==4, by(txcode) legend(label(1 "Ipilimumab") label(2 "Nivolumab") label(3 "Ipilimumab + Nivolumab")) title("CheckMate 067")
stphplot if studycode==5, by(txcode) legend(label(1 "Ipilimumab") label(2 "Nivolumab + Ipilimumab")) title("CheckMate 069")
stphplot if studycode==6, by(txcode) legend(label(1 "Vemurafenib") label(2 "Vemurafenib + Cobimetinib")) title("COBRIM")
stphplot if studycode==7, by(txcode) legend(label(1 "Dabrafenib") label(2 "Dabrafenib + Trametinib")) title("COMBI-d")
stphplot if studycode==8, by(txcode) legend(label(1 "Dabrafenib + Trametinib") label(2 "Vemurafenib")) title("COMBI-v")
stphplot if studycode==9, by(txcode) legend(label(1 "Ipilimumab") label(2 "Ipilimumab + Sargramostin")) title("Hodi 2014")
stphplot if studycode==10, by(txcode) legend(label(1 "Ipilimumab") label(2 "Pemborlizumab")) title("Keynote 006")
stphplot if studycode==11, by(txcode) legend(label(1 "Dacarbazine") label(2 "Tremelimumab")) title("Ribas 2013")
stphplot if studycode==12, by(txcode) legend(label(1 "Dacarbazine") label(2 "Ipilimumab + Dacarbazine")) title("Robert 2011")
stphplot if studycode==13, by(txcode) legend(label(1 "Dacarbazine") label(2 "Selumetinib + Dacarbazine")) title("Robert 2013")