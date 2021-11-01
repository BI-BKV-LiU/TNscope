#!/usr/bin/awk -f

#@include "genes_of_interest.awk"

BEGIN{
    OFS="\t";
	GoI = "GoI"

	goi["AKT1"] = GoI;
	goi["ALK"] = GoI;
	goi["BRAF"] = GoI;
	goi["CDKN2A"] = GoI;
	goi["CTNNB1"] = GoI;
	goi["DDR2"] = GoI;
	goi["EGFR"] = GoI;
	goi["ERBB2"] = GoI;
	goi["ESR1"] = GoI;
	goi["FGFR1"] = GoI;
	goi["FGFR2"] = GoI;
	goi["FGFR3"] = GoI;
	goi["FLT3"] = GoI;
	goi["GNA11"] = GoI;
	goi["GNAQ"] = GoI;
	goi["HRAS"] = GoI;
	goi["IDH1"] = GoI;
	goi["IDH2"] = GoI;
	goi["KIT"] = GoI;
	goi["KRAS"] = GoI;
	goi["MAP2K1"] = GoI;
	goi["MAP2K2"] = GoI;
	goi["MET"] = GoI;
	goi["MLH1"] = GoI;
	goi["MSH2"] = GoI;
	goi["MSH6"] = GoI;
	goi["PMS2"] = GoI;
	goi["NRAS"] = GoI;
	goi["NTRK1"] = GoI;
	goi["NTRK2"] = GoI;
	goi["NTRK3"] = GoI;
	goi["PDGFRA"] = GoI;
	goi["PIK3CA"] = GoI;
	goi["POLE"] = GoI;
	goi["PTEN"] = GoI;
	goi["RAF1"] = GoI;
	goi["RICTOR"] = GoI;
	goi["ROS1"] = GoI;
	goi["RET"] = GoI;
	goi["SMAD4"] = GoI;
	goi["STK11"] = GoI;
	goi["TERT"] = GoI;
	goi["BRCA1"] = GoI;
	goi["BRCA2"] = GoI;
	goi["APC"] = GoI;
	goi["PTCH1"] = GoI;
	goi["SUFU"] = GoI;
	goi["SMO"] = GoI;
	goi["MYC"] = GoI;
	goi["MYCN"] = GoI;
	goi["MDM2"] = GoI;
	goi["CDK4"] = GoI;
	goi["SMARCB1"] = GoI;
	goi["GNAS"] = GoI;
	goi["TSC2"] = GoI;
	goi["PALB2"] = GoI;
	goi["TP53"] = GoI;
	goi["ATM"] = GoI;
	goi["CHEK2"] = GoI;
	goi["RAD51C"] = GoI;
	goi["RAD51D"] = GoI;
	goi["BRIP1"] = GoI;
	goi["BAP1"] = GoI;
	goi["MUTYH"] = GoI;
	goi["DPYD"] = GoI;
	goi["TPMT"] = GoI;
	goi["UGT1A1"] = GoI;
	goi["CYP2D6"] = GoI;
}
{
	chr = $1
	start = $2
	end = $3
	symbol = $4

	if (symbol in goi){
		$0 = $0"_"goi[symbol]
		print $0
	} else {
		print $0
	}
}
