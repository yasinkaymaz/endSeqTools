#!/bin/bash

GTFfile=$1
UTRoutFile=${GTFfile%.gtf}_UTR.gtf

#GeneFeatureOutFile=$2
#Plus strand

awk '{
	if($3=="UTR" && ($7=="+" || $7=="-") && ($22 == "\"KNOWN\";") && ($20=="\"protein_coding\";" || $20=="\"snRNA\";" || $20=="\"lincRNA\";" || $20=="\"antisense\";") ){
		print $0
	}
}' $GTFfile |sort -k1,1 -k4,4g > $UTRoutFile




