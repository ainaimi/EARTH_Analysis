options nofmterr;
%include 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_formats_macros\macros\macros.sas';
libname outcome 'S:\EOME\HauserEARTH\SQS2\SAS_files\master_datasets\2020';
libname EMR 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_emr';
%include 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_formats_macros\formats\Ruth BQ Formats.txt';
%include 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_formats_macros\formats\Ruth Female FQ Formats.txt';
%include 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_formats_macros\formats\Ruth moreBQ Formats.txt';
libname HERCULES 'E:\Emory\HERCULES\Mixtures Project';
run;

******************************************************************************************************************************************************;
******************2. Female CYCLE DATA****************************************************************************************************************;
******************************************************************************************************************************************************;

data femaleoutcomes; set outcome.mghwomen_cyclelevel_100120;
if cyclenum=. then cyclenum=99; 
if id="e3000" then delete; /*because no AFC outcomes*/
if id="d2428" then delete; /*because no AFC outcomes*/
if id="x2714" then delete; /*because no AFC outcomes*/
run; 
proc sort data=femaleoutcomes; by id cyclenum; run; 
data femaleoutcomes2; retain id proc_LuprStartDate FSHDay3_Date1; set femaleoutcomes; keep id cyclenum VISITDATE wbq_currstat; run;


/**************3. created a dataset "undermed" to identify women under medication**************************************/
data undermed; set femaleoutcomes; where proc_LuprStopDate ne .; keep id  proc_LuprStartDate proc_LuprStopDate; run; 
proc sort data=undermed; by ID proc_LuprStartDate; run;
data undermed;
	set undermed;
	by ID proc_LuprStartDate;
	if first.ID then Timepoint=1;
	else Timepoint+1;
	run;

proc transpose data=undermed
			   out=undermed1 (drop=_NAME_ drop=_LABEL_) prefix=proc_LuprStartDate;
	by id;
	id  timepoint;
	var proc_LuprStartDate;
run;
proc transpose data=undermed
			   out=undermed2 (drop=_NAME_ drop=_LABEL_) prefix=proc_LuprStopDate;
	by id;
	id  timepoint;
	var proc_LuprStopDate;
run;
data undermed_final;
merge undermed1 undermed2;
by id;
run;

/************5 import AFC dataset**************************************/
data allantral; set EMR.fhorm_afcleftright_071420; 
where FHORM_AFScanAvail=1;
AFScanAvail=FHORM_AFScanAvail;
AFScanDate=FHORM_AFScanDate;
AFC_L=FHORM_AntFollCountLeft;
AFC_R=FHORM_AntFollCountRight; 
drop FHORM_AFScanAvail  FHORM_AFScanDate  FHORM_AntFollCountLeft FHORM_AntFollCountRight;
run; 

proc sort data=allantral ; by id AFScanDate; run;
proc sort data=undermed_final; by id;run; 
data AFCmed; 
merge allantral(in=a)  undermed_final;
by id;
if a;
format AFScanDate MMDDYY8.;
run;

/**********6. identify the duplicate records***********/
proc sort data=AFCmed; by id AFScanDate; run; 
data afcmed_dup; set AFCmed;
format AFScanDate MMDDYY8.;
duplicate + 1;
  by id AFScanDate;
  if first.AFScanDate then duplicate = 1;
run;

data afcmed_nodup;retain id AFScanDate ; set afcmed_dup;
where duplicate=1;
run; /*1022 available scans, no duplicates*/

proc sort data=afcmed_nodup; by id;
data women_only; set afcmed_nodup; by id; if first.id; run;
/*869 unique women*/

/***********7  kick out AFC data under medication
Irene:If any of the patients whose AFC was done on lupron have another AFC that was done within 3 months not on lupron, you could consider keeping those.
************/
data afc_nomed; set afcmed_nodup;
format AFScanDate MMDDYY8.;
format proc_LuprStartDate1 MMDDYY8.;
format proc_LuprStopDate1 MMDDYY8.;
format proc_LuprStartDate2 MMDDYY8.;
format proc_LuprStopDate2 MMDDYY8.;
format proc_LuprStartDate3 MMDDYY8.;
format proc_LuprStopDate3 MMDDYY8.;
format proc_LuprStartDate4 MMDDYY8.;
format proc_LuprStopDate4 MMDDYY8.;
format proc_LuprStartDate5 MMDDYY8.;
format proc_LuprStopDate5 MMDDYY8.;
format proc_LuprStartDate6 MMDDYY8.;
format proc_LuprStopDate6 MMDDYY8.;
AFC_medstart1=AFScanDate-proc_LuprStartDate1;
AFC_medstop1=AFScanDate-proc_LuprStopDate1;
AFC_medstart2=AFScanDate-proc_LuprStartDate2;
AFC_medstop2=AFScanDate-proc_LuprStopDate2;
AFC_medstart3=AFScanDate-proc_LuprStartDate3;
AFC_medstop3=AFScanDate-proc_LuprStopDate3;
AFC_medstart4=AFScanDate-proc_LuprStartDate4;
AFC_medstop4=AFScanDate-proc_LuprStopDate4;
AFC_medstart5=AFScanDate-proc_LuprStartDate5;
AFC_medstop5=AFScanDate-proc_LuprStopDate5;
AFC_medstart5=AFScanDate-proc_LuprStartDate6;
AFC_medstop5=AFScanDate-proc_LuprStopDate6;
if 0<AFC_medstart1<30  or -1<AFC_medstop1<30 or
0<AFC_medstart2<30  or -1<AFC_medstop2<30  or
0<AFC_medstart3<30  or -1<AFC_medstop3<30 or
0<AFC_medstart4<30  or -1<AFC_medstop4<30 or
0<AFC_medstart5<30  or -1<AFC_medstop5<30 or
0<AFC_medstart6<30  or -1<AFC_medstop6<30 then luprafc=1; else luprafc=0;
run;

data afc_nomed_2;
set AFC_nomed;
if AFScanAvail =1 then AFS_yes=1; else AFS_yes=0; 
if AFC_L=77 then PCOS_L=1; else PCOS_L=0;
if AFC_R=77 then PCOS_R=1; else PCOS_R=0;
if AFC_L=78 then diffsee_L=1;else diffsee_L=0;
if AFC_R=78 then diffsee_R=1;else diffsee_R=0;
if AFC_L=79 then noovary_L=1;else noovary_L=0;
if AFC_R=79 then noovary_R=1;else noovary_R=0;
if AFC_L in (.,99) then afcmiss_L=1; else afcmiss_L=0;
if AFC_R in (.,99) then afcmiss_R=1; else afcmiss_R=0;
AFC_t=AFC_L+AFC_R;
run; 

data afc_nomed_3; set afc_nomed_2; if luprafc=1 then delete;run; /*980 AF Scans*/

*Information on exclusions for Supplemental Figure 1*;
proc freq data=afc_nomed_3; table noovary_L*noovary_R; run;
proc freq data=afc_nomed_3; table PCOS_L*PCOS_R; run;
proc freq data=afc_nomed_3; table diffsee_L*diffsee_R; run;
proc freq data=afc_nomed_3; table afcmiss_L*afcmiss_R; run;

/*******8. kick out AFC =PCOS, difficult to visualized or oopherectomy or missing*/
data afc_nomed_4; set afc_nomed_3; 
if noovary_L=1 or noovary_R=1 then delete; 
if PCOS_L=1 or PCOS_R=1 or diffsee_L=1 or diffsee_R=1  then delete;
if afcmiss_L=1 or afcmiss_R=1 then delete; 
run; /*882 scans*/

proc sort data=afc_nomed_4; by id AFScanDate;
data afc_nomed_4_N; set afc_nomed_4; by id AFScanDate; if first.id; 
keep id AFScanDate AFC_L AFC_R AFC_T; run; /*775 women*/

proc sort data=femaleoutcomes; by id;
data alldata;
set femaleoutcomes;
by id;
if first.id;
if smokstat=1 then eversmok=0; else if smokstat=2 then eversmok=1;  else if smokstat=3 then eversmok=1;
if educ=. then educ1=1; *Anyone missing education is set to the most popular category*;
else if educ=1 or educ=2 or educ=3 or educ=4 then educ1=0; else if educ=5 then educ1=1; else educ1=2; 
enroll_year=year(visitdate);
if wbq_everpregnant=1 then gravid=1; else gravid=0;
if wbq_previnftexam=1 then inftexam=1; else inftexam=0;
if wbq_inftrt=1 then infttrt=1; else infttrt=0;
if WBQ_CYCLELENGTH=888 or WBQ_CYCLELENGTH=. then abnormalcycle=9; 
else if WBQ_CYCLELENGTH<24 or WBQ_CYCLELENGTH>38 then abnormalcycle=1; else abnormalcycle=0;
run;

data covariates;
set alldata;
keep id idnum idnew visitdate enroll_year bmi age educ1 smokstat eversmok visitdate wbq_brthdate e2_triggerresult races 
infdx1 white sartnew2 previousIUI previousIVF gravid inftexam infttrt FSHDay3_Date1 FSHDay3Result1 
WBQ_CYCLELENGTH WBQ_CYCLEREG abnormalcycle;
run;

proc sort data=covariates; by id;
proc sort data=afc_nomed_4_N; by id;
data AFC_final;
merge afc_nomed_4_N(in=b) covariates;
by id;
if b;
run;

proc univariate data=AFC_final;
var AFC_t;
histogram;
run;
/*775 women with one AFC scan (no repeated scans)*/


libname pht 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_phthalates';
proc sort data = pht.phthalates121517; by id; run;
proc sort data=AFC_final; by id;
data AFC_pht;
merge AFC_final(in=b) pht.phthalates121517;
by id;
daydiff=AFScanDate-sampledate;
if daydiff=. then delete;
daydiff_abs=abs(daydiff); 
run;

proc sort data=AFC_pht; by id daydiff_abs;
data AFC_pht_N; set AFC_pht; by id; 
dehp=(mehp*(1/278.34))+(mehhp*(1/294.34))+(meohp*(1/292.33))+(mecpp*(1/308.33));
sgratio_pht=((1.015-1)/(urine_sg-1));
phtsampledate=sampledate; 
urine_sg_pht=urine_sg;
format phtsampledate DATE9.;
if first.id; 
run; 
data pht_AFC; set AFC_pht_N; keep id phtsampledate urine_sg_pht sgratio_pht 
dehp mehp mehhp meohp mecpp MEP1 MBP mCPP miBP mBZP1 mCOP mCNP;
run; 


/*683 women have at least one phthalate measurement
  only women with phthalates measured prior to AFC ---> 218 women
  only women with phthalates measured prior to and up to a month after AFC ----> 371 women
  only women with phthalates measured within a month prior to or after AFC ----> 309 women;*/

libname bpa 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_phenols';
proc contents data=bpa.phenols013018; run;
proc means data=bpa.phenols013018 n nmiss min q1 median q3 max; 
var bpa bpf bps bp_3 b_pb m_pb p_pb tcs; run;
proc sort data = bpa.phenols013018; by id; run;
proc sort data=AFC_final; by id;
data AFC_bpa;
merge AFC_final(in=b) bpa.phenols013018;
by id;
daydiff=AFScanDate-sampledate;
if daydiff=. then delete;
daydiff_abs=abs(daydiff); 
run;

proc sort data=AFC_bpa; by id daydiff_abs;
data AFC_bpa_N; set AFC_bpa; by id daydiff_abs; 
sgratio_bpa=((1.015-1)/(urine_sg-1));
bpasampledate=sampledate; 
urine_sg_bpa=urine_sg;
format bpasampledate DATE9.;
if first.id; 
run; 
data bpa_AFC; set AFC_bpa_N; keep id bpasampledate urine_sg_bpa sgratio_bpa 
bpa bpf bps bp_3 b_pb m_pb p_pb tcs;
run; 
/*683 women have at least one bpa measurement*;*/


libname opfr 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_flame_retardants'; 
proc contents data=opfr.pfr_metabolites_170301; run;
proc means data=opfr.pfr_metabolites_170301 n nmiss min q1 median q3 max; 
var BDCIPP DPHP ipPPP; run;
proc sort data=opfr.pfr_metabolites_170301; by id; run;
proc sort data=AFC_final; by id;
data AFC_opfr;
merge AFC_final(in=b) opfr.pfr_metabolites_170301;
by id;
daydiff=AFScanDate-sampledate;
if daydiff=. then delete;
daydiff_abs=abs(daydiff); 
run;

proc sort data=AFC_opfr; by id daydiff_abs;
data AFC_opfr_N; set AFC_opfr; by id;
sgratio_opfr=((1.015-1)/(sg-1));
opfrsampledate=sampledate; 
urine_sg_opfr=sg;
format opfrsampledate DATE9.;
if first.id; 
run; 
data opfr_AFC; set AFC_opfr_N; keep id opfrsampledate urine_sg_opfr sgratio_opfr 
BDCIPP DPHP ipPPP;
run; 
/*196 women have at least one OPFR measurement*;*/

libname hg 'S:\EOME\HauserEARTH\SQS2\SAS_files\sysfiles_metals'; 
proc contents data=hg.hair_hg_021821; run;
proc means data=hg.hair_hg_021821 n nmiss min q1 median q3 max; var hg; run;
proc sort data=hg.hair_hg_021821; by id; run;
proc sort data=AFC_final; by id;
data AFC_hg;
merge AFC_final(in=b) hg.hair_hg_021821;
by id;
daydiff=AFScanDate-sampledate;
if daydiff=. then delete;
daydiff_abs=abs(daydiff); 
run;
proc sort data=AFC_hg; by id daydiff_abs;
data AFC_hg_N; set AFC_hg; by id daydiff_abs; if first.id; hgsampledate=sampledate; format hgsampledate DATE9.; run; 
/*700 women have at least one hair Hg measurement*;*/

data hg_AFC; set AFC_hg_N; keep id hgsampledate hg; run;


********************MERGING ALL OF THE CHEMICALS INTO ONE DATASET WITH AFC (one line per woman)**********************;
proc sort data=pht_AFC; by id; proc sort data=bpa_AFC; by id; proc sort data=opfr_AFC; by id; proc sort data=hg_AFC; by id; 
proc sort data=AFC_final; by id;
data HERCULES.AFC_mixtures;
merge pht_AFC bpa_AFC opfr_AFC hg_AFC AFC_final;
by id;
run;

proc corr data=AFC_mixtures spearman;
var dehp mehp mehhp meohp mecpp MEP1 MBP mCPP miBP mBZP1 mCOP mCNP
bpa bpf bps bp_3 b_pb m_pb p_pb tcs
BDCIPP DPHP ipPPP hg;
run;
