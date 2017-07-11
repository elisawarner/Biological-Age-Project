LIBNAME nhanes "M:\NHANES III\Original Dataset\";

PROC SORT DATA = NHANES.PARTI;
	BY AGE;
RUN;

DATA nhanes.PARTI;
	SET nhanes.PARTI;
		LABEL SEQN = "SEQN"
			  BUP = "Urea Nitrogen"
			  TCP = "Total cholesterol"
			  CRP = "C-reactive protein"
			  CEP = "Creatinine"
			  APPSI = "Alkaline phosphatase"
			  AMP = "Albumin"
			  GHP = "Glycated hemoglobin"
			  AGE = "Age"
			  samp_wt = "Sample Weight"
			  CMVOD = "Cytomegalovirus OD"
			  FEV = "Forced expiratory volume"
			  WBC = "White blood cell count"
			  SBP = "Systolic blood pressure"
			  FVC = "Forced Vital Capacity"
			  FEV_FVC = "Forced Vital Capacity Ratio"
			  FEMALE = "Female? 0=No 1=Yes"
			  ;
RUN;

*Split age into categories;
DATA PARTI;
	SET nhanes.PARTI;
	IF AGE > 29 and AGE < 60 THEN age_cat = 0;
	ELSE IF AGE >= 60 and AGE < 76 THEN age_cat = 1;
	ELSE age_cat = 2;
RUN;

PROC FREQ DATA = PARTI;
	TABLES age_cat;
RUN;

PROC PRINT DATA = PARTI; 
	WHERE age_cat = 2;
	ID SEQN;
	VAR AGE;
RUN;

*Double check means;
PROC MEANS DATA = PARTI;
	VAR AGE CRP CEP GHP AMP TCP CMVOD APPSI FEV BUP SBP;
	WHERE AGE >= 30 and AGE <= 75;
	WEIGHT SAMP_WT;
*	CLASS age_cat;
RUN;

PROC MEANS DATA = nhanes.PARTI;
	VAR AGE;
	WHERE AGE < 76;
	WEIGHT samp_wt;
RUN;

/****** MLR METHOD ******/

DATA FINAL;
	SET NHANES.PARTI;
	WHERE AGE >= 30 AND AGE <= 75;
RUN;

*Regression Females (Removed CMV data);
ods graphics on;
PROC GLM DATA = final ORDER=FORMATTED PLOTS=(DIAGNOSTICS RESIDUALS);
	MODEL AGE = CRP CEP GHP AMP TCP CMVOD APPSI FEV BUP SBP / solution clparm;
	WEIGHT samp_wt;
	TITLE "Regression for Females";
	WHERE FEMALE = 1;
RUN;
QUIT;
ods graphics off;

*Regression Males;
ods graphics on;
PROC GLM DATA = final ORDER=FORMATTED PLOTS=(DIAGNOSTICS RESIDUALS);
	MODEL AGE = CRP CEP GHP AMP TCP CMVOD APPSI FEV BUP SBP / solution clparm;
	TITLE "Regression for Males";
	WEIGHT samp_wt;
	WHERE FEMALE = 0;
RUN;
QUIT;
ods graphics off;

*Create BA variable (BA MLR);
DATA final;
	SET final;
	IF FEMALE = 1
		THEN BAM = 33.89-(1.232*CRP)-(0.542*CEP)+(0.187*GHP)-(1.549*AMP)+(0.047*TCP)-(0.606*CMVOD)+(0.005*APPSI)-(0.008*FEV)+(0.529*BUP)+(0.196*SBP);
	ELSE BAM = 54.52-(0.289*CRP)+(0.810*CEP)+(0.634*GHP)-(6.357*AMP)+(0.016*TCP)+(1.096*CMVOD)-(0.0097*APPSI)-(0.006*FEV)+(0.508*BUP)+(0.206*SBP);
RUN;

*Get average BA and CA vars to standardize;
PROC MEANS DATA = final;
	VAR BAM; OUTPUT OUT = mean1 MEAN = avg_BA; 
	VAR AGE;
RUN;

*Create standardized BA;
*In Levine's dataset, I did not standardize because my mean was within .2 of the CA sample mean;
DATA final;
	SET final;
	BAZscore = (BA - 51.2723050)/10.2470791;
	BAstandardized = (BAZscore * 13.6218507)+49.8870806;
RUN;

PROC UNIVARIATE DATA = final;
	VAR BA HSAGEIR BAstandardized;
	histogram;
	TITLE "";
RUN;

PROC SGSCATTER DATA = final;
	PLOT BA * HSAGEIR;
RUN;


/****** KLEMERA AND DOUBAL METHOD ******/

*Men;

DATA finalkm;
	SET nhanes.partI;
	WHERE AGE >= 30 and AGE <= 75 and FEMALE = 0;
RUN; 

/*
*FINDING VALUES FOR VARIABLES (Run SLR);
*Macro for j=CA;
%MACRO SLR (j);
	ODS GRAPHICS ON;
	PROC GLM DATA=finalkm ORDER=FORMATTED PLOTS=(DIAGNOSTICS RESIDUALS);
		MODEL &j = AGE / solution clparm;
		WEIGHT samp_wt;
		TITLE &J;
	RUN;
	QUIT;
	ODS GRAPHICS OFF;
%MEND SLR;

%SLR(CRP)
%SLR(CEP)
%SLR(GHP)
%SLR(AMP)
%SLR(TCP)
%SLR(CMVOD)
%SLR(APPSI)
%SLR(FEV)
%SLR(BUP)
%SLR(SBP)*/

DATA finalkm;
	SET finalkm;

	Agemax = 75;
	Agemin = 30;
	n = 4359;
	m = 10; *number of covariates;

	CRP_k = 0.0057474554; *k slope;
	CRP_s = sqrt(0.298914); *s MSE;
	CRP_q = 0.0838117260; *q Intercept;
	CRP_x = (CRP-CRP_q)*(CRP_k/(CRP_s**2)); *numerator;
	CRP_d = ((CRP_k / CRP_s)**2); *denominator;
	CRP_r = sqrt(0.021917); *r from r2;
	CRP_rn = (CRP_r**2)/(sqrt(1-(CRP_r**2))); *numerator of rchar;
	CRP_rd = (CRP_r)/(sqrt(1-(CRP_r**2))); *denominator of rchar;

	CEP_k = 0.003244968;
	CEP_s = sqrt(0.0762434);
	CEP_q = 1.026886232;
	CEP_x = (CEP-CEP_q)*(CEP_k/(CEP_s**2));
	CEP_d = ((CEP_k / CEP_s)**2);
	CEP_r = sqrt(0.027241);
	CEP_rn = (CEP_r**2)/(sqrt(1-(CEP_r**2)));
	CEP_rd = (CEP_r)/(sqrt(1-(CEP_r**2)));

	GHP_k = 0.016103234;
	GHP_s = sqrt(0.981220);
	GHP_q = 4.694422950;
	GHP_x = (GHP-GHP_q)*(GHP_k/(GHP_s**2));
	GHP_d = ((GHP_k / GHP_s)**2);
	GHP_r = sqrt(0.050862);
	GHP_rn = (GHP_r**2)/(sqrt(1-(GHP_r**2)));
	GHP_rd = (GHP_r)/(sqrt(1-(GHP_r**2)));

	AMP_k = -0.008042133;
	AMP_s = sqrt(0.1275978);
	AMP_q = 4.635353910;
	AMP_x = (AMP-AMP_q)*(AMP_k/(AMP_s**2));
	AMP_d = ((AMP_k / AMP_s)**2);
	AMP_r = sqrt(0.093200);
	AMP_rn = (AMP_r**2)/(sqrt(1-(AMP_r**2)));
	AMP_rd = (AMP_r)/(sqrt(1-(AMP_r**2)));

	TCP_k = 0.4548463;
	TCP_s = sqrt(1929.243);
	TCP_q = 187.0878515;
	TCP_x = (TCP-TCP_q)*(TCP_k/(TCP_s**2));
	TCP_d = ((TCP_k / TCP_s)**2);
	TCP_r = sqrt(0.021282);
	TCP_rn = (TCP_r**2)/(sqrt(1-(TCP_r**2)));
	TCP_rd = (TCP_r)/(sqrt(1-(TCP_r**2)));

	CMVOD_k = 0.0244357162;
	CMVOD_s = sqrt(1.578067);
	CMVOD_q = 0.4608126456;
	CMVOD_x = (CMVOD-CMVOD_q)*(CMVOD_k/(CMVOD_s**2));
	CMVOD_d = ((CMVOD_k / CMVOD_s)**2);
	CMVOD_r = sqrt(0.071257);
	CMVOD_rn = (CMVOD_r**2)/(sqrt(1-(CMVOD_r**2)));
	CMVOD_rd = (CMVOD_r)/(sqrt(1-(CMVOD_r**2)));

	APPSI_k = 0.29606655;
	APPSI_s = sqrt(1075.323);
	APPSI_q = 69.64017103;
	APPSI_x = (APPSI-APPSI_q)*(APPSI_k/(APPSI_s**2));
	APPSI_d = ((APPSI_k / APPSI_s)**2);
	APPSI_r = sqrt(0.016260);
	APPSI_rn = (APPSI_r**2)/(sqrt(1-(APPSI_r**2)));
	APPSI_rd = (APPSI_r)/(sqrt(1-(APPSI_r**2)));

	FEV_k = -41.305355;
	FEV_s = sqrt(592796);
	FEV_q = 5573.361947;
	FEV_x = (FEV-FEV_q)*(FEV_k/(FEV_s**2));
	FEV_d = ((FEV_k / FEV_s)**2);
	FEV_r = sqrt(0.368526);
	FEV_rn = (FEV_r**2)/(sqrt(1-(FEV_r**2)));
	FEV_rd = (FEV_r)/(sqrt(1-(FEV_r**2)));

	BUP_k = 0.09116738;
	BUP_s = sqrt(26.7675);
	BUP_q = 11.11701117;
	BUP_x = (BUP-BUP_q)*(BUP_k/(BUP_s**2));
	BUP_d = ((BUP_k / BUP_s)**2);
	BUP_r = sqrt(0.059232);
	BUP_rn = (BUP_r**2)/(sqrt(1-(BUP_r**2)));
	BUP_rd = (BUP_r)/(sqrt(1-(BUP_r**2)));

	SBP_k = 0.5381768;
	SBP_s = sqrt(242.725);
	SBP_q = 100.1292243;
	SBP_x = (SBP-SBP_q)*(SBP_k/(SBP_s**2));
	SBP_d = ((SBP_k / SBP_s)**2);
	SBP_r = sqrt(0.194820);
	SBP_rn = (SBP_r**2)/(sqrt(1-(SBP_r**2)));
	SBP_rd = (SBP_r)/(sqrt(1-(SBP_r**2)));

	BANumerator = CEP_x + CRP_x + GHP_x + AMP_x + TCP_x + CMVOD_x + APPSI_x + FEV_x + BUP_x + SBP_x;
	BADenominator = CEP_d + CRP_d  + GHP_d + AMP_d + TCP_d + CMVOD_d + APPSI_d + FEV_d + BUP_d + SBP_d;
	BA = BANumerator / BADenominator;
	RCHARup = CRP_rn + CEP_rn + GHP_rn + AMP_rn + TCP_rn + CMVOD_rn + APPSI_rn + FEV_rn + BUP_rn + SBP_rn;
	RCHARdown = CRP_rd + CEP_rd + GHP_rd + AMP_rd + TCP_rd + CMVOD_rd + APPSI_rd + FEV_rd + BUP_rd + SBP_rd;
	RCHAR = RCHARup / RCHARdown;
	diffn = (BA - Age);
RUN;

PROC PRINT DATA = finalkm;
	var diffn;
	sum diffn;
RUN;

DATA finalkm;
	SET finalkm;
	S2BAp = (diffn-(16130.15/n))**2; *prelim S2BA, val is sum of diffn;
RUN;

PROC SGSCATTER DATA = finalkm;
	PLOT BA * AGE;
RUN;

PROC PRINT DATA = finalkm;
	var S2BAp;
	sum S2BAP;
RUN;

*Calculating S2(BA);
DATA finalkm;
	SET finalkm;
	S2BA = (924635.88/n) - ( ((1-(RCHAR**2))/(RCHAR**2))*(((Agemax-Agemin)**2)/(12 * m))); *var is sum S2BAP;
RUN;

*Final BA calculations (BA adjusted);
DATA finalkm;
	SET finalkm;

	BARNumerator = CEP_x + CRP_x + GHP_x + AMP_x + TCP_x + CMVOD_x + APPSI_x + FEV_x + BUP_x + SBP_x + (m*(Age/S2BA));
	BARDenominator = CEP_d + CRP_d  + GHP_d + AMP_d + TCP_d + CMVOD_d + APPSI_d + FEV_d + BUP_d + SBP_d + (m/S2BA);
	BAR = BARNumerator / BARDenominator;
RUN;

PROC MEANS DATA = finalkm;
	VAR BAR;
	VAR AGE;
	Title "Biological Age Males";
RUN;

*Women;

DATA finalkw;
	SET nhanes.partI;
	WHERE AGE >= 30 and AGE <= 75 and FEMALE = 1;
RUN; 

/*
*FINDING VALUES FOR VARIABLES (Run SLR);
*Macro for j=CA;
%MACRO SLR (j);
	ODS GRAPHICS ON;
	PROC GLM DATA=finalkw ORDER=FORMATTED PLOTS=(DIAGNOSTICS RESIDUALS);
		MODEL &j = AGE / solution clparm;
		WEIGHT samp_wt;
		TITLE &J;
	RUN;
	QUIT;
	ODS GRAPHICS OFF;
%MEND SLR;

%SLR(CRP)
%SLR(CEP)
%SLR(GHP)
%SLR(AMP)
%SLR(TCP)
%SLR(CMVOD)
%SLR(APPSI)
%SLR(FEV)
%SLR(BUP)
%SLR(SBP)*/

DATA finalkw;
	SET finalkw;

	Agemax = 75;
	Agemin = 30;
	n = 4942;
	m = 10; *number of covariates;

	CRP_k = 0.0056035735; *k slope;
	CRP_s = sqrt(0.559497); *s MSE;
	CRP_q = 0.2075161376; *q Intercept;
	CRP_x = (CRP-CRP_q)*(CRP_k/(CRP_s**2)); *numerator;
	CRP_d = ((CRP_k / CRP_s)**2); *denominator;
	CRP_r = sqrt(0.010789); *r from r2;
	CRP_rn = (CRP_r**2)/(sqrt(1-(CRP_r**2)));
	CRP_rd = (CRP_r)/(sqrt(1-(CRP_r**2))); *denominator of r char;

	CEP_k = 0.0031710591;
	CEP_s = sqrt(0.0655090);
	CEP_q = 0.8175042832;
	CEP_x = (CEP-CEP_q)*(CEP_k/(CEP_s**2));
	CEP_d = ((CEP_k / CEP_s)**2);
	CEP_r = sqrt(0.028967);
	CEP_rn = (CEP_r**2)/(sqrt(1-(CEP_r**2)));
	CEP_rd = (CEP_r)/(sqrt(1-(CEP_r**2)));

	GHP_k = 0.023305405;
	GHP_s = sqrt(1.145501);
	GHP_q = 4.275135665;
	GHP_x = (GHP-GHP_q)*(GHP_k/(GHP_s**2));
	GHP_d = ((GHP_k / GHP_s)**2);
	GHP_r = sqrt(0.084372);
	GHP_rn = (GHP_r**2)/(sqrt(1-(GHP_r**2)));
	GHP_rd = (GHP_r)/(sqrt(1-(GHP_r**2)));

	AMP_k = -0.003314656;
	AMP_s = sqrt(0.1297046);
	AMP_q = 4.236146884;
	AMP_x = (AMP-AMP_q)*(AMP_k/(AMP_s**2));
	AMP_d = ((AMP_k / AMP_s)**2);
	AMP_r = sqrt(0.016195);
	AMP_rn = (AMP_r**2)/(sqrt(1-(AMP_r**2)));
	AMP_rd = (AMP_r)/(sqrt(1-(AMP_r**2)));

	TCP_k = 1.4537342;
	TCP_s = sqrt(1878.71);
	TCP_q = 140.6943597;
	TCP_x = (TCP-TCP_q)*(TCP_k/(TCP_s**2));
	TCP_d = ((TCP_k / TCP_s)**2);
	TCP_r = sqrt(0.179394);
	TCP_rn = (TCP_r**2)/(sqrt(1-(TCP_r**2)));
	TCP_rd = (TCP_r)/(sqrt(1-(TCP_r**2)));

	CMVOD_k = 0.0218676205;
	CMVOD_s = sqrt(1.499649);
	CMVOD_q = 0.8962358207;
	CMVOD_x = (CMVOD-CMVOD_q)*(CMVOD_k/(CMVOD_s**2));
	CMVOD_d = ((CMVOD_k / CMVOD_s)**2);
	CMVOD_r = sqrt(0.058353);
	CMVOD_rn = (CMVOD_r**2)/(sqrt(1-(CMVOD_r**2)));
	CMVOD_rd = (CMVOD_r)/(sqrt(1-(CMVOD_r**2)));

	APPSI_k = 0.65327953;
	APPSI_s = sqrt(962.063);
	APPSI_q = 48.70406255;
	APPSI_x = (APPSI-APPSI_q)*(APPSI_k/(APPSI_s**2));
	APPSI_d = ((APPSI_k / APPSI_s)**2);
	APPSI_r = sqrt(0.079367);
	APPSI_rn = (APPSI_r**2)/(sqrt(1-(APPSI_r**2)));
	APPSI_rd = (APPSI_r)/(sqrt(1-(APPSI_r**2)));

	FEV_k = -32.674667;
	FEV_s = sqrt(287892);
	FEV_q = 4179.144341;
	FEV_x = (FEV-FEV_q)*(FEV_k/(FEV_s**2));
	FEV_d = ((FEV_k / FEV_s)**2);
	FEV_r = sqrt(0.418841);
	FEV_rn = (FEV_r**2)/(sqrt(1-(FEV_r**2)));
	FEV_rd = (FEV_r)/(sqrt(1-(FEV_r**2)));

	BUP_k = 0.133564868;
	BUP_s = sqrt(21.5894);
	BUP_q = 7.028488448;
	BUP_x = (BUP-BUP_q)*(BUP_k/(BUP_s**2));
	BUP_d = ((BUP_k / BUP_s)**2);
	BUP_r = sqrt(0.138366);
	BUP_rn = (BUP_r**2)/(sqrt(1-(BUP_r**2)));
	BUP_rd = (BUP_r)/(sqrt(1-(BUP_r**2)));

	SBP_k = 0.80892684;
	SBP_s = sqrt(260.636);
	SBP_q = 82.24622001;
	SBP_x = (SBP-SBP_q)*(SBP_k/(SBP_s**2));
	SBP_d = ((SBP_k / SBP_s)**2);
	SBP_r = sqrt(0.327919);
	SBP_rn = (SBP_r**2)/(sqrt(1-(SBP_r**2)));
	SBP_rd = (SBP_r)/(sqrt(1-(SBP_r**2)));

	BANumerator = CEP_x + CRP_x + GHP_x + AMP_x + TCP_x + CMVOD_x + APPSI_x + FEV_x + BUP_x + SBP_x;
	BADenominator = CEP_d + CRP_d  + GHP_d + AMP_d + TCP_d + CMVOD_d + APPSI_d + FEV_d + BUP_d + SBP_d;
	BA = BANumerator / BADenominator;
	RCHARup = CRP_rn + CEP_rn + GHP_rn + AMP_rn + TCP_rn + CMVOD_rn + APPSI_rn + FEV_rn + BUP_rn + SBP_rn;
	RCHARdown = CRP_rd + CEP_rd + GHP_rd + AMP_rd + TCP_rd + CMVOD_rd + APPSI_rd + FEV_rd + BUP_rd + SBP_rd;
	RCHAR = RCHARup / RCHARdown;
	diffn = (BA - Age)/n;
RUN;
/*
PROC PRINT DATA = finalkw;
	var diffn;
	sum diffn;
RUN;*/

DATA finalkw;
	SET finalkw;
	S2BAp = ((BA-Age)-3.27401)**2; *prelim S2BA, val is sum of diffn;
RUN;
/*
PROC SGSCATTER DATA = finalkw;
	PLOT BA * AGE;
RUN;

PROC PRINT DATA = finalkw;
	var S2BAp;
	sum S2BAP;
RUN;*/

*Calculating S2(BA);
DATA finalkw;
	SET finalkw;
	S2BA = (756409.98/n) - ( ((1-RCHAR)/RCHAR)*(((Agemax-Agemin)**2)/(12 * m))); *var is sum S2BAP;
RUN;

*Final BA calculations (BA adjusted);
DATA finalkw;
	SET finalkw;

	BARNumerator = CEP_x + CRP_x + GHP_x + AMP_x + TCP_x + CMVOD_x + APPSI_x + FEV_x + BUP_x + SBP_x + (m*(Age/S2BA));
	BARDenominator = CEP_d + CRP_d  + GHP_d + AMP_d + TCP_d + CMVOD_d + APPSI_d + FEV_d + BUP_d + SBP_d + (m/S2BA);
	BAR = BARNumerator / BARDenominator;
RUN;

PROC MEANS DATA = finalkw;
	VAR BAR;
	VAR AGE;
	Label "Biological Age Females";
RUN;

DATA finalkw;
	SET finalkw;
	RENAME BAR = BAW;
RUN;

/* MERGE */

PROC SORT DATA = final;
	BY SEQN;
RUN;

PROC SORT DATA = finalkm;
	BY SEQN;
RUN;

PROC SORT DATA = finalkw;
	BY SEQN;
RUN;

DATA nhanes.final;
	MERGE final finalkm finalkw;
	BY SEQN;
	KEEP SEQN BAM BAR BAW AGE SAMP_WT ;
RUN;

DATA nhanes.final;
	SET nhanes.final;
	BAK = BAR;
	If BAK = . THEN BAK = BAW;
	LABEL BAM = 'Biological Age MLR'
		  BAK = 'Biological Age K&D';
RUN;

PROC MEANS DATA = nhanes.final;
	var Age BAM BAK;
	WEIGHT SAMP_WT;
RUN;


*Double checking stuff;
PROC MEANS DATA = finalkm min max mean;
	var S2BA;
	var rchar;
	Title 'KDM Males';
RUN;

PROC MEANS DATA = finalkw min max mean;
	var S2BA;
	var rchar;
	Title 'KDM Females';
	WEIGHT SAMP_WT;
RUN;

PROC SORT DATA = nhanes.final;
	BY AGE;
RUN;

ODS GRAPHICS ON;
PROC GLM DATA=nhanes.final ORDER=FORMATTED PLOTS(MAXPOINT=10000)=(DIAGNOSTICS RESIDUALS);
	MODEL AGE = BAK / solution;
*	WEIGHT samp_wt;
	TITLE "AGE vs BA";
	OUTPUT OUT=check lcl=low ucl=high stdi=stdi stdp=stdp stdr=stdr;
RUN;
QUIT;
ODS GRAPHICS OFF;

DATA check;
	SET check;
	diff = high-low;
RUN;

PROC SGSCATTER DATA = check;
	PLOT diff*age stdi*age stdp*Age;
RUN;

DATA check;
	SET check;
	diff2 = BAK - Age;
RUN;

PROC SGSCATTER DATA = check;
	PLOT diff2*age;
RUN;

/**** MORTALITY DATA ****/
PROC SORT DATA = check;
	BY SEQN;
RUN;

DATA doublef;
	MERGE check (KEEP = MORTSTAT SEQN) nhanes.final (IN=A);
	BY SEQN;
	IF A;
	IF MORTSTAT = . THEN DELETE;
RUN;


*Chronological Age c=.828;
ods graphics on;
PROC LOGISTIC DATA = doublef DESCENDING PLOTS(MAXPOINTS=10000 ONLY) = (ROC EFFECT) ROCOPTIONS(NODETAILS);
	MODEL MORTSTAT = Age / NOFIT;
	ROC 'Chronological Age' Age;
	ROC 'Chance';
	ROCCONTRAST REFERENCE('Chance') /ESTIMATE;
	WEIGHT SAMP_WT;
RUN;
QUIT;
ods graphics off;*/

/* MLR */
ods graphics on;
PROC LOGISTIC DATA = doublef DESCENDING PLOTS(MAXPOINTS=10000 ONLY) = (ROC EFFECT) ROCOPTIONS(NODETAILS);
	MODEL MORTSTAT = BAM Age/ NOFIT;
	ROC 'BA (MLR)' BAM;
	ROC 'Chronological Age' Age;
	ROC 'BA (KDM)' BAK;
	ROCCONTRAST REFERENCE('Chronological Age') /ESTIMATE;
	WEIGHT Samp_wt;
RUN;
QUIT;
ods graphics off;

/* KDM */
ods graphics on;
PROC LOGISTIC DATA = doublef DESCENDING PLOTS(MAXPOINTS=10000 ONLY) = (ROC EFFECT) ROCOPTIONS(NODETAILS);
	MODEL MORTSTAT = BAK Age/ NOFIT;
	ROC 'BA (KDM)' BAK;
	ROC 'Chronological Age' Age;
	ROCCONTRAST REFERENCE('Chronological Age') /ESTIMATE;
	WEIGHT Samp_wt;
RUN;
QUIT;
ods graphics off;

/* CA */
ods graphics on;
PROC LOGISTIC DATA = doublef DESCENDING PLOTS(MAXPOINTS=10000 ONLY) = (ROC EFFECT) ROCOPTIONS(NODETAILS);
	MODEL MORTSTAT = Age/ NOFIT;
	ROC 'CA' Age;
	ROC 'Chronological Age' Age;
	ROCCONTRAST REFERENCE('Chronological Age') /ESTIMATE;
	WEIGHT Samp_wt;
RUN;
QUIT;
ods graphics off;

PROC SGSCATTER DATA = doublef;
	PLOT Samp_wt*Age;
RUN;
