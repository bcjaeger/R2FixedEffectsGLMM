
/*  ------PROGRAM DESCRIPTION:  Inputs, Processes, Outputs-------- 	*/
/*																	*/
/*    Inputs:    Simulated longitudinal data						*/
/*    Processes: Computes Rsq_beta for all outcomes					*/
/*    Outputs:   A dataset containing the Rsq estimates				*/	
/*																	*/
/*  ---------------------PROGRAMMER HISTORY-----------------------	*/
/*																	*/
/*    Programmer(s)  Date(s)    Brief Description of Modifications.	*/
/*    Byron Jaeger   28JUL2015  Initial Program						*/
/*    Byron Jaeger   05NOV2015  Updated for R2 paper				*/
/*																	*/
/*  -------------------------------------------------------------- 	*/


ods listing; options ls=75 ps=50 PAGENO=1;
%LET PATH1=C:\Users\Byron\Desktop\School\Dissertation\SAS Code\;		* Path name									;
libname mylib "&PATH1\GLMM R2"; 										* Set libname								;
%inc "&PATH1\GLMM R2\Glimmix_R2_V3.sas"; 								* The Glimmix Macro							;
/* Read in simulated longitudinal data */

%MACRO simdatR2(iter); 
%do i=1 %to &iter;
	PROC IMPORT OUT= WORK.LongDat&i 
            	DATAFILE= "C:\Users\Byron\Desktop\School\Dissertation\R Code\Data\LongDat&i" 
            	DBMS=CSV REPLACE;
     	GETNAMES=YES;
     	DATAROW=2; 
	RUN;
%end; 
%MEND simdatR2;

%LET Pred1 = X1i obstime;								* Predictors for Reduced Model				;
%LET Pred2 = X1i obstime X1i*obstime;					* Predictors for Full Model 				;
%LET Cov_1 = int;										* Covariance structure 1					;
%LET Cov_2 = int obstime;								* Covariance structure 2					;

/* 
 *  create a macro variable for the data set names. 				
 *  the variable memname contains the list of data sets.			
 *  make sure that you only get the distinct data set names 		
 *  then write that list to the macro variable dsnames
 */

%LET nsims = 1000;
%simdatR2(iter=1000);

proc sql noprint;
	select distinct(memname) into :dsnames separated by ' '
	from dictionary.columns
	where libname="WORK";
quit;
%put &dsnames;  /* this shows you the list of the data set names in the log */

%macro R2_SIM(Out_Loc =);

%let q=1;								/* give the loop a starting point 								   */
%do %until (%scan(&dsnames,&q," ")=);	/* run this loop until after the last item in the list is reached  */
%let currentds=%scan(&dsnames,&q," ");	/* assign the current dataset name as the macro variable currentds */

*Poisson Outcomes;
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Pois_Pred1_Cov1_&q,_ID=id,_OUTCOME=yij_pois,_PREDICTORS=&Pred1,_CLASS=X1i,_RANDOM=&Cov_1,_ERROR=Poisson,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Pois_Pred2_Cov1_&q,_ID=id,_OUTCOME=yij_pois,_PREDICTORS=&Pred2,_CLASS=X1i,_RANDOM=&Cov_1,_ERROR=Poisson,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Pois_Pred1_Cov2_&q,_ID=id,_OUTCOME=yij_pois,_PREDICTORS=&Pred1,_CLASS=X1i,_RANDOM=&Cov_2,_ERROR=Poisson,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Pois_Pred2_Cov2_&q,_ID=id,_OUTCOME=yij_pois,_PREDICTORS=&Pred2,_CLASS=X1i,_RANDOM=&Cov_2,_ERROR=Poisson,SINTAX_ALL=);

*Binomial Outcomes;
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Bin_Pred1_Cov1_&q,_ID=id,_OUTCOME=yij_Bin,_PREDICTORS=&Pred1,_CLASS=X1i,_RANDOM=&Cov_1,_ERROR=Binomial,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Bin_Pred2_Cov1_&q,_ID=id,_OUTCOME=yij_Bin,_PREDICTORS=&Pred2,_CLASS=X1i,_RANDOM=&Cov_1,_ERROR=Binomial,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Bin_Pred1_Cov2_&q,_ID=id,_OUTCOME=yij_Bin,_PREDICTORS=&Pred1,_CLASS=X1i,_RANDOM=&Cov_2,_ERROR=Binomial,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Bin_Pred2_Cov2_&q,_ID=id,_OUTCOME=yij_Bin,_PREDICTORS=&Pred2,_CLASS=X1i,_RANDOM=&Cov_2,_ERROR=Binomial,SINTAX_ALL=);

*Normal Outcomes;
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Norm_Pred1_Cov1_&q,_ID=id,_OUTCOME=yij,_PREDICTORS=&Pred1,_CLASS=X1i,_RANDOM=&Cov_1,_ERROR=Normal,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Norm_Pred2_Cov1_&q,_ID=id,_OUTCOME=yij,_PREDICTORS=&Pred2,_CLASS=X1i,_RANDOM=&Cov_1,_ERROR=Normal,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Norm_Pred1_Cov2_&q,_ID=id,_OUTCOME=yij,_PREDICTORS=&Pred1,_CLASS=X1i,_RANDOM=&Cov_2,_ERROR=Normal,SINTAX_ALL=);
%Glimmix_R2 (_DATA_IN=&currentds,_OUTNAME=Norm_Pred2_Cov2_&q,_ID=id,_OUTCOME=yij,_PREDICTORS=&Pred2,_CLASS=X1i,_RANDOM=&Cov_2,_ERROR=Normal,SINTAX_ALL=);

%let q=%eval(&q+1);  /* increase q by 1 */
%end;

Data R2_Pois_Pred1_Cov1; set R2_Pois_Pred1_Cov1_1 - R2_Pois_Pred1_Cov1_&nsims; Outcome = 'Poisson'; Pred = 'Reduced Model'; Cov = 'Covariance 1'; run;
Data R2_Pois_Pred2_Cov1; set R2_Pois_Pred2_Cov1_1 - R2_Pois_Pred2_Cov1_&nsims; Outcome = 'Poisson'; Pred = 'Full Model'; Cov = 'Covariance 1'; run;
Data R2_Pois_Pred1_Cov2; set R2_Pois_Pred1_Cov2_1 - R2_Pois_Pred1_Cov2_&nsims; Outcome = 'Poisson'; Pred = 'Reduced Model'; Cov = 'Covariance 2'; run;
Data R2_Pois_Pred2_Cov2; set R2_Pois_Pred2_Cov2_1 - R2_Pois_Pred2_Cov2_&nsims; Outcome = 'Poisson'; Pred = 'Full Model'; Cov = 'Covariance 2'; run;

Data R2_Bin_Pred1_Cov1; set R2_Bin_Pred1_Cov1_1 - R2_Bin_Pred1_Cov1_&nsims; Outcome = 'Binomial'; Pred = 'Reduced Model'; Cov = 'Covariance 1'; run;
Data R2_Bin_Pred2_Cov1; set R2_Bin_Pred2_Cov1_1 - R2_Bin_Pred2_Cov1_&nsims; Outcome = 'Binomial'; Pred = 'Full Model'; Cov = 'Covariance 1'; run;
Data R2_Bin_Pred1_Cov2; set R2_Bin_Pred1_Cov2_1 - R2_Bin_Pred1_Cov2_&nsims; Outcome = 'Binomial'; Pred = 'Reduced Model'; Cov = 'Covariance 2'; run;
Data R2_Bin_Pred2_Cov2; set R2_Bin_Pred2_Cov2_1 - R2_Bin_Pred2_Cov2_&nsims; Outcome = 'Binomial'; Pred = 'Full Model'; Cov = 'Covariance 2'; run;

Data R2_Norm_Pred1_Cov1; set R2_Norm_Pred1_Cov1_1 - R2_Norm_Pred1_Cov1_&nsims; Outcome = 'Normal'; Pred = 'Reduced Model'; Cov = 'Covariance 1'; run;
Data R2_Norm_Pred2_Cov1; set R2_Norm_Pred2_Cov1_1 - R2_Norm_Pred2_Cov1_&nsims; Outcome = 'Normal'; Pred = 'Full Model'; Cov = 'Covariance 1'; run;
Data R2_Norm_Pred1_Cov2; set R2_Norm_Pred1_Cov2_1 - R2_Norm_Pred1_Cov2_&nsims; Outcome = 'Normal'; Pred = 'Reduced Model'; Cov = 'Covariance 2'; run;
Data R2_Norm_Pred2_Cov2; set R2_Norm_Pred2_Cov2_1 - R2_Norm_Pred2_Cov2_&nsims; Outcome = 'Normal'; Pred = 'Full Model'; Cov = 'Covariance 2'; run;

Data R2_Pois; set R2_Pois_Pred1_Cov1 R2_Pois_Pred2_Cov1 R2_Pois_Pred1_Cov2 R2_Pois_Pred2_Cov2; run;
Data R2_Bin; set R2_Bin_Pred1_Cov1 R2_Bin_Pred2_Cov1 R2_Bin_Pred1_Cov2 R2_Bin_Pred2_Cov2; run;
Data R2_Norm; set R2_Norm_Pred1_Cov1 R2_Norm_Pred2_Cov1 R2_Norm_Pred1_Cov2 R2_Norm_Pred2_Cov2; run;


Data R2_out; Set R2_Pois R2_Bin R2_Norm; drop label ProbF SumSq Effect; run;

Data Dat; set R2_out; where NumDF > 1; run;

PROC EXPORT DATA= Dat 
            OUTFILE= &Out_Loc. 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

%mend;
ods results=off; %R2_SIM(Out_Loc = "C:\Users\Byron\Desktop\School\Dissertation\R Code\Data\R2_SIM_13OCT2015.csv");


