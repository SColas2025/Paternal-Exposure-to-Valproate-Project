***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring
*
* This code makes effect estimation for NDD (including ASD) using crude Cox model
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*           Perform cox proportional hazards regression model
*           Calculate dfbetas and exclude influential data
*			Re-estimation of the Cox model
*			Check proportionality assumption	
*			Table creation
***************************************************************************************************************************;

***************************************************************************************************************************;
*Step 1: Read in dataset;            
****************************************************************************************************************************;
* Read in "mydata" dataset, described in "01_data_description.docx";
data in_data; 
	set dta.Mydata;
run;

***************************************************************************************************************************;
*Step2: Perform cox proportional hazards regression model	
***************************************************************************************************************************;
 data cox; 
	set in_data;

    * Select cohort; 
	if C1_COMP=1;

	* Follow up time (in years);
	* EXIT_DT takes into account date of first diagnosis of NDD including ASD; 
	fu_time=(EXIT_DT-INDEX_DT)/365; 

	label fu_time='Follow-up time (years)'; 
run; 

* Remove all formats, so that it is easier to define reference groups for cox regression;
proc datasets lib=work memtype=data noprint;
	modify cox;
	attrib _all_ format=;
run;
quit;

* Cox regression;
* covs(aggregate) = correct the standard errors using sandwich estimators. 
Some offspring share the same father and therefore the correction is necessary.;  
proc phreg data=cox covs(aggregate); 
    class PATEXPVALP(ref='0'); * reference = lamotrigine/levetiracetam; 
	model fu_time*O_NDDALL(0)=PATEXPVALP/rl;   * O_NDDALL = 0 = censored, O_NDDALL = 1 = NOT censored; 
	id patid; 
	output out=diagnosis dfbeta=dfbet;
run;  

***************************************************************************************************************************;
*Step3: Calculate dfbetas and exclude influential data
***************************************************************************************************************************;
* Outliers are defined as dfbetas outside the range +/-2/sqrt(n), where n is the number of observations; 
proc sql;
	select 2/sqrt(count(*)) into:threshold from cox;  
quit; 

* 1 = outlier, outside the range, 0=No outlier, inside the range;  
data cox_; 
	set diagnosis; 
	if abs(dfbet)>&threshold then outlier=1; * outside the range +/-2/sqrt(n); 
	else if abs(dfbet)<=&threshold then outlier=0; * inside the range; 
run; 

***************************************************************************************************************************;
*Step4: Re-estimation of the Cox model
***************************************************************************************************************************;
* Cox regression;
* covs(aggregate) = correct the standard errors using sandwich estimators. 
Some offspring share the same father and therefore the correction is necessary.;  
ods output ParameterEstimates=Est CensoredSummary=summary; 
proc phreg data=cox_ covs(aggregate); 
    class PATEXPVALP(ref='0'); * reference = lamotrigine/levetiracetam; 
	model fu_time*O_NDDALL(0)=PATEXPVALP/rl;   * O_NDDALL = 0 = censored, O_NDDALL = 1 = NOT censored; 
	id patid;  
	output out=residuals RESSCH=res_PATEXPVALP;
	where outlier=0;  *exclude outliers; 
run;  

* Extract results for table creation;
data to_report_; 
	set est; 
    keep HazardRatio HRLowerCL HRUpperCL ProbChiSq HR CI pval var;

	length HR CI pval $20.;

	*if not missing(total) then N_total_obs=put(total, 8.);
	if not missing(HazardRatio) then HR=put(HazardRatio, 10.2);

	* confidence interval; 
	if not missing(HRLowerCL) then do; 
		if not missing(HRUpperCL) then do; 
		CI='(' || trim(left(put(HRLowerCL,10.2))) || ' , ' || trim(left(put(HRUpperCL,10.2))) || ')';
		end;
		else CI='(' || trim(left(put(HRLowerCL,10.2))) || ' , ' || '-' || ')';
    end;
    else do;
		if not missing(HRUpperCL) then do;
			CI='(' || '-' || ' , ' || trim(left(put(HRUpperCL,10.2))) || ')'; 
        end;  
    end;

    if not missing(ProbChiSq) then pval=put(ProbChiSq, pvalue9.4); 

	array cols N_total_obs HR CI pval;
	do over cols; 
		if missing(cols) then cols='-';
	end;

	var="Paternal exposure: valproate vs lamotrigine/levetiracetam"; 

run;

proc sql; 
     select count(*) into:total_obs from cox_;
     select count(*) into:normal_obs from cox_ where outlier=0; 
quit; 

data to_report;
	set to_report_;
	N_total_obs=put(&total_obs.,8.);
	N_normal_obs=put(&normal_obs., 8.);
	pct_normal=cats(put((&normal_obs./&total_obs.)*100, 10.2), "^{}");
run; 

***************************************************************************************************************************;
*Step 5: Check proportionality assumption	
***************************************************************************************************************************;
* 1. Test the proportionality assumption with a statistical test; 
data event; 
	set residuals;
    if O_NDDALL=1; 
run; 

proc rank data=event out=ranked ties=mean;
	var fu_time;
    ranks timerank; 
run;
 
* Test schoenfeld residuals have no correlation with the rank of the survival time;
* If the p-value>=0.05, it means there is no correlation between the residuals and the survival time. 
The proportionality assumption holds;  
proc corr data=ranked;
	var res_PATEXPVALP;
	with timerank;
run; 

* 2. Residual plots; 
data residuals_2; 	
	set residuals; 
	label res_PATEXPVALP='Schoenfeld residual for valproate exposure';
run;

* Schoenfeld residuals are plotted against time and if the plot shows 
a non-random pattern over time, there is evidence of violation of the proportional hazard assumption.; 
%macro res_check(vars=res_PATEXPVALP); 

	%do i = 1 %to %sysfunc(countw(&vars.));
	 %let var =  %upcase(%scan(&vars., &i.));

	 
        %let track_date=%sysfunc(today(),date9.);
	 	ods listing image_dpi=300;
		ods listing gpath="&output.\02 Figures\";
		ods graphics / imagename="&var._plot_table34 &track_date." imagefmt=jpg;
	    title "Schoenfeld residual plot";
		proc sgplot data=residuals_2;
			scatter y=&var. x=fu_time/ name="Schoenfeld_CHECK";
			loess y=&var. x=fu_time/SMOOTH=.8 DEGREE=1 INTERPOLATION= CUBIC ;
			keylegend "Schoenfeld_CHECK"/ title=' ' Position=bottom;
		run;
		ods listing close; 		
	%end;
%mend; 
%res_check;

* 3. log-log plots.;	
* Log-log plots are another visual diagnostic tool: in this case, -log(-log(Survival(t)) is 
plotted against log(t) for each exposure group. this plot should give parallel lines.; 
%let track_date=%sysfunc(today(),date9.); 
ods listing image_dpi=300;
ods listing gpath="&output.\02 Figures\";
ods graphics / imagename="loglog_table34 &track_date." imagefmt=jpg;
proc lifetest data=cox_ (where=(outlier=0)) plot=loglogs; * (where=(outlier=0)) exclude outliers;  
  time fu_time*O_NDDALL(0);
  strata PATEXPVALP;
run; 
ods listing close;

***************************************************************************************************************************;
*Step 6: Table creation	
***************************************************************************************************************************;
** Set title **;
options noquotelenmax;
%let title="Table 34. Effect estimation for NDD using crude Cox model; primary outcome";

** Set footnotes **;
%let fn1="a) Influential subjects will be identified using the dfbetas for the main exposure coefficient."; 

title &title.;
footnote1 j=l &fn1.;

%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table34 &track_date..xlsx" options(SHEET_NAME="Table34"  embedded_titles='yes' embedded_footnotes='yes');
ods escapechar="^";
proc report data=to_report  split = '^';
	column  
             ("Variable" "" var) 
			 ("Total N" "N" N_total_obs)
			 ("Number of subjects ^ included in the model ^ (after excluding influential subjects) a" ("N" N_normal_obs) ("%" pct_normal) )
             ("Model estimates" ("HR" HR) ("95% CI" CI) ("P-value" pval));
    
	define var    / display ""
                      style(header)={vjust=middle just=center width= 8 cm}
					  style(column)={just=l cellwidth = 8 cm};
    define  N_total_obs	/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
    define 	N_normal_obs	/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define 	pct_normal	/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define HR		/ display "" 
					  style(column)={vjust=middle just=center width= 2 cm}
                      style(column)={just=c cellwidth = 2 cm tagattr="format:0.00"};
	define CI		/ display "" 
					  style(column)={vjust=middle just=center width= 2 cm}
                      style(column)={just=c cellwidth = 2 cm};
	define pval		/ display "" 
					  style(column)={vjust=middle just=center width= 2 cm}
                      style(column)={just=c cellwidth = 2 cm tagattr="format:0.0000"};
run;
ods excel close;


**********************************************************END*****************************************************************;

