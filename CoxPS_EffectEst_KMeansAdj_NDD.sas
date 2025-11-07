***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring
*
* This code makes effect estimation for NDD (including ASD) using PS-weighted Cox model adjusted for K-means exposure cluster
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*			Choose PS model and select cohort
*           Perform cox proportional hazards regression model
*			Check proportionality assumption	
*			Table creation
***************************************************************************************************************************;


***************************************************************************************************************************;
*Step 1: Read in dataset;            
****************************************************************************************************************************;
* Read in "ps_weight_logistic" dataset, described in "01_data_description.docx";
data ps_data; 
	set dta.ps_weight_logistic;
run;

***************************************************************************************************************************;
*Step2: Choose PS model and select cohort 	
***************************************************************************************************************************;
* Choose a dataset for the finally chosen PS model;  
data cox; 
	set ps_data; * PS model with logistic regression is chosen for the primary outcome analysis, due to best balance obtained; 

    * select cohort AND ps-weighted offspring; 
	if C1_COMP=1 and not missing(weight);
	rename weight=sw; * stabilized weights;

	* Follow up time (in years);
	* EXIT_DT takes into account date of first diagnosis of NDD including ASD; 
	fu_time=(EXIT_DT-INDEX_DT)/365;

    * 'missing' categories are coded as . missing. 
	These missing values won't be included in the cox model;	
    if M_SMOKE_12MLB=9 then M_SMOKE_12MLB=.;
	if M_SMOKE_PREG=9 then M_SMOKE_PREG=.;
 
    label fu_time='Follow-up time (years)'; 
run; 
* All patients should have outlier=0 (not outlier);

* Remove all formats, so that it is easier to define reference groups for logistic regression;
proc datasets lib=work memtype=data noprint;
	modify cox;
	attrib _all_ format=;
run;
quit;

***************************************************************************************************************************;
*Step 3: Perform cox proportional hazards regression model	
***************************************************************************************************************************;
* Inclusion of confounders retained in the PS that were still unbalanced after weighting will be considered for inclusion; 
* Check table 30 to see which confounders are imbalanced. %candidates = imbalanced confounders; 
* The imbalanced confounders should be assciated with both exposure and outcome;
* Including kmeans variable;  
%let candidates=				/* No confounder meets the criteria;*/
/* k-means variable */ F_KMN_3MLB; 

* Macro to select reference groups in logistic regression; 
* Select the lowest category as the reference group;
* Set 26-30 as reference for mother's age, 31-35 for father's age;
* 26-30 years old is the most populated group for mothers. 30-35 years old is most populated for fathers.; 
%macro get_ref(vars=, ds=); 

%if %length(&vars)>0 %then %do;

	%do i=1 %to %sysfunc(countw(&vars.));
    %let var=%upcase(%scan(&vars., &i.));
    %put &var.; 

    %if "&var."="M_AGECATINDEX" %then %do;
		%let var_cat=&var. (ref='3');
	%end; 

    %else %if "&var."="F_AGECATINDEX" %then %do;
		%let var_cat=&var. (ref='4');
	%end;
	%else %do; 
		* get the lowest category; 
		proc sql noprint; 
			select cat("(ref='", min(&var.), "')") 
        	into:low_cat 
        	from &ds.
            where not missing(&var.); 
    	quit;

		* variable name + lowest category; 
		%let var_cat=&var. &low_cat;
    %end;
 
	* accumulating var_cat; 
    %let refs=&refs &var_cat;
		 
	%end; 

%end; 

%mend;

* reference groups; 
%let refs=; 
%get_ref(vars=&candidates., ds=cox);
%put &refs;

*note: the hazard ratio for PATEXPVALP will NOT be calculated and SAS will give a note.: 
NOTE: Hazard ratios that cannot be conveniently calculated or displayed are
      set to missing in the ParameterEstimates table. Use the HAZARDRATIO
      statement to compute the needed hazard ratios.;  

* a related SAS document (https://url.de.m.mimecastprotect.com/s/HCj2CEqYoWtV1ojjvIyHjt77GW0?domain=support.sas.com)
the lack of hazard ratios for effects which are involved in interaction terms occurs 
because their hazard ratios depend on the level of the other variable(s) in the interaction.; 

* covs(aggregate) = correct the standard errors using sandwich estimators. 
Some offspring share the same father and therefore the correction is necessary.; 
ods output ParameterEstimates=param CensoredSummary=summary  HazardRatios=Hazard; 
proc phreg data=cox namelen=32 covs(aggregate); 
    class PATEXPVALP(ref='0') &refs.; * reference = lamotrigine/levetiracetam; 
	model fu_time*O_NDDALL(0)=PATEXPVALP F_KMN_3MLB PATEXPVALP*F_KMN_3MLB &candidates. /rl;   * O_NDDALL = 0 = censored, O_NDDALL = 1 = NOT censored; 
    id patid; 
    weight sw; * PS weighting, use stabilized weights; 
    hazardratio 'Effect of valproate across k-means clusters' PATEXPVALP / at (F_KMN_3MLB=all) diff=ref ;
    output out=residuals RESSCH=res1-res100;  * set 100 residual variable names. There would be less than 100 variables in reality, 
	but 100 is used here to have enough amount of names;
run;  

* Extract variable names and category; 
data param_2; 
	set param; 
	Variable=upcase(Parameter);
	if variable="PATEXPVALP*F_KMN_3MLB" then category=input(ClassVal1, 8.); * group category;
    else category=input(ClassVal0, 8.);  
	proc sort; by variable category; 
run; 

* Extract categories of interaction terms; 
data hazard_2; 
	set hazard;
    * extract kmeans category; 
	if not missing(Description) then do; 
		category=input(substr(Description,length(Description),1),8.);
		ClassVal0=put(category, 1.);
    end; 
	* assign a variable name; 
    variable="PATEXPVALP*F_KMN_3MLB_2"; 
run;  

* Add the two datasets vertically;
data param_3;
    length variable label $100; 
	set param_2
        hazard_2;
run;  

* Create an empty dataset with nvar, _line and labels for table creation;
* If interactions are used in the model, add interaction terms to comm; 
data comm;
	length variable $ 200 level $ 200;
	infile datalines delimiter='+';
	input variable category level nvar _line;
	datalines;
    PATEXPVALP + 1 + Paternal exposure: valproate vs lamotrigine/levetiracetam + 0 + 0
    F_KMN_3MLB +  + K-means exposure cluster: + 0 + 1
	F_KMN_3MLB + 2 + K-means exposure cluster 2 + 0 + 2
	F_KMN_3MLB + 3 + K-means exposure cluster 3 + 0 + 3
	F_KMN_3MLB + 4 + K-means exposure cluster 4 + 0 + 4
	PATEXPVALP*F_KMN_3MLB +  + Paternal exposure * cluster: + 0 + 5
	PATEXPVALP*F_KMN_3MLB + 2 + Valproate * cluster 2 + 0 + 6
	PATEXPVALP*F_KMN_3MLB + 3 + Valproate * cluster 3 + 0 + 7
	PATEXPVALP*F_KMN_3MLB + 4 + Valproate * cluster 4 + 0 + 8
	PATEXPVALP*F_KMN_3MLB_2 + + Effect of valproate across K-means cluster: + 0 + 9
    PATEXPVALP*F_KMN_3MLB_2 + 1 + Valproate vs lamotrigine/levetiracetam in cluster 1 + 0 + 10
	PATEXPVALP*F_KMN_3MLB_2 + 2 + Valproate vs lamotrigine/levetiracetam in cluster 2 + 0 + 11
	PATEXPVALP*F_KMN_3MLB_2 + 3 + Valproate vs lamotrigine/levetiracetam in cluster 3 + 0 + 12
    PATEXPVALP*F_KMN_3MLB_2 + 4 + Valproate vs lamotrigine/levetiracetam in cluster 4 + 0 + 13
	GENDER + + Gender + 1 + 0  
	GENDER + 2 + Female + 1 + 3  
    O_CONCMV_1YR + 1 + Congenital CMV + 2 + 1
    O_CONRUBE_1YR + 1 + Congenital rubella + 3 + 1
    O_FAS_EXIT + 1 + Foetal alcohol syndrome + 4 + 1
    O_FRAGILEX_EXIT + 1  + Fragile X syndrome + 5 + 1
    O_LEJEUNE_EXIT + 1 + Lejeune/cri du chat syndrome + 6 + 1
    O_TUBSIS_EXIT + 1 + Tuberous sclerosis + 7 + 1 
    M_AGECATINDEX +  + Mother's age (categorical) + 8 + 0
    M_AGECATINDEX + 1 + <=20 years + 8 + 2
    M_AGECATINDEX + 2 + 21-25 + 8 + 3
	M_AGECATINDEX + 4 + 31-35 + 8 + 5
 	M_AGECATINDEX + 5 + 36-40 + 8 + 6
	M_AGECATINDEX + 6 + >40 + 8 + 7
	M_AFFECTDIS_INDEX + 1 + Maternal affective disorder + 9 + 1
	M_DIAB_INDEX + 1 + Maternal diabetes + 10 + 1
    M_GESTDIAB_PREG + 1 + Maternal gestational diabetes + 11 + 1 
    M_NEURODIS_INDEX  + 1 + Maternal Neurotic disorder + 12 + 1 
    M_SCHIZO_INDEX  + 1 + Maternal schizophrenia, schizotypal and delusional disorders + 13 + 1
    M_OBES_12MLB  + 1 + Maternal obesity + 14 + 1 
    M_CMV_PREG  + 1 + Maternal CMV + 15 + 1 
    M_RUBELLA_PREG   + 1 + Maternal rubella + 16 + 1  
    M_ALCABUSE_12MLB   + 1 + Maternal alcohol abuse prior to LMP2 + 17 + 1 
    M_ALCABUSE_PREG  + 1 + Maternal alcohol abuse during pregnancy + 18 + 1 
    M_SUBSABUSE_12MLB   + 1 + Maternal substance abuse prior to LMP2 + 19 + 1  
    M_SUBSABUSE_PREG  + 1 + Maternal substance abuse during pregnancy + 20 + 1  
    M_SMOKE_12MLB  +  + Maternal smoking prior to LMP2 + 21 + 0 
	M_SMOKE_12MLB  + 1 + Yes + 21 + 3
    M_SMOKE_PREG  + + Maternal smoking during pregnancy  + 22 + 0
	M_SMOKE_PREG  + 1 + Yes + 22 + 3 
	M_PPICAT_3MLB + + Maternal polypharmacy index prior to LMP2 (categorical) + 23 + 0  
    M_PPICAT_3MLB + 1 + 1-4 + 23 + 3
	M_PPICAT_3MLB + 2 + 5-10 + 23 + 4
	M_PPICAT_3MLB + 3 + >10 + 23 + 5
    M_PPICAT_PREG + + Maternal polypharmacy index during pregnancy (categorical) + 24 + 0
	M_PPICAT_PREG + 1 + 1-4 + 24 + 3
	M_PPICAT_PREG + 2 + 5-10 + 24 + 4
	M_PPICAT_PREG + 3 + >10 + 24 + 5
	M_VALPSYCH_12MLB + 1 + Maternal concomitant medications associated with valproate-indicated psychiatric conditions prior to LMP2 - mothers with at least one prescription + 25 + 1
    M_VALPSYCH_PREG  + 1 + Maternal concomitant medications associated with valproate-indicated psychiatric conditions during pregnancy - mothers with at least one prescription + 26 + 1
    M_NEUROPSYCH_12MLB  + 1 + Maternal concomitant medications associated with neuropsychiatric adverse events prior to LMP2 - mothers with at least one prescription + 27 + 1 
    M_NEUROPSYCH_PREG + 1 + Maternal concomitant medications associated with neuropsychiatric adverse events during pregnancy - mothers with at least one prescription + 28 + 1 
    F_AFFDISEXBM_LMP + 1 + Paternal affective disorder  + 29 + 1
    F_BIPOLAR_LMP + 1 + Paternal bipolar affective disorder + 30 +  1
    F_MANIA_LMP + 1 + Paternal mania + 31 + 1
    F_NEURODIS_LMP + 1 + Paternal neurotic disorder + 32 + 1 
    F_SCHIZO_LMP + 1 + Paternal schizophrenia, schizotypal and delusional disorders + 33 + 1 
    F_SUBSABUSE_12MLB + 1 + Paternal substance abuse + 34 + 1 
    F_PPICAT_3MLB + + Paternal polypharmacy index (categorical) + 35 + 0
	F_PPICAT_3MLB + 1 + 1-4 + 35 + 3
	F_PPICAT_3MLB + 2 + 5-10 + 35 + 4
	F_PPICAT_3MLB + 3 + >10 + 35 + 5
    F_VALPSYCH_12MLB + 1 + Paternal concomitant medications associated with valproate-indicated psychiatric conditions  - fathers with at least one prescription + 36 + 1
    F_NEUROPSYCH_12MLB + 1 + Paternal concomitant medications associated with neuropsychiatric adverse events - fathers with at least one prescription + 37 + 1
    F_AGECATINDEX +  + Father's age (categorical) + 38 + 0
	F_AGECATINDEX + 1 + <=20 years + 38 + 2
    F_AGECATINDEX + 2 + 21-25 + 38 + 3
	F_AGECATINDEX + 3 + 26-30 + 38 + 4
 	F_AGECATINDEX + 5 + 36-40 + 38 + 6
	F_AGECATINDEX + 6 + >40 + 38 + 7
    F_YRCATCONCEPT +  + Year of offspring conception + 39 + 0
	F_YRCATCONCEPT + 3 + 2011-2015 + 39 + 3
	F_YRCATCONCEPT + 4 + 2016-2019 + 39 + 4
   ; 
   proc sort; by variable category; 
run; 

* Add variable labels; 
data param_4; 
	merge comm (in=a)
	      param_3 (in=b);
    by variable category;
	proc sort; by nvar _line;
run;

* Choose variables included in the model;
proc sql; 
	create table param_5 as 
    select *
    from param_4
    where variable in (select variable from param_3); 
quit; 

* BY variable is missing in merging, and it is OK. This is just adding "Total" variable to PARAM dataset; 
data to_report; 
	merge summary 
          param_5;

    * remove non-existent k-means categories;
    if not missing(category) and missing(ClassVal0) then delete;
  
    * hazard ratios for the interaction term;
    if variable="PATEXPVALP*F_KMN_3MLB_2" then do;
	    HRLowerCL=RobustWaldLower;
        HRUpperCL=RobustWaldUpper;
	end; 

	length total_char HR CI pval $20.;

    if not missing(total) then total_char=put(total, best8.);
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

	array cols total_char HR CI pval;
	do over cols; 
		if missing(cols) and not missing(category) then cols='-';
	end; 
	* add this code to keep two decimals;
	if HR not in ('', '-') then HR = cats(HR, '^{}');

	keep total variable level total_char HR CI pval nvar _line;

run;

***************************************************************************************************************************;
*Step 4: Check proportionality assumption	
***************************************************************************************************************************;
* get residual variable names;
proc contents data=residuals out=residual_contents noprint; run;
proc sql; select name into:res_name separated by ' ' from residual_contents where name contains 'res'; quit;  

%put &res_name;

* 1. Test the proportionality assumption with a statistical test;
* For each approach, the proportionality assumption will be tested 
using statistical tests and graphical diagnostics based on the scaled Schoenfeld residuals; 
data event; 
	set residuals;
    if O_NDDALL=1; 
run; 

proc rank data=event out=ranked ties=mean;
	var fu_time;
    ranks timerank; 
run;
 
* Test schoenfeld residuals have no correlation with the rank of the survival time;
* if the p-value>=0.05, it means there is no correlation between the residuals and the survival time. 
The proportionality assumption holds;  
proc corr data=ranked;
	var &res_name; * residuals; 
	with timerank;
run; 

* 2. Residual plots;
* Schoenfeld residuals are plotted against time and if the plot shows 
a non-random pattern over time, there is evidence of violation of the proportional hazard assumption.; 
%macro res_check(vars=&res_name.); 

	%do i = 1 %to %sysfunc(countw(&vars.));
	 %let var =  %upcase(%scan(&vars., &i.));

        %let track_date=%sysfunc(today(),date9.);
	 	ods listing image_dpi=300;
		ods listing gpath="&output.\02 Figures\";
		ods graphics / imagename="&var._plot_table36 &track_date." imagefmt=jpg;
	    title "Schoenfeld residual plot";
		proc sgplot data=residuals;
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
ods graphics / imagename="loglog_table36 &track_date." imagefmt=jpg;
proc lifetest data=to_cox plot=loglogs notable;
  time fu_time*O_NDDALL(0);
  strata PATEXPVALP;
run; 
ods listing close; 

***************************************************************************************************************************;
*Step 5: Table creation	
***************************************************************************************************************************;
** Set title **;
options noquotelenmax;
%let title="Table 36. Effect estimation for NDD using PS-weighted Cox model adjusted for K-means exposure cluster; primary outcome";

** Set footnotes **;
title &title.;

%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table36 &track_date..xlsx" options(SHEET_NAME="Table36"  embedded_titles='yes' embedded_footnotes='yes');
ods escapechar="^";
proc report data=to_report  split = '*';
	column  
             ("Variable" "" level) 
             ("Total N" "" total_char)
             ("Model estimates" ("HR" HR) ("95% CI" CI) ("P-value" pval));
    
	define level    / display ""
                      style(header)={vjust=middle just=center width= 8 cm}
					  style(column)={just=l cellwidth = 8 cm};
	define total_char		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define HR		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm tagattr="format:0.00"};
	define CI		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define pval		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm tagattr="format:0.0000"};

	compute level;
		if find(level, "Paternal exposure *", 'i')>=1 or find(level, "K-means exposure cluster:", 'i')>=1 or
		   find(level, "Effect of", 'i')>=1 or 
           find(level,"Gender",'i') >= 1 or find(level,"age",'i') >= 1 or find(level,"Smoking",'i') >= 1 or
           find(level,"polypharmacy",'i') >= 1 or find(level," offspring conception",'i') >= 1 then do;
            call define (_col_,"style","style={fontweight=bold background=lightgrey}");
		end; 
    endcomp;

run;
ods excel close;



**********************************************************END*****************************************************************;

