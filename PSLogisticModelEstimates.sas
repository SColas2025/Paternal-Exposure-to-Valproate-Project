***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring
*
* This code prints and summarises variable estimates from logistic regression propensity score model
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*           Perform PS model with logistic regression and get odds ratios and p-value	
*           Extract results from PS model and create table
***************************************************************************************************************************;

***************************************************************************************************************************;
*Step 1: Read in dataset;            
****************************************************************************************************************************;
* Read in "PS_weight_logistic" dataset, described in "01_data_description.docx";
data PS_weight_logistic; 
	set dta.PS_weight_logistic;
run;


***************************************************************************************************************************;
*Step2: Perform PS model with logistic regression and get odds ratios and p-value	
***************************************************************************************************************************;
* Check the log file of PS_weight_logistic_v1.0.sas;
* See which covariates are selected as &candidates_revised_2. and add those selected covariates below; 
%let candidates_revised_2=
GENDER  M_AFFECTDIS_INDEX   M_GESTDIAB_PREG M_NEURODIS_INDEX   M_OBES_12MLB           
M_SMOKE_12MLB M_SMOKE_PREG M_VALPSYCH_12MLB M_VALPSYCH_PREG M_NEUROPSYCH_PREG F_AFFDISEXBM_LMP 
F_BIPOLAR_LMP F_MANIA_LMP F_NEURODIS_LMP F_SCHIZO_LMP F_SUBSABUSE_12MLB F_PPICAT_3MLB 
F_NEUROPSYCH_12MLB F_YRCATCONCEPT;


* Macro to select reference groups in logistic regression; 
* Select the lowest category as the reference group;
* Set 26-30 as reference for mother's age, 31-35 for father's age;
* 26-30 years old is the most populated group for mothers. 30-35 years old is most populated for fathers.; 
%macro get_ref(vars=, ds=); 

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

%mend;

* outlier=0 should go into the logistic model; 
proc freq data=PS_weight_logistic; 
	table outlier/list missing; 
run;

data non_outlier;
	set PS_weight_logistic;
	if outlier=0;
run;

* new reference groups; 
%let refs=; 
%get_ref(vars=&candidates_revised_2., ds=non_outlier); 

* Check reference groups for logistic regression; 
%put &refs; 

* PS model without outliers: logistic regression; 
* Valproate = event (1), lamotrigine/levetiracetam = reference (0);  
ods output OddsRatios=Odds ParameterEstimates=Param; 
proc logistic data=non_outlier namelen=32;
	class &refs. /param=ref;  
	model PATEXPVALP(event='1') = &candidates_revised_2.; 
run;

***************************************************************************************************************************;
*Step 3: Extract results from PS model and create table	
***************************************************************************************************************************;
* Get variable name, group category, reference category;
data odds_2;
	set odds;
 
	* variable name; 
	if not missing(Effect) then do; 
		Variable=upcase(strip(scan(Effect,1,' ')));

		* category and reference; 
		do i=1 to countw(Effect, " ");
			if scan(Effect, i, " ")="vs" then do; 
			category=input(scan(Effect, i-1, " "), 8.);
			reference=input(scan(Effect, i+1, " "), 8.);
			end; 
		end;
	end; 
run;  

* The following notes are expected and can be ignored.; 
*NOTE: Numeric values have been converted to character
      values at the places given by: (Line):(Column).
      7591:21
*NOTE: Variable ClassVal1 is uninitialized.; 
data param_2; 
	set param; 
	Variable=upcase(Variable);
	category=input(ClassVal0, 8.); * group category; 
	category2=input(ClassVal1, 8.); * if there are interaction terms, this variable gives second group categories; 
run; 

* Merge ORs and p-values;
proc sql; 
	create table odds_pval as
    select upcase(a.Variable) as variable, 
	       a.category,
           a.category2, 
		   b.reference,
		   b.OddsRatioEst, 
		   b.LowerCL, 
           b.UpperCL, 
           a.ProbChiSq
    from param_2 a
    full join odds_2 b
	on a.variable=b.variable and a.category=b.category
    order by a.variable, a.category; 
quit;  

proc sort data=odds_pval; by variable category; run;  

* Create a empty dataset with nvar, _line and labels for table creation;
* If interactions are used in the model, add interaction terms to comm; 
data comm;
	length variable $ 200 level $ 200;
	infile datalines delimiter='+';
	input variable category level nvar _line;
	datalines;
	+ + Offspring risk factors/confounders + 0.5 + 0 
	GENDER + + Gender ^{super b} + 1 + 0  
	GENDER + 1 + Male  + 1 + 2 
	GENDER + 2 + Female + 1 + 3  
    O_CONCMV_1YR + 1 + Congenital CMV ^{super c} + 2 + 1
    O_CONRUBE_1YR + 1 + Congenital rubella ^{super c} + 3 + 1
    O_FAS_EXIT + 1 + Foetal alcohol syndrome ^{super c} + 4 + 1
    O_FRAGILEX_EXIT + 1  + Fragile X syndrome ^{super c} + 5 + 1
    O_LEJEUNE_EXIT + 1 + Lejeune/cri du chat syndrome ^{super c} + 6 + 1
    O_TUBSIS_EXIT + 1 + Tuberous sclerosis ^{super c} + 7 + 1
    + + Maternal risk factors/confounders + 7.5 + 0 
    M_AGECATINDEX +  + Mother's age ^{super b} (categorical) + 8 + 0
    M_AGECATINDEX + 1 + <=20 years + 8 + 2
    M_AGECATINDEX + 2 + 21-25 + 8 + 3
	M_AGECATINDEX + 3 + 26-30 + 8 + 4
	M_AGECATINDEX + 4 + 31-35 + 8 + 5
 	M_AGECATINDEX + 5 + 36-40 + 8 + 6
	M_AGECATINDEX + 6 + >40 + 8 + 7
	M_AFFECTDIS_INDEX + 1 + Affective disorder ^{super d} + 9 + 1
	M_DIAB_INDEX + 1 + Diabetes ^{super d} + 10 + 1
    M_GESTDIAB_PREG + 1 + Gestational diabetes ^{super e} + 11 + 1 
    M_NEURODIS_INDEX  + 1 + Neurotic disorder ^{super d} + 12 + 1 
    M_SCHIZO_INDEX  + 1 + Schizophrenia, schizotypal and delusional disorders ^{super d} + 13 + 1
    M_OBES_12MLB  + 1 + Obesity ^{super f} + 14 + 1 
    M_CMV_PREG  + 1 + CMV ^{super f} + 15 + 1 
    M_RUBELLA_PREG   + 1 + Rubella ^{super f} + 16 + 1  
    M_ALCABUSE_12MLB   + 1 + Alcohol abuse prior to LMP2 ^{super f} + 17 + 1 
    M_ALCABUSE_PREG  + 1 + Alcohol abuse during pregnancy ^{super e} + 18 + 1 
    M_SUBSABUSE_12MLB   + 1 + Substance abuse prior to LMP2 ^{super f} + 19 + 1  
    M_SUBSABUSE_PREG  + 1 + Substance abuse during pregnancy ^{super e} + 20 + 1  
    M_SMOKE_12MLB  +  + Smoking prior to LMP2 ^{super f} + 21 + 0 
	M_SMOKE_12MLB  + 0 + No + 21 + 2 
	M_SMOKE_12MLB  + 1 + Yes + 21 + 3
    M_SMOKE_PREG  + + Smoking during pregnancy ^{super e}  + 22 + 0
    M_SMOKE_PREG  + 0 + No + 22 + 2 
	M_SMOKE_PREG  + 1 + Yes + 22 + 3 
	M_PPICAT_3MLB + + Maternal polypharmacy index prior to LMP2 ^{super h} (categorical) + 23 + 0  
    M_PPICAT_3MLB + 0 + 0 + 23 + 2 
    M_PPICAT_3MLB + 1 + 1-4 + 23 + 3
	M_PPICAT_3MLB + 2 + 5-10 + 23 + 4
	M_PPICAT_3MLB + 3 + >10 + 23 + 5
    M_PPICAT_PREG + + Maternal polypharmacy index during pregnancy ^{super e} (categorical) + 24 + 0
    M_PPICAT_PREG + 0 + 0 + 24 + 2
	M_PPICAT_PREG + 1 + 1-4 + 24 + 3
	M_PPICAT_PREG + 2 + 5-10 + 24 + 4
	M_PPICAT_PREG + 3 + >10 + 24 + 5
	M_VALPSYCH_12MLB + 1 + Concomitant medications associated with valproate-indicated psychiatric conditions prior to LMP2 ^{super f} - mothers with at least one prescription + 25 + 1
    M_VALPSYCH_PREG  + 1 + Concomitant medications associated with valproate-indicated psychiatric conditions during pregnancy ^{super e} - mothers with at least one prescription + 26 + 1
    M_NEUROPSYCH_12MLB  + 1 + Concomitant medications associated with neuropsychiatric adverse events prior to LMP2 ^{super f} - mothers with at least one prescription + 27 + 1 
    M_NEUROPSYCH_PREG + 1 + Concomitant medications associated with neuropsychiatric adverse events during pregnancy ^{super e} - mothers with at least one prescription + 28 + 1
	+ + Paternal risk factors/confounders + 28.5 + 0 
    F_AFFDISEXBM_LMP + 1 + Affective disorder ^{super d,g} + 29 + 1
    F_BIPOLAR_LMP + 1 + Bipolar affective disorder ^{super d} + 30 +  1
    F_MANIA_LMP + 1 + Mania ^{super d} + 31 + 1
    F_NEURODIS_LMP + 1 + Neurotic disorder ^{super d} + 32 + 1 
    F_SCHIZO_LMP + 1 + Schizophrenia, schizotypal and delusional disorders ^{super d} + 33 + 1 
    F_SUBSABUSE_12MLB + 1 + Substance abuse ^{super f} + 34 + 1 
    F_PPICAT_3MLB + + Paternal polypharmacy index ^{super h} (categorical) + 35 + 0
    F_PPICAT_3MLB + 0 + 0 + 35 + 2
	F_PPICAT_3MLB + 1 + 1-4 + 35 + 3
	F_PPICAT_3MLB + 2 + 5-10 + 35 + 4
	F_PPICAT_3MLB + 3 + >10 + 35 + 5
    F_VALPSYCH_12MLB + 1 + Concomitant medications associated with valproate-indicated psychiatric conditions ^{super f} - fathers with at least one prescription + 36 + 1
    F_NEUROPSYCH_12MLB + 1 + Concomitant medications associated with neuropsychiatric adverse events ^{super f} - fathers with at least one prescription + 37 + 1
    F_AGECATINDEX +  + Father's age ^{super b} (categorical) + 38 + 0
    F_AGECATINDEX + 1 + <=20 years + 38 + 2
    F_AGECATINDEX + 2 + 21-25 + 38 + 3
	F_AGECATINDEX + 3 + 26-30 + 38 + 4
	F_AGECATINDEX + 4 + 31-35 + 38 + 5
 	F_AGECATINDEX + 5 + 36-40 + 38 + 6
	F_AGECATINDEX + 6 + >40 + 38 + 7
    F_YRCATCONCEPT +  + Year of offspring conception ^{super i,j} + 39 + 0
    F_YRCATCONCEPT + 2 + 2006-2010 + 39 + 2
	F_YRCATCONCEPT + 3 + 2011-2015 + 39 + 3
	F_YRCATCONCEPT + 4 + 2016-2019 + 39 + 4
   ; 
   proc sort; by variable category; 
run; 

data odds_pval_2; 
	merge comm (in=a)
	      odds_pval (in=b);
    by variable category;
	length odds CI pval $20;

    if not missing(OddsRatioEst) and UpperCL<999 then do;  * the limit to the upper limit: 999, upper limit more than 999 is not reasonable;
		odds=strip(put(OddsRatioEst, 10.2));
        CI=trim(left(put(LowerCL,10.2))) || ' - ' || trim(left(put(UpperCL,10.2)));
        pval=strip(put(ProbChiSq, pvalue9.4));
	end; 

    drop OddsRatioEst ProbChiSq LowerCL UpperCL; 
	array vars odds CI pval;

	* put reference and '-'; 
	do over vars;
	    * reference for mother/father's age; 
        if (nvar=8 and _line=4) or (nvar=38 and _line=5) then odds='Reference';
		else if nvar ne 8 and nvar ne 38 and _line=2 then odds='Reference';
    	if _line>=1 and missing(vars) then vars='-';
	end;

	* delete intercept;
    if variable='INTERCEPT' then delete;  
	proc sort; by nvar _line;
run;

* Choose variables included in PS weighting;
proc sql; 
	create table to_report_ as 
    select *
    from odds_pval_2
    where variable in (select variable from odds_pval) or nvar in (0.5, 7.5, 28.5); 
quit; 

data to_report;
	set to_report_;

	* add the following code snippet for printing 2 decimals;
	if odds not in ('-', 'Reference', '') then odds = cats(odds, '^{}');
run;

* Create table;
options noquotelenmax;
** Set title **;
%let title="Table 26. Variable estimates from logistic regression propensity score model; primary outcome";

** Set footnotes **;
%let fn1="a) Candidate covariates will be considered to enter the PS model if associated with the study outcome based on univariate analyses. Additionally, two-way interactions will be included in the PS model if identified as clinically meaningful.";
%let fn2="b) at index (childbirth)";
%let fn3="c) between index and exit date";
%let fn4="d) all available data prior to index date";
%let fn5="e) during pregnancy (from LMP2 until index date)";
%let fn6="f) 12-months lookback from LMP2";
%let fn7="g) excluding bipolar affective disorder and mania";
%let fn8="h) 3-months lookback from LMP2";
%let fn9="i) at mother s LMP2 ^{newline} j) calendar years will be grouped in each country according to the length of the study period";

title &title.;
footnote2 j=l &fn1.;
footnote3 j=l &fn2.;
footnote4 j=l &fn3.;
footnote5 j=l &fn4.;
footnote6 j=l &fn5.;
footnote7 j=l &fn6.;
footnote8 j=l &fn7.;
footnote9 j=l &fn8.;
footnote10 j=l &fn9.;

%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table26 &track_date..xlsx" options(SHEET_NAME="Table26"  embedded_titles='yes' embedded_footnotes='yes');
ods escapechar="^";
proc report data=to_report  split = '*';
	column (('NDD' ('Variable (or interaction)^{super a}' level))
            ('Estimate' 
             (("OR" odds) 
              ("95% CI" CI) 
              ("P-value" pval)  
              )));
    
	define  level   / display ""
                      style(header)={vjust=middle just=center width= 8 cm}
					  style(column)={just=l cellwidth = 8 cm};
	define 	odds	/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm tagattr="format:0.00"};
	define 	CI	/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define 	pval	/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm tagattr="format:0.0000"}; * printing 4 decimal places; 
	compute level;
		if level in ('Offspring risk factors/confounders', 'Maternal risk factors/confounders',
                     'Paternal risk factors/confounders') then do;
            call define (_col_,"style","style={fontweight=bold background=lightgrey}");
		end; 
		if find(level,"Gender",'i') >= 1 or find(level,"age",'i') >= 1 or find(level,"Smoking",'i') >= 1 or
           find(level,"polypharmacy",'i') >= 1 or find(level," offspring conception",'i') >= 1 then do;
            call define (_col_,"style","style={fontweight=bold}");
        end;
    endcomp;

run;
ods excel close;

**********************************************************END*****************************************************************;


