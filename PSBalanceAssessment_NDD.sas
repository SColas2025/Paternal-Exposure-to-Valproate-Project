***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring
*
* This code makes assessment of balance of risk factors/confounders after PS weighting (PS scores obtained using 
* logistic regression) for the primary outcome (NDD including ASD)
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*           Assess absolute standardized difference	
*           Assess variance ratio
*			Table creation
***************************************************************************************************************************;

***************************************************************************************************************************;
*Step 1: Read in dataset;            
****************************************************************************************************************************;
* Read in "PS_weight_logistic" dataset, described in "01_data_description.docx";
data PS_weight_logistic; 
	set dta.PS_weight_logistic;
run;

***************************************************************************************************************************;
*Step 2: Assess absolute standardized difference 		
***************************************************************************************************************************;
* Select weighted patients; 
data weighted;
	set PS_weight_logistic;

	* patients not missing weights are the ones included in the PS model and having no outliers; 
	if not missing(weight);

	* rename stabilized weights; 
	rename weight=sw;

	if PATEXPVALP=0 then PS_group='Lamotrigine/Levetiracetam';
	else if PATEXPVALP=1 then PS_group='Valproate'; 
	proc sort; by PS_group; 
run;
* All patients should have outlier=0; 

* Remove format;
proc datasets lib=work memtype=data noprint;
	modify weighted;
	attrib _all_ format=;
run;
quit;
 
* Make an empty shell for table creation;
data label;
	length level $ 150 index1 $ 30 balance1 $ 30 index2 $ 30 balance2 $ 30 nvar 8;
	infile datalines delimiter='+';
	input level index1 balance1 index2 balance2 nvar;
	datalines;
	Offspring risk factors/confounders +  +  +  +  + 0.5 
	Gender ^{super c} +  +  +  +  + 1
	Congenital CMV ^{super d} +  +  +  +  + 2
	Congenital rubella ^{super d} +  +  +  +  + 3
	Foetal alcohol syndrome ^{super d} +  +  +  +  + 4
	Fragile X syndrome ^{super d} +  +  +  +  + 5
	Lejeune/cri du chat syndrome ^{super d} +  +  +  +  + 6
	Tuberous sclerosis ^{super d} +  +  +  +  + 7
	Maternal risk factors/confounders +  +  +  +  + 7.5 
	Mother's age ^{super c} (categorical) +  +  +  +  + 8
	Affective disorder ^{super e} +  +  +  +  + 9
	Diabetes ^{super e} +  +  +  +  + 10
	Gestational diabetes ^{super f} +  +  +  +  + 11
	Neurotic disorder ^{super e} +  +  +  +  + 12
	Schizophrenia, schizotypal and delusional disorders ^{super e} +  +  +  +  + 13
	Obesity ^{super g} +  +  +  +  + 14
	CMV ^{super g} +  +  +  +  + 15
	Rubella ^{super g} +  +  +  +  + 16
	Alcohol abuse prior to LMP2 ^{super g} +  +  +  +  + 17
	Alcohol abuse during pregnancy ^{super f} +  +  +  +  + 18
	Substance abuse prior to LMP2 ^{super g} +  +  +  +  + 19
	Substance abuse during pregnancy ^{super f} +  +  +  +  + 20
	Smoking prior to LMP2 ^{super g} +  +  +  +  + 21
	Smoking during pregnancy ^{super f} +  +  +  +  + 22
	Maternal polypharmacy index prior to LMP2 ^{super i} (categorical) +  +  +  +  + 23
	Maternal polypharmacy index during pregnancy ^{super f} (categorical) +  +  +  +  + 24
	Concomitant medications associated with valproate-indicated psychiatric conditions prior to LMP2 ^{super g} - mothers with at least one prescription +  +  +  +  + 25
	Concomitant medications associated with valproate-indicated psychiatric conditions during pregnancy ^{super f} - mothers with at least one prescription +  +  +  +  + 26
	Concomitant medications associated with neuropsychiatric adverse events prior to LMP2 ^{super g} - mothers with at least one prescription +  +  +  +  + 27
	Concomitant medications associated with neuropsychiatric adverse events during pregnancy ^{super f} - mothers with at least one prescription +  +  +  +  + 28
	Paternal risk factors/confounders +  +  +  +  + 28.5
	Affective disorder ^{super e,h} +  +  +  +  + 29
	Bipolar affective disorder ^{super e} +  +  +  +  + 30
	Mania ^{super e} +  +  +  +  + 31
	Neurotic disorder ^{super e} +  +  +  +  + 32
	Schizophrenia, schizotypal and delusional disorders ^{super e} +  +  +  +  + 33
	Substance abuse ^{super g} +  +  +  +  + 34
	Paternal polypharmacy index ^{super i} (categorical) +  +  +  +  + 35
	Concomitant medications associated with valproate-indicated psychiatric conditions ^{super g} - fathers with at least one prescription +  +  +  +  + 36
	Concomitant medications associated with neuropsychiatric adverse events ^{super g} - fathers with at least one prescription +  +  +  +  + 37
	Father's age ^{super c} (categorical) +  +  +  +  + 38
	Year of offspring conception ^{super j} +  +  +  +  + 39
	;
run;


* Make an empty binomial proportion;
* Print this when either valproate/comparator group does not have value = 1 in a binary variable;
data empty_freq; 
	length PS_group $30.;
    infile datalines delimiter=',';
	input PS_group _bin_;
	datalines; 
	Lamotrigine/Levetiracetam, 0 
    Valproate, 0
	;
run;

%let vars=/* offspring factors */ 
            gender  O_CONCMV_1YR  O_CONRUBE_1YR  O_FAS_EXIT  O_FRAGILEX_EXIT  O_LEJEUNE_EXIT  O_TUBSIS_EXIT
            /* maternal factors */ 
            M_AGECATINDEX  M_AFFECTDIS_INDEX  M_DIAB_INDEX  M_GESTDIAB_PREG  M_NEURODIS_INDEX  M_SCHIZO_INDEX 
            M_OBES_12MLB  M_CMV_PREG  M_RUBELLA_PREG  M_ALCABUSE_12MLB  M_ALCABUSE_PREG  M_SUBSABUSE_12MLB  M_SUBSABUSE_PREG 
            M_SMOKE_12MLB  M_SMOKE_PREG  M_PPICAT_3MLB  M_PPICAT_PREG
            M_VALPSYCH_12MLB  M_VALPSYCH_PREG  M_NEUROPSYCH_12MLB  M_NEUROPSYCH_PREG
			/* paternal factors */
            F_AFFDISEXBM_LMP   F_BIPOLAR_LMP    F_MANIA_LMP  F_NEURODIS_LMP  F_SCHIZO_LMP  F_SUBSABUSE_12MLB
            F_PPICAT_3MLB    F_VALPSYCH_12MLB    F_NEUROPSYCH_12MLB  F_AGECATINDEX  F_YRCATCONCEPT; 

* type = 1 : binary variable, 2 = categorical variable with >2 levels;
%let types=1 1 1 1 1 1 1 
           2 1 1 1 1 1 
           1 1 1 1 1 1 1 
           1 1 2 2 
           1 1 1 1 
           1 1 1 1 1 1 
           2 1 1 2 2; 
* Macro to calculate the standardized difference;
options minoperator mindelimiter=',';
%macro check_stand_diff; 
    %do i = 1 %to %sysfunc(countw(&vars.));
    %let var=%upcase(%scan(&vars.,&i.)); 
    %let type=%scan(&types.,&i.);
 
    %put &var.;
    * balance assessment for a binary variable;
    %if &type=1 %then %do;

	    proc sql noprint; 
			select count(distinct &var.) into:N_cats from weighted; * the number of variable levels; 
			select max(&var.) into:max_level_val from weighted where PATEXPVALP=1; * maximum variable level value in valproate group; 
			select max(&var.) into:max_level_comp from weighted where PATEXPVALP=0; * maximum variable level value in comparator group; 
	    quit; 

      %if "&var"="GENDER" %then %do;  * Gender variable has female=2 as reference; 
		proc freq data=weighted noprint;
			by PS_group; 
   			tables &var. / binomial(level='2' CL=wald); /* get the binomial proportion of gender = female */
			weight sw; * PS weighting: get weighted proportion; 
    		output out=freq_&var. binomial;
		run;
      %end; 
	  %else %if &N_cats=2 and &max_level_val=1 and &max_level_comp=1 %then %do; * if 2 levels exist in both valproate and comparator group, then get the binomial proportion;  

        proc freq data=weighted noprint;
			by PS_group; 
   			tables &var. / binomial(level='1' CL=wald); /* get the binomial proportion of variable = yes */
			weight sw; * PS weighting: get weighted proportion;
    		output out=freq_&var. binomial;
		run;

	  %end; 
	  %else %if &N_cats=2 and &max_level_val=1 and &max_level_comp=0 %then %do; * if the valproate group has 2 levels, but the comparator group has 1 level only;

        %put &var. = 1 does not exist in the comparator group.;
		* calculate binomial proportion only for the valproate group;
        proc freq data=weighted noprint;
            by PS_group; 
   			tables &var. / binomial(level='1' CL=wald); /* level = '1' = YES. Get the proportion of "YES" */
			weight sw; * PS weighting: get weighted proportion;
   			output out=Freq_&var. binomial;
			where PATEXPVALP=1; /* valproate group only */
		run;

		* Create proportion = 0 for comparator; 
		data Freq_&var.;
		    length PS_group $25.;
			set Freq_&var.; 
			output; 
    		if _n_=1 then do;
			call missing(of _all_);
			PS_group='Lamotrigine/Levetiracetam';
			_BIN_=0;
			output; 
			end; 
		run; 

	  %end; 
	  %else %if &N_cats=2 and &max_level_val=0 and &max_level_comp=1 %then %do; * if the comparator group has 2 levels, but the valproate group has 1 level only;
 
	    %put &var. = 1 does not exist in the valproate group.;
		* calculate binomial proportion only for the comparator group;
        proc freq data=weighted noprint;
            by PS_group; 
   			tables &var. / binomial(level='1' CL=wald); /* level = '1' = YES. Get the proportion of "YES" */
			weight sw; * PS weighting: get weighted proportion;
   			output out=Freq_&var. binomial;
			where PATEXPVALP=0; /* comparator group only */
		run;

		* create proportion = 0 for comparator; 
		data Freq_&var.;
		    length PS_group $25.;
			set Freq_&var.; 
			output; 
    		if _n_=1 then do;
			call missing(of _all_);
			PS_group='Valproate';
			_BIN_=0;
			output; 
			end; 
		run; 

	  %end;
      %else %if &N_cats=1 %then %do;  * if the binary variable has only one level in the weighted data;

	  %put &var. has only 1 level in the weighted patient data. The standardized difference is not calculated.;
		data freq_&var.;
			set empty_freq;
		run; 
	  %end; 


	  * Take the proportion of valproate and comparator groups; 
        data _NULL_;
			set freq_&var.;
			if PS_group='Valproate' then call symputx("p_val", _bin_);
			if PS_group='Lamotrigine/Levetiracetam' then call symputx("p_comp", _bin_);
		run;
		
		* the numerator and denominator of the standardized difference; 
		%let numerator = %sysevalf(&p_val-&p_comp);
		%let denominator = %sysfunc(sqrt(%sysevalf(&p_val*(1-&p_val)+&p_comp*(1-&p_comp))/2)); 
		
			* calculate the absolute standardized difference;
			%if &denominator>0 %then %let stand_diff = %sysfunc(abs(%sysevalf(&numerator/&denominator))); 
			%else %if &denominator=0 %then %do; * division by 0 is not possible; 
		        %put Note: The denominator is 0 for the standardized difference of &var.. It is not possible to calculate the standardized difference.;  
				%let stand_diff = .; * division by 0 is not possible;
		    %end; 
	%end;


    * balance assessment for a categorical variable with >2 levels;
   %if &type=2 %then %do;

	   data freq_&var.;
			set empty_freq; 
	   run; 

	   * Mahalanobis distance;
	   proc DISCRIM data=weighted 
	 		distance anova MANOVA CROSSLISTERR;
			class PATEXPVALP;
			var &var.;
			weight sw; * PS weighting: use stabilized weights; 
			ods output  DistGeneralized=MD_&var.;
	   run;

	   * rename variable 0 to _0 and 1 to _1. SAS 9.04.01M7 creates variable names 0 and 1 and these should be renamed; 
	   options validvarname=any dkrocond=nowarn;
	   data MD_&var.;
	   		set MD_&var.;
			rename '0'n=_0 '1'n=_1;
	   run;
       options validvarname=V7 dkrocond=WARN; * set the options back to default;

	   * extract Mahalanobis distance and name it "the standardized difference"; 
       data _NULL_;
			set MD_&var. (rename=(_0=psm));
			if FromPATEXPVALP="1" then call symputx("stand_diff", psm);
	   run; 

   %end; 
	
		* print the standardized difference and balance assessment;
		data label;
			set label; 
			if nvar=&i. then do; 
				if &type.=1 then index1=put(&stand_diff, 10.2);  * standardized difference for a binary variable;
				else if &type.=2 then index1=put(&stand_diff, 10.2)|| "*";  * Mahalanobis distance, add a footnote that this is mahalanobis;

			* (applied to binary variables) if the absolulte standardized difference is below 0.1, then it is a good balance; 
		    if strip(index1)='.' then balance1='- ^{super **}'; * add a footnote explaining why the standardized difference is not calculated; 
		    else if input(index1, 10.2)<=0.1 then balance1='Yes';
		    else if input(index1, 10.2)>0.1 then balance1='No';  
            end; 
		run;   
  %end;  
%mend; 
%check_stand_diff;  

***************************************************************************************************************************;
*Step 3: Assess variance ratio 		
***************************************************************************************************************************;
* In this macro, the variance ratio is calculated;
%let vars=/* offspring factors */ 
            gender  O_CONCMV_1YR  O_CONRUBE_1YR  O_FAS_EXIT  O_FRAGILEX_EXIT  O_LEJEUNE_EXIT  O_TUBSIS_EXIT
            /* maternal factors */ 
            M_AGECATINDEX  M_AFFECTDIS_INDEX  M_DIAB_INDEX  M_GESTDIAB_PREG  M_NEURODIS_INDEX  M_SCHIZO_INDEX 
            M_OBES_12MLB  M_CMV_PREG  M_RUBELLA_PREG  M_ALCABUSE_12MLB  M_ALCABUSE_PREG  M_SUBSABUSE_12MLB  M_SUBSABUSE_PREG 
            M_SMOKE_12MLB  M_SMOKE_PREG  M_PPICAT_3MLB  M_PPICAT_PREG
            M_VALPSYCH_12MLB  M_VALPSYCH_PREG  M_NEUROPSYCH_12MLB  M_NEUROPSYCH_PREG
			/* paternal factors */
            F_AFFDISEXBM_LMP   F_BIPOLAR_LMP    F_MANIA_LMP  F_NEURODIS_LMP  F_SCHIZO_LMP  F_SUBSABUSE_12MLB
            F_PPICAT_3MLB    F_VALPSYCH_12MLB    F_NEUROPSYCH_12MLB  F_AGECATINDEX  F_YRCATCONCEPT;
%macro check_var_ratio;

	%do i = 1 %to %sysfunc(countw(&vars.));
    %let var=%scan(&vars.,&i.); 

	proc means data=weighted noprint;
		by PS_group;
		var &var.;
		weight sw; * PS weighting: get weighted variance; 
		output out=var_&var. var=var;
	run;

    * variance of a covariate in the valproate/reference group;
	data _NULL_;
		set var_&var.;
		if PS_group='Valproate' then call symputx("v_val", var); * v_val = numerator; 
		if PS_group='Lamotrigine/Levetiracetam' then call symputx("v_comp", var); * v_comp = denominator; 
	run;

	* calculate the variance ratio;
	* variance ratio = the ratio of the sample variance of each covariate for the weighted valproate and comparator group; 
	%if &v_comp>0 and &v_val>0 %then %do; 
		%let var_ratio = %sysevalf(&v_val/&v_comp);
	%end; 
	%else %if &v_comp=0 %then %do; * cannot calculate the variance ratio when the denominator=0; 
		%put Note: the variance of &var. is 0 in comparator group. It is not possible to calculate the variance ratio.; 
    	%let var_ratio = .; 
    %end;  
    %else %if &v_val=0 %then %do; * cannot calculate the variance ratio when the denominator=0; 
		%put Note: the variance of &var. is 0 in valproate group. It is not possible to calculate the variance ratio.; 
    	%let var_ratio = .;   
    %end; 
	
	* print the variable ratio and balance assessment;
	data label;
		set label; 
		if nvar=&i. then do; 
		index2=put(&var_ratio, 10.2);
		
		* if the variance ratio is between 0 and 2, it is a good balance;  
		if strip(index2)='.' then balance2='- ^{super ***}';
		else if 0<=input(index2, 10.2)<=2 then balance2='Yes';
		else if input(index2, 10.2)>2 then balance2='No';  
		
    	end; 
	run;

%end; 
%mend;
%check_var_ratio; 		
  
* Mask standardized differences and variance ratios which are not calculated;
* This can happen when a comorbidity / medication has 0 hits;   
data to_report; 
	set label; 
	array cols index1 balance1 index2 balance2; 
	do over cols;
		if nvar not in (0.5, 7.5, 28.5) then do; 
			if missing(cols) or strip(cols)='.' then cols='-';
		end;  
	end; 
	* add this code to keep two decimals;
	if index1 not in ('-', '') then index1 = cats(index1, '^{}');
run;  

***************************************************************************************************************************;
*Step 4: Table creation	
***************************************************************************************************************************;
** Set title **;
options noquotelenmax;
%let title="Table 30. Balance of risk factors/confounders after PS weighting (PS scores obtained using logistic regression); primary outcome";
%let fn1="a) absolute standardized difference below 0.1";
%let fn2="b) variance ratio between 0 and 2";
%let fn3="c) at index (childbirth)";
%let fn4="d) between index and exit date";
%let fn5="e) all available data prior to index date";
%let fn6="f) during pregnancy (from LMP2 until index date)";
%let fn7="g) 12-months lookback from LMP2";
%let fn8="h) excluding bipolar affective disorder and mania";
%let fn9="i) 3-months lookback from LMP2";
%let fn10="j) at mother's LMP2
^{newline} * Mahalanobis distance is calculated for categorical variables with more than 2 levels.
^{newline} ** The standardized difference is not calculated if a binary variable has only 1 category level in the weighted patient data.
^{newline} *** The variance ratio is not calculated if a variable has only 1 category level in one of valproate and comparator groups (the denominator of the variance ratio is 0)."; 


title &title.;
footnote1 j=l &fn1.;
footnote2 j=l &fn2.;
footnote3 j=l &fn3.;
footnote4 j=l &fn4.;
footnote5 j=l &fn5.;
footnote6 j=l &fn6.;
footnote7 j=l &fn7.;
footnote8 j=l &fn8.;
footnote9 j=l &fn9.;
footnote10 j=l &fn10.;

%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table30 &track_date..xlsx" options(SHEET_NAME="Table30"  embedded_titles='yes' embedded_footnotes='yes');
ods escapechar="^";
proc report data=to_report  split = '*';
	column ( 
             ("NDD" level) 
             ("Absolute standardized difference" index1) 
             ("Balanced achieved ^{super a}" balance1)
             ("Variance ratio * (valproate vs lamotrigine/levetiracetam)" index2) 
             ("Balanced achieved ^{super b}" balance2) 
           );
    
	define level    / display "" 
                      style(header)={vjust=middle just=center width= 8 cm}
					  style(column)={just=l cellwidth = 8 cm};
	define index1		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm tagattr="Format:0.00"}; * using the format to create 2 decimal places; 
	define balance1		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define index2		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm tagattr="Format:0.00"}; * using the format to create 2 decimal places; 
	define balance2		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm}; 

	compute level;
		if level in ("Offspring risk factors/confounders", "Maternal risk factors/confounders", "Paternal risk factors/confounders") then do;
            call define (_col_,"style","style={fontweight=bold background=lightgrey}");
        end;
		
    endcomp;

run;
ods excel close;
  
**********************************************************END*****************************************************************;


