***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring
*
* This code makes variable estimates from logistic regression propensity score model
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*           Assessment of association between the potential confounders and risk factors and the outcome
*           Perform logistic regression PS model, Cook's distance, excluding influential data
*			Re-estimate PS model (logistic regression)
*			PS weighting
***************************************************************************************************************************;


***************************************************************************************************************************;
*Step 1: Read in dataset;            
****************************************************************************************************************************;
* Read in "mydata" dataset, described in "01_data_description.docx";
data in_data; 
	set dta.mydata;
run;

***************************************************************************************************************************;
*Step 2: Assessment of association between the potential confounders and risk factors and the outcome	
***************************************************************************************************************************;
* Select cohort for the analyses; 
data PS_data;
	set in_data; 

	if C1_COMP=1; 

	* do not include 'missing' category of smoking in the analysis;	
    if M_SMOKE_12MLB=9 then M_SMOKE_12MLB=.;
	if M_SMOKE_PREG=9 then M_SMOKE_PREG=.; 
run; 

* Remove all formats, so that it is easier to define reference groups for logistic regression PS model;
proc datasets lib=work memtype=data noprint;
	modify PS_data;
	attrib _all_ format=;
run;
quit;

* Select candidate covariates; 
* Candidate covariates will be considered to enter the model if the odds ratio in a single association 
with the outcome is higher than 1.1 or lower than 0.9
for categorical variables with more than two categories (i.e. more than one OR estimated) 
the Wald test will be evaluated to assess whether the variable is significantly associated with the outcome.;

/* 
Variables that are evaluated:

21 Risk factors for mother
-	Age (index)
-	Obesity (12 months look back from LMP2)
-	Smoking (12 months look back from LMP2)
-	Smoking during pregnancy
-	Substance abuse (12 months look back from LMP2)
-	Substance abuse during pregnancy
-	Alcohol abuse (12 months look back from LMP2)
-	Alcohol abuse during pregnancy
-	Schizophrenia, schizotypal and delusional disorders (ever)
-	Affective Disorder (ever)
-	Neurotic Disorder (ever)   
-	Rubella (during pregnancy)
-	CMV (during pregnancy)
-	Diabetes (ever)           
-	Gestational Diabetes (during pregnancy)
-	Any concomitant medications associated with valproate-indicated psychiatric conditions (12 months look back from LMP2)
-	Any concomitant medications associated with valproate-indicated psychiatric conditions during pregnancy
-	Any concomitant medications associated with neuropsychiatric adverse effects (12 months look back from LMP2)
-	Any concomitant medications associated with neuropsychiatric adverse effects during pregnancy
-	Maternal polypharmacy index (3 months look back from LMP2)
-	Maternal polypharmacy index during pregnancy

11 risk factors for father
-	Substance Abuse (12 months look back from LMP2) 
-	Affective Disorders (excluding bipolar and mania) (ever) 
-	Schizophrenia, schizotypal and delusional disorders (ever) 
-	Neurotic Disorder (ever)   
-	Any concomitant medications associated with valproate-indicated psychiatric conditions (12 months look back from LMP2)
-	Any concomitant medications associated with neuropsychiatric adverse effects (12 months look back from LMP2)
-	Age (index)
-	Bipolar Affective Disorder (ever)
-	Mania (ever)
-	Year of offspring conception (mother s LMP2)
-	Paternal polypharmacy index (3 months look back from LMP2)

7 risk factors for offspring
-	Gender
-	Fragile X Syndrome
-	Lejeune/ cri du chat syndrome
-	Tuberous Sclerosis
-	Congenital CMV
-	Congenital Rubella
-	Foetal Alcohol syndrome
 
*/

%let ds=PS_data; 
%let covars= /* offspring factors */ 
            gender  O_CONCMV_1YR  O_CONRUBE_1YR  O_FAS_EXIT  O_FRAGILEX_EXIT  O_LEJEUNE_EXIT  O_TUBSIS_EXIT
            /* maternal factors */  
            M_AGECATINDEX  M_AFFECTDIS_INDEX  M_DIAB_INDEX  M_GESTDIAB_PREG  M_NEURODIS_INDEX  M_SCHIZO_INDEX 
            M_OBES_12MLB  M_CMV_PREG  M_RUBELLA_PREG  M_ALCABUSE_12MLB  M_ALCABUSE_PREG  M_SUBSABUSE_12MLB  M_SUBSABUSE_PREG 
            M_SMOKE_12MLB  M_SMOKE_PREG  M_PPICAT_3MLB  M_PPICAT_PREG
            M_VALPSYCH_12MLB  M_VALPSYCH_PREG  M_NEUROPSYCH_12MLB  M_NEUROPSYCH_PREG

			/* paternal factors */
            F_AFFDISEXBM_LMP   F_BIPOLAR_LMP    F_MANIA_LMP  F_NEURODIS_LMP  F_SCHIZO_LMP  F_SUBSABUSE_12MLB
            F_PPICAT_3MLB    F_VALPSYCH_12MLB    F_NEUROPSYCH_12MLB  F_AGECATINDEX  F_YRCATCONCEPT; 
%let candidates=; 
%macro check_odds();

	%do i=1 %to %sysfunc(countw(&covars.));
        %let var=%upcase(%scan(&covars., &i.));
		%put &var.; 

		* count how many categories exist for a variable; 
		proc sql noprint; 
			select count (distinct &var.) into:N_cat
			from &ds.; 
		quit; 

        * variables with more than 2 categories; 
		%if &N_cat>=3 %then %do;
 
		    * the Wald test; 
			proc logistic data=&ds.;
				class &var.(ref=first)/param=ref; *set lowest level as reference;
				model O_NDDALL(event='1')=&var./ridging=ABSOLUTE;
				ods output GlobalTests=GlobalTests&i.;
			run;		

			proc sql noprint; select ProbChiSq format=8.4 into:wald_pval from GlobalTests&i. where test='Wald'; quit;

			* if the p-value is less than 0.05 (significant), add this variable as a candidate; 
			%if &wald_pval<0.05 %then %do; 
				%let candidates=&candidates &var.;
			%end; 

		%end; 
		%else %if &N_cat=2 %then %do; 

		    * Odds ratio; 
			proc logistic data=&ds.;
				class &var.(ref=first)/param=ref; *set lowest level as reference;
				model O_NDDALL(event='1')=&var./ridging=ABSOLUTE;
				ods output OddsRatios=singleOdds&i.;
			run;

			
            proc sql noprint; select abs(OddsRatioEst-1) into:OR_abs from singleOdds&i.; quit; 
			
            * if the OR is >1.1 or <0.9, add this variable as a candidate;	
            %if &OR_abs>0.1 %then %do; 
				%let candidates=&candidates &var.;
			%end; 

		%end; 
		%else %if 0<=&N_cat<=1 %then %do; 
			%put &var. has less than 2 categories available. It is not possible to perform logistic regression PS model.;  
		%end; 

	%end; 
%mend; 
%check_odds();

* Check candidate covariates for logistic regression PS model; 
%put &candidates; 

* Select reference groups in logistic regression PS model; 
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

* Keep patient id and candidates;
data PSmodel; 
	set &ds.;
	keep OFFID MATID PATID PATEXPVALP PATEXPCOMP &candidates.;
run;

***************************************************************************************************************************;
*Step 3: Perform logistic regression PS model, Cook's distance, excluding influential data
***************************************************************************************************************************;
/*To identify outlying values of variables, such as gestational age and weight, the Cook s distance will be calculated for 
each offspring after fitting the PS model. Offspring that have a Cook s distance greater than four times the mean will be 
classified as influential and excluded from the model (and subsequent PS weighting), before re-estimating the model. This 
assessment will be repeated for each PS model*/

* Detect quasi-complete separation; 
%macro check_sep(covar=, data=);

%global newvar;
%let newvar=&covar.;

%do i=1 %to %sysfunc(countw(&covar.));
        %let var=%scan(&covar.,&i., ' ');
		%put &var.;

		proc sql;
         	create table cat_valp as select distinct &var. from &data. where PATEXPVALP=1 and not missing(&var.); 
			create table cat_comp as select distinct &var. from &data. where PATEXPVALP=0 and not missing(&var.); 
		quit; 

	    proc sql noprint; 
		    select count(*) into:N_cat_valp from cat_valp;
			select count(*) into:N_cat_comp from cat_comp;
        quit; 

		proc compare base=cat_valp
             compare=cat_comp noprint;
        run;

    %if &N_cat_valp=1 %then %do;
		%put &var. has only one level in valproate group and is removed from covariates.;
		* remove the variable from candidates;
		%let newvar = %sysfunc(tranwrd(&newvar., &var., %str()));
    %end;

	%if &N_cat_comp=1 %then %do;
		%put &var. has only one level in comparator group and is removed from covariates.;
		* remove the variable from candidates;
		%let newvar = %sysfunc(tranwrd(&newvar., &var., %str()));
    %end;

	 * check if there is an overlap in categories between valproate and comparator group; 
	* &sysinfo=0 means cat_valp = cat_comp. The categories in valproate and comparator groups are the same.; 
    %if "&sysinfo" ne "0" %then %do;
		%put Quasi-complete separation detected: &var. predicts paternal exposure perfectly. This variable is removed from covariates.;
		%let newvar = %sysfunc(tranwrd(&newvar., &var., %str()));
    %end;

	proc delete data=cat_valp cat_comp; run;
 
%end; 
	    
%mend;

%check_sep(covar=&candidates., data=PSmodel);

%let candidates_revised_1=&newvar.;
%symdel newvar;

* compare old and new candidate variables; 
%put &candidates.; 
* variables which can cause quasi-complete separation are removed;  
%put &candidates_revised_1.;

* Initial logistic regression PS model;   
* Calculate cook's distance; 
* if there is still quasi-complete separation, the following warnings can happen 
WARNING: The negative of the Hessian is not positive definite. The convergence is questionable.
WARNING: The procedure is continuing but the validity of the model fit is questionable.
WARNING: The specified model did not converge.; 
proc genmod data=PSmodel;
	class PATEXPVALP; 
	model PATEXPVALP(event='1') = &candidates_revised_1./dist=bin link=logit;
	output out=stdres cooksd=cookd;
run; 

* Outliers are defined as cook's distance greater than 4 times the mean; 
proc sql;
	select mean(cookd)*4 into:threshold from stdres;  
quit; 

* Estimate propensity scores and cook's distance;
data PSmodel_2;
	set stdres;
 
	* 0 = include in the PS model, 1 = outlier, exclude from the PS model; 
	if not missing(cookd) then do; 
		if cookd>&threshold then exclude=1;
		else if cookd<=&threshold then exclude=0;
	end; 

	* when an observation is missing some variable values, proc genmod does not use this observation, and no cookd printed; 
	if missing(cookd) then exclude=.;
run;

* Observations are missing smoking variables; 
* They are excluded from proc genmod, no cook's distance for them; 
data miss; 
	set stdres;
    if missing(M_SMOKE_PREG); 
run; 

***************************************************************************************************************************;
*Step 4: Re-estimate PS model (logistic regression)	
***************************************************************************************************************************;
* Check for quasi-complete separation;
%check_sep(covar=&candidates_revised_1., data=PSmodel_2 (where=(exclude=0)) );

%let candidates_revised_2=&newvar.;
%symdel newvar;

* compare old and new candidate variables; 
%put &candidates_revised_1.;
* variables which can cause quasi-complete separation are removed;  
%put &candidates_revised_2.;


* new reference groups; 
%let refs=; 
%get_ref(vars=&candidates_revised_2., ds=PSmodel_2 (where=(exclude=0)) ); 

* check reference groups for logistic regression PS model; 
%put &refs; 

* Re-estimation of PS model: logistic regression with candidate covariates; 
* Valproate = event (1), lamotrigine/levetiracetam = reference (0);  
ods output OddsRatios=Odds ParameterEstimates=Param; 
proc logistic data=PSmodel_2 namelen=32;
	class &refs. /param=ref;  
	model PATEXPVALP(event='1') = &candidates_revised_2.; 
	output out=propen prob=prob_yes_valp;  
	where exclude=0;  * important: exclude outliers; 
run;

* Histogram of propensity scores;
* valproate and comparator groups should overlap;  
proc sgplot data=propen;
	histogram prob_yes_valp / group=PATEXPVALP transparency=0.5;      
	density prob_yes_valp / type=kernel group=PATEXPVALP; 
run;

***************************************************************************************************************************;
*Step 5: PS weighting
***************************************************************************************************************************;
/*PS weighting is used instead of matching for increased generalisability and to avoid the exclusion of patients from 
the adjusted analyses due to lack of matches; this is deemed as particularly important in situations where the outcome 
of interest is relatively rare or infrequent. 
Stabilised weights will be calculated as follows: denoting paternal exposure by Z (with Z=1 indicating paternal exposure to 
valproate and Z=0 indicating paternal exposure to lamotrigine/levetiracetam) and the propensity score as e= Pr(Z=1|X) (X is 
the set of covariates included in the PS model), stabilised weights will be obtained as:
w=Pr?(Z=1)/e  (if Z=1)
w=Pr?(Z=0)/(1-e)  (if Z=0)
With Pr(Z=1) and Pr(Z=0) representing the marginal probability of paternal exposure to valproate (Z=1) or 
lamotrigine/levetiracetam (Z=0), obtained as the proportion of offspring in the two exposure groups in the population of interest. 
Stabilised weights will be preferred to unstabilised weights as PS scores close to 0 or 1 often lead to very large unstabilised weights.
*/


* Pr(Z=1) = the marginal probability of paternal exposure to valproate (Z=1);
* This is proportion of valproate exposure, after excluding outliers (i.e., P_val);
proc sql; 
	select sum(PATEXPVALP)/count(*) into:P_val
	from propen;
quit; 
	
* Stabilized weights;
data ps_weight; 
	set propen; 

	* valproate group, if Z=1; 
	if PATEXPVALP=1 then weight=&P_val./prob_yes_valp; * prob_yes_valp = propensity score (e= Pr(Z=1|X)); 
	* lamotrigine/levetiracetam, if Z=0; 
	else if PATEXPVALP=0 then weight=(1-&P_val.)/(1-prob_yes_valp); 
run; 

* Merge cook's distance, propensity score and weights;
* c1_comp = total patients, before outlier exclusion;
* propen = patients after outlier exclusion;   
proc sql; 
	create table PS_weight_logistic as
	select a.*, b.cookd, c.prob_yes_valp, c.weight
	from PS_data a
	left join stdres b
	on a.OFFID=b.OFFID and a.MATID=b.MATID and a.PATID=b.PATID
    left join ps_weight c
    on a.OFFID=c.OFFID and a.MATID=c.MATID and a.PATID=c.PATID; 
quit; 

* All patients with ALL covariates;
data dta.PS_weight_logistic;
	set PS_weight_logistic; 

	* 1 = outlier, 0 = no outlier, . = excluded from logistic regression PS model due to missingness;
	if not missing(cookd) then do; 
		if cookd>&threshold then outlier=1;
		else if cookd<=&threshold then outlier=0;
	end;
    if missing(cookd) then outlier=.; 

	label outlier="1=outlier (excluded from PS weighting), 0=not outlier, missing=missingness in covariates";
run;

**********************************************************END*****************************************************************;


