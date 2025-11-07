***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring 
*
* This code summarises univariate association between all pre-specified maternal risk factors/confounders, 
* and the primary outcome (binary variable indicating if the offspring was diagnosed with NDD including ASD during follow-up)
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*			Association between maternal risk factors/confounders and NDD	
*			Table creation
***************************************************************************************************************************;

***************************************************************************************************************************;
*Step 1: Read in dataset;            
****************************************************************************************************************************;
* Read in "mydata" dataset, described in "01_data_description.docx";
data in_data; 
	set dta.mydata;
run;

* Read in table grid containing summary data, described in "01_data_description.docx";
data sum_data;
	set dta.summary_maternal;
run;

***************************************************************************************************************************;
*Step 2: Association between maternal risk factors/confounders and NDD	
***************************************************************************************************************************;
* Create input dataset for association; 
data mother_part2; 
     set in_data;

	 *select cohort;
    if C1_COMP=1;

	 * create new smoking variables. For them 'missing' category will be excluded from statistical tests;
     M_SMOKE_12MLB_2=M_SMOKE_12MLB;
     M_SMOKE_PREG_2=M_SMOKE_PREG;

	 * replace 9 with .; 
     if M_SMOKE_12MLB_2=9 then M_SMOKE_12MLB_2=.;
	 if M_SMOKE_PREG_2=9 then M_SMOKE_PREG_2=.;

    label   /* mother: age*/
           	M_AGECATINDEX="Mother's age ^{super a} (categorical)"

		   /* mother: comorbidities */
		  	M_AFFECTDIS_INDEX='Affective disorder ^{super b}'
          	M_DIAB_INDEX='Diabetes ^{super b}'
		  	M_GESTDIAB_PREG='Gestational diabetes ^{super c}' 
		  	M_NEURODIS_INDEX='Neurotic disorder ^{super b}'
          	M_SCHIZO_INDEX='Schizophrenia, schizotypal and delusional disorders ^{super b}'
		  	M_OBES_12MLB='Obesity ^{super d}'
          	M_CMV_PREG='CMV ^{super c}'
          	M_RUBELLA_PREG='Rubella ^{super c}'
		 
 			/* mother: lifestyle */
          	M_ALCABUSE_12MLB='Alcohol abuse prior to LMP2 ^{super d}' 
          	M_ALCABUSE_PREG='Alcohol abuse during pregnancy ^{super c}' 
          	M_SUBSABUSE_12MLB='Substance abuse prior to LMP2 ^{super d}'  
          	M_SUBSABUSE_PREG='Substance abuse during pregnancy ^{super c}' 

			/* smoking */
		  	M_SMOKE_12MLB='Smoking prior to LMP2 ^{super e}' 
          	M_SMOKE_PREG='Smoking during pregnancy ^{super c}'

			/* mother: PPI */
		 	M_PPICAT_3MLB='Maternal polypharmacy index prior to LMP2 ^{super e} (categorical)'
			M_PPICAT_PREG='Maternal polypharmacy index during pregnancy ^{super c} (categorical)'

			/* mother: concomitant medication */
          	M_VALPSYCH_12MLB='Concomitant medications associated with valproate-indicated psychiatric conditions prior to LMP2 ^{super d} - mothers with at least one prescription' 
          	M_VALPSYCH_PREG='Concomitant medications associated with valproate-indicated psychiatric conditions during pregnancy ^{super c} - mothers with at least 1 prescription'
          	M_NEUROPSYCH_12MLB='Concomitant medications associated with neuropsychiatric adverse events prior to LMP2 ^{super d} - mothers with at least one prescription' 
          	M_NEUROPSYCH_PREG='Concomitant medications associated with neuropsychiatric adverse events during pregnancy ^{super c} - mothers with at least one prescription';

run;  

* Remove all formats, so that it is easier to define reference groups for logistic regression;
proc datasets lib=work memtype=data noprint;
      modify mother_part2;
      attrib _all_ format=;
run;
quit;

* Create empty oddsratio;
data OR_empty;
    length oddsci $40.; 
	oddsci='-'; 
	OddsRatioEst=.;
    LowerCL=.; 
    UpperCL=.;
run;

* Create empty wald test;
data WT_empty;
    length test $16. tp $40.;
    test="Wald"; 
	tp='-';
    ChiSq=.;
    ProbChiSq=.; 
run;

* Perform logistic regression for estimating the associations;
*There is a warning about Ridging has failed to improve the loglikelihood for one variable in the mockup data,if such warning present
	when running on real data, this can be fixed by using "INEST= option" or "RIDGING=NON"; 
%let devar=M_AGECATINDEX M_AFFECTDIS_INDEX M_DIAB_INDEX M_GESTDIAB_PREG M_NEURODIS_INDEX
			M_SCHIZO_INDEX M_OBES_12MLB M_CMV_PREG M_RUBELLA_PREG M_ALCABUSE_12MLB
			M_ALCABUSE_PREG M_SUBSABUSE_12MLB M_SUBSABUSE_PREG M_SMOKE_12MLB_2 M_SMOKE_PREG_2
			M_PPICAT_3MLB M_PPICAT_PREG M_VALPSYCH_12MLB M_VALPSYCH_PREG M_NEUROPSYCH_12MLB 
			M_NEUROPSYCH_PREG; /*Maternal risk factors/confounders included*/
%let ds=mother_part2;
%macro Logistic;
	%do i = 1 %to %sysfunc(countw(&devar.));
	%let var=%scan(&devar.,&i.);

	* the following warning may occur but can be ignored: 
         WARNING: Output 'GlobalTests' was not created.  Make sure that
         the output object name, label, or path is spelled
         correctly.  Also, verify that the appropriate
         procedure options are used to produce the requested
         output object.  For example, verify that the NOPRINT
         option is not used.;
	* WARNING: There is possibly a quasi-complete separation of data points. The maximum likelihood
         estimate may not exist.
	* WARNING: The LOGISTIC procedure continues in spite of the above warning. Results shown are based
         on the last maximum likelihood iteration. Validity of the model fit is questionable.; 
	proc sql noprint;
      select count(*) into:no_event from &ds. where not missing(&var.) and O_NDDALL=1; 
      select count(*) into:no_nonevent from &ds. where not missing(&var.) and O_NDDALL=0; quit;

	* If no events or non-events, use the empty oddsratio;
	%if &no_event=0 or &no_nonevent=0 %then %do;
	data OddsRatios&i.;
		    length Effect $200.; 
			set OR_empty;
            Effect="&var.";
			nvar=&i.; 
 
            if Effect not in ("M_AGECATINDEX", "M_PPICAT_PREG", "M_PPICAT_3MLB") then do; 
            	category=1;
			end; 
	run; 
	%end;
	%else %if "&var"="M_AGECATINDEX" %then %do;
	* Perform logistic regression for age;
    proc logistic data=&ds.;
		class &var.(ref='3')/param=ref; *set age 26-30 as reference;
		model O_NDDALL(event='1')=&var./ridging=ABSOLUTE ;
		ods output OddsRatios=OddsRatios&i. GlobalTests= GlobalTests&i.;
	run;
	%end; 
	%else %do; 
	* Perform logistic regression for other variables;
	proc logistic data=&ds.;
		class &var.(ref=first)/param=ref; *set lowest level as reference;
		model O_NDDALL(event='1')=&var./ridging=ABSOLUTE ;
		ods output OddsRatios=OddsRatios&i. GlobalTests= GlobalTests&i.;
	run;
	%end; 

	proc sql noprint; select count(*) into:N_odds from OddsRatios&i.; quit;  

	* if oddsratios has 0 observations, use the empty oddsratio; 
	%if &N_odds=0 %then %do;
		data OddsRatios&i.;
		    length variable $200.; 
			set OR_empty;
            Variable="&var.";
			nvar=&i.; 
 
            if Variable not in ("M_AGECATINDEX", "M_PPICAT_PREG", "M_PPICAT_3MLB") then do; 
            	category=1;
			end; 
		run; 
    %end;
    %else %do;
		data OddsRatios&i.; 
			set OddsRatios&i.;
            nvar=&i.;  
		run; 
    %end;  

	* if GlobalTests exists;
	%if %sysfunc(exist(GlobalTests&i.)) %then %do;
    %end;
    %else %do; * if GlobalTests does not exist, use the empty wald test; 
		data GlobalTests&i.;
			set WT_empty;
		run; 
    %end; 

	%end;
%mend Logistic;
%Logistic;

* Merge all OddsRatios from all risk factors/confounders;
data odds;
    length Effect $200 oddsci $40;
	set Oddsratios1-Oddsratios21;
    if not missing(OddsRatioEst) and UpperCL<999 then do; * the limit to the upper limit: 999, upper limit more than 999 is not reasonable; 
		oddsci=cats(put(OddsRatioEst,10.2),'(',put(LowerCL,10.2),',',put(UpperCL,10.2),')'); 
	end; 
run;  

* Get variable name, group category, reference category;
data odds_2;
    length variable $ 200; 
	set odds; 
	* variable name; 
	if not missing(Effect) then do; 
		Variable=strip(scan(Effect,1,' '));

		* category and reference; 
		do i=1 to countw(Effect, " ");
			if scan(Effect, i, " ")="vs" then do; 
			category=input(scan(Effect, i-1, " "), 8.);
			reference=input(scan(Effect, i+1, " "), 8.);
			end; 
		end;
	end; 
    proc sort; by nvar category;
run; 

* Create an empty table shell with category groups for table creation;
data comm2;
    length level $ 200;
	infile datalines delimiter='+';
	input category level nvar;
	datalines;
	1 + <=20 years + 1  
    2 + 21-25 + 1 
    4 + 31-35 + 1  
	5 + 36-40 + 1  
	6 + >40 + 1  
    1 + Yes + 2 
    1 + Yes + 3 
    1 + Yes + 4 
    1 + Yes + 5 
    1 + Yes + 6 
    1 + Yes + 7 
    1 + Yes + 8 
    1 + Yes + 9 
    1 + Yes + 10  
    1 + Yes + 11 
    1 + Yes + 12 
    1 + Yes + 13 
    1 + Yes + 14
    9 + Missing + 14 
    1 + Yes + 15 
    9 + Missing + 15   
    1 + 1-4 + 16 
    2 + 5-10 + 16 
    3 + >10 + 16  
    1 + 1-4 + 17   
    2 + 5-10 + 17
    3 + >10 + 17 
    1 + Yes + 18 
    1 + Yes + 19 
    1 + Yes + 20 
    1 + Yes + 21 
	; 
	proc sort; by nvar category;
run; 

* Merge odds ration with category names;
data OR_all;
	merge comm2 (in=a)
	      odds_2 (in=b);
    by nvar category;
	if a; 
	if missing(oddsci) then oddsci='-';
	keep level oddsci nvar category;
	proc sort; by nvar category;
run;

* Merge all Waldtests;
data WT_all(keep=tp level nvar _line);
	length level $ 200 tp $40;
	set Globaltests1-Globaltests21;
	where test="Wald";
	if not missing(ChiSq) then do; 
		tp=cats(put(ChiSq,10.2),' (',put(ProbChiSq,pvalue9.4), ')');
	end; 
	infile datalines delimiter='+';
	input level nvar _line;
	datalines;
	Wald test + 1 + 9
	Yes + 2 + 3 
	Yes + 3 + 3
	Yes + 4 + 3
	Yes + 5 + 3
	Yes + 6 + 3
	Yes + 7 + 3
	Yes + 8 + 3
	Yes + 9 + 3
	Yes + 10 + 3
	Yes + 11 + 3
	Yes + 12 + 3
	Yes + 13 + 3
	Wald test without 'Missing' category + 14 + 5
	Wald test without 'Missing' category + 15 + 5
	Wald test + 16 + 6
	Wald test + 17 + 6
	Yes + 18 +3
	Yes + 19 +3
	Yes + 20 +3
    Yes + 21 +3
;
run;

* Merge all results;
proc sort data=sum_data; by nvar level;
proc sort data=Or_all; by nvar level;
proc sort data=WT_all; by nvar level; run;

data To_report(drop=c1 c2 c3);
	merge sum_data(in=a)Or_all(in=b) WT_all(in=c);
	by nvar level;
	if a; 
	
	*Set reference and -;
	if (nvar=1 and _line=4) or (nvar in (14,15) and _line=2) or (nvar in (16,17) and _line=2)
		then oddsci="Reference";
		else if (nvar=1 and _line=9) or (nvar in (14,15) and _line in (5,6)) or (nvar in (16,17) and _line=6)
		then oddsci="-";

	if _line>1 and missing(tp) then tp="-";

	if level="No" then do;
		if strip(pct1) in ('100.00%', '0.00%') then oddsci='-'; else oddsci='Reference'; 
		tp="-";
	end;

	proc sort; by nvar _line; 
run;

***************************************************************************************************************************;
*Step 3: Table creation		
***************************************************************************************************************************;
* Get bold labels; 
proc contents data=mother_part2 out=contents noprint; run;
proc sql; 
	select cats('"',LABEL,'"')into: labels separated by ' ' 
	from contents where not missing(LABEL);
quit;
%put &labels;

options noquotelenmax;
** Set title **;
%let title="Table 24. Association between potential maternal risk factors/confounders and NDD; primary outcome";

** Set footnotes **;
%let fn1="Legend: The overall column represents the number and percentage of offspring with each characteristic (percentage is calculated over the total number of offspring). The columns Event/Non-Event represent the number and percentage of events and non-events (NDD) in each subgroup defined by the characteristic (percentage is calculated over the number of offspring with the characteristic, i.e. row percentage). The association between each characteristics and the outcome is tested by fitting a logistic regression model, and the odds ratios (OR) with 95% confidence intervals (CI) and likelihood ratio (Wald) test are reported.";
%let fn2="a) at index (childbirth)";
%let fn3="b) all available data prior to index date";
%let fn4="c) during pregnancy (from LMP2 until index date)";
%let fn5="d) 12-months lookback from LMP2";
%let fn6="e) 3-months lookback from LMP2";

title &title.;
footnote1 j=l &fn1.;
footnote2 j=l &fn2.;
footnote3 j=l &fn3.;
footnote4 j=l &fn4.;
footnote5 j=l &fn5.;
footnote6 j=l &fn6.;

ods listing close;
%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table24 &track_date..xlsx" options(SHEET_NAME="Table24"  embedded_titles='yes' embedded_footnotes='yes');
ods escapechar="^";
proc report data=to_report  split = '*';
	column (
            ("NDD" (" " level)) 
             ("Overall" ("N" k1) ("%" pct1)) 
             ("Event" ("N" k2) ("%" pct2))
             ("Non-event" ("N" k3) ("%" pct3)) 
             ("Association" ("OR (95% CI)" oddsci) ("Test statistics (p-value)" tp)) 
             );
    
	define level    / display ""
                      style(header)={vjust=middle just=center width= 8 cm}
					  style(column)={just=l cellwidth = 8 cm};
	define k1		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define pct1		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define k2		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define pct2		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm}; 
	define k3		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	define pct3		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	
	define oddsci		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};

	define tp		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	
	compute level;
		if level in ("Maternal risk factors/confounders") then do;
            call define (_col_,"style","style={fontweight=bold background=lightgray}");
        end;
		if level in (&labels.) then do;
            call define (_col_,"style","style={fontweight=bold}");
        end;
    endcomp;

run;
ods excel close;

**********************************************************END*****************************************************************;
