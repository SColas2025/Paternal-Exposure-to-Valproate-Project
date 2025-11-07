***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring
*
* This code summarises univariate association between all pre-specified paternal risk factors/confounders, 
* and the primary outcome (binary variable indicating if the offspring was diagnosed with NDD including ASD during follow-up)
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*			Association between paternal risk factors/confounders and NDD	
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
	set dta.summary_paternal;
run;

***************************************************************************************************************************;
*Step 2: Association between paternal risk factors/confounders and NDD	
***************************************************************************************************************************;
* Create input dataset for association; 
data father_part2; 
     set in_data;
	 *select cohort;
	 if C1_COMP=1;

    label  /* comorbidities */
			F_AFFDISEXBM_LMP='Affective disorder excluding bipolar affective disorder and mania ^{super a}' 
			F_BIPOLAR_LMP='Bipolar affective disorder ^{super a}' 
			F_MANIA_LMP='Mania ^{super a}'
			F_NEURODIS_LMP='Neurotic disorder ^{super a}' 
			F_SCHIZO_LMP='Schizophrenia, schizotypal and delusional disorders ^{super a}'

			/* lifestyle */
	      	F_SUBSABUSE_12MLB='Substance abuse ^{super c}'

			/* PPI */
          	F_PPICAT_3MLB='Paternal polypharmacy index ^{super d} (categorical)'

			/* father: concomitant medication */
          	F_VALPSYCH_12MLB='Concomitant medications associated with valproate-indicated psychiatric conditions ^{super c} - fathers with at least one prescription'   
          	F_NEUROPSYCH_12MLB='Concomitant medications associated with neuropsychiatric adverse events ^{super c} - fathers with at least one prescription'

			/* father: age*/
           	F_AGECATINDEX="Father's age ^{super e} (categorical)" 

			/*Year of offspring conception*/
          	F_YRCATCONCEPT='Year of offspring conception ^{super f,g}';
run; 

* Remove all formats, so that it is possible to extract category & reference for odds ratios; 
proc datasets lib=work memtype=data noprint;
	modify father_part2;
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
%let devar= F_AFFDISEXBM_LMP F_BIPOLAR_LMP F_MANIA_LMP F_NEURODIS_LMP F_SCHIZO_LMP
			F_SUBSABUSE_12MLB F_PPICAT_3MLB F_VALPSYCH_12MLB F_NEUROPSYCH_12MLB
			F_AGECATINDEX F_YRCATCONCEPT; /*Paternal risk factors/confounders included*/
%let ds=father_part2;
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
     %if &no_event=0 or &no_nonevent=0  %then %do;
     data OddsRatios&i.;
		    length Effect $200.; 
			set OR_empty;
            Effect="&var.";
            nvar=&i.; 

            if Effect not in ("F_PPICAT_6MLB", "F_AGECATINDEX", "F_YRCATCONCEPT") then do; 
            	category=1;
			end; 
	  run; 
    %end;
	%else %if "&var"="F_AGECATINDEX" %then %do;
	* Perform logistic regression for age;
    proc logistic data=&ds.;
		class &var.(ref='4')/param=ref; *set age 31-35 as reference;
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

            if Variable not in ("F_PPICAT_3MLB", "F_AGECATINDEX", "F_YRCATCONCEPT") then do; 
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
	set Oddsratios1-Oddsratios11;
    if not missing(OddsRatioEst) and UpperCL<999 then do;  * the limit to the upper limit: 999, upper limit more than 999 is not reasonable;
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
	1 + Yes + 1
	1 + Yes + 2
	1 + Yes + 3
	1 + Yes + 4
    1 + Yes + 5
	1 + Yes + 6 
    1 + 1-4 + 7   
    2 + 5-10 + 7  
    3 + >10 + 7
    1 + Yes + 8
    1 + Yes + 9 
    1 + <=20 years + 10
    2 + 21-25 + 10  
    3 + 26-30 + 10  
    5 + 36-40 + 10 
    6 + >40 + 10  
    3 + 2011-2015 + 11 
	4 + 2016-2019 + 11
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
	keep level category oddsci nvar;
	proc sort; by nvar category; 
run;

* Merge all Waldtests;
data WT_all(keep=tp level nvar _line);
	length level $ 200 tp $40;
	set Globaltests1-Globaltests11; 
	where test="Wald";
	if not missing(ChiSq) then do;
		tp=cats(put(ChiSq,10.2),' (',put(ProbChiSq,pvalue9.4), ')');
	end; 
	infile datalines delimiter='+';
	input level nvar _line;
	datalines;
	Yes + 1 + 3
    Yes + 2 + 3 
	Yes + 3 + 3
	Yes + 4 + 3
	Yes + 5 + 3
    Yes + 6 + 3
	Wald test + 7 + 6
	Yes + 8 + 3
	Yes + 9 + 3
	Wald test + 10 + 9
	Wald test + 11 + 5
    ;
run;

* Merge all results;
proc sort data=sum_data; by nvar level;
proc sort data=Or_all; by nvar level;
proc sort data=WT_all; by nvar level;
run;

data To_report(drop=c1 c2 c3);
	merge sum_data(in=a )Or_all(in=b) WT_all(in=c);
	by nvar level;
	if a;
	
	*Set reference and -;
	if (nvar=7 and _line=2) or (nvar=10 and _line=5) or (nvar=11 and _line=2) then oddsci="Reference";
    if (nvar in (7) and _line in (2,3,4,5)) or
       (nvar in (10) and _line in (2,3,4,5,6,7)) or (nvar in (11) and _line in (3,4))
		then tp="-"; 

    * put '-' in missing p-value; 
	if (nvar=7 and _line=6) or (nvar=10 and _line=9) or (nvar=11 and _line=7) then oddsci="-";
	if _line>1 and missing(tp) then tp='-';

	* if odds ratio is masked, then mask tp, too; 
    if category=1 and oddsci='-' then tp='-'; 

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
proc contents data=father_part2 out=contents noprint; run;
proc sql; 
	select cats('"',LABEL,'"')into: labels separated by ' ' 
	from contents where not missing(LABEL);
quit;
%put &labels;

options noquotelenmax;
** Set title **;
%let title="Table 25. Association between potential paternal risk factors/confounders and NDD; primary outcome";

** Set footnotes **;
%let fn1="Legend: The overall column represents the number and percentage of offspring with each characteristic (percentage is calculated over the total number of offspring). The columns Event/Non-Event represent the number and percentage of events and non-events (NDD) in each subgroup defined by the characteristic (percentage is calculated over the number of offspring with the characteristic, i.e. row percentage). The association between each characteristics and the outcome is tested by fitting a logistic regression model, and the odds ratios (OR) with 95% confidence intervals (CI) and likelihood ratio (Wald) test are reported.";
%let fn2="a) all available data prior to index date (childbirth)";
%let fn4="c) 12-months lookback from LMP2";
%let fn5="d) 3-months lookback from LMP2";
%let fn6="e) at index (date of childbirth)";
%let fn7="f) at mother s LMP2";
%let fn8="g) calendar years will be grouped in each country according to the length of the study period";

title &title.;
footnote1 j=l &fn1.;
footnote2 j=l &fn2.;
footnote4 j=l &fn4.;
footnote5 j=l &fn5.;
footnote6 j=l &fn6.;
footnote7 j=l &fn7.;
footnote8 j=l &fn8.;

ods listing close;
%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table25 &track_date..xlsx" options(SHEET_NAME="Table25"  embedded_titles='yes' embedded_footnotes='yes');
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
		if level in ("Paternal risk factors/confounders") then do;
            call define (_col_,"style","style={fontweight=bold background=lightgray}");
        end;
		if level in (&labels.) then do;
            call define (_col_,"style","style={fontweight=bold}");
        end;
    endcomp;

run;
ods excel close;

**********************************************************END*****************************************************************;
