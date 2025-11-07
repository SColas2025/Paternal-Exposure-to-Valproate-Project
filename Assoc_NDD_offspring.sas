***************************************************************************************************************************;
*
* SAS code for publication
* Article: Paternal Valproate Exposure and Neurodevelopmental Disorder and Congenital Malformation Risk in Offspring 
*
* This code summarises univariate association between all pre-specified offspring risk factors/confounders, 
* and the primary outcome (binary variable indicating if the offspring was diagnosed with NDD including ASD during follow-up)
*
* By: IQVIA
*
* Overview: 
*           Read in dataset
*			Association between offspring risk factors/confounders and NDD	
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
	set dta.summary_offspring;
run;

***************************************************************************************************************************;
*Step 2: Association between offspring risk factors/confounders and NDD	
***************************************************************************************************************************;
* Create input dataset for association; 
data offspring_part2; 
     set in_data;

	*select cohort;
    if C1_COMP=1; 

	label  
           gender='Gender ^{super a}'
	       O_CONCMV_1YR='Congenital CMV ^{super b}'
           O_CONRUBE_1YR='Congenital rubella ^{super b}'
           O_FAS_EXIT='Foetal alcohol syndrome ^{super b}'
           O_FRAGILEX_EXIT='Fragile X syndrome ^{super b}'
           O_LEJEUNE_EXIT='Lejeune/cri du chat syndrome ^{super b}'
           O_TUBSIS_EXIT='Tuberous sclerosis ^{super b}';

run;  

* Create empty oddsratio;
data OR_empty;
    length oddsci $40.; 
	oddsci='-'; 
run;

* Create empty wald test;
data WT_empty;
    length test $16. tp $40.;
    test="Wald"; 
	tp='-'; 
run;

* Perform logistic regression for estimating the associations;
* NOTE: for some variables, the following warning may occur and can be ignored: 
 WARNING: Output 'GlobalTests' was not created.  Make sure that
         the output object name, label, or path is spelled
         correctly.  Also, verify that the appropriate
         procedure options are used to produce the requested
         output object.  For example, verify that the NOPRINT
         option is not used.; 
%let devar=gender O_CONCMV_1YR O_CONRUBE_1YR O_FAS_EXIT O_FRAGILEX_EXIT O_LEJEUNE_EXIT O_TUBSIS_EXIT;/*Offspring risk factors/confounders included*/
%let ds=offspring_part2;
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

%if &no_event=0 or &no_nonevent=0 %then %do;
	data OddsRatios&i.;
		    length Effect $200.; 
			set OR_empty;
            Effect="&var.";
			end; 
	run; 
%end;
%else %if &i=1 %then %do;
	* Perform logistic regression for gender;
	proc logistic data=&ds. ;
		class &var.(ref=last)/param=ref; *set male as reference; 
		model O_NDDALL(event='1')=&var./ridging=ABSOLUTE ;
		ods output OddsRatios=OddsRatios&i. GlobalTests= GlobalTests&i.;
	run;

%end;
%else %do;
	* Perform logistic regression for other variables;
	proc logistic data=&ds. ;
		class &var.(ref=first)/param=ref; *set lowest level as reference;
		model O_NDDALL(event='1')=&var./ridging=ABSOLUTE ;
		ods output OddsRatios=OddsRatios&i. GlobalTests= GlobalTests&i.;
	run;
%end;

proc sql noprint; select count(*) into:N_odds from OddsRatios&i.; quit;  

	* if oddsratios has 0 observations, use the empty oddsratio; 
	%if &N_odds=0 %then %do;
		data OddsRatios&i.;
			set OR_empty;
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
data OR_all(keep=oddsci level nvar _line);
	length Effect $ 40. oddsci $ 40. level $ 150;
	set Oddsratios1-Oddsratios7;
	if not missing(OddsRatioEst) and UpperCL<999 then do; * the limit to the upper limit: 999, upper limit more than 999 is not reasonable; 
		oddsci=cats(put(OddsRatioEst,10.2),' (',put(LowerCL,10.2),',',put(UpperCL,10.2),')');
    end; 

	infile datalines delimiter='+';
	input level nvar _line;
	datalines;	
	Female + 1 + 3 
	Congenital CMV ^{super b} + 2 + 3 
	Congenital rubella ^{super b} + 3 + 3 
    Foetal alcohol syndrome ^{super b} + 4 + 3 
    Fragile X syndrome ^{super b} + 5 + 3 
    Lejeune/cri du chat syndrome ^{super b} + 6 + 3 
    Tuberous sclerosis ^{super b} + 7 + 3 
	
run;

* Merge all Waldtests from all risk factors/confounders;
data WT_all(keep=tp level nvar _line);
	length tp $ 40 level $ 150;
	set Globaltests1-Globaltests7;
	if test="Wald";
	if not missing(ChiSq) then do; 
		tp=cats(put(ChiSq,10.2),' (',put(ProbChiSq,pvalue9.4),')');
	end; 
	
	infile datalines delimiter='+';
	input level nvar _line;
	datalines;	
	Wald test + 1 + 5 
	Congenital CMV ^{super b} + 2 + 3 
	Congenital rubella ^{super b} + 3 + 3 
    Foetal alcohol syndrome ^{super b} + 4 + 3 
    Fragile X syndrome ^{super b} + 5 + 3 + 3
    Lejeune/cri du chat syndrome ^{super b} + 6 + 3 
    Tuberous sclerosis ^{super b} + 7 + 3 
	;
run;

* Merge all results;
proc sort data=OR_all; by nvar level;
proc sort data=WT_all; by nvar level;
run;
data _To_report;
	merge OR_all(in=a) WT_all(in=b);
	by nvar level;
	if a or b;
	rename oddsci=oddsratios tp=waldtest;
	proc sort; by nvar _line; 
run;

proc sql;
	create table To_report_all as 
	select a.level, a.nvar, a._line, a.k1, a.pct1, a.k2, a.pct2, a.k3, a.pct3, 
			b.oddsratios,b.waldtest
	from sum_data a
	left join _To_report b
	on a.nvar=b.nvar and a._line=b._line
	order by a.nvar, a._line;
quit;

data To_report;
	set To_report_all;

	if (nvar=1 and _line=2) then oddsratios='Reference';
		else if (nvar=1 and _line=4 and missing(oddsratios)) then oddsratios='-';
		else if (nvar=1 and _line=5) then oddsratios='-';

	if nvar=1 and _line in (2,3,4) then waldtest="-";

	if level="No" then do;
		if strip(pct1) in ('100.00%', '0.00%') then oddsratios='-'; else oddsratios='Reference'; 
		waldtest="-";
	end;
 
run;

***************************************************************************************************************************;
*Step 3: Table creation		
***************************************************************************************************************************;
* Get bold labels; 
proc contents data=offspring_part2 out=contents noprint; run;
proc sql; 
	select cats('"',LABEL,'"')into: labels separated by ' ' 
	from contents where not missing(LABEL);
quit;
%put &labels;

** Set title **;
%let title="Table 23. Association between potential offspring risk factors/confounders and NDD; primary outcome";

** Set footnotes **;
options noquotelenmax;
%let fn1="Legend: The overall column represents the number and percentage of offspring with each characteristic (percentage is calculated over the total number of offspring). The columns Event/Non-Event represent the number and percentage of events and non-events (NDD) in each subgroup defined by the characteristic (percentage is calculated over the number of offspring with the characteristic, i.e. row percentage). The association between each characteristics and the outcome is tested by fitting a logistic regression model, and the odds ratios (OR) with 95% confidence intervals (CI) and likelihood ratio (Wald) test are reported.";
%let fn2="a) at index (childbirth)";
%let fn3="b) between index and exit date";

title &title.;
footnote1 j=l &fn1.;
footnote2 j=l &fn2.;
footnote3 j=l &fn3.;

ods listing close;
%let track_date=%sysfunc(today(),date9.);
ods excel file = "&output.\01 Tables\Table23 &track_date..xlsx" options(SHEET_NAME="Table23"  embedded_titles='yes' embedded_footnotes='yes');
ods escapechar="^";
proc report data=to_report  split = '*';
	column (
             ("NDD" (" " level)) 
             ("Overall" ("N" k1) ("%" pct1)) 
             ("Event" ("N" k2) ("%" pct2))
             ("Non-event" ("N" k3) ("%" pct3)) 
             ("Association" ("OR (95% CI)" oddsratios) ("Test statistics (p-value)" waldtest))  
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
	
	define oddsratios		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};

	define waldtest		/ display "" 
					  style(column)={vjust=middle just=center width= 3 cm}
                      style(column)={just=c cellwidth = 3 cm};
	
	compute level;
		if level in ( "Offspring risk factors/confounders") then do;
            call define (_col_,"style","style={fontweight=bold background=lightgray}");
        end;
		if level in (&labels.) then do;
            call define (_col_,"style","style={fontweight=bold}");
        end;
    endcomp;

run;
ods excel close;


**********************************************************END*****************************************************************;
