/* Below are assumption and variable inputs */
		
		proc datasets lib=work kill;quit;		
			%let Admission_Rate=1.0; /*setting to 1 as our total population is reduced to 3% of total to mimic 3% admit rate*/
			%LET HOSP_LOS=7; /*Hospital Length of Stay*/
			%let DiagnosedRate=1.0; /* 1.0 is no adjustment for missed diagnosis*/
			%let ICUPercent =.25; /*Percent of admitted needing ICU*/
			%let KnownAdmits=5.0; /*Number of total admits at start of curve*/
/* population of 21 counties arround cuyahoga= 4390484 *.03*(.23+.29)	=68473		 */
/*N (Population)*/%let Population=68490;/*this is CC market share and UH market share* 3% admit *total population of 21 county */
			%let MarketSharePercent=1.0/*Market share percent - this is no longer used so set to 1.0*/
			%let RecoveryDays=14; /* number of days remaining in Infectious column*/
			%let doublingtime=4.25; /*Not used in fitting algorithm in favor of calculated R(0) and R(t)*/
			%let DAY_ZERO='13MAR2020'd;
			%let ECMO_RATE=.03; /*Percent of admissions that require ECMO*/
			%let ECMO_LOS=6; /*Length of stay for those using ECMO*/
			%LET VENT_RATE=.22; /*percent of admitted needing ventilation*/
			%LET VENT_LOS=10; /*Length of stay for those on invasive ventilator*/
			%LET DIAL_RATE=.03; /*percent of admissions needing dialysis*/
			%LET DIAL_LOS=11; /* length of stay for those utilizing dialysis*/
			%LET ICU_LOS=17; /*length of stay for those in ICU*/
			%LET Incubationperiod=0; 
			%LET initrecovered=0;
			%LET FATALITYRATE=0;
			
			
			
			%LET HOSP_RATE = %SYSEVALF(&Admission_Rate. * &DiagnosedRate.); /*rate of those infected being admitted (Set to 1.0 for admit fitting)*/
				%LET ICU_RATE = %SYSEVALF(&ICUPercent. * &DiagnosedRate.);
				%LET VENT_RATE = %SYSEVALF(&VentPErcent. * &DiagnosedRate.);
			* calculated parameters used in models;
/* 			%let I=5; /*Static entry admission number to allow fitting  */
				%LET I = %SYSEVALF(&KnownAdmits. / 
											&MarketSharePercent. / 
											(&Admission_Rate. * &DiagnosedRate.));
/* the following is specific to CCF coding and included prior to the %EasyRun Macro */
	libname DL_RA teradata server=tdprod1 database=abc;
	libname DL_COV teradata server=tdprod1 database=def;
	
	PROC IMPORT DATAFILE="Filepath/DailyAdmissions.csv"
		DBMS=CSV
		OUT=DailyAdmits
		REPLACE;
		GETNAMES=YES;
RUN;
	data work.dailyAdmitCumSum;
	set work.dailyadmits;
	CUMULATIVE_CASE_COUNT +UHCCFMerge_Admits;
	keep CUMULATIVE_CASE_COUNT date;
	where date>'12MAR2020'd and Date<'17NOV2020'd;
	;run;
data STORE.FIT_INPUT;set work.dailyAdmitCumSum;
			TIME+1;
			run;
/*Stock SIER fit*/

			/*BELOW is new version as of 7.27.20*/
PROC MODEL DATA = STORE.FIT_INPUT OUTMODEL=SEIRMOD_I NOPRINT; 
					/* Parameters of interest */
					PARMS /*R0 3.0*/ /*I0 4*/ RI5 -1  /*DI '28MAR2020'd*/;
/* 					BOUNDS 1 <= R0 <= 5; */
/* R(0)=4.15 and the rest are R(t) values adjusted over time as model fit was determined and variable could be removed for static number */
/*RI6 is being fit with a R0 value and used to calculate future SIR curve.*/					
					RESTRICT /*R0 +*/ 4.2-2.8-0.25-.1709 +.353 + -.3385 + RI6 > 0;
/* 					restrict DI>'26MAR2020'd; */
/* 					Restrict DI2>'15MAY2020'd; */
					/* Fixed values */
					N = &Population.;
					INF = &RecoveryDays.;
/* 					SIGMA = &SIGMA.; */
					STEP = CDF('NORMAL',DATE, '30MAR2020'd, 1);
					STEP2 = CDF('NORMAL',DATE, '06APR2020'd, 1);
					STEP3 = CDF('NORMAL',DATE, '22MAY2020'd, 1);
					STEP4 = CDF('NORMAL',DATE, '24JUN2020'd, 1);
					STEP5 = CDF('NORMAL',DATE, '20JUL2020'd, 1);
					STEP6 = CDF('NORMAL',DATE, '09OCT2020'd, 1);
					
/* 					STEP6= CDF('NORMAL',DATE, DI6,1); */
					/* Differential equations */
					GAMMA = 1 / INF;
					BETA = (/*R0+*/4.2 + -2.8*STEP + -0.25*STEP2 -.1709*STEP3 + .353*STEP4 + -.3385*STEP5+ RI6*STEP6) * GAMMA / N;
					/* Differential equations */
					/* a. Decrease in healthy susceptible persons through infections: number of encounters of (S,I)*TransmissionProb*/
					DERT.S_N = -BETA * S_N * I_N;
					/* b. inflow from a. -Decrease in Exposed: alpha*e "promotion" inflow from E->I;*/
/* 					DERT.E_N = BETA * S_N * I_N - SIGMA * E_N; */
					/* c. inflow from b. - outflow through recovery or death during illness*/
					DERT.I_N =  BETA * S_N * I_N - GAMMA * I_N;
					/* d. Recovered and death humans through "promotion" inflow from c.*/
					DERT.R_N = GAMMA * I_N;
					CUMULATIVE_CASE_COUNT = I_N + R_N;
					/* Fit the data */
					FIT CUMULATIVE_CASE_COUNT INIT=(S_N=&Population. I_N=&I. R_N=0) / TIME=TIME DYNAMIC OUTPREDICT OUTACTUAL OUT=FIT_PRED LTEBOUND=1E-10 OUTEST=FIT_PARMS
						;
					OUTVARS S_N /*E_N*/ I_N R_N;
				QUIT;
				

	DATA FIT_PRED;
					SET FIT_PRED;
					LABEL CUMULATIVE_CASE_COUNT='Cumulative Incidence';
					FORMAT ModelType $30. DATE DATE9.; 
/* 					DATE = &FIRST_CASE. + TIME - 1; */
					ModelType="TMODEL - SEIR - FIT";
				run;
				DATA FIT_PARMS;
					SET FIT_PARMS;
					FORMAT ModelType $30.; 
					ModelType="TMODEL - SEIR - FIT";
					ScenarioNameUnique=cats("&Scenario.",' (',ScenarioIndex,'-',"&SYSUSERID.",'-',"&ScenarioSource.",')');
				run;

			/*Capture basline R0, date of Intervention effect, R0 after intervention*/
				/*This one below is used for the static beginning parameters*/



				
/*FINAL DATA STEPS  */
			DATA DINIT(Label="Initial Conditions of Simulation"); 
					FORMAT DATE DATE9.; 
					DO TIME = 0 TO 465; 
						S_N = &Population.;
/* 						E_N = 0; */
						I_N = &I. / 1;
						R_N = 0;
						
						DATE = '13MAR2020'd + TIME;
						OUTPUT; 
					END; 
				RUN;
	
			
			proc model data=DINIT model=SEIRMOD_I;
			solve S_N /*E_N*/ I_N R_N / OUT = TMODEL_SEIR_FIT_I;
			;quit;
/*This one will be for splitting scenarios and oscilation*/

		DATA Store.TMODEL_SEIR_FIT_I;
					FORMAT ModelType $30.DATE ADMIT_DATE DATE9. Scenarioname $30. ScenarioNameUnique $100.;
					ModelType="TMODEL - SEIR - FIT";
					LABEL HOSPITAL_OCCUPANCY="Hospital Occupancy" ICU_OCCUPANCY="ICU Occupancy" VENT_OCCUPANCY="Ventilator Utilization"
						ECMO_OCCUPANCY="ECMO Utilization" DIAL_OCCUPANCY="Dialysis Utilization";
					RETAIN LAG_S LAG_I LAG_R LAG_N CUMULATIVE_SUM_HOSP CUMULATIVE_SUM_ICU CUMULATIVE_SUM_VENT CUMULATIVE_SUM_ECMO CUMULATIVE_SUM_DIAL Cumulative_sum_fatality
						CUMULATIVE_SUM_MARKET_HOSP CUMULATIVE_SUM_MARKET_ICU CUMULATIVE_SUM_MARKET_VENT CUMULATIVE_SUM_MARKET_ECMO CUMULATIVE_SUM_MARKET_DIAL cumulative_Sum_Market_Fatality;
					LAG_S = S_N; 
					LAG_E = E_N; 
					LAG_I = I_N; 
					LAG_R = R_N; 
					LAG_N = N; 
					SET TMODEL_SEIR_FIT_I(RENAME=(TIME=DAY) DROP=_ERRORS_ _MODE_ _TYPE_);
					N = SUM(S_N, E_N, I_N, R_N);
					SCALE = LAG_N / N;
				/* START: Common Post-Processing Across each Model Type and Approach */
					NEWINFECTED=LAG&IncubationPeriod(SUM(LAG(SUM(S_N,E_N)),-1*SUM(S_N,E_N)));
					IF NEWINFECTED < 0 THEN NEWINFECTED=0;
	HOSP = NEWINFECTED * &HOSP_RATE. * &MarketSharePercent.; 
	Hosp=Hosp+0;
					ICU = NEWINFECTED * &ICU_RATE. * &MarketSharePercent. * &HOSP_RATE.;
					VENT = NEWINFECTED * &VENT_RATE. * &MarketSharePercent. * &HOSP_RATE.;
					ECMO = NEWINFECTED * &ECMO_RATE. * &MarketSharePercent. * &HOSP_RATE.;
					DIAL = NEWINFECTED * &DIAL_RATE. * &MarketSharePercent. * &HOSP_RATE.;
					Fatality = NEWINFECTED * &FatalityRate * &MarketSharePercent. * &HOSP_RATE.;
					MARKET_HOSP = NEWINFECTED * &HOSP_RATE.;
					MARKET_ICU = NEWINFECTED * &ICU_RATE. * &HOSP_RATE.;
					MARKET_VENT = NEWINFECTED * &VENT_RATE. * &HOSP_RATE.;
					MARKET_ECMO = NEWINFECTED * &ECMO_RATE. * &HOSP_RATE.;
					MARKET_DIAL = NEWINFECTED * &DIAL_RATE. * &HOSP_RATE.;
					Market_Fatality = NEWINFECTED * &FatalityRate. * &HOSP_RATE.;
					CUMULATIVE_SUM_HOSP + HOSP;
					CUMULATIVE_SUM_ICU + ICU;
					CUMULATIVE_SUM_VENT + VENT;
					CUMULATIVE_SUM_ECMO + ECMO;
					CUMULATIVE_SUM_DIAL + DIAL;
					Cumulative_sum_fatality + Fatality;
					CUMULATIVE_SUM_MARKET_HOSP + MARKET_HOSP;
					CUMULATIVE_SUM_MARKET_ICU + MARKET_ICU;
					CUMULATIVE_SUM_MARKET_VENT + MARKET_VENT;
					CUMULATIVE_SUM_MARKET_ECMO + MARKET_ECMO;
					CUMULATIVE_SUM_MARKET_DIAL + MARKET_DIAL;
					cumulative_Sum_Market_Fatality + Market_Fatality;
					CUMADMITLAGGED=ROUND(LAG&HOSP_LOS.(CUMULATIVE_SUM_HOSP),1) ;
					CUMICULAGGED=ROUND(LAG&ICU_LOS.(CUMULATIVE_SUM_ICU),1) ;
					CUMVENTLAGGED=ROUND(LAG&VENT_LOS.(CUMULATIVE_SUM_VENT),1) ;
					CUMECMOLAGGED=ROUND(LAG&ECMO_LOS.(CUMULATIVE_SUM_ECMO),1) ;
					CUMDIALLAGGED=ROUND(LAG&DIAL_LOS.(CUMULATIVE_SUM_DIAL),1) ;
					CUMMARKETADMITLAG=ROUND(LAG&HOSP_LOS.(CUMULATIVE_SUM_MARKET_HOSP));
					CUMMARKETICULAG=ROUND(LAG&ICU_LOS.(CUMULATIVE_SUM_MARKET_ICU));
					CUMMARKETVENTLAG=ROUND(LAG&VENT_LOS.(CUMULATIVE_SUM_MARKET_VENT));
					CUMMARKETECMOLAG=ROUND(LAG&ECMO_LOS.(CUMULATIVE_SUM_MARKET_ECMO));
					CUMMARKETDIALLAG=ROUND(LAG&DIAL_LOS.(CUMULATIVE_SUM_MARKET_DIAL));
					ARRAY FIXINGDOT _NUMERIC_;
					DO OVER FIXINGDOT;
						IF FIXINGDOT=. THEN FIXINGDOT=0;
					END;
					HOSPITAL_OCCUPANCY= ROUND(CUMULATIVE_SUM_HOSP-CUMADMITLAGGED,1) /*REMOVING STARTING CENSUS CONSTANT*/;
					ICU_OCCUPANCY= ROUND(CUMULATIVE_SUM_ICU-CUMICULAGGED,1);
					VENT_OCCUPANCY= ROUND(CUMULATIVE_SUM_VENT-CUMVENTLAGGED,1);
					ECMO_OCCUPANCY= ROUND(CUMULATIVE_SUM_ECMO-CUMECMOLAGGED,1);
					DIAL_OCCUPANCY= ROUND(CUMULATIVE_SUM_DIAL-CUMDIALLAGGED,1);
					Deceased_Today = Fatality;
					Total_Deaths = Cumulative_sum_fatality;
					MedSurgOccupancy=Hospital_Occupancy-ICU_Occupancy;
					MARKET_HOSPITAL_OCCUPANCY= ROUND(CUMULATIVE_SUM_MARKET_HOSP-CUMMARKETADMITLAG,1);
					MARKET_ICU_OCCUPANCY= ROUND(CUMULATIVE_SUM_MARKET_ICU-CUMMARKETICULAG,1);
					MARKET_VENT_OCCUPANCY= ROUND(CUMULATIVE_SUM_MARKET_VENT-CUMMARKETVENTLAG,1);
					MARKET_ECMO_OCCUPANCY= ROUND(CUMULATIVE_SUM_MARKET_ECMO-CUMMARKETECMOLAG,1);
					MARKET_DIAL_OCCUPANCY= ROUND(CUMULATIVE_SUM_MARKET_DIAL-CUMMARKETDIALLAG,1);	
					Market_Deceased_Today = Market_Fatality;
					Market_Total_Deaths = cumulative_Sum_Market_Fatality;
					Market_MEdSurg_Occupancy=Market_Hospital_Occupancy-MArket_ICU_Occupancy;
					DATE = &DAY_ZERO. + round(DAY,1);
					ADMIT_DATE = SUM(DATE, &IncubationPeriod.);
					LABEL
						ADMIT_DATE = "Date of Admission"
						DATE = "Date of Infection"
						DAY = "Day of Pandemic"
						HOSP = "New Hospitalized Patients"
						HOSPITAL_OCCUPANCY = "Current Hospitalized Census"
						MARKET_HOSP = "New Region Hospitalized Patients"
						MARKET_HOSPITAL_OCCUPANCY = "Current Region Hospitalized Census"
						ICU = "New Hospital ICU Patients"
						ICU_OCCUPANCY = "Current Hospital ICU Census"
						MARKET_ICU = "New Region ICU Patients"
						MARKET_ICU_OCCUPANCY = "Current Region ICU Census"
						MedSurgOccupancy = "Current Hospital Medical and Surgical Census (non-ICU)"
						Market_MedSurg_Occupancy = "Current Region Medical and Surgical Census (non-ICU)"
						VENT = "New Hospital Ventilator Patients"
						VENT_OCCUPANCY = "Current Hospital Ventilator Patients"
						MARKET_VENT = "New Region Ventilator Patients"
						MARKET_VENT_OCCUPANCY = "Current Region Ventilator Patients"
						DIAL = "New Hospital Dialysis Patients"
						DIAL_OCCUPANCY = "Current Hospital Dialysis Patients"
						MARKET_DIAL = "New Region Dialysis Patients"
						MARKET_DIAL_OCCUPANCY = "Current Region Dialysis Patients"
						ECMO = "New Hospital ECMO Patients"
						ECMO_OCCUPANCY = "Current Hospital ECMO Patients"
						MARKET_ECMO = "New Region ECMO Patients"
						MARKET_ECMO_OCCUPANCY = "Current Region ECMO Patients"
						Deceased_Today = "New Hospital Mortality: Fatality=Deceased_Today"
						Fatality = "New Hospital Mortality: Fatality=Deceased_Today"
						Total_Deaths = "Cumulative Hospital Mortality"
						Market_Deceased_Today = "New Region Mortality"
						Market_Fatality = "New Region Mortality"
						Market_Total_Deaths = "Cumulative Region Mortality"
						N = "Region Population"
						S_N = "Current Susceptible Population"
						E_N = "Current Exposed Population"
						I_N = "Current Infected Population"
						R_N = "Current Recovered Population"
						NEWINFECTED = "New Infected Population"
						ModelType = "Model Type Used to Generate Scenario"
						SCALE = "Ratio of Previous Day Population to Current Day Population"
						ScenarioIndex = "Scenario ID: Order"
						ScenarioSource = "Scenario ID: Source (BATCH or UI)"
						ScenarioUser = "Scenario ID: User who created Scenario"
						ScenarioNameUnique = "Unique Scenario Name"
						Scenarioname = "Scenario Name"
						;
				/* END: Common Post-Processing Across each Model Type and Approach */
					DROP LAG:  ;
					/*Hosp=Hosp+5.53;*/ /*STATIC ADDITION */
				RUN;