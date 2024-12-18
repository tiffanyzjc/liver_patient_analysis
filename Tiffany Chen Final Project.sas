/* The following options are described in HomeworkOptionsAnnotated.sas in compass.
   Please refer to that file to determine which settings you wish to use or modify
   for your report.
*/
*ods html close;
*options nodate nonumber leftmargin=1in rightmargin=1in;
*title;
*ods escapechar="~";
*ods graphics on / width=4in height=3in;
*ods rtf file='FileNameWithFullPath.rtf'
        nogtitle startpage=no;
*ods noproctitle;

/* The raw data in Indian Liver Patient Dataset (ILPD).csv is from

   https://archive.ics.uci.edu/ml/datasets/ILPD+(Indian+Liver+Patient+Dataset) 

   Dua, D. and Graff, C. (2019). 
   UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. 
   Irvine, CA: University of California, School of Information and Computer Science. 
*/
proc import datafile="/home/u63746690/datasets/Indian Liver Patient Dataset (ILPD).csv"
	out= liver
	dbms = csv
	replace;
	getnames=no;
run;
/* after importing, rename the variables to match the data description */
data liver;
	set liver;
	Age=VAR1; Gender=VAR2; TB=VAR3;	DB=VAR4; Alkphos=VAR5;
	Alamine=VAR6; Aspartate=VAR7; TP=VAR8; ALB=VAR9; AGRatio=VAR10;
	LiverPatient=VAR11;
	*if VAR11=1 then LiverPatient='Yes';
		*Else LiverPatient='No';
	drop VAR1--VAR11;
run;

proc sort data=liver;
	by LiverPatient;
run;

*Part 1: Descriptive Analysis;
* - univariate results by LiverPatient status;
* -- LiverPatient = Yes: age mean = 46.15, tb median = 1.4, db median = 0.5, alkphos median = 229, alamine median = 41, aspartate median = 52.5, tp mean = 6.46, alb mean = 3.06, agratio mean = 0.9;
* -- LiverPatient = No: age mean = 41.24, tb median = 0.8, db median = 0.2, alkphos median = 186, alamine median = 27, aspartate median = 29, tp mean = 6.54, alb mean = 3.34, agratio mean = 1.03;
* -- Except for TP with LiverPatient = No, tests for normality of all other variables indicate non-normal distributions;
proc univariate data=liver normaltest;
  var Age TB--AGRatio;
  histogram Age TB--AGRatio /normal;
  probplot Age TB--AGRatio;
  by LiverPatient; 
  ods select Moments BasicMeasures Histogram ProbPlot TestsForNormality;
run;
* - wilcoxon rank sum test to determine whether there are significant differences in measurements for LiverPatients = Yes/No;
* -- significantly different values in LiverPatients = Yes/No for all measurements except TP;
proc npar1way data=liver wilcoxon;
  class LiverPatient;
  var Age TB--AGRatio;
  ods exclude KruskalWallisTest;
run;
* - correlation between LiverPatient status and measurements;
* -- similarly, there is a significant correlation between LiverPatient status and all measurements except TP;
proc corr data=liver pearson spearman;
  var Age TB--AGRatio LiverPatient;
  ods select PearsonCorr SpearmanCorr;
run;
* - scatterplot visualizations for correlation between all variables with each other, grouped by LiverPatient status;
proc sgscatter data=liver;
  title "Scatterplot";
  matrix Age TB--AGRatio / group=LiverPatient;
run;
* - scatterplot visualizations for variables with strongest correlation to LiverPatient status;
proc sgplot data=liver;
  scatter y=LiverPatient x=Aspartate / group=LiverPatient;
run;
proc sgplot data=liver;
  scatter y=LiverPatient x=TB / group=LiverPatient;
run;
proc sgplot data=liver;
  scatter y=LiverPatient x=DB / group=LiverPatient;
run;
proc sgplot data=liver;
  scatter y=LiverPatient x=Alamine / group=LiverPatient;
run;

*Part 2: Modeling LiverPatient with clinical measurements, gender, and age;
* - logistic regression;
* -- best model: backward selection, LiverPatient = Age + DB + Alamine + TP + ALB;
proc logistic data=liver;
	class Gender/param=ref;
	model LiverPatient = Age--AGRatio/
		selection=backward;
	ods select OddsRatios ParameterEstimates 
		GlobalTests ModelInfo FitStatistics;
run;

proc logistic data=liver;
	class Gender/param=ref;
	model LiverPatient = Age--AGRatio/
		selection=forward;
	ods select OddsRatios ParameterEstimates 
		GlobalTests ModelInfo FitStatistics;
run;

proc logistic data=liver;
	class Gender/param=ref;
	model LiverPatient = Age--AGRatio/
		selection=stepwise;
	ods select OddsRatios ParameterEstimates 
		GlobalTests ModelInfo FitStatistics;
run;
* -- diagnostics for final model;
proc logistic data=liver plots=effect;
	class Gender/param=ref;
	model LiverPatient = Age--AGRatio/
		selection=backward;
    output predprobs=individual out=liverout;
run;
* -- compare levels observed with levels predicted;
proc freq data=liverout;
    tables LiverPatient*_into_/nopercent norow nocol;
run;

*Part 3: Modeling Total Proteins with other clinical measurements, gender, and age for those who are liver patients;
* - dataset including only observations for liver patients;
data patients;
  set liver;
  if LiverPatient=1;
run;
* - gamma model (response yi > 0);
proc genmod data=patients;
	class Gender;
	model TP = Age Gender TB DB Alkphos Alamine Aspartate ALB AGRatio / dist=gamma 
		link=log type1 type3;	
	output out=gammares pred=presp_n stdreschi=presids
		stdresdev= dresids;
	ods select ModelInfo ModelFit ParameterEstimates Type1 Type3;	
run;
* -- the Type 3 analysis suggests removing Alkphos from the model (largest p-value);
proc genmod data=patients;
	class Gender;
	model TP = Age Gender TB DB Alamine Aspartate ALB AGRatio / dist=gamma 
		link=log type1 type3;	
	output out=gammares pred=presp_n stdreschi=presids
		stdresdev= dresids;
	ods select ModelInfo ModelFit ParameterEstimates Type1 Type3;	
run;
* -- the Type 3 analysis suggests removing Gender from the model (largest p-value);
proc genmod data=patients;
	class Gender;
	model TP = Age TB DB Alamine Aspartate ALB AGRatio / dist=gamma 
		link=log type1 type3;	
	output out=gammares pred=presp_n stdreschi=presids
		stdresdev= dresids;
	ods select ModelInfo ModelFit ParameterEstimates Type1 Type3;	
run;
* -- the Type 3 analysis suggests removing Age from the model (largest p-value);
proc genmod data=patients;
	class Gender;
	model TP = TB DB Alamine Aspartate ALB AGRatio / dist=gamma 
		link=log type1 type3;	
	output out=gammares pred=presp_n stdreschi=presids
		stdresdev= dresids;
	ods select ModelInfo ModelFit ParameterEstimates Type1 Type3;	
run;
* -- the Type 3 analysis suggests removing TB from the model (largest p-value);
proc genmod data=patients;
	class Gender;
	model TP = DB Alamine Aspartate ALB AGRatio / dist=gamma 
		link=log type1 type3;	
	output out=gammares pred=presp_n stdreschi=presids
		stdresdev= dresids;
	ods select ModelInfo ModelFit ParameterEstimates Type1 Type3;	
run;
* -- final model: TP = DB + Alamine + Aspartate + ALB + AGRatio;
* -- plot diagnostic residuals: overall residuals are relatively flat with the exception of a few outliers. the model is relatively good;
proc genmod data=patients plots=(stdreschi stdresdev);
	model TP = DB Alamine Aspartate ALB AGRatio / dist=gamma 
		link=log type1 type3;
	ods select ModelInfo DiagnosticPlot;	
run;

*Part 4: Grouping observations based on clinical measurements and age, comparing the groupings with LiverPatient status
* - remove outlier observations;
proc univariate data=liver;
  var Age TB--AGRatio;
  ods select ExtremeObs;
run;
* -- remove obs 560, 582, 560, 355, 377, 356, 397, 287, 541, 84, 83, 137, 378, 325, 368 which had outlying values for >= 2 variables;
data liver2;
  set liver;
  id = _N_;
  if id not in(560, 582, 355, 377, 356, 397, 287, 541, 84, 83, 137, 378, 325, 368);
run;
* - cluster;
proc cluster data=liver2 method=average outtree=liveravg ccc pseudo print=15 plots=all;
	var Age: TB--AGRatio;
   	copy LiverPatient;
	ods select ClusterHistory Dendrogram CccPsfAndPsTSqPlot;
run;
* -- 6 clusters is a relatively good choice based on the pseudo t^2 and pseudo F statistics;
proc tree data=liveravg noprint ncl=6 out=clusters;
   copy Age: TB--AGRatio LiverPatient:;
run;
* -- most observations for both liver patients and non liver patients were grouped into the same cluster (cluster 1), but some liver patients were clustered into other groups based on clinical measurements and age;
proc freq data=clusters;
  tables cluster*LiverPatient/ nopercent norow nocol;
run;
* - compare means between clusters;
* -- cluster 1 vs clusters 2 and 3: liver patients and non liver patients generally shared similar values for clinical measurements and age. However, some liver patients tended to have higher TB, DB, Alamine, and Aspartate, along with lower ALB and AGRatio;
proc sort data=clusters;
 by cluster;
run;
proc means data=clusters;
 var Age TB--AGRatio;
 by cluster;
run;
* - PCA to extract the most prominent features in the clusters;
* -- based on the scree plot, three features is sufficient to accurately characterize the clusters, and will describe 67.5% of the variance in the clustering;
* -- princomp1: the first most prominent feature determining the clustering is a comparison between ALB vs TB and DB measurements;
* -- princomp2: comparison between Age vs Alamine, TP, and ALB;
* -- princomp3: comparison between Alamine and Aspartate vs TB and DB;
proc princomp data=clusters n=6 out=pcout;
  var Age TB--AGRatio;
run
proc sgplot data=pcout;
  scatter y=prin1 x=prin2 / markerchar=LiverPatient;
run;
proc sgplot data=pcout;
  scatter y=prin1 x=prin3 / markerchar=LiverPatient;
run;

*Part 5: Classifying four groups based on gender and LiverPatient status (F/Yes, F/No, M/Yes, M/No)
* - combine LiverPatient and Gender variables into a new patient variable;
data livergen;
  set liver;
  if LiverPatient = 1 then Liver='Y';
           			  else Liver='N';
  if Gender = 'Female' then Gender='F';
           			   else Gender='M';
  patient=Liver||Gender;
  drop LiverPatient;
run;
* - test for equal covariance using chi-square test;
* -- F statistic is significant at the 0.1 level, so we assume unequal covariance and proceed with QDA;
proc sort data=livergen;
	by patient;
run;
proc discrim data=livergen method=normal pool=test
       testout=plotpQ testoutd=plotdQ
	   short noclassify crosslisterr;
   class patient;
   var Age TB--AGRatio;
run;
* - QDA;
* -- find discriminant and classification count errors;
* -- the manova test statistics are significant at the 0.05 level (p <0.0001), indicating that the age and clinical measurement variables may be able to discriminate between liver patient statuses/genders;
* -- the classification table indicates that about 20.41% of non liver patient females, 43.97% of non liver patient males, 20.88% of liver patient females, and 43.05% of liver patient males were correctly classified;
* -- the overall classification error rate is 67.93%, meaning that 67.93% of all observations are estimated to be misclassified into liver patient status/gender groups by the model;
proc discrim data=livergen crossvalidate manova;
   class patient;
   var Age TB--AGRatio;
   ods select ChiSq MultStat ClassifiedCrossVal ErrorCrossVal;
run;
* - stepwise selection to find the most significant variables;
* -- variables selected: DB, Alkphos, Age, Alamine;
proc stepdisc data=livergen sle=.05 sls=.05;
   	class patient;
   	var Age TB--AGRatio;
	ods select Summary;
run;
* - QDA with the selected variables;
* -- the manova test statistics are significant at the 0.05 level (p <0.0001), indicating that the age and clinical measurement variables may be able to discriminate between liver patient statuses/genders;
* -- the classification table indicates that about 80% of non liver patient females, 29.91% of non liver patient males, 10.87% of liver patient females, and 35.80% of liver patient males were correctly classified;
* -- the overall classification error rate is 65.52%, meaning that 65.52% of all observations are estimated to be misclassified into liver patient status/gender groups by the model;
* -- the reduced model performs slightly better overall, as it has a lower overall error rate.;
* -- liver patient males, non liver patient males, and liver patient females tended to be overclassified into the non liver patient female group. This may indicate that all four groups share similarities in Age, DB, Alkphos, and Alamine measurements.;
* -- male liver patients were nearly just as likely to be classified as non liver patient males as they were to be classified as non liver patient females. This may indicate similarities in Age, DB, Alkphos, and Alamine measurements across all males in general - regardless of their liver patient status.;
proc discrim data=livergen pool=test crossvalidate;
  	class patient;
  	var Age DB Alkphos Alamine;
   	priors proportional;
   	ods select Chisq ClassifiedCrossVal ErrorCrossVal;
run;



























