//__________________________________________________//
AliEbyEMultiplicityFluctuationAnalysis *GetAnalysisMFObject() {
 
   gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/macro/ConfigureEbyEBaseObject.C");
   AliEbyEEventSelector *base = 0;
   base =  GetEbyEAnalysisBaseObject();
   
   AliEbyEMultiplicityFluctuationAnalysis *analysis = new AliEbyEMultiplicityFluctuationAnalysis();
   
   analysis->SetBaseAnalysis(base);
   return analysis;
}

