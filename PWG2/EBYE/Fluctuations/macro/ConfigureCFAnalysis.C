//__________________________________________________//
AliEbyEChargeFluctuationAnalysis *GetAnalysisCFObject() {
 
   gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/macro/ConfigureEbyEBaseObject.C");
   AliEbyEEventSelector *base = 0;
   base =  GetEbyEAnalysisBaseObject();
   
   AliEbyEChargeFluctuationAnalysis *analysis = new AliEbyEChargeFluctuationAnalysis();
   
   analysis->SetBaseAnalysis(base);
   return analysis;
}

