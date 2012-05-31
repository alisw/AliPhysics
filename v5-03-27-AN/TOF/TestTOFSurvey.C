void TestTOFSurvey() 
{
  /////////////////////////////////////////////////////////////////////////
  //
  // Simple Test macro for Conversion of Survey Data into AlignObj for TOF
  // Author: S.Arcelli
  // root[1] AliTOFtestProb("match-6KevPYTHIA-0.2T.root")
  /////////////////////////////////////////////////////////////////////////
 AliTOFAlignment *al = new AliTOFAlignment();
 Float_t mis[6]={0.,0.,5.,10.,0.,0.}; //+5 cm shift z and 10 deg rot in phi
 al->BuildGeomForSurvey();
 al->InsertMisAlignment(mis);
 al->AlignFromSurvey();
}
