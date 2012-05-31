void MakeVZEROTimeSlewingEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the time slewing OCDB object
  // x = ADC charge / TDC threshold
  TF1 *slew = new TF1("TimeSlewing","[0]*TMath::Power(x,[1])",1,1024);
  slew->SetParameter(0,1.05e1);
  slew->SetParameter(1,-5.0e-1);

  /*
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.04508e+01   1.25434e-02   1.62336e-04  -1.01124e-02
   2  p1          -4.85664e-01   1.26751e-03   5.74823e-06  -2.65093e-01
   3  p2           7.79978e+01   7.05932e-03   4.05090e-05  -3.55601e-02
 FCN=7919.32 FROM MIGRAD    STATUS=CONVERGED     104 CALLS         105 TOTAL
                     EDM=4.59683e-11    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.03716e+01   8.01903e-03   8.58839e-05   2.80998e-03
   2  p1          -4.90768e-01   1.00097e-03   3.29876e-06   9.17350e-02
   3  p2           7.80244e+01   5.57426e-03   3.72049e-05   1.45924e-02
  */

  /*
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.06394e+01   8.78531e-03  -0.00000e+00   2.64914e-06
   2  p1          -5.00000e-01     fixed    
   3  p2           7.80541e+01   2.10785e-03   0.00000e+00   3.36855e-07

  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           1.04238e+01   5.26216e-03   8.98026e-05   1.59005e-05
   2  p1          -5.00000e-01     fixed    
   3  p2           7.80744e+01   1.37139e-03   3.72288e-05   2.78012e-05
  */
	
  TObjString str("VZERO Time-slewing correction");

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time-slewing correction used in reconstruction and MC simulation");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/TimeSlewing",0,AliCDBRunRange::Infinity());

  storLoc->Put(slew, id, md);

  storLoc->Delete();
  delete md;

}
