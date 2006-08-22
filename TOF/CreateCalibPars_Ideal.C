void CreateCalibPars_Ideal(){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliTOFGeometry *geom = new AliTOFGeometryV5(); 
  AliTOFcalib *tofcalib = new AliTOFcalib(geom);
  AliTOFCal *tofCal= new AliTOFCal(geom);
  tofCal->CreateArray();//"empty" channels as a default for ideal (par,delay=0)
  TH1F *hToT= new TH1F(); //"empty" ToT histo as a default for ideal 
  tofcalib->WriteSimParOnCDB("TOF/Calib",0,0,tofCal,hToT);
}


