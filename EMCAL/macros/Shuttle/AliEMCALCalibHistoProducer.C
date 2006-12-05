void AliEMCALCalibHistoProducer(const char* file="EMCAL.raw.root") 
{
  // Script is to be run at DAQ computers (LDC, GDC or HLT);
  // it fills the histograms with amplitudes per channel.
  // These histograms will be processed in the Shuttle preprocessor
  // e.g., to extract mean amplitudes needed for relative channel calibration
  //
  // This example assumes that the input data is supplied from the
  // raw data file in the ROOT format.
  //
  //  Gustavo Conesa Balbastre, December 2006 
  // Load EMCAL shuttle library as it is not linked to aliroot
  gSystem->Load("libEMCALshuttle");

  AliRawReaderRoot* rf = new AliRawReaderRoot(file);
  AliEMCALCalibHistoProducer hp(rf);
//   hp.SetSMInstalled(10,kFALSE); //Supermodule not installed
//   hp.SetSMInstalled(11,kFALSE);
  hp.Run();
}
