void AliPHOSCalibHistoProducer(const char* file="2006run2211.root") 
{
  // Script is to be run at DAQ computers (LDC, GDC, HLT or MOOD);
  // it fills the histograms with amplitudes per channel.
  // These histograms will be processed in the Shuttle preprocessor
  // e.g., to extract mean amplitudes needed for relative channel calibration
  //
  // This example assumes that the input data is supplied from the
  // raw data file in the ROOT format.
  //
  // Author: Boris Polichtchouk, 4 October 2006

  // Load PHOS shuttle library as it is not linked to aliroot
  gSystem->Load("libPHOSshuttle");

  AliRawReaderRoot* rf = new AliRawReaderRoot(file);
  AliPHOSCalibHistoProducer hp(rf);
  hp.Run();
}
