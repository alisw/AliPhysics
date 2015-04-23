void MakeFMDCDBObjects()
{
  AliFMDCalibFaker calibFaker(AliFMDCalibFaker::kAll, "local://$ALICE_SOURCE/OCDB");
  calibFaker.SetRunRange(0, AliCDBRunRange::Infinity());
  calibFaker.Exec(); //FMD/Calib/AltroMap Dead Pedestal PulseGain RecoParam SampleRate StripRange ZeroSuppression
}
//FMD/Align/Data
