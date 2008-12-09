//____________________________________________________________________
//
// $Id$
//
// Make fake alignment data.
//
/** @file    MakeCalibration.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:58:27 2006
    @brief   Make fake calibration data.
*/
/** Make fake calibration data 
    @ingroup simple_script
 */
Float_t
AdcPerMip2Gain(Int_t adc) 
{
  return 2.9;// / adc * AliFMDParameters::Instance()->GetEdepMip();
}

void
MakeCalibration(const char* base="local://$ALICE_ROOT")
{
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage(base);

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libFMDanalysis.so");
  gSystem->Load("libFMDutil.so");
  AliFMDCalibFaker f(AliFMDCalibFaker::kAll, 0);
  f.SetRunRange(0,999999);
  f.SetGainSeed(AdcPerMip2Gain(60)); // From astrid test beam 
  f.SetThresholdFactor(3);
  f.SetPedestalRange(80,130); // From ASTRID test-beam
  f.SetDeadChance(0);
  f.SetZeroThreshold(0);
  f.SetStripRange(0, 127);
  f.SetRate(4);
  f.Exec();
}
//____________________________________________________________________
//
// EOF
//
