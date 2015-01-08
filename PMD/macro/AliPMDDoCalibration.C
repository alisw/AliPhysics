//
// This macro does the calibration by finding the isolated cells
//
void AliPMDDoCalibration()
{
  AliPMDCalibrator *calibrator = new AliPMDCalibrator();
  calibrator->Init();
  calibrator->CalculateIsoCell();
  calibrator->Store();
}
