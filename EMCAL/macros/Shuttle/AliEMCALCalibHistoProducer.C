///
/// \file  AliEMCALCalibHistoProducer.C
/// \ingroup EMCAL_Shuttle
/// \brief Make Preprocessor OCDB configuration
///
/// Obsolete??
/// Script is to be run at DAQ computers (LDC, GDC or HLT);
/// it fills the histograms with amplitudes per channel.
/// These histograms will be processed in the Shuttle preprocessor
/// e.g., to extract mean amplitudes needed for relative channel calibration
///
/// This example assumes that the input data is supplied from the
/// raw data file in the ROOT format.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS), 2006
///

///
/// Main method
///
void AliEMCALCalibHistoProducer(const char* file="EMCAL.raw.root") 
{
  AliRawReaderRoot* rf = new AliRawReaderRoot(file);
 
  AliEMCALCalibHistoProducer hp(rf);

//   hp.SetSMInstalled(10,kFALSE); //Supermodule not installed
//   hp.SetSMInstalled(11,kFALSE);
 
  hp.Run();
}
