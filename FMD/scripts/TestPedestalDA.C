//____________________________________________________________________
//
// $Id: TestRawIO.C 13249 2006-03-24 16:09:36Z cholm $
//
// Test of AliFMDPedestalDA
//
/** @ingroup simple_script
 */
#ifndef __CINT__
# include <TSystem.h>
# include <TStopwatch.h>
# include <AliCDBManager.h>
# include <AliRawReader.h>
# include <AliFMDPedestalDA.h>
# include <AliFMDParameters.h>
#endif
void TestPedestalDA(Char_t* fileName="data.raw", Int_t runNumber=1, 
		    Bool_t  oldFormat=kTRUE, Bool_t diagnostics=kFALSE)
{
  
  //This script runs the pedestal DA using the class AliFMDPedestalDA
#ifdef __CINT__
  // Load utility library
  gSystem->Load("libFMDutil");
#endif

  // Set-up CDB interface 
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(runNumber);
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Set debug level
  AliLog::SetModuleDebugLevel("FMD", 1);

  // Set-up paramters 
  AliFMDParameters* params = AliFMDParameters::Instance();
  params->Init();
  params->UseCompleteHeader(oldFormat);
  
  // Set-up raw readers 
  AliRawReader *reader = AliRawReader::Create(fileName);


  // Set-up timer 
  TStopwatch timer;
  timer.Start();

  Bool_t append = false;
  // Make and run DA
  AliFMDPedestalDA pedestalDA;
  pedestalDA.SetSaveDiagnostics(diagnostics);
  pedestalDA.Run(reader, append);

  // Stop and print summary
  timer.Stop();
  timer.Print();
}
//
// EOF
//
