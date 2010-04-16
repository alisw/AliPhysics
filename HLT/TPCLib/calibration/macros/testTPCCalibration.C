// $Id$
/*
 * Example macro to run the HLT tracker embedded into AliRoot reconstruction. 
 * The reconstruction is done from the TPC raw data. The calibration components 
 * have been added to the standard HLT TPC reconstruction.
 *
 * They take input from TPC-clusters, TPC-globalmerger and GLOBAL-esd-converter. 
 *
 * Usage:
 * <pre>
 * 
 * Usage: aliroot -b -q -l testTPCCalibration.C'("file", "cdb", minEvent, maxEvent, option)'
 *
 * Examples:
 *     testTPCCalibration.C'("raw://run12345", minEvent, MaxEvent,"calibtime")'
 *     testTPCCalibration.C'("raw.root", "local://$PWD", minEvent, MaxEvent, "calibtime")'
 * 
 * Defaults
 *     cdb="raw://"  -> takes OCDB from GRID
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *     option="task" -> loads the analysis task that calls both the drift velocity and the gain calibration
 * 
 * </pre>
 *
 * To see usage examples run:
 * 
 * aliroot -q testTPCCalibration.C
 *
 * The argument option should be 'calibtime' for now, as the rest of the options 
 * call components that are not working yet.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 * 
 * @ingroup alihlt_tpc
 * @author Kalliopi.Kanaki@ift.uib.no
 */

void testTPCCalibration(const char* filename, const char *cdbURI, int minEvent=-1, int maxEvent=-1, const char* option="task"){

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");  
  
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "Please delete file galice.root or run at a different place." << endl;
    return;
  }

  if (!filename) {
    cerr << "Please specify input file name." << endl;
    return;
  }

  /////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
 
  AliHLTSystem *gHLT = AliHLTPluginBase::GetInstance();
  
  TString calibInput;
  if(calibInput.Length()>0) calibInput+=" ";
  calibInput+="GLOBAL-esd-converter TPC-clusters TPC-globalmerger";
    
  TString calibOption = option;
  if     (calibOption.CompareTo("task")==0)          AliHLTConfiguration calibtimeconf("TPCcalib",      "TPCCalibration",   calibInput.Data(), "");
  else if(calibOption.CompareTo("calibtime")==0)     AliHLTConfiguration calibtimeconf("calibTime",     "TPCCalibTime",     calibInput.Data(), "");
  else if(calibOption.CompareTo("calibtimegain")==0) AliHLTConfiguration calibtimeconf("calibTimeGain", "TPCCalibTimeGain", calibInput.Data(), "");
  else 
    {
      cerr << "\nPlease specify an option for the calibration component you want to run, e.g. aliroot -b -q testTPCCalibration.C'(\"raw.root\",\"calibtime\")' \n"<< endl;
      return;
    }


  /////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off
 
  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri  = cdbURI;
  TString strfile = filename;
  if(struri.BeginsWith("raw://") ||  strfile.Contains("://") && !strfile.Contains("local://")){
     TGrid::Connect("alien");
  }
 
  // Set the CDB storage location
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  
  AliReconstruction rec;
  rec.SetInput(filename);
  rec.SetRunVertexFinder(kFALSE);
  //rec.SetWriteESDfriend(kTRUE);
 
  if(minEvent>=0 || maxEvent>minEvent){
     if(minEvent<0) minEvent=0;
     if(maxEvent<minEvent) maxEvent=minEvent;
     rec.SetEventRange(minEvent,maxEvent);
  } 
  
  if     (calibOption.CompareTo("task")==0)          rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=TPCcalib"); 
  if     (calibOption.CompareTo("calibtime")==0)     rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=calibTime"); 
  else if(calibOption.CompareTo("calibtimegain")==0) rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=calibTimeGain");
  
  rec.SetRunPlaneEff(kFALSE);    
  // switch off cleanESD
  rec.SetCleanESD(kFALSE);
 
  rec.SetRunReconstruction("HLT");
  rec.SetRunQA(":");

  AliLog::Flush();
  rec.Run();
}

void testTPCCalibration(const char *filename, int minEvent=-1, int maxEvent=-1, const char* option="task"){
 
  testTPCCalibration(filename, "raw://", minEvent, maxEvent, option);
}

void testTPCCalibration(){

  cout << "testTPCCalibration: Run AliRoot reconstruction and TPC drift velocity calibration locally" << endl;
  cout << " Usage: aliroot -b -q -l testTPCCalibration.C'(\"file\", \"cdb\", minEvent, maxEvent, option)'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     testTPCCalibration.C'(\"raw://run12345\", minEvent, MaxEvent,\"calibtime\")'" << endl;
  cout << "     testTPCCalibration.C'(\"raw.root\", \"local://$PWD\", minEvent, MaxEvent, \"calibtime\")'" << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     cdb=\"raw://\"  -> takes OCDB from GRID" << endl;
  cout << "     minEvent=-1   -> no lower event selection" << endl;
  cout << "     maxEvent=-1   -> no upper event selection" << endl;
  cout << "     option=\"task\" -> loads the analysis task that calls both the drift velocity and the gain calibration" << endl;
}





// void testTPCCalibration(const char* input="./", const char* option="task"){
// 
//   gSystem->Load("libANALYSIS");
//   gSystem->Load("libTPCcalib");  
//   
//   if(!gSystem->AccessPathName("galice.root")){
//     cerr << "Please delete file galice.root or run at a different place." << endl;
//     return;
//   }
// 
//   if (!input) {
//     cerr << "Please specify input or run without arguments." << endl;
//     return;
//   }
// 
//   /////////////////////////////////////////////////////////////////////////
//   //
//   // init the HLT system in order to define the analysis chain below
//   //
//  
//   AliHLTSystem *gHLT=AliHLTPluginBase::GetInstance();
//   
//   TString calibInput;
//   if(calibInput.Length()>0) calibInput+=" ";
//   calibInput+="GLOBAL-esd-converter TPC-globalmerger TPC-clusters";
//     
//   TString calibOption = option;
//   if     (calibOption.CompareTo("task")==0)          AliHLTConfiguration calibtimeconf("TPCcalib",      "TPCCalibration",   calibInput.Data(), "");
//   else if(calibOption.CompareTo("calibtime")==0)     AliHLTConfiguration calibtimeconf("calibTime",     "TPCCalibTime",     calibInput.Data(), "");
//   else if(calibOption.CompareTo("calibtimegain")==0) AliHLTConfiguration calibtimeconf("calibTimeGain", "TPCCalibTimeGain", calibInput.Data(), "");
//   else 
//     {
//       cerr << "\nPlease specify an option for the calibration component you want to run, e.g. aliroot -b -q testTPCCalibration.C'(\"raw.root\",\"calibtime\")' \n"<< endl;
//       return;
//     }
// 
// 
//   /////////////////////////////////////////////////////////////////////////
//   //
//   // Init and run the reconstruction
//   // All but HLT reconstruction is switched off
//   //
//   AliReconstruction rec;
//   rec.SetInput(input);
//   rec.SetRunVertexFinder(kFALSE);
//   rec.SetRunLocalReconstruction("HLT");
//   rec.SetRunTracking("");
//   rec.SetLoadAlignFromCDB(0);
//   rec.SetRunQA(":");
//   rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//   rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
// 
// 
//   // NOTE: FillESD is a step in the AliReconstruction sequence and has
//   // nothing to do with the fact that this macro writes ESD output
//   // HLT processes the HLTOUT during FillESD and extracts data which
//   // has already been prepared. This step is currently not necessary for
//   // this macro
//   rec.SetFillESD("");
//   if     (calibOption.CompareTo("task")==0)          rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=TPCcalib"); 
//   if     (calibOption.CompareTo("calibtime")==0)     rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=calibTime"); 
//   else if(calibOption.CompareTo("calibtimegain")==0) rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTTPC.so libAliHLTGlobal.so libAliHLTTPCCalibration.so loglevel=0x7c chains=calibTimeGain");
// 
//   rec.Run();
// }
