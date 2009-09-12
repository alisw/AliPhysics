// $Id: $
/*
 * Example macro to run TPC calibration interface embedded into AliRoot reconstruction. 
 * The input used is unpacked from the HLTOUT of recorded data, in which the merged
 * tracks, the ESD and the TPC clusters have been added. The seed maker and the calibration 
 * components follow.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q calib-from-hltout.C'("raw.root","task"     )'     | tee calib-from-hltout.log
 *   aliroot -b -q calib-from-hltout.C'("raw.root","calibtime")'     | tee calib-from-hltout.log
 *   aliroot -b -q calib-from-hltout.C'("raw.root","calibtimegain")' | tee calib-from-hltout.log
 * </pre>
 *
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way. Keep in mind that you should use the proper default and
 * specific storage otherwise the reconstruction will break.
 *
 * @ingroup alihlt_tpc
 * @author Kalliopi.Kanaki@ift.uib.no
 */

void calib_from_hltout(const char* input="./", const char* option="task"){

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");  
  
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "Please delete file galice.root or run at a different place." << endl;
    return;
  }

  if (!input) {
    cerr << "Please specify input or run without arguments." << endl;
    return;
  }

  /////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
 
  AliHLTSystem *gHLT = AliHLTPluginBase::GetInstance();
  if(!gHLT){
    cerr << "fatal error: cannot get HLT instance" << endl;
  }

  AliHLTConfiguration publisher("hltout", "AliHLTOUTPublisher", NULL, "");
  
  AliHLTConfiguration seedconf("seeds", "TPCCalibSeedMaker", "hltout","");
  
  TString calibOption = option;
  if     (calibOption.CompareTo("task")==0)          AliHLTConfiguration calibtimeconf("TPCcalib",      "TPCCalibration",  "hltout seeds" ,"");
  else if(calibOption.CompareTo("calibtime")==0)     AliHLTConfiguration calibtimeconf("calibTime",     "TPCCalibTime",    "hltout seeds", "");
  else if(calibOption.CompareTo("calibtimegain")==0) AliHLTConfiguration calibtimeconf("calibTimeGain", "TPCCalibTimeGain","hltout seeds", "");
  else 
    {
      cerr << "\nPlease specify an option for the calibration component you want to run." << endl;
      return;
    }


  /////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  //rec.SetDefaultStorage("local:///home/kanaki/calibComp/libraryTest/OCDB/");
  //rec.SetSpecificStorage("GRP/GRP/Data", "local:///home/kanaki/calibComp/libraryTest/OCDB/");


  // NOTE: FillESD is a step in the AliReconstruction sequence and has
  // nothing to do with the fact that this macro writes ESD output
  // HLT processes the HLTOUT during FillESD and extracts data which
  // has already been prepared. This step is currently not necessary for
  // this macro
  rec.SetFillESD("");
  
  if     (calibOption.CompareTo("task")==0)          rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTGlobal.so libAliHLTTPC.so libAliHLTTPCCalibration.so loglevel=0x7c chains=TPCcalib"); 
  if     (calibOption.CompareTo("calibtime")==0)     rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTGlobal.so libAliHLTTPC.so libAliHLTTPCCalibration.so loglevel=0x7c chains=calibTime"); 
  else if(calibOption.CompareTo("calibtimegain")==0) rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTGlobal.so libAliHLTTPC.so libAliHLTTPCCalibration.so loglevel=0x7c chains=calibTimeGain");
  
  //rec.SetEventRange(0,20);
  rec.Run();
}
