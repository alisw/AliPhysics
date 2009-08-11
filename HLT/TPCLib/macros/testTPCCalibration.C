// $Id$
/*
 * Example macro to run the HLT tracker embedded into AliRoot reconstruction. 
 * The reconstruction is done from the TPC raw data. The seed maker and the 
 * calibration components have been added to the standard HLT TPC reconstruction.
 *
 * Its output is contained in TPC-globalmerger and all the TPC clusters necessary
 * to run the calibration are in TPC-clusters.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q testTPCCalibration.C'("./")' | tee testTPCCalibration.log
 *   aliroot -b -q testTPCCalibration.C'("./","calibtime")' | tee testTPCCalibration.log
 *   aliroot -b -q testTPCCalibration.C'("./","calibtimegain")' | tee testTPCCalibration.log
 * </pre>
 *
 * The macro assumes raw data to be available in the rawx folders, either
 * simulated or real data. If this is the only input specified, the macro
 * will load the calibration analysis task. If the second argument is specified,
 * then the individual calibration classes can be tested.
 *
 * <pre>
 *   aliroot -b -q testTPCCalibration.C'("raw.root","calibtime")'
 * </pre>
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Kalliopi.Kanaki@ift.uib.no
 */

void testTPCCalibration(const char* input="./", const char* option="task"){

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
 
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
  
  TString seedMakerInput;
  if(seedMakerInput.Length()>0) seedMakerInput+=" ";
  seedMakerInput+="TPC-globalmerger";
  seedMakerInput+=" ";
  seedMakerInput+="TPC-clusters";

  AliHLTConfiguration seedconf("seeds", "TPCCalibSeedMaker", seedMakerInput.Data(), "");

  TString calibInput;
  if(calibInput.Length()>0) calibInput+=" ";
  calibInput+="TPC-esd-converter";
  calibInput+=" ";
  calibInput+="seeds";
    
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
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");

  // NOTE: FillESD is a step in the AliReconstruction sequence and has
  // nothing to do with the fact that this macro writes ESD output
  // HLT processes the HLTOUT during FillESD and extracts data which
  // has already been prepared. This step is currently not necessary for
  // this macro
  rec.SetFillESD("");
  if     (calibOption.CompareTo("task")==0)          rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=TPCcalib"); 
  if     (calibOption.CompareTo("calibtime")==0)     rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=calibTime"); 
  else if(calibOption.CompareTo("calibtimegain")==0) rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=calibTimeGain");

  rec.Run();
}
