// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/**
 * @file   transformFCFtoSFC.C
 * @author Kalliopi.Kanaki@ift.uib.no
 * @date
 * @brief  Test macro for the AliHLTTPCHWClusterTransformComponent.cxx
 *
 *
 * This macro tests the component that does the tranformation from the HW
 * cluster format to the software cluster structure.
 * For the macro to work, the user needs to use the FilePublisher component
 * to make the input files available. The data specification has to be given
 * and the input files have to be present.
 *
 * In order for the reconstruction to work properly at the absence of raw/ folders,
 * a directory has to be specified as the input argument of the function
 * and an empty raw0 folder will be created if not present.
 *
 * The macro looks for a folder called FCFFiles by default, where the produced FCF files
 * should be (format TPC_ddlnumber.fcf or .bin). Argument number 2 can be set to look in a different
 * folder if wanted.
 *
 * In addition, since the $ALICE_ROOT/OCDB/GRP/GRP/Data entry has been removed as obsolete,
 * the user needs to produce 1 simulated event, in order to create a proper GRP entry in the
 * local folder, which will be used by the line
 *
 * rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
 *
 * The component runs without arguments, except for the case when we want to differentiate
 * between the FCF and the SFC output. Then the argument -change-dataId will change the
 * data Id of the FCF output.
 *
 */
void transformFCFtoSCF(const char* input="./",const char* dirName="./"){
  
  gSystem->Exec("rm galice.root");  
  
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "Please delete file galice.root or run at a different place." << endl;
    return;
  }

  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  
  TString inputStr;
  inputStr.Form("%s",input);
  //  if(inputStr.)

  if(inputStr.CompareTo("./") == 0){
    if(gSystem->AccessPathName("./raw0")){
      cout<<"No raw folder found, making dummy raw0 directory...."<<endl;
      gSystem->Exec("mkdir raw0");
    }
    inputStr.Form("%s/",gSystem->pwd());
  }

  if(!inputStr.EndsWith("/")){
    inputStr+="/";
  }

  TString dir;
  dir.Form("%s%s/raw0",inputStr.Data(),dirName);

  if(!dir.EndsWith("/")){
    dir+="/";
  }

  if(gSystem->AccessPathName(dir)){
    cerr << "Input directory does not exist: "<<dir.Data() << endl;
    return;
  }
  

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // 
  // init the HLT system in order to define the analysis chain below
  //
  gSystem->Load("libHLTrec.so");
  AliHLTSystem *gHLT = AliHLTReconstructorBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
    
  int iMinSlice =  0; 
  int iMaxSlice = 35;
  int iMinPart  =  0;
  int iMaxPart  =  5;

  TString dumpOutput;
  TString FCFInput;
  TString mergerInput;
  TString allclusters;

  for(int slice=iMinSlice; slice<=iMaxSlice; slice++){ 
      
      TString trackerInput;
      for(int part=iMinPart; part<=iMaxPart; part++){

          int ddlno=768;
          if (part>1) ddlno+=72+4*slice+(part-2);
          else ddlno+=2*slice+part;

          TString file;
          file.Form("raw0/TPC_%d.bin",ddlno);
          
	  if(gSystem->AccessPathName(file)){
      	    cerr << "Input file does not exist: "<< file.Data() << endl;
      	    continue;
          }
	  	  
          TString publisher, fcf, dumpout;
          publisher.Form("DP_%02d_%d", slice, part);
          TString argument;
          argument.Form("-datatype 'HWCLUST1' 'TPC '  -datafile %s -dataspec 0x%02x%02x%02x%02x", file.Data(), slice, slice, part, part);
        	
          AliHLTConfiguration pubconf(publisher.Data(), "FilePublisher", NULL, argument.Data());
          fcf.Form("FCF_%02d_%d", slice, part);
                        
          AliHLTConfiguration hwconf(fcf.Data(), "TPCHWClusterTransform", publisher.Data(), "");
          
	  if(trackerInput.Length()>0) trackerInput+=" ";
          trackerInput+=fcf;
	  if(allclusters.Length()>0) allclusters+=" ";
          allclusters+=fcf;
	  
          dumpout.Form("DUMP_%02d_%d", slice, part);
          
	  TString argDump, dFile;
          dFile.Form("FCF_%d.dump",ddlno);
          argDump.Form("-directory FCFClusterDump -subdir=raw -datafile %s -specfmt= -blcknofmt= -idfmt= -skip-datatype", dFile.Data());
          
	  AliHLTConfiguration clusDumpconf(dumpout.Data(), "TPCClusterDump", fcf.Data(), argDump.Data());               
         
          if(dumpOutput.Length()>0) dumpOutput+=",";
          dumpOutput+=dumpout;
      }
      
      TString tracker;
      tracker.Form("TR_%02d", slice);
      AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "");
      
      if(mergerInput.Length()>0) mergerInput+=" ";
      mergerInput+=tracker;
  }
  
  AliHLTConfiguration mergerconf("globalmerger","TPCCAGlobalMerger",mergerInput.Data(),"");
  AliHLTConfiguration esdconf("ESD","GlobalEsdConverter","globalmerger","");
  AliHLTConfiguration sink("esdfile", "EsdCollector", "ESD", "-directory FCF-hlt-tpc-esd"); 

  TString histoInput; 
  if(histoInput.Length()>0) histoInput+=" ";
  histoInput+=allclusters;
  histoInput+=" ";
  histoInput+="globalmerger";
  
  AliHLTConfiguration histconf("histo","TPCTrackHisto",histoInput.Data(),"");
  
  AliHLTConfiguration cfcompconf("comparison","TPCCFComparison",allclusters.Data(),"");
  
  AliHLTConfiguration rfwconf("RFW","ROOTFileWriter","histo","-datafile FCF_trackhisto -overwrite -concatenate-events");

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off

  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");  
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  TString option;
  option.Form("libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTGlobal.so loglevel=0x7c chains=comparison,RFW,esdfile,%s",dumpOutput.Data());
  rec.SetOption("HLT", option);
  rec.Run();


}
