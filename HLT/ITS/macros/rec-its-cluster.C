// $Id$
/*
 * Example macro to run the Offline ITS Cluster Finding in HLT.
 * The reconstruction is done from the ITS raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-its-cluster.C | tee rec-its-cluster.log
 *   aliroot -b -q rec-its-cluster.C'("./","spd")' | tee rec-its-cluster.log
 * </pre>
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-its-cluster.C'("input.root")'
 * </pre>
 *
 * In the second parameter you can change the clusterfinders to run only
 * SPD, SDD, SSD or All:
 *    - SPD runs only Silicon Pixels
 *    - SDD runs only Silicon Drift
 *    - SSD runs only Silicon Stips
 *    - All will run the full ITS
 *
 * All is default. The chain ends in a histogram component. 
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_its
 * @author st05886@ift.uib.no
 */

void rec_its_cluster(const char* input="./", char* opt="All")
{

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }
  
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
  int verbosity=1; // 1 higher verbosity, 0 reduced output
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Setting up which output to give
  //
  TString option;
  option+="libAliHLTUtil.so libAliHLTRCU.so libAliHLTITS.so libAliHLTSample.so ";
  option+="loglevel=";
  option+=verbosity>0?"0x7c":"0x79";
  option+=" chains=";
  Bool_t runspd=kFALSE, runsdd=kFALSE, runssd=kFALSE;
  TString allArgs=opt;
  TString argument;
  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries(); i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
 
      if (argument.CompareTo("spd", TString::kIgnoreCase)==0) {
	runspd=kTRUE;
	continue;
      } 
      if (argument.CompareTo("sdd", TString::kIgnoreCase)==0) {
	runsdd=kTRUE;
	continue;
      }
      if (argument.CompareTo("ssd",TString::kIgnoreCase)==0) {
	runssd=kTRUE;
	continue;
      }
      if (argument.CompareTo("all",TString::kIgnoreCase)==0) {
	runspd=kTRUE;
	runsdd=kTRUE;
	runssd=kTRUE;
	continue;
      }
      else {
	break;
      }
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
    
  //The spec starts from 0x1 in SPD, SDD and SSD. So 0x1 is ddl 0 for SPD, 0x10 is ddl 1, and so on
  //in SDD 0x1 is ddl 256, 0x10 is ddl 257, and so on. This means that the spec has to be set to 0x1 
  //before the loops over the clusterfinder

  int minddl=0;          //min ddl number for SPD
  int maxddl=19;         //max ddl number for SPD
  int spec=0x1;          //spec for ddl's
  int ddlno=0;
  TString dummyInput="";

  if(runspd){
    for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
      TString arg, publisher, cf;
      
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISPD ' -dataspec 0x%08x %s",ddlno, spec, (verbosity>0?" -verbose":""));
      publisher.Form("DP_%d", ddlno);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      cf.Form("CF_%d",ddlno);
      AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSPD", publisher.Data(), "");
      
      if (dummyInput.Length()>0) dummyInput+=" ";
      dummyInput+=cf;
      
      spec=spec<<1;
    }
  }
  
  if(runsdd){
    minddl=256;        //min ddl number for SDD    
    maxddl=279;        //max ddl number for SDD
    spec=0x1;          //spec for ddl's
    
    for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
      TString arg, publisher, cf;
      
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISDD ' -dataspec 0x%08x %s",ddlno, spec, (verbosity>0?" -verbose":"")); 
      publisher.Form("DP_%d", ddlno);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      cf.Form("CF_%d",ddlno);
      AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSDD", publisher.Data(), "");
      
      if (dummyInput.Length()>0) dummyInput+=" ";
      dummyInput+=cf;
      
      spec=spec<<1;
    }
  }

  if(runssd){
    minddl=512;      //min ddl number for SSD     
    maxddl=527;      //max ddl number for SSD
    spec=0x1;        //spec for ddl's
    
    for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
      TString arg, publisher, cf;
      
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISSD ' -dataspec 0x%08x %s",ddlno, spec, (verbosity>0?" -verbose":""));
      publisher.Form("DP_%d", ddlno);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      cf.Form("CF_%d",ddlno);
      AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSSD", publisher.Data(), "");
      
      if (dummyInput.Length()>0) dummyInput+=" ";
      dummyInput+=cf;
      
      spec=spec<<1;
    }
  }

  AliHLTConfiguration cfconf("clusterHisto","ITSClusterHisto",dummyInput.Data(),"-pushback-period=10");
  AliHLTConfiguration fwconf("histFile","ROOTFileWriter", "clusterHisto","-datafile ClusterHisto -concatenate-events -overwrite");

  option+="histFile";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
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
  rec.SetOption("HLT", option);
  rec.Run();
}



