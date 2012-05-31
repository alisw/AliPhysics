// $Id$
/*
 * Example macro to run the ITS tracker with the TPC reconstruction.
 * The reconstruction is done from the TPC and ITS raw data.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-tpc-its.C | tee rec-tpc-its.log
 *   aliroot -b -q rec-tpc-its.C'("./","spd")' | tee rec-tpc-its.log
 * </pre>
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-tpc-its.C'("input.root")'
 * </pre>
 * 
 * In the second parameter you can change the clusterfinders to run only
 * SPD, SDD, SSD or All:
 *    - SPD runs only Silicon Pixels
 *    - SDD runs only Silicon Drift
 *    - SSD runs only Silicon Stips
 *    - All will run the full ITS. This is default
 *
 * In the first section, an analysis chain is defined. The scale of the
 * chain can be defined by choosing the range of sectors and partitions.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Gaute.Ovrebekk@ift.uib.no
 */
void rec_tpc_its(const char* input="./", char* opt="All")
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
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Setting up which output to give
  //
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTITS.so libAliHLTGlobal.so loglevel=0x7c chains=";
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
	cerr << "Unknown argument" << endl;
	break;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  int minddl=0;          //min ddl number for SPD
  int maxddl=19;         //max ddl number for SPD
  int spec=0x1;          //spec for ddl's
  int ddlno=0;
  TString cfout="";

  if(runspd){
    for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
      TString arg, publisher, cf;
      
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISPD ' -dataspec 0x%08x -verbose",ddlno, spec);
      publisher.Form("DP_%d", ddlno);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      cf.Form("CF_%d",ddlno);
      AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSPD", publisher.Data(), "");
      
      if (cfout.Length()>0) cfout+=" ";
      cfout+=cf;
      
      spec=spec<<1;
    }
  }
  
  if(runsdd){
    minddl=256;        //min ddl number for SDD    
    maxddl=279;        //max ddl number for SDD
    spec=0x1;          //spec for ddl's
    
    for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
      TString arg, publisher, cf;
      
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISDD ' -dataspec 0x%08x -verbose",ddlno, spec); 
      publisher.Form("DP_%d", ddlno);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      cf.Form("CF_%d",ddlno);
      AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSDD", publisher.Data(), "");
      
      if (cfout.Length()>0) cfout+=" ";
      cfout+=cf;
      
      spec=spec<<1;
    }
  }

  if(runssd){
    minddl=512;      //min ddl number for SSD     
    maxddl=527;      //max ddl number for SSD
    spec=0x1;        //spec for ddl's
    
    for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
      TString arg, publisher, cf;
      
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISSD ' -dataspec 0x%08x -verbose",ddlno, spec);
      publisher.Form("DP_%d", ddlno);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      cf.Form("CF_%d",ddlno);
      AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSSD", publisher.Data(), "");
      
      if (cfout.Length()>0) cfout+=" ";
      cfout+=cf;
      
      spec=spec<<1;
    }
  }

  TString ITSinput = "TPC-globalmerger ";
  ITSinput += cfout;
  AliHLTConfiguration itstrackerconf("itstracker","ITSTracker",ITSinput.Data(),""); 
  //option+="itstracker";
  
  AliHLTConfiguration globalConverter("globalConverter", "GlobalEsdConverter"   , "TPC-globalmerger itstracker", "");
  AliHLTConfiguration sink("esdfile", "EsdCollector"   , "globalConverter", "-directory hlt-tpc-esd");
  option+="esdfile";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd())); 
  rec.SetOption("HLT", option);
  rec.Run();
}
