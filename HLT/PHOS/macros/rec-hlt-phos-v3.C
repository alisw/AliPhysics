
#include "AliHLTPHOSConstants.h"

//void rec_hlt_phos(const char* input="./")//, char* opt="decoder ESD")
void rec_hlt_phos()//, char* opt="decoder ESD")
{
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }
 
//   if (!input) {
//     cerr << "please specify input or run without arguments" << endl;
//     return;
//   }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
   
  int moduleStart = 2;
  int moduleEnd = 4;
  int rcuStart = 0;
  int rcuEnd = 3;
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTPHOS.so loglevel=0x7f chains=ESD_WRITER";
  //TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTPHOS.so loglevel=0x7f chains=PHS-EM";
  TString emInput;
  
  for (int module = moduleStart; module <= moduleEnd; module++) 
    {
      TString clInput;
      
      for(int rcu = rcuStart; rcu <= rcuEnd; rcu++) 
	{
	  TString arg, publisher, ra, dm;
	  // raw data publisher components
	  publisher.Form("PHS-RP_%02d_%d", module, rcu);
	  arg.Form("-minid %d -datatype 'DDL_RAW ' 'PHOS'  -dataspec 0x%x ", 1792 + module*(PhosHLTConst::NRCUSPERMODULE) + rcu, 0x1 << (module*PhosHLTConst::NRCUSPERMODULE + rcu));
	  cout << arg << endl;
	  AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
	  
	  // Raw analyzer
	  arg = "";
	  ra.Form("PHS-RA_%02d_%d", module, rcu);
	  AliHLTConfiguration rawConf(ra.Data(), "PhosRawCrudev3", publisher.Data(), arg.Data());
	  
	  // digit maker components
	  dm.Form("PHS-DM_%02d_%d", module, rcu);
	  arg="";
	  arg.Form("-sethighgainfactor 0.005 -setlowgainfactor 0.08 -setdigitthresholds 0.005 0.002");
	  cout << arg << endl;
	  AliHLTConfiguration dmConf(dm.Data(), "PhosDigitMaker", ra.Data(), arg.Data());

	  if(clInput.Length() > 0) clInput += " ";
	  clInput+=dm;
	}
      cout << clInput << endl;
      TString arg, cl, ca;

      cl.Form("PHS-CL_%02d", module);
      arg = "";
      arg.Form("-digitthreshold 0.005 -recpointthreshold 0.1 -modulemode");
      cout << arg << endl;
      AliHLTConfiguration clConf(cl.Data(), "PhosClusterizer", clInput.Data(), arg.Data());
	
      ca.Form("PHS-CA_%02d", module);
      arg = "";
      AliHLTConfiguration caConf(ca.Data(), "PhosClusterAnalyser", cl.Data(), arg.Data());

      if(emInput.Length() > 0) emInput += " ";
      emInput += ca;
    }
      
  TString arg, em;
  
  // tracker finder components
  em.Form("PHS-EM");
  arg = "";

  AliHLTConfiguration emConf(em.Data(), "PhosEsdEntriesMaker", emInput.Data(), " ");



  //  AliHLTConfiguration writerConf("ESD_WRITER", "PhosEsdCaloClusterWriter", em.Data(), "-filename test.root -writemodulo 5");

  //  AliHLTConfiguration esdWriter("ESD_WRITER", "ROOTFileWriter"   , "PHS-EM", "PhosEsdEntriesMaker", "-datafile ESDClusters -concatenate-events -overwrite");
  AliHLTConfiguration esdWriter("ESD_WRITER", "ROOTFileWriter"   , em.Data(), "-datafile ESDClusters -concatenate-events -overwrite");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off 
  //
  AliReconstruction rec;
  //  rec.SetInput(input);
  rec.SetInput("./");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT PHOS");
  rec.SetRunTracking(":");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");

  // NOTE: FillESD is a step in the AliReconstruction sequence and has
  // nothing to do with the fact that this macro writes ESD output
  // HLT processes the HLTOUT during FillESD and extracts data which
  // has already been prepared. This step is currently not necessary for
  // this macro
  //  rec.SetRunLocalReconstruction("PHOS") ;

  rec.SetFillESD("PHOS");
    rec.SetOption("HLT", option);
  //  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTPHOS.so loglevel=0x7f chains=ESD_WRITER" )
  rec.Run();
}
