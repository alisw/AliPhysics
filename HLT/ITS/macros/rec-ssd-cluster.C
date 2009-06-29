// $Id$
//aliroot -b -q rec-spd-cluster.C | tee rec-spd-cluster.log

void rec_ssd_cluster(const char* input="./", char* opt="")
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
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTITS.so libAliHLTSample.so loglevel=0x7c chains=";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  
  int minddl=512;
  int maxddl=527;
  int spec=0x1;
  int ddlno=0;

  TString dummyInput="";
  for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
    TString arg, publisher, cf;
    
    arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISSD ' -dataspec 0x%08x -verbose",ddlno, spec);
    publisher.Form("DP_%d", ddlno);
    AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
    
    cf.Form("CF_%d",ddlno);
    AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSSD", publisher.Data(), "");

    if (dummyInput.Length()>0) dummyInput+=" ";
    dummyInput+=cf;

    spec=spec<<1;
  }

  //add dummy
  //AliHLTConfiguration dummyconf("dummy", "Dummy", dummyInput.Data(), "-output_percentage 0");
  //option+="dummy";

  AliHLTConfiguration cfconf("clusterHisto","ITSClusterHisto",dummyInput.Data(),"");
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



