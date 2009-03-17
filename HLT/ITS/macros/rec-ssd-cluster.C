
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
  gSystem->Load("libHLTrec.so");
  AliHLTSystem* gHLT=AliHLTReconstructorBase::GetInstance();
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Setting up which output to give
  //
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTITS.so libAliHLTSample.so loglevel=0x7c chains=";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  
  int minddl=0x00000000;
  int maxddl=0x00004000;

  int ddl=0;
  int ddlno=0;

  TString dummyInput="";
  for(int ddl=minddl;ddl<=maxddl;){
    TString arg, publisher, cf;
    //arg.Form("-minid %d -datatype 'DDL_RAW ' 'ITS '  -dataspec 0x%02x%02x%02x%02x -verbose", ddl, 00, 00, 00, 00);   
    //arg.Form("-detector ITSSPD -skipempty -datatype 'DDL_RAW ' 'ITS ' -verbose");
    //arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);
    //arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISPD ' -dataspec 0x%08x -verbose",ddlno, ddl);
    //arg.Form("-detector ITSSPD -datatype 'DDL_RAW ' 'ISPD ' -skipempty -dataspec 0x%08x -verbose",ddl);
    arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISSD ' -dataspec 0x%08x -verbose",ddlno, ddl);
    publisher.Form("DP_%d", ddl);
    AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
    
    cf.Form("CF_%d",ddl);
    AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSSD", publisher.Data(), "");

    if (dummyInput.Length()>0) dummyInput+=" ";
    dummyInput+=cf;

    ddlno++;
    if(ddl==0x0000000){ddl++;}else{ddl = ddl << 1;}
  }

  //add dummy
  AliHLTConfiguration dummyconf("dummy", "Dummy", dummyInput.Data(), "-output_percentage 0");

  option+="dummy";

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



