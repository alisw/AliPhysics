void addInput( TString &input, const TString &add )
{
  if (input.Length()>0) input+=" ";
  input+=add;
}
 
void rec_hlt_tpc_its(const char* input="./")
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
 
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTGlobal.so loglevel=0x7c chains=";

  TString histoInput;
  TString esdInput;
  TString Bz = "";
  //Bz = "-solenoidBz -5";

  // TPC 
  {
    TString tpcMergerInput;

    int iMinSlice=0;
    int iMaxSlice=35;
    int iMinPart=0;
    int iMaxPart=5;

    for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
      TString trackerInput;
      for (int part=iMinPart; part<=iMaxPart; part++) {
	TString arg, publisher, cf;
	TString clusterHistoOutput;
	// raw data publisher components
	int ddlno=768;
	if (part>1) ddlno+=72+4*slice+(part-2);
	else ddlno+=2*slice+part;
	arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);
	publisher.Form("DP_%02d_%d", slice, part);
	new AliHLTConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
	
	cf.Form("CF_%02d_%d", slice, part);
	new AliHLTConfiguration(cf.Data(), "TPCClusterFinder32Bit", publisher.Data(), "-release-memory");
	
	addInput( trackerInput, cf);
      }
      TString tracker;
      tracker.Form("TR_%02d", slice);
      new AliHLTConfiguration(tracker.Data(), "TPCCATracker", trackerInput.Data(), Bz.Data());
      
      addInput( tpcMergerInput, tracker);
    }

    new AliHLTConfiguration("tpcMerger","TPCCAGlobalMerger",tpcMergerInput.Data(), Bz.Data());
    addInput( esdInput, "tpcMerger");
  }

  // ITS 
  {
    Bool_t runspd=1, runsdd=1, runssd=1;
    TString cfout="";

    if(runspd){
      int minddl=0;          //min ddl number for SPD
      int maxddl=19;         //max ddl number for SPD
      int spec=0x1;          //spec for ddl's
      for(int ddlno=minddl;ddlno<=maxddl;ddlno++){  
	TString arg, publisher, cf;	
	publisher.Form("DP_%d", ddlno);
	cf.Form("CF_%d",ddlno);
	arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISPD ' -dataspec 0x%08x -verbose",ddlno, spec);
	new AliHLTConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());	
	new AliHLTConfiguration (cf.Data(), "ITSClusterFinderSPD", publisher.Data(), "");	
	addInput( cfout, cf);
	spec=spec<<1;
      }       
      AliHLTConfiguration itsVtx("itsVtx", "ITSVertexerSPD"   , cfout.Data(), Bz.Data());
      addInput( esdInput, "itsVtx");
      addInput( histoInput, "itsVtx");
    }
  
    if(runsdd){
      int minddl=256;        //min ddl number for SDD    
      int maxddl=279;        //max ddl number for SDD
      int spec=0x1;          //spec for ddl's    
      for(int ddlno=minddl;ddlno<=maxddl;ddlno++){  
	TString arg, publisher, cf;      
	publisher.Form("DP_%d", ddlno);
	cf.Form("CF_%d",ddlno);
	arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISDD ' -dataspec 0x%08x -verbose",ddlno, spec); 
	AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());	
	AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSDD", publisher.Data(), "");	
	addInput( cfout, cf);
	spec=spec<<1;
      }
    }

    if(runssd){
      minddl=512;      //min ddl number for SSD     
      maxddl=527;      //max ddl number for SSD
      spec=0x1;        //spec for ddl's      
      for(int ddlno=minddl;ddlno<=maxddl;ddlno++){  
	TString arg, publisher, cf;      
	publisher.Form("DP_%d", ddlno);
	cf.Form("CF_%d",ddlno);
	arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISSD ' -dataspec 0x%08x -verbose",ddlno, spec);
	AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());      
	AliHLTConfiguration cfconf(cf.Data(), "ITSClusterFinderSSD", publisher.Data(), "");      
	addInput( cfout, cf);      
	spec=spec<<1;
      }
    }

    addInput( cfout, "tpcMerger");      
    AliHLTConfiguration its("itsTracker", "ITSTracker"   , cfout.Data(), Bz.Data());
    addInput( esdInput, "itsTracker");
  }

  
  new AliHLTConfiguration ("esd-converter", "GlobalEsdConverter", esdInput.Data(), Bz.Data());
  new AliHLTConfiguration ("globalVtx", "GlobalVertexer"   , "esd-converter", "");
  new AliHLTConfiguration ("v0Histo", "V0Histo"   , "globalVtx", "");

  addInput( histoInput, "globalVtx v0Histo");  

  new AliHLTConfiguration("histFilter", "BlockFilter", histoInput.Data(), "-datatype 'ROOTHIST' '    '");
  new AliHLTConfiguration("histoRootFile", "ROOTFileWriter", "histFilter", "-datafile hltHisto -concatenate-events -overwrite");    
  
  new AliHLTConfiguration("sink-esd-file", "EsdCollector" , "globalVtx", "");

  option+="histoRootFile,sink-esd-file";


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;


  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // rec.SetDefaultStorage("local://$HOME/HCDB"); 
   rec.SetSpecificStorage("GRP/GRP/Data",
  			 Form("local://%s",gSystem->pwd()));

  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("HLT");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetOption("HLT",option);
  // switch off cleanESD
  rec.SetCleanESD(kFALSE);
  //rec.SetEventRange(0, 0);
  rec.Run();
}
