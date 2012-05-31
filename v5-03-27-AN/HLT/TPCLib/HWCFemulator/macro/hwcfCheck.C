// $Id: rec-hlt-tpc.C 36001 2009-10-26 17:08:56Z kaamodt $
/*
 * Example macro to run the HLT Conformal mapping tracker embedded into
 * AliRoot reconstruction. The reconstruction is done from the TPC raw
 * data.
 *
 * Usage:
 * <pre>
 *   rm galice.root & aliroot -b -q hwcfCheck.C | tee rec-hlt-tpc.log
 *   rm galice.root & aliroot -b -q hwcfCheck.C'("./")' | tee rec-hlt-tpc.log
 * </pre>
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   rm galice.root & aliroot -b -q hwcfCheck.C'("./")'
 * </pre>
 *
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 */

TString GetFileList( char *basedir, char *filename )
{
  TString s = "";
  int evnt=0;
  do{
    TString dir;
    dir.Form("%sraw%d",basedir,evnt);
    if( gSystem->AccessPathName(dir.Data()) ) break;
    TString add;
    add.Form(" -datafile %s/%s", dir.Data(), filename);
    if( evnt>0 ) s+=" -nextevent";
    s+=add;
    evnt++;
  }while(1);
  return s;
}

void hwcfCheck(const char* basedir="./")
{
  
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }

  if (!basedir) {
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
  int clusterFinderType=0; // 0 = v3; 1 = decoder; 2 = packed (offline v1)
  bool bUseCA=true;   // use the CA tracker and merger
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  int iMinSlice=0;
  int iMaxSlice=0;
  int iMinPart=0;
  int iMaxPart=0;
  TString clusters;
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {    
    for (int part=iMinPart; part<=iMaxPart; part++) {
 
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;     
 
      // raw data publisher 
      
      TString arg, publisher, cf, publisherHW;

      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
            

      // FPGA cluster finder emulator

      cf.Form("CF_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(cf.Data(), "TPCHWClusterFinderEmulator", publisher.Data(), "");
            
      if(clusters.Length()>0) clusters+=" ";
      clusters+=cf;
      
      TString fname;
      fname.Form("HWCF_%d.ddl",ddlno);
      TString flist = GetFileList(basedir,fname.Data());

      cout<<"\n\nFileList for FPGA ddl "<<ddlno<<":"<<endl;
      cout<<flist.Data()<<endl;
      cout<<endl;

      arg.Form("-datatype 'HWCLUST1' 'TPC ' %s -dataspec 0x%02x%02x%02x%02x", flist.Data(), slice, slice, part, part);

      publisherHW.Form("DPHW_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisherHW.Data(), "FilePublisher", "", arg.Data());
          
      if(clusters.Length()>0) clusters+=" ";
      clusters+=publisherHW;
    }
  }

  
  AliHLTConfiguration pubconf("ConsistencyControl", "TPCHWCFConsistenyControl", clusters.Data(), "");
  
  option+="ConsistencyControl";


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off
  //
  AliReconstruction rec;
  rec.SetInput(basedir);
  rec.SetRunVertexFinder(0);
  rec.SetRunVertexFinder(0);
  rec.SetRunVertexFinderTracks(0);
  rec.SetRunCascadeFinder(0);
  rec.SetRunMultFinder(0);
  rec.SetRunReconstruction("HLT");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");   
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  rec.SetOption("HLT", option);
  rec.Run();
}
