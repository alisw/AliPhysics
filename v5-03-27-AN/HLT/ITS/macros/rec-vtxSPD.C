// $Id: rec-tpc-its.C 35299 2009-10-07 09:07:29Z kkanaki $
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
 * In the first section, an analysis chain is defined. The scale of the
 * chain can be defined by choosing the range of sectors and partitions.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 * @author Gaute.Ovrebekk@ift.uib.no
 */
void rec_vtxSPD(const char* input="./")
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


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  int minddl=0;          //min ddl number for SPD
  int maxddl=19;         //max ddl number for SPD
  int spec=0x1;          //spec for ddl's
  int ddlno=0;
  TString cfout="";

  for(ddlno=minddl;ddlno<=maxddl;ddlno++){  
    TString arg, publisher, cf;
    
    arg.Form("-minid %d -datatype 'DDL_RAW ' 'ISPD ' -dataspec 0x%08x -verbose",ddlno, spec);
    publisher.Form("DP_%d", ddlno);
    new AliHLTConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
    
    cf.Form("CF_%d",ddlno);
    new AliHLTConfiguration(cf.Data(), "ITSClusterFinderSPD", publisher.Data(), "");
    
    if (cfout.Length()>0) cfout+=" ";
    cfout+=cf;
    
    spec=spec<<1;
  }

  
  AliHLTConfiguration itsvtx("itsvtx", "ITSVertexerSPD"   , cfout.Data(), "");

  AliHLTConfiguration rootFileWriterClusters("historootfile", "ROOTFileWriter", "itsvtx" , "-datafile ITSHistograms -concatenate-events -overwrite");


  AliHLTConfiguration globalConverter("globalConverter", "GlobalEsdConverter"   , "itsvtx", "");
  AliHLTConfiguration sink("esdfile", "EsdCollector"   , "globalConverter", "");
  option+="historootfile,esdfile ";

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
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd())); 
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  rec.SetRunQA(":");
  rec.SetOption("HLT", option);
  //rec.SetEventRange(0,0);
rec.Run();
}
