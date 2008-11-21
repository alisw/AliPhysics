// $Id$
/**
 * Configuration for the rec-krypton-hist.C macro
 */
{
  int iMinSlice=21; 
  int iMaxSlice=21;
  int iMinPart=1;
  int iMaxPart=1;

  TString arg, publisher, noise, noiseInput, krypton,kryptonAP, clusHist,clusDump, rootFileWriter, rootFileWriter2, histoInput1, histoInput2, kelly, clusterDumpInput, clusterDumpInput2,activepads, activePadsInput;
  
  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    for (int part=iMinPart; part<=iMaxPart; part++) {
      
      // raw data publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);

      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL, arg.Data());
      
      krypton.Form("KF_%02d_%d", slice, part);
      AliHLTConfiguration kryptonconf(krypton.Data(), "TPCKryptonClusterFinder", publisher.Data(), "");
      
      activepads.Form("AP_%02d_%d", slice, part);
      activePadsInput+=publisher.Data();
      activePadsInput+=" ";
      activePadsInput+=krypton.Data();
      AliHLTConfiguration activepadsconf(activepads.Data(), "AltroChannelSelector", activePadsInput.Data(),"");
      activePadsInput="";
      kryptonAP.Form("KFAP_%02d_%d", slice, part);
      AliHLTConfiguration kryptonconf2(kryptonAP.Data(), "TPCKryptonClusterFinder", activepads.Data(), "");

      if(histoInput1.Length()>0) histoInput1+=" ";
      histoInput1+=krypton;

      if(clusterDumpInput.Length()>0) clusterDumpInput+=" ";
      clusterDumpInput+=krypton;

      if(clusterDumpInput2.Length()>0) clusterDumpInput2+=" ";
      clusterDumpInput2+=kryptonAP;


    } // end of partition loop

    clusHist.Form("HG_%02d", slice);
    AliHLTConfiguration clusHistconf(clusHist.Data(), "TPCClusterHisto", histoInput1.Data(), "");
    histoInput1="";
    
    if (histoInput2.Length()>0) histoInput2+=" ";
    histoInput2+=clusHist;

 } // end of slice loop

  TString argdump;
  argdump.Form("-directory /home/kenneth/SimpleComponentWrapper/aliroot_configs/100PerEvent/clusterDump");
  cout<<"clusterDumpInput    "<<clusterDumpInput<<endl;
  AliHLTConfiguration clusDumpconf("dumpOut", "TPCClusterDump", clusterDumpInput.Data(),"-directory-clusterdump /home/kenneth/SimpleComponentWrapper/aliroot_configs/100PerEvent/clusterDump/");

  TString argdump2;
  argdump2.Form("-directory-clusterdump /home/kenneth/SimpleComponentWrapper/aliroot_configs/100PerEvent/clusterDumpAP/");

  AliHLTConfiguration clusDumpconf2("dumpOut2", "TPCClusterDump", clusterDumpInput2.Data(), argdump2.Data());

  if (histoInput2.Length()>0) histoInput2+=" ";
  histoInput2+="dumpOut";

  if (histoInput2.Length()>0) histoInput2+=" ";
  histoInput2+="dumpOut2";

  AliHLTConfiguration histconf("hist", "TPCHistogramHandler", histoInput2.Data(),"-sum-krypton-histograms 1");
  
  
  //root file writer component
  rootFileWriter.Form("-datafile %s","partitionHist.root");
  AliHLTConfiguration rootFileconf("sink1", "ROOTFileWriter", "hist", rootFileWriter.Data());

}
