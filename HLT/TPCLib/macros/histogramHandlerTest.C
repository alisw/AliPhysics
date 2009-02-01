
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors Kenneth Aamodt <kenneth.aamodt@cern.ch>                *
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
 * @file   histogramHandlerTest.C
 * @author Kenneth.Aamodt@cern.ch
 * @date   
 * @brief  Test macro for histogram handler
 *
 * To test the histogram handler the AliHLTTPCClusterHistoComponent is used as input to the 
 * histogram handler. There are two clusterfinders running, each of these sends its output to 
 * one clusterhisto component each, which again sends its data to 2 histogram handlers.
 * One of the Histogram handlers get data from one of the cluster histo components, while the
 * other get the output from both.
 * The content of the two rootfiles can now be compared to eachother.
 * The histograms in histogramHandlerFile2... should now have double the amount of entries, and the
 * height in the y-axis should be 2 times that of the histogramHandlerFile1...
 */

void histogramHandlerTest(const char* input="./"){
  
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // 
  // init the HLT system in order to define the analysis chain below
  //
  gSystem->Load("libHLTrec.so");
  AliHLTSystem* gHLT=AliHLTReconstructorBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
    
  int iMinSlice=0; 
  int iMaxSlice=35;
  int iMinPart=0;
  int iMaxPart=5;

  TString histogramHandlerInput1;
  TString histogramHandlerInput2;
  TString histogramHandlerOutput1;
  TString histogramHandlerOutput2;
  TString rootFileWriterInput1;
  TString rootFileWriterInput2;

  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString clusterFinderOutput1;
      TString clusterFinderOutput2;
      TString clusterHistoInput1;
      TString clusterHistoInput2;
      TString clusterHistoOutput1;
      TString clusterHistoOutput2;
      TString arg, publisher;
      // digit publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", ddlno, slice, slice, part, part);

      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

      // first clusterfinder
      clusterFinderOutput1.Form("CF1_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(clusterFinderOutput1.Data(), "TPCClusterFinderDecoder", publisher.Data(), "-timebins 446");//-timebins set to simulated data
      if (clusterHistoInput1.Length()>0) clusterHistoInput1+=" ";
      clusterHistoInput1+=clusterFinderOutput1;

      // second clusterfinder
      clusterFinderOutput2.Form("CF2_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(clusterFinderOutput2.Data(), "TPCClusterFinderDecoder", publisher.Data(), "-timebins 446");//-timebins set to simulated data
      if (clusterHistoInput2.Length()>0) clusterHistoInput2+=" ";
      clusterHistoInput2+=clusterFinderOutput2;
    

      // first cluster histo component
      clusterHistoOutput1.Form("CH1_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(clusterHistoOutput1.Data(), "TPCClusterHisto", clusterHistoInput1.Data(), "");
      if (histogramHandlerInput1.Length()>0) histogramHandlerInput1+=" ";
      histogramHandlerInput1+=clusterHistoOutput1;

      //second cluster histo component
      clusterHistoOutput2.Form("CH2_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(clusterHistoOutput2.Data(), "TPCClusterHisto", clusterHistoInput2.Data(), "");
      if (histogramHandlerInput1.Length()>0) histogramHandlerInput1+=" ";
      histogramHandlerInput1+=clusterHistoOutput2;
      if (histogramHandlerInput2.Length()>0) histogramHandlerInput2+=" ";
      histogramHandlerInput2+=clusterHistoOutput2;

    }
  }
  
  // first histogram handler component
  histogramHandlerOutput1.Form("HH1_%02d_%d", slice, part);
  AliHLTConfiguration cfconf(histogramHandlerOutput1.Data(), "TPCHistogramHandler", histogramHandlerInput1.Data(), "-use-general");
  if (rootFileWriterInput1.Length()>0) rootFileWriterInput1+=" ";
  rootFileWriterInput1+=histogramHandlerOutput1;
  
  // second histogram handler component
  histogramHandlerOutput2.Form("HH2_%02d_%d", slice, part);
  AliHLTConfiguration cfconf(histogramHandlerOutput2.Data(), "TPCHistogramHandler", histogramHandlerInput2.Data(), "-use-general");
  if (rootFileWriterInput2.Length()>0) rootFileWriterInput2+=" ";
  rootFileWriterInput2+=histogramHandlerOutput2;
  
  AliHLTConfiguration rootFileWriter1("RootFileWriter1", "ROOTFileWriter", rootFileWriterInput1.Data() , "-datafile histogramHandlerFile1");

  AliHLTConfiguration rootFileWriter2("RootFileWriter2", "ROOTFileWriter", rootFileWriterInput2.Data() , "-datafile histogramHandlerFile2");


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  //rec.SetFillESD("HLT");
  rec.SetFillESD("");
  rec.SetFillTriggerESD(false);
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=RootFileWriter1,RootFileWriter2");
  rec.Run();
}
