
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>           *
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
 * @file   noiseMapTest.C
 * @author Kalliopi.Kanaki@ift.uib.no
 * @date   
 * @brief  Test macro for noise map component in combination with histogram handler

 * Usage:
 * <pre>
 *   aliroot -b -q noiseMapTest.C | tee noiseMapTest.log
 *   aliroot -b -q noiseMapTest.C'("./")' | tee noiseMapTest.log
 * </pre>
 *
 * The macro assumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q noiseMapTest.C'("input.root")'
 * </pre>
 *
 * In the first section an analysis chain is defined. The scale of the
 * chain can be defined by choosing the range of sectors and partitions.
 *
 * The reconstruction is steered by the AliReconstruction object in the
 * usual way.
 *
 * @ingroup alihlt_tpc
 */

void noiseMapHistoHandler(const char* input="./"){
  
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // 
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
    
  int iMinSlice = 35; 
  int iMaxSlice = 35;
  int iMinPart  = 0;
  int iMaxPart  = 5;

  TString rootFileWriter, histoInput;

  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    for (int part=iMinPart; part<=iMaxPart; part++) {
     
      TString argument, publisher, noise;
     
      // raw data publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      argument.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x", ddlno, slice, slice, part, part);
      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL, argument.Data());
      
      // noise filling
      noise.Form("NM_%02d_%d", slice, part);
      AliHLTConfiguration noiseconf(noise.Data(), "TPCNoiseMap", publisher.Data(), "");
      
      
      if (histoInput.Length()>0) histoInput+=" ";
      histoInput+=noise;
 
    }
  }
  
  //cout << "histogram handler input:   " << histoInput.Data() << endl;
  
  AliHLTConfiguration histconf("hist", "TPCHistogramHandler", histoInput.Data(),"-use-general -ignore-specification");
  //AliHLTConfiguration histconf("hist", "TPCHistogramHandler", histoInput.Data(),"-sum-noise-histograms");
  
  //root file writer component
  
  //AliHLTConfiguration rootFileconf("sink1", "ROOTFileWriter", noise.Data(), "-datafile partition");
  AliHLTConfiguration rootFileconf("sink1", "ROOTFileWriter", "hist", "-datafile TPC -concatenate-events -concatenate-blocks");
  
  
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
  rec.SetOption("HLT", "libAliHLTTPC.so loglevel=0x7c chains=sink1");
  rec.Run();
}
