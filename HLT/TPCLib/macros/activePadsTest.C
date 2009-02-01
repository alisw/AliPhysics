
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
 * @file   activePadsTest.C
 * @author Kenneth.Aamodt@cern.ch
 * @date   
 * @brief  Test macro of the active pad selection
 *
 * The macro reads the simulated TPC digits from the RunLoader and
 * prints the digit info to stout.
 *
 * The macro reads ddl data and apply zerosuppression on it. The zerosuppression
 * component then ships out a list of pads in addition to the zerosuppressed data.
 * This list of pads is here used in combination with the original raw data to
 * remove all pads with no data surviving the zerosuppression with the altro 
 * selection component. This is to study undershoot of signals etc. The output of
 * altro selection component is the NON zerosuppressed data of all the pads which
 * there should be a signal in. To test this the macro runs in several steps. 
 * Point 3 and 4 has to be done in this way to avoid confusion in what ddl data is 
 * read by the component. 
 * 1. read the data.
 * 2. dump it to file. (will be compared to the zerosuppressed data later)
 * 3. apply zerosuppression producing activepads list
 * 4. apply zerosuppression producing zerosuppressed data.
 * 5. dump the data from 4. to file
 * 6. send the active pad list together with the original data to the altro selection component
 * 7. dump the reduced data to file
 * 8. again apply zero suppression to this reduced data
 * 9. dump the result from zero suppression from 8.(on the reduced data)
 * 10. compare the outputs from 2 and 7 at see that they are different (if not; this is not a good test) (diff dump1/ev.... dump2/ev...)
 * 11. compare the outputs from 5. and 9. to see that they are alike. (diff dump3/ev... dump4/ev..)
 *
 * NB: 10. and 11. must by now be done manually.... will change soon 
 * aliroot -b -q activePadsTest.C
 *
 */

void activePadsTest(const char* input="./"){
  
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

//   AliCDBManager* pManager=AliCDBManager::Instance(NULL, 0);
//   pManager->SetRun(0);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
    
  int iMinSlice=8; 
  int iMaxSlice=8;
  int iMinPart=2;
  int iMaxPart=2;
  TString digitDumpInput_OriginalData;
  TString digitDumpInput_AfterAPSelection;
  TString digitDumpInput_ZeroSuppressedOriginalData;
  TString digitDumpInput_ZerosuppressedAPSelectionData;
  TString dumpHwAddressInput;

  for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
    for (int part=iMinPart; part<=iMaxPart; part++) {
      TString arg, publisher, inputAltroSelection, zsactivepadlist, zsoriginaldata, altroChannelSelector, zsACSdata, zsACSdataOutput;
      // digit publisher components
      int ddlno=768;
      if (part>1) ddlno+=72+4*slice+(part-2);
      else ddlno+=2*slice+part;
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);

      publisher.Form("DP_%02d_%d", slice, part);
      AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      if(digitDumpInput_OriginalData.Length() >0) digitDumpInput_OriginalData += " ";
      digitDumpInput_OriginalData += publisher;
      if(inputAltroSelection.Length() >0) inputAltroSelection += " ";
      inputAltroSelection += publisher;

      // zero suppression component (3.)
      zsactivepadlist.Form("ZSAP_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(zsactivepadlist.Data(), "TPCZeroSuppression", publisher.Data(), "signal-threshold 2 start-timebin 70 end-timebin 900 -skip-sending-data -send-hw-list");
      if (inputAltroSelection.Length()>0) inputAltroSelection+=" ";
      inputAltroSelection+=zsactivepadlist;
      if (dumpHwAddressInput.Length()>0) dumpHwAddressInput+=" ";
      dumpHwAddressInput+=zsactivepadlist;
    
      // zero suppression component (4.)
      zsoriginaldata.Form("ZSDDL_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(zsoriginaldata.Data(), "TPCZeroSuppression", publisher.Data(), "signal-threshold 2 start-timebin 70 end-timebin 900");
      if(digitDumpInput_ZeroSuppressedOriginalData.Length() >0) digitDumpInput_ZeroSuppressedOriginalData += " ";
      digitDumpInput_ZeroSuppressedOriginalData += zsoriginaldata;

      //altro channel selector
      altroChannelSelector.Form("ACS_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(altroChannelSelector.Data(), "AltroChannelSelector", inputAltroSelection.Data(), "");
      if(zsACSdata.Length() >0) zsACSdata += " ";
      zsACSdata+= altroChannelSelector;
      if(digitDumpInput_AfterAPSelection.Length() >0) digitDumpInput_AfterAPSelection += " ";
      digitDumpInput_AfterAPSelection += altroChannelSelector;

      // zero suppression component (8.)
      zsACSdataOutput.Form("ZSACS_%02d_%d", slice, part);
      AliHLTConfiguration cfconf(zsACSdataOutput.Data(), "TPCZeroSuppression", zsACSdata.Data(), "signal-threshold 2 start-timebin 70 end-timebin 900");
      if(digitDumpInput_ZerosuppressedAPSelectionData.Length() >0) digitDumpInput_ZerosuppressedAPSelectionData += " ";
      digitDumpInput_ZerosuppressedAPSelectionData += zsACSdataOutput;
      
    }
  }
  //Dumping data (2.)
  AliHLTConfiguration dump1("Dump1", "TPCDigitDump", digitDumpInput_OriginalData.Data() , "-digitreader decoder -directory dump1 -unsorted -concatenate-blocks");
  //Dumping data (7.)
  AliHLTConfiguration dump2("Dump2", "TPCDigitDump", digitDumpInput_AfterAPSelection.Data() , "-digitreader decoder -directory dump2 -unsorted -concatenate-blocks");
  //Dumping data (5.)
  AliHLTConfiguration dump3("Dump3", "TPCDigitDump", digitDumpInput_ZeroSuppressedOriginalData.Data() , "-digitreader decoder -directory dump3 -unsorted -concatenate-blocks");
  //Dumping data (9.)
  AliHLTConfiguration dump4("Dump4", "TPCDigitDump", digitDumpInput_ZerosuppressedAPSelectionData.Data() , "-digitreader decoder -directory dump4 -unsorted -concatenate-blocks");
  //Dumping hw lists
  AliHLTConfiguration dump5("Dump5", "FileWriter", dumpHwAddressInput.Data() , "-directory hwlists");


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
  rec.SetOption("HLT", "libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so loglevel=0x7c chains=Dump1,Dump2,Dump3,Dump4,Dump5");
  rec.Run();
}
