// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTReadoutListDumpComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Base class for writer components to store data in a ROOT file

                                                                          */

#include "AliHLTReadoutListDumpComponent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTCTPData.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTReadoutListDumpComponent)

AliHLTReadoutListDumpComponent::AliHLTReadoutListDumpComponent()
  : AliHLTFileWriter()
  , fMode(AliHLTReadoutListDumpComponent::kModeBinaryList)
  , fBitsHisto(NULL)
  , fBitsVsCTP(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTReadoutListDumpComponent::~AliHLTReadoutListDumpComponent()
{
  // see header file for class documentation
}

int AliHLTReadoutListDumpComponent::InitWriter()
{
  // see header file for class documentation
  int iResult=0;
  fBitsHisto=CreateReadoutListHistogram();
  fBitsVsCTP=CreateReadoutListVsCTPHistogram();
  if (!fBitsHisto || !fBitsVsCTP) return -ENOMEM;

  return iResult;
}

int AliHLTReadoutListDumpComponent::CloseWriter()
{
  // see header file for class documentation
  TString filename=GetDirectory();
  if (!filename.IsNull() && !filename.EndsWith("/")) filename+="/";
  filename+=fBitsHisto->GetName();
  filename+="_";
  filename+=GetRunNo();
  filename+=".root";

  TFile out(filename, "RECREATE");
  fBitsHisto->Write();
  fBitsVsCTP->Write();
  out.Close();

  delete fBitsHisto;
  fBitsHisto=NULL;
  delete fBitsVsCTP;
  fBitsVsCTP=NULL;
  return 0;
}

int AliHLTReadoutListDumpComponent::DumpEvent( const AliHLTComponentEventData& /*evtData*/,
					    const AliHLTComponentBlockData* /*blocks*/, 
					    AliHLTComponentTriggerData& trigData )
{
  // see header file for class documentation
  int iResult=0;
  if (!IsDataEvent()) return 0;

  if (fMode==AliHLTReadoutListDumpComponent::kModeBinaryList) {
    const AliHLTComponentDataType hltrdlstdt=AliHLTComponentDataTypeInitializer("HLTRDLST", "HLT ");
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
	 pBlock && iResult>=0;
	 pBlock=GetNextInputBlock()) {
      if (pBlock->fDataType!=hltrdlstdt) continue;
      if (pBlock->fSize==sizeof(AliHLTEventDDL) or pBlock->fSize==sizeof(AliHLTEventDDLV0)) {
	HLTDebug("Filling histograms from binary buffer");
	AliHLTReadoutList readoutlist(*reinterpret_cast<AliHLTEventDDL*>(pBlock->fPtr));
	FillReadoutListHistogram(fBitsHisto, &readoutlist);
	FillReadoutListVsCTP(fBitsVsCTP, &readoutlist, &trigData);
      } else {
	HLTError("HLTRDLST size missmatch: %d, expected %d or %d", pBlock->fSize, sizeof(AliHLTEventDDL), sizeof(AliHLTEventDDLV0));
      }
    }
  } else if (fMode==AliHLTReadoutListDumpComponent::kModeHLTDecision) {
    for (const TObject* pObject=GetFirstInputObject();
	 pObject && iResult>=0;
	 pObject=GetNextInputObject()) {
      const AliHLTTriggerDecision* pDecision=dynamic_cast<const AliHLTTriggerDecision*>(pObject);
      if (!pDecision) continue;

      AliHLTReadoutList list=pDecision->ReadoutList();
      HLTDebug("Filling histograms from HLT decision object");
      FillReadoutListHistogram(fBitsHisto, &list);
      FillReadoutListVsCTP(fBitsVsCTP, &list, &trigData);
    }
  } else {
    HLTError("invalid mode %d", fMode);
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTReadoutListDumpComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=-EINVAL;
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -binary
  if (argument.CompareTo("-binary")==0) {
    fMode=AliHLTReadoutListDumpComponent::kModeBinaryList;
    return 1;
  }

  // -decision
  if (argument.CompareTo("-decision")==0) {
    fMode=AliHLTReadoutListDumpComponent::kModeHLTDecision;
    return 1;
  }

  return iResult;
}

TH1I* AliHLTReadoutListDumpComponent::CreateReadoutListHistogram()
{
  // see header file for class documentation
  int bins=gkAliHLTDDLListSize*sizeof(AliHLTUInt32_t)*8;
  TH1I* histo=new TH1I("HLTRDLST","HLT readout list", bins, 0, bins);
  return histo;
}

TH2I* AliHLTReadoutListDumpComponent::CreateReadoutListVsCTPHistogram()
{
  // see header file for class documentation
  int bins=gkAliHLTDDLListSize*sizeof(AliHLTUInt32_t)*8;
  TH2I* histo=new TH2I("HLTRDLSTvsCTP","HLT readout list vs. CTP trigger", bins, 0, bins, gkNCTPTriggerClasses, 0, gkNCTPTriggerClasses);
  return histo;
}

int AliHLTReadoutListDumpComponent::FillReadoutListHistogram(TH1I* histo, const AliHLTReadoutList* list)
{
  // see header file for class documentation
  if (!histo || !list) return -EINVAL;
  if (list->BufferSize()!=sizeof(AliHLTEventDDL)) return -EBADF;
  if (list->Buffer()->fCount!=(unsigned)gkAliHLTDDLListSize) return -EBADF;

  for (int word=0; word<gkAliHLTDDLListSize; word++) {
    for (unsigned bit=0; bit<sizeof(AliHLTUInt32_t)*8; bit++) {
      if (list->Buffer()->fList[word]&0x1<<bit) histo->Fill(word*sizeof(AliHLTUInt32_t)*8+bit);
    }
  }
  
  return 0;
}

int AliHLTReadoutListDumpComponent::FillReadoutListVsCTP(TH2I* histo, const AliHLTReadoutList* list, const AliHLTComponentTriggerData* trigData)
{
  // see header file for class documentation
  if (!histo || !list || !trigData) return -EINVAL;
  if (list->BufferSize()!=sizeof(AliHLTEventDDL)) return -EBADF;
  if (list->Buffer()->fCount!=(unsigned)gkAliHLTDDLListSize) return -EBADF;

  AliHLTUInt64_t triggerMask=AliHLTCTPData::ActiveTriggers(*trigData);
  AliHLTUInt64_t bit0=0x1;
  for (int word=0; word<gkAliHLTDDLListSize; word++) {
    for (unsigned bit=0; bit<sizeof(AliHLTUInt32_t)*8; bit++) {
      if (list->Buffer()->fList[word]&0x1<<bit) {
	for (int trigger=0; trigger<gkNCTPTriggerClasses; trigger++) {
	  if ((triggerMask&(bit0<<trigger))!=0) {
	    histo->Fill(word*sizeof(AliHLTUInt32_t)*8+bit, trigger);
	  }
	}
      }
    }
  }

  return 0;
}
