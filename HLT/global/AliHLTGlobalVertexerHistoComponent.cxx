//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timur.Pocheptsov@cern.ch                              *
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

/** @file   AliHLTGlobalVertexerComponent.h
    @author Timur Pocheptsov
    @brief  Component for monitor primary vertex
*/

#include "AliESDVertex.h"
#include "AliESDEvent.h"

#include "AliHLTGlobalVertexerHistoComponent.h"

ClassImp(AliHLTGlobalVertexerHistoComponent)

AliHLTGlobalVertexerHistoComponent::AliHLTGlobalVertexerHistoComponent()
  : 
  fPrimaryXY(),
  fPrimaryZX(),
  fPrimaryZY(),
  fSPDVertexXY(),
  fSPDVertexX(),
  fSPDVertexY(),
  fSPDVertexZ()
{
  //Default ctor.
  fPrimaryXY.SetName("primVertexXY");
  fPrimaryXY.SetTitle("HLT: Primary vertex distribution in XY");
  fPrimaryXY.SetMarkerStyle(8);
  fPrimaryXY.SetMarkerSize(0.4);
  fPrimaryXY.SetYTitle("Y [cm]");
  fPrimaryXY.SetXTitle("X [cm]");
  fPrimaryXY.SetStats(0);
  fPrimaryXY.SetBit(TH1::kCanRebin);

  fPrimaryZX.SetName("primVertexZX");
  fPrimaryZX.SetTitle("HLT: Primary vertex distribution in ZX");
  fPrimaryZX.SetMarkerStyle(8);
  fPrimaryZX.SetMarkerSize(0.4);
  fPrimaryZX.SetYTitle("X [cm]");
  fPrimaryZX.SetXTitle("Z [cm]");
  fPrimaryZX.SetStats(0);
  fPrimaryZX.SetBit(TH1::kCanRebin);

  fPrimaryZY.SetName("primVertexZY");
  fPrimaryZY.SetTitle("HLT: Primary vertex distribution in ZY");
  fPrimaryZY.SetMarkerStyle(8);
  fPrimaryZY.SetMarkerSize(0.4);
  fPrimaryZY.SetYTitle("Y [cm]");
  fPrimaryZY.SetXTitle("Z [cm]");
  fPrimaryZY.SetStats(0);
  fPrimaryZY.SetBit(TH1::kCanRebin);


  fSPDVertexXY.SetName("spdVertexXY");
  fSPDVertexXY.SetTitle("HLT: SPDVertex vertex distribution in XY");
  fSPDVertexXY.SetMarkerStyle(8);
  fSPDVertexXY.SetMarkerSize(0.4);
  fSPDVertexXY.SetYTitle("Y [cm]");
  fSPDVertexXY.SetXTitle("X [cm]");
  fSPDVertexXY.SetStats(0);
  fSPDVertexXY.SetBit(TH1::kCanRebin);

  fSPDVertexX.SetName("spdVertexX");
  fSPDVertexX.SetTitle("HLT: SPDVertex vertex distribution in X");
  fSPDVertexX.SetMarkerStyle(8);
  fSPDVertexX.SetMarkerSize(0.4);
  fSPDVertexX.SetXTitle("X [cm]");
  fSPDVertexX.SetStats(0);
  fSPDVertexX.SetBit(TH1::kCanRebin);

  fSPDVertexY.SetName("spdVertexY");
  fSPDVertexY.SetTitle("HLT: SPDVertex vertex distribution in Y");
  fSPDVertexY.SetMarkerStyle(8);
  fSPDVertexY.SetMarkerSize(0.4);
  fSPDVertexY.SetXTitle("Y [cm]");
  fSPDVertexY.SetStats(0);
  fSPDVertexY.SetBit(TH1::kCanRebin);
 
  fSPDVertexZ.SetName("spdVertexZ");
  fSPDVertexZ.SetTitle("HLT: SPDVertex vertex distribution in Z");
  fSPDVertexZ.SetMarkerStyle(8);
  fSPDVertexZ.SetMarkerSize(0.4);
  fSPDVertexZ.SetXTitle("Z [cm]");
  fSPDVertexZ.SetStats(0);
  fSPDVertexZ.SetBit(TH1::kCanRebin);
}

const char* AliHLTGlobalVertexerHistoComponent::GetComponentID()
{
  //Unique component id.
  return "GlobalVertexerHisto";
}

void AliHLTGlobalVertexerHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  //
  list.clear();
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS );
}



AliHLTComponentDataType AliHLTGlobalVertexerHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTGlobalVertexerHistoComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeHistogram|kAliHLTDataOriginOut);
  tgtList.push_back(kAliHLTDataTypeHistogram|kAliHLTDataOriginITSSPD);
 return tgtList.size();
}


void AliHLTGlobalVertexerHistoComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //8000 is a temporary number. Find the way to
  //estimate 3 serialized TH2F. 
  constBase = GetOutputDataSize();
  inputMultiplier = 1.;
}

AliHLTComponent* AliHLTGlobalVertexerHistoComponent::Spawn()
{
  //Should be const, I think, but it's not.
  return new AliHLTGlobalVertexerHistoComponent;
}

int AliHLTGlobalVertexerHistoComponent::DoInit(int /*argc*/, const char** /*argv*/)
{
  //Clear all bin contents and statistics.
  fPrimaryXY.Reset();
  fPrimaryZX.Reset();
  fPrimaryZY.Reset();
  fSPDVertexXY.Reset();
  fSPDVertexX.Reset();
  fSPDVertexY.Reset();
  fSPDVertexZ.Reset();
  //Set bin numbers and axis ranges.
  fPrimaryXY.SetBins(60,  -2.,  2., 60, -2., 2.);
  fPrimaryZX.SetBins(60, -15., 15., 60, -2., 2.);
  fPrimaryZY.SetBins(60, -10., 10., 60, -2., 2.);
  fSPDVertexXY.SetBins(100,  -2.,  2., 100, -2., 2.);
  fSPDVertexX.SetBins(100,  -2.,  2.);
  fSPDVertexY.SetBins(100,  -2.,  2.);
  fSPDVertexZ.SetBins(100,  -15.,  15.);

  return 0;
}

int AliHLTGlobalVertexerHistoComponent::DoDeinit()
{
  //Nothing to delete or "do de init" yet.
  return 0;
}

int AliHLTGlobalVertexerHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
                                                const AliHLTComponentBlockData* /*blocks*/,
                                                AliHLTComponentTriggerData& /*trigData*/,
                                                AliHLTUInt8_t* /*outputPtr*/,
                                                AliHLTUInt32_t& size,
                                                AliHLTComponentBlockDataList& /*outputBlocks*/)
{
  //Fill histogramms.
  if (GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR))
    return 0;

  if (GetOutputDataSize() > size) {
    Logging(kHLTLogFatal, "HLT::AliHLTGlobalVertexerHistoComponent::DoEvent", "Too much data", 
            "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
            GetOutputDataSize(), size);
    return -ENOSPC;
  }

  int iResult = 0;

  bool wasITS = 0;
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS); iter != NULL; iter = GetNextInputObject() ) {
    if( AliESDVertex *vertex = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) ) ){
      if( vertex && vertex->GetNContributors() >= 3 ){
	fSPDVertexXY.Fill(vertex->GetX(), vertex->GetY());
	fSPDVertexX.Fill(vertex->GetX() );
	fSPDVertexY.Fill(vertex->GetY() );
	fSPDVertexZ.Fill(vertex->GetZ() );
	wasITS = 1;
      }
    } else {
      HLTError("ITS SPD vertex object is corrupted");
      iResult = -EINVAL;    
    }
  }

  for (const TObject* iter = GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut); iter; iter = GetNextInputObject() ) {
    if (AliESDEvent* event = dynamic_cast<AliESDEvent*>((TObject*)iter)) {
      event->GetStdContent();
      const AliESDVertex* vertex = event->GetPrimaryVertexTracks();
      if (vertex && vertex->GetNContributors() >= 3) {
        fPrimaryXY.Fill(vertex->GetX(), vertex->GetY());
        fPrimaryZX.Fill(vertex->GetZ(), vertex->GetX());
        fPrimaryZY.Fill(vertex->GetZ(), vertex->GetY());
      }
      if( !wasITS ){
	vertex = event->GetPrimaryVertexSPD();
	if( vertex && vertex->GetNContributors() >= 3 ){
	  fSPDVertexXY.Fill(vertex->GetX(), vertex->GetY());
	  fSPDVertexX.Fill(vertex->GetX() );
	  fSPDVertexY.Fill(vertex->GetY() );
	  fSPDVertexZ.Fill(vertex->GetZ() );	  
	}
      }
    } else {
      HLTError("ESD event object is corrupted");
      iResult = -EINVAL;    
    }
  }
 
  PushBack(&fPrimaryXY, kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, 0);
  PushBack(&fPrimaryZX, kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, 0);
  PushBack(&fPrimaryZY, kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, 0);

  PushBack(&fSPDVertexXY, kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, 0);
  PushBack(&fSPDVertexX, kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, 0);
  PushBack(&fSPDVertexY, kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, 0);
  PushBack(&fSPDVertexZ, kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, 0);
  
  return iResult; 
}

int AliHLTGlobalVertexerHistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  //Do nothing at the moment, but at least log (base class version).
  return AliHLTComponent::Reconfigure(cdbEntry, chainId);
}

unsigned long AliHLTGlobalVertexerHistoComponent::GetOutputDataSize()const
{
  //8000 is a temporary number. Find the way to
  //estimate serialized TH2F. 
  return 20000;
}
