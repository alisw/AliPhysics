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
#include "TTimeStamp.h"
#include "TSystem.h"

#include "AliHLTGlobalVertexerHistoComponent.h"

ClassImp(AliHLTGlobalVertexerHistoComponent)

AliHLTGlobalVertexerHistoComponent::AliHLTGlobalVertexerHistoComponent()
  : 
  fUID(0),
  fRefreshPeriod(1000),
  fFillSecond(0),
  fFillSecondSPD(0)
{
  //Default ctor.
  for( int i=0; i<2; i++ ){
    fPrimaryXY[i].SetName("primVertexXY");
    fPrimaryXY[i].SetTitle("HLT: Primary vertex distribution in XY");
    fPrimaryXY[i].SetMarkerStyle(8);
    fPrimaryXY[i].SetMarkerSize(0.4);
    fPrimaryXY[i].SetYTitle("Y [cm]");
    fPrimaryXY[i].SetXTitle("X [cm]");
    //fPrimaryXY[i].SetStats(0);
    //fPrimaryXY[i].SetBit(TH1::kCanRebin);
    
    fPrimaryX[i].SetName("primVertexX");
    fPrimaryX[i].SetTitle("HLT: Primary vertex distribution in X");
    fPrimaryX[i].SetMarkerStyle(8);
    fPrimaryX[i].SetMarkerSize(0.4);
    fPrimaryX[i].SetXTitle("X [cm]");
    //fPrimaryX[i].SetStats(0);
    //fPrimaryX[i].SetBit(TH1::kCanRebin);
    
    fPrimaryY[i].SetName("primVertexY");
    fPrimaryY[i].SetTitle("HLT: Primary vertex distribution in Y");
    fPrimaryY[i].SetMarkerStyle(8);
    fPrimaryY[i].SetMarkerSize(0.4);
    fPrimaryY[i].SetXTitle("Y [cm]");
    //fPrimaryY[i].SetStats(0);
    //fPrimaryX[i].SetBit(TH1::kCanRebin);
    
    fPrimaryZ[i].SetName("primVertexZ");
    fPrimaryZ[i].SetTitle("HLT: Primary vertex distribution in Z");
    fPrimaryZ[i].SetMarkerStyle(8);
    fPrimaryZ[i].SetMarkerSize(0.4);
    fPrimaryZ[i].SetXTitle("Z [cm]");
    //fPrimaryZ[i].SetStats(0);
    //fPrimaryX[i].SetBit(TH1::kCanRebin);
    
    
    fSPDVertexXY[i].SetName("spdVertexXY");
    fSPDVertexXY[i].SetTitle("HLT: SPD vertex distribution in XY");
    fSPDVertexXY[i].SetMarkerStyle(8);
    fSPDVertexXY[i].SetMarkerSize(0.4);
    fSPDVertexXY[i].SetYTitle("Y [cm]");
    fSPDVertexXY[i].SetXTitle("X [cm]");
    //fSPDVertexXY[i].SetStats(0);
    //fSPDVertexXY[i].SetBit(TH1::kCanRebin);
    
    fSPDVertexX[i].SetName("spdVertexX");
    fSPDVertexX[i].SetTitle("HLT: SPD vertex distribution in X");
    fSPDVertexX[i].SetMarkerStyle(8);
    fSPDVertexX[i].SetMarkerSize(0.4);
    fSPDVertexX[i].SetXTitle("X [cm]");
    //fSPDVertexX[i].SetStats(0);
    //fSPDVertexX[i].SetBit(TH1::kCanRebin);
    
    fSPDVertexY[i].SetName("spdVertexY");
    fSPDVertexY[i].SetTitle("HLT: SPD vertex distribution in Y");
    fSPDVertexY[i].SetMarkerStyle(8);
    fSPDVertexY[i].SetMarkerSize(0.4);
    fSPDVertexY[i].SetXTitle("Y [cm]");
    //fSPDVertexY[i].SetStats(0);
    //fSPDVertexY[i].SetBit(TH1::kCanRebin);
    
    fSPDVertexZ[i].SetName("spdVertexZ");
    fSPDVertexZ[i].SetTitle("HLT: SPD vertex distribution in Z");
    fSPDVertexZ[i].SetMarkerStyle(8);
    fSPDVertexZ[i].SetMarkerSize(0.4);
    fSPDVertexZ[i].SetXTitle("Z [cm]");
    //fSPDVertexZ[i].SetStats(0);
    //fSPDVertexZ[i].SetBit(TH1::kCanRebin);
  }
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
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS );
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut );
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);
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

void AliHLTGlobalVertexerHistoComponent::SetDefaultConfiguration()
{
  // Set default configuration for the CA tracker component
  // Some parameters can be later overwritten from the OCDB

  fRefreshPeriod = 1000;
}

int AliHLTGlobalVertexerHistoComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the CA tracker component from the string

  int iResult = 0;
  if ( !arguments ) return iResult;

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );

  int nArgs =  pTokens ? pTokens->GetEntries() : 0;

  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    if ( argument.IsNull() ) continue;

    if ( argument.CompareTo( "-refresh" ) == 0 ) {
      if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
      fRefreshPeriod = ( ( TObjString* )pTokens->At( i ) )->GetString().Atoi();
      HLTInfo( "N events for refresh of the run vertex is set to: %d", fRefreshPeriod );
      continue;
    }
    

    HLTError( "Unknown option \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }
  delete pTokens;

  if ( bMissingParam ) {
    HLTError( "Specifier missed for parameter \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }

  return iResult;
}


int AliHLTGlobalVertexerHistoComponent::DoInit(int argc, const char** argv)
{
  //Clear all bin contents and statistics.

  SetDefaultConfiguration();

  
  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  int ret = ReadConfigurationString( arguments.Data() );
  for( int i=0; i<2; i++ ){
    fPrimaryXY[i].Reset();
    fPrimaryX[i].Reset();
    fPrimaryY[i].Reset();
    fPrimaryZ[i].Reset();
    fSPDVertexXY[i].Reset();
    fSPDVertexX[i].Reset();
    fSPDVertexY[i].Reset();
    fSPDVertexZ[i].Reset();
    //Set bin numbers and axis ranges[i].
    fPrimaryXY[i].SetBins(1000,  -1.,  1., 1000, -1., 1.);
    fPrimaryX[i].SetBins( 1000, -1., 1.);
    fPrimaryY[i].SetBins( 1000, -1., 1.);
    fPrimaryZ[i].SetBins( 1000, -30., 30.); 
    fSPDVertexXY[i].SetBins(1000,  -1.,  1., 1000, -1., 1.);
    fSPDVertexX[i].SetBins(1000,  -1.,  1.);
    fSPDVertexY[i].SetBins(1000,  -1.,  1.);
    fSPDVertexZ[i].SetBins(1000,  -30.,  30.);
  }
  fFillSecond = 0;
  fFillSecondSPD = 0;
  fUID = 0;
  return ret;
}

int AliHLTGlobalVertexerHistoComponent::DoDeinit()
{
  //Nothing to delete or "do de init" yet.
  fUID = 0;
  return 0;
}

int AliHLTGlobalVertexerHistoComponent::DoEvent(const AliHLTComponentEventData& evtData,
                                                const AliHLTComponentBlockData* /*blocks*/,
                                                AliHLTComponentTriggerData& /*trigData*/,
                                                AliHLTUInt8_t* /*outputPtr*/,
                                                AliHLTUInt32_t& size,
                                                AliHLTComponentBlockDataList& /*outputBlocks*/)
{
  //Fill histogramms.
  if (GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR))
    return 0;

  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
    //cout<<"\nSet id to "<<fUID<<endl;
  }

  if (GetOutputDataSize() > size) {
    Logging(kHLTLogFatal, "HLT::AliHLTGlobalVertexerHistoComponent::DoEvent", "Too much data", 
            "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
            GetOutputDataSize(), size);
    return -ENOSPC;
  }

  if( fRefreshPeriod>0 ){
    if( fPrimaryXY[0].GetEntries()+1>=fRefreshPeriod/2 ) fFillSecond = 1;
    if( fSPDVertexXY[0].GetEntries()+1>=fRefreshPeriod/2 ) fFillSecondSPD = 1;
  }

  int iResult = 0;

  const AliESDVertex *vertexITS = 0;
  const AliESDVertex *vertexGlobal = 0;


  {
    const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS);
    if( iter != NULL  ) {
      if( !( vertexITS = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) ) ) ){    
	HLTError("ITS SPD vertex object is corrupted");
	iResult = -EINVAL;    
      }
    }
  }

  {
    const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginOut);
    if( iter != NULL ){
      if( !( vertexGlobal = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) ) ) ){
	HLTError("Global vertex object is corrupted");
	iResult = -EINVAL;    
      }
    }
  }

  if( !vertexITS || !vertexGlobal ){

    for (const TObject* iter = GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginOut); iter; iter = GetNextInputObject() ) {
      if ( AliESDEvent* event = dynamic_cast<AliESDEvent*>((TObject*)iter) ) {
	event->GetStdContent();
	if( !vertexGlobal ) vertexGlobal = event->GetPrimaryVertexTracks();
	if( !vertexITS ) vertexITS = event->GetPrimaryVertexSPD();      
      } else {
	HLTError("ESD event object is corrupted");
	iResult = -EINVAL;    
      }
    }
  }
  
  if (vertexGlobal && vertexGlobal->GetNContributors() >= 5) {
    for( int i=0; i<=fFillSecond; i++ ){
      if(  fRefreshPeriod>0 && fPrimaryXY[i].GetEntries()>=fRefreshPeriod ){
	fPrimaryXY[i].Reset();
	fPrimaryX[i].Reset();
	fPrimaryY[i].Reset();
	fPrimaryZ[i].Reset();      
      }
      fPrimaryXY[i].Fill(vertexGlobal->GetX(), vertexGlobal->GetY());
      fPrimaryX[i].Fill(vertexGlobal->GetX());
      fPrimaryY[i].Fill(vertexGlobal->GetY());
      fPrimaryZ[i].Fill(vertexGlobal->GetZ());
    }
  }

  if( vertexITS && vertexITS->GetNContributors() >= 5 ){
    for( int i=0; i<=fFillSecondSPD; i++ ){
      if(  fRefreshPeriod>0 && fSPDVertexXY[i].GetEntries()>=fRefreshPeriod ){
	fSPDVertexXY[i].Reset();
	fSPDVertexX[i].Reset();
	fSPDVertexY[i].Reset();
	fSPDVertexZ[i].Reset();
      }
      fSPDVertexXY[i].Fill(vertexITS->GetX(), vertexITS->GetY());
      fSPDVertexX[i].Fill(vertexITS->GetX() );
      fSPDVertexY[i].Fill(vertexITS->GetY() );
      fSPDVertexZ[i].Fill(vertexITS->GetZ() );	  
    }
  }
  
  int i = ( fPrimaryXY[1].GetEntries() > fPrimaryXY[0].GetEntries() );

  PushBack(&fPrimaryXY[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, fUID);
  PushBack(&fPrimaryZ[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, fUID);
  PushBack(&fPrimaryX[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, fUID);
  PushBack(&fPrimaryY[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginOut, fUID);

  i = ( fSPDVertexXY[1].GetEntries() > fSPDVertexXY[0].GetEntries() );

  //cout<<"bla NEntr     = "<<fPrimaryXY[0].GetEntries()<<" / "<<fPrimaryXY[1].GetEntries()<<endl;
  //cout<<"bla NEntr SPD = "<<fSPDVertexXY[0].GetEntries()<<" / "<<fSPDVertexXY[1].GetEntries()<<endl;
  
  PushBack(&fSPDVertexXY[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, fUID);
  PushBack(&fSPDVertexZ[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, fUID);
  PushBack(&fSPDVertexX[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, fUID);
  PushBack(&fSPDVertexY[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginITSSPD, fUID);
  
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
