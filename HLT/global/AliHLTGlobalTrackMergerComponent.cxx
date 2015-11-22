// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jacek Otwinowski <Jacek.Otwinowski@gsi.de>            *
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

/** @file   AliHLTGlobalTrackMergerComponent.cxx
    @author Jacek Otwinowski
    @date   
    @brief  HLT global track merger component.
*/

#include <climits>
#include <cassert>
#include <cstdlib>
#include <cerrno>

#include "AliESDEvent.h"
#include "AliTracker.h"
#include "AliHLTGlobalTrackMerger.h"
#include "AliHLTGlobalTrackMergerComponent.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalTrackMergerComponent)

//_____________________________________________________________________________
AliHLTGlobalTrackMergerComponent::AliHLTGlobalTrackMergerComponent() : AliHLTGlobalEsdConverterComponent(), 
  fGlobalTrackMerger(0)
{
}

//_____________________________________________________________________________
AliHLTGlobalTrackMergerComponent::~AliHLTGlobalTrackMergerComponent()
{
  // see header file for class documentation
}

//_____________________________________________________________________________
const char* AliHLTGlobalTrackMergerComponent::GetComponentID()
{
  // see header file for class documentation
  return "GlobalTrackMerger";
}

//_____________________________________________________________________________
void AliHLTGlobalTrackMergerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back( kAliHLTDataTypeTrack );
}

//_____________________________________________________________________________
// AliHLTComponentDataType AliHLTGlobalTrackMergerComponent::GetOutputDataType()
// {
//   // see header file for class documentation
//   return kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC;
// }

//_____________________________________________________________________________
// void AliHLTGlobalTrackMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
// {
//   // see header file for class documentation
//   // XXX TODO: Find more realistic values.
//   constBase = 20000;
//   inputMultiplier = 1.0;
// }

//_____________________________________________________________________________
AliHLTComponent* AliHLTGlobalTrackMergerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalTrackMergerComponent;
}

//_____________________________________________________________________________
void AliHLTGlobalTrackMergerComponent::SetMergerParameters(Double_t maxy,Double_t maxz,Double_t maxkappa,Double_t maxpsi,Double_t maxtgl)
{
  // see header file for class documentation
  fGlobalTrackMerger->SetParameter( maxy, maxz, maxkappa, maxpsi, maxtgl );
}

//_____________________________________________________________________________
int AliHLTGlobalTrackMergerComponent::DoInit( int /*argc*/, const char** /*argv*/ )
{
  // see header file for class documentation
  const char* argv = "-notree";
  AliHLTGlobalEsdConverterComponent::DoInit(1,&argv);

  int iResult = 0;
  
  // Init merger
  fGlobalTrackMerger = new AliHLTGlobalTrackMerger();

  if (!fGlobalTrackMerger ) {
     HLTError("failed creating internal objects");
     iResult=-ENOMEM;
     return iResult;
  }

  SetMergerParameters();

  return iResult;
}

//_____________________________________________________________________________
int AliHLTGlobalTrackMergerComponent::DoDeinit()
{
  // see header file for class documentation
  if(fGlobalTrackMerger) delete fGlobalTrackMerger; fGlobalTrackMerger =0;
  AliHLTGlobalEsdConverterComponent::DoDeinit();
  return 0;
}

//_____________________________________________________________________________
int AliHLTGlobalTrackMergerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, 
					      AliHLTComponentTriggerData& /*trigData*/)
{
  //
  // global track merger function
  // takes TRD and TPC tracks and merges them
  //
  HLTInfo("DoEvent processing data");

  // see header file for class documentation
  int iResult=0;

  if(!fGlobalTrackMerger || !fESD) {
    HLTError("component not initialized");
    iResult=-ENOMEM;
    return iResult;
  }

  fESD->Reset(); 
  fESD->SetMagneticField(fSolenoidBz);
  //fESD->SetMagneticField(AliTracker::GetBz());

  if ((iResult=ProcessBlocks(NULL, fESD, NULL))<0) return iResult;
  if(!(fESD->GetNumberOfTracks()>0)) return iResult;

   // merge tracks
   Bool_t isMerged = fGlobalTrackMerger->Merge(fESD);
   if(!isMerged) {
     HLTInfo("No merged tracks");
   }

   // try to propagate all tracks to DCA to primary vertex
   fGlobalTrackMerger->PropagateTracksToDCA(fESD);

   // calculate specification
   // AliHLTUInt32_t iSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, 0, 5 );
   // HLTInfo("minSlice %d, maxSlice %d", minSlice, maxSlice);

   // send output data
   //PushBack(fESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC, iSpecification);
   PushBack(fESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginOut);

   // clean ESD event content
   fESD->Reset();
  
return iResult;
}

	
