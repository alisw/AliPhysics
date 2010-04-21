/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTRDEsdWriterComponent.cxx
    @author Mateusz Ploskon
    @date   
    @brief  Writer component to store tracks of the HLT TRD

                                                                          */
#include "AliHLTTRDEsdWriterComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"
#include "AliHLTCTPData.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTRDtrackV1.h"
#include "TTree.h"
#include "TFile.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDEsdWriterComponent)

AliHLTTRDEsdWriterComponent::AliHLTTRDEsdWriterComponent()
:AliHLTRootFileWriterComponent()
  ,fTree(NULL)
  ,fFrTree(NULL)
  ,fESD(NULL)
  ,fESDfriend(NULL)
  ,fFile(NULL)
  ,fFrFile(NULL)
  ,fTracksArray(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTRDEsdWriterComponent::AliHLTTRDEsdWriterComponent(const AliHLTTRDEsdWriterComponent&)
  :AliHLTRootFileWriterComponent()
  ,fTree(NULL)
  ,fFrTree(NULL)
  ,fESD(NULL)
  ,fESDfriend(NULL)
  ,fFile(NULL)
  ,fFrFile(NULL)
  ,fTracksArray(NULL)
{
}

void AliHLTTRDEsdWriterComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data  
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTTRDDefinitions::fgkTracksDataType );
  list.push_back( AliHLTTRDDefinitions::fgkHiLvlTracksDataType );
}

AliHLTTRDEsdWriterComponent::~AliHLTTRDEsdWriterComponent()
{
  // see header file for class documentation
}

int AliHLTTRDEsdWriterComponent::InitWriter()
{
  // see header file for class documentation
  
  fFile = new TFile("AliHLTTRDESDs.root", "recreate");
  fESD = new AliESDEvent;
  fESD->CreateStdContent();
  fTree = new TTree("esdTree", "Tree with HLT::TRD ESD objects");
  fESD->WriteToTree(fTree);
  fFrFile = new TFile("AliHLTTRDESDfriends.root", "recreate");
  fESDfriend = new AliESDfriend();
  fFrTree = new TTree("esdFriendTree", "Tree with HLT::TRD ESD Friend objects");
  fFrTree->Branch("ESDfriend.","AliESDfriend", &fESDfriend);
  fESD->AddObject(fESDfriend);
  fFile->cd();
  fTree->GetUserInfo()->Add(fESD);
  fTracksArray = new TClonesArray("AliTRDtrackV1");

  SetupCTPData();

  return 0;
}

int AliHLTTRDEsdWriterComponent::CloseWriter()
{
  // see header file for class documentation

  //fTree->Print();
  fFile->cd();
  fTree->Write(fTree->GetName(),TObject::kOverwrite);
  fFile->Write();
  fFile->Close();
  fFrFile->cd();
  fFrTree->Write(fFrTree->GetName(),TObject::kOverwrite);
  fFrFile->Write();
  fFrFile->Close();
  delete fFile; fFile=0;
  delete fFrFile; fFrFile=0;
  //delete fTree;
  delete fTracksArray; fTracksArray=0;

  return AliHLTRootFileWriterComponent::CloseWriter();
}

int AliHLTTRDEsdWriterComponent::DumpEvent( const AliHLTComponentEventData& /*evtData*/,
					    const AliHLTComponentBlockData* /*blocks*/, 
					    AliHLTComponentTriggerData& trigData )
{
  TClonesArray* TCAarray[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t usedEntries = 0;
  Int_t blockOrObject = 0;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(AliHLTTRDDefinitions::fgkTracksDataType); pBlock; pBlock=GetNextInputBlock()) 
    {
      TCAarray[0] = fTracksArray;
      AliHLTTRDUtils::ReadTracks(TCAarray[0], pBlock->fPtr, pBlock->fSize);
      usedEntries = 1;
      blockOrObject = -1;
    }

  for(const TObject *iter = GetFirstInputObject(AliHLTTRDDefinitions::fgkHiLvlTracksDataType); iter; iter = GetNextInputObject()) 
    {
      if(blockOrObject<0){
	HLTError("You may not mix high level and low level!");
	return -1;
      }

      TCAarray[usedEntries] = dynamic_cast<TClonesArray*>(const_cast<TObject*>(iter));
      if(!TCAarray[usedEntries])continue;
      usedEntries++;
      blockOrObject = 1;
    }

  if(!blockOrObject)
    return 0;

  fESD->Reset(); 
  fESD->SetMagneticField(GetBz());
  fESD->SetRunNumber(GetRunNo());
  fESD->SetPeriodNumber(GetPeriodNumber());
  fESD->SetOrbitNumber(GetOrbitNumber());
  fESD->SetBunchCrossNumber(GetBunchCrossNumber());
  fESD->SetTimeStamp(GetTimeStamp());
  fESD->SetEventType(7);

  const AliHLTCTPData* pCTPData=CTPData();
  if (pCTPData) {
    AliHLTUInt64_t mask=pCTPData->ActiveTriggers(trigData);
    for (int index=0; index<gkNCTPTriggerClasses; index++) {
      if ((mask&((AliHLTUInt64_t)0x1<<index)) == 0) continue;
      fESD->SetTriggerClass(pCTPData->Name(index), index);
    }
    fESD->SetTriggerMask(mask);
  }
  
  for(int i=0; i<usedEntries; i++){
    const TClonesArray* inArr = TCAarray[i];
    for(int ii=0; ii<inArr->GetEntriesFast(); ii++){
      AliTRDtrackV1* inV1 = (AliTRDtrackV1*)inArr->UncheckedAt(ii);
      AliESDtrack *esdTrack = new AliESDtrack();
      esdTrack->UpdateTrackParams(inV1, AliESDtrack::kTRDout);
      esdTrack->SetLabel(inV1->GetLabel());
      inV1->UpdateESDtrack(esdTrack);
      AliTRDtrackV1 *calibTrack = new AliTRDtrackV1(*inV1);
      calibTrack->SetOwner();
      esdTrack->AddCalibObject(calibTrack);
      fESD->AddTrack(esdTrack);
      delete esdTrack;
    }
  }
  
  fESD->GetESDfriend(fESDfriend);
  Int_t nb = fTree->Fill();
  HLTInfo("Tree filled with %i bytes\n", nb);  
  nb = fFrTree->Fill();
  HLTInfo("FrTree filled with %i bytes\n", nb);
  fESD->Reset();
  fESDfriend->~AliESDfriend();
  new (fESDfriend) AliESDfriend();

  if(blockOrObject<0){
    TCAarray[0]->Delete();
  }

  return 0;
}
