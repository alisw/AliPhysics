/*************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                           CutQA                                      //
//                                                                      //
/*
   Allow to monitor how many tracks,pair,events pass the selection criterion 
   in any of the cuts added to the corresponding filters. All you need to 
   add to your config is the following:

   dielectron->SetCutQA();


*/
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliDielectronCutQA.h"

#include <TList.h>
#include <TCollection.h>

#include "AliDielectronCutGroup.h"
#include "AliAnalysisCuts.h"

#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliDielectronPair.h"


ClassImp(AliDielectronCutQA)


AliDielectronCutQA::AliDielectronCutQA() :
  TNamed(),
  fQAHistArray()
{
  //
  // Default constructor
  //
  for(Int_t itype=0; itype<kNtypes; itype++) {
    fCutQA[itype]=0x0;
    fNCuts[itype]=1;
    for(Int_t i=0; i<20; i++) {
      fCutNames[i][itype]="";
    }
  }
  fTypeKeys[kTrack] = "Track";
  fTypeKeys[kPair]  = "Pair";
  fTypeKeys[kEvent] = "Event";
  fQAHistArray.SetOwner();
}

//_____________________________________________________________________
AliDielectronCutQA::AliDielectronCutQA(const char* name, const char* title) :
  TNamed(name, title),
  fQAHistArray()
{
  //
  // Named Constructor
  //
  for(Int_t itype=0; itype<kNtypes; itype++) {
    fCutQA[itype]=0x0;
    fNCuts[itype]=1;
    for(Int_t i=0; i<20; i++) {
      fCutNames[i][itype]="";
    }
  }
  fTypeKeys[kTrack] = "Track";
  fTypeKeys[kPair]  = "Pair";
  fTypeKeys[kEvent] = "Event";
  fQAHistArray.SetOwner();
}

//_____________________________________________________________________
AliDielectronCutQA::~AliDielectronCutQA() 
{
  //
  //Default Destructor
  //
  fQAHistArray.Delete();
}

//_____________________________________________________________________
void AliDielectronCutQA::Init()
{

  fQAHistArray.SetName(Form("%s",GetName()));

  // loop over all types
  for(Int_t itype=0; itype<kNtypes; itype++) {
    //    printf("\n type: %d\n",itype);
    fCutNames[0][itype]="no cuts";

    // create histogram based on added cuts
    fCutQA[itype] = new TH1F(fTypeKeys[itype],
			     Form("%sQA;cuts;# passed %ss",fTypeKeys[itype],fTypeKeys[itype]),
			     fNCuts[itype],0,fNCuts[itype]);
    // loop over all cuts
    for(Int_t i=0; i<fNCuts[itype]; i++) {
      fCutQA[itype]->GetXaxis()->SetBinLabel(i+1,fCutNames[i][itype]);
      //      printf(" %s \n",fCutNames[i][itype]);
    }
    fQAHistArray.AddLast(fCutQA[itype]);
  }

}

//_____________________________________________________________________
void AliDielectronCutQA::AddTrackFilter(AliAnalysisFilter *trackFilter)
{
  //
  // add track filter cuts to the qa histogram
  //


  TIter listIterator(trackFilter->GetCuts());
  while (AliAnalysisCuts *thisCut = (AliAnalysisCuts*) listIterator()) {
    Bool_t addCut=kTRUE;

    // add new cut class to the array
    if(addCut) {
      fCutNames[fNCuts[kTrack]][kTrack]=thisCut->GetTitle();
      //      printf("add cut %s to %d \n",thisCut->GetTitle(),fNCuts[kTrack]);
      fNCuts[kTrack]++;
    }

  } // pair filter loop

}


//_____________________________________________________________________
void AliDielectronCutQA::AddPairFilter(AliAnalysisFilter *pairFilter)
{
  //
  // add track filter cuts to the qa histogram
  //


  TIter listIterator(pairFilter->GetCuts());
  while (AliAnalysisCuts *thisCut = (AliAnalysisCuts*) listIterator()) {
    Bool_t addCut=kTRUE;

    // add new cut class to the array
    if(addCut) {
      fCutNames[fNCuts[kPair]][kPair]=thisCut->GetTitle();
      //  printf("add cut %s to %d \n",thisCut->GetTitle(),fNCuts[kPair]);
      fNCuts[kPair]++;
    }

  } // trk filter loop

}

//_____________________________________________________________________
void AliDielectronCutQA::AddEventFilter(AliAnalysisFilter *eventFilter)
{
  //
  // add track filter cuts to the qa histogram
  //


  TIter listIterator(eventFilter->GetCuts());
  while (AliAnalysisCuts *thisCut = (AliAnalysisCuts*) listIterator()) {
    Bool_t addCut=kTRUE;

    // add new cut class to the array
    if(addCut) {
      fCutNames[fNCuts[kEvent]][kEvent]=thisCut->GetTitle();
      //      printf("add cut %s to %d \n",thisCut->GetTitle(),fNCuts[kEvent]);
      fNCuts[kEvent]++;
    }

  } // trk filter loop

}

//_____________________________________________________________________
void AliDielectronCutQA::Fill(UInt_t mask, TObject *obj)
{
  //
  // fill the corresponding step in the qa histogram
  //

  UInt_t idx = GetObjIndex(obj);

  Int_t cutstep=1;
  for (Int_t iCut=0; iCut<fNCuts[idx]-1;++iCut) {
    //    UInt_t cutMask=1<<iCut;         // for each cut
    UInt_t cutMask=(1<<(iCut+1))-1; // increasing cut match

    if ((mask&cutMask)==cutMask) {
      fCutQA[idx]->Fill(cutstep);
      ++cutstep;
    }

  }

}

//_____________________________________________________________________
void AliDielectronCutQA::FillAll(TObject *obj)
{
  //
  // fill the corresponding step in the qa histogram
  //

  UInt_t idx = GetObjIndex(obj);
  fCutQA[idx]->Fill(0);

}

//______________________________________________________________________
UInt_t AliDielectronCutQA::GetObjIndex(TObject *obj)
{
  //
  // return the corresponding idex
  //
  //  printf("INFO: object type is a %s \n", obj->IsA()->GetName());
  if(obj->InheritsFrom(AliDielectronPair::Class())) return kPair;
  if(obj->InheritsFrom(AliVTrack::Class())        ) return kTrack;
  if(obj->InheritsFrom(AliVParticle::Class())     ) return kTrack;
  if(obj->InheritsFrom(AliVEvent::Class())        ) return kEvent;
  printf("FATAL: object type %s not yet supported, please let the author know\n", obj->IsA()->GetName());
  return -1;
  //TODO complete

}




