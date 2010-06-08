/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id: $ */

//__________________________________________
// Event mixing class. A 5-dimensinal pool fPool is maintained:
// type_idx(fixed), vertex(fixed), centrality(fixed), 
// and event(dynamic), particle(dynamic)
// fixed dimensions are fixed at task initialization time (via AliJetCorrelSelector)
// event&particle are allow to float during runtime (via linked lists TList=event & CorrelList_t=particle)
//-- Author: Paul Constantin

#include "AliJetCorrelMixer.h"

using namespace std;

ClassImp(AliJetCorrelMixer)

AliJetCorrelMixer::AliJetCorrelMixer() :
  fSelector(NULL), fMaker(NULL), fWriter(NULL),
  fTriggEvnt(NULL), fAssocEvnt(NULL),
  fAssocIter(NULL), fTriggIter(NULL) {
  // constructor
}

AliJetCorrelMixer::~AliJetCorrelMixer(){
  // destructor
  CleanPool();
}

void AliJetCorrelMixer::Init(AliJetCorrelSelector * const s, AliJetCorrelMaker * const m, AliJetCorrelWriter * const w){
  // initialization method
  fSelector = s;
  fMaker = m;
  fWriter = w;
  
  for(UInt_t vBin=0; vBin<fSelector->NoOfBins(t_vert); vBin++)
    for(UInt_t cBin=0; cBin<fSelector->NoOfBins(t_cent); cBin++)
      for(UInt_t ia=0; ia<fMaker->NoOfAssoc(); ia++)
	fPool[ia][vBin][cBin] = new TList;
}

void AliJetCorrelMixer::FillPool(CorrelList_t *partList, UInt_t pIdx, UInt_t vBin, UInt_t cBin){
  // pool filling method
  if(partList->Size()<1) return;
  UInt_t pSize = fPool[pIdx][vBin][cBin]->GetSize();
  // when pool depth is reached, pop pool before new event push (keep const depth)
  if(pSize>=fSelector->PoolDepth()) fPool[pIdx][vBin][cBin]->RemoveFirst();
  // incoming list is cleared at end-of-event; hence, store in pool a deep copy:
  fPool[pIdx][vBin][cBin]->AddLast(partList->DeepCopy());
}

void AliJetCorrelMixer::Mix(UInt_t vBin, UInt_t cBin, UInt_t ia, UInt_t ic) {
  // rolling buffer mixing method
  TListIter* iterPool=(TListIter*)fPool[ia][vBin][cBin]->MakeIterator();
  if(!fTriggEvnt) {std::cerr<<"AliJetCorrelMixer::Mix - ERROR: Trigger list not set!"<<std::endl; exit(-1);}
  fTriggIter = fTriggEvnt->Head();
  while((fAssocEvnt=(CorrelList_t*)iterPool->Next())){
    fAssocIter = fAssocEvnt->Head();  
    if(fTriggEvnt->EvtID()==fAssocEvnt->EvtID()) continue; // don't mix same event!
    
    while(!fTriggIter.HasEnded()){
      while(!fAssocIter.HasEnded()){
	fWriter->FillCorrelations(1,ic,cBin,vBin,fTriggIter.Data(),fAssocIter.Data()); // trigg first!
	fAssocIter.Move();
      } // loop over associated particles
      fAssocIter = fAssocEvnt->Head(); // reset associated particle iterator to list head
      fTriggIter.Move();
    } // loop over trigger particles
    fTriggIter = fTriggEvnt->Head(); // reset trigger particle iterator to list head
  } // loop over associated pool
  delete iterPool;
}

void AliJetCorrelMixer::CleanPool(){
  // pool cleaning
  UInt_t size = fMaker->NoOfAssoc();
  for(UInt_t k=0; k<size; k++){
    for(UInt_t vBin=0; vBin<fSelector->NoOfBins(t_vert); vBin++)
      for(UInt_t cBin=0; cBin<fSelector->NoOfBins(t_cent); cBin++)
	fPool[k][vBin][cBin]->Delete(); // Remove all list objects AND delete all heap based objects
  }
}

void AliJetCorrelMixer::ShowSummary(UInt_t pIdx, UInt_t vBin, UInt_t cBin) const {
  // pool printout method
  UInt_t totalPoolSize=0;
  TListIter* iter=(TListIter*)fPool[pIdx][vBin][cBin]->MakeIterator();
  CorrelList_t* partList;
  while((partList=(CorrelList_t*)iter->Next())) totalPoolSize += partList->Size();
  delete iter;
  std::cout<<"Pool["<<vBin<<"]["<<cBin<<"]: nevt="<<fPool[pIdx][vBin][cBin]->GetSize()<<" npart="<<totalPoolSize<<std::endl;
}
