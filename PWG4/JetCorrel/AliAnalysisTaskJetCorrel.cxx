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
// Main class for two-particle correlations.
// Calls AliJetCorrelSelector and AliJetCorrelMaker for setup, then
// AliJetCorrelReader for ESD/AOD input reading into CorrelList_t lists, then
// AliJetCorrelMixer for event mixing and AliJetCorrelWriter for output histos
//-- Author: Paul Constantin

#include "AliAnalysisTaskJetCorrel.h"

using namespace std;

ClassImp(AliAnalysisTaskJetCorrel)

AliAnalysisTaskJetCorrel::AliAnalysisTaskJetCorrel() : 
  AliAnalysisTaskSE("JetCorrelTask"), fjcESD(NULL), fOutputContainer(NULL),
  fSelector(new AliJetCorrelSelector), fNumCorrel(0), fNumTrigg(0), fNumAssoc(0), fNumEvts(0),
  fMaker(new AliJetCorrelMaker), fWriter(new AliJetCorrelWriter), fReader(new AliJetCorrelReader),
  fMixer(new AliJetCorrelMixer), fTriggList(NULL), fAssocList(NULL) {
  // default constructor
}

AliAnalysisTaskJetCorrel::AliAnalysisTaskJetCorrel(AliJetCorrelSelector *s) : 
  AliAnalysisTaskSE("JetCorrelTask"), fjcESD(NULL), fOutputContainer(NULL),
  fSelector(s), fNumCorrel(0), fNumTrigg(0), fNumAssoc(0), fNumEvts(0),
  fMaker(new AliJetCorrelMaker), fWriter(new AliJetCorrelWriter), fReader(new AliJetCorrelReader),
  fMixer(new AliJetCorrelMixer), fTriggList(NULL), fAssocList(NULL) {
  // constructor
  fNumCorrel = fSelector->NoOfCorrel();
  if(!fMaker->Init(fNumCorrel,fSelector->CorrelTypes()))
    {std::cerr<<"AliJetCorrelMaker initialization failed. Bye!"<<std::endl; exit(-1);}
  fNumTrigg  = fMaker->NoOfTrigg();
  fNumAssoc  = fMaker->NoOfAssoc();
//  fMaker->Show();

  fWriter->Init(fSelector,fMaker);
  fReader->Init(fSelector,fWriter);
  fMixer->Init(fSelector,fMaker,fWriter);

  fTriggList = new CorrelList_t[fNumTrigg];
  fAssocList = new CorrelList_t[fNumAssoc];

  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
}

AliAnalysisTaskJetCorrel::~AliAnalysisTaskJetCorrel(){
  // destructor
  if(fTriggList) {delete [] fTriggList;}
  if(fAssocList) {delete [] fAssocList;}
  if(fOutputContainer) {fOutputContainer->Clear(); delete fOutputContainer;}
  if(fMaker) delete fMaker;
  if(fReader) delete fReader;
  if(fWriter) delete fWriter;
  if(fMixer) delete fMixer;
  fNumEvts=0;
}

void AliAnalysisTaskJetCorrel::ConnectInputData(Option_t *) {
  // connects to input data stream
  TTree* tree = dynamic_cast<TTree*>(GetInputData(0));
  if(tree){
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(esdH){
      fjcESD = dynamic_cast<AliESDEvent*>(esdH->GetEvent());
      if(!fjcESD) {std::cerr<<"AliAnalysisTaskJetCorrel::ConnectInputData - ERROR: no event"<<std::endl; exit(-1);}
    } else {std::cerr<<"AliAnalysisTaskJetCorrel::ConnectInputData - ERROR: no ESD Input"<<std::endl; exit(-1);}
  } else {std::cerr<<"AliAnalysisTaskJetCorrel::ConnectInputData - ERROR: no input tree"<<std::endl; exit(-1);}
}

void AliAnalysisTaskJetCorrel::CreateOutputObjects(){  
  // call writer object methods for histogram booking inside the output container list
  // OpenFile(0);
  if(!fOutputContainer) {
    fOutputContainer = new TList();
    fOutputContainer->SetName("JetCorrelHistos");
  }
  fWriter->CreateGeneric(fOutputContainer);
  if(fSelector->GenQA()) fWriter->CreateQA(fOutputContainer);
  fWriter->CreateCorrelations(fOutputContainer);
}

void AliAnalysisTaskJetCorrel::Exec(Option_t */*option*/){
  // get the event and pass it to the data reader object
  //fjcESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fjcESD)
    {std::cerr<<"AliAnalysisTaskJetCorrel::Exec() - ERROR: Cannot get event "<<fNumEvts<<std::endl; exit(-1);}
  fReader->SetEvent(fjcESD);

  // get global event pars and apply global cuts
  //  if(!fSelector->SelectedEvtTrigger(fjcESD)) return; // using AliPhysicsSelection
  Float_t cent = fReader->GetMultiplicity(); // use multiplicity in p-p
  Float_t zvtx = fReader->GetVertex();
  Int_t cBin = fSelector->GetBin(t_cent,cent);
  Int_t vBin = fSelector->GetBin(t_vert,zvtx);
  if(cBin<0 || vBin<0 || fReader->VtxOutPipe()) return; // event fails centrality or vertex selection
  fWriter->FillGlobal(cent,zvtx);
  fNumEvts++;

  // loop over correlations
  for(UInt_t iCor=0; iCor<fNumCorrel; iCor++){
    UInt_t idxTrigg = fMaker->IdxTrigg(iCor);
    UInt_t idxAssoc = fMaker->IdxAssoc(iCor);
    Bool_t aFilled = fAssocList[idxAssoc].Filled();

    fTriggList[idxTrigg].Label(fMaker->TriggType(iCor),fNumEvts);
    fAssocList[idxAssoc].Label(fMaker->AssocType(iCor),fNumEvts);
    fReader->FillLists(&fTriggList[idxTrigg],&fAssocList[idxAssoc]); // trigger list first!

    UInt_t nTriggs = fTriggList[idxTrigg].Size();
    UInt_t nAssocs = fAssocList[idxAssoc].Size();    
//     std::cout<<" Correl:"<<fMaker->Descriptor(iCor)<<"("<<idxTrigg<<"/"<<idxAssoc<<")"
// 	     <<" triggs:"<<nTriggs<<" assocs:"<<nAssocs<<std::endl;

    if(nTriggs<1 && nAssocs<1) continue;
    fWriter->FillSingleHistos(cBin, &fTriggList[idxTrigg], idxTrigg, &fAssocList[idxAssoc], idxAssoc);
    
    if(!aFilled) fMixer->FillPool(&fAssocList[idxAssoc], idxAssoc, vBin, cBin);

    if(nTriggs>0){
      CrossCorrelate(&fTriggList[idxTrigg], &fAssocList[idxAssoc], cBin, vBin, iCor); // same-event correlation
      fMixer->CurrTrigList(&fTriggList[idxTrigg]);
      fMixer->Mix(vBin, cBin, idxAssoc, iCor); // mixed-event correlation
    }
    
    PostData(0, fOutputContainer);
  } // loop over correlations

  // clear the lists
  for(UInt_t it=0; it<fNumTrigg; it++) fTriggList[it].Reset();
  for(UInt_t ia=0; ia<fNumAssoc; ia++) fAssocList[ia].Reset();
}

void AliAnalysisTaskJetCorrel::Terminate(Option_t */*option*/){
  // clean pool, print stats
  fMixer->CleanPool();
//   std::cout<<"CorrelParticle="<<sizeof(CorrelParticle_t)
// 	   <<" CorrelTrack="<<sizeof(CorrelTrack_t)
// 	   <<" CorrelRecoParent="<<sizeof(CorrelRecoParent_t)<<std::endl;
//   fWriter->ShowStats();
}

void AliAnalysisTaskJetCorrel::CrossCorrelate(CorrelList_t * const TriggList, CorrelList_t * const AssocList,
					      UInt_t cBin, UInt_t vBin, UInt_t iCor){
  // do the same event histogram filling
  if(TriggList->Size()<1 || AssocList->Size()<1) return;
  CorrelListIter_t iterTrigg = TriggList->Head();
  CorrelListIter_t iterAssoc = AssocList->Head();
  while(!iterTrigg.HasEnded()){
    while(!iterAssoc.HasEnded()){
      fWriter->FillCorrelations(0,iCor,cBin,vBin,iterTrigg.Data(),iterAssoc.Data()); // trigg first!
      iterAssoc.Move();
    }
    iterAssoc = AssocList->Head(); // reset associated particle iterator to list head
    iterTrigg.Move();
  }
  iterTrigg = TriggList->Head(); // reset trigger particle iterator to list head
}
