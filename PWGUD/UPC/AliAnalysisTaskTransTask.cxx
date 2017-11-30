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

/* AliAnaysisTaskTransTask
 *
 * simple task to study CTRUE events to estimate number of UPC events lost due to vetoes
 *
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskTransTask.h"

class AliAnalysisTaskTransTask;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskTransTask) // classimp: necessary for root

AliAnalysisTaskTransTask::AliAnalysisTaskTransTask() : AliAnalysisTaskSE(), 
  fAOD(0), fOutputList(0), fAnaTree(0), fRunNum(0), fTracklets(0), fCtrue(0),
  fL0inputs(0), fL1inputs(0), fZem1Energy(0), fZem2Energy(0),
  fZNCEnergy(0), fZNAEnergy(0), fZPCEnergy(0), fZPAEnergy(0),
  fV0ADecision(-10), fV0CDecision(-10), fADADecision(-10), fADCDecision(-10), 
  fIR1Map(0),fIR2Map(0),fCounter(0)

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskTransTask::AliAnalysisTaskTransTask(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0), fOutputList(0), fAnaTree(0), fRunNum(0), fTracklets(0), fCtrue(0),
  fL0inputs(0), fL1inputs(0), fZem1Energy(0), fZem2Energy(0),						       
  fZNCEnergy(0), fZNAEnergy(0), fZPCEnergy(0), fZPAEnergy(0),
  fV0ADecision(-10), fV0CDecision(-10), fADADecision(-10), fADCDecision(-10), 
  fIR1Map(0),fIR2Map(0),fCounter(0)

{
  for (Int_t i=0;i<4;i++) fZNATDC[i]=fZNCTDC[i]=fZPATDC[i]=fZPCTDC[i]=0;
    // constructor
    DefineInput(0, TChain::Class());   
    DefineOutput(1, TTree::Class());
    DefineOutput(2, TList::Class());   
                                        
}

//_____________________________________________________________________________
AliAnalysisTaskTransTask::~AliAnalysisTaskTransTask()
{
  // destructor
  if(fOutputList) {delete fOutputList;}
  if(fAnaTree) {delete fAnaTree;}
  if(fCounter) {delete fCounter;}
}

//_____________________________________________________________________________
void AliAnalysisTaskTransTask::UserCreateOutputObjects()
{
  // define list of histos
  fOutputList = new TList();  
  fOutputList->SetOwner(kTRUE);
  // define hist and add it to the list
  fCounter = new TH1I("fCounter", "fCounter", 100, 1, 101);
  fOutputList->Add(fCounter);
  // define tree
  fAnaTree = new TTree("fAnaTree", "fAnaTree");
  fAnaTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  fAnaTree ->Branch("fTracklets", &fTracklets, "fTracklets/I");	
  fAnaTree ->Branch("fCtrue", &fCtrue, "fCtrue/I");
  fAnaTree ->Branch("fC1zed", &fC1zed, "fC1zed/I");  
  fAnaTree ->Branch("fL0inputs", &fL0inputs, "L0inputs/i");
  fAnaTree ->Branch("fL1inputs", &fL1inputs, "L1inputs/i");      
  fAnaTree ->Branch("fZem1Energy", &fZem1Energy, "fZem1Energy/D");
  fAnaTree ->Branch("fZem2Energy", &fZem2Energy, "fZem2Energy/D");  
  fAnaTree ->Branch("fZNCEnergy", &fZNCEnergy, "fZNCEnergy/D");  
  fAnaTree ->Branch("fZNAEnergy", &fZNAEnergy, "fZNAEnergy/D");
  fAnaTree ->Branch("fZPCEnergy", &fZPCEnergy, "fZPCEnergy/D");
  fAnaTree ->Branch("fZPAEnergy", &fZPAEnergy, "fZPAEnergy/D");  
  fAnaTree ->Branch("fZNATDC", &fZNATDC[0], "fZNATDC[4]/D");
  fAnaTree ->Branch("fZNCTDC", &fZNCTDC[0], "fZNCTDC[4]/D");  
  fAnaTree ->Branch("fZPATDC", &fZPATDC[0], "fZPATDC[4]/D");
  fAnaTree ->Branch("fZPCTDC", &fZPCTDC[0], "fZPCTDC[4]/D");  
  fAnaTree ->Branch("fV0ADecision", &fV0ADecision, "fV0ADecision/I");
  fAnaTree ->Branch("fV0CDecision", &fV0CDecision, "fV0CDecision/I");  
  fAnaTree ->Branch("fADADecision", &fADADecision, "fADADecision/I");
  fAnaTree ->Branch("fADCDecision", &fADCDecision, "fADCDecision/I");
  fAnaTree ->Branch("fIR1Map", &fIR1Map);
  fAnaTree ->Branch("fIR2Map", &fIR2Map);  

 
  // post output
  PostData(1, fAnaTree);
  PostData(2, fOutputList);          
}

//_____________________________________________________________________________
void AliAnalysisTaskTransTask::UserExec(Option_t *)
{
  fCounter->Fill(1);
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {PostData(2, fOutputList); return;}
  fCounter->Fill(2);    

  // trigger
  fCtrue = -1;
  TString trigger = fAOD->GetFiredTriggerClasses();
  if (trigger.Contains("CTRUE-B")) fCtrue = 1;
  if (trigger.Contains("CTRUE-A")) fCtrue = 2;
  if (trigger.Contains("CTRUE-C")) fCtrue = 3;
  if (trigger.Contains("CTRUE-E")) fCtrue = 4;
  fC1zed = -1;
  if (trigger.Contains("C1ZED-B")) fC1zed = 1;
  if (trigger.Contains("C1ZED-A")) fC1zed = 2;
  if (trigger.Contains("C1ZED-C")) fC1zed = 3;
  if (trigger.Contains("C1ZED-E")) fC1zed = 4;
  if (fCtrue == -1 && fC1zed == -1) {PostData(2, fOutputList); return;}
  fCounter->Fill(3);    

  // tracks
  // if (fAOD->GetNumberOfTracks() > 0) {PostData(2, fOutputList); return;}
  fCounter->Fill(4);    

  // event info
  fRunNum = fAOD ->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();
  
  // trigger inputs
  fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();
  
  //Past-future protection maps
  fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();
  
  //ZDC
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {PostData(2, fOutputList); return;}
  fCounter->Fill(5);    

  fZem1Energy = dataZDC->GetZEM1Energy();
  fZem2Energy = dataZDC->GetZEM2Energy();  
  fZNAEnergy = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy = dataZDC->GetZPCTowerEnergy()[0];

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);
  
  // V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {PostData(2, fOutputList); return;}
  fCounter->Fill(6);    

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(!dataAD) {PostData(2, fOutputList); return;}
  fCounter->Fill(7);    

  fADADecision = dataAD->GetADADecision();
  fADCDecision = dataAD->GetADCDecision();

  // fill the tree
  fAnaTree->Fill();
  
  PostData(1, fAnaTree);
  PostData(2, fOutputList);          
}

//_____________________________________________________________________________
void AliAnalysisTaskTransTask::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
