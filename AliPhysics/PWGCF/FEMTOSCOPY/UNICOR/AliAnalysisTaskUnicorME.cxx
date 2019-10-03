/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
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

//Author: Dariusz Miskowiec 2007

//=============================================================================
// unicor analysis task
//=============================================================================
#include "AliMultiEventInputHandler.h"
#include "AliESDEvent.h"
#include "AliUnicorAnalGlobal.h"
#include "AliUnicorAnalSingle.h"
#include "AliUnicorAnalCorrel.h"
#include "AliUnicorAnalPtfluc.h"
#include "AliUnicorAnalHighpt.h"
#include "AliUnicorEventAliceESD.h"
#include "AliAnalysisTaskUnicorME.h"

ClassImp(AliAnalysisTaskUnicorME)

//=============================================================================
AliAnalysisTaskUnicorME::AliAnalysisTaskUnicorME(const char *name) : 
  AliAnalysisTaskME(name), 
  fEv0(0),
  fEv1(0),
  fOutputList(0)
{
  // constructor

  fEv0 = new AliUnicorEventAliceESD();
  fEv1 = new AliUnicorEventAliceESD();
  DefineOutput(1, TList::Class());
}
//=============================================================================
void AliAnalysisTaskUnicorME::UserCreateOutputObjects() 
{
  // executed once on each worker 

  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->Add(new AliUnicorAnalGlobal("dag"));
  fOutputList->Add(new AliUnicorAnalSingle("all",fEv0->Etamin(),fEv0->Etamax(),0));
  fOutputList->Add(new AliUnicorAnalSingle("pim",fEv0->Etamin(),fEv0->Etamax(),-211));
  fOutputList->Add(new AliUnicorAnalSingle("pip",fEv0->Etamin(),fEv0->Etamax(), 211));
  AliUnicorAnalCorrel::AnalysisFrame frame = AliUnicorAnalCorrel::kLCMS;
  fOutputList->Add(new AliUnicorAnalCorrel("cnn",fEv0->Etamin(),fEv0->Etamax(),-211,-211, frame));
  fOutputList->Add(new AliUnicorAnalCorrel("cpp",fEv0->Etamin(),fEv0->Etamax(), 211, 211, frame));
  fOutputList->Add(new AliUnicorAnalCorrel("cnp",fEv0->Etamin(),fEv0->Etamax(),-211, 211, frame));
  fOutputList->Add(new AliUnicorAnalPtfluc("ptf",0,0));
  fOutputList->Add(new AliUnicorAnalHighpt("hpt",fEv0->Etamin(),fEv0->Etamax(),0,0));
}
//=============================================================================
void AliAnalysisTaskUnicorME::UserExec(Option_t */*option*/)
{
  // process one event pair

  static int nic0=0;
  static int nic1=0;
  static int nic2=0;
  static int nic3=0;

  nic0++;

  if (fInputHandler->GetBufferSize() < 2) return;
  nic1++;
  AliESDEvent *esd0 = dynamic_cast<AliESDEvent*>(GetEvent(0));
  AliESDEvent *esd1 = dynamic_cast<AliESDEvent*>(GetEvent(1));

  if (!esd0) return;
  nic2++;
  //  if (esd0->GetEventType() != 7) return; // physics event
  fEv0->SetESD(esd0);
  if (!fEv0->Good()) return;
  nic3++;
  ((AliUnicorAnalGlobal *) fOutputList->At(0))->Process(fEv0);
  ((AliUnicorAnalSingle *) fOutputList->At(1))->Process(fEv0);
  ((AliUnicorAnalSingle *) fOutputList->At(2))->Process(fEv0);
  ((AliUnicorAnalSingle *) fOutputList->At(3))->Process(fEv0);
  ((AliUnicorAnalCorrel *) fOutputList->At(4))->Process(0,fEv0,fEv0,0);
  ((AliUnicorAnalCorrel *) fOutputList->At(4))->Process(2,fEv0,fEv0,TMath::DegToRad()*180);
  ((AliUnicorAnalCorrel *) fOutputList->At(5))->Process(0,fEv0,fEv0,0);
  ((AliUnicorAnalCorrel *) fOutputList->At(5))->Process(2,fEv0,fEv0,TMath::DegToRad()*180);
  ((AliUnicorAnalCorrel *) fOutputList->At(6))->Process(0,fEv0,fEv0,0);
  ((AliUnicorAnalCorrel *) fOutputList->At(6))->Process(2,fEv0,fEv0,TMath::DegToRad()*180);
  ((AliUnicorAnalPtfluc *) fOutputList->At(7))->Process(0,fEv0,fEv0);
  ((AliUnicorAnalHighpt *) fOutputList->At(8))->Process(fEv0,fEv0);

  if (!esd1) return;
  //  if (esd1->GetEventType() != 7) return; // physics event
  printf("esd0 nr %3d mult %3d     esd1 nr %3d mult %3d\n",
	 esd0->GetEventNumberInFile(), esd0->GetNumberOfTracks(), 
	 esd1->GetEventNumberInFile(), esd1->GetNumberOfTracks());
  fEv1->SetESD(esd1);
  if (!fEv1->Good()) return;
  ((AliUnicorAnalCorrel *) fOutputList->At(4))->Process(1,fEv0,fEv1,0);
  ((AliUnicorAnalCorrel *) fOutputList->At(5))->Process(1,fEv0,fEv1,0);
  ((AliUnicorAnalCorrel *) fOutputList->At(6))->Process(1,fEv0,fEv1,0);
  ((AliUnicorAnalPtfluc *) fOutputList->At(7))->Process(1,fEv0,fEv1);
  ((AliUnicorAnalHighpt *) fOutputList->At(8))->Process(fEv0,fEv1);

  printf("counts: %6d %6d %6d %6d\n",nic0,nic1,nic2,nic3);
  PostData(1, fOutputList);
} 
//=============================================================================
void AliAnalysisTaskUnicorME::Terminate(Option_t */*option*/)
{
  // terminate

  printf("terminate \n");
  TList *outputlist = (TList*) GetOutputData(1);
  int n = outputlist->GetEntries();
  if (n) ((AliUnicorAnal *) outputlist->At(0))->Save("unicor-result.root","recreate");
  for (int i=1; i<n; i++) ((AliUnicorAnal *) outputlist->At(i))->Save("unicor-result.root");
}
//=============================================================================
