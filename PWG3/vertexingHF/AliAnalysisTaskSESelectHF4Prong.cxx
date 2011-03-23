/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the selection of heavy flavor
// decay candidates and creation a stand-alone AOD for 
// 4prong D0 decay.
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
//         F.Colamaria, fabio.colamaria@ba.infn.it
/////////////////////////////////////////////////////////////


#include "Riostream.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TBits.h" 
#include "TNtuple.h"

#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSESelectHF4Prong.h"
#include "AliAODPidHF.h"
#include "AliRDHFCuts.h"

ClassImp(AliAnalysisTaskSESelectHF4Prong)


//________________________________________________________________________
AliAnalysisTaskSESelectHF4Prong::AliAnalysisTaskSESelectHF4Prong():
AliAnalysisTaskSE(),
fVerticesHFTClArr(0),
fCharm4ProngTClArr(0),
fSelected(0),
fOutput(0),
fOutput2(0),
fOutput3(0),
fOutput4(0),
fOutput5(0),
fOutputC(0),
fhInvMassD0Sum10MevBin1(0),
fhInvMassD0barSum10MevBin1(0),
fhInvMassSumAll10MevBin1(0),
fhInvMassD0Sum5MevBin1(0),
fhInvMassD0barSum5MevBin1(0),
fhInvMassSumAll5MevBin1(0),
fhInvMassD0Sum10MevBin2(0),
fhInvMassD0barSum10MevBin2(0),
fhInvMassSumAll10MevBin2(0),
fhInvMassD0Sum5MevBin2(0),
fhInvMassD0barSum5MevBin2(0),
fhInvMassSumAll5MevBin2(0),
fhInvMassD0Sum10MevBin3(0),
fhInvMassD0barSum10MevBin3(0),
fhInvMassSumAll10MevBin3(0),
fhInvMassD0Sum5MevBin3(0),
fhInvMassD0barSum5MevBin3(0),
fhInvMassSumAll5MevBin3(0),
fhInvMassD0Sum10MevBin4(0),
fhInvMassD0barSum10MevBin4(0),
fhInvMassSumAll10MevBin4(0),
fhInvMassD0Sum5MevBin4(0),
fhInvMassD0barSum5MevBin4(0),
fhInvMassSumAll5MevBin4(0),
fhInvMassD0Sum10MevBin5(0),
fhInvMassD0barSum10MevBin5(0),
fhInvMassSumAll10MevBin5(0),
fhInvMassD0Sum5MevBin5(0),
fhInvMassD0barSum5MevBin5(0),
fhInvMassSumAll5MevBin5(0),
fhInvMassMultipleOnlyBin1(0),
fhInvMassMultipleOnlyBin2(0),
fhInvMassMultipleOnlyBin3(0),
fhInvMassMultipleOnlyBin4(0),
fhInvMassMultipleOnlyBin5(0),
fScatterP4PID(0),
fPtVsY(0),
fPtVsYAll(0),
fEventCounter(0),
fCutDCA(0),
fCutDCA3(0),
fCutDCA2(0),
fCutDCA5(0),
fCutVertexDist2(0),
fCutVertexDist3(0),
fCutVertexDist4(0),
fCutCosinePoint(0),
fCutPt(0),
fCutY(0),
fPIDSel(0),
fPIDSelBin1(0),
fPIDSelBin2(0),
fPIDSelBin3(0),
fPIDSelBin4(0),
fPIDSelBin5(0),
fMultipleHyps(0),
fMultipleHypsType(0),
fPtSel(0),
fCuts(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSESelectHF4Prong::AliAnalysisTaskSESelectHF4Prong(const char *name,AliRDHFCutsD0toKpipipi* cuts):
AliAnalysisTaskSE(name),
fVerticesHFTClArr(0),
fCharm4ProngTClArr(0),
fSelected(0),
fOutput(0),
fOutput2(0),
fOutput3(0),
fOutput4(0),
fOutput5(0),
fOutputC(0),
fhInvMassD0Sum10MevBin1(0),
fhInvMassD0barSum10MevBin1(0),
fhInvMassSumAll10MevBin1(0),
fhInvMassD0Sum5MevBin1(0),
fhInvMassD0barSum5MevBin1(0),
fhInvMassSumAll5MevBin1(0),
fhInvMassD0Sum10MevBin2(0),
fhInvMassD0barSum10MevBin2(0),
fhInvMassSumAll10MevBin2(0),
fhInvMassD0Sum5MevBin2(0),
fhInvMassD0barSum5MevBin2(0),
fhInvMassSumAll5MevBin2(0),
fhInvMassD0Sum10MevBin3(0),
fhInvMassD0barSum10MevBin3(0),
fhInvMassSumAll10MevBin3(0),
fhInvMassD0Sum5MevBin3(0),
fhInvMassD0barSum5MevBin3(0),
fhInvMassSumAll5MevBin3(0),
fhInvMassD0Sum10MevBin4(0),
fhInvMassD0barSum10MevBin4(0),
fhInvMassSumAll10MevBin4(0),
fhInvMassD0Sum5MevBin4(0),
fhInvMassD0barSum5MevBin4(0),
fhInvMassSumAll5MevBin4(0),
fhInvMassD0Sum10MevBin5(0),
fhInvMassD0barSum10MevBin5(0),
fhInvMassSumAll10MevBin5(0),
fhInvMassD0Sum5MevBin5(0),
fhInvMassD0barSum5MevBin5(0),
fhInvMassSumAll5MevBin5(0),
fhInvMassMultipleOnlyBin1(0),
fhInvMassMultipleOnlyBin2(0),
fhInvMassMultipleOnlyBin3(0),
fhInvMassMultipleOnlyBin4(0),
fhInvMassMultipleOnlyBin5(0),
fScatterP4PID(0),
fPtVsY(0),
fPtVsYAll(0),
fEventCounter(0),
fCutDCA(0),
fCutDCA3(0),
fCutDCA2(0),
fCutDCA5(0),
fCutVertexDist2(0),
fCutVertexDist3(0),
fCutVertexDist4(0),
fCutCosinePoint(0),
fCutPt(0),
fCutY(0),
fPIDSel(0),
fPIDSelBin1(0),
fPIDSelBin2(0),
fPIDSelBin3(0),
fPIDSelBin4(0),
fPIDSelBin5(0),
fMultipleHyps(0),
fMultipleHypsType(0),
fPtSel(0),
fCuts(0)
{
  // Standard constructor
   
  fCuts=cuts;

  // Input slot #0 works with an Ntuple
//  DefineInput(0, TTree::Class());

  // Output slot #0 writes into a TTree container
  // Output slots #1-6 writes into a TList container
  DefineOutput(0, TTree::Class()); //default
  DefineOutput(1, TList::Class()); //histos inv. mass bin1
  DefineOutput(2, TList::Class()); //histos inv. mass bin2
  DefineOutput(3, TList::Class()); //histos inv. mass bin3
  DefineOutput(4, TList::Class()); //histos inv. mass bin4
  DefineOutput(5, TList::Class()); //histos inv. mass bin5
  DefineOutput(6, TList::Class()); //histos of cuts
  DefineOutput(7, AliRDHFCutsD0toKpipipi::Class()); //cuts
}

//________________________________________________________________________
AliAnalysisTaskSESelectHF4Prong::~AliAnalysisTaskSESelectHF4Prong()
{
  // Destructor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  } 
  if (fOutput2) {
    delete fOutput2;
    fOutput2 = 0;
  }
  if (fOutput3) {
    delete fOutput3;
    fOutput3 = 0;
  } 
  if (fOutput4) {
    delete fOutput4;
    fOutput4 = 0;
  }
  if (fOutput5) {
    delete fOutput5;
    fOutput5 = 0;
  }
  if (fOutputC) {
    delete fOutputC;
    fOutputC = 0;
  }   
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }

}  

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSESelectHF4Prong::Init() \n");

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSESelectHF4Prong::UserCreateOutputObjects() \n");

  fVerticesHFTClArr = new TClonesArray("AliAODVertex", 0);
  fVerticesHFTClArr->SetName("VerticesHF");
  AddAODBranch("TClonesArray", &fVerticesHFTClArr);

  fCharm4ProngTClArr = new TClonesArray("AliAODRecoDecayHF4Prong", 0);
  fCharm4ProngTClArr->SetName("Charm4Prong");
  AddAODBranch("TClonesArray", &fCharm4ProngTClArr);

  fOutput = new TList();
  fOutput->SetOwner();

  fOutput2 = new TList();
  fOutput2->SetOwner();

  fOutput3 = new TList();
  fOutput3->SetOwner();

  fOutput4 = new TList();
  fOutput4->SetOwner();

  fOutput5 = new TList();
  fOutput5->SetOwner();

  fOutputC = new TList();
  fOutputC->SetOwner();

  fhInvMassD0Sum10MevBin1 = new TH1F("fhInvMassD0Sum10MevBin1", "D0 invariant mass Bin1 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0Sum10MevBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum10MevBin1->SetMinimum(0);
  fOutput->Add(fhInvMassD0Sum10MevBin1);

  fhInvMassD0barSum10MevBin1 = new TH1F("fhInvMassD0barSum10MevBin1", "D0bar invariant mass Bin1 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0barSum10MevBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum10MevBin1->SetMinimum(0);
  fOutput->Add(fhInvMassD0barSum10MevBin1);

  fhInvMassSumAll10MevBin1 = new TH1F("fhInvMassSumAll10MevBin1", "D0/D0bar invariant mass Bin1 (good hyps); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassSumAll10MevBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll10MevBin1->SetMinimum(0);
  fOutput->Add(fhInvMassSumAll10MevBin1);

  fhInvMassD0Sum5MevBin1 = new TH1F("fhInvMassD0Sum5MevBin1", "D0 invariant mass Bin1 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0Sum5MevBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum5MevBin1->SetMinimum(0);
  fOutput->Add(fhInvMassD0Sum5MevBin1);

  fhInvMassD0barSum5MevBin1 = new TH1F("fhInvMassD0barSum5MevBin1", "D0bar invariant mass Bin1 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0barSum5MevBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum5MevBin1->SetMinimum(0);
  fOutput->Add(fhInvMassD0barSum5MevBin1);

  fhInvMassSumAll5MevBin1 = new TH1F("fhInvMassSumAll5MevBin1", "D0/D0bar invariant mass Bin1 (good hyps); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassSumAll5MevBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll5MevBin1->SetMinimum(0);
  fOutput->Add(fhInvMassSumAll5MevBin1);

  fhInvMassD0Sum10MevBin2 = new TH1F("fhInvMassD0Sum10MevBin2", "D0 invariant mass Bin2 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0Sum10MevBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum10MevBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassD0Sum10MevBin2);

  fhInvMassD0barSum10MevBin2 = new TH1F("fhInvMassD0barSum10MevBin2", "D0bar invariant mass Bin2 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0barSum10MevBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum10MevBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassD0barSum10MevBin2);

  fhInvMassSumAll10MevBin2 = new TH1F("fhInvMassSumAll10MevBin2", "D0/D0bar invariant mass Bin2 (good hyps); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassSumAll10MevBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll10MevBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassSumAll10MevBin2);

  fhInvMassD0Sum5MevBin2 = new TH1F("fhInvMassD0Sum5MevBin2", "D0 invariant mass Bin2 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0Sum5MevBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum5MevBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassD0Sum5MevBin2);

  fhInvMassD0barSum5MevBin2 = new TH1F("fhInvMassD0barSum5MevBin2", "D0bar invariant mass Bin2 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0barSum5MevBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum5MevBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassD0barSum5MevBin2);

  fhInvMassSumAll5MevBin2 = new TH1F("fhInvMassSumAll5MevBin2", "D0/D0bar invariant mass Bin2 (good hyps); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassSumAll5MevBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll5MevBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassSumAll5MevBin2);

  fhInvMassD0Sum10MevBin3 = new TH1F("fhInvMassD0Sum10MevBin3", "D0 invariant mass Bin3 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0Sum10MevBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum10MevBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassD0Sum10MevBin3);

  fhInvMassD0barSum10MevBin3 = new TH1F("fhInvMassD0barSum10MevBin3", "D0bar invariant mass Bin3 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0barSum10MevBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum10MevBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassD0barSum10MevBin3);

  fhInvMassSumAll10MevBin3 = new TH1F("fhInvMassSumAll10MevBin3", "D0/D0bar invariant mass Bin3 (good hyps); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassSumAll10MevBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll10MevBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassSumAll10MevBin3);

  fhInvMassD0Sum5MevBin3 = new TH1F("fhInvMassD0Sum5MevBin3", "D0 invariant mass Bin3 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0Sum5MevBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum5MevBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassD0Sum5MevBin3);

  fhInvMassD0barSum5MevBin3 = new TH1F("fhInvMassD0barSum5MevBin3", "D0bar invariant mass Bin3 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0barSum5MevBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum5MevBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassD0barSum5MevBin3);

  fhInvMassSumAll5MevBin3 = new TH1F("fhInvMassSumAll5MevBin3", "D0/D0bar invariant mass Bin3 (good hyps); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassSumAll5MevBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll5MevBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassSumAll5MevBin3);

  fhInvMassD0Sum10MevBin4 = new TH1F("fhInvMassD0Sum10MevBin4", "D0 invariant mass Bin4 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0Sum10MevBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum10MevBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassD0Sum10MevBin4);

  fhInvMassD0barSum10MevBin4 = new TH1F("fhInvMassD0barSum10MevBin4", "D0bar invariant mass Bin4 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0barSum10MevBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum10MevBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassD0barSum10MevBin4);

  fhInvMassSumAll10MevBin4 = new TH1F("fhInvMassSumAll10MevBin4", "D0/D0bar invariant mass Bin4 (good hyps); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassSumAll10MevBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll10MevBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassSumAll10MevBin4);

  fhInvMassD0Sum5MevBin4 = new TH1F("fhInvMassD0Sum5MevBin4", "D0 invariant mass Bin4 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0Sum5MevBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum5MevBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassD0Sum5MevBin4);

  fhInvMassD0barSum5MevBin4 = new TH1F("fhInvMassD0barSum5MevBin4", "D0bar invariant mass Bin4 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0barSum5MevBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum5MevBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassD0barSum5MevBin4);

  fhInvMassSumAll5MevBin4 = new TH1F("fhInvMassSumAll5MevBin4", "D0/D0bar invariant mass Bin4 (good hyps); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassSumAll5MevBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll5MevBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassSumAll5MevBin4);

  fhInvMassD0Sum10MevBin5 = new TH1F("fhInvMassD0Sum10MevBin5", "D0 invariant mass Bin5 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0Sum10MevBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum10MevBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassD0Sum10MevBin5);

  fhInvMassD0barSum10MevBin5 = new TH1F("fhInvMassD0barSum10MevBin5", "D0bar invariant mass Bin5 (good hyp); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassD0barSum10MevBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum10MevBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassD0barSum10MevBin5);

  fhInvMassSumAll10MevBin5 = new TH1F("fhInvMassSumAll10MevBin5", "D0/D0bar invariant mass Bin5 (good hyps); Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2);
  fhInvMassSumAll10MevBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll10MevBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassSumAll10MevBin5);

  fhInvMassD0Sum5MevBin5 = new TH1F("fhInvMassD0Sum5MevBin5", "D0 invariant mass Bin5 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0Sum5MevBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0Sum5MevBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassD0Sum5MevBin5);

  fhInvMassD0barSum5MevBin5 = new TH1F("fhInvMassD0barSum5MevBin5", "D0bar invariant mass Bin5 (good hyp); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassD0barSum5MevBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassD0barSum5MevBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassD0barSum5MevBin5);

  fhInvMassSumAll5MevBin5 = new TH1F("fhInvMassSumAll5MevBin5", "D0/D0bar invariant mass Bin5 (good hyps); Inv. mass [GeV]; Entries/5 MeV",120,1.6,2.2);
  fhInvMassSumAll5MevBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassSumAll5MevBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassSumAll5MevBin5);

  fhInvMassMultipleOnlyBin1 = new TH1F("fhInvMassMultipleOnlyBin1", "D0/D0bar invariant mass Bin1 (good hyps) - Multple hyps accepted only; Inv. mass [GeV]; Entries/10 MeV",120,1.6,2.2);
  fhInvMassMultipleOnlyBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassMultipleOnlyBin1->SetMinimum(0);
  fOutput->Add(fhInvMassMultipleOnlyBin1);

  fhInvMassMultipleOnlyBin2 = new TH1F("fhInvMassMultipleOnlyBin2", "D0/D0bar invariant mass Bin2 (good hyps) - Multple hyps accepted only; Inv. mass [GeV]; Entries/10 MeV",120,1.6,2.2);
  fhInvMassMultipleOnlyBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassMultipleOnlyBin2->SetMinimum(0);
  fOutput2->Add(fhInvMassMultipleOnlyBin2);

  fhInvMassMultipleOnlyBin3 = new TH1F("fhInvMassMultipleOnlyBin3", "D0/D0bar invariant mass Bin3 (good hyps) - Multple hyps accepted only; Inv. mass [GeV]; Entries/10 MeV",120,1.6,2.2);
  fhInvMassMultipleOnlyBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassMultipleOnlyBin3->SetMinimum(0);
  fOutput3->Add(fhInvMassMultipleOnlyBin3);

  fhInvMassMultipleOnlyBin4 = new TH1F("fhInvMassMultipleOnlyBin4", "D0/D0bar invariant mass Bin4 (good hyps) - Multple hyps accepted only; Inv. mass [GeV]; Entries/10 MeV",120,1.6,2.2);
  fhInvMassMultipleOnlyBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassMultipleOnlyBin4->SetMinimum(0);
  fOutput4->Add(fhInvMassMultipleOnlyBin4);

  fhInvMassMultipleOnlyBin5 = new TH1F("fhInvMassMultipleOnlyBin5", "D0/D0bar invariant mass Bin5 (good hyps) - Multple hyps accepted only; Inv. mass [GeV]; Entries/10 MeV",120,1.6,2.2);
  fhInvMassMultipleOnlyBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhInvMassMultipleOnlyBin5->SetMinimum(0);
  fOutput5->Add(fhInvMassMultipleOnlyBin5);

  fScatterP4PID = new TH2F("fScatterP4PID", "Transverse momentum of K vs l-s Pi (D0 + D0bar); Pt of K [GeV/c]; Pt of Pi [GeV/c]",500,0.,5.,500,0.,5.);
  fScatterP4PID->SetMinimum(0);
  fOutput->Add(fScatterP4PID);

  fPtVsY = new TH2F("fPtVsY", "Pt vs Y PPR Sel. Candidates; Pt [GeV/c]; Y",250,0.,25.,300,-3.,3.);
  fPtVsY->SetMinimum(0);
  fOutputC->Add(fPtVsY);

  fPtVsYAll = new TH2F("fPtVsYAll", "Pt vs Y All Candidates; Pt [GeV/c]; Y",250,0.,25.,300,-3.,3.);
  fPtVsYAll->SetMinimum(0);
  fOutputC->Add(fPtVsYAll);

  fEventCounter = new TH1F("fEventCounter", "N° of total events; NA; Events",1,0.,1.);
  fEventCounter->SetMinimum(0);
  fOutputC->Add(fEventCounter);

  fCutDCA = new TH1F("fCutDCA", "DCA of candidate (couple); DCA [cm]; Entries/micron",500,0.,0.05);
  fCutDCA->SetMinimum(0);
  fOutputC->Add(fCutDCA);

  fCutDCA3 = new TH1F("fCutDCA3", "DCA of candidate (trips); DCA [cm]; Entries/micron",500,0.,0.05);
  fCutDCA3->SetMinimum(0);
  fOutputC->Add(fCutDCA3);

  fCutDCA2 = new TH1F("fCutDCA2", "DCA of candidate (quads1); DCA [cm]; Entries/micron",500,0.,0.05);
  fCutDCA2->SetMinimum(0);
  fOutputC->Add(fCutDCA2);

  fCutDCA5 = new TH1F("fCutDCA5", "DCA of candidate (quads2); DCA [cm]; Entries/micron",500,0.,0.05);
  fCutDCA5->SetMinimum(0);
  fOutputC->Add(fCutDCA5);

  fCutVertexDist2 = new TH1F("fCutVertexDist2", "Distance Vtx doubl.-Primary Vtx; Distance [cm]; Entries/15 micron",500,0.,0.75);
  fCutVertexDist2->SetMinimum(0);
  fOutputC->Add(fCutVertexDist2);

  fCutVertexDist3 = new TH1F("fCutVertexDist3", "Distance Vtx trips-Primary Vtx; Distance [cm]; Entries/10 micron",500,0.,0.5);
  fCutVertexDist3->SetMinimum(0);
  fOutputC->Add(fCutVertexDist3);

  fCutVertexDist4 = new TH1F("fCutVertexDist4", "Distance Vtx quads-Primary Vtx; Distance [cm]; Entries/5 micron",500,0.,0.25);
  fCutVertexDist4->SetMinimum(0);
  fOutputC->Add(fCutVertexDist4);

  fCutCosinePoint = new TH1F("fCutCosinePoint", "Cosine of angle of pointing; Cos(Thetapt.); Entries/10^(-3)",250,0.75,1.);
  fCutCosinePoint->SetMinimum(0);
  fOutputC->Add(fCutCosinePoint);

  fCutPt = new TH1F("fCutPt", "Pt of candidate D0; Pt [GeV/c]; Entries/5 MeV",3000,0.,15.);
  fCutPt->SetMinimum(0);
  fOutputC->Add(fCutPt);

  fCutY = new TH1F("fCutY", "Y of candidate D0; Pt [GeV/c]; Entries/5 MeV",900,-9.,9.);
  fCutY->SetMinimum(0);
  fOutputC->Add(fCutY);

  fPIDSel = new TH1F("fPIDSel", "Ratio of D0 selected by PID for Correction",3,0.,3.);
  fPIDSel->SetMinimum(0);
  fPIDSel->GetXaxis()->SetBinLabel(1,"D0allhyp All");
  fPIDSel->GetXaxis()->SetBinLabel(2,"D0allhyp PID");
  fPIDSel->GetXaxis()->SetBinLabel(3,"D0allhyp PID (hypok)");

  fPIDSelBin1 = new TH1F("fPIDSelBin1", "Ratio of D0 selected by PID for Correction",3,0.,3.);
  fPIDSelBin1->SetMinimum(0);
  fPIDSelBin1->GetXaxis()->SetBinLabel(1,"D0allhyp All");
  fPIDSelBin1->GetXaxis()->SetBinLabel(2,"D0allhyp PID");
  fPIDSelBin1->GetXaxis()->SetBinLabel(3,"D0allhyp PID (hypok)");

  fPIDSelBin2 = new TH1F("fPIDSelBin1", "Ratio of D0 selected by PID for Correction",3,0.,3.);
  fPIDSelBin2->SetMinimum(0);
  fPIDSelBin2->GetXaxis()->SetBinLabel(1,"D0allhyp All");
  fPIDSelBin2->GetXaxis()->SetBinLabel(2,"D0allhyp PID");
  fPIDSelBin2->GetXaxis()->SetBinLabel(3,"D0allhyp PID (hypok)");

  fPIDSelBin3 = new TH1F("fPIDSelBin1", "Ratio of D0 selected by PID for Correction",3,0.,3.);
  fPIDSelBin3->SetMinimum(0);
  fPIDSelBin3->GetXaxis()->SetBinLabel(1,"D0allhyp All");
  fPIDSelBin3->GetXaxis()->SetBinLabel(2,"D0allhyp PID");
  fPIDSelBin3->GetXaxis()->SetBinLabel(3,"D0allhyp PID (hypok)");

  fPIDSelBin4 = new TH1F("fPIDSelBin1", "Ratio of D0 selected by PID for Correction",3,0.,3.);
  fPIDSelBin4->SetMinimum(0);
  fPIDSelBin4->GetXaxis()->SetBinLabel(1,"D0allhyp All");
  fPIDSelBin4->GetXaxis()->SetBinLabel(2,"D0allhyp PID");
  fPIDSelBin4->GetXaxis()->SetBinLabel(3,"D0allhyp PID (hypok)");

  fPIDSelBin5 = new TH1F("fPIDSelBin1", "Ratio of D0 selected by PID for Correction",3,0.,3.);
  fPIDSelBin5->SetMinimum(0);
  fPIDSelBin5->GetXaxis()->SetBinLabel(1,"D0allhyp All");
  fPIDSelBin5->GetXaxis()->SetBinLabel(2,"D0allhyp PID");
  fPIDSelBin5->GetXaxis()->SetBinLabel(3,"D0allhyp PID (hypok)");

  fMultipleHyps = new TH1F("fMultipleHyps", "N. of hyp. accepted for each candidate (accounted N. times)",8,0.,8.);
  fMultipleHyps->SetMinimum(0);
  fMultipleHyps->GetXaxis()->SetBinLabel(1,"1 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(2,"2 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(3,"3 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(4,"4 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(5,"1 (PPR+PID)");
  fMultipleHyps->GetXaxis()->SetBinLabel(6,"2 (PPR+PID)");
  fMultipleHyps->GetXaxis()->SetBinLabel(7,"3 (PPR+PID)");
  fMultipleHyps->GetXaxis()->SetBinLabel(8,"4 (PPR+PID)");

  fMultipleHypsType = new TH1F("fMultipleHypsType", "Type of hyp. accepted for each candidate",8,0.,8.);
  fMultipleHypsType->SetMinimum(0);
  fMultipleHypsType->GetXaxis()->SetBinLabel(1,"D0");
  fMultipleHypsType->GetXaxis()->SetBinLabel(2,"D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(3,"2D0");
  fMultipleHypsType->GetXaxis()->SetBinLabel(4,"2D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(5,"D0+D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(6,"2D0+D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(7,"D0+2D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(8,"2D0+2D0bar");

  fOutputC->Add(fMultipleHyps);
  fOutputC->Add(fMultipleHypsType);

  fOutputC->Add(fPIDSel);
  fOutput->Add(fPIDSelBin1);
  fOutput2->Add(fPIDSelBin2);
  fOutput3->Add(fPIDSelBin3);
  fOutput4->Add(fPIDSelBin4);
  fOutput5->Add(fPIDSelBin5);

  fPtSel = new TH1F("fPtSel", "Pt of candidates accepted; Pt [GeV/c]; Entries/10 MeV",2000,0.,20.);
  fPtSel->SetMinimum(0);
  fOutputC->Add(fPtSel);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates selection and histograms
  
  AliAODEvent *aodIn = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *inputArrayCharm4Prong = 0;

  if(!aodIn && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodIn = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      // load D0 candidates                                                   
      inputArrayCharm4Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
    }
  } else {
    // load D0 candidates                                                   
    inputArrayCharm4Prong=(TClonesArray*)aodIn->GetList()->FindObject("Charm4Prong");
  }

  if(!inputArrayCharm4Prong) {
    printf("AliAnalysisTaskSESelectHF4Prong::UserExec: D0to3Kpi branch not found!\n");
    return;
  }

  //print event info
//  aodIn->GetHeader()->Print();
  
  //Event counter ++
  fEventCounter->Fill(0);

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodIn->GetPrimaryVertex()) return;

  // primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodIn->GetPrimaryVertex();
//  vtx1->Print();

  // make trkIDtoEntry register (temporary)
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<aodIn->GetNumberOfTracks();it++) {
    AliAODTrack *track = aodIn->GetTrack(it);
    trkIDtoEntry[track->GetID()]=it;
  }

  Int_t iOutVerticesHF=0,iOutCharm4Prong=0;
  fVerticesHFTClArr->Delete();
  iOutVerticesHF = fVerticesHFTClArr->GetEntriesFast();
  TClonesArray &verticesHFRef = *fVerticesHFTClArr;
  fCharm4ProngTClArr->Delete();
  iOutCharm4Prong = fCharm4ProngTClArr->GetEntriesFast();
  TClonesArray &aodCharm4ProngRef = *fCharm4ProngTClArr;

  // loop over D0->K3pi candidates
  Int_t nInCharm4Prong = inputArrayCharm4Prong->GetEntriesFast();
  printf("Number of D0->K3pi: %d\n",nInCharm4Prong);

  for (Int_t iCharm4Prong = 0; iCharm4Prong < nInCharm4Prong; iCharm4Prong++) {
    AliAODRecoDecayHF4Prong *dIn = (AliAODRecoDecayHF4Prong*)inputArrayCharm4Prong->UncheckedAt(iCharm4Prong);
    Bool_t unsetvtx=kFALSE;

    if(!dIn->GetOwnPrimaryVtx()) {
      dIn->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
   
   //fill histos of cuts
   Double_t dca = dIn->GetDCA();
   Double_t dist2 = dIn->GetDist12toPrim();
   Double_t dist3 = dIn->GetDist3toPrim();
   Double_t dist4 = dIn->GetDist4toPrim();
   Double_t cosine = dIn->CosPointingAngle();
   Double_t ptPart = dIn->Pt();
   Double_t yPart = dIn->YD0();
   Double_t dcatrip = dIn->GetDCA(3);
   Double_t dcaquad1 = dIn->GetDCA(2);
   Double_t dcaquad2 = dIn->GetDCA(5);

   Double_t ptBinH[6] = {2.5,4.5,6.0,8.0,12.0,25.0};  //bin i has pt between values i and i+1

   fCutDCA->Fill(dca);
   fCutDCA3->Fill(dcatrip);
   fCutDCA2->Fill(dcaquad1);
   fCutDCA5->Fill(dcaquad2);
   fCutVertexDist2->Fill(dist2);
   fCutVertexDist3->Fill(dist3);
   fCutVertexDist4->Fill(dist4);
   fCutCosinePoint->Fill(cosine);
   fCutPt->Fill(ptPart);
   fCutY->Fill(yPart);
   fPtVsYAll->Fill(ptPart,yPart);

   //flags initialization
   fSelected = fCuts->IsSelected(dIn,AliRDHFCuts::kCandidate);
   Int_t selD01 = fCuts->D01Selected(dIn,AliRDHFCuts::kCandidate);
   Int_t selD02 = fCuts->D02Selected(dIn,AliRDHFCuts::kCandidate);  
   Int_t selD0bar1 = fCuts->D0bar1Selected(dIn,AliRDHFCuts::kCandidate);
   Int_t selD0bar2 = fCuts->D0bar2Selected(dIn,AliRDHFCuts::kCandidate);
   Int_t flagAccLim = 1;

   //Limited Acceptance
   if(ptPart > 5.) {
     if (TMath::Abs(yPart) > 0.8) flagAccLim = 0;
   } 
   else {
     Double_t maxFiducialY = -0.2/15*ptPart*ptPart+1.9/15*ptPart+0.5; 
     Double_t minFiducialY = 0.2/15*ptPart*ptPart-1.9/15*ptPart-0.5;		
     if (yPart < minFiducialY || yPart > maxFiducialY) flagAccLim = 0;;
   }

   //number of CANDIDATES (regardless of hypotheses) passing PPR
   if(fSelected==1||fSelected==2||fSelected==3) {
      fPtVsY->Fill(ptPart,yPart);
      fPIDSel->Fill(0);
        if (ptPart >= ptBinH[0] && ptPart < ptBinH[1]) 	    fPIDSelBin1->Fill(0);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fPIDSelBin2->Fill(0);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fPIDSelBin3->Fill(0);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fPIDSelBin4->Fill(0);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fPIDSelBin5->Fill(0);
        }
      

   PostData(6,fOutputC);
   
   //selection
   if((fSelected==1||fSelected==2||fSelected==3) && flagAccLim == 1) {  
      // get daughter AOD tracks
      AliAODTrack *trk0 = (AliAODTrack*)dIn->GetDaughter(0);
      AliAODTrack *trk1 = (AliAODTrack*)dIn->GetDaughter(1);
      AliAODTrack *trk2 = (AliAODTrack*)dIn->GetDaughter(2);
      AliAODTrack *trk3 = (AliAODTrack*)dIn->GetDaughter(3);
      if(!trk0 || !trk1 || !trk2 || !trk3) {
	trk0=aodIn->GetTrack(trkIDtoEntry[dIn->GetProngID(0)]);
	trk1=aodIn->GetTrack(trkIDtoEntry[dIn->GetProngID(1)]);
	trk2=aodIn->GetTrack(trkIDtoEntry[dIn->GetProngID(2)]);
	trk3=aodIn->GetTrack(trkIDtoEntry[dIn->GetProngID(3)]);
      }
      printf("pt of positive track #1: %f\n",trk0->Pt());
      printf("pt of negative track #1: %f\n",trk1->Pt());
      printf("pt of positive track #2: %f\n",trk2->Pt());
      printf("pt of negative track #2: %f\n",trk3->Pt());
      
      dIn->InvMassD0(fmassD0);
      dIn->InvMassD0bar(fmassD0bar);

      //fill histos (combining selection form cuts & PID (with rho information))
      Int_t hypD01 = 0, hypD02 = 0, hypD0bar1 = 0, hypD0bar2 = 0;
      Int_t pid1 = 0, pid2 = 0, pidbar1 = 0, pidbar2 = 0;
      Int_t pidSelection = fCuts->IsSelectedFromPID(dIn, &pid1, &pid2, &pidbar1, &pidbar2);

   //number of CANDIDATES (regardless of hypotheses) passing PPR + PID - PAY ATTENTION: hypoth. for PID and PPR may not be the same!
   if (pidSelection > 0) {
        fPIDSel->Fill(1);
        if (ptPart >= ptBinH[0] && ptPart < ptBinH[1]) 	    fPIDSelBin1->Fill(1);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fPIDSelBin2->Fill(1);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fPIDSelBin3->Fill(1);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fPIDSelBin4->Fill(1);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fPIDSelBin5->Fill(1);
      }
   
	   //number of hypoteses accepted per candidate after PPR
	   if(selD01+selD02+selD0bar1+selD0bar2 == 1) fMultipleHyps->Fill(0);
	   if(selD01+selD02+selD0bar1+selD0bar2 == 2) fMultipleHyps->Fill(1);
	   if(selD01+selD02+selD0bar1+selD0bar2 == 3) fMultipleHyps->Fill(2);
	   if(selD01+selD02+selD0bar1+selD0bar2 == 4) fMultipleHyps->Fill(3);

   //combine PPD + PID cuts
   if (selD01 == 1 && pid1==1) hypD01 = 1;
   if (selD02 == 1 && pid2==1) hypD02 = 1;
   if (selD0bar1 == 1 && pidbar1==1) hypD0bar1 = 1;
   if (selD0bar2 == 1 && pidbar2==1) hypD0bar2 = 1;

   //number of CANDIDATES (regardless of hypotheses) passing PPR + PID - PAY ATTENTION: hypoth. for PID and PPR must match (at least one)!
   if (hypD01 == 1 || hypD02 == 1 || hypD0bar1 == 1 || hypD0bar2 == 1) {
        fPIDSel->Fill(2);
        if (ptPart >= ptBinH[0] && ptPart < ptBinH[1]) 	    fPIDSelBin1->Fill(2);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fPIDSelBin2->Fill(2);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fPIDSelBin3->Fill(2);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fPIDSelBin4->Fill(2);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fPIDSelBin5->Fill(2);
      }

	   //number of hypoteses accepted per candidate after PPR and PID
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 1) fMultipleHyps->Fill(4);
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 2) fMultipleHyps->Fill(5);
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 3) fMultipleHyps->Fill(6);
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 4) fMultipleHyps->Fill(7);

	   //type of hypoteses accepted per candidate after PPR and PID
	   if(hypD01+hypD02 == 1 && hypD0bar1+hypD0bar2 == 0) fMultipleHypsType->Fill(0);
	   if(hypD01+hypD02 == 0 && hypD0bar1+hypD0bar2 == 1) fMultipleHypsType->Fill(1);
	   if(hypD01+hypD02 == 2 && hypD0bar1+hypD0bar2 == 0) fMultipleHypsType->Fill(2);
	   if(hypD01+hypD02 == 0 && hypD0bar1+hypD0bar2 == 2) fMultipleHypsType->Fill(3);
	   if(hypD01+hypD02 == 1 && hypD0bar1+hypD0bar2 == 1) fMultipleHypsType->Fill(4);
	   if(hypD01+hypD02 == 2 && hypD0bar1+hypD0bar2 == 1) fMultipleHypsType->Fill(5);
	   if(hypD01+hypD02 == 1 && hypD0bar1+hypD0bar2 == 2) fMultipleHypsType->Fill(6);
	   if(hypD01+hypD02 == 2 && hypD0bar1+hypD0bar2 == 2) fMultipleHypsType->Fill(7);


//
// This section is auxiliary: (multiple hypotheses histos)
//
 if (hypD01+hypD02+hypD0bar1+hypD0bar2 == 2 || hypD01+hypD02+hypD0bar1+hypD0bar2 == 3 || hypD01+hypD02+hypD0bar1+hypD0bar2 == 4) 
 {
    
   if (ptPart > ptBinH[0]) {
      // D01 hyp.
      if(hypD01==1) {
        if (ptPart < ptBinH[1]) 	 		    fhInvMassMultipleOnlyBin1->Fill(fmassD0[0]);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fhInvMassMultipleOnlyBin2->Fill(fmassD0[0]);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fhInvMassMultipleOnlyBin3->Fill(fmassD0[0]);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fhInvMassMultipleOnlyBin4->Fill(fmassD0[0]);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fhInvMassMultipleOnlyBin5->Fill(fmassD0[0]);
       }
      // D02 hyp.
      if(hypD02==1) {
        if (ptPart < ptBinH[1]) 	 		    fhInvMassMultipleOnlyBin1->Fill(fmassD0[1]);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fhInvMassMultipleOnlyBin2->Fill(fmassD0[1]);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fhInvMassMultipleOnlyBin3->Fill(fmassD0[1]);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fhInvMassMultipleOnlyBin4->Fill(fmassD0[1]);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fhInvMassMultipleOnlyBin5->Fill(fmassD0[1]);
       }
      // D0bar1 hyp.
      if(hypD0bar1==1) {
	if (ptPart < ptBinH[1]) 	 	    	    fhInvMassMultipleOnlyBin1->Fill(fmassD0bar[0]);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fhInvMassMultipleOnlyBin2->Fill(fmassD0bar[0]);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fhInvMassMultipleOnlyBin3->Fill(fmassD0bar[0]);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fhInvMassMultipleOnlyBin4->Fill(fmassD0bar[0]);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fhInvMassMultipleOnlyBin5->Fill(fmassD0bar[0]); 
       } 
      // D0bar2 hyp.
      if(hypD0bar2==1) {
        if (ptPart < ptBinH[1]) 	 		    fhInvMassMultipleOnlyBin1->Fill(fmassD0bar[1]);
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) fhInvMassMultipleOnlyBin2->Fill(fmassD0bar[1]);
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) fhInvMassMultipleOnlyBin3->Fill(fmassD0bar[1]);
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) fhInvMassMultipleOnlyBin4->Fill(fmassD0bar[1]);
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5]) fhInvMassMultipleOnlyBin5->Fill(fmassD0bar[1]);
       }
     }

   }
// 
//end of auxiliary section
//


   //All histos are filled if Pt of candidate is greater than minimum of first bin (in this way: bin1+bin2+...binN = whole)
   if (ptPart > ptBinH[0]) {
    
      // D01 hyp.
      if(hypD01==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk1->Pt(),trk3->Pt());
        if (ptPart < ptBinH[1]) 	 		 {fhInvMassD0Sum10MevBin1->Fill(fmassD0[0]);
			         	 	  	  fhInvMassD0Sum5MevBin1->Fill(fmassD0[0]);
	                                 		  fhInvMassSumAll10MevBin1->Fill(fmassD0[0]);
        	                         		  fhInvMassSumAll5MevBin1->Fill(fmassD0[0]);}
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) {fhInvMassD0Sum10MevBin2->Fill(fmassD0[0]);
				      	 	  	     fhInvMassD0Sum5MevBin2->Fill(fmassD0[0]);
	                          	     		     fhInvMassSumAll10MevBin2->Fill(fmassD0[0]);
        	              	            		     fhInvMassSumAll5MevBin2->Fill(fmassD0[0]);}
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) {fhInvMassD0Sum10MevBin3->Fill(fmassD0[0]);
			   	       	  	  	     fhInvMassD0Sum5MevBin3->Fill(fmassD0[0]);
                          		      	  	     fhInvMassSumAll10MevBin3->Fill(fmassD0[0]);
                          	      		  	     fhInvMassSumAll5MevBin3->Fill(fmassD0[0]);}
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) {fhInvMassD0Sum10MevBin4->Fill(fmassD0[0]);
			   	       	  	  	     fhInvMassD0Sum5MevBin4->Fill(fmassD0[0]);
                          		      	  	     fhInvMassSumAll10MevBin4->Fill(fmassD0[0]);
                          	      		  	     fhInvMassSumAll5MevBin4->Fill(fmassD0[0]);}
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5])	{fhInvMassD0Sum10MevBin5->Fill(fmassD0[0]);
				         			fhInvMassD0Sum5MevBin5->Fill(fmassD0[0]);
			                 			fhInvMassSumAll10MevBin5->Fill(fmassD0[0]);
			                 			fhInvMassSumAll5MevBin5->Fill(fmassD0[0]);} 
     }
      // D02 hyp.
      if(hypD02==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk3->Pt(),trk1->Pt());
        if (ptPart < ptBinH[1]) 	 		 {fhInvMassD0Sum10MevBin1->Fill(fmassD0[1]);
			         	 	  	  fhInvMassD0Sum5MevBin1->Fill(fmassD0[1]);
	                                 		  fhInvMassSumAll10MevBin1->Fill(fmassD0[1]);
        	                         		  fhInvMassSumAll5MevBin1->Fill(fmassD0[1]);}
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) {fhInvMassD0Sum10MevBin2->Fill(fmassD0[1]);
				      	 	  	     fhInvMassD0Sum5MevBin2->Fill(fmassD0[1]);
	                          	     		     fhInvMassSumAll10MevBin2->Fill(fmassD0[1]);
        	              	            		     fhInvMassSumAll5MevBin2->Fill(fmassD0[1]);}
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) {fhInvMassD0Sum10MevBin3->Fill(fmassD0[1]);
			   	       	  	  	     fhInvMassD0Sum5MevBin3->Fill(fmassD0[1]);
                          		      	  	     fhInvMassSumAll10MevBin3->Fill(fmassD0[1]);
                          	      		  	     fhInvMassSumAll5MevBin3->Fill(fmassD0[1]);}
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) {fhInvMassD0Sum10MevBin4->Fill(fmassD0[1]);
			   	       	  	  	     fhInvMassD0Sum5MevBin4->Fill(fmassD0[1]);
                          		      	  	     fhInvMassSumAll10MevBin4->Fill(fmassD0[1]);
                          	      		  	     fhInvMassSumAll5MevBin4->Fill(fmassD0[1]);}
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5])	{fhInvMassD0Sum10MevBin5->Fill(fmassD0[1]);
				        			 fhInvMassD0Sum5MevBin5->Fill(fmassD0[1]);
			                 			 fhInvMassSumAll10MevBin5->Fill(fmassD0[1]);
			                			 fhInvMassSumAll5MevBin5->Fill(fmassD0[1]);}
     }
      // D0bar1 hyp.
      if(hypD0bar1==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk0->Pt(),trk2->Pt());
        if (ptPart < ptBinH[1]) 	 		 {fhInvMassD0barSum10MevBin1->Fill(fmassD0bar[0]);
			         	 	  	  fhInvMassD0barSum5MevBin1->Fill(fmassD0bar[0]);
                             		    		  fhInvMassSumAll10MevBin1->Fill(fmassD0bar[0]);
                                 			  fhInvMassSumAll5MevBin1->Fill(fmassD0bar[0]);}
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) {fhInvMassD0barSum10MevBin2->Fill(fmassD0bar[0]);
				      	 	  	     fhInvMassD0barSum5MevBin2->Fill(fmassD0bar[0]);
                           	     	 	  	     fhInvMassSumAll10MevBin2->Fill(fmassD0bar[0]);
                       	            	 		     fhInvMassSumAll5MevBin2->Fill(fmassD0bar[0]);}
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) {fhInvMassD0barSum10MevBin3->Fill(fmassD0bar[0]);
			   	         	 	     fhInvMassD0barSum5MevBin3->Fill(fmassD0bar[0]);
                          	      	 		     fhInvMassSumAll10MevBin3->Fill(fmassD0bar[0]);
                          	      	 		     fhInvMassSumAll5MevBin3->Fill(fmassD0bar[0]);}
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) {fhInvMassD0barSum10MevBin4->Fill(fmassD0bar[0]);
			   	         	 	     fhInvMassD0barSum5MevBin4->Fill(fmassD0bar[0]);
                          	      	 		     fhInvMassSumAll10MevBin4->Fill(fmassD0bar[0]);
                          	      	 		     fhInvMassSumAll5MevBin4->Fill(fmassD0bar[0]);}
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5])	{fhInvMassD0barSum10MevBin5->Fill(fmassD0bar[0]);
				         			fhInvMassD0barSum5MevBin5->Fill(fmassD0bar[0]);
			                 			fhInvMassSumAll10MevBin5->Fill(fmassD0bar[0]);
			                 			fhInvMassSumAll5MevBin5->Fill(fmassD0bar[0]);}   
     } 
      // D0bar2 hyp.
      if(hypD0bar2==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk2->Pt(),trk0->Pt());
        if (ptPart < ptBinH[1]) 	 		 {fhInvMassD0barSum10MevBin1->Fill(fmassD0bar[1]);
			         	 	  	  fhInvMassD0barSum5MevBin1->Fill(fmassD0bar[1]);
                             		    		  fhInvMassSumAll10MevBin1->Fill(fmassD0bar[1]);
                                 			  fhInvMassSumAll5MevBin1->Fill(fmassD0bar[1]);}
        else if (ptPart >= ptBinH[1] && ptPart < ptBinH[2]) {fhInvMassD0barSum10MevBin2->Fill(fmassD0bar[1]);
				      	 	  	     fhInvMassD0barSum5MevBin2->Fill(fmassD0bar[1]);
                           	     	 	  	     fhInvMassSumAll10MevBin2->Fill(fmassD0bar[1]);
                       	            	 		     fhInvMassSumAll5MevBin2->Fill(fmassD0bar[1]);}
        else if (ptPart >= ptBinH[2] && ptPart < ptBinH[3]) {fhInvMassD0barSum10MevBin3->Fill(fmassD0bar[1]);
			   	         	 	     fhInvMassD0barSum5MevBin3->Fill(fmassD0bar[1]);
                          	      	 		     fhInvMassSumAll10MevBin3->Fill(fmassD0bar[1]);
                          	      	 		     fhInvMassSumAll5MevBin3->Fill(fmassD0bar[1]);}
        else if (ptPart >= ptBinH[3] && ptPart < ptBinH[4]) {fhInvMassD0barSum10MevBin4->Fill(fmassD0bar[1]);
			   	         	 	     fhInvMassD0barSum5MevBin4->Fill(fmassD0bar[1]);
                          	      	 		     fhInvMassSumAll10MevBin4->Fill(fmassD0bar[1]);
                          	      	 		     fhInvMassSumAll5MevBin4->Fill(fmassD0bar[1]);}
        else if (ptPart >= ptBinH[4] && ptPart < ptBinH[5])	{fhInvMassD0barSum10MevBin5->Fill(fmassD0bar[1]);
				         			fhInvMassD0barSum5MevBin5->Fill(fmassD0bar[1]);
			                 			fhInvMassSumAll10MevBin5->Fill(fmassD0bar[1]);
			                 			fhInvMassSumAll5MevBin5->Fill(fmassD0bar[1]);} 
     }
   }

      PostData(1,fOutput);
      PostData(2,fOutput2);
      PostData(3,fOutput3);
      PostData(4,fOutput4);
      PostData(5,fOutput5);

      // HERE ONE COULD RECALCULATE THE VERTEX USING THE KF PACKAGE

      // clone candidate for output AOD
      if(hypD01||hypD02||hypD0bar1||hypD0bar2) {
      AliAODVertex *v = new(verticesHFRef[iOutVerticesHF++]) 
	AliAODVertex(*(dIn->GetSecondaryVtx()));
      AliAODRecoDecayHF4Prong *dOut=new(aodCharm4ProngRef[iOutCharm4Prong++]) 
	AliAODRecoDecayHF4Prong(*dIn);
      dOut->SetSecondaryVtx(v);
      dOut->SetOwnPrimaryVtx((AliAODVertex*)((dIn->GetOwnPrimaryVtx())->Clone()));
      v->SetParent(dOut); 
      }
    } //end of selection loop (starts with fSelected == 1, 2 or 3)
    if(unsetvtx) dIn->UnsetOwnPrimaryVtx();
  } // end loop on D0->K3pi

  printf("Number of selected D0->K3pi: %d\n",iOutCharm4Prong);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSESelectHF4Prong: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fOutput2 = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutput2) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fOutput3 = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutput3) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fOutput4 = dynamic_cast<TList*> (GetOutputData(4));
  if (!fOutput4) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  if (!fOutput5) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fOutputC = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputC) {
    printf("ERROR: fOutputC not available\n");
    return;
  }

  fhInvMassD0Sum10MevBin1 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum10MevBin1"));
  fhInvMassD0barSum10MevBin1 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum10MevBin1"));
  fhInvMassSumAll10MevBin1 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll10MevBin1"));
  fhInvMassD0Sum5MevBin1 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum5MevBin1"));
  fhInvMassD0barSum5MevBin1 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum5MevBin1"));
  fhInvMassSumAll5MevBin1 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll5MevBin1"));
  fhInvMassD0Sum10MevBin2 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum10MevBin2"));
  fhInvMassD0barSum10MevBin2 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum10MevBin2"));
  fhInvMassSumAll10MevBin2 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll10MevBin2"));
  fhInvMassD0Sum5MevBin2 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum5MevBin2"));
  fhInvMassD0barSum5MevBin2 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum5MevBin2"));
  fhInvMassSumAll5MevBin2 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll5MevBin2"));
  fhInvMassD0Sum10MevBin3 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum10MevBin3"));
  fhInvMassD0barSum10MevBin3 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum10MevBin3"));
  fhInvMassSumAll10MevBin3 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll10MevBin3"));
  fhInvMassD0Sum5MevBin3 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum5MevBin3"));
  fhInvMassD0barSum5MevBin3 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum5MevBin3"));
  fhInvMassSumAll5MevBin3 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll5MevBin3"));
  fhInvMassD0Sum10MevBin4 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum10MevBin4"));
  fhInvMassD0barSum10MevBin4 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum10MevBin4"));
  fhInvMassSumAll10MevBin4 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll10MevBin4"));
  fhInvMassD0Sum5MevBin4 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum5MevBin4"));
  fhInvMassD0barSum5MevBin4 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum5MevBin4"));
  fhInvMassSumAll5MevBin4 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll5MevBin4"));
  fhInvMassD0Sum10MevBin5 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum10MevBin5"));
  fhInvMassD0barSum10MevBin5 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum10MevBin5"));
  fhInvMassSumAll10MevBin5 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll10MevBin5"));
  fhInvMassD0Sum5MevBin5 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0Sum5MevBin5"));
  fhInvMassD0barSum5MevBin5 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassD0barSum5MevBin5"));
  fhInvMassSumAll5MevBin5 = dynamic_cast<TH1F*>(fOutput->FindObject("fhInvMassSumAll5MevBin5"));

  fScatterP4PID = dynamic_cast<TH2F*>(fOutput->FindObject("fScatterP4PID"));
  fPtVsY = dynamic_cast<TH2F*>(fOutputC->FindObject("fPtVsY"));
  fPtVsYAll = dynamic_cast<TH2F*>(fOutputC->FindObject("fPtVsYAll"));

  fEventCounter = dynamic_cast<TH1F*>(fOutputC->FindObject("fEventCounter"));

  fCutDCA = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutDCA"));
  fCutDCA3 = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutDCA3"));
  fCutDCA2 = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutDCA2"));
  fCutDCA5 = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutDCA5"));
  fCutVertexDist2 = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutVertexDist2"));
  fCutVertexDist3 = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutVertexDist3"));
  fCutVertexDist4 = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutVertexDist4"));
  fCutCosinePoint = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutCosinePoint"));

  fCutPt = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutPt"));
  fCutY = dynamic_cast<TH1F*>(fOutputC->FindObject("fCutY"));
  fPIDSel = dynamic_cast<TH1F*>(fOutputC->FindObject("fPIDSel"));
  fPIDSelBin1 = dynamic_cast<TH1F*>(fOutputC->FindObject("fPIDSelBin1"));
  fPIDSelBin2 = dynamic_cast<TH1F*>(fOutputC->FindObject("fPIDSelBin2"));
  fPIDSelBin3 = dynamic_cast<TH1F*>(fOutputC->FindObject("fPIDSelBin3"));
  fPIDSelBin4 = dynamic_cast<TH1F*>(fOutputC->FindObject("fPIDSelBin4"));
  fPIDSelBin5 = dynamic_cast<TH1F*>(fOutputC->FindObject("fPIDSelBin5"));
}

