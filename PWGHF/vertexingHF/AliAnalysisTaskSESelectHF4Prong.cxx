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
fMCTruth(0),
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
fhReflBin1(0),
fhReflBin2(0),
fhReflBin3(0),
fhReflBin4(0),
fhReflBin5(0),
fhReflD0Bin1(0),
fhReflD0Bin2(0),
fhReflD0Bin3(0),
fhReflD0Bin4(0),
fhReflD0Bin5(0),
fhReflD0barBin1(0),
fhReflD0barBin2(0),
fhReflD0barBin3(0),
fhReflD0barBin4(0),
fhReflD0barBin5(0),
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

//_______________________________________________________________________
AliAnalysisTaskSESelectHF4Prong::AliAnalysisTaskSESelectHF4Prong(const char *name,AliRDHFCutsD0toKpipipi* cuts):
AliAnalysisTaskSE(name),
fVerticesHFTClArr(0),
fCharm4ProngTClArr(0),
fSelected(0),
fMCTruth(0),
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
fhReflBin1(0),
fhReflBin2(0),
fhReflBin3(0),
fhReflBin4(0),
fhReflBin5(0),
fhReflD0Bin1(0),
fhReflD0Bin2(0),
fhReflD0Bin3(0),
fhReflD0Bin4(0),
fhReflD0Bin5(0),
fhReflD0barBin1(0),
fhReflD0barBin2(0),
fhReflD0barBin3(0),
fhReflD0barBin4(0),
fhReflD0barBin5(0),
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
  PrintPtBinHandMCFlag();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSESelectHF4Prong::UserCreateOutputObjects() \n");

  if(fDebug > 1) PrintPtBinHandMCFlag();

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

  fhReflBin1 = new TH2F("fhReflBin1", "Invariant Mass Histogram for reflections; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,8,0.,8.);
  fhReflBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhReflBin1->SetMinimum(0);
  fOutput->Add(fhReflBin1);

  fhReflBin2 = new TH2F("fhReflBin2", "Invariant Mass Histogram for reflections; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,8,0.,8.);
  fhReflBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhReflBin2->SetMinimum(0);
  fOutput2->Add(fhReflBin2);

  fhReflBin3 = new TH2F("fhReflBin3", "Invariant Mass Histogram for reflections; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,8,0.,8.);
  fhReflBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhReflBin3->SetMinimum(0);
  fOutput3->Add(fhReflBin3);

  fhReflBin4 = new TH2F("fhReflBin4", "Invariant Mass Histogram for reflections; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,8,0.,8.);
  fhReflBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhReflBin4->SetMinimum(0);
  fOutput4->Add(fhReflBin4);

  fhReflBin5 = new TH2F("fhReflBin5", "Invariant Mass Histogram for reflections; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,8,0.,8.);
  fhReflBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhReflBin5->SetMinimum(0);
  fOutput5->Add(fhReflBin5);

  fhReflD0Bin1 = new TH2F("fhReflD0Bin1", "Invariant Mass Histogram for reflections - D0 hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0Bin1->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0Bin1->SetMinimum(0);
  fOutput->Add(fhReflD0Bin1);

  fhReflD0Bin2 = new TH2F("fhReflD0Bin2", "Invariant Mass Histogram for reflections - D0 hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0Bin2->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0Bin2->SetMinimum(0);
  fOutput2->Add(fhReflD0Bin2);

  fhReflD0Bin3 = new TH2F("fhReflD0Bin3", "Invariant Mass Histogram for reflections - D0 hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0Bin3->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0Bin3->SetMinimum(0);
  fOutput3->Add(fhReflD0Bin3);

  fhReflD0Bin4 = new TH2F("fhReflD0Bin4", "Invariant Mass Histogram for reflections - D0 hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0Bin4->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0Bin4->SetMinimum(0);
  fOutput4->Add(fhReflD0Bin4);

  fhReflD0Bin5 = new TH2F("fhReflD0Bin5", "Invariant Mass Histogram for reflections - D0 hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0Bin5->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0Bin5->SetMinimum(0);
  fOutput5->Add(fhReflD0Bin5);

  fhReflD0barBin1 = new TH2F("fhReflD0barBin1", "Invariant Mass Histogram for reflections - D0bar hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0barBin1->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0barBin1->SetMinimum(0);
  fOutput->Add(fhReflD0barBin1);

  fhReflD0barBin2 = new TH2F("fhReflD0barBin2", "Invariant Mass Histogram for reflections - D0bar hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0barBin2->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0barBin2->SetMinimum(0);
  fOutput2->Add(fhReflD0barBin2);

  fhReflD0barBin3 = new TH2F("fhReflD0barBin3", "Invariant Mass Histogram for reflections - D0bar hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0barBin3->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0barBin3->SetMinimum(0);
  fOutput3->Add(fhReflD0barBin3);

  fhReflD0barBin4 = new TH2F("fhReflD0barBin4", "Invariant Mass Histogram for reflections - D0bar hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0barBin4->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0barBin4->SetMinimum(0);
  fOutput4->Add(fhReflD0barBin4);

  fhReflD0barBin5 = new TH2F("fhReflD0barBin5", "Invariant Mass Histogram for reflections - D0bar hyp.; Inv. mass [GeV]; Entries/10 MeV",60,1.6,2.2,7,0.,7.);
  fhReflD0barBin5->Sumw2(); //Create structure to store sum of squares of weights
  fhReflD0barBin5->SetMinimum(0);
  fOutput5->Add(fhReflD0barBin5);

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

  fMultipleHyps = new TH2F("fMultipleHyps", "N. of hyp. accepted for each candidate (accounted N. times)",8,0.,8.,5,0.,5.);
  fMultipleHyps->SetMinimum(0);
  fMultipleHyps->GetXaxis()->SetBinLabel(1,"1 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(2,"2 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(3,"3 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(4,"4 (PPR)");
  fMultipleHyps->GetXaxis()->SetBinLabel(5,"1 (PPR+PID)");
  fMultipleHyps->GetXaxis()->SetBinLabel(6,"2 (PPR+PID)");
  fMultipleHyps->GetXaxis()->SetBinLabel(7,"3 (PPR+PID)");
  fMultipleHyps->GetXaxis()->SetBinLabel(8,"4 (PPR+PID)");
  fMultipleHyps->GetYaxis()->SetBinLabel(1,"PtBin 1");
  fMultipleHyps->GetYaxis()->SetBinLabel(2,"PtBin 2");
  fMultipleHyps->GetYaxis()->SetBinLabel(3,"PtBin 3");
  fMultipleHyps->GetYaxis()->SetBinLabel(4,"PtBin 4");
  fMultipleHyps->GetYaxis()->SetBinLabel(5,"PtBin 5");

  fMultipleHypsType = new TH2F("fMultipleHypsType", "Type of hyp. accepted for each candidate",8,0.,8.,5,0.,5.);
  fMultipleHypsType->SetMinimum(0);
  fMultipleHypsType->GetXaxis()->SetBinLabel(1,"D0");
  fMultipleHypsType->GetXaxis()->SetBinLabel(2,"D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(3,"2D0");
  fMultipleHypsType->GetXaxis()->SetBinLabel(4,"2D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(5,"D0+D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(6,"2D0+D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(7,"D0+2D0bar");
  fMultipleHypsType->GetXaxis()->SetBinLabel(8,"2D0+2D0bar");
  fMultipleHypsType->GetYaxis()->SetBinLabel(1,"PtBin 1");
  fMultipleHypsType->GetYaxis()->SetBinLabel(2,"PtBin 2");
  fMultipleHypsType->GetYaxis()->SetBinLabel(3,"PtBin 3");
  fMultipleHypsType->GetYaxis()->SetBinLabel(4,"PtBin 4");
  fMultipleHypsType->GetYaxis()->SetBinLabel(5,"PtBin 5");

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
        if (ptPart >= fPtBinH[0] && ptPart < fPtBinH[1])      fPIDSelBin1->Fill(0);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fPIDSelBin2->Fill(0);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fPIDSelBin3->Fill(0);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fPIDSelBin4->Fill(0);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fPIDSelBin5->Fill(0);
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
      Double_t binPt = -1;
      Int_t hypD01 = 0, hypD02 = 0, hypD0bar1 = 0, hypD0bar2 = 0;
      Int_t pid1 = 0, pid2 = 0, pidbar1 = 0, pidbar2 = 0;
      Int_t pidSelection = fCuts->IsSelectedFromPID(dIn, &pid1, &pid2, &pidbar1, &pidbar2);

   //number of CANDIDATES (regardless of hypotheses) passing PPR + PID - PAY ATTENTION: hypoth. for PID and PPR may not be the same!
   if (pidSelection > 0) {
        fPIDSel->Fill(1);
        if (ptPart >= fPtBinH[0] && ptPart < fPtBinH[1])      {fPIDSelBin1->Fill(1); binPt = 0;}
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) {fPIDSelBin2->Fill(1); binPt = 1;}
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) {fPIDSelBin3->Fill(1); binPt = 2;}
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) {fPIDSelBin4->Fill(1); binPt = 3;}
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) {fPIDSelBin5->Fill(1); binPt = 4;}
      }
   
	   //number of hypoteses accepted per candidate after PPR
	   if(selD01+selD02+selD0bar1+selD0bar2 == 1) fMultipleHyps->Fill((Double_t)0,binPt);
	   if(selD01+selD02+selD0bar1+selD0bar2 == 2) fMultipleHyps->Fill((Double_t)1,binPt);
	   if(selD01+selD02+selD0bar1+selD0bar2 == 3) fMultipleHyps->Fill((Double_t)2,binPt);
	   if(selD01+selD02+selD0bar1+selD0bar2 == 4) fMultipleHyps->Fill((Double_t)3,binPt);

   //combine PPD + PID cuts
   if (selD01 == 1 && pid1==1) hypD01 = 1;
   if (selD02 == 1 && pid2==1) hypD02 = 1;
   if (selD0bar1 == 1 && pidbar1==1) hypD0bar1 = 1;
   if (selD0bar2 == 1 && pidbar2==1) hypD0bar2 = 1;

   //number of CANDIDATES (regardless of hypotheses) passing PPR + PID - PAY ATTENTION: hypoth. for PID and PPR must match (at least one)!
   if (hypD01 == 1 || hypD02 == 1 || hypD0bar1 == 1 || hypD0bar2 == 1) {
        fPIDSel->Fill(2);
        if (ptPart >= fPtBinH[0] && ptPart < fPtBinH[1])      {fPIDSelBin1->Fill(2); binPt = 0;}
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) {fPIDSelBin2->Fill(2); binPt = 1;}
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) {fPIDSelBin3->Fill(2); binPt = 2;}
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) {fPIDSelBin4->Fill(2); binPt = 3;}
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) {fPIDSelBin5->Fill(2); binPt = 4;}
      }

	   //number of hypoteses accepted per candidate after PPR and PID
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 1) fMultipleHyps->Fill((Double_t)4,binPt);
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 2) fMultipleHyps->Fill((Double_t)5,binPt);
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 3) fMultipleHyps->Fill((Double_t)6,binPt);
	   if(hypD01+hypD02+hypD0bar1+hypD0bar2 == 4) fMultipleHyps->Fill((Double_t)7,binPt);

	   //type of hypoteses accepted per candidate after PPR and PID
	   if(hypD01+hypD02 == 1 && hypD0bar1+hypD0bar2 == 0) fMultipleHypsType->Fill((Double_t)0,binPt);
	   if(hypD01+hypD02 == 0 && hypD0bar1+hypD0bar2 == 1) fMultipleHypsType->Fill((Double_t)1,binPt);
	   if(hypD01+hypD02 == 2 && hypD0bar1+hypD0bar2 == 0) fMultipleHypsType->Fill((Double_t)2,binPt);
	   if(hypD01+hypD02 == 0 && hypD0bar1+hypD0bar2 == 2) fMultipleHypsType->Fill((Double_t)3,binPt);
	   if(hypD01+hypD02 == 1 && hypD0bar1+hypD0bar2 == 1) fMultipleHypsType->Fill((Double_t)4,binPt);
	   if(hypD01+hypD02 == 2 && hypD0bar1+hypD0bar2 == 1) fMultipleHypsType->Fill((Double_t)5,binPt);
	   if(hypD01+hypD02 == 1 && hypD0bar1+hypD0bar2 == 2) fMultipleHypsType->Fill((Double_t)6,binPt);
	   if(hypD01+hypD02 == 2 && hypD0bar1+hypD0bar2 == 2) fMultipleHypsType->Fill((Double_t)7,binPt);

   //Call function for reflection analysis
   if (hypD01+hypD02+hypD0bar1+hypD0bar2 > 0) AnalysisReflection(aodIn, dIn, hypD01, hypD02, hypD0bar1, hypD0bar2);

   //All histos are filled if Pt of candidate is greater than minimum of first bin (in this way: bin1+bin2+...binN = whole)
   if (ptPart > fPtBinH[0]) {
    
      // D01 hyp.
      if(hypD01==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk1->Pt(),trk3->Pt());
        if (ptPart < fPtBinH[1]) 	 		  {fhInvMassD0Sum10MevBin1->Fill(fmassD0[0]);
			         	 	  	  fhInvMassD0Sum5MevBin1->Fill(fmassD0[0]);
	                                 		  fhInvMassSumAll10MevBin1->Fill(fmassD0[0]);
        	                         		  fhInvMassSumAll5MevBin1->Fill(fmassD0[0]);}
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) {fhInvMassD0Sum10MevBin2->Fill(fmassD0[0]);
				      	 	  	      fhInvMassD0Sum5MevBin2->Fill(fmassD0[0]);
	                          	     		      fhInvMassSumAll10MevBin2->Fill(fmassD0[0]);
        	              	            		      fhInvMassSumAll5MevBin2->Fill(fmassD0[0]);}
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) {fhInvMassD0Sum10MevBin3->Fill(fmassD0[0]);
			   	       	  	  	      fhInvMassD0Sum5MevBin3->Fill(fmassD0[0]);
                          		      	  	      fhInvMassSumAll10MevBin3->Fill(fmassD0[0]);
                          	      		  	      fhInvMassSumAll5MevBin3->Fill(fmassD0[0]);}
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) {fhInvMassD0Sum10MevBin4->Fill(fmassD0[0]);
			   	       	  	  	      fhInvMassD0Sum5MevBin4->Fill(fmassD0[0]);
                          		      	  	      fhInvMassSumAll10MevBin4->Fill(fmassD0[0]);
                          	      		  	      fhInvMassSumAll5MevBin4->Fill(fmassD0[0]);}
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5])	{fhInvMassD0Sum10MevBin5->Fill(fmassD0[0]);
				         			fhInvMassD0Sum5MevBin5->Fill(fmassD0[0]);
			                 			fhInvMassSumAll10MevBin5->Fill(fmassD0[0]);
			                 			fhInvMassSumAll5MevBin5->Fill(fmassD0[0]);} 
     }
      // D02 hyp.
      if(hypD02==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk3->Pt(),trk1->Pt());
        if (ptPart < fPtBinH[1]) 	 		  {fhInvMassD0Sum10MevBin1->Fill(fmassD0[1]);
			         	 	  	  fhInvMassD0Sum5MevBin1->Fill(fmassD0[1]);
	                                 		  fhInvMassSumAll10MevBin1->Fill(fmassD0[1]);
        	                         		  fhInvMassSumAll5MevBin1->Fill(fmassD0[1]);}
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) {fhInvMassD0Sum10MevBin2->Fill(fmassD0[1]);
				      	 	  	      fhInvMassD0Sum5MevBin2->Fill(fmassD0[1]);
	                          	     		      fhInvMassSumAll10MevBin2->Fill(fmassD0[1]);
        	              	            		      fhInvMassSumAll5MevBin2->Fill(fmassD0[1]);}
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) {fhInvMassD0Sum10MevBin3->Fill(fmassD0[1]);
			   	       	  	  	      fhInvMassD0Sum5MevBin3->Fill(fmassD0[1]);
                          		      	  	      fhInvMassSumAll10MevBin3->Fill(fmassD0[1]);
                          	      		  	      fhInvMassSumAll5MevBin3->Fill(fmassD0[1]);}
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) {fhInvMassD0Sum10MevBin4->Fill(fmassD0[1]);
			   	       	  	  	      fhInvMassD0Sum5MevBin4->Fill(fmassD0[1]);
                          		      	  	      fhInvMassSumAll10MevBin4->Fill(fmassD0[1]);
                          	      		  	      fhInvMassSumAll5MevBin4->Fill(fmassD0[1]);}
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5])	{fhInvMassD0Sum10MevBin5->Fill(fmassD0[1]);
				        			 fhInvMassD0Sum5MevBin5->Fill(fmassD0[1]);
			                 			 fhInvMassSumAll10MevBin5->Fill(fmassD0[1]);
			                			 fhInvMassSumAll5MevBin5->Fill(fmassD0[1]);}
     }
      // D0bar1 hyp.
      if(hypD0bar1==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk0->Pt(),trk2->Pt());
        if (ptPart < fPtBinH[1]) 	 		  {fhInvMassD0barSum10MevBin1->Fill(fmassD0bar[0]);
			         	 	  	  fhInvMassD0barSum5MevBin1->Fill(fmassD0bar[0]);
                             		    		  fhInvMassSumAll10MevBin1->Fill(fmassD0bar[0]);
                                 			  fhInvMassSumAll5MevBin1->Fill(fmassD0bar[0]);}
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) {fhInvMassD0barSum10MevBin2->Fill(fmassD0bar[0]);
				      	 	  	      fhInvMassD0barSum5MevBin2->Fill(fmassD0bar[0]);
                           	     	 	  	      fhInvMassSumAll10MevBin2->Fill(fmassD0bar[0]);
                       	            	 		      fhInvMassSumAll5MevBin2->Fill(fmassD0bar[0]);}
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) {fhInvMassD0barSum10MevBin3->Fill(fmassD0bar[0]);
			   	         	 	      fhInvMassD0barSum5MevBin3->Fill(fmassD0bar[0]);
                          	      	 		      fhInvMassSumAll10MevBin3->Fill(fmassD0bar[0]);
                          	      	 		      fhInvMassSumAll5MevBin3->Fill(fmassD0bar[0]);}
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) {fhInvMassD0barSum10MevBin4->Fill(fmassD0bar[0]);
			   	         	 	      fhInvMassD0barSum5MevBin4->Fill(fmassD0bar[0]);
                          	      	 		      fhInvMassSumAll10MevBin4->Fill(fmassD0bar[0]);
                          	      	 		      fhInvMassSumAll5MevBin4->Fill(fmassD0bar[0]);}
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5])	{fhInvMassD0barSum10MevBin5->Fill(fmassD0bar[0]);
				         			fhInvMassD0barSum5MevBin5->Fill(fmassD0bar[0]);
			                 			fhInvMassSumAll10MevBin5->Fill(fmassD0bar[0]);
			                 			fhInvMassSumAll5MevBin5->Fill(fmassD0bar[0]);}   
     } 
      // D0bar2 hyp.
      if(hypD0bar2==1) {
	fPtSel->Fill(ptPart);
        fScatterP4PID->Fill(trk2->Pt(),trk0->Pt());
        if (ptPart < fPtBinH[1]) 	 		  {fhInvMassD0barSum10MevBin1->Fill(fmassD0bar[1]);
			         	 	  	  fhInvMassD0barSum5MevBin1->Fill(fmassD0bar[1]);
                             		    		  fhInvMassSumAll10MevBin1->Fill(fmassD0bar[1]);
                                 			  fhInvMassSumAll5MevBin1->Fill(fmassD0bar[1]);}
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) {fhInvMassD0barSum10MevBin2->Fill(fmassD0bar[1]);
				      	 	  	      fhInvMassD0barSum5MevBin2->Fill(fmassD0bar[1]);
                           	     	 	  	      fhInvMassSumAll10MevBin2->Fill(fmassD0bar[1]);
                       	            	 		      fhInvMassSumAll5MevBin2->Fill(fmassD0bar[1]);}
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) {fhInvMassD0barSum10MevBin3->Fill(fmassD0bar[1]);
			   	         	 	      fhInvMassD0barSum5MevBin3->Fill(fmassD0bar[1]);
                          	      	 		      fhInvMassSumAll10MevBin3->Fill(fmassD0bar[1]);
                          	      	 		      fhInvMassSumAll5MevBin3->Fill(fmassD0bar[1]);}
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) {fhInvMassD0barSum10MevBin4->Fill(fmassD0bar[1]);
			   	         	 	      fhInvMassD0barSum5MevBin4->Fill(fmassD0bar[1]);
                          	      	 		      fhInvMassSumAll10MevBin4->Fill(fmassD0bar[1]);
                          	      	 		      fhInvMassSumAll5MevBin4->Fill(fmassD0bar[1]);}
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) 	{fhInvMassD0barSum10MevBin5->Fill(fmassD0bar[1]);
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
    }
    if(unsetvtx) dIn->UnsetOwnPrimaryVtx();
  } // end loop on D0->K3pi

  printf("Number of selected D0->K3pi: %d\n",iOutCharm4Prong);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::AnalysisReflection(AliAODEvent* aodIn, AliAODRecoDecayHF4Prong* d, Int_t hypD01, Int_t hypD02, Int_t hypD0bar1, Int_t hypD0bar2)
{
  /*
  ---STUDY OF REFLECTIONS ON CANDIDATE IN ANALYSIS---
  Layers of TH2F fhRefl:
  0 = all candidates, single hyps selected only;
  1 = all candidates, multi hyps selected only;
  2 = true D0/D0bar only, all hypotheses selected;
  3 = true D0/D0bar only, single hyps selected only; -> I assume that the unique selected hypotesis is the right one!
  4 = true D0/D0bar only, multi hyps selected only;
  5 = true D0/D0bar only, multi TRUE hyps selected only;
  6 = true D0/D0bar only, multi FAKE hyps selected only;
  7 = false D0/D0bar only (background)
  Layers of TH2F fhReflD0 (idem for D0bar, inverting particles):
  0 = true D0 only, true hypothesis for both single and multi hyps cases;
  1 = true D0 only, true hypothesis for single hyps case only; -> I assume that the unique selected hypotesis is the right one!
  2 = true D0 only, true hypothesis for multi hyps case only;
  3 = true D0 only, other D0 wrong hypothesis (multi case, obviously);
  4 = true D0 only, D0bar1 wrong hypothesis (multi case, obviously);
  5 = true D0 only, D0bar2 wrong hypothesis (multi case, obviously);
  6 = true D0 only, D0bar1 + D0bar2 wrong hypotheses
  */

  Int_t flagMult = hypD01+hypD02+hypD0bar1+hypD0bar2; //single or multi hyps selected
  Int_t flagLayer1 = -1, flagLayer2 = -1, flagLayer3 = -1, flagLayer4 = -1; //to select layers in fhRefl

  if (flagMult == 1) {flagLayer1 = 0; flagLayer2 = 0; flagLayer3 = 0; flagLayer4 = 0;}
  else if (flagMult == 2 || flagMult == 3 || flagMult == 4) {flagLayer1 = 1; flagLayer2 = 1; flagLayer3 = 1; flagLayer4 = 1;}

  FillReflHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, flagLayer1, flagLayer2, flagLayer3, flagLayer4);  //Fill layers 0 and 1

  if (fMCTruth==0) return;

  //start of MC Truth phase
  TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodIn->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) {
    AliError("Could not find Monte-Carlo in AOD");
    return;}

  Int_t pdgCand = 421;
  Int_t pdgDgCharm4Prong[4]={321,211,211,211}; //pdg of daughters

  Int_t mcLabel = d->MatchToMC(pdgCand,mcArray,4,pdgDgCharm4Prong); //selection of true or false candidate (regardless of hypothesis) through MCtruth
  printf("MatchToMC = %d\n",mcLabel); 

  if (mcLabel==-1) { //fill layer 7 (background)
  FillReflHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 7, 7, 7, 7);
  return;}

  FillReflHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 2, 2, 2, 2); //fill layer 2 (true D0/D0bar)

  Int_t truthHyp = StudyMCTruth(mcArray, d); //calls function which studies which hypothesis is true for candidate

  if (flagMult == 1) { //fill layer 3 (true D0/D0bar - single hyps only)
    FillReflHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 3, 3, 3, 3);
    switch(truthHyp) { //fills fhReflD0 and fhReflD0bar with layers containing all true D0 or D0bars, for single and all hypotheses (layers 1, 0)
      case(1): FillReflD0Histos(d, hypD01, 0, 0, 0, 0, 0, 0, 0);
               FillReflD0Histos(d, hypD01, 0, 0, 0, 1, 1, 1, 1); //here I discard the cases in which only a wrong hyp is selected (very very few)
               break;
      case(2): FillReflD0Histos(d, 0, hypD02, 0, 0, 0, 0, 0, 0);
               FillReflD0Histos(d, 0, hypD02, 0, 0, 1, 1, 1, 1);
               break;
      case(3): FillReflD0barHistos(d, 0, 0, hypD0bar1, 0, 0, 0, 0, 0); 
	       FillReflD0barHistos(d, 0, 0, hypD0bar1, 0, 1, 1, 1, 1);
	       break;
      case(4): FillReflD0barHistos(d, 0, 0, 0, hypD0bar2, 0, 0, 0, 0); 
	       FillReflD0barHistos(d, 0, 0, 0, hypD0bar2, 1, 1, 1, 1);
	       break;
    }
  }
  else { 
    flagLayer1 = 6; flagLayer2 = 6; flagLayer3 = 6; flagLayer4 = 6; 
    switch(truthHyp) { //fills fhReflD0 and fhReflD0bar with layers containing all true D0 or D0bars, for multi and all hypotheses (layers 2, 0)
      case(1): FillReflD0Histos(d, hypD01, 0, 0, 0, 0, 0, 0, 0);
               FillReflD0Histos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 2, 3, 4, 5);
               FillReflD0Histos(d, 0, 0, hypD0bar1, hypD0bar2, 2, 3, 6, 6); //merge of opposite particle hyps (D0bar) in layer 6
	       flagLayer1 = 5;
               break;
      case(2): FillReflD0Histos(d, 0, hypD02, 0, 0, 0, 0, 0, 0);
               FillReflD0Histos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 3, 2, 4, 5);
               FillReflD0Histos(d, 0, 0, hypD0bar1, hypD0bar2, 0, 0, 6, 6); //merge of opposite particle hyps (D0bar) in layer 6
	       flagLayer2 = 5;
               break;
      case(3): FillReflD0barHistos(d, 0, 0, hypD0bar1, 0, 0, 0, 0, 0); 
	       FillReflD0barHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 4, 5, 2, 3);
               FillReflD0barHistos(d, hypD01, hypD02, 0, 0, 6, 6, 0, 0); //merge of opposite particle hyps (D0) in layer 6
	       flagLayer3 = 5;
	       break;
      case(4): FillReflD0barHistos(d, 0, 0, 0, hypD0bar2, 0, 0, 0, 0); 
	       FillReflD0barHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 4, 5, 3, 2);
               FillReflD0barHistos(d, hypD01, hypD02, 0, 0, 6, 6, 0, 0); //merge of opposite particle hyps (D0) in layer 6
	       flagLayer4 = 5;
	       break;
    }
    FillReflHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, flagLayer1, flagLayer2, flagLayer3, flagLayer4); //fill layers 5 and 6 (true and false hyps for multi)
    FillReflHistos(d, hypD01, hypD02, hypD0bar1, hypD0bar2, 4, 4, 4, 4); //fill layer 4 (true D0/D0bar - multi hyps only)
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskSESelectHF4Prong::StudyMCTruth(TClonesArray* mcArray, AliAODRecoDecayHF4Prong* d)
{
  /* 
  ---STUDY OF MCTRUTH ON CANDIDATE IN ANALYSIS---
  Flag Truth (output):
  0 = problems in daughter tracks found
  1 = candidate is D01 (piKpipi)
  2 = candidate is D02 (pipipiK)
  3 = candidate is D0bar1 (Kpipipi)
  4 = candidate is D0bar2 (pipiKpi)
  */

  Int_t truthHyp = 0;

  AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
  AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
  AliAODTrack *trk2 = (AliAODTrack*)d->GetDaughter(2);
  AliAODTrack *trk3 = (AliAODTrack*)d->GetDaughter(3);
  Int_t labels[4];
  Int_t pdg[4];
  labels[0] = trk0->GetLabel();
  labels[1] = trk1->GetLabel();
  labels[2] = trk2->GetLabel();
  labels[3] = trk3->GetLabel();
  if (labels[0]<=0 || labels[1]<=0 || labels[2]<=0 || labels[3]<=0) {AliWarning("Negative Label for daughter, skipping"); return truthHyp;}
  AliAODMCParticle* mc0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labels[0]));
  AliAODMCParticle* mc1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labels[1]));
  AliAODMCParticle* mc2 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labels[2]));
  AliAODMCParticle* mc3 = dynamic_cast<AliAODMCParticle*>(mcArray->At(labels[3]));
  if (!mc0 || !mc1 || !mc2 || !mc3) {AliWarning("At least one Daughter Particle not found in tree, skipping"); return truthHyp;}
  pdg[0] = TMath::Abs(mc0->GetPdgCode());
  pdg[1] = TMath::Abs(mc1->GetPdgCode());
  pdg[2] = TMath::Abs(mc2->GetPdgCode());
  pdg[3] = TMath::Abs(mc3->GetPdgCode());
  if (pdg[0]==211 && pdg[1]==321 && pdg[2]==211 && pdg[3]==211) truthHyp = 1;
  if (pdg[0]==211 && pdg[1]==211 && pdg[2]==211 && pdg[3]==321) truthHyp = 2;
  if (pdg[0]==321 && pdg[1]==211 && pdg[2]==211 && pdg[3]==211) truthHyp = 3;
  if (pdg[0]==211 && pdg[1]==211 && pdg[2]==321 && pdg[3]==211) truthHyp = 4;

  return truthHyp;
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::FillReflHistos(AliAODRecoDecayHF4Prong* d, Int_t hypD01, Int_t hypD02, Int_t hypD0bar1, Int_t hypD0bar2, Int_t flagLayer1, Int_t flagLayer2, Int_t flagLayer3, Int_t flagLayer4)
{

  Double_t ptPart = d->Pt();

    if (ptPart > fPtBinH[0]) {
      // D01 hyp.
      if(hypD01==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflBin1->Fill(fmassD0[0],(double)flagLayer1);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflBin2->Fill(fmassD0[0],(double)flagLayer1);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflBin3->Fill(fmassD0[0],(double)flagLayer1);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflBin4->Fill(fmassD0[0],(double)flagLayer1);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflBin5->Fill(fmassD0[0],(double)flagLayer1);
      }
      // D02 hyp.
      if(hypD02==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflBin1->Fill(fmassD0[1],(double)flagLayer2);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflBin2->Fill(fmassD0[1],(double)flagLayer2);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflBin3->Fill(fmassD0[1],(double)flagLayer2);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflBin4->Fill(fmassD0[1],(double)flagLayer2);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflBin5->Fill(fmassD0[1],(double)flagLayer2);
      }
      // D0bar1 hyp.
      if(hypD0bar1==1) {
	if (ptPart < fPtBinH[1]) 	 	    	      fhReflBin1->Fill(fmassD0bar[0],(double)flagLayer3);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflBin2->Fill(fmassD0bar[0],(double)flagLayer3);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflBin3->Fill(fmassD0bar[0],(double)flagLayer3);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflBin4->Fill(fmassD0bar[0],(double)flagLayer3);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflBin5->Fill(fmassD0bar[0],(double)flagLayer3);
      } 
      // D0bar2 hyp.
      if(hypD0bar2==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflBin1->Fill(fmassD0bar[1],(double)flagLayer4);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflBin2->Fill(fmassD0bar[1],(double)flagLayer4);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflBin3->Fill(fmassD0bar[1],(double)flagLayer4);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflBin4->Fill(fmassD0bar[1],(double)flagLayer4);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflBin5->Fill(fmassD0bar[1],(double)flagLayer4);
      }
    }

}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::FillReflD0Histos(AliAODRecoDecayHF4Prong* d, Int_t hypD01, Int_t hypD02, Int_t hypD0bar1, Int_t hypD0bar2, Int_t flagLayforD01, Int_t flagLayforD02, Int_t flagLayforD03, Int_t flagLayforD04)
{

  Double_t ptPart = d->Pt();

    if (ptPart > fPtBinH[0]) {
      // D01 hyp.
      if(hypD01==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflD0Bin1->Fill(fmassD0[0],(double)flagLayforD01);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0Bin2->Fill(fmassD0[0],(double)flagLayforD01);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0Bin3->Fill(fmassD0[0],(double)flagLayforD01);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0Bin4->Fill(fmassD0[0],(double)flagLayforD01);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0Bin5->Fill(fmassD0[0],(double)flagLayforD01);
      }
      // D02 hyp.
      if(hypD02==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflD0Bin1->Fill(fmassD0[1],(double)flagLayforD02);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0Bin2->Fill(fmassD0[1],(double)flagLayforD02);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0Bin3->Fill(fmassD0[1],(double)flagLayforD02);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0Bin4->Fill(fmassD0[1],(double)flagLayforD02);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0Bin5->Fill(fmassD0[1],(double)flagLayforD02);
      }
      // D0bar1 hyp.
      if(hypD0bar1==1) {
	if (ptPart < fPtBinH[1]) 	 	    	      fhReflD0Bin1->Fill(fmassD0bar[0],(double)flagLayforD03);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0Bin2->Fill(fmassD0bar[0],(double)flagLayforD03);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0Bin3->Fill(fmassD0bar[0],(double)flagLayforD03);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0Bin4->Fill(fmassD0bar[0],(double)flagLayforD03);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0Bin5->Fill(fmassD0bar[0],(double)flagLayforD03);
      } 
      // D0bar2 hyp.
      if(hypD0bar2==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflD0Bin1->Fill(fmassD0bar[1],(double)flagLayforD04);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0Bin2->Fill(fmassD0bar[1],(double)flagLayforD04);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0Bin3->Fill(fmassD0bar[1],(double)flagLayforD04);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0Bin4->Fill(fmassD0bar[1],(double)flagLayforD04);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0Bin5->Fill(fmassD0bar[1],(double)flagLayforD04);
      }
    }

}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::FillReflD0barHistos(AliAODRecoDecayHF4Prong* d, Int_t hypD01, Int_t hypD02, Int_t hypD0bar1, Int_t hypD0bar2, Int_t flagLayforD0bar1, Int_t flagLayforD0bar2, Int_t flagLayforD0bar3, Int_t flagLayforD0bar4)
{

  Double_t ptPart = d->Pt();

    if (ptPart > fPtBinH[0]) {
      // D01 hyp.
      if(hypD01==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflD0barBin1->Fill(fmassD0[0],(double)flagLayforD0bar1);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0barBin2->Fill(fmassD0[0],(double)flagLayforD0bar1);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0barBin3->Fill(fmassD0[0],(double)flagLayforD0bar1);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0barBin4->Fill(fmassD0[0],(double)flagLayforD0bar1);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0barBin5->Fill(fmassD0[0],(double)flagLayforD0bar1);
      }
      // D02 hyp.
      if(hypD02==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflD0barBin1->Fill(fmassD0[1],(double)flagLayforD0bar2);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0barBin2->Fill(fmassD0[1],(double)flagLayforD0bar2);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0barBin3->Fill(fmassD0[1],(double)flagLayforD0bar2);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0barBin4->Fill(fmassD0[1],(double)flagLayforD0bar2);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0barBin5->Fill(fmassD0[1],(double)flagLayforD0bar2);
      }
      // D0bar1 hyp.
      if(hypD0bar1==1) {
	if (ptPart < fPtBinH[1]) 	 	    	      fhReflD0barBin1->Fill(fmassD0bar[0],(double)flagLayforD0bar3);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0barBin2->Fill(fmassD0bar[0],(double)flagLayforD0bar3);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0barBin3->Fill(fmassD0bar[0],(double)flagLayforD0bar3);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0barBin4->Fill(fmassD0bar[0],(double)flagLayforD0bar3);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0barBin5->Fill(fmassD0bar[0],(double)flagLayforD0bar3);
      } 
      // D0bar2 hyp.
      if(hypD0bar2==1) {
        if (ptPart < fPtBinH[1]) 	 		      fhReflD0barBin1->Fill(fmassD0bar[1],(double)flagLayforD0bar4);
        else if (ptPart >= fPtBinH[1] && ptPart < fPtBinH[2]) fhReflD0barBin2->Fill(fmassD0bar[1],(double)flagLayforD0bar4);
        else if (ptPart >= fPtBinH[2] && ptPart < fPtBinH[3]) fhReflD0barBin3->Fill(fmassD0bar[1],(double)flagLayforD0bar4);
        else if (ptPart >= fPtBinH[3] && ptPart < fPtBinH[4]) fhReflD0barBin4->Fill(fmassD0bar[1],(double)flagLayforD0bar4);
        else if (ptPart >= fPtBinH[4] && ptPart < fPtBinH[5]) fhReflD0barBin5->Fill(fmassD0bar[1],(double)flagLayforD0bar4);
      }
    }

}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::SetPtBinH(Double_t* ptlimits)
{
  for(int i=0; i<6; i++) {fPtBinH[i]=ptlimits[i];}
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::PrintPtBinHandMCFlag()
{
  printf("PtBin limits---------\n");
  for (int i=0; i<5; i++) {
    printf("Bin %d = %.1f to %.1f\n",i+1,fPtBinH[i],fPtBinH[i+1]);
  }
  printf("MC Truth = %d\n",fMCTruth);
  printf("---------------------\n");
}

//________________________________________________________________________
void AliAnalysisTaskSESelectHF4Prong::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //

/*
  Double_t entries[10] = {0,0,0,0,0,0,0,0,0,0};
  for(int i=1;i<=fhReflD0Bin1->GetNbinsX();i++) {
    for(int j=1;j<=fhReflD0Bin1->GetNbinsY();j++) {
      entries[0] += fhReflD0Bin1->GetBinContent(i,j);
      entries[1] += fhReflD0Bin2->GetBinContent(i,j);
      entries[2] += fhReflD0Bin3->GetBinContent(i,j);
      entries[3] += fhReflD0Bin4->GetBinContent(i,j);
      entries[4] += fhReflD0Bin5->GetBinContent(i,j);
      entries[5] += fhReflD0barBin1->GetBinContent(i,j);
      entries[6] += fhReflD0barBin2->GetBinContent(i,j);
      entries[7] += fhReflD0barBin3->GetBinContent(i,j);
      entries[8] += fhReflD0barBin4->GetBinContent(i,j);
      entries[9] += fhReflD0barBin5->GetBinContent(i,j);
    }
  }
  fhReflD0Bin1->SetEntries(entries[0]);
  fhReflD0Bin2->SetEntries(entries[1]);  
  fhReflD0Bin3->SetEntries(entries[2]);
  fhReflD0Bin4->SetEntries(entries[3]);
  fhReflD0Bin5->SetEntries(entries[4]);
  fhReflD0barBin1->SetEntries(entries[5]);
  fhReflD0barBin2->SetEntries(entries[6]);
  fhReflD0barBin3->SetEntries(entries[7]);
  fhReflD0barBin4->SetEntries(entries[8]);
  fhReflD0barBin5->SetEntries(entries[9]);*/

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
  fOutput5 = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutput5) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fOutputC = dynamic_cast<TList*> (GetOutputData(6));
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

