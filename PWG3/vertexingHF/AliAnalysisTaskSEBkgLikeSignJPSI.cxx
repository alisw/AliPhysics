/**************************************************************************
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

///////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for reading both reconstructed JPSI -> ee candidates
// and like sign pairs and for drawing corresponding distributions
//
// Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
///////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEBkgLikeSignJPSI.h"

ClassImp(AliAnalysisTaskSEBkgLikeSignJPSI)

//________________________________________________________________________
AliAnalysisTaskSEBkgLikeSignJPSI::AliAnalysisTaskSEBkgLikeSignJPSI():
AliAnalysisTaskSE(),
fOutput(0), 
fHistMassJPSI(0),
fHistMassLS(0),
fHistCtsJPSI(0),           
fHistCtsLS(0),
fHistCtsLSpos(0),
fHistCtsLSneg(0),
fHistCPtaJPSI(0),          
fHistCPtaLS(0),
fHistd0d0JPSI(0),          
fHistd0d0LS(0),
fHistDCAJPSI(0),           
fHistDCALS(0),
fVHF(0),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)
{
  //
  // Default constructor
  //
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSEBkgLikeSignJPSI::AliAnalysisTaskSEBkgLikeSignJPSI(const char *name):
AliAnalysisTaskSE(name),
fOutput(0),
fHistMassJPSI(0),
fHistMassLS(0),
fHistCtsJPSI(0),
fHistCtsLS(0),
fHistCtsLSpos(0),
fHistCtsLSneg(0),
fHistCPtaJPSI(0),
fHistCPtaLS(0),
fHistd0d0JPSI(0),
fHistd0d0LS(0),
fHistDCAJPSI(0),
fHistDCALS(0),
fVHF(0),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)
{
  //
  // Default constructor
  //
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSEBkgLikeSignJPSI::~AliAnalysisTaskSEBkgLikeSignJPSI()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }

}
//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignJPSI::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEBkgLikeSignJPSI::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  fVHF->PrintStatus();

  return;

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignJPSI::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEBkgLikeSignJPSI::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();

  fHistMassJPSI = new TH1F("fHistMassJPSI", "J/#Psi invariant mass; M [GeV]; Entries",200,2.8,3.25);
  fHistMassJPSI->Sumw2();
  fHistMassJPSI->SetMinimum(0);
  fOutput->Add(fHistMassJPSI);

  fHistMassLS = new TH1F("fHistMassLS", "Like sign pairs invariant mass; M [GeV]; Entries",200,2.8,3.25);
  fHistMassLS->Sumw2();
  fHistMassLS->SetMinimum(0);
  fOutput->Add(fHistMassLS);

  fHistCtsJPSI = new TH1F("fHistCtsJPSI", "J/#Psi cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsJPSI->Sumw2();
  fHistCtsJPSI->SetMinimum(0);
  fOutput->Add(fHistCtsJPSI);

  fHistCtsLS = new TH1F("fHistCtsLS", "Like sign pairs cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsLS->Sumw2();
  fHistCtsLS->SetMinimum(0);
  fOutput->Add(fHistCtsLS);

  fHistCtsLSpos = new TH1F("fHistCtsLSpos", "Like sign ++ pairs cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsLSpos->Sumw2();
  fHistCtsLSpos->SetMinimum(0);
  fOutput->Add(fHistCtsLSpos);

  fHistCtsLSneg = new TH1F("fHistCtsLSneg", "Like sign -- pairs cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsLSneg->Sumw2();
  fHistCtsLSneg->SetMinimum(0);
  fOutput->Add(fHistCtsLSneg);

  fHistCPtaJPSI = new TH1F("fHistCPtaJPSI", "J/#Psi cosine of pointing angle; Cos#Theta_{point}; Entries",200,-1.,1.);
  fHistCPtaJPSI->Sumw2();
  fHistCPtaJPSI->SetMinimum(0);
  fOutput->Add(fHistCPtaJPSI);

  fHistCPtaLS = new TH1F("fHistCPtaLS", "Like sign pairs cosine of pointing angle; Cos#Theta_{point}; Entries",200,-1.,1.);
  fHistCPtaLS->Sumw2();
  fHistCPtaLS->SetMinimum(0);
  fOutput->Add(fHistCPtaLS);

  fHistd0d0JPSI = new TH1F("fHistd0d0JPSI", "J/#Psi product of impact parameters; d0xd0 [#mu m^{2}]; Entries",200,-100000.,100000.);
  fHistd0d0JPSI->Sumw2(); 
  fHistd0d0JPSI->SetMinimum(0);
  fOutput->Add(fHistd0d0JPSI);

  fHistd0d0LS = new TH1F("fHistd0d0LS", "Like sign pairs product of impact parameters; d0xd0 [#mu m^{2}]; Entries",200,-100000.,100000.);
  fHistd0d0LS->Sumw2();
  fHistd0d0LS->SetMinimum(0);
  fOutput->Add(fHistd0d0LS);

  fHistDCAJPSI = new TH1F("fHistDCAJPSI", "J/#Psi distance of closest approach; dca [10^{2}#mu m]; Entries",100,0.,5.);
  fHistDCAJPSI->Sumw2(); 
  fHistDCAJPSI->SetMinimum(0);
  fOutput->Add(fHistDCAJPSI);

  fHistDCALS = new TH1F("fHistDCALS", "Like sign pairs distance of closest approach; dca [10^{2}#mu m]; Entries",100,0.,5.);
  fHistDCALS->Sumw2(); 
  fHistDCALS->SetMinimum(0);
  fOutput->Add(fHistDCALS);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignJPSI::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  // load heavy flavour vertices
  TClonesArray *arrayVerticesHF =
    (TClonesArray*)aod->GetList()->FindObject("VerticesHF");
  if(!arrayVerticesHF) {
    printf("AliAnalysisTaskSEBkgLikeSignJPSI::UserExec: VerticesHF branch not found!\n");
    return;
  }

  // load JPSI->ee candidates                                                   
  TClonesArray *arrayJPSItoEle =
    (TClonesArray*)aod->GetList()->FindObject("JPSItoEle");
  if(!arrayJPSItoEle) {
    printf("AliAnalysisTaskSEBkgLikeSignJPSI::UserExec: JPSItoEle branch not found!\n");
    return;
  }

  // load like sign candidates
  TClonesArray *arrayLikeSign =
    (TClonesArray*)aod->GetList()->FindObject("LikeSign2Prong");
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSEBkgLikeSignJPSI::UserExec: LikeSign2Prong branch not found!\n");
    return;
  }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  // make trkIDtoEntry register (temporary)
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
    AliAODTrack *track = aod->GetTrack(it);
    trkIDtoEntry[track->GetID()]=it;
  }

  // loop over Like sign candidates
  Int_t nPosPairs=0,nNegPairs=0;
  Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
  printf("+++\n+++Number of like sign pairs ---> %d \n+++\n", nLikeSign);

  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
        d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
        unsetvtx=kTRUE;
    }
    Int_t okBtoJPSIls=0;
    if(d->SelectBtoJPSI(fVHF->GetBtoJPSICuts(),okBtoJPSIls)) {
       fHistMassLS->Fill(d->InvMassJPSIee());
       fHistCPtaLS->Fill(d->CosPointingAngle());
       fHistd0d0LS->Fill(1e8*d->Prodd0d0());
       fHistCtsLS->Fill(d->CosThetaStarJPSI());
       fHistDCALS->Fill(100*d->GetDCA());
       PostData(1,fOutput);
       AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
       if(!trk0) {
          trk0=aod->GetTrack(trkIDtoEntry[d->GetProngID(0)]);
          printf("references to standard AOD not available \n");
       }
       if((trk0->Charge())==1) {
          nPosPairs++;
          fHistCtsLSpos->Fill(d->CosThetaStarJPSI());
          PostData(1,fOutput);
        } else {
          nNegPairs++;
          fHistCtsLSneg->Fill(d->CosThetaStarJPSI());
          PostData(1,fOutput);
        }
      PostData(1,fOutput);
      }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  printf("------------ N. of positive pairs in Event ----- %d \n", nPosPairs);
  printf("------------ N. of negative pairs in Event ----- %d \n", nNegPairs);

  fTotPosPairs += nPosPairs;
  fTotNegPairs += nNegPairs;

  // loop over JPSI candidates
  Int_t nBtoJpsiToEle = arrayJPSItoEle->GetEntriesFast();
  printf("Number of like JPSI -> ee candidates ---> %d \n", nBtoJpsiToEle);

  for (Int_t iBtoJpsiToEle = 0; iBtoJpsiToEle < nBtoJpsiToEle; iBtoJpsiToEle++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayJPSItoEle->UncheckedAt(iBtoJpsiToEle);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    Int_t okBtoJPSI=0;
    if(d->SelectBtoJPSI(fVHF->GetBtoJPSICuts(),okBtoJPSI)) {
      fHistMassJPSI->Fill(d->InvMassJPSIee());
      fHistCtsJPSI->Fill(d->CosThetaStarJPSI());
      fHistd0d0JPSI->Fill(1e8*d->Prodd0d0());
      fHistCPtaJPSI->Fill(d->CosPointingAngle());
      fHistDCAJPSI->Fill(100*d->GetDCA());
      PostData(1,fOutput);
    }
   if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignJPSI::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEBkgLikeSignJPSI: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs*fTotNegPairs); 

  fHistMassJPSI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMassJPSI"));
  fHistMassLS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMassLS"));
  fHistCtsJPSI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsJPSI"));
  fHistCtsLS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLS"));
  fHistCtsLSpos = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLSpos"));
  fHistCtsLSneg = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLSneg"));
  fHistCPtaJPSI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCPtaJPSI"));
  fHistCPtaLS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCPtaLS"));
  fHistd0d0JPSI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistd0d0JPSI"));
  fHistd0d0LS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistd0d0LS"));
  fHistDCAJPSI = dynamic_cast<TH1F*>(fOutput->FindObject("fHistDCAJPSI"));
  fHistDCALS = dynamic_cast<TH1F*>(fOutput->FindObject("fHistDCALS"));

  fHistMassLS->Scale((1/fLsNormalization)*fHistMassLS->GetEntries());
  fHistCtsLS->Scale((1/fLsNormalization)*fHistCtsLS->GetEntries());
  fHistCtsLSpos->Scale((1/fLsNormalization)*fHistCtsLSpos->GetEntries());
  fHistCtsLSneg->Scale((1/fLsNormalization)*fHistCtsLSneg->GetEntries());
  fHistCPtaLS->Scale((1/fLsNormalization)*fHistCPtaLS->GetEntries());
  fHistd0d0LS->Scale((1/fLsNormalization)*fHistd0d0LS->GetEntries());
  fHistDCALS->Scale((1/fLsNormalization)*fHistDCALS->GetEntries());

  return;
}
