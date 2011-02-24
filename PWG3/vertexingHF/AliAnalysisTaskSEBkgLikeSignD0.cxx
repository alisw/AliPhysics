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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for reading both reconstructed D0->Kpi candidates
// and like sign pairs and for drawing corresponding distributions
//
// Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
///////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TH1F.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEBkgLikeSignD0.h"

ClassImp(AliAnalysisTaskSEBkgLikeSignD0)

//________________________________________________________________________
AliAnalysisTaskSEBkgLikeSignD0::AliAnalysisTaskSEBkgLikeSignD0():
AliAnalysisTaskSE(),
fOutput(0), 
fHistMassD0(0),
fHistMassLS(0),
fHistCtsD0(0),           
fHistCtsLS(0),
fHistCtsLSpos(0),
fHistCtsLSneg(0),
fHistCPtaD0(0),          
fHistCPtaLS(0),
fHistd0d0D0(0),          
fHistd0d0LS(0),
fHistDCAD0(0),           
fHistDCALS(0),
fVHF(0),
fNentries(0),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)
{
  //
  // Default constructor
  //
}

//________________________________________________________________________
AliAnalysisTaskSEBkgLikeSignD0::AliAnalysisTaskSEBkgLikeSignD0(const char *name):
AliAnalysisTaskSE(name),
fOutput(0),
fHistMassD0(0),
fHistMassLS(0),
fHistCtsD0(0),
fHistCtsLS(0),
fHistCtsLSpos(0),
fHistCtsLSneg(0),
fHistCPtaD0(0),
fHistCPtaLS(0),
fHistd0d0D0(0),
fHistd0d0LS(0),
fHistDCAD0(0),
fHistDCALS(0),
fVHF(0),
fNentries(0),
fTotPosPairs(0),
fTotNegPairs(0),
fLsNormalization(1.)
{
  //
  // Standard constructor
  //
  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TH1F container
  DefineOutput(2,TH1F::Class());  //My private output

}

//________________________________________________________________________
AliAnalysisTaskSEBkgLikeSignD0::~AliAnalysisTaskSEBkgLikeSignD0()
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

  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }

}
//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignD0::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEBkgLikeSignD0::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  //fVHF->SetD0toKpiCuts(0.7,0.02,0.8,0.7,0.7,1,1,-0.00025,0.8);
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignD0::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEBkgLikeSignD0::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();

  fHistMassD0 = new TH1F("fHistMassD0", "D0 invariant mass; M [GeV]; Entries",200,1.765,1.965);
  fHistMassD0->Sumw2();
  fHistMassD0->SetMinimum(0);
  fOutput->Add(fHistMassD0);

  fHistMassLS = new TH1F("fHistMassLS", "Like sign pairs invariant mass; M [GeV]; Entries",200,1.765,1.965);
  fHistMassLS->Sumw2();
  fHistMassLS->SetMinimum(0);
  fOutput->Add(fHistMassLS);

  fHistCtsD0 = new TH1F("fHistCtsD0", "D0 cosine of decay angle; Cos#Theta*; Entries",200,-1.,1.);
  fHistCtsD0->Sumw2();
  fHistCtsD0->SetMinimum(0);
  fOutput->Add(fHistCtsD0);

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

  fHistCPtaD0 = new TH1F("fHistCPtaD0", "D0 cosine of pointing angle; Cos#Theta_{point}; Entries",200,0,1.);
  fHistCPtaD0->Sumw2();
  fHistCPtaD0->SetMinimum(0);
  fOutput->Add(fHistCPtaD0);

  fHistCPtaLS = new TH1F("fHistCPtaLS", "Like sign pairs cosine of pointing angle; Cos#Theta_{point}; Entries",200,0,1.);
  fHistCPtaLS->Sumw2();
  fHistCPtaLS->SetMinimum(0);
  fOutput->Add(fHistCPtaLS);

  fHistd0d0D0 = new TH1F("fHistd0d0D0", "D0 product of impact parameters; d0xd0 [#mu m^{2}]; Entries",200,-100000.,100000.);
  fHistd0d0D0->Sumw2(); 
  fHistd0d0D0->SetMinimum(0);
  fOutput->Add(fHistd0d0D0);

  fHistd0d0LS = new TH1F("fHistd0d0LS", "Like sign pairs product of impact parameters; d0xd0 [#mu m^{2}]; Entries",200,-100000.,100000.);
  fHistd0d0LS->Sumw2();
  fHistd0d0LS->SetMinimum(0);
  fOutput->Add(fHistd0d0LS);

  fHistDCAD0 = new TH1F("fHistDCAD0", "D0 distance of closest approach; dca [10^{2}#mu m]; Entries",100,0.,5.);
  fHistDCAD0->Sumw2(); 
  fHistDCAD0->SetMinimum(0);
  fOutput->Add(fHistDCAD0);

  fHistDCALS = new TH1F("fHistDCALS", "Like sign pairs distance of closest approach; dca [10^{2}#mu m]; Entries",100,0.,5.);
  fHistDCALS->Sumw2(); 
  fHistDCALS->SetMinimum(0);
  fOutput->Add(fHistDCALS);

  fNentries=new TH1F("nentriesLS", "Look at the number of entries! it is = to the  number of AODs", 2,1.,2.);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignD0::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *arrayD0toKpi = 0;
  TClonesArray *arrayLikeSign = 0;

  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      // load D0 candidates                                                   
      arrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      // load like sign candidates
      arrayLikeSign=(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign2Prong");
    }
  } else if(aod) {
    // load D0 candidates                                                   
    arrayD0toKpi=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    // load like sign candidates
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign2Prong");
  }


  if(!aod || !arrayD0toKpi) {
    printf("AliAnalysisTaskSEBkgLikeSignD0::UserExec: D0toKpi branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSEBkgLikeSignD0::UserExec: LikeSign2Prong branch not found!\n");
    return;
  }

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;


  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  // make trkIDtoEntry register (temporary)
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
    AliAODTrack *track = aod->GetTrack(it);
    trkIDtoEntry[track->GetID()]=it;
  }

  //histogram filled with 1 for every AOD
  fNentries->Fill(1);
  PostData(2,fNentries);

  // loop over Like sign candidates
  Int_t nPosPairs=0,nNegPairs=0;
  Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
  if(fDebug>1) printf("+++\n+++Number of like sign pairs ---> %d \n+++\n", nLikeSign);

  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
        d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
        unsetvtx=kTRUE;
    }
    /*
    Int_t okD0ls=0; Int_t okD0barls=0;
    //if(d->SelectD0(fVHF->GetD0toKpiCuts(),okD0ls,okD0barls)) {
    if(d) {
       AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
       if(!trk0) {
          trk0=aod->GetTrack(trkIDtoEntry[d->GetProngID(0)]);
          printf("references to standard AOD not available \n");
       }
       if(okD0ls) fHistMassLS->Fill(d->InvMassD0());
       if(okD0barls) fHistMassLS->Fill(d->InvMassD0bar());
       fHistCPtaLS->Fill(d->CosPointingAngle());
       fHistd0d0LS->Fill(1e8*d->Prodd0d0());
       if(okD0ls) fHistCtsLS->Fill(d->CosThetaStarD0());
       if(okD0barls) fHistCtsLS->Fill(d->CosThetaStarD0bar());
       fHistDCALS->Fill(100*d->GetDCA());
       //PostData(1,fOutput);
       if((trk0->Charge())==1) {
          nPosPairs++;
          fHistCtsLSpos->Fill(d->CosThetaStarD0());
          //PostData(1,fOutput);
        } else {
          nNegPairs++;
          fHistCtsLSneg->Fill(d->CosThetaStarD0());
          //PostData(1,fOutput);
        }
      PostData(1,fOutput);
      }
    */
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  if(fDebug>1) printf("------------ N. of positive pairs in Event ----- %d \n", nPosPairs);
  if(fDebug>1) printf("------------ N. of negative pairs in Event ----- %d \n", nNegPairs);

  fTotPosPairs += nPosPairs;
  fTotNegPairs += nNegPairs;

  // loop over D0 candidates
  Int_t nD0toKpi = arrayD0toKpi->GetEntriesFast();
  if(fDebug>1) printf("Number of like D0 -> Kpi candidates ---> %d \n", nD0toKpi);

  for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    //Int_t okD0=0; Int_t okD0bar=0;
    //if(d->SelectD0(fVHF->GetD0toKpiCuts(),okD0,okD0bar)) {
    fHistMassD0->Fill(d->InvMassD0());
    fHistMassD0->Fill(d->InvMassD0bar());
    fHistCtsD0->Fill(d->CosThetaStarD0());
    fHistCtsD0->Fill(d->CosThetaStarD0bar());
    fHistd0d0D0->Fill(1e8*d->Prodd0d0());
    fHistCPtaD0->Fill(d->CosPointingAngle());
    fHistDCAD0->Fill(100*d->GetDCA());
    PostData(1,fOutput);
    
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEBkgLikeSignD0::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEBkgLikeSignD0: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs*fTotNegPairs); 

  fHistMassD0   = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMassD0"));
  fHistMassLS   = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMassLS"));
  fHistCtsD0    = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsD0"));
  fHistCtsLS    = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLS"));
  fHistCtsLSpos = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLSpos"));
  fHistCtsLSneg = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCtsLSneg"));
  fHistCPtaD0   = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCPtaD0"));
  fHistCPtaLS   = dynamic_cast<TH1F*>(fOutput->FindObject("fHistCPtaLS"));
  fHistd0d0D0   = dynamic_cast<TH1F*>(fOutput->FindObject("fHistd0d0D0"));
  fHistd0d0LS   = dynamic_cast<TH1F*>(fOutput->FindObject("fHistd0d0LS"));
  fHistDCAD0    = dynamic_cast<TH1F*>(fOutput->FindObject("fHistDCAD0"));
  fHistDCALS    = dynamic_cast<TH1F*>(fOutput->FindObject("fHistDCALS"));

  if(fLsNormalization>0) {
    if(fHistMassLS) fHistMassLS->Scale((1/fLsNormalization)*fHistMassLS->GetEntries());
    if(fHistCtsLS) fHistCtsLS->Scale((1/fLsNormalization)*fHistCtsLS->GetEntries());
    if(fHistCtsLSpos) fHistCtsLSpos->Scale((1/fLsNormalization)*fHistCtsLSpos->GetEntries());
    if(fHistCtsLSneg) fHistCtsLSneg->Scale((1/fLsNormalization)*fHistCtsLSneg->GetEntries());
    if(fHistCPtaLS) fHistCPtaLS->Scale((1/fLsNormalization)*fHistCPtaLS->GetEntries());
    if(fHistd0d0LS) fHistd0d0LS->Scale((1/fLsNormalization)*fHistd0d0LS->GetEntries());
    if(fHistDCALS) fHistDCALS->Scale((1/fLsNormalization)*fHistDCALS->GetEntries());
  }

  return;
}
