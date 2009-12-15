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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for D* candidates invariant mass histogram
// and comparison with the MC truth.
//
// Authors: Y.Wang, yifei@physi.uni-heidelberg.de
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TH1F.h>


#include "AliPID.h"
#include "AliTPCpidESD.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDStar.h"

ClassImp(AliAnalysisTaskSEDStar)


//________________________________________________________________________
AliAnalysisTaskSEDStar::AliAnalysisTaskSEDStar():
AliAnalysisTaskSE(),
fOutput(0), 
fhistMass(0),
fhistSgn(0),
fhistBkg(0),
fReadMC(0),
fPID(1),
fNSigma(3),
fVHF(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDStar::AliAnalysisTaskSEDStar(const char *name):
AliAnalysisTaskSE(name),
fOutput(0), 
fhistMass(0),
fhistSgn(0),
fhistBkg(0),
fReadMC(0),
fPID(1),
fNSigma(3),
fVHF(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSEDStar::~AliAnalysisTaskSEDStar()
{
  // Destructor

  if (fhistMass){
    delete fhistMass;
    fhistMass=0;
  }

  if (fhistSgn){
    delete fhistSgn;
    fhistSgn=0;
  }

  if (fhistBkg){
    delete fhistBkg;
    fhistBkg=0;
  }

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
void AliAnalysisTaskSEDStar::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEDStar::Init() \n");

//  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  // set dedidcated cuts
  //fVHF->SetD0fromDstarCuts(0.7,0.03,0.8,0.06,0.06,0.05,0.05,-0.0002,0.6); //2.p-p vertex reconstructed
  fVHF->SetD0fromDstarCuts(0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.);
  fVHF->SetDstarCuts(0.3, 0.1, 0.05, 100000000000.0, 0.5);
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDStar::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDStar::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();

    fhistMass = new TH1F("fhistMass","D^{*+}-D^{0} invariant mass; M [GeV]; Entries",200,0,0.3);
    fhistSgn = new TH1F("fhistSgn", "D^{*+}-D^{0} invariant mass - MC; M [GeV]; Entries",200,0,0.3);
    fhistBkg = new TH1F("fhistBkg", "D^{*+}-Background invariant mass - MC; M [GeV]; Entries",200,0,0.3);
 
  fOutput->Add(fhistMass);
  fOutput->Add(fhistSgn);
  fOutput->Add(fhistBkg);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDStar::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
 
  // load D*->D0 pi candidates                                                 
  TClonesArray *inputArray=0;
  
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
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
    }
  } else {
    inputArray=(TClonesArray*)aod->GetList()->FindObject("Dstar");
  }


  if(!inputArray) {
    printf("AliAnalysisTaskSEDStar::UserExec: input branch not found!\n");
    return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //vtx1->Print();
  
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSEDStar::UserExec: MC particles branch not found!\n");
      return;
    }
  
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDStar::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  // loop over D*->D0 pi candidates
  Int_t nInDstar = inputArray->GetEntriesFast();
  if(fDebug > 1) printf("Number of D*->D0 pi: %d\n",nInDstar);
///////
  
  for (Int_t iDstar = 0; iDstar < nInDstar; iDstar++) {
    AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*)inputArray->UncheckedAt(iDstar);

    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      d->Get2Prong()->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }
    
    if(d->SelectDstar(fVHF->GetDstarCuts(),fVHF->GetD0fromDstarCuts(),kTRUE)) {//selected
      Double_t invmassDelta = d->DeltaInvMass();
      //TVector3 p3Trk0(d->PxProng(0),d->PyProng(0),d->PzProng(0)); // pi_s
      //TVector3 p3Trk1(d->PxProng(1),d->PyProng(1),d->PzProng(1)); // D0
      //Double_t CosOpenAngle = p3Trk0.Dot(p3Trk1)/(p3Trk0.Mag()*p3Trk1.Mag());

      //PID of D0 daughters
      AliAODTrack *pos = (AliAODTrack*)d->Get2Prong()->GetDaughter(0);
      AliAODTrack *neg = (AliAODTrack*)d->Get2Prong()->GetDaughter(1);

      if (fPID) {
        if(fDebug > 1) printf("AnalysisTaskSEDStar::TPCPIDon \n");
        if(fDebug > 1) printf("AnalysisTaskSEDStar::NSigmaTPC: %f\n", fNSigma);

        if (d->Charge()>0){
          if(!SelectTPCPID(pos, 2, fNSigma)) continue;//pion+
          if(!SelectTPCPID(neg, 3, fNSigma)) continue;//kaon-
        }else{
          if(!SelectTPCPID(pos, 3, fNSigma)) continue;//kaon+
          if(!SelectTPCPID(neg, 2, fNSigma)) continue;//pion-
        }
      }

      if(fReadMC) {
      Int_t labDstar = d->MatchToMC(413,421,mcArray); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
      if(labDstar>=0) {
        //AliAODMCParticle *partDstar = (AliAODMCParticle*)mcArray->At(labDstar);
        //AliAODMCParticle *partPis = (AliAODMCParticle*)mcArray->At(partDstar->GetDaughter(1));
        //AliAODMCParticle *partD0 = (AliAODMCParticle*)mcArray->At(partDstar->GetDaughter(0));
        //AliAODMCParticle *partD0daughter = (AliAODMCParticle*)mcArray->At(partD0->GetDaughter(0));
        fhistSgn->Fill(invmassDelta);

      }
      else {//background
        fhistBkg->Fill(invmassDelta);
      }
      }
      //no MC info, just cut selection
      fhistMass->Fill(invmassDelta);

    } //else cout<<"NOT SELECTED"<<endl;

    if(unsetvtx) {
      d->UnsetOwnPrimaryVtx();
      d->Get2Prong()->UnsetOwnPrimaryVtx();
    }
  }  

  // Post the data
  PostData(1,fOutput);


  return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSEDStar::SelectTPCPID(AliAODTrack *trk, Int_t pid, Double_t nsig){//pid(0-4): {e,mu,pi,K,p}
  Bool_t flag=kTRUE;
  const Double_t mip=50.0, Res=0.07;
  if ((trk->GetStatus()&AliESDtrack::kTPCpid )==0) return flag;
  AliAODPid *detpid = trk->GetDetPid();
  Double_t dedx = detpid->GetTPCsignal()/mip;
  Double_t mass = AliPID::ParticleMass(pid);
  AliTPCpidESD tpcpid;
  Double_t mean = tpcpid.Bethe(trk->P()/mass);
  Double_t nsigma = (dedx-mean)/(Res*mean);
  if (TMath::Abs(nsigma)>nsig) flag=kFALSE;
  return flag;
}

//________________________________________________________________________
void AliAnalysisTaskSEDStar::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDStar: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  
  fhistMass = dynamic_cast<TH1F*>(fOutput->FindObject("fhistMass"));
  fhistSgn = dynamic_cast<TH1F*>(fOutput->FindObject("fhistSgn"));
  fhistBkg = dynamic_cast<TH1F*>(fOutput->FindObject("fhistBkg"));

  return;
}

