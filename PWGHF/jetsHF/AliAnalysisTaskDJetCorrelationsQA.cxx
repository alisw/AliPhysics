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
//
//
//  Analysis Task for MC and data QA
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// S. Antônio (University of São Paulo) antonio.silva@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL)  XMZhang@lbl.gov
// B. Trzeciak (Utrecht University) barbara.antonina.trzeciak@cern.ch
// J. Kvapil (University og Bitmingham) jakub.kvapil@cern.ch
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include "TROOT.h"
#include <THnSparse.h>
#include <TSystem.h>
#include <TObjectTable.h>
#include "AliMultSelection.h"

#include "AliAnalysisTaskDJetCorrelationsQA.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliParticleContainer.h"
#include "AliEmcalParticle.h"
#include "AliLocalRhoParameter.h"
#include "AliAnalysisTaskLocalRho.h"

#include "AliVertexingHFUtils.h"

ClassImp(AliAnalysisTaskDJetCorrelationsQA)


//_______________________________________________________________________________

AliAnalysisTaskDJetCorrelationsQA::AliAnalysisTaskDJetCorrelationsQA() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskFlavourJetCorrelations",kTRUE),
fUseMCInfo(kTRUE),
fIsPPData(kTRUE),
fIsPbPbData(kFALSE),
fMultiplicityEstimator(kNtrk10),
fUseReco(kTRUE),
fUsePythia(kFALSE),
fBuildRM(kFALSE),
fBuildRMEff(kFALSE),
fCandidateType(),
fCorrelationMethod(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(nullptr),
fMinMass(),
fMaxMass(),
fCandidateArray(nullptr),
fSideBandArray(nullptr),
fAnalyseDBkg(kFALSE),
fNAxesBigSparse(9),
fUseCandArray(kFALSE),
fUseSBArray(kFALSE),
fhstat(),
fhCentDjet(),
fhPtJetTrks(),
fhPhiJetTrks(),
fhEtaJetTrks(),
fhPtJet(),
fhPhiJet(),
fhEtaJet(),
fhInvMassptD(),
fhDiffSideBand(),
fhInvMassptDbg(),
fhPtPion(),
fhsDphiz(),
fResponseMatrix(),
fhPhiJetPtJet_incl_MC(),
fhPhiJetTrksPtJet_incl_MC(),
fhEtaJetTrksPtJet_incl_MC(),
fhEtaJetPtJet_incl_MC(),
fhAreaJetPtJet_incl_MC(),
fhJetTrksPtJet_incl_MC(),
fhPhiJetPtJet_incl_Reco(),
fhPhiJetTrksPtJet_incl_Reco(),
fhEtaJetTrksPtJet_incl_Reco(),
fhEtaJetPtJet_incl_Reco(),
fhAreaJetPtJet_incl_Reco(),
fhJetTrksPtJet_incl_Reco(),
fhPhiJetPtJet_Djet_MC(),
fhPhiJetTrksPtJet_Djet_MC(),
fhEtaJetTrksPtJet_Djet_MC(),
fhEtaJetPtJet_Djet_MC(),
fhAreaJetPtJet_Djet_MC(),
fhJetTrksPtJet_Djet_MC(),
fhPhiJetPtJet_Djet_Reco(),
fhPhiJetTrksPtJet_Djet_Reco(),
fhEtaJetTrksPtJet_Djet_Reco(),
fhEtaJetPtJet_Djet_Reco(),
fhAreaJetPtJet_Djet_Reco(),
fhJetTrksPtJet_Djet_Reco()

{
   //
   // Default ctor
}

//_______________________________________________________________________________

AliAnalysisTaskDJetCorrelationsQA::AliAnalysisTaskDJetCorrelationsQA(const Char_t* name, AliRDHFCuts* cuts,ECandidateType candtype) :
AliAnalysisTaskEmcalJet(name,kTRUE),
fUseMCInfo(kTRUE),
fIsPPData(kTRUE),
fIsPbPbData(kFALSE),
fMultiplicityEstimator(kNtrk10),
fUseReco(kTRUE),
fUsePythia(kFALSE),
fBuildRM(kFALSE),
fBuildRMEff(kFALSE),
fCandidateType(),
fCorrelationMethod(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(nullptr),
fMinMass(),
fMaxMass(),
fCandidateArray(nullptr),
fSideBandArray(nullptr),
fAnalyseDBkg(kFALSE),
fNAxesBigSparse(9),
fUseCandArray(kFALSE),
fUseSBArray(kFALSE),
fhstat(),
fhCentDjet(),
fhPtJetTrks(),
fhPhiJetTrks(),
fhEtaJetTrks(),
fhPtJet(),
fhPhiJet(),
fhEtaJet(),
fhInvMassptD(),
fhDiffSideBand(),
fhInvMassptDbg(),
fhPtPion(),
fhsDphiz(),
fResponseMatrix(),
fhPhiJetPtJet_incl_MC(),
fhPhiJetTrksPtJet_incl_MC(),
fhEtaJetTrksPtJet_incl_MC(),
fhEtaJetPtJet_incl_MC(),
fhAreaJetPtJet_incl_MC(),
fhJetTrksPtJet_incl_MC(),
fhPhiJetPtJet_incl_Reco(),
fhPhiJetTrksPtJet_incl_Reco(),
fhEtaJetTrksPtJet_incl_Reco(),
fhEtaJetPtJet_incl_Reco(),
fhAreaJetPtJet_incl_Reco(),
fhJetTrksPtJet_incl_Reco(),
fhPhiJetPtJet_Djet_MC(),
fhPhiJetTrksPtJet_Djet_MC(),
fhEtaJetTrksPtJet_Djet_MC(),
fhEtaJetPtJet_Djet_MC(),
fhAreaJetPtJet_Djet_MC(),
fhJetTrksPtJet_Djet_MC(),
fhPhiJetPtJet_Djet_Reco(),
fhPhiJetTrksPtJet_Djet_Reco(),
fhEtaJetTrksPtJet_Djet_Reco(),
fhEtaJetPtJet_Djet_Reco(),
fhAreaJetPtJet_Djet_Reco(),
fhJetTrksPtJet_Djet_Reco()
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //

   Info("AliAnalysisTaskFlavourJetCorrelations","Calling Constructor");
   fCuts=cuts;
   fCandidateType=candtype;
   const Int_t nptbins=fCuts->GetNPtBins();
   Float_t defaultSigmaD013[20]={0.012, 0.012, 0.012, 0.015, 0.015,0.018,0.018,0.020,0.020,0.030,0.030,0.037,0.040,0.040,0.040,0.040,0.040,0.040,0.040,0.040};

   switch(fCandidateType){
   case 0:
      fPDGmother=421;
      fNProngs=2;
      fPDGdaughters[0]=211;//pi
      fPDGdaughters[1]=321;//K
      fPDGdaughters[2]=0; //empty
      fPDGdaughters[3]=0; //empty
      fBranchName="D0toKpi";
      break;
   case 1:
      fPDGmother=413;
      fNProngs=3;
      fPDGdaughters[1]=211;//pi soft
      fPDGdaughters[0]=421;//D0
      fPDGdaughters[2]=211;//pi fromD0
      fPDGdaughters[3]=321; // K from D0
      fBranchName="Dstar";

      if(nptbins<20){
      	 for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]= defaultSigmaD013[ipt];
      } else {
      	 AliFatal(Form("Default sigma D0 not enough for %d pt bins, use SetSigmaD0ForDStar to set them",nptbins));
      }
      break;
   default:
      printf("%d not accepted!!\n",fCandidateType);
      break;
   }

   if(fCandidateType==kD0toKpi)SetMassLimits(0.15,fPDGmother);
   if(fCandidateType==kDstartoKpipi) SetMassLimits(0.035, fPDGmother);
   if(fUseCandArray) DefineInput(1, TClonesArray::Class());
   if(fUseSBArray) DefineInput(2, TClonesArray::Class());

      DefineOutput(1,TList::Class()); // histos
      DefineOutput(2,AliRDHFCuts::Class()); // my cuts

}

//_______________________________________________________________________________

AliAnalysisTaskDJetCorrelationsQA::~AliAnalysisTaskDJetCorrelationsQA() {
   //
   // destructor
   //

   Info("~AliAnalysisTaskFlavourJetCorrelations","Calling Destructor");

   delete fCuts;

}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::Init(){
   //
   // Initialization
   //

   if(fDebug > 1) printf("AnalysisTaskRecoJetCorrelations::Init() \n");

   switch(fCandidateType){
   case 0:
      {
      	 AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      	 copyfCuts->SetName("AnalysisCutsDzero");
      	 // Post the data
      	 PostData(2,copyfCuts);
      }
      break;
   case 1:
      {
      	 AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      	 copyfCuts->SetName("AnalysisCutsDStar");
      	 // Post the cuts
      	 PostData(2,copyfCuts);
      }
      break;
   default:
      return;
   }

   return;
}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::UserCreateOutputObjects() {
   // output
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
   AliAnalysisTaskEmcal::UserCreateOutputObjects();

   // define histograms
   // the TList fOutput is already defined in  AliAnalysisTaskEmcal::UserCreateOutputObjects()
   DefineHistoForAnalysis();
   PostData(1,fOutput);

   return;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskDJetCorrelationsQA::Run()
{
   // user exec from AliAnalysisTaskEmcal is used

   // Load the event
   AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
  if (matchingAODdeltaAODlevel<=0) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return kFALSE;
  }

   TClonesArray *arrayDStartoD0pi=nullptr;

   if (!aodEvent && AODEvent() && IsStandardAOD()) {

      // In case there is an AOD handler writing a standard AOD, use the AOD
      // event in memory rather than the input (ESD) event.
      aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());

      // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
      // have to taken from the AOD event hold by the AliAODExtension
      AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
      if(aodHandler->GetExtensions()) {
      	 AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      	 AliAODEvent *aodFromExt = ext->GetAOD();
      	 arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
      }
   } else if(aodEvent){
      arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject(fBranchName.Data());
   }

   if (!arrayDStartoD0pi) {
      AliInfo(Form("Could not find array %s, skipping the event",fBranchName.Data()));
      //  return;
   } else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast()));

   TClonesArray* mcArray = nullptr;
   if (fUseMCInfo) { //not used at the moment,uncomment return if you use
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
      	 printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
      }
   }

    //D meson candidates. Also background if is MC

    if(fUseCandArray)
    {
        fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
        if (!fCandidateArray) return kFALSE;
        for(Int_t icand=0; icand<fCandidateArray->GetEntriesFast(); icand++)
        {
            fhstat->Fill(2);
        }
    }
    if(fUseSBArray)
    {
        fSideBandArray = dynamic_cast<TClonesArray*>(GetInputData(2));
        if (!fSideBandArray) return kFALSE;
    }

    fhstat->Fill(0);

   // fix for temporary bug in ESDfilter
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return kFALSE;

    //Event selection
    Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
    TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
    if(!iseventselected) return kFALSE;

    fhstat->Fill(1);

    Float_t nTracklets = 0;
    Int_t nTrackletsEta10 = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
    Int_t nTrackletsEta16 = static_cast<Int_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.6,1.6));
    nTracklets = static_cast<Float_t>(nTrackletsEta10);

    // multiplicity estimator with VZERO
    Float_t vzeroMult=0;
    AliAODVZERO *vzeroAOD = (AliAODVZERO*)aodEvent->GetVZEROData();
    if(vzeroAOD) vzeroMult = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();
    Float_t multiplicity = nTracklets; // set to the Ntracklet estimator
    if(fMultiplicityEstimator==kVZERO) { multiplicity = vzeroMult; }

    fhMultiplicity->Fill(multiplicity);

    Float_t lPercentile = -1.;
    if(!fUseMCInfo && !fIsPPData) {
        lPercentile = fCuts->GetCentrality(aodEvent);
        fhCentDjet->Fill(lPercentile);
    }


// for MC response matrix of efficiency studies, fMultCand option only
if(fUseMCInfo && fBuildRMEff){

    AliJetContainer* mcjets = nullptr;

    if(!fAnalyseDBkg) mcjets = GetJetContainer(1);
    else mcjets = GetJetContainer(2);
    if(!mcjets) return kFALSE;

    AliParticleContainer *MCParticlesCont = mcjets->GetParticleContainer();

    mcjets->ResetCurrentID();
    AliEmcalJet* jet = nullptr;

    while ((jet = mcjets->GetNextJet()))
    {
        UInt_t rejectionReason = 0;
        Bool_t OKjet = mcjets->AcceptJet(jet, rejectionReason);
        if(!OKjet) {
            fhstat->Fill(5);
            continue;
        }
        //remove pure ghost jets
        if(jet->Pt() < 1E-9) continue;

        fhstat->Fill(3); //Jet accepted
        fhPhiJet->Fill(jet->Phi());
        fhEtaJet->Fill(jet->Eta());
        fhPtJet->Fill(jet->Pt());

        fhPhiJetPtJet_incl_MC->Fill(jet->Pt(),jet->Phi());
        fhEtaJetPtJet_incl_MC->Fill(jet->Pt(),jet->Eta());
        fhAreaJetPtJet_incl_MC->Fill(jet->Pt(),jet->Area());
        Int_t ntrjet=  jet->GetNumberOfTracks();
        fhJetTrksPtJet_incl_MC->Fill(jet->Pt(),ntrjet);

        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliAODMCParticle* jetTrk=(AliAODMCParticle*)jet->TrackAt(itrk,MCParticlesCont->GetArray());
            if (!jetTrk) continue;
            fhPtJetTrks->Fill(jetTrk->Pt());
            fhPhiJetTrks->Fill(jetTrk->Phi());
            fhEtaJetTrks->Fill(jetTrk->Eta());
            fhPhiJetTrksPtJet_incl_MC->Fill(jet->Pt(),jetTrk->Phi());
            fhEtaJetTrksPtJet_incl_MC->Fill(jet->Pt(),jetTrk->Eta());
        } //end loop on jet tracks
    } // end loop on mc jets

    //loop over rec jets
    AliJetContainer* recojets = GetJetContainer(0);

    AliParticleContainer *recoParticlesCont = recojets->GetParticleContainer();
    recojets->ResetCurrentID();
    AliEmcalJet* recojet = nullptr;
    while ((recojet = recojets->GetNextJet()))
    {
        UInt_t rejectionReason = 0;
        Bool_t OKjet = recojets->AcceptJet(recojet, rejectionReason);
        if(!OKjet) {
            continue;
        }
        //remove pure ghost jets
        if(recojet->Pt() < 1E-9) continue;

        fhPhiJetPtJet_incl_Reco->Fill(recojet->Pt(),recojet->Phi());
        fhEtaJetPtJet_incl_Reco->Fill(recojet->Pt(),recojet->Eta());
        fhAreaJetPtJet_incl_Reco->Fill(recojet->Pt(),recojet->Area());

        Int_t ntrjet=  recojet->GetNumberOfTracks();
        fhJetTrksPtJet_incl_Reco->Fill(recojet->Pt(),ntrjet);

        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliAODMCParticle* jetTrk=(AliAODMCParticle*)recojet->TrackAt(itrk,recoParticlesCont->GetArray());
            if (!jetTrk) continue;
            fhPhiJetTrksPtJet_incl_Reco->Fill(recojet->Pt(),jetTrk->Phi());
            fhEtaJetTrksPtJet_incl_Reco->Fill(recojet->Pt(),jetTrk->Eta());
        } //end loop on jet tracks
    } // end loop on reco jets

    // Get HF accepted MC jet
    AliEmcalJet* MCjet = nullptr;
    FindMCJet(MCjet);
    if(!MCjet) return kFALSE;
    //if( TMath::Abs(MCjet->Eta()) > (0.9 - mcjets->GetJetRadius()) ) return kFALSE;

    if(fCorrelationMethod==kConstituent)
    {
            if(fBuildRMEff==kTRUE) CreateMCResponseMatrix(MCjet, aodEvent);
    }
    /* the other method not enabled for now
     * else if(fCorrelationMethod==kAngular)
    {
        if(fCandidateArray->GetEntriesFast()>0) AngularCorrelationMethod(kFALSE,aodEvent);
        if(fAnalyseDBkg==kTRUE && fSideBandArray->GetEntriesFast()>0) AngularCorrelationMethod(kTRUE,aodEvent);
    }*/

}
else {

    AliJetContainer* JetCont = GetJetContainer(0);
    if(!JetCont) return kFALSE;
    AliParticleContainer *ParticlesCont = JetCont->GetParticleContainer();

    AliJetContainer* JetContSB = nullptr;
    AliParticleContainer *ParticlesContSB = nullptr;
    if(fAnalyseDBkg)
    {
        JetContSB = GetJetContainer(1);
        if(!JetContSB) return kFALSE;
        ParticlesContSB = JetContSB->GetParticleContainer();
    }

   //Distribution of all particles in the event
   Int_t ntrarr=ParticlesCont->GetNParticles();
   for(Int_t i=0;i<ntrarr;i++)
   {
       AliVParticle* jetTrk= ParticlesCont->GetParticle(i);
       if (!jetTrk) continue;
       AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(jetTrk);
       if (emcpart) jetTrk = emcpart->GetTrack();
       fhPtJetTrks->Fill(jetTrk->Pt());
       fhPhiJetTrks->Fill(jetTrk->Phi());
       fhEtaJetTrks->Fill(jetTrk->Eta());
   }

    JetCont->ResetCurrentID();
    AliEmcalJet* jet = nullptr;
    //loop over reco
    while ((jet = JetCont->GetNextJet()))
    {
        UInt_t rejectionReason = 0;
        Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
        if(!OKjet) {
            fhstat->Fill(5);
            continue;
        }

        //remove pure ghost jets
        if(jet->Pt() < 1E-9) continue;

        Double_t JetPtCorr = 0;
        if(fUseMCInfo && fUsePythia){
            JetPtCorr = jet->Pt();
        }
        else {
            JetPtCorr = jet->Pt() - jet->Area()*JetCont->GetRhoVal(); //background subtraction
            if(fLocalRho)
            {
                JetPtCorr = jet->Pt() - jet->Area()*fLocalRho->GetLocalVal(jet->Phi(),JetCont->GetJetRadius(),JetCont->GetRhoVal()); //
            }
        }

        fhstat->Fill(3); //Jet accepted
        fhPhiJet->Fill(jet->Phi());
        fhEtaJet->Fill(jet->Eta());
        fhPtJet->Fill(JetPtCorr);

        fhPhiJetPtJet_incl_Reco->Fill(jet->Pt(),jet->Phi());
        fhEtaJetPtJet_incl_Reco->Fill(jet->Pt(),jet->Eta());
        fhAreaJetPtJet_incl_Reco->Fill(jet->Pt(),jet->Area());
        Int_t ntrjet=  jet->GetNumberOfTracks();
        fhJetTrksPtJet_incl_Reco->Fill(jet->Pt(),ntrjet);

        for(Int_t itrk=0;itrk<ntrjet;itrk++){
            AliAODMCParticle* jetTrk=(AliAODMCParticle*)jet->TrackAt(itrk,ParticlesCont->GetArray());
            if (!jetTrk) continue;
            fhPhiJetTrksPtJet_incl_Reco->Fill(jet->Pt(),jetTrk->Phi());
            fhEtaJetTrksPtJet_incl_Reco->Fill(jet->Pt(),jetTrk->Eta());
        } //end loop on jet tracks

    }

    if(ParticlesCont->GetNParticles()>0) fhstat->Fill(2);

    if(fCorrelationMethod==kConstituent)
    {
        if(ParticlesCont->GetNParticles()>0) ConstituentCorrelationMethod(kFALSE,aodEvent);
        if(fAnalyseDBkg==kTRUE && ParticlesContSB->GetNParticles()>0) ConstituentCorrelationMethod(kTRUE,aodEvent);
    }

    else if(fCorrelationMethod==kAngular)
    {
        if(fCandidateArray->GetEntriesFast()>0) AngularCorrelationMethod(kFALSE,aodEvent);
        if(fAnalyseDBkg==kTRUE && fSideBandArray->GetEntriesFast()>0) AngularCorrelationMethod(kTRUE,aodEvent);
    }
}


   PostData(1,fOutput);
   return kTRUE;
}
void AliAnalysisTaskDJetCorrelationsQA::ConstituentCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent)
{

    AliJetContainer* JetCont = nullptr;

    if(!IsBkg) JetCont = GetJetContainer(0);
    else JetCont = GetJetContainer(1);

    Double_t rho = 0;
    if(!JetCont->GetRhoName().IsNull()) rho = JetCont->GetRhoVal();

    AliEmcalJet *jet;

    GetHFJet(jet,IsBkg);

    if(jet)
    {
        if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetCont->GetJetRadius(),JetCont->GetRhoVal());
        FillDJetHistograms(jet,rho,IsBkg,aodEvent);
    }

}
void AliAnalysisTaskDJetCorrelationsQA::AngularCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent)
{
    AliJetContainer* JetCont = GetJetContainer(0);

    Int_t ncand = 0;
    if(!IsBkg) ncand = fCandidateArray->GetEntriesFast();
    else ncand = fSideBandArray->GetEntriesFast();

    Double_t rho = 0;
    if(!JetCont->GetRhoName().IsNull()) rho = JetCont->GetRhoVal();

    for(Int_t icand = 0; icand<ncand; icand++)
    {
        AliVParticle* charm = nullptr;
        if(!IsBkg) charm = (AliVParticle*)fCandidateArray->At(icand);
        else charm = (AliVParticle*)fSideBandArray->At(icand);
        if(!charm) continue;

        Int_t JetTag = AliEmcalJet::kD0;
        if (fCandidateType == kDstartoKpipi) JetTag = AliEmcalJet::kDStar;
        //loop over jets
        JetCont->ResetCurrentID();
        AliEmcalJet* jet = nullptr;
        while ((jet = JetCont->GetNextJet()))
        {
            UInt_t rejectionReason = 0;
            Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
            if(!OKjet) continue;

            if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetCont->GetJetRadius(),JetCont->GetRhoVal());

            if(DeltaR(jet,charm,rho)<JetCont->GetJetRadius() && CheckDeltaR(jet,charm)<JetCont->GetJetRadius())
            {
                if(!IsBkg) jet->AddFlavourTag(JetTag);
                jet->AddFlavourTrack(charm);
                FillDJetHistograms(jet,rho,IsBkg,aodEvent);
            }
        }
    }

}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::CreateResponseMatrix(AliEmcalJet* jet)
{
    AliJetContainer* JetContRec = GetJetContainer(0);

    Double_t rho = 0;
    if(!fUsePythia) rho = JetContRec->GetRhoVal();

    AliEmcalJet* MCjet = nullptr;

    FindMCJet(MCjet);

    if(!jet) AliDebug(2, "No Reconstructed Level Jet Found!");
    else if(!MCjet) AliDebug(2, "No Generated Level Jet Found!");
    else
    {
        if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetContRec->GetJetRadius(),JetContRec->GetRhoVal());
        AliVParticle *Drec = jet->GetFlavourTrack();
        AliVParticle *Dgen = MCjet->GetFlavourTrack();
        Double_t zRec = Z(Drec,jet,rho);
        Double_t zGen = Z(Dgen,MCjet,0);
        Double_t JetPtRec = jet->Pt() - rho*jet->Area();
        Double_t JetPtGen = MCjet->Pt();
        Double_t JetYRec = jet->Y();
        Double_t JetYGen = MCjet->Y();
        Double_t JetEtaRec = jet->Eta();
        Double_t JetEtaGen = MCjet->Eta();
        Double_t JetPhiRec = jet->Phi();
        Double_t JetPhiGen = MCjet->Phi();
        Double_t JetnTrkRec = jet->GetNumberOfTracks();
        Double_t JetnTrkGen = MCjet->GetNumberOfTracks();
        Double_t JetAreaRec = jet->GetNumberOfTracks();
        Double_t JetAreaGen = MCjet->GetNumberOfTracks();
        AliAODMCParticle* DTrk=(AliAODMCParticle*)Dgen;
        Double_t DYGen = DTrk->Y();
        Double_t DPtGen = DTrk->Pt();
        Double_t DEtaGen = DTrk->Eta();
        Double_t DPhiGen = DTrk->Phi();
        Double_t DYRec = -999;
        Double_t DPtRec = -999;
        Double_t DEtaRec = -999;
        Double_t DPhiRec = -999;

        Double_t pTRes = JetPtGen ? ((JetPtRec - JetPtGen) / JetPtGen) : -999;
        Double_t zRes = zGen ? ((zRec - zGen) / zGen) : -999;

        Int_t pdgMeson = 413;
        if (fCandidateType == kD0toKpi) pdgMeson = 421;

        fhPhiJetPtJet_Djet_Reco->Fill(JetPtRec,jet->Phi());
        fhEtaJetPtJet_Djet_Reco->Fill(JetPtRec,jet->Eta());
        fhAreaJetPtJet_Djet_Reco->Fill(JetPtRec,jet->Area());
        fhJetTrksPtJet_Djet_Reco->Fill(JetPtRec,jet->GetNumberOfTracks());

        fhPhiJetPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->Phi());
        fhEtaJetPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->Eta());
        fhAreaJetPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->Area());
        fhJetTrksPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->GetNumberOfTracks());

        for(Int_t i=0;i<jet->GetNumberOfTracks();i++)
        {
            AliVParticle* jetTrk= jet->Track(i);
            if (!jetTrk) continue;
            fhPhiJetTrksPtJet_Djet_Reco->Fill(JetPtRec,jetTrk->Phi());
            fhEtaJetTrksPtJet_Djet_Reco->Fill(JetPtRec,jetTrk->Eta());
        }
        for(Int_t i=0;i<MCjet->GetNumberOfTracks();i++)
        {
            AliVParticle* jetTrk= MCjet->Track(i);
            if (!jetTrk) continue;
            fhPhiJetTrksPtJet_Djet_MC->Fill(MCjet->Pt(),jetTrk->Phi());
            fhEtaJetTrksPtJet_Djet_MC->Fill(MCjet->Pt(),jetTrk->Eta());
        }

        if(fCandidateType==kD0toKpi)
        {
            AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Drec;
            DYRec = dzero->Y(pdgMeson);
            DPtRec = dzero->Pt();
            DEtaRec = dzero->Eta();
            DPhiRec = dzero->Phi();

        }
        else if(fCandidateType==kDstartoKpipi)
        {
            AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Drec;
            DYRec = dstar->Y(pdgMeson);
            DPtRec = dstar->Pt();
            DEtaRec = dstar->Eta();
            DPhiRec = dstar->Phi();

        }

        Double_t fillRM[25] = {zRec,DPtRec,DYRec,DEtaRec,DPhiRec,JetPtRec,JetYRec,JetEtaRec,JetPhiRec,JetAreaRec,JetnTrkRec,zGen,DPtGen,DYGen,DEtaGen,DPhiGen,JetPtGen,JetYGen,JetEtaGen,JetPhiGen,JetAreaGen,JetnTrkGen,pTRes,zRes,-1};
        fResponseMatrix->Fill(fillRM,1.);
    }
}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::CreateMCResponseMatrix(AliEmcalJet* MCjet, AliAODEvent* aodEvent)
{
    if(!MCjet) AliDebug(2, "No Generated Level Jet Found!");

    Float_t nTracklets = static_cast<Float_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
    // multiplicity estimator with VZERO
    Float_t vzeroMult=0;
    AliAODVZERO *vzeroAOD = (AliAODVZERO*)aodEvent->GetVZEROData();
    if(vzeroAOD) vzeroMult = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();
    Float_t multiplicity = nTracklets; // set to the Ntracklet estimator
    if(fMultiplicityEstimator==kVZERO) { multiplicity = vzeroMult; }

    AliVParticle *Dgen = MCjet->GetFlavourTrack();
    AliAODMCParticle* DTrk = (AliAODMCParticle*)Dgen;
    Int_t pdg = DTrk->GetPdgCode();

    Double_t phiGen = MCjet->Phi();
    Double_t nTrkGen = MCjet->GetNumberOfTracks();

    Double_t zGen = Z(Dgen,MCjet,0);
    Double_t DPtGen = Dgen->Pt();
    Double_t DYGen = DTrk->Y();
    Double_t DEtaGen = Dgen->Eta();
    Double_t DPhiGen = Dgen->Phi();
    Double_t JetPtGen = MCjet->Pt();
    Double_t JetYGen = MCjet->Y();
    Double_t JetEtaGen = MCjet->Eta();
    Double_t JetPhiGen = MCjet->Phi();
    Double_t JetAreaGen = MCjet->Area();
    Double_t JetnTrkGen = MCjet->GetNumberOfTracks();

    Double_t zRec = -999;
    Double_t DPtRec = -999;
    Double_t DYRec = -999;
    Double_t DEtaRec = -999;
    Double_t DPhiRec = -999;
    Double_t JetPtRec = -999;
    Double_t JetYRec = -999;
    Double_t JetEtaRec = -999;
    Double_t JetPhiRec = -999;
    Double_t JetAreaRec = -999;
    Double_t JetnTrkRec = -999;

    Double_t pTRes = -999;
    Double_t zRes = -999;

    AliJetContainer* JetContRec = GetJetContainer(0);
    if(JetContRec){
        AliEmcalJet* jet;
        GetHFJet(jet,kFALSE);

        if (jet){
            AliVParticle *Drec = jet->GetFlavourTrack();

            Double_t rho = 0;
            if(!fUsePythia){
                rho = JetContRec->GetRhoVal();
                if(fLocalRho) rho = fLocalRho->GetLocalVal(jet->Phi(),JetContRec->GetJetRadius(),JetContRec->GetRhoVal());
            }
            zRec = 0;
            if(rho>0) zRec = Z(Drec,jet,rho);
            else zRec = Z(Drec,jet);

            DPtRec = Drec->Pt();
            DYRec = Drec->Y();
            DEtaRec = Drec->Eta();
            DPhiRec = Drec->Phi();
            JetPtRec = jet->Pt() - rho*jet->Area();
            JetYRec = jet->Y();
            JetEtaRec = jet->Eta();
            JetPhiRec = jet->Phi();
            JetAreaRec = jet->Area();
            JetnTrkRec = jet->GetNumberOfTracks();

            pTRes = JetPtGen ? ((JetPtRec - JetPtGen) / JetPtGen) : -999;
            zRes = zGen ? ((zRec - zGen) / zGen) : -999;

            Int_t pdgMeson = 413;
            if (fCandidateType == kD0toKpi) pdgMeson = 421;

             if(fCandidateType==kD0toKpi)
            {
                AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Drec;
                DYRec = dzero->Y(pdgMeson);

            }
            else if(fCandidateType==kDstartoKpipi)
            {
                AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Drec;
                DYRec = dstar->Y(pdgMeson);

            }

            Bool_t bDInEMCalAcc=InEMCalAcceptance(Drec);
            Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

            fhPhiJetPtJet_Djet_Reco->Fill(JetPtRec,jet->Phi());
            fhEtaJetPtJet_Djet_Reco->Fill(JetPtRec,jet->Eta());
            fhAreaJetPtJet_Djet_Reco->Fill(JetPtRec,jet->Area());
            fhJetTrksPtJet_Djet_Reco->Fill(JetPtRec,jet->GetNumberOfTracks());

            fhPhiJetPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->Phi());
            fhEtaJetPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->Eta());
            fhAreaJetPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->Area());
            fhJetTrksPtJet_Djet_MC->Fill(MCjet->Pt(),MCjet->GetNumberOfTracks());

            for(Int_t i=0;i<jet->GetNumberOfTracks();i++)
            {
                AliVParticle* jetTrk= jet->Track(i);
                if (!jetTrk) continue;
                fhPhiJetTrksPtJet_Djet_Reco->Fill(JetPtRec,jetTrk->Phi());
                fhEtaJetTrksPtJet_Djet_Reco->Fill(JetPtRec,jetTrk->Eta());
            }
            for(Int_t i=0;i<MCjet->GetNumberOfTracks();i++)
            {
                AliVParticle* jetTrk= MCjet->Track(i);
                if (!jetTrk) continue;
                fhPhiJetTrksPtJet_Djet_MC->Fill(MCjet->Pt(),jetTrk->Phi());
                fhEtaJetTrksPtJet_Djet_MC->Fill(MCjet->Pt(),jetTrk->Eta());
            }

            if(fCandidateType==kD0toKpi)
            {
                AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Drec;
                FillHistogramsD0JetCorr(dzero,zRec,Drec->Pt(),JetPtRec,JetEtaRec,kFALSE,bDInEMCalAcc,bJetInEMCalAcc,aodEvent,pdg,JetPhiRec,JetnTrkRec,JetYRec,JetAreaRec);
            }
            else if(fCandidateType==kDstartoKpipi)
            {
                AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Drec;
                FillHistogramsDstarJetCorr(dstar,zRec,Drec->Pt(),JetPtRec,JetEtaRec,kFALSE,bDInEMCalAcc,bJetInEMCalAcc);
            }

        } // if HF reco jet
    } // if jet cont reco

    Double_t fillRM[25] = {zRec,DPtRec,DYRec,DEtaRec,DPhiRec,JetPtRec,JetYRec,JetEtaRec,JetPhiRec,JetAreaRec,JetnTrkRec,zGen,DPtGen,DYGen,DEtaGen,DPhiGen,JetPtGen,JetYGen,JetEtaGen,JetPhiGen,JetAreaGen,JetnTrkGen,pTRes,zRes,multiplicity};
    fResponseMatrix->Fill(fillRM,1.);

}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::FillDJetHistograms(AliEmcalJet* jet, Double_t rho, Bool_t IsBkg, AliAODEvent* aodEvent)
{


    AliVParticle *Dmeson = jet->GetFlavourTrack(0);
    Double_t JetPtCorr = jet->Pt() - rho*jet->Area();
    Double_t JetEtaRec = jet->Eta();
    Double_t JetYRec = jet->Y();
    Double_t JetAreaRec = jet->Area();
    Double_t JetPhiRec = jet->Phi();
    Double_t JetNtrkRec = jet->GetNumberOfTracks();
    Double_t z = 0;
    if(rho>0) z = Z(Dmeson,jet,rho);
    else z = Z(Dmeson,jet);

    if(IsBkg==kFALSE && fBuildRM==kTRUE) CreateResponseMatrix(jet);
    else{
        fhPhiJetPtJet_Djet_Reco->Fill(JetPtCorr,jet->Phi());
        fhEtaJetPtJet_Djet_Reco->Fill(JetPtCorr,jet->Eta());
        fhAreaJetPtJet_Djet_Reco->Fill(JetPtCorr,jet->Area());
        fhJetTrksPtJet_Djet_Reco->Fill(JetPtCorr,jet->GetNumberOfTracks());

        for(Int_t i=0;i<jet->GetNumberOfTracks();i++)
        {
            AliVParticle* jetTrk= jet->Track(i);
            if (!jetTrk) continue;
            fhPhiJetTrksPtJet_Djet_Reco->Fill(JetPtCorr,jetTrk->Phi());
            fhEtaJetTrksPtJet_Djet_Reco->Fill(JetPtCorr,jetTrk->Eta());
        }
    }

    Bool_t bDInEMCalAcc=InEMCalAcceptance(Dmeson);
    Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

    if(fCandidateType==kD0toKpi)
    {
        AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Dmeson;
        FillHistogramsD0JetCorr(dzero,z,Dmeson->Pt(),JetPtCorr,JetEtaRec,IsBkg,bDInEMCalAcc,bJetInEMCalAcc,aodEvent,-999,JetPhiRec,JetNtrkRec,JetYRec,JetAreaRec);
    }
    else if(fCandidateType==kDstartoKpipi)
    {
        AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Dmeson;
        FillHistogramsDstarJetCorr(dstar,z,Dmeson->Pt(),JetPtCorr,JetEtaRec,IsBkg,bDInEMCalAcc,bJetInEMCalAcc);
    }


}

void AliAnalysisTaskDJetCorrelationsQA::GetHFJet(AliEmcalJet*& jet, Bool_t IsBkg)
{
    AliJetContainer* JetCont = nullptr;

    if(!IsBkg) JetCont = GetJetContainer(0);
    else JetCont = GetJetContainer(1);

    AliParticleContainer *ParticlesCont = JetCont->GetParticleContainer();

    Bool_t JetIsHF = kFALSE;

    JetCont->ResetCurrentID();
    while ((jet = JetCont->GetNextJet()))
    {
        UInt_t rejectionReason = 0;
        Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
        if(!OKjet) continue;

        Int_t JetTag = AliEmcalJet::kD0;
        TString recoDecayClassName("AliAODRecoDecayHF2Prong");
        if (fCandidateType == kDstartoKpipi)
        {
            JetTag = AliEmcalJet::kDStar;
            recoDecayClassName = "AliAODRecoCascadeHF";
        }
        //loop on jet particles
        Int_t ntrjet=  jet->GetNumberOfTracks();
        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliVParticle* jetTrk=jet->TrackAt(itrk,ParticlesCont->GetArray());
            if(!jetTrk) continue;
            AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(jetTrk);
            if(emcpart) jetTrk = emcpart->GetTrack();
            if(strcmp(jetTrk->ClassName(),recoDecayClassName)==0)
            {
                JetIsHF = kTRUE;
                if(!IsBkg) jet->AddFlavourTag(JetTag);
                jet->AddFlavourTrack(jetTrk);
                break;
            }
        } //end loop on jet tracks
        if(JetIsHF) break;
    } //end of jet loop

    if(!JetIsHF) jet = nullptr;

}
//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::FindMCJet(AliEmcalJet*& mcjet)
{
    Bool_t HFMCjet = kFALSE;

    AliJetContainer* mcjets = nullptr;

    if(!fAnalyseDBkg) mcjets = GetJetContainer(1);
    else mcjets = GetJetContainer(2);

    AliParticleContainer *ParticlesCont = mcjets->GetParticleContainer();
    mcjets->ResetCurrentID();

    Int_t njet=0;

    while ((mcjet = mcjets->GetNextAcceptJet()))
    {
        njet++;
        //loop on jet particles
        Int_t ntrjet=  mcjet->GetNumberOfTracks();

        for(Int_t itrk=0;itrk<ntrjet;itrk++)
        {
            AliAODMCParticle* jetTrk=(AliAODMCParticle*)mcjet->TrackAt(itrk,ParticlesCont->GetArray());

            if(TMath::Abs(jetTrk->GetPdgCode())==fPDGmother)
            {
                HFMCjet=kTRUE;
                mcjet->AddFlavourTrack(jetTrk);
                break;
            }
        } //end loop on jet tracks
        if(HFMCjet==kTRUE) break;
    }

    if(!HFMCjet) mcjet = nullptr;
}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::Terminate(Option_t*)
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   Info("Terminate"," terminate");
   AliAnalysisTaskSE::Terminate();

   fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1));
   if (!fOutput) {
      printf("ERROR: fOutput not available\n");
      return;
   }
}

//_______________________________________________________________________________

void  AliAnalysisTaskDJetCorrelationsQA::SetMassLimits(Double_t range, Int_t pdg){
   Float_t mass=0;
   Int_t abspdg=TMath::Abs(pdg);

   mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
   // compute the Delta mass for the D*
   if(fCandidateType==kDstartoKpipi){
      Float_t mass1=0;
      mass1=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      mass = mass-mass1;

      fMinMass = mass-range;
      fMaxMass = mass+range;
  }
  else {
      fMinMass = 1.5;
      fMaxMass = 2.2;
  }

   AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
   if (fMinMass<0 || fMaxMass<=0 || fMaxMass<fMinMass) AliFatal("Wrong mass limits!\n");
}

//_______________________________________________________________________________

void  AliAnalysisTaskDJetCorrelationsQA::SetMassLimits(Double_t lowlimit, Double_t uplimit){
   if(uplimit>lowlimit)
   {
      fMinMass = lowlimit;
      fMaxMass = uplimit;
   }
   else{
      printf("Error! Lower limit larger than upper limit!\n");
      fMinMass = uplimit - uplimit*0.2;
      fMaxMass = uplimit;
   }
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskDJetCorrelationsQA::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
   if(nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   if(!width) return kFALSE;
   for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
   return kTRUE;
}

//_______________________________________________________________________________

Double_t AliAnalysisTaskDJetCorrelationsQA::Z(AliVParticle* part,AliEmcalJet* jet, Double_t rho) const{

    Double_t p[3],pj[3];
    Bool_t okpp=part->PxPyPz(p);
    Bool_t okpj=jet->PxPyPz(pj);

    if(!okpp || !okpj)
    {
        printf("Problems getting momenta\n");
        return -999;
    }

    //Background Subtraction
    //It corrects the each component of the jet momentum for Z calculation

    pj[0] = jet->Px() - jet->Area()*(rho*TMath::Cos(jet->AreaPhi()));
    pj[1] = jet->Py() - jet->Area()*(rho*TMath::Sin(jet->AreaPhi()));
    pj[2] = jet->Pz() - jet->Area()*(rho*TMath::SinH(jet->AreaEta()));

    return Z(p,pj);
}
//_______________________________________________________________________________

Double_t AliAnalysisTaskDJetCorrelationsQA::Z(AliVParticle* part,AliEmcalJet* jet) const{

    Double_t p[3],pj[3];
    Bool_t okpp=part->PxPyPz(p);
    Bool_t okpj=jet->PxPyPz(pj);

    if(!okpp || !okpj)
    {
        printf("Problems getting momenta\n");
        return -999;
    }

    return Z(p,pj);
}
//_______________________________________________________________________________
Double_t AliAnalysisTaskDJetCorrelationsQA::Z(Double_t* p, Double_t *pj) const{

   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1]+pj[2]*pj[2];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet2);
   return z;
}


//_______________________________________________________________________________
Double_t AliAnalysisTaskDJetCorrelationsQA::ZT(Double_t* p, Double_t *pj) const{

   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1])/(pjet2);
   return z;
}

//_______________________________________________________________________________

Bool_t  AliAnalysisTaskDJetCorrelationsQA::DefineHistoForAnalysis(){

   // Statistics
   Int_t nbins=8;
   if(fUseMCInfo) nbins+=2;

   fhstat=new TH1I("hstat","Statistics",nbins,-0.5,nbins-0.5);
   fhstat->GetXaxis()->SetBinLabel(1,"N ev anal");
   fhstat->GetXaxis()->SetBinLabel(2,"N ev sel");
   fhstat->GetXaxis()->SetBinLabel(3,"N cand sel");
   fhstat->GetXaxis()->SetBinLabel(4,"N jets");
   fhstat->GetXaxis()->SetBinLabel(5,"N cand in jet");
   fhstat->GetXaxis()->SetBinLabel(6,"N jet rej");
   fhstat->GetXaxis()->SetBinLabel(7,"N cand sel & !jet");
   fhstat->GetXaxis()->SetBinLabel(8,"N jets & cand");
   if(fUseMCInfo) {
    fhstat->GetXaxis()->SetBinLabel(3,"N Signal sel & jet");
    fhstat->GetXaxis()->SetBinLabel(5,"N Signal in jet");
    fhstat->GetXaxis()->SetBinLabel(9,"N Bkg sel & jet");
    fhstat->GetXaxis()->SetBinLabel(10,"N Bkg in jet");
   }
   fhstat->SetNdivisions(1);
   fOutput->Add(fhstat);


   const Int_t nbinsmass=300;
   const Int_t nbinspttrack=500;
   const Int_t nbinsptjet=200;
   const Int_t nbinsptD=100;
   const Int_t nbinsphi=200;
   const Int_t nbinseta=100;
   const Int_t nbinsArea=10;
   const Int_t nbinNtracks=50;

   //binning for THnSparse
   const Int_t nbinsSpsmass=280;
   const Int_t nbinsSpsptjet=200;
   const Int_t nbinsSpsptD=100;
   const Int_t nbinsSpsz=160;
   const Int_t nbinsSpsy=150;

   Int_t nbinsCent=100;
   if(fIsPPData || fUseMCInfo) nbinsCent = 1;
   Int_t nbinsMult = 500;

   const Float_t pttracklims[2]={0.,200.};
   const Float_t ptjetlims[2]={-50.,150.};
   const Float_t ptDlims[2]={0.,50.};
   const Float_t zlims[2]={-1.2,2};
   const Float_t philims[2]={0.,6.3};
   const Float_t etalims[2]={-1.5,1.5};
   const Float_t nTracksLims[2]={-0.5,49.5};

   Float_t centLims[2] = {0,100};
   Float_t multLims[2] = {-0.5,499.5};
   if(fIsPbPbData) multLims[1] = 9999.5;

   fhCentDjet=new TH1F("hCentDjet","Centrality",nbinsCent,centLims[0],centLims[1]);
   fOutput->Add(fhCentDjet);

   fhMultiplicity=new TH1F("fhMultiplicity","Multiplicity",nbinsMult,multLims[0],multLims[1]);
   fOutput->Add(fhMultiplicity);

   // jet related fistograms

   fhPhiJetTrks    = new TH1F("hPhiJetTrks","Jet tracks #phi distribution; #phi (rad)",  nbinsphi,philims[0],philims[1]);
   fhPhiJetTrks->Sumw2();
   fhEtaJetTrks    = new TH1F("hEtaJetTrks","Jet tracks #eta distribution; #eta",  nbinseta,etalims[0],etalims[1]);
   fhEtaJetTrks->Sumw2();
   fhPtJetTrks     = new TH1F("hPtJetTrks",  "Jet tracks Pt distribution; p_{T} (GeV/c)",nbinspttrack,pttracklims[0],pttracklims[1]);
   fhPtJetTrks->Sumw2();

   fhPhiJet    = new TH1F("hPhiJet","Jet #phi distribution; #phi (rad)",  nbinsphi,philims[0],philims[1]);
   fhPhiJet->Sumw2();
   fhEtaJet    = new TH1F("hEtaJet","Jet #eta distribution; #eta", nbinseta,etalims[0],etalims[1]);
   fhEtaJet->Sumw2();
   fhPtJet      = new TH1F("hPtJet",  "Jet Pt distribution; p_{T} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
   fhPtJet->Sumw2();

   fhPhiJetPtJet_incl_MC  = new TH2F("fhPhiJetPtJet_incl_MC","Jet #phi distribution vs Jet p_{T} MC; Jet p_{T};#varphi^{jet} (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhPhiJetTrksPtJet_incl_MC = new TH2F("fhPhiJetTrksPtJet_incl_MC","Jet tracks #phi distribution vs Jet p_{T} MC; Jet p_{T};#phi (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhEtaJetPtJet_incl_MC  = new TH2F("fhEtaJetPtJet_incl_MC","Jet #eta distribution vs Jet p_{T} MC; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhEtaJetTrksPtJet_incl_MC = new TH2F("fhEtaJetTrksPtJet_incl_MC","Jet tracks #eta distribution vs Jet p_{T} MC; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhAreaJetPtJet_incl_MC= new TH2F("fhAreaJetPtJet_incl_MC","Jet area distribution vs Jet p_{T} MC; Jet p_{T};Area",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsArea,0,1);
   fhJetTrksPtJet_incl_MC= new TH2F("fhJetTrksPtJet_incl_MC","Jet N tracks distribution vs Jet p_{T} MC; Jet p_{T};N tracks",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinNtracks,nTracksLims[0],nTracksLims[1]);
   fhPhiJetPtJet_incl_Reco  = new TH2F("fhPhiJetPtJet_incl_Reco","Jet #phi distribution vs Jet p_{T} rec; Jet p_{T};#phi (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhPhiJetTrksPtJet_incl_Reco = new TH2F("fhPhiJetTrksPtJet_incl_Reco","Jet tracks #phi distribution vs Jet p_{T} rec; Jet p_{T};#phi (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhEtaJetPtJet_incl_Reco  = new TH2F("fhEtaJetPtJet_incl_Reco","Jet #eta distribution vs Jet p_{T} rec; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhEtaJetTrksPtJet_incl_Reco = new TH2F("fhEtaJetTrksPtJet_incl_Reco","Jet tracks #eta distribution vs Jet p_{T} rec; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhAreaJetPtJet_incl_Reco= new TH2F("fhAreaJetPtJet_incl_Reco","Jet area distribution vs Jet p_{T} rec; Jet p_{T};Area",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsArea,0,1);
   fhJetTrksPtJet_incl_Reco= new TH2F("fhJetTrksPtJet_incl_Reco","Jet N tracks distribution vs Jet p_{T} rec; Jet p_{T};N tracks",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinNtracks,nTracksLims[0],nTracksLims[1]);

   fhPhiJetPtJet_Djet_MC  = new TH2F("fhPhiJetPtJet_Djet_MC","D Jet #phi distribution vs Jet p_{T} MC; Jet p_{T};#varphi^{jet} (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhPhiJetTrksPtJet_Djet_MC = new TH2F("fhPhiJetTrksPtJet_Djet_MC","D Jet tracks #phi distribution vs Jet p_{T} MC; Jet p_{T};#phi (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhEtaJetPtJet_Djet_MC  = new TH2F("fhEtaJetPtJet_Djet_MC","D Jet #eta distribution vs Jet p_{T} MC; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhEtaJetTrksPtJet_Djet_MC = new TH2F("fhEtaJetTrksPtJet_Djet_MC","D Jet tracks #eta distribution vs Jet p_{T} MC; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhAreaJetPtJet_Djet_MC= new TH2F("fhAreaJetPtJet_Djet_MC","D Jet area distribution vs Jet p_{T} MC; Jet p_{T};Area",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsArea,0,1);
   fhJetTrksPtJet_Djet_MC= new TH2F("fhJetTrksPtJet_Djet_MC","D Jet N tracks distribution vs Jet p_{T} MC; Jet p_{T};N tracks",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinNtracks,nTracksLims[0],nTracksLims[1]);
   fhPhiJetPtJet_Djet_Reco  = new TH2F("fhPhiJetPtJet_Djet_Reco","D Jet #phi distribution vs Jet p_{T} rec; Jet p_{T};#phi (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhPhiJetTrksPtJet_Djet_Reco = new TH2F("fhPhiJetTrksPtJet_Djet_Reco","D Jet tracks #phi distribution vs Jet p_{T} rec; Jet p_{T};#phi (rad)",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsphi,philims[0],philims[1]);
   fhEtaJetPtJet_Djet_Reco  = new TH2F("fhEtaJetPtJet_Djet_Reco","D Jet #eta distribution vs Jet p_{T} rec; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhEtaJetTrksPtJet_Djet_Reco = new TH2F("fhEtaJetTrksPtJet_Djet_Reco","D Jet tracks #eta distribution vs Jet p_{T} rec; Jet p_{T};#eta",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinseta,etalims[0],etalims[1]);
   fhAreaJetPtJet_Djet_Reco= new TH2F("fhAreaJetPtJet_Djet_Reco","D Jet area distribution vs Jet p_{T} rec; Jet p_{T};Area",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinsArea,0,1);
   fhJetTrksPtJet_Djet_Reco= new TH2F("fhJetTrksPtJet_Djet_Reco","D Jet N tracks distribution vs Jet p_{T} rec; Jet p_{T};N tracks",nbinsptjet,ptjetlims[0],ptjetlims[1], nbinNtracks,nTracksLims[0],nTracksLims[1]);


   fhPhiJetPtJet_incl_MC->Sumw2();
   fhPhiJetTrksPtJet_incl_MC->Sumw2();
   fhEtaJetPtJet_incl_MC->Sumw2();
   fhEtaJetTrksPtJet_incl_MC ->Sumw2();
   fhAreaJetPtJet_incl_MC->Sumw2();
   fhJetTrksPtJet_incl_MC->Sumw2();
   fhPhiJetPtJet_incl_Reco->Sumw2();
   fhPhiJetTrksPtJet_incl_Reco->Sumw2();
   fhEtaJetPtJet_incl_Reco->Sumw2();
   fhEtaJetTrksPtJet_incl_Reco->Sumw2();
   fhAreaJetPtJet_incl_Reco->Sumw2();
   fhJetTrksPtJet_incl_Reco->Sumw2();

   fhPhiJetPtJet_Djet_MC->Sumw2();
   fhPhiJetTrksPtJet_Djet_MC->Sumw2();
   fhEtaJetPtJet_Djet_MC->Sumw2();
   fhEtaJetTrksPtJet_Djet_MC ->Sumw2();
   fhAreaJetPtJet_Djet_MC->Sumw2();
   fhJetTrksPtJet_Djet_MC->Sumw2();
   fhPhiJetPtJet_Djet_Reco->Sumw2();
   fhPhiJetTrksPtJet_Djet_Reco->Sumw2();
   fhEtaJetPtJet_Djet_Reco->Sumw2();
   fhEtaJetTrksPtJet_Djet_Reco->Sumw2();
   fhAreaJetPtJet_Djet_Reco->Sumw2();
   fhJetTrksPtJet_Djet_Reco->Sumw2();

   fOutput->Add(fhPhiJetTrks);
   fOutput->Add(fhEtaJetTrks);
   fOutput->Add(fhPtJetTrks);
   fOutput->Add(fhPhiJet);
   fOutput->Add(fhEtaJet);
   fOutput->Add(fhPtJet);

   fOutput->Add(fhPhiJetPtJet_incl_MC);
   fOutput->Add(fhPhiJetTrksPtJet_incl_MC);
   fOutput->Add(fhEtaJetPtJet_incl_MC);
   fOutput->Add(fhEtaJetTrksPtJet_incl_MC);
   fOutput->Add(fhAreaJetPtJet_incl_MC);
   fOutput->Add(fhJetTrksPtJet_incl_MC);
   fOutput->Add(fhPhiJetPtJet_incl_Reco);
   fOutput->Add(fhPhiJetTrksPtJet_incl_Reco);
   fOutput->Add(fhEtaJetPtJet_incl_Reco);
   fOutput->Add(fhEtaJetTrksPtJet_incl_Reco);
   fOutput->Add(fhAreaJetPtJet_incl_Reco);
   fOutput->Add(fhJetTrksPtJet_incl_Reco);

   fOutput->Add(fhPhiJetPtJet_Djet_MC);
   fOutput->Add(fhPhiJetTrksPtJet_Djet_MC);
   fOutput->Add(fhEtaJetPtJet_Djet_MC);
   fOutput->Add(fhEtaJetTrksPtJet_Djet_MC);
   fOutput->Add(fhAreaJetPtJet_Djet_MC);
   fOutput->Add(fhJetTrksPtJet_Djet_MC);
   fOutput->Add(fhPhiJetPtJet_Djet_Reco);
   fOutput->Add(fhPhiJetTrksPtJet_Djet_Reco);
   fOutput->Add(fhEtaJetPtJet_Djet_Reco);
   fOutput->Add(fhEtaJetTrksPtJet_Djet_Reco);
   fOutput->Add(fhAreaJetPtJet_Djet_Reco);
   fOutput->Add(fhJetTrksPtJet_Djet_Reco);

      if(fCandidateType==kDstartoKpipi)
      {
	       if(fAnalyseDBkg){
      	    fhDiffSideBand = new TH2F("hDiffSideBand","M(kpipi)-M(kpi) Side Band Background",nbinsmass,fMinMass,fMaxMass,nbinsptD, ptDlims[0],ptDlims[1]);
      	    fhDiffSideBand->SetStats(kTRUE);
      	    fhDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV");
      	    fhDiffSideBand->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	    fhDiffSideBand->Sumw2();
      	    fOutput->Add(fhDiffSideBand);
      	 }

      	 fhPtPion = new TH1F("hPtPion","Primary pions candidates pt ",500,0,10);
      	 fhPtPion->SetStats(kTRUE);
      	 fhPtPion->GetXaxis()->SetTitle("GeV/c");
      	 fhPtPion->GetYaxis()->SetTitle("Entries");
      	 fhPtPion->Sumw2();
      	 fOutput->Add(fhPtPion);

      }


      fhInvMassptD = new TH2F("hInvMassptD","D (Delta R < Rjet) invariant mass distribution p_{T}^{j} > threshold",nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      fhInvMassptD->SetStats(kTRUE);
      fhInvMassptD->GetXaxis()->SetTitle("mass (GeV)");
      fhInvMassptD->GetYaxis()->SetTitle("#it{p}_{t}^{D} (GeV/c)");
      fhInvMassptD->Sumw2();

      fOutput->Add(fhInvMassptD);

      if(fUseMCInfo){
	       fhInvMassptDbg = new TH2F("hInvMassptDbg","Background D (Delta R < Rjet) invariant mass distribution p_{T}^{j} > threshold",nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      	 fhInvMassptDbg->GetXaxis()->SetTitle(Form("%s (GeV)",(fCandidateType==kDstartoKpipi) ? "M(Kpipi)" : "M(Kpi)"));
      	 fhInvMassptDbg->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	 fhInvMassptDbg->Sumw2();
      	 fOutput->Add(fhInvMassptDbg);

      }

    fhsDphiz=nullptr; //definition below according to the switches

    if(fUseMCInfo){
        if(fCandidateType==kDstartoKpipi){
          AliInfo("Creating a 9 axes container (MB background candidates)");
          const Int_t nAxis=9;
          const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,nbinsSpsy,nbinsSpsy, 2, 2, 2};
          const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,etalims[0],etalims[0], -0.5,-0.5,-0.5};
          const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,etalims[1],etalims[1], 1.5, 1.5 , 1.5};
          fNAxesBigSparse=nAxis;
          fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass., y^{D}, #eta^{jet}, Bkg?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);

        }
        else{
          const Int_t nAxis=15;
          const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,nbinsSpsy,nbinsSpsy,nbinsSpsmass,nbinsSpsmass, 2, 2, 2,nbinsphi,nbinNtracks,nbinsMult,nbinsCent};
          const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,etalims[0],etalims[0],fMinMass,fMinMass, -0.5,-0.5,-0.5,philims[0],nTracksLims[0],multLims[0],centLims[0]};
          const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,etalims[1],etalims[1],fMaxMass,fMaxMass, 1.5, 1.5 , 1.5,philims[1],nTracksLims[1],multLims[1],centLims[1]};
          fNAxesBigSparse=nAxis;
          fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass., y^{D}, #eta^{jet}, mass_{true}, mass_{refl}, Bkg?, D in EMCal acc?, jet in EMCal acc?, #varphi^{jet}, Ntracks^{jet},mult,centr", nAxis, nbinsSparse, minSparse, maxSparse);
        }

    } else{
        AliInfo("Creating a 14 axes container");
        const Int_t nAxis=14;
        const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptD,nbinsSpsy,nbinseta,nbinsphi,nbinsSpsptjet,nbinsSpsy,nbinseta,nbinsphi,nbinsArea,nbinNtracks,nbinsSpsmass,nbinsMult,nbinsCent};
        const Double_t minSparse[nAxis]={zlims[0],ptDlims[0],etalims[0],etalims[0],philims[0],ptjetlims[0],etalims[0],etalims[0],philims[0],0,nTracksLims[0],fMinMass,multLims[0],centLims[0]};
        const Double_t maxSparse[nAxis]={zlims[1],ptDlims[1],etalims[1],etalims[1],philims[1],ptjetlims[1],etalims[1],etalims[1],philims[1],1,nTracksLims[1],fMaxMass,multLims[1],centLims[1]};
        fNAxesBigSparse=nAxis;
        fhsDphiz=new THnSparseF("hsDphiz","z, p_{T}^{D}, y^{D}, #eta^{D}, #varphi^{D},p_{T}^{jet},y^{jet}, #eta^{jet}, #varphi^{jet}, Area^{jet}, Ntracks^{jet}, mass^{D}, mult, cent", nAxis, nbinsSparse, minSparse, maxSparse);
    }

    if(!fhsDphiz) AliFatal("No THnSparse created");
    fhsDphiz->Sumw2();

    fOutput->Add(fhsDphiz);

    if(fBuildRM == kTRUE || fBuildRMEff == kTRUE)
    {
        const Int_t nRMAxis=25;
        const Int_t RMnbins[nRMAxis] = {nbinsSpsz,nbinsSpsptD,nbinsSpsy,nbinseta,nbinsphi,nbinsSpsptjet,nbinsSpsy,nbinseta,nbinsphi,nbinsArea,nbinNtracks,nbinsSpsz,nbinsSpsptD,nbinsSpsy,nbinseta,nbinsphi,nbinsSpsptjet,nbinsSpsy,nbinseta,nbinsphi,nbinsArea,nbinNtracks,100,100,nbinsMult};
        const Double_t RMmin[nRMAxis]= {zlims[0],ptDlims[0],etalims[0],etalims[0],philims[0],ptjetlims[0],etalims[0],etalims[0],philims[0],0,nTracksLims[0],zlims[0],ptDlims[0],etalims[0],etalims[0],philims[0],ptjetlims[0],etalims[0],etalims[0],philims[0],0,nTracksLims[0],-1,-1,multLims[0]};
        const Double_t RMmax[nRMAxis]= {zlims[1],ptDlims[1],etalims[1],etalims[1],philims[1],ptjetlims[1],etalims[1],etalims[1],philims[1],1,nTracksLims[1],zlims[1],ptDlims[1],etalims[1],etalims[1],philims[1],ptjetlims[1],etalims[1],etalims[1],philims[1],1,nTracksLims[1],1,1,multLims[1]};
        fResponseMatrix = new THnSparseF("ResponseMatrix","z, p_{T}^{D}, y^{D}, #eta^{D}, #varphi^{D},p_{T}^{jet},y^{jet}, #eta^{jet}, #varphi^{jet}, Area^{jet}, Ntracks^{jet}: Rec and Gen, p_{T}^{res}, z^{res},mult",nRMAxis,RMnbins,RMmin,RMmax);
        fResponseMatrix->Sumw2();
        fOutput->Add(fResponseMatrix);
    }

   PostData(1, fOutput);

   return kTRUE;
}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t z, Double_t ptD, Double_t ptj, Double_t jetEta, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent, Int_t pdgTrue,Double_t JetPhi,Double_t JetNTracks,Double_t JetY,Double_t JetArea){


    Float_t nTracklets = static_cast<Float_t>(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
    // multiplicity estimator with VZERO
    Float_t vzeroMult=0;
    AliAODVZERO *vzeroAOD = (AliAODVZERO*)aodEvent->GetVZEROData();
    if(vzeroAOD) vzeroMult = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();
    Float_t multiplicity = nTracklets; // set to the Ntracklet estimator
    if(fMultiplicityEstimator==kVZERO) { multiplicity = vzeroMult; }

    Float_t lPercentile = -1.;
    if(!fIsPPData) {
        lPercentile = fCuts->GetCentrality(aodEvent);
    }

   Double_t masses[2]={0.,0.};
   Int_t pdgdaughtersD0[2]={211,321};//pi,K
   Int_t pdgdaughtersD0bar[2]={321,211};//K,pi
   Int_t pdgMeson = 413;
   if (fCandidateType == kD0toKpi) pdgMeson = 421;

   masses[0]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
   masses[1]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar

   Double_t *point= nullptr;

   if(!fUseMCInfo)
   {
      point=new Double_t[14];
      point[0]=z;
      point[1]=ptD;
      point[2]=candidate->Y(pdgMeson);
      point[3]=candidate->Eta();
      point[4]=candidate->Phi();
      point[5]=ptj;
      point[6]=JetY;
      point[7]=jetEta;
      point[8]=JetPhi;
      point[9]=JetArea;
      point[10]=JetNTracks;
      point[11]=masses[0];
      point[12]=static_cast<Double_t>(multiplicity);
      point[13]=static_cast<Double_t>(lPercentile);
   }
   else
   {
      point=new Double_t[15];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=candidate->Y(pdgMeson);
      point[5]=jetEta;
      point[6]=-999;
      point[7]=-999;
      point[8]=static_cast<Double_t>(IsBkg ? 1 : 0);
      point[9]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[10]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
      point[11]=JetPhi;
      point[12]=JetNTracks;
      point[13]=static_cast<Double_t>(multiplicity);
      point[14]=static_cast<Double_t>(lPercentile);

   }
    Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);

    if(!fUseMCInfo) {
      if(isselected==1 || isselected==3)
      {
          fhInvMassptD->Fill(masses[0],ptD);
          point[11]=masses[0];
          fhsDphiz->Fill(point,1.);
      }
      if(isselected>=2)
      {
          fhInvMassptD->Fill(masses[1],ptD);
          point[11]=masses[1];
          fhsDphiz->Fill(point,1.);
      }
    }
    else {
      if(isselected==1 || isselected==3)
      {
          fhInvMassptD->Fill(masses[0],ptD);
          point[3]=masses[0];

          if(isselected==1) {
            point[6]=masses[0];
            point[7]=-999;
          }
          else if(isselected==3 && pdgTrue==421){
            point[6]=masses[0];
            point[7]=masses[1];
          }
          else if(isselected==3 && pdgTrue==-421){
            point[6]=masses[1];
            point[7]=masses[0];
          }
          fhsDphiz->Fill(point,1.);
      }
      if(isselected>=2)
      {
          fhInvMassptD->Fill(masses[1],ptD);
          point[3]=masses[1];
          if(isselected==2) {
            point[6]=masses[1];
            point[7]=-999;
          }
          else if(isselected==3){
            point[6]=-999;
            point[7]=-999;
          }
          fhsDphiz->Fill(point,1.);
      }

    }

    delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar,  Double_t z, Double_t ptD, Double_t ptj, Double_t jetEta, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){

    AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
    Double_t deltamass= dstar->DeltaInvMass();
    Int_t pdgMeson = 413;
    if (fCandidateType == kD0toKpi) pdgMeson = 421;

    fhPtPion->Fill(softpi->Pt());

    Double_t *point= nullptr;
    if(!fAnalyseDBkg)
    {
        point=new Double_t[6];
        point[0]=z;
        point[1]=ptj;
        point[2]=ptD;
        point[3]=deltamass;
        point[4]=dstar->Y(pdgMeson);
        point[5]=jetEta;
    }
    else
    {
        point=new Double_t[9];
        point[0]=z;
        point[1]=ptj;
        point[2]=ptD;
        point[3]=deltamass;
        point[4]=dstar->Y(pdgMeson);
        point[5]=jetEta;
        point[6]=static_cast<Double_t>(IsBkg ? 1 : 0);
        point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
        point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
    }

    if(!point){
        AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
        return;
    }


    fhInvMassptD->Fill(deltamass,ptD);
    fhsDphiz->Fill(point,1.);
    delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskDJetCorrelationsQA::FillHistogramsMCGenDJetCorr(Double_t z,Double_t ptD,Double_t ptjet, Double_t yD, Double_t jetEta, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){

    Double_t pdgmass=0;
    if(fCandidateType==kD0toKpi) pdgmass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    if(fCandidateType==kDstartoKpipi) pdgmass=TDatabasePDG::Instance()->GetParticle(413)->Mass()-TDatabasePDG::Instance()->GetParticle(421)->Mass();

    Double_t *point=nullptr;

    if(fNAxesBigSparse==6){
        point=new Double_t[6];
        point[0]=z;
        point[1]=ptjet;
        point[2]=ptD;
        point[3]=pdgmass;
        point[4]=yD;
        point[5]=jetEta;
    }
    if(fNAxesBigSparse==9){
        point=new Double_t[9];
        point[0]=z;
        point[1]=ptjet;
        point[2]=ptD;
        point[3]=pdgmass;
        point[4]=yD;
        point[5]=jetEta;
        point[6]=1;
        point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
        point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
    }

    if(!point){
        AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
        return;
    }


    point[3]=pdgmass;
    fhsDphiz->Fill(point,1.);

    delete[] point;
}

//_______________________________________________________________________________

Float_t AliAnalysisTaskDJetCorrelationsQA::DeltaR(AliEmcalJet *p1, AliVParticle *p2, Double_t rho) const {
   //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
   //It recalculates the eta-phi values if it was asked for background subtraction of the jets
   if(!p1 || !p2) return -1;

    Double_t phi1=p1->Phi(),eta1=p1->Eta();

    if (rho>0)
    {
            Double_t pj[3];
            Bool_t okpj=p1->PxPyPz(pj);
            if(!okpj){
                printf("Problems getting momenta\n");
                return -999;
            }
            pj[0] = p1->Px() - p1->Area()*(rho*TMath::Cos(p1->AreaPhi()));
            pj[1] = p1->Py() - p1->Area()*(rho*TMath::Sin(p1->AreaPhi()));
            pj[2] = p1->Pz() - p1->Area()*(rho*TMath::SinH(p1->AreaEta()));
            //Image of the function Arccos(px/pt) where pt = sqrt(px*px+py*py) is:
            //[0,pi]    if py > 0 and
            //[pi,2pi]  if py < 0
            phi1 = TMath::ACos(pj[0]/TMath::Sqrt(pj[0]*pj[0]+pj[1]*pj[1]));
            if(pj[1]<0) phi1 = 2*TMath::Pi() - phi1;
            eta1 = TMath::ASinH(pj[2]/TMath::Sqrt(pj[0]*pj[0]+pj[1]*pj[1]));
    }

   Double_t phi2 = p2->Phi(),eta2 = p2->Eta() ;

   Double_t dPhi=TMath::Abs(phi1-phi2);
   if(dPhi>TMath::Pi()) dPhi = (2*TMath::Pi())-dPhi;


   Double_t dEta=eta1-eta2;
   Double_t deltaR=TMath::Sqrt(dEta*dEta + dPhi*dPhi );
   return deltaR;

}
//_______________________________________________________________________________
Float_t AliAnalysisTaskDJetCorrelationsQA::CheckDeltaR(AliEmcalJet *p1, AliVParticle *p2) const {
    //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
    if(!p1 || !p2) return -1;

    Double_t phi1=p1->Phi(),eta1=p1->Eta();

    Double_t phi2 = p2->Phi(),eta2 = p2->Eta() ;

    Double_t dPhi=TMath::Abs(phi1-phi2);
    if(dPhi>TMath::Pi()) dPhi = (2*TMath::Pi())-dPhi;


    Double_t dEta=eta1-eta2;
    Double_t deltaR=TMath::Sqrt(dEta*dEta + dPhi*dPhi );
    return deltaR;

}

//_______________________________________________________________________________

Int_t AliAnalysisTaskDJetCorrelationsQA::IsDzeroSideBand(AliAODRecoCascadeHF *candDstar){

   Double_t ptD=candDstar->Pt();
   Int_t bin = fCuts->PtBin(ptD);
   if (bin < 0){
      // /PWGHF/vertexingHF/AliRDHFCuts::PtBin(Double_t) const may return values below zero depending on config.
      bin = 9999; // void result code for coverity (bin later in the code non-zero) - will coverity pick up on this?
      return -1;  // out of bounds
   }

   Double_t invM=candDstar->InvMassD0();
   Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();

   Float_t fourSigmal= mPDGD0-4.*fSigmaD0[bin] , sixSigmal= mPDGD0-8.*fSigmaD0[bin];
   Float_t fourSigmar= mPDGD0+4.*fSigmaD0[bin] , sixSigmar= mPDGD0+8.*fSigmaD0[bin];

   if((invM>=sixSigmal && invM<fourSigmal) || (invM>fourSigmar && invM<=sixSigmar)) return 1;
   else return 0;

}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskDJetCorrelationsQA::InEMCalAcceptance(AliVParticle *vpart){
   //check eta phi of a VParticle: return true if it is in the EMCal acceptance, false otherwise

   Double_t phiEMCal[2]={1.405,3.135},etaEMCal[2]={-0.7,0.7};
   Bool_t binEMCal=kTRUE;
   Double_t phi=vpart->Phi(), eta=vpart->Eta();
   if(phi<phiEMCal[0] || phi>phiEMCal[1]) binEMCal=kFALSE;
   if(eta<etaEMCal[0] || eta>etaEMCal[1]) binEMCal=kFALSE;
   return binEMCal;


}
