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
//  Analysis Taks for reconstructed particle correlation
//  (first implementation done for D mesons) with jets
//  (use the so called Emcal framework)

//  Migrated from PWGJE/FlavourJetTasks to PWGHF/jetsHF
//
//  Dev task. Copied from the PWGHF/jetsHF/AliAnalysisTaskDJetCorrelations
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// S. Antônio (University of São Paulo) antonio.silva@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL)  XMZhang@lbl.gov
// B. Trzeciak (Utrecht University) barbara.antonina.trzeciak@cern.ch
// A. Mohanty (Utrecht University) auro.mohanty@cern.ch
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include "TROOT.h"
#include <THnSparse.h>
#include <TSystem.h>
#include <TObjectTable.h>
#include "AliMultSelection.h"

#include "AliATDJetCorrDev.h"
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

ClassImp(AliATDJetCorrDev)


//_______________________________________________________________________________
AliATDJetCorrDev::AliATDJetCorrDev() :
AliAnalysisTaskEmcalJet("AliATDJetCorrDev",kTRUE),
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
fJetRadius(0.4),
fRecluster(kFALSE),
fhLPThetaEnergy(),
fMinLnOneByTheta(),
fMaxLnOneByTheta(),
fMinEnergy(),
fMaxEnergy(),
fhDmesonOrNot(),
fMinDelMass(),
fMaxDelMass()
{
   //
   // Default ctor
}

//_______________________________________________________________________________
AliATDJetCorrDev::AliATDJetCorrDev(const Char_t* name, AliRDHFCuts* cuts, ECandidateType candtype) :
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
fJetRadius(0.4),
fRecluster(kFALSE),
fhLPThetaEnergy(),
fMinLnOneByTheta(),
fMaxLnOneByTheta(),
fMinEnergy(),
fMaxEnergy(),
fhDmesonOrNot(),
fMinDelMass(),
fMaxDelMass()
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //

   Info("AliATDJetCorrDev","Calling Constructor");
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
AliATDJetCorrDev::~AliATDJetCorrDev() {
   //
   // destructor
   //

   Info("~AliATDJetCorrDev","Calling Destructor");

   delete fCuts;

}

//_______________________________________________________________________________
void AliATDJetCorrDev::Init(){
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
void AliATDJetCorrDev::UserCreateOutputObjects() {
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
Bool_t AliATDJetCorrDev::Run()
{
    // user exec from AliAnalysisTaskEmcal is used

    // Load the event
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

    // check #1
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<=0) {
        // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
        return kFALSE;
    }

    TClonesArray *arrayDmesontoDaughters=nullptr;

    // 1. using AOD: Finding the above TClonesArray* arrayDmesontoDaughters
    if (!aodEvent && AODEvent() && IsStandardAOD()) {

        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());

        // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
        // have to be taken from the AOD event held by the AliAODExtension
        AliAODHandler* aodHandler = (AliAODHandler*)
        ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if(aodHandler->GetExtensions()) {
            AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
            AliAODEvent *aodFromExt = ext->GetAOD();
            arrayDmesontoDaughters=(TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
        }
    }
    else if(aodEvent){
        arrayDmesontoDaughters=(TClonesArray*)aodEvent->GetList()->FindObject(fBranchName.Data());
    }

    // check #2: if the TClonesArray not found.
    if (!arrayDmesontoDaughters) {
        AliInfo(Form("Could not find array %s, skipping the event",fBranchName.Data()));
        //  return;
    }
    else AliDebug(2, Form("Found %d vertices",arrayDmesontoDaughters->GetEntriesFast()));

    // 2. array for MC
    TClonesArray* mcArray = nullptr;
    if (fUseMCInfo) { //not used at the moment,uncomment return if you use
        mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        if (!mcArray) {
            printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
        }
    }

    //D meson candidates. Also background if is MC

    // 3. candidate array. default, false.
    if(fUseCandArray)
    {
        fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
        if (!fCandidateArray) return kFALSE;
        for(Int_t icand=0; icand<fCandidateArray->GetEntriesFast(); icand++)
        {
            fhstat->Fill(2);
        }
    }
    // 4. side band array. default, false.
    if(fUseSBArray)
    {
        fSideBandArray = dynamic_cast<TClonesArray*>(GetInputData(2));
        if (!fSideBandArray) return kFALSE;
    }

    fhstat->Fill(0); // number of events analysed.

    // 5. bug fix
    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return kFALSE;

    //Event selection
    Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
    TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses(); //unused variable firedTriggerClasses
    if(!iseventselected) return kFALSE;

    fhstat->Fill(1); // number of events selected.
    Float_t lPercentile = -1.;
    if(!fUseMCInfo && !fIsPPData) {
        lPercentile = fCuts->GetCentrality(aodEvent);
        fhCentDjet->Fill(lPercentile);
    }

    // 6. 
    // for MC response matrix of efficiency studies, fMultCand option only
    if(fUseMCInfo && fBuildRMEff){
    //----------------------------
    
        AliJetContainer* mcjets = nullptr;
        if(!fAnalyseDBkg) mcjets = GetJetContainer(1);
        else mcjets = GetJetContainer(2);
        if(!mcjets) return kFALSE;
    
        AliParticleContainer *MCParticlesCont = mcjets->GetParticleContainer();
    
        mcjets->ResetCurrentID();
        AliEmcalJet* jet=nullptr;
    
        // while loop over MC jets
        // a) ghost jets skipped while iterating over all jets from jet container
        // b) tracks of each jet have their properties stored in histograms
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
    
            // total number of tracks in a jet
            Int_t ntrjet =  jet->GetNumberOfTracks();
    
            // for loop on jet tracks
            for(Int_t itrk = 0; itrk < ntrjet; itrk++)
            {
                AliAODMCParticle* jetTrk=(AliAODMCParticle*)jet->TrackAt(itrk,MCParticlesCont->GetArray());
                if (!jetTrk) continue;
                fhPtJetTrks->Fill(jetTrk->Pt());
                fhPhiJetTrks->Fill(jetTrk->Phi());
                fhEtaJetTrks->Fill(jetTrk->Eta());
            }//end loop on jet tracks
    
        }// end while loop over MC jets
    
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
    
    // end if(fUseMCInfo && fBuildRMEff) for MC response matrix of efficiency studies, fMultCand option only
    }
    //----------------------------
    // 6. (else) 
    //
    else {
        AliJetContainer* JetCont = GetJetContainer(0); // contains all the jets
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

        JetCont->ResetCurrentID(); // resets internal jet iterator
        AliEmcalJet* jet=nullptr;

        // 6... looping over jets 
        while ((jet = JetCont->GetNextJet()))
        {
            UInt_t rejectionReason = 0;
            Bool_t OKjet = JetCont->AcceptJet(jet, rejectionReason);
            if(!OKjet) {
                fhstat->Fill(5); // jets rejected
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
        }

        if(ParticlesCont->GetNParticles()>0) fhstat->Fill(2); // ??

        // Finally, find the jets with a D meson inside
        if(fCorrelationMethod==kConstituent) // constituent method finds if D is a constituent of the jet
        {
            if(ParticlesCont->GetNParticles()>0) ConstituentCorrelationMethod(kFALSE,aodEvent);
            if(fAnalyseDBkg==kTRUE && ParticlesContSB->GetNParticles()>0) ConstituentCorrelationMethod(kTRUE,aodEvent);
        }

        else if(fCorrelationMethod==kAngular) // angular method finds if D is within the cone of the jet
        {
            if(fCandidateArray->GetEntriesFast()>0) AngularCorrelationMethod(kFALSE,aodEvent);
            if(fAnalyseDBkg==kTRUE && fSideBandArray->GetEntriesFast()>0) AngularCorrelationMethod(kTRUE,aodEvent);
        }
    }
    // end if(fUseMCInfo && fBuildRMEff) else
    // not for MC response matrix of efficiency studies, fMultCand option only
    // but for data, or MC without response matrix studies


   PostData(1,fOutput);
   return kTRUE;
}
//_______________________________________________________________________________
void AliATDJetCorrDev::ConstituentCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent)
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
    if(jet && fRecluster)
    {
        // Pseudo code:
        //  1. Recluster the jet
        //  2. Start declustering
        //  2.1. Get two subjets of the jet
        //  2.2. Get and fill in a histogram, the theta and E of the soft subjet/radiator as, you follow the hard subjet
        //
        // ClusterSequence runs the jet clustering. see https://github.com/alisw/fastjet/blob/master/fastjet/example/01-basic.cc line 71
        fastjet::ClusterSequence* cs = Recluster(jet); // used to give warning:did you mean 'fastjet::contrib::Recluster'? (FixIt) 
        if (cs)
        {
            std::vector<fastjet::PseudoJet> recl_jets = sorted_by_pt( cs->inclusive_jets() );
            if( recl_jets.size() > 0 )
            {
               DeclusterTheJet( recl_jets[0], jet);            // 2. Declustering the jet
               FillLundPlane(4.0, 5.2);                        // 3. Fill the Lund Plane
            }
        }
        delete cs;
        
        // The FastJet::contrib way
        //fastjet::PseudoJet recl_jet = ReclusteredJet(jet);  // 1. getting the reclustered jet
        //// the AliEmcal jet is provided to compare the last hard jet of recl_jet with D meson from AliEmcal jet
        //DeclusterTheJet(recl_jet, jet);                     // 2. Declustering the jet
        //FillLundPlane(4.0, 5.2);
    }

}
//_______________________________________________________________________________
void AliATDJetCorrDev::AngularCorrelationMethod(Bool_t IsBkg, AliAODEvent* aodEvent)
{
    AliJetContainer* JetCont = GetJetContainer(0);

    Int_t ncand = 0;
    if(!IsBkg) ncand = fCandidateArray->GetEntriesFast();
    else ncand = fSideBandArray->GetEntriesFast();

    Double_t rho = 0;
    if(!JetCont->GetRhoName().IsNull()) rho = JetCont->GetRhoVal();

    for(Int_t icand = 0; icand<ncand; icand++)
    {
        AliVParticle* charm=nullptr;
        if(!IsBkg) charm = (AliVParticle*)fCandidateArray->At(icand);
        else charm = (AliVParticle*)fSideBandArray->At(icand);
        if(!charm) continue;

        Int_t JetTag = AliEmcalJet::kD0;
        if (fCandidateType == kDstartoKpipi) JetTag = AliEmcalJet::kDStar;
        //loop over jets
        JetCont->ResetCurrentID();
        AliEmcalJet* jet=nullptr;
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
void AliATDJetCorrDev::CreateResponseMatrix(AliEmcalJet* jet)
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
        Double_t etaRec = jet->Eta();
        Double_t etaGen = MCjet->Eta();
        Double_t DYRec = -999;
        AliAODMCParticle* DTrk=(AliAODMCParticle*)Dgen;
        Double_t DYGen = DTrk->Y();
        Double_t pTRes = JetPtGen ? ((JetPtRec - JetPtGen) / JetPtGen) : -999;
        Double_t zRes = zGen ? ((zRec - zGen) / zGen) : -999;

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

        Double_t fillRM[13] = {zRec,JetPtRec,Drec->Pt(),DYRec,etaRec,zGen,JetPtGen,Dgen->Pt(),DYGen,etaGen,pTRes,zRes,-1};
        fResponseMatrix->Fill(fillRM,1.);
    }
}

//_______________________________________________________________________________
void AliATDJetCorrDev::CreateMCResponseMatrix(AliEmcalJet* MCjet, AliAODEvent* aodEvent)
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
    Double_t zGen = Z(Dgen,MCjet,0);
    Double_t JetPtGen = MCjet->Pt();
    Double_t JetEtaGen = MCjet->Eta();
    Double_t DPtGen = Dgen->Pt();
    //Double_t JetnTrkGen = MCjet->GetNumberOfTracks();//unused variable
    AliAODMCParticle* DTrk = (AliAODMCParticle*)Dgen;
    Double_t DYGen = DTrk->Y();
    Int_t pdg = DTrk->GetPdgCode();

    Double_t zRec = -999;
    Double_t JetPtRec = -999;
    Double_t JetEtaRec = -999;
    Double_t DPtRec = -999;
    Double_t DYRec = -999;
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
            JetPtRec = jet->Pt() - rho*jet->Area();
            JetEtaRec = jet->Eta();
            DPtRec = Drec->Pt();
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

            if(fCandidateType==kD0toKpi)
            {
                AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Drec;
                FillHistogramsD0JetCorr(dzero,zRec,Drec->Pt(),JetPtRec,JetEtaRec,JetnTrkRec,kFALSE,bDInEMCalAcc,bJetInEMCalAcc,aodEvent,pdg);
            }
            else if(fCandidateType==kDstartoKpipi)
            {
                AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Drec;
                FillHistogramsDstarJetCorr(dstar,zRec,Drec->Pt(),JetPtRec,JetEtaRec,kFALSE,bDInEMCalAcc,bJetInEMCalAcc);
            }

        } // if HF reco jet
    } // if jet cont reco

    Double_t fillRM[13] = {zRec,JetPtRec,DPtRec,DYRec,JetEtaRec,zGen,JetPtGen,DPtGen,DYGen,JetEtaGen,pTRes,zRes,multiplicity};
    fResponseMatrix->Fill(fillRM,1.);

}

//_______________________________________________________________________________
void AliATDJetCorrDev::FillDJetHistograms(AliEmcalJet* jet, Double_t rho, Bool_t IsBkg, AliAODEvent* aodEvent)
{

    AliVParticle *Dmeson = jet->GetFlavourTrack(0);
    Double_t JetPtCorr = jet->Pt() - rho*jet->Area();
    Double_t JetEtaRec = jet->Eta();
    Double_t JetnTrkRec = jet->GetNumberOfTracks();
    Double_t z = 0;
    if(rho>0) z = Z(Dmeson,jet,rho);
    else z = Z(Dmeson,jet);

    if(IsBkg==kFALSE && fBuildRM==kTRUE) CreateResponseMatrix(jet);

    Bool_t bDInEMCalAcc=InEMCalAcceptance(Dmeson);
    Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

    if(fCandidateType==kD0toKpi)
    {
        AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)Dmeson;
        FillHistogramsD0JetCorr(dzero,z,Dmeson->Pt(),JetPtCorr,JetEtaRec,JetnTrkRec,IsBkg,bDInEMCalAcc,bJetInEMCalAcc,aodEvent,-999);
    }
    else if(fCandidateType==kDstartoKpipi)
    {
        AliAODRecoCascadeHF* dstar=(AliAODRecoCascadeHF*)Dmeson;
        FillHistogramsDstarJetCorr(dstar,z,Dmeson->Pt(),JetPtCorr,JetEtaRec,IsBkg,bDInEMCalAcc,bJetInEMCalAcc);
    }


}

//_______________________________________________________________________________
void AliATDJetCorrDev::GetHFJet(AliEmcalJet*& jet, Bool_t IsBkg)
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
void AliATDJetCorrDev::FindMCJet(AliEmcalJet*& mcjet)
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
    if(!HFMCjet) mcjet=nullptr;
}

//_______________________________________________________________________________
void AliATDJetCorrDev::Terminate(Option_t*)
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
void  AliATDJetCorrDev::SetMassLimits(Double_t range, Int_t pdg){
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
void  AliATDJetCorrDev::SetMassLimits(Double_t lowlimit, Double_t uplimit){
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
Bool_t AliATDJetCorrDev::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
   if(nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   if(!width) return kFALSE;
   for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
   return kTRUE;
}

//_______________________________________________________________________________
Double_t AliATDJetCorrDev::Z(AliVParticle* part,AliEmcalJet* jet, Double_t rho) const{

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
Double_t AliATDJetCorrDev::Z(AliVParticle* part,AliEmcalJet* jet) const{

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
Double_t AliATDJetCorrDev::Z(Double_t* p, Double_t *pj) const{
    //p is momentum array of D meson
    //pj is momentum array of jet

    Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1]+pj[2]*pj[2];
    Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet2);
    return z;
}


//_______________________________________________________________________________
Double_t AliATDJetCorrDev::ZT(Double_t* p, Double_t *pj) const{

    Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1];
    Double_t z=(p[0]*pj[0]+p[1]*pj[1])/(pjet2);
    return z;
}

//_______________________________________________________________________________
Bool_t  AliATDJetCorrDev::DefineHistoForAnalysis(){

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

    fhCentDjet=new TH1F("hCentDjet","Centrality",100,0,100);
    fOutput->Add(fhCentDjet);

    const Int_t nbinsmass=300;
    const Int_t nbinspttrack=500;
    const Int_t nbinsptjet=200;
    const Int_t nbinsptD=100;
    const Int_t nbinsphi=200;
    const Int_t nbinseta=100;

    //binning for THnSparse
    const Int_t nbinsSpsmass=280;
    const Int_t nbinsSpsptjet=200;
    const Int_t nbinsSpsptD=100;
    const Int_t nbinsSpsz=160;
    const Int_t nbinsSpsy=150;
    const Int_t nbinNtracks=50;

    //binning for Lund Plane
    const Int_t nbinstheta=300;
    const Int_t nbinsenergy=300;
    const Int_t nbinsdelmass=300;
    //bin ranges for reclustering histograms
    fMinDelMass = -2.0;
    fMaxDelMass =  2.0;
    fMinLnOneByTheta = 0;
    fMaxLnOneByTheta = 100;
    fMinEnergy = 0;
    fMaxEnergy = 100;

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
    fhPtJet     = new TH1F("hPtJet",  "Jet Pt distribution; p_{T} (GeV/c)",nbinsptjet,ptjetlims[0],ptjetlims[1]);
    fhPtJet->Sumw2();

    // reclutering/Lund plane related histograms
    fhLPThetaEnergy = new TH2F("hLPThetaEnergy","Lund Plane in theta, energy of the radiator",nbinstheta,fMinLnOneByTheta,fMaxLnOneByTheta,nbinsenergy,fMinEnergy,fMaxEnergy);
    fhLPThetaEnergy->SetStats(kTRUE);
    fhLPThetaEnergy->GetYaxis()->SetTitle("ln (1/#theta)");
    fhLPThetaEnergy->GetXaxis()->SetTitle("E");
    fhLPThetaEnergy->Sumw2();

    fhDmesonOrNot   = new TH1F("fhDmesonOrNot","Inv mass (D meson - hard track)",  nbinsdelmass,fMinDelMass,fMaxDelMass);
    fhDmesonOrNot->SetStats(kTRUE);

    fOutput->Add(fhLPThetaEnergy);
    fOutput->Add(fhDmesonOrNot);

    fOutput->Add(fhPhiJetTrks);
    fOutput->Add(fhEtaJetTrks);
    fOutput->Add(fhPtJetTrks);
    fOutput->Add(fhPhiJet);
    fOutput->Add(fhEtaJet);
    fOutput->Add(fhPtJet);

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
            const Int_t nAxis=14;
            const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,nbinsSpsy,nbinsSpsy,nbinsSpsmass,nbinsSpsmass, nbinNtracks,nbinsMult,nbinsCent,2, 2, 2};
            const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,etalims[0],etalims[0],fMinMass,fMinMass,nTracksLims[0],multLims[0],centLims[0], -0.5,-0.5,-0.5};
            const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,etalims[1],etalims[1],fMaxMass,fMaxMass,nTracksLims[1],multLims[1],centLims[1], 1.5, 1.5 , 1.5};
            fNAxesBigSparse=nAxis;
            fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass., y^{D}, #eta^{jet}, mass_{true}, mass_{refl}, nTrk, mult, cent, Bkg?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);
        }

    }
    else{
        AliInfo("Creating a 9 axes container");
        const Int_t nAxis=9;
        const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,nbinsSpsy,nbinsSpsy,nbinNtracks,nbinsMult,nbinsCent};
        const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,etalims[0],etalims[0],nTracksLims[0],multLims[0],centLims[0]};
        const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,etalims[1],etalims[1],nTracksLims[1],multLims[1],centLims[1]};
        fNAxesBigSparse=nAxis;
        fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass., y^{D}, #eta^{jet}, nTrk, mult, cent,", nAxis, nbinsSparse, minSparse, maxSparse);
    }

    if(!fhsDphiz) AliFatal("No THnSparse created");
    fhsDphiz->Sumw2();

    fOutput->Add(fhsDphiz);

    if(fBuildRM == kTRUE || fBuildRMEff == kTRUE)
    {

        const Int_t nRMAxis=13;
        const Int_t RMnbins[nRMAxis] = {nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsy,nbinsSpsy,nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsy,nbinsSpsy,100,100,nbinsMult};
        const Double_t RMmin[nRMAxis]={zlims[0],ptjetlims[0],ptDlims[0],etalims[0],etalims[0],zlims[0],ptjetlims[0],ptDlims[0],etalims[0],etalims[0],-1,-1,multLims[0]};
        const Double_t RMmax[nRMAxis]={zlims[1],ptjetlims[1],ptDlims[1],etalims[1],etalims[1],zlims[1],ptjetlims[1],ptDlims[1],etalims[1],etalims[1],1,1,multLims[1]};
        fResponseMatrix = new THnSparseF("ResponseMatrix","z, p_{T}^{jet}, p_{T}^{D}, y^{D}, #eta^{jet}: Rec and Gen, p_{T}^{res}, z^{res}",nRMAxis,RMnbins,RMmin,RMmax);
        fResponseMatrix->Sumw2();
        fOutput->Add(fResponseMatrix);
    }

    PostData(1, fOutput);

    return kTRUE;
}

//_______________________________________________________________________________
void AliATDJetCorrDev::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t z, Double_t ptD, Double_t ptj, Double_t jetEta, Double_t nJetTrack, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent, Int_t pdgTrue){


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

   Double_t *point=nullptr;

   if(!fUseMCInfo)
   {
      point=new Double_t[9];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=candidate->Y(pdgMeson);
      point[5]=jetEta;
      point[6]=nJetTrack;
      point[7]=multiplicity;
      point[8]=lPercentile;
   }
   else
   {
      point=new Double_t[14];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=candidate->Y(pdgMeson);
      point[5]=jetEta;
      point[6]=-999;
      point[7]=-999;
      point[8]=nJetTrack;
      point[9]=multiplicity;
      point[10]=lPercentile;
      point[11]=static_cast<Double_t>(IsBkg ? 1 : 0);
      point[12]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[13]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);

   }
    Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);

    if(!fUseMCInfo) {
      if(isselected==1 || isselected==3)
      {
          fhInvMassptD->Fill(masses[0],ptD);
          point[3]=masses[0];
          fhsDphiz->Fill(point,1.);
      }
      if(isselected>=2)
      {
          fhInvMassptD->Fill(masses[1],ptD);
          point[3]=masses[1];
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
void AliATDJetCorrDev::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar,  Double_t z, Double_t ptD, Double_t ptj, Double_t jetEta, Bool_t IsBkg, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){

    AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
    Double_t deltamass= dstar->DeltaInvMass();
    Int_t pdgMeson = 413;
    if (fCandidateType == kD0toKpi) pdgMeson = 421;

    fhPtPion->Fill(softpi->Pt());

    Double_t *point=nullptr;
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
void AliATDJetCorrDev::FillHistogramsMCGenDJetCorr(Double_t z,Double_t ptD,Double_t ptjet, Double_t yD, Double_t jetEta, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){

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
Float_t AliATDJetCorrDev::DeltaR(AliEmcalJet *p1, AliVParticle *p2, Double_t rho) const {
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
Float_t AliATDJetCorrDev::CheckDeltaR(AliEmcalJet *p1, AliVParticle *p2) const {
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
Int_t AliATDJetCorrDev::IsDzeroSideBand(AliAODRecoCascadeHF *candDstar){

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
Bool_t AliATDJetCorrDev::InEMCalAcceptance(AliVParticle *vpart){
   //check eta phi of a VParticle: return true if it is in the EMCal acceptance, false otherwise

   Double_t phiEMCal[2]={1.405,3.135},etaEMCal[2]={-0.7,0.7};
   Bool_t binEMCal=kTRUE;
   Double_t phi=vpart->Phi(), eta=vpart->Eta();
   if(phi<phiEMCal[0] || phi>phiEMCal[1]) binEMCal=kFALSE;
   if(eta<etaEMCal[0] || eta>etaEMCal[1]) binEMCal=kFALSE;
   return binEMCal;

}

////_______________________________________________________________________________
//fastjet::PseudoJet AliATDJetCorrDev::ReclusteredJet(const AliEmcalJet* jet){
//    fastjet::Recluster recluster_ca_inf(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R); // the function that reclusters with the given algorithm and jet-radius (max possible in this case)
//    fastjet::PseudoJet rec_jet_ca_inf;// = recluster_ca_inf(jet); // the given jet is reclustered now, and stored as a PseudoJet
//    return (rec_jet_ca_inf);
//}
//    
//_______________________________________________________________________________
fastjet::ClusterSequence* AliATDJetCorrDev::Recluster(const AliEmcalJet* jet){
    // reclustering a fastjet jet with CA algorithm
    // ref(arxiv/1111/6097, pg. 13, section 3.3: fastjet::ClusterSequence)
    // ClusterSequence object is created to run the jet clustering
    // it needs:    1. const std::vector<L> & input_particles,
    //              2. const JetDefintion & jet_def
    // where input_particles is the vector of initial particles of any type ( PseudoJet, HepLorentzVector etc.) that can be used to initialize a PseudoJet
    // and jet_def contains the full specification of the clustering.
    
    std::vector<fastjet::PseudoJet> input_particles;
    UShort_t ntracks = jet->GetNumberOfTracks();
    for (Int_t j = 0; j < ntracks; j++)
    {
        input_particles.push_back( fastjet::PseudoJet( jet->Track(j)->Px(), jet->Track(j)->Py(), jet->Track(j)->Pz(), jet->Track(j)->E()) );
    }
    fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, fJetRadius, fastjet::E_scheme, fastjet::Best);
    return (new fastjet::ClusterSequence(input_particles, jet_def) );
}

//_______________________________________________________________________________
void AliATDJetCorrDev::DeclusterTheJet(fastjet::PseudoJet fj_jet, AliEmcalJet* ali_jet)
{
    fastjet::PseudoJet jj = fj_jet;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;

    while(jj.has_parents(j1, j2))
    {
        if(j1.perp2() < j2.perp2()) std::swap(j1,j2); // j1 is assumed to be harder subjet
                                                    // if j1 is not harder of the two, swap it with j2
                                                    // so that j1 is now the harder subjet
                                                    //
        Double_t del_R = j1.delta_R(j2);            // find the angular distance between the two
                                                    // subjets j1 and j2
                                                    //
        Double_t rad_E = j2.E();                    // the radiator energy
        // Fill the 2D histogram with del_R and del_E
        FillLundPlane(del_R, rad_E);
        jj = j1;
    }
    // check stats here if the final jet `jj' is a D meson or not
    // put the stats into a histogram
    // // AliVParticle* hardTrack = ... maybe cast pseudojet to alivparticle and then get the flavour track of this jet. compare these two tracks
    // // or compare their invariant masses, or their transverse momenta
    // Doublt_t hardDMass = jj.perp() - Dmeson->Pt();
    //AliVParticle *Dmeson = jet->GetFlavourTrack(0);
    //
    fastjet::PseudoJet hardTrk=jj;
    AliVParticle* Dmeson = ali_jet->GetFlavourTrack(0);
    Double_t deltaM = (Dmeson->M() - hardTrk.m());
    fhDmesonOrNot->Fill(deltaM);

}

//_______________________________________________________________________________
void AliATDJetCorrDev::FillLundPlane(Double_t sj_deltaR, Double_t sj_energy)
{
    fhLPThetaEnergy->Fill(TMath::Log(1.0/sj_deltaR), sj_energy);
    
}


