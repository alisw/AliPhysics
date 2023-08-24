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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to acquire MC-level predictions in general
// for several LF-related particle species. First deployed to deal with
// the Pb-Pb 5 TeV strangeness analysis. Adapted to several other use cases
// afterwards, including prompt/non-prompt HF
//
// Please report any bugs, complaints, suggestions to:
// --- david.dobrigkeit.chinellato@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TDatabasePDG.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
#include "AliPWG0Helper.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskMCPredictionsEE.h"
////
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVertexingHFUtils.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictionsEE)

AliAnalysisTaskMCPredictionsEE::AliAnalysisTaskMCPredictionsEE()
: AliAnalysisTaskSE(),
fkSelectINELgtZERO(kTRUE),
fTree(0),
fEtaThreshold(8.),
fEnergyThreshold(0.),
fEtaBarrel(0.5),
fNMPI(0),
fMC_NPart(0),
fMC_NColl(0),
fMC_b(0.),
fInelGT0(0),
fNchEta(0),
fNLambdaEta(0),
fNXiEta(0),
fNAntiXiEta(0),
fNOmegaEta(0),
fNPiEta(0),
fNPi0Eta(0),
fNKchEta(0),
fNK0Eta(0),
fEffEnergy(0)
{
  DefineOutput(1, TTree::Class());
}

AliAnalysisTaskMCPredictionsEE::AliAnalysisTaskMCPredictionsEE(const char *name)
: AliAnalysisTaskSE(name),
fkSelectINELgtZERO(kTRUE),
fTree(0),
fEtaThreshold(8.),
fEnergyThreshold(0.),
fEtaBarrel(0.5),
fNMPI(0),
fMC_NPart(0),
fMC_NColl(0),
fMC_b(0.),
fInelGT0(0),
fNchEta(0),
fNLambdaEta(0),
fNXiEta(0),
fNAntiXiEta(0),
fNOmegaEta(0),
fNPiEta(0),
fNPi0Eta(0),
fNKchEta(0),
fNK0Eta(0),
fEffEnergy(0)
{
  DefineOutput(1, TTree::Class());
}


AliAnalysisTaskMCPredictionsEE::~AliAnalysisTaskMCPredictionsEE()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictionsEE::UserCreateOutputObjects()
{
  // Called once
  fDB = TDatabasePDG::Instance();

  // create tree
  fTree = new TTree("fTree","fTree");
  fTree->Branch("nmpi", &fNMPI, "nmpi/I");
  //fMC_NPart
  //fMC_NColl
  //fMC_b
  fTree->Branch("nchEta", &fNchEta, "nchEta/I");
  fTree->Branch("nLambdaEta", &fNLambdaEta, "nLambdaEta/I");
  fTree->Branch("nXiEta", &fNXiEta, "nXiEta/I");
  fTree->Branch("nAntiXiEta", &fNAntiXiEta, "nAntiXiEta/I");
  fTree->Branch("nOmegaEta", &fNOmegaEta, "nOmegaEta/I");
  fTree->Branch("nPiEta", &fNPiEta, "nPiEta/I");
  fTree->Branch("nPi0Eta", &fNPi0Eta, "nPi0Eta/I");
  fTree->Branch("nKchEta", &fNKchEta, "nKchEta/I");
  fTree->Branch("nK0Eta", &fNK0Eta, "nK0Eta/I");
  fTree->Branch("effEnergy",&fEffEnergy, "effEnergy/F");

  //TTree: Normal
  PostData(1, fTree);
  
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskMCPredictionsEE::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  
  
  // Connect to the InputEvent
  // After these lines, we should have an ESD/AOD event + the number of V0s in it.
  
  // Appropriate for ESD analysis!
  
  lMCevent = MCEvent();
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  
  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  
  //------------------------------------------------
  // Multiplicity Information Acquistion
  //------------------------------------------------
  
  //Monte Carlo Level information !
  //--------- GENERATED NUMBER OF CHARGED PARTICLES
  // ---> Variable Definition
  
  int dNdeta = 0;
  bool lEvSel_INELgtZEROStackPrimaries = kFALSE;

  // reset variables
  fNchEta = 0;
  fNLambdaEta = 0;
  fNXiEta = 0;
  fNAntiXiEta = 0;
  fNOmegaEta = 0;
  fNPiEta = 0;
  fNPi0Eta = 0;
  fNKchEta = 0;
  fNK0Eta = 0;
  fEffEnergy = 0.;
  
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  {   // This is the begining of the loop on tracks
    TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
    if(!particleOne) continue;
    if(!particleOne->GetPDG()) continue;
    Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
    if(TMath::Abs(lThisCharge)<0.001) continue;
    if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
    
    //    Double_t gpt = particleOne -> Pt();
    Double_t geta = particleOne -> Eta();
    
    if( TMath::Abs(geta) < 0.5 ) dNdeta++;
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZEROStackPrimaries = kTRUE;
  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------
  
  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZEROStackPrimaries && fkSelectINELgtZERO ) return;
    
  //------------------------------------------------
  // Acquire information on Npart, Ncoll, b
  //------------------------------------------------
  
  //Npart and Ncoll information
  AliGenHijingEventHeader* hHijing=0;
  AliGenDPMjetEventHeader* dpmHeader=0;
  AliGenEventHeader* mcGenH = lMCevent->GenEventHeader();
  
  Int_t fMC_NPart = -1;
  Int_t fMC_NColl = -1;
  Float_t fMC_b = -1;
  Int_t fMC_NMPI = -1;
  
  if (mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())){
    AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (mcGenH);
    if(fMcPythiaHeader){
      fNMPI = fMcPythiaHeader->GetNMPI();
    }
  }
  
  //DPMJet/HIJING info if available
  if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class()))
    hHijing = (AliGenHijingEventHeader*)mcGenH;
  else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
    TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
    hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing"));
    if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing pPb_0"));
    if (!hHijing) hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("Hijing_0"));
  }
  else if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
    dpmHeader = (AliGenDPMjetEventHeader*)mcGenH;
  }
  if(hHijing)   {
    fMC_NPart = hHijing->ProjectileParticipants()+hHijing->TargetParticipants();
    fMC_NColl = hHijing->NN()+hHijing->NNw()+hHijing->NwN()+hHijing->NwNw();
  }
  if(dpmHeader) {
    fMC_NPart =dpmHeader->ProjectileParticipants()+dpmHeader->TargetParticipants();
    fMC_NColl =dpmHeader->NN()+dpmHeader->NNw()+dpmHeader->NwN()+dpmHeader->NwNw();
  }
  
  //check EPOS info, if available
  if ( IsEPOSLHC() ){
    AliGenHepMCEventHeader *lHepMCHeader = 0x0;
    if (mcGenH->InheritsFrom(AliGenHepMCEventHeader::Class()))
      lHepMCHeader = (AliGenHepMCEventHeader*)mcGenH;
    
    if (lHepMCHeader ){
      fMC_NPart = lHepMCHeader->Npart_proj()+lHepMCHeader->Npart_targ();
      fMC_NColl = lHepMCHeader->N_Nwounded_collisions() +
      lHepMCHeader->Nwounded_N_collisions() +
      lHepMCHeader->Nwounded_Nwounded_collisions();
      
      fMC_b = lHepMCHeader->impact_parameter();
    }
  }
  
  Int_t lThisPDG  = 0;
  Double_t lThisRap  = 0;
  Double_t lThisPt   = 0;
  Bool_t lIsPhysicalPrimary = kFALSE;

  loopMC(lMCevent); 

  fTree->Fill();
  
  // Post output data.
  PostData(1, fTree);
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictionsEE::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fTree = dynamic_cast<TTree*> (GetOutputData(0));
  if (!fTree) {
    Printf("ERROR: fTree not available");
    return;
  }

  system("touch ok.job");
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictionsEE::IsHijing() const {
  //Function to check if this is Hijing MC
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())){
      //Option 1: Just Hijing
      lReturnValue = kTRUE;
    } else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
      //Option 2: cocktail involving Hijing
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      TIter next(headers);
      while (const TObject *obj=next()){
        //Look for an object inheriting from the hijing header class
        if ( obj->InheritsFrom(AliGenHijingEventHeader::Class()) ){ lReturnValue = kTRUE; }
      }
    }
  }
  return lReturnValue;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictionsEE::IsDPMJet() const {
  //Function to check if this is DPMJet
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    if (mcGenH->InheritsFrom(AliGenDPMjetEventHeader::Class())) {
      //DPMJet Header is there!
      lReturnValue = kTRUE;
    }
  }
  return lReturnValue;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictionsEE::IsEPOSLHC() const {
  //Function to check if this is DPMJet
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    //A bit uncivilized, but hey, if it works...
    TString lHeaderTitle = mcGenH->GetName();
    if (lHeaderTitle.Contains("EPOSLHC")) {
      //This header has "EPOS" in its title!
      lReturnValue = kTRUE;
    }
  }
  return lReturnValue;
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictionsEE::loopMC(AliMCEvent *mcEvent){
  // # of reco tracks and primaries from MC
  Int_t ntra = mcEvent->GetNumberOfTracks();

  //  fMaxChargePt = 0;
  fEffEnergy = 0;
  
  // reset variables
  //  fP_cand_leadA=0;
  //  fN_cand_leadA=0;
  //  fP_cand_leadC=0;
  //  fN_cand_leadC=0;

  //  fNch=0,
  fNchEta=0,
    //fNchEtaA=0,fNchEtaC=0,
    fNLambdaEta=0,fNXiEta=0,fNAntiXiEta=0,
    //fNXiEtaFrag=0,fNXiEtaUp=0,fNXiEtaDown=0,
    fNOmegaEta=0,fNPiEta=0,fNPi0Eta=0,fNKchEta=0,fNK0Eta=0;
  //  fSumPtLambdaEta=fSumPtXiEta=fSumPtOmegaEta=fSumPtPiEta=0;
  //fEnergyEta=0;
  
  // loop MC particles
  Int_t nPrim = mcEvent->Stack()->GetNprimary();
  for (Int_t i = 0; i < nPrim; i++){
    TParticle *part = mcEvent->Stack()->Particle(i);
    if (!AliPWG0Helper::IsPrimaryCharged(part, nPrim)) continue;
    Double_t eta = part->Eta();
    if (fabs(eta) < 1.0) fInelGT0 = true;
  }

  Float_t px,py,pz,pt;
  for(Int_t i = 0; i < ntra ;i++){
    // get particle from stack
    AliMCParticle *MCpart = (AliMCParticle *) mcEvent->GetTrack(i);
    TParticle *part = MCpart->Particle(); 
    
    Int_t status = (part->GetStatusCode() == 1);
    
    // get momentum
    px = part->Px();
    py = part->Py();
    pz = part->Pz();
    pt = sqrt(px*px + py*py);
    
    Int_t charge = 0;

    if(fDB->GetParticle(part->GetPdgCode())) charge = Int_t(fDB->GetParticle(part->GetPdgCode())->Charge());
    else continue; // skip particles with undefined pdg

    if(charge) {
      status = AliPWG0Helper::IsPrimaryCharged(part, nPrim); // official definition of charged primary
    }
    else{
      if (part->GetFirstDaughter() != -1 && part->GetFirstDaughter() < nPrim) status = 0;
    }

    if(status && charge){ // charged particles
      //      fNch++;

      if(TMath::Abs(part->Eta())<fEtaBarrel){
	//	fEnergyEta+=part->Energy();
	fNchEta++;
	//        if(part->Eta() < 0) fNchEtaA++;
	//        else fNchEtaC++;

	//	if(pt > fMaxChargePt) fMaxChargePt = pt;
	
	if(TMath::Abs(part->GetPdgCode()) == 211){ // charged pions
	  fNPiEta++;
	  //	  fSumPtPiEta += pt;
	}
        else if(TMath::Abs(part->GetPdgCode()) == 321){ // neutral kaons
          fNKchEta++;
        }
      }
    }
    else if(! charge){
      if(TMath::Abs(part->Eta())<fEtaBarrel){
        if(TMath::Abs(part->GetPdgCode()) == 111){ // neutral pions
          fNPi0Eta++;
        }
        else if(TMath::Abs(part->GetPdgCode()) == 311){ // K0s neutral kaons
          fNK0Eta++;
        }
      }
    }

    if(TMath::Abs(part->Eta())<fEtaBarrel){
      if(TMath::Abs(part->GetPdgCode()) == 3122){ // lambda
	fNLambdaEta++;
	//	fSumPtLambdaEta += pt;
      }
      if(TMath::Abs(part->GetPdgCode()) == 3312){ // Xi
        if(part->GetPdgCode() > 0){
	  //          fPtXiEta[fNXiEta] = pt;
          fNXiEta++;
        }
        else{
	  //          fPtXiEta[fNAntiXiEta] = pt;
          fNAntiXiEta++;
        }
	//	fSumPtXiEta += pt;

        // look if it comes from fragmenetation
	/*
        Int_t imoth = part->GetFirstMother();
        AliMCParticle *partMCM;
        TParticle *partM;
        while(imoth >= 0 && fNXiEtaFrag < 100){
          partMCM = (AliMCParticle *) mcEvent->GetTrack(imoth);  
          partM = partMCM->Particle();
          if(TMath::Abs(partM->GetPdgCode()) == 3){
            fPtXiEtaFrag[fNXiEtaFrag] = pt;
            fNXiEtaFrag++;
            imoth = -1;
          }
          else if(TMath::Abs(partM->GetPdgCode()) == 1){
            fPtXiEtaUp[fNXiEtaUp] = pt;
            fNXiEtaUp++;
            imoth = -1;
          }
          else if(TMath::Abs(partM->GetPdgCode()) == 2){
            fPtXiEtaDown[fNXiEtaDown] = pt;
            fNXiEtaDown++;
            imoth = -1;
          }
          else{
            imoth = partM->GetFirstMother();
          }
	}
	*/
      }
      if(TMath::Abs(part->GetPdgCode()) == 3334){ // Omega
	fNOmegaEta++;
	//	fSumPtOmegaEta += pt;
     }
    }
    
    // keep only stable particles at forward pseudorapidity
    if(status!=1 || TMath::Abs(part->Eta())<fEtaThreshold || part->Energy() < fEnergyThreshold) continue;
    
    Int_t imoth = part->GetFirstMother();
    
    Int_t statusM = 0;
    Int_t pdgM = 0;
    
    AliMCParticle *partMCM;
    if(imoth>7){ // first 8 usually refered to initial incident protons
      partMCM = (AliMCParticle *) mcEvent->GetTrack(imoth);  
      TParticle *partM = partMCM->Particle();
      
      statusM = partM->GetStatusCode();
      
      if(statusM==11 &&  TMath::Abs(partM->GetPdgCode())>99) pdgM = partM->GetPdgCode();
      else statusM = 0;
    }
    
    if(statusM != 11) imoth = -1; // only decay accepted to define mother

    fEffEnergy -= part->Energy();

    /*
    if((part->Eta())>0.0) { // positive Z (A side)
      if(charge && fP_cand_leadA < fgkDim){ // charged candidate
	fE_p_cand_leadA[fP_cand_leadA]=part->Energy();
	fPdg_cand_leadP2[fP_cand_leadA]=part->GetPdgCode();
	fPdgM_cand_leadP2[fP_cand_leadA]=pdgM;
	fLabel_cand_P2[fP_cand_leadA]=i;
	fLabelM_cand_P2[fP_cand_leadA]=imoth;
	fMomCP2x[fP_cand_leadA] = px;
	fMomCP2y[fP_cand_leadA] = py;
	fMomCP2z[fP_cand_leadA] = pz;
	fP_cand_leadA++;
      } 
      else if(!charge && fN_cand_leadA < fgkDim){ // neutral candidate
	fE_n_cand_leadA[fN_cand_leadA]=part->Energy();
	fPdg_cand_leadN2[fN_cand_leadA]=part->GetPdgCode(); 
	fPdgM_cand_leadN2[fN_cand_leadA]=pdgM;
	fLabel_cand_N2[fN_cand_leadA]=i;
	fLabelM_cand_N2[fN_cand_leadA]=imoth;
	fMomCN2x[fN_cand_leadA] = px;
	fMomCN2y[fN_cand_leadA] = py;
	fMomCN2z[fN_cand_leadA] = pz;
	fN_cand_leadA++;
      }
    }
    else{ // negative Z (C side)
      if(charge && fP_cand_leadC < fgkDim){ // charged candidate
	fE_p_cand_leadC[fP_cand_leadC]=part->Energy(); 
	fPdg_cand_leadP1[fP_cand_leadC]=part->GetPdgCode(); 
	fPdgM_cand_leadP1[fP_cand_leadC]=pdgM;
	fLabel_cand_P1[fP_cand_leadC]=i; 
	fLabelM_cand_P1[fP_cand_leadC]=imoth;
	fMomCP1x[fP_cand_leadC] = px;
	fMomCP1y[fP_cand_leadC] = py;
	fMomCP1z[fP_cand_leadC] = pz;
	fP_cand_leadC++;
      }
      else if(!charge && fN_cand_leadC < fgkDim){ // neutral candidate
	fE_n_cand_leadC[fN_cand_leadC]=part->Energy(); 
	fPdg_cand_leadN1[fN_cand_leadC]=part->GetPdgCode(); 
	fPdgM_cand_leadN1[fN_cand_leadC]=pdgM;
	fLabel_cand_N1[fN_cand_leadC]=i; 
	fLabelM_cand_N1[fN_cand_leadC]=imoth;
	fMomCN1x[fN_cand_leadC] = px;
	fMomCN1y[fN_cand_leadC] = py;
	fMomCN1z[fN_cand_leadC] = pz;
	fN_cand_leadC++;
      }
    }
    */
  } 
}
