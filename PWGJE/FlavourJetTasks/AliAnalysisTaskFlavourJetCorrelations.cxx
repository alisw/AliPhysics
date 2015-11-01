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
//
//-----------------------------------------------------------------------
// Authors:
// C. Bianchin (Utrecht University) chiara.bianchin@cern.ch
// S. Antônio (University of São Paulo) antonio.silva@cern.ch
// A. Grelli (Utrecht University) a.grelli@uu.nl
// X. Zhang (LBNL)  XMZhang@lbl.gov
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include "TROOT.h"
#include <TH3F.h>
#include <THnSparse.h>
#include <TSystem.h>
#include <TObjectTable.h>

#include "AliAnalysisTaskFlavourJetCorrelations.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliPicoTrack.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRhoParameter.h"
#include "AliParticleContainer.h"

ClassImp(AliAnalysisTaskFlavourJetCorrelations)


//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations() :
AliAnalysisTaskEmcalJet("",kTRUE),
fUseMCInfo(kTRUE), 
fUseReco(kTRUE),
fCandidateType(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(0),
fMinMass(),
fMaxMass(),  
fJetArrName(0),
fTrackArrName(0),
fCandArrName(0),
fJetRadius(0),
fCandidateArray(0),
fSideBandArray(0),
fPmissing(),
fPtJet(0),
fRhoValue(0),
fIsDInJet(0),
fTypeDInJet(2),
fTrackArr(0),
fSwitchOnSB(0),
fSwitchOnPhiAxis(0),
fSwitchOnOutOfConeAxis(0),
fSwitchOnSparses(1),
fNAxesBigSparse(9),
fJetCont(0),
fTrackCont(0),
fClusterCont(0),
fhstat(),
fhCentDjet(),
fhPtJetTrks(),
fhPhiJetTrks(),
fhEtaJetTrks(),
fhPtJet(),
fhPhiJet(),
fhEtaJet(),
fhNtrArr(),
fhNJetPerEv(),
fhdeltaRJetTracks(),
fhInvMassptD(),
fhDiffSideBand(),
fhInvMassptDbg(),
fhPtPion(),
fhControlDInJ(),
fhIDddaugh(),
fhIDddaughOut(),
fhIDjetTracks(),
fhDRdaughOut(),
fhzDinjet(),
fhmissingp(),
fhMissPi(),
fhDeltaPtJet(),
fhRelDeltaPtJet(),
fhzDT(),
fhDeltaRD(),
fhDeltaRptDptj(),
fhDeltaRptDptjB(),
fhsDphiz()

{
   //
   // Default ctor
   //
   //SetMakeGeneralHistograms(kTRUE)(),
   
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelations::AliAnalysisTaskFlavourJetCorrelations(const Char_t* name, AliRDHFCuts* cuts,ECandidateType candtype) :
AliAnalysisTaskEmcalJet(name,kTRUE),
fUseMCInfo(kTRUE),
fUseReco(kTRUE),  
fCandidateType(),
fPDGmother(),
fNProngs(),
fPDGdaughters(),
fBranchName(),
fCuts(0),
fMinMass(),
fMaxMass(),  
fJetArrName(0),
fTrackArrName(0),
fCandArrName(0),
fJetRadius(0),
fCandidateArray(0),
fSideBandArray(0),
fPmissing(),
fPtJet(0),
fRhoValue(0),
fIsDInJet(0),
fTypeDInJet(2),
fTrackArr(0),
fSwitchOnSB(0),
fSwitchOnPhiAxis(0),
fSwitchOnOutOfConeAxis(0),
fSwitchOnSparses(1),
fNAxesBigSparse(9),
fJetCont(0),
fTrackCont(0),
fClusterCont(0),
fhstat(),
fhCentDjet(),
fhPtJetTrks(),
fhPhiJetTrks(),
fhEtaJetTrks(),
fhPtJet(),
fhPhiJet(),
fhEtaJet(),
fhNtrArr(),
fhNJetPerEv(),
fhdeltaRJetTracks(),
fhInvMassptD(),
fhDiffSideBand(),
fhInvMassptDbg(),
fhPtPion(),
fhControlDInJ(),
fhIDddaugh(),
fhIDddaughOut(),
fhIDjetTracks(),
fhDRdaughOut(),
fhzDinjet(),
fhmissingp(),
fhMissPi(),
fhDeltaPtJet(),
fhRelDeltaPtJet(),
fhzDT(),
fhDeltaRD(),
fhDeltaRptDptj(),
fhDeltaRptDptjB(),
fhsDphiz()
{
   //
   // Constructor. Initialization of Inputs and Outputs
   //
   
   Info("AliAnalysisTaskFlavourJetCorrelations","Calling Constructor");
   fCuts=cuts;
   fCandidateType=candtype;
   const Int_t nptbins=fCuts->GetNPtBins();
   Float_t defaultSigmaD013[13]={0.012, 0.012, 0.012, 0.015, 0.015,0.018,0.018,0.020,0.020,0.030,0.030,0.037,0.040};
   
   switch(fCandidateType){
   case 0:
      fPDGmother=421;
      fNProngs=2;
      fPDGdaughters[0]=211;//pi 
      fPDGdaughters[1]=321;//K
      fPDGdaughters[2]=0; //empty
      fPDGdaughters[3]=0; //empty
      fBranchName="D0toKpi";
      fCandArrName="D0";
      break;
   case 1: 
      fPDGmother=413;
      fNProngs=3;
      fPDGdaughters[1]=211;//pi soft
      fPDGdaughters[0]=421;//D0
      fPDGdaughters[2]=211;//pi fromD0
      fPDGdaughters[3]=321; // K from D0
      fBranchName="Dstar";
      fCandArrName="Dstar";
      
      if(nptbins<=13){
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
   if(fCandidateType==kDstartoKpipi) SetMassLimits(0.015, fPDGmother);
   
      DefineInput(1, TClonesArray::Class());
      DefineInput(2, TClonesArray::Class());
      
      DefineOutput(1,TList::Class()); // histos
      DefineOutput(2,AliRDHFCuts::Class()); // my cuts
   
}

//_______________________________________________________________________________

AliAnalysisTaskFlavourJetCorrelations::~AliAnalysisTaskFlavourJetCorrelations() {
   //
   // destructor
   //
   
   Info("~AliAnalysisTaskFlavourJetCorrelations","Calling Destructor");  
   
   delete fCuts;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::Init(){
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

void AliAnalysisTaskFlavourJetCorrelations::UserCreateOutputObjects() { 
   // output 
   Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
   AliAnalysisTaskEmcal::UserCreateOutputObjects();
   
   fJetCont = GetJetContainer(0);
   if(fJetCont){
      fTrackCont =   fJetCont->GetParticleContainer();
      fClusterCont = fJetCont->GetClusterContainer();
   }

   
   // define histograms
   // the TList fOutput is already defined in  AliAnalysisTaskEmcal::UserCreateOutputObjects()
   DefineHistoForAnalysis();
   PostData(1,fOutput);
   
   return;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::Run()
{
   // user exec from AliAnalysisTaskEmcal is used
    
   // Load the event
   AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
   
   TClonesArray *arrayDStartoD0pi=0;
   
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
   
   TClonesArray* mcArray = 0x0;
   if (fUseMCInfo) { //not used at the moment,uncomment return if you use
      mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) {
      	 printf("AliAnalysisTaskSEDStarSpectra::UserExec: MC particles not found!\n");
      	 //return kFALSE;
      }
   }
   
   //retrieve jets
   //this is a duplication of fTrackCont, but is is used in the loop of line 598 and changing it needs a thorough test 
   fTrackArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackArrName));
   //clusArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("CaloClustersCorr"));
   //fJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrName));
   //fJetContainer=GetJetContainer(0);
   //if(!fJetContainer) {
   //   AliError("Jet Container 0 not found");
   //   return kFALSE;
   //}
   fJetRadius=fJetCont->GetJetRadius();
   
   if(!fTrackArr){
      AliInfo(Form("Could not find the track array\n"));
      return kFALSE;
   }
   
    
   fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
   if (!fCandidateArray) return kFALSE;
   if ((fCandidateType==1 && fSwitchOnSB) || fUseMCInfo) {
      fSideBandArray = dynamic_cast<TClonesArray*>(GetInputData(2));
      if (!fSideBandArray) return kFALSE;
   }
   //Printf("ncandidates found %d",fCandidateArray->GetEntriesFast());
   
   fhstat->Fill(0);
   
   // fix for temporary bug in ESDfilter 
   // the AODs with null vertex pointer didn't pass the PhysSel
   if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return kFALSE;
   
   //Event selection
   Bool_t iseventselected=fCuts->IsEventSelected(aodEvent);
   TString firedTriggerClasses=((AliAODEvent*)aodEvent)->GetFiredTriggerClasses();
   if(!iseventselected) return kFALSE;
   
   fhstat->Fill(1);
   fhCentDjet->Fill(fCent);
   
   //retrieve charm candidates selected
   Int_t candidates = fCandidateArray->GetEntriesFast();
   Int_t njets=fJetCont->GetNJets();
    //Background Subtraction for jets
    //If there's no background subtraction rhoval=0 and momentum correction is simply not took into account
    if (!fJetCont->GetRhoName().IsNull()) fRhoValue = fJetCont->GetRhoVal();
   
   Int_t ntrarr=fTrackCont->GetNParticles();
   fhNtrArr->Fill(ntrarr);
   
   for(Int_t i=0;i<ntrarr;i++){
      AliAODTrack *jtrack=static_cast<AliAODTrack*>(fTrackCont->GetParticle(i));
      //reusing histograms
      if(!jtrack || jtrack->IsMuonTrack()) continue; // added check on IsMuonTrack because in some cases the DCA() gives crazy  YAtDCA values that issue floating point exception 
      fhPtJetTrks->Fill(jtrack->Pt());
      fhPhiJetTrks->Fill(jtrack->Phi());
      fhEtaJetTrks->Fill(jtrack->Eta());

   }
   
   Int_t cntjet=0;
   //loop over jets
   for (Int_t iJets = 0; iJets<njets; iJets++) {
      fPmissing[0]=0;
      fPmissing[1]=0;
      fPmissing[2]=0;
      
      //Printf("Jet N %d",iJets);
      AliEmcalJet* jet = (AliEmcalJet*)fJetCont->GetJet(iJets);
      //Printf("Corr task Accept Jet");
      if(!AcceptJet(jet)) {
      	 fhstat->Fill(5);
      	 continue;
      }
      
      //jets variables
      fPtJet = jet->Pt() - jet->Area()*fRhoValue; //It takes into account the background subtraction
      Double_t origPtJet=fPtJet;
      
      //used for normalization
      fhstat->Fill(3);
      cntjet++;
      // fill energy, eta and phi of the jet
      fhPhiJet ->Fill(jet->Phi());
      fhEtaJet ->Fill(jet->Eta());
      fhPtJet  ->Fill(fPtJet);
      
      //loop on jet particles
      Int_t ntrjet=  jet->GetNumberOfTracks(); 
      Double_t sumPtTracks=0,sumPzTracks=0;
      Int_t zg1jtrk=0;
      for(Int_t itrk=0;itrk<ntrjet;itrk++){
      	 
      	 AliPicoTrack* jetTrk=(AliPicoTrack*)jet->TrackAt(itrk,fTrackArr);     
      	 fhdeltaRJetTracks->Fill(DeltaR(jet,jetTrk));
      	 Double_t ztrackj=Z(jetTrk,jet);
     	 sumPtTracks+=jetTrk->Pt(); 
     	 sumPzTracks+=jetTrk->Pz(); 
     	 if(ztrackj>1){
     	    zg1jtrk++;
     	 }
     	 
     	 
     	 
      }//end loop on jet tracks
      


      	 //Printf("N candidates %d ", candidates);
      	 for(Int_t ic = 0; ic < candidates; ic++) {
      	    fhstat->Fill(7);
      	    // D* candidates
      	    AliVParticle* charm=0x0;
      	    charm=(AliVParticle*)fCandidateArray->At(ic);
      	    if(!charm) continue;
      	    AliAODRecoDecayHF *charmdecay=(AliAODRecoDecayHF*) charm;
      	    fIsDInJet=IsDInJet(jet, charmdecay, kTRUE);
      	    if (fIsDInJet) FlagFlavour(jet);
      	    if (jet->TestFlavourTag(AliEmcalJet::kDStar) || jet->TestFlavourTag(AliEmcalJet::kD0)) fhstat->Fill(4);
      	    
      	    //Note: the z component of the jet momentum comes from the eta-phi direction of the jet particles, it is not calculated from the z component of the tracks since, as default, the scheme used for jet reco is the pt-scheme which sums the scalar component, not the vectors. Addind the D daughter momentum component by componet as done here is not 100% correct, but the difference is small, for fairly collimated particles.

      	    Double_t pjet[3];
      	    jet->PxPyPz(pjet);
             //It corrects the jet momentum if it was asked for jet background subtraction
             pjet[0] = jet->Px() - jet->Area()*(fRhoValue*TMath::Cos(jet->AreaPhi()));
             pjet[1] = jet->Py() - jet->Area()*(fRhoValue*TMath::Sin(jet->AreaPhi()));
             pjet[2] = jet->Pz() - jet->Area()*(fRhoValue*TMath::SinH(jet->AreaEta()));
             
            //It corrects the jet momentum due to D daughters out of the jet
      	    RecalculateMomentum(pjet,fPmissing);      	          	    
      	    fPtJet=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]); //recalculated jet pt
      	    if((jet->Pt()-jet->Area()*fRhoValue)<0) fPtJet = -fPtJet;
      	    
      	    //debugging histograms
      	    Double_t pmissing=TMath::Sqrt(fPmissing[0]*fPmissing[0]+fPmissing[1]*fPmissing[1]+fPmissing[2]*fPmissing[2]); //recalculated jet total momentum
      	    for(Int_t k=0;k<3;k++) fhMissPi[k]->Fill(fPmissing[k]);
      	    
      	    fhmissingp->Fill(pmissing);
      	    Double_t ptdiff=fPtJet-origPtJet;
      	    fhDeltaPtJet->Fill(ptdiff);
      	    fhRelDeltaPtJet->Fill(ptdiff/origPtJet);
      	    
      	    FillHistogramsRecoJetCorr(charm, jet, aodEvent);
      	    
      	 }//end loop on candidates
      	 
      	 //retrieve side band background candidates for Dstar
      	 Int_t nsbcand = 0; if(fSideBandArray) nsbcand=fSideBandArray->GetEntriesFast();
      	 
      	 for(Int_t ib=0;ib<nsbcand;ib++){
      	    if(fCandidateType==kDstartoKpipi && !fUseMCInfo){
      	       AliAODRecoCascadeHF *sbcand=(AliAODRecoCascadeHF*)fSideBandArray->At(ib);
      	       if(!sbcand) continue;
      	        
      	       fIsDInJet=IsDInJet(jet, sbcand,kFALSE);
      	       Double_t pjet[3];
      	       jet->PxPyPz(pjet);
                //It corrects the jet momentum if it was asked for jet background subtraction
                pjet[0] = jet->Px() - jet->Area()*(fRhoValue*TMath::Cos(jet->AreaPhi()));
                pjet[1] = jet->Py() - jet->Area()*(fRhoValue*TMath::Sin(jet->AreaPhi()));
                pjet[2] = jet->Pz() - jet->Area()*(fRhoValue*TMath::SinH(jet->AreaEta()));
                
               //It corrects the jet momentum due to D daughters out of the jet
      	       RecalculateMomentum(pjet,fPmissing);      	          	    
      	       fPtJet=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]); //recalculated jet pt
      	       if((jet->Pt()-jet->Area()*fRhoValue)<0) fPtJet = -fPtJet;
     	       SideBandBackground(sbcand,jet);
     	       
      	    }
      	    if(fUseMCInfo){
      	       
      	       AliAODRecoDecayHF* charmbg = 0x0;
      	       charmbg=(AliAODRecoDecayHF*)fSideBandArray->At(ib);
      	       if(!charmbg) continue;
      	       fhstat->Fill(8);
      	       fIsDInJet=IsDInJet(jet, charmbg,kFALSE);
      	       if (fIsDInJet) {
      	       	  FlagFlavour(jet); //this are backgroud HF jets, but flagged as signal at the moment. Can use the bkg flavour flag in the future. This info is not stored now a part in the jet
      	       	  fhstat->Fill(9);
      	       }
      	       Double_t pjet[3];
      	       jet->PxPyPz(pjet);
      	       //It corrects the jet momentum if it was asked for jet background subtraction
      	       pjet[0] = jet->Px() - jet->Area()*(fRhoValue*TMath::Cos(jet->AreaPhi()));
      	       pjet[1] = jet->Py() - jet->Area()*(fRhoValue*TMath::Sin(jet->AreaPhi()));
      	       pjet[2] = jet->Pz() - jet->Area()*(fRhoValue*TMath::SinH(jet->AreaEta()));
                
               //It corrects the jet momentum due to D daughters out of the jet
      	       RecalculateMomentum(pjet,fPmissing);      	          	    
      	       fPtJet=TMath::Sqrt(pjet[0]*pjet[0]+pjet[1]*pjet[1]); //recalculated jet pt
      	       if((jet->Pt()-jet->Area()*fRhoValue)<0) fPtJet = -fPtJet;
      	       MCBackground(charmbg,jet);
      	    }
      	 }

   } // end of jet loop
   
   fhNJetPerEv->Fill(cntjet);
   PostData(1,fOutput);
   return kTRUE;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::Terminate(Option_t*)
{    
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
   Info("Terminate"," terminate");
   AliAnalysisTaskSE::Terminate();
   
   fOutput = dynamic_cast<TList*> (GetOutputData(1));
   if (!fOutput) {     
      printf("ERROR: fOutput not available\n");
      return;
   }
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourJetCorrelations::SetMassLimits(Double_t range, Int_t pdg){
   Float_t mass=0;
   Int_t abspdg=TMath::Abs(pdg);
   
   mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
   // compute the Delta mass for the D*
   if(fCandidateType==kDstartoKpipi){
      Float_t mass1=0;
      mass1=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      mass = mass-mass1;
   }
   
   fMinMass = mass-range;
   fMaxMass = mass+range;
   
   AliInfo(Form("Setting mass limits to %f, %f",fMinMass,fMaxMass));
   if (fMinMass<0 || fMaxMass<=0 || fMaxMass<fMinMass) AliFatal("Wrong mass limits!\n");
}

//_______________________________________________________________________________

void  AliAnalysisTaskFlavourJetCorrelations::SetMassLimits(Double_t lowlimit, Double_t uplimit){
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

Bool_t AliAnalysisTaskFlavourJetCorrelations::SetD0WidthForDStar(Int_t nptbins,Float_t *width){
   if(nptbins>30) {
      AliInfo("Maximum number of bins allowed is 30!");
      return kFALSE;
   }
   if(!width) return kFALSE;
   for(Int_t ipt=0;ipt<nptbins;ipt++) fSigmaD0[ipt]=width[ipt];
   return kTRUE;
}

//_______________________________________________________________________________

Double_t AliAnalysisTaskFlavourJetCorrelations::Z(AliVParticle* part,AliEmcalJet* jet, Bool_t transverse) const{
   if(!part || !jet){
      printf("Particle or jet do not exist!\n");
      return -99;
   }
   Double_t p[3],pj[3];
   Bool_t okpp=part->PxPyPz(p);
   Bool_t okpj=jet->PxPyPz(pj);
    
    //Background Subtraction
    //It corrects the each component of the jet momentum for Z calculation
 
    pj[0] = jet->Px() - jet->Area()*(fRhoValue*TMath::Cos(jet->AreaPhi()));
    pj[1] = jet->Py() - jet->Area()*(fRhoValue*TMath::Sin(jet->AreaPhi()));
    pj[2] = jet->Pz() - jet->Area()*(fRhoValue*TMath::SinH(jet->AreaEta()));

    
    
   if(!okpp || !okpj){
      printf("Problems getting momenta\n");
      return -999;
   }
   
   RecalculateMomentum(pj, fPmissing);
   if(transverse) return ZT(p,pj);
   else
   return Z(p,pj);
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelations::Z(Double_t* p, Double_t *pj) const{
   
   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1]+pj[2]*pj[2];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1]+p[2]*pj[2])/(pjet2);
   return z;
}


//_______________________________________________________________________________
Double_t AliAnalysisTaskFlavourJetCorrelations::ZT(Double_t* p, Double_t *pj) const{
   
   Double_t pjet2=pj[0]*pj[0]+pj[1]*pj[1];
   Double_t z=(p[0]*pj[0]+p[1]*pj[1])/(pjet2);
   return z;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::RecalculateMomentum(Double_t* pj, const Double_t* pmissing) const {

   pj[0]+=pmissing[0];
   pj[1]+=pmissing[1];
   pj[2]+=pmissing[2];

}

//_______________________________________________________________________________

Bool_t  AliAnalysisTaskFlavourJetCorrelations::DefineHistoForAnalysis(){
   
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
   Int_t nbinsptjet=150;
   if(!fJetCont->GetRhoName().IsNull()) nbinsptjet=200;
   const Int_t nbinsptD=100;
   const Int_t nbinsz=144;
   const Int_t nbinsphi=200;
   const Int_t nbinseta=100;
   
   //binning for THnSparse
   const Int_t nbinsSpsmass=60;
   const Int_t nbinsSpsptjet=200;
   const Int_t nbinsSpsptD=50;
   const Int_t nbinsSpsz=144;
   const Int_t nbinsSpsphi=100;
    
   const Float_t pttracklims[2]={0.,200.};
   Float_t ptjetlims[2]={0.,150.};
   if(!fJetCont->GetRhoName().IsNull()) ptjetlims[0]=-50.;
   const Float_t ptDlims[2]={0.,50.};
   const Float_t zlims[2]={0.,1.2};
   const Float_t philims[2]={0.,6.3};
   const Float_t etalims[2]={-1.5,1.5};
   
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
   
   fhdeltaRJetTracks=new TH1F("hdeltaRJetTracks","Delta R of tracks in the jets",200, 0.,10.);
   fhdeltaRJetTracks->Sumw2();
   
   fhNtrArr= new TH1F("hNtrArr", "Number of tracks in the array of jets; number of tracks",500,0,1000);
   fhNtrArr->Sumw2();
   
   fhNJetPerEv=new TH1F("hNJetPerEv","Number of jets used per event; number of jets/ev",10,-0.5,9.5);
   fhNJetPerEv->Sumw2();
   
   
   fOutput->Add(fhPhiJetTrks);
   fOutput->Add(fhEtaJetTrks);
   fOutput->Add(fhPtJetTrks);
   fOutput->Add(fhPhiJet);
   fOutput->Add(fhEtaJet);
   fOutput->Add(fhPtJet);
   fOutput->Add(fhdeltaRJetTracks);
   fOutput->Add(fhNtrArr);
   fOutput->Add(fhNJetPerEv);
   


      
      //debugging histograms
      fhControlDInJ=new TH1I("hControlDInJ","Checks D in Jet",8, -0.5,7.5);
      fhControlDInJ->GetXaxis()->SetBinLabel(1,"DR In,1 daugh out");
      fhControlDInJ->GetXaxis()->SetBinLabel(2,"DR In,2 daugh out");
      fhControlDInJ->GetXaxis()->SetBinLabel(3,"DR In,3 daugh out");
      fhControlDInJ->GetXaxis()->SetBinLabel(4,"Tot tracks non matched");
      fhControlDInJ->GetXaxis()->SetBinLabel(5,"D0 daug missing");
      fhControlDInJ->GetXaxis()->SetBinLabel(6,"soft pi missing");
      fhControlDInJ->GetXaxis()->SetBinLabel(7,"IDprong diff GetDau");
      fhControlDInJ->GetXaxis()->SetBinLabel(8,"still z>1 (cand)");
      
      fhControlDInJ->SetNdivisions(1);
      fhControlDInJ->GetXaxis()->SetLabelSize(0.05);
      fOutput->Add(fhControlDInJ);
      
      fhmissingp=new TH1F("hmissingp", "Distribution of missing momentum (modulo);|p_{missing}|", 200, 0.,20.);
      fOutput->Add(fhmissingp);
      
      fhMissPi=new TH1F*[3];
      for(Int_t k=0;k<3;k++) {
      	 fhMissPi[k]=new TH1F(Form("hMissP%d",k), Form("Missing p component {%d};p_{%d}",k,k), 400, -10.,10.);
      	 fOutput->Add(fhMissPi[k]);
      }
      fhDeltaPtJet=new TH1F("hDeltaPtJet", "Difference between the jet pt before and after correction;p_{T}^{jet,recalc}-p_{T}^{jet,orig}", 200, 0.,20.);
      
      fOutput->Add(fhDeltaPtJet);
      fhRelDeltaPtJet=new TH1F("hRelDeltaPtJet", "Difference between the jet pt before and after correction/ original pt;(p_{T}^{jet,recalc}-p_{T}^{jet,orig})/p_{T}^{jet,orig}", 200, 0.,1.);
      fOutput->Add(fhRelDeltaPtJet);
      
      fhzDinjet=new TH1F("hzDinjet","Z of candidates with daughters in jet;z",nbinsz,zlims[0],zlims[1]);
      fOutput->Add(fhzDinjet);
      
      
      //calculate frag func with pt (simply ptD(or track)\cdot pt jet /ptjet^2)
      fhzDT=new TH1F("hzDT", Form("Z of D %s in jet in transverse components;(p_{T}^{D} dot p_{T}^{jet})/p_{T}^{jet}^{2} ", fUseMCInfo ? "(S+B)" : ""),nbinsz,zlims[0],zlims[1]);
      fOutput->Add(fhzDT);
      
      fhIDddaugh=new TH1I("hIDddaugh", "ID of daughters;ID", 2001,-1000,1000);
      fOutput->Add(fhIDddaugh);
      fhIDddaughOut=new TH1I("hIDddaughOut", "ID of daughters not found in jet;ID", 2001,-1000,1000);
      fOutput->Add(fhIDddaughOut);
      fhIDjetTracks=new TH1I("hIDjetTracks", "ID of jet tracks;ID", 2001,-1000,1000);
      fOutput->Add(fhIDjetTracks);
      
      fhDRdaughOut=new TH1F("hDRdaughOut", "#Delta R of daughters not belonging to the jet tracks (D in jet);#Delta R",200, 0.,2.);
      fOutput->Add(fhDRdaughOut);
      
      
      if(fCandidateType==kDstartoKpipi) 
      {
      	 if(fSwitchOnSB){
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
       
      
      fhDeltaRD=new TH1F("hDeltaRD",Form("#Delta R distribution of D candidates %s selected;#Delta R", fUseMCInfo ? "(S+B)" : ""),200, 0.,10.);
      fhDeltaRD->Sumw2();
      fOutput->Add(fhDeltaRD);
      
      fhDeltaRptDptj=new TH3F("hDeltaRptDptj",Form("#Delta R distribution of D candidates %s selected;#Delta R;#it{p}_{T}^{D} (GeV/c);#it{p}_{T}^{jet} (GeV/c)", fUseMCInfo ? "(S)" : ""),100, 0.,5.,nbinsptD, ptDlims[0],ptDlims[1],nbinsptjet,ptjetlims[0],ptjetlims[1]);
      fhDeltaRptDptj->Sumw2();
      fOutput->Add(fhDeltaRptDptj);
      
      if(fUseMCInfo){
      	 fhDeltaRptDptjB=new TH3F("hDeltaRptDptjB",Form("#Delta R distribution of D candidates (B) selected;#Delta R;#it{p}_{T}^{D} (GeV/c);#it{p}_{T}^{jet} (GeV/c)"),100, 0.,5.,nbinsptD, ptDlims[0],ptDlims[1],nbinsptjet,ptjetlims[0],ptjetlims[1]);
      	 fhDeltaRptDptjB->Sumw2();
      	 fOutput->Add(fhDeltaRptDptjB);
      }
      
      //background (side bands for the Dstar and like sign for D0)
      AliJetContainer *jetCont=GetJetContainer(0);
      if(!jetCont){
      	 Printf("Container 0 not found, try with name %s", fJetArrName.Data());
      	 jetCont=GetJetContainer(fJetArrName);
      	 if(!jetCont) Printf("NOT FOUND AGAIN");
      	 return kFALSE;
      }
      Printf("CONTAINER NAME IS %s", jetCont->GetArrayName().Data());
      //fJetRadius=jetCont->GetJetRadius();
      fhInvMassptD = new TH2F("hInvMassptD",Form("D (Delta R < %.1f) invariant mass distribution p_{T}^{j} > threshold",fJetRadius),nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      fhInvMassptD->SetStats(kTRUE);
      fhInvMassptD->GetXaxis()->SetTitle("mass (GeV)");
      fhInvMassptD->GetYaxis()->SetTitle("#it{p}_{t}^{D} (GeV/c)");
      fhInvMassptD->Sumw2();
      
      fOutput->Add(fhInvMassptD);
      
      if(fUseMCInfo){
      	 fhInvMassptDbg = new TH2F("hInvMassptDbg",Form("Background D (Delta R < %.1f) invariant mass distribution p_{T}^{j} > threshold",fJetRadius),nbinsmass,fMinMass,fMaxMass,nbinsptD,ptDlims[0],ptDlims[1]);
      	 fhInvMassptDbg->GetXaxis()->SetTitle(Form("%s (GeV)",(fCandidateType==kDstartoKpipi) ? "M(Kpipi)" : "M(Kpi)"));
      	 fhInvMassptDbg->GetYaxis()->SetTitle("p_{t}^{D} (GeV/c)");
      	 fhInvMassptDbg->Sumw2();
      	 fOutput->Add(fhInvMassptDbg);
      	 
      }
      
      if(fSwitchOnSparses){
      	 Double_t pi=TMath::Pi(), philims2[2];
      	 philims2[0]=-pi*1./2.;
      	 philims2[1]=pi*3./2.;
      	 fhsDphiz=0x0; //definition below according to the switches
      	 
      	 if(fSwitchOnSB && fSwitchOnPhiAxis && fSwitchOnOutOfConeAxis){
      	    AliInfo("Creating a 9 axes container: This might cause large memory usage");
      	    const Int_t nAxis=9;   
      	    const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsphi,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,2, 2, 2, 2};
      	    const Double_t minSparse[nAxis]={zlims[0],philims2[0],ptjetlims[0],ptDlims[0],fMinMass,-0.5, -0.5,-0.5,-0.5};
      	    const Double_t maxSparse[nAxis]={zlims[1],philims2[1],ptjetlims[1],ptDlims[1],fMaxMass, 1.5, 1.5, 1.5 , 1.5};
      	    fNAxesBigSparse=nAxis;
      	    
      	    fhsDphiz=new THnSparseF("hsDphiz","Z and #Delta#phi vs p_{T}^{jet}, p_{T}^{D}, mass. SB? D within the jet cone?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);
      	 }
      	 
      	 if(!fSwitchOnPhiAxis || !fSwitchOnOutOfConeAxis || !fSwitchOnSB){
      	    fSwitchOnPhiAxis=0;
      	    fSwitchOnOutOfConeAxis=0;
      	    fSwitchOnSB=0;
      	    if(fUseMCInfo){
      	       AliInfo("Creating a 7 axes container (MB background candidates)");
      	       const Int_t nAxis=7;   
      	       const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,2, 2, 2};
      	       const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass, -0.5,-0.5,-0.5};
      	       const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass, 1.5, 1.5 , 1.5};
      	       fNAxesBigSparse=nAxis;
      	       fhsDphiz=new THnSparseF("hsDphiz","Z vs p_{T}^{jet}, p_{T}^{D}, mass. Bkg?, D in EMCal acc?, jet in EMCal acc?", nAxis, nbinsSparse, minSparse, maxSparse);
      	       
      	    } else{
      	       AliInfo("Creating a 4 axes container");
      	       const Int_t nAxis=5;
      	       const Int_t nbinsSparse[nAxis]={nbinsSpsz,nbinsSpsptjet,nbinsSpsptD,nbinsSpsmass,63};
      	       const Double_t minSparse[nAxis]={zlims[0],ptjetlims[0],ptDlims[0],fMinMass,0};
      	       const Double_t maxSparse[nAxis]={zlims[1],ptjetlims[1],ptDlims[1],fMaxMass,6.3};
      	       fNAxesBigSparse=nAxis;      	 
      	       
      	       fhsDphiz=new THnSparseF("hsDphiz","Z, p_{T}^{jet}, p_{T}^{D}, mass., Jet phi", nAxis, nbinsSparse, minSparse, maxSparse);
      	    }
      	 }
      	 if(!fhsDphiz) AliFatal("No THnSparse created");
      	 fhsDphiz->Sumw2();
      	 
      	 fOutput->Add(fhsDphiz);
      }

   PostData(1, fOutput);
   
   return kTRUE; 
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsRecoJetCorr(AliVParticle* candidate, AliEmcalJet *jet,  AliAODEvent* aodEvent){
   
   Double_t ptD=candidate->Pt();
   Double_t deltaR=DeltaR(jet,candidate);
   Double_t phiD=candidate->Phi();
   Double_t deltaphi = jet->Phi()-phiD;
   if(deltaphi<=-(TMath::Pi())/2.) deltaphi = deltaphi+2.*(TMath::Pi());
   if(deltaphi>(3.*(TMath::Pi()))/2.) deltaphi = deltaphi-2.*(TMath::Pi());
   Double_t z=Z(candidate,jet);
   /*
   if(z>1) {
      ((TH1I*)fOutput->FindObject("hControlDInJ"))->Fill(7);
      Double_t pmissing=TMath::Sqrt(fPmissing[0]*fPmissing[0]+fPmissing[1]*fPmissing[1]+fPmissing[2]*fPmissing[2]);
      
      for(Int_t k=0;k<3;k++) ((TH1F*)fOutput->FindObject(Form("hMissP%d",k)))->Fill(fPmissing[k]);
      
      ((TH1F*)fOutput->FindObject("hmissingp"))->Fill(pmissing);
      Double_t ptdiff=fPtJet-jet->Pt();
      ((TH1F*)fOutput->FindObject("hDeltaPtJet"))->Fill(ptdiff);
      ((TH1F*)fOutput->FindObject("hRelDeltaPtJet"))->Fill(ptdiff/jet->Pt());

      
   }
   */
   if(fIsDInJet)fhzDT->Fill(Z(candidate,jet,kTRUE));
      
   fhDeltaRD->Fill(deltaR);
   fhDeltaRptDptj->Fill(deltaR,ptD,fPtJet);
   
   Bool_t bDInEMCalAcc=InEMCalAcceptance(candidate);
   Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);
   if(fUseReco){
      if(fCandidateType==kD0toKpi) {
      	 AliAODRecoDecayHF* dzero=(AliAODRecoDecayHF*)candidate;
      	 
      	 FillHistogramsD0JetCorr(dzero,deltaphi,z,ptD,fPtJet,jet->Phi(),bDInEMCalAcc,bJetInEMCalAcc, aodEvent);
      	 
      }
      
      if(fCandidateType==kDstartoKpipi) {
      	 AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)candidate;
      	 FillHistogramsDstarJetCorr(dstar,deltaphi,z,ptD,fPtJet,jet->Phi(),bDInEMCalAcc,bJetInEMCalAcc);
      	 
      }
   } else{ //generated level
      //AliInfo("Non reco");
      FillHistogramsMCGenDJetCorr(deltaphi,z,ptD,fPtJet,jet->Phi(),bDInEMCalAcc,bJetInEMCalAcc);
   }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t dPhi, Double_t z, Double_t ptD, Double_t ptj, Double_t phij, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc, AliAODEvent* aodEvent){


   Double_t masses[2]={0.,0.};
   Int_t pdgdaughtersD0[2]={211,321};//pi,K 
   Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
   
   masses[0]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
   masses[1]=candidate->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
   
   //TH3F* hPtJetWithD=(TH3F*)fOutput->FindObject("hPtJetWithD");
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=dPhi;
      point[2]=ptj;
      point[3]=ptD;
      point[4]=masses[0];
      point[5]=0;
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==5){
      point=new Double_t[5];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=phij;
      
}
  if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=masses[0];
      point[4]=0;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);

   }
   if(!point){
      AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
      return;
   }
   
   //Printf("Candidate in FillHistogramsD0JetCorr IsA %s", (candidate->IsA())->GetName());   
   Int_t isselected=fCuts->IsSelected(candidate,AliRDHFCuts::kAll,aodEvent);
   if(isselected==1 || isselected==3) {
      
      //if(fIsDInJet) hPtJetWithD->Fill(ptj,masses[0],ptD);
      
      FillMassHistograms(masses[0], ptD);
      if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);
   }
   if(isselected>=2) {
      //if(fIsDInJet) hPtJetWithD->Fill(ptj,masses[1],ptD);
      
      FillMassHistograms(masses[1], ptD);
      if(fNAxesBigSparse==9) point[4]=masses[1];
      if(fNAxesBigSparse==5 || fNAxesBigSparse==7) point[3]=masses[1];
      if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);
   }
   delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t dPhi,  Double_t z, Double_t ptD, Double_t ptj, Double_t phij, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){
  //dPhi and z not used at the moment,but will be (re)added

   AliAODTrack *softpi = (AliAODTrack*)dstar->GetBachelor();
   Double_t deltamass= dstar->DeltaInvMass(); 
   //Double_t massD0= dstar->InvMassD0();
   
   
   fhPtPion->Fill(softpi->Pt());
   
   //TH3F* hPtJetWithD=(TH3F*)fOutput->FindObject("hPtJetWithD");
   Int_t isSB=0;//IsDzeroSideBand(dstar);
   
   //Double_t point[]={z,dPhi,ptj,ptD,deltamass,static_cast<Double_t>(isSB), static_cast<Double_t>(fIsDInJet ? 1 : 0),bDInEMCalAcc,bJetInEMCalAcc};
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=dPhi;
      point[2]=ptj;
      point[3]=ptD;
      point[4]=deltamass;
      point[5]=static_cast<Double_t>(isSB);
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==5){
      point=new Double_t[5];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=deltamass;
      point[4]=phij;
   }
   if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=ptj;
      point[2]=ptD;
      point[3]=deltamass;
      point[4]=0;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }

   if(!point){
      AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
      return;
   }

   //if(fIsDInJet) hPtJetWithD->Fill(ptj,deltamass,ptD);
   
   FillMassHistograms(deltamass, ptD);
   if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);
   delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillHistogramsMCGenDJetCorr(Double_t dPhi,Double_t z,Double_t ptD,Double_t ptjet, Double_t phij, Bool_t bDInEMCalAcc, Bool_t bJetInEMCalAcc){
   
   Double_t pdgmass=0;
   if(fCandidateType==kD0toKpi) pdgmass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
   if(fCandidateType==kDstartoKpipi) pdgmass=TDatabasePDG::Instance()->GetParticle(413)->Mass()-TDatabasePDG::Instance()->GetParticle(421)->Mass();
   //TH3F* hPtJetWithD=(TH3F*)fOutput->FindObject("hPtJetWithD");
   
   //Double_t point[9]={z,dPhi,ptjet,ptD,pdgmass,0, static_cast<Double_t>(fIsDInJet ? 1 : 0),bDInEMCalAcc,bJetInEMCalAcc};
   
   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=dPhi;
      point[2]=ptjet;
      point[3]=ptD;
      point[4]=pdgmass;
      point[5]=0;
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(fNAxesBigSparse==5){
      point=new Double_t[5];
      point[0]=z;
      point[1]=ptjet;
      point[2]=ptD;
      point[3]=pdgmass;
      point[4]=phij;
   }
      if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=ptjet;
      point[2]=ptD;
      point[3]=pdgmass;
      point[4]=1;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }

   if(!point){
      AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
      return;
   }

   
   if(fNAxesBigSparse==9) point[4]=pdgmass;
   if(fNAxesBigSparse==5 || fNAxesBigSparse==7) point[3]=pdgmass;
   if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);
   //if(fIsDInJet) {
   //  hPtJetWithD->Fill(ptjet,pdgmass,ptD); // candidates within a jet
   //}
   delete[] point;
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FillMassHistograms(Double_t mass,Double_t ptD){
   
   if(fIsDInJet) {
      fhInvMassptD->Fill(mass,ptD);
   }
}

//________________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::FlagFlavour(AliEmcalJet *jet){
   
   AliEmcalJet::EFlavourTag tag=AliEmcalJet::kDStar;
   if (fCandidateType==kD0toKpi) tag=AliEmcalJet::kD0;
   if (fIsDInJet) jet->AddFlavourTag(tag);
   
   return;
   
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::SideBandBackground(AliAODRecoCascadeHF *candDstar, AliEmcalJet *jet){
   
   //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
   // (expected detector resolution) on the left and right frm the D0 mass. Each band
   //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
   
   //TH3F* hPtJetWithDsb=(TH3F*)fOutput->FindObject("hPtJetWithDsb");
   
   Bool_t bDInEMCalAcc=InEMCalAcceptance(candDstar);  
   Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);
   
   Double_t deltaM=candDstar->DeltaInvMass(); 
   //Printf("Inv mass = %f between %f and %f or %f and %f?",invM, sixSigmal,fourSigmal,fourSigmar,sixSigmar);
   Double_t z=Z(candDstar,jet);
   Double_t ptD=candDstar->Pt();

   Double_t dPhi=jet->Phi()-candDstar->Phi();

   if(dPhi<=-(TMath::Pi())/2) dPhi = dPhi+2*(TMath::Pi());
   if(dPhi>(3*(TMath::Pi()))/2) dPhi = dPhi-2*(TMath::Pi());
   
   Int_t isSideBand=1;//IsDzeroSideBand(candDstar);
   //if no SB analysis we should not enter here, so no need of checking the number of axes
   Double_t point[9]={z,dPhi,fPtJet,ptD,deltaM,static_cast<Double_t>(isSideBand), static_cast<Double_t>(fIsDInJet ? 1 : 0),static_cast<Double_t>(bDInEMCalAcc? 1 : 0),static_cast<Double_t>(bJetInEMCalAcc? 1 : 0)};
   //fill the background histogram with the side bands or when looking at MC with the real background
   if(isSideBand==1){
      fhDiffSideBand->Fill(deltaM,ptD); // M(Kpipi)-M(Kpi) side band background    
      //hdeltaPhiDjaB->Fill(deltaM,ptD,dPhi);
      if(fSwitchOnSparses) fhsDphiz->Fill(point,1.);
      //if(fIsDInJet){  // evaluate in the near side	
      //	 hPtJetWithDsb->Fill(fPtJet,deltaM,ptD);
      //}
   }
}

//_______________________________________________________________________________

void AliAnalysisTaskFlavourJetCorrelations::MCBackground(AliAODRecoDecayHF *candbg,AliEmcalJet* jet){
   
   //need mass, deltaR, pt jet, ptD
   //two cases: D0 and Dstar
   
   Int_t isselected=fCuts->IsSelected(candbg,AliRDHFCuts::kAll, AODEvent());
   if(!isselected) return;
   
   Double_t ptD=candbg->Pt();
   Double_t deltaR=DeltaR(jet,candbg);
   Double_t phiD=candbg->Phi();
   Double_t deltaphi = jet->Phi()-phiD;
   if(deltaphi<=-(TMath::Pi())/2.) deltaphi = deltaphi+2.*(TMath::Pi());
   if(deltaphi>(3.*(TMath::Pi()))/2.) deltaphi = deltaphi-2.*(TMath::Pi());
   Double_t z=Z(candbg,jet);

   if(fIsDInJet) fhzDT->Fill(Z(candbg,jet,kTRUE));
   
   
   
   
   fhDeltaRD->Fill(deltaR);
   fhDeltaRptDptjB->Fill(deltaR,ptD,fPtJet);

   Bool_t bDInEMCalAcc=InEMCalAcceptance(candbg);
   Bool_t bJetInEMCalAcc=InEMCalAcceptance(jet);

   //TH3F* hPtJetWithDsb=(TH3F*)fOutput->FindObject("hPtJetWithDsb");

   Double_t *point=0x0;
   if(fNAxesBigSparse==9){
      point=new Double_t[9];
      point[0]=z;
      point[1]=deltaphi;
      point[2]=fPtJet;
      point[3]=ptD;
      point[4]=-999; //set below
      point[5]=1;
      point[6]=static_cast<Double_t>(fIsDInJet ? 1 : 0);
      point[7]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[8]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }

   if(fNAxesBigSparse==7){
      point=new Double_t[7];
      point[0]=z;
      point[1]=fPtJet;
      point[2]=ptD;
      point[3]=-999; //set below
      point[4]=1;
      point[5]=static_cast<Double_t>(bDInEMCalAcc ? 1 : 0);
      point[6]=static_cast<Double_t>(bJetInEMCalAcc ? 1 : 0);
   }
   if(!point){
      AliError(Form("Numer of THnSparse entries %d not valid", fNAxesBigSparse));
      return;
   }

   if(fCandidateType==kDstartoKpipi){
      AliAODRecoCascadeHF* dstarbg = (AliAODRecoCascadeHF*)candbg;
      Double_t deltaM=dstarbg->DeltaInvMass();
      fhInvMassptDbg->Fill(deltaM,ptD);
      //if(fIsDInJet) hPtJetWithDsb->Fill(fPtJet,deltaM,ptD);
      if(fNAxesBigSparse==9) point[4]=deltaM;
      if(fNAxesBigSparse==5 || fNAxesBigSparse==7) point[3]=deltaM;
      if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);      
   }
   
   if(fCandidateType==kD0toKpi){
      Double_t masses[2]={0.,0.};
      Int_t pdgdaughtersD0[2]={211,321};//pi,K 
      Int_t pdgdaughtersD0bar[2]={321,211};//K,pi 
      
      masses[0]=candbg->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0); //D0
      masses[1]=candbg->InvMass(fNProngs,(UInt_t*)pdgdaughtersD0bar); //D0bar
      
      
      if(isselected==1 || isselected==3) {
      	 //if(fIsDInJet) hPtJetWithDsb->Fill(fPtJet,masses[0],ptD);
      	 fhInvMassptDbg->Fill(masses[0],ptD);
      	 if(fNAxesBigSparse==9) point[4]=masses[0];
      	 if(fNAxesBigSparse==5 || fNAxesBigSparse==7) point[3]=masses[0];
      	 if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);
     }
      if(isselected>=2) {
      	 //if(fIsDInJet) hPtJetWithDsb->Fill(fPtJet,masses[1],ptD);
      	 fhInvMassptDbg->Fill(masses[1],ptD);
      	 if(fNAxesBigSparse==9) point[4]=masses[1];
      	 if(fNAxesBigSparse==5 || fNAxesBigSparse==7) point[3]=masses[1];
      	 if(fSwitchOnSparses && (fSwitchOnOutOfConeAxis || fIsDInJet)) fhsDphiz->Fill(point,1.);
      	 
      }
      
      
   }
   delete[] point;
}

//_______________________________________________________________________________

Float_t AliAnalysisTaskFlavourJetCorrelations::DeltaR(AliEmcalJet *p1, AliVParticle *p2) const {
   //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
   //It recalculates the eta-phi values if it was asked for background subtraction of the jets
   if(!p1 || !p2) return -1;
    
    Double_t phi1=p1->Phi(),eta1=p1->Eta();
    
    //It subtracts the backgroud of jets if it was asked for it.

    if (!fJetCont->GetRhoName().IsNull()) {
            Double_t pj[3];
            Bool_t okpj=p1->PxPyPz(pj);
            if(!okpj){
                printf("Problems getting momenta\n");
                return -999;
            }
            pj[0] = p1->Px() - p1->Area()*(fRhoValue*TMath::Cos(p1->AreaPhi()));
            pj[1] = p1->Py() - p1->Area()*(fRhoValue*TMath::Sin(p1->AreaPhi()));
            pj[2] = p1->Pz() - p1->Area()*(fRhoValue*TMath::SinH(p1->AreaEta()));
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
Float_t AliAnalysisTaskFlavourJetCorrelations::CheckDeltaR(AliEmcalJet *p1, AliVParticle *p2) const {
    //Calculate DeltaR between p1 and p2: DeltaR=sqrt(Delataphi^2+DeltaEta^2)
    //It recalculates the eta-phi values if it was asked for background subtraction of the jets
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

Int_t AliAnalysisTaskFlavourJetCorrelations::IsDzeroSideBand(AliAODRecoCascadeHF *candDstar){
   
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

Bool_t AliAnalysisTaskFlavourJetCorrelations::AreDaughtersInJet(AliEmcalJet *thejet, AliAODRecoDecayHF* charm, Int_t*& daughOutOfJetID, AliAODTrack**& dtrks, Bool_t fillH){
   //returns true/false and the array of particles not found among jet constituents
   
   
   Int_t daughtersID[3];
   Int_t ndaugh=0;
   Int_t dmatchedID[3]={0,0,0};
   Int_t countmatches=0;
   if (fCandidateType==kDstartoKpipi){
      AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)charm;
      AliAODRecoDecayHF2Prong* dzero=dstar->Get2Prong();
      daughtersID[0]=dzero->GetProngID(0);
      daughtersID[1]=dzero->GetProngID(1);
      daughtersID[2]=dstar->GetBachelor()->GetID();
      ndaugh=3;
     
      dtrks=new AliAODTrack*[3];
      dtrks[0]=(AliAODTrack*)dzero->GetDaughter(0);
      dtrks[1]=(AliAODTrack*)dzero->GetDaughter(1);
      dtrks[2]=(AliAODTrack*)dstar->GetBachelor();
  
      //check
      if(fillH){
      	 if(daughtersID[0]!=((AliAODTrack*)dtrks[0])->GetID() || daughtersID[1]!=((AliAODTrack*)dtrks[1])->GetID())  fhControlDInJ->Fill(6);
      	 
      	 fhIDddaugh->Fill(daughtersID[0]);
      	 fhIDddaugh->Fill(daughtersID[1]);
      	 fhIDddaugh->Fill(daughtersID[2]);
      	 
      }
      //Printf("ID daughters %d, %d, %d",daughtersID[0], daughtersID[1], daughtersID[2]);
   }
   
   if (fCandidateType==kD0toKpi){
      daughtersID[0]=charm->GetProngID(0);
      daughtersID[1]=charm->GetProngID(1);
      ndaugh=2;
      if(fillH){
      	 fhIDddaugh->Fill(daughtersID[0]);
      	 fhIDddaugh->Fill(daughtersID[1]);
      }
      dtrks=new AliAODTrack*[2];
      dtrks[0]=(AliAODTrack*)charm->GetDaughter(0);
      dtrks[1]=(AliAODTrack*)charm->GetDaughter(1);

   }
   
   const Int_t cndaugh=ndaugh;
   daughOutOfJetID=new Int_t[cndaugh];

   Int_t nchtrk=thejet->GetNumberOfTracks();
   for(Int_t itk=0;itk<nchtrk;itk++){
      AliVTrack* tkjet=((AliPicoTrack*)thejet->TrackAt(itk,fTrackArr))->GetTrack();
      if(!tkjet) continue;
      Int_t idtkjet=tkjet->GetID();
      if(fillH) fhIDjetTracks->Fill(idtkjet);
      for(Int_t id=0;id<ndaugh;id++){
      	 if(idtkjet==daughtersID[id]) {
      	    dmatchedID[id]=daughtersID[id]; //daughter which matches a track in the jet
      	    countmatches++;
      	    
      	 }
      	 
      	 if(countmatches==ndaugh) break;
      }
   }
   //reverse: include in the array the ID of daughters not matching

   for(Int_t id=0;id<ndaugh;id++){
      if(dmatchedID[id]==0) {
      	 daughOutOfJetID[id]=daughtersID[id];
      	 if(fillH) fhIDddaughOut->Fill(daughOutOfJetID[id]);
      }
      else daughOutOfJetID[id]=0;
   }
   if(fillH){
      if((ndaugh-countmatches) == 1) fhControlDInJ->Fill(0);
      if((ndaugh-countmatches) == 2) fhControlDInJ->Fill(1);
      if((ndaugh-countmatches) == 3) fhControlDInJ->Fill(2);
   }
   if(ndaugh!=countmatches){
      return kFALSE;
   }
   
   return kTRUE;
}

//_______________________________________________________________________________

Bool_t AliAnalysisTaskFlavourJetCorrelations::IsDInJet(AliEmcalJet *thejet, AliAODRecoDecayHF* charm, Bool_t fillH){
   
   //check the conditions type:
   //type 0 : DeltaR < jet radius, don't check daughters (and don't correct pt jet) 
   //type 1 : DeltaR < jet radius and check for all daughters among jet tracks, don't correct ptjet
   //type 2 (default) : DeltaR < jet radius and check for all daughters among jet tracks, if not present, correct the ptjet
   //type 3 (under development) :  DeltaR < jet radius and check for all daughters among jet tracks, if not present, correct the ptjet usign the pt-scheme 
  
   fPmissing[0]=0; 
   fPmissing[1]=0;
   fPmissing[2]=0;
   
   Bool_t testDeltaR=kFALSE;
    
   if(!fJetCont->GetRhoName().IsNull())
   {
        if(DeltaR(thejet,charm)<fJetRadius && CheckDeltaR(thejet,charm)<fJetRadius) testDeltaR=kTRUE;
   }
   else
   {
       if(DeltaR(thejet,charm)<fJetRadius) testDeltaR=kTRUE;
   }

   //for type 3 this check should be redone after the modification of the jet axis
   
   Int_t* daughOutOfJet;
   AliAODTrack** charmDaugh;
   Bool_t testDaugh=AreDaughtersInJet(thejet, charm, daughOutOfJet,charmDaugh,fillH);
   if(testDaugh && testDeltaR) {
      //Z of candidates with daughters in the jet
      fhzDinjet->Fill(Z(charm,thejet));
      
   }
   
   TVector3 thejetv;    //(x=0,y=0,z=0)
   thejetv.SetPtEtaPhi(thejet->Pt(),thejet->Eta(),thejet->Phi());
   TVector3 newjet(thejetv);
   TVector3 correction; //(x=0,y=0,z=0)
   TVector3 charmcand;  //(x=0,y=0,z=0)
   charmcand.SetPtEtaPhi(charm->Pt(),charm->Eta(),charm->Phi());
   
   if(!testDaugh && testDeltaR && fTypeDInJet>=2){
      
      Int_t ndaugh=3;
      if(fCandidateType==kD0toKpi) ndaugh=2;
      Int_t nOut=ndaugh;
      
      for(Int_t id=0;id<ndaugh;id++){
      	 if(daughOutOfJet[id]!=0){ //non-matched daughter
      	    //get track and its momentum
      	    nOut--;
      	    if(fillH) {
      	       fhControlDInJ->Fill(3);
      	       if(id<2) fhControlDInJ->Fill(4);
      	       if(id==2)fhControlDInJ->Fill(5);
      	       fhDRdaughOut->Fill(DeltaR(thejet, charmDaugh[id]));
      	    }
      	    if(fTypeDInJet==2){
      	       
      	       newjet.SetX(newjet(0)+charmDaugh[id]->Px()); 
      	       newjet.SetY(newjet(1)+charmDaugh[id]->Py());
      	       newjet.SetZ(newjet(2)+charmDaugh[id]->Pz());
      	    }
      	    if(fTypeDInJet==3){
      	       
      	       Double_t ptdaug  = charmDaugh[id]->Pt();
      	       Double_t ptjet   = newjet.Perp();
      	       Double_t ptn     = ptjet+ptdaug;
      	       Double_t phidaug = charmDaugh[id]->Phi();
      	       Double_t phijet  = newjet.Phi();
      	       Double_t phin    = (phijet*ptjet+phidaug*ptdaug)/(ptjet+ptdaug);
      	       
      	       Double_t etadaug = charmDaugh[id]->Eta();
      	       Double_t etajet  = newjet.Eta();
      	       Double_t etan    = (etajet*ptjet+etadaug*ptdaug)/(ptjet+ptdaug);
      	       
      	       newjet.SetPtEtaPhi(ptn,etan,phin);
 
      	    }
      	 }      
      }
      
      correction=newjet-thejetv;
      fPmissing[0]=correction(0);
      fPmissing[1]=correction(1);
      fPmissing[2]=correction(2);

      //now the D is within the jet
      testDaugh=kTRUE;
      if(fTypeDInJet==3){ //recalc R(D,jet)
      	 Double_t dRnew=newjet.DeltaR(charmcand);
      	 if(dRnew<fJetRadius) testDeltaR=kTRUE;
      	 else testDeltaR=kFALSE;
      }
   }
   delete[] daughOutOfJet;
   delete[] charmDaugh;
   
   Bool_t result=0;
   switch(fTypeDInJet){
   case 0:
      result=testDeltaR;
      break;
   case 1:
      result=testDeltaR && testDaugh;
      break;
   case 2: //this case defines fPmissing
      result=testDeltaR && testDaugh; 
      break;
   case 3: //this case defines fPmissing and recalculate R(jet,D) with the updated jet axis
      result=testDeltaR && testDaugh; 
      break;
      
   default:
      AliInfo("Selection type not specified, use 1");
      result=testDeltaR && testDaugh;
      break;
   }
   return result;
}

Bool_t AliAnalysisTaskFlavourJetCorrelations::InEMCalAcceptance(AliVParticle *vpart){
   //check eta phi of a VParticle: return true if it is in the EMCal acceptance, false otherwise
   
   Double_t phiEMCal[2]={1.405,3.135},etaEMCal[2]={-0.7,0.7};
   Bool_t binEMCal=kTRUE;
   Double_t phi=vpart->Phi(), eta=vpart->Eta();
   if(phi<phiEMCal[0] || phi>phiEMCal[1]) binEMCal=kFALSE;
   if(eta<etaEMCal[0] || eta>etaEMCal[1]) binEMCal=kFALSE;
   return binEMCal;


}
