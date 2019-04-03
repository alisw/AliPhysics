#ifndef ALIANALYSISTASKSE_H

#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TArrayF.h>
#include <TArrayD.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliEmcalList.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliParticleContainer.h"
#include "AliInputEventHandler.h"
#include "math.h"
#endif

#include <time.h>
#include <TRandom3.h>
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliPicoTrack.h>
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAnalysisUtils.h"
#include "AliRhoParameter.h"
#include "TVector3.h"
#include "AliVVertex.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskJetPP.h"
#include "AliHeader.h" 
#include "AliRunLoader.h"  
#include "AliVVZERO.h" 
#include "AliVZDC.h"
#include "AliExternalTrackParam.h"
using namespace std;

// Analyis of jets in pp collisions 
// Author Peter Pribeli   (12.Sep. 2017)

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetPP)
/// \endcond

///Default constructor
AliAnalysisTaskJetPP::AliAnalysisTaskJetPP(): 
AliAnalysisTaskEmcalJet("AliAnalysisTaskJetPP", kTRUE),  
  fUseDefaultVertexCut(1), fUsePileUpCut(1), fIsMC(0), 
 fSignalJetRadius(0.4), 
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fASPDCvsTCut(65.), fBSPDCvsTCut(4.), fMinJetArea(0.0),  
fCentralityType("V0A"), fHelperClass(0), fInitializedLocal(0),
fZVertexCut(10.0), fhJetPt(0x0), fhCuts(0x0), fhTrackPt(0x0),fhJetConstituentPt(0x0),fhJetEtaPt(0x0),fhAktJetEtaPhi(0x0),fhKtJetEtaPhi(0x0),fhJetAreaPt(0x0),fhJetPhiPt(0x0),fhZVertex(0x0),fhYVertex(0x0),fhXVertex(0x0),fhJetPtRho(0x0),fhJetPtConeRho(0x0),fhJetPtCMSRho(0x0),fhRho(0x0),fhConeRho(0x0),fhCMSRho(0x0),fhTrackEtaPt(0x0),fhGenTrackEtaPt(0x0),fhTrackPhiPt(0x0),fhTrackEtaPhi(0x0),fhKTJetPt(0x0),fhZVertexBC(0x0), fhRemx(0x0), fhPrimGenTrkPt(0x0), fhGenJetPt(0x0), fhRecTrkPt(0x0), fhFakeTrkPt(0x0),fhMult(0x0),fhTrackPhiCG(0x0),fhTrackPhiTPCG(0x0),fhAtimesRhoMedian(0x0),fhAtimesRhoCone(0x0),fhAtimesRhoCMS(0x0)
{
  //Arrays initiation
   for(int i=0;i<2;i++){
      fhInvPtQVsPhi[i] = NULL;   
      fhInvPtQVsEta[i] = NULL; 
      fhInvPtQVsPhiASide[i] = NULL; 
      fhInvPtQVsPhiCSide[i] = NULL; 
      fhSigmaPtOverPtVsPt[i] = NULL; 
   }
}

/// Standard constructor
AliAnalysisTaskJetPP::AliAnalysisTaskJetPP(const char *name) : 
AliAnalysisTaskEmcalJet(name,kTRUE),  
  fUseDefaultVertexCut(1), fUsePileUpCut(1), fIsMC(0), 
 fSignalJetRadius(0.4), 
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fASPDCvsTCut(65.), fBSPDCvsTCut(4.), fMinJetArea(0.0),  
fCentralityType("V0A"), fHelperClass(0), fInitializedLocal(0),
fZVertexCut(10.0),fhJetPt(0x0),fhCuts(0x0), fhTrackPt(0x0),fhJetConstituentPt(0x0),fhJetEtaPt(0x0),fhAktJetEtaPhi(0x0),fhKtJetEtaPhi(0x0),fhJetAreaPt(0x0),fhJetPhiPt(0x0),fhZVertex(0x0),fhYVertex(0x0),fhXVertex(0x0),fhJetPtRho(0x0),fhJetPtConeRho(0x0),fhJetPtCMSRho(0x0),fhRho(0x0),fhConeRho(0x0),fhCMSRho(0x0),fhTrackEtaPt(0x0),fhTrackPhiPt(0x0),fhTrackEtaPhi(0x0),fhGenTrackEtaPt(0x0),fhKTJetPt(0x0),fhZVertexBC(0x0), fhRemx(0x0), fhPrimGenTrkPt(0x0), fhGenJetPt(0x0), fhRecTrkPt(0x0), fhFakeTrkPt(0x0),fhMult(0x0),fhTrackPhiCG(0x0),fhTrackPhiTPCG(0x0),fhAtimesRhoMedian(0x0),fhAtimesRhoCone(0x0),fhAtimesRhoCMS(0x0)

{
  //Arrays initiation
   for(int i=0;i<2;i++){
      fhInvPtQVsPhi[i] = NULL;   
      fhInvPtQVsEta[i] = NULL; 
      fhInvPtQVsPhiASide[i] = NULL; 
      fhInvPtQVsPhiCSide[i] = NULL; 
      fhSigmaPtOverPtVsPt[i] = NULL; 
   }


   DefineOutput(1, AliEmcalList::Class());
}
  
//________________________________________________________________________
Bool_t AliAnalysisTaskJetPP::IsEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION RECONSTRUCTED DATA

   if(!event) return kFALSE;
 
   //___________________________________________________
   //BEFORE VERTEX CUT
   fhZVertexBC->Fill(event->GetPrimaryVertex()->GetZ()); 
   if(fUseDefaultVertexCut){
      if(!fHelperClass || !fHelperClass->IsVertexSelected2013pA(event)){
         return kFALSE;
      }
      if(fHelperClass && fHelperClass->IsVertexSelected2013pA(event)) fhCuts->Fill(1.5);//events that passed the vertex cut
   }else{

      if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > fZVertexCut){ 
                 return kFALSE;
      }
      else fhCuts->Fill(1.5);//events that passed the vertex cut
   }
   fhZVertex->Fill(event->GetPrimaryVertex()->GetZ());
   fhXVertex->Fill(event->GetPrimaryVertex()->GetX());
   fhYVertex->Fill(event->GetPrimaryVertex()->GetY());

   //___________________________________________________
   //AFTER VERTEX CUT
   return kTRUE;
}
//PILEUP CUT
Bool_t AliAnalysisTaskJetPP::IsSPDClusterVsTrackletBG(AliVEvent *event){ //NEW PILEUP
   if(fUsePileUpCut){
      Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
      Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
      Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
      if (nClustersLayer0 + nClustersLayer1 > fASPDCvsTCut + nTracklets*fBSPDCvsTCut){
         return kTRUE;
      }
   }
   fhCuts->Fill(2.5);
   return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPP::IsTrackInAcceptance(AliVParticle* track, Bool_t isprimary){
   // CHECK THE TRACK PT AND ETA RANGE 
   if(!track) return kFALSE;

   if(isprimary) {
      if(!track->Charge()) return kFALSE;
      if(!(static_cast<AliAODMCParticle*>(track))->IsPhysicalPrimary()) return kFALSE;
  }
   
   if(TMath::Abs(track->Eta()) <= fTrackEtaWindow){ //APPLY TRACK ETA CUT
      if(track->Pt() >= fMinTrackPt){   //APPLY TRACK MINIMUM PT CUT
         return kTRUE;
      }
   }
   return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPP::IsSignalJetInAcceptance(AliEmcalJet *jet, Bool_t suppressGhost){   
   //select jets in acceptance 
   if(!jet) return kFALSE;
   if(TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow){
      if(jet->Area() >= fMinJetArea){
         if(suppressGhost){
            if(jet->Pt() >= fMinTrackPt) return kTRUE;
         }else{
            return kTRUE;
         }
      }
   }  
   return kFALSE;
}

/// Perform steps needed to initialize the analysis.
void AliAnalysisTaskJetPP::ExecOnceLocal(){
   // Initialization of jet containers done in  AliAnalysisTaskEmcalJet::ExecOnce()
   //Read arrays of jets and tracks
   fInitializedLocal = kTRUE; 

   // Initialize helper class (for vertex selection & pile up correction)
   fHelperClass = new AliAnalysisUtils();
   fHelperClass->SetCutOnZVertexSPD(kFALSE); // kFALSE: no cut; kTRUE: |zvtx-SPD - zvtx-TPC|<0.5cm

   return;
}

/// Overloads base class method. Fills the output histograms
/// return kTRUE if successful
Bool_t AliAnalysisTaskJetPP::FillHistograms(){
   if(!InputEvent()){
      AliError("??? Event pointer == 0 ???");
      return kFALSE;
   }
   fhCuts->Fill(0.5);
   //Execute only once:  Get tracks, jets from arrays if not already given 
   if(!fInitializedLocal) ExecOnceLocal(); 

   //_________________________________________________________________
   //  FILL EVENT STATISTICS

   //Select events (vertex, pile-up,...) 
   if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec
   if(IsSPDClusterVsTrackletBG(InputEvent())) return kFALSE;


   // JET+TRACK CONTAINERS
   AliJetContainer *jetContRecAKT = NULL; //AKTjet container from reconstruced tracks
   AliJetContainer *jetContRecKT = NULL; //KTjet container from reconstruced tracks
   AliJetContainer *jetContGenAKT  = NULL; //AKT jet container from generated MC tracks
   AliJetContainer *jetContGenKT  = NULL; //KT jet container from generated MC tracks

   AliEmcalJet  *jetRec = NULL;
   AliEmcalJet  *jetGen = NULL;//MC

   AliParticleContainer *trkContRec = NULL; //track array of real reconstructed tracks 
   AliParticleContainer *parContGen = NULL; //track array of generated MC tracks 

   AliVParticle *constTrackRec = NULL; //rec jet constituent
   AliVParticle *constTrackGen = NULL; //gen jet constituent
   //_________________________________________________________
   //READ JET TRACK CONTAINERS
   jetContRecAKT  = GetJetContainer(kContainerOne); //AKT jet container
   jetContRecKT   = GetJetContainer(kContainerTwo); //KT  jet container
   if(fIsMC){
      jetContGenAKT  = GetJetContainer(kContainerThree); //MC AKT jet container
      jetContGenKT   = GetJetContainer(kContainerFour); //MC KT  jet container
   } 

   trkContRec     = GetParticleContainer(kContainerOne); //reconstructed hybrid tracks 
   if(fIsMC) parContGen     = GetParticleContainer(kContainerTwo); //MC particles 

   Int_t fMultCounter;
   Double_t xyz[50];
   Double_t pxpypz[50];
   Double_t cv[21];
   Int_t itrkq;

   AliAODTrack *constTrackRecAod=NULL;

// Lop over reconstructed trcks
   if(trkContRec){
      fMultCounter=0;
      for (auto trackIterator : trkContRec->accepted_momentum()){
         constTrackRecAod = (AliAODTrack*)trackIterator.second;
         if(!constTrackRecAod->IsHybridGlobalConstrainedGlobal()) continue;
         if(!constTrackRecAod) continue;
         if(!IsTrackInAcceptance(constTrackRecAod)) continue;
            Double_t phi = Convert(constTrackRecAod->Phi());        
	         fhTrackPt->Fill(constTrackRecAod->Pt());
            fhTrackEtaPt->Fill(constTrackRecAod->Eta(),constTrackRecAod->Pt()); 
            fhTrackPhiPt->Fill(phi,constTrackRecAod->Pt()); 
            fhTrackEtaPhi->Fill(constTrackRecAod->Eta(),phi);            

            if(constTrackRecAod->IsGlobalConstrained()){
               fhTrackPhiCG->Fill(constTrackRecAod->Pt(), phi); //global constrained
            }else{
               fhTrackPhiTPCG->Fill(constTrackRecAod->Pt(), phi); //complementary
            }
            itrkq = (constTrackRecAod->Charge()<0) ? 0 : 1;

            fhInvPtQVsPhi[itrkq]->Fill(phi, 1.0/constTrackRecAod->Pt());
            fhInvPtQVsEta[itrkq]->Fill(constTrackRecAod->Eta(), 1.0/constTrackRecAod->Pt());

            if(constTrackRecAod->Eta()>0){
               fhInvPtQVsPhiASide[itrkq]->Fill(phi, 1.0/constTrackRecAod->Pt());
            }else{
               fhInvPtQVsPhiCSide[itrkq]->Fill(phi, 1.0/constTrackRecAod->Pt());
            }
            memset(cv, 0, sizeof(Double_t) * 21); //cleanup arrays
            memset(pxpypz, 0, sizeof(Double_t) * 50);
            memset(xyz, 0, sizeof(Double_t) * 50);
            constTrackRecAod->GetXYZ(xyz);
            constTrackRecAod->GetPxPyPz(pxpypz);
            constTrackRecAod->GetCovarianceXYZPxPyPz(cv);
            AliExternalTrackParam  par(xyz, pxpypz, cv, constTrackRecAod->Charge());
            fhSigmaPtOverPtVsPt[itrkq]->Fill(constTrackRecAod->Pt(), TMath::Abs(sqrt(par.GetSigma1Pt2())/par.GetSigned1Pt()));
            fMultCounter++;
      }
      fhMult->Fill(fMultCounter);
   }
   //_________________________________________________   
   // Background
   Double_t rho    = EstimateBgKT(jetContRecKT);
   fhRho->Fill(rho); 

   Double_t conerho = EstimateLocalBg(jetContRecAKT,trkContRec);
   fhConeRho->Fill(conerho); 

   Double_t cmsrho = EstimateBgKTCMS(jetContRecKT);
   fhCMSRho->Fill(cmsrho);
   //_________________________________________________
   // Loop over AKT jets
   if(jetContRecAKT){ 
      jetContRecAKT->ResetCurrentID();
      for(auto jetIterator : jetContRecAKT->accepted_momentum()) {//loop over reconstructed jets
         jetRec=(AliEmcalJet*)jetIterator.second;
         if(!jetRec) continue;
         if(!IsSignalJetInAcceptance(jetRec)) continue; //check jet eta and pt
         Double_t phi = Convert(jetRec->Phi());
         fhJetPt->Fill(jetRec->Pt()); //fill AKT jet pT
         fhJetPtRho->Fill(jetRec->Pt()-rho*jetRec->Area()); //fill AKT jet pT without bkg
         fhJetPtConeRho->Fill(jetRec->Pt()-conerho*jetRec->Area()); //fill AKT jet pT without local bkg
         fhJetPtCMSRho->Fill(jetRec->Pt()-cmsrho*jetRec->Area()); //fill AKT jet pT without CMS bkg
         fhJetEtaPt->Fill(jetRec->Eta(),jetRec->Pt());
         fhAktJetEtaPhi->Fill(jetRec->Eta(),phi);
         fhJetPhiPt->Fill(phi,jetRec->Pt());
         fhJetAreaPt->Fill(jetRec->Area(),jetRec->Pt());
         fhAtimesRhoMedian->Fill(rho*jetRec->Area());
         fhAtimesRhoCone->Fill(conerho*jetRec->Area());
         fhAtimesRhoCMS->Fill(cmsrho*jetRec->Area());

         for(Int_t iq=0; iq < jetRec->GetNumberOfTracks(); iq++){
            constTrackRec = static_cast<AliVParticle*> (jetRec->TrackAt(iq, trkContRec->GetArray())); //matched rec and emb tracks
            if(!constTrackRec) continue;
	       fhJetConstituentPt->Fill(constTrackRec->Pt());
 
         }
      }
   }
  // Loop over KT jets
  if(jetContRecKT){
   for(auto jetIterator : jetContRecKT->accepted_momentum()) {
      jetRec=(AliEmcalJet*)jetIterator.second;
	   if(!jetRec) continue;
	   if(!IsSignalJetInAcceptance(jetRec)) continue;

	   fhKTJetPt->Fill(jetRec->Pt());
        }
   }

   // Single particle efficiency and contamination
   //primary particles spectrum
   if(parContGen){
      for(auto trackIterator : parContGen->accepted_momentum()) {
         constTrackGen=(AliVParticle*)trackIterator.second;
         if(!constTrackGen) continue;
         if(IsTrackInAcceptance(constTrackGen, kTRUE)){
            fhPrimGenTrkPt->Fill(constTrackGen->Pt()); //pT spectrum of generator level particles eta-pt, phi-pt, eta-phi histogram
            fhGenTrackEtaPt->Fill(constTrackGen->Eta(),constTrackGen->Pt()); //pT spectrum of generator level particles eta-pt, phi-pt, eta-phi histogram
         }
      }
   }


   //reconstructed primary + secondary particles
   Bool_t bRecPrim = kFALSE; //tags the reconstructed primary particles
   if(trkContRec && parContGen){
      for(auto trackIterator : trkContRec->accepted_momentum()) {
         constTrackRecAod=(AliAODTrack*)trackIterator.second;
         if(!constTrackRecAod->IsHybridGlobalConstrainedGlobal()) continue;
         if(!constTrackRecAod) continue;
         if(!IsTrackInAcceptance(constTrackRecAod, kFALSE)) continue; //reconstructed level tracks
         bRecPrim = kFALSE; //reset matching flag

         parContGen->ResetCurrentID();
         for(auto trackIterator : parContGen->accepted_momentum()) {
            constTrackGen=(AliAODTrack*)trackIterator.second;
            if(!constTrackGen) continue;
            if(!IsTrackInAcceptance(constTrackGen, kTRUE)) continue; //gen level physical primary
            if(TMath::Abs(constTrackRecAod->GetLabel()) == TMath::Abs(constTrackGen->GetLabel())){
               fhRecTrkPt->Fill(constTrackGen->Pt());
               bRecPrim = kTRUE; //matched
               break;
            }
         }
         if(!bRecPrim){
            fhFakeTrkPt->Fill(constTrackRecAod->Pt());//Pt of fake tracks
         }
      }
   }
   // Response matrix
   Double_t ptGenCorr; //GEN jet pt corrected for rho
   Double_t ptRecCorr; //REC jet pt corrected for rho

   //Response matrix normalization - spectrum of all generator level jets in acceptance
   if(jetContGenAKT){
      for(auto jetIterator : jetContGenAKT->accepted_momentum()) {
         jetGen=(AliEmcalJet*)jetIterator.second;
         if(!jetGen) continue;
         if(!IsSignalJetInAcceptance(jetGen,kTRUE)) continue; //cuts on eta, pT ,are
         fhGenJetPt->Fill(jetGen->Pt()); //Pt spectrum of MC jets 
      }
   }
      //Find closest gen level+rec level  jets
   if(jetContRecAKT){ 
      for(auto jetIterator : jetContRecAKT->accepted_momentum()) { 
         jetRec=(AliEmcalJet*)jetIterator.second;
         if(!jetRec) continue;
         if(!IsSignalJetInAcceptance(jetRec,kTRUE)) continue; //cuts on eta, pT ,area

         jetGen = 0x0;
         jetGen = jetRec->ClosestJet(); 

         if(!jetGen){ //did not find matching generator level jet
            continue;     
         }
         fhRemx->Fill(jetRec->Pt(),jetGen->Pt());//response matrix
      }
   }
     

   return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskJetPP::Terminate(Option_t *){
   //Treminate 
   PostData(1, fOutput);

   // Mandatory
   fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1)); // '1' refers to the output slot
   if(!fOutput) {
      printf("ERROR: Output list not available\n");
      return;
   }
}

//________________________________________________________________________
AliAnalysisTaskJetPP::~AliAnalysisTaskJetPP(){
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   delete fHelperClass;
 
} 

/// Overloads base class method. Creates output objects
void AliAnalysisTaskJetPP::UserCreateOutputObjects(){
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TString name;

   //Histograms
   fhJetPt = new TH1D("fhJetPt","Pt spectrum of the AKT jets",320,-20.0,300.0);
   fhJetPtRho = new TH1D("fhJetPtRho","Pt spectrum of the AKT jets - kt bg",320,-20.0,300.0);
   fhJetPtConeRho = new TH1D("fhJetPtConeRho","p_{T} spectrum of the a-k_{T} jets - local kt bg",320,-20.0,300.0);
   fhJetPtCMSRho = new TH1D("fhJetPtCMSRho","p_{T} spectrum of the a-k_{T} jets - CMS bg",320,-20.0,300.0);
   fhAtimesRhoMedian = new TH1D("fhAtimesRhoMedian","a-k_{T} background times area",320,-20.0,300.0);
   fhAtimesRhoCone = new TH1D("fhAtimesRhoCone","cone background times area",320,-20.0,300.0);
   fhAtimesRhoCMS = new TH1D("fhAtimesRhoCMS","CMS a-k_{T} background times area",320,-20.0,300.0);
   fhJetConstituentPt = new TH1D("fhJetConstituentPt","Pt spectrum of the a-k_{T} jet constituents",300,0.0,300.0);
   fhTrackPt = new TH1D("fhTrackPt","p_{T} spectrum of tracks",300,0.0,300.0);
   fhJetAreaPt = new TH2D("fhJetAreaPt","Jet area vs pt",50,0.,2.,300,0.0,300.0);
   fhJetEtaPt = new TH2D("fhJetEtaPt","#eta-p_{T} distribution of jets",50,-fSignalJetEtaWindow,fSignalJetEtaWindow,300,0.0,300.0);
   fhAktJetEtaPhi = new TH2D("fhAktJetEtaPhi","Eta-Phi distribution of a-k_{T} jets",50,-fSignalJetEtaWindow,fSignalJetEtaWindow,50,-TMath::Pi(),TMath::Pi());
   fhKtJetEtaPhi = new TH2D("fhKtJetEtaPhi","#eta-#phi distribution of k_{T} jets",50,-fSignalJetEtaWindow,fSignalJetEtaWindow,50,-TMath::Pi(),TMath::Pi());
   fhTrackEtaPhi = new TH2D("fhTrackEtaPhi","#eta-#phi distribution of tracks",50,-fTrackEtaWindow,fTrackEtaWindow,50,-TMath::Pi(),TMath::Pi());
   fhTrackEtaPt = new TH2D("fhTrackEtaPt","#eta-p_{T} distribution of tracks",50,-fTrackEtaWindow,fTrackEtaWindow,300,0.0,300.0);
   fhGenTrackEtaPt = new TH2D("fhGenTrackEtaPt","#eta-p_{T} distribution of tracks",50,-fTrackEtaWindow,fTrackEtaWindow,300,0.0,300.0);
   fhJetPhiPt = new TH2D("fhJetPhiPt","#phi-p_{T} distribution of jets",50,-TMath::Pi(),TMath::Pi(),300,0.0,300.0);
   fhTrackPhiPt = new TH2D("fhTrackPhiPt","#phi-p_{T} distribution of tracks",50,-TMath::Pi(),TMath::Pi(),300,0.0,300.0);
   fhZVertex = new TH1D("fhZVertex","Z vertex",300,-15.0,15.0);
   fhZVertexBC = new TH1D("fhZVertexBC","Z vertex before the cut",200,-50.0,50.0);
   fhXVertex = new TH1D("fhXVertex","X vertex",300,-15.0,15.0);
   fhYVertex = new TH1D("fhYVertex","Y vertex",300,-15.0,15.0);
   fhCuts = new TH1D("fhCuts","Cuts statistics",3,0,3);
   fhRho = new TH1D("fhRho","k_{T} jet background",500,0,50.0);
   fhConeRho = new TH1D("fhConeRho","Local a-k_{T} jet background",500,0,50.0);
   fhCMSRho = new TH1D("fhCMSRho","CMS-style a-k_{T} jet background",500,0,50.0);
   fhKTJetPt = new TH1D("fhKTJetPt","KT jets p_{T} spectrum",300,0.0,300.0);
   fhPrimGenTrkPt = new TH1D("fhPrimGenTrkPt","MC track p_{T}",300,0.0,300.0);
   fhGenJetPt = new TH1D("fhGenJetPt","MC jet p_{T}",300,0.0,300.0);
   fhRecTrkPt = new TH1D("fhRecTrkPt","Correctly rec. track p_{T}",300,0.0,300.0);
   fhFakeTrkPt = new TH1D("fhFakeTrkPt","Fake track p_{T}",300,0.0,300.0);
   fhRemx = new TH2D("fhRemx","Response matrix",300,0.0,300.0,300,0.0,300.0);
   fhMult = new TH1D("fhMult","Reconstructed track multiplicity",300,0,300.0);

   fhTrackPhiCG = new TH2D("fhTrackPhiCG","Global tracks azimuth vs p_{T} ",300,0,300,50,-TMath::Pi(),TMath::Pi());
   fhTrackPhiTPCG = new TH2D("fhTrackPhiTPCG","Complementary tracks azimuth vs p_{T)",300,0,300.0,50,-TMath::Pi(),TMath::Pi());
   for(int i=0;i<2;i++){
      name = Form("fhSigmaPtOverPtVsPt%i",i); 
      fhSigmaPtOverPtVsPt[i] = new TH2D(name.Data(),"#sigma p_{T}/p_{T} vs p_{T}",100,0.0,100.0,300,0.0,3.0);
      fOutput->Add(fhSigmaPtOverPtVsPt[i]);
      name = Form("fhInvPtQVsPhi%i",i); 
      fhInvPtQVsPhi[i] = new TH2D(name.Data(),"Negative track azimuth vs 1/pt",50,-TMath::Pi(),TMath::Pi(),300,0.0,0.25);
      fOutput->Add(fhInvPtQVsPhi[i]);
      name = Form("fhInvPtQVsEta%i",i); 
      fhInvPtQVsEta[i] = new TH2D(name.Data(),"Negative track azimuth vs 1/pt",50,-fTrackEtaWindow,fTrackEtaWindow,300,0.0,0.25);
      fOutput->Add(fhInvPtQVsEta[i]); 
      name = Form("fhInvPtQVsPhiASide%i",i); 
      fhInvPtQVsPhiASide[i] = new TH2D(name.Data(),"Negative phi vs 1/pt A-eta",50,-TMath::Pi(),TMath::Pi(),300,0.0,0.25);
      fOutput->Add(fhInvPtQVsPhiASide[i]);
      name = Form("fhInvPtQVsPhiCSide%i",i); 
      fhInvPtQVsPhiCSide[i] = new TH2D(name.Data(),"Negative phi vs 1/pt C-eta",50,-TMath::Pi(),TMath::Pi(),300,0.0,0.25);
      fOutput->Add(fhInvPtQVsPhiCSide[i]); 
  }

   fOutput->Add(fhJetPt);
   fOutput->Add(fhJetPtRho);
   fOutput->Add(fhJetPtConeRho);
   fOutput->Add(fhJetPtCMSRho);
   fOutput->Add(fhAtimesRhoMedian);
   fOutput->Add(fhAtimesRhoCone);
   fOutput->Add(fhAtimesRhoCMS);
   fOutput->Add(fhJetAreaPt);

   fOutput->Add(fhJetConstituentPt);
   fOutput->Add(fhTrackPt);
   fOutput->Add(fhJetEtaPt);
   fOutput->Add(fhTrackEtaPt);
   fOutput->Add(fhJetPhiPt);

   fOutput->Add(fhTrackPhiPt);
   fOutput->Add(fhZVertex);
   fOutput->Add(fhXVertex);
   fOutput->Add(fhYVertex);
   fOutput->Add(fhCuts);

   fOutput->Add(fhRho);
   fOutput->Add(fhConeRho);
   fOutput->Add(fhCMSRho);
   fOutput->Add(fhKTJetPt);
   fOutput->Add(fhZVertexBC);
   fOutput->Add(fhPrimGenTrkPt);

   fOutput->Add(fhGenJetPt);
   fOutput->Add(fhRecTrkPt);
   fOutput->Add(fhFakeTrkPt);
   fOutput->Add(fhRemx);
   fOutput->Add(fhMult);
   fOutput->Add(fhGenTrackEtaPt);

   fOutput->Add(fhTrackPhiCG);
   fOutput->Add(fhTrackPhiTPCG);
   fOutput->Add(fhAktJetEtaPhi);
   fOutput->Add(fhKtJetEtaPhi);
   fOutput->Add(fhTrackEtaPhi);
   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutput->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
      if(hn){
         hn->Sumw2();
      }
   }
   TH1::AddDirectory(oldStatus);


   PostData(1, fOutput);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetPP::RetrieveEventObjects() {
   //
   // retrieve event objects
   //
   if(!AliAnalysisTaskEmcalJet::RetrieveEventObjects())  return kFALSE;
 
   return kTRUE;
}
/// Run the nalysis
Bool_t AliAnalysisTaskJetPP::Run()
{
   
   return kTRUE;
}

/// This function estimates the underlying event background based on the standard area based approach.
/// @param jetCont The jet container
/// @return rhoCone Background density
Double_t AliAnalysisTaskJetPP::EstimateBgKT(AliJetContainer *jetCont){
   Double_t rhoCone    = 0.0;
 
   if(!jetCont) return rhoCone;   

   AliEmcalJet*  jet        = NULL;
   static Double_t rhovec[999];
   Int_t nJetAcc = 0;
   Double_t jetpt;

   for(auto jetIterator : jetCont->accepted_momentum()) {
      jet=(AliEmcalJet*)jetIterator.second;
      if(!jet) continue;
      if(!IsSignalJetInAcceptance(jet,kFALSE)) continue; //check jet eta and pt

      jetpt = jet->Pt(); 
      if(jetpt <0.005) jetpt = 0.; //set pt of ghost jets identical to zero
      rhovec[nJetAcc] = jetpt/jet->Area();
      nJetAcc++;
   }

   if(nJetAcc>0){
      rhoCone = TMath::Median(nJetAcc, rhovec); //take the median of pTjet/Ajet
   }
 
  return rhoCone; 
}
/// This function estimates the udnerlying event background based on two perpendicula cones to the leading jet.
/// @param jetCont The jet container
/// @param trkCont The track container
/// @return rhoCone Background density
Double_t AliAnalysisTaskJetPP::EstimateLocalBg(AliJetContainer *jetCont,AliParticleContainer *trkCont){
   Double_t rhoCone    = 0.0;
 
   if(!jetCont) return rhoCone;   

   AliEmcalJet*  jet        = NULL;
   AliVParticle* track      = NULL;
   Double_t trackpt=0;
   Double_t tracketa=0;
   Double_t trackphi=0;
   Double_t coner = 0.4;//the radius of the perpendicular cone
   Double_t coner2 = coner*coner;
   Double_t highjetpt=-1;
   Double_t highjeteta=0;
   Double_t highjetphi=0;
   Double_t perpphi1=0;
   Double_t perpphi2=0;
   Double_t bgpt=0;

   Double_t X=0;
   Double_t Y=0;
   Double_t Z=0;

   for(auto jetIterator : jetCont->accepted_momentum()) {
      jet=(AliEmcalJet*)jetIterator.second;
      if(!jet) continue;
      if(!IsSignalJetInAcceptance(jet)) continue; //check jet eta and pt
      if(jet->Pt()>highjetpt){ 
	 highjetpt = jet->Pt(); 
	 highjeteta = jet->Eta();
	 highjetphi = Convert(jet->Phi());
      }
   }
   if(highjetpt<0) return 0.;
   perpphi1=Convert(highjetphi+TMath::Pi()/2);//the two perpendicular angles
   perpphi2=Convert(highjetphi-TMath::Pi()/2);

   for(auto trackIterator : trkCont->accepted_momentum()) {
      track=(AliVParticle*)trackIterator.second;
      if(!((AliAODTrack*)track)->IsHybridGlobalConstrainedGlobal()) continue; 
      if(!track) continue;
      if(!IsTrackInAcceptance(track)) continue;

      trackpt = track->Pt();
      tracketa = track->Eta();
      trackphi = Convert(track->Phi());
      X=tracketa-highjeteta;
      Y=Convert(trackphi-perpphi1);
      Z=Convert(trackphi-perpphi2);     

      if((X*X+Y*Y) < coner2 || (X*X+Z*Z) < coner2){
         bgpt=bgpt+trackpt;
      }
   }
   rhoCone = bgpt/(2*TMath::Pi()*coner2); //calculate the pT density
    
   return rhoCone; 
}
/// This function estimates the underlying event background based on the improved CMS method.
/// @param jetCont The jet container
/// @return rhoCMS Background density
Double_t AliAnalysisTaskJetPP::EstimateBgKTCMS(AliJetContainer *jetCont){
   Double_t rhoCMS    = 0.0;
 
   if(!jetCont) return rhoCMS;   

   AliEmcalJet*  jet        = NULL;
   static Double_t rhovec[999];
   Int_t nJetAcc = 0;
   Double_t jetpt = 0.0;
   Double_t physJetArea = 0.0;
   Double_t allJetArea = 0.0; 

   for(auto jetIterator : jetCont->accepted_momentum()) {
      jet=(AliEmcalJet*)jetIterator.second;
      if(!jet) continue;
      if(!IsSignalJetInAcceptance(jet,kFALSE)) continue; //check jet eta and pt

      if(jet->Pt() > 0.01){//Select only physical jets 
         physJetArea += jet->Area(); //Area of physical jets 
         jetpt = jet->Pt();          
         rhovec[nJetAcc] = jetpt/jet->Area();
         nJetAcc++;
      }
      allJetArea += jet->Area();//Area of all jets
   }

   if(nJetAcc>0 && allJetArea!=0){
      rhoCMS = TMath::Median(nJetAcc, rhovec)*physJetArea/allJetArea; //Calculate the CMS background
   }
 
  return rhoCMS; 
}
/// This function convertes the range of pseudorapidity to the interval (-Pi,Pi)
/// @param input Azimuth
/// @return phi Converted azimuth
Double_t AliAnalysisTaskJetPP::Convert(Double_t input){//Converts an angle to the range (-Pi,Pi)
   Double_t phi;
   phi = TMath::ATan2(TMath::Sin(input),TMath::Cos(input));
   return phi;
}
