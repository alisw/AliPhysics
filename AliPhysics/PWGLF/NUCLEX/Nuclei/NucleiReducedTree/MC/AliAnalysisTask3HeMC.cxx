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
//
// Task for the analysis of the production of 3He as a function of multiplicity
//
// Author:
//  Sebastian Hornung <Sebastian.Hornung@cern.ch>

// Root header
#include <TChain.h>
#include <TList.h>

//AliRoot header
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"

#include "TH1F.h" // for dummy
#include "TH2D.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"

//AliPhysics header
#include "AliAnalysisTask3HeMC.h"

//____________________________________________________________
AliAnalysisTask3HeMC::AliAnalysisTask3HeMC():
AliAnalysisTaskSE(),
kAODanalysis(kFALSE),
kHasMCdata(kFALSE),
kpp(kFALSE),
kpPb(kFALSE),
kPbPb(kFALSE),
fEtaMin(-0.9),
fEtaMax(0.9),
fMaxDCAxy(2),
fMaxDCAz(3),
fRejectKinkMother(kFALSE),
Multiplicity(-1),
MultiplicityPercentile(-1),
kTPCnSigmaCut(-5.0),
fHistMultEstZDep(NULL),
fHistTPCdEdxRigidity(NULL),
fHistTPCsigHe3(NULL),
fHistTPCsigHe4(NULL),
fHistTPCdEdxSigmaHe3(NULL),
fHistTPCdEdxSigmaHe4(NULL),
fHistPtTrueRecHe3(NULL),
fHistPtTrueRecHe4(NULL),
fHistTrueHe3(NULL),
fHistTrueHe4(NULL),
fEventCuts(),
fPIDResponse(NULL),
fUtils(NULL),
fOutput(NULL),
fQAList(NULL),
fAODMCHeader(NULL),
fAODArrayMCParticles(NULL)
{
}

//____________________________________________________________
AliAnalysisTask3HeMC::AliAnalysisTask3HeMC(const char *name):
AliAnalysisTaskSE(name),
kAODanalysis(kFALSE),
kHasMCdata(kFALSE),
kpp(kFALSE),
kpPb(kFALSE),
kPbPb(kFALSE),
fEtaMin(-0.9),
fEtaMax(0.9),
fMaxDCAxy(2),
fMaxDCAz(3),
fRejectKinkMother(kFALSE),
Multiplicity(-1),
MultiplicityPercentile(-1),
kTPCnSigmaCut(-5.0),
fHistMultEstZDep(NULL),
fHistMultValueVsPercientile(NULL),
fHistTPCdEdxRigidity(NULL),
fHistTPCsigHe3(NULL),
fHistTPCsigHe4(NULL),
fHistTPCdEdxSigmaHe3(NULL),
fHistTPCdEdxSigmaHe4(NULL),
fHistPtTrueRecHe3(NULL),
fHistPtTrueRecHe4(NULL),
fHistTrueHe3(NULL),
fHistTrueHe4(NULL),
fEventCuts(),
fPIDResponse(NULL),
fUtils(NULL),
fOutput(NULL),
fQAList(NULL),
fAODMCHeader(NULL),
fAODArrayMCParticles(NULL)
{
   
   fUtils = new AliAnalysisUtils();
   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
   DefineOutput(2, TList::Class());
}

//____________________________________________________________
AliAnalysisTask3HeMC::~AliAnalysisTask3HeMC(){
   //
   // Destructor
   //
   
   // Delete output objects only if we are not running in PROOF mode because otherwise this produces a crash during merging
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(mgr && mgr->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
      if(fOutput) delete fOutput;
      if(fQAList) delete fQAList;
      if(fUtils) delete fUtils;
   }
}

//___________________________________________________
void AliAnalysisTask3HeMC::UserCreateOutputObjects(){ //Objects that are output are initialized in the function
   AliDebug(4, "Creating Output Objects");
   
   // Make lists for Output
   if(!fOutput) fOutput = new TList();
   fOutput->SetOwner();
   
   if(!fQAList) fQAList = new TList();
   fQAList->SetOwner();
   
   //   // Automatic determination of the analysis mode
   //   AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   //   if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
   //      SetAODAnalysis();
   //   } else {
   //      SetESDAnalysis();
   //      if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) // possible for AOD too?
   //         SetHasMCData();
   //   }
   // Analysis mode
   AliInfo(Form("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD"));
   AliInfo(Form("MC Data available %s\n", HasMCData() ? "Yes" : "No"));
   
   fEventCuts.AddQAplotsToList(fQAList); /// Add event selection QA plots
   
   // QA histograms
   fHistMultEstZDep = new TH2D("fHistMultEstZDep","Multiplicity vs z_{vertex}", 100,0, 100, 201, -10, 10);
   fHistMultEstZDep->GetXaxis()->SetTitle("Multiplicity");
   fHistMultEstZDep->GetYaxis()->SetTitle("z_{Prim. vertex}");
   fHistMultEstZDep->Sumw2();
   fQAList->Add(fHistMultEstZDep);

   fHistMultValueVsPercientile = new TH2D("fHistMultValueVsPercientile","Multiplicity vs z_{vertex}",100, 0, 1000, 100, 0, 100);
   fHistMultValueVsPercientile->GetXaxis()->SetTitle("Multiplicity");
   fHistMultValueVsPercientile->GetYaxis()->SetTitle("Multiplicity Percentile");
   fHistMultValueVsPercientile->Sumw2();
   fQAList->Add(fHistMultValueVsPercientile);

   fHistTPCdEdxRigidity = new TH2D("fHistTPCdEdxRigidity","TPC dE/dx distribution",100,0.1,10,150, 0, 1500);
   fHistTPCdEdxRigidity->GetXaxis()->SetTitle("#it{p/z} (GeV/#it{c})");
   fHistTPCdEdxRigidity->GetYaxis()->SetTitle("TPC d#it{E}/d#it{x} (a.u.)");
   fHistTPCdEdxRigidity->Sumw2();
   fQAList->Add(fHistTPCdEdxRigidity);
   
   fHistTPCsigHe3 = new TProfile2D("fHistTPCsigHe3", "hist to check cut on TPC sigma for He3", 100, 0.1, 6.0, 1000, 0.0, 1000.);
   fQAList->Add(fHistTPCsigHe3);
   fHistTPCsigHe4 = new TProfile2D("fHistTPCsigHe4", "hist to check cut on TPC sigma for He4", 100, 0.1, 6.0, 1000, 0.0, 1000.);
   fQAList->Add(fHistTPCsigHe4);
   
   //Result
   const Int_t nDim=11;
   const Int_t nBinMultiplicity=100;
   const Int_t nBinRigidity=50;
   const Int_t nBinTPCdEdx=80;
   const Int_t nBinPt=15;
   const Int_t nBinCharge=2;
   const Int_t nBinDCAxy=120;
   const Int_t nBinDCAz=300;
   const Int_t nBinRapidity=20;
   const Int_t nPtBins=15;
   const Int_t nITSBins=6;
   const Int_t nTPCclusterBins=4;

   //										mult  Rig	dE/dx	pT		Charge	DCAxy	DCAz	y	source	ITShits	TPCcluster
   Double_t BinEdge_min[nDim] = {0,		0,		-4,	0.45,	-2,	-0.6,	-1.5,	-1,	-0.5,		0.5,	59.9};
   Double_t BinEdge_max[nDim] = {100,	5,		4,		7.95,	2,		0.6,	1.5,	1,		2.5,		6.5,	79.9};
   Int_t nBin[nDim] = {nBinMultiplicity,nBinRigidity,nBinTPCdEdx,nBinPt,nBinCharge, nBinDCAxy, nBinDCAz, nBinRapidity, 3,nITSBins,nTPCclusterBins};
   fHistTPCdEdxSigmaHe3 =new THnSparseD("fHistTPCdEdxSigmaHe3","TPC dE/dx distribution [sigma_{3He}]",nDim,nBin,BinEdge_min,BinEdge_max);
   fHistTPCdEdxSigmaHe3->GetAxis(0)->SetTitle("Multiplicity");
   fHistTPCdEdxSigmaHe3->GetAxis(1)->SetTitle("#it{p/z} (GeV/#it{c})");
   fHistTPCdEdxSigmaHe3->GetAxis(2)->SetTitle("TPC d#it{E}/d#it{x} (#sigma_{3He})");
   fHistTPCdEdxSigmaHe3->GetAxis(3)->SetTitle("#it{p}_{T}^{corr} (GeV/#it{c})");
   fHistTPCdEdxSigmaHe3->GetAxis(4)->SetTitle("(Anit-)Particle");
   fHistTPCdEdxSigmaHe3->GetAxis(5)->SetTitle("DCA in xy (cm)");
   fHistTPCdEdxSigmaHe3->GetAxis(6)->SetTitle("DCA in z (cm)");
   fHistTPCdEdxSigmaHe3->GetAxis(7)->SetTitle("Rapidity (lab frame)");
   fHistTPCdEdxSigmaHe3->GetAxis(8)->SetTitle("Source");
   fHistTPCdEdxSigmaHe3->GetAxis(9)->SetTitle("ITS hits");
   fHistTPCdEdxSigmaHe3->GetAxis(10)->SetTitle("TPC cluster");
   fHistTPCdEdxSigmaHe3->Sumw2();
   
   fHistTPCdEdxSigmaHe4 = (THnSparseD*) fHistTPCdEdxSigmaHe3->Clone("fHistTPCdEdxSigmaHe4");
   fHistTPCdEdxSigmaHe4->SetTitle("TPC dE/dx distribution [sigma_{4He}]");
   fHistTPCdEdxSigmaHe4->GetAxis(2)->SetTitle("TPC d#it{E}/d#it{x} (#sigma_{4He})");
   
   fOutput->Add(fHistTPCdEdxSigmaHe3);
   fOutput->Add(fHistTPCdEdxSigmaHe4);
   
   if(HasMCData()){
      fHistPtTrueRecHe3 = new TH2D("fHistPtTrueRecHe3", "(#it{p}_{T,true}- #it{p}_{T,rec}) vs #it{p}_{T,rec}", 100, 0.1, 10.1, 120, -3, 3);
      fHistPtTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{T,rec} (GeV/#it{c})");
      fHistPtTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{T,true}- #it{p}_{T,rec}) (GeV/#it{c})");
      fOutput->Add(fHistPtTrueRecHe3);

      fHistPtTrueRecHe4 = (TH2D*) fHistPtTrueRecHe3->Clone("fHistPtTrueRecHe4");
      fHistPtTrueRecHe4->SetName("fHistPtTrueRecHe4");
      fOutput->Add(fHistPtTrueRecHe4);

      const Int_t nDimMC=5;
      //                           mult        pTtrue    Source        Charge      y
      Double_t BinEdgeMC_min[nDim] = {0, 			0.45,	  	-0.5, 	   	-2,	      -1};
      Double_t BinEdgeMC_max[nDim] = {100,		7.95,	  	2.5,   	  		2,		       1};
      Int_t nBinMC[nDimMC] = {nBinMultiplicity,nBinPt,3,nBinCharge, nBinRapidity};
      
      fHistTrueHe3 = new THnSparseD("fHistTrueHe3", "MC primary He3 paritcle", nDimMC, nBinMC, BinEdgeMC_min, BinEdgeMC_max);
      fHistTrueHe3->GetAxis(0)->SetTitle("Multiplicity");
      fHistTrueHe3->GetAxis(1)->SetTitle("#it{p}_{T}^{corr} (GeV/#it{c})");
      fHistTrueHe3->GetAxis(2)->SetTitle("Source");
      fHistTrueHe3->GetAxis(3)->SetTitle("(Anit-)Particle");
      fHistTrueHe3->GetAxis(4)->SetTitle("Rapidity (lab frame)");
      
      fHistTrueHe4 = (THnSparseD*) fHistTrueHe3->Clone("fHistTrueHe4");
      fHistTrueHe4->SetTitle("MC primary He4 paritcle");
      
      fOutput->Add(fHistTrueHe3);
      fOutput->Add(fHistTrueHe4);
   }

   PostData(1,fOutput);
   PostData(2,fQAList);
}

//___________________________________________________
void AliAnalysisTask3HeMC::UserExec(Option_t *){ //called for each event
   
   AliDebug(4, "Starting Single Event Analysis");
   if(!fInputEvent){
      AliError("Reconstructed Event not available");
      return;
   }
   
   if(IsESDanalysis() && HasMCData()){
      AliDebug(3, Form("MC Event: %p", fMCEvent));
      if(!fMCEvent){
         AliError("No MC Event, but MC Data required");
         return;
      }
      // Protect against missing MC trees
      AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if(!mcH){
         AliError("No MC Event Handler available");
         return;
      }
      if(!mcH->InitOk()) return;
      if(!mcH->TreeK()) return;
      if(!mcH->TreeTR()) return;
   }
   
   if(IsAODanalysis() && HasMCData()){
      //       take MC info
      fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
      if(!fAODMCHeader){
         AliError("No AliAODMCHeader");
         return;
      }
      fAODArrayMCParticles = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!fAODArrayMCParticles){
         AliError("No AOD MC particles");
         return;
      }
   }
   
   //PID response
   fPIDResponse = fInputHandler->GetPIDResponse();
   if(!fPIDResponse) {
      AliError("No PID Response found");
      return;
   }
   
   //Event Cut
   if (!fEventCuts.AcceptEvent(fInputEvent)) {
      PostData(2, fQAList);
      return;
   }
   
   // Multiplicity / centrality framework
   AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
   if(!MultSelection->IsEventSelected()) return;
   if (MultSelection){
      AliMultEstimator* MultiplicityEstimator = NULL;

      // Percentile
      MultiplicityPercentile = (Double_t) MultSelection->GetMultiplicityPercentile("V0A"); // SPDTracklets V0A V0M
      // Value
      MultiplicityEstimator = MultSelection->GetEstimator("V0A"); // Multiplicity estimater based on SPD tracklets
      Multiplicity = (Double_t) MultiplicityEstimator->GetValue(); // returns multiplicity only for events used in calibration

      fHistMultValueVsPercientile->Fill(Multiplicity, MultiplicityPercentile);
   }else{
      //If this happens, re-check if AliMultSelectionTask ran before your task!
      AliInfo("AliMultSelection object not found!");
      return;
   }
   Double_t zVertex = fInputEvent->GetPrimaryVertex()->GetZ();
   fHistMultEstZDep->Fill(MultiplicityPercentile, zVertex);
   
   // Run analysis for AOD
   if (IsAODanalysis()) {
      ProcessAOD();
   }
   
   PostData(1,fOutput);
   PostData(2,fQAList);
}

//____________________________________________________________
void AliAnalysisTask3HeMC::Terminate(Option_t *){
   //
   // Terminate not implemented at the moment
   //
}

//___________________________________________________
void AliAnalysisTask3HeMC::ProcessAOD(){ // Loop over tracks and perform analysis
   AliDebug(3, "Processing AOD Event");
   
   AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
   if(!fAOD){
      AliError("AOD Event required for AOD Analysis");
      return;
   }
   
   Int_t fNumberOfVertices = fAOD->GetNumberOfVertices();
   Double_t fListOfmotherkink[fNumberOfVertices];
   Int_t fNumberOfMotherkink = 0;
   for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++){
      AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
      if(!aodvertex) continue;
      if(aodvertex->GetType()==AliAODVertex::kKink){
         AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
         if(!mother) continue;
         Int_t idmother = mother->GetID();
         fListOfmotherkink[fNumberOfMotherkink] = idmother;
         fNumberOfMotherkink++;
      }
   }

   //loop over MC particles
   if (HasMCData()) {
      for (Int_t iPartMC=0; iPartMC <= fAODArrayMCParticles->GetLast(); iPartMC++) {
         AliAODMCParticle* mcpart = (AliAODMCParticle*) fAODArrayMCParticles->At(iPartMC);
         if (!mcpart) {
            AliError("No MC particle found");
            continue;
         } else {
            Int_t ParticleType = 0;
            if (TMath::Abs(mcpart->GetPdgCode())== 1000020030) ParticleType = 1; //He3
            if (TMath::Abs(mcpart->GetPdgCode())== 1000020040) ParticleType = 2; //He4
            //            if (TMath::Abs(mcpart->GetPdgCode())== 1000010030) ParticleType = 3; //triton
            if (ParticleType == 0) continue;

            Double_t TruePt = mcpart->Pt();
            Double_t Charge = (Double_t) mcpart->Charge() / (Double_t) 3.; // Charge is given in fractional charge (like quarks) 3 -> 1e
            Double_t RapidityMC = mcpart->Y();
            // Correction to CMS rapidity
            RapidityMC = RapidityMC - 0.465;

            Double_t Source = -1;
            if (mcpart->IsPhysicalPrimary()){
               Source=0;
            } else if (mcpart->IsSecondaryFromMaterial()){
               Source=1;
            } else if (mcpart->IsSecondaryFromWeakDecay()){
               Source=2;
            }

            Double_t FillSparseMC[] = {MultiplicityPercentile,TruePt,Source,Charge/2,RapidityMC};
            if (ParticleType==1) {
               fHistTrueHe3->Fill(FillSparseMC);
            } else if (ParticleType ==2){
               fHistTrueHe4->Fill(FillSparseMC);
            }
         }
      }
   }

   Int_t nTracks(fAOD->GetNumberOfTracks());
   for(Int_t iTrack=0; iTrack < nTracks; iTrack++) { // loop over all the tracks
      
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
      if(!track) continue;

      Double_t DCAxy = GetDCAxy(track);
      Double_t DCAz = GetDCAz(track);

      Int_t ParticleType = 0; // needed since He3 and He4 are not separated in the TPC for high pT
      Double_t TruePt = -1;
      Double_t Source = -1;

      if (HasMCData()) {
         AliAODMCParticle* mcpart = (AliAODMCParticle*) fAODArrayMCParticles->At(TMath::Abs(track->GetLabel()));
         if (!mcpart) {
            AliError("No MC particle found");
            continue;
         } else {
            if (TMath::Abs(mcpart->GetPdgCode())== 1000020030) ParticleType = 1; //He3
            if (TMath::Abs(mcpart->GetPdgCode())== 1000020040) ParticleType = 2; //He4
            //            if (TMath::Abs(mcpart->GetPdgCode())== 1000010030) ParticleType = 3; //triton
            if (ParticleType == 0) continue;

            TruePt = mcpart->Pt();

            if (mcpart->IsPhysicalPrimary()){
               Source=0;
            } else if (mcpart->IsSecondaryFromMaterial()){
               Source=1;
            } else if (mcpart->IsSecondaryFromWeakDecay()){
               Source=2;
            }
         }
      }

      // Geomatrical acceptance
      if((track->Eta() < fEtaMin) || (track->Eta() > fEtaMax)) continue;

      // check AOD filter bit
      if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

      // Reject kink mothers
      if(fRejectKinkMother){
         Bool_t kinkmotherpass = kTRUE;
         for(Int_t kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++){
            if(track->GetID() == fListOfmotherkink[kinkmother]){
               kinkmotherpass = kFALSE;
               continue;
            }
         }
         if(!kinkmotherpass) continue;
      }

      // additional track cuts
      if(TMath::Abs(DCAxy) > fMaxDCAxy) continue;
      if(TMath::Abs(DCAz) > fMaxDCAz) continue;

      if (track->GetTPCNcls() < 60) continue;
      if(track->GetITSNcls() < 1) continue;

      Double_t ITShits = (Double_t) track->GetITSNcls();
      Double_t TPCcluster = (Double_t) track->GetTPCNcls();


      // Get particle properties
      Double_t TrackPt = 2*track->Pt(); // use function for pT correction
      Double_t Charge = (Double_t) track->Charge(); // not the physical charge is always -1 or 1; more like sign of charge
      Double_t Rigidity = TMath::Abs(track->P());
      Double_t dEdx = track->GetTPCsignal();

      // Fill TPC QA histogram
      fHistTPCdEdxRigidity->Fill(Rigidity,dEdx);

      //Fill TPC nSigma histogram He3
      Double_t dEdxSigmaTriton= fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
      Double_t dEdxSigmaHe3= fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
      Double_t dEdxSigmaHe4= fPIDResponse->NumberOfSigmasTPC(track, AliPID::kAlpha);

      // basic PID cut
      Double_t Rapidity = 0;

      // For particles with charge 2
      Double_t FillSparse[] = {MultiplicityPercentile,Rigidity,dEdxSigmaHe3,TrackPt,Charge,DCAxy,DCAz,Rapidity,Source, ITShits, TPCcluster};

      if (dEdxSigmaTriton > kTPCnSigmaCut){

         Rapidity = CalculateRapidity(track->Px(), track->Py(), track->Pz(), "3He");

         //Fill 3He candidates
         if (ParticleType==1){
            if (track->GetTPCNcls() >= 70 && track->GetITSNcls() > 1 && HasMCData() && TMath::Abs(dEdxSigmaHe3) < 3){
               if (Rapidity < 0 && Rapidity > -1. ) {
                  fHistPtTrueRecHe3->Fill(TrackPt, (TruePt-TrackPt));
               }
            }

            fHistTPCsigHe3->Fill(TrackPt,   dEdx, dEdxSigmaHe3);
            FillSparse[2] = dEdxSigmaHe3;
            FillSparse[7] = CalculateRapidity(track->Px(), track->Py(), track->Pz(), "3He");
            FillSparse[3] = 2*track->Pt() + 0.00068 - 0.3641*TMath::Exp(-0.2386*TMath::Power(2*track->Pt(),3.0));
            fHistTPCdEdxSigmaHe3->Fill(FillSparse);
         }

         //Fill 4He candidates
         if (ParticleType==2){
            if (track->GetTPCNcls() >= 70 && track->GetITSNcls() > 1 && HasMCData() && TMath::Abs(dEdxSigmaHe4) < 3){
               if (Rapidity < 0 && Rapidity > -1. ) {
                  fHistPtTrueRecHe4->Fill(TrackPt, (TruePt-TrackPt));
               }
            }

            fHistTPCsigHe4->Fill(TrackPt,   dEdx, dEdxSigmaHe4);
            FillSparse[2] = dEdxSigmaHe4;
            FillSparse[7] = CalculateRapidity(track->Px(), track->Py(), track->Pz(), "4He");
            FillSparse[3] = 2*track->Pt(); // Correction not yet available
            fHistTPCdEdxSigmaHe4->Fill(FillSparse);
         }
      }
   }
}

//________________________________________________________________________
void AliAnalysisTask3HeMC::BinLogAxis(TAxis *axis) {
   //
   // Method for the correct logarithmic binning of histograms
   //
   int bins = axis->GetNbins();

   Double_t from = axis->GetXmin();
   Double_t to = axis->GetXmax();
   Double_t *newBins = new Double_t[bins + 1];

   newBins[0] = from;
   Double_t factor = pow(to/from, 1./bins);

   for (int i = 1; i <= bins; i++) {
      newBins[i] = factor * newBins[i-1];
   }
   axis->Set(bins, newBins);
   delete [] newBins;
}

Double_t AliAnalysisTask3HeMC::GetDCAxy (AliAODTrack *track)  {

   AliAODEvent *fAODevent = dynamic_cast<AliAODEvent *>(fInputEvent);

   //   Double_t v[3];
   //   Double_t pos[3];
   //
   //   AliAODVertex *vertex = fAODevent->GetPrimaryVertex();
   //   vertex->GetXYZ(v);
   //   track->GetXYZ(pos);
   //
   //   Float_t DCAx = pos[0] - v[0];
   //   Float_t DCAy = pos[1] - v[1];
   //   Float_t DCAxy = TMath::Sqrt(DCAx*DCAx + DCAy*DCAy);

   Double_t impactParameter[2];
   Double_t covarianceMatrix[3];
   if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

   Double_t DCAxy = impactParameter[0];

   return DCAxy;
}

Double_t AliAnalysisTask3HeMC::GetDCAz (AliAODTrack *track)  {

   AliAODEvent *fAODevent = dynamic_cast<AliAODEvent *>(fInputEvent);

   //   Double_t v[3];
   //   Double_t pos[3];
   //
   //   AliAODVertex *vertex = fAODevent->GetPrimaryVertex();
   //   vertex->GetXYZ(v);
   //   track->GetXYZ(pos);
   //
   //   Float_t DCAz = TMath::Abs(pos[2] - v[2]);

   Double_t impactParameter[2];
   Double_t covarianceMatrix[3];
   if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

   Double_t DCAz = impactParameter[1];

   return DCAz;
}
//______________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTask3HeMC::CalculateRapidity(Double_t px, Double_t py, Double_t pz, TString type){
   // calculation of the rapidity
   Double_t p = TMath::Sqrt(px*px + py*py+ pz*pz);
   Double_t energy = 0;
   if (type == "4He") {
      energy = TMath::Sqrt(p*2*p*2 + AliPID::ParticleMass(AliPID::kAlpha)*AliPID::ParticleMass(AliPID::kAlpha));
   } else {
      energy = TMath::Sqrt(p*2*p*2 + AliPID::ParticleMass(AliPID::kHe3)*AliPID::ParticleMass(AliPID::kHe3));
   }
   Double_t rap = 0.5*TMath::Log((energy + pz*2)/(energy - pz*2));
   // Correction to CMS rapidity
   rap = rap - 0.465;

   return rap;
}
