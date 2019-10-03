/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                 Dealing with Wide pT Window Modified to ESDs            //
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliHelperPID.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliTracker.h"


using std::endl;
using std::cout;

#include "AliEbyEPidTTask.h"

ClassImp(AliEbyEPidTTask)

//-----------------------------------------------------------------------
AliEbyEPidTTask::AliEbyEPidTTask( const char *name ) : 
AliAnalysisTaskSE( name ), 
  fArrayMC(NULL),
  fESDtrackCuts(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),

  fThnList(NULL), 
  fAODtrackCutBit(128),
  fHelperPID(0x0),
  fEventCounter(NULL), 
  fPidCont(0x0),
  fHistCent(0x0),  
  fHistPt00(0x0),
  fHistPt10(0x0),
  fHistPt20(0x0),
  fHistPt30(0x0),
  fHistPt01(0x0),
  fHistPt11(0x0),
  fHistPt21(0x0),
  fHistPt31(0x0),

  fVxMax(3.), 
  fVyMax(3.), 
  fVzMax(10.), 
  fPtMin(0.2),   
  fPtMax(3.), 
  fEtaMin(-1.), 
  fEtaMax(1.),  
  fDcaXy(10.),
  fDcaZ(10.),  
  
  fIsMC(kFALSE),
  fIsAOD(kFALSE),
  fDebug(kFALSE),
  fIsQa(kFALSE),
  fIsTrig(kFALSE),

  fRunNumber(0),
  fNumberOfTracks(15000),
  fNumberOfTracksM(15000),
  fNTracks(0) {
   
  DefineOutput(1, TList::Class()); //! Connect Outpput....
  DefineOutput(2, TTree::Class()); //! Connect Outpput....
}

AliEbyEPidTTask::~AliEbyEPidTTask() {
  //!   Cleaning up
  if (fThnList)   delete fThnList;
  if (fHelperPID) delete fHelperPID;
  if (fPidCont)    delete fPidCont;
}

//---------------------------------------------------------------------------------
void AliEbyEPidTTask::UserCreateOutputObjects() {
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);

 if (!fIsAOD) {
    if(!fESDtrackCuts)
      fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    else 
      Printf(" >>>>  User Track Cuts <<<< ");
    fESDtrackCuts->Print();
    
    Printf(" >>> DCAxy in TC [%8.4f:%8.4f]", 
	   fESDtrackCuts->GetMinDCAToVertexXY(), fESDtrackCuts->GetMaxDCAToVertexXY());
    Printf(" >>> DCAz in TC  [%8.4f:%8.4f]", 
	   fESDtrackCuts->GetMinDCAToVertexZ(), fESDtrackCuts->GetMaxDCAToVertexZ());
	       
    Float_t r1,r2;
    fESDtrackCuts->GetPtRange(r1,r2);
    Printf(" >>> Pt in TC  [%10.4f:%10.4f]",r1,r2);

    fESDtrackCuts->GetRapRange(r1,r2);
    Printf(" >>> Rap in TC [%10.4f:%10.4f]",r1,r2);

    fESDtrackCuts->GetEtaRange(r1,r2);
    Printf(" >>> Eta in TC [%10.4f:%10.4f]",r1,r2);
  }     
 
 const Int_t NptBins = 52;
  Double_t ptBin[NptBins + 1] = {0.05, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
 
 
  fHistPt00 = new TH2F("fHistPtChargeMinus","Centrality vs p_{T} : N_{ch};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  fHistPt10 = new TH2F("fHistPtPionMinus","Centrality vs p_{T} : N_{#pi};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  fHistPt20 = new TH2F("fHistPtKaonMinus","Centrality vs p_{T} : N_{k};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  fHistPt30 = new TH2F("fHistPtProtonMinus","Centrality vs p_{T} : N_{p};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  
  fHistPt01 = new TH2F("fHistPtChargePlus","Centrality vs p_{T} : N_{ch};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  fHistPt11 = new TH2F("fHistPtPionPlus","Centrality vs p_{T} : N_{#pi};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  fHistPt21 = new TH2F("fHistPtKaonPlus","Centrality vs p_{T} : N_{k};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  fHistPt31 = new TH2F("fHistPtProtonPlus","Centrality vs p_{T} : N_{p};#it{Bin};p_{T}", 20,-0.5,19.5,NptBins, ptBin);
  
  fThnList->Add(fHistPt00);
  fThnList->Add(fHistPt10);
  fThnList->Add(fHistPt20);
  fThnList->Add(fHistPt30);

  fThnList->Add(fHistPt01);
  fThnList->Add(fHistPt11);
  fThnList->Add(fHistPt21);
  fThnList->Add(fHistPt31);

  fHistCent = new TH1F("fHistCentPid","Centrality",20,-0.5,19.5);			 
  fThnList->Add(fHistCent);
  

  fEventCounter = new TH1D("fEventCounter","EventCounter", 100, -0.5,99.5);
  fThnList->Add(fEventCounter);
  
   TDirectory *owd = gDirectory;
  OpenFile(1);
  fPidCont = new TTree("Event","fPidCont B");
  owd->cd();

  fPidCont->Branch("fRunNumber", &fRunNumber,  "fRunNumber/I");
  fPidCont->Branch("cent",fCentrality,"fCentrality[6]/F");
  if (fIsTrig) fPidCont->Branch("Trigger",fTrigMask,  "fTrigMask[5]/I");
  fPidCont->Branch("vertex",fVtx,"fVtx[3]/F");
  fPidCont->Branch("fNumberOfTracks", &fNumberOfTracks,"fNumberOfTracks/I");
  fPidCont->Branch("fTrackPt",fTrackPt,"fTrackPt[fNumberOfTracks]/F");

  fPidCont->Branch("fTrackPhi",fTrackPhi,"fTrackPhi[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackEta",fTrackEta,"fTrackEta[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackDxy",fTrackDxy,"fTrackDxy[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackDz",fTrackDz,"fTrackDz[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackPid",fTrackPid,"fTrackPid[fNumberOfTracks]/I");
  fPidCont->Branch("fTrackTpcNcl",fTrackTpcNcl,"fTrackTpcNcl[fNumberOfTracks]/I");
  fPidCont->Branch("fTrackCnDf",fTrackCnDf,"fTrackCnDf[fNumberOfTracks]/F");
  fPidCont->Branch("fTrackBit",fTrackBit,"fTrackBit[fNumberOfTracks]/I");
  
  if (fHelperPID) {
    fThnList->Add(new TList);
    TList *list =  static_cast<TList*>(fThnList->Last());
    list->SetName("HelperPID");
    list->SetOwner(kTRUE);
    TList *ll = (TList*)fHelperPID->GetOutputList();
    for (Int_t ikey = 0; ikey < ll->GetEntries(); ikey++) {
      list->Add(ll->At(ikey));
    }
  }

  
  if(fIsMC) {
    fPidCont->Branch("fTrackLabel",fTrackLabel,"fTrackLabel[fNumberOfTracks]/I");
    fPidCont->Branch("fPidStat",fPidStat,"fPidStat[fNumberOfTracks]/I");
    
    fPidCont->Branch("fNumberOfTracksM", &fNumberOfTracksM,"fNumberOfTracksM/I");
    fPidCont->Branch("fTrackPtM",fTrackPtM,"fTrackPtM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackLabelM",fTrackLabelM,"fTrackLabelM[fNumberOfTracksM]/I");
    fPidCont->Branch("fTrackPhiM",fTrackPhiM,"fTrackPhiM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackEtaM",fTrackEtaM,"fTrackEtaM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackPidM",fTrackPidM,"fTrackPidM[fNumberOfTracksM]/I");
  }

  PostData(1, fThnList);
  PostData(2, fPidCont);  
}

//----------------------------------------------------------------------------------
void AliEbyEPidTTask::UserExec( Option_t * ){
  fEventCounter->Fill(1);
  AliVEvent *event = InputEvent();
  if (!event) return;
  AliInputEventHandler* fInputEventHandler = static_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fInputEventHandler) return;
  
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if(!vertex) return;
  
  Bool_t vtest = kFALSE;
  Double32_t fCov[6];
  vertex->GetCovarianceMatrix(fCov);
  if(vertex->GetNContributors() > 0) {
    if(fCov[5] != 0) {
      vtest = kTRUE;
    }
  }
  if(!vtest)return;
  
  if(TMath::Abs(vertex->GetX()) > fVxMax) return;
  fVtx[0] = vertex->GetX();
  if(TMath::Abs(vertex->GetY()) > fVyMax) return;
  fVtx[1] = vertex->GetY();
  if(TMath::Abs(vertex->GetZ()) > fVzMax) return;
  fVtx[2] = vertex->GetZ();
  
  AliCentrality *centrality = event->GetCentrality();
  
  if (centrality->GetQuality() != 0) return;
  
  fRunNumber = event->GetRunNumber();
  if (fIsTrig) {
    fTrigMask[0] = 0;  
    if ((fInputEventHandler->IsEventSelected() & AliVEvent::kMB))          
      fTrigMask[0] = 1;
    fTrigMask[1] = 0;  
    if ((fInputEventHandler->IsEventSelected() & AliVEvent::kCentral))     
      fTrigMask[1] = 1;
    fTrigMask[2] = 0;  
    if ((fInputEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) 
      fTrigMask[2] = 1;
    fTrigMask[3] = 0;  
    if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      
      fTrigMask[3] = 1;
    fTrigMask[4] = 0;  
    if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      
      fTrigMask[4] = 1;
  } else {
    if (!(fInputEventHandler->IsEventSelected() & AliVEvent::kMB))  return;        
  }
  
  // AliAODHeader *header = (AliAODHeader*) event->GetHeader();
  // if (!header) return;

  fCentrality[0] = centrality->GetCentralityPercentile("V0M");
  fCentrality[1] = centrality->GetCentralityPercentile("CL1");
  fCentrality[2] = centrality->GetCentralityPercentile("TRK");
    
  fCentrality[3] = event->GetVZEROData()->GetMTotV0A();
  fCentrality[4] = event->GetVZEROData()->GetMTotV0C();
  fCentrality[5] = event->GetZDCData()->GetZDCParticipants();
 
  fHistCent->Fill(fCentrality[0]);

  // fCentrality[3] = header->GetCentralityP()->GetCentralityPercentile("V0M");
  // fCentrality[4] = header->GetCentralityP()->GetCentralityPercentile("V0A");
  // fCentrality[5] = header->GetCentralityP()->GetCentralityPercentile("V0C");
  // fCentrality[6] = header->GetCentralityP()->GetCentralityPercentile("FMD");
  // fCentrality[7] = header->GetCentralityP()->GetCentralityPercentile("V0MvsFMD");

  
  //  Printf("%f %f %f  %f %f  %f %f: %f %f", fCentrality[0],fCentrality[1],
  //	  fCentrality[2], fCentrality[3],fCentrality[4], fCentrality[5],fCentrality[6], fCentrality[7],fVtx[2]);
  

  //---------- Initiate MC
  if (fIsMC) {
    fMCEvent = NULL;
    fEventCounter->Fill(8);
    if (fIsAOD) {
      fArrayMC = NULL;
      fArrayMC = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!"); 
    } else {
      fMCEvent = MCEvent();
      if (!fMCEvent) {
	Printf("ERROR: Could not retrieve MC event");
	return;
      }
      fMCStack = fMCEvent->Stack();
      if (!fMCStack) {
	Printf("ERROR: Could not retrieve MC stack");
	return;
      }
    }
  }
    //----------

  fEventCounter->Fill(3);
  fNTracks  = event->GetNumberOfTracks();  
  
  Int_t iTracks = 0; 
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(event->GetTrack(idxTrack)); 
    if(!AcceptTrackL(track)) continue;
    
   
      
    Int_t bit = -1;
    Float_t dca[2], cov[3], ndf; // 
     if (track->InheritsFrom("AliESDtrack")) {
      (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
    } else if (track->InheritsFrom("AliAODTrack")) {
       AliAODTrack* clone = dynamic_cast<AliAODTrack*>(track->Clone("trk_clone"));
       if(clone->TestBit(AliAODTrack::kIsDCA)){
       	 dca[0] = clone->DCA();
	 dca[1] = clone->ZAtDCA();
       } else {
	 Double_t dcaa[2], cova[3];
	 Double_t gBzKg = dynamic_cast<AliAODEvent*>(InputEvent())->GetMagneticField();
	 Bool_t propagate = clone->PropagateToDCA(vertex,gBzKg,1000000.,dcaa,cova);
	 if (!propagate) { dca[0] = -999; dca[1] = -999;}
	 dca[0] = Float_t(dcaa[0]);
	 dca[1] = Float_t(dcaa[1]);
       }
       ndf = clone->Chi2perNDF();
           
       //	      Double_t v[3], pos[3];
       //      vertex->GetXYZ(v);
       //     track->GetXYZ(pos);
       //     
       //     Double_t DCAX = pos[0] - v[0];
       //     Double_t DCAY = pos[1] - v[1];
       //     Double_t DCAZ = pos[2] - v[2];
       //     
       //     Double_t DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));
       //     Printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f", 
       //     dca[0], dca[1], DCAXY,  DCAZ, clone->DCA(), clone->ZAtDCA());
       // bit = (Int_t)clone->GetFilterMap();   
       
       // if (clone->IsTPCOnly() && clone->IsHybridTPCConstrainedGlobal()) bit = 2;
       // else if (clone->IsHybridTPCConstrainedGlobal()) bit = 3;
       // else if (clone->IsTPCOnly()) bit = 1;
       // else 

	 bit = (Int_t)clone->GetFilterMap();    

       delete clone;  
     }
     
     // Printf("%f %f %f",dca[0], dca[1], ndf);
     //  if ( TMath::Abs(dca[0]) > fDcaXy ) continue;
     // if  ( TMath::Abs(dca[1]) > fDcaZ )  continue;

    Int_t icharge = track->Charge() < 0 ? -1 : 1;
    Int_t a = fHelperPID->GetParticleSpecies(track,kTRUE);

    Int_t b = -999;
    if (a == 0 )      b = 1;
    else if (a == 1 ) b = 2;
    else if (a == 2 ) b = 3;
    else              b = 4;    

    Float_t lPt = (Float_t)track->Pt();

    if (icharge < 0) {
      fHistPt00->Fill(fCentrality[0],lPt); 
      if (b == 1) fHistPt10->Fill(fCentrality[0],lPt); 
      else if (b == 2)fHistPt20->Fill(fCentrality[0],lPt); 
      else if (b == 3)fHistPt30->Fill(fCentrality[0],lPt); 
    } else if (icharge>0) {
      fHistPt01->Fill(fCentrality[0],lPt); 
      if (b == 1) fHistPt11->Fill(fCentrality[0],lPt); 
      else if (b == 2)fHistPt21->Fill(fCentrality[0],lPt); 
      else if (b == 3)fHistPt31->Fill(fCentrality[0],lPt); 
    }

    //    cout << bit <<endl;

    fTrackBit[iTracks] = bit;

    fTrackCnDf[iTracks]   =  ndf;
    fTrackTpcNcl[iTracks] = track->GetTPCClusterInfo(2,1);

    fTrackPt[iTracks]  = (Float_t)track->Pt();
    fTrackPhi[iTracks] = (Float_t)track->Phi();
    fTrackEta[iTracks] = (Float_t)track->Eta();
    fTrackDxy[iTracks] = dca[0];
    fTrackDz[iTracks]  = dca[1];
    fTrackPid[iTracks] = icharge*b;

    //----------------------------------------
    if (fIsMC) {
      Int_t label  = TMath::Abs(track->GetLabel()); 
      fTrackLabel[iTracks] = label;

      Bool_t isPhysicalPrimary        = 0;
      Bool_t isSecondaryFromWeakDecay = 0;
      Bool_t isSecondaryFromMaterial  = 0;
      AliVParticle* particle = NULL;
      if (track->InheritsFrom("AliESDtrack")) {
	particle = static_cast<AliVParticle*>(fMCEvent->GetTrack(label));
	if (!particle) return;
	isPhysicalPrimary        = fMCStack->IsPhysicalPrimary(label);
	isSecondaryFromWeakDecay = fMCStack->IsSecondaryFromWeakDecay(label);
	isSecondaryFromMaterial  = fMCStack->IsSecondaryFromMaterial(label);
      } else {
	particle                 =  static_cast<AliVParticle*>(fArrayMC->At(label));
	isPhysicalPrimary        =  (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
	isSecondaryFromWeakDecay =  (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
	isSecondaryFromMaterial  =  (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
      }

      fPidStat[iTracks] = -1;
      if (isPhysicalPrimary) {
	if (b != 4) {
	  if (particle->PdgCode() == (track->Charge()*GetPDG(b))) fPidStat[iTracks] = 1;
	  else  fPidStat[iTracks]  = 2;
	}
	else  fPidStat[iTracks] = 3;
      } else if(isSecondaryFromWeakDecay)fPidStat[iTracks] = 5;
      else if (isSecondaryFromMaterial) fPidStat[iTracks] = 6;
      else fPidStat[iTracks] = 7;
    }
    //========================
    /*Printf("%6d %10.5f %10.5f %10.5f %10.5f %10.5f %2d %10.6f %4d %3d %10d",iTracks, 
	   fTrackPt[iTracks],fTrackPhi[iTracks],fTrackEta[iTracks],
	   fTrackDxy[iTracks],fTrackDz[iTracks],fTrackPid[iTracks], 
	   fTrackCnDf[iTracks], fTrackTpcNcl[iTracks], fPidStat[iTracks], fTrackLabel[iTracks]);
    */
    iTracks++;
  }
  fNumberOfTracks = iTracks;
  fEventCounter->Fill(7);
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    Int_t mTracks = 0;
    fEventCounter->Fill(8);
    if (fIsAOD) {
      // fArrayMC = NULL;
      // fArrayMC = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
      // if (!fArrayMC)
      //    AliFatal("No array of MC particles found !!!"); 
      
      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;

	if (!particle->IsPhysicalPrimary()) continue;
	if (!AcceptTrackLMC((AliVParticle*)particle)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? -1 : 1;
	Int_t iPid    = -999;  

	if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	else    iPid = 4;

	fTrackPtM[mTracks]    = (Float_t)particle->Pt();
	fTrackPhiM[mTracks]   = (Float_t)particle->Phi();
	fTrackEtaM[mTracks]   = (Float_t)particle->Eta();
	fTrackPidM[mTracks]   = icharge*iPid;
      	fTrackLabelM[mTracks] = idxMC;
	/*
	  Printf("%6d  %10.5f  %10.5f  %10.5f  %4d  %10d",mTracks, 
	  fTrackPtM[mTracks],fTrackPhiM[mTracks],fTrackEtaM[mTracks],fTrackPidM[mTracks], 
	  fTrackLabelM[mTracks]);
	
	*/
	mTracks++;
      }
      fEventCounter->Fill(9);
      fNumberOfTracksM = mTracks;
    } else  {

      //  AliMCEvent* mcEvent = MCEvent();
      // if (!mcEvent) {
      //	Printf("ERROR: Could not retrieve MC event");
      //	return;
      // }
      // AliStack* stack = mcEvent->Stack();
      // if (!stack) {
      //	Printf("ERROR: Could not retrieve MC stack");
      //	return;
      // }

      fEventCounter->Fill(10);
      for (Int_t idxMC = 0; idxMC < fMCStack->GetNprimary(); ++idxMC) {
	AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	if (!particle) 
	  continue;
	if(!fMCStack->IsPhysicalPrimary(idxMC))  continue;
	
	if (!AcceptTrackLMC(particle)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? -1 : 1;
	
	Int_t iPid = -999;  
	if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	else  iPid = 4;
	
	fTrackPtM[mTracks]  = (Float_t)particle->Pt();
	fTrackPhiM[mTracks] = (Float_t)particle->Phi();
	fTrackEtaM[mTracks] = (Float_t)particle->Eta();
	fTrackPidM[mTracks] = icharge*iPid;
	   
	fTrackLabelM[mTracks] = idxMC;
	mTracks++;
      }
      fEventCounter->Fill(11);
      fNumberOfTracksM = mTracks;
    }
  }
  fEventCounter->Fill(12);
  
  fEventCounter->Fill(5);
  fPidCont->Fill();  
  PostData(1, fThnList); 
  PostData(2, fPidCont);
}

//___________________________________________________________
Bool_t AliEbyEPidTTask::AcceptTrackL(AliVTrack *track) const {
  if (!track) 
    return kFALSE; 
  if (track->Charge() == 0) 
    return kFALSE; 
    
  if (fIsAOD) {  // AOD
    AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
    if (!trackAOD) {
      AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
      return kFALSE; 
    }
    if (!trackAOD->TestFilterBit(fAODtrackCutBit))
      return kFALSE;
  } else {      // ESDs
    if(!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))  return kFALSE;
  }
  
  if(track->Pt() < fPtMin || track->Pt() > fPtMax )  return kFALSE; 
  if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE; 
  return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyEPidTTask::AcceptTrackLMC(AliVParticle *particle) const {
  if(!particle) return kFALSE;
  if (particle->Charge() == 0.0) return kFALSE;
  if (particle->Pt() < fPtMin || particle->Pt() > fPtMax) return kFALSE;
  if (TMath::Abs(particle->Eta()) > fEtaMax) return kFALSE;

  return kTRUE;
}

//___________________________________________________________
void AliEbyEPidTTask::Terminate( Option_t * ){
  Info("AliEbyEPidTTask"," Task Successfully finished");
}
