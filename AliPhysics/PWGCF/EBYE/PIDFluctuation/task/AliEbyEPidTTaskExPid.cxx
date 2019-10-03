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
//         AliEbyE Analysis for producing a tree with complete PID info    //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
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
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"

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

#include "AliEbyEPidTTaskExPid.h"

ClassImp(AliEbyEPidTTaskExPid)

//-----------------------------------------------------------------------
AliEbyEPidTTaskExPid::AliEbyEPidTTaskExPid( const char *name ) : 
AliAnalysisTaskSE( name ), 
  fArrayMC(NULL),
  fESDtrackCuts(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fPIDResponse(0x0),
  fAODtrackCutBit(128),
  fPidCont(0x0),

  fPtMin(0.2),   
  fPtMax(4.), 
  fEtaMin(-1.5), 
  fEtaMax(1.5),  
  
  fIsMC(kFALSE),
  fIsAOD(kFALSE),
  
  fRunNumber(0),
  fNumberOfTracks(15000),
  fNumberOfTracksM(15000),
  fNTracks(0) {
   
  DefineOutput(1, TTree::Class()); //! Connect Outpput....
}

AliEbyEPidTTaskExPid::~AliEbyEPidTTaskExPid() {
  if (fPIDResponse) delete fPIDResponse;
  if (fPidCont)   delete fPidCont;
}

const Int_t AliEbyEPidTTaskExPid::fgkPidEx[] = {0,2,3,4};

//---------------------------------------------------------------------------------
void AliEbyEPidTTaskExPid::UserCreateOutputObjects() {
  
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
  
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fPidCont = new TTree("Event","fPidCont B");
  owd->cd();
  
  fPidCont->Branch("fRunNumber", &fRunNumber,  "fRunNumber/I");
  fPidCont->Branch("fCentrality",fCentrality,"fCentrality[6]/F");
  fPidCont->Branch("Vertex",fVtx,"fVtx[3]/F");
  fPidCont->Branch("fNumberOfTracks", &fNumberOfTracks,"fNumberOfTracks/I");
  fPidCont->Branch("fPParam",fPParam,"fPParam[fNumberOfTracks][15]/F");
  fPidCont->Branch("fTParam",fTParam,"fTParam[fNumberOfTracks][11]/F");
  
  if(fIsMC) {
    Printf(" >>>>>>>>>>>>>>>>>> Update container for MC Events >>>>>>>>>>>>>>>>>>");

    fPidCont->Branch("fTrackLabel",fTrackLabel,"fTrackLabel[fNumberOfTracks]/I");
    fPidCont->Branch("fPidStat",fPidStat,"fPidStat[fNumberOfTracks]/I");
    
    fPidCont->Branch("fNumberOfTracksM", &fNumberOfTracksM,"fNumberOfTracksM/I");
    fPidCont->Branch("fTrackPtM",fTrackPtM,"fTrackPtM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackLabelM",fTrackLabelM,"fTrackLabelM[fNumberOfTracksM]/I");
    
    fPidCont->Branch("fTrackPhiM",fTrackPhiM,"fTrackPhiM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackEtaM",fTrackEtaM,"fTrackEtaM[fNumberOfTracksM]/F");
    fPidCont->Branch("fTrackPidM",fTrackPidM,"fTrackPidM[fNumberOfTracksM]/I");
  }
  
  PostData(1, fPidCont);  
}

//----------------------------------------------------------------------------------
void AliEbyEPidTTaskExPid::UserExec( Option_t * ){
  AliVEvent *event = InputEvent();
  if (!event) return;
  AliInputEventHandler* fInputEventHandler = static_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fInputEventHandler) return;
  
  fPIDResponse = fInputEventHandler->GetPIDResponse();
  if(!fPIDResponse) return;

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
  
 
  fVtx[0] = vertex->GetX();
  fVtx[1] = vertex->GetY();
  fVtx[2] = vertex->GetZ();
  
  AliCentrality *centrality = event->GetCentrality();
  
  if (centrality->GetQuality() != 0) return;
  
  fRunNumber = event->GetRunNumber();
  if (!(fInputEventHandler->IsEventSelected() & AliVEvent::kMB))  return;        
  
  // AliAODHeader *header = (AliAODHeader*) event->GetHeader();
  // if (!header) return;

  fCentrality[0] = centrality->GetCentralityPercentile("V0M");
  fCentrality[1] = centrality->GetCentralityPercentile("CL1");
  fCentrality[2] = centrality->GetCentralityPercentile("TRK");
  
  fCentrality[3] = event->GetVZEROData()->GetMTotV0A();
  fCentrality[4] = event->GetVZEROData()->GetMTotV0C();
  fCentrality[5] = event->GetZDCData()->GetZDCParticipants();
  
  //---------- Initiate MC
  if (fIsMC) {
    fMCEvent = NULL;
   
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
     /* Part from Hans
   // Bigger loop has bad precision, we're nearly one centimeter too far, 
     // go back in small steps.
     while (shiftedRadiusSquared>RSquaredWanted){
	// Propagate a mm inwards
	x-=.1;
	if(!etp.PropagateTo(x,bfield)){
	  // Propagation failed but we're already with a
	  // cm precision at R=1.25m so we only break the 
	  // inner loop
	  break;
	}
	// Get the global position
	etp.GetXYZ(xyz);
	// Calculate shifted radius, squared
	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
			     }
     // We reached R=1.25m with a precission of a cm to a mm,
     // set the spatial position
     fXSftR125[0]=xyz[0]-priVtx[0];
     fXSftR125[1]=xyz[1]-priVtx[1];
     fXSftR125[2]=xyz[2]-priVtx[2];
     // Done
     return;
     } // End of if roughly reached radius
   } // End of coarse propagation loop
     */
     fPParam[iTracks][0] = (Float_t)track->GetITSsignal();
     fPParam[iTracks][1] = (Float_t)track->GetTPCsignal();
     fPParam[iTracks][2] = (Float_t)track->GetTOFsignal();
     for (Int_t i = 0; i < 4; i++) {
       fPParam[iTracks][i+3]  = (Float_t)fPIDResponse->NumberOfSigmasITS(track, (AliPID::EParticleType)fgkPidEx[i]);
       fPParam[iTracks][i+7]  = (Float_t)fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)fgkPidEx[i]);
       fPParam[iTracks][i+11] = (Float_t)fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)fgkPidEx[i]);
     }
     
     fTParam[iTracks][0] = (Float_t)track->Charge();
     fTParam[iTracks][1] = (Float_t)track->GetStatus();
     fTParam[iTracks][2] = (Float_t)track->GetTPCNclsF();
     fTParam[iTracks][3] = (Float_t)track->P();
     fTParam[iTracks][4] = (Float_t)track->Pt();
     fTParam[iTracks][5] = (Float_t)ndf;
     fTParam[iTracks][6] = (Float_t)track->GetTPCClusterInfo(2,1);
     fTParam[iTracks][7] = (Float_t)track->Phi();
     fTParam[iTracks][8] = (Float_t)track->Eta();
     fTParam[iTracks][9] = dca[0];
     fTParam[iTracks][10] = dca[1];
    
     /*
     Printf(" %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f", 
	    fPParam[iTracks][0], fPParam[iTracks][1], fPParam[iTracks][2], 
	    fPParam[iTracks][3], fPParam[iTracks][4], fPParam[iTracks][5], 
	    fPParam[iTracks][6], fPParam[iTracks][7], fPParam[iTracks][8],  
	    fPParam[iTracks][9], fPParam[iTracks][10], fPParam[iTracks][11],
	    fPParam[iTracks][12], fPParam[iTracks][13], fPParam[iTracks][14], fTParam[iTracks][0]);
     */
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
      if (isPhysicalPrimary) fPidStat[iTracks] = 0;
      else if   (isSecondaryFromWeakDecay)           fPidStat[iTracks] = 1; 
      else if   (isSecondaryFromMaterial)            fPidStat[iTracks] = 2; 
      else {
	if (TMath::Abs(particle->PdgCode()) ==   11)      fPidStat[iTracks] = 3; 
	else if (TMath::Abs(particle->PdgCode()) ==   13) fPidStat[iTracks] = 4; 
	else if (TMath::Abs(particle->PdgCode()) ==  211) fPidStat[iTracks] = 5; 
	else if (TMath::Abs(particle->PdgCode()) ==  321) fPidStat[iTracks] = 6; 
	else if (TMath::Abs(particle->PdgCode()) == 2212) fPidStat[iTracks] = 7; 
	else                                              fPidStat[iTracks] = 8; 
      }
    }
    iTracks++;
  }
  fNumberOfTracks = iTracks;
 
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    Int_t mTracks = 0;
    if (fIsAOD) {
      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;
	
	if (!particle->IsPhysicalPrimary()) continue;
	if (!AcceptTrackLMC((AliVParticle*)particle)) continue;
	
	fTrackPtM[mTracks]    = (Float_t)particle->Pt();
	fTrackPhiM[mTracks]   = (Float_t)particle->Phi();
	fTrackEtaM[mTracks]   = (Float_t)particle->Eta();
	//fTrackPidM[mTracks]   = particle->PdgCode();
	fTrackLabelM[mTracks] = idxMC;
	mTracks++;
      }
      fNumberOfTracksM = mTracks;
    } else  {
      for (Int_t idxMC = 0; idxMC < fMCStack->GetNprimary(); ++idxMC) {
		AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	if (!particle) 
	  continue;
	if(!fMCStack->IsPhysicalPrimary(idxMC))  continue;
	
	if (!AcceptTrackLMC(particle)) continue;
	fTrackPtM[mTracks]  = (Float_t)particle->Pt();
	fTrackPhiM[mTracks] = (Float_t)particle->Phi();
	fTrackEtaM[mTracks] = (Float_t)particle->Eta();
	//	fTrackPidM[mTracks] = particle->PdgCode();
	
	fTrackLabelM[mTracks] = idxMC;
	mTracks++;
      }
      
      fNumberOfTracksM = mTracks;
    }
  }
  fPidCont->Fill();  
  PostData(1, fPidCont);
}

//___________________________________________________________
Bool_t AliEbyEPidTTaskExPid::AcceptTrackL(AliVTrack *track) const {
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
Bool_t AliEbyEPidTTaskExPid::AcceptTrackLMC(AliVParticle *particle) const {
  if(!particle) return kFALSE;
  if (particle->Charge() == 0.0) return kFALSE;
  if (particle->Pt() < fPtMin || particle->Pt() > fPtMax) return kFALSE;
  if (TMath::Abs(particle->Eta()) > fEtaMax) return kFALSE;

  return kTRUE;
}

//___________________________________________________________
void AliEbyEPidTTaskExPid::Terminate( Option_t * ){
  Info("AliEbyEPidTTaskExPid"," Task Successfully finished");
}
