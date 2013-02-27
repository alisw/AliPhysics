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

// Class for heavy-flavour electron  with EMCal triggered events
// Author: Shingo Sakai


#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

#include "AliAnalysisTaskHFECal.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliEMCALTrack.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliStack.h"

#include "AliCentrality.h"

ClassImp(AliAnalysisTaskHFECal)
//________________________________________________________________________
AliAnalysisTaskHFECal::AliAnalysisTaskHFECal(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fMC(0)
  ,stack(0)
  ,fGeom(0)
  ,fOutputList(0)
  ,fqahist(1) 
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fmcData(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.01)	
  ,fNoEvents(0)
  ,fEMCAccE(0)
  ,hEMCAccE(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fIncpT(0)	
  ,fIncpTM20(0)	
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fInvmassLSmc(0)		
  ,fInvmassULSmc(0)		
  ,fInvmassLSmc0(0)		
  ,fInvmassLSmc1(0)		
  ,fInvmassLSmc2(0)		
  ,fInvmassLSmc3(0)		
  ,fInvmassULSmc0(0)		
  ,fInvmassULSmc1(0)		
  ,fInvmassULSmc2(0)		
  ,fInvmassULSmc3(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fPhotoElecPt(0)
  ,fPhoElecPt(0)
  ,fPhoElecPtM20(0)
  ,fSameElecPt(0)
  ,fSameElecPtM20(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fEleInfo(0)
  /*
  ,fClsEBftTrigCut(0)	 
  ,fClsEAftTrigCut(0)	 
  ,fClsEAftTrigCut1(0)	 
  ,fClsEAftTrigCut2(0)	 
  ,fClsEAftTrigCut3(0)	 
  ,fClsEAftTrigCut4(0)	 
  ,fClsETime(0)       
  ,fClsETime1(0)       
  ,fTrigTimes(0)
  ,fCellCheck(0)
  */
  ,fInputHFEMC(0)
  ,fInputAlle(0)
  ,fIncpTMChfe(0)	
  ,fIncpTMChfeAll(0)	
  ,fIncpTMCM20hfe(0)	
  ,fIncpTMCM20hfeAll(0)	
  ,fIncpTMCM20hfeCheck(0)	
  ,fInputHFEMC_weight(0)
  ,fIncpTMCM20hfeCheck_weight(0)
  ,fIncpTMCpho(0)	
  ,fIncpTMCM20pho(0)	
  ,fPhoElecPtMC(0)
  ,fPhoElecPtMCM20(0)
  ,fSameElecPtMC(0)
  ,fSameElecPtMCM20(0)
  ,fIncpTMCM20pho_pi0e(0)	
  ,fPhoElecPtMCM20_pi0e(0)
  ,fSameElecPtMCM20_pi0e(0)
  ,fIncpTMCM20pho_eta(0)	
  ,fPhoElecPtMCM20_eta(0)
  ,fSameElecPtMCM20_eta(0)
  ,fIncpTMCpho_pi0e_TPC(0)	
  ,fPhoElecPtMC_pi0e_TPC(0)
  ,fSameElecPtMC_pi0e_TPC(0)
  ,CheckNclust(0)
  ,CheckNits(0)
  ,Hpi0pTcheck(0)
  ,HETApTcheck(0)
  ,HphopTcheck(0)
  ,HDpTcheck(0)
  ,HBpTcheck(0)
  ,fpTCheck(0)
  ,fMomDtoE(0) 
  ,fLabelCheck(0)
  ,fgeoFake(0)
{
  //Named constructor
  
  fPID = new AliHFEpid("hfePid");
  fTrackCuts = new AliESDtrackCuts();
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHFECal::AliAnalysisTaskHFECal() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskHFECal")
  ,fESD(0)
  ,fMC(0)
  ,stack(0)
  ,fGeom(0)
  ,fOutputList(0)
  ,fqahist(1)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
  ,fmcData(kFALSE)
  ,fVz(0.0)
  ,fCFM(0)	
  ,fPID(0)       
  ,fPIDqa(0)	       
  ,fOpeningAngleCut(0.1)
  ,fInvmassCut(0.01)	
  ,fNoEvents(0)
  ,fEMCAccE(0)
  ,hEMCAccE(0)
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fIncpT(0)	 
  ,fIncpTM20(0)	 
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fInvmassLSmc(0)		
  ,fInvmassULSmc(0)		
  ,fInvmassLSmc0(0)		
  ,fInvmassLSmc1(0)		
  ,fInvmassLSmc2(0)		
  ,fInvmassLSmc3(0)		
  ,fInvmassULSmc0(0)		
  ,fInvmassULSmc1(0)		
  ,fInvmassULSmc2(0)		
  ,fInvmassULSmc3(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fPhotoElecPt(0)
  ,fPhoElecPt(0)
  ,fPhoElecPtM20(0)
  ,fSameElecPt(0)
  ,fSameElecPtM20(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fEleInfo(0)
  /*
  ,fClsEBftTrigCut(0)	 
  ,fClsEAftTrigCut(0)	 
  ,fClsEAftTrigCut1(0)	 
  ,fClsEAftTrigCut2(0)	 
  ,fClsEAftTrigCut3(0)	 
  ,fClsEAftTrigCut4(0)	 
  ,fClsETime(0)       
  ,fClsETime1(0)       
  ,fTrigTimes(0)
  ,fCellCheck(0)
  */
  ,fInputHFEMC(0)
  ,fInputAlle(0)
  ,fIncpTMChfe(0)	
  ,fIncpTMChfeAll(0)	
  ,fIncpTMCM20hfe(0)	
  ,fIncpTMCM20hfeAll(0)	
  ,fIncpTMCM20hfeCheck(0)	
  ,fInputHFEMC_weight(0)
  ,fIncpTMCM20hfeCheck_weight(0)
  ,fIncpTMCpho(0)	
  ,fIncpTMCM20pho(0)	
  ,fPhoElecPtMC(0)
  ,fPhoElecPtMCM20(0)
  ,fSameElecPtMC(0)
  ,fSameElecPtMCM20(0)
  ,fIncpTMCM20pho_pi0e(0)	
  ,fPhoElecPtMCM20_pi0e(0)
  ,fSameElecPtMCM20_pi0e(0)
  ,fIncpTMCM20pho_eta(0)	
  ,fPhoElecPtMCM20_eta(0)
  ,fSameElecPtMCM20_eta(0)
  ,fIncpTMCpho_pi0e_TPC(0)	
  ,fPhoElecPtMC_pi0e_TPC(0)
  ,fSameElecPtMC_pi0e_TPC(0)
  ,CheckNclust(0)
  ,CheckNits(0)
  ,Hpi0pTcheck(0)
  ,HETApTcheck(0)
  ,HphopTcheck(0)
  ,HDpTcheck(0)
  ,HBpTcheck(0)
  ,fpTCheck(0)
  ,fMomDtoE(0)
  ,fLabelCheck(0)
  ,fgeoFake(0)
{
	//Default constructor
	fPID = new AliHFEpid("hfePid");

	fTrackCuts = new AliESDtrackCuts();
	
	// Constructor
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 id reserved by the base class for AOD
	// Output slot #1 writes into a TH1 container
	// DefineOutput(1, TH1I::Class());
	DefineOutput(1, TList::Class());
	//DefineOutput(3, TTree::Class());
}
//_________________________________________

AliAnalysisTaskHFECal::~AliAnalysisTaskHFECal()
{
  //Destructor 
  
  delete fOutputList;
  delete fGeom;
  delete fPID;
  delete fCuts;
  delete fCFM;
  delete fPIDqa;
  delete fTrackCuts;
}
//_________________________________________

void AliAnalysisTaskHFECal::UserExec(Option_t*)
{
  //Main loop
  //Called for each event
  
  // create pointer to event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }
  
  if(!fPID->IsInitialized()){ 
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
    fPID->InitializePID(fESD->GetRunNumber());
  }
 
  if(fmcData)fMC = MCEvent();
  //AliStack* stack = NULL;
  if(fmcData && fMC)stack = fMC->Stack();

  Float_t cent = -1.;
  AliCentrality *centrality = fESD->GetCentrality(); 
  cent = centrality->GetCentralityPercentile("V0M");

  //---- fill MC track info
  if(fmcData && fMC)
    {
    Int_t nParticles = stack->GetNtrack();
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = stack->Particle(iParticle);
      int fPDG = particle->GetPdgCode(); 
      double mcZvertex = fMC->GetPrimaryVertex()->GetZ();
      double pTMC = particle->Pt();
      double proR = particle->R();
      double etaMC = particle->Eta();
      if(fabs(etaMC)>0.6)continue;
      Bool_t mcInDtoE= kFALSE;
      Bool_t mcInBtoE= kFALSE;

      Bool_t MChijing = fMC->IsFromBGEvent(iParticle);
      //if(!MChijing)printf("not MC hijing");
      int iHijing = 1;
      if(!MChijing)iHijing = 0;
      double mcphoinfo[5];
      mcphoinfo[0] = cent;
      mcphoinfo[1] = pTMC;
      mcphoinfo[2] = iHijing;
      //if(fPDG==111)Hpi0pTcheck->Fill(pTMC,iHijing);
      //if(fPDG==221)HETApTcheck->Fill(pTMC,iHijing);
      if(fPDG==111)Hpi0pTcheck->Fill(mcphoinfo);
      if(fPDG==221)HETApTcheck->Fill(mcphoinfo);
      if(fabs(fPDG)==411 || fabs(fPDG)==413 || fabs(fPDG)==421 || fabs(fPDG)==423 || fabs(fPDG)==431)HDpTcheck->Fill(pTMC,iHijing);
      if(fabs(fPDG)==511 || fabs(fPDG)==513 || fabs(fPDG)==521 || fabs(fPDG)==523 || fabs(fPDG)==531)HBpTcheck->Fill(pTMC,iHijing);

      if(particle->GetFirstMother()>-1 && fPDG==22)
        {
	 int parentPID = stack->Particle(particle->GetFirstMother())->GetPdgCode();  
         if(parentPID==111 || parentPID==221)HphopTcheck->Fill(pTMC,iHijing); // pi0->g & eta->g
        } 

      if(particle->GetFirstMother()>-1 && fabs(fPDG)==11)
        {
	    int parentPID = stack->Particle(particle->GetFirstMother())->GetPdgCode();  
	    double pTMCparent = stack->Particle(particle->GetFirstMother())->Pt();  
            if((fabs(parentPID)==411 || fabs(parentPID)==413 || fabs(parentPID)==421 || fabs(parentPID)==423 || fabs(parentPID)==431)&& fabs(fPDG)==11)mcInDtoE = kTRUE;
            if((fabs(parentPID)==511 || fabs(parentPID)==513 || fabs(parentPID)==521 || fabs(parentPID)==523 || fabs(parentPID)==531)&& fabs(fPDG)==11)mcInBtoE = kTRUE;
            if((mcInBtoE || mcInDtoE) && fabs(mcZvertex)<10.0)
               {
                fInputHFEMC->Fill(cent,pTMC);
                double mcinfo[5];
                mcinfo[0] = cent;
                mcinfo[1] = pTMC;
                mcinfo[2] = 0.0;
                mcinfo[3] = iHijing;
                mcinfo[4] = pTMCparent;
                fInputHFEMC_weight->Fill(mcinfo);
               }
         }


         if(proR<7.0 && fabs(fPDG)==11)fInputAlle->Fill(cent,pTMC);

      }
    } 

  fNoEvents->Fill(0);

  Int_t fNOtrks =  fESD->GetNumberOfTracks();
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();
  
  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(1);
  
  if(fNOtrks<1) return;
  fNoEvents->Fill(2);
  
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
  }
  
  fPID->SetPIDResponse(pidResponse);
  
  fCFM->SetRecEventInfo(fESD);
  
  //Float_t cent = -1.;
  //AliCentrality *centrality = fESD->GetCentrality(); 
  //cent = centrality->GetCentralityPercentile("V0M");
  fCent->Fill(cent);
  
  //if(cent>90.) return;
	
 // Calorimeter info.
 
   FindTriggerClusters();

  // make EMCAL array 
  for(Int_t iCluster=0; iCluster<fESD->GetNumberOfCaloClusters(); iCluster++)
     {
      AliESDCaloCluster *clust = fESD->GetCaloCluster(iCluster);
      if(clust && clust->IsEMCAL())
        {
         double clustE = clust->E();
         float  emcx[3]; // cluster pos
         clust->GetPosition(emcx);
         TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
         double emcphi = clustpos.Phi(); 
         double emceta = clustpos.Eta();
         double calInfo[5];
         calInfo[0] = emcphi; calInfo[1] = emceta; calInfo[2] = clustE; calInfo[3] = cent; calInfo[4] = clust->Chi2(); 
         //fEMCAccE->Fill(calInfo); 
         //if(clustE>3.0)fEMCAccE->Fill(calInfo); 
         //if(fqahist==1 && clustE>1.5)fEMCAccE->Fill(calInfo); 
         hEMCAccE->Fill(cent,clustE); 
        }
   }

  // Track loop 
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
   
    int parentlabel = 99999;
    int parentPID = 99999;
    int grand_parentlabel = 99999;
    int grand_parentPID = 99999;
    Bool_t mcPho = kFALSE;
    Bool_t mcDtoE= kFALSE;
    Bool_t mcBtoE= kFALSE;
    Bool_t mcOrgPi0 = kFALSE;
    Bool_t mcOrgEta = kFALSE;
    double mcele = -1.;
    double mcpT = 0.0;
    double mcMompT = 0.0;
    //double mcGrandMompT = 0.0;
    double mcWeight = -10.0;

    int iHijing = 1;
    int mcLabel = -1;

    if(fmcData && fMC && stack)
      {
       Int_t label = TMath::Abs(track->GetLabel());
       mcLabel = track->GetLabel();
       
       //if(mcLabel>-1)
       {
      
	       Bool_t MChijing = fMC->IsFromBGEvent(label);
	       if(!MChijing)iHijing = 0;

	       TParticle* particle = stack->Particle(label);
	       int mcpid = particle->GetPdgCode();
	       mcpT = particle->Pt();
	       //printf("MCpid = %d",mcpid);
	       if(particle->GetFirstMother()>-1)
	       {
		       //int parentlabel = particle->GetFirstMother();
		       //int parentPID = stack->Particle(particle->GetFirstMother())->GetPdgCode();
		       //mcMompT = stack->Particle(particle->GetFirstMother())->Pt();
		       FindMother(particle, parentlabel, parentPID);
		       mcMompT = stack->Particle(parentlabel)->Pt();
		       if((parentPID==22 || parentPID==111 || parentPID==221)&& fabs(mcpid)==11)mcPho = kTRUE;
		       if((fabs(parentPID)==411 || fabs(parentPID)==413 || fabs(parentPID)==421 || fabs(parentPID)==423 || fabs(parentPID)==431)&& fabs(mcpid)==11)mcDtoE = kTRUE;
		       if((fabs(parentPID)==511 || fabs(parentPID)==513 || fabs(parentPID)==521 || fabs(parentPID)==523 || fabs(parentPID)==531)&& fabs(mcpid)==11)mcBtoE = kTRUE;

		       // make D->e pT correlation
		       if(mcDtoE)fMomDtoE->Fill(mcpT,mcMompT); 

		       //cout << "check PID = " << parentPID << endl;
		       //cout << "check pho = " << mcPho << endl;
		       //cout << "check D or B = " << mcDtoE << endl;
		       // pi->e (Dalitz)
		       if(parentPID==111 && fabs(mcpid)==11 && mcMompT>0.0)
		       {
			       //cout << "find pi0->e " <<  endl;
			       mcOrgPi0 = kTRUE;
			       mcWeight = GetMCweight(mcMompT); 
		       }
		       // eta->e (Dalitz)
		       if(parentPID==221 && fabs(mcpid)==11 && mcMompT>0.0)
		       {
			       //cout << "find Eta->e " <<  endl;
			       mcOrgEta = kTRUE;
			       mcWeight = GetMCweightEta(mcMompT); 
		       }

		       // access grand parent 
		       TParticle* particle_parent = stack->Particle(parentlabel); // get parent pointer
		       //if(particle_parent->GetFirstMother()>-1 && parentPID==22 && fabs(mcpid)==11) // get grand parent g->e
		       if(particle_parent->GetFirstMother()>-1 && (parentPID==22 || parentPID==111) && fabs(mcpid)==11) // get grand parent g->e
		       {
			       //int grand_parentPID = stack->Particle(particle_parent->GetFirstMother())->GetPdgCode();
			       //double pTtmp = stack->Particle(particle_parent->GetFirstMother())->Pt();
			       FindMother(particle_parent, grand_parentlabel, grand_parentPID);
			       double mcGrandpT = stack->Particle(grand_parentlabel)->Pt();
			       if(grand_parentPID==111 && mcGrandpT>0.0)
			       {
				       // check eta->pi0 decay !
				       int grand2_parentlabel = 99999; int grand2_parentPID = 99999;
				       TParticle* particle_grand = stack->Particle(grand_parentlabel); // get parent pointer
				       FindMother(particle_grand, grand2_parentlabel, grand2_parentPID);
				       if(grand2_parentPID==221)
				       {
					       //cout << "find Eta->e " <<  endl;
					       double mcGrandpT2 = stack->Particle(grand2_parentlabel)->Pt();
					       mcOrgEta = kTRUE;
					       mcWeight = GetMCweight(mcGrandpT2);  
					       mcMompT = mcGrandpT2; 
				       }
				       else
				       {
					       //cout << "find pi0->e " <<  endl;
					       mcOrgPi0 = kTRUE;
					       mcWeight = GetMCweight(mcGrandpT);  
					       mcMompT = mcGrandpT; 
				       }
			       }

			       if(grand_parentPID==221 && mcGrandpT>0.0)
			       {
				       //cout << "find Eta->e " <<  endl;
				       mcOrgEta = kTRUE;
				       mcOrgPi0 = kFALSE;
				       mcWeight = GetMCweightEta(mcGrandpT); 
				       mcMompT = mcGrandpT; 
			       }
		       }
	       }

	       if(fabs(mcpid)==11 && mcDtoE)mcele= 1.; 
	       if(fabs(mcpid)==11 && mcBtoE)mcele= 2.; 
	       if(fabs(mcpid)==11 && mcPho)mcele= 3.; 

       } // end of mcLabel>-1

      } // end of MC info.
 
    //cout << "Pi0 = " << mcOrgPi0 << " ; Eta = " << mcOrgEta << endl; 
    //printf("weight = %f\n",mcWeight);

    if(TMath::Abs(track->Eta())>0.6) continue;
    if(TMath::Abs(track->Pt()<2.0)) continue;
    
    fTrackPtBefTrkCuts->Fill(track->Pt());		
    // RecKine: ITSTPC cuts  
    if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    
    //RecKink
    if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
      if(track->GetKinkIndex(0) != 0) continue;
    } 
    
    // RecPrim
    if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
    
    // HFEcuts: ITS layers cuts
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
    
    // HFE cuts: TPC PID cleanup
    if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;

    int nTPCcl = track->GetTPCNcls();
    //int nTPCclF = track->GetTPCNclsF(); // warnings
    int nITS = track->GetNcls(0);
    
    fTrackPtAftTrkCuts->Fill(track->Pt());		
    
    Double_t mom = -999., eop=-999., pt = -999., dEdx=-999., fTPCnSigma=-10, phi=-999., eta=-999.;
    pt = track->Pt();
    if(pt<2.0)continue;
    
    // Track extrapolation
    
    Int_t charge = track->Charge();
    fTrkpt->Fill(pt);
    mom = track->P();
    phi = track->Phi();
    eta = track->Eta();
    dEdx = track->GetTPCsignal();
    fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
    if(mcLabel==-1) // nsigma mean correction
      {
       double mean_corr = NsigCorr(cent);
       printf("correction %f\n",mean_corr);
       fTPCnSigma -= mean_corr;
      }    

        double ncells = -1.0;
        double m20 = -1.0;
        double m02 = -1.0;
        double disp = -1.0;
        double rmatch = -1.0;  
        double nmatch = -1.0;
        double oppstatus = 0.0;

    Bool_t fFlagPhotonicElec = kFALSE;
    Bool_t fFlagConvinatElec = kFALSE;

    Int_t clsId = track->GetEMCALcluster();
    if (clsId>0){
      AliESDCaloCluster *clust = fESD->GetCaloCluster(clsId);
      if(clust && clust->IsEMCAL()){

	        double clustE = clust->E();
                eop = clustE/fabs(mom);
                //double clustT = clust->GetTOF();
                ncells = clust->GetNCells();
                m02 = clust->GetM02();
                m20 = clust->GetM20();
                disp = clust->GetDispersion();
		double delphi = clust->GetTrackDx(); 
		double deleta = clust->GetTrackDz(); 
		rmatch = sqrt(pow(delphi,2)+pow(deleta,2));
		nmatch = clust->GetNTracksMatched();

		if(fTPCnSigma>-1.5 && fTPCnSigma<3.0)
		{
		  SelectPhotonicElectron(iTracks,cent,track,fFlagPhotonicElec,fFlagConvinatElec,fTPCnSigma,m20,eop,mcele,mcWeight,iHijing,mcOrgPi0,mcOrgEta);
		}
		if(fFlagPhotonicElec)oppstatus = 1.0;
		if(fFlagConvinatElec)oppstatus = 2.0;
		if(fFlagPhotonicElec && fFlagConvinatElec)oppstatus = 3.0;

		  double valdedx[16];
		  valdedx[0] = pt; valdedx[1] = nITS; valdedx[2] = phi; valdedx[3] = eta; valdedx[4] = fTPCnSigma;
		  //valdedx[5] = eop; valdedx[6] = rmatch; valdedx[7] = ncells,  valdedx[8] = nTPCclF; valdedx[9] = m20; valdedx[10] = mcpT;
		  valdedx[5] = eop; valdedx[6] = rmatch; valdedx[7] = ncells,  valdedx[8] = nmatch; valdedx[9] = m20; valdedx[10] = mcpT;
		  valdedx[11] = cent; valdedx[12] = charge; valdedx[13] = oppstatus; valdedx[14] = nTPCcl;
                  valdedx[15] = mcele;
                  if(fqahist==1)fEleInfo->Fill(valdedx);
                 

      }
    }
     
    if(nITS<2.5)continue;
    if(nTPCcl<100)continue;
   
    CheckNclust->Fill(nTPCcl); 
    CheckNits->Fill(nITS); 

    fdEdxBef->Fill(mom,fTPCnSigma);
    fTPCnsigma->Fill(mom,fTPCnSigma);
    if(fTPCnSigma >= -1.0 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,eop);

    Int_t pidpassed = 1;
 
    // check reco eff. with TPC

    double phoval[5];
    phoval[0] = cent;
    phoval[1] = pt;
    phoval[2] = fTPCnSigma;
    phoval[3] = iHijing;
    phoval[4] = mcMompT;
   
    if((fTPCnSigma >= -1.0 && fTPCnSigma <= 3) && mcele>0 && mcPho && mcOrgPi0)
      {
        if(iHijing==1)mcWeight = 1.0; 
        fIncpTMCpho_pi0e_TPC->Fill(phoval,mcWeight);    
        if(fFlagPhotonicElec) fPhoElecPtMC_pi0e_TPC->Fill(phoval,mcWeight);
        if(fFlagConvinatElec) fSameElecPtMC_pi0e_TPC->Fill(phoval,mcWeight);
      }

    //--- track accepted
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    double binct = 10.5;
    if((0.0< cent) && (cent<5.0)) binct = 0.5;
    if((5.0< cent) && (cent<10.0)) binct = 1.5;
    if((10.0< cent) && (cent<20.0)) binct = 2.5;
    if((20.0< cent) && (cent<30.0)) binct = 3.5;
    if((30.0< cent) && (cent<40.0)) binct = 4.5;
    if((40.0< cent) && (cent<50.0)) binct = 5.5;
    if((50.0< cent) && (cent<60.0)) binct = 6.5;
    if((60.0< cent) && (cent<70.0)) binct = 7.5;
    if((70.0< cent) && (cent<80.0)) binct = 8.5;
    if((80.0< cent) && (cent<90.0)) binct = 9.5;
    if((90.0< cent) && (cent<100.0)) binct = 10.5;

    hfetrack.SetCentrality((int)binct); //added
    hfetrack.SetPbPb();
    if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;

    if(pidpassed==0) continue;
    
    fTrkEovPAft->Fill(pt,eop);
    fdEdxAft->Fill(mom,fTPCnSigma);

    // Fill real data
    fIncpT->Fill(cent,pt);    
    if(fFlagPhotonicElec) fPhoElecPt->Fill(cent,pt);
    if(fFlagConvinatElec) fSameElecPt->Fill(cent,pt);

    if(m20>0.0 && m20<0.3)
      { 
       fIncpTM20->Fill(cent,pt);    
       if(fFlagPhotonicElec) fPhoElecPtM20->Fill(cent,pt);
       if(fFlagConvinatElec) fSameElecPtM20->Fill(cent,pt);
     }
    
    // MC
    // check label for electron candidiates

    int idlabel = 1;
    if(mcLabel==0)idlabel = 0;
    fLabelCheck->Fill(pt,idlabel);
    if(mcLabel==0)fgeoFake->Fill(phi,eta);

    if(mcele>0) // select MC electrons
      {

          fIncpTMChfeAll->Fill(cent,pt);    
          if(m20>0.0 && m20<0.3)fIncpTMCM20hfeAll->Fill(cent,pt);    

       if(mcBtoE || mcDtoE) // select B->e & D->e
         {
          fIncpTMChfe->Fill(cent,pt);    
          //if(m20>0.0 && m20<0.3)fIncpTMCM20hfe->Fill(cent,pt);    
          //if(m20>0.0 && m20<0.3)fIncpTMCM20hfeCheck->Fill(cent,mcpT);    
          if(m20>0.0 && m20<0.3)
            {
                cout << "MC label = " << mcLabel << endl;
                fIncpTMCM20hfe->Fill(cent,pt);    
                fIncpTMCM20hfeCheck->Fill(cent,mcpT);    
                fIncpTMCM20hfeCheck_weight->Fill(phoval);    
            }
         }
      
       if(mcPho) // select photonic electrons
        {

         fIncpTMCpho->Fill(phoval);    
         if(fFlagPhotonicElec) fPhoElecPtMC->Fill(phoval);
         if(fFlagConvinatElec) fSameElecPtMC->Fill(phoval);

         if(m20>0.0 && m20<0.3) 
           {
            fIncpTMCM20pho->Fill(phoval);    
            if(fFlagPhotonicElec) fPhoElecPtMCM20->Fill(phoval);
            if(fFlagConvinatElec) fSameElecPtMCM20->Fill(phoval);
            // pi0->g->e
            if(mcWeight>-1)
              {
               if(iHijing==1)mcWeight = 1.0; 
               if(mcOrgPi0)
                 {
                  fIncpTMCM20pho_pi0e->Fill(phoval,mcWeight);    
                  if(fFlagPhotonicElec) fPhoElecPtMCM20_pi0e->Fill(phoval,mcWeight);
                  if(fFlagConvinatElec) fSameElecPtMCM20_pi0e->Fill(phoval,mcWeight);
                  //fIncpTMCM20pho_pi0e->Fill(phoval);   // v5-04-02-AN & v5-04-06-AN 
                  //if(fFlagPhotonicElec) fPhoElecPtMCM20_pi0e->Fill(phoval);
                  //if(fFlagConvinatElec) fSameElecPtMCM20_pi0e->Fill(phoval);
                 }
               if(mcOrgEta)
                 {
                  fIncpTMCM20pho_eta->Fill(phoval,mcWeight);    
                  if(fFlagPhotonicElec) fPhoElecPtMCM20_eta->Fill(phoval,mcWeight);
                  if(fFlagConvinatElec) fSameElecPtMCM20_eta->Fill(phoval,mcWeight);
                  //fIncpTMCM20pho_eta->Fill(phoval);    
                  //if(fFlagPhotonicElec) fPhoElecPtMCM20_eta->Fill(phoval);
                  //if(fFlagConvinatElec) fSameElecPtMCM20_eta->Fill(phoval);
                 }
               // --- eta
              }
           }
        } 
      } 
 }
 PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskHFECal::UserCreateOutputObjects()
{
  //--- Check MC
 
  //Bool_t mcData = kFALSE;
  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())
    {
     fmcData = kTRUE;
     printf("+++++ MC Data available");
    }
  if(fmcData)
    {
     printf("++++++++= MC analysis \n");
    }
  else
   {
     printf("++++++++ real data analysis \n");
   }

   printf("+++++++ QA hist %d",fqahist);

  //---- Geometry
  fGeom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

  //--------Initialize PID
  //fPID->SetHasMCData(kFALSE);
  fPID->SetHasMCData(fmcData);
  if(!fPID->GetNumberOfPIDdetectors()) 
    {
      fPID->AddDetector("TPC", 0);
      fPID->AddDetector("EMCAL", 1);
    }
  
  Double_t params[4];
  const char *cutmodel;
  cutmodel = "pol0";
  params[0] = -1.0; //sigma min
  double maxnSig = 3.0;
  if(fmcData)
    {
     params[0] = -5.0; //sigma min
     maxnSig = 5.0; 
    } 

  for(Int_t a=0;a<11;a++)fPID->ConfigureTPCcentralityCut(a,cutmodel,params,maxnSig);

  fPID->SortDetectors(); 
  fPIDqa = new AliHFEpidQAmanager();
  fPIDqa->Initialize(fPID);
 
  //------- fcut --------------  
  fCuts = new AliHFEcuts();
  fCuts->CreateStandardCuts();
  //fCuts->SetMinNClustersTPC(100);
  fCuts->SetMinNClustersTPC(90);
  fCuts->SetMinRatioTPCclusters(0.6);
  fCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  //fCuts->SetMinNClustersITS(3);
  fCuts->SetMinNClustersITS(2);
  fCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  fCuts->SetCheckITSLayerStatus(kFALSE);
  fCuts->SetVertexRange(10.);
  fCuts->SetTOFPIDStep(kFALSE);
  fCuts->SetPtRange(2, 50);
  fCuts->SetMaxImpactParam(3.,3.);

  //--------Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
  fCFM->SetNStepParticle(kNcutSteps);
  for(Int_t istep = 0; istep < kNcutSteps; istep++)
    fCFM->SetParticleCutsList(istep, NULL);
  
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  fCuts->Initialize(fCFM);
  
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->Add(fPIDqa->MakeList("PIDQA"));
  
  fNoEvents = new TH1F("fNoEvents","",4,-0.5,3.5) ;
  fOutputList->Add(fNoEvents);
  
  Int_t binsE[5] =    {250, 100,  1000, 200,   10};
  Double_t xminE[5] = {1.0,  -1,   0.0,   0, -0.5}; 
  Double_t xmaxE[5] = {3.5,   1, 100.0, 100,  9.5}; 
  fEMCAccE = new THnSparseD("fEMCAccE","EMC acceptance & E;#eta;#phi;Energy;Centrality;trugCondition;",5,binsE,xminE,xmaxE);
  if(fqahist==1)fOutputList->Add(fEMCAccE);

  hEMCAccE = new TH2F("hEMCAccE","Cluster Energy",200,0,100,100,0,20);
  fOutputList->Add(hEMCAccE);

  fTrkpt = new TH1F("fTrkpt","track pt",100,0,50);
  fOutputList->Add(fTrkpt);
  
  fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",100,0,50);
  fOutputList->Add(fTrackPtBefTrkCuts);
  
  fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",100,0,50);
  fOutputList->Add(fTrackPtAftTrkCuts);
  
  fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma",100,0,50,200,-10,10);
  fOutputList->Add(fTPCnsigma);
  
  fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",100,0,50,100,0,2);
  fOutputList->Add(fTrkEovPBef);
  
  fTrkEovPAft = new TH2F("fTrkEovPAft","track E/p after HFE pid",100,0,50,100,0,2);
  fOutputList->Add(fTrkEovPAft);
  
  fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",100,0,50,200,-10,10);
  fOutputList->Add(fdEdxBef);
  
  fdEdxAft = new TH2F("fdEdxAft","track dEdx vs p after HFE pid",100,0,50,200,-10,10);
  fOutputList->Add(fdEdxAft);
  
  fIncpT = new TH2F("fIncpT","HFE pid electro vs. centrality",200,0,100,100,0,50);
  fOutputList->Add(fIncpT);

  fIncpTM20 = new TH2F("fIncpTM20","HFE pid electro vs. centrality with M20",200,0,100,100,0,50);
  fOutputList->Add(fIncpTM20);
  
  Int_t nBinspho[9] =  { 200, 100, 500, 12,   50,    4, 200,   8, 100};
  Double_t minpho[9] = {  0.,  0.,  0., -2.5,  0, -0.5,   0,-1.5,   0};   
  Double_t maxpho[9] = {100., 50., 0.5, 3.5,   1,  3.5,   2, 6.5,  50};   

  fInvmassLS = new THnSparseD("fInvmassLS", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2); nSigma; angle; m20cut; eop; Mcele;", 9, nBinspho,minpho, maxpho);
  if(fqahist==1)fOutputList->Add(fInvmassLS);
  
  fInvmassULS = new THnSparseD("fInvmassULS", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2); nSigma; angle; m20cut; eop; MCele", 9, nBinspho,minpho, maxpho);
  if(fqahist==1)fOutputList->Add(fInvmassULS);
  
  fInvmassLSmc = new THnSparseD("fInvmassLSmc", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2); nSigma; angle; m20cut; eop; Mcele;", 9, nBinspho,minpho, maxpho);
  if(fqahist==1)fOutputList->Add(fInvmassLSmc);
  
  fInvmassULSmc = new THnSparseD("fInvmassULSmc", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2); nSigma; angle; m20cut; eop; MCele", 9, nBinspho,minpho, maxpho);
  if(fqahist==1)fOutputList->Add(fInvmassULSmc);

  fInvmassLSmc0 = new TH2D("fInvmassLSmc0", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassLSmc0->Sumw2();
  fOutputList->Add(fInvmassLSmc0);
  
  fInvmassLSmc1 = new TH2D("fInvmassLSmc1", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassLSmc1->Sumw2();
  fOutputList->Add(fInvmassLSmc1);

  fInvmassLSmc2 = new TH2D("fInvmassLSmc2", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassLSmc2->Sumw2();
  fOutputList->Add(fInvmassLSmc2);

  fInvmassLSmc3 = new TH2D("fInvmassLSmc3", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassLSmc3->Sumw2();
  fOutputList->Add(fInvmassLSmc3);

  fInvmassULSmc0 = new TH2D("fInvmassULSmc0", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassULSmc0->Sumw2();
  fOutputList->Add(fInvmassULSmc0);

  fInvmassULSmc1 = new TH2D("fInvmassULSmc1", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassULSmc1->Sumw2();
  fOutputList->Add(fInvmassULSmc1);

  fInvmassULSmc2 = new TH2D("fInvmassULSmc2", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassULSmc2->Sumw2();
  fOutputList->Add(fInvmassULSmc2);

  fInvmassULSmc3 = new TH2D("fInvmassULSmc3", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2)",20,0,20,500,0,0.5 );
  fInvmassULSmc3->Sumw2();
  fOutputList->Add(fInvmassULSmc3);

  fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleLS);
  
  fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleULS);
  
  fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",100,0,50);
  fOutputList->Add(fPhotoElecPt);
  
  fPhoElecPt = new TH2F("fPhoElecPt", "Pho-inclusive electron pt",200,0,100,100,0,50);
  fOutputList->Add(fPhoElecPt);
  
  fPhoElecPtM20 = new TH2F("fPhoElecPtM20", "Pho-inclusive electron pt with M20",200,0,100,100,0,50);
  fOutputList->Add(fPhoElecPtM20);

  fSameElecPt = new TH2F("fSameElecPt", "Same-inclusive electron pt",200,0,100,100,0,50);
  fOutputList->Add(fSameElecPt);

  fSameElecPtM20 = new TH2F("fSameElecPtM20", "Same-inclusive electron pt with M20",200,0,100,100,0,50);
  fOutputList->Add(fSameElecPtM20);

  fCent = new TH1F("fCent","Centrality",200,0,100) ;
  fOutputList->Add(fCent);
 
  // Make common binning
  const Double_t kMinP = 0.;
  const Double_t kMaxP = 50.;
  //const Double_t kTPCSigMim = 40.;
  //const Double_t kTPCSigMax = 140.;

  // 1st histogram: TPC dEdx with/without EMCAL (p, pT, TPC Signal, phi, eta,  Sig,  e/p,  ,match, cell, M02, M20, Disp, Centrality, select)
  Int_t nBins[16] =  {  250,    10,   60,    20,   100,  300,  50,   40,   200, 200,  250, 200,    3,    5, 100,    8};
  Double_t min[16] = {kMinP,  -0.5, 1.0,  -1.0,  -6.0,    0,   0,    0,   0.0, 0.0,  0.0,   0, -1.5, -0.5,  80, -1.5};
  Double_t max[16] = {kMaxP,   9.5, 4.0,   1.0,   4.0,  3.0, 0.1,   40,   200, 2.0, 50.0, 100,  1.5,  4.5, 180,  6.5};
  fEleInfo = new THnSparseD("fEleInfo", "Electron Info; pT [GeV/c]; TPC signal;phi;eta;nSig; E/p;Rmatch;Ncell;clsF;M20;mcpT;Centrality;charge;opp;same;trigCond;MCele", 16, nBins, min, max);
  if(fqahist==1)fOutputList->Add(fEleInfo);

  //<---  Trigger info
  /*
  fClsEBftTrigCut = new TH1F("fClsEBftTrigCut","cluster E before trigger selection",1000,0,100);
  fOutputList->Add(fClsEBftTrigCut);

  fClsEAftTrigCut = new TH1F("fClsEAftTrigCut","cluster E if cls has 0 trigcut channel",1000,0,100);
  fOutputList->Add(fClsEAftTrigCut);

  fClsEAftTrigCut1 = new TH1F("fClsEAftTrigCut1","cluster E if cls with trig channel",1000,0,100);
  fOutputList->Add(fClsEAftTrigCut1);

  fClsEAftTrigCut2 = new TH1F("fClsEAftTrigCut2","cluster E if cls with trigcut channel",1000,0,100);
  fOutputList->Add(fClsEAftTrigCut2);

  fClsEAftTrigCut3 = new TH1F("fClsEAftTrigCut3","cluster E if cls with trigcut channel + nCell>Ecorrect",1000,0,100);
  fOutputList->Add(fClsEAftTrigCut3);

  fClsEAftTrigCut4 = new TH1F("fClsEAftTrigCut4","cluster E if cls with trigcut channel + nCell>Ecorrect + cls time cut",1000,0,100);
  fOutputList->Add(fClsEAftTrigCut4);

  fClsETime = new TH2F("fClsETime", "Cls time vs E; E; time;",1000,0,100,1000,-0.0000002,0.0000009);
  fOutputList->Add(fClsETime);

  fClsETime1 = new TH2F("fClsETime1", "Cls time vs E if cls contains trigger channel; E; time;",1000,0,100,1000,-0.0000002,0.0000009);
  fOutputList->Add(fClsETime1);

  fTrigTimes = new TH1F("fTrigTimes", "Trigger time; time; N;",25,0,25);
  fOutputList->Add(fTrigTimes);

  fCellCheck = new TH2F("fCellCheck", "Cell vs E; E GeV; Cell ID",10,6,26,12000,0,12000);
  fOutputList->Add(fCellCheck);
  */
  //<---------- MC

  fInputHFEMC = new TH2F("fInputHFEMC","Input MC HFE pid electro vs. centrality",200,0,100,100,0,50);
  fOutputList->Add(fInputHFEMC);

  fInputAlle = new TH2F("fInputAlle","Input MC electro vs. centrality",200,0,100,100,0,50);
  fOutputList->Add(fInputAlle);

  fIncpTMChfe = new TH2F("fIncpTMChfe","MC HFE pid electro vs. centrality",200,0,100,100,0,50);
  fOutputList->Add(fIncpTMChfe);

  fIncpTMChfeAll = new TH2F("fIncpTMChfeAll","MC Alle pid electro vs. centrality",200,0,100,100,0,50);
  fOutputList->Add(fIncpTMChfeAll);

  fIncpTMCM20hfe = new TH2F("fIncpTMCM20hfe","MC HFE pid electro vs. centrality with M20",200,0,100,100,0,50);
  fOutputList->Add(fIncpTMCM20hfe);

  fIncpTMCM20hfeAll = new TH2F("fIncpTMCM20hfeAll","MC Alle pid electro vs. centrality with M20",200,0,100,100,0,50);
  fOutputList->Add(fIncpTMCM20hfeAll);

  fIncpTMCM20hfeCheck = new TH2F("fIncpTMCM20hfeCheck","MC HFE pid electro vs. centrality with M20 Check",200,0,100,100,0,50);
  fOutputList->Add(fIncpTMCM20hfeCheck);

  Int_t nBinspho2[5] =  { 200, 100,    7,    3, 700};
  Double_t minpho2[5] = {  0.,  0., -2.5, -0.5, 0.};   
  Double_t maxpho2[5] = {100., 50.,  4.5,  2.5, 70.};   
  
  fInputHFEMC_weight = new THnSparseD("fInputHFEMC_weight", "MC HFE electron pt",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fInputHFEMC_weight);
  
  fIncpTMCM20hfeCheck_weight = new THnSparseD("fIncpTMCM20hfeCheck_weight", "HFE electron pt with M20",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fIncpTMCM20hfeCheck_weight);
  
  fIncpTMCpho = new THnSparseD("fIncpTMCpho","MC Pho pid electro vs. centrality",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fIncpTMCpho);

  fIncpTMCM20pho = new THnSparseD("fIncpTMCM20pho","MC Pho pid electro vs. centrality with M20",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fIncpTMCM20pho);

  fPhoElecPtMC = new THnSparseD("fPhoElecPtMC", "MC Pho-inclusive electron pt",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fPhoElecPtMC);
  
  fPhoElecPtMCM20 = new THnSparseD("fPhoElecPtMCM20", "MC Pho-inclusive electron pt with M20",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fPhoElecPtMCM20);

  fSameElecPtMC = new THnSparseD("fSameElecPtMC", "MC Same-inclusive electron pt",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fSameElecPtMC);

  fSameElecPtMCM20 = new THnSparseD("fSameElecPtMCM20", "MC Same-inclusive electron pt with M20",5,nBinspho2,minpho2,maxpho2);
  fOutputList->Add(fSameElecPtMCM20);

  fIncpTMCM20pho_pi0e = new THnSparseD("fIncpTMCM20pho_pi0e","MC Pho pi0->e pid electro vs. centrality with M20",5,nBinspho2,minpho2,maxpho2);
  fIncpTMCM20pho_pi0e->Sumw2();
  fOutputList->Add(fIncpTMCM20pho_pi0e);

  fPhoElecPtMCM20_pi0e = new THnSparseD("fPhoElecPtMCM20_pi0e", "MC Pho-inclusive electron pt with M20 pi0->e",5,nBinspho2,minpho2,maxpho2);
  fPhoElecPtMCM20_pi0e->Sumw2();
  fOutputList->Add(fPhoElecPtMCM20_pi0e);

  fSameElecPtMCM20_pi0e = new THnSparseD("fSameElecPtMCM20_pi0e", "MC Same-inclusive electron pt pi0->e",5,nBinspho2,minpho2,maxpho2);
  fSameElecPtMCM20_pi0e->Sumw2();
  fOutputList->Add(fSameElecPtMCM20_pi0e);
 // 
  fIncpTMCM20pho_eta = new THnSparseD("fIncpTMCM20pho_eta","MC Pho pi0->e pid electro vs. centrality with M20",5,nBinspho2,minpho2,maxpho2);
  fIncpTMCM20pho_eta->Sumw2();
  fOutputList->Add(fIncpTMCM20pho_eta);

  fPhoElecPtMCM20_eta = new THnSparseD("fPhoElecPtMCM20_eta", "MC Pho-inclusive electron pt with M20 pi0->e",5,nBinspho2,minpho2,maxpho2);
  fPhoElecPtMCM20_eta->Sumw2();
  fOutputList->Add(fPhoElecPtMCM20_eta);

  fSameElecPtMCM20_eta = new THnSparseD("fSameElecPtMCM20_eta", "MC Same-inclusive electron pt pi0->e",5,nBinspho2,minpho2,maxpho2);
  fSameElecPtMCM20_eta->Sumw2();
  fOutputList->Add(fSameElecPtMCM20_eta);
  // ------------
  fIncpTMCpho_pi0e_TPC = new THnSparseD("fIncpTMCpho_pi0e_TPC","MC Pho pi0->e pid electro vs. centrality with M20",5,nBinspho2,minpho2,maxpho2);
  fIncpTMCpho_pi0e_TPC->Sumw2();
  fOutputList->Add(fIncpTMCpho_pi0e_TPC);

  fPhoElecPtMC_pi0e_TPC = new THnSparseD("fPhoElecPtMC_pi0e_TPC", "MC Pho-inclusive electron pt with  pi0->e",5,nBinspho2,minpho2,maxpho2);
  fPhoElecPtMC_pi0e_TPC->Sumw2();
  fOutputList->Add(fPhoElecPtMC_pi0e_TPC);

  fSameElecPtMC_pi0e_TPC = new THnSparseD("fSameElecPtMC_pi0e_TPC", "MC Same-inclusive electron pt pi0->e",5,nBinspho2,minpho2,maxpho2);
  fSameElecPtMC_pi0e_TPC->Sumw2();
  fOutputList->Add(fSameElecPtMC_pi0e_TPC);
  //-------------

  CheckNclust = new TH1D("CheckNclust","cluster check",200,0,200);
  fOutputList->Add(CheckNclust);

  CheckNits = new TH1D("CheckNits","ITS cluster check",8,-0.5,7.5);
  fOutputList->Add(CheckNits);
  /*
  Hpi0pTcheck = new TH2D("Hpi0pTcheck","Pi0 pT from Hijing",100,0,50,3,-0.5,2.5);
  fOutputList->Add(Hpi0pTcheck);

  HETApTcheck = new TH2D("HETApTcheck","Eta pT from Hijing",100,0,50,3,-0.5,2.5);
  fOutputList->Add(HETApTcheck);
  */

  Int_t nBinspho3[3] =  { 200, 100, 3};
  Double_t minpho3[3] = {  0.,  0., -0.5};   
  Double_t maxpho3[3] = {100., 50., 2.5};   

  Hpi0pTcheck = new THnSparseD("Hpi0pTcheck","Pi0 pT from Hijing",3,nBinspho3,minpho3,maxpho3);
  fOutputList->Add(Hpi0pTcheck);

  HETApTcheck = new THnSparseD("HETApTcheck","Eta pT from Hijing",3,nBinspho3,minpho3,maxpho3);
  fOutputList->Add(HETApTcheck);
  //--
  HphopTcheck = new TH2D("HphopTcheck","Pho pT from Hijing",100,0,50,3,-0.5,2.5);
  fOutputList->Add(HphopTcheck);
  //
  HDpTcheck = new TH2D("HDpTcheck","D pT from Hijing",100,0,50,3,-0.5,2.5);
  fOutputList->Add(HDpTcheck);

  HBpTcheck = new TH2D("HBpTcheck","B pT from Hijing",100,0,50,3,-0.5,2.5);
  fOutputList->Add(HBpTcheck);
  //

  fpTCheck = new TH1D("fpTCheck","pT check",500,0,50);
  fOutputList->Add(fpTCheck);

  fMomDtoE = new TH2D("fMomDtoE","D->E pT correlations;e p_{T} GeV/c;D p_{T} GeV/c",400,0,40,400,0,40);
  fOutputList->Add(fMomDtoE);

  fLabelCheck = new TH2D("fLabelCheck","MC label",50,0,50,5,-1.5,3.5);
  fOutputList->Add(fLabelCheck);

  fgeoFake = new TH2D("fgeoFake","Label==0 eta and phi",628,0,6.28,200,-1,1);
  fOutputList->Add(fgeoFake);

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHFECal::Terminate(Option_t *)
{
  // Info("Terminate");
	AliAnalysisTaskSE::Terminate();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHFECal::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
  // Check single track cuts for a given cut step
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  return kTRUE;
}
//_________________________________________
//void AliAnalysisTaskHFECal::SelectPhotonicElectron(Int_t itrack, Double_t cent, AliESDtrack *track, Bool_t &fFlagPhotonicElec)
//void AliAnalysisTaskHFECal::SelectPhotonicElectron(Int_t itrack, Double_t cent, AliESDtrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec, Double_t nSig)
void AliAnalysisTaskHFECal::SelectPhotonicElectron(Int_t itrack, Double_t cent, AliESDtrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec, Double_t nSig, Double_t shower, Double_t ep, Double_t mce, Double_t w, Int_t ibgevent, Bool_t tagpi0, Bool_t tageta)
{
  //Identify non-heavy flavour electrons using Invariant mass method
  
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetRequireITSRefit(kTRUE);
  fTrackCuts->SetEtaRange(-0.9,0.9);
  //fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts->SetMinNClustersTPC(90);
  
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  double ptEle = track->Pt();  //add
  if(ibgevent==0 && w > 0.0)
     {
      fpTCheck->Fill(ptEle,w);
     }

  Bool_t flagPhotonicElec = kFALSE;
  Bool_t flagConvinatElec = kFALSE;
  
  int p1 = 0;
  if(mce==3)
     {
       Int_t label = TMath::Abs(track->GetLabel());
       TParticle* particle = stack->Particle(label);
       p1 = particle->GetFirstMother();
     }

  //for(Int_t jTracks = itrack+1; jTracks<fESD->GetNumberOfTracks(); jTracks++){
  for(Int_t jTracks = 0; jTracks<fESD->GetNumberOfTracks(); jTracks++){
    AliESDtrack* trackAsso = fESD->GetTrack(jTracks);
    if (!trackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }
    if(itrack==jTracks)continue;    
    int jbgevent = 0;    

    int p2 = 0;
    if(mce==3)
    {
      Int_t label2 = TMath::Abs(trackAsso->GetLabel());
      TParticle* particle2 = stack->Particle(label2);
      Bool_t MChijing_ass = fMC->IsFromBGEvent(label2);
      if(MChijing_ass)jbgevent =1;
      if(particle2->GetFirstMother()>-1)
         p2 = particle2->GetFirstMother();
    }

    Double_t dEdxAsso = -999., ptPrim=-999., ptAsso=-999., openingAngle = -999.;
    Double_t mass=999., width = -999;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    //ptPrim = track->Pt();
    ptPrim = ptEle;

    dEdxAsso = trackAsso->GetTPCsignal();
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();
    

    if(ptAsso <0.5) continue;
    if(!fTrackCuts->AcceptTrack(trackAsso)) continue;
    if(dEdxAsso <65 || dEdxAsso>100) continue; //11a pass1
    
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
 
    //printf("chargeAsso = %d\n",chargeAsso);
    //printf("charge = %d\n",charge);
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;
    
    //printf("fFlagLS = %d\n",fFlagLS);
    //printf("fFlagULS = %d\n",fFlagULS);
    //printf("\n");

    AliKFParticle ge1(*track, fPDGe1);
    AliKFParticle ge2(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);
    
    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
    
    //  for v5-03-70-AN comment
    AliKFVertex primV(*pVtx);
    primV += recg;
    recg.SetProductionVertex(primV);
    
    recg.SetMassConstraint(0,0.0001);
    
    openingAngle = ge1.GetAngle(ge2);
    if(fFlagLS) fOpeningAngleLS->Fill(openingAngle);
    if(fFlagULS) fOpeningAngleULS->Fill(openingAngle);
    
    
    recg.GetMass(mass,width);
    
    double ishower = 0;
    if(shower>0.0 && shower<0.3)ishower = 1;

    double phoinfo[9];
    phoinfo[0] = cent;
    phoinfo[1] = ptPrim;
    phoinfo[2] = mass;
    phoinfo[3] = nSig;
    phoinfo[4] = openingAngle;
    phoinfo[5] = ishower;
    phoinfo[6] = ep;
    phoinfo[7] = mce;
    phoinfo[8] = ptAsso;

    if(fFlagLS) fInvmassLS->Fill(phoinfo);
    if(fFlagULS) fInvmassULS->Fill(phoinfo);
    if(fFlagLS && ibgevent==0 && jbgevent==0) fInvmassLSmc->Fill(phoinfo,w);
    if(fFlagULS && ibgevent==0 && jbgevent==0)
       {
         fInvmassULSmc->Fill(phoinfo,w);
       }
    //printf("fInvmassCut %f\n",fInvmassCut);
    //printf("openingAngle %f\n",fOpeningAngleCut);

    if(openingAngle > fOpeningAngleCut) continue;
    
    // for real data  
    //printf("mce =%f\n",mce);
    if(mce<-0.5) // mce==-1. is real
       {
         //printf("Real data\n");
	 if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
	 //if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec && (p1==p2)){ <--- only MC train (55,56) v5-03-68-AN & 69 for check
	       flagPhotonicElec = kTRUE;
	      }
	 if(mass<fInvmassCut && fFlagLS && !flagConvinatElec){
	       flagConvinatElec = kTRUE;
	      }
        }
    // for MC data  
    else
       {
         //printf("MC data\n");

         if(w>0.0)
           {
           //cout << "tagpi0 = " << tagpi0 << " ; tageta = " << tageta << endl;
           if(fFlagLS && ibgevent==0 && jbgevent==0 && tagpi0) fInvmassLSmc0->Fill(ptPrim,mass,w);
           if(fFlagULS && ibgevent==0 && jbgevent==0 && tagpi0) fInvmassULSmc0->Fill(ptPrim,mass,w);
           if(fFlagLS && ibgevent==0 && jbgevent==0 && tageta) fInvmassLSmc1->Fill(ptPrim,mass,w);
           if(fFlagULS && ibgevent==0 && jbgevent==0 && tageta) fInvmassULSmc1->Fill(ptPrim,mass,w);
           if(fFlagLS && ibgevent==0 && jbgevent==0 && (p1==p2) && tagpi0) fInvmassLSmc2->Fill(ptPrim,mass,w);
           if(fFlagULS && ibgevent==0 && jbgevent==0 && (p1==p2) && tagpi0) fInvmassULSmc2->Fill(ptPrim,mass,w);
           if(fFlagLS && ibgevent==0 && jbgevent==0 && (p1==p2) && tageta) fInvmassLSmc3->Fill(ptPrim,mass,w);
           if(fFlagULS && ibgevent==0 && jbgevent==0 && (p1==p2) && tageta) fInvmassULSmc3->Fill(ptPrim,mass,w);
          }
	 if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec && (ibgevent==jbgevent)){
	 //if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec && (p1==p2)){ <--- only MC train (55,56) v5-03-68-AN & 69 for check
	       flagPhotonicElec = kTRUE;
	      }
	 if(mass<fInvmassCut && fFlagLS && !flagConvinatElec && (ibgevent==jbgevent)){
	       flagConvinatElec = kTRUE;
	      }
        }

  }
  fFlagPhotonicElec = flagPhotonicElec;
  fFlagConvinatElec = flagConvinatElec;
  
}
//-------------------------------------------

void AliAnalysisTaskHFECal::FindMother(TParticle* part, int &label, int &pid)
{
 //int label = 99999;
 //int pid = 99999;

 if(part->GetFirstMother()>-1)
   {
    label = part->GetFirstMother();
    pid = stack->Particle(label)->GetPdgCode();
   }
   //cout << "Find Mother : label = " << label << " ; pid" << pid << endl;
}

double AliAnalysisTaskHFECal::GetMCweight(double mcPi0pT)
{
        double weight = 1.0;

	if(mcPi0pT>0.0 && mcPi0pT<5.0)
	{
		weight = 0.323*mcPi0pT/(TMath::Exp(-1.6+0.767*mcPi0pT+0.0285*mcPi0pT*mcPi0pT));
	}
	else
	{
		weight = 115.0/(0.718*mcPi0pT*TMath::Power(mcPi0pT,3.65));
	}
  return weight;
}

double AliAnalysisTaskHFECal::GetMCweightEta(double mcEtapT)
{
  double weight = 1.0;

  weight = 223.3/TMath::Power((TMath::Exp(-0.17*mcEtapT-0.0322*mcEtapT*mcEtapT)+mcEtapT/1.69),5.65);
  return weight;
}


//_________________________________________
void AliAnalysisTaskHFECal::FindTriggerClusters()
{
  //cout << "finding trigger patch" << endl; 
  // constants
  const int nModuleCols = 2;
  const int nModuleRows = 5;
  const int nColsFeeModule = 48;
  const int nRowsFeeModule = 24;
  const int nColsFaltroModule = 24;
  const int nRowsFaltroModule = 12;
  //const int faltroWidthMax = 20;

  // part 1, trigger extraction -------------------------------------
  Int_t globCol, globRow;
  //Int_t ntimes=0, nTrigChannel=0, nTrigChannelCut=0, trigInCut=0;
  Int_t ntimes=0, nTrigChannel=0, nTrigChannelCut=0;

  //Int_t trigtimes[faltroWidthMax];
  Double_t cellTime[nColsFeeModule*nModuleCols][nRowsFeeModule*nModuleRows];
  Double_t cellEnergy[nColsFeeModule*nModuleCols][nRowsFeeModule*nModuleRows];
  //Double_t fTrigCutLow = 6;
  //Double_t fTrigCutHigh = 10;
  Double_t fTimeCutLow =  469e-09;
  Double_t fTimeCutHigh = 715e-09;

  AliESDCaloTrigger * fCaloTrigger = fESD->GetCaloTrigger( "EMCAL" );

  // erase trigger maps
  for(Int_t i = 0; i < nColsFaltroModule*nModuleCols; i++ )
  {
    for(Int_t j = 0; j < nRowsFaltroModule*nModuleRows; j++ )
    {
      ftriggersCut[i][j] = 0;
      ftriggers[i][j] = 0;
      ftriggersTime[i][j] = 0;
    }
  }

  Int_t iglobCol=0, iglobRow=0;
  // go through triggers
  if( fCaloTrigger->GetEntries() > 0 )
  {
    // needs reset
    fCaloTrigger->Reset();
    while( fCaloTrigger->Next() )
    {
      fCaloTrigger->GetPosition( globCol, globRow );
      fCaloTrigger->GetNL0Times( ntimes );
      /*
      // no L0s
      if( ntimes < 1 )   continue;
      // get precise timings
      fCaloTrigger->GetL0Times( trigtimes );
      trigInCut = 0;
      for(Int_t i = 0; i < ntimes; i++ )
      {
        // save the first trigger time in channel
        if( i == 0 || triggersTime[globCol][globRow] > trigtimes[i] )
          triggersTime[globCol][globRow] = trigtimes[i];
        //printf("trigger times: %d\n",trigtimes[i]);
        // check if in cut
        if(trigtimes[i] > fTrigCutLow && trigtimes[i] < fTrigCutHigh )
          trigInCut = 1;

        fTrigTimes->Fill(trigtimes[i]);
      }
      */
   
      //L1 analysis from AliAnalysisTaskEMCALTriggerQA
      Int_t bit = 0;
      fCaloTrigger->GetTriggerBits(bit);
      //cout << "bit = " << bit << endl;

      Int_t ts = 0;
      fCaloTrigger->GetL1TimeSum(ts);
      //cout << "ts = " << ts << endl;
      if (ts > 0)ftriggers[globCol][globRow] = 1;
      // number of triggered channels in event
      nTrigChannel++;
      // ... inside cut
      if(ts>0 && (bit >> 6 & 0x1))
      {
        iglobCol = globCol;
        iglobRow = globRow;
        nTrigChannelCut++;
        //cout << "ts cut = " << ts << endl;
        //cout << "globCol = " << globCol << endl;
        //cout << "globRow = " << globRow << endl;
        ftriggersCut[globCol][globRow] = 1;
      }

    } // calo trigger entries
  } // has calo trigger entries

  // part 2 go through the clusters here -----------------------------------
  //cout << " part 2 go through the clusters here ----------------------------------- " << endl; 
  Int_t nCluster=0, nCell=0, iCell=0, gCell=0;
  Short_t cellAddr, nSACell;
  Int_t mclabel;
  //Int_t nSACell, iSACell, mclabel;
  Int_t iSACell;
  Double_t cellAmp=0, cellTimeT=0, clusterTime=0, efrac=0;
  Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta, gphi, geta, feta, fphi;

  TRefArray *fCaloClusters = new TRefArray();
  fESD->GetEMCALClusters( fCaloClusters ); 
  nCluster = fCaloClusters->GetEntries();


  // save all cells times if there are clusters  
  if( nCluster > 0 ){
    // erase time map
    for(Int_t i = 0; i < nColsFeeModule*nModuleCols; i++ ){ 
      for(Int_t j = 0; j < nRowsFeeModule*nModuleRows; j++ ){
        cellTime[i][j] = 0.;
        cellEnergy[i][j] = 0.;
      }
    }

    // get cells
    AliESDCaloCells *fCaloCells = fESD->GetEMCALCells();
    //AliVCaloCells fCaloCells = fESD->GetEMCALCells();
    nSACell = fCaloCells->GetNumberOfCells();
    for(iSACell = 0; iSACell < nSACell; iSACell++ ){ 
      // get the cell info *fCal
      fCaloCells->GetCell( iSACell, cellAddr, cellAmp, cellTimeT , mclabel, efrac);

      // get cell position 
      fGeom->GetCellIndex( cellAddr, nSupMod, nModule, nIphi, nIeta ); 
      fGeom->GetCellPhiEtaIndexInSModule( nSupMod,nModule, nIphi, nIeta, iphi, ieta);

      // convert co global phi eta  
      gphi = iphi + nRowsFeeModule*(nSupMod/2);
      geta = ieta + nColsFeeModule*(nSupMod%2);

      // save cell time and energy
      cellTime[geta][gphi] = cellTimeT;
      cellEnergy[geta][gphi] = cellAmp;

    }
  }

  Int_t nClusterTrig, nClusterTrigCut;
  UShort_t *cellAddrs;
  Double_t clsE=-999, clsEta=-999, clsPhi=-999;
  Float_t clsPos[3] = {0.,0.,0.};

  for(Int_t icl=0; icl<fESD->GetNumberOfCaloClusters(); icl++)
  {
    AliESDCaloCluster *cluster = fESD->GetCaloCluster(icl);
    if(!cluster || !cluster->IsEMCAL()) continue;

    // get cluster cells
    nCell = cluster->GetNCells();

    // get cluster energy
    clsE = cluster->E();

    // get cluster position
    cluster->GetPosition(clsPos);
    TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);
    clsEta = clsPosVec.Eta();
    clsPhi = clsPosVec.Phi();

    // get the cell addresses
    cellAddrs = cluster->GetCellsAbsId();

    // check if the cluster contains cell, that was marked as triggered
    nClusterTrig = 0;
    nClusterTrigCut = 0;

    // loop the cells to check, if cluser in acceptance
    // any cluster with a cell outside acceptance is not considered
    for( iCell = 0; iCell < nCell; iCell++ )
    {
     // check hot cell
     //if(clsE>6.0)fCellCheck->Fill(clsE,cellAddrs[iCell]); 

      // get cell position
      fGeom->GetCellIndex( cellAddrs[iCell], nSupMod, nModule, nIphi, nIeta );
      fGeom->GetCellPhiEtaIndexInSModule( nSupMod,nModule, nIphi, nIeta, iphi, ieta);

      // convert co global phi eta
      gphi = iphi + nRowsFeeModule*(nSupMod/2);
      geta = ieta + nColsFeeModule*(nSupMod%2);

      if( cellTime[geta][gphi] > 0. ){ 
        clusterTime += cellTime[geta][gphi];
        gCell++;
      }

      // get corresponding FALTRO
      fphi = gphi / 2;
      feta = geta / 2;

        //cout << "fphi = " << fphi << endl;
        //cout << "feta = " << feta << endl;

      // try to match with a triggered
      if( ftriggers[feta][fphi]==1)
      {  nClusterTrig++;
      }
      if( ftriggersCut[feta][fphi]==1)
      { nClusterTrigCut++;
      }

      //cout << "nClusterTrigCut : " << nClusterTrigCut << endl;
     
    } // cells


    if( gCell > 0 ) 
      clusterTime = clusterTime / (Double_t)gCell;
    // fix the reconstriction code time 100ns jumps
    if( fESD->GetBunchCrossNumber() % 4 < 2 )
      clusterTime -= 0.0000001;

    //fClsETime->Fill(clsE,clusterTime);
    //fClsEBftTrigCut->Fill(clsE);

    if(nClusterTrig>0){
      //fClsETime1->Fill(clsE,clusterTime);
    }

    if(nClusterTrig>0){
      cluster->SetChi2(1);
      //fClsEAftTrigCut1->Fill(clsE);                                               
    }

    if(nClusterTrigCut>0){
      cluster->SetChi2(2);
      //fClsEAftTrigCut2->Fill(clsE);
    }

    if(nClusterTrigCut>0 && ( nCell > (1 + clsE / 3)))
    {
      cluster->SetChi2(3);
      //fClsEAftTrigCut3->Fill(clsE);
    }

    if(nClusterTrigCut>0 && (nCell > (1 + clsE / 3) )&&( clusterTime > fTimeCutLow && clusterTime < fTimeCutHigh ))
    {
      // cluster->SetChi2(4);
      //fClsEAftTrigCut4->Fill(clsE);
    }
    if(nClusterTrigCut<1)
    {
      cluster->SetChi2(0);

      //fClsEAftTrigCut->Fill(clsE);
    }

  } // clusters
}


double AliAnalysisTaskHFECal::NsigCorr(float cent)
{
 double shift = 0.0;
 if(cent>=20 && cent<30)shift = 0.156;
 if(cent>=30 && cent<40)shift = 0.316;
 if(cent>=40 && cent<50)shift = 0.336;
 if(cent>=50 && cent<70)shift = 0.440;
 if(cent>=70 && cent<90)shift = 0.534;
 return shift;
}




