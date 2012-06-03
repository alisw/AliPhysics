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
  ,fGeom(0)
  ,fOutputList(0)
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
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fIncpT(0)	
  ,fIncpTM20(0)	
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
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
  ,fInputHFEMC(0)
  ,fIncpTMChfe(0)	
  ,fIncpTMCM20hfe(0)	
  ,fIncpTMCpho(0)	
  ,fIncpTMCM20pho(0)	
  ,fPhoElecPtMC(0)
  ,fPhoElecPtMCM20(0)
  ,fSameElecPtMC(0)
  ,fSameElecPtMCM20(0)
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
  ,fGeom(0)
  ,fOutputList(0)
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
  ,fTrkpt(0)
  ,fTrkEovPBef(0)	 
  ,fTrkEovPAft(0)	 
  ,fdEdxBef(0)	 
  ,fdEdxAft(0)	 
  ,fIncpT(0)	 
  ,fIncpTM20(0)	 
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
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
  ,fInputHFEMC(0)
  ,fIncpTMChfe(0)	
  ,fIncpTMCM20hfe(0)	
  ,fIncpTMCpho(0)	
  ,fIncpTMCM20pho(0)	
  ,fPhoElecPtMC(0)
  ,fPhoElecPtMCM20(0)
  ,fSameElecPtMC(0)
  ,fSameElecPtMCM20(0)
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
  AliStack* stack = NULL;
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
      Bool_t mcInDtoE= kFALSE;
      Bool_t mcInBtoE= kFALSE;

      if(particle->GetFirstMother()>-1 && fabs(fPDG)==11)
        {
	    int parentPID = stack->Particle(particle->GetFirstMother())->GetPdgCode();  
            if((fabs(parentPID)==411 || fabs(parentPID)==413 || fabs(parentPID)==421 || fabs(parentPID)==423 || fabs(parentPID)==431)&& fabs(fPDG)==11)mcInDtoE = kTRUE;
            if((fabs(parentPID)==511 || fabs(parentPID)==513 || fabs(parentPID)==521 || fabs(parentPID)==523 || fabs(parentPID)==531)&& fabs(fPDG)==11)mcInBtoE = kTRUE;
            if((mcInBtoE || mcInDtoE) && fabs(mcZvertex)<10.0)fInputHFEMC->Fill(cent,pTMC);
         }

      }
    } 

  fNoEvents->Fill(0);

  Int_t fNOtrks =  fESD->GetNumberOfTracks();
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();
  
  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(1);
  
  if(fNOtrks<2) return;
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
  
  if(cent>90.) return;
	
 // Calorimeter info.
 
   FindTriggerClusters();

  // make EMCAL array 
  for(Int_t iCluster=0; iCluster<fESD->GetNumberOfCaloClusters(); iCluster++)
     {
      AliESDCaloCluster *clust = fESD->GetCaloCluster(iCluster);
      if(clust->IsEMCAL())
        {
         double clustE = clust->E();
         float  emcx[3]; // cluster pos
         clust->GetPosition(emcx);
         TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
         double emcphi = clustpos.Phi(); 
         double emceta = clustpos.Eta();
         double calInfo[5];
         calInfo[0] = emcphi; calInfo[1] = emceta; calInfo[2] = clustE; calInfo[3] = cent; calInfo[4] = clust->Chi2(); 
         if(clustE>1.5)fEMCAccE->Fill(calInfo); 
        }
   }

  // Track loop 
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
   
    Bool_t mcPho = kFALSE;
    Bool_t mcDtoE= kFALSE;
    Bool_t mcBtoE= kFALSE;
    double mcele = -1.;
    if(fmcData && fMC && stack)
      {
       Int_t label = track->GetLabel();
       TParticle* particle = stack->Particle(label);
       int mcpid = particle->GetPdgCode();
       printf("MCpid = %d",mcpid);
       if(particle->GetFirstMother()>-1)
         {
          int parentPID = stack->Particle(particle->GetFirstMother())->GetPdgCode();
          if((parentPID==22 || parentPID==111 || parentPID==221)&& fabs(mcpid)==11)mcPho = kTRUE;
          if((fabs(parentPID)==411 || fabs(parentPID)==413 || fabs(parentPID)==421 || fabs(parentPID)==423 || fabs(parentPID)==431)&& fabs(mcpid)==11)mcDtoE = kTRUE;
          if((fabs(parentPID)==511 || fabs(parentPID)==513 || fabs(parentPID)==521 || fabs(parentPID)==523 || fabs(parentPID)==531)&& fabs(mcpid)==11)mcBtoE = kTRUE;
         } 
       if(fabs(mcpid)==11)mcele= 1.; 
      }
 
    if(TMath::Abs(track->Eta())>0.7) continue;
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

    
    fTrackPtAftTrkCuts->Fill(track->Pt());		
    
    Double_t mom = -999., eop=-999., pt = -999., dEdx=-999., fTPCnSigma=-10, phi=-999., eta=-999.;
    //pt = track->Pt();
    //if(pt<2.0)continue;
    
    // Track extrapolation
    
    Int_t charge = track->Charge();
    fTrkpt->Fill(pt);
    mom = track->P();
    phi = track->Phi();
    eta = track->Eta();
    dEdx = track->GetTPCsignal();
    fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
    
        double ncells = -1.0;
        double m20 = -1.0;
        double m02 = -1.0;
        double disp = -1.0;
        double rmatch = -1.0;  
        double nmatch = -1.0;
        double oppstatus = 0.0;
        double samestatus = 0.0;

    Bool_t fFlagPhotonicElec = kFALSE;
    Bool_t fFlagConvinatElec = kFALSE;
    SelectPhotonicElectron(iTracks,cent,track,fFlagPhotonicElec,fFlagConvinatElec);
    if(fFlagPhotonicElec)oppstatus = 1.0;
    if(fFlagConvinatElec)samestatus = 1.0;

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

		  double valdedx[17];
		  valdedx[0] = pt; valdedx[1] = dEdx; valdedx[2] = phi; valdedx[3] = eta; valdedx[4] = fTPCnSigma;
		  valdedx[5] = eop; valdedx[6] = rmatch; valdedx[7] = ncells,  valdedx[8] = m02; valdedx[9] = m20; valdedx[10] = disp;
		  valdedx[11] = cent; valdedx[12] = charge; valdedx[13] = oppstatus; valdedx[14] = samestatus; valdedx[15] = clust->Chi2();
                  valdedx[16] = mcele;
                  fEleInfo->Fill(valdedx);
                 

      }
    }
        
    fdEdxBef->Fill(mom,fTPCnSigma);
    fTPCnsigma->Fill(mom,fTPCnSigma);
    if(fTPCnSigma >= -1.0 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,eop);

    Int_t pidpassed = 1;
    

    //--- track accepted
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    int centf = (int)cent;
    hfetrack.SetCentrality(centf); //added
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
    if(mcele==1)
      {
       if(mcBtoE || mcDtoE)
         {
          fIncpTMChfe->Fill(cent,pt);    
          if(m20>0.0 && m20<0.3)fIncpTMCM20hfe->Fill(cent,pt);    
         }
       if(mcPho)
        {
         fIncpTMCpho->Fill(cent,pt);    
         if(fFlagPhotonicElec) fPhoElecPtMC->Fill(cent,pt);
         if(fFlagConvinatElec) fSameElecPtMC->Fill(cent,pt);

         if(m20>0.0 && m20<0.3) 
           {
            fIncpTMCM20pho->Fill(cent,pt);    
            if(fFlagPhotonicElec) fPhoElecPtMCM20->Fill(cent,pt);
            if(fFlagConvinatElec) fSameElecPtMCM20->Fill(cent,pt);
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
  
  fPID->SortDetectors(); 
  fPIDqa = new AliHFEpidQAmanager();
  fPIDqa->Initialize(fPID);
  
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
  
  Int_t binsE[5] =    {250, 100,  1000, 100,   10};
  Double_t xminE[5] = {1.0,  -1,   0.0,   0, -0.5}; 
  Double_t xmaxE[5] = {3.5,   1, 100.0, 100,  9.5}; 
  fEMCAccE = new THnSparseD("fEMCAccE","EMC acceptance & E;#eta;#phi;Energy;Centrality;trugCondition;",5,binsE,xminE,xmaxE);
  fOutputList->Add(fEMCAccE);

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
  
  fIncpT = new TH2F("fIncpT","HFE pid electro vs. centrality",100,0,100,100,0,50);
  fOutputList->Add(fIncpT);

  fIncpTM20 = new TH2F("fIncpTM20","HFE pid electro vs. centrality with M20",100,0,100,100,0,50);
  fOutputList->Add(fIncpTM20);

  Int_t nBinspho[3] =  { 100, 100, 500};
  Double_t minpho[3] = {  0.,  0., 0.};   
  Double_t maxpho[3] = {100., 50., 0.5};   

  fInvmassLS = new THnSparseD("fInvmassLS", "Inv mass of LS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2);", 3, nBinspho,minpho, maxpho);
  fOutputList->Add(fInvmassLS);
  
  fInvmassULS = new THnSparseD("fInvmassULS", "Inv mass of ULS (e,e); cent; p_{T} (GeV/c); mass(GeV/c^2);", 3, nBinspho,minpho, maxpho);
  fOutputList->Add(fInvmassULS);
  
  fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleLS);
  
  fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
  fOutputList->Add(fOpeningAngleULS);
  
  fPhotoElecPt = new TH1F("fPhotoElecPt", "photonic electron pt",100,0,50);
  fOutputList->Add(fPhotoElecPt);
  
  fPhoElecPt = new TH2F("fPhoElecPt", "Pho-inclusive electron pt",100,0,100,100,0,50);
  fOutputList->Add(fPhoElecPt);
  
  fPhoElecPtM20 = new TH2F("fPhoElecPtM20", "Pho-inclusive electron pt with M20",100,0,100,100,0,50);
  fOutputList->Add(fPhoElecPtM20);

  fSameElecPt = new TH2F("fSameElecPt", "Same-inclusive electron pt",100,0,100,100,0,50);
  fOutputList->Add(fSameElecPt);

  fSameElecPtM20 = new TH2F("fSameElecPtM20", "Same-inclusive electron pt with M20",100,0,100,100,0,50);
  fOutputList->Add(fSameElecPtM20);

  fCent = new TH1F("fCent","Centrality",100,0,100) ;
  fOutputList->Add(fCent);
 
  // Make common binning
  const Double_t kMinP = 2.;
  const Double_t kMaxP = 50.;
  const Double_t kTPCSigMim = 40.;
  const Double_t kTPCSigMax = 140.;

  // 1st histogram: TPC dEdx with/without EMCAL (p, pT, TPC Signal, phi, eta,  Sig,  e/p,  ,match, cell, M02, M20, Disp, Centrality, select)
  Int_t nBins[17] =  {  480,        200,   60,    20,   600,  300, 100,   40,   200, 200, 200, 100,    3,    3,    3,   10,    3};
  Double_t min[17] = {kMinP,  kTPCSigMim, 1.0,  -1.0,  -8.0,    0,   0,    0,   0.0, 0.0, 0.0,   0, -1.5, -0.5, -0.5, -0.5, -1.5};
  Double_t max[17] = {kMaxP,  kTPCSigMax, 4.0,   1.0,   4.0,  3.0, 0.1,   40,   2.0, 2.0, 2.0, 100,  1.5,  2.5,  2.5,  9.5,  1.5};
  fEleInfo = new THnSparseD("fEleInfo", "Electron Info; pT [GeV/c]; TPC signal;phi;eta;nSig; E/p;Rmatch;Ncell;M02;M20;Disp;Centrality;charge;opp;same;trigCond;", 17, nBins, min, max);
  fOutputList->Add(fEleInfo);

  //<---  Trigger info

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

  //<---------- MC

  fInputHFEMC = new TH2F("fInputHFEMC","Input MC HFE pid electro vs. centrality",100,0,100,100,0,50);
  fOutputList->Add(fInputHFEMC);

  fIncpTMChfe = new TH2F("fIncpTMChfe","MC HFE pid electro vs. centrality",100,0,100,100,0,50);
  fOutputList->Add(fIncpTMChfe);

  fIncpTMCM20hfe = new TH2F("fIncpTMCM20hfe","MC HFE pid electro vs. centrality with M20",100,0,100,100,0,50);
  fOutputList->Add(fIncpTMCM20hfe);

  fIncpTMCpho = new TH2F("fIncpTMCpho","MC Pho pid electro vs. centrality",100,0,100,100,0,50);
  fOutputList->Add(fIncpTMCpho);

  fIncpTMCM20pho = new TH2F("fIncpTMCM20pho","MC Pho pid electro vs. centrality with M20",100,0,100,100,0,50);
  fOutputList->Add(fIncpTMCM20pho);

  fPhoElecPtMC = new TH2F("fPhoElecPtMC", "MC Pho-inclusive electron pt",100,0,100,100,0,50);
  fOutputList->Add(fPhoElecPtMC);
  
  fPhoElecPtMCM20 = new TH2F("fPhoElecPtMCM20", "MC Pho-inclusive electron pt with M20",100,0,100,100,0,50);
  fOutputList->Add(fPhoElecPtM20);

  fSameElecPtMC = new TH2F("fSameElecPtMC", "MC Same-inclusive electron pt",100,0,100,100,0,50);
  fOutputList->Add(fSameElecPtMC);

  fSameElecPtMCM20 = new TH2F("fSameElecPtMCM20", "MC Same-inclusive electron pt with M20",100,0,100,100,0,50);
  fOutputList->Add(fSameElecPtMCM20);


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
void AliAnalysisTaskHFECal::SelectPhotonicElectron(Int_t itrack, Double_t cent, AliESDtrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec)
{
  //Identify non-heavy flavour electrons using Invariant mass method
  
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  //fTrackCuts->SetRequireTPCRefit(kTRUE);
  //fTrackCuts->SetEtaRange(-0.7,0.7);
  //fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fTrackCuts->SetMinNClustersTPC(70);
  
  const AliESDVertex *pVtx = fESD->GetPrimaryVertex();
  
  Bool_t flagPhotonicElec = kFALSE;
  Bool_t flagConvinatElec = kFALSE;
  
  //for(Int_t jTracks = itrack+1; jTracks<fESD->GetNumberOfTracks(); jTracks++){
  for(Int_t jTracks = 0; jTracks<fESD->GetNumberOfTracks(); jTracks++){
    AliESDtrack* trackAsso = fESD->GetTrack(jTracks);
    if (!trackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }
    if(itrack==jTracks)continue;    

    Double_t dEdxAsso = -999., ptPrim=-999., ptAsso=-999., openingAngle = -999.;
    Double_t mass=999., width = -999;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    ptPrim = track->Pt();

    dEdxAsso = trackAsso->GetTPCsignal();
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();
    
    //if(ptAsso <0.3) continue;
    if(ptAsso <0.5) continue;
    if(!fTrackCuts->AcceptTrack(trackAsso)) continue;
    if(dEdxAsso <65 || dEdxAsso>100) continue; //11a pass1
    
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
    
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;
    
    AliKFParticle ge1(*track, fPDGe1);
    AliKFParticle ge2(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);
    
    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
    
    AliKFVertex primV(*pVtx);
    primV += recg;
    recg.SetProductionVertex(primV);
    
    recg.SetMassConstraint(0,0.0001);
    
    openingAngle = ge1.GetAngle(ge2);
    if(fFlagLS) fOpeningAngleLS->Fill(openingAngle);
    if(fFlagULS) fOpeningAngleULS->Fill(openingAngle);
    
    if(openingAngle > fOpeningAngleCut) continue;
    
    recg.GetMass(mass,width);
    
    double phoinfo[3];
    phoinfo[0] = cent;
    phoinfo[1] = ptPrim;
    phoinfo[2] = mass;

    if(fFlagLS) fInvmassLS->Fill(phoinfo);
    if(fFlagULS) fInvmassULS->Fill(phoinfo);
          
    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
    if(mass<fInvmassCut && fFlagLS && !flagConvinatElec){
      flagConvinatElec = kTRUE;
    }
    
  }
  fFlagPhotonicElec = flagPhotonicElec;
  fFlagConvinatElec = flagConvinatElec;
  
}


//_________________________________________
void AliAnalysisTaskHFECal::FindTriggerClusters()
{
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

      Int_t ts = 0;
      fCaloTrigger->GetL1TimeSum(ts);
      if (ts > 0)ftriggers[globCol][globRow] = 1;
      // number of triggered channels in event
      nTrigChannel++;
      // ... inside cut
      if(ts>0 && (bit >> 6 & 0x1))
      {
        iglobCol = globCol;
        iglobRow = globRow;
        nTrigChannelCut++;
        ftriggersCut[globCol][globRow] = 1;
      }

    } // calo trigger entries
  } // has calo trigger entries

  // part 2 go through the clusters here -----------------------------------
  Int_t nCluster=0, nCell=0, iCell=0, gCell=0;
  Short_t cellAddr, nSACell, mclabel;
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
     if(clsE>6.0)fCellCheck->Fill(clsE,cellAddrs[iCell]); 

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

      // try to match with a triggered
      if( ftriggers[feta][fphi]==1)
      {  nClusterTrig++;
      }
      if( ftriggersCut[feta][fphi]==1)
      { nClusterTrigCut++;
      }

    } // cells


    if( gCell > 0 ) 
      clusterTime = clusterTime / (Double_t)gCell;
    // fix the reconstriction code time 100ns jumps
    if( fESD->GetBunchCrossNumber() % 4 < 2 )
      clusterTime -= 0.0000001;

    fClsETime->Fill(clsE,clusterTime);
    fClsEBftTrigCut->Fill(clsE);

    if(nClusterTrig>0){
      fClsETime1->Fill(clsE,clusterTime);
    }

    if(nClusterTrig>0){
      cluster->SetChi2(1);
      fClsEAftTrigCut1->Fill(clsE);                                               
    }

    if(nClusterTrigCut>0){
      cluster->SetChi2(2);
      fClsEAftTrigCut2->Fill(clsE);
    }

    if(nClusterTrigCut>0 && ( nCell > (1 + clsE / 3)))
    {
      cluster->SetChi2(3);
      fClsEAftTrigCut3->Fill(clsE);
    }

    if(nClusterTrigCut>0 && (nCell > (1 + clsE / 3) )&&( clusterTime > fTimeCutLow && clusterTime < fTimeCutHigh ))
    {
      // cluster->SetChi2(4);
      fClsEAftTrigCut4->Fill(clsE);
    }
    if(nClusterTrigCut<1)
    {
      cluster->SetChi2(0);

      fClsEAftTrigCut->Fill(clsE);
    }

  } // clusters
}







