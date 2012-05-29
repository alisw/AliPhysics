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

#include "AliCentrality.h"

ClassImp(AliAnalysisTaskHFECal)
//________________________________________________________________________
AliAnalysisTaskHFECal::AliAnalysisTaskHFECal(const char *name) 
  : AliAnalysisTaskSE(name)
  ,fESD(0)
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
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
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fPhotoElecPt(0)
  ,fPhoElecPt(0)
  ,fSameElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fEleInfo(0)
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
  ,fOutputList(0)
  ,fTrackCuts(0)
  ,fCuts(0)
  ,fIdentifiedAsOutInz(kFALSE)
  ,fPassTheEventCut(kFALSE)
  ,fRejectKinkMother(kFALSE)
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
  ,fInvmassLS(0)		
  ,fInvmassULS(0)		
  ,fOpeningAngleLS(0)	
  ,fOpeningAngleULS(0)	
  ,fPhotoElecPt(0)
  ,fPhoElecPt(0)
  ,fSameElecPt(0)
  ,fTrackPtBefTrkCuts(0)	 
  ,fTrackPtAftTrkCuts(0)	 	  
  ,fTPCnsigma(0)
  ,fCent(0)
  ,fEleInfo(0)
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
  
  Float_t cent = -1.;
  AliCentrality *centrality = fESD->GetCentrality(); 
  cent = centrality->GetCentralityPercentile("V0M");
  fCent->Fill(cent);
  
  if(cent>90.) return;
	
 // Calorimeter info.
 
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
         double calInfo[4];
         calInfo[0] = emcphi; calInfo[1] = emceta; calInfo[2] = clustE; calInfo[3] = cent; 
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

		  double valdedx[16];
		  valdedx[0] = mom; valdedx[1] = pt; valdedx[2] = dEdx; valdedx[3] = phi; valdedx[4] = eta; valdedx[5] = fTPCnSigma;
		  valdedx[6] = eop; valdedx[7] = rmatch; valdedx[8] = ncells,  valdedx[9] = m02; valdedx[10] = m20; valdedx[11] = disp;
		  valdedx[12] = cent; valdedx[13] = charge; valdedx[14] = oppstatus; valdedx[15] = samestatus;
                  fEleInfo->Fill(valdedx);
                 

      }
    }
        
    fdEdxBef->Fill(mom,dEdx);
    fTPCnsigma->Fill(mom,fTPCnSigma);
    if(fTPCnSigma >= -1.0 && fTPCnSigma <= 3)fTrkEovPBef->Fill(pt,eop);

    Int_t pidpassed = 1;
    

    //--- track accepted
    AliHFEpidObject hfetrack;
    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    hfetrack.SetRecTrack(track);
    hfetrack.SetPbPb();
    if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;

    if(pidpassed==0) continue;
    
    fTrkEovPAft->Fill(pt,eop);
    fdEdxAft->Fill(mom,dEdx);
    fIncpT->Fill(cent,pt);    

    
    if(fFlagPhotonicElec) fPhoElecPt->Fill(cent,pt);
    if(fFlagConvinatElec) fSameElecPt->Fill(cent,pt);
 }
 PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskHFECal::UserCreateOutputObjects()
{
  //--------Initialize PID
  fPID->SetHasMCData(kFALSE);
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
  
  Int_t binsE[4] =    {250, 100,  1000, 100};
  Double_t xminE[4] = {1.0,  -1,   0.0,   0}; 
  Double_t xmaxE[4] = {3.5,   1, 100.0, 100}; 
  fEMCAccE = new THnSparseD("fEMCAccE","EMC acceptance & E;#eta;#phi;Energy;Centrality",4,binsE,xminE,xmaxE);
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
  
  fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",100,0,50,150,0,150);
  fOutputList->Add(fdEdxBef);
  
  fdEdxAft = new TH2F("fdEdxAft","track dEdx vs p after HFE pid",100,0,50,150,0,150);
  fOutputList->Add(fdEdxAft);
  
  fIncpT = new TH2F("fIncpT","HFE pid electro vs. centrality",100,0,100,100,0,50);
  fOutputList->Add(fIncpT);


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
  
  fSameElecPt = new TH2F("fSameElecPt", "Same-inclusive electron pt",100,0,100,100,0,50);
  fOutputList->Add(fSameElecPt);

  fCent = new TH1F("fCent","Centrality",100,0,100) ;
  fOutputList->Add(fCent);
 
  // Make common binning
  const Double_t kMinP = 2.;
  const Double_t kMaxP = 50.;
  const Double_t kTPCSigMim = 40.;
  const Double_t kTPCSigMax = 140.;

  // 1st histogram: TPC dEdx with/without EMCAL (p, pT, TPC Signal, phi, eta,  Sig,  e/p,  ,match, cell, M02, M20, Disp, Centrality, select)
  Int_t nBins[16] =  {  480,   480,        200,   60,    20,   600,  300, 100,   40,   200, 200, 200, 100,    3,    3,    3};
  Double_t min[16] = {kMinP, kMinP,  kTPCSigMim, 1.0,  -1.0,  -8.0,    0,   0,    0,   0.0, 0.0, 0.0,   0, -1.5, -0.5, -0.5};
  Double_t max[16] = {kMaxP, kMaxP,  kTPCSigMax, 4.0,   1.0,   4.0,  3.0, 0.1,   40,   2.0, 2.0, 2.0, 100,  1.5,  2.5,  2.5};
  fEleInfo = new THnSparseD("fEleInfo", "Electron Info; p [GeV/c]; pT [GeV/c]; TPC signal;phi;eta;nSig; E/p;Rmatch;Ncell;M02;M20;Disp; Centrality; charge", 16, nBins, min, max);
  fOutputList->Add(fEleInfo);

//_________________________________________________________
 
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
  
  //fTrackCuts->SetAcceptKinkDaughters(kFALSE);
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

