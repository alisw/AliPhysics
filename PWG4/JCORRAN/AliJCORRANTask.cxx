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

//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TChain.h>
#include <TVector3.h> 
#include <TFile.h>
#include <TH1.h> 
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TFormula.h>
#include <TString.h>
#include <TRefArray.h>
#include <TNtuple.h>

#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"

#include "AliJCORRANTask.h" 
#include "AliAnalysisManager.h"
#include "AliESDEvent.h" 
#include "AliMCEvent.h" 
#include "AliStack.h" 
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliESDCaloCluster.h" 
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeoUtils.h"
#include "AliESDtrackCuts.h"

#include "AliPhJTrackList.h"
#include "AliPhJMCTrackList.h"
#include "AliPhJPhotonList.h"
#include "AliPhJHeaderList.h"
#include "AliJTrack.h"
#include "AliJPhoton.h"
#include "AliJHeader.h"
#include "AliAODTracklets.h"
#include "AliMultiplicity.h"
#include "JConst.h"
#include "AliESDRun.h"
#include "AliJRunHeader.h"


//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask() :   
  AliAnalysisTaskSE("PWG4JCORRAN"),
  fInputFormat(0),
  fEsdTrackCuts(0), 
  fDownscaling(1),
  fLowerCutOnLPMom(0),
  fLowerCutOnLeadingCaloClusterE(0), 
  fLowerCutOnCaloClusterE(0.2),
  fIsRealOrMC(0),
  fAODName("jcorran.root"),
  fTrackList(0x0),
  fMCTrackList(0x0),
  fPhotonList(0x0),
  fHeaderList(0x0),
  fAliRunHeader(0x0),
  fPHOSGeom(0x0),
  fEMCALGeom(0x0)
{
  //Default constructor
  fInputFormat = "ESD";

  for(Int_t i=0;i<kRangeTriggerTableAlice;i++)   fActiveTriggers[i]=" ";
  for(Int_t i=0;i<kRangeTriggerTableJCorran;i++) fTriggerTableJCorran[i]=" ";
  
  DefineInput (0, TChain::Class());
  
}

//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const char *name, TString inputformat) : 
  AliAnalysisTaskSE(name), 
  fInputFormat(inputformat),  
  fEsdTrackCuts(0),    // to be set by setters in AddAliJCORRANTask macro
  fDownscaling(1),
  fLowerCutOnLPMom(0),
  fLowerCutOnLeadingCaloClusterE(0),
  fLowerCutOnCaloClusterE(0.2),
  fIsRealOrMC(0),
  fAODName("jcorran.root"),
  fTrackList(0x0),
  fMCTrackList(0x0),
  fPhotonList(0x0),
  fHeaderList(0x0),
  fAliRunHeader(0x0),
  fPHOSGeom(0x0),
  fEMCALGeom(0x0)
{
  // Constructor.
  if(fDebug > 5) cout << "---- AliJCORRANTask Constructor ----"<<endl;

  for(Int_t i=0;i<kRangeTriggerTableAlice;i++)   fActiveTriggers[i]=" ";
  for(Int_t i=0;i<kRangeTriggerTableJCorran;i++) fTriggerTableJCorran[i]=" ";
  
  DefineInput (0, TChain::Class());

 

}

//____________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const AliJCORRANTask& ap) :
  AliAnalysisTaskSE(ap.GetName()), 
  fInputFormat(ap.fInputFormat),
  fEsdTrackCuts(ap.fEsdTrackCuts), 
  fDownscaling(ap.fDownscaling),   
  fLowerCutOnLPMom(ap.fLowerCutOnLPMom),  
  fLowerCutOnLeadingCaloClusterE(ap.fLowerCutOnLeadingCaloClusterE),  
  fLowerCutOnCaloClusterE(ap.fLowerCutOnCaloClusterE),
  fIsRealOrMC(ap.fIsRealOrMC),
  fAODName(ap.fAODName),
  fTrackList(ap.fTrackList),
  fMCTrackList(ap.fMCTrackList),
  fPhotonList(ap.fPhotonList),
  fHeaderList(ap.fHeaderList),
  fAliRunHeader(ap.fAliRunHeader),
  fPHOSGeom(ap.fPHOSGeom),
  fEMCALGeom(ap.fEMCALGeom)
{ 
  // cpy ctor
}

//_____________________________________________________________________________
AliJCORRANTask& AliJCORRANTask::operator= (const AliJCORRANTask& ap)
{
// assignment operator

  this->~AliJCORRANTask();
  new(this) AliJCORRANTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJCORRANTask::~AliJCORRANTask()
{
  // destructor 
  
  delete fTrackList;
  if(fMCTrackList) delete fMCTrackList;
  delete fPhotonList;
  delete fHeaderList;
  delete fAliRunHeader;
  if(fPHOSGeom)  delete fPHOSGeom;
  if(fEMCALGeom) delete fEMCALGeom;

}

//________________________________________________________________________
void AliJCORRANTask::UserCreateOutputObjects()
{  

  fTrackList    = new AliPhJTrackList(kALICE);
  fMCTrackList  = new AliPhJMCTrackList(kALICE);
  fPhotonList   = new AliPhJPhotonList(kALICE);
  fHeaderList   = new AliPhJHeaderList(kALICE);
  
  fAliRunHeader = new AliJRunHeader();

  fPHOSGeom  = new AliPHOSGeoUtils("PHOSgeo") ;
  fEMCALGeom = new AliEMCALGeoUtils("EMCAL_COMPLETE");
 
  // create the jcorran output deltaAOD
  //if(fDebug > 5) cout << "AliJCORRANTask UserCreateOutputObjects----------------------"<<endl;
  
  if(fDebug > 1) printf("AliJCORRANTask::UserCreateOutPutData() \n");
  if(!AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }   


  AddAODBranch("AliPhJTrackList", &fTrackList,fAODName.Data());
  AddAODBranch("AliPhJPhotonList", &fPhotonList,fAODName.Data());
  AddAODBranch("AliPhJHeaderList", &fHeaderList,fAODName.Data()); 

  if(fIsRealOrMC){ 
    AddAODBranch("AliPhJMCTrackList", &fMCTrackList);
  }
  
}

//______________________________________________________________________________
void AliJCORRANTask::UserExec(Option_t */*option*/) 
{
  // Processing of one event
  //if(fDebug > 5) cout << "------- AliJCORRANTask Exec-------"<<endl;
  if(!((Entry()-1)%100)) 
      AliInfo(Form(" Processing event # %lld",  Entry())); 
  Bool_t storeEvent = kFALSE;//based on offline trigger decide whetehr to store the event or not 
  if(fIsRealOrMC){
    storeEvent = kTRUE; //store all MC events
  }else{ //when we are processing real events store only selected events
    if(StoreDownscaledMinBiasEvent() || 
      ContainsESDHighPtTrack()   || 
      ContainsESDHighECaloClusters()){
        storeEvent = kTRUE; 
    }
  }

  fTrackList->Reset();
  fMCTrackList->Reset();
  fPhotonList->Reset();
  fHeaderList->Reset();
 

  static int runId=-1;
 
  if(fInputFormat=="ESD"){
   //  if(fDebug > 5) cout <<"--------- Reading ESD --------"<< endl; 
     AliESDEvent* esd = (AliESDEvent*)InputEvent();

     AliMCEvent* mcEvent = NULL; 
     if(fIsRealOrMC)  mcEvent = MCEvent();


     //=========== FILL AND STORE RUN HEADER   (ONLY ONCE PER RUN) =============
     if(esd->GetRunNumber() != runId){ //new run has started
 
       const AliESDRun* esdRun = esd->GetESDRun();

       for(Int_t triggerBit=0; triggerBit<kRangeTriggerTableAlice; triggerBit++){
         fActiveTriggers[triggerBit] = esdRun->GetTriggerClass(triggerBit);
       }
       runId = esd->GetRunNumber();

       //Polarity of magnetic field in L3 solenoid
       Short_t l3MgFieldPolarity=0; // (LHC convention: +current -> +Bz)
       if(esd->GetCurrentL3() >0) l3MgFieldPolarity =  1;
       if(esd->GetCurrentL3() <0) l3MgFieldPolarity = -1;

       //Create JCorran trigger mask mapping
       fTriggerTableJCorran[kMinBiasTriggerBitJCorran]="Minimum Bias";//minimum bias occupies first trigger bit
      
       //=========== Fill Run header object ===============
       fAliRunHeader->SetRunNumber(runId);
       fAliRunHeader->SetActiveTriggersAlice(fActiveTriggers);
       fAliRunHeader->SetL3Field(l3MgFieldPolarity, esd->GetMagneticField());
       fAliRunHeader->SetActiveTriggersJCorran(fTriggerTableJCorran,kRangeTriggerTableJCorran);

       //Store Run header
       (OutputTree()->GetUserInfo())->Add(fAliRunHeader);  //FK// 
    }

    if(storeEvent){ 
      //-------------- reset all the arrays -------------
         //store event only when it is downscaled min bias
      // or contais high pt hadron
      // or contains high energy cluster in EMCAL or PHOS
      ReadESDTracks(esd);
      ReadESDCaloClusters(esd);
      ReadESDHeader(esd);
      if(fIsRealOrMC) ReadMCTracks(mcEvent);
    }
  }else if( fInputFormat == "AODout" || fInputFormat == "AODin") {
  
    AliAODEvent* aod = NULL;
    if(fInputFormat == "AODout"){ // reading from AOD output handler
      aod = AODEvent();
    }else if(fInputFormat == "AODin"){ // reading from AOD input handler
      aod = (AliAODEvent*)InputEvent();
    }

    if(storeEvent){ 
      //-------------- reset all the arrays -------------
         
      ReadAODTracks(aod);
      ReadAODCaloClusters(aod);
      ReadAODHeader(aod);
    }
    
  }else{
    cout << "Error: Not correct InputDataFormat especified " << endl;
    return;
  }
}

//______________________________________________________________________________
void AliJCORRANTask::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJCORRANTask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  OutputTree()->Print(); 
  if(fInputFormat == "AODout") ReadFilter(); // change it to save this info also from AODin !!!! 

  ((AliJRunHeader *) (OutputTree()->GetUserInfo())->First())->PrintOut();  

  cout<<"PWG4JCORRAN Analysis DONE !!"<<endl; 


}

//______________________________________________________________________________
void AliJCORRANTask::ReadESDTracks(const AliESDEvent * esd)
{

    // Read the AliESDtrack and fill the list of AliJTrack containers
    Int_t nt = esd->GetNumberOfTracks();
   // if(fDebug > 5) cout << "ESD::NumberOfTracks = " << nt << endl;
    Float_t ptot, pt, eta;
    TVector3 p3(0,0,0);
    Short_t ntrk = 0;
    Double_t pid[10];

    //loop over tracks
    for(Int_t it = 0; it < nt; it++) { 

        AliESDtrack *track = esd->GetTrack(it);
        if(! fEsdTrackCuts->IsSelected(track)) continue; //apply quality selection criteria

        UInt_t status = track->GetStatus();
	    
	Float_t impactXY, impactZ;
	track->GetImpactParameters(impactXY,impactZ);
	//------------ T P C ------------
        Float_t tpcDca[2], tpcCov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
	track->GetImpactParametersTPC(tpcDca,tpcCov);

        Float_t tpcDedx = track->GetTPCsignal();
        Int_t   tpcNcls = track->GetTPCNcls();

	Int_t nClust         = track->GetTPCclusters(0);
	Int_t nFindableClust = track->GetTPCNclsF();
	Float_t tpcChi2PerCluster = 0.;
	if(nClust>0.) tpcChi2PerCluster = track->GetTPCchi2()/Float_t(nClust);
	Float_t tpcClustPerFindClust = 0.;
	if(nFindableClust>0.) tpcClustPerFindClust = Float_t(nClust)/nFindableClust;
        //--------------------------------

        Double_t mom[3];
        track->GetPxPyPz(mom);
        p3.SetXYZ(mom[0],mom[1],mom[2]);
        ptot = track->GetP();
        pt = p3.Pt();
        eta = p3.Eta();
        track->GetESDpid(pid);
       
        Double_t extCov[15];
        track->GetExternalCovariance(extCov);

        Double_t extDiaCov[5];
        extDiaCov[0]=extCov[0];
        extDiaCov[1]=extCov[2];
        extDiaCov[2]=extCov[5];
        extDiaCov[3]=extCov[9];
        extDiaCov[4]=extCov[14];

      //  Int_t itsLabel = track->GetITSLabel(); //FK//
      //  Int_t tpcLabel = track->GetTPCLabel(); //FK//   

        //create a new AliJTrack and fill the track info
	fTrackList->AddAliJTrack(ntrk);
	AliJTrack *ctrack = fTrackList->GetAliJTrack(ntrk);

        ctrack->SetPtot(ptot);
        ctrack->SetPt(pt);
        ctrack->SetTheta(p3.Theta());
        ctrack->SetPhi(p3.Phi());
        ctrack->SetPID(pid);
        ctrack->SetFlavor(kNone);//kHadron);
	ctrack->SetCharge(track->Charge());
        ctrack->ConvertAliPID();
        ctrack->SetEta(eta);

	Int_t itof = track->GetTOFcluster(); // index of the assigned TOF cluster
	if(itof>=0) ctrack->SetTof(track->GetTOFsignal());

	ctrack->SetTPCdEdx(tpcDedx);
	ctrack->SetTPCnClust(tpcNcls);
	ctrack->SetTPCDCAXY(tpcDca[0]);
	ctrack->SetTPCDCAZ(tpcDca[1]);
	ctrack->SetTPCClustPerFindClust(tpcClustPerFindClust);
	ctrack->SetTPCChi2PerClust(tpcChi2PerCluster);
	ctrack->SetImapactXY(impactXY);
	ctrack->SetImapactZ(impactZ);

        ctrack->SetKinkIndex(track->GetKinkIndex(0));
        ctrack->SetStatus(status);
        ctrack->SetExternalDiaCovariance(extDiaCov);

      //  ctrack->SetITSLabel(itsLabel);//FK//
      //  ctrack->SetTPCLabel(tpcLabel);//FK//


	fTrackList->SetNTracks(++ntrk);

     } // end tracks loop
}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODTracks(const AliAODEvent * aod)
{
    // Read the AliAODtrack and fill the list of AliJTrack containers
    Int_t nt = aod->GetNumberOfTracks();
    if(fDebug < 5) cout << "AOD::NumberOfTracks = " << nt << endl;
    Short_t ntrk = 0;
    //loop over tracks
    for(Int_t it = 0; it < nt; it++) { 

        AliAODTrack *track = aod->GetTrack(it);
        if(track->GetType() != AliAODTrack::kPrimary) continue; // only primaries 
        if(! fEsdTrackCuts->IsSelected(track)) continue; //apply loose selection criteria
        //create a new AliJTrack and fill the track info
	fTrackList->AddAliJTrack(ntrk);
        AliJTrack *ctrack = fTrackList->GetAliJTrack(ntrk);

        ctrack->SetPtot(track->P());
        ctrack->SetPt(track->Pt());
        ctrack->SetTheta(track->Theta());
        ctrack->SetPhi(track->Phi());
        ctrack->SetEta(track->Eta());
        ctrack->SetPID((Double_t*)track->PID());
        ctrack->SetFlavor(kNone); //kHadron);
        ctrack->SetCharge(track->Charge());
        ctrack->SetChi2perNDF(track->Chi2perNDF());
        ctrack->SetChi2Trig(track->GetChi2MatchTrigger());
        ctrack->SetRecFlags(track->GetFlags());
      
       
       // ctrack->SetITSLabel(track->GetLabel());//FK//?
      //  ctrack->SetTPCLabel(track->GetLabel());//FK//?
	fTrackList->SetNTracks(++ntrk);

     } // end tracks loop
}

//______________________________________________________________________________
void AliJCORRANTask::ReadESDCaloClusters(const AliESDEvent* esd)
{
  // Read the AliESDCaloClusters and fill the list of AliJPhoton containers
  Short_t nPhotons = 0;
  Int_t numberOfCaloClusters = esd->GetNumberOfCaloClusters() ;
  if(fDebug < 5) cout << "ESD::number of ESD caloclusters " << numberOfCaloClusters << endl;
  // loop over all the Calo Clusters
  for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {
    AliESDCaloCluster *caloCluster = esd->GetCaloCluster(icluster) ;
    if(!caloCluster) continue;
    if(caloCluster->GetTrackMatched()==-1){
      if(caloCluster->E()<fLowerCutOnCaloClusterE) continue;                  //FK//
      // we will not implement any PID cut here      
      fPhotonList->AddAliJPhoton(nPhotons);
      AliJPhoton *pht = fPhotonList->GetAliJPhoton(nPhotons);
      pht->SetFlavor(kPhoton);
      pht->SetE(caloCluster->E());
      pht->SetChi2(caloCluster->GetClusterChi2());
      pht->SetPID(caloCluster->GetPid());
      Float_t pos[3]; caloCluster->GetPosition(pos) ;
      pht->SetX(pos[0]);
      pht->SetY(pos[1]);
      pht->SetZ(pos[2]);
      pht->SetPhi(atan2(pos[1],pos[0]));
      pht->SetTheta(atan2(sqrt(pos[0]*pos[0]+pos[1]*pos[1]),pos[2]));
      pht->SetPt(caloCluster->E()*sin(atan2(sqrt(pos[0]*pos[0]+pos[1]*pos[1]),pos[2])));
      pht->SetCharge(0);
      if(caloCluster->IsEMCAL()) pht->SetCaloType(AliJPhoton::kEMCALCalo);
      if(caloCluster->IsPHOS())  pht->SetCaloType(AliJPhoton::kPHOSCalo);
      pht->SetDistToBadChannel(caloCluster->GetDistanceToBadChannel());
      pht->SetDispersion(caloCluster->GetClusterDisp());
      pht->SetM20(caloCluster->GetM20());
      pht->SetM02(caloCluster->GetM02());
      pht->SetEmcCpvDist(caloCluster->GetEmcCpvDistance());
      pht->SetNCells(caloCluster->GetNCells());
      pht->SetCellsAbsId(caloCluster->GetCellsAbsId());
      Int_t imoduleID = GetSuperModuleNumber(caloCluster->IsEMCAL(), caloCluster->GetCellAbsId(0));
      pht->SetSuperModuleID(imoduleID);
      
      fPhotonList->SetNPhotons(++nPhotons);
    } // end if 
  } //PHOS and EMCAL clusters

}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODCaloClusters(const AliAODEvent* aod)
{
  // Read the AliAODCaloClusters and fill the list of AliJPhoton containers
  Int_t numberOfCaloClusters = aod->GetNCaloClusters() ;
  if(fDebug < 5) cout << "AOD::number of ESD caloclusters " << numberOfCaloClusters << endl;
  Short_t nPhotons = 0;
  for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {
       AliAODCaloCluster * caloCluster = aod->GetCaloCluster(icluster) ;
       if(!caloCluster) continue; 
       if(caloCluster->GetNTracksMatched() > 0) continue;
      // we will not implement any PID cut here      
      fPhotonList->AddAliJPhoton(nPhotons);
      AliJPhoton *pht = fPhotonList->GetAliJPhoton(nPhotons);
      
      pht->SetE(caloCluster->E());
      pht->SetFlavor(kPhoton);
      pht->SetChi2(caloCluster->Chi2());
      pht->SetPID((Double_t*)caloCluster->PID());
      Float_t pos[3]; caloCluster->GetPosition(pos) ;
      pht->SetX(pos[0]);
      pht->SetY(pos[1]);
      pht->SetZ(pos[2]);
      pht->SetPhi(atan2(pos[1],pos[0]));
      pht->SetTheta(atan2(sqrt(pos[0]*pos[1]+pos[1]*pos[1]),pos[2]));
      pht->SetPt(caloCluster->E()*sin(atan2(sqrt(pos[0]*pos[0]+pos[1]*pos[1]),pos[2])));
      pht->SetCharge(0);
      if(caloCluster->IsEMCALCluster()) pht->SetCaloType(AliJPhoton::kEMCALCalo);
      if(caloCluster->IsPHOSCluster())  pht->SetCaloType(AliJPhoton::kPHOSCalo);
      pht->SetDistToBadChannel(caloCluster->GetDistToBadChannel());
      pht->SetDispersion(caloCluster->GetDispersion());
      pht->SetM20(caloCluster->GetM20());
      pht->SetM02(caloCluster->GetM02());
      pht->SetEmcCpvDist(caloCluster->GetEmcCpvDistance());
      pht->SetNCells(int(caloCluster->GetNCells()));
      pht->SetCellsAbsId(caloCluster->GetCellsAbsId());
      Int_t imoduleID = GetSuperModuleNumber(caloCluster->IsEMCALCluster(), caloCluster->GetCellAbsId(0));
      pht->SetSuperModuleID(imoduleID);
      
      fPhotonList->SetNPhotons(++nPhotons);
  
  } // clusters
}

//______________________________________________________________________________
void AliJCORRANTask::ReadESDHeader(const AliESDEvent *esd)
{
    // Read the AliESDEvent and fill the list of AliJHeader containers
    Short_t nHeaders = 0;
    //create a header and fill it
    fHeaderList->AddAliJHeader(nHeaders);
    AliJHeader *hdr = fHeaderList->GetAliJHeader(nHeaders);
	
    AliMultiplicity *fSPDMult =(AliMultiplicity *) esd->GetMultiplicity();
    hdr->SetSPDTrackletMult(fSPDMult->GetNumberOfTracklets());

    hdr->SetEventID( esd->GetEventNumberInFile());

    hdr->SetTriggerMaskAlice(esd->GetTriggerMask()); //ULong64_t
    hdr->SetTriggerMaskJCorran(ConvertTriggerMask(/*esd->GetTriggerMask()*/)); //UInt_t

    const AliESDVertex * vtxESD = esd->GetVertex();
    hdr->SetZVertex(vtxESD->GetZv());
    hdr->SetZVertexErr(vtxESD->GetZRes());
    
    fHeaderList->SetNHeaders(++nHeaders);
}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODHeader(const AliAODEvent *aod)
{
  //read AOD event header
  Short_t nHeaders = 0;
  //create a header and fill it
  fHeaderList->AddAliJHeader(nHeaders);
  AliJHeader *hdr = fHeaderList->GetAliJHeader(nHeaders);
	 
  //load aod event header
  AliAODHeader * aodh = aod->GetHeader();

  hdr->SetCentrality(int(aodh->GetCentrality())); 
  
  hdr->SetTriggerMaskAlice(aodh->GetTriggerMask()); //ULong64_t
  hdr->SetTriggerMaskJCorran(ConvertTriggerMask(/*esd->GetTriggerMask()*/)); //UInt_t
  hdr->SetEventType(aodh->GetEventType());
	
  fHeaderList->SetNHeaders(++nHeaders);
}

//______________________________________________________________________________
Int_t AliJCORRANTask::GetSuperModuleNumber(bool isemcal, Int_t absId)
{
  //get super module number 
  if(isemcal){
    return fEMCALGeom->GetSuperModuleNumber(absId) ;
  } else {
    Int_t    relId[4];
    if ( absId >= 0) {
      fPHOSGeom->AbsToRelNumbering(absId,relId);
      return relId[0]-1; 
    } else return -1;
  }//PHOS

  return -1;
}

//_____________________________________________________________________________

UInt_t AliJCORRANTask::ConvertTriggerMask(/*Long64_t alicetriggermask*/){
 //convert alice trigger mask to jcorran trigger mask
  UInt_t triggerMaskJC=0;
 
  Bool_t isSelected = kTRUE;
  isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  if(isSelected){ //tag collision candidates
    triggerMaskJC += (1<<kMinBiasTriggerBitJCorran); 
  }

  return triggerMaskJC;
}


//______________________________________________________________________________
bool AliJCORRANTask::StoreDownscaledMinBiasEvent(){
  
  bool isThisEventToBeStored = kFALSE;
  static Long_t evtNumber=0;

  if( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() ){
    //collision candidate

    if(evtNumber% fDownscaling ==0){ isThisEventToBeStored = kTRUE;  } //store every 50th collision candidate event
    evtNumber++;
  }
  return isThisEventToBeStored;
}

//______________________________________________________________________________
bool AliJCORRANTask::ContainsESDHighPtTrack(){

  bool isThisEventToBeStored = kFALSE; //initialize return value

  if(fInputFormat=="ESD"){

    AliESDEvent* esd = NULL; 
    esd = (AliESDEvent*)InputEvent();

    Int_t nt = esd->GetNumberOfTracks();

    for(Int_t it = 0; it < nt; it++) {
      AliESDtrack *track = esd->GetTrack(it);
      //Does event contain high pt particle above thereshold in GeV 
      if(track->Pt() > fLowerCutOnLPMom && fEsdTrackCuts->IsSelected(track)){
        isThisEventToBeStored = kTRUE;
        break; 
      }
    }
  }else{

    AliAODEvent* aod=NULL;
    if(fInputFormat == "AODout"){ // reading from AOD output handler
      aod = AODEvent();
    }else if(fInputFormat == "AODin"){ // reading from AOD input handler
      aod = (AliAODEvent*)InputEvent();
    }

    Int_t nt = aod->GetNumberOfTracks();

    for(Int_t it = 0; it < nt; it++) {
      AliAODTrack *track = aod->GetTrack(it);
      //Does event contain high pt particle above threshold in GeV 
      if(track->Pt() > fLowerCutOnLPMom && IsSelectedAODTrack(track)){
        isThisEventToBeStored = kTRUE;
        break; 
      }
    }
  }

  return isThisEventToBeStored;
}

//______________________________________________________________________________
bool AliJCORRANTask::ContainsESDHighECaloClusters(){
  bool isThisEventToBeStored = kFALSE; //initialize return value

  if(fInputFormat=="ESD"){

    AliESDEvent* esd = NULL; 
    esd = (AliESDEvent*)InputEvent();

    Int_t numberOfCaloClusters = esd->GetNumberOfCaloClusters() ;
    // loop over all the Calo Clusters
    for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {
      AliESDCaloCluster *caloCluster = esd->GetCaloCluster(icluster) ;
      if(!caloCluster) continue;
      if(caloCluster->GetTrackMatched()==-1){
        //sotre calo clusters above 1 GeV
        if( caloCluster->E() > fLowerCutOnLeadingCaloClusterE){
          isThisEventToBeStored = kTRUE;
          break; 
        }
      }
    }
  }else{

    AliAODEvent* aod=NULL;
    if(fInputFormat == "AODout"){ // reading from AOD output handler
      aod = AODEvent();
    }else if(fInputFormat == "AODin"){ // reading from AOD input handler
      aod = (AliAODEvent*)InputEvent();
    }

    Int_t numberOfCaloClusters = aod->GetNCaloClusters() ;
    // loop over all the Calo Clusters
    for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {
      AliAODCaloCluster *caloCluster = aod->GetCaloCluster(icluster) ;
      if(!caloCluster) continue;
      if(caloCluster->GetNTracksMatched() > 0) continue;
      //sotre calo clusters above 1 GeV
      if( caloCluster->E() > fLowerCutOnLeadingCaloClusterE){
        isThisEventToBeStored = kTRUE;
        break; 
      }
    }
  }

  return isThisEventToBeStored;
}


//______________________________________________________________________________

void AliJCORRANTask::ReadMCTracks(AliMCEvent *fMC)
{
  //AliGenEventHeader* genHeader = fMC->GenEventHeader();
  //AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
         
  //Double_t ptHard = 0;
  //Double_t nTrials = 1; // Trials for MC trigger weigth for real data

  //nTrials = pythiaGenHeader->Trials();
  //ptHard  = pythiaGenHeader->GetPtHard();

  AliStack *stack = fMC->Stack();
  Int_t np        = fMC->GetNumberOfTracks();
  //Int_t nprim = stack->GetNtrack();
  //  if(np!=nprim) cout << "GetNumberOfTracks = "<< np <<"\t, stack = "<< nprim << endl;
 
  Short_t ntrack = 0;

  for(Int_t itrack = 0; itrack < np; itrack++){
    AliMCParticle *track = (AliMCParticle*) fMC->GetTrack(itrack);
    if(!track){
      Printf("ERROR: Could not receive track %d", itrack);
      continue;
    }

    Bool_t isPrimary = stack->IsPhysicalPrimary(itrack);

    //create a new JMCTrack and fill the track info
    fMCTrackList->AddJMCTrack(ntrack);
    AliJMCTrack *ctrack = fMCTrackList->GetTrack(ntrack);

    TParticle *partStack = stack->Particle(itrack);
    Int_t   pdg  = partStack->GetPdgCode();
    Float_t engy = partStack->Energy();
    Float_t pt   = partStack->Pt();
    Float_t ptot = partStack->P();
    Float_t eta  = partStack->Eta();
    Float_t theta  = partStack->Theta();
    Float_t phi    = atan2(sin( partStack->Phi()), cos(partStack->Phi()));
    Short_t ch     = (Short_t) partStack->GetPDG()->Charge();
    Int_t label    = track->GetLabel();
    Int_t   status = partStack->GetStatusCode();

    ctrack->SetLabel(label);
    ctrack->SetPdgCode(pdg);
    ctrack->SetPt(pt);
    ctrack->SetTheta(theta);
    ctrack->SetEta(eta);
    ctrack->SetPhi(phi);
    ctrack->SetE(engy);
    ctrack->SetCharge(ch);
    ctrack->SetPtot(ptot);
    ctrack->SetStatusCode(status);
    ctrack->SetIsPrimary(isPrimary);

    ctrack->SetProductionVertex(partStack->Vx(),partStack->Vy(),partStack->Vz());

    //ctrack->SetPtHard(ptHard);
      
    //bool isInc = (status ==  1 && icode ==  22); //Inclusive
    bool ispi0 = (status == 11 && pdg == 111); //kPizero
    bool isDgamma = (status == 6 || status == 7) && pdg == 22; // Direct photon
    bool inPHOS  = (ispi0||isDgamma)&&fabs(eta)<0.12; 
    bool inEMCAL = (ispi0||isDgamma)&&fabs(eta)<0.7; 
    bool inTPC   = fabs(eta)<0.9; 
    ctrack->SetMother(0,partStack->GetFirstMother());
    ctrack->SetMother(1,partStack->GetSecondMother());
    ctrack->SetDaughter(0,partStack->GetFirstDaughter());
    ctrack->SetDaughter(1,partStack->GetLastDaughter());
    ctrack->SetIsInPHOS(inPHOS);
    ctrack->SetIsInEMCAL(inEMCAL);
    ctrack->SetIsInTPC(inTPC);

    fMCTrackList->SetNTracks(++ntrack);
  }// loop for al primary tracks
}

//______________________________________________________________________________

bool AliJCORRANTask::IsSelectedAODTrack(AliAODTrack   *track){

  if(fIsRealOrMC &&  track->GetType() != AliAODTrack::kPrimary) return kFALSE; // only primaries 
  if(fEsdTrackCuts->GetMinNClusterTPC() > track->GetTPCNcls()) return kFALSE;
  if(fEsdTrackCuts->GetRequireTPCRefit()  && ((track->GetStatus() & AliJTrack::kTPCrefit) == 0)) return kFALSE;
  if(fEsdTrackCuts->GetRequireITSRefit()  && ((track->GetStatus() & AliJTrack::kITSrefit) == 0)) return kFALSE;

  return kTRUE;
}

//______________________________________________________________________________






