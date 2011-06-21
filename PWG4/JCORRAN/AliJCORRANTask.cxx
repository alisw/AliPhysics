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
#include "AliPHOSGeoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliESDtrackCuts.h"
#include "AliAODVertex.h" 
#include "AliAODTracklets.h" 
#include "AliAODPid.h" 

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
  fInputFormat("ESD"),
  fEsdTrackCuts(0x0), 
  fDownscaling(1),
  fLowerCutOnLPMom(0),
  fLowerCutOnLeadingCaloClusterE(0),
  fLowerCutOnCaloClusterE(0.2),
  fIsRealOrMC(0),
  fAODName("jcorran.root"),
  f1CutMaxDCAToVertexXYPtDep(0x0),
  fTrackList(0x0),
  fMCTrackList(0x0),
  fPhotonList(0x0),
  fHeaderList(0x0),
  fAliRunHeader(0x0),
  fPHOSGeom(0x0)
{
  //Default constructor
  for(Int_t i=0;i<kRangeTriggerTableAlice;i++)   fActiveTriggers[i]=" ";
  for(Int_t i=0;i<kRangeTriggerTableJCorran;i++) fTriggerTableJCorran[i]=" ";
  
  DefineInput (0, TChain::Class());
}

//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const char *name, TString inputformat):
  AliAnalysisTaskSE(name), 
  fInputFormat(inputformat),  
  fEsdTrackCuts(0x0), 
  fDownscaling(1),
  fLowerCutOnLPMom(0),
  fLowerCutOnLeadingCaloClusterE(0),
  fLowerCutOnCaloClusterE(0.2),
  fIsRealOrMC(0),
  fAODName("jcorran.root"),
  f1CutMaxDCAToVertexXYPtDep(0x0),
  fTrackList(0x0),
  fMCTrackList(0x0),
  fPhotonList(0x0),
  fHeaderList(0x0),
  fAliRunHeader(0x0),
  fPHOSGeom(0x0)
{
  // Constructor
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
  f1CutMaxDCAToVertexXYPtDep(ap.f1CutMaxDCAToVertexXYPtDep), 
  fTrackList(ap.fTrackList),
  fMCTrackList(ap.fMCTrackList),
  fPhotonList(ap.fPhotonList),
  fHeaderList(ap.fHeaderList),
  fAliRunHeader(ap.fAliRunHeader),
  fPHOSGeom(ap.fPHOSGeom)
{ 
  // cpy ctor
  for(int k=0; k < kRangeTriggerTableAlice; k++)
    fActiveTriggers[k] = ap.fActiveTriggers[k];
  
  for(int j=0; j < kRangeTriggerTableJCorran; j++)
    fTriggerTableJCorran[j] = ap.fTriggerTableJCorran[j];
}

//_____________________________________________________________________________
AliJCORRANTask& AliJCORRANTask::operator = (const AliJCORRANTask& ap)
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
  if(f1CutMaxDCAToVertexXYPtDep) delete f1CutMaxDCAToVertexXYPtDep;
  if(fTrackList)    delete fTrackList;
  if(fMCTrackList)  delete fMCTrackList;
  if(fPhotonList)   delete fPhotonList;
  if(fHeaderList)   delete fHeaderList;
  if(fAliRunHeader) delete fAliRunHeader;
  if(fPHOSGeom)     delete fPHOSGeom  ;
  GetEMCALGeoUtils(kTRUE);
}

//________________________________________________________________________
void AliJCORRANTask::UserCreateOutputObjects()
{  
  // create the jcorran outputs objects
  fTrackList    = new AliPhJTrackList(kALICE);
  if(fIsRealOrMC) fMCTrackList  = new AliPhJMCTrackList(kALICE);
  fPhotonList   = new AliPhJPhotonList(kALICE);
  fHeaderList   = new AliPhJHeaderList(kALICE);

  fAliRunHeader = new AliJRunHeader();

  fPHOSGeom  = new AliPHOSGeoUtils("PHOSgeo") ;

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
    AddAODBranch("AliPhJMCTrackList", &fMCTrackList, fAODName.Data());
  }
  
}

//______________________________________________________________________________
void AliJCORRANTask::UserExec(Option_t* /*option*/) 
{
  // Processing of one event
  //if(fDebug > 5) cout << "------- AliJCORRANTask Exec-------"<<endl;
  if(!((Entry()-1)%100)) 
      AliInfo(Form(" Processing event # %lld",  Entry())); 


  Bool_t storeEvent = kFALSE;//based on offline trigger decide whether to store the event or not 
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
  if(fIsRealOrMC) fMCTrackList->Reset();
  fPhotonList->Reset();
  fHeaderList->Reset();

  static int runId=-1;
 
  if(fInputFormat=="ESD"){   //   Reading ESD  
    AliESDEvent* esd = (AliESDEvent*)InputEvent();
    if(!esd) return;
                       //cout << storeEvent<<"    "<<esd->GetFiredTriggerClasses() << endl;
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

       //Create internal JCorran trigger mask.  Mapping between trigger and trigger bit
       fTriggerTableJCorran[kMinBiasTriggerBitJCorran]="Minimum Bias";//minimum bias occupies first trigger bit
       fTriggerTableJCorran[kHighMultTriggerBitJCorran]="High Multiplicity";//high multiplicity trigger => second trigger bit
      
       //=========== Fill Run header object ===============
       fAliRunHeader->SetRunNumber(runId);
       fAliRunHeader->SetActiveTriggersAlice(fActiveTriggers);
       fAliRunHeader->SetL3Field(l3MgFieldPolarity, esd->GetMagneticField());
       fAliRunHeader->SetActiveTriggersJCorran(fTriggerTableJCorran,kRangeTriggerTableJCorran);

      //Store Run header
      ( OutputTree()->GetUserInfo())->Add(fAliRunHeader); 
    }//end RunHeader

    if(storeEvent){ //FK//
      //-------------- reset all the arrays -------------
      //store event only when it is downscaled min bias
      // or contais high pt hadron
      // or contains high energy cluster in EMCAL or PHOS
      ReadESDTracks(esd);
      ReadESDCaloClusters(esd);
      ReadESDHeader(esd);
      if(fIsRealOrMC) ReadMCTracks(mcEvent);
    }
  }else if( fInputFormat == "AOD") {
  
    AliAODEvent* aod = AODEvent();
    if(!aod) return; 
   
    //=========== FILL AND STORE RUN HEADER   (ONLY ONCE PER RUN) =============
    if(aod->GetRunNumber() != runId){ //new run has started
 
      runId = aod->GetRunNumber();
      // trigger names???

      //Polarity of magnetic field in L3 solenoid
      Short_t l3MgFieldPolarity=0; // (LHC convention: +current -> +Bz)
      if( aod->GetMagneticField() > 0 ) l3MgFieldPolarity =  1;
      if( aod->GetMagneticField() < 0 ) l3MgFieldPolarity = -1;

      fAliRunHeader->SetRunNumber(runId);
      fAliRunHeader->SetL3Field(l3MgFieldPolarity, aod->GetMagneticField());
      //Store Run header
      (OutputTree()->GetUserInfo())->Add(fAliRunHeader);  

    }

    ReadAODTracks(aod);
    ReadAODCaloClusters(aod);
    ReadAODHeader(aod);
    
  }else{
    cout << "Error: Not correct InputDataFormat especified " << endl;
    return;
  }


  if(fTrackList->GetNTracks() > 0 || 
    fPhotonList->GetNPhotons() >0 || 
    ( fIsRealOrMC && fMCTrackList->GetNTracks()>0)){

    AliAODHandler* outputHandler = 
    (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    outputHandler->SetFillAOD(kTRUE);
  }
}

//______________________________________________________________________________
void AliJCORRANTask::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialization") ; 

  TString formula(fEsdTrackCuts->GetMaxDCAToVertexXYPtDep());
  if(formula.Length()>0){ // momentum dep DCA cut for AOD
    formula.ReplaceAll("pt","x");
    if(f1CutMaxDCAToVertexXYPtDep) delete f1CutMaxDCAToVertexXYPtDep; 
    f1CutMaxDCAToVertexXYPtDep = new TFormula("f1CutMaxDCAToVertexXYPtDep",formula.Data());
  }
}

//______________________________________________________________________________
void AliJCORRANTask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
 // OutputTree()->Print();

  // ((AliJRunHeader *) (( OutputTree()->GetUserInfo())->First()))->PrintOut();  

  cout<<"PWG4JCORRAN Analysis DONE !!"<<endl; 
}

//______________________________________________________________________________
void AliJCORRANTask::ReadESDTracks(const AliESDEvent * esd)
{
    // Read the AliESDtrack and fill the list of AliJTrack containers
    Int_t nt = esd->GetNumberOfTracks();
   // if(fDebug < 5) cout << "ESD::NumberOfTracks = " << nt << endl;
    Short_t ntrk = 0;
    Double_t pid[10];
    
    //loop over tracks
    for(Int_t it = 0; it < nt; it++) { 

        AliESDtrack *track = esd->GetTrack(it);
        if(! fEsdTrackCuts->IsSelected(track)) continue; // apply track selection criteria

	//------------ T P C ------------
        Float_t tpcDca[2], tpcCov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
	track->GetImpactParametersTPC(tpcDca,tpcCov);

	Int_t nClust         = track->GetTPCclusters(0);
	Int_t nFindableClust = track->GetTPCNclsF();
	Float_t tpcChi2PerCluster = 0.;
	if(nClust>0.) tpcChi2PerCluster = track->GetTPCchi2()/Float_t(nClust);

	Float_t tpcClustPerFindClust = 0.;
	if(nFindableClust>0.) tpcClustPerFindClust = Float_t(nClust)/nFindableClust;
        //--------------------------------
       
             //create a new AliJTrack and fill the track info
	fTrackList->AddAliJTrack(ntrk);
	AliJTrack *ctrack = fTrackList->GetAliJTrack(ntrk);

        ctrack->SetPtot(track->P());//here
        ctrack->SetPt(track->Pt());
        ctrack->SetTheta(track->Theta());
        ctrack->SetPhi(atan2(sin(track->Phi()), cos(track->Phi())));
        track->GetESDpid(pid);
        ctrack->SetPID(pid);

        ctrack->SetFlavor(kNone);
	ctrack->SetCharge(track->Charge());
        ctrack->ConvertAliPID();
        ctrack->SetEta(track->Eta());

	Int_t itof = track->GetTOFcluster(); // index of the assigned TOF cluster
	if(itof>=0) ctrack->SetTof(track->GetTOFsignal());

	ctrack->SetTPCdEdx(track->GetTPCsignal());
	ctrack->SetTPCnClust(track->GetTPCNcls());
	ctrack->SetTPCDCAXY(tpcDca[0]);
	ctrack->SetTPCDCAZ(tpcDca[1]);
	ctrack->SetTPCClustPerFindClust(tpcClustPerFindClust);
	ctrack->SetTPCChi2PerClust(tpcChi2PerCluster);

	Float_t impactXY, impactZ; //get track impact parameters
	track->GetImpactParameters(impactXY,impactZ);
	ctrack->SetImapactXY(impactXY);
	ctrack->SetImapactZ(impactZ);

        ctrack->SetKinkIndex(track->GetKinkIndex(0));
        ctrack->SetStatus(track->GetStatus()); // ULong_t 

        if(fIsRealOrMC){
          ctrack->SetITSLabel(track->GetITSLabel());
          ctrack->SetTPCLabel(track->GetTPCLabel());
        }

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
        if(!AcceptAODTrack(track)) continue; 
        //if(! fEsdTrackCuts->IsSelected(track)) continue; //apply loose selection criteria
       //FK//if(track->GetType() != AliAODTrack::kPrimary) continue; // only primaries 
     
        //create a new AliJTrack and fill the track info
	fTrackList->AddAliJTrack(ntrk);
        AliJTrack *ctrack = fTrackList->GetAliJTrack(ntrk);

        ctrack->SetPtot(track->P());
        ctrack->SetPt(track->Pt());
        ctrack->SetTheta(track->Theta());
        //ctrack->SetPhi(track->Phi());
        ctrack->SetPhi(atan2(sin(track->Phi()), cos(track->Phi())));
        ctrack->SetPID((Double_t*)track->PID());
        ctrack->SetFlavor(kNone);
        ctrack->SetCharge(track->Charge());
        ctrack->SetEta(track->Eta());
     
        AliAODPid* pidAOD = track->GetDetPid();
        if(pidAOD){ 
          ctrack->SetTof(pidAOD->GetTOFsignal());
          ctrack->SetTPCdEdx(pidAOD->GetTPCsignal());
        } 
 
        Double_t impactDCA[3];
        if( track->GetPosition(impactDCA)){
          ctrack->SetImapactXY(sqrt(impactDCA[0]*impactDCA[0] + impactDCA[1]*impactDCA[1]));//impactXY);
          ctrack->SetImapactZ(impactDCA[2]);//impactZ);          
        }       

        ctrack->SetChi2perNDF(track->Chi2perNDF());
        ctrack->SetTPCnClust(track->GetTPCNcls());
        ctrack->SetChi2Trig(track->GetChi2MatchTrigger());
     
        ctrack->SetStatus(track->GetStatus());//
        //ctrack->SetRecFlags(track->GetFlags());//?status

        if(fIsRealOrMC){
          Int_t  label = track->GetLabel();
          ctrack->SetITSLabel(label);
          ctrack->SetTPCLabel(label);       
        }

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
    if(caloCluster->GetNTracksMatched()<=0){
      if(caloCluster->E()<fLowerCutOnCaloClusterE) continue;   
 
      // we will not implement any PID cut here      
      fPhotonList->AddAliJPhoton(nPhotons);
      AliJPhoton *pht = fPhotonList->GetAliJPhoton(nPhotons);
      pht->SetFlavor(kNone);//kPhoton);
      pht->SetE(caloCluster->E());
      pht->SetChi2(caloCluster->Chi2());
      pht->SetPID(caloCluster->GetPID());
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
      pht->SetDispersion(caloCluster->GetDispersion());
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
  Int_t numberOfCaloClusters = aod->GetNumberOfCaloClusters();  
  if(fDebug < 5) cout << "AOD::number of ESD caloclusters " << numberOfCaloClusters << endl;
  Short_t nPhotons = 0;

  for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {

      AliAODCaloCluster * caloCluster = aod->GetCaloCluster(icluster) ;
      if(!caloCluster) continue; 
      if(caloCluster->GetNTracksMatched() > 0) continue;
      if(caloCluster->E() < fLowerCutOnCaloClusterE) continue; 
      // we will not implement any PID cut here      
      fPhotonList->AddAliJPhoton(nPhotons);

      AliJPhoton *pht = fPhotonList->GetAliJPhoton(nPhotons);
      pht->SetFlavor(kNone);
      pht->SetE(caloCluster->E());
      pht->SetChi2(caloCluster->Chi2());
      pht->SetPID((Double_t*)caloCluster->GetPID());
      Float_t pos[3]; caloCluster->GetPosition(pos);
      pht->SetX(pos[0]);
      pht->SetY(pos[1]);
      pht->SetZ(pos[2]);
      pht->SetPhi(atan2(pos[1],pos[0]));
      pht->SetTheta(atan2(sqrt(pos[0]*pos[1]+pos[1]*pos[1]),pos[2]));
      pht->SetPt(caloCluster->E()*sin(atan2(sqrt(pos[0]*pos[0]+pos[1]*pos[1]),pos[2])));
      pht->SetCharge(0);
      if(caloCluster->IsEMCAL()) pht->SetCaloType(AliJPhoton::kEMCALCalo);
      if(caloCluster->IsPHOS())  pht->SetCaloType(AliJPhoton::kPHOSCalo);
      pht->SetDistToBadChannel(caloCluster->GetDistanceToBadChannel());
      pht->SetDispersion(caloCluster->GetDispersion());
      pht->SetM20(caloCluster->GetM20());
      pht->SetM02(caloCluster->GetM02());
      pht->SetEmcCpvDist(caloCluster->GetEmcCpvDistance());
      pht->SetNCells(int(caloCluster->GetNCells()));
      pht->SetCellsAbsId(caloCluster->GetCellsAbsId());
      Int_t imoduleID = GetSuperModuleNumber(caloCluster->IsEMCAL(), caloCluster->GetCellAbsId(0));
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
    if(fSPDMult) hdr->SetSPDTrackletMult(fSPDMult->GetNumberOfTracklets());

    hdr->SetEventID( esd->GetEventNumberInFile());

    hdr->SetTriggerMaskAlice(esd->GetTriggerMask()); //ULong64_t
    hdr->SetTriggerMaskJCorran(ConvertTriggerMask()); //UInt_t

    const AliESDVertex * vtxESD = esd->GetVertex();
    if(vtxESD){
      hdr->SetZVertex(vtxESD->GetZv());
      hdr->SetZVertexErr(vtxESD->GetZRes());
    }
    hdr->SetEventType(esd->GetEventType());
    
    fHeaderList->SetNHeaders(++nHeaders);
}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODHeader(const AliAODEvent *aod)
{  
   //Read the AliAODEvent and fill the list of AliJHeader containers
   static int eventID = 0; //FK//?? dummy indexing of events (I cannot see how to get evet ID from AOD) 

    //read AOD event header
    Short_t nHeaders = 0;

    //create a header and fill it
    fHeaderList->AddAliJHeader(nHeaders);
    AliJHeader *hdr = fHeaderList->GetAliJHeader(nHeaders);
	
    const AliAODTracklets *trackletsSPD = aod->GetTracklets();
    if(trackletsSPD){
      hdr->SetSPDTrackletMult(trackletsSPD->GetNumberOfTracklets());
    }
    hdr->SetEventID( eventID++ );  //FK//?? I cannot see how to get evet ID from AOD	
    const AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
    if(vtxAOD){
      hdr->SetZVertex((float) vtxAOD->GetZ());
      float sigmaVtx[3];
      vtxAOD->GetSigmaXYZ(sigmaVtx);
      hdr->SetZVertexErr((float) sqrt(sigmaVtx[2]));//TMath::Sqrt(fCovZZ)
    }

    //load aod event header
    AliAODHeader * aodh = aod->GetHeader();
    if(aodh){
      hdr->SetCentrality(int(aodh->GetCentrality())); 
      hdr->SetEventType(aodh->GetEventType());
      hdr->SetTriggerMaskAlice(aodh->GetTriggerMask()); //ULong64_t
    }
    
    hdr->SetTriggerMaskJCorran(ConvertTriggerMask()); //UInt_t
   
    fHeaderList->SetNHeaders(++nHeaders);
}

//______________________________________________________________________________
Int_t AliJCORRANTask::GetSuperModuleNumber(bool isemcal, Int_t absId)
{
  //get super module number 
  if(isemcal){
    return GetEMCALGeoUtils()->GetSuperModuleNumber(absId) ;
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

UInt_t AliJCORRANTask::ConvertTriggerMask(){

  //convert alice trigger mask to jcorran trigger mask
  UInt_t triggerMaskJC=0;
  if(!fIsRealOrMC){  //REAL data 
    if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
       ->IsEventSelected() & AliVEvent::kMB){
      // minimum bias TBit 0 
      triggerMaskJC += (1<<kMinBiasTriggerBitJCorran); 
    }

    if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
       ->IsEventSelected() & AliVEvent::kHighMult){
      //high multiplicity trigger TBit 1 
      triggerMaskJC += (1<<kHighMultTriggerBitJCorran);
    }
  }else{
    triggerMaskJC=1; //MC data, at the moment all events filled as MB
  }

  return triggerMaskJC;
}


//______________________________________________________________________________
bool AliJCORRANTask::StoreDownscaledMinBiasEvent(){
  //Decide whether to downscale this MinBiasEvent and store it 
  bool isThisEventToBeStored = kFALSE;
  static Long_t evtNumber=0;

  if( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() ){
    //collision candidate

    if(evtNumber % fDownscaling ==0){ isThisEventToBeStored = kTRUE; } //store every Xth collision candidate event
    evtNumber++;
  }
  return isThisEventToBeStored;
}
//______________________________________________________________________________
bool AliJCORRANTask::ContainsESDHighPtTrack(){
  //If there was an identified high pT particle above threshold reutrn flag to store this event 
  bool isThisEventToBeStored = kFALSE;

  if(fInputFormat=="ESD"){  // E S D

    AliESDEvent* esd = (AliESDEvent*)InputEvent();
    if(!esd) return kFALSE;
    Int_t nt = esd->GetNumberOfTracks();

    for(Int_t it = 0; it < nt; it++){
      AliESDtrack *track = esd->GetTrack(it);
      //Does event contain high pt particle above thereshold in GeV 
      if(track->Pt() > fLowerCutOnLPMom && fEsdTrackCuts->IsSelected(track)){
        isThisEventToBeStored = kTRUE;
        break;
      }
    }
  }else{  // A O D
    AliAODEvent* aod=(AliAODEvent*) InputEvent(); 
    if(!aod) return kFALSE;
    Int_t nt = aod->GetNumberOfTracks();

    for(Int_t it = 0; it < nt; it++) {
      AliAODTrack *track = aod->GetTrack(it);
      //Does event contain high pt particle above threshold in GeV 
      if(track->Pt() > fLowerCutOnLPMom && AcceptAODTrack(track)){
        isThisEventToBeStored = kTRUE;      
        break;
      }
    }
  }

  return isThisEventToBeStored;
}

//______________________________________________________________________________
bool AliJCORRANTask::ContainsESDHighECaloClusters(){
  //Check whether there was in the event high E calo clustre and renturn a flag whether to store this event. 
  bool isThisEventToBeStored = kFALSE;

  if(fInputFormat=="ESD"){

    AliESDEvent* esd = (AliESDEvent*)InputEvent();
    if(!esd) return kFALSE;
    Int_t numberOfCaloClusters = esd->GetNumberOfCaloClusters() ;
    // loop over all the Calo Clusters
    for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {
      AliESDCaloCluster *caloCluster = esd->GetCaloCluster(icluster) ;
      if(!caloCluster) continue;
      if(caloCluster->GetNTracksMatched() ==-1){
        //sotre calo clusters above 1 GeV
        if( caloCluster->E() > fLowerCutOnLeadingCaloClusterE){
          isThisEventToBeStored = kTRUE;
          break;
        }
      }
    }
  }else{

    AliAODEvent* aod = (AliAODEvent*)InputEvent();
    if(!aod) return kFALSE;
    Int_t numberOfCaloClusters = aod->GetNumberOfCaloClusters();
    // loop over all the Calo Clusters
    for(Int_t icluster = 0 ; icluster < numberOfCaloClusters ; icluster++) {
      AliAODCaloCluster *caloCluster = aod->GetCaloCluster(icluster) ;
      if(!caloCluster) continue;
      if(caloCluster->GetNTracksMatched() == -1){ 
        //sotre calo clusters above 1 GeV
        if( caloCluster->E() > fLowerCutOnLeadingCaloClusterE){
          isThisEventToBeStored = kTRUE;
          break;
        }
      }
    }
  }

  return isThisEventToBeStored;
}


//--------------------------------------------------------------------

void AliJCORRANTask::ReadMCTracks(AliMCEvent *fMC){
  //store MC information from AliStack
  AliStack *stack = fMC->Stack();
  Int_t np    = fMC->GetNumberOfTracks();

  //  AliGenEventHeader* genHeader = fMC->GenEventHeader();
  //  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  //  Double_t ptHard = 0;
  //  Double_t nTrials = 1; // Trials for MC trigger weigth for real data
  //  nTrials = pythiaGenHeader->Trials();
  //  ptHard  = pythiaGenHeader->GetPtHard();
  //  Int_t nprim = stack->GetNtrack();
 
  Short_t ntrack = 0;

  for(Int_t iTrack = 0; iTrack < np; iTrack++){
    AliMCParticle *track = (AliMCParticle*) fMC->GetTrack(iTrack);
    if(!track){
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    Bool_t isPrimary = fMC->Stack()->IsPhysicalPrimary(iTrack);
    if(isPrimary){
      //create a new JMCTrack and fill the track info
      fMCTrackList->AddJMCTrack(ntrack);
      AliJMCTrack *ctrack = fMCTrackList->GetTrack(ntrack);
  
      TParticle *partStack = stack->Particle(iTrack);
      Int_t   pdg  = partStack->GetPdgCode();
      Float_t engy = partStack->Energy();
      Float_t pt   = partStack->Pt();
      Float_t ptot = partStack->P();
      Float_t eta  = partStack->Eta();
      Float_t theta  = partStack->Theta();
      Float_t phi    = atan2(sin( partStack->Phi()), cos(partStack->Phi()));
      Short_t ch     = (Short_t) partStack->GetPDG()->Charge();
      Int_t label    = track->GetLabel();

      ctrack->SetLabel(label);
      ctrack->SetPdgCode(pdg);
      ctrack->SetPt(pt);
      ctrack->SetTheta(theta);
      ctrack->SetEta(eta);
      ctrack->SetPhi(phi);
      ctrack->SetE(engy);
      ctrack->SetCharge(ch);
      ctrack->SetPtot(ptot);
      ctrack->SetIsPrimary(isPrimary);

      ctrack->SetProductionVertex(partStack->Vx(),partStack->Vy(),partStack->Vz());
   /*
      Int_t   status  = partStack->GetStatusCode();
      ctrack->SetStatusCode(status);

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
*/
      fMCTrackList->SetNTracks(++ntrack);
    }// loop for al primary tracks
  }
}

//--------------------------------------------------------------------
bool AliJCORRANTask::AcceptAODTrack(AliAODTrack* aodTrack){
  //This function mimics for the AliAODTracks object the AliESDtrackCut function IsSelected
  //Cuts are taken from fEsdTrackCuts object 
  if(fEsdTrackCuts->GetMinNClusterTPC() > aodTrack->GetTPCNcls()) return kFALSE;

  //if(fEsdTrackCuts->GetMaxChi2PerClusterTPC() <    );//<-------- how to check? 
  //      ctrack->SetChi2perNDF(track->Chi2perNDF());

  //         C h e c k    r e f i t

  if(fEsdTrackCuts->GetRequireTPCRefit()  && 
    ((aodTrack->GetStatus() & AliJTrack::kTPCrefit) == 0)) return kFALSE;
  if(fEsdTrackCuts->GetRequireITSRefit()  && 
    ((aodTrack->GetStatus() & AliJTrack::kITSrefit) == 0)) return kFALSE;


   //            C u t s   o n    D C A
   Float_t impactDCA[3];
   if( aodTrack->GetPosition(impactDCA)){
     if((fEsdTrackCuts->GetMaxDCAToVertexXY()>0) && 
       (fEsdTrackCuts->GetMaxDCAToVertexXY() < sqrt(impactDCA[0]*impactDCA[0] + impactDCA[1]*impactDCA[1]))) return kFALSE;
     if((fEsdTrackCuts->GetMaxDCAToVertexZ()>0) &&
       (fEsdTrackCuts->GetMaxDCAToVertexZ() < TMath::Abs(impactDCA[2]))) return kFALSE;
   } else return kFALSE; 
 
   if(f1CutMaxDCAToVertexXYPtDep){ 
     if(f1CutMaxDCAToVertexXYPtDep->Eval(aodTrack->Pt()) < sqrt(impactDCA[0]*impactDCA[0] + impactDCA[1]*impactDCA[1])) return kFALSE; 
   }


  // if(fEsdTrackCuts->GetAcceptKinkDaughters()) //<--------how to check ?
  //  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);


  return kTRUE; 
}

//--------------------------------------------------------------------
AliEMCALGeometry * AliJCORRANTask::GetEMCALGeoUtils (bool doDelete){
  //include EMCAL singleton modul
  static AliEMCALGeometry* emcalgeo = 0x0;
  if( ! emcalgeo ){
    emcalgeo = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
  }
  if( emcalgeo && doDelete ){ //FK// !emcalgeo
    delete emcalgeo;
    emcalgeo=0x0;
  }
  return emcalgeo;
}



