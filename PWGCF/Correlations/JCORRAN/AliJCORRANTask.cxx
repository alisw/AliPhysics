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
#include <TVectorT.h> 
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
#include "AliAODVertex.h" 
#include "AliAODTracklets.h" 
#include "AliAODPid.h" 
#include "AliESDUtils.h"
//#include "AliESDVZERO.h" 
#include "AliCentrality.h" 
#include "AliAODTracklets.h"
#include "AliMultiplicity.h"
#include "AliJConst.h"
#include "AliESDRun.h"
#include "AliESDVZERO.h"
#include "AliExternalTrackParam.h"
//== EMCAL
#include "AliESDCaloCluster.h"
#include "AliEMCALGeometry.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALPIDUtils.h"

#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJPhoton.h"
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"

//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask() :   
    AliAnalysisTaskSE("PWG4JCORRAN"),
    fRunType("LHC10h"), // enable filling EP info
    fInputFormat("ESD"),
    fEsdTrackCuts(0x0), 
    fESDFilter(0x0), 
    fIsRealOrMC(), 
    fAODName("jcorran.root"),
    fStoreEventPlaneSource(false), 
    fStoreTPCTrack(false), 
    fOADBPath(),
    fTrackList(0),
    fMCTrackList(0x0),
    fPhotonList(0x0),
    fHeaderList(0x0),
    fRunInfoList(0x0),
    fPIDesd(0x0),
    fPIDResponse(0x0),
    fPIDCombined(0x0),
    fVZEROData(0x0), 
    fTZEROData(0x0), 
    //fFMDData(0x0), 
    fZDCData(0x0), 
    fAliRunHeader(0x0),
    fEMCALGeoUtils(0x0),    
    fPHOSGeom(0x0)
{
  //Default constructor
  for(Int_t i=0;i<kRangeTriggerTableAlice;i++)   fActiveTriggers[i]=" ";
  for(Int_t i=0;i<kRangeTriggerTableJCorran;i++) fTriggerTableJCorran[i]=" ";

  fIsRealOrMC.ResizeTo(1);
  fIsRealOrMC[0]=0;

  DefineInput (0, TChain::Class());
  DefineInput (1, TList::Class());
  DefineOutput (1, TList::Class());
  //  DefineOutput (2, TList::Class());
}

//______________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
    fRunType("LHC10h"), // enable filling EP info
    fInputFormat(inputformat),  
    fEsdTrackCuts(0x0), 
    fESDFilter(0x0), 
    fIsRealOrMC(), 
    fAODName("jcorran.root"),
    fStoreEventPlaneSource(false), 
    fStoreTPCTrack(false), 
    fOADBPath(),
    fTrackList(0),
    fMCTrackList(0x0),
    fPhotonList(0x0),
    fHeaderList(0x0),
    fRunInfoList(0x0),
    fPIDesd(0x0),
    fPIDResponse(0x0),
    fPIDCombined(0x0),
    fVZEROData(0x0), 
    fTZEROData(0x0), 
    //fFMDData(0x0), 
    fZDCData(0x0), 
    fAliRunHeader(0x0),
    fEMCALGeoUtils(0x0),    
    fPHOSGeom(0x0)
{
  // Constructor
  if(fDebug > 5) cout << "---- AliJCORRANTask Constructor ----"<<endl;

  for(Int_t i=0;i<kRangeTriggerTableAlice;i++)   fActiveTriggers[i]=" ";
  for(Int_t i=0;i<kRangeTriggerTableJCorran;i++) fTriggerTableJCorran[i]=" ";

  fIsRealOrMC.ResizeTo(1);

  fIsRealOrMC[0]=0;

  DefineInput (0, TChain::Class());
  //  DefineInput (1, TList::Class());
  DefineOutput (1, TList::Class());
  //  DefineOutput (2, TList::Class());
}

//____________________________________________________________________________
AliJCORRANTask::AliJCORRANTask(const AliJCORRANTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
    fRunType(ap.fRunType), 
    fInputFormat(ap.fInputFormat),
    fEsdTrackCuts(ap.fEsdTrackCuts), 
    fESDFilter(ap.fESDFilter), 
    fIsRealOrMC(ap.fIsRealOrMC), 
    fAODName(ap.fAODName), 
    fStoreEventPlaneSource(ap.fStoreEventPlaneSource), 
    fStoreTPCTrack(ap.fStoreTPCTrack), 
    fOADBPath(ap.fOADBPath),
    fTrackList(ap.fTrackList),
    fMCTrackList(ap.fMCTrackList),
    fPhotonList(ap.fPhotonList),
    fHeaderList(ap.fHeaderList),
    fRunInfoList(ap.fRunInfoList),
    fPIDesd(ap.fPIDesd),
    fPIDResponse(ap.fPIDResponse),
    fPIDCombined(ap.fPIDCombined),
    fVZEROData(ap.fVZEROData), 
    fTZEROData(ap.fTZEROData), 
    //fFMDData(ap.fFMDData), 
    fZDCData(ap.fZDCData), 
    fAliRunHeader(ap.fAliRunHeader),
    fEMCALGeoUtils(ap.fEMCALGeoUtils),    
    fPHOSGeom(ap.fPHOSGeom)
{ 
  // cpy ctor
  for(int k=0; k < kRangeTriggerTableAlice; k++)
    fActiveTriggers[k] = ap.fActiveTriggers[k];

  for(int j=0; j < kRangeTriggerTableJCorran; j++)
    fTriggerTableJCorran[j] = ap.fTriggerTableJCorran[j];

  fIsRealOrMC.ResizeTo(1);

  fIsRealOrMC[0]                    = ap.fIsRealOrMC[0];

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
  delete fTrackList;
  delete fMCTrackList;
  delete fPhotonList;
  delete fHeaderList;
  delete fAliRunHeader;
  delete fRunInfoList;
  delete fPIDesd;
  delete fOADBPath;
  delete fPIDResponse;
  delete fPIDCombined;
  delete fEMCALGeoUtils;
  delete fPHOSGeom;
  delete fVZEROData;
  delete fTZEROData;
  delete fZDCData;
  //  delete fFMDData;


}

//________________________________________________________________________

void AliJCORRANTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJCORRANTask::UserCreateOutPutData() \n");

  //=== Get AnalysisManager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if(!man->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }
  man->RegisterExtraFile(fAODName.Data());

  //=== Other Objects
  fEMCALGeoUtils = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE");
  fPHOSGeom = new AliPHOSGeoUtils();

  //=== Set Tree and TClonesArray
  //== TRACKS
  AddListAODBranch("AliJTrackList", "AliJTrack", &fTrackList, 1000);
  AddListAODBranch("AliJPhotonList", "AliJPhoton", &fPhotonList, 1000);
  if((bool)fIsRealOrMC[0]) 
    AddListAODBranch("AliJTMCrackList", "AliJMCTrack", &fMCTrackList, 1000);
  //== Event Header
  AddListAODBranch("AliJEventHeaderList", "AliJEventHeader", &fHeaderList, 1000);
  //== RUN HEADER
  fAliRunHeader = new AliJRunHeader();
  fRunInfoList	= new TList();
  fRunInfoList->SetName("RunInfoList");
  fRunInfoList->SetOwner();
  fRunInfoList->Clear();
  fRunInfoList->Add(fAliRunHeader);
  //== EventPlane SRC
  if( fStoreEventPlaneSource ){
    fVZEROData = new AliESDVZERO;
    fTZEROData = new AliESDTZERO;
    fZDCData   = new AliESDZDC;
    AddAODBranch("AliESDVZERO", &fVZEROData,  fAODName.Data());
    AddAODBranch("AliESDTZERO", &fTZEROData,  fAODName.Data());
    AddAODBranch("AliESDZDC",   &fZDCData,    fAODName.Data());
  }
  //== PID
  fPIDesd = new AliESDpid;
  fPIDCombined = new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  fPIDResponse = ((AliInputEventHandler*) (man->GetInputEventHandler()))->GetPIDResponse();
  fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
  if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());

  cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;
  PostData(1,fRunInfoList);

}

//______________________________________________________________________________
void AliJCORRANTask::UserExec(Option_t* /*option*/) 
{

  // Processing of one event
  if(fDebug > 5) cout << "------- AliJCORRANTask Exec-------"<<endl;
  if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry())); 

  //=== Init Variables
  fTrackList->Clear();
  if((bool)fIsRealOrMC[0]) fMCTrackList->Clear();
  fPhotonList->Clear();
  fHeaderList->Clear();

  //=== COMMON for ESD and AOD
  static int runId=-1;
  AliVEvent *event = InputEvent();
  if(!event) return;
  AliMCEvent* mcEvent = NULL;
  if((bool)fIsRealOrMC[0])  mcEvent = MCEvent();
  // RUN Header
  if(event->GetRunNumber() != runId){ //new run has started
    runId = event->GetRunNumber();
    //Polarity of magnetic field in L3 solenoid
    Short_t l3MgFieldPolarity=0; // (LHC convention: +current -> +Bz)
    //Create internal JCorran trigger mask.  Mapping between trigger and trigger bit
    fTriggerTableJCorran[kMinBiasTriggerBitJCorran]="Minimum Bias";//minimum bias occupies first trigger bit
    fTriggerTableJCorran[kHighMultTriggerBitJCorran]="High Multiplicity";//high multiplicity trigger => second trigger bit
    //=========== Fill Run header object ===============
    fAliRunHeader->SetRunNumber(runId);
    fAliRunHeader->SetActiveTriggersJCorran(fTriggerTableJCorran,kRangeTriggerTableJCorran);
    SetAliceTriggerDef(fAliRunHeader);//TODO for AOD
    SetAliceFilterMapDef(fAliRunHeader);// TODO for AOD
    //FOR ESD
    if(fInputFormat=="ESD"){
      AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
      if(!esd) return;
      if(esd->GetCurrentL3() >0) l3MgFieldPolarity =  1;
      if(esd->GetCurrentL3() <0) l3MgFieldPolarity = -1;
      fAliRunHeader->SetL3Field(l3MgFieldPolarity, esd->GetMagneticField());
      const AliESDRun* esdRun = esd->GetESDRun();
      for(Int_t triggerBit=0; triggerBit<kRangeTriggerTableAlice; triggerBit++){
        fActiveTriggers[triggerBit] = esdRun->GetTriggerClass(triggerBit);
      }
      fAliRunHeader->SetActiveTriggersAlice(fActiveTriggers);
    }
    fRunInfoList->Add(fAliRunHeader);
    cout << "Add(fAliRunHeader) is done =============" << endl;
  }

  //=== If ESD or AOD
  if(fInputFormat=="ESD"){   //Reading ESD  
    if(fDebug > 5) cout << "\t------- Start ESD "<<endl;
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if(!esd) return;
    ReadESDHeader(esd);
    ReadESDTracks(esd);
    //ReadESDCaloClusters(esd);
    if((bool)fIsRealOrMC[0]) ReadMCTracks(mcEvent);
  }else if( fInputFormat == "AOD") {
    if(fDebug > 5) cout << "\t------- Start AOD "<<endl;
    AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
    if(!aod) return;
    ReadAODHeader(aod);
    ReadAODTracks(aod);
    //ReadAODCaloClusters(aod);
  }else{
    cout << "Error: Not correct InputDataFormat especified " << endl;
    return;
  }

  //=== TODO : need this?
  AliAODHandler* outputHandler = 
      (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  outputHandler->SetFillAOD(kTRUE);
  outputHandler->SetFillExtension(kTRUE);
  PostData(1,fRunInfoList);

  if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJCORRANTask::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialization") ; 

  TString formula(fEsdTrackCuts->GetMaxDCAToVertexXYPtDep());
  if(formula.Length()>0){ // momentum dep DCA cut for AOD
    formula.ReplaceAll("pt","x");
  }
}

//______________________________________________________________________________
void AliJCORRANTask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  cout<<"PWG4JCORRAN Analysis DONE !!"<<endl; 
  // Printout fRunInfoList here
  fRunInfoList = dynamic_cast<TList*> (GetOutputData(1));
  if(fRunInfoList)
  {
    fAliRunHeader = dynamic_cast<AliJRunHeader*> (fRunInfoList->FindObject("AliJRunHeader"));
    if(fAliRunHeader) {fAliRunHeader->Print();}
  }
  else
  {
    cout << "WARNING : Run Information List is empty" << endl;
  }

}

//______________________________________________________________________________
void AliJCORRANTask::ReadESDTracks(AliESDEvent * esd)
  //void AliJCORRANTask::ReadESDTracks(const AliESDEvent * esd)
{
  // Read the AliESDtrack and fill the list of AliJTrack containers
  Int_t nt = esd->GetNumberOfTracks();
  if(fDebug > 5) cout << "ESD::NumberOfTracks = " << nt << endl;
  Short_t ntrk = 0;

  //loop over tracks
  for(Int_t it = 0; it < nt; it++) { 

    AliESDtrack *track = esd->GetTrack(it);
    if( !track ) continue;
    if(! fEsdTrackCuts->IsSelected(track)) continue; // apply track selection criteria
    UInt_t filterMap = fESDFilter->IsSelected( track );
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
    AliJTrack * ctrack = new( (*fTrackList)[fTrackList->GetEntriesFast()] ) AliJTrack;
    ctrack->SetPxPyPzE(track->Px(), track->Py(), track->Pz(), 0 );
    ctrack->SetTPCdEdx( track->GetTPCsignal()  );
    if( fStoreTPCTrack ){
      AliESDtrack *tpcTrack = AliESDtrackCuts::GetTPCOnlyTrack( esd,  it );
      if( !tpcTrack ) continue;
      ctrack->SetTPCTrack( tpcTrack->Px(),  tpcTrack->Py(),  tpcTrack->Pz());
    }
    ReadESDPID( track, ctrack );
    ctrack->SetParticleType(kNone);
    ctrack->SetCharge(track->Charge());
    ctrack->SetFilterMap( filterMap );

    if(fDebug > 5 && track->P()>1 ) cout << "P = " << track->P() << endl;

    ++ntrk;
  } // end tracks loop
}

//_________________________________________________________________________________-
void AliJCORRANTask::ReadESDPID(AliESDtrack *track, AliJTrack *ctrack)
{
  // Probability calculation for PID
  Double_t probTPC[AliPID::kSPECIES]={0.};
  Double_t probTOF[AliPID::kSPECIES]={0.};
  Double_t probTPCTOF[AliPID::kSPECIES]={0.};
  //
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
  UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
  if (detUsed != 0) {  // TPC is available
    for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec) {
      Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)ispec);
      ctrack->SetTPCsigma(AliJTrack::AliJTrkPID(ispec), nSigmaTPC);
    }

    // compute priors for TPC+TOF, even if we ask just TOF for PID
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
    detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
    if (detUsed != 0) {  // TOF is available
      for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec) {
        Double_t nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)ispec);
        ctrack->SetTOFsigma(AliJTrack::AliJTrkPID(ispec), nSigmaTOF);
      }
    }

    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
    detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
  }

  for (int ip=0 ; ip < (AliPID::kSPECIES); ip++) {
    ctrack->SetPID(AliJTrack::AliJTrkPID(ip), probTOF[ip], AliJTrack::kTOF);
    ctrack->SetPID(AliJTrack::AliJTrkPID(ip), probTPC[ip], AliJTrack::kTPC);
    ctrack->SetPID(AliJTrack::AliJTrkPID(ip), probTPCTOF[ip], AliJTrack::kTPCTOF);
  }

  // TOF beta and expected beta
  Float_t beta = -1;  Float_t minP = 0.2; Float_t minTPCsignal = 40; Float_t minTOFLength = 365.; Float_t minTOFsignal = 12000;
  Double_t inttimes[10];  Float_t betaTh[10]; Double_t dEdxTh[10];
  track->GetIntegratedTimes(inttimes);
  for(int i=0;i<10;i++) {
    betaTh[i] = -1.; // initialize expected beta = -1
    dEdxTh[i] = -1.; // initialize expected dEdx = -1
  }
  if(track->P() > minP && track->GetTPCsignal() > minTPCsignal && (track->GetStatus() & AliESDtrack::kTOFout) 
     && (track->GetStatus() & AliESDtrack::kTIME) && (track->GetIntegratedLength() > minTOFLength)&& track->GetTOFsignal() > minTOFsignal) {
    Double_t consCal = 33.3564095198152043; 
    beta = track->GetIntegratedLength() / (track->GetTOFsignal() - fPIDesd->GetTOFResponse().GetStartTime(track->P())) * consCal;
    for(int i=0; i<10; i++) {
      betaTh[i] = track->GetIntegratedLength() / ( inttimes[i] ) * consCal;
    }
  }
  ctrack->SetTOFbeta(beta);

  Double_t ptpc[3];
  track->GetInnerPxPyPz(ptpc);
  Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);
  for(int ip=0; ip < (AliJTrack::kNAliJTrkPID); ip++) {
    ctrack->SetExpectedTOFbeta(AliJTrack::AliJTrkPID(ip), betaTh[ip]);
    // expected dEdx
    dEdxTh[ip] = fPIDesd->GetTPCResponse().GetExpectedSignal(momtpc, AliPID::EParticleType(ip));
    ctrack->SetExpectedTPCdEdx(AliJTrack::AliJTrkPID(ip), dEdxTh[ip]);
  }



}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODTracks(const AliAODEvent * aod)
{
  // Read the AliAODtrack and fill the list of AliJTrack containers
  Int_t nt = aod->GetNumberOfTracks();
  if(fDebug > 5) cout << "AOD::NumberOfTracks = " << nt << endl;
  cout << "AOD::NumberOfTracks = " << nt << endl;
  //Short_t ntrk = 0;
  //loop over tracks
  for(Int_t it = 0; it < nt; it++) { 

    AliAODTrack *track = aod->GetTrack(it);
    //if(!AcceptAODTrack(track)) continue; 
    //if(! fEsdTrackCuts->IsSelected(track)) continue; //apply loose selection criteria
    //FK//if(track->GetType() != AliAODTrack::kPrimary) continue; // only primaries 

    AliJTrack * ctrack = new( (*fTrackList)[fTrackList->GetEntriesFast()] ) AliJTrack;
    ctrack->SetPxPyPzE(track->Px(), track->Py(), track->Pz(), 0 );
    //TODO if( fStoreTPCTrack )
    ctrack->SetParticleType(kNone);
    ctrack->SetCharge(track->Charge());
    ctrack->SetStatus(track->GetStatus());//
    ctrack->SetFlags( track->GetFlags() );
    ctrack->SetLabel( track->GetLabel() );
    //FilterMap
    UInt_t filterMap=0;
    for( unsigned int i=0;i<sizeof(filterMap)*8;i++ ){
      filterMap |= track->TestFilterBit( 1UL<<i ); 
    }
    ctrack->SetFilterMap( filterMap );

    //PID TODO
    double const * pid = track->PID();
    ctrack->SetPID(AliJTrack::kElectronAliJ,pid[AliAODTrack::kElectron],AliJTrack::kTOF);
    ctrack->SetPID(AliJTrack::kMuonAliJ,    pid[AliAODTrack::kMuon],    AliJTrack::kTOF);
    ctrack->SetPID(AliJTrack::kPionAliJ,    pid[AliAODTrack::kPion],    AliJTrack::kTOF);
    ctrack->SetPID(AliJTrack::kKaonAliJ,    pid[AliAODTrack::kKaon],    AliJTrack::kTOF);
    ctrack->SetPID(AliJTrack::kProtonAliJ,  pid[AliAODTrack::kProton],  AliJTrack::kTOF);
    //TPC
    ctrack->SetTPCnClust(track->GetTPCNcls());
    ctrack->SetTPCdEdx( track->GetTPCsignal()  );

    if((bool) fIsRealOrMC[0]){
      //Int_t  label = track->GetLabel();
      //ctrack->SetITSLabel(label);
      //ctrack->SetTPCLabel(label);       
    }


    if(fDebug > 5 && track->P()>1 ) cout << "P = " << track->P() << endl;
  } // end tracks loop
}

//______________________________________________________________________________
void AliJCORRANTask::ReadESDCaloClusters(const AliESDEvent* esd)
{
  //AliVEvent * event = InputEvent();
  AliVEvent * event = (AliVEvent*)esd;
  TRefArray *caloClustersArr=new TRefArray();
  event->GetEMCALClusters(caloClustersArr);
  Double_t v[3];
  event->GetPrimaryVertex()->GetXYZ(v);

  AliEMCALRecoUtils * fRecoUtils = new AliEMCALRecoUtils;

  //const Int_t kNumberOfEMCALClusters =caloClustersArr->GetEntries();
  Int_t numberOfCaloClusters = caloClustersArr->GetEntries() ;
  if(fDebug > 5) cout << "ESD::number of ESD caloclusters " << numberOfCaloClusters << endl;

  AliVCaloCells *emCells =event->GetEMCALCells();

  int nPhotons = 0;
  for(Int_t icluster=0; icluster<numberOfCaloClusters; icluster++){
    AliVCluster *c1 = (AliVCluster*) caloClustersArr->At(icluster);
    if( !c1 ) continue;
    //== remove bad channels
    if(fRecoUtils->ClusterContainsBadChannel(fEMCALGeoUtils, c1->GetCellsAbsId(), c1->GetNCells())) continue;
    //== check energy and position
    if(fRecoUtils->IsRecalibrationOn()){
      fRecoUtils->RecalibrateClusterEnergy(fEMCALGeoUtils, c1, emCells);
      fRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeoUtils, emCells, c1);
      fRecoUtils->RecalculateClusterPID(c1);
    }
    //== correct non linearity
    c1->SetE(fRecoUtils->CorrectClusterEnergyLinearity(c1));

    //== corrected clusters
    if(c1->E() < 0.8 || c1->E() > 30) continue; //TODO
    //fRecoUtils->GetMaxEnergyCell(fEMCALGeo, emCells, c1, absID1, iSM, ieta1, iphi1, shared); 

    AliJPhoton *pht = new( (*fPhotonList)[nPhotons++] ) AliJPhoton;
    pht->SetParticleType(kNone);//kPhoton);
    pht->SetChi2(c1->Chi2());
    pht->SetPID(c1->GetPID());
    Float_t pos[3];
    TLorentzVector p1;
    c1->GetPosition(pos);
    c1->GetMomentum(p1, v);
    //TODO
    pht->SetPositionX(pos[0]);
    pht->SetPositionY(pos[1]);
    pht->SetPositionZ(pos[2]);
    pht->SetPxPyPzE( p1.Px(), p1.Py(), p1.Pz(), p1.E());
    pht->SetTrackDx( c1->GetTrackDx() );
    pht->SetTrackDz( c1->GetTrackDz() );
    pht->SetCharge(0);
    if(c1->IsEMCAL()) pht->SetCaloType(AliJPhoton::kEMCALCalo);
    if(c1->IsPHOS())  pht->SetCaloType(AliJPhoton::kPHOSCalo);
    pht->SetDistToBadChannel(c1->GetDistanceToBadChannel());
    pht->SetDispersion(c1->GetDispersion());
    pht->SetM20(c1->GetM20());
    pht->SetM02(c1->GetM02());
    pht->SetEmcCpvDist(c1->GetEmcCpvDistance());
    pht->SetNCells(c1->GetNCells());
    pht->SetCellsAmplitudeFraction( c1->GetCellsAmplitudeFraction() );
    pht->SetCellsAbsId(c1->GetCellsAbsId());
    Int_t imoduleID = GetSuperModuleNumber(c1->IsEMCAL(), c1->GetCellAbsId(0));
    pht->SetSuperModuleID(imoduleID);
  }
  delete fRecoUtils;
  delete caloClustersArr;

}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODCaloClusters(const AliAODEvent* aod)
{
  if( !aod ) aod=0;
  // Read the AliAODCaloClusters and fill the list of AliJPhoton containers
}

//______________________________________________________________________________
AliJEventHeader* AliJCORRANTask::ReadCommonHeader(AliVEvent *event){
  //Read the AliVEvent and fill the list of AliJEventHeader containers
  //create a header and fill it
  AliJEventHeader *hdr = new( (*fHeaderList)[fHeaderList->GetEntriesFast()] ) AliJEventHeader;

  AliVVZERO *v0 = event->GetVZEROData();
  if( v0 ) hdr->SetV0Mult(v0->GetMTotV0A() + v0->GetMTotV0C());
  // Get Centrality as a percent from 0% to 100%
  AliCentrality *cent = event->GetCentrality();
  if( cent ){
    hdr->SetCentrality( cent->GetCentralityPercentile("V0M"));
    hdr->SetCentralityArray(AliJEventHeader::kcV0M, cent->GetCentralityPercentile("V0M"));
    hdr->SetCentralityArray(AliJEventHeader::kcFMD, cent->GetCentralityPercentile("FMD"));
    hdr->SetCentralityArray(AliJEventHeader::kcTRK, cent->GetCentralityPercentile("TRK"));
    hdr->SetCentralityArray(AliJEventHeader::kcTKL, cent->GetCentralityPercentile("TKL"));
    hdr->SetCentralityArray(AliJEventHeader::kcCL0, cent->GetCentralityPercentile("CL0"));
    hdr->SetCentralityArray(AliJEventHeader::kcCL1, cent->GetCentralityPercentile("CL1"));
    hdr->SetCentralityArray(AliJEventHeader::kcV0MvsFMD, cent->GetCentralityPercentile("V0MvsFMD"));
    hdr->SetCentralityArray(AliJEventHeader::kcTKLvsV0, cent->GetCentralityPercentile("TKLvsV0"));
    hdr->SetCentralityArray(AliJEventHeader::kcZEMvsZDC, cent->GetCentralityPercentile("ZEMvsZDC"));
  }
  hdr->SetTriggerMaskAlice(event->GetTriggerMask()); //ULong64_t
  hdr->SetTriggerMaskJCorran(ConvertTriggerMask()); //UInt_t
  hdr->SetEventType(event->GetEventType());
  int ncontributors = 0;
  const AliVVertex * vtxESD = event->GetPrimaryVertex();
  if(vtxESD){
    hdr->SetXVertex(vtxESD->GetX()); //FK// EFF
    hdr->SetYVertex(vtxESD->GetY()); //FK// EFF
    hdr->SetZVertex(vtxESD->GetZ());
    //hdr->SetZVertexErr(vtxESD->GetZRes());
    double covMat[6];
    vtxESD->GetCovarianceMatrix(covMat);
    hdr->SetZVertexErr(TMath::Sqrt(covMat[5])); // GetZRes := TMath::Sqrt(fCovZZ)
    ncontributors = vtxESD->GetNContributors(); // get number of contributors to vertex 
    hdr->SetVtxMult( vtxESD->GetNContributors() );
  }else{
    hdr->SetZVertex(9999);
    hdr->SetZVertexErr(9999);
  }
  hdr->SetVtxMult(ncontributors); //FK// EFF contrib to vertex
  return hdr;
}
//______________________________________________________________________________
void AliJCORRANTask::ReadESDHeader(AliESDEvent *esd)
{
  // Read the AliESDEvent and fill the list of AliJEventHeader containers
  if(!esd) return;
  AliESDUtils::RefitESDVertexTracks( esd ); // TODO only for LHC11a right?
  AliJEventHeader *hdr = ReadCommonHeader( esd );
  //create a header and fill it
  AliMultiplicity *fSPDMult =(AliMultiplicity *) esd->GetMultiplicity();
  if(fSPDMult) hdr->SetSPDTrackletMult(fSPDMult->GetNumberOfTracklets());

  //TODO  Store Detector data
  if( fStoreEventPlaneSource ){
    *fVZEROData = *esd->GetVZEROData();
    *fTZEROData = AliESDTZERO(*esd->GetESDTZERO());
    *fZDCData  = *esd->GetESDZDC();
  }
  hdr->SetEventID( esd->GetEventNumberInFile());
  const AliESDVertex * vtxESD = esd->GetPrimaryVertex();
  if( vtxESD->GetStatus() == 0 ) hdr->SetVtxMult( 0 );
  // if fNcontributes > 0 then status is always true. do we need this?
}

//______________________________________________________________________________
void AliJCORRANTask::ReadAODHeader(AliAODEvent *aod)
{  
  //Read the AliAODEvent and fill the list of AliJEventHeader containers
  AliJEventHeader *hdr = ReadCommonHeader( aod );

  const AliAODTracklets *trackletsSPD = aod->GetTracklets();
  if(trackletsSPD){
    hdr->SetSPDTrackletMult(trackletsSPD->GetNumberOfTracklets());
  }
  //TODO hdr->SetEventID( esd->GetEventNumberInFile());
}

//______________________________________________________________________________
Int_t AliJCORRANTask::GetSuperModuleNumber(bool isemcal, Int_t absId)
{
  //get super module number 
  if(isemcal){
    //    return GetEMCALGeoUtils()->GetSuperModuleNumber(absId) ;
    return fEMCALGeoUtils->GetSuperModuleNumber(absId) ;

  } else {
    Int_t    relId[4];
    if ( absId >= 0) {
      fPHOSGeom->AbsToRelNumbering(absId,relId);
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

  return triggerMaskJC;
}


//______________________________________________________________________________

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
      AliJMCTrack *ctrack = new( (*fMCTrackList)[ntrack++] ) AliJMCTrack;;

      TParticle *partStack = stack->Particle(iTrack);
      Int_t   pdg  = partStack->GetPdgCode();
      // BS unused : Float_t engy = partStack->Energy();
      // BS unused : Float_t pt   = partStack->Pt();
      // BS unused : Float_t ptot = partStack->P();
      // BS unused : Float_t eta  = partStack->Eta();
      // BS unused : Float_t theta  = partStack->Theta();
      // BS unused : Float_t phi    = atan2(sin( partStack->Phi()), cos(partStack->Phi()));
      // BS unused : Short_t ch     = (Short_t) partStack->GetPDG()->Charge();
      Int_t label    = track->GetLabel();

      ctrack->SetLabel(label);
      ctrack->SetPdgCode(pdg);
      ctrack->SetPxPyPzE( partStack->Px(), partStack->Py(), partStack->Pz(), partStack->Energy());
      //ctrack->SetCharge(ch);
      ctrack->SetFlag(AliJMCTrack::kPrimary, isPrimary);

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

  /*
     if(fEsdTrackCuts->GetRequireTPCRefit()  && 
     ((aodTrack->GetStatus() & AliJTrack::kTPCrefit) == 0)) return kFALSE;
     if(fEsdTrackCuts->GetRequireITSRefit()  && 
     ((aodTrack->GetStatus() & AliJTrack::kITSrefit) == 0)) return kFALSE;
     */

  //            C u t s   o n    D C A
  Float_t impactDCA[3];
  if( aodTrack->GetPosition(impactDCA)){
    if((fEsdTrackCuts->GetMaxDCAToVertexXY()>0) && 
       (fEsdTrackCuts->GetMaxDCAToVertexXY() < sqrt(impactDCA[0]*impactDCA[0] + impactDCA[1]*impactDCA[1]))) return kFALSE;
    if((fEsdTrackCuts->GetMaxDCAToVertexZ()>0) &&
       (fEsdTrackCuts->GetMaxDCAToVertexZ() < TMath::Abs(impactDCA[2]))) return kFALSE;
  } else return kFALSE; 



  // if(fEsdTrackCuts->GetAcceptKinkDaughters()) //<--------how to check ?
  //  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);


  return kTRUE; 
}

bool AliJCORRANTask::SetAliceTriggerDef(AliJRunHeader *RunInfo){
  RunInfo->AddAliceTriggerDef( "kMB", AliVEvent::kMB );
  if( fRunType == "LHC10h" )
  {
    RunInfo->AddAliceTriggerDef( "kINT7", AliVEvent::kINT7 );
    RunInfo->AddAliceTriggerDef( "kMUON", AliVEvent::kMUON );
    RunInfo->AddAliceTriggerDef( "kHighMult", AliVEvent::kHighMult );
    RunInfo->AddAliceTriggerDef( "kEMC1", AliVEvent::kEMC1 );
    RunInfo->AddAliceTriggerDef( "kCINT5", AliVEvent::kCINT5 );
    RunInfo->AddAliceTriggerDef( "kCMUS5", AliVEvent::kCMUS5 );
    RunInfo->AddAliceTriggerDef( "kMUSH7", AliVEvent::kMUSH7 );
    RunInfo->AddAliceTriggerDef( "kMUL7", AliVEvent::kMUL7 );
    RunInfo->AddAliceTriggerDef( "kMUU7", AliVEvent::kMUU7 );
    RunInfo->AddAliceTriggerDef( "kEMC7", AliVEvent::kEMC7 );
    RunInfo->AddAliceTriggerDef( "kMUS7", AliVEvent::kMUS7 );
    RunInfo->AddAliceTriggerDef( "kPHI1", AliVEvent::kPHI1 );
    RunInfo->AddAliceTriggerDef( "kPHI7", AliVEvent::kPHI7 );
    RunInfo->AddAliceTriggerDef( "kUserDefined", AliVEvent::kUserDefined );
    RunInfo->AddAliceTriggerDef( "kFastOnly", AliVEvent::kFastOnly );
    RunInfo->AddAliceTriggerDef( "kAny", AliVEvent::kAny );
    RunInfo->AddAliceTriggerDef( "kAnyINT", AliVEvent::kAnyINT );
  }
  else{
    // Default
    RunInfo->AddAliceTriggerDef( "kINT7", AliVEvent::kINT7 );
    RunInfo->AddAliceTriggerDef( "kMUON", AliVEvent::kMUON );
    RunInfo->AddAliceTriggerDef( "kHighMult", AliVEvent::kHighMult );
    RunInfo->AddAliceTriggerDef( "kEMC1", AliVEvent::kEMC1 );
    RunInfo->AddAliceTriggerDef( "kCINT5", AliVEvent::kCINT5 );
    RunInfo->AddAliceTriggerDef( "kCMUS5", AliVEvent::kCMUS5 );
    RunInfo->AddAliceTriggerDef( "kMUSH7", AliVEvent::kMUSH7 );
    RunInfo->AddAliceTriggerDef( "kMUL7", AliVEvent::kMUL7 );
    RunInfo->AddAliceTriggerDef( "kMUU7", AliVEvent::kMUU7 );
    RunInfo->AddAliceTriggerDef( "kEMC7", AliVEvent::kEMC7 );
    RunInfo->AddAliceTriggerDef( "kMUS7", AliVEvent::kMUS7 );
    RunInfo->AddAliceTriggerDef( "kPHI1", AliVEvent::kPHI1 );
    RunInfo->AddAliceTriggerDef( "kPHI7", AliVEvent::kPHI7 );
    RunInfo->AddAliceTriggerDef( "kUserDefined", AliVEvent::kUserDefined );
    RunInfo->AddAliceTriggerDef( "kFastOnly", AliVEvent::kFastOnly );
    RunInfo->AddAliceTriggerDef( "kAny", AliVEvent::kAny );
    RunInfo->AddAliceTriggerDef( "kAnyINT", AliVEvent::kAnyINT );
  }
  return true;
}

bool AliJCORRANTask::SetAliceFilterMapDef(AliJRunHeader *RunInfo) {
  if( fRunType == "LHC10h" )
  {
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsL",BIT(0));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsITsa",BIT(1));
    RunInfo->AddAliceFilterMapDef("ItsStrong",BIT(2));
    RunInfo->AddAliceFilterMapDef("ElectronID",BIT(3));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsH",BIT(4));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsH2",BIT(5));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsH3",BIT(6));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsTPCOnly",BIT(7));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsRaa",BIT(8));
  }
  else
  {
    // Default
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsL",BIT(0));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsITsa",BIT(1));
    RunInfo->AddAliceFilterMapDef("ItsStrong",BIT(2));
    RunInfo->AddAliceFilterMapDef("ElectronID",BIT(3));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsH",BIT(4));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsH2",BIT(5));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsH3",BIT(6));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsTPCOnly",BIT(7));
    RunInfo->AddAliceFilterMapDef("EsdTrackCutsRaa",BIT(8));
  }
  return true;
}

void AliJCORRANTask::PrintOut() {
  AliJRunHeader * RunInfo = fAliRunHeader;
  cout << "===== TriggerDef =====" << endl;
  cout << RunInfo->GetAliceTriggerDef("kMB") << endl;
  cout << RunInfo->GetAliceTriggerDef("kINT7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kMUON") << endl;
  cout << RunInfo->GetAliceTriggerDef("kHighMult") << endl;
  cout << RunInfo->GetAliceTriggerDef("kEMC1") << endl;
  cout << RunInfo->GetAliceTriggerDef("kCINT5") << endl;
  cout << RunInfo->GetAliceTriggerDef("kCMUS5") << endl;
  cout << RunInfo->GetAliceTriggerDef("kMUSH7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kMUL7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kMUU7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kEMC7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kMUS7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kPHI1") << endl;
  cout << RunInfo->GetAliceTriggerDef("kPHI7") << endl;
  cout << RunInfo->GetAliceTriggerDef("kUserDefined") << endl;
  cout << RunInfo->GetAliceTriggerDef("kFastOnly") << endl;
  cout << RunInfo->GetAliceTriggerDef("kAny") << endl;
  cout << RunInfo->GetAliceTriggerDef("kAnyINT") << endl;
  cout << "===== FilterMapDef =====" << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsL") << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsITsa") << endl;
  cout << RunInfo->GetAliceFilterMapDef("ItsStrong") << endl;
  cout << RunInfo->GetAliceFilterMapDef("ElectronID") << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsH") << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsH2") << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsH3") << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsTPCOnly") << endl;
  cout << RunInfo->GetAliceFilterMapDef("EsdTrackCutsRaa") << endl;
}


//********************************************
//    UTILS
//********************************************
void AliJCORRANTask::AddListAODBranch(const char* aname, const char* cname, TClonesArray **obj, int nlist){
  *obj = new TClonesArray(cname, nlist);
  (*obj)->SetName(aname);
  AddAODBranch("TClonesArray", obj, fAODName.Data() );
}
