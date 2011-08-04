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

//==============================================================================
// AliHMPIDAnalysysTask - Class representing a basic analysis tool of HMPID data
// A set of histograms is created.
//==============================================================================
//
// By means of AliHMPIDPerformanceTask.C macro it is possible to use this class
// to perform the analysis on local data, on data on alien using local machine
// and on CAF.

#ifndef AliHMPIDPERFORMANCETASK_CXX
#define AliHMPIDPERFORMANCETASK_CXX

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDfriend.h"
#include "AliHMPIDCluster.h"
#include "AliESDfriendTrack.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliHMPIDPerformanceTask.h"
#include "AliCentrality.h"

ClassImp(AliHMPIDPerformanceTask)

//__________________________________________________________________________
AliHMPIDPerformanceTask::AliHMPIDPerformanceTask() :
  fESD(0x0),fESDfriend(0x0),fEsdTrackCuts(0x0),fClu(0x0),fMC(0x0),fEsdVtx(0x0),fCentrality(0x0),fEsdTrack(0x0),fEsdFriendTrack(0x0),
  fCalibObject(0x0),  fUseMC(kTRUE),  fHmpHistList(0x0),  fHmpNevents(0x0),fHmpNevPerTrigClass(0x0),
  fGlobalEventNumber(0),
  //tree stuff
  fvType(0),
  fvFiredTriggerClasses(0x0), fvRunNumber(0), fvBunchCrossNumber(0),fvOrbitNumber(0),  fvPeriodNumber(0), fvMagField(0),
  fvVertexX(0),fvVertexY(0),fvVertexZ(0),fvVertexNContributors(0),  fvPesd(0),fvPhmpMag(0),fvCentrality(0),
  fvAcceptedTracks(0),fvRefMultTpc(0),  fvHmpChi2(0),fvHmpCluIndx(0),fvHmpMipX(0),fvHmpMipY(0),fvHmpMipQ(0),fvHmpMipNPhots(0),fvHmpSignal(0),
  fvHmpTrkX(0),fvHmpTrkY(0),fvHmpTrkTheta(0),fvHmpTrkPhi(0),  fvHmpCluQ(0), fvHmpCluX(0),fvHmpCluY(0),fvHmpCluCh(0),fvHmpCluSize(0),fvHmpCluBox(0),fvHmpCluStatus(0),
  fvEsdTrackAccepted(0),  fvKinkIndex(0),fvTofSignal(0),
  
  fTree(0x0)
{
  //Default ctor
  //
   for (Int_t i=0; i<3; i++) fvPhmp[i]=0;
   for (Int_t i=0; i<5; i++) fvHmpPid[i] = 0;
}

//___________________________________________________________________________
AliHMPIDPerformanceTask::AliHMPIDPerformanceTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0x0), fESDfriend(0x0), fEsdTrackCuts(0x0), fClu(0x0),fMC(0x0), fEsdVtx(0x0),fCentrality(0x0),
  fEsdTrack(0x0),fEsdFriendTrack(0x0),  fCalibObject(0x0),   fUseMC(kTRUE),
  fHmpHistList(0x0),  fHmpNevents(0x0),fHmpNevPerTrigClass(0x0),
  fGlobalEventNumber(0),  fvType(0),  fvFiredTriggerClasses(0x0), fvRunNumber(0), fvBunchCrossNumber(0),fvOrbitNumber(0),
  fvPeriodNumber(0),fvMagField(0),  fvVertexX(0),fvVertexY(0),fvVertexZ(0),fvVertexNContributors(0), 
  fvPesd(0),fvPhmpMag(0),fvCentrality(0),  fvAcceptedTracks(0),fvRefMultTpc(0),
  fvHmpChi2(0),fvHmpCluIndx(0),fvHmpMipX(0),fvHmpMipY(0),fvHmpMipQ(0),fvHmpMipNPhots(0),fvHmpSignal(0),
  fvHmpTrkX(0),fvHmpTrkY(0),fvHmpTrkTheta(0),fvHmpTrkPhi(0),  fvHmpCluQ(0), fvHmpCluX(0),fvHmpCluY(0),fvHmpCluCh(0),fvHmpCluSize(0),fvHmpCluBox(0),fvHmpCluStatus(0),
  fvEsdTrackAccepted(0),    fvKinkIndex(0),fvTofSignal(0),
  
  fTree(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  for (Int_t i=0; i<3; i++) fvPhmp[i]=0;
  for (Int_t i=0; i<5; i++) fvHmpPid[i] = 0;
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}

//___________________________________________________________________________
AliHMPIDPerformanceTask& AliHMPIDPerformanceTask::operator=(const AliHMPIDPerformanceTask& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fESD             = c.fESD;
    fESDfriend       = c.fESDfriend;
    fEsdTrackCuts    = c.fEsdTrackCuts;
    fClu             = c.fClu;
    fMC              = c.fMC;
    fEsdVtx			 = c.fEsdVtx;
    fCentrality       = c.fCentrality;
    fEsdTrack         = c.fEsdTrack;
    fEsdFriendTrack   =c.fEsdFriendTrack;
    fCalibObject      = c.fCalibObject;
          
    fUseMC           = c.fUseMC;
    fHmpHistList     = c.fHmpHistList;
    fHmpNevents      = c.fHmpNevents;
    fHmpNevPerTrigClass = c.fHmpNevPerTrigClass;
    fGlobalEventNumber = c.fGlobalEventNumber;
    
    fvFiredTriggerClasses =c.fvFiredTriggerClasses;
    fvRunNumber = c.fvRunNumber;
    fvBunchCrossNumber = c.fvBunchCrossNumber;
    fvOrbitNumber = c.fvOrbitNumber;
    fvPeriodNumber = c.fvPeriodNumber;
    fvMagField    =c.fvMagField;    
    fvVertexX     = c.fvVertexX;
    fvVertexY     = c.fvVertexY;
    fvVertexZ     = c.fvVertexZ;
    fvVertexNContributors = c.fvVertexNContributors;
    fvPesd = c.fvPesd;
    fvPhmpMag = c.fvPhmpMag;
    fvCentrality = c.fvCentrality;
    fvAcceptedTracks = c.fvAcceptedTracks;
    fvRefMultTpc = c.fvRefMultTpc;
    
    fvHmpChi2 = c.fvHmpChi2;
    fvHmpCluIndx = c.fvHmpCluIndx;
    fvHmpMipX = c.fvHmpMipX;
    fvHmpMipY = c.fvHmpMipY;
    fvHmpMipQ = c.fvHmpMipQ;
    fvHmpMipNPhots = c.fvHmpMipNPhots;
    fvHmpSignal = c.fvHmpSignal;
    fvHmpTrkX = c.fvHmpTrkX;
    fvHmpTrkY = c.fvHmpTrkY;
    fvHmpTrkTheta = c.fvHmpTrkTheta;
    fvHmpTrkPhi = c.fvHmpTrkPhi;
    
    fvHmpCluQ = c.fvHmpCluQ;
    fvHmpCluX = c.fvHmpCluX;
    fvHmpCluY = c.fvHmpCluY;
    fvHmpCluCh = c.fvHmpCluCh;
    fvHmpCluSize = c.fvHmpCluSize;
    fvHmpCluBox = c.fvHmpCluBox;
    fvHmpCluStatus = c.fvHmpCluStatus;
    fvEsdTrackAccepted = c.fvEsdTrackAccepted;           
    
    fvKinkIndex = c.fvKinkIndex;
    fvTofSignal = c.fvTofSignal;
    
    fTree            = c.fTree;

    for (Int_t i=0; i<3; i++) fvPhmp[i] = c.fvPhmp[i];
    for (Int_t i=0; i<5; i++) fvHmpPid[i] = c.fvHmpPid[i];
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDPerformanceTask::AliHMPIDPerformanceTask(const AliHMPIDPerformanceTask& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fESDfriend(c.fESDfriend),fEsdTrackCuts(c.fEsdTrackCuts),fClu(c.fClu),
  fMC(c.fMC),  fEsdVtx(c.fEsdVtx),    fCentrality(c.fCentrality),fEsdTrack(c.fEsdTrack),
  fEsdFriendTrack(c.fEsdFriendTrack),
  fCalibObject(c.fCalibObject),
    
  fUseMC(c.fUseMC),
  fHmpHistList(c.fHmpHistList),
  fHmpNevents(c.fHmpNevents),
  fHmpNevPerTrigClass(c.fHmpNevPerTrigClass),
  fGlobalEventNumber(c.fGlobalEventNumber),  
  fvFiredTriggerClasses(c.fvFiredTriggerClasses),
  fvRunNumber(c.fvRunNumber),
  fvBunchCrossNumber(c.fvBunchCrossNumber),
  fvOrbitNumber(c.fvOrbitNumber),
  fvPeriodNumber(c.fvPeriodNumber),
  fvMagField(c.fvMagField),  
  fvVertexX(c.fvVertexX),
  fvVertexY(c.fvVertexY),
  fvVertexZ(c.fvVertexZ),
  fvVertexNContributors(c.fvVertexNContributors),
  fvPesd(c.fvPesd),
  fvPhmpMag(c.fvPhmpMag), 
  fvCentrality(c.fvCentrality),
  fvAcceptedTracks(c.fvAcceptedTracks),
  fvRefMultTpc(c.fvRefMultTpc),
  fvHmpChi2 ( c.fvHmpChi2),
  fvHmpCluIndx ( c.fvHmpCluIndx),
  fvHmpMipX ( c.fvHmpMipX),
  fvHmpMipY ( c.fvHmpMipY),
  fvHmpMipQ ( c.fvHmpMipQ),
  fvHmpMipNPhots ( c.fvHmpMipNPhots),
  fvHmpSignal ( c.fvHmpSignal),
  fvHmpTrkX ( c.fvHmpTrkX),
  fvHmpTrkY ( c.fvHmpTrkY),
  fvHmpTrkTheta ( c.fvHmpTrkTheta),
  fvHmpTrkPhi ( c.fvHmpTrkPhi),
  
  fvHmpCluQ ( c.fvHmpCluQ),
  fvHmpCluX ( c.fvHmpCluX),
  fvHmpCluY ( c.fvHmpCluY),
  fvHmpCluCh ( c.fvHmpCluCh),
  fvHmpCluSize ( c.fvHmpCluSize),
  fvHmpCluBox ( c.fvHmpCluBox),
  fvHmpCluStatus ( c.fvHmpCluStatus),
  
  fvEsdTrackAccepted(c.fvEsdTrackAccepted),
  fvKinkIndex(c.fvKinkIndex),
  fvTofSignal(c.fvTofSignal),
  
  fTree(c.fTree)
{
  //
  // Copy Constructor
  //
  for (Int_t i=0; i<3; i++) fvPhmp[i] = c.fvPhmp[i];
  for (Int_t i=0; i<5; i++) fvHmpPid[i] = c.fvHmpPid[i];
}

//___________________________________________________________________________
AliHMPIDPerformanceTask::~AliHMPIDPerformanceTask() {
  //
  //destructor
  //
  Info("~AliHMPIDPerformanceTask","Calling Destructor");
  if (fHmpHistList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHmpHistList;
}

//___________________________________________________________________________
void AliHMPIDPerformanceTask::ConnectInputData(Option_t *)
{
  // Connect ESD here
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
  } else
    fESD = esdH->GetEvent();
  
  if (fUseMC){
    // Connect MC
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcH) {
      AliDebug(2,Form("ERROR: Could not get MCEventHandler"));
      fUseMC = kFALSE;
    } else
      fMC = mcH->MCEvent();
    if (!fMC) AliDebug(2,Form("ERROR: Could not get MCEvent"));
  }
  
  //Init ESD track cuts
  fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE); 
  
}
//___________________________________________________________________________
void AliHMPIDPerformanceTask::UserExec(Option_t *)
{  
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  fHmpNevents->Fill(0);
  
  //fEsdVtx = fESD->GetPrimaryVertex();
  
  /* check which trigger is fired */          fHmpNevPerTrigClass->Fill(0);
  if(fESD->IsTriggerClassFired("kMB"))        fHmpNevPerTrigClass->Fill(1);
  if(fESD->IsTriggerClassFired("kINT7"))      fHmpNevPerTrigClass->Fill(2);
  if(fESD->IsTriggerClassFired("kCINT5"))     fHmpNevPerTrigClass->Fill(3);
  if(fESD->IsTriggerClassFired("kFastOnly"))  fHmpNevPerTrigClass->Fill(4);
  if(fESD->IsTriggerClassFired("kHighMult"))  fHmpNevPerTrigClass->Fill(5);
  if(fESD->IsTriggerClassFired("kAnyINT"))    fHmpNevPerTrigClass->Fill(6);
  
  fvFiredTriggerClasses = Form("%s",(fESD->GetFiredTriggerClasses()).Data());
  
  //___ check vertex
  if(!fESD->GetPrimaryVertex() || !fESD->GetPrimaryVertex() || fESD->GetPrimaryVertex()->GetNContributors()<1) {     // SPD vertex
    if(!fESD->GetPrimaryVertexSPD()) return;
    if(!fESD->GetPrimaryVertexSPD()->GetStatus()) return;
    if(fESD->GetPrimaryVertexSPD()->GetNContributors()<1) return; // no good vertex, skip event  
  }
  
  fHmpNevents->Fill(1);
  if(fESD->GetPrimaryVertex()) {
    fvVertexX = fESD->GetPrimaryVertex()->GetX();
    fvVertexY = fESD->GetPrimaryVertex()->GetY();
    fvVertexZ = fESD->GetPrimaryVertex()->GetZ();
    fvVertexNContributors = fESD->GetPrimaryVertex()->GetNContributors();
    fvMagField = fESD->GetMagneticField();
  }
  
  fCentrality = fESD->GetCentrality();
  if(fCentrality)
    {
      fvCentrality = fCentrality->GetCentralityPercentile("V0M");
    }
  
  // Exec function
  // Loop over tracks and call  Process function
  
  
  if(fEsdTrackCuts) 
    {
      fvAcceptedTracks  =  fEsdTrackCuts->CountAcceptedTracks(fESD);
     fvRefMultTpc         =  fEsdTrackCuts->GetReferenceMultiplicity(fESD, kTRUE);
    }
  
  
  
  fESDfriend=(AliESDfriend*)(((AliESDEvent*)fESD)->FindListObject("AliESDfriend"));
  
  //  if (!fESDfriend) { Printf("ERROR: fESDfriend not available");   }
  
  Double_t ktol = 0.001;
  
  fvRunNumber = fESD->GetRunNumber();
  fvBunchCrossNumber = fESD->GetBunchCrossNumber();
  fvOrbitNumber = fESD->GetOrbitNumber();
  fvPeriodNumber = fESD->GetPeriodNumber();
  
  fGlobalEventNumber =   (ULong64_t)fESD->GetBunchCrossNumber()     + (ULong64_t)fESD->GetOrbitNumber()*3564+     (ULong64_t)fESD->GetPeriodNumber()*16777215*3564;
  
  fEsdTrackCuts -> CountAcceptedTracks(fESD);
  
  
  for (Int_t iTrack=0; iTrack<fESD->GetNumberOfTracks(); iTrack++) { // tracks loop    
    
    fEsdTrack=fESD->GetTrack(iTrack);
    if(!fEsdTrack) continue;
    
    if(fEsdTrackCuts && fEsdTrackCuts->AcceptTrack(fEsdTrack) == kTRUE ) fvEsdTrackAccepted = 1;
    else fvEsdTrackAccepted = 0;
    
    
    if(Equal(fEsdTrack->GetHMPIDsignal(),-20.,ktol)) continue;
    if(fEsdTrack->GetHMPIDcluIdx() < 0) continue;
    
    fvPesd = fEsdTrack->P();

    //fEsdTrack->GetOuterHmpPxPyPz(fvPhmp);     
    fvPhmpMag = TMath::Sqrt(fvPhmp[0]*fvPhmp[0] +  fvPhmp[1]*fvPhmp[1] + fvPhmp[2]*fvPhmp[2] );

    fvTofSignal = fEsdTrack->GetTOFsignal();
    
    fEsdFriendTrack=fESDfriend->GetTrack(iTrack);
    
    //___ fill HMP stuff
    fvHmpChi2 = fEsdTrack -> GetHMPIDchi2();
    fvHmpCluIndx = fEsdTrack -> GetHMPIDcluIdx();
    fEsdTrack->GetHMPIDmip(fvHmpMipX,fvHmpMipY,fvHmpMipQ,fvHmpMipNPhots);
    
    //oid	GetHMPIDpid(Double_t* p) const
    fvHmpSignal= fEsdTrack->GetHMPIDsignal();
    fEsdTrack->GetHMPIDtrk(fvHmpTrkX,fvHmpTrkY,fvHmpTrkTheta, fvHmpTrkPhi);
    fEsdTrack->GetHMPIDpid(fvHmpPid);  
    
    fvKinkIndex = fEsdTrack->GetKinkIndex(0);    
    
    fvType = 0;    fTree->Fill();
    if(fEsdFriendTrack) {
      for(Int_t j=0;(    fCalibObject=fEsdFriendTrack->GetCalibObject(j));++j){ // calib object loop
	if((fClu=dynamic_cast<AliHMPIDCluster*>(fCalibObject))){ // if HMPID calib object
	  fvHmpCluQ = fClu->Q();
	  fvHmpCluX = fClu->X();
	  fvHmpCluY = fClu->Y();
	  fvHmpCluCh = fClu->Ch();         
	  fvHmpCluSize = fClu->Size();         
	  fvHmpCluBox = fClu->Box();         
	  fvHmpCluStatus = fClu->Status();
	  
	  fvType = 1;  fTree->Fill();
	  
	} // if HMPID calib object 
      } // calib object loop 
    }//there is fEsdFriendTrack
    
    
    
  } // tracks loop
  
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHmpHistList);
  PostData(2,fTree);
}
//___________________________________________________________________________
void AliHMPIDPerformanceTask::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  Info("Terminate"," ");
  
  if (!fUseMC) return;
  
  fHmpHistList = dynamic_cast<TList*> (GetOutputData(1));
  
  if (!fHmpHistList) {
    AliError("Histogram List is not available");
    return;
  }
  
  AliAnalysisTaskSE::Terminate();
  
}

//___________________________________________________________________________
void AliHMPIDPerformanceTask::UserCreateOutputObjects() {
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //

  //slot #1
//   OpenFile(1);
   fHmpHistList = new TList();
   fHmpHistList->SetOwner();

   fHmpNevents = new TH1F("fHmpNevents","Number of events",2,0,2);
   fHmpHistList->Add(fHmpNevents);

   fHmpNevPerTrigClass = new TH1I("fHmpNevPerTrigClass","Number of event per trigger class;Trigger class;Events",10,0,10);
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(1,"All coll.");
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(2,"MB");
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(3,"INT7");
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(4,"CINT5");
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(5,"FastOnly");
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(6,"HighMult");
   fHmpNevPerTrigClass->GetXaxis()->SetBinLabel(7,"AnyINT");
   fHmpHistList->Add(fHmpNevPerTrigClass);
 
   
 //   OpenFile(2);
   fTree = new TTree("TreePhot","Tree with photon clusters information");
   fTree->Branch("fvFiredTriggerClasses",fvFiredTriggerClasses,"fvFiredTriggerClasses/C");
   fTree->Branch("fvRunNumber",&fvRunNumber,"fvRunNumber/I");
   fTree->Branch("fvBunchCrossNumber",&fvBunchCrossNumber,"fvBunchCrossNumber/I");
   fTree->Branch("fvOrbitNumber",&fvOrbitNumber,"fvOrbitNumber/I");
   fTree->Branch("fvPeriodNumber",&fvPeriodNumber,"fvPeriodNumber/I");
   fTree->Branch("fvMagField",&fvMagField,"fvMagField/D");
   fTree->Branch("fvVertexX",&fvVertexX,"fvVertexX/D");
   fTree->Branch("fvVertexY",&fvVertexY,"fvVertexY/D");
   fTree->Branch("fvVertexZ",&fvVertexZ,"fvVertexZ/D");
   fTree->Branch("fvVertexNContributors",&fvVertexNContributors,"fvVertexNContributors/I");
   fTree->Branch("fvCentrality",&fvCentrality,"fvCentrality/D");
   fTree->Branch("fvAcceptedTracks",&fvAcceptedTracks,"fvAcceptedTracks/I");
   fTree->Branch("fvRefMultTpc",&fvRefMultTpc,"fvRefMultTpc/I");
   fTree->Branch("fvHmpChi2",&fvHmpChi2,"fvHmpChi2/D");
   fTree->Branch("fvHmpCluIndx",&fvHmpCluIndx,"fvHmpCluIndx/I");
   fTree->Branch("fvHmpMipX",&fvHmpMipX,"fvHmpMipX/F");
   fTree->Branch("fvHmpMipY",&fvHmpMipY,"fvHmpMipY/F");
   fTree->Branch("fvHmpMipQ",&fvHmpMipQ,"fvHmpMipQ/I");
   fTree->Branch("fvHmpMipNPhots",&fvHmpMipNPhots,"fvHmpMipNPhots/I");
   fTree->Branch("fvHmpSignal",&fvHmpSignal,"fvHmpSignal/D");
   fTree->Branch("fvHmpTrkX",&fvHmpTrkX,"fvHmpTrkX/F");
   fTree->Branch("fvHmpTrkY",&fvHmpTrkY,"fvHmpTrkY/F");
   fTree->Branch("fvHmpTrkTheta",&fvHmpTrkTheta,"fvHmpTrkTheta/F");
   fTree->Branch("fvHmpTrkPhi",&fvHmpTrkPhi,"fvHmpTrkPhi/F");
   fTree->Branch("fvHmpPid",fvHmpPid,"fvHmpPid/D");
   fTree->Branch("fvType",&fvType,"fvType/I");
   fTree->Branch("fvHmpCluQ",&fvHmpCluQ,"fvHmpCluQ/D");
   fTree->Branch("fvHmpCluX",&fvHmpCluX,"fvHmpCluX/D");
   fTree->Branch("fvHmpCluY",&fvHmpCluY,"fvHmpCluY/D");
   fTree->Branch("fvHmpCluCh",&fvHmpCluCh,"fvHmpCluCh/I");
   fTree->Branch("fvHmpCluSize",&fvHmpCluSize,"fvHmpCluSize/I");
   fTree->Branch("fvHmpCluBox",&fvHmpCluBox,"fvHmpCluBox/I");
   fTree->Branch("fvHmpCluStatus",&fvHmpCluStatus,"fvHmpCluStatus/I");
   fTree->Branch("fvEsdTrackAccepted",&fvEsdTrackAccepted,"fvEsdTrackAccepted/I");
   fTree->Branch("fvKinkIndex",&fvKinkIndex,"fvKinkIndex/I");
   fTree->Branch("fvTofSignal",&fvTofSignal,"fvTofSignal/D");  
    
   PostData(1,fHmpHistList);
   PostData(2,fTree);
}
//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDPerformanceTask::Equal(Double_t x, Double_t y, Double_t tolerance)
{
 return abs(x - y) <= tolerance ;
}
//____________________________________________________________________________________________________________________________________
   
#endif
