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
// AliHMPIDTaskQA - Class representing a quality check tool of HMPID 
// A set of histograms is created.
//==============================================================================


#ifndef AliHMPIDTASKQA_CXX
#define AliHMPIDTASKQA_CXX


#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliHMPIDTaskQA.h"

ClassImp(AliHMPIDTaskQA)

//__________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA() :
  fESD(0x0),fMC(0x0),fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpNevents(0x0),
  fZvertex(0x0),
  fTrackCuts(0x0),
  fTrackFilter(0x0),
  fTree(0x0)
{
  //
  //Default ctor
  //
  for (Int_t i=0; i<23; i++) fVar[i]=0;
}

//___________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0x0), fMC(0x0), fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpNevents(0x0),
  fZvertex(0x0),
  fTrackCuts(0x0),
  fTrackFilter(0x0),
  fTree(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  for (Int_t i=0; i<23; i++) fVar[i]=0;

  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}

//___________________________________________________________________________
AliHMPIDTaskQA& AliHMPIDTaskQA::operator=(const AliHMPIDTaskQA& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fESD             = c.fESD;
    fMC              = c.fMC;
    fUseMC           = c.fUseMC;
    fHmpHistList     = c.fHmpHistList;
    fHmpNevents      = c.fHmpNevents;
    fZvertex         = c.fZvertex;
    fTrackCuts       = c.fTrackCuts;
    fTrackFilter     = c.fTrackFilter;
    fTree            = c.fTree;
    for (Int_t i=0; i<23; i++) fVar[i]=c.fVar[i];
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDTaskQA::AliHMPIDTaskQA(const AliHMPIDTaskQA& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fMC(c.fMC),fUseMC(c.fUseMC),
  fHmpHistList(c.fHmpHistList),
  fHmpNevents(c.fHmpNevents),
  fZvertex(c.fZvertex),
  fTrackCuts(c.fTrackCuts),
  fTrackFilter(c.fTrackFilter),
  fTree(c.fTree)
{
  //
  // Copy Constructor
  //
  for (Int_t i=0; i<23; i++) fVar[i]=c.fVar[i];
}
 
//___________________________________________________________________________
AliHMPIDTaskQA::~AliHMPIDTaskQA() {
  //
  //destructor
  //
  Info("~AliHMPIDTaskQA","Calling Destructor");
  if (fHmpHistList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHmpHistList;
}

//___________________________________________________________________________
void AliHMPIDTaskQA::ConnectInputData(Option_t *)
{
  // Connect ESD here

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
  } else
    fESD = (AliESDEvent*)esdH->GetEvent();

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
  
  fTrackCuts = new AliESDtrackCuts("fTrackCuts", "Standard");
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetMinNClustersTPC(70);  // changed default value 80 -> 70 (mtangaro)
  fTrackCuts->SetMaxChi2PerClusterTPC(4);
  fTrackCuts->SetMaxDCAToVertexXY(3);
  fTrackCuts->SetMaxDCAToVertexZ(3);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetRequireITSRefit(kTRUE);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  fTrackFilter->AddCuts(fTrackCuts);   
    
}
//***************************************************************************************************************************************************************************
void AliHMPIDTaskQA::UserExec(Option_t *)
{
  
  fHmpNevents->Fill(0);

  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  
  if(!vertex || !vertex->GetStatus() || vertex->GetNContributors()<1) {
    
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(!vertex) return;
    if(!vertex->GetStatus()) return;
    if(vertex->GetNContributors()<1) return; // no good vertex, skip event
    
  }
  
  fHmpNevents->Fill(1);
  
  Double_t vtxPos[3] = {999, 999, 999};
  if(vertex) vertex->GetXYZ(vtxPos);
  fZvertex->Fill(vtxPos[2]);

  // vertex cut
  if (TMath::Abs(vtxPos[2]) > 10.) return;
  
  AliESDtrack *track = 0;

  Double_t ktol = 0.001;

  //
  // Main loop function, executed on Event basis
  //
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    track = fESD->GetTrack(iTrack);
    if(!track) continue;
    
    // pT min
    if ( track->Pt() < 0.2 ) continue;
    
    // TPC clust
    if( track->GetNcls(1) < 70 ) continue;

    // DCA cut
    Float_t b[2];
    Float_t bCov[3];
    track->GetImpactParameters(b,bCov);  
    if( TMath::Sqrt(b[0]*b[0] + b[1]*b[1]) > 2 ) continue;
    
    // ITS-TPC refit
    if( !track->IsOn(AliESDtrack::kITSrefit) ) continue;
    if( !track->IsOn(AliESDtrack::kTPCrefit) ) continue;
    
    //HMPID cuts 
    
    if(Equal(track->GetHMPIDsignal(),-20.,ktol)) continue;
    if(track->GetHMPIDcluIdx() < 0) continue;

    Int_t q, nph;
    Float_t x, y;
    Float_t xpc, ypc, th, ph;
    track->GetHMPIDmip(x,y,q,nph);
    track->GetHMPIDtrk(xpc,ypc,th,ph);
    
    if(Equal(x,0.,ktol) && Equal(y,0.,ktol) && Equal(xpc,0.,ktol) && Equal(ypc,0.,ktol)) continue;

    Double_t pHmp[3] = {0}, pHmp3 = 0;
    if (track->GetOuterHmpPxPyPz(pHmp)) pHmp3 = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]);
    
    fVar[0] = track->GetHMPIDcluIdx()/1000000;
    fVar[1] = pHmp3;
    fVar[2] = (Float_t)track->P();
    fVar[3] = xpc;
    fVar[4] = ypc;
    fVar[5] = x;
    fVar[6] = y;
    fVar[7] = (Float_t)track->GetHMPIDsignal();
    fVar[8] = q;
    fVar[9] = th;
    fVar[10] = ph;
    fVar[11] = (Float_t)track->GetSign();
    fVar[12] = (Float_t)nph;
    fVar[13] = track->GetHMPIDcluIdx()%1000000/1000; // cluster size
    fVar[14] = (Float_t)track->Eta();
    fVar[15] = (Float_t)track->Phi();
    fVar[16] = (Float_t)pHmp[0];
    fVar[17] = (Float_t)pHmp[1];
    fVar[18] = (Float_t)pHmp[2];
    fVar[19] = (Float_t)track->Px();
    fVar[20] = (Float_t)track->Py();
    fVar[21] = (Float_t)track->Pz();
    fVar[22] = (Float_t)track->GetHMPIDchi2();
       
    fTree->Fill();
   
  }//track loop

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHmpHistList);
  PostData(2,fTree);
}
//___________________________________________________________________________
void AliHMPIDTaskQA::Terminate(Option_t*)
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
void AliHMPIDTaskQA::UserCreateOutputObjects() {
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

   fZvertex = new TH1F("fZvertex","Z primary vertex distribution",4000,-20,20);
   fHmpHistList->Add(fZvertex);
   
//   OpenFile(2);
   fTree = new TTree("Tree","Tree with data");
   fTree->Branch("Chamber",&fVar[0]);
   fTree->Branch("pHmp3",&fVar[1]);
   fTree->Branch("P",&fVar[2]);
   fTree->Branch("Xpc",&fVar[3]);
   fTree->Branch("Ypc",&fVar[4]);
   fTree->Branch("X",&fVar[5]);
   fTree->Branch("Y",&fVar[6]);
   fTree->Branch("HMPIDsignal",&fVar[7]);
   fTree->Branch("Charge",&fVar[8]);
   fTree->Branch("Theta",&fVar[9]);
   fTree->Branch("Phi",&fVar[10]);
   fTree->Branch("Sign",&fVar[11]);
   fTree->Branch("NumPhotons",&fVar[12]);
   fTree->Branch("ClustSize",&fVar[13]);
   fTree->Branch("Eta",&fVar[14]);
   fTree->Branch("PhiTrack",&fVar[15]);
   fTree->Branch("pHmpX",&fVar[16]);
   fTree->Branch("pHmpY",&fVar[17]);
   fTree->Branch("pHmpZ",&fVar[18]);
   fTree->Branch("Px",&fVar[19]);
   fTree->Branch("Py",&fVar[20]);
   fTree->Branch("Pz",&fVar[21]);
   fTree->Branch("HmpSigma",&fVar[22]);
   
   PostData(1,fHmpHistList);
   PostData(2,fTree);
}

//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDTaskQA::Equal(Double_t x, Double_t y, Double_t tolerance)
{
  return TMath::Abs(x - y) <= tolerance ;
}
   
#endif
