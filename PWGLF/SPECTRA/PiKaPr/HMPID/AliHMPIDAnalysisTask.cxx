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
// By means of AliHMPIDAnalysisTask.C macro it is possible to use this class
// to perform the analysis on local data, on data on alien using local machine
// and on CAF.

#ifndef AliHMPIDAnalysisTASK_CXX
#define AliHMPIDAnalysisTASK_CXX


#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "AliPIDResponse.h"
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
#include "AliHMPIDAnalysisTask.h"

ClassImp(AliHMPIDAnalysisTask)

//__________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask() :
  fESD(0x0),fMC(0x0),fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpNevents(0x0),
  fZvertex(0x0),
  fPIDResponse(0x0),
  fTrackCuts(0x0),
  fTrackFilter(0x0),
  fTree(0x0)
{
  //
  //Default ctor
  //
  for (Int_t i=0; i<69; i++) fVar[i]=0;
}

//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0x0), fMC(0x0), fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpNevents(0x0),
  fZvertex(0x0),
  fPIDResponse(0x0),
  fTrackCuts(0x0),
  fTrackFilter(0x0),
  fTree(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  for (Int_t i=0; i<69; i++) fVar[i]=0;

  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}

//___________________________________________________________________________
AliHMPIDAnalysisTask& AliHMPIDAnalysisTask::operator=(const AliHMPIDAnalysisTask& c)
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
    fPIDResponse     = c.fPIDResponse;
    fTrackCuts       = c.fTrackCuts;
    fTrackFilter     = c.fTrackFilter;
    fTree            = c.fTree;
    for (Int_t i=0; i<69; i++) fVar[i]=c.fVar[i];
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fMC(c.fMC),fUseMC(c.fUseMC),
  fHmpHistList(c.fHmpHistList),
  fHmpNevents(c.fHmpNevents),
  fZvertex(c.fZvertex),
  fPIDResponse(c.fPIDResponse),  
  fTrackCuts(c.fTrackCuts),
  fTrackFilter(c.fTrackFilter),
  fTree(c.fTree)
{
  //
  // Copy Constructor
  //
  for (Int_t i=0; i<69; i++) fVar[i]=c.fVar[i];
}
 
//___________________________________________________________________________
AliHMPIDAnalysisTask::~AliHMPIDAnalysisTask() {
  //
  //destructor
  //
  Info("~AliHMPIDAnalysisTask","Calling Destructor");
  if (fHmpHistList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHmpHistList;
}

//___________________________________________________________________________
void AliHMPIDAnalysisTask::ConnectInputData(Option_t *)
{
  // Connect ESD here

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
  } else
    fESD = esdH->GetEvent();
    fPIDResponse = esdH->GetPIDResponse();

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
  fTrackCuts->SetMinNClustersTPC(80);
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
void AliHMPIDAnalysisTask::UserExec(Option_t *)
{
  
  fHmpNevents->Fill(0);

  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  
  if(!vertex || !vertex->GetStatus() || vertex->GetNContributors()<1) {
    
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(!vertex) {
      PostData(1,fHmpHistList);
      PostData(2,fTree);
      return;
    }
    if(!vertex->GetStatus()){
      PostData(1,fHmpHistList);
      PostData(2,fTree);
      return;
    }  
    if(vertex->GetNContributors()<1){
      PostData(1,fHmpHistList);
      PostData(2,fTree);
      return;
    }    
  }
  
  Double_t vtxPos[3] = {999, 999, 999};
  if(vertex) vertex->GetXYZ(vtxPos);

  fHmpNevents->Fill(1);
  
  fZvertex->Fill(vtxPos[2]);
      
  if(!fPIDResponse){
      PostData(1,fHmpHistList);
      PostData(2,fTree);
      return;
   }
      
  Double_t probsHMP[5], probsTPC[5], probsTOF[5];
  Float_t  nSigmasTOF[5], nSigmasTPC[5];
      
  AliESDtrack *track = 0;

  Double_t ktol = 0.001;
  Double_t r[3], rin[3], rout[3];

  //
  // Main loop function, executed on Event basis
  //
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    track = fESD->GetTrack(iTrack);
    if(!track) continue;
    
    //if(!fTrackFilter->IsSelected(track)) continue; 
    //if(!track->GetStatus() || !AliESDtrack::kTOFout || !AliESDtrack::kTIME || track->GetTOFsignal()<12000 || track->GetTOFsignal()>200000 || track->GetIntegratedLength()<365) continue; 
     
    track->GetXYZ(r);
    track->GetInnerXYZ(rin);
    track->GetOuterXYZ(rout);

    if(Equal(track->GetHMPIDsignal(),-20.,ktol)) continue;
    if(track->GetHMPIDcluIdx() < 0) continue;

    Int_t q, nph;
    Float_t x, y;
    Float_t xpc, ypc, th, ph;
    track->GetHMPIDmip(x,y,q,nph);
    track->GetHMPIDtrk(xpc,ypc,th,ph);

    Int_t ITSrefit = track->IsOn(AliESDtrack::kITSrefit);
    Int_t TPCrefit = track->IsOn(AliESDtrack::kTPCrefit);
    
    if(Equal(x,0.,ktol) && Equal(y,0.,ktol) && Equal(xpc,0.,ktol) && Equal(ypc,0.,ktol)) continue;

    Double_t pHmp[3] = {0}, pHmp3 = 0;
    if (track->GetOuterHmpPxPyPz(pHmp)) pHmp3 = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]);

    Float_t b[2];
    Float_t bCov[3];
    track->GetImpactParameters(b,bCov);    
    
    track->GetHMPIDpid(probsHMP);
    fPIDResponse->ComputeTOFProbability(track,AliPID::kSPECIES,probsTOF);
    fPIDResponse->ComputeTPCProbability(track,AliPID::kSPECIES,probsTPC);
    
    for(Int_t ispecie=0; ispecie<AliPID::kSPECIES; ispecie++){
      nSigmasTOF[ispecie] = fPIDResponse->NumberOfSigmasTOF((AliVTrack*)track,(AliPID::EParticleType)ispecie);
      nSigmasTPC[ispecie] = fPIDResponse->NumberOfSigmasTPC((AliVTrack*)track,(AliPID::EParticleType)ispecie);      
    }
        
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
    fVar[13] = (Float_t)track->GetNcls(1);
    fVar[14] = (Float_t)probsHMP[0];
    fVar[15] = (Float_t)probsHMP[1];
    fVar[16] = (Float_t)probsHMP[2];
    fVar[17] = (Float_t)probsHMP[3];
    fVar[18] = (Float_t)probsHMP[4];
    fVar[19] = (Float_t)track->GetTOFsignal();
    fVar[20] = (Float_t)track->GetKinkIndex(0);
    fVar[21] = (Float_t)track->Xv();
    fVar[22] = (Float_t)track->Yv();
    fVar[23] = (Float_t)track->Zv();
    fVar[24] = (Float_t)track->GetTPCchi2();
    fVar[25] = b[0];
    fVar[26] = b[1];
    fVar[27] = track->GetHMPIDcluIdx()%1000000/1000;
    fVar[28] = vtxPos[0];
    fVar[29] = vtxPos[1];
    fVar[30] = vtxPos[2];
    fVar[31] = (Float_t)ITSrefit;
    fVar[32] = (Float_t)TPCrefit;
    fVar[33] = (Float_t)track->Eta();
    fVar[34] = (Float_t)r[0];
    fVar[35] = (Float_t)r[1];
    fVar[36] = (Float_t)r[2];
    fVar[37] = (Float_t)rout[0];           
    fVar[38] = (Float_t)rout[1];
    fVar[39] = (Float_t)rout[2];
    fVar[40] = (Float_t)pHmp[0];
    fVar[41] = (Float_t)pHmp[1];
    fVar[42] = (Float_t)pHmp[2];
    fVar[43] = (Float_t)track->Px();
    fVar[44] = (Float_t)track->Py();
    fVar[45] = (Float_t)track->Pz();
    fVar[46] = (Float_t)track->GetTPCClusterInfo(2,1);
    fVar[47] = (Float_t)track->GetTPCNclsF();
    fVar[48] = nSigmasTOF[0];
    fVar[49] = nSigmasTOF[1];
    fVar[50] = nSigmasTOF[2];
    fVar[51] = nSigmasTOF[3];
    fVar[52] = nSigmasTOF[4];
    fVar[53] = nSigmasTPC[0];
    fVar[54] = nSigmasTPC[1];
    fVar[55] = nSigmasTPC[2];
    fVar[56] = nSigmasTPC[3];
    fVar[57] = nSigmasTPC[4];
    fVar[58] = (Float_t)probsTOF[0];
    fVar[59] = (Float_t)probsTOF[1];
    fVar[60] = (Float_t)probsTOF[2];
    fVar[61] = (Float_t)probsTOF[3];
    fVar[62] = (Float_t)probsTOF[4];
    fVar[63] = (Float_t)probsTPC[0];
    fVar[64] = (Float_t)probsTPC[1];
    fVar[65] = (Float_t)probsTPC[2];
    fVar[66] = (Float_t)probsTPC[3];
    fVar[67] = (Float_t)probsTPC[4];
    fVar[68] = (Float_t)track->GetHMPIDchi2();
       
    fTree->Fill();
   
  }//track loop

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHmpHistList);
  PostData(2,fTree);
}
//___________________________________________________________________________
void AliHMPIDAnalysisTask::Terminate(Option_t*)
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
void AliHMPIDAnalysisTask::UserCreateOutputObjects() {
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
   fTree->Branch("NumTPCclust",&fVar[13]);
   fTree->Branch("ProbHMP0",&fVar[14]);
   fTree->Branch("ProbHMP1",&fVar[15]);
   fTree->Branch("ProbHMP2",&fVar[16]);
   fTree->Branch("ProbHMP3",&fVar[17]);
   fTree->Branch("ProbHMP4",&fVar[18]);
   fTree->Branch("TOFsignal",&fVar[19]);
   fTree->Branch("KinkIndex",&fVar[20]);
   fTree->Branch("Xv",&fVar[21]);
   fTree->Branch("Yv",&fVar[22]);
   fTree->Branch("Zv",&fVar[23]);
   fTree->Branch("TPCchi2",&fVar[24]);
   fTree->Branch("b0",&fVar[25]);
   fTree->Branch("b1",&fVar[26]);
   fTree->Branch("ClustSize",&fVar[27]);
   fTree->Branch("PrimVertexX",&fVar[28]);
   fTree->Branch("PrimVertexY",&fVar[29]);
   fTree->Branch("PrimVertexZ",&fVar[30]);
   fTree->Branch("ITSrefit",&fVar[31]);
   fTree->Branch("TPCrefit",&fVar[32]);   
   fTree->Branch("Eta",&fVar[33]);
   fTree->Branch("xTrack",&fVar[34]);
   fTree->Branch("yTrack",&fVar[35]);
   fTree->Branch("zTrack",&fVar[36]);
   fTree->Branch("xOuterTrack",&fVar[37]);
   fTree->Branch("yOuterTrack",&fVar[38]);
   fTree->Branch("zOuterTrack",&fVar[39]);
   fTree->Branch("pHmpX",&fVar[40]);
   fTree->Branch("pHmpY",&fVar[41]);
   fTree->Branch("pHmpZ",&fVar[42]);
   fTree->Branch("Px",&fVar[43]);
   fTree->Branch("Py",&fVar[44]);
   fTree->Branch("Pz",&fVar[45]);
   fTree->Branch("nCrossedRowsTPC",&fVar[46]);
   fTree->Branch("findableClustTPC",&fVar[47]);
   fTree->Branch("NSigmasTOF0",&fVar[48]);   
   fTree->Branch("NSigmasTOF1",&fVar[49]);   
   fTree->Branch("NSigmasTOF2",&fVar[50]);   
   fTree->Branch("NSigmasTOF3",&fVar[51]);   
   fTree->Branch("NSigmasTOF4",&fVar[52]);   
   fTree->Branch("NSigmasTPC0",&fVar[53]);   
   fTree->Branch("NSigmasTPC1",&fVar[54]);   
   fTree->Branch("NSigmasTPC2",&fVar[55]);   
   fTree->Branch("NSigmasTPC3",&fVar[56]);   
   fTree->Branch("NSigmasTPC4",&fVar[57]);   
   fTree->Branch("ProbTOF0",&fVar[58]);
   fTree->Branch("ProbTOF1",&fVar[59]);
   fTree->Branch("ProbTOF2",&fVar[60]);
   fTree->Branch("ProbTOF3",&fVar[61]);
   fTree->Branch("ProbTOF4",&fVar[62]);
   fTree->Branch("ProbTPC0",&fVar[63]);
   fTree->Branch("ProbTPC1",&fVar[64]);
   fTree->Branch("ProbTPC2",&fVar[65]);
   fTree->Branch("ProbTPC3",&fVar[66]);
   fTree->Branch("ProbTPC4",&fVar[67]);
   fTree->Branch("HmpSigma",&fVar[68]);
   
   PostData(1,fHmpHistList);
   PostData(2,fTree);
}

//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDAnalysisTask::Equal(Double_t x, Double_t y, Double_t tolerance)
{
 return abs(x - y) <= tolerance ;
}
   
#endif
