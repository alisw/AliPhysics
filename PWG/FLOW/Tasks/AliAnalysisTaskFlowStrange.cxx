/*************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use,copy,modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee,provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowStrange:
// Analysis task to select K0/Lambda candidates for flow analysis.
// Uses one AliESDtrackCuts object for both daughters and
// QA histograms to monitor the reconstruction.
// Authors: Cristian Ivan (civan@cern.ch)
//          Carlos Perez (cperez@cern.ch)
//          Pawel Debski (pdebski@cern.ch)
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TVector3.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"

#include "AliAODMCParticle.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "TMath.h"
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowBayesianPID.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowVector.h"

#include "AliAnalysisTaskFlowStrange.h"

ClassImp(AliAnalysisTaskFlowStrange)

//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange() :
  AliAnalysisTaskSE(),
  fPIDResponse(NULL),
  fBayesianPID(NULL),
  fDebug(kFALSE),
  fUseEventSelection(kTRUE),
  fDoQA(kFALSE),
  fDoExtraQA(kFALSE),
  fRunOnpA(kFALSE),
  fRunOnpp(kFALSE),
  fCalib(NULL),
  fPsi2(0.0),
  fSpecie(0),
  fMCmatch(-1),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fCutsEvent(NULL),
  fCutsRFPTPC(NULL),
  fCutsRFPVZE(NULL),
  fCutsPOI(NULL),
  fCutsDau(NULL),
  fFlowEventTPC(NULL),
  fFlowEventVZE(NULL),
  fCandidates(NULL),
  fQAList(NULL)
{
  //ctor
  for (Int_t i=0; i!=11; ++i)
    fV0Cuts[i] = 0;
}
//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange(const char *name,
						       AliFlowEventCuts *cutsEvent,
						       AliFlowTrackCuts *cutsRFPTPC,
						       AliFlowTrackCuts *cutsRFPVZE,
						       AliESDtrackCuts *cutsDau) :
  AliAnalysisTaskSE(name),
  fPIDResponse(NULL),
  fBayesianPID(NULL),
  fDebug(kFALSE),
  fUseEventSelection(kTRUE),
  fDoQA(kFALSE),
  fDoExtraQA(kFALSE),
  fRunOnpA(kFALSE),
  fRunOnpp(kFALSE),
  fCalib(NULL),
  fPsi2(0.0),
  fSpecie(0),
  fMCmatch(-1),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fCutsEvent(cutsEvent),
  fCutsRFPTPC(cutsRFPTPC),
  fCutsRFPVZE(cutsRFPVZE),
  fCutsPOI(NULL),
  fCutsDau(cutsDau),
  fFlowEventTPC(NULL),
  fFlowEventVZE(NULL),
  fCandidates(NULL),
  fQAList(NULL)
{
  //ctor
  for (Int_t i=0; i!=11; ++i)
    fV0Cuts[i] = 0;

  DefineInput( 0,TChain::Class());
  DefineOutput(1,AliFlowEventSimple::Class()); // TPC object
  DefineOutput(2,AliFlowEventSimple::Class()); // VZE object
  DefineOutput(3,TList::Class());
}
//=======================================================================
AliAnalysisTaskFlowStrange::~AliAnalysisTaskFlowStrange()
{
  //dtor
  if (fQAList)       delete fQAList;
  if (fFlowEventTPC) delete fFlowEventTPC;
  if (fFlowEventVZE) delete fFlowEventVZE;
  if (fCandidates)   delete fCandidates;
  if (fCutsDau)      delete fCutsDau;
  if (fCutsPOI)      delete fCutsPOI;
  if (fCutsRFPTPC)   delete fCutsRFPTPC;
  if (fCutsRFPVZE)   delete fCutsRFPVZE;
  if (fCalib)        delete fCalib;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserCreateOutputObjects()
{
  //UserCreateOutputObjects
  fQAList=new TList();
  fQAList->SetOwner();
  AddQAEvents();
  AddQACandidates();

  AliFlowCommonConstants *cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(1); cc->SetMultMin(0);   cc->SetMultMax(1);
  cc->SetNbinsPt(120); cc->SetPtMin(0.0);   cc->SetPtMax(12.0);
  cc->SetNbinsPhi(1);  cc->SetPhiMin(0.0);  cc->SetPhiMax(TMath::TwoPi());
  cc->SetNbinsEta(18);  cc->SetEtaMin(-0.9); cc->SetEtaMax(+0.9);
  cc->SetNbinsQ(1);    cc->SetQMin(0.0);    cc->SetQMax(1.0);
  cc->SetNbinsMass(fMassBins);
  cc->SetMassMin(fMinMass);
  cc->SetMassMax(fMaxMass);

  fCutsPOI = new AliFlowTrackCuts("dumb_cuts");
  fCutsPOI->SetParamType( fCutsRFPTPC->GetParamType() );
  fCutsPOI->SetPtRange(+1.0,-1.0);
  fCutsPOI->SetEtaRange(+1.0,-1.0);

  fBayesianPID = new AliFlowBayesianPID();
  fBayesianPID->SetNewTrackParam();
  fFlowEventTPC = new AliFlowEvent(3000);
  fFlowEventVZE = new AliFlowEvent(500);
  fCandidates = new TObjArray(100);
  fCandidates->SetOwner();

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  PostData(1,fFlowEventTPC);
  PostData(2,fFlowEventVZE);
  PostData(3,fQAList);

}
//=======================================================================
void AliAnalysisTaskFlowStrange::Exec(Option_t* option)
{
  //bypassing ::exec (needed because of AMPT)
  if(fMCmatch==-1) // executes EVENT in data analysis
    AliAnalysisTaskSE::Exec(option);
  else // executes EVENT in monteCarlo
    AliAnalysisTaskFlowStrange::MyUserExec(option);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserExec(Option_t *option)
{
  // dummy user exec
  if(fRunOnpA) { // temporal extra cuts for pA (2BE REMOVED!)
    AliAODEvent *aod=dynamic_cast<AliAODEvent*>(InputEvent());
    if(!aod) return;
    if(aod->GetHeader()->GetEventNumberESDFile() == 0) return; //rejecting first chunk
     // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAVertexSelectionStudies
    const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
    if (!trkVtx || trkVtx->GetNContributors()<=0) return;
    TString vtxTtl = trkVtx->GetTitle();
    if (!vtxTtl.Contains("VertexerTracks")) return;
    Float_t zvtx = trkVtx->GetZ();
    const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
    if (spdVtx->GetNContributors()<=0) return;
    TString vtxTyp = spdVtx->GetTitle();
    Double_t cov[6]={0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
    if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
    if (TMath::Abs(zvtx) > 10) return;
  }

  AliAnalysisTaskFlowStrange::MyUserExec(option);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQAEvents()
{
  // function to add event qa
  TList *tQAEvents=new TList();
  tQAEvents->SetName("Event");
  tQAEvents->SetOwner();

  TH1D *tEvent = new TH1D("Events","Number of Events",3,0,3); tQAEvents->Add(tEvent);
  tEvent->GetXaxis()->SetBinLabel(1,"reached");
  tEvent->GetXaxis()->SetBinLabel(2,"selected");
  tEvent->GetXaxis()->SetBinLabel(3,"unexpected");

  TH2D *tH2D;
  TH1D *tTPCRFP = new TH1D("RFPTPC","TPC Reference Flow Particles;multiplicity",3000,0,3000); tQAEvents->Add(tTPCRFP);
  TH1D *tVZERFP = new TH1D("RFPVZE","VZERO Reference Flow Particles;multiplicity",3000,0,30000); tQAEvents->Add(tVZERFP);
  tH2D = new TH2D("TPCPhiEta","TPC RFP;Phi;Eta",100,0,TMath::TwoPi(),100,-1.0,+1.0); tQAEvents->Add( tH2D );
  tH2D = new TH2D("VZEPhiEta","VZE RFP;Phi;Eta",20,0,TMath::TwoPi(),40,-4.0,+6.0); tQAEvents->Add( tH2D );
  TH1D *tPOI = new TH1D("POI","POIs;multiplicity",500,0,500); tQAEvents->Add( tPOI );
  if(fDoQA) {
    printf("QA enabled\n");
    tH2D = new TH2D("VTXZ","VTXZ;Global||SPD;SPD",60,-25,+25,60,-25,+25); tQAEvents->Add( tH2D );
    TH3D *tH3D = new TH3D("EVPLANE","EVPLANE;TPC;V0A;V0C",72,0,TMath::Pi(),72,0,TMath::Pi(),72,0,TMath::Pi()); tQAEvents->Add( tH3D );
    tH2D = new TH2D("VTXZSEL","VTXZ SEL;Global||SPD;SPD",40,-10,+10,40,-10,+10); tQAEvents->Add( tH2D );
    tH3D = new TH3D("PRIMVERTEX","PRIMVERTEX;#sigma_{x};#sigma_{y};#sigma_{z}",100,0,5e-3,100,0,5e-3,100,0,8e-3); tQAEvents->Add( tH3D );
    double dArrayPt[25] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
			   1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
			   2.5, 3.0, 3.5, 4.0, 5.0 };
    tH2D = new TH2D("TPCPtQx","TPCPtQx;pT;Qx/M",24,dArrayPt,100,-0.2,+0.2); tQAEvents->Add( tH2D );
    tH2D = new TH2D("TPCPtQy","TPCPtQy;pT;Qy/M",24,dArrayPt,100,-0.2,+0.2); tQAEvents->Add( tH2D );
    tH2D = new TH2D("TPCPtEta","TPCPtEta;pT;Eta/M",24,dArrayPt,100,-0.3,+0.3); tQAEvents->Add( tH2D );
    tH2D = new TH2D("TPCQxQy","TPCQxQy;Qx/M;Qy/M",100,-0.3,+0.3,100,-0.3,+0.3); tQAEvents->Add( tH2D );
    tH2D = new TH2D("VZEQxQy","VZEQxQy;Qx/M;Qy/M",100,-0.3,+0.3,100,-0.3,+0.3); tQAEvents->Add( tH2D );
  }
  TProfile *tCuts = new TProfile("Cuts","Analysis Cuts",11,0,11);
  tCuts->Fill(0.5,fV0Cuts[0],1); tCuts->GetXaxis()->SetBinLabel(1,"dl");
  tCuts->Fill(1.5,fV0Cuts[1],1); tCuts->GetXaxis()->SetBinLabel(2,"dca");
  tCuts->Fill(2.5,fV0Cuts[2],1); tCuts->GetXaxis()->SetBinLabel(3,"ctp");
  tCuts->Fill(3.5,fV0Cuts[3],1); tCuts->GetXaxis()->SetBinLabel(4,"d0");
  tCuts->Fill(4.5,fV0Cuts[4],1); tCuts->GetXaxis()->SetBinLabel(5,"d0xd0");
  tCuts->Fill(5.5,fV0Cuts[5],1); tCuts->GetXaxis()->SetBinLabel(6,"qt");
  tCuts->Fill(6.5,fV0Cuts[6],1); tCuts->GetXaxis()->SetBinLabel(7,"min eta");
  tCuts->Fill(7.5,fV0Cuts[7],1); tCuts->GetXaxis()->SetBinLabel(8,"max eta");
  tCuts->Fill(8.5,fV0Cuts[8],1); tCuts->GetXaxis()->SetBinLabel(9,"PID sigmas");
  tCuts->Fill(9.5,fV0Cuts[9],1); tCuts->GetXaxis()->SetBinLabel(10,"ctau");
  tCuts->Fill(10.5,fV0Cuts[10],1); tCuts->GetXaxis()->SetBinLabel(11,"dlxy");
  tQAEvents->Add(tCuts);
  fQAList->Add(tQAEvents);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQACandidates()
{
  // function to add histogramming for candidates
  TList *tList;
  TH1D *tH1D;
  TH2D *tH2D;
  TH3D *tH3D;

  Int_t nMass = 88;
  Double_t dMinMass = 0.412, dMaxMass=0.588;
  switch(fSpecie) {
  case(1):
    nMass = 92;
    dMinMass = 1.075;
    dMaxMass = 1.167;
    break;
  case(90):
    nMass = 100;
    dMinMass = 0.0;
    dMaxMass = 1.0;
    break;
  }
  double dArrayPt[29] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
			 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
			 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10., 12. };
  tList=new TList();
  tList->SetName("Candidates");
  tList->SetOwner();
  tH2D = new TH2D("V0MASS","V0MASS;pT;Mass",28,dArrayPt,nMass,dMinMass,dMaxMass); tList->Add(tH2D);
  tH2D = new TH2D("V0PhiEta","V0PhiEta;Phi;Eta",100,0,TMath::TwoPi(),100,-1.0,+1.0); tList->Add(tH2D);
  fQAList->Add(tList);

  if(fDoExtraQA) {
    printf("Extra QA enabled\n");
    tList = new TList(); tList->SetOwner(); tList->SetName("QACutsBefore_IP");
    tH3D = new TH3D("BefVOL", "VOLUME;Phi;Eta;Pt [GeV]",       63, 0.0,+6.3, 40,-1.0,+1.0, 60,0,12); tList->Add( tH3D );
    tH3D = new TH3D("BefRAP", "RAPIDITY;y;Pt [GeV]",           40,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefDLZ", "DLZ;[cm];Pt [GeV]",             50,-5.0,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefDLXY","DLXY;[cm];Pt [GeV]",           100,0.00,100., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefCTau","CTau;[cm];Pt [GeV]",           250,0.00,250., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefDCA", "DCA;[cm];Pt [GeV]",             50,0.00,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefCTP", "CTP;;Pt [GeV]",                 80,0.99,1.00, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefD0",  "D0;[cm];Pt [GeV]",              50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefD0D0","D0D0;[cm^{2}];Pt [GeV]",        50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefAP",  "AP;#alpha;q_{t}[GeV];Pt [GeV]", 80,-1.0,+1.0, 90,0,0.3, 60,0,12); tList->Add( tH3D );
    fQAList->Add(tList);
    tList = new TList(); tList->SetOwner(); tList->SetName("QACutsBefore_OP");
    tH3D = new TH3D("BefVOL", "VOLUME;Phi;Eta;Pt [GeV]",       63, 0.0,+6.3, 40,-1.0,+1.0, 60,0,12); tList->Add( tH3D );
    tH3D = new TH3D("BefRAP", "RAPIDITY;y;Pt [GeV]",           40,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefDLZ", "DLZ;[cm];Pt [GeV]",             50,-5.0,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefDLXY","DLXY;[cm];Pt [GeV]",           100,0.00,100., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefCTau","CTau;[cm];Pt [GeV]",           250,0.00,250., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefDCA", "DCA;[cm];Pt [GeV]",             50,0.00,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefCTP", "CTP;;Pt [GeV]",                 50,0.99,1.00, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefD0",  "D0;[cm];Pt [GeV]",              50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefD0D0","D0D0;[cm^{2}];Pt [GeV]",        50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("BefAP",  "AP;#alpha;q_{t}[GeV];Pt [GeV]", 80,-1.0,+1.0, 90,0,0.3, 60,0,12); tList->Add( tH3D );
    fQAList->Add(tList);
    tList = new TList(); tList->SetOwner(); tList->SetName("QACutsAfter_IP");
    tH3D = new TH3D("AftVOL", "VOLUME;Phi;Eta;Pt [GeV]",       63, 0.0,+6.3, 40,-1.0,+1.0, 60,0,12); tList->Add( tH3D );
    tH3D = new TH3D("AftRAP", "RAPIDITY;y;Pt [GeV]",           40,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftDLZ", "DLZ;[cm];Pt [GeV]",             50,-5.0,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftDLXY","DLXY;[cm];Pt [GeV]",           100,0.00,100., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftCTau","CTau;[cm];Pt [GeV]",           250,0.00,250., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftDCA", "DCA;[cm];Pt [GeV]",             50,0.00,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftCTP", "CTP;;Pt [GeV]",                 50,0.99,1.00, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftD0",  "D0;[cm];Pt [GeV]",              50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftD0D0","D0D0;[cm^{2}];Pt [GeV]",        50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftAP",  "AP;#alpha;q_{t}[GeV];Pt [GeV]", 80,-1.0,+1.0, 90,0,0.3, 60,0,12); tList->Add( tH3D );
    fQAList->Add(tList);
    tList = new TList(); tList->SetOwner(); tList->SetName("QACutsAfter_OP");
    tH3D = new TH3D("AftVOL", "VOLUME;Phi;Eta;Pt [GeV]",       63, 0.0,+6.3, 40,-1.0,+1.0, 60,0,12); tList->Add( tH3D );
    tH3D = new TH3D("AftRAP", "RAPIDITY;y;Pt [GeV]",           40,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftDLZ", "DLZ;[cm];Pt [GeV]",             50,-5.0,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftDLXY","DLXY;[cm];Pt [GeV]",           100,0.00,100., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftCTau","CTau;[cm];Pt [GeV]",           250,0.00,250., 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftDCA", "DCA;[cm];Pt [GeV]",             50,0.00,+5.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftCTP", "CTP;;Pt [GeV]",                 50,0.99,1.00, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftD0",  "D0;[cm];Pt [GeV]",              50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftD0D0","D0D0;[cm^{2}];Pt [GeV]",        50,-1.0,+1.0, 60,0,12, nMass, dMinMass, dMaxMass); tList->Add( tH3D );
    tH3D = new TH3D("AftAP",  "AP;#alpha;q_{t}[GeV];Pt [GeV]", 80,-1.0,+1.0, 90,0,0.3, 60,0,12); tList->Add( tH3D );
    fQAList->Add(tList);
  }

  if(fMCmatch>0) { // only for MCanalysis
    printf("MC mode enabled\n");
    tList = new TList(); tList->SetOwner(); tList->SetName("QAMC");
    tH3D = new TH3D("MCPDG","MC PDGm;p_{T} (GeV);Mass [GeV]",24,0,12,nMass,dMinMass,dMaxMass,30,0,30); tList->Add( tH3D );
    tH3D->GetZaxis()->SetBinLabel( 1,"NONE");
    tH3D->GetZaxis()->SetBinLabel( 2,"NONE chked");
    tH3D->GetZaxis()->SetBinLabel( 3,"NONPDG");
    tH3D->GetZaxis()->SetBinLabel( 4,"NONPDG chked"); // 3
    tH3D->GetZaxis()->SetBinLabel( 5,"Truths");
    tH3D->GetZaxis()->SetBinLabel( 6,"Truths chked"); // 1
    tH3D->GetZaxis()->SetBinLabel( 7,"K_L0 chked");
    tH3D->GetZaxis()->SetBinLabel( 8,"K_S0 chked");
    tH3D->GetZaxis()->SetBinLabel( 9,"Lambda0 chked");
    tH3D->GetZaxis()->SetBinLabel(10,"Lambda0_bar chked");
    tH3D->GetZaxis()->SetBinLabel(11,"phi chked");
    tH3D->GetZaxis()->SetBinLabel(12,"rho0 chked");
    tH3D->GetZaxis()->SetBinLabel(13,"omega chked");
    tH3D->GetZaxis()->SetBinLabel(14,"f_0 chked");
    tH3D->GetZaxis()->SetBinLabel(15,"e- chked");
    tH3D->GetZaxis()->SetBinLabel(16,"e+ chked");
    tH3D->GetZaxis()->SetBinLabel(17,"pi+ chked");
    tH3D->GetZaxis()->SetBinLabel(18,"pi- chked");
    tH3D->GetZaxis()->SetBinLabel(19,"mu+ chked");
    tH3D->GetZaxis()->SetBinLabel(20,"mu- chked");
    tH3D->GetZaxis()->SetBinLabel(21,"K+ chked");
    tH3D->GetZaxis()->SetBinLabel(22,"K- chked");
    tH3D->GetZaxis()->SetBinLabel(23,"proton chked");
    tH3D->GetZaxis()->SetBinLabel(24,"antiproton chked");
  
    tH3D = new TH3D("MCMOTHER","MC MOTHER;truth p_{T} (GeV); mother p_{T} (GeV)",60,0,12,60,0,12,50,0,50); tList->Add( tH3D );
    tH3D->GetZaxis()->SetBinLabel( 1,"NONPDG chked");
    tH3D->GetZaxis()->SetBinLabel( 2,"PRIMARY chked");
    tH3D->GetZaxis()->SetBinLabel( 3,"Xi0 chked");
    tH3D->GetZaxis()->SetBinLabel( 4,"Xi0_bar chked");
    tH3D->GetZaxis()->SetBinLabel( 5,"Xi- chked");
    tH3D->GetZaxis()->SetBinLabel( 6,"Xi-_bar chked");
    tH3D->GetZaxis()->SetBinLabel( 7,"Omega- chked");
    tH3D->GetZaxis()->SetBinLabel( 8,"Omega+ chked");
    tH3D->GetZaxis()->SetBinLabel( 9,"K+ chked");
    tH3D->GetZaxis()->SetBinLabel(10,"K- chked");
    tH3D->GetZaxis()->SetBinLabel(11,"Lambda0 chked");
    tH3D->GetZaxis()->SetBinLabel(12,"Lambda0_bar chked");
    tH3D->GetZaxis()->SetBinLabel(13,"K0 chked");
    tH3D->GetZaxis()->SetBinLabel(14,"K0_bar chked");
    tH3D->GetZaxis()->SetBinLabel(15,"K_S0 chked");
    tH3D->GetZaxis()->SetBinLabel(16,"K_L0 chked");
    tH3D->GetZaxis()->SetBinLabel(17,"phi chked");
    tH3D->GetZaxis()->SetBinLabel(18,"D+ chked");
    tH3D->GetZaxis()->SetBinLabel(19,"D- chked");
    tH3D->GetZaxis()->SetBinLabel(20,"D0 chked");
    tH3D->GetZaxis()->SetBinLabel(21,"D0_bar chked");
    tH3D->GetZaxis()->SetBinLabel(22,"D_s+ chked");
    tH3D->GetZaxis()->SetBinLabel(23,"D_s- chked");
    tH3D->GetZaxis()->SetBinLabel(24,"Lambda_c+ chked");
    tH3D->GetZaxis()->SetBinLabel(25,"Lambda_c- chked");
    tH3D->GetZaxis()->SetBinLabel(26,"Sigma*- chked");
    tH3D->GetZaxis()->SetBinLabel(27,"Sigma*-_bar chked");
    tH3D->GetZaxis()->SetBinLabel(28,"Sigma*+ chked");
    tH3D->GetZaxis()->SetBinLabel(29,"Sigma*+_bar chked");
    tH3D->GetZaxis()->SetBinLabel(30,"Sigma*0 chked");
    tH3D->GetZaxis()->SetBinLabel(31,"Sigma*0_bar chked");
    tH3D->GetZaxis()->SetBinLabel(32,"Sigma- chked");
    tH3D->GetZaxis()->SetBinLabel(33,"Sigma+_bar chked");
    tH3D->GetZaxis()->SetBinLabel(34,"Sigma0 chked");
    tH3D->GetZaxis()->SetBinLabel(35,"Sigma0_bar chked");
    tH3D->GetZaxis()->SetBinLabel(36,"NONE chked");
    tH1D = new TH1D("MCDAUGHTERS","MC DAUGHTERS",10,0,10); tList->Add( tH1D );
    tH2D = new TH2D("MCRADS","MC RADS",2,0,2,200,0,10); tList->Add( tH2D );
    fQAList->Add(tList);
  } // only for MCanalysis

}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyUserExec(Option_t *)
{
  // user exec
  AliAODEvent *tAOD=dynamic_cast<AliAODEvent*>(InputEvent());
  Bool_t acceptEvent=kFALSE;
  fCandidates->SetLast(-1);
  if(tAOD) {
    ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("Events"))->Fill(0);
    Double_t tVtxZ = tAOD->GetPrimaryVertex()->GetZ();
    Double_t tSPDVtxZ = tAOD->GetPrimaryVertexSPD()->GetZ();
    if( fDoQA )
      ((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
    Bool_t tESelection = TMath::Abs(tVtxZ-tSPDVtxZ) < 0.5;
    Bool_t tDSelection = kTRUE;
    if(fUseEventSelection) {
      tDSelection = fCutsEvent->IsSelected(tAOD);
    } else {
      if(TMath::Abs(tVtxZ)>10.0) tDSelection = kFALSE; // Cut on VtxZ mandatory!
    }
    if(tDSelection&&tESelection) {
      acceptEvent=kTRUE;
      ReadFromAODv0(tAOD);
      if( fDoQA )
	((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("VTXZSEL"))->Fill( tVtxZ, tSPDVtxZ );
    }
  }
  if(!acceptEvent) return;
  // QA filling
  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("Events"))->Fill(1);
  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("RFPTPC"))->Fill( fFlowEventTPC->GetNumberOfRPs() );
  Double_t mult=0;
  for(Int_t i=0;i!=fFlowEventVZE->GetNumberOfRPs();++i) {
    AliFlowTrackSimple *pTrack = fFlowEventVZE->GetTrack(i);
    mult += pTrack->Weight();
  }
  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("RFPVZE"))->Fill( mult );
  ((TH1D*)((TList*)fQAList->FindObject("Event"))->FindObject("POI"))->Fill( fCandidates->GetEntriesFast() );
  AddCandidates();

  PostData(1,fFlowEventTPC);
  PostData(2,fFlowEventVZE);
  PostData(3,fQAList);
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddCandidates()
{
  // adds candidates to flow events (untaging if necessary)
  if(fDebug) printf("I received %d candidates\n",fCandidates->GetEntriesFast());
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand = dynamic_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if(!cand) continue;
    if(fDebug) printf(" >Checking at candidate %d with %d daughters: mass %f\n",
		      iCand,cand->GetNDaughters(),cand->Mass());
    // untagging ===>
    for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
      if(fDebug) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
      for(int iRPs=0; iRPs!=fFlowEventTPC->NumberOfTracks(); ++iRPs ) {
        AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEventTPC->GetTrack( iRPs ));
        if (!iRP) continue;
        if( !iRP->InRPSelection() ) continue;
        if( cand->GetIDDaughter(iDau) == iRP->GetID() ) {
          if(fDebug) printf(" was in RP set");
          iRP->SetForRPSelection(kFALSE);
          fFlowEventTPC->SetNumberOfRPs( fFlowEventTPC->GetNumberOfRPs() -1 );
        }
      }
      if(fDebug) printf("\n");
    }
    // <=== untagging
    cand->SetForPOISelection(kTRUE);
    fFlowEventTPC->InsertTrack( ((AliFlowTrack*) cand) );
    fFlowEventVZE->InsertTrack( ((AliFlowTrack*) cand) );
  }
  if(fDebug) printf("TPCevent %d | VZEevent %d\n",
		    fFlowEventTPC->NumberOfTracks(),
		    fFlowEventVZE->NumberOfTracks() );
}

//=======================================================================
void AliAnalysisTaskFlowStrange::ChargedParticleAnalysis(AliAODEvent *tAOD)
{
  if(fMCmatch>0) fBayesianPID->SetMC(kTRUE);
  fBayesianPID->SetDetResponse(tAOD,999);
  for(int iRPs=0; iRPs!=fFlowEventTPC->NumberOfTracks(); ++iRPs ) {
    AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEventTPC->GetTrack( iRPs ));
    if(!iRP) continue;
    Int_t ntracks = tAOD->GetNTracks();
    AliAODTrack *iT = NULL;
    for(Int_t it=0; it!=ntracks; ++it) {
      iT = (AliAODTrack*) tAOD->GetTrack(it);
      if(iT->GetID() == iRP->GetID()) break;
      iT = NULL;
    }
    if(!iT) continue;
    fBayesianPID->ComputeProb(iT,tAOD);
    Float_t *prob = fBayesianPID->GetProb();
    Float_t tofMismProb = fBayesianPID->GetTOFMismProb();
    if( fBayesianPID->GetCurrentMask(1) && tofMismProb<0.5 ) {
      iRP->SetMass( prob[3] );
      iRP->SetForPOISelection(kTRUE);
      iRP->SetForRPSelection(kFALSE);
      fFlowEventVZE->InsertTrack( ((AliFlowTrack*) iRP) );
      iRP->SetForRPSelection(kTRUE);
      // === MATCHED TO MC ===>>
      if(fMCmatch>0) {
      	TClonesArray* mcArray = dynamic_cast<TClonesArray*>(tAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      	if(mcArray) {
      	  TString sPDG="NONE";
      	  TString sMOTHER = "NONE";
      	  Double_t ptTruth=-1, ptMom=-1;
      	  AliAODMCParticle *iTMC = dynamic_cast<AliAODMCParticle*>(mcArray->At( TMath::Abs(iT->GetLabel()) ));
          if(iTMC) {
            sPDG="NONPDG";
            ptTruth = iTMC->Pt();
      	    Int_t iPDG = iTMC->GetPdgCode();
            TDatabasePDG *pdgDatabase = TDatabasePDG::Instance();
            if(pdgDatabase->GetParticle(iPDG))
              sPDG = (pdgDatabase->GetParticle(iPDG))->GetName();
            if(iTMC->GetMother()>=0) {
              AliAODMCParticle *iTMom = dynamic_cast<AliAODMCParticle*>(mcArray->At(iTMC->GetMother()));
              if(iTMom) {
                ptMom = iTMom->Pt();
                sMOTHER="NONPDG";
                Int_t iMomPDG = iTMom->GetPdgCode();
                pdgDatabase = TDatabasePDG::Instance();
                if(pdgDatabase->GetParticle(iMomPDG))
                sMOTHER = (pdgDatabase->GetParticle(iMomPDG))->GetName();
              }
            } else {
              sMOTHER="PRIMARY";
            }
          }
          Double_t dMASS=prob[3];
          if(dMASS>1) dMASS=1;
          if(dMASS<0) dMASS=0;
          if(iT->GetLabel()>=0) {
            sPDG = Form("%s chked",sPDG.Data());
            sMOTHER = Form("%s chked",sMOTHER.Data());
          }
          ((TH3D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCPDG"))->Fill(iT->Pt(),dMASS,sPDG.Data(),1);
          if(dMASS>=0.90)
            ((TH3D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCMOTHER"))->Fill(ptTruth,ptMom,sMOTHER.Data(),1);
        }
      }
      // <<=== MATCHED TO MC ===
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromAODv0(AliAODEvent *tAOD)
{
  fCutsRFPTPC->SetEvent(tAOD,MCEvent());
  fCutsRFPVZE->SetEvent(tAOD,MCEvent());
  fCutsPOI->SetEvent(tAOD,MCEvent());
  fFlowEventTPC->Fill(fCutsRFPTPC,fCutsPOI);
  fFlowEventVZE->Fill(fCutsRFPVZE,fCutsPOI);
  fPsi2 = (fFlowEventVZE->GetQ()).Phi()/2;

  for(Int_t i=0; i!=fFlowEventTPC->NumberOfTracks(); i++) {
    AliFlowTrackSimple *pTrack = (AliFlowTrackSimple*) fFlowEventTPC->GetTrack(i);
    if(!pTrack) continue;
    if(!pTrack->InRPSelection()) continue;
    ((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("TPCPhiEta"))->Fill( pTrack->Phi(), pTrack->Eta(), pTrack->Weight() );
  }
  for(Int_t i=0; i!=fFlowEventVZE->NumberOfTracks(); i++) {
    AliFlowTrackSimple *pTrack = (AliFlowTrackSimple*) fFlowEventVZE->GetTrack(i);
    if(!pTrack) continue;
    if(!pTrack->InRPSelection()) continue;
    ((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("VZEPhiEta"))->Fill( pTrack->Phi(), pTrack->Eta(), pTrack->Weight() );
  }

  if(fDoQA) {
    AliFlowVector Qs[2];
    fFlowEventVZE->TagSubeventsInEta(-5,1,1,+5);
    fFlowEventVZE->Get2Qsub(Qs,2);
    Double_t dEPV0C = Qs[1].Phi()/2;
    Double_t dEPV0A = Qs[0].Phi()/2;
    Double_t dEPTPC = (fFlowEventTPC->GetQ()).Phi()/2;
    ((TH3D*)((TList*)fQAList->FindObject("Event"))->FindObject("EVPLANE"))->Fill( dEPTPC, dEPV0A, dEPV0C );
    Double_t dVZEQx = (fFlowEventVZE->GetQ(2)).Px() / (fFlowEventVZE->GetQ(2)).GetMult();
    Double_t dVZEQy = (fFlowEventVZE->GetQ(2)).Py() / (fFlowEventVZE->GetQ(2)).GetMult();
    ((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("VZEQxQy"))->Fill( dVZEQx, dVZEQy );
    ((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("TPCQxQy"))->Fill( (fFlowEventTPC->GetQ(2)).Px() / (fFlowEventTPC->GetQ(2)).GetMult(),
										  (fFlowEventTPC->GetQ(2)).Py() / (fFlowEventTPC->GetQ(2)).GetMult() );
    // TPC Q
    const int ngr=24;
    double dArrayPt[25] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
			   1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
			   2.5, 3.0, 3.5, 4.0, 5.0 };
    double dTPCPt[ngr];
    double dTPCQx[ngr];
    double dTPCQy[ngr];
    double dTPCEta[ngr];
    double dTPCM[ngr];
    for(int i=0; i!=ngr; ++i) {
      dTPCPt[i] = 0;
      dTPCQx[i] = 0;
      dTPCQy[i] = 0;
      dTPCEta[i] = 0;
      dTPCM[i] = 0;
    }
    for(Int_t i=0; i!=fFlowEventTPC->NumberOfTracks(); i++) {
      AliFlowTrackSimple *pTrack = (AliFlowTrackSimple*) fFlowEventTPC->GetTrack(i);
      if(!pTrack) continue;
      if(!pTrack->InRPSelection()) continue;
      Double_t dPt  = pTrack->Pt();
      int npt=-1;
      for(int pt=0; pt!=ngr; ++pt)
	if( (dPt > dArrayPt[pt])&&(dPt < dArrayPt[pt+1]) ) {
	  npt = pt;
	  break;
	}
      if(npt<0) continue;
      Double_t dPhi = pTrack->Phi();
      Double_t dWeight = pTrack->Weight();
      dTPCPt[npt] += dWeight*dPt;
      dTPCQx[npt] += dWeight*TMath::Cos(2*dPhi);
      dTPCQy[npt] += dWeight*TMath::Sin(2*dPhi);
      dTPCEta[npt] += dWeight*pTrack->Eta();
      dTPCM[npt] += dWeight;
    }
    for(int i=0; i!=ngr; ++i)
      if( dTPCM[i]>0 ) {
	((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("TPCPtQx"))->Fill( dTPCPt[i]/dTPCM[i], dTPCQx[i]/dTPCM[i] );
	((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("TPCPtQy"))->Fill( dTPCPt[i]/dTPCM[i], dTPCQy[i]/dTPCM[i] );
	((TH2D*)((TList*)fQAList->FindObject("Event"))->FindObject("TPCPtEta"))->Fill( dTPCPt[i]/dTPCM[i], dTPCEta[i]/dTPCM[i] );
      }
    // End of TPC Q
    const AliAODVertex* trkVtx = tAOD->GetPrimaryVertex();
    Double_t cov[6]={0};
    trkVtx->GetCovarianceMatrix(cov);
    ((TH3D*)((TList*)fQAList->FindObject("Event"))->FindObject("PRIMVERTEX"))->Fill( TMath::Sqrt( cov[0] ),
										     TMath::Sqrt( cov[2] ),
										     TMath::Sqrt( cov[5] ) );
  }

  if(fSpecie>80) {
    ChargedParticleAnalysis(tAOD);
    return;
  }
  Int_t nV0s = tAOD->GetNumberOfV0s();
  AliAODv0 *myV0;
  Double_t dMASS=0.0;
  for (Int_t i=0; i!=nV0s; ++i) {
    myV0 = (AliAODv0*) tAOD->GetV0(i);
    if(!myV0) continue;
    if(myV0->Pt()<0.1) continue; // skipping low momentum
    Int_t pass = PassesAODCuts(myV0,tAOD);
    if(pass==0) continue;
    if(fSpecie==0) {
      dMASS = myV0->MassK0Short();
    } else {
      dMASS = myV0->MassLambda();
      if(pass==2) dMASS = myV0->MassAntiLambda();
    }
    MakeTrack(dMASS, myV0->Pt(), myV0->Phi(), myV0->Eta(),
              ((AliAODTrack*) myV0->GetDaughter(0))->GetID(),
              ((AliAODTrack*) myV0->GetDaughter(1))->GetID());
  }
  return;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD)
{
  if (myV0->GetOnFlyStatus() ) return 0;
  //the following is needed in order to evualuate track-quality
  AliAODTrack *iT, *jT;
  AliAODVertex *vV0s = myV0->GetSecondaryVtx();
  Double_t pos[3],cov[6];
  vV0s->GetXYZ(pos);
  vV0s->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  // TESTING CHARGE
  int iPos, iNeg;
  iT=(AliAODTrack*) myV0->GetDaughter(0);
  if(iT->Charge()>0) {
    iPos = 0; iNeg = 1;
  } else {
    iPos = 1; iNeg = 0;
  }
  // END OF TEST

  iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive
  AliESDtrack ieT( iT );
  ieT.SetTPCClusterMap( iT->GetTPCClusterMap() );
  ieT.SetTPCSharedMap( iT->GetTPCSharedMap() );
  ieT.SetTPCPointsF( iT->GetTPCNclsF() );
  ieT.RelateToVertex(&vESD, tAOD->GetMagneticField(), 100);
  if (!fCutsDau->IsSelected( &ieT ) ) return 0;

  jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative
  AliESDtrack jeT( jT );
  jeT.SetTPCClusterMap( jT->GetTPCClusterMap() );
  jeT.SetTPCSharedMap( jT->GetTPCSharedMap() );
  jeT.SetTPCPointsF( jT->GetTPCNclsF() );
  jeT.RelateToVertex(&vESD, tAOD->GetMagneticField(), 100);
  if (!fCutsDau->IsSelected( &jeT ) ) return 0;

  Double_t pvertex[3];
  pvertex[0]=tAOD->GetPrimaryVertex()->GetX();
  pvertex[1]=tAOD->GetPrimaryVertex()->GetY();
  pvertex[2]=tAOD->GetPrimaryVertex()->GetZ();
  Double_t dRAP;
  if(fSpecie==0)
    dRAP=myV0->RapK0Short();
  else
    dRAP=myV0->RapLambda();
  Double_t dDL=myV0->DecayLengthV0( pvertex );
  Double_t dDLZ=myV0->DecayVertexV0Z()-pvertex[2];
  Double_t dDLXY=myV0->RadiusV0();
  Double_t dDCA=myV0->DcaV0Daughters();
  Double_t dCTP=myV0->CosPointingAngle( pvertex );
  Double_t dD0P=ieT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
  Double_t dD0M=jeT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
  Double_t dD0D0=dD0P*dD0M;
  Double_t dQT=myV0->PtArmV0();
  Double_t dALPHA=myV0->AlphaV0(); // AlphaV0 -> AODRecoDecat::Alpha -> return 1.-2./(1.+QlProng(0)/QlProng(1));
  if(myV0->ChargeProng(iPos)<0) dALPHA = -dALPHA; // protects for a change in convention
  Double_t dPT=myV0->Pt();
  Double_t dETA=myV0->Eta();
  Int_t passes = 1;
  if(fSpecie==1&&dALPHA<0) passes = 2; // antilambda
  Double_t dMASS = myV0->MassK0Short();
  if(fSpecie==1) {
    if(passes==2) dMASS = myV0->MassAntiLambda();
    else dMASS = myV0->MassLambda();
  }
  Double_t dCT = dDL*dMASS/myV0->P();
  Double_t ctaucut = 2.68;
  if(fSpecie==1) ctaucut = 7.89;
  ctaucut *= fV0Cuts[9];
  if(dDL<fV0Cuts[0]) passes = 0;
  if(dDCA >fV0Cuts[1]) passes = 0;
  if(dCTP <fV0Cuts[2]) passes = 0;
  if(TMath::Abs(dD0P) <fV0Cuts[3]) passes = 0;
  if(TMath::Abs(dD0M) <fV0Cuts[3]) passes = 0;
  if(dD0D0>fV0Cuts[4]) passes = 0;
  if(dETA <fV0Cuts[6]) passes = 0;
  if(dETA >fV0Cuts[7]) passes = 0;
  if(fSpecie==0) if(dQT<+fV0Cuts[5]*dALPHA) passes = 0;
  if(fSpecie==0) if(dQT<-fV0Cuts[5]*dALPHA) passes = 0;
  if(dCT>ctaucut) passes = 0;
  if(dDLXY<fV0Cuts[10]) passes = 0;

  Double_t dPHI = myV0->Phi();
  Double_t dDPHI = dPHI - fPsi2;
  if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
  if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
  TString sIPOP = "IP";
  if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) )
    sIPOP = "OP";
  if(fDoExtraQA) {
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefRAP")) ->Fill(dRAP,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefDLZ")) ->Fill(dDLZ,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefDLXY"))->Fill(dDLXY, dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefCTau"))->Fill(dCT,   dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefDCA")) ->Fill(dDCA,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefCTP")) ->Fill(dCTP,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefD0"))  ->Fill(dD0M,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefD0"))  ->Fill(dD0P,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefD0D0"))->Fill(dD0D0, dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefAP"))  ->Fill(dALPHA,dQT,dPT);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsBefore_%s",sIPOP.Data())))->FindObject("BefVOL")) ->Fill(dPHI, dETA,dPT);
  }

  if(passes&&fV0Cuts[8]) {
    switch(fSpecie) {
    case 0: // K0 PID
      if( (jT->GetTPCmomentum()<15) &&
        (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kPion))>fV0Cuts[8]) )
        passes = 0;
      if( (iT->GetTPCmomentum()<15) &&
        (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kPion))>fV0Cuts[8]) )
        passes = 0;
      break;
    case 1: // Lambda PID  i==pos j ==neg
      if(passes==1) {
        if( (iT->GetTPCmomentum()<15) &&
          (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kProton))>fV0Cuts[8]) )
          passes = 0;
        if( (jT->GetTPCmomentum()<15) &&
          (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kPion))>fV0Cuts[8]) )
          passes = 0;
      }
      if(passes==2) {
        if( (iT->GetTPCmomentum()<15) &&
          (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(iT,AliPID::kPion))>fV0Cuts[8]) )
          passes = 0;
        if( (jT->GetTPCmomentum()<15) &&
          (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(jT,AliPID::kProton))>fV0Cuts[8]) )
          passes = 0;
      }
      break;
    }
  }
  // === MATCHED TO MC ===>>
  // AliAODRecoDecay::MatchToMC()
  if((passes>0)&&(fMCmatch>0)) {
    TClonesArray* mcArray = dynamic_cast<TClonesArray*>(tAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    Int_t mcpasses=0;
    if(mcArray) {
      TString sPDG="NONE";
      TString sMOTHER = "NONE";
      Double_t ptTruth=-1, ptMom=-1;
      Int_t abc=4; // 0:reject   1:truth   2:another decay   3:nonv0 t2   4:nonv0 t1
      AliAODMCParticle *iTMC = dynamic_cast<AliAODMCParticle*>(mcArray->At( TMath::Abs(iT->GetLabel()) ));
      AliAODMCParticle *jTMC = dynamic_cast<AliAODMCParticle*>(mcArray->At( TMath::Abs(jT->GetLabel()) ));
      if((iTMC)&&(jTMC)) {
        if((iTMC->GetMother()>=0)&&(jTMC->GetMother()>=0)) {
          AliAODMCParticle *iTMom = dynamic_cast<AliAODMCParticle*>(mcArray->At(iTMC->GetMother()));
          AliAODMCParticle *jTMom = dynamic_cast<AliAODMCParticle*>(mcArray->At(jTMC->GetMother()));
          if((iTMom)&&(jTMom)) {
            sPDG="NONPDG";
            abc=3;
            if(iTMC->GetMother()==jTMC->GetMother()) {
              Int_t iMomPDG = iTMom->GetPdgCode();
              abc=2;
              TDatabasePDG *pdgDatabase = TDatabasePDG::Instance();
              if(pdgDatabase->GetParticle(iMomPDG))
                sPDG = (pdgDatabase->GetParticle(iMomPDG))->GetName();
              Int_t pdgcode;
              if(fSpecie) pdgcode = 3122; else pdgcode = 310;
              if(TMath::Abs(iMomPDG)==pdgcode) {
                Double_t pxSumDgs = iTMC->Px()+jTMC->Px();
                Double_t pySumDgs = iTMC->Py()+jTMC->Py();
                Double_t pzSumDgs = iTMC->Pz()+jTMC->Pz();
                Double_t pxMother = iTMom->Px();
                Double_t pyMother = iTMom->Py();
                Double_t pzMother = iTMom->Pz();
                sPDG = "Truths";
                Double_t rads = TMath::Sqrt( iTMC->Xv()*iTMC->Xv() + iTMC->Yv()*iTMC->Yv() );
                ((TH1D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCDAUGHTERS"))->Fill(iTMom->GetNDaughters());
                if((TMath::Abs(pxMother-pxSumDgs)/(TMath::Abs(pxMother)+1.e-13)) < 0.00001 &&
                  (TMath::Abs(pyMother-pySumDgs)/(TMath::Abs(pyMother)+1.e-13)) < 0.00001 &&
                  (TMath::Abs(pzMother-pzSumDgs)/(TMath::Abs(pzMother)+1.e-13)) < 0.00001) {
                  ((TH2D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCRADS"))->Fill(1.0,rads);
		  abc=1;
		  if(iTMom->GetMother()>=0) {
		    AliAODMCParticle *iGrandMa = dynamic_cast<AliAODMCParticle*>(mcArray->At(iTMom->GetMother()));
		    if(iGrandMa) {
		      Int_t iGrandMaPDG = iGrandMa->GetPdgCode();
		      ptMom = iGrandMa->Pt();
		      ptTruth = iTMom->Pt();
		      if(pdgDatabase->GetParticle(iGrandMaPDG))
			sMOTHER = (pdgDatabase->GetParticle(iGrandMaPDG))->GetName();
		    }
                  } else {
                    sMOTHER = "PRIMARY";
                  }
                } else {
                  ((TH2D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCRADS"))->Fill(0.0,rads);
                }
              }
            }
          }
        }
      }
      Double_t dMCMASS;
      if(fSpecie==0) {
        dMCMASS = myV0->MassK0Short();
        if(dMCMASS>0.588) dMCMASS=0.587999;
        if(dMCMASS<0.412) dMCMASS=0.412001;
      } else {
        dMCMASS = myV0->MassLambda();
        if(passes==2) dMCMASS = myV0->MassAntiLambda();
        if(dMCMASS>1.167) dMCMASS=1.1669;
        if(dMCMASS<1.075) dMCMASS=1.0751;
      }
      if((iT->GetLabel()>=0)&&(jT->GetLabel()>=0)) {
        sPDG = Form("%s chked",sPDG.Data());
        sMOTHER = Form("%s chked",sMOTHER.Data());
        if(fMCmatch==abc) mcpasses=1;
      }
      ((TH3D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCPDG"))->Fill(dPT,dMCMASS,sPDG.Data(),1);
      ((TH3D*)((TList*)fQAList->FindObject("QAMC"))->FindObject("MCMOTHER"))->Fill(ptTruth,ptMom,sMOTHER.Data(),1);
    }//mcArray
    passes=mcpasses;
  }
  // <<=== MATCHED TO MC ===
  if(passes&&fDoExtraQA) {
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftRAP")) ->Fill(dRAP,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftDLZ")) ->Fill(dDLZ,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftDLXY"))->Fill(dDLXY, dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftCTau"))->Fill(dCT,   dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftDCA")) ->Fill(dDCA,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftCTP")) ->Fill(dCTP,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftD0"))  ->Fill(dD0M,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftD0"))  ->Fill(dD0P,  dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftD0D0"))->Fill(dD0D0, dPT, dMASS);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftAP"))  ->Fill(dALPHA,dQT,dPT);
    ((TH3D*)((TList*)fQAList->FindObject(Form("QACutsAfter_%s",sIPOP.Data())))->FindObject("AftVOL")) ->Fill(dPHI, dETA,dPT);
  }
  if(passes&&fDoQA) {
    ((TH2D*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0MASS"))->Fill( dPT, dMASS );
    ((TH2D*)((TList*)fQAList->FindObject("Candidates"))->FindObject("V0PhiEta"))->Fill( dPHI, dETA );
  }
  return passes;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Terminate(Option_t *)
{
  //terminate
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeTrack( Double_t mass, Double_t pt, Double_t phi,
					    Double_t eta, Int_t iid, Int_t jid )
{
  // create track for flow tasks
  if(fCandidates->GetLast()+1>=fCandidates->GetSize()) {
    fCandidates->Expand( 2*fCandidates->GetSize() );
  }
  Bool_t overwrite = kTRUE;
  AliFlowCandidateTrack *oTrack = (static_cast<AliFlowCandidateTrack*> (fCandidates->At( fCandidates->GetLast()+1 )));
  if( !oTrack ) { // creates new
    oTrack = new AliFlowCandidateTrack();
    overwrite = kFALSE;
  } else { // overwrites
    oTrack->ClearMe();
  }
  oTrack->SetMass(mass);
  oTrack->SetPt(pt);
  oTrack->SetPhi(phi);
  oTrack->SetEta(eta);
  oTrack->AddDaughter(iid);
  oTrack->AddDaughter(jid);
  oTrack->SetForPOISelection(kTRUE);
  oTrack->SetForRPSelection(kFALSE);
  if(overwrite) {
    fCandidates->SetLast( fCandidates->GetLast()+1 );
  } else {
    fCandidates->AddLast(oTrack);
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass)
{
  // setter for mass bins
  fMassBins = massBins;
  fMinMass = minMass;
  fMaxMass = maxMass;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SetCuts(Double_t cuts[11]) {
  // defines cuts to be used
  // fV0Cuts[11] dl3d dca ctp d0 d0d0 qt minEta maxEta PID ct dlxy
  for(int i=0; i!=11; ++i) fV0Cuts[i] = cuts[i];
}

