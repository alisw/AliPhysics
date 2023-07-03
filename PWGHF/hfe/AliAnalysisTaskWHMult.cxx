/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 * Edited by Taiga Kawaguchi                                              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskWHMult
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskWHMult.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "TFile.h"
#include "TRandom.h"
#include "TGrid.h"

class AliAnalysisTaskWHMult;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskWHMult) // classimp: necessary for root

AliAnalysisTaskWHMult::AliAnalysisTaskWHMult() : AliAnalysisTaskSE(), 
  fAOD(0),
  //fMCEvent(0),
  fVevent(0),
  fPIDResponse(0),
  fMCTrackpart(0),
  fMCarray(0),
  fMCparticle(0),
  fOutputList(0),
  fHistPt(0),
  pVertex(0),
  cutVer(0),
  EtavsPhi(0),
  TPCSig(0),
  Cent(0),
  fzvtx_Ntrkl(0),
  fzvtx_Ntrkl_cal(0),
  fNtrklNch(0),
  fPDG(0),
  fMPDG(0),
  fFPDG(0),
  fHistClustE(0),
  fHistClustEMatch(0),
  fHistClustLongAxis(0),
  fHistClustLongAxisE(0),
  fHistNsigmaP(0),
  fHistMCNsigmaP(0),
  fPtEoverPE(0),
  fPtEoverPMCE(0),
  fPtEoverPEGeo(0),
  fHistMCClsLAE(0),
  fPtEoverPH(0),
  fPtEoverPHGeo(0),
  fREisolation(0),
  fREisoW(0),
  fREisoHF(0),
  fREisoWhpt(0),
  fREisoHFhpt(0),
  fREisoWvhpt(0),
  fREisoHFvhpt(0),
  EleTraDelPhi(0),
  fHistPt_We(0),
  fHistPt_HFe(0),
  fPt_maxtrack_W(0),
  fNtrkl_PtOfMaxTrk_W(0),
  fHistPt_We_Ntrkl(0),
  fW_true(0),
  fW_false(0),
  EleTraDelPhi_fullrange(0),
  fEMCEG1(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskWHMult::AliAnalysisTaskWHMult(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),
  //fMCEvent(0),
  fVevent(0),
  fPIDResponse(0),
  fMCTrackpart(0),
  fMCarray(0),
  fMCparticle(0),
  fOutputList(0),
  fHistPt(0),
  pVertex(0),
  cutVer(0),
  EtavsPhi(0),
  TPCSig(0),
  Cent(0),
  fzvtx_Ntrkl(0),
  fzvtx_Ntrkl_cal(0),
  fNtrklNch(0),
  fPDG(0),
  fMPDG(0),
  fFPDG(0),
  fHistClustE(0),
  fHistClustEMatch(0),
  fHistClustLongAxis(0),
  fHistClustLongAxisE(0),
  fHistNsigmaP(0),
  fHistMCNsigmaP(0),
  fPtEoverPE(0),
  fPtEoverPMCE(0),
  fPtEoverPEGeo(0),
  fHistMCClsLAE(0),
  fPtEoverPH(0),
  fPtEoverPHGeo(0),
  fREisolation(0),
  fREisoW(0),
  fREisoHF(0),
  fREisoWhpt(0),
  fREisoHFhpt(0),
  fREisoWvhpt(0),
  fREisoHFvhpt(0),
  EleTraDelPhi(0),
  fHistPt_We(0),
  fHistPt_HFe(0),
  fPt_maxtrack_W(0),
  fNtrkl_PtOfMaxTrk_W(0),
  fHistPt_We_Ntrkl(0),
  fW_true(0),
  fW_false(0),
  EleTraDelPhi_fullrange(0),
  fEMCEG1(kFALSE)
{
  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                      // this chain is created by the analysis manager, so no need to worry about it, 
                                      // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                      // you can add more output objects by calling DefineOutput(2, classname::Class())
                                      // if you add more output objects, make sure to call PostData for all of them, and to
                                      // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskWHMult::~AliAnalysisTaskWHMult()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskWHMult::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fHistPt = new TH1F("fHistPt", "fHistPt", 200, 0, 200);

  pVertex = new TH1F("pVertex", "pVertex", 100, -20, 20);       // my histgram
  pVertex->GetXaxis()->SetTitle("vertex (collision point on z-axis)");
  pVertex->GetYaxis()->SetTitle("events");

  cutVer = new TH1F("cutVer", "cutVer", 100, -20, 20);
  cutVer->GetXaxis()->SetTitle("vertex (collision point on z-axis)");
  cutVer->GetYaxis()->SetTitle("events");

  EtavsPhi = new TH2F("EtavsPhi","#eta vs #phi", 1000, -1.5, 1.5, 1000, -1, 7);
  EtavsPhi->GetXaxis()->SetTitle("#eta");
  EtavsPhi->GetYaxis()->SetTitle("#phi (rad)");

  TPCSig = new TH2F("TPCSig","TPC signal", 5000, 0.1, 30, 10000, 0, 10000);
  TPCSig->GetXaxis()->SetTitle("P (GeV/c)");
  TPCSig->GetYaxis()->SetTitle("dE/dx");
  TPCSig->SetMarkerStyle(7);

  Cent = new TH1F("Cent","centrality", 120, -10, 110);
  Cent->GetXaxis()->SetTitle("centrality (%)");
  Cent->GetYaxis()->SetTitle("counts");

  fzvtx_Ntrkl = new TH2F("fzvtx_Ntrkl","vertexZ vs Number of Tracklets; vertex Z (cm); Number of Tracklets",100,-10,10,300,0,300);
  fzvtx_Ntrkl_cal=new TH2F("fzvtx_Ntrkl_cal","vertexZ vs Number of Tracklets after calibration;vertex Z (cm);Number of Tracklets",100,-10,10,300,0,300);
  fNtrklNch = new TH2F("fNtrklNch","N_{tracklets} vs N_{ch}; N_{tracklets}; N_{ch}",200,0,200,200,0,200);

  fPDG = new TH1F("fPDG","pdg code", 2000, -1000, 1000);
  fMPDG = new TH1F("fMPDG","mother pdg code",2000,-1000,1000);
  fFPDG = new TH1F("fFPDG","father pdg code",2000,-1000,1000);

  fHistClustE = new TH1F("fHistClustE","Cluster Energy ;E (GeV) ;Entries",100,0,100);
  fHistClustEMatch = new TH1F("fHistClustEMatch","Cluster Energy after track matching (-1 < n#sigma < 3);E (GeV) ;Entries",2000,0,100);

  fHistClustLongAxis = new TH1F("fHistClustLongAxis","Events in M02;M02;Entries",100,0,2);
  fHistClustLongAxisE = new TH1F("fHistClustLongAxisE","Events in M02 (-1 < n#sigma < 3) (0.7 < E/P < 1.5) (P_{T} > 1.5);M02;Entries",100,0,2);

  fHistNsigmaP = new TH2F("fHistNsigmaP","n#sigma vs P ;P (GeV/c) ;n#sigma",1000,0,30,1000,-10,10);
  fHistMCNsigmaP = new TH2F("fHistMCNsigmaP","MC n#sigma vs P ;P (GeV/c) ;n#sigma",1000,0,30,1000,-10,10);

  fPtEoverPE = new TH2F("fPtvsEoverPE","P_{T} vs E/P (-1 < n#sigma < 3) ;P_{T} (GeV/c) ;E/P",1000,0,100,50,0,3);
  fPtEoverPMCE = new TH2F("fPtEoverPMCE","MC Events P_{T} vs E/P (-1 < n#sigma < 3) ;P_{T} (GeV/c) ;E/P",1000,0,100,50,0,3);
  fPtEoverPEGeo = new TH2F("fPtEoverPEGeo","P_{T} vs E/P (-1 < n#sigma < 3) (0.1 < M02 < 0.6) ;P_{T} (GeV/c) ;E/P",1000,0,100,50,0,3);
  fHistMCClsLAE = new TH1F("fHistMCClsLAE","MC Events in M02 (-1 < n#sigma < 3);M02;Entries",100,0,2);
  fPtEoverPH = new TH2F("fPtvsEoverPH","P_{T} vs E/P (n#sigma < -3) ;P_{T} (GeV/c) ;E/P",1000,0,100,50,0,3);
  fPtEoverPHGeo = new TH2F("fPtEoverPHGeo","P_{T} vs E/P (n#sigma < -3) (0.1 < M02 < 0.6) ;P_{T} (GeV/c) ;E/P",1000,0,100,50,0,3);

  fREisolation = new TH1F("fREisolation","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}",100,0,1);
  fREisolation->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisolation->GetYaxis()->SetTitle("Entries");

  fREisoW = new TH1F("fREisoW","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (W #rightarrow e) (P_{T} > 10GeV)",100,0,1);
  fREisoW->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisoW->GetYaxis()->SetTitle("Entries");
  fREisoHF = new TH1F("fREisoHF","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (b,c #rightarrow e) (P_{T} > 10GeV)",100,0,1);
  fREisoHF->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisoHF->GetYaxis()->SetTitle("Entries");
  fREisoWhpt = new TH1F("fREisoWhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (W #rightarrow e) (P_{T} > 20GeV)",100,0,1);
  fREisoWhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisoWhpt->GetYaxis()->SetTitle("Entries");
  fREisoHFhpt = new TH1F("fREisoHFhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (b,c #rightarrow e) (P_{T} > 20GeV)",100,0,1);
  fREisoHFhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisoHFhpt->GetYaxis()->SetTitle("Entries");
  fREisoWvhpt = new TH1F("fREisoWvhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (W #rightarrow e) (P_{T} > 30GeV)",100,0,1);
  fREisoWvhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisoWvhpt->GetYaxis()->SetTitle("Entries");
  fREisoHFvhpt = new TH1F("fREisoHFvhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (b,c #rightarrow e) (P_{T} > 30GeV)",100,0,1);
  fREisoHFvhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREisoHFvhpt->GetYaxis()->SetTitle("Entries");

  EleTraDelPhi = new TH2F("EleTraDelPhi","#Delta #phi = #phi_{track} - #phi_{e #leftarrow W}",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
  EleTraDelPhi->GetXaxis()->SetTitle("#Delta #phi (rad)");
  EleTraDelPhi->GetYaxis()->SetTitle("track P_{T}(GeV/c)");

  fHistPt_We = new TH1F("fHistPt_We","P_{T, e #leftarrow W}",200,0,200);
  fHistPt_We->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt_We->GetYaxis()->SetTitle("Counts");
  fHistPt_HFe = new TH1F("fHistPt_HFe","P_{T, e #leftarrow b,c}",200,0,200);
  fHistPt_HFe->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt_HFe->GetYaxis()->SetTitle("Counts");

  fPt_maxtrack_W = new TH1F("fPt_maxtrack_W","P_{T, trk} / P_{T, e #leftarrow W}", 200, 0, 2);
  fPt_maxtrack_W->GetXaxis()->SetTitle("P_{T, trk} / P_{T, e #leftarrow W}");
  fPt_maxtrack_W->GetYaxis()->SetTitle("Counts");
  fNtrkl_PtOfMaxTrk_W = new TH2F("fNtrkl_PtOfMaxTrk_W","N_{tracklets} vs P_{T, trk} / P_{T, e #leftarrow W}", 200, 0, 200, 200, 0, 2);
  fNtrkl_PtOfMaxTrk_W->GetXaxis()->SetTitle("N_{tracklets}");
  fNtrkl_PtOfMaxTrk_W->GetYaxis()->SetTitle("P_{T, trk} / P_{T, e #leftarrow W}");

  fHistPt_We_Ntrkl = new TH2F("fHistPt_We_Ntrkl","P_{T,e#leftarrowW} vs N_{tracklets}",200,0,200,200,0,200);
  fHistPt_We_Ntrkl->GetXaxis()->SetTitle("P_{T,e#leftarrowW} (GeV/c)");
  fHistPt_We_Ntrkl->GetYaxis()->SetTitle("N_{tracklets}");

  fW_true = new TH1F("fW_true","The true of e#leftarrowW; P_{T} (GeV/c); Entries",200,0,200);
  fW_false = new TH1F("fW_false","The false of e#leftarrowW; P_{T} (GeV/c); Entries",200,0,200);

  EleTraDelPhi_fullrange = new TH1F("EleTraDelPhi_fullrange","#Delta #phi full range; #Delta #phi; track",200,-7,7);

  fOutputList->Add(fHistPt);
  fOutputList->Add(pVertex);
  fOutputList->Add(cutVer);
  fOutputList->Add(EtavsPhi);
  fOutputList->Add(TPCSig);
  fOutputList->Add(Cent);
  fOutputList->Add(fzvtx_Ntrkl);
  fOutputList->Add(fzvtx_Ntrkl_cal);
  fOutputList->Add(fNtrklNch);
  fOutputList->Add(fPDG);
  fOutputList->Add(fMPDG);
  fOutputList->Add(fFPDG);

  fOutputList->Add(fHistClustE);
  fOutputList->Add(fHistClustEMatch);
  fOutputList->Add(fHistClustLongAxis);
  fOutputList->Add(fHistClustLongAxisE);

  fOutputList->Add(fHistNsigmaP);
  fOutputList->Add(fHistMCNsigmaP);
  fOutputList->Add(fPtEoverPE);
  fOutputList->Add(fPtEoverPMCE);
  fOutputList->Add(fPtEoverPEGeo);
  fOutputList->Add(fHistMCClsLAE);
  fOutputList->Add(fPtEoverPH);
  fOutputList->Add(fPtEoverPHGeo);
  fOutputList->Add(fREisolation);
  fOutputList->Add(fREisoW);
  fOutputList->Add(fREisoHF);
  fOutputList->Add(fREisoWhpt);
  fOutputList->Add(fREisoHFhpt);
  fOutputList->Add(fREisoWvhpt);
  fOutputList->Add(fREisoHFvhpt);
  fOutputList->Add(EleTraDelPhi);
  fOutputList->Add(fHistPt_We);
  fOutputList->Add(fHistPt_HFe);
  fOutputList->Add(fPt_maxtrack_W);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W);
  fOutputList->Add(fHistPt_We_Ntrkl);
  fOutputList->Add(fW_true);
  fOutputList->Add(fW_false);

  fOutputList->Add(EleTraDelPhi_fullrange);

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskWHMult::UserExec(Option_t *)
{
  //!!!!!!!!!!!//
  //this function is called once for each event
  //!!!!!!!!!!!//

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  if(!fAOD) return;

  if(fAOD)fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent) {
	  printf("ERROR: fVEvent not available\n");
	  return;
  }   

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
  }

  //fMCEvent = MCEvent();

  if(!gGrid){
    cout << "no Grid connection, connecting to the Grid ..." << endl;
    TGrid::Connect("alien//");
  }

  //==== trigger selection ======

  TString firedTrigger;
  TString TriggerEG1("EG1");

  fVevent->GetFiredTriggerClasses();
  if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
  if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}



  //=====reference objects file=====
  TFile* RefData=TFile::Open("alien:///alice/cern.ch/user/t/takawagu/Data.root");
  TTree* tree=(TTree*)RefData->Get("tree");
  //TTree* treee=(TTree*)RefData->Get("treee");

  Int_t iTracks(fAOD->GetNumberOfTracks());

  Int_t Nclust = fAOD->GetNumberOfCaloClusters();

  Double_t vertexZ = fAOD->GetPrimaryVertex()->GetZ();
  pVertex->Fill(vertexZ);
  if(abs(vertexZ)<=10.){
    cutVer->Fill(vertexZ);
  }
  const AliVVertex *pVtx = fAOD->GetPrimaryVertex();

  Float_t centrality(0);
  AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
  if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");	//get centrality
  Cent->Fill(centrality);

//---------------SPD tracklets---------------
  Int_t nTracklets = 0;
  Double_t nAcc = 0;		//number of tracklets in etarange
  Double_t corr_nAcc = 0;	//corrected with vertexZ
  Double_t etaRange = 1.0;

  AliAODTracklets* tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
  nTracklets = tracklets->GetNumberOfTracklets();
  for (Int_t nn=0;nn<nTracklets;nn++) {
    Double_t trkltheta = tracklets->GetTheta(nn);
    Double_t trkleta = -TMath::Log(TMath::Tan(trkltheta/2.0));
    if (TMath::Abs(trkleta)<etaRange) nAcc++;
  }
  fzvtx_Ntrkl->Fill(vertexZ,nAcc);

  //=====vertex Z correction=====
  if(abs(vertexZ)<=10.){
    Int_t BinOfvertexZ = fzvtx_Ntrkl->GetXaxis()->FindBin(vertexZ);
    corr_nAcc = GetCorrectedNtrackletsD(tree, nAcc, BinOfvertexZ);
    fzvtx_Ntrkl_cal->Fill(vertexZ,corr_nAcc);
  }

  Int_t Nch = 0;
  if(fMCarray){
    Nch = CountNch();
    if(abs(vertexZ)<=10.) fNtrklNch->Fill(corr_nAcc,Nch);
  }

//---------------track loop---------------
  for(Int_t i(0); i < iTracks; i++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));

    ////////////////////////
    // Get MC information //
    ////////////////////////
    Int_t ilabel = TMath::Abs(track->GetLabel());	// track label(charge) at MC event
    Int_t pdg = -999;					// initialize pdg code
    Double_t pid_ele = 0.0;				// initialize electron PID
    Double_t pTmom = -1.0;				// initialize Mother's pT
    Int_t pidM = -1;					// initialize Mother's PDG
    Int_t ilabelM = -1;					// initialize Mother's label(charge)
    Double_t pTdad = -1.0;				// initialize Father's pT
    Int_t pidF = -1;					// initialize Father's PDG
    Int_t ilabelF = -1;					// initialize Father's label(charge)
    Double_t pTGMom = -1.0;				// initialize Grand Mother's pT
    Int_t pidGM = -1;					// initialize Grand Mother's PDG
    Int_t ilabelGM = -1.0;				// initialize Grand Mother's label(charge)
    Bool_t pidW = 0;					// initialize W boson PID
    Bool_t pidHF = 0;					// initialize HF PID

    if (ilabel>0 && fMCarray) {
      fMCTrackpart = (AliAODMCParticle*) fMCarray->At(ilabel);	// define fMCarray: container of particles & connect data and MC by lable
      pdg = fMCTrackpart->GetPdgCode();				// get pdg code
      if (TMath::Abs(pdg) == 11) pid_ele = 1.0;			// find electron: pid_ele = 1.0
      //find electron mother's information from FindMother class
      if (pid_ele == 1.0) {
        FindMother(fMCTrackpart, ilabelM, pidM, pTmom);
        FindFather(fMCTrackpart, ilabelF, pidF, pTdad);
        //find W
        if (TMath::Abs(pidF) == 24) pidW = 1;
        //find HF
        if (TMath::Abs(pidM) == 411 || TMath::Abs(pidM) == 421 || TMath::Abs(pidM) == 413 || TMath::Abs(pidM) == 423 || TMath::Abs(pidM) == 431 || TMath::Abs(pidM) == 433 || TMath::Abs(pidM) == 511 || TMath::Abs(pidM) == 521 || TMath::Abs(pidM) == 513 || TMath::Abs(pidM) == 523 || TMath::Abs(pidM) == 531 || TMath::Abs(pidM) == 533) pidHF = 1;
      }
    }
    fPDG->Fill(pdg);
    fMPDG->Fill(pidM);
    fFPDG->Fill(pidF);
    //////////////////////////////////////
    // Found electron Mother infomation //
    //////////////////////////////////////

//---------------FilterBit 1---------------
    if(!track || !track->TestFilterBit(1)) continue;
    fHistPt->Fill(track->Pt());

  // Comparison with TPC analysis
    Double_t DCA[2] = {-999.,-999.}, covar[3];
    /////////////////////////
    //      track cut      //
    /////////////////////////
  //===== 1. TPC and ITS refit cut =====
    if (!(track->GetStatus()&AliAODTrack::kITSrefit) || !(track->GetStatus()&AliAODTrack::kTPCrefit)) continue;
  //===== 2. AOD filter bit required =====
    if (!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;		//minimum cut
  //===== 3. TPC cluster cut =====
    if (track->GetTPCNcls() < 80) continue;
  //===== 4. ITS cluster cut =====
    if (track->GetITSNcls() < 2) continue;
  //===== 5. SPD hit cut =====
    if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
  //===== 6. Eta cut =====
    if (track->Eta() > 0.9 || track->Eta() < -0.9) continue;
  //===== 7. DCA cut =====
    if (track -> PropagateToDCA(pVtx,fAOD -> GetMagneticField(),20.,DCA,covar))
    {
      if (TMath::Abs(DCA[0]) > 2.4 || TMath::Abs(DCA[1]) > 3.2) continue;
    }
  //===== 8. chi2 cut =====
    Double_t TPCchi2NDF = track->GetTPCchi2perNDF();
    Double_t ITSchi2 = track->GetITSchi2();
    if ((ITSchi2 >= 25) || (TPCchi2NDF >= 4)) continue;
  //===== 9. NCrossedRow cut =====
    if (track->GetTPCCrossedRows() < 100) continue;


//---------------vertex within 10cm---------------> event selection
    if(abs(vertexZ)<=10.)
    {
      double eta = -log(TMath::Tan((track->Theta())/2.));
      EtavsPhi->Fill(eta,track->Phi());

      TPCSig->Fill(track->P(), track->GetTPCsignal());

      Double_t fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
      fHistNsigmaP->Fill(track->P(),fTPCnSigma);
      if (pid_ele == 1.0) fHistMCNsigmaP->Fill(track->P(), fTPCnSigma);

//---------------track matching---------------
      Int_t EMCalIndex = -1;
      EMCalIndex = track->GetEMCALcluster();			//EMCal index of HIT tracks
      if (track->Eta() > 0.6 || track->Eta() < -0.6) continue;	//For EMCal acceptance
      AliVCluster*clustMatch = 0x0;
      if (EMCalIndex >= 0)
      {
        clustMatch = (AliVCluster*)fAOD->GetCaloCluster(EMCalIndex);    //information of cluster matching with tracks
        Double_t clustEmatch = clustMatch->E();				//information of cluster energy matching with tracks
        fHistClustEMatch->Fill(clustEmatch);

        Double_t clustLongE = clustMatch->GetM02();
        fHistClustLongAxis->Fill(clustLongE);

        Double_t EoverP = clustEmatch/track->P();

//---------------nsigma cut (electron)---------------
        if (fTPCnSigma > -1. && fTPCnSigma < 3.)
        {
          fPtEoverPE->Fill(track->Pt(),EoverP);
          if (TMath::Abs(pdg) == 11) fPtEoverPMCE->Fill(track->Pt(),EoverP);

//---------------EoverP cut (electron)---------------
          if (EoverP > 0.7 && EoverP < 1.5)
          {
            if (TMath::Abs(pdg) == 11) fHistMCClsLAE->Fill(clustLongE);
            fHistClustLongAxisE->Fill(clustLongE);

//---------------M02 cut (electron)---------------THIS IS MOSTLY ELECTRON
            if (clustLongE > 0.1 && clustLongE < 0.6)
            {
              fPtEoverPEGeo->Fill(track->Pt(),EoverP);
              Float_t showerx[3];
              clustMatch->GetPosition(showerx);
              TVector3 Corepos(showerx[0],showerx[1],showerx[2]);
              Double_t Corephi = Corepos.Phi();
              Double_t Coreeta = Corepos.Eta();
              Double_t eleE = clustMatch->E();
              Double_t RsumE = 0;
//---------------Isolation cut EMCal Cluster loop start---------------
              for(Int_t icl=0;icl<Nclust;icl++){
                AliVCluster*clust = 0x0;
                clust = (AliVCluster*)fAOD->GetCaloCluster(icl);
                if(clust && clust->IsEMCAL()){ 
                  Double_t clustE = clust->E();
                  fHistClustE->Fill(clustE);
                  Float_t aroclsx[3];
                  clust->GetPosition(aroclsx);
                  TVector3 aroClsPos(aroclsx[0],aroclsx[1],aroclsx[2]);
                  Double_t aroClsphi = aroClsPos.Phi();
                  Double_t aroClseta = aroClsPos.Eta();
                  Double_t R = sqrt(pow(Corephi - aroClsphi,2.0)+pow(Coreeta - aroClseta,2.0));
                  if (R < 0.3) {
                    RsumE = RsumE + clustE;
                  }
                }
              }		//cluster loop end
              Double_t Eiso = (RsumE - eleE)/eleE;
              fREisolation->Fill(Eiso);
              //=====MC Data=====
              if (track->Pt() > 10.) {
                if (pidW == 1) fREisoW->Fill(Eiso);
                if (pidHF == 1) fREisoHF->Fill(Eiso);
              }
              if (track->Pt() > 20.) {
                if (pidW == 1) fREisoWhpt->Fill(Eiso);
                if (pidHF == 1) fREisoHFhpt->Fill(Eiso);
              }
              if (track->Pt() > 30.) {
                if (pidW == 1) fREisoWvhpt->Fill(Eiso);
                if (pidHF == 1) fREisoHFvhpt->Fill(Eiso);
              }
              //=================
//---------------isolation cut-----------------------
              if (Eiso >= 0.1) fHistPt_HFe->Fill(track->Pt());
              if (Eiso >= 0. && Eiso <= 0.05) {
                fHistPt_We->Fill(track->Pt());
//---------------Pt cut---------------
                if (track->Pt() > 10.)		//Pt identify W and HF
                {
                  if (pidW == 1) fW_true->Fill(track->Pt());
                  if (pidW != 1) fW_false->Fill(track->Pt());
                  //////////////////
                  // e <- W EVENT //
                  //////////////////
//---------------another track loop start------------------
                  Int_t MaxPtTrackNum = 0;
                  if (MaxPtTrackNum == i) MaxPtTrackNum = 1;
                  for(Int_t j(0); j < iTracks; j++) {
                    AliAODTrack* anotrack = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
                    if(!anotrack || !anotrack->TestFilterBit(1)) continue;
                    if (j == i) continue;	//reject e<-W candidate
                    
                    /////////////////////////
                    //      track cut      //
                    /////////////////////////
                    //===== 1. TPC and ITS refit cut =====
                    if (!(anotrack->GetStatus()&AliAODTrack::kITSrefit) || !(anotrack->GetStatus()&AliAODTrack::kTPCrefit)) continue;
                    //===== 2. AOD filter bit required =====
                    if (!anotrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
                    //===== 3. TPC cluster cut =====
                    if (anotrack->GetTPCNcls() < 80) continue;
                    //===== 4. ITS cluster cut =====
                    if (anotrack->GetITSNcls() < 2) continue;
                    //===== 5. SPD hit cut =====
                    if (!(anotrack->HasPointOnITSLayer(0) || anotrack->HasPointOnITSLayer(1))) continue;
                    //===== 6. Eta cut =====
                    if (anotrack->Eta() > 0.9 || anotrack->Eta() < -0.9) continue;
                    //===== 7. DCA cut =====
                    Double_t anoDCA[2] = {-999.,-999.}, anocovar[3];
                    if (anotrack -> PropagateToDCA(pVtx,fAOD -> GetMagneticField(),20.,anoDCA,anocovar))
                    {
                      if (TMath::Abs(anoDCA[0]) > 2.4 || TMath::Abs(anoDCA[1]) > 3.2) continue;
                    }
                    //===== 8. chi2 cut =====
                    Double_t anoITSchi2 = anotrack->GetITSchi2();
                    Double_t anoTPCchi2NDF = anotrack->GetTPCchi2perNDF();
                    if ((anoITSchi2 >= 25) || (anoTPCchi2NDF >= 4)) continue;
                    //===== 9. NCrossedRow cut =====
                    if (anotrack->GetTPCCrossedRows() < 100) continue;

                    Double_t anoeta = -log(TMath::Tan((anotrack->Theta())/2.));
                    Double_t DelPhi = (anotrack->Phi()) - (track->Phi());
                    EleTraDelPhi_fullrange->Fill(DelPhi);
                    if (DelPhi < -1.*TMath::Pi()/3.) DelPhi = DelPhi + 2.*TMath::Pi();
                    if (DelPhi > 5.*TMath::Pi()/3.) DelPhi = DelPhi - 2.*TMath::Pi();
                    EleTraDelPhi->Fill(DelPhi,anotrack->Pt());

                    //=====Highest Pt Track=====
                    if (DelPhi >= 5.*TMath::Pi()/6. && DelPhi <= 7.*TMath::Pi()/6.) {
                      AliAODTrack* Maxtrack = static_cast<AliAODTrack*>(fAOD->GetTrack(MaxPtTrackNum));
                      Double_t MaxDelPhi = (Maxtrack->Phi()) - (track->Phi());
                      if (MaxDelPhi < -1.*TMath::Pi()/3.) MaxDelPhi = MaxDelPhi + 2.*TMath::Pi();
                      if (MaxDelPhi > 5.*TMath::Pi()/3.) MaxDelPhi = MaxDelPhi - 2.*TMath::Pi();
                      if (MaxDelPhi < 5.*TMath::Pi()/6. || MaxDelPhi > 7.*TMath::Pi()/6.) MaxPtTrackNum = j;
                      if (anotrack->Pt() >= Maxtrack->Pt()) MaxPtTrackNum = j;
                    }
                  }	//other track loop end

                  //=====Highest Pt Track / e<-W Pt Track
                  AliAODTrack* MaxPtTrk = static_cast<AliAODTrack*>(fAOD->GetTrack(MaxPtTrackNum));
                  fPt_maxtrack_W->Fill(MaxPtTrk->Pt()/track->Pt());
                  fNtrkl_PtOfMaxTrk_W->Fill(corr_nAcc,MaxPtTrk->Pt()/track->Pt());

                  fHistPt_We_Ntrkl->Fill(track->Pt(),corr_nAcc);

////////still=== e <- W ===////////
                }	//Pt cut
              }		//isolation cut
            }		//M02 cut
          }		//EoverP cut
        }		//nsigma cut

//---------------nsigma cut (hadron)---------------
        if (fTPCnSigma < -3.)
        {
          fPtEoverPH->Fill(track->Pt(),EoverP);
//---------------M02 cut---------------
          if (clustLongE > 0.1 && clustLongE < 0.6) fPtEoverPHGeo->Fill(track->Pt(),EoverP);
        }	//hadron filter end
      }		//track match end
    }		//Zvertex cut end
  }		//track loop end

  PostData(1, fOutputList);
}
//----------MC Mother's Particle Information----------
void AliAnalysisTaskWHMult::FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom)
{
  if (part->GetMother() > -1) {
    label = part->GetMother();
    AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
    //part = (AliAODMCParticle*)fMCarray->At(label);
    pid = partM->GetPdgCode();
    ptmom = partM->Pt();
  }
  //cout << "Find Mother : label = " << label << " ; pid" << pid << endl;
}
//----------MC Mother's W Information. Father is Origin Mather----------
void AliAnalysisTaskWHMult::FindFather(AliAODMCParticle* part, int &label, int &pid, double &ptdad)
{
  while (part->GetMother() > 0) {
    label = part->GetMother();
    //AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
    part = (AliAODMCParticle*)fMCarray->At(label);
    pid = part->GetPdgCode();
    ptdad = part->Pt();
  }
  //cout << "Find Mother : label = " << label << " ; pid" << pid << endl;
}
//----------vertexZ correction for the number of tracklets----------
Double_t AliAnalysisTaskWHMult::GetCorrectedNtrackletsD(TTree* tree, Double_t uncorrectedNacc, Int_t BinOfvertexZ)
{
  //if(BinOfvertexZ<1 || BinOfvertexZ > 100) return uncorrectedNacc;
  if(!tree) return uncorrectedNacc;

  Int_t TreeEntry = BinOfvertexZ - 1;
  Double_t clb_vZ_Ntl = 0.;
  tree->SetBranchAddress("vertexZ_Ntrkl",&clb_vZ_Ntl);
  tree->GetEntry(TreeEntry);

  Double_t deltaM = 0;
  deltaM = uncorrectedNacc*(clb_vZ_Ntl - 1.);

  Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));
  if (correctedNacc<0) correctedNacc = 0;

  return correctedNacc;
}
//----------count number of charged particles generated MC----------
Int_t AliAnalysisTaskWHMult::CountNch()
{
  Int_t Nch = 0;
  for (int iMCpart = 0;iMCpart < fMCarray->GetEntriesFast();iMCpart++) {
    fMCparticle = (AliAODMCParticle*)fMCarray->At(iMCpart);		//call all generated particles in fMCarray of each event
    Int_t chargetrue = fMCparticle->Charge();				//get charge
    Double_t pdgEta = fMCparticle->Eta();				//acceptance
    Bool_t PrimaryParticle = fMCparticle->IsPhysicalPrimary();		//choose primary particle
    if(chargetrue!=0){
      if(TMath::Abs(pdgEta)<1.0){
        if(PrimaryParticle){
          Nch++;
        }
      }
    }
  }
  return Nch;
}

void AliAnalysisTaskWHMult::Terminate(Option_t *)
{
  // terminate
}
