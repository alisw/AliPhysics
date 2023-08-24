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
  fVevent(0),
  fPIDResponse(0),
  fMCTrackpart(0),
  fMCarray(0),
  fMCparticle(0),
  fOutputList(0),
  tree(0),
  fNevents(0),
  fHistPt(0),
  pVertex_all(0),
  pVertex(0),
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
  fREisolation(),
  fREiso_MCW(0),
  fREiso_MCHF(0),
  fREiso_MCWhpt(0),
  fREiso_MCHFhpt(0),
  fREiso_MCWvhpt(0),
  fREiso_MCHFvhpt(0),
  fdPhi_trkW_Pt(),
  fdPhi_trkHF_Pt(),
  fdPhi_trkW_Pt_hpt(),
  fdPhi_trkHF_Pt_hpt(),
  fdPhi_trkW_ePt(),
  fdPhi_trkHF_ePt(),
  fHistPt_We(),
  fHistPt_HFe(),
  fPt_maxtrack_W(),
  fPt_maxtrack_W_hpt(),
  fNtrkl_PtOfMaxTrk_W(),
  fNtrkl_PtOfMaxTrk_W_hpt(),
  fHistPt_We_Ntrkl(),
  fdPhi_trkW_full(),
  fdPhi_trkHF_full(),
  fdPhi_trkW_full_hpt(),
  fdPhi_trkHF_full_hpt(),
  fNtrkl_ClustE(0),
  fEMCEG1(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskWHMult::AliAnalysisTaskWHMult(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),
  fVevent(0),
  fPIDResponse(0),
  fMCTrackpart(0),
  fMCarray(0),
  fMCparticle(0),
  fOutputList(0),
  tree(0),
  fNevents(0),
  fHistPt(0),
  pVertex_all(0),
  pVertex(0),
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
  fREisolation(),
  fREiso_MCW(0),
  fREiso_MCHF(0),
  fREiso_MCWhpt(0),
  fREiso_MCHFhpt(0),
  fREiso_MCWvhpt(0),
  fREiso_MCHFvhpt(0),
  fdPhi_trkW_Pt(),
  fdPhi_trkHF_Pt(),
  fdPhi_trkW_Pt_hpt(),
  fdPhi_trkHF_Pt_hpt(),
  fdPhi_trkW_ePt(),
  fdPhi_trkHF_ePt(),
  fHistPt_We(),
  fHistPt_HFe(),
  fPt_maxtrack_W(),
  fPt_maxtrack_W_hpt(),
  fNtrkl_PtOfMaxTrk_W(),
  fNtrkl_PtOfMaxTrk_W_hpt(),
  fHistPt_We_Ntrkl(),
  fdPhi_trkW_full(),
  fdPhi_trkHF_full(),
  fdPhi_trkW_full_hpt(),
  fdPhi_trkHF_full_hpt(),
  fNtrkl_ClustE(0),
  fEMCEG1(kFALSE)
{
  // constructor
  DefineInput(0, TChain::Class());

  DefineOutput(1, TList::Class()); 
}
//_____________________________________________________________________________
AliAnalysisTaskWHMult::~AliAnalysisTaskWHMult()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskWHMult::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fNevents = new TH1F("fNevents","Number of events",8,-0.5,7.5);

  fHistPt = new TH1F("fHistPt", "fHistPt", 200, 0, 200);

  pVertex_all = new TH1F("pVertex_all", "pVertex_all", 100, -20, 20);
  pVertex_all->GetXaxis()->SetTitle("Z vertex (collision point on z-axis)");
  pVertex_all->GetYaxis()->SetTitle("events");

  pVertex = new TH1F("pVertex", "pVertex", 100, -20, 20);
  pVertex->GetXaxis()->SetTitle("Z vertex (collision point on z-axis)");
  pVertex->GetYaxis()->SetTitle("events");

  EtavsPhi = new TH2F("EtavsPhi","#eta vs #phi", 500, -1.5, 1.5, 500, -1, 7);
  EtavsPhi->GetXaxis()->SetTitle("#eta");
  EtavsPhi->GetYaxis()->SetTitle("#phi (rad)");

  TPCSig = new TH2F("TPCSig","TPC signal", 5000, 0.1, 30, 10000, 0, 10000);
  TPCSig->GetXaxis()->SetTitle("p (GeV/c)");
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
  fHistClustEMatch = new TH1F("fHistClustEMatch","Cluster Energy after track matching ;E (GeV) ;Entries",2000,0,100);

  fHistClustLongAxis = new TH1F("fHistClustLongAxis","Events in M02;M02;Entries",100,0,2);
  fHistClustLongAxisE = new TH1F("fHistClustLongAxisE","Events in M02 (-1 < n#sigma < 3) (0.7 < E/p < 1.5) (P_{T} > 1.5);M02;Entries",100,0,2);

  fHistNsigmaP = new TH2F("fHistNsigmaP","n#sigma vs p ;p (GeV/c) ;n#sigma",1000,0,30,1000,-10,10);
  fHistMCNsigmaP = new TH2F("fHistMCNsigmaP","MC n#sigma vs p ;p (GeV/c) ;n#sigma",1000,0,30,1000,-10,10);

  fPtEoverPE = new TH2F("fPtvsEoverPE","p_{T} vs E/p (-1 < n#sigma < 3) ;p_{T} (GeV/c) ;E/p",1000,0,100,100,0,3);
  fPtEoverPMCE = new TH2F("fPtEoverPMCE","MC Events p_{T} vs E/p (-1 < n#sigma < 3) ;p_{T} (GeV/c) ;E/p",1000,0,100,100,0,3);
  fPtEoverPEGeo = new TH2F("fPtEoverPEGeo","p_{T} vs E/p (-1 < n#sigma < 3) (0.1 < M02 < 0.6) ;p_{T} (GeV/c) ;E/p",1000,0,100,100,0,3);
  fHistMCClsLAE = new TH1F("fHistMCClsLAE","MC Events in M02 (-1 < n#sigma < 3);M02;Entries",100,0,2);
  fPtEoverPH = new TH2F("fPtvsEoverPH","p_{T} vs E/p (n#sigma < -3) ;p_{T} (GeV/c) ;E/p",1000,0,100,100,0,3);
  fPtEoverPHGeo = new TH2F("fPtEoverPHGeo","p_{T} vs E/p (n#sigma < -3) (0.1 < M02 < 0.6) ;p_{T} (GeV/c) ;E/p",1000,0,100,100,0,3);

  for (Int_t isoR=0;isoR<3;isoR++) {
    fREisolation[isoR] = new TH1F(Form("fREisolation_%d",isoR),Form("(#sum_{R<0.%d}E_{shower}-E_{electron})/E_{electron}",3+isoR),100,0,1);
    fREisolation[isoR]->GetXaxis()->SetTitle(Form("(#sum_{R<0.%d}E_{shower}-E_{electron})/E_{electron}",3+isoR));
    fREisolation[isoR]->GetYaxis()->SetTitle("Entries");
  }

  fREiso_MCW = new TH1F("fREiso_MCW","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (e #leftarrow W) (p_{T} > 10GeV)",100,0,1);
  fREiso_MCW->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREiso_MCW->GetYaxis()->SetTitle("Entries");
  fREiso_MCHF = new TH1F("fREiso_MCHF","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (e #leftarrow b,c) (p_{T} > 10GeV)",100,0,1);
  fREiso_MCHF->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREiso_MCHF->GetYaxis()->SetTitle("Entries");
  fREiso_MCWhpt = new TH1F("fREiso_MCWhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (e #leftarrow W) (p_{T} > 20GeV)",100,0,1);
  fREiso_MCWhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREiso_MCWhpt->GetYaxis()->SetTitle("Entries");
  fREiso_MCHFhpt = new TH1F("fREiso_MCHFhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (e #leftarrow b,c) (p_{T} > 20GeV)",100,0,1);
  fREiso_MCHFhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREiso_MCHFhpt->GetYaxis()->SetTitle("Entries");
  fREiso_MCWvhpt = new TH1F("fREiso_MCWvhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (e #leftarrow W) (p_{T} > 30GeV)",100,0,1);
  fREiso_MCWvhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREiso_MCWvhpt->GetYaxis()->SetTitle("Entries");
  fREiso_MCHFvhpt = new TH1F("fREiso_MCHFvhpt","(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron} (e #leftarrow b,c) (p_{T} > 30GeV)",100,0,1);
  fREiso_MCHFvhpt->GetXaxis()->SetTitle("(#sum_{R<0.3}E_{shower}-E_{electron})/E_{electron}");
  fREiso_MCHFvhpt->GetYaxis()->SetTitle("Entries");

  for (Int_t isoR=0;isoR<3;isoR++) {
    fdPhi_trkW_Pt[isoR] = new TH2F(Form("fdPhi_trkW_Pt_%d",isoR),"",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
    fdPhi_trkW_Pt[isoR]->SetTitle(Form("#Delta #phi = #phi_{trk}-#phi_{can} (Eiso_{R<0.%d})",3+isoR));
    fdPhi_trkW_Pt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkW_Pt[isoR]->GetYaxis()->SetTitle("p_{T,trk}(GeV/c)");
    fdPhi_trkHF_Pt[isoR] = new TH2F(Form("fdPhi_trkHF_Pt_%d",isoR),"",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
    fdPhi_trkHF_Pt[isoR]->SetTitle(Form("#Delta #phi = #phi_{trk}-#phi_{HFcan} (Eiso_{R<0.%d})",3+isoR));
    fdPhi_trkHF_Pt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkHF_Pt[isoR]->GetYaxis()->SetTitle("p_{T,trk}(GeV/c)");
    fdPhi_trkW_Pt_hpt[isoR] = new TH2F(Form("fdPhi_trkW_Pt_hpt_%d",isoR),"",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
    fdPhi_trkW_Pt_hpt[isoR]->SetTitle(Form("#Delta #phi = #phi_{trk}-#phi_{can} (Eiso_{R<0.%d}) (p_{T,can}>30GeV/c)",3+isoR));
    fdPhi_trkW_Pt_hpt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkW_Pt_hpt[isoR]->GetYaxis()->SetTitle("p_{T,trk}(GeV/c)");
    fdPhi_trkHF_Pt_hpt[isoR] = new TH2F(Form("fdPhi_trkHF_Pt_hpt_%d",isoR),"",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
    fdPhi_trkHF_Pt_hpt[isoR]->SetTitle(Form("#Delta #phi = #phi_{trk}-#phi_{HFcan} (Eiso_{R<0.%d}) (p_{T,HFcan}>30GeV/c)",3+isoR));
    fdPhi_trkHF_Pt_hpt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkHF_Pt_hpt[isoR]->GetYaxis()->SetTitle("p_{T,trk}(GeV/c)");

    fdPhi_trkW_ePt[isoR] = new TH2F(Form("fdPhi_trkW_ePt_%d",isoR),"",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
    fdPhi_trkW_ePt[isoR]->SetTitle(Form("#Delta #phi = #phi_{trk}-#phi_{can} (Eiso_{R<0.%d})",3+isoR));
    fdPhi_trkW_ePt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkW_ePt[isoR]->GetYaxis()->SetTitle("p_{T,can}(GeV/c)");
    fdPhi_trkHF_ePt[isoR] = new TH2F(Form("fdPhi_trkHF_ePt_%d",isoR),"",50,-TMath::Pi()/3.,5.*TMath::Pi()/3.,1000,0,100);
    fdPhi_trkHF_ePt[isoR]->SetTitle(Form("#Delta #phi = #phi_{trk}-#phi_{HFcan} (Eiso_{R<0.%d})",3+isoR));
    fdPhi_trkHF_ePt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkHF_ePt[isoR]->GetYaxis()->SetTitle("p_{T,HFcan}(GeV/c)");

    fHistPt_We[isoR] = new TH1F(Form("fHistPt_We_%d",isoR),Form("p_{T,can} (Eiso_{R<0.%d})",3+isoR),200,0,200);
    fHistPt_We[isoR]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPt_We[isoR]->GetYaxis()->SetTitle("Counts");
    fHistPt_HFe[isoR] = new TH1F(Form("fHistPt_HFe_%d",isoR),Form("p_{T,HFcan} (Eiso_{R<0.%d})",3+isoR),200,0,200);
    fHistPt_HFe[isoR]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPt_HFe[isoR]->GetYaxis()->SetTitle("Counts");

    fPt_maxtrack_W[isoR] = new TH1F(Form("fPt_maxtrack_W_%d",isoR),Form("p_{T,trk}/p_{T,can} (Eiso_{R<0.%d})",3+isoR),200,0,2);
    fPt_maxtrack_W[isoR]->GetXaxis()->SetTitle("p_{T,trk}/p_{T,can}");
    fPt_maxtrack_W[isoR]->GetYaxis()->SetTitle("Counts");
    fPt_maxtrack_W_hpt[isoR] = new TH1F(Form("fPt_maxtrack_W_hpt_%d",isoR),Form("p_{T,trk}/p_{T,can} (Eiso_{R<0.%d}) (p_{T,can}>30GeV/c)",3+isoR),200,0,2);
    fPt_maxtrack_W_hpt[isoR]->GetXaxis()->SetTitle("p_{T,trk}/p_{T,can}");
    fPt_maxtrack_W_hpt[isoR]->GetYaxis()->SetTitle("Counts");

    fNtrkl_PtOfMaxTrk_W[isoR] = new TH2F(Form("fNtrkl_PtOfMaxTrk_W_%d",isoR),"",200,0,200,200,0,2);
    fNtrkl_PtOfMaxTrk_W[isoR]->SetTitle(Form("N_{tracklets} vs p_{T,trk}/p_{T,can} (Eiso_{R<0.%d})",3+isoR));
    fNtrkl_PtOfMaxTrk_W[isoR]->GetXaxis()->SetTitle("N_{tracklets}");
    fNtrkl_PtOfMaxTrk_W[isoR]->GetYaxis()->SetTitle("p_{T,trk}/p_{T,can}");
    fNtrkl_PtOfMaxTrk_W_hpt[isoR] = new TH2F(Form("fNtrkl_PtOfMaxTrk_W_hpt_%d",isoR),"",200,0,200,200,0,2);
    fNtrkl_PtOfMaxTrk_W_hpt[isoR]->SetTitle(Form("N_{tracklets} vs p_{T,trk}/p_{T,can} (Eiso_{R<0.%d}) (p_{T,can}>30GeV/c)",3+isoR));
    fNtrkl_PtOfMaxTrk_W_hpt[isoR]->GetXaxis()->SetTitle("N_{tracklets}");
    fNtrkl_PtOfMaxTrk_W_hpt[isoR]->GetYaxis()->SetTitle("p_{T,trk}/p_{T,can}");

    fHistPt_We_Ntrkl[isoR] = new TH2F(Form("fHistPt_We_Ntrkl_%d",isoR),Form("p_{T,can} vs N_{tracklets} (Eiso_{R<0.%d})",3+isoR),200,0,200,200,0,200);
    fHistPt_We_Ntrkl[isoR]->GetXaxis()->SetTitle("p_{T,can} (GeV/c)");
    fHistPt_We_Ntrkl[isoR]->GetYaxis()->SetTitle("N_{tracklets}");

    fdPhi_trkW_full[isoR] = new TH2F(Form("fdPhi_trkW_Pt_full_%d",isoR),Form("raw #Delta #phi (Eiso_{R<0.%d})",3+isoR),200,-7,7,200,0,200);
    fdPhi_trkW_full[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkW_full[isoR]->GetYaxis()->SetTitle("p_{T,trk} (GeV/c)");
    fdPhi_trkHF_full[isoR] = new TH2F(Form("fdPhi_trkHF_Pt_full_%d",isoR),Form("raw #Delta #phi (Eiso_{R<0.%d})",3+isoR),200,-7,7,200,0,200);
    fdPhi_trkHF_full[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkHF_full[isoR]->GetYaxis()->SetTitle("p_{T,trk} (GeV/c)");
    fdPhi_trkW_full_hpt[isoR] = new TH2F(Form("fdPhi_trkW_Pt_full_hpt_%d",isoR),"",200,-7,7,200,0,200);
    fdPhi_trkW_full_hpt[isoR]->SetTitle(Form("raw #Delta #phi (Eiso_{R<0.%d}) (p_{T,can}>30GeV/c)",3+isoR));
    fdPhi_trkW_full_hpt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkW_full_hpt[isoR]->GetYaxis()->SetTitle("p_{T,trk} (GeV/c)");
    fdPhi_trkHF_full_hpt[isoR] = new TH2F(Form("fdPhi_trkHF_Pt_full_hpt_%d",isoR),"",200,-7,7,200,0,200);
    fdPhi_trkHF_full_hpt[isoR]->SetTitle(Form("raw #Delta #phi (Eiso_{R<0.%d}) (p_{T,HFcan}>30GeV/c)",3+isoR));
    fdPhi_trkHF_full_hpt[isoR]->GetXaxis()->SetTitle("#Delta #phi (rad)");
    fdPhi_trkHF_full_hpt[isoR]->GetYaxis()->SetTitle("p_{T,trk} (GeV/c)");
  }

  fNtrkl_ClustE = new TH2F("fNtrkl_ClustE","N_{tracklets} vs cluster energy; N_{tracklets}; E (GeV)",200,0,200,100,0,50);


  fOutputList->Add(fNevents);
  fOutputList->Add(fHistPt);
  fOutputList->Add(pVertex_all);
  fOutputList->Add(pVertex);
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
  fOutputList->Add(fREisolation[0]);
  fOutputList->Add(fREisolation[1]);
  fOutputList->Add(fREisolation[2]);
  fOutputList->Add(fREiso_MCW);
  fOutputList->Add(fREiso_MCHF);
  fOutputList->Add(fREiso_MCWhpt);
  fOutputList->Add(fREiso_MCHFhpt);
  fOutputList->Add(fREiso_MCWvhpt);
  fOutputList->Add(fREiso_MCHFvhpt);
  fOutputList->Add(fdPhi_trkW_Pt[0]);
  fOutputList->Add(fdPhi_trkW_Pt[1]);
  fOutputList->Add(fdPhi_trkW_Pt[2]);
  fOutputList->Add(fdPhi_trkHF_Pt[0]);
  fOutputList->Add(fdPhi_trkHF_Pt[1]);
  fOutputList->Add(fdPhi_trkHF_Pt[2]);
  fOutputList->Add(fdPhi_trkW_Pt_hpt[0]);
  fOutputList->Add(fdPhi_trkW_Pt_hpt[1]);
  fOutputList->Add(fdPhi_trkW_Pt_hpt[2]);
  fOutputList->Add(fdPhi_trkHF_Pt_hpt[0]);
  fOutputList->Add(fdPhi_trkHF_Pt_hpt[1]);
  fOutputList->Add(fdPhi_trkHF_Pt_hpt[2]);
  fOutputList->Add(fdPhi_trkW_ePt[0]);
  fOutputList->Add(fdPhi_trkW_ePt[1]);
  fOutputList->Add(fdPhi_trkW_ePt[2]);
  fOutputList->Add(fdPhi_trkHF_ePt[0]);
  fOutputList->Add(fdPhi_trkHF_ePt[1]);
  fOutputList->Add(fdPhi_trkHF_ePt[2]);
  fOutputList->Add(fHistPt_We[0]);
  fOutputList->Add(fHistPt_We[1]);
  fOutputList->Add(fHistPt_We[2]);
  fOutputList->Add(fHistPt_HFe[0]);
  fOutputList->Add(fHistPt_HFe[1]);
  fOutputList->Add(fHistPt_HFe[2]);
  fOutputList->Add(fPt_maxtrack_W[0]);
  fOutputList->Add(fPt_maxtrack_W[1]);
  fOutputList->Add(fPt_maxtrack_W[2]);
  fOutputList->Add(fPt_maxtrack_W_hpt[0]);
  fOutputList->Add(fPt_maxtrack_W_hpt[1]);
  fOutputList->Add(fPt_maxtrack_W_hpt[2]);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W[0]);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W[1]);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W[2]);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W_hpt[0]);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W_hpt[1]);
  fOutputList->Add(fNtrkl_PtOfMaxTrk_W_hpt[2]);
  fOutputList->Add(fHistPt_We_Ntrkl[0]);
  fOutputList->Add(fHistPt_We_Ntrkl[1]);
  fOutputList->Add(fHistPt_We_Ntrkl[2]);
  fOutputList->Add(fdPhi_trkW_full[0]);
  fOutputList->Add(fdPhi_trkW_full[1]);
  fOutputList->Add(fdPhi_trkW_full[2]);
  fOutputList->Add(fdPhi_trkHF_full[0]);
  fOutputList->Add(fdPhi_trkHF_full[1]);
  fOutputList->Add(fdPhi_trkHF_full[2]);
  fOutputList->Add(fdPhi_trkW_full_hpt[0]);
  fOutputList->Add(fdPhi_trkW_full_hpt[1]);
  fOutputList->Add(fdPhi_trkW_full_hpt[2]);
  fOutputList->Add(fdPhi_trkHF_full_hpt[0]);
  fOutputList->Add(fdPhi_trkHF_full_hpt[1]);
  fOutputList->Add(fdPhi_trkHF_full_hpt[2]);
  fOutputList->Add(fNtrkl_ClustE);

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
  if(!tree)
    {
     cout << "no tree for correcting mult" << endl; 
     cout << "get from Grid" << endl; 
     TFile* RefData=TFile::Open("alien:///alice/cern.ch/user/t/takawagu/Data.root");
     tree=(TTree*)RefData->Get("tree");
     //TTree* treee=(TTree*)RefData->Get("treee");
    }
 
  Int_t iTracks(fAOD->GetNumberOfTracks());

  Int_t Nclust = fAOD->GetNumberOfCaloClusters();

  //===== Global Vtx =====
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t NcontV = pVtx->GetNContributors();
  Double_t Zvertex = pVtx->GetZ();
  pVertex_all->Fill(Zvertex);
  //===== SPD Vtx =====
  const AliVVertex *pVtxSPD = fVevent->GetPrimaryVertexSPD();
  Double_t ZvertexSPD = pVtxSPD->GetZ();
  Double_t NcontVSPD = pVtxSPD->GetNContributors();
  Double_t cov[6] = {0};
  pVtxSPD->GetCovarianceMatrix(cov);

  Float_t centrality(0);
  AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
  if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");	//get centrality
  Cent->Fill(centrality);

    ////////////////////////
    //  Event Selection   //
    ////////////////////////
  fNevents->Fill(0);
  // ===== 1. remove pile up events =====
  if (fVevent->IsPileupFromSPDInMultBins()) return;
  fNevents->Fill(1);
  // ===== 2. Global contributors cut =====
  if (NcontV < 2) return;
  fNevents->Fill(2);
  // ===== 3. SPD contributors cut =====
  if (NcontVSPD < 1) return;
  fNevents->Fill(3);
  // ===== 4. select events where SPD and primary vertex match =====
  if (TMath::Abs(ZvertexSPD - Zvertex) > 0.5) return;
  fNevents->Fill(4);
  // ===== 5. SPD vertex resolution cut =====
  if (TMath::Sqrt(cov[5]) > 0.25) return;
  fNevents->Fill(5);
  // ===== 6. Z vertex position cut =====
  if (TMath::Abs(Zvertex) > 10.) return;
  fNevents->Fill(6);
  pVertex->Fill(Zvertex);

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
  fzvtx_Ntrkl->Fill(Zvertex,nAcc);

  //=====vertex Z correction=====
  Int_t BinOfvertexZ = fzvtx_Ntrkl->GetXaxis()->FindBin(Zvertex);
  corr_nAcc = GetCorrectedNtrackletsD(tree, nAcc, BinOfvertexZ, fAOD);
  fzvtx_Ntrkl_cal->Fill(Zvertex,corr_nAcc);

  Int_t Nch = 0;
  if(fMCarray){
    Nch = CountNch();
    fNtrklNch->Fill(corr_nAcc,Nch);
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


    double eta = -log(TMath::Tan((track->Theta())/2.));
    EtavsPhi->Fill(eta,track->Phi());

    TPCSig->Fill(track->P(), track->GetTPCsignal());

    Double_t fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    fHistNsigmaP->Fill(track->P(),fTPCnSigma);
    if (pid_ele == 1.0) fHistMCNsigmaP->Fill(track->P(), fTPCnSigma);

//---------------track matching---------------
    Int_t EMCalIndex = -1;
    EMCalIndex = track->GetEMCALcluster();
    if (track->Eta() > 0.6 || track->Eta() < -0.6) continue;
    AliVCluster*clustMatch = 0x0;
    if (EMCalIndex >= 0)
    {
      clustMatch = (AliVCluster*)fAOD->GetCaloCluster(EMCalIndex);
      Double_t clustEmatch = clustMatch->E();
      fHistClustEMatch->Fill(clustEmatch);
      fNtrkl_ClustE->Fill(corr_nAcc,clustEmatch);

      Double_t clustLongE = clustMatch->GetM02();
      fHistClustLongAxis->Fill(clustLongE);

      Double_t EoverP = clustEmatch/track->P();

//---------------nsigma cut (electron)---------------
      if (fTPCnSigma > -1. && fTPCnSigma < 3.)
      {
        fPtEoverPE->Fill(track->Pt(),EoverP);
        if (TMath::Abs(pdg) == 11) fPtEoverPMCE->Fill(track->Pt(),EoverP);

        //to compare all tracks and electron M02
        if (EoverP > 0.7 && EoverP < 1.5)
        {
          if (TMath::Abs(pdg) == 11) fHistMCClsLAE->Fill(clustLongE);
          fHistClustLongAxisE->Fill(clustLongE);
        }

//---------------M02 cut (electron)---------------
        if (clustLongE > 0.1 && clustLongE < 0.6)
        {
          fPtEoverPEGeo->Fill(track->Pt(),EoverP);

//---------------EoverP cut (electron)---------------
          //////////////
          // ELECTRON //
          //////////////
          if (EoverP > 0.7 && EoverP < 1.5)
          {
            Float_t showerx[3];
            clustMatch->GetPosition(showerx);
            TVector3 Corepos(showerx[0],showerx[1],showerx[2]);
            Double_t Corephi = Corepos.Phi();
            Double_t Coreeta = Corepos.Eta();
            Double_t eleE = clustMatch->E();
            Double_t RsumE[3] = {0};
//---------------EMCal Cluster loop for Isolation cut start---------------
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
                if (R < 0.3) RsumE[0] = RsumE[0] + clustE;
                if (R < 0.4) RsumE[1] = RsumE[1] + clustE;
                if (R < 0.5) RsumE[2] = RsumE[2] + clustE;
              }
            }
            Double_t Eiso[3] = {0};
            for (Int_t isoR=0;isoR<3;isoR++) {
              Eiso[isoR] = (RsumE[isoR] - eleE)/eleE;
              fREisolation[isoR]->Fill(Eiso[isoR]);
            }
            //=====MC Data=====
            if (track->Pt() > 10.) {
              if (pidW == 1) fREiso_MCW->Fill(Eiso[0]);
              if (pidHF == 1) fREiso_MCHF->Fill(Eiso[0]);
            }
            if (track->Pt() > 20.) {
              if (pidW == 1) fREiso_MCWhpt->Fill(Eiso[0]);
              if (pidHF == 1) fREiso_MCHFhpt->Fill(Eiso[0]);
            }
            if (track->Pt() > 30.) {
              if (pidW == 1) fREiso_MCWvhpt->Fill(Eiso[0]);
              if (pidHF == 1) fREiso_MCHFvhpt->Fill(Eiso[0]);
            }
            //=================
            for (Int_t isoR=0;isoR<3;isoR++) {
              if (Eiso[isoR] >= 0.1) fHistPt_HFe[isoR]->Fill(track->Pt());
              if (Eiso[isoR] <= 0.05) fHistPt_We[isoR]->Fill(track->Pt());

//---------------1st Pt cut---------------
              if (track->Pt() > 10.)
              {
//---------------another track loop start---------------
                Int_t MaxPtTrackNum = 0;
                if (MaxPtTrackNum == i) MaxPtTrackNum = 1;
                for(Int_t j(0); j < iTracks; j++) {
                  AliAODTrack* anotrack = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
                  if(!anotrack || !anotrack->TestFilterBit(1)) continue;
                  if (j == i) continue;     //reject e<-W candidate

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
                  Double_t dPhi = anotrack->Phi() - track->Phi();

                  if (Eiso[isoR] >= 0.05) fdPhi_trkW_full[isoR]->Fill(dPhi,anotrack->Pt());
                  if (Eiso[isoR] <= 0.1) fdPhi_trkHF_full[isoR]->Fill(dPhi,anotrack->Pt());
                  if (track->Pt() > 30.) {
                    if (Eiso[isoR] >= 0.05) fdPhi_trkW_full_hpt[isoR]->Fill(dPhi,anotrack->Pt());
                    if (Eiso[isoR] <= 0.1) fdPhi_trkHF_full_hpt[isoR]->Fill(dPhi,anotrack->Pt());
                  }
                  //=== change dPhi range ===
                  if (dPhi < -1.*TMath::Pi()/3.) dPhi = dPhi + 2.*TMath::Pi();
                  if (dPhi > 5.*TMath::Pi()/3.) dPhi = dPhi - 2.*TMath::Pi();

                  if (Eiso[isoR] <= 0.05) {
                    fdPhi_trkW_Pt[isoR]->Fill(dPhi,anotrack->Pt());
                    fdPhi_trkW_ePt[isoR]->Fill(dPhi,track->Pt());
                    if (track->Pt() > 30.) fdPhi_trkW_Pt_hpt[isoR]->Fill(dPhi,anotrack->Pt());
                  }
                  if (Eiso[isoR] >= 0.1) {
                    fdPhi_trkHF_Pt[isoR]->Fill(dPhi,anotrack->Pt());
                    fdPhi_trkHF_ePt[isoR]->Fill(dPhi,track->Pt());
                    if (track->Pt() > 30.) fdPhi_trkHF_Pt_hpt[isoR]->Fill(dPhi,anotrack->Pt());
                  }

                  //===== Highest Pt track selection =====
                  if (dPhi >= 5.*TMath::Pi()/6. && dPhi <= 7.*TMath::Pi()/6.) {
                    AliAODTrack* Maxtrack = static_cast<AliAODTrack*>(fAOD->GetTrack(MaxPtTrackNum));
                    Double_t dPhiMaxTrk = (Maxtrack->Phi()) - (track->Phi());
                    //=== change dPhi range ===
                    if (dPhiMaxTrk < -1.*TMath::Pi()/3.) dPhiMaxTrk = dPhiMaxTrk + 2.*TMath::Pi();
                    if (dPhiMaxTrk > 5.*TMath::Pi()/3.) dPhiMaxTrk = dPhiMaxTrk - 2.*TMath::Pi();
                    //=== select back-to-back track ===
                    if (dPhiMaxTrk < 5.*TMath::Pi()/6. || dPhiMaxTrk > 7.*TMath::Pi()/6.) MaxPtTrackNum = j;
                    //=== higher Pt track ===
                    if (anotrack->Pt() >= Maxtrack->Pt()) MaxPtTrackNum = j;
                  }
                }   //other track loop end

                //===== Highest Pt Track / e<-W Pt Track =====
                AliAODTrack* MaxPtTrk = static_cast<AliAODTrack*>(fAOD->GetTrack(MaxPtTrackNum));

//---------------e<-W isolation cut---------------
                //////////////////////
                // e <- W candidate //
                //////////////////////
                if (Eiso[isoR] <= 0.05) {
                  fPt_maxtrack_W[isoR]->Fill(MaxPtTrk->Pt()/track->Pt());
                  fNtrkl_PtOfMaxTrk_W[isoR]->Fill(corr_nAcc,MaxPtTrk->Pt()/track->Pt());
                  fHistPt_We_Ntrkl[isoR]->Fill(track->Pt(),corr_nAcc);
//---------------2nd Pt cut----------
                  if (track->Pt() > 30.) {
                    fPt_maxtrack_W_hpt[isoR]->Fill(MaxPtTrk->Pt()/track->Pt());
                    fNtrkl_PtOfMaxTrk_W_hpt[isoR]->Fill(corr_nAcc,MaxPtTrk->Pt()/track->Pt());
                  }
                }       //isolation cut end
              }         //standard Pt cut end
            }           //varied R of Eiso loop end
          }             //M02 cut
        }               //EoverP cut
      }                 //nsigma cut

//---------------nsigma cut (hadron)---------------
      if (fTPCnSigma < -3.)
      {
        fPtEoverPH->Fill(track->Pt(),EoverP);
//---------------M02 cut---------------
        if (clustLongE > 0.1 && clustLongE < 0.6) fPtEoverPHGeo->Fill(track->Pt(),EoverP);
      }		//hadron filter end
    }		//track match end
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
Double_t AliAnalysisTaskWHMult::GetCorrectedNtrackletsD(TTree* tree, Double_t uncorrectedNacc, Int_t BinOfvertexZ, AliAODEvent* fAOD)
{
  //if(BinOfvertexZ<1 || BinOfvertexZ > 100) return uncorrectedNacc;
  if(!tree) return uncorrectedNacc;

  Int_t TreeEntry = BinOfvertexZ - 1;
  Double_t corr_ZvNtl = 0.;
  Double_t corr_ZvNtl_data_a = 0.;
  Double_t corr_ZvNtl_data_b = 0.;
  Double_t corr_ZvNtl_data_c = 0.;
  Double_t corr_ZvNtl_MC_a = 0.;
  Double_t corr_ZvNtl_MC_b = 0.;
  Double_t corr_ZvNtl_MC_c = 0.;

  Int_t runNo = fAOD->GetRunNumber();

  if (!fMCarray) {
    //===== data 2016 =====
    if (runNo >= 256941 && runNo <= 264347) tree->SetBranchAddress("vertexZ_Ntrkl_data_a",&corr_ZvNtl_data_a);  //2016
    //===== data 2017 =====
    if (runNo >= 270581 && runNo <= 282704) tree->SetBranchAddress("vertexZ_Ntrkl_data_b",&corr_ZvNtl_data_b);  //2017
    //===== data 2018 =====
    if (runNo >= 285008 && runNo <= 294925) tree->SetBranchAddress("vertexZ_Ntrkl_data_c",&corr_ZvNtl_data_c);  //2018
  }
  if (fMCarray) {
    //===== MC 2016 =====
    if (runNo >= 256941 && runNo <= 264347) tree->SetBranchAddress("vertexZ_Ntrkl_MC_a",&corr_ZvNtl_MC_a);      //2016
    //===== MC 2017 =====
    if (runNo >= 270581 && runNo <= 282704) tree->SetBranchAddress("vertexZ_Ntrkl_MC_b",&corr_ZvNtl_MC_b);      //2017
    //===== MC 2018 =====
    if (runNo >= 285008 && runNo <= 294925) tree->SetBranchAddress("vertexZ_Ntrkl_MC_c",&corr_ZvNtl_MC_c);      //2018
  }

  tree->GetEntry(TreeEntry);

  if (!fMCarray) {
    //===== data 2016 =====
    if (runNo >= 256941 && runNo <= 264347) corr_ZvNtl = corr_ZvNtl_data_a;
    //===== data 2017 =====
    if (runNo >= 270581 && runNo <= 282704) corr_ZvNtl = corr_ZvNtl_data_b;
    //===== data 2018 =====
    if (runNo >= 285008 && runNo <= 294925) corr_ZvNtl = corr_ZvNtl_data_c;
  }
  if (fMCarray) {
    //===== MC 2016 =====
    if (runNo >= 256941 && runNo <= 264347) corr_ZvNtl = corr_ZvNtl_MC_a;
    //===== MC 2017 =====
    if (runNo >= 270581 && runNo <= 282704) corr_ZvNtl = corr_ZvNtl_MC_b;
    //===== MC 2018 =====
    if (runNo >= 285008 && runNo <= 294925) corr_ZvNtl = corr_ZvNtl_MC_c;
  }

  Double_t deltaM = 0;
  deltaM = uncorrectedNacc*(corr_ZvNtl - 1.);

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
