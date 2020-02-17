/**************************************************************************
 * Copyright(c) 2007-2019, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Analysis task to produce Ds candidates mass spectra           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
// Mantainers: F. Catalano, fabio.catalano@cern.ch               //
//             F. Grosa, fabrizio.grosa@cern.ch                  //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliHFMLResponseDstoKKpi.h"

#include "AliAnalysisTaskSEDs.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEDs);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs() : AliAnalysisTaskSE()                                   
{
  /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDs::AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi *analysiscuts, Bool_t createMLtree) : AliAnalysisTaskSE(name),
                                                                                                                     fCreateMLtree(createMLtree)
{                                                   
  /// Standard constructor
  SetAnalysisCuts(analysiscuts);
  Int_t nptbins = fAnalysisCuts->GetNPtBins();
  Float_t *ptlim = fAnalysisCuts->GetPtBinLimits();
  SetPtBins(nptbins, ptlim);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, AliNormalizationCounter::Class());
  if(fCreateMLtree)
    DefineOutput(4, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEDs::~AliAnalysisTaskSEDs()
{
  // Destructor
  if (fOutput && !fOutput->IsOwner())
  {
    delete fHistNEvents;
    delete fHistAllV0multNTPCout;
    delete fHistSelV0multNTPCout;
    delete fCosPHist3D;
    delete fCosPxyHist3D;
    delete fDLenHist3D;
    delete fDLenxyHist3D;
    delete fNDLenxyHist3D;
    delete fSigVertHist3D;
    delete fDCAHist3D;
    delete fNormIPHist3D;
    delete fCosPiDsHist3D;
    delete fCosPiKPhiHist3D;
    delete fPtProng0Hist3D;
    delete fPtProng1Hist3D;
    delete fPtProng2Hist3D;

    for (Int_t iHist = 0; iHist < 4; iHist++)
    {
      delete fChanHist[iHist];
    }
    for (Int_t iHist = 0; iHist < 4 * fNPtBins; iHist++)
    {
      delete fMassHist[iHist];
      delete fMassHistPhi[iHist];
      delete fMassHistK0st[iHist];
      delete fCosPHist[iHist];
      delete fDLenHist[iHist];
      delete fSumd02Hist[iHist];
      delete fSigVertHist[iHist];
      delete fPtMaxHist[iHist];
      delete fPtCandHist[iHist];
      delete fDCAHist[iHist];
      delete fPtProng0Hist[iHist];
      delete fPtProng1Hist[iHist];
      delete fPtProng2Hist[iHist];
      delete fDalitz[iHist];
      delete fDalitzPhi[iHist];
      delete fDalitzK0st[iHist];
    }
    for (Int_t iPt = 0; iPt < fNPtBins; iPt++)
    {
      delete fMassHistKK[iPt];
      delete fMassHistKpi[iPt];
      delete fMassHistKKVsKKpi[iPt];
      delete fMassHistKpiVsKKpi[iPt];
      delete fMassRotBkgHistPhi[iPt];
      delete fMassLSBkgHistPhi[iPt];
      delete fMassRSBkgHistPhi[iPt];
    }
    delete fPtVsMass;
    delete fPtVsMassPhi;
    delete fPtVsMassK0st;
    delete fYVsPt;
    delete fYVsPtSig;
    for (Int_t iHist = 0; iHist < 3; iHist++)
    {
      delete fHistCentrality[iHist];
      delete fHistCentralityMult[iHist];
    }

    delete fnSparse;
    if (fFillImpParSparse)
    {
      delete fImpParSparse;
      for (Int_t iHist = 0; iHist < 4; iHist++)
        delete fImpParSparseMC[iHist];
    }
    for (Int_t iHist = 0; iHist < 5; iHist++)
    {
      delete fnSparseMC[iHist];
    }
    if (fFillSparseDplus)
    {
      for (Int_t iHist = 0; iHist < 4; iHist++)
          delete fnSparseMCDplus[iHist];
    }

    if(fApplyML && fEnablePIDMLSparses)
    {
      for (Int_t iHist=0; iHist<2; iHist++)
      {
          delete fnSparseNsigmaPIDVsML[iHist];
      }
    }
  }
  delete fOutput;
  delete fListCuts;
  delete fCounter;
  delete fAnalysisCuts;

  if(fApplyML && fMLResponse)
    delete fMLResponse;

  if(fCreateMLtree && fMLhandler)
    delete fMLhandler;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::Init()
{
  /// Initialization

  if (fDebug > 1)
    printf("AnalysisTaskSEDs::Init() \n");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("CutObjects");

  AliRDHFCutsDstoKKpi *analysis = new AliRDHFCutsDstoKKpi(*fAnalysisCuts);
  analysis->SetName("AnalysisCuts");

  fListCuts->Add(analysis);
  PostData(2, fListCuts);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserCreateOutputObjects()
{
  /// Create the output container
  //
  if (fDebug > 1)
    printf("AnalysisTaskSEDs::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "number of events ", 16, -0.5, 15.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
  fHistNEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents Mismatched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(4, "nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(5, "n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(6, "n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(7, "n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(8, "n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(9, "n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(10, "n. rejected for pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(11, "no. of out centrality events");
  fHistNEvents->GetXaxis()->SetBinLabel(12, "no. of 3 prong candidates");
  fHistNEvents->GetXaxis()->SetBinLabel(13, "no. of Ds after filtering cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(14, "no. of Ds after selection cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(15, "no. of not on-the-fly rec Ds");
  fHistNEvents->GetXaxis()->SetBinLabel(16, "no. of Ds rejected by preselect");

  fHistNEvents->GetXaxis()->SetNdivisions(1, kFALSE);

  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistCentrality[0] = new TH1F("hCentr", "centrality", 10000, 0., 100.);
  fHistCentrality[1] = new TH1F("hCentr(selectedCent)", "centrality(selectedCent)", 10000, 0., 100.);
  fHistCentrality[2] = new TH1F("hCentr(OutofCent)", "centrality(OutofCent)", 10000, 0., 100.);
  fHistCentralityMult[0] = new TH2F("hCentrMult", "centrality vs mult", 100, 0.5, 30000.5, 40, 0., 100.);
  fHistCentralityMult[1] = new TH2F("hCentrMult(selectedCent)", "centrality vs mult(selectedCent)", 100, 0.5, 30000.5, 40, 0., 100.);
  fHistCentralityMult[2] = new TH2F("hCentrMult(OutofCent)", "centrality vs mult(OutofCent)", 100, 0.5, 30000.5, 40, 0., 100.);
  for (Int_t iHist = 0; iHist < 3; iHist++)
  {
    fOutput->Add(fHistCentrality[iHist]);
    fOutput->Add(fHistCentralityMult[iHist]);
  }
  if (fDoCutV0multTPCout)
  {
    fHistAllV0multNTPCout = new TH2F("HistAllV0multNTPCout", "V0mult vs # TPCout (all) ;V0mult ;# TPCout", 1000, 0., 40000, 1000, 0, 30000);
    fHistSelV0multNTPCout = new TH2F("HistSelV0multNTPCout", "V0mult vs # TPCout (sel) ;V0mult ;# TPCout", 1000, 0., 40000, 1000, 0, 30000);
    fOutput->Add(fHistAllV0multNTPCout);
    fOutput->Add(fHistSelV0multNTPCout);
  }

  Double_t massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();

  Int_t nInvMassBins = (Int_t)(fMassRange / fMassBinSize + 0.5);
  if (nInvMassBins % 2 == 1)
    nInvMassBins++;
  Double_t minMass = massDs - 0.5 * nInvMassBins * fMassBinSize;
  Double_t maxMass = massDs + 0.5 * nInvMassBins * fMassBinSize;
  fminMass = minMass;
  fmaxMass = maxMass;

  TString hisname;
  TString htype;
  Int_t index;
  for (Int_t iType = 0; iType < 4; iType++)
  {
    for (Int_t iPt = 0; iPt < fNPtBins; iPt++)
    {
      if (iType == 0)
      {
        htype = "All";
        index = GetHistoIndex(iPt);
      }
      else if (iType == 1)
      {
        htype = "Sig";
        index = GetSignalHistoIndex(iPt);
      }
      else if (iType == 2)
      {
        htype = "Bkg";
        index = GetBackgroundHistoIndex(iPt);
      }
      else
      {
        htype = "ReflSig";
        index = GetReflSignalHistoIndex(iPt);
      }
      hisname.Form("hMass%sPt%d", htype.Data(), iPt);
      fMassHist[index] = new TH1F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass);
      fMassHist[index]->Sumw2();
      hisname.Form("hMass%sPt%dphi", htype.Data(), iPt);
      fMassHistPhi[index] = new TH1F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass);
      fMassHistPhi[index]->Sumw2();
      hisname.Form("hMass%sPt%dk0st", htype.Data(), iPt);
      fMassHistK0st[index] = new TH1F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass);
      fMassHistK0st[index]->Sumw2();
      hisname.Form("hCosP%sPt%d", htype.Data(), iPt);
      fCosPHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.5, 1.);
      hisname.Form("hCosPxy%sPt%d", htype.Data(), iPt);
      fCosPxyHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.5, 1.);
      hisname.Form("hDLen%sPt%d", htype.Data(), iPt);
      fDLenHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 0.5);
      hisname.Form("hDLenxy%sPt%d", htype.Data(), iPt);
      fDLenxyHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 0.5);
      hisname.Form("hNDLenxy%sPt%d", htype.Data(), iPt);
      fNDLenxyHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 11.);
      hisname.Form("hSumd02%sPt%d", htype.Data(), iPt);
      fSumd02Hist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 1.);
      hisname.Form("hSigVert%sPt%d", htype.Data(), iPt);
      fSigVertHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 0.1);
      hisname.Form("hPtMax%sPt%d", htype.Data(), iPt);
      fPtMaxHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.5, 20.);
      hisname.Form("hPtCand%sPt%d", htype.Data(), iPt);
      fPtCandHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.5, 20.);
      hisname.Form("hDCA%sPt%d", htype.Data(), iPt);
      fDCAHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 0.1);
      hisname.Form("hNormIP%sPt%d", htype.Data(), iPt);
      fNormIPHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 6.);
      hisname.Form("hCosPiDs%sPt%d", htype.Data(), iPt);
      fCosPiDsHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.5, 1.);
      hisname.Form("hCosPiKPhi%sPt%d", htype.Data(), iPt);
      fCosPiKPhiHist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0., 0.5);
      hisname.Form("hPtProng0%sPt%d", htype.Data(), iPt);
      fPtProng0Hist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.0, 20.);
      hisname.Form("hPtProng1%sPt%d", htype.Data(), iPt);
      fPtProng1Hist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.0, 20.);
      hisname.Form("hPtProng2%sPt%d", htype.Data(), iPt);
      fPtProng2Hist[index] = new TH1F(hisname.Data(), hisname.Data(), 100, 0.0, 20.);
      hisname.Form("hDalitz%sPt%d", htype.Data(), iPt);
      fDalitz[index] = new TH2F(hisname.Data(), hisname.Data(), 100, 0., 2., 100, 0., 2.);
      hisname.Form("hDalitz%sPt%dphi", htype.Data(), iPt);
      fDalitzPhi[index] = new TH2F(hisname.Data(), hisname.Data(), 100, 0., 2., 100, 0., 2.);
      hisname.Form("hDalitz%sPt%dk0st", htype.Data(), iPt);
      fDalitzK0st[index] = new TH2F(hisname.Data(), hisname.Data(), 100, 0., 2., 100, 0., 2.);
    }
  }

  for (Int_t iHist = 0; iHist < 4 * fNPtBins; iHist++)
  {
    fOutput->Add(fMassHist[iHist]);
    fOutput->Add(fMassHistPhi[iHist]);
    fOutput->Add(fMassHistK0st[iHist]);
    fOutput->Add(fPtCandHist[iHist]);
    if (fDoCutVarHistos)
    {
      fOutput->Add(fCosPHist[iHist]);
      fOutput->Add(fCosPxyHist[iHist]);
      fOutput->Add(fDLenHist[iHist]);
      fOutput->Add(fDLenxyHist[iHist]);
      fOutput->Add(fNDLenxyHist[iHist]);
      fOutput->Add(fSumd02Hist[iHist]);
      fOutput->Add(fSigVertHist[iHist]);
      fOutput->Add(fPtMaxHist[iHist]);
      fOutput->Add(fDCAHist[iHist]);
      fOutput->Add(fNormIPHist[iHist]);
      fOutput->Add(fCosPiDsHist[iHist]);
      fOutput->Add(fCosPiKPhiHist[iHist]);
      fOutput->Add(fPtProng0Hist[iHist]);
      fOutput->Add(fPtProng1Hist[iHist]);
      fOutput->Add(fPtProng2Hist[iHist]);
      fOutput->Add(fDalitz[iHist]);
      fOutput->Add(fDalitzPhi[iHist]);
      fOutput->Add(fDalitzK0st[iHist]);
    }
  }

  fChanHist[0] = new TH1F("hChanAll", "KKpi and piKK candidates", 64, -0.5, 63.5);
  fChanHist[1] = new TH1F("hChanSig", "KKpi and piKK candidates", 64, -0.5, 63.5);
  fChanHist[2] = new TH1F("hChanBkg", "KKpi and piKK candidates", 64, -0.5, 63.5);
  fChanHist[3] = new TH1F("hChanReflSig", "KKpi and piKK candidates", 64, -0.5, 63.5);
  for (Int_t iHist = 0; iHist < 4; iHist++)
  {
    fChanHist[iHist]->SetMinimum(0);
    fOutput->Add(fChanHist[iHist]);
  }

  fCosPHist3D = new TH3F("fCosPHist3D", "CosP vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0.5, 1.);
  fCosPxyHist3D = new TH3F("fCosPxyHist3D", "CosPxy vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0.5, 1.);
  fDLenHist3D = new TH3F("fDLenHist3D", "DLen vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 0.5);
  fDLenxyHist3D = new TH3F("fDLenxyHist3D", "DLenxy vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 0.5);
  fNDLenxyHist3D = new TH3F("fNDLenxyHist3D", "NDLenxy vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 11.);
  fSigVertHist3D = new TH3F("fSigVertHist3D", "SigVert vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 0.1);
  fDCAHist3D = new TH3F("fDCAHist3D", "DCA vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 0.1);
  fNormIPHist3D = new TH3F("fNormIPHist3D", "nIP vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 6.);
  fCosPiDsHist3D = new TH3F("fCosPiDsHist3D", "CosPiDs vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0.5, 1.);
  fCosPiKPhiHist3D = new TH3F("fCosPiKPhiHist3D", "CosPiKPhi vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0., 0.5);
  fPtProng0Hist3D = new TH3F("fPtProng0Hist3D", "Pt prong0 vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0.0, 20.);
  fPtProng1Hist3D = new TH3F("fPtProng1Hist3D", "Pt prong1 vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0.0, 20.);
  fPtProng2Hist3D = new TH3F("fPtProng2Hist3D", "Pt prong2 vs Ds mass", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins], 0., fPtLimits[fNPtBins], 100, 0.0, 20.);

  if (!fReadMC && fDoCutVarHistos)
  {
    fOutput->Add(fCosPHist3D);
    fOutput->Add(fCosPxyHist3D);
    fOutput->Add(fDLenHist3D);
    fOutput->Add(fDLenxyHist3D);
    fOutput->Add(fNDLenxyHist3D);
    fOutput->Add(fSigVertHist3D);
    fOutput->Add(fDCAHist3D);
    fOutput->Add(fNormIPHist3D);
    fOutput->Add(fCosPiDsHist3D);
    fOutput->Add(fCosPiKPhiHist3D);
    fOutput->Add(fPtProng0Hist3D);
    fOutput->Add(fPtProng1Hist3D);
    fOutput->Add(fPtProng2Hist3D);
  }

  fPtVsMass = new TH2F("hPtVsMass", "PtVsMass (prod. cuts)", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins] * 2, 0., fPtLimits[fNPtBins]);
  fPtVsMassPhi = new TH2F("hPtVsMassPhi", "PtVsMass (phi selection)", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins] * 10, 0., fPtLimits[fNPtBins]);
  fPtVsMassK0st = new TH2F("hPtVsMassK0st", "PtVsMass (K0* selection)", nInvMassBins, minMass, maxMass, (Int_t)fPtLimits[fNPtBins] * 10, 0., fPtLimits[fNPtBins]);
  fYVsPt = new TH2F("hYVsPt", "YvsPt (prod. cuts)", (Int_t)fPtLimits[fNPtBins] * 2, 0., fPtLimits[fNPtBins], 80, -2., 2.);
  fYVsPtSig = new TH2F("hYVsPtSig", "YvsPt (MC, only sig., prod. cuts)", (Int_t)fPtLimits[fNPtBins] * 2, 0., fPtLimits[fNPtBins], 80, -2., 2.);

  for (Int_t iPt = 0; iPt < fNPtBins; iPt++)
  {
    hisname.Form("hMassKKPt%d", iPt);
    fMassHistKK[iPt] = new TH1F(hisname.Data(), hisname.Data(), 200, 0.95, 1.35);
    fMassHistKK[iPt]->Sumw2();
    fOutput->Add(fMassHistKK[iPt]);
    hisname.Form("hMassKpiPt%d", iPt);
    fMassHistKpi[iPt] = new TH1F(hisname.Data(), hisname.Data(), 200, 0.7, 1.1);
    fMassHistKpi[iPt]->Sumw2();
    fOutput->Add(fMassHistKpi[iPt]);
    hisname.Form("hMassKKVsKKpiPt%d", iPt);
    fMassHistKKVsKKpi[iPt] = new TH2F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass, 200, 0.95, 1.15);
    fOutput->Add(fMassHistKKVsKKpi[iPt]);
    hisname.Form("hMassKpiVsKKpiPt%d", iPt);
    fMassHistKpiVsKKpi[iPt] = new TH2F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass, 200, 0.7, 1.1);
    fOutput->Add(fMassHistKpiVsKKpi[iPt]);
    if (fDoRotBkg)
    {
      hisname.Form("hMassAllPt%dphi_RotBkg", iPt);
      fMassRotBkgHistPhi[iPt] = new TH1F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass);
      fMassRotBkgHistPhi[iPt]->Sumw2();
      fOutput->Add(fMassRotBkgHistPhi[iPt]);
    }
    if (fDoBkgPhiSB)
    {
      hisname.Form("fMassLSBkgHistPhiPt%d", iPt);
      fMassLSBkgHistPhi[iPt] = new TH1F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass);
      fMassLSBkgHistPhi[iPt]->Sumw2();
      fOutput->Add(fMassLSBkgHistPhi[iPt]);
      hisname.Form("fMassRSBkgHistPhiPt%d", iPt);
      fMassRSBkgHistPhi[iPt] = new TH1F(hisname.Data(), hisname.Data(), nInvMassBins, minMass, maxMass);
      fMassRSBkgHistPhi[iPt]->Sumw2();
      fOutput->Add(fMassRSBkgHistPhi[iPt]);
    }
  }

  fOutput->Add(fPtVsMass);
  fOutput->Add(fPtVsMassPhi);
  fOutput->Add(fPtVsMassK0st);
  fOutput->Add(fYVsPt);
  fOutput->Add(fYVsPtSig);

  //Sparses for Cut variation
  if (fFillSparse)
    CreateCutVarsAndEffSparses();

  //Sparses for Impact parameter fits
  if (fFillImpParSparse)
    CreateImpactParameterSparses();

  //Counter for Normalization
  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();
  PostData(3, fCounter);

  //Loading of ML models
  if(fApplyML) {
    fMLResponse = new AliHFMLResponseDstoKKpi("DstoKKpiMLResponse", "DstoKKpiMLResponse", fConfigPath.Data());
    fMLResponse->MLResponseInit();

    if(fEnablePIDMLSparses)
      CreatePIDMLSparses();
  }

  //Create ML tree
  if(fCreateMLtree) {
    OpenFile(4);
    fMLhandler = new AliHFMLVarHandlerDstoKKpi(fPIDopt, AliHFMLVarHandlerDstoKKpi::kDeltaMassKKPhi);
    fMLhandler->SetAddSingleTrackVars(fAddSingleTrackVar);
    if(fReadMC && fFillOnlySignal)
      fMLhandler->SetFillOnlySignal();
    fMLtree = fMLhandler->BuildTree("treeMLDs", "treeMLDs");
    fMLtree->SetMaxVirtualSize(1.e+8);
    PostData(4, fMLtree);
  }

  //Set seed of gRandom
  if(fCreateMLtree && fEnableEvtSampling)
    gRandom->SetSeed(fSeedSampling);
  
  PostData(1, fOutput);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::UserExec(Option_t * /*option*/)
{
  /// Ds selection for current event, fill mass histos and selection variable histo
  /// separate signal and backgound if fReadMC is activated

  if(fCreateMLtree && fEnableEvtSampling && gRandom->Rndm() > fFracEvtToKeep)
    return;
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent *>(InputEvent());

  fHistNEvents->Fill(0); // all events
  if (fAODProtection >= 0)
  {
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
    {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(2);
      return;
    }
    fHistNEvents->Fill(1);
  }

  TClonesArray *array3Prong = 0;
  if (!aod && AODEvent() && IsStandardAOD())
  {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent *>(AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler *aodHandler = (AliAODHandler *)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if (aodHandler->GetExtensions())
    {
      AliAODExtension *ext = (AliAODExtension *)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong = (TClonesArray *)aodFromExt->GetList()->FindObject("Charm3Prong");
    }
  }
  else if (aod)
  {
    array3Prong = (TClonesArray *)aod->GetList()->FindObject("Charm3Prong");
  }

  if (!aod || !array3Prong)
  {
    printf("AliAnalysisTaskSEDs::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField()) < 0.001) return;

  fHistNEvents->Fill(3); // count event
  // Post the data already here
  PostData(1, fOutput);

  fCounter->StoreEvent(aod, fAnalysisCuts, fReadMC);

  Bool_t isEvSel = fAnalysisCuts->IsEventSelected(aod);
  Float_t ntracks = aod->GetNumberOfTracks();
  Float_t evCentr = fAnalysisCuts->GetCentrality(aod);

  fHistCentrality[0]->Fill(evCentr);
  fHistCentralityMult[0]->Fill(ntracks, evCentr);
  if (fAnalysisCuts->IsEventRejectedDueToTrigger())
    fHistNEvents->Fill(5);
  if (fAnalysisCuts->IsEventRejectedDueToNotRecoVertex())
    fHistNEvents->Fill(6);
  if (fAnalysisCuts->IsEventRejectedDueToVertexContributors())
    fHistNEvents->Fill(7);
  if (fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion())
    fHistNEvents->Fill(8);
  if (fAnalysisCuts->IsEventRejectedDueToPileup())
    fHistNEvents->Fill(9);
  if (fAnalysisCuts->IsEventRejectedDueToCentrality())
  {
    fHistNEvents->Fill(10);
    fHistCentrality[2]->Fill(evCentr);
    fHistCentralityMult[2]->Fill(ntracks, evCentr);
  }
  
  // no physics selection can be applied for upgrade studies
  if(fSystem == kUpgr && fAnalysisCuts->IsEventRejectedDuePhysicsSelection())
    isEvSel = kTRUE;
  
  Int_t runNumber = aod->GetRunNumber();

  TClonesArray *arrayMC = 0;
  AliAODMCHeader *mcHeader = 0;

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex *)aod->GetPrimaryVertex();
  //    vtx1->Print();

  // load MC particles
  if (fReadMC)
  {

    arrayMC = (TClonesArray *)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!arrayMC)
    {
      printf("AliAnalysisTaskSEDs::UserExec: MC particles branch not found!\n");
      return;
    }

    // load MC header
    mcHeader = (AliAODMCHeader *)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader)
    {
      printf("AliAnalysisTaskSEDs::UserExec: MC header branch not found!\n");
      return;
    }
  }

  if (fReadMC && fFillSparse)
  {
    if (aod->GetTriggerMask() == 0 && (runNumber >= 195344 && runNumber <= 195677))
      // protection for events with empty trigger mask in p-Pb
      return;
    if (fAnalysisCuts->GetUseCentrality() > 0 && fAnalysisCuts->IsEventSelectedInCentrality(aod) != 0)
      // events not passing the centrality selection can be removed immediately.
      return;
    FillMCGenAccHistos(arrayMC, mcHeader);
  }

  if (!isEvSel) return;

  fHistNEvents->Fill(4);
  fHistCentrality[1]->Fill(evCentr);
  fHistCentralityMult[1]->Fill(ntracks, evCentr);

  if (fDoCutV0multTPCout)
  {
    //cut on V0mult. vs #tracks kTPCout
    Float_t V0mult = 0.;
    AliAODVZERO *aodVZERO = (AliAODVZERO *)aod->GetVZEROData();
    if (aodVZERO)
    {
      for (int ich = 0; ich < 64; ich++)
        V0mult += aodVZERO->GetMultiplicity(ich);
    }
    Int_t nTPCout = 0;
    for (Int_t i = 0; i < aod->GetNumberOfTracks(); ++i)
    {
      AliVTrack *track = aod->GetTrack(i);
      if (!track)
        continue;
      if ((track->GetStatus() & AliVTrack::kTPCout))
        nTPCout++;
    }
    fHistAllV0multNTPCout->Fill(V0mult, nTPCout);
    if (nTPCout > (0.32 * V0mult + 750))
      return;
    else
      fHistSelV0multNTPCout->Fill(V0mult, nTPCout);
  }

  Int_t n3Prong = array3Prong->GetEntriesFast();
  if (fDebug > 1)
    printf("Number of Ds->KKpi: %d\n", n3Prong);

  Int_t pdgDstoKKpi[3] = {321, 321, 211};
  Int_t nSelected = 0;
  Int_t nFiltered = 0;
  Double_t massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();

  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF vHF = AliAnalysisVertexingHF();

  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++)
  {

    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong *)array3Prong->UncheckedAt(i3Prong);
    fHistNEvents->Fill(11);

    if (fUseSelectionBit && !(d->HasSelectionBit(AliRDHFCuts::kDsCuts))) continue;
    nFiltered++;
    fHistNEvents->Fill(12);

    TObjArray arrTracks(3);
    for(Int_t ipr=0;ipr<3;ipr++){
      AliAODTrack *tr=vHF.GetProng(aod,d,ipr);
      arrTracks.AddAt(tr,ipr);
    }
    if(!fAnalysisCuts->PreSelect(arrTracks)){
      fHistNEvents->Fill(15);
      continue;
    }

    if (!(vHF.FillRecoCand(aod, d)))
    {                         ////Fill the data members of the candidate only if they are empty.
      fHistNEvents->Fill(14); //monitor how often this fails
      continue;
    }

    Bool_t unsetvtx = kFALSE;
    if (!d->GetOwnPrimaryVtx())
    {
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx = kTRUE;
      // NOTE: the ow primary vertex should be unset, otherwise there is a memory leak
      // Pay attention if you use continue inside this loop!!!
    }

    Bool_t recVtx = kFALSE;
    AliAODVertex *origownvtx = nullptr;

    Double_t ptCand = d->Pt();
    Int_t iPtBin = TMath::BinarySearch(fNPtBins, fPtLimits, (Float_t)ptCand);
    Double_t rapid = d->YDs();
    fYVsPt->Fill(ptCand, rapid);
    Bool_t isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(ptCand, rapid);
    Double_t massKK_KKpi = 0.;
    Double_t massKK_piKK = 0.;
    Double_t massKp = 0;
    Double_t masspK = 0;
    Double_t invMass_KKpi = 0.;
    Double_t invMass_piKK = 0.;

    if (fCreateMLtree && fEnableCandSampling) // apply sampling in pt
    { 
      Double_t pseudoRand = ptCand * 1000. - (long)(ptCand * 1000);
      if (pseudoRand > fFracCandToKeep && ptCand < fMaxCandPtSampling)
      {
        if (unsetvtx)
          d->UnsetOwnPrimaryVtx();
        continue;
      }
    }

    if (!isFidAcc) // check if candidate is in acceptance
    {
      if (unsetvtx)
        d->UnsetOwnPrimaryVtx();
      continue;
    }

    Int_t retCodeAnalysisCuts = fAnalysisCuts->IsSelected(d, AliRDHFCuts::kAll, aod);
    Int_t retCodeNoRes = retCodeAnalysisCuts;
    Bool_t origRes = fAnalysisCuts->IsCutOnResonancesApplied();
    if (origRes)
    {
      fAnalysisCuts->ApplyCutOnResonances(kFALSE);
      retCodeNoRes = fAnalysisCuts->IsSelected(d, AliRDHFCuts::kAll, aod);
      fAnalysisCuts->ApplyCutOnResonances(origRes);
    }

    if (retCodeNoRes & 1)
    { //KKpi
      massKK_KKpi = d->InvMass2Prongs(0, 1, 321, 321);
      massKp = d->InvMass2Prongs(1, 2, 321, 211);
      invMass_KKpi = d->InvMassDsKKpi();
      fMassHistKK[iPtBin]->Fill(massKK_KKpi);
      fMassHistKpi[iPtBin]->Fill(massKp);
      fMassHistKKVsKKpi[iPtBin]->Fill(invMass_KKpi, massKK_KKpi);
      fMassHistKpiVsKKpi[iPtBin]->Fill(invMass_KKpi, massKp);
    }
    if (retCodeNoRes & 2)
    { //piKK
      massKK_piKK = d->InvMass2Prongs(1, 2, 321, 321);
      masspK = d->InvMass2Prongs(0, 1, 211, 321);
      invMass_piKK = d->InvMassDspiKK();
      fMassHistKK[iPtBin]->Fill(massKK_piKK);
      fMassHistKpi[iPtBin]->Fill(masspK);
      fMassHistKKVsKKpi[iPtBin]->Fill(invMass_piKK, massKK_piKK);
      fMassHistKpiVsKKpi[iPtBin]->Fill(invMass_piKK, masspK);
    }

    Int_t isKKpi = retCodeAnalysisCuts & 1;
    Int_t ispiKK = retCodeAnalysisCuts & 2;
    Int_t isPhiKKpi = retCodeAnalysisCuts & 4;
    Int_t isPhipiKK = retCodeAnalysisCuts & 8;
    Int_t isK0starKKpi = retCodeAnalysisCuts & 16;
    Int_t isK0starpiKK = retCodeAnalysisCuts & 32;

    if (retCodeAnalysisCuts <= 0)
    {
      if (unsetvtx)
        d->UnsetOwnPrimaryVtx();
      continue;
    }

    if (fAnalysisCuts->GetIsPrimaryWithoutDaughters())
    {
      if (d->GetOwnPrimaryVtx())
        origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
      if (fAnalysisCuts->RecalcOwnPrimaryVtx(d, aod))
        recVtx = kTRUE;
      else
        fAnalysisCuts->CleanOwnPrimaryVtx(d, aod, origownvtx);
    }

    fHistNEvents->Fill(13);
    nSelected++;

    Int_t index = GetHistoIndex(iPtBin);
    fPtCandHist[index]->Fill(ptCand);

    Double_t weightKKpi = 1.;
    Double_t weightpiKK = 1.;
    if (fAnalysisCuts->GetPidOption() == AliRDHFCutsDstoKKpi::kBayesianWeights)
    {
      weightKKpi = fAnalysisCuts->GetWeightForKKpi();
      weightpiKK = fAnalysisCuts->GetWeightForpiKK();
      if (weightKKpi > 1. || weightKKpi < 0.)
        weightKKpi = 0.;
      if (weightpiKK > 1. || weightpiKK < 0.)
        weightpiKK = 0.;
    }

    fChanHist[0]->Fill(retCodeAnalysisCuts);

    const Int_t nProng = 3;
    Int_t indexMCKKpi = -1;
    Int_t indexMCpiKK = -1;
    Int_t labDs = -1;
    Int_t labDplus = -1;
    Int_t pdgCode0 = -999;

    AliAODMCParticle *partDs = nullptr;
    Int_t orig = 0;
    Int_t orig_HijingCheck = 0;
    Bool_t isCandInjected = kFALSE;
    Float_t trueImpParDsFromB = 99999.;

    if (fReadMC)
    {
      labDs = d->MatchToMC(431, arrayMC, nProng, pdgDstoKKpi);
      if (labDs >= 0)
      {
        partDs = (AliAODMCParticle*)arrayMC->At(labDs);
        Int_t labDau0 = ((AliAODTrack *)d->GetDaughter(0))->GetLabel();
        AliAODMCParticle *p = (AliAODMCParticle *)arrayMC->UncheckedAt(TMath::Abs(labDau0));
        pdgCode0 = TMath::Abs(p->GetPdgCode());

        if (isKKpi)
        {
          if (pdgCode0 == 321)
          {
            indexMCKKpi = GetSignalHistoIndex(iPtBin);
            fYVsPtSig->Fill(ptCand, rapid);
            fChanHist[1]->Fill(retCodeAnalysisCuts);
          }
          else
          {
            indexMCKKpi = GetReflSignalHistoIndex(iPtBin);
            fChanHist[3]->Fill(retCodeAnalysisCuts);
          }
        }
        if (ispiKK)
        {
          if (pdgCode0 == 211)
          {
            indexMCpiKK = GetSignalHistoIndex(iPtBin);
            fYVsPtSig->Fill(ptCand, rapid);
            fChanHist[1]->Fill(retCodeAnalysisCuts);
          }
          else
          {
            indexMCpiKK = GetReflSignalHistoIndex(iPtBin);
            fChanHist[3]->Fill(retCodeAnalysisCuts);
          }
        }
      }
      else
      {
        if(fKeepOnlyBkgFromHIJING) {
            isCandInjected = AliVertexingHFUtils::IsCandidateInjected(d, mcHeader, arrayMC);
        }
        if(!isCandInjected) {
            indexMCpiKK = GetBackgroundHistoIndex(iPtBin);
            indexMCKKpi = GetBackgroundHistoIndex(iPtBin);
            fChanHist[2]->Fill(retCodeAnalysisCuts);
            orig=6;
        }

        labDplus = d->MatchToMC(411, arrayMC, nProng, pdgDstoKKpi);
        if (labDplus >= 0)
        {
          partDs = (AliAODMCParticle*)arrayMC->At(labDplus);
          Int_t labDau0 = ((AliAODTrack *)d->GetDaughter(0))->GetLabel();
          AliAODMCParticle *p = (AliAODMCParticle *)arrayMC->UncheckedAt(TMath::Abs(labDau0));
          pdgCode0 = TMath::Abs(p->GetPdgCode());
        }
      }
      if(partDs){
        orig = AliVertexingHFUtils::CheckOrigin(arrayMC,partDs,kTRUE);
        orig_HijingCheck = AliVertexingHFUtils::CheckOrigin(arrayMC,partDs,kFALSE);
      }
       
    }

    if (isKKpi)
    {
      if (fDoRotBkg && TMath::Abs(massKK_KKpi - massPhi) <= fMaxDeltaPhiMass4Rot)
        GenerateRotBkg(d, 1, iPtBin);

      fMassHist[index]->Fill(invMass_KKpi, weightKKpi);
      fPtVsMass->Fill(invMass_KKpi, ptCand, weightKKpi);

      if (fDoBkgPhiSB && (0.010 < TMath::Abs(massKK_KKpi - massPhi)) && (TMath::Abs(massKK_KKpi - massPhi) < 0.030))
      {
        if (massKK_KKpi < massPhi)
          fMassLSBkgHistPhi[iPtBin]->Fill(invMass_KKpi);
        else
          fMassRSBkgHistPhi[iPtBin]->Fill(invMass_KKpi);
      }

      if (isPhiKKpi)
      {
        fMassHistPhi[index]->Fill(invMass_KKpi, weightKKpi);
        fPtVsMassPhi->Fill(invMass_KKpi, ptCand, weightKKpi);
      }
      if (isK0starKKpi)
      {
        fMassHistK0st[index]->Fill(invMass_KKpi, weightKKpi);
        fPtVsMassK0st->Fill(invMass_KKpi, ptCand, weightKKpi);
      }
      if (fReadMC && indexMCKKpi != -1)
      {
        fMassHist[indexMCKKpi]->Fill(invMass_KKpi, weightKKpi);
        if (isPhiKKpi)
          fMassHistPhi[indexMCKKpi]->Fill(invMass_KKpi, weightKKpi);
        if (isK0starKKpi)
          fMassHistK0st[indexMCKKpi]->Fill(invMass_KKpi, weightKKpi);
      }
    }
    if (ispiKK)
    {
      if (fDoRotBkg && TMath::Abs(massKK_piKK - massPhi) <= fMaxDeltaPhiMass4Rot)
        GenerateRotBkg(d, 2, iPtBin);

      fMassHist[index]->Fill(invMass_piKK, weightpiKK);
      fPtVsMass->Fill(invMass_piKK, ptCand, weightpiKK);

      if (fDoBkgPhiSB && (0.010 < TMath::Abs(massKK_piKK - massPhi)) && (TMath::Abs(massKK_piKK - massPhi) < 0.030))
      {
        if (massKK_piKK < massPhi)
          fMassLSBkgHistPhi[iPtBin]->Fill(invMass_piKK);
        else
          fMassRSBkgHistPhi[iPtBin]->Fill(invMass_piKK);
      }

      if (isPhipiKK)
      {
        fMassHistPhi[index]->Fill(invMass_piKK, weightpiKK);
        fPtVsMassPhi->Fill(invMass_piKK, ptCand, weightpiKK);
      }
      if (isK0starpiKK)
      {
        fMassHistK0st[index]->Fill(invMass_piKK, weightpiKK);
        fPtVsMassK0st->Fill(invMass_piKK, ptCand, weightpiKK);
      }
      if (fReadMC && indexMCpiKK != -1)
      {
        fMassHist[indexMCpiKK]->Fill(invMass_piKK, weightpiKK);
        if (isPhipiKK)
          fMassHistPhi[indexMCpiKK]->Fill(invMass_piKK, weightpiKK);
        if (isK0starpiKK)
          fMassHistK0st[indexMCpiKK]->Fill(invMass_piKK, weightpiKK);
      }
    }

    ///////////////////// CODE FOR NSPARSES /////////////////////////
    if (fFillSparse && (isPhiKKpi || isPhipiKK))
    {
      Double_t deltaMassKK = -999.;
      Double_t dlen = d->DecayLength();
      Double_t dlenxy = d->DecayLengthXY();
      Double_t normdlxy = d->NormalizedDecayLengthXY();
      Double_t cosp = d->CosPointingAngle();
      Double_t cospxy = d->CosPointingAngleXY();
      Double_t sigvert = d->GetSigmaVert();
      Double_t cosPiDs = -99.;
      Double_t cosPiKPhiNoabs = -99.;
      Double_t cosPiKPhi = -99.;
      Double_t normIP = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(d, aod->GetMagneticField());
      Double_t absimpparxy = TMath::Abs(d->ImpParXY());
      //variables for ML application
      AliAODPidHF *Pid_HF = nullptr;
      Double_t modelPred = -1.;
      bool isMLsel = true;
      if(fApplyML)
        Pid_HF = fAnalysisCuts->GetPidHF();

      if (isPhiKKpi)
      {
        deltaMassKK = TMath::Abs(massKK_KKpi-massPhi);
        cosPiDs = d->CosPiDsLabFrameKKpi();
        cosPiKPhi = d->CosPiKPhiRFrameKKpi();
        cosPiKPhiNoabs = cosPiKPhi * cosPiKPhi * cosPiKPhi;
        cosPiKPhi = TMath::Abs(cosPiKPhiNoabs);

        if(fApplyML)
        {
          isMLsel = fMLResponse->IsSelected(modelPred, d, aod->GetMagneticField(), Pid_HF, 0);

          if(fEnablePIDMLSparses && (!fReadMC || (indexMCKKpi == GetSignalHistoIndex(iPtBin) && orig == 4)))
          {
            Double_t var4nSparsePID[knVarPID] = {ptCand, modelPred,
                                                fMLResponse->GetVariable("nsigTPC_Pi_0"), fMLResponse->GetVariable("nsigTPC_K_0"),
                                                fMLResponse->GetVariable("nsigTOF_Pi_0"), fMLResponse->GetVariable("nsigTOF_K_0"),
                                                fMLResponse->GetVariable("nsigTPC_Pi_1"), fMLResponse->GetVariable("nsigTPC_K_1"),
                                                fMLResponse->GetVariable("nsigTOF_Pi_1"), fMLResponse->GetVariable("nsigTOF_K_1"),
                                                fMLResponse->GetVariable("nsigTPC_Pi_2"), fMLResponse->GetVariable("nsigTPC_K_2"),
                                                fMLResponse->GetVariable("nsigTOF_Pi_2"), fMLResponse->GetVariable("nsigTOF_K_2")};

            Double_t var4nSparsePIDcomb[knVarPIDcomb] = {ptCand, modelPred,
                                                        fMLResponse->GetVariable("nsigComb_Pi_0"), fMLResponse->GetVariable("nsigComb_K_0"),
                                                        fMLResponse->GetVariable("nsigComb_Pi_1"), fMLResponse->GetVariable("nsigComb_K_1"),
                                                        fMLResponse->GetVariable("nsigComb_Pi_2"), fMLResponse->GetVariable("nsigComb_K_2")};

            fnSparseNsigmaPIDVsML[0]->Fill(var4nSparsePID);
            fnSparseNsigmaPIDVsML[1]->Fill(var4nSparsePIDcomb);
          }
        }

        Double_t var4nSparse[knVarForSparse] = {invMass_KKpi, ptCand, deltaMassKK * 1000, dlen * 1000, dlenxy * 1000,
                                                normdlxy, cosp * 100, cospxy * 100, sigvert * 1000, cosPiDs * 10, cosPiKPhi * 10,
                                                TMath::Abs(normIP), absimpparxy * 10000, modelPred};

        if(!fApplyML || isMLsel)
        {
          if (!fReadMC)
          {
            fnSparse->Fill(var4nSparse);
          }
          else
          {
            if (indexMCKKpi == GetSignalHistoIndex(iPtBin))
            {
              if (orig == 4) fnSparseMC[2]->Fill(var4nSparse);
              else if (orig == 5) fnSparseMC[3]->Fill(var4nSparse);
            }
            else if(indexMCKKpi == GetBackgroundHistoIndex(iPtBin) && fFillBkgSparse)
            {
                fnSparseMC[4]->Fill(var4nSparse);
            }
            else if (fFillSparseDplus && labDplus >= 0 && pdgCode0 == 321)
            {
              if (orig == 4) fnSparseMCDplus[2]->Fill(var4nSparse);
              else if (orig == 5) fnSparseMCDplus[3]->Fill(var4nSparse);
            }
          }
        }
      }
      if (isPhipiKK)
      {
        deltaMassKK = TMath::Abs(massKK_piKK - massPhi);
        cosPiDs = d->CosPiDsLabFramepiKK();
        cosPiKPhi = d->CosPiKPhiRFramepiKK();
        cosPiKPhiNoabs = cosPiKPhi * cosPiKPhi * cosPiKPhi;
        cosPiKPhi = TMath::Abs(cosPiKPhiNoabs);

        if(fApplyML)
        {
          isMLsel = fMLResponse->IsSelected(modelPred, d, aod->GetMagneticField(), Pid_HF, 1);
          if(fEnablePIDMLSparses && (!fReadMC || (indexMCpiKK == GetSignalHistoIndex(iPtBin) && orig==4)))
          {
            Double_t var4nSparsePID[knVarPID] = {ptCand, modelPred,
                                                fMLResponse->GetVariable("nsigTPC_Pi_0"), fMLResponse->GetVariable("nsigTPC_K_0"),
                                                fMLResponse->GetVariable("nsigTOF_Pi_0"), fMLResponse->GetVariable("nsigTOF_K_0"),
                                                fMLResponse->GetVariable("nsigTPC_Pi_1"), fMLResponse->GetVariable("nsigTPC_K_1"),
                                                fMLResponse->GetVariable("nsigTOF_Pi_1"), fMLResponse->GetVariable("nsigTOF_K_1"),
                                                fMLResponse->GetVariable("nsigTPC_Pi_2"), fMLResponse->GetVariable("nsigTPC_K_2"),
                                                fMLResponse->GetVariable("nsigTOF_Pi_2"), fMLResponse->GetVariable("nsigTOF_K_2")};

            Double_t var4nSparsePIDcomb[knVarPIDcomb] = {ptCand, modelPred,
                                                        fMLResponse->GetVariable("nsigComb_Pi_0"), fMLResponse->GetVariable("nsigComb_K_0"),
                                                        fMLResponse->GetVariable("nsigComb_Pi_1"), fMLResponse->GetVariable("nsigComb_K_1"),
                                                        fMLResponse->GetVariable("nsigComb_Pi_2"), fMLResponse->GetVariable("nsigComb_K_2")};

            fnSparseNsigmaPIDVsML[0]->Fill(var4nSparsePID);
            fnSparseNsigmaPIDVsML[1]->Fill(var4nSparsePIDcomb);
          }
        }

        Double_t var4nSparse[knVarForSparse] = {invMass_piKK, ptCand, deltaMassKK * 1000, dlen * 1000, dlenxy * 1000,
                                                normdlxy, cosp * 100, cospxy * 100, sigvert * 1000, cosPiDs * 10, cosPiKPhi * 10,
                                                TMath::Abs(normIP), absimpparxy * 10000, modelPred};

        if(!fApplyML || isMLsel)
        {
          if (!fReadMC)
          {
            fnSparse->Fill(var4nSparse);
          }
          else
          {
            if (indexMCpiKK == GetSignalHistoIndex(iPtBin))
            {
              if (orig == 4) fnSparseMC[2]->Fill(var4nSparse);
              else if (orig == 5) fnSparseMC[3]->Fill(var4nSparse);
            }
            else if(indexMCpiKK == GetBackgroundHistoIndex(iPtBin) && fFillBkgSparse)
            {
                fnSparseMC[4]->Fill(var4nSparse);
            }
            else if (fFillSparseDplus && labDplus >= 0 && pdgCode0 == 211)
            {
              if (orig == 4) fnSparseMCDplus[2]->Fill(var4nSparse);
              else if (orig == 5) fnSparseMCDplus[3]->Fill(var4nSparse);
            }
          }
        }
      }
    }

    if(fFillImpParSparse && (isPhiKKpi || isPhipiKK))
    {
      Double_t impParxy = d->ImpParXY() * 10000.;
      Double_t array4ImpPar[3] = {invMass_KKpi, ptCand, impParxy};
      if (isPhiKKpi)
      {
        if (!fReadMC)
            fImpParSparse->Fill(array4ImpPar);
        else
        {
          if (orig == 4 && indexMCKKpi == GetSignalHistoIndex(iPtBin))
            fImpParSparseMC[0]->Fill(array4ImpPar);
          else if (orig == 5 && indexMCKKpi == GetSignalHistoIndex(iPtBin))
          {
            fImpParSparseMC[1]->Fill(array4ImpPar);
            trueImpParDsFromB = GetTrueImpactParameterDstoPhiPi(mcHeader, arrayMC, partDs) * 10000;
            Double_t array4ImpParTrueB[3] = {invMass_KKpi, ptCand, trueImpParDsFromB};
            fImpParSparseMC[2]->Fill(array4ImpParTrueB);
          }
          else if (indexMCKKpi == GetBackgroundHistoIndex(iPtBin) && fFillBkgSparse)
            fImpParSparseMC[3]->Fill(array4ImpPar);
        }
      }
      if (isPhipiKK)
      {
        if (!fReadMC)
          fImpParSparse->Fill(array4ImpPar);
        else
        {
          if (orig == 4 && indexMCpiKK == GetSignalHistoIndex(iPtBin))
            fImpParSparseMC[0]->Fill(array4ImpPar);
          else if (orig == 5 && indexMCpiKK == GetSignalHistoIndex(iPtBin))
          {
            fImpParSparseMC[1]->Fill(array4ImpPar);
            trueImpParDsFromB = GetTrueImpactParameterDstoPhiPi(mcHeader, arrayMC, partDs) * 10000;
            Double_t array4ImpParTrueB[3] = {invMass_piKK, ptCand, trueImpParDsFromB};
            fImpParSparseMC[2]->Fill(array4ImpParTrueB);
          }
          else if (indexMCpiKK == GetBackgroundHistoIndex(iPtBin) && fFillBkgSparse)
            fImpParSparseMC[3]->Fill(array4ImpPar);
        }
      }
    }

    ///////////////////CODE FOR ML TREE/////////////////////////////
    if (fCreateMLtree && (isPhiKKpi || isPhipiKK))
    {
      AliAODPidHF *Pid_HF = fAnalysisCuts->GetPidHF();

      if (isPhiKKpi)
      {
        bool issignal = kFALSE;
        bool isbkg = kFALSE;
        bool isprompt = kFALSE;
        bool isFD = kFALSE;
        bool isrefl = kFALSE;
        bool isSignalWoQuark = kFALSE;

        if(fReadMC) {
          if(labDs >= 0) {
            if(orig == 4 || orig_HijingCheck==4)
              isprompt = kTRUE;
            else if(orig == 5|| orig_HijingCheck==5)
              isFD = kTRUE;

            if(orig >= 4 && orig_HijingCheck >= 4) {
              if(pdgCode0 == 321)
                issignal = kTRUE;
              else if(pdgCode0 == 211)
                isrefl = kTRUE;
            }
            if(orig < 4 && orig_HijingCheck >=4)
              isSignalWoQuark = kTRUE;
          } else {
            if(!isCandInjected)
              isbkg = kTRUE;
            if(labDplus >= 0)
            {
              fMLhandler->SetIsDplustoKKpi(kTRUE);
              if(orig == 4)
                isprompt = kTRUE;
              else if(orig == 5)
                isFD = kTRUE;
            }
          }
        }

        fMLhandler->SetCandidateType(issignal, isbkg, isprompt, isFD, isrefl);
        fMLhandler->SetlsSignalWoQuark(isSignalWoQuark);
        fMLhandler->SetVariables(d, aod->GetMagneticField(), AliHFMLVarHandlerDstoKKpi::kKKpi, Pid_HF);
        if(!(fReadMC && !issignal && !isbkg && !isprompt && !isFD && !isrefl))
          fMLhandler->FillTree();
      }

      if (isPhipiKK)
      {
        bool issignal = kFALSE;
        bool isbkg = kFALSE;
        bool isprompt = kFALSE;
        bool isFD = kFALSE;
        bool isrefl = kFALSE;

        if(fReadMC) {
          if(labDs >= 0) {
            if(orig == 4 || orig_HijingCheck==4)
              isprompt = kTRUE;
            else if(orig == 5 || orig_HijingCheck==5)
              isFD = kTRUE;

            if(orig >= 4 && orig_HijingCheck >=4) {
              if(pdgCode0 == 211)
                issignal = kTRUE;
              else if(pdgCode0 == 321)
                isrefl = kTRUE;
            }
            if(orig < 4 && orig_HijingCheck >= 4)
              isSignalWoQuark = kTRUE;
          } else {
            if(!isCandInjected)
              isbkg = kTRUE;
            if(labDplus >= 0)
            {
              fMLhandler->SetIsDplustoKKpi(kTRUE);
              if(orig == 4)
                isprompt = kTRUE;
              else if(orig == 5)
                isFD = kTRUE;
            }
          }
        }

        fMLhandler->SetCandidateType(issignal, isbkg, isprompt, isFD, isrefl);
        fMLhandler->SetlsSignalWoQuark(isSignalWoQuark);
        fMLhandler->SetVariables(d, aod->GetMagneticField(), AliHFMLVarHandlerDstoKKpi::kpiKK, Pid_HF);
        if(!(fReadMC && !issignal && !isbkg && !isprompt && !isFD && !isrefl))
          fMLhandler->FillTree();
      }
    }

    ////////////////////////////////////////////////////////////////
    if (fDoCutVarHistos)
    {
      Double_t dlen = d->DecayLength();
      Double_t dlenxy = d->DecayLengthXY();
      Double_t normdlxy = d->NormalizedDecayLengthXY();
      Double_t cosp = d->CosPointingAngle();
      Double_t cospxy = d->CosPointingAngleXY();
      Double_t pt0 = d->PtProng(0);
      Double_t pt1 = d->PtProng(1);
      Double_t pt2 = d->PtProng(2);
      Double_t sigvert = d->GetSigmaVert();
      Double_t cosPiDs = -99.;
      Double_t cosPiKPhi = -99.;
      Double_t normIP = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(d, aod->GetMagneticField());
      Double_t sumD02 = d->Getd0Prong(0) * d->Getd0Prong(0) + d->Getd0Prong(1) * d->Getd0Prong(1) + d->Getd0Prong(2) * d->Getd0Prong(2);
      Double_t dca = d->GetDCA();
      Double_t ptmax = 0;
      for (Int_t i = 0; i < 3; i++)
        if (d->PtProng(i) > ptmax)
          ptmax = d->PtProng(i);

      fCosPHist[index]->Fill(cosp);
      fCosPxyHist[index]->Fill(cospxy);
      fDLenHist[index]->Fill(dlen);
      fDLenxyHist[index]->Fill(dlenxy);
      fNDLenxyHist[index]->Fill(normdlxy);
      fSigVertHist[index]->Fill(sigvert);
      fSumd02Hist[index]->Fill(sumD02);
      fPtMaxHist[index]->Fill(ptmax);
      fDCAHist[index]->Fill(dca);
      fNormIPHist[index]->Fill(normIP);
      fCosPiDsHist[index]->Fill(cosPiDs);
      fCosPiKPhiHist[index]->Fill(cosPiKPhi);
      fPtProng0Hist[index]->Fill(pt0);
      fPtProng1Hist[index]->Fill(pt1);
      fPtProng2Hist[index]->Fill(pt2);
      if (isKKpi)
      {
        cosPiDs = d->CosPiDsLabFrameKKpi();
        cosPiKPhi = d->CosPiKPhiRFrameKKpi();
        cosPiKPhi = TMath::Abs(cosPiKPhi * cosPiKPhi * cosPiKPhi);

        if (!fReadMC)
        {
          fCosPHist3D->Fill(invMass_KKpi, ptCand, cosp);
          fCosPxyHist3D->Fill(invMass_KKpi, ptCand, cospxy);
          fDLenHist3D->Fill(invMass_KKpi, ptCand, dlen);
          fDLenxyHist3D->Fill(invMass_KKpi, ptCand, dlenxy);
          fNDLenxyHist3D->Fill(invMass_KKpi, ptCand, normdlxy);
          fSigVertHist3D->Fill(invMass_KKpi, ptCand, sigvert);
          fDCAHist3D->Fill(invMass_KKpi, ptCand, dca);
          fNormIPHist3D->Fill(invMass_KKpi, ptCand, normIP);
          fCosPiDsHist3D->Fill(invMass_KKpi, ptCand, cosPiDs);
          fCosPiKPhiHist3D->Fill(invMass_KKpi, ptCand, cosPiKPhi);
          fPtProng0Hist3D->Fill(invMass_KKpi, ptCand, pt0);
          fPtProng1Hist3D->Fill(invMass_KKpi, ptCand, pt1);
          fPtProng2Hist3D->Fill(invMass_KKpi, ptCand, pt2);
        }
        fDalitz[index]->Fill(massKK_KKpi, massKp);
        if (isPhiKKpi)
          fDalitzPhi[index]->Fill(massKK_KKpi, massKp);
        if (isK0starKKpi)
          fDalitzK0st[index]->Fill(massKK_KKpi, massKp);
        if (fReadMC && indexMCKKpi != -1)
        {
          fDalitz[indexMCKKpi]->Fill(massKK_KKpi, massKp);
          if (isPhiKKpi)
            fDalitzPhi[indexMCKKpi]->Fill(massKK_KKpi, massKp);
          if (isK0starKKpi)
            fDalitzK0st[indexMCKKpi]->Fill(massKK_KKpi, massKp);
          fCosPHist[indexMCKKpi]->Fill(cosp);
          fCosPxyHist[indexMCKKpi]->Fill(cospxy);
          fDLenHist[indexMCKKpi]->Fill(dlen);
          fDLenxyHist[indexMCKKpi]->Fill(dlenxy);
          fNDLenxyHist[indexMCKKpi]->Fill(normdlxy);
          fSigVertHist[indexMCKKpi]->Fill(sigvert);
          fSumd02Hist[indexMCKKpi]->Fill(sumD02);
          fPtMaxHist[indexMCKKpi]->Fill(ptmax);
          fPtCandHist[indexMCKKpi]->Fill(ptCand);
          fDCAHist[indexMCKKpi]->Fill(dca);
          fNormIPHist[indexMCKKpi]->Fill(normIP);
          fCosPiDsHist[indexMCKKpi]->Fill(cosPiDs);
          fCosPiKPhiHist[indexMCKKpi]->Fill(cosPiKPhi);
          fPtProng0Hist[indexMCKKpi]->Fill(pt0);
          fPtProng1Hist[indexMCKKpi]->Fill(pt1);
          fPtProng2Hist[indexMCKKpi]->Fill(pt2);
        }
      }
      if (ispiKK)
      {
        cosPiDs = d->CosPiDsLabFramepiKK();
        cosPiKPhi = d->CosPiKPhiRFramepiKK();
        cosPiKPhi = TMath::Abs(cosPiKPhi * cosPiKPhi * cosPiKPhi);

        if (!fReadMC)
        {
          fCosPHist3D->Fill(invMass_piKK, ptCand, cosp);
          fCosPxyHist3D->Fill(invMass_piKK, ptCand, cospxy);
          fDLenHist3D->Fill(invMass_piKK, ptCand, dlen);
          fDLenxyHist3D->Fill(invMass_piKK, ptCand, dlenxy);
          fNDLenxyHist3D->Fill(invMass_piKK, ptCand, normdlxy);
          fSigVertHist3D->Fill(invMass_piKK, ptCand, sigvert);
          fDCAHist3D->Fill(invMass_piKK, ptCand, dca);
          fNormIPHist3D->Fill(invMass_piKK, ptCand, normIP);
          fCosPiDsHist3D->Fill(invMass_piKK, ptCand, cosPiDs);
          fCosPiKPhiHist3D->Fill(invMass_piKK, ptCand, cosPiKPhi);
          fPtProng0Hist3D->Fill(invMass_piKK, ptCand, pt0);
          fPtProng1Hist3D->Fill(invMass_piKK, ptCand, pt1);
          fPtProng2Hist3D->Fill(invMass_piKK, ptCand, pt2);
        }
        fDalitz[index]->Fill(massKK_piKK, masspK);
        if (isPhipiKK)
          fDalitzPhi[index]->Fill(massKK_piKK, masspK);
        if (isK0starpiKK)
          fDalitzK0st[index]->Fill(massKK_piKK, masspK);

        if (fReadMC && indexMCpiKK != -1)
        {
          fDalitz[indexMCpiKK]->Fill(massKK_piKK, masspK);
          if (isPhipiKK)
            fDalitzPhi[indexMCpiKK]->Fill(massKK_piKK, masspK);
          if (isK0starpiKK)
            fDalitzK0st[indexMCpiKK]->Fill(massKK_piKK, masspK);
          fCosPHist[indexMCpiKK]->Fill(cosp);
          fCosPxyHist[indexMCpiKK]->Fill(cospxy);
          fDLenHist[indexMCpiKK]->Fill(dlen);
          fDLenxyHist[indexMCpiKK]->Fill(dlenxy);
          fNDLenxyHist[indexMCpiKK]->Fill(normdlxy);
          fSigVertHist[indexMCpiKK]->Fill(sigvert);
          fSumd02Hist[indexMCpiKK]->Fill(sumD02);
          fPtMaxHist[indexMCpiKK]->Fill(ptmax);
          fPtCandHist[indexMCpiKK]->Fill(ptCand);
          fDCAHist[indexMCpiKK]->Fill(dca);
          fNormIPHist[indexMCpiKK]->Fill(normIP);
          fCosPiDsHist[indexMCpiKK]->Fill(cosPiDs);
          fCosPiKPhiHist[indexMCpiKK]->Fill(cosPiKPhi);
          fPtProng0Hist[indexMCpiKK]->Fill(pt0);
          fPtProng1Hist[indexMCpiKK]->Fill(pt1);
          fPtProng2Hist[indexMCpiKK]->Fill(pt2);
        }
      }
    }

    if (unsetvtx)
      d->UnsetOwnPrimaryVtx();
    if (recVtx)
      fAnalysisCuts->CleanOwnPrimaryVtx(d, aod, origownvtx);
  }

  fCounter->StoreCandidates(aod, nFiltered, kTRUE);
  fCounter->StoreCandidates(aod, nSelected, kFALSE);

  PostData(1, fOutput);
  PostData(3, fCounter);

  return;
}

//_________________________________________________________________
void AliAnalysisTaskSEDs::Terminate(Option_t * /*option*/)
{
  /// Terminate analysis
  //
  if (fDebug > 1)
    printf("AnalysisTaskSEDs: Terminate() \n");
  fOutput = dynamic_cast<TList *>(GetOutputData(1));
  if (!fOutput)
  {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F *>(fOutput->FindObject("hNEvents"));
  if (fHistNEvents)
  {
    printf("Number of analyzed events = %d\n", (Int_t)fHistNEvents->GetBinContent(2));
  }
  else
  {
    printf("ERROR: fHistNEvents not available\n");
    return;
  }
  return;
}

//_________________________________________________________________
void AliAnalysisTaskSEDs::FillMCGenAccHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader)
{
  /// Fill MC histos for cuts study
  ///    - at GenLimAccStep and AccStep (if fFillAcceptanceLevel=kFALSE)
  ///    - at AccStep (if fFillAcceptanceLevel=kTRUE)

  Int_t nProng = 3;
  Double_t zMCVertex = mcHeader->GetVtxZ(); //vertex MC
  if (TMath::Abs(zMCVertex) <= fAnalysisCuts->GetMaxVtxZ())
  {
    for (Int_t iPart = 0; iPart < arrayMC->GetEntriesFast(); iPart++)
    {

      AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(iPart));

      if (TMath::Abs(mcPart->GetPdgCode()) == 431)
      {
        Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, kTRUE); //Prompt = 4, FeedDown = 5

        Int_t deca = 0;
        Bool_t isGoodDecay = kFALSE;
        Int_t labDau[3] = {-1, -1, -1};
        Bool_t isFidAcc = kFALSE;
        Bool_t isDaugInAcc = kFALSE;

        deca = AliVertexingHFUtils::CheckDsDecay(arrayMC, mcPart, labDau);
        if (deca == 1)
          isGoodDecay = kTRUE; // == 1 -> Phi pi -> kkpi

        if (labDau[0] == -1)
          continue; //protection against unfilled array of labels

        if (isGoodDecay)
        {
          Double_t pt = mcPart->Pt();
          Double_t rapid = mcPart->Y();
          isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(pt, rapid);
          isDaugInAcc = CheckDaugAcc(arrayMC, nProng, labDau);

          if ((fFillAcceptanceLevel && isFidAcc && isDaugInAcc) || (!fFillAcceptanceLevel && TMath::Abs(rapid)<0.5))
          {
            Double_t var4nSparseAcc[knVarForSparseAcc] = {pt, rapid * 10};            
            if (orig == 4)
              fnSparseMC[0]->Fill(var4nSparseAcc);
            if (orig == 5)
              fnSparseMC[1]->Fill(var4nSparseAcc);
          }
        }
      }
      if (fFillSparseDplus && TMath::Abs(mcPart->GetPdgCode()) == 411)
      {
        Int_t orig = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, kTRUE); //Prompt = 4, FeedDown = 5

        Int_t deca = 0;
        Bool_t isGoodDecay = kFALSE;
        Int_t labDau[3] = {-1, -1, -1};
        Bool_t isFidAcc = kFALSE;
        Bool_t isDaugInAcc = kFALSE;

        deca = AliVertexingHFUtils::CheckDplusKKpiDecay(arrayMC, mcPart, labDau);
        if (deca == 1)
          isGoodDecay = kTRUE; // == 1 -> Phi pi -> kkpi

        if (labDau[0] == -1)
          continue; //protection against unfilled array of labels

        if (isGoodDecay)
        {
          Double_t pt = mcPart->Pt();
          Double_t rapid = mcPart->Y();
          isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(pt, rapid);
          isDaugInAcc = CheckDaugAcc(arrayMC, nProng, labDau);

          if ((fFillAcceptanceLevel && isFidAcc && isDaugInAcc) || (!fFillAcceptanceLevel && TMath::Abs(rapid)<0.5))
          {
            Double_t var4nSparseAcc[knVarForSparseAcc] = {pt, rapid * 10};
            if (orig == 4)
              fnSparseMCDplus[0]->Fill(var4nSparseAcc);
            if (orig == 5)
              fnSparseMCDplus[1]->Fill(var4nSparseAcc);
          }
        }
      }
    }
  }
}

//_________________________________________________________________
Bool_t AliAnalysisTaskSEDs::CheckDaugAcc(TClonesArray *arrayMC, Int_t nProng, Int_t *labDau)
{
  /// check if the decay products are in the good eta and pt range

  for (Int_t iProng = 0; iProng < nProng; iProng++)
  {
    AliAODMCParticle *mcPartDaughter = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDau[iProng]));
    if (!mcPartDaughter)
    {
      return kFALSE;
    }
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > 0.9 || pt < 0.1)
    {
      return kFALSE;
    }
  }
  return kTRUE;
}

//_________________________________________________________________
void AliAnalysisTaskSEDs::GenerateRotBkg(AliAODRecoDecayHF3Prong *d, Int_t dec, Int_t iPtBin)
{

  const Int_t nprongs = 3;
  Double_t PxProng[nprongs], PyProng[nprongs], PzProng[nprongs], P2Prong[nprongs], mProng[nprongs];
  Double_t Px, Py, Pz, P2;
  UInt_t pdg[3] = {321, 321, 211};
  int idPion = 2;
  if (dec == 2)
  {
    pdg[0] = 211;
    pdg[2] = 321;
    idPion = 0;
  }

  for (Int_t ip = 0; ip < nprongs; ip++)
  {
    PxProng[ip] = d->PxProng(ip);
    PyProng[ip] = d->PxProng(ip);
    PzProng[ip] = d->PzProng(ip);
    P2Prong[ip] = d->P2Prong(ip);
    mProng[ip] = TDatabasePDG::Instance()->GetParticle(pdg[ip])->Mass();
  }

  for (Int_t i = 0; i < 9; i++)
  { //9 rotations implemented for the pion track around Pi
    Px = 0.;
    Py = 0.;
    Pz = 0.;

    Double_t phirot = TMath::Pi() * (5 / 6. + 1 / 27. * i);

    PxProng[idPion] = PxProng[idPion] * TMath::Cos(phirot) - PyProng[idPion] * TMath::Sin(phirot);
    PyProng[idPion] = PxProng[idPion] * TMath::Sin(phirot) + PyProng[idPion] * TMath::Cos(phirot);

    for (Int_t j = 0; j < nprongs; j++)
    {
      Px += PxProng[j];
      Py += PyProng[j];
      Pz += PzProng[j];
    }
    P2 = Px * Px + Py * Py + Pz * Pz;

    Double_t energysum = 0.;
    for (Int_t j = 0; j < nprongs; j++)
    {
      energysum += TMath::Sqrt(mProng[j] * mProng[j] + P2Prong[j]);
    }
    Double_t mass = TMath::Sqrt(energysum * energysum - P2);
    if ((fminMass <= mass) && (mass < fmaxMass))
      fMassRotBkgHistPhi[iPtBin]->Fill(mass);
  }
}

//_________________________________________________________________________
void AliAnalysisTaskSEDs::CreateCutVarsAndEffSparses()
{

  Double_t massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  Int_t nInvMassBins = (Int_t)(0.7 / fMassBinSize + 0.5);
  Double_t minMass = massDs - 0.5 * nInvMassBins * fMassBinSize;
  Double_t maxMass = massDs + 0.5 * nInvMassBins * fMassBinSize;

  Int_t nSparseAxes = knVarForSparse;
  if(!fApplyML)
    nSparseAxes--;

  Int_t nPtBins = (Int_t)fPtLimits[fNPtBins];
  if(fUseFinPtBinsForSparse)
    nPtBins = nPtBins*10;

  Int_t nBinsReco[knVarForSparse];
  Double_t xminReco[knVarForSparse];
  Double_t xmaxReco[knVarForSparse];
  TString axis[knVarForSparse] = {"invMassDsAllPhi", "#it{p}_{T}", "#Delta Mass(KK)", "dlen", "dlen_{xy}", "normdl_{xy}", "cosP", "cosP_{xy}",
                                  "sigVert", "cosPiDs", "|cosPiKPhi^{3}|", "normIP", "ImpPar_{xy}", "ML model output"};

  if (fSystem == kpp)
  {
    std::vector<Int_t> nBinsRecoVec = {nInvMassBins, nPtBins, 30, 20, 20, 20, 20, 20, 14, 6, 6, 12, 30, fNMLBins};
    std::vector<Double_t> xminRecoVec = {minMass, 0., 0., 0., 0., 0., 90., 90., 0., 7., 0., 0., 0., fMLOutputMin};
    std::vector<Double_t> xmaxRecoVec = {maxMass, fPtLimits[fNPtBins], 15., 100., 100., 10., 100., 100., 70., 10., 3., 6., 300., fMLOutputMax};
    std::copy(nBinsRecoVec.begin(),nBinsRecoVec.end(),nBinsReco);
    std::copy(xminRecoVec.begin(),xminRecoVec.end(),xminReco);
    std::copy(xmaxRecoVec.begin(),xmaxRecoVec.end(),xmaxReco);
  }
  else if (fSystem == kPbPb)
  {
    nInvMassBins = (Int_t)(0.45 / fMassBinSize + 0.5);
    minMass = massDs - 0.5 * nInvMassBins * fMassBinSize;
    maxMass = massDs + 0.5 * nInvMassBins * fMassBinSize;
    std::vector<Int_t> nBinsRecoVec = {nInvMassBins, nPtBins, 15, 10, 10, 10, 10, 10, 14, 6, 6, 12, 30, fNMLBins};
    std::vector<Double_t> xminRecoVec = {minMass, 0., 0., 0., 0., 0., 95., 95., 0., 7., 0., 0., 0., fMLOutputMin};
    std::vector<Double_t> xmaxRecoVec = {maxMass, fPtLimits[fNPtBins], 15., 100., 100., 10., 100., 100., 70., 10., 3., 6., 300., fMLOutputMax};
    std::copy(nBinsRecoVec.begin(),nBinsRecoVec.end(),nBinsReco);
    std::copy(xminRecoVec.begin(),xminRecoVec.end(),xminReco);
    std::copy(xmaxRecoVec.begin(),xmaxRecoVec.end(),xmaxReco);
  }
  else if (fSystem == kUpgr)
  {
    std::vector<Int_t> nBinsRecoVec = {nInvMassBins, nPtBins, 40, 120, 120, 50, 60, 60, 30, 12, 12, 20, 100, fNMLBins};
    std::vector<Double_t> xminRecoVec = {minMass, 0., 0., 0., 0., 0., 0.97, 0.97, 0., 0.7, 0., 0., 0., fMLOutputMin};
    std::vector<Double_t> xmaxRecoVec = {maxMass, fPtLimits[fNPtBins], 20., 1200., 1200., 25., 1., 1., 150., 1., 0.3, 5., 50., fMLOutputMax};
    std::copy(nBinsRecoVec.begin(),nBinsRecoVec.end(),nBinsReco);
    std::copy(xminRecoVec.begin(),xminRecoVec.end(),xminReco);
    std::copy(xmaxRecoVec.begin(),xmaxRecoVec.end(),xmaxReco);
  }

  Int_t nBinsAcc[knVarForSparseAcc] = {nPtBins, 20};
  Double_t xminAcc[knVarForSparseAcc] = {0., -10.};
  Double_t xmaxAcc[knVarForSparseAcc] = {fPtLimits[fNPtBins], 10.};

  if (fReadMC)
  {
    TString label[3] = {"fromC", "fromB", "bkg"};
    Int_t nSparseReco = 2;
    if(fFillBkgSparse)
      nSparseReco = 3;

    for (Int_t iHist = 0; iHist < 2; iHist++)
    {
      TString titleSparse = Form("MC nSparse (%s)- %s", fFillAcceptanceLevel ? "Acc.Step" : "Gen.Acc.Step", label[iHist].Data());
      fnSparseMC[iHist] = new THnSparseF(Form("fnSparseAcc_%s", label[iHist].Data()), titleSparse.Data(), knVarForSparseAcc, nBinsAcc, xminAcc, xmaxAcc);
      fnSparseMC[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/c)");
      fnSparseMC[iHist]->GetAxis(1)->SetTitle("#it{y}");
      fOutput->Add(fnSparseMC[iHist]);

      //Dplus
      if (fFillSparseDplus)
      {
        titleSparse = Form("MC nSparse D^{+} (%s)- %s", fFillAcceptanceLevel ? "Acc.Step" : "Gen.Acc.Step", label[iHist].Data());
        fnSparseMCDplus[iHist] = new THnSparseF(Form("fnSparseAccDplus_%s", label[iHist].Data()), titleSparse.Data(), knVarForSparseAcc, nBinsAcc, xminAcc, xmaxAcc);
        fnSparseMCDplus[iHist]->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/c)");
        fnSparseMCDplus[iHist]->GetAxis(1)->SetTitle("#it{y}");
        fOutput->Add(fnSparseMCDplus[iHist]);
      }
    }
    for (Int_t iHist = 2; iHist < nSparseReco + 2; iHist++)
    {
      fnSparseMC[iHist] = new THnSparseF(Form("fnSparseReco_%s", label[iHist - 2].Data()), Form("MC nSparse (Reco Step)- %s", label[iHist - 2].Data()), nSparseAxes, nBinsReco, xminReco, xmaxReco);
      for (Int_t iAxis = 0; iAxis < nSparseAxes; iAxis++)
      {
        fnSparseMC[iHist]->GetAxis(iAxis)->SetTitle(Form("%s", axis[iAxis].Data()));
      }
      fOutput->Add(fnSparseMC[iHist]);

      //Dplus
      if (fFillSparseDplus && iHist<4)
      {
        fnSparseMCDplus[iHist] = new THnSparseF(Form("fnSparseRecoDplus_%s", label[iHist - 2].Data()), Form("MC nSparse D^{+} (Reco Step)- %s", label[iHist - 2].Data()), nSparseAxes, nBinsReco, xminReco, xmaxReco);
        for (Int_t iAxis = 0; iAxis < nSparseAxes; iAxis++)
        {
          fnSparseMCDplus[iHist]->GetAxis(iAxis)->SetTitle(Form("%s", axis[iAxis].Data()));
        }
        fOutput->Add(fnSparseMCDplus[iHist]);
      }
    }
  } //end MC
  else
  {
    fnSparse = new THnSparseF("fnSparse", "nSparse", nSparseAxes, nBinsReco, xminReco, xmaxReco);
    for (Int_t iAxis = 0; iAxis < nSparseAxes; iAxis++)
    {
      fnSparse->GetAxis(iAxis)->SetTitle(Form("%s", axis[iAxis].Data()));
    }
    fOutput->Add(fnSparse);
  }
}

//_________________________________________________________________________
void AliAnalysisTaskSEDs::CreateImpactParameterSparses()
{

  /// Histos for impact parameter study
  Double_t massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  Int_t nInvMassBins = (Int_t)(0.7 / fMassBinSize + 0.5);
  Double_t minMass = massDs - 0.5 * nInvMassBins * fMassBinSize;
  Double_t maxMass = massDs + 0.5 * nInvMassBins * fMassBinSize;

  //dimensions for THnSparse
  TString axTit[kVarForImpPar] = {"M_{K#pi#pi} (GeV/c^{2})", "#it{p}_{T} (GeV/c)", "Imp Par (#mum)"};

  Int_t nbins[kVarForImpPar] = {nInvMassBins, (Int_t)fPtLimits[fNPtBins] * 2, 1000};
  Double_t xmin[kVarForImpPar] = {minMass, 0., -1000};
  Double_t xmax[kVarForImpPar] = {maxMass, fPtLimits[fNPtBins], 1000};

  //mass, pt, imppar
  fImpParSparse = new THnSparseF("hMassPtImpParAll", "Mass vs. pt vs. imppar - All", kVarForImpPar, nbins, xmin, xmax);
  fImpParSparseMC[0] = new THnSparseF("hMassPtImpParPrompt", "Mass vs. pt vs. imppar - promptD", kVarForImpPar, nbins, xmin, xmax);
  fImpParSparseMC[1] = new THnSparseF("hMassPtImpParBfeed", "Mass vs. pt vs. imppar - DfromB", kVarForImpPar, nbins, xmin, xmax);
  fImpParSparseMC[2] = new THnSparseF("hMassPtImpParTrueBfeed", "Mass vs. pt vs. true imppar -DfromB", kVarForImpPar, nbins, xmin, xmax);
  fImpParSparseMC[3] = new THnSparseF("hMassPtImpParBkg", "Mass vs. pt vs. imppar - backgr.", kVarForImpPar, nbins, xmin, xmax);

  for(Int_t iax=0; iax<kVarForImpPar; iax++){
    fImpParSparse->GetAxis(iax)->SetTitle(axTit[iax].Data());
    for(Int_t ih=0; ih<4; ih++) fImpParSparseMC[ih]->GetAxis(iax)->SetTitle(axTit[iax].Data());
  }


  if (!fReadMC)
    fOutput->Add(fImpParSparse);
  else
  {
    for (Int_t iSparse = 0; iSparse < 4; iSparse++)
    {
      fOutput->Add(fImpParSparseMC[iSparse]);
    }
  }
}

//_________________________________________________________________________
void AliAnalysisTaskSEDs::CreatePIDMLSparses()
{
  Int_t nPtBins = (Int_t)fPtLimits[fNPtBins];
  if(fUseFinPtBinsForSparse)
     nPtBins = nPtBins*10;
  Int_t nPIDbins = 80;
  Double_t PIDmin = -20.;
  Double_t PIDmax = 20.;
  Double_t PIDcombMin = 0.;
  Double_t PIDcombMax = 40.;

  TString PIDvarnames[knVarPID] = {"#it{p}_{T}", "ML model output", "n#sigmaTPCPi_0", "n#sigmaTPCK_0", "n#sigmaTOFPi_0", "n#sigmaTOFK_0",
                                         "n#sigmaTPCPi_1", "n#sigmaTPCK_1", "n#sigmaTOFPi_1", "n#sigmaTOFK_1", "n#sigmaTPCPi_2", "n#sigmaTPCK_2",
                                         "n#sigmaTOFPi_2", "n#sigmaTOFK_2"};
  Int_t nBinsPID[knVarPID] = {nPtBins, fNMLBins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins,
                                    nPIDbins, nPIDbins};
  Double_t xminPID[knVarPID] = {0., fMLOutputMin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin, PIDmin};
  Double_t xmaxPID[knVarPID] = {fPtLimits[fNPtBins], fMLOutputMax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax, PIDmax};

  TString PIDvarnamesComb[knVarPIDcomb] = {"#it{p}_{T}", "ML model output", "n#sigmaCombPi_0", "n#sigmaCombK_0", "n#sigmaCombPi_1", "n#sigmaCombK_1",
                                                 "n#sigmaCombPi_2", "n#sigmaCombK_2"};
  Int_t nBinsPIDcomb[knVarPIDcomb] = {nPtBins, fNMLBins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins, nPIDbins};
  Double_t xminPIDcomb[knVarPIDcomb] = {0., fMLOutputMin, PIDcombMin, PIDcombMin, PIDcombMin, PIDcombMin, PIDcombMin, PIDcombMin};
  Double_t xmaxPIDcomb[knVarPIDcomb] = {fPtLimits[fNPtBins], fMLOutputMax, PIDcombMax, PIDcombMax, PIDcombMax, PIDcombMax, PIDcombMax, PIDcombMax};

  fnSparseNsigmaPIDVsML[0] = new THnSparseF("fnSparsePID", "nSparsePID", knVarPID, nBinsPID, xminPID, xmaxPID);
  for (Int_t iAxis = 0; iAxis < knVarPID; iAxis++)
    fnSparseNsigmaPIDVsML[0]->GetAxis(iAxis)->SetTitle(Form("%s", PIDvarnames[iAxis].Data()));
  fOutput->Add(fnSparseNsigmaPIDVsML[0]);

  fnSparseNsigmaPIDVsML[1] = new THnSparseF("fnSparsePIDcomb", "nSparsePIDcomb", knVarPIDcomb, nBinsPIDcomb, xminPIDcomb, xmaxPIDcomb);
  for (Int_t iAxis = 0; iAxis < knVarPIDcomb; iAxis++)
    fnSparseNsigmaPIDVsML[1]->GetAxis(iAxis)->SetTitle(Form("%s", PIDvarnamesComb[iAxis].Data()));
  fOutput->Add(fnSparseNsigmaPIDVsML[1]);
}

//_________________________________________________________________________________________________
Float_t AliAnalysisTaskSEDs::GetTrueImpactParameterDstoPhiPi(const AliAODMCHeader *mcHeader, TClonesArray *arrayMC, const AliAODMCParticle *partDs) const
{
  /// true impact parameter calculation

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partDs->XvYvZv(origD);
  Short_t charge = partDs->Charge();
  Double_t pXdauTrue[3], pYdauTrue[3], pZdauTrue[3];
  for (Int_t iDau = 0; iDau < 3; iDau++)
  {
    pXdauTrue[iDau] = 0.;
    pYdauTrue[iDau] = 0.;
    pZdauTrue[iDau] = 0.;
  }

  Int_t nDau = partDs->GetNDaughters();
  Int_t labelFirstDau = partDs->GetDaughterLabel(0);
  if (nDau == 2)
  {
    Int_t theDau = 0;
    for (Int_t iDau = 0; iDau < 2; iDau++)
    {
      Int_t ind = labelFirstDau + iDau;
      AliAODMCParticle *part = dynamic_cast<AliAODMCParticle *>(arrayMC->At(ind));
      if (!part)
      {
        AliError("Daughter particle not found in MC array");
        return 99999.;
      }
      Int_t pdgCode = TMath::Abs(part->GetPdgCode());
      if (pdgCode == 211)
      {
        pXdauTrue[theDau] = part->Px();
        pYdauTrue[theDau] = part->Py();
        pZdauTrue[theDau] = part->Pz();
        ++theDau;
      }
      else
      {
        Int_t nDauRes = part->GetNDaughters();
        if (nDauRes == 2)
        {
          Int_t labelFirstDauRes = part->GetDaughterLabel(0);
          for (Int_t iDauRes = 0; iDauRes < 2; iDauRes++)
          {
            Int_t indDR = labelFirstDauRes + iDauRes;
            AliAODMCParticle *partDR = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indDR));
            if (!partDR)
            {
              AliError("Daughter particle not found in MC array");
              return 99999.;
            }

            Int_t pdgCodeDR = TMath::Abs(partDR->GetPdgCode());
            if (pdgCodeDR == 321)
            {
              pXdauTrue[theDau] = partDR->Px();
              pYdauTrue[theDau] = partDR->Py();
              pZdauTrue[theDau] = partDR->Pz();
              ++theDau;
            }
          }
        }
      }
    }
  }
  else
  {
    AliError("Wrong number of decay prongs");
    return 99999.;
  }

  Double_t d0dummy[3] = {0., 0., 0.};
  AliAODRecoDecayHF aodDsMC(vtxTrue, origD, 3, charge, pXdauTrue, pYdauTrue, pZdauTrue, d0dummy);
  return aodDsMC.ImpParXY();
}

//________________________________________________________________________
void AliAnalysisTaskSEDs::SetPtBins(Int_t n, Float_t *lim)
{
  /// define pt bins for analysis
  if (n > kMaxPtBins)
  {
    printf("Max. number of Pt bins = %d\n", kMaxPtBins);
    fNPtBins = kMaxPtBins;
    fPtLimits[0] = 0.;
    fPtLimits[1] = 1.;
    fPtLimits[2] = 3.;
    fPtLimits[3] = 5.;
    fPtLimits[4] = 10.;
    for (Int_t iPt = 5; iPt < kMaxPtBins + iPt; iPt++)
      fPtLimits[iPt] = 999.;
  }
  else
  {
    fNPtBins = n;
    for (Int_t iPt = 0; iPt < fNPtBins + 1; iPt++)
      fPtLimits[iPt] = lim[iPt];
    for (Int_t iPt = fNPtBins + 1; iPt < kMaxPtBins + 1; iPt++)
      fPtLimits[iPt] = 999.;
  }
  if (fDebug > 1)
  {
    printf("Number of Pt bins = %d\n", fNPtBins);
    for (Int_t iPt = 0; iPt < fNPtBins; iPt++)
      printf(" Bin%d = %8.2f-%8.2f\n", iPt, fPtLimits[iPt], fPtLimits[iPt + 1]);
  }
}
