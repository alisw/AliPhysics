#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGrid.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliAODVZERO.h"
#include "AliAnalysisTaskPi0EtaV2.h"
#include "AliConvEventCuts.h"
#include "AliAODConversionMother.h"
#include "TProfile.h"
#include "AliMultSelection.h"
#include "TProfile2D.h"
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

class AliAnalysisTaskPi0EtaV2;

using namespace std;

ClassImp(AliAnalysisTaskPi0EtaV2)

    //________________________________________________________________________
    AliAnalysisTaskPi0EtaV2::AliAnalysisTaskPi0EtaV2() : AliAnalysisTaskSE(),
                                                         fInputEvent(NULL),
                                                         fnCuts(0),
                                                         fiCut(0),
                                                         fTrainConfig(0),
                                                         fPeriod(""),
                                                         fOutputAODBranchName(""),
                                                         fOutputBGBranchName(""),
                                                         fEventCutArray(NULL),
                                                         fOutputContainer(NULL),
                                                         fCutFolder(NULL),
                                                         fESDList(NULL),
                                                         fQAList(NULL),
                                                         fHistoMotherInvMassPtPhiV0A(NULL),
                                                         fHistoMotherInvMassPtPhiV0C(NULL),
                                                         fHistoMotherBackInvMassPtdPhiV0A(NULL),
                                                         fHistoMotherBackInvMassPtdPhiV0C(NULL),
                                                         fEventCount(NULL),
                                                         fHistoMotherInvMassPtV0CInPlane(NULL),
                                                         fHistoMotherInvMassPtV0AInPlane(NULL),
                                                         fHistoMotherInvMassPtV0COutPlane(NULL),
                                                         fHistoMotherInvMassPtV0AOutPlane(NULL),
                                                         fHistoMotherBackInvMassPtV0CInPlane(NULL),
                                                         fHistoMotherBackInvMassPtV0AInPlane(NULL),
                                                         fHistoMotherBackInvMassPtV0COutPlane(NULL),
                                                         fHistoMotherBackInvMassPtV0AOutPlane(NULL),
                                                         fHistoMotherInvMassPtV0CCos2phi(NULL),
                                                         fHistoMotherInvMassPtV0ACos2phi(NULL),
                                                         fHistoMotherBackInvMassPtV0CCos2phi(NULL),
                                                         fHistoMotherBackInvMassPtV0ACos2phi(NULL),
                                                         runNumList(NULL),
                                                         fhistPhi(NULL),
                                                         fhistPhiBG(NULL),
                                                         fHistRunNumBin(NULL),
                                                         runNum(0),
                                                         oldRunNum(0),
                                                         centSPD1(-999),
                                                         IsVZEROCalibOn(kTRUE),
                                                         IsQAVZERO(kTRUE),
                                                         IsUseInOutPlane(kTRUE),
                                                         fListVZEROCalib(NULL),
                                                         fVZEROCalibFile(NULL),
                                                         fPsi2V0C(-999),
                                                         fPsi2V0A(-999),
                                                         fHist2DPsi2V0CCent(NULL),
                                                         fHist2DPsi2V0ACent(NULL),
                                                         hMultV0(NULL),
                                                         contMult(NULL),
                                                         contQxncm(NULL),
                                                         contQyncm(NULL),
                                                         contQxnam(NULL),
                                                         contQynam(NULL),
                                                         fHCorrectV0ChWeghts(NULL),
                                                         fProfileV0CQxCentGE(NULL),
                                                         fProfileV0CQyCentGE(NULL),
                                                         fProfileV0CQxVtxGE(NULL),
                                                         fProfileV0CQyVtxGE(NULL),
                                                         fHist2CalibPsi2V0CCentGE(NULL),
                                                         fProfileV0AQxCentGE(NULL),
                                                         fProfileV0AQyCentGE(NULL),
                                                         fProfileV0AQxVtxGE(NULL),
                                                         fProfileV0AQyVtxGE(NULL),
                                                         fHist2CalibPsi2V0ACentGE(NULL),
                                                         fProfileV0CQxCentRC(NULL),
                                                         fProfileV0CQyCentRC(NULL),
                                                         fProfileV0CQxVtxRC(NULL),
                                                         fProfileV0CQyVtxRC(NULL),
                                                         fHist2CalibPsi2V0CCentRC(NULL),
                                                         fProfileV0AQxCentRC(NULL),
                                                         fProfileV0AQyCentRC(NULL),
                                                         fProfileV0AQxVtxRC(NULL),
                                                         fProfileV0AQyVtxRC(NULL),
                                                         fHist2CalibPsi2V0ACentRC(NULL),
                                                         fHist2V0Res(NULL)
{
    for (int i = 0; i < 2; ++i)
    {
        hQx2mV0[i] = NULL;
        hQy2mV0[i] = NULL;
    }
}

//________________________________________________________________________
AliAnalysisTaskPi0EtaV2::AliAnalysisTaskPi0EtaV2(const char *name) : AliAnalysisTaskSE(name),
                                                                     fInputEvent(NULL),
                                                                     fnCuts(0),
                                                                     fiCut(0),
                                                                     fTrainConfig(0),
                                                                     fPeriod(""),
                                                                     fOutputAODBranchName(""),
                                                                     fOutputBGBranchName(""),
                                                                     fEventCutArray(NULL),
                                                                     fOutputContainer(NULL),
                                                                     fCutFolder(NULL),
                                                                     fESDList(NULL),
                                                                     fQAList(NULL),
                                                                     fHistoMotherInvMassPtPhiV0A(NULL),
                                                                     fHistoMotherInvMassPtPhiV0C(NULL),
                                                                     fHistoMotherBackInvMassPtdPhiV0A(NULL),
                                                                     fHistoMotherBackInvMassPtdPhiV0C(NULL),
                                                                     fEventCount(NULL),
                                                                     fHistoMotherInvMassPtV0CInPlane(NULL),
                                                                     fHistoMotherInvMassPtV0AInPlane(NULL),
                                                                     fHistoMotherInvMassPtV0COutPlane(NULL),
                                                                     fHistoMotherInvMassPtV0AOutPlane(NULL),
                                                                     fHistoMotherBackInvMassPtV0CInPlane(NULL),
                                                                     fHistoMotherBackInvMassPtV0AInPlane(NULL),
                                                                     fHistoMotherBackInvMassPtV0COutPlane(NULL),
                                                                     fHistoMotherBackInvMassPtV0AOutPlane(NULL),
                                                                     fHistoMotherInvMassPtV0CCos2phi(NULL),
                                                                     fHistoMotherInvMassPtV0ACos2phi(NULL),
                                                                     fHistoMotherBackInvMassPtV0CCos2phi(NULL),
                                                                     fHistoMotherBackInvMassPtV0ACos2phi(NULL),
                                                                     runNumList(NULL),
                                                                     fhistPhi(NULL),
                                                                     fhistPhiBG(NULL),
                                                                     fHistRunNumBin(NULL),
                                                                     runNum(0),
                                                                     oldRunNum(0),
                                                                     centSPD1(-999),
                                                                     IsVZEROCalibOn(kTRUE),
                                                                     IsQAVZERO(kTRUE),
                                                                     IsUseInOutPlane(kTRUE),
                                                                     fListVZEROCalib(NULL),
                                                                     fVZEROCalibFile(NULL),
                                                                     fPsi2V0C(-999),
                                                                     fPsi2V0A(-999),
                                                                     fHist2DPsi2V0CCent(NULL),
                                                                     fHist2DPsi2V0ACent(NULL),
                                                                     hMultV0(NULL),
                                                                     contMult(NULL),
                                                                     contQxncm(NULL),
                                                                     contQyncm(NULL),
                                                                     contQxnam(NULL),
                                                                     contQynam(NULL),
                                                                     fHCorrectV0ChWeghts(NULL),
                                                                     fProfileV0CQxCentGE(NULL),
                                                                     fProfileV0CQyCentGE(NULL),
                                                                     fProfileV0CQxVtxGE(NULL),
                                                                     fProfileV0CQyVtxGE(NULL),
                                                                     fHist2CalibPsi2V0CCentGE(NULL),
                                                                     fProfileV0AQxCentGE(NULL),
                                                                     fProfileV0AQyCentGE(NULL),
                                                                     fProfileV0AQxVtxGE(NULL),
                                                                     fProfileV0AQyVtxGE(NULL),
                                                                     fHist2CalibPsi2V0ACentGE(NULL),
                                                                     fProfileV0CQxCentRC(NULL),
                                                                     fProfileV0CQyCentRC(NULL),
                                                                     fProfileV0CQxVtxRC(NULL),
                                                                     fProfileV0CQyVtxRC(NULL),
                                                                     fHist2CalibPsi2V0CCentRC(NULL),
                                                                     fProfileV0AQxCentRC(NULL),
                                                                     fProfileV0AQyCentRC(NULL),
                                                                     fProfileV0AQxVtxRC(NULL),
                                                                     fProfileV0AQyVtxRC(NULL),
                                                                     fHist2CalibPsi2V0ACentRC(NULL),
                                                                     fHist2V0Res(NULL)
{
    for (int i = 0; i < 2; ++i)
    {
        hQx2mV0[i] = NULL;
        hQy2mV0[i] = NULL;
    }

    DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                     // this chain is created by the analysis manager, so no need to worry about it,
                                     // it does its work automatically
    DefineOutput(1, TList::Class()); // define the ouptut slots
}

//________________________________________________________________________
AliAnalysisTaskPi0EtaV2::~AliAnalysisTaskPi0EtaV2()
{
    if (fOutputContainer)
        delete fOutputContainer;
}
//________________________________________________________________________
void AliAnalysisTaskPi0EtaV2::UserCreateOutputObjects()
{
    // create a new TList that OWNS its objects
    if (fOutputContainer != NULL)
    {
        delete fOutputContainer;
        fOutputContainer = NULL;
    }
    if (fOutputContainer == NULL)
    {
        fOutputContainer = new TList();
        fOutputContainer->SetOwner(kTRUE);
    }
    fCutFolder = new TList *[fnCuts];
    fESDList = new TList *[fnCuts];
    fQAList = new TList *[fnCuts];
    fHistoMotherInvMassPtPhiV0C = new TH3F *[fnCuts];
    fHistoMotherInvMassPtPhiV0A = new TH3F *[fnCuts];
    fHistoMotherBackInvMassPtdPhiV0C = new TH3F *[fnCuts];
    fHistoMotherBackInvMassPtdPhiV0A = new TH3F *[fnCuts];
    fHistoMotherInvMassPtV0CInPlane = new TH2F *[fnCuts];
    fHistoMotherInvMassPtV0AInPlane = new TH2F *[fnCuts];
    fHistoMotherInvMassPtV0COutPlane = new TH2F *[fnCuts];
    fHistoMotherInvMassPtV0AOutPlane = new TH2F *[fnCuts];
    fHistoMotherBackInvMassPtV0CInPlane = new TH2F *[fnCuts];
    fHistoMotherBackInvMassPtV0AInPlane = new TH2F *[fnCuts];
    fHistoMotherBackInvMassPtV0COutPlane = new TH2F *[fnCuts];
    fHistoMotherBackInvMassPtV0AOutPlane = new TH2F *[fnCuts];
    fHistoMotherInvMassPtV0CCos2phi = new TProfile2D *[fnCuts];
    fHistoMotherInvMassPtV0ACos2phi = new TProfile2D *[fnCuts];
    fHistoMotherBackInvMassPtV0CCos2phi = new TProfile2D *[fnCuts];
    fHistoMotherBackInvMassPtV0ACos2phi = new TProfile2D *[fnCuts];
    fEventCount = new TH1D *[fnCuts];
    fHist2DPsi2V0CCent = new TH2D *[fnCuts];
    fHist2DPsi2V0ACent = new TH2D *[fnCuts];
    fProfileV0CQxCentGE = new TProfile *[fnCuts];
    fProfileV0CQyCentGE = new TProfile *[fnCuts];
    fProfileV0CQxVtxGE = new TProfile *[fnCuts];
    fProfileV0CQyVtxGE = new TProfile *[fnCuts];
    fHist2CalibPsi2V0CCentGE = new TH2D *[fnCuts];
    fProfileV0AQxCentGE = new TProfile *[fnCuts];
    fProfileV0AQyCentGE = new TProfile *[fnCuts];
    fProfileV0AQxVtxGE = new TProfile *[fnCuts];
    fProfileV0AQyVtxGE = new TProfile *[fnCuts];
    fHist2CalibPsi2V0ACentGE = new TH2D *[fnCuts];
    fProfileV0CQxCentRC = new TProfile *[fnCuts];
    fProfileV0CQyCentRC = new TProfile *[fnCuts];
    fProfileV0CQxVtxRC = new TProfile *[fnCuts];
    fProfileV0CQyVtxRC = new TProfile *[fnCuts];
    fHist2CalibPsi2V0CCentRC = new TH2D *[fnCuts];
    fProfileV0AQxCentRC = new TProfile *[fnCuts];
    fProfileV0AQyCentRC = new TProfile *[fnCuts];
    fProfileV0AQxVtxRC = new TProfile *[fnCuts];
    fProfileV0AQyVtxRC = new TProfile *[fnCuts];
    fHist2CalibPsi2V0ACentRC = new TH2D *[fnCuts];
    fHist2V0Res = new TProfile *[fnCuts];

    Int_t nBinsPt = 108;
    Double_t minPt = 0.;
    Double_t maxPt = 40;
    Double_t *arrPtBinning = new Double_t[1200];
    Int_t nBinsMinv = 800;
    Double_t arrMinvBin[nBinsMinv + 1];
    Float_t maxMinv = 0.8;
    Double_t minPhi = 0.;
    Double_t maxPhi = TMath::Pi();
    Int_t nPhiBins = 6;
    Double_t arrPhiBin[nPhiBins + 1];
    for (Int_t i = 0; i < nPhiBins + 1; i++) // phiBin[0,pi]
    {
        arrPhiBin[i] = minPhi + ((maxPhi - minPhi) / nPhiBins) * i;
    }
    for (Int_t i = 0; i < nBinsMinv + 1; i++) // MassBin
    {
        arrMinvBin[i] = (maxMinv / nBinsMinv) * i;
    }

    for (Int_t i = 0; i < nBinsPt + 1; i++) // PtBin
    {
        if (i < 1)
            arrPtBinning[i] = 0.5 * i;
        else if (i < 56)
            arrPtBinning[i] = 0.5 + 0.1 * (i - 1);
        else if (i < 80)
            arrPtBinning[i] = 6. + 0.25 * (i - 56);
        else if (i < 108)
            arrPtBinning[i] = 12. + 1.0 * (i - 80);
        else
            arrPtBinning[i] = maxPt;
    }
    for (Int_t iCut = 0; iCut < fnCuts; iCut++)
    {
        fCutFolder[iCut] = new TList();
        fCutFolder[iCut]->SetName(Form("Cent Bin -%i", iCut));
        fCutFolder[iCut]->SetOwner(kTRUE);
        fOutputContainer->Add(fCutFolder[iCut]);
        fESDList[iCut] = new TList();
        fESDList[iCut]->SetName("ESD histograms");
        fESDList[iCut]->SetOwner(kTRUE);
        fQAList[iCut] = new TList();
        fQAList[iCut]->SetName("VZERO QA");
        fQAList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fESDList[iCut]);
        fCutFolder[iCut]->Add(fQAList[iCut]);
        fHistoMotherInvMassPtPhiV0A[iCut] = new TH3F("ESD_Mother_InvMass_Pt_PhiV0A", "ESD_Mother_InvMass_Pt_PhiV0A", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherInvMassPtPhiV0A[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtPhiV0A[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherInvMassPtPhiV0A[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherInvMassPtPhiV0A[iCut]);
        fHistoMotherInvMassPtPhiV0C[iCut] = new TH3F("ESD_Mother_InvMass_Pt_PhiV0C", "ESD_Mother_InvMass_Pt_PhiV0C", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherInvMassPtPhiV0C[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtPhiV0C[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherInvMassPtPhiV0C[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherInvMassPtPhiV0C[iCut]);
        fHistoMotherBackInvMassPtdPhiV0A[iCut] = new TH3F("ESD_Background_InvMass_Pt_PhiV0A", "ESD_Background_InvMass_Pt_PhiV0A", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherBackInvMassPtdPhiV0A[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtdPhiV0A[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherBackInvMassPtdPhiV0A[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtdPhiV0A[iCut]);
        fHistoMotherBackInvMassPtdPhiV0C[iCut] = new TH3F("ESD_Background_InvMass_Pt_PhiV0C", "ESD_Background_InvMass_Pt_PhiV0C", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherBackInvMassPtdPhiV0C[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtdPhiV0C[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherBackInvMassPtdPhiV0C[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtdPhiV0C[iCut]);
        fHistoMotherInvMassPtV0CInPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0C_InPlane", "ESD_Mother_InvMass_Pt_V0C_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0CInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0CInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0CInPlane[iCut]);
        fHistoMotherInvMassPtV0AInPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0A_InPlane", "ESD_Mother_InvMass_Pt_V0A_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0AInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0AInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0AInPlane[iCut]);
        fHistoMotherInvMassPtV0COutPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0C_OutPlane", "ESD_Mother_InvMass_Pt_V0C_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0COutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0COutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0COutPlane[iCut]);
        fHistoMotherInvMassPtV0AOutPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0A_OutPlane", "ESD_Mother_InvMass_Pt_V0A_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0AOutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0AOutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0AOutPlane[iCut]);
        fHistoMotherBackInvMassPtV0CInPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0C_InPlane", "ESD_Background_InvMass_Pt_V0C_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0CInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0CInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0CInPlane[iCut]);
        fHistoMotherBackInvMassPtV0AInPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0A_InPlane", "ESD_Background_InvMass_Pt_V0A_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0AInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0AInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0AInPlane[iCut]);
        fHistoMotherBackInvMassPtV0COutPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0C_OutPlane", "ESD_Background_InvMass_Pt_V0C_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0COutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0COutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0COutPlane[iCut]);
        fHistoMotherBackInvMassPtV0AOutPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0A_OutPlane", "ESD_Background_InvMass_Pt_V0A_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0AOutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0AOutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0AOutPlane[iCut]);
        fHistoMotherInvMassPtV0CCos2phi[iCut] = new TProfile2D("InvMassPtV0CCos2phi", "InvMassPtV0CCos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0ACos2phi[iCut] = new TProfile2D("InvMassPtV0ACos2phi", "InvMassPtV0ACos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0CCos2phi[iCut] = new TProfile2D("BackInvMassPtV0CCos2phi", "BackInvMassPtV0CCos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0ACos2phi[iCut] = new TProfile2D("BackInvMassPtV0ACos2phi", "BackInvMassPtV0ACos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0CCos2phi[iCut]);
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0ACos2phi[iCut]);
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0CCos2phi[iCut]);
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0ACos2phi[iCut]);

        fEventCount[iCut] = new TH1D("EventCount", "EventCount", 100, 0, 100);
        fESDList[iCut]->Add(fEventCount[iCut]);
        // VZERO
        fHist2DPsi2V0CCent[iCut] = new TH2D("fHist2DPsi2V0CCent", "", 20, 0, 100, 100, -2 * TMath::Pi(), 2 * TMath::Pi());
        fHist2DPsi2V0ACent[iCut] = new TH2D("fHist2DPsi2V0ACent", "", 20, 0, 100, 100, -2 * TMath::Pi(), 2 * TMath::Pi());
        fProfileV0CQxCentGE[iCut] = new TProfile("fProfileV0CQxCentGE", "", 100, 0, 100.);
        fProfileV0CQyCentGE[iCut] = new TProfile("fProfileV0CQyCentGE", "", 100, 0, 100.);
        fProfileV0CQxVtxGE[iCut] = new TProfile("fProfileV0CQxVzGE", "", 20, -10, 10);
        fProfileV0CQyVtxGE[iCut] = new TProfile("fProfileV0CQyVzGE", "", 20, -10, 10);
        fHist2CalibPsi2V0CCentGE[iCut] = new TH2D("fHist2CalibPsi2V0CCentGE", "", 20, 0, 100, 50, 0, TMath::Pi());
        fProfileV0AQxCentGE[iCut] = new TProfile("fProfileV0AQxCentGE", "", 100, 0, 100.);
        fProfileV0AQyCentGE[iCut] = new TProfile("fProfileV0AQyCentGE", "", 100, 0, 100.);
        fProfileV0AQxVtxGE[iCut] = new TProfile("fProfileV0AQxVzGE", "", 20, -10, 10);
        fProfileV0AQyVtxGE[iCut] = new TProfile("fProfileV0AQyVzGE", "", 20, -10, 10);
        fHist2CalibPsi2V0ACentGE[iCut] = new TH2D("fHist2CalibPsi2V0ACentGE", "", 20, 0, 100, 50, 0, TMath::Pi());
        fProfileV0CQxCentRC[iCut] = new TProfile("fProfileV0CQxCentRC", "", 100, 0, 100.);
        fProfileV0CQyCentRC[iCut] = new TProfile("fProfileV0CQyCentRC", "", 100, 0, 100.);
        fProfileV0CQxVtxRC[iCut] = new TProfile("fProfileV0CQxVzRC", "", 20, -10, 10);
        fProfileV0CQyVtxRC[iCut] = new TProfile("fProfileV0CQyVzRC", "", 20, -10, 10);
        fHist2CalibPsi2V0CCentRC[iCut] = new TH2D("fHist2CalibPsi2V0CCentRC", "", 20, 0, 100, 50, 0, TMath::Pi());
        fProfileV0AQxCentRC[iCut] = new TProfile("fProfileV0AQxCentRC", "", 100, 0, 100.);
        fProfileV0AQyCentRC[iCut] = new TProfile("fProfileV0AQyCentRC", "", 100, 0, 100.);
        fProfileV0AQxVtxRC[iCut] = new TProfile("fProfileV0AQxVzRC", "", 20, -10, 10);
        fProfileV0AQyVtxRC[iCut] = new TProfile("fProfileV0AQyVzRC", "", 20, -10, 10);
        fHist2CalibPsi2V0ACentRC[iCut] = new TH2D("fHist2CalibPsi2V0ACentRC", "", 20, 0, 100, 50, 0, TMath::Pi());
        fHist2V0Res[iCut] = new TProfile("fHist2V0Res", "", 4, 0, 4);
        fESDList[iCut]->Add(fHist2DPsi2V0CCent[iCut]);
        fESDList[iCut]->Add(fHist2DPsi2V0ACent[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxVtxGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyVtxGE[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0CCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxVtxGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyVtxGE[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0ACentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxVtxRC[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyVtxRC[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0CCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxVtxRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyVtxRC[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0ACentRC[iCut]);
        fESDList[iCut]->Add(fHist2V0Res[iCut]);
    }
    ////runNumber hist
    ////runNumber for 18qr
    TString runNumListString[214] = {"296623", "296622", "296621", "296619", "296618", "296616", "296615", "296594", "296553", "296552",
                                     "296551", "296550", "296548", "296547", "296516", "296512", "296511", "296510", "296509", "296472",
                                     "296433", "296424", "296423", "296420", "296419", "296415", "296414", "296383", "296381", "296380",
                                     "296379", "296378", "296377", "296376", "296375", "296312", "296309", "296304", "296303", "296280",
                                     "296279", "296273", "296270", "296269", "296247", "296246", "296244", "296243", "296242", "296241",
                                     "296240", "296198", "296197", "296196", "296195", "296194", "296192", "296191", "296143", "296142",
                                     "296135", "296134", "296133", "296132", "296123", "296074", "296066", "296065", "296063", "296062",
                                     "296060", "296016", "295942", "295941", "295937", "295936", "295913", "295910", "295909", "295861",
                                     "295860", "295859", "295856", "295855", "295854", "295853", "295831", "295829", "295826", "295825",
                                     "295822", "295819", "295818", "295816", "295791", "295788", "295786", "295763", "295762", "295759",
                                     "295758", "295755", "295754", "295725", "295723", "295721", "295719", "295718", "295717", "295714",
                                     "295712", "295676", "295675", "295673", "295668", "295667", "295666", "295615", "295612", "295611",
                                     "295610", "295589", "295588", "295586", "295585", "297595", "297590", "297588", "297558", "297544", "297542", "297541", "297540", "297537", "297512",
                                     "297483", "297479", "297452", "297451", "297450", "297446", "297442", "297441", "297415", "297414",
                                     "297413", "297406", "297405", "297380", "297379", "297372", "297367", "297366", "297363", "297336",
                                     "297335", "297333", "297332", "297317", "297311", "297310", "297278", "297222", "297221", "297218",
                                     "297196", "297195", "297193", "297133", "297132", "297129", "297128", "297124", "297123", "297119",
                                     "297118", "297117", "297085", "297035", "297031", "296966", "296941", "296938", "296935", "296934",
                                     "296932", "296931", "296930", "296903", "296900", "296899", "296894", "296852", "296851", "296850",
                                     "296848", "296839", "296838", "296836", "296835", "296799", "296794", "296793", "296790", "296787",
                                     "296786", "296785", "296784", "296781", "296752", "296694", "296693", "296691", "296690"};

    runNumList = new std::map<int, int>;
    for (int i = 0; i < 214; i++)
        runNumList->insert(std::pair<int, int>(runNumListString[i].Atoi(), i + 1));
    fHistRunNumBin = new TH1I("runNumBin", "", (int)runNumList->size(), 1, (int)runNumList->size() + 1);
    std::map<int, int>::iterator iter;
    for (auto runNum1 : *runNumList)
        fHistRunNumBin->GetXaxis()->SetBinLabel(runNum1.second, Form("%i", runNum1.first));
    fQAList[0]->Add(fHistRunNumBin);

    fhistPhi = new TH1F *[214];
    fhistPhiBG = new TH1F *[214];
    for (std::map<int, int>::reverse_iterator r_it = runNumList->rbegin(); r_it != runNumList->rend(); ++r_it)
    {
        fhistPhi[r_it->second - 1] = new TH1F(Form("fhistPhi%i", r_it->first), Form("fhistPhi%i", r_it->first), 200, -TMath::TwoPi(), TMath::TwoPi());
        fhistPhiBG[r_it->second - 1] = new TH1F(Form("fhistPhi_BG%i", r_it->first), Form("fhistPhi_BG%i", r_it->first), 200, -TMath::TwoPi(), TMath::TwoPi());
        fQAList[0]->Add(fhistPhi[r_it->second - 1]);
        fQAList[0]->Add(fhistPhiBG[r_it->second - 1]);
    }
    ////////////////////////
    // VZERO
    ////////////////////////
    if (!gGrid)
    {
        TGrid::Connect("alien://");
    }

    if (IsVZEROCalibOn)
    {
        if (fPeriod.EqualTo("LHC15o"))
        {
            fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/VZEROCalibFile15o.root", "READ");
            fListVZEROCalib = dynamic_cast<TList *>(fVZEROCalibFile->Get("VZEROCalibList"));
            if (fListVZEROCalib)
            {
                // V0C Qx Mean
                // V0C Qy Mean
                // V0A Qx Mean
                // V0A Qy Mean
                contQxncm = (AliOADBContainer *)fListVZEROCalib->FindObject(Form("fqxc%im", 2));
                contQyncm = (AliOADBContainer *)fListVZEROCalib->FindObject(Form("fqyc%im", 2));
                contQxnam = (AliOADBContainer *)fListVZEROCalib->FindObject(Form("fqxa%im", 2));
                contQynam = (AliOADBContainer *)fListVZEROCalib->FindObject(Form("fqya%im", 2));
                contMult = (AliOADBContainer *)fListVZEROCalib->FindObject("hMultV0BefCorPfpx");
            }
            else
                std::cout << "!!!!!!!!!!!!!!!VZERO List not Found!!!!!!!!!!!!!!!" << std::endl;
        }
        if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r"))
        {
            fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/j/jwan/CalibFile/calibq2V0C18qrP3.root", "READ");
            if (fVZEROCalibFile)
            {
                // V0C Qx Mean
                // V0C Qy Mean
                // V0A Qx Mean
                // V0A Qy Mean
                contQxncm = (AliOADBContainer *)fVZEROCalibFile->GetObjectChecked("fqxc2m", "AliOADBContainer");
                contQyncm = (AliOADBContainer *)fVZEROCalibFile->GetObjectChecked("fqyc2m", "AliOADBContainer");
                contQxnam = (AliOADBContainer *)fVZEROCalibFile->GetObjectChecked("fqxa2m", "AliOADBContainer");
                contQynam = (AliOADBContainer *)fVZEROCalibFile->GetObjectChecked("fqya2m", "AliOADBContainer");
            }
            else
                std::cout << "!!!!!!!!!!!!!!!VZERO File not Found!!!!!!!!!!!!!!!" << std::endl;
        }
    }

    // add the list to our output file
    PostData(1, fOutputContainer); //
}
//________________________________________________________________________
void AliAnalysisTaskPi0EtaV2::UserExec(Option_t *)
{

    for (Int_t iCut = 0; iCut < fnCuts; iCut++)
    {
        fiCut = iCut;
        fInputEvent = InputEvent();
        // AliMultSelection *fMultSel = (AliMultSelection *)InputEvent()->FindListObject("MultSelection");
        // centSPD1 = fMultSel->GetMultiplicityPercentile("CL1");
        centSPD1 = ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetCentrality(fInputEvent);
        int centBin = 999; //  depend on GammaCalo trainconfig
        if (centSPD1 >= 0 && centSPD1 < 10)
        {
            centBin = 0;
        }
        else if (centSPD1 >= 10 && centSPD1 < 30)
        {
            centBin = 1;
        }
        else if (centSPD1 >= 30 && centSPD1 < 50)
        {
            centBin = 2;
        }
        else if (centSPD1 >= 50 && centSPD1 < 90)
        {
            centBin = 3;
        }
        if (iCut != centBin)
            continue;
        fOutputAODBranchName = Form("pi0Array_%i", fTrainConfig);
        fOutputBGBranchName = Form("pi0BackgroundArray_%i", fTrainConfig);
        //    cout << "fOutputAODBranchName===" << fOutputAODBranchName.Data() << endl;
        TClonesArray *arrClustersPi0 = NULL;
        TClonesArray *arrClustersBg = NULL;
        arrClustersPi0 = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(fOutputAODBranchName.Data()));
        arrClustersBg = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(fOutputBGBranchName.Data()));
        if (!arrClustersPi0 || !arrClustersBg)
        {
            cout << "No Clusters found in this event" << endl;
            return;
        }
        int nclusPi0 = arrClustersPi0->GetEntriesFast();
        int nclusBg = arrClustersBg->GetEntriesFast();
        //    if (nclusPi0 == 0)
        //        return;
        //    cout << "nclusPi0==" << nclusPi0 << endl;
        // VZERO Plane
        runNum = fInputEvent->GetRunNumber();
        if (IsVZEROCalibOn)
        {
            if (runNum != oldRunNum)
            {
                if (!LoadCalibHistForThisRun())
                {
                    //            cout << "=====faild to LoadCalibHistForThisRun=====" << endl;
                    return;
                }
                oldRunNum = runNum;
            }
            if (!GetVZEROPlane())
            {
                //         cout << "==============faild to GetVZEROPlane===========" << endl;
                return;
            }
        }
        fHist2V0Res[iCut]->Fill(centBin + 0.5, cos(2 * (fPsi2V0A - fPsi2V0C)));
        fHist2DPsi2V0CCent[iCut]->Fill(centSPD1, fPsi2V0C);
        fHist2DPsi2V0ACent[iCut]->Fill(centSPD1, fPsi2V0A);
        fEventCount[iCut]->Fill(centSPD1);
        if (runNumList->find(runNum) == runNumList->end())
            return;
        int rumNumbin = runNumList->at(runNum) - 1;
        fHistRunNumBin->Fill(rumNumbin + 1);
        // Same Event
        for (Long_t i = 0; i < nclusPi0; i++)
        {
            AliAODConversionMother *pi0cand = NULL;
            pi0cand = new AliAODConversionMother(*(AliAODConversionMother *)arrClustersPi0->At(i));
            fHistoMotherInvMassPtV0CCos2phi[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), cos(2 * (pi0cand->Phi() - fPsi2V0C)));
            fHistoMotherInvMassPtV0ACos2phi[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), cos(2 * (pi0cand->Phi() - fPsi2V0A)));
            if (pi0cand->M() > 0.1 && pi0cand->M() < 0.15)
            {
                fhistPhi[rumNumbin]->Fill(pi0cand->Phi());
            }
            else
            {
                fhistPhiBG[rumNumbin]->Fill(pi0cand->Phi());
            }
            double dphiV0A = TVector2::Phi_0_2pi(pi0cand->Phi() - fPsi2V0A);
            if (dphiV0A > TMath::Pi())
            {
                dphiV0A -= TMath::Pi();
            }
            double dphiV0C = TVector2::Phi_0_2pi(pi0cand->Phi() - fPsi2V0C);
            if (dphiV0C > TMath::Pi())
            {
                dphiV0C -= TMath::Pi();
            }
            if (!IsUseInOutPlane)
            {
                fHistoMotherInvMassPtPhiV0C[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), dphiV0C, pi0cand->GetWeight());
                fHistoMotherInvMassPtPhiV0A[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), dphiV0A, pi0cand->GetWeight());
            }
            else
            {
                if (dphiV0A < 0.25 * TMath::Pi() || (dphiV0A > 0.75 * TMath::Pi() && dphiV0A < 1.25 * TMath::Pi()) || dphiV0A > 1.75 * TMath::Pi())
                {
                    fHistoMotherInvMassPtV0AInPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
                }
                else
                {
                    fHistoMotherInvMassPtV0AOutPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
                }

                if (dphiV0C < 0.25 * TMath::Pi() || (dphiV0C > 0.75 * TMath::Pi() && dphiV0C < 1.25 * TMath::Pi()) || dphiV0C > 1.75 * TMath::Pi())
                {
                    fHistoMotherInvMassPtV0CInPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
                }
                else
                {
                    fHistoMotherInvMassPtV0COutPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
                }
            }
            delete pi0cand;
            pi0cand = 0x0;
        }
        // Background
        for (Long_t i = 0; i < nclusBg; i++)
        {
            AliAODConversionMother *Bgcand = NULL;
            Bgcand = new AliAODConversionMother(*(AliAODConversionMother *)arrClustersBg->At(i));
            // Use cos(2*(psi-Psi))
            fHistoMotherBackInvMassPtV0CCos2phi[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), cos(2 * (Bgcand->Phi() - fPsi2V0C)));
            fHistoMotherBackInvMassPtV0ACos2phi[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), cos(2 * (Bgcand->Phi() - fPsi2V0A)));
            //  fhistPhiBG[fiCut]->Fill(Bgcand->Phi());
            double dphiV0A = TVector2::Phi_0_2pi(Bgcand->Phi() - fPsi2V0A);
            if (dphiV0A > TMath::Pi())
            {
                dphiV0A -= TMath::Pi();
            }
            double dphiV0C = TVector2::Phi_0_2pi(Bgcand->Phi() - fPsi2V0C);
            if (dphiV0C > TMath::Pi())
            {
                dphiV0C -= TMath::Pi();
            }
            if (!IsUseInOutPlane)
            {
                fHistoMotherBackInvMassPtdPhiV0C[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), dphiV0C, Bgcand->GetWeight());
                fHistoMotherBackInvMassPtdPhiV0A[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), dphiV0A, Bgcand->GetWeight());
            }
            else
            {
                if (dphiV0A < 0.25 * TMath::Pi() || (dphiV0A > 0.75 * TMath::Pi() && dphiV0A < 1.25 * TMath::Pi()) || dphiV0A > 1.75 * TMath::Pi())
                {
                    fHistoMotherBackInvMassPtV0AInPlane[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), Bgcand->GetWeight());
                }
                else
                {
                    fHistoMotherBackInvMassPtV0AOutPlane[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), Bgcand->GetWeight());
                }

                if (dphiV0C < 0.25 * TMath::Pi() || (dphiV0C > 0.75 * TMath::Pi() && dphiV0C < 1.25 * TMath::Pi()) || dphiV0C > 1.75 * TMath::Pi())
                {
                    fHistoMotherBackInvMassPtV0CInPlane[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), Bgcand->GetWeight());
                }
                else
                {
                    fHistoMotherBackInvMassPtV0COutPlane[fiCut]->Fill(Bgcand->M(), Bgcand->Pt(), Bgcand->GetWeight());
                }
            }

            delete Bgcand;
            Bgcand = 0x0;
        }
    }
    PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskPi0EtaV2::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}

//________________________________________________________________________
bool AliAnalysisTaskPi0EtaV2::LoadCalibHistForThisRun()
{
    if (fPeriod.EqualTo("LHC15o"))
    {
        if (IsVZEROCalibOn)
        {
            //  hMultV0->Reset();
            //  for (int i = 0; i < 2; i++)
            //  {
            //      hQx2mV0[i]->Reset();
            //      hQy2mV0[i]->Reset();
            //  }
            if (!contQxncm || !contQyncm || !contQxnam || !contQynam)
                return false;
            hMultV0 = ((TH1D *)contMult->GetObject(runNum));
            hQx2mV0[0] = ((TH1D *)contQxncm->GetObject(runNum));
            hQy2mV0[0] = ((TH1D *)contQyncm->GetObject(runNum));
            hQx2mV0[1] = ((TH1D *)contQxnam->GetObject(runNum));
            hQy2mV0[1] = ((TH1D *)contQynam->GetObject(runNum));
            if (!hMultV0)
                return false;
            if (!hQx2mV0[0])
                return false;
            if (!hQy2mV0[0])
                return false;
            if (!hQx2mV0[1])
                return false;
            if (!hQy2mV0[1])
                return false;
        }
    }
    if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r"))
    {
        //   for (int i = 0; i < 2; i++)
        //   {
        //       hQx2mV0[i]->Reset();
        //       hQy2mV0[i]->Reset();
        //   }
        if (!contQxncm || !contQyncm || !contQxnam || !contQynam)
        {
            cout << "contQncm not found" << endl;
            return false;
        }
        //  hQx2mV0[0] = ((TH1D *)contQxncm->GetObject(runNum));
        //  hQy2mV0[0] = ((TH1D *)contQyncm->GetObject(runNum));
        //  hQx2mV0[1] = ((TH1D *)contQxnam->GetObject(runNum));
        //  hQy2mV0[1] = ((TH1D *)contQynam->GetObject(runNum));

        hQx2mV0[0] = ((TH1D *)contQxncm->GetDefaultObject(Form("hV0QxMeanCRun%d", runNum)));
        hQy2mV0[0] = ((TH1D *)contQyncm->GetDefaultObject(Form("hV0QyMeanCRun%d", runNum)));
        hQx2mV0[1] = ((TH1D *)contQxnam->GetDefaultObject(Form("hV0QxMeanARun%d", runNum)));
        hQy2mV0[1] = ((TH1D *)contQynam->GetDefaultObject(Form("hV0QyMeanARun%d", runNum)));
        if (!hQx2mV0[0] || !hQy2mV0[0] || !hQx2mV0[1] || !hQy2mV0[1])
        {
            //     cout << "=======hQ2mV0 not found======" << endl;
            return false;
        }
        //    fHCorrectV0ChWeghts->Reset();
        fHCorrectV0ChWeghts = (TH2F *)fVZEROCalibFile->GetObjectChecked(Form("hWgtV0ChannelsvsVzRun%d", runNum), "TH2F");
        if (!fHCorrectV0ChWeghts)
        {
            //       cout << "=======fHCorrectV0ChWeghts not found======" << endl;
            return false;
        }
    }
    return true;
}

//________________________________________________________________________
bool AliAnalysisTaskPi0EtaV2::GetVZEROPlane()
{
    double multV0Ch[64] = {0};
    double V0XMean[3] = {0};
    double V0YMean[3] = {0};

    // [0]: M; [1]: C; [2]: A;
    double qxGE[3] = {0}, qyGE[3] = {0};
    double qxRC[3] = {0}, qyRC[3] = {0};
    double multRingGE[3] = {0};
    double psi2GE[3] = {0};
    double psi2RC[3] = {0};

    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    if (fPeriod.EqualTo("LHC15o"))
        for (int iCh = 0; iCh < 64; ++iCh)
            multV0Ch[iCh] = hMultV0->GetBinContent(iCh + 1);
    int iCentSPD = (int)centSPD1;
    if (iCentSPD >= 90)
        return false;
    V0XMean[0] = -999.;
    V0YMean[0] = -999.;
    for (int i = 0; i < 2; ++i)
    { // [1]: C; [2]: A;
        V0XMean[i + 1] = hQx2mV0[i]->GetBinContent(iCentSPD + 1);
        V0YMean[i + 1] = hQy2mV0[i]->GetBinContent(iCentSPD + 1);
    }

    for (int iCh = 0; iCh < 64; ++iCh)
    {
        double phi = TMath::Pi() / 8. + TMath::Pi() / 4. * (iCh % 8);
        double multCh = 0.;
        AliAODVZERO *aodV0 = (AliAODVZERO *)fInputEvent->GetVZEROData();
        multCh = aodV0->GetMultiplicity(iCh);
        if (iCh < 32)
        { // C
            double multChGEC = -1;
            if (fPeriod.EqualTo("LHC15o"))
            {
                if (iCh < 8)
                    multChGEC = multCh / multV0Ch[iCh] * multV0Ch[0];
                else if (iCh >= 8 && iCh < 16)
                    multChGEC = multCh / multV0Ch[iCh] * multV0Ch[8];
                else if (iCh >= 16 && iCh < 24)
                    multChGEC = multCh / multV0Ch[iCh] * multV0Ch[16];
                else if (iCh >= 24 && iCh < 32)
                    multChGEC = multCh / multV0Ch[iCh] * multV0Ch[24];
            }
            if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r"))
            {
                int ibinV0 = fHCorrectV0ChWeghts->FindBin(vertex[2], iCh);
                double V0chGE = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
                multChGEC = multCh * V0chGE;
            }
            if (multChGEC < 0)
                continue;

            // for V0C GE
            qxGE[1] += multChGEC * TMath::Cos(2 * phi);
            qyGE[1] += multChGEC * TMath::Sin(2 * phi);
            multRingGE[1] += multChGEC;
        }
        else if (iCh >= 32 && iCh < 64)
        { // A
            double multChGEA = -1;
            if (fPeriod.EqualTo("LHC15o"))
            {
                if (iCh >= 32 && iCh < 40)
                    multChGEA = multCh / multV0Ch[iCh] * multV0Ch[32];
                else if (iCh >= 40 && iCh < 48)
                    multChGEA = multCh / multV0Ch[iCh] * multV0Ch[40];
                else if (iCh >= 48 && iCh < 56)
                    multChGEA = multCh / multV0Ch[iCh] * multV0Ch[48];
                else if (iCh >= 56 && iCh < 64)
                    multChGEA = multCh / multV0Ch[iCh] * multV0Ch[56];
            }
            if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r"))
            {
                int ibinV0 = fHCorrectV0ChWeghts->FindBin(vertex[2], iCh);
                double V0chGE = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
                multChGEA = multCh * V0chGE;
            }
            if (multChGEA < 0)
                continue;
            // for V0A GE
            qxGE[2] += multChGEA * TMath::Cos(2 * phi);
            qyGE[2] += multChGEA * TMath::Sin(2 * phi);
            multRingGE[2] += multChGEA;
        }
    }
    if (multRingGE[1] < 1.e-6 || multRingGE[2] < 1.e-6)
    {
        return false;
    }

    // VZERO GE Plane
    for (int i = 1; i < 3; i++)
    {
        psi2GE[i] = GetEventPlane(qxGE[i], qyGE[i], 2.);
        if (TMath::IsNaN(psi2GE[i]))
        {
            return false;
        }
    }
    // VZERO Recenter
    for (int i = 1; i < 3; ++i)
    {
        double qxMean = V0XMean[i];
        double qyMean = V0YMean[i];
        if (TMath::IsNaN(qxMean) || TMath::IsNaN(qyMean))
            continue;
        if (qyMean < -900 || qxMean < -900)
            continue;
        // For 10 h, we've stored the qx/y of V0M, and they cannot been found in A.Dorbin's calib file for 15o period!
        qxRC[i] = qxGE[i] - qxMean;
        qyRC[i] = qyGE[i] - qyMean;
        psi2RC[i] = GetEventPlane(qxRC[i], qyRC[i], 2.);
        if (TMath::IsNaN(psi2RC[i]))
        {
            return false;
        }
    }
    // VZERO QA //0GE 1RC
    if (IsQAVZERO)
    {
        // V0C
        fProfileV0CQxCentGE[fiCut]->Fill(centSPD1, qxGE[1]);
        fProfileV0CQyCentGE[fiCut]->Fill(centSPD1, qyGE[1]);
        fProfileV0CQxVtxGE[fiCut]->Fill(vertex[2], qxGE[1]);
        fProfileV0CQyVtxGE[fiCut]->Fill(vertex[2], qyGE[1]);
        fHist2CalibPsi2V0CCentGE[fiCut]->Fill(centSPD1, psi2GE[1]);

        fProfileV0CQxCentRC[fiCut]->Fill(centSPD1, qxRC[1]);
        fProfileV0CQyCentRC[fiCut]->Fill(centSPD1, qyRC[1]);
        fProfileV0CQxVtxRC[fiCut]->Fill(vertex[2], qxRC[1]);
        fProfileV0CQyVtxRC[fiCut]->Fill(vertex[2], qyRC[1]);
        fHist2CalibPsi2V0CCentRC[fiCut]->Fill(centSPD1, psi2RC[1]);
        // V0A
        fProfileV0AQxCentGE[fiCut]->Fill(centSPD1, qxGE[2]);
        fProfileV0AQyCentGE[fiCut]->Fill(centSPD1, qyGE[2]);
        fProfileV0AQxVtxGE[fiCut]->Fill(vertex[2], qxGE[2]);
        fProfileV0AQyVtxGE[fiCut]->Fill(vertex[2], qyGE[2]);
        fHist2CalibPsi2V0ACentGE[fiCut]->Fill(centSPD1, psi2GE[2]);

        fProfileV0AQxCentRC[fiCut]->Fill(centSPD1, qxRC[2]);
        fProfileV0AQyCentRC[fiCut]->Fill(centSPD1, qyRC[2]);
        fProfileV0AQxVtxRC[fiCut]->Fill(vertex[2], qxRC[2]);
        fProfileV0AQyVtxRC[fiCut]->Fill(vertex[2], qyRC[2]);
        fHist2CalibPsi2V0ACentRC[fiCut]->Fill(centSPD1, psi2RC[2]);
    }
    fPsi2V0C = psi2RC[1];
    fPsi2V0A = psi2RC[2];
    return true;
}
//________________________________________________________________________
double AliAnalysisTaskPi0EtaV2::GetEventPlane(double qx, double qy, double harmonic)
{
    double psi = (1. / harmonic) * TMath::ATan2(qy, qx);
    if (psi < 0)
        return psi += TMath::TwoPi() / harmonic;
    else
        return psi;
}
