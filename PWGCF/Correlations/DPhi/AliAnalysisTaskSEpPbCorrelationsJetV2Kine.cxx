#include <TH1D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TParticle.h>
#include <TList.h>

#include "iostream"
#include "AliLog.h"
#include "AliVParticle.h"
#include "AliTHn.h"

#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventPoolManager.h"

#include "AliAnalysisTaskSEpPbCorrelationsJetV2Kine.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSEpPbCorrelationsMCLEGOYS.h"
#include "AliAnalysisTaskSEpPbCorrelationsJetV2.h"

using namespace std;

ClassImp(AliAnalysisTaskSEpPbCorrelationsJetV2Kine)
ClassImp(AliAssociatedTrackYSLEGOMC)
ClassImp(AliTrigAssoPairST)

//_____________________________________________________________________________
AliAnalysisTaskSEpPbCorrelationsJetV2Kine::AliAnalysisTaskSEpPbCorrelationsJetV2Kine() :
AliAnalysisTaskSE(),
fListOutput(nullptr),
fPoolMgr1(0x0),
fPoolMgr2(0x0),
fCen1(0.),
fCen2(0.),
fPtMin(0.2),
fPtMax(10.),
fMode("TPCTPC"),
fAsscoptCut(0.5),
fCentrality(0.),
fNbinsAssocPt(1),
//fNbinsZvtx(1),
fNbinsPtTrig(1),
fPtTrigAxis(0x0),
fVtxZ(0x0),
fPtAssocAxis(0x0),
fHistTPCTPC_SS(0),
fHistTPCTPC_Mixed_SS(0),
fHistTPCTrig(0),
fHistTPCTPCFMDA(0),
fHistTPCTPCFMDA_Mixed(0),
fHistTPCTPCFMDATrig(0),
fHistTPCFMDA(0),
fHistTPCFMDA_Mixed(0),
fHistTPCFMDC(0),
fHistTPCFMDC_Mixed(0),
fHistFMDAFMDC(0),
fHistFMDAFMDC_Mixed(0),
fHistFMDTrig(0),
fT(1.)
{
//
//  Default constructor
//
}

//_____________________________________________________________________________
AliAnalysisTaskSEpPbCorrelationsJetV2Kine::AliAnalysisTaskSEpPbCorrelationsJetV2Kine(const char *name) :
AliAnalysisTaskSE(name),
fListOutput(nullptr),
fPoolMgr1(0x0),
fPoolMgr2(0x0),
fCen1(0.),
fCen2(0.),
fPtMin(0.2),
fPtMax(10.),
fMode("TPCTPC"),
fAsscoptCut(0.5),
fCentrality(0.),
fNbinsAssocPt(1),
//fNbinsZvtx(1),
fNbinsPtTrig(1),
fPtTrigAxis(0x0),
fPtAssocAxis(0x0),
fVtxZ(0x0),
fHistTPCTPC_SS(0),
fHistTPCTPC_Mixed_SS(0),
fHistTPCTrig(0),
fHistTPCTPCFMDA(0),
fHistTPCTPCFMDA_Mixed(0),
fHistTPCTPCFMDATrig(0),
fHistTPCFMDA(0),
fHistTPCFMDA_Mixed(0),
fHistTPCFMDC(0),
fHistTPCFMDC_Mixed(0),
fHistFMDAFMDC(0),
fHistFMDAFMDC_Mixed(0),
fHistFMDTrig(0),
fT(1.)
{
//
//  Constructor
//
  DefineInput(1, TH1D::Class());
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEpPbCorrelationsJetV2Kine::~AliAnalysisTaskSEpPbCorrelationsJetV2Kine()
{
//
//  Default destructor
//
  if (fListOutput) { delete fListOutput; fListOutput = nullptr; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::UserCreateOutputObjects()
{
//
//  void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::UserCreateOutputObjects()
//

  if (fListOutput) {
    delete fListOutput;
    fListOutput = nullptr;
  }

  fListOutput = new TList();
  fListOutput->SetOwner(kTRUE);
  fListOutput->Add(new TH1D("hTrials", "", 5, 0.5, 5.5));
  fListOutput->Add(new TH1D("hImp", "", 200, 0., 20.));
  fListOutput->Add(new TH1D("hRP",  "", 360, 0., TMath::TwoPi()));

  fListOutput->Add(new TH2D("hMu",  "", 200, 0., 50., 150,-4.,-2.5));
  fListOutput->Add(new TH2D("hMu_All",  "", 200, 0., 50., 1000,-10.,10));
  fListOutput->Add(new TH2D("hImp_Charge",  "", 1000, 0., 1000.,300, 0., 30.));
  fListOutput->Add(new TH1D("hChargeV0A",  "", 1000, 0., 1000.));
  fListOutput->Add(new TH1D("hCentrality",  "", 100, 0., 100.));
  fListOutput->Add(new TH2D("hCentrality_MidCharged",  "", 100, 0., 100., 1000, 0., 1000.));
  
  fListOutput->Add(new TH2D("hEta_Phi",  "", 80, -2.5, 5.5, 72, 0 ,2.*TMath::Pi()));

  fListOutput->Add(new TH2D("hPt_Cen","", 40, 0., 20., 100, 0., 100.));
  fListOutput->Add(new TH2D("hPt_Cen_mid","", 40, 0., 20., 100, 0., 100.));
 /*
  TObject *p(nullptr);
  TIter next(fListOutput);
  while ((p = next())) {
    auto h(dynamic_cast<TH1*>(p));
    if (!h) continue;
    if (h->InheritsFrom("TProfile")) continue;
    h->Sumw2();
  }
 */
//=============================== For TPC-TPC
 const Int_t nbins_dPhidEtaPt[] = {fNbinsPtTrig, fNbinsAssocPt, 32, 72};
 const Int_t nVar = sizeof(nbins_dPhidEtaPt) / sizeof(Int_t);

 const TArrayD *aBin_Trig(fPtTrigAxis->GetXbins());
 const Int_t nBin_Trig(aBin_Trig->GetSize() - 1);
 Double_t dBin_Trig[nBin_Trig+1];

 const TArrayD *aBin_Asso(fPtAssocAxis->GetXbins());
 const Int_t nBin_Asso(aBin_Asso->GetSize() - 1);
 Double_t dBin_Asso[nBin_Asso+1];

 for (Int_t i=0; i<=nBin_Trig; ++i) dBin_Trig[i] = (*aBin_Trig)[i];
 for (Int_t i=0; i<=nBin_Asso; ++i) dBin_Asso[i] = (*aBin_Asso)[i];
 
 fHistTPCTPC_SS = new AliTHn("fHistTPCTPC_SS", "fHistTPCTPC_SS", 1, nVar, nbins_dPhidEtaPt);

 fHistTPCTPC_SS->SetBinLimits(0, dBin_Trig);
 fHistTPCTPC_SS->SetBinLimits(1, dBin_Asso);
 fHistTPCTPC_SS->SetBinLimits(2, -1.6, 1.6);
 fHistTPCTPC_SS->SetBinLimits(3, -0.5*TMath::Pi(), 1.5*TMath::Pi());

 fHistTPCTPC_SS->SetVarTitle(0, "leading p_{T} GeV/c");
 fHistTPCTPC_SS->SetVarTitle(1, "asso p_{T} GeV/c");
 fHistTPCTPC_SS->SetVarTitle(2, "#Delta#eta");
 fHistTPCTPC_SS->SetVarTitle(3, "#Delta#phi");
 
 fHistTPCTPC_Mixed_SS = new AliTHn("fHistTPCTPC_Mixed_SS", "fHistTPCTPC_Mixed_SS", 1, nVar, nbins_dPhidEtaPt);

 fHistTPCTPC_Mixed_SS->SetBinLimits(0, dBin_Trig);
 fHistTPCTPC_Mixed_SS->SetBinLimits(1, dBin_Asso);
 fHistTPCTPC_Mixed_SS->SetBinLimits(2, -1.6, 1.6);
 fHistTPCTPC_Mixed_SS->SetBinLimits(3, -0.5*TMath::Pi(), 1.5*TMath::Pi());

 fHistTPCTPC_Mixed_SS->SetVarTitle(0, "leading p_{T} GeV/c");
 fHistTPCTPC_Mixed_SS->SetVarTitle(1, "asso p_{T} GeV/c");
 fHistTPCTPC_Mixed_SS->SetVarTitle(2, "#Delta#eta");
 fHistTPCTPC_Mixed_SS->SetVarTitle(3, "#Delta#phi");

 Double_t dCen[] = {0.,10.,60.,80.,100.};
 Int_t nCen = sizeof(dCen) / sizeof(Double_t) - 1;

 //const Int_t nbins_dTrig[] = {fNbinsPtTrig, nCen};
 const Int_t nbins_dTrig[] = {fNbinsPtTrig};
 const Int_t nVar_Trig = sizeof(nbins_dTrig) / sizeof(Int_t);

 fHistTPCTrig = new AliTHn("fHistTPCTrig", "fHistTPCTrig", 1, nVar_Trig, nbins_dTrig);
 fHistTPCTrig->SetBinLimits(0, dBin_Trig);
 //fHistTPCTrig->SetBinLimits(1, dCen);

 fHistTPCTrig->SetVarTitle(0, "leading p_{T} GeV/c");
// fHistTPCTrig->SetVarTitle(1, "Centrality");

 if(fMode == "TPCTPC" || fMode == "TPCTPCFMDA") 
 {
  fListOutput->Add(fHistTPCTPC_SS);
  fListOutput->Add(fHistTPCTPC_Mixed_SS);
  fListOutput->Add(fHistTPCTrig);
 }

//=============================== For TPC-FMDA, TPC-FMDC, FMDA-FMDC
 const Double_t binning_detaTPCFMDA[] = {-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8};
 const Int_t ndetatpcfmda= sizeof(binning_detaTPCFMDA)/sizeof(Double_t) - 1;

 const  Double_t binning_detaTPCFMDC[] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4};
 const Int_t ndetatpcfmdc= sizeof(binning_detaTPCFMDC)/sizeof(Double_t) - 1;

 const Double_t binning_detaFMDAFMDC[] = {3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4};
 const Int_t ndetafmdafmdc= sizeof(binning_detaFMDAFMDC)/sizeof(Double_t) - 1;

 //const Int_t nbins_dTPCFMDA[] = {fNbinsPtTrig, ndetatpcfmda, 72, nCen};
 const Int_t nbins_dTPCFMDA[] = {fNbinsPtTrig, ndetatpcfmda, 72};
 const Int_t nVar_TPCFMDA  = sizeof(nbins_dTPCFMDA)  / sizeof(Int_t); 

 //const Int_t nbins_dTPCFMDC[] = {fNbinsPtTrig, ndetatpcfmdc, 72, nCen};
 const Int_t nbins_dTPCFMDC[] = {fNbinsPtTrig, ndetatpcfmdc, 72};
 const Int_t nVar_TPCFMDC  = sizeof(nbins_dTPCFMDC)  / sizeof(Int_t);

 //const Int_t nbins_dFMDAFMDC[] = {ndetafmdafmdc, 72, nCen};
 const Int_t nbins_dFMDAFMDC[] = {ndetafmdafmdc, 72};
 const Int_t nVar_FMDAFMDC = sizeof(nbins_dFMDAFMDC) / sizeof(Int_t);

 fHistTPCFMDA  = new AliTHn("fHistTPCFMDA", "fHistTPCFMDA",  1, nVar_TPCFMDA , nbins_dTPCFMDA);
 fHistTPCFMDA->SetBinLimits(0, dBin_Trig);
 fHistTPCFMDA->SetBinLimits(1, binning_detaTPCFMDA);
 fHistTPCFMDA->SetBinLimits(2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 //fHistTPCFMDA->SetBinLimits(3, dCen);

 fHistTPCFMDA->SetVarTitle(0,"p_{T}_{TPC-trig} GeV/c");
 fHistTPCFMDA->SetVarTitle(1,"#Delta#eta_{TPC-FMDA}");
 fHistTPCFMDA->SetVarTitle(2,"#Delta#phi_{TPC-FMDA}");
 //fHistTPCFMDA->SetVarTitle(3,"Centrality"); 

 fHistTPCFMDA_Mixed  = new AliTHn("fHistTPCFMDA_Mixed", "fHistTPCFMDA_Mixed",  1, nVar_TPCFMDA , nbins_dTPCFMDA);
 fHistTPCFMDA_Mixed->SetBinLimits(0, dBin_Trig);
 fHistTPCFMDA_Mixed->SetBinLimits(1, binning_detaTPCFMDA);
 fHistTPCFMDA_Mixed->SetBinLimits(2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 //fHistTPCFMDA_Mixed->SetBinLimits(3, dCen);

 fHistTPCFMDA_Mixed->SetVarTitle(0,"p_{T}_{TPC-trig} GeV/c");
 fHistTPCFMDA_Mixed->SetVarTitle(1,"#Delta#eta_{TPC-FMDA}");
 fHistTPCFMDA_Mixed->SetVarTitle(2,"#Delta#phi_{TPC-FMDA}");
 //fHistTPCFMDA_Mixed->SetVarTitle(3,"Centrality");

 if(fMode == "TPCFMDA")
 {
  fListOutput->Add(fHistTPCFMDA);
  fListOutput->Add(fHistTPCFMDA_Mixed);
  fListOutput->Add(fHistTPCTrig);
 }


 fHistTPCFMDC  = new AliTHn("fHistTPCFMDC", "fHistTPCFMDC",  1, nVar_TPCFMDC , nbins_dTPCFMDC);
 fHistTPCFMDC->SetBinLimits(0, dBin_Trig);
 fHistTPCFMDC->SetBinLimits(1, binning_detaTPCFMDC);
 fHistTPCFMDC->SetBinLimits(2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
// fHistTPCFMDC->SetBinLimits(3, dCen);

 fHistTPCFMDC->SetVarTitle(0,"p_{T}_{TPC-trig} GeV/c");
 fHistTPCFMDC->SetVarTitle(1,"#Delta#eta_{TPC-FMDA}");
 fHistTPCFMDC->SetVarTitle(2,"#Delta#phi_{TPC-FMDA}");
// fHistTPCFMDC->SetVarTitle(3,"Centrality"); 

 fHistTPCFMDC_Mixed  = new AliTHn("fHistTPCFMDC_Mixed", "fHistTPCFMDC_Mixed",  1, nVar_TPCFMDC , nbins_dTPCFMDC);
 fHistTPCFMDC_Mixed->SetBinLimits(0, dBin_Trig);
 fHistTPCFMDC_Mixed->SetBinLimits(1, binning_detaTPCFMDC);
 fHistTPCFMDC_Mixed->SetBinLimits(2, -0.5*TMath::Pi(), 1.5*TMath::Pi());
// fHistTPCFMDC_Mixed->SetBinLimits(3, dCen);

 fHistTPCFMDC_Mixed->SetVarTitle(0,"p_{T}_{TPC-trig} GeV/c");
 fHistTPCFMDC_Mixed->SetVarTitle(1,"#Delta#eta_{TPC-FMDA}");
 fHistTPCFMDC_Mixed->SetVarTitle(2,"#Delta#phi_{TPC-FMDA}");
// fHistTPCFMDC_Mixed->SetVarTitle(3,"Centrality");

 if(fMode == "TPCFMDC")
 {
  fListOutput->Add(fHistTPCFMDC);
  fListOutput->Add(fHistTPCFMDC_Mixed);
  fListOutput->Add(fHistTPCTrig);
 }

 fHistFMDAFMDC = new AliTHn("fHistFMDAFMDC","fHistFMDAFMDC", 1, nVar_FMDAFMDC, nbins_dFMDAFMDC);
 fHistFMDAFMDC->SetBinLimits(0, binning_detaFMDAFMDC);
 fHistFMDAFMDC->SetBinLimits(1, -0.5*TMath::Pi(), 1.5*TMath::Pi());
// fHistFMDAFMDC->SetBinLimits(2, dCen);

 fHistFMDAFMDC->SetVarTitle(0,"#Delta#eta_{FMDA-FMDC}");
 fHistFMDAFMDC->SetVarTitle(1,"#Delta#phi_{FMDA-FMDC}");
// fHistFMDAFMDC->SetVarTitle(2,"Centrality");

 fHistFMDAFMDC_Mixed = new AliTHn("fHistFMDAFMDC_Mixed","fHistFMDAFMDC_Mixed", 1, nVar_FMDAFMDC, nbins_dFMDAFMDC);
 fHistFMDAFMDC_Mixed->SetBinLimits(0, binning_detaFMDAFMDC);
 fHistFMDAFMDC_Mixed->SetBinLimits(1, -0.5*TMath::Pi(), 1.5*TMath::Pi());
// fHistFMDAFMDC_Mixed->SetBinLimits(2, dCen);

 fHistFMDAFMDC_Mixed->SetVarTitle(0,"#Delta#eta_{FMDA-FMDC}");
 fHistFMDAFMDC_Mixed->SetVarTitle(1,"#Delta#phi_{FMDA-FMDC}");
// fHistFMDAFMDC_Mixed->SetVarTitle(2,"Centrality");

 const Int_t nbins_dFMDTrig[] = {nCen};
 const Int_t nVar_FMDTrig = sizeof(nbins_dFMDTrig) / sizeof(Int_t);

 fHistFMDTrig = new AliTHn("fHistFMDTrig","fHistFMDTrig", 1, nVar_FMDTrig, nbins_dFMDTrig); 
 fHistFMDTrig->SetBinLimits(0, dCen);
 fHistFMDTrig->SetVarTitle(0, "Centrality");

 if(fMode == "FMDAFMDC")
 {
  fListOutput->Add(fHistFMDAFMDC);
  fListOutput->Add(fHistFMDAFMDC_Mixed);
  fListOutput->Add(fHistFMDTrig);
 }


//=============================== For TPC-TPC-FMDA

 const Int_t nbins_dTPCTPCFMDA[] = {16, fNbinsPtTrig, ndetatpcfmda, 18, 24};
 const Int_t nVar_TPCTPCFMDA = sizeof(nbins_dTPCTPCFMDA) / sizeof(Int_t);

 fHistTPCTPCFMDA = new AliTHn("fHistTPCTPCFMDA","fHistTPCTPCFMDA", 1, nVar_TPCTPCFMDA, nbins_dTPCTPCFMDA);
 fHistTPCTPCFMDA->SetBinLimits(0, -1.6,1.6);
 fHistTPCTPCFMDA->SetBinLimits(1, dBin_Trig);
 fHistTPCTPCFMDA->SetBinLimits(2, binning_detaTPCFMDA);
 fHistTPCTPCFMDA->SetBinLimits(3, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistTPCTPCFMDA->SetBinLimits(4, -0.5*TMath::Pi(), 1.5*TMath::Pi());
// fHistTPCTPCFMDA->SetBinLimits(5, dBin_Asso);
  
 fHistTPCTPCFMDA->SetVarTitle(0,"#Delta#eta_{TPC-TPC}");
 fHistTPCTPCFMDA->SetVarTitle(1,"p_{T}_{TPC-trig} GeV/c");
 fHistTPCTPCFMDA->SetVarTitle(2,"#Delta#eta_{TPC-FMD}");
 fHistTPCTPCFMDA->SetVarTitle(3,"#Delta#phi_{TPC-FMD}");
 fHistTPCTPCFMDA->SetVarTitle(4,"#Delta#phi_{TPC-TPC}");
// fHistTPCTPCFMDA->SetVarTitle(5,"p_{T}_{TPC-assoc} GeV/c");

 fHistTPCTPCFMDA_Mixed = new AliTHn("fHistTPCTPCFMDA_Mixed","fHistTPCTPCFMDA_Mixed", 1, nVar_TPCTPCFMDA, nbins_dTPCTPCFMDA);
 fHistTPCTPCFMDA_Mixed->SetBinLimits(0, -1.6,1.6);
 fHistTPCTPCFMDA_Mixed->SetBinLimits(1, dBin_Trig);
 fHistTPCTPCFMDA_Mixed->SetBinLimits(2, binning_detaTPCFMDA);
 fHistTPCTPCFMDA_Mixed->SetBinLimits(3, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistTPCTPCFMDA_Mixed->SetBinLimits(4, -0.5*TMath::Pi(), 1.5*TMath::Pi());
  
 fHistTPCTPCFMDA_Mixed->SetVarTitle(0,"#Delta#eta_{TPC-TPC}");
 fHistTPCTPCFMDA_Mixed->SetVarTitle(1,"p_{T}_{TPC-trig} GeV/c");
 fHistTPCTPCFMDA_Mixed->SetVarTitle(2,"#Delta#eta_{TPC-FMD}");
 fHistTPCTPCFMDA_Mixed->SetVarTitle(3,"#Delta#phi_{TPC-FMD}");
 fHistTPCTPCFMDA_Mixed->SetVarTitle(4,"#Delta#phi_{TPC-TPC}");


 const Int_t nbins_dTPCTPCFMDA_Trig[] = {fNbinsPtTrig, 16, 24};
 const Int_t nVar_TPCTPCFMDA_Trig = sizeof(nbins_dTPCTPCFMDA_Trig) / sizeof(Int_t);
 fHistTPCTPCFMDATrig = new AliTHn("fHistTPCTPCFMDATrig","fHistTPCTPCFMDATrig", 1, nVar_TPCTPCFMDA_Trig, nbins_dTPCTPCFMDA_Trig);
 fHistTPCTPCFMDATrig->SetBinLimits(0, dBin_Trig); 
 fHistTPCTPCFMDATrig->SetBinLimits(1, -1.6 , 1.6); 
 fHistTPCTPCFMDATrig->SetBinLimits(2, -0.5*TMath::Pi(),1.5*TMath::Pi()); 
 
 fHistTPCTPCFMDATrig->SetVarTitle(0, "leading p_{T} GeV/c");
 fHistTPCTPCFMDATrig->SetVarTitle(1, "#Delta#eta");
 fHistTPCTPCFMDATrig->SetVarTitle(2, "#Delta#phi");

 if(fMode == "TPCTPCFMDA")
 {
  fListOutput->Add(fHistTPCTPCFMDA);
  fListOutput->Add(fHistTPCTPCFMDA_Mixed);
  fListOutput->Add(fHistTPCTPCFMDATrig);
 }

//=====================================================================


 const Int_t nzvtx = 1;
 Double_t zvtxbins[nzvtx+1] = {-10.,10.};

 const Int_t nCentralityBins  = 15;
 Double_t centBins[] = {0.,1.,2.,3.,4.,5., 10.,20.,30.,40.,50.,60.,70.,80.,90.,100.1};

 fPoolMgr1 = new AliEventPoolManager(20000, 50000, nCentralityBins, centBins, nzvtx, zvtxbins);
 if (!fPoolMgr1) return;
 fPoolMgr1->SetTargetValues(50000, 0.1, 5);


 fPoolMgr2 = new AliEventPoolManager(2000, 50000, nCentralityBins, centBins, nzvtx, zvtxbins);
 if (!fPoolMgr2) return;
 fPoolMgr2->SetTargetValues(50000, 0.1, 5);

//=============================================================================

  PostData(1, fListOutput);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::UserExec(Option_t *)
{
//
//  void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::UserExec(Option_t *)
//

  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(dynamic_cast<AliVEvent*>(MCEvent()));
  if(!ev) { AliFatal("MC event not found!"); }

  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!");  }


  const AliVVertex *spdVtx = ev->GetPrimaryVertex();
  fVtxZ = spdVtx->GetZ();
  Int_t nMCAllTracks = ev->GetNumberOfTracks();


//============================================================================= Random number
// if(rand()%4 == 1) return;    
//============================================================================  For Multiplicity
  Int_t nV0A = 0;
  Int_t nCharged_mid = 0;
  for (auto i=0; i<ev->GetNumberOfTracks(); ++i)
  {
   AliMCParticle *mcTrack = (AliMCParticle *)ev->GetTrack(i);
   if (!mcTrack) {
      Error("ReadEventAODMC", "Could not receive particle %d", i);
      continue;
   }
   Int_t pdgabs=TMath::Abs(mcTrack->PdgCode());
   Bool_t TrIsPrim=mcTrack->IsPhysicalPrimary();
   Bool_t TrCharge=mcTrack->Charge()!=0;

   if (!TrCharge)        continue;
   if (!TrIsPrim)           continue;
   if (mcTrack->Pt() < 0.001 || mcTrack->Pt() > 50.) continue;
   if (pdgabs!=211 && pdgabs!=321 && pdgabs!=2212) continue; //only charged pi+K+p
   if (pdgabs==9902210) return; //no diffractive protons

   const auto dEtapp(mcTrack->Eta());

   if(dEtapp>2.8 && dEtapp<5.1)
   {
    ++nV0A;
   }
   if(dEtapp>-0.5 && dEtapp<0.5)
   {
    ++nCharged_mid;
   } 
  }

  (static_cast<TH1D*>(fListOutput->FindObject("hChargeV0A")))->Fill(nV0A);

//Centrality Selection
 TH1D *h_Charge = dynamic_cast<TH1D*>(GetInputData(1));
 Double_t dMul_All = h_Charge->Integral(1,-1);
 fCentrality = h_Charge->Integral(h_Charge->FindBin(nV0A),-1) / dMul_All * 100;
 (static_cast<TH1D*>(fListOutput->FindObject("hCentrality")))->Fill(fCentrality);
 (static_cast<TH2D*>(fListOutput->FindObject("hCentrality_MidCharged")))->Fill(fCentrality,nCharged_mid);
 if(fMode == "Centrality") return;

 for (auto i=0; i<ev->GetNumberOfTracks(); ++i)
 {
  AliMCParticle *mcTrack = (AliMCParticle *)ev->GetTrack(i);
   if (!mcTrack) {
      Error("ReadEventAODMC", "Could not receive particle %d", i);
      continue;
   }
  Int_t pdgabs=TMath::Abs(mcTrack->PdgCode());
  Bool_t TrIsPrim=mcTrack->IsPhysicalPrimary();
  Bool_t TrCharge=mcTrack->Charge()!=0;

  if (!TrCharge)        continue;
  if (!TrIsPrim)           continue;
  if (mcTrack->Pt() < 0.001 || mcTrack->Pt() > 50.) continue;
  if (pdgabs!=211 && pdgabs!=321 && pdgabs!=2212) continue; //only charged pi+K+p
  if (pdgabs==9902210) return; //no diffractive protons 

  const auto dEtapp(mcTrack->Eta());
  const auto dEtappt(mcTrack->Pt());
 
  (static_cast<TH2D*>(fListOutput->FindObject("hPt_Cen")))->Fill(dEtappt,fCentrality);
  if(dEtapp>-0.8 && dEtapp<0.8) (static_cast<TH2D*>(fListOutput->FindObject("hPt_Cen_mid")))->Fill(dEtappt,fCentrality);
 } 

 if(fMode == "RCP") return;
 if(fCentrality <fCen1 || fCentrality>fCen2) return;

//Construct the container for "TPC" and "FMD"
 TObjArray *selectedTracksTPC1 = new TObjArray;
 selectedTracksTPC1->SetOwner(kTRUE);

 TObjArray *selectedTracksTPC2 = new TObjArray;
 selectedTracksTPC2->SetOwner(kTRUE);

 TObjArray *selectedTracksFMDA = new TObjArray;
 selectedTracksFMDA->SetOwner(kTRUE);

 TObjArray *selectedTracksFMDC = new TObjArray;
 selectedTracksFMDC->SetOwner(kTRUE);

 TObjArray *selected_TPC_Pairs = new TObjArray; 
 selected_TPC_Pairs->SetOwner(kTRUE);

 Double_t mcTrackEta=-999.;
 Double_t mcTrackPt=-999.;
 Double_t mcTrackPhi=-999.;
 Bool_t TrIsPrim=kFALSE; 
 Bool_t TrCharge=kFALSE; 

 for (auto i=0; i<ev->GetNumberOfTracks(); ++i)
 {
  AliMCParticle *mcTrack = (AliMCParticle *)ev->GetTrack(i);
  if (!mcTrack) {
     Error("ReadEventAODMC", "Could not receive particle %d", i);
     continue;
  }

  TrIsPrim=mcTrack->IsPhysicalPrimary();
  TrCharge=mcTrack->Charge()!=0;
  if(!TrCharge)        continue;
  if(!TrIsPrim)           continue;
 
  mcTrackEta = mcTrack->Eta();
  mcTrackPt  = mcTrack->Pt();
  mcTrackPhi = mcTrack->Phi(); 

  if(mcTrackEta>-0.8 && mcTrackEta<0.8)
  {
   if(mcTrackPt>fPtMax) continue;
   if(mcTrackPt<fPtMin) continue;
   selectedTracksTPC1->Add(new AliAssociatedTrackYSLEGOMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
   selectedTracksTPC2->Add(new AliAssociatedTrackYSLEGOMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1)); 
   (static_cast<TH1D*>(fListOutput->FindObject("hEta_Phi")))->Fill(mcTrackEta,mcTrackPhi);
  }

  if(mcTrackEta>1.7 && mcTrackEta<4.9)
  {
   if(mcTrackPt>fPtMax) continue;
   if(mcTrackPt<fPtMin) continue;
   selectedTracksFMDA->Add(new AliAssociatedTrackYSLEGOMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
   (static_cast<TH1D*>(fListOutput->FindObject("hEta_Phi")))->Fill(mcTrackEta,mcTrackPhi);
  }

  if(mcTrackEta>-3.4 && mcTrackEta<-1.7)
  {
   if(mcTrackPt>fPtMax) continue;
   if(mcTrackPt<fPtMin) continue;
   selectedTracksFMDC->Add(new AliAssociatedTrackYSLEGOMC(mcTrack->Charge(),mcTrackEta,mcTrack->Phi(),mcTrack->Pt(),mcTrack->GetLabel(),-999,-999,0, 1));
   (static_cast<TH1D*>(fListOutput->FindObject("hEta_Phi")))->Fill(mcTrackEta,mcTrackPhi);
  }
 }

 if(fMode == "TPCTPC")
 {  
  FillCorrelationTracksTPCTPC(selectedTracksTPC1,selectedTracksTPC2,fHistTPCTrig,fHistTPCTPC_SS,selected_TPC_Pairs,fCentrality);
  FillCorrelationTracksTPCTPCMixed(fCentrality, selectedTracksTPC1, selectedTracksTPC2, fHistTPCTPC_Mixed_SS);
 }

 if(fMode == "TPCFMDA")
 {
  FillCorrelationTracksTPCTPC(selectedTracksTPC1,selectedTracksFMDA,fHistTPCTrig,fHistTPCFMDA,selected_TPC_Pairs,fCentrality);
  FillCorrelationTracksTPCTPCMixed(fCentrality, selectedTracksTPC1, selectedTracksFMDA, fHistTPCFMDA_Mixed);
 }

 if(fMode == "TPCFMDC")
 {
  FillCorrelationTracksTPCTPC(selectedTracksTPC1,selectedTracksFMDC,fHistTPCTrig,fHistTPCFMDC,selected_TPC_Pairs,fCentrality);
  FillCorrelationTracksTPCTPCMixed(fCentrality, selectedTracksTPC1, selectedTracksFMDC, fHistTPCFMDC_Mixed);
 }

 if(fMode == "FMDAFMDC")
 {
  FillCorrelationTracksTPCTPC(selectedTracksFMDA,selectedTracksFMDC,fHistFMDTrig,fHistFMDAFMDC,selected_TPC_Pairs,fCentrality);
  FillCorrelationTracksTPCTPCMixed(fCentrality, selectedTracksFMDA, selectedTracksFMDC, fHistFMDAFMDC_Mixed);
 }


 if(fMode == "TPCTPCFMDA")
 { 
  FillCorrelationTracksTPCTPC(selectedTracksTPC1,selectedTracksTPC2,fHistTPCTrig,fHistTPCTPC_SS,selected_TPC_Pairs, fCentrality);
  FillCorrelationTracksTPCTPCFMDA(selected_TPC_Pairs,selectedTracksFMDA,fHistTPCTPCFMDATrig,fHistTPCTPCFMDA); 
  FillCorrelationTracksTPCTPCFMDA_Mixed(fCentrality,selected_TPC_Pairs,selectedTracksFMDA,fHistTPCTPCFMDA_Mixed);
 }

 selectedTracksTPC1->Clear();
 delete selectedTracksTPC1;
 selectedTracksTPC2->Clear();
 delete selectedTracksTPC2;
 selectedTracksFMDA->Clear();
 delete selectedTracksFMDA;
 selectedTracksFMDC->Clear();
 delete selectedTracksFMDC;
 selected_TPC_Pairs->Clear();
 delete selected_TPC_Pairs;

 ++fT;
 //cout << "Event " << fT <<endl;
  return;
}

void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::FillCorrelationTracksTPCTPC(TObjArray *triggerArray, TObjArray *assoArray, AliTHn *triggerHist, AliTHn *correlationHist, TObjArray *selected_TPC_Pairs, Double_t dcentrality) {
  
 if (!triggerHist || !correlationHist)    return;

 //========= For TPC-TPC
 if(fMode == "TPCTPC" || fMode == "TPCTPCFMDA")
 {
  Double_t binscontTrigTPCTPC[1];
  Double_t binscontTPCTPC[4] = {0.}; 

  for(Int_t i = 0; i < triggerArray->GetEntriesFast(); i++)
  {
   AliAssociatedTrackYSLEGOMC *trigger = (AliAssociatedTrackYSLEGOMC *)triggerArray->At(i);
   if (!trigger)    continue;
   Int_t trigID = trigger->GetID();
   Double_t triggerPt = trigger->Pt();
   Double_t triggerEta = trigger->Eta();
   Double_t triggerPhi = trigger->Phi();

   binscontTrigTPCTPC[0] = triggerPt;

   triggerHist->Fill(binscontTrigTPCTPC, 0);

   for(Int_t j = 0; j < assoArray->GetEntriesFast(); j++)
   {
    AliAssociatedTrackYSLEGOMC *associate =   (AliAssociatedTrackYSLEGOMC*)assoArray->At(j);
    if (!associate)        continue;
   
    if (triggerPt < associate->Pt())          continue;
    if (trigID == associate->GetID())         continue;
    if((trigger->Charge())*(associate->Charge())<0) continue;
   
    Double_t dTPC_Pairs_Eta = triggerEta-associate->Eta();
    Double_t dTPC_Pairs_phi = RangePhi(triggerPhi-associate->Phi());
    //cout << dTPC_Pairs_phi << endl;
    binscontTPCTPC[0] = triggerPt;
    binscontTPCTPC[1] = associate->Pt();
    binscontTPCTPC[2] = dTPC_Pairs_Eta;
    binscontTPCTPC[3] = dTPC_Pairs_phi;

    correlationHist->Fill(binscontTPCTPC, 0);
   
    selected_TPC_Pairs->Add(new AliTrigAssoPairST(trigger->Charge(), triggerEta, triggerPhi, triggerPt, associate->Pt(), trigger->GetID(), -999, -999, 0, 1, dTPC_Pairs_Eta,dTPC_Pairs_phi));
   }
  }
 }

 if(fMode == "TPCFMDA" || fMode == "TPCFMDC")
 {
  Double_t binscontTrigTPCFMD[1];
  Double_t binscontTPCFMD[3] = {0.}; 
  for(Int_t i = 0; i < triggerArray->GetEntriesFast(); i++)
  {
   AliAssociatedTrackYSLEGOMC *trigger = (AliAssociatedTrackYSLEGOMC *)triggerArray->At(i);
   if (!trigger)    continue;
   Int_t trigID = trigger->GetID();
   Double_t triggerPt = trigger->Pt();
   Double_t triggerEta = trigger->Eta();
   Double_t triggerPhi = trigger->Phi();
 
   binscontTrigTPCFMD[0] = triggerPt;
   triggerHist->Fill(binscontTrigTPCFMD, 0);
  
   for(Int_t j = 0; j < assoArray->GetEntriesFast(); j++)
   {
    AliAssociatedTrackYSLEGOMC *associate =   (AliAssociatedTrackYSLEGOMC*)assoArray->At(j);
    if (!associate)        continue;

    Double_t dEta_TPCFMD = triggerEta-associate->Eta();
    Double_t dPhi_TPCFMD = RangePhi(triggerPhi-associate->Phi());
   //cout << dTPC_Pairs_phi << endl;
    binscontTPCFMD[0] = triggerPt;
    binscontTPCFMD[1] = dEta_TPCFMD;
    binscontTPCFMD[2] = dPhi_TPCFMD;
    correlationHist->Fill(binscontTPCFMD, 0);
   }
  }
 }

 if(fMode == "FMDAFMDC")
 {
  Double_t binscontTrigFMDAFMDC[1];
  Double_t binscontFMDAFMDC[2] = {0.}; 
  for(Int_t i = 0; i < triggerArray->GetEntriesFast(); i++)
  {
   AliAssociatedTrackYSLEGOMC *trigger = (AliAssociatedTrackYSLEGOMC *)triggerArray->At(i);
   if (!trigger)    continue;
   Double_t triggerEta = trigger->Eta();
   Double_t triggerPhi = trigger->Phi();
   binscontTrigFMDAFMDC[0] = dcentrality;
   fHistFMDTrig->Fill(binscontTrigFMDAFMDC, 0);  
   for(Int_t j = 0; j < assoArray->GetEntriesFast(); j++)
   {
    AliAssociatedTrackYSLEGOMC *associate =   (AliAssociatedTrackYSLEGOMC*)assoArray->At(j);
    if (!associate)        continue;
    Double_t dEta_FMDAFMDC = triggerEta-associate->Eta();
    Double_t dPhi_FMDAFMDC = RangePhi(triggerPhi-associate->Phi());
    binscontFMDAFMDC[0] = dEta_FMDAFMDC;
    binscontFMDAFMDC[1] = dPhi_FMDAFMDC;
  //  binscontFMDAFMDC[2] = dcentrality;
    correlationHist->Fill(binscontFMDAFMDC, 0);
   }
  }
 }
}


void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::FillCorrelationTracksTPCTPCMixed(Double_t dcentrality, TObjArray *triggerArray, TObjArray *associateArray, AliTHn *MixedCorrelation) {

 // check if mixed event pool is ready
 AliEventPool* pool = fPoolMgr1->GetEventPool(dcentrality,0.);
 if (!pool){
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", dcentrality,
                  0.));
 }
  
 Bool_t poolReady = (pool->IsReady() || pool->NTracksInPool() > 5000 || pool->GetCurrentNEvents() > 5);
 Double_t binscont[4];
 Double_t binscontTPCFMD[3];
 Double_t binscontFMDAFMDC[2];

 if(poolReady)
 {
  for(Int_t j = 0; j < triggerArray->GetEntriesFast(); ++j) {
  
   AliAssociatedTrackYSLEGOMC *trigger = (AliAssociatedTrackYSLEGOMC*)triggerArray->At(j);
   Double_t triggerPt = trigger->Pt();
   Double_t triggerEta = trigger->Eta();
   Double_t triggerPhi = trigger->Phi();
 
   Int_t ptBin = fPtTrigAxis->FindBin(triggerPt);
   if (ptBin<1 || ptBin>fNbinsPtTrig) continue;  
    
   for(Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) {
    
    TObjArray *mixEvents = pool->GetEvent(jMix);
    for (Int_t jTrk=0; jTrk<mixEvents->GetEntriesFast(); jTrk++) {

     AliAssociatedTrackYSLEGOMC *associate = (AliAssociatedTrackYSLEGOMC*)mixEvents->At(jTrk);
     if((trigger->Charge())*(associate->Charge())<0) continue;

      Int_t assocPtBin = fPtAssocAxis->FindBin(associate->Pt());
      if (assocPtBin<1 || assocPtBin>fNbinsAssocPt) continue;

      Double_t dphi = RangePhi(triggerPhi - associate->Phi());

      if(fMode == "TPCTPC")
      {
       binscont[0] = triggerPt;
       binscont[1] = associate->Pt();
       binscont[2] = triggerEta - associate->Eta();
       binscont[3] = dphi;
       MixedCorrelation->Fill(binscont,0);
      //fHistTPCTPC_Mixed_SS->Fill(binscont,0);
      }
    
      if(fMode == "TPCFMDA" || fMode == "TPCFMDC")  
      {
        binscontTPCFMD[0] = triggerPt;
        binscontTPCFMD[1] = triggerEta - associate->Eta();
        binscontTPCFMD[2] = dphi;
        //binscontTPCFMD[3] = dcentrality;
        MixedCorrelation->Fill(binscontTPCFMD,0);
      }

      if(fMode == "FMDAFMDC")
      { 
       binscontFMDAFMDC[0] = triggerEta - associate->Eta();
       binscontFMDAFMDC[1] = dphi;
       //binscontFMDAFMDC[2] = dcentrality;
       MixedCorrelation->Fill(binscontFMDAFMDC,0);
      }
    }
   }
  }
 }
  TObjArray* tracksClone=CloneTrack(associateArray);
  //TObjArray* tracksClone=CloneTrack(triggerArray);
  pool->UpdatePool(tracksClone);
}

void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::FillCorrelationTracksTPCTPCFMDA(TObjArray *triggerArray, TObjArray *assoArray, AliTHn *triggerHist, AliTHn *correlationHist)
{
 if (!triggerHist || !correlationHist)    return;
 
 Double_t binscontTrigTPCTPCFMDA[3] = {0.};
 Double_t binscontTPCTPCFMDA[5] = {0.};
 
 for(Int_t i = 0; i < triggerArray->GetEntriesFast(); i++)
 {
  AliTrigAssoPairST *trigger_TPC_Pair = (AliTrigAssoPairST*)triggerArray->At(i);
  if(!trigger_TPC_Pair) continue;
  if(trigger_TPC_Pair->Pt_Asso()< fAsscoptCut) continue;
  
  binscontTrigTPCTPCFMDA[0] = trigger_TPC_Pair->Pt();
  binscontTrigTPCTPCFMDA[1] = trigger_TPC_Pair->Getdeta_pairs();
  binscontTrigTPCTPCFMDA[2] = trigger_TPC_Pair->Getdphi_pairs();

  triggerHist->Fill(binscontTrigTPCTPCFMDA, 0);
  for(Int_t k = 0; k < assoArray->GetEntriesFast(); k++)
  {
   binscontTPCTPCFMDA[0] = trigger_TPC_Pair->Getdeta_pairs(); // TPC eta1-eta2  
   binscontTPCTPCFMDA[1] = trigger_TPC_Pair->Pt();; // trigger TPC pt
   AliAssociatedTrackYSLEGOMC* associate = (AliAssociatedTrackYSLEGOMC*) assoArray->At(k);
   if(!associate) continue;
   Double_t dFMD_Eta = associate->Eta();
   binscontTPCTPCFMDA[2] = trigger_TPC_Pair->Eta() - dFMD_Eta; // eta1-eta3
   binscontTPCTPCFMDA[3] = RangePhi(trigger_TPC_Pair->Phi()-associate->Phi()); // phi1-phi3
   binscontTPCTPCFMDA[4] = trigger_TPC_Pair->Getdphi_pairs(); // TPC phi1-phi2 

   correlationHist->Fill(binscontTPCTPCFMDA, 0);
  }
 }
}


void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::FillCorrelationTracksTPCTPCFMDA_Mixed(Double_t centrality, TObjArray *triggerArray, TObjArray *AssociateArray, AliTHn *correlationHist)
{
 if (!correlationHist)    return;
 
 Double_t binscont[5] = {0.};
 AliEventPool *pool = fPoolMgr2->GetEventPool(centrality, 0.);

 if (!pool){
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality,
                  0.));
 }

  Bool_t poolReady = (pool->IsReady() || pool->NTracksInPool() > 5000 || pool->GetCurrentNEvents() > 5);

  if(pool->IsReady())
  {
   Int_t nMix = pool->GetCurrentNEvents();
   for(Int_t j2 = 0; j2 < triggerArray->GetEntriesFast(); j2++)
   {
    AliTrigAssoPairST *trigger_TPC_Pair = (AliTrigAssoPairST*)triggerArray->At(j2);
    if (!trigger_TPC_Pair) continue;
    if (trigger_TPC_Pair->Pt_Asso()< fAsscoptCut) continue;
    
    binscont[0] = trigger_TPC_Pair->Getdeta_pairs(); // TPC eta1-eta2  
    binscont[1] = trigger_TPC_Pair->Pt(); // trigger TPC pt
    for (Int_t jMix = 0; jMix < nMix; jMix++) {
     TObjArray *mixEvents = pool->GetEvent(jMix);
     for(Int_t k = 0; k < mixEvents->GetEntriesFast(); k++)
     {
      AliAssociatedTrackYSLEGOMC* associate = (AliAssociatedTrackYSLEGOMC*)  mixEvents->At(k);
      if(!associate) continue;
      Double_t dFMD_Eta = associate->Eta();
      binscont[2] = trigger_TPC_Pair->Eta() - dFMD_Eta; // eta1-eta3
      binscont[3] = RangePhi(trigger_TPC_Pair->Phi()-associate->Phi()); // phi1-phi3
      binscont[4] = trigger_TPC_Pair->Getdphi_pairs(); // TPC phi1-phi2
  //    binscont[5] = trigger_TPC_Pair->Pt_Asso(); // associate TPC pt 
     correlationHist->Fill(binscont, 0);
     }
    }
   }
  }

  TObjArray* tracksClone=CloneTrack(AssociateArray);
  pool->UpdatePool(tracksClone);
}




void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::SetTrigPtBinning(Int_t nBins, Double_t *limits) {

  fNbinsPtTrig = nBins;
  fPtTrigAxis  = new TAxis(fNbinsPtTrig, limits);

}

void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::SetAssocPtBinning(Int_t nBins, Double_t *limits) {
/*
  if (nBins>fNMaxBinsAssocPt) {
    AliInfo(Form("WARNING : only %d assoc pt bins (out of the %d proposed) will be considered",fNMaxBinsAssocPt,nBins));
    nBins = fNMaxBinsAssocPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one assoc pt bin must be considered");
    nBins = 1;
  }
*/
  fNbinsAssocPt = nBins;
  fPtAssocAxis  = new TAxis(fNbinsAssocPt, limits);

}

Double_t AliAnalysisTaskSEpPbCorrelationsJetV2Kine::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}

void AliAnalysisTaskSEpPbCorrelationsJetV2Kine::Terminate(Option_t *) {

  if (fPoolMgr1)    delete fPoolMgr1;
  if (fPoolMgr2)    delete fPoolMgr2;

  fListOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fListOutput) {
    printf("ERROR: Output list not available\n");
    return;
  }
/*
  fOutputList1 = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputList1) {
    printf("ERROR: Output list not available\n");
    return;
  }
*/
}

TObjArray* AliAnalysisTaskSEpPbCorrelationsJetV2Kine::CloneTrack(TObjArray*selectedTrackArray){
  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i = 0; i < selectedTrackArray->GetEntriesFast(); i++) {
    AliAssociatedTrackYSLEGOMC *particle =  (AliAssociatedTrackYSLEGOMC *)selectedTrackArray->At(i);
    tracksClone->Add(new AliAssociatedTrackYSLEGOMC(particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
                                              particle->GetID(), particle->GetIDFirstDaughter(),
                                              particle->GetIDSecondDaughter(), particle->WhichCandidate(),
                                              particle->Multiplicity()));
  }

  return tracksClone;
}

