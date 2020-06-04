#include "AliAnalysisTaskNucleiYield.h"

#include <array>
#include <string>
#include <algorithm>

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TRandom3.h>

// ALIROOT includes
#include "AdditionalFunctions.h"
#include "AliAnalysisManager.h"
#include "AliAODHeader.h"
#include "AliCentrality.h"
#include "AliPWGFunc.h"
#include "AliPDG.h"
#include "AliMultSelection.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"

#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"

using TMath::TwoPi;
using std::string;

///\cond CLASSIMP
ClassImp(AliAnalysisTaskNucleiYield);
///\endcond

namespace {

double TOFsignal(double *x, double *par) {
  double &norm = par[0];
  double &mean = par[1];
  double &sigma = par[2];
  double &tail = par[3];

  if (x[0] <= (tail + mean))
    return norm * TMath::Gaus(x[0], mean, sigma);
  else
    return norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean) / (sigma * sigma));
}

void SetupTRD2013(TF1* neg[4], TF1* pos[4]) { 
  const double fgkPhiParamPos[4][4] = {
      {1.38984e+00, -2.10187e+01, 5.81724e-02, 1.91938e+01},
      {2.02372e+00, -2.44456e+00, 8.99000e-01, 9.22399e-01},
      {4.21954e+00, -2.56555e+01, 4.17557e-02, 2.40301e+01},
      {5.17499e+00, -2.69241e+00, 6.97167e-01, 1.25974e+00}};
    const double fgkPhiParamNeg[4][4] = {
      {2.81984e+00, -1.81497e-01, -2.03494e+00, 2.64148e-01},
      {5.79322e+00, -5.44966e-02, -1.10803e+00, 1.29737e+00},
      {5.60000e+00, -2.06000e-01, -1.97130e+00, 2.67181e-01},
      {9.72180e+00, -4.35801e-02, -1.14550e+00, 1.49160e+00}};
  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    for (int iParam = 0; iParam < 4; ++iParam)
    {
      neg[iFunction]->SetParameter(iParam, fgkPhiParamNeg[iFunction][iParam]);
      pos[iFunction]->SetParameter(iParam, fgkPhiParamPos[iFunction][iParam]);
    }
  }
}

}

bool AliAnalysisTaskNucleiYield::IsInTRD(float pt, float phi, float sign) {
    bool withTRD[2]{
        phi < fTRDboundariesNeg[0]->Eval(pt) ||
            (phi > fTRDboundariesNeg[1]->Eval(pt) && phi < fTRDboundariesNeg[2]->Eval(pt)) ||
            phi > fTRDboundariesNeg[3]->Eval(pt),
        phi < fTRDboundariesPos[0]->Eval(pt) ||
            (phi > fTRDboundariesPos[1]->Eval(pt) && phi < fTRDboundariesPos[2]->Eval(pt)) ||
            phi > fTRDboundariesPos[3]->Eval(pt)};
    bool positive = sign > 0;
    return withTRD[positive];
}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskNucleiYield::AliAnalysisTaskNucleiYield(TString taskname)
  :AliAnalysisTaskSE(taskname.Data())
   ,fEventCut{false}
   ,fFilterBit{BIT(4)}
   ,fPropagateTracks{false}
   ,fPtCorrectionA{3}
   ,fPtCorrectionM{3}
   ,fOptionalTOFcleanup{-1.}
   ,fCurrentFileName{""}
   ,fTOFfunction{nullptr}
   ,fList{nullptr}
   ,fRTree{nullptr}
   ,fSTree{nullptr}
   ,fCutVec{}
   ,fPDG{0}
   ,fPDGMass{0}
   ,fPDGMassOverZ{0}
   ,fCharge{1.f}
   ,fIsMC{false}
   ,fFillOnlyEventHistos{false}
   ,fPID{nullptr}
   ,fMagField{0.f}
   ,fCentrality{-1.f}
   ,fDCAzLimit{10.}
   ,fDCAzNbins{400}
   ,fSigmaLimit{6.}
   ,fSigmaNbins{240}
   ,fTOFSigmaLimit{12.}
   ,fTOFSigmaNbins{240}
   ,fTOFlowBoundary{-2.4}
   ,fTOFhighBoundary{3.6}
   ,fTOFnBins{75}
   ,fDisableITSatHighPt{100.f}
   ,fDisableTPCpidAtHighPt{100.f}
   ,fEnablePtCorrection{false}
   ,fRequireITSrecPoints{2u}
   ,fRequireTPCrecPoints{0u}
   ,fRequireITSsignal{0u}
   ,fRequireSDDrecPoints{0u}
   ,fRequireSPDrecPoints{1u}
   ,fRequireTPCsignal{70u}
   ,fRequireEtaMin{-0.8f}
   ,fRequireEtaMax{0.8f}
   ,fRequireYmin{-0.5f}
   ,fRequireYmax{0.5f}
   ,fRequireMaxChi2{4.f}
   ,fRequireMaxDCAxy{0.12f}
   ,fRequireMaxDCAz{1.f}
   ,fRequireTPCpidSigmas{3.f}
   ,fRequireITSpidSigmas{-1.f}
   ,fRequireTOFpidSigmas{-1.f}
   ,fRequireMinEnergyLoss{0.}
   ,fRequireDeadZoneWidth{0.}
   ,fRequireCutGeoNcrNclLength{0.}
   ,fRequireCutGeoNcrNclGeom1Pt{0.}
   ,fCutGeoNcrNclFractionNcr{0.}
   ,fCutGeoNcrNclFractionNcl{0.}
   ,fRequireVetoSPD{false}
   ,fRequireMaxMomentum{-1.}
   ,fFixForLHC14a6{false}
   ,fRequireTPCfoundFraction{0.}
   ,fPtShapeFunction{kNoPtShape}
   ,fPtShapeMaximum{0.f}
   ,fITSelectronRejectionSigma{-1.}
   ,fBeamRapidity{0.f}
   ,fEstimator{0}
   ,fEnableFlattening{false}
   ,fSaveTrees{false}
   ,fRecNucleus{0.,0.,0.,0.,0.,0.,0,0,0}
   ,fSimNucleus{0.,0.,0,0,0}
   ,fParticle{AliPID::kUnknown}
   ,fCentBins{0}
   ,fDCABins{0}
   ,fPtBins{0}
   ,fCustomTPCpid{0}
   ,fFlatteningProbs{0}
   ,fPtShapeParams{0}
   ,fFunctCollection{nullptr}
   ,fPtShape{nullptr}
   ,fNormalisationHist{nullptr}
   ,fProduction{nullptr}
   ,fReconstructed{{nullptr}}
   ,fTotal{nullptr}
   ,fPtCorrection{nullptr}
   ,fPcorrectionTPC{nullptr}
   ,fDCAPrimary{{nullptr}}
   ,fDCASecondary{{nullptr}}
   ,fDCASecondaryWeak{{nullptr}}
   ,fTOFsignal{nullptr}
   ,fTOFT0FillSignal{nullptr}
   ,fTOFNoT0FillSignal{nullptr}
   ,fTPCcounts{nullptr}
   ,fMultDistributionTPC{nullptr}
   ,fMultDistributionTOF{nullptr}
   ,fTOFnSigma{nullptr}
   ,fTOFT0FillNsigma{nullptr}
   ,fTOFNoT0FillNsigma{nullptr}
   ,fTPCsignalTpl{nullptr}
   ,fTPCbackgroundTpl{nullptr}
   ,fDCAxy{{nullptr}}
   ,fDCAz{{nullptr}}
   ,fHist2Phi{nullptr}
   ,fTRDboundariesPos{nullptr}
   ,fTRDboundariesNeg{nullptr}
   ,fTRDvintage{0}
   ,fTRDin{false}
   ,fNanoPIDindexTPC{-1}
   ,fNanoPIDindexTOF{-1}
   ,fRefMult{-1}
   {
     gRandom->SetSeed(0); //TODO: provide a simple method to avoid "complete randomness"
     Float_t aCorrection[3] = {-2.10154e-03,-4.53472e-01,-3.01246e+00};
     Float_t mCorrection[3] = {-2.00277e-03,-4.93461e-01,-3.05463e+00};
     fPtCorrectionA.Set(3, aCorrection);
     fPtCorrectionM.Set(3, mCorrection);
     DefineInput(0, TChain::Class());
     DefineOutput(1, TList::Class());
     DefineOutput(2, TTree::Class());
     DefineOutput(3, TTree::Class());
   }

/// Standard destructor
///
AliAnalysisTaskNucleiYield::~AliAnalysisTaskNucleiYield(){
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fList) delete fList;
  if (fRTree) delete fRTree;
  if (fSTree) delete fSTree;
  if (fTOFfunction) delete fTOFfunction;
  if (fFunctCollection) delete fFunctCollection;
  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    if (fTRDboundariesPos[iFunction])
      delete fTRDboundariesPos[iFunction];
    if (fTRDboundariesNeg[iFunction])
      delete fTRDboundariesNeg[iFunction];
  }
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskNucleiYield::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);

  const Int_t nPtBins = fPtBins.GetSize() - 1;
  const Int_t nCentBins = fCentBins.GetSize() - 1;
  const Int_t nDCAbins = fDCABins.GetSize() - 1;
  const float *pTbins = fPtBins.GetArray();
  const float *centBins = fCentBins.GetArray();
  double doubleCentBins[nCentBins+1];
  std::copy(centBins,centBins+nCentBins+1,doubleCentBins);
  const float *dcaBins = fDCABins.GetArray();

  char   letter[2] = {'A','M'};
  string tpctof[2] = {"TPC","TOF"};
  string tpctofMC[2] = {"TPC","TPC_TOF"};

  if (fPtShapeFunction != kNoPtShape)
    fFunctCollection = new AliPWGFunc;
  switch (fPtShapeFunction) {
    case kBlastWaveShape:
      fPtShape = fFunctCollection->GetBGBW(fPDGMass, fPtShapeParams[0], fPtShapeParams[1], fPtShapeParams[2], 1.);
    case kTsallisShape:
      fPtShape = LevyTsallis("nuclei_levytsallis", fPDGMass, fPtShapeParams[0], fPtShapeParams[1], 1.);
  }
  if (fPtShape)
    fPtShapeMaximum = fPtShape->GetMaximum(0,10,1.e-10,10000);

  if (fIsMC) {
    fProduction = new TH1F("fProduction",";#it{p} (GeV/#it{c});Entries",100,-10,10);
    fList->Add(fProduction);

    for (int iC = 0; iC < 2; ++iC) {
      fTotal[iC] = new TH2F(Form("f%cTotal",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); Counts",
          nCentBins,centBins,nPtBins,pTbins);
      fPtCorrection[iC] = new TH2F(Form("f%cPtCorrection",letter[iC]),
          ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{MC}-#it{p}_{T}^{rec} (GeV/#it{c});Entries",
          160,0.4,6.,80,-1.,1.);
      fPcorrectionTPC[iC] = new TH2F(Form("f%cPcorrectionTPC",letter[iC]),
          ";#it{p}^{rec} (GeV/#it{c});#it{p}^{MC}-#it{p}^{rec} (GeV/#it{c});Entries",
          160,0.4,6.,80,-1.,1.);
      fList->Add(fTotal[iC]);
      fList->Add(fPtCorrection[iC]);
      fList->Add(fPcorrectionTPC[iC]);
      for (int iT = 0; iT < 2; ++iT) {
        fReconstructed[iT][iC] = new TH2F(Form("f%cITS_%s",letter[iC],tpctofMC[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); Counts",
            nCentBins,centBins,nPtBins,pTbins);
        fDCAPrimary[iT][iC] = new TH3F(Form("f%cDCAPrimary%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fDCASecondary[iT][iC] = new TH3F(Form("f%cDCASecondary%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fDCASecondaryWeak[iT][iC] = new TH3F(Form("f%cDCASecondaryWeak%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fList->Add(fReconstructed[iT][iC]);
        fList->Add(fDCAPrimary[iT][iC]);
        fList->Add(fDCASecondary[iT][iC]);
        fList->Add(fDCASecondaryWeak[iT][iC]);
      }
    }
  } else {

    float tofBins[fTOFnBins + 1];
    const float deltaTOF = (fTOFhighBoundary - fTOFlowBoundary) / fTOFnBins;
    for (int i = 0; i <= fTOFnBins; ++i)
      tofBins[i] = i * deltaTOF + fTOFlowBoundary;
    float dcazBins[fDCAzNbins + 1];
    const float deltaDCAz = 2.f * fDCAzLimit / fDCAzNbins;
    for (int i = 0; i <= fDCAzNbins; ++i)
      dcazBins[i] = i * deltaDCAz - fDCAzLimit;
    float sigmaBins[fSigmaNbins + 1];
    const float deltaSigma = 2.f * fSigmaLimit / fSigmaNbins;
    for (int i = 0; i <= fSigmaNbins; ++i)
      sigmaBins[i] = i * deltaSigma - fSigmaLimit;
    float TOFSigmaBins[fTOFSigmaNbins + 1];
    const float deltaTOFSigma = 2.f * fTOFSigmaLimit / fTOFSigmaNbins;
    for (int i = 0; i <= fTOFSigmaNbins; ++i)
      TOFSigmaBins[i] = i * deltaTOFSigma - fTOFSigmaLimit;

    float nSigmasBins[51];
    float multBins[51];
    for (int i = 0; i <= 50; ++i) {
      nSigmasBins[i] = -5. + i * 0.2;
      multBins[i] = i * 2.;
    }
  
    fMultDistributionTPC = new TH3F("fMultDistributionTPC",";Reference Multiplicity;#it{p}_{T} (Gev/#it{c});TPC n#sigma", 50, multBins, nPtBins,pTbins, 50, nSigmasBins);
    fMultDistributionTOF = new TH3F("fMultDistributionTOF",";Reference Multiplicity;#it{p}_{T} (Gev/#it{c});TOF n#sigma", 50, multBins, nPtBins,pTbins, 50, nSigmasBins);
    fList->Add(fMultDistributionTPC);
    fList->Add(fMultDistributionTOF);

    for (int iC = 0; iC < 2; ++iC) {
      fTOFsignal[iC] = new TH3F(Form("f%cTOFsignal",letter[iC]),
          ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{m}^{2}-m_{PDG}^{2} (GeV/#it{c}^{2})^{2}",
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fTOFT0FillSignal[iC] = new TH3F(Form("f%cTOFT0FillSignal",letter[iC]),
          ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{m}^{2}-m_{PDG}^{2} (GeV/#it{c}^{2})^{2}",
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fTOFNoT0FillSignal[iC] = new TH3F(Form("f%cTOFNoT0FillSignal",letter[iC]),
          ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{m}^{2}-m_{PDG}^{2} (GeV/#it{c}^{2})^{2}",
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fTOFnSigma[iC] = new TH3F(Form("f%cTOFnSigma",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,fTOFSigmaNbins,TOFSigmaBins);
      fTOFT0FillNsigma[iC] = new TH3F(Form("f%cTOFT0FillNsigma",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,fTOFSigmaNbins,TOFSigmaBins);
      fTOFNoT0FillNsigma[iC] = new TH3F(Form("f%cTOFNoT0FillNsigma",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,fTOFSigmaNbins,TOFSigmaBins);
      fTPCcounts[iC] = new TH3F(Form("f%cTPCcounts",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,fSigmaNbins,sigmaBins);
      fTPCsignalTpl[iC] = new TH3F(Form("f%cTPCsignalTpl",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,fSigmaNbins,sigmaBins);
      fTPCbackgroundTpl[iC] = new TH3F(Form("f%cTPCbackgroundTpl",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,fSigmaNbins,sigmaBins);
      fHist2Phi[iC] = new TH2F(Form("fHist2Phi%c", letter[iC]), Form("%c; #Phi (rad) ;#it{p}_{T} (Gev/#it{c});", letter[iC]), 100, 0, TMath::TwoPi(), 100, 0, 7);

      fList->Add(fTOFsignal[iC]);
      fList->Add(fTOFT0FillSignal[iC]);
      fList->Add(fTOFNoT0FillSignal[iC]);
      fList->Add(fTOFnSigma[iC]);
      fList->Add(fTOFT0FillNsigma[iC]);
      fList->Add(fTOFNoT0FillNsigma[iC]);
      fList->Add(fTPCcounts[iC]);
      fList->Add(fTPCsignalTpl[iC]);
      fList->Add(fTPCbackgroundTpl[iC]);
      fList->Add(fHist2Phi[iC]);

      for (int iT = 0; iT < 2; ++iT) {
        fDCAxy[iT][iC] = new TH3F(Form("f%cDCAxy%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it[c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fDCAz[iT][iC] = new TH3F(Form("f%cDCAz%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)",
            nCentBins,centBins,nPtBins,pTbins,fDCAzNbins,dcazBins);
        fList->Add(fDCAxy[iT][iC]);
        fList->Add(fDCAz[iT][iC]);
      }
    }
  }

  std::array<std::string,4> norm_labels = {
    "No cuts",
    "Event selection",
    "Vertex reconstruction and quality",
    "Vertex position"
  };
  fNormalisationHist = new TH2F("fNormalisationHist",";Centrality (%%);",nCentBins,doubleCentBins,norm_labels.size(),-.5,norm_labels.size() - 0.5);
  for (size_t iB = 1; iB <= norm_labels.size(); iB++) fNormalisationHist->GetYaxis()->SetBinLabel(iB,norm_labels[iB-1].data());
  fList->Add(fNormalisationHist);

  fTOFfunction = new TF1("fTOFfunction", TOFsignal, -2440., 2440., 4);
  if (fTOFfunctionPars.GetSize() == 4)
    fTOFfunction->SetParameters(fTOFfunctionPars.GetArray());

  AliPDG::AddParticlesToPdgDataBase();
  PostData(1,fList);

  if (fSaveTrees) {
    OpenFile(1);
    fRTree = new TTree("RTree", "Reconstructed nuclei");
    fRTree->Branch("RLightNucleus", &fRecNucleus);
    PostData(2, fRTree);

    if (fIsMC) {
      fSTree = new TTree("STree", "Simulated nuclei");
      fSTree->Branch("SLightNucleus", &fSimNucleus);
      PostData(3, fSTree);
    }
  }

  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    fTRDboundariesNeg[iFunction] = new TF1(Form("fNeg%i", iFunction), "[0]-exp([1]*pow(x,[2])+[3])", 0.2, 10);
    fTRDboundariesPos[iFunction] = new TF1(Form("fPos%i", iFunction), "[0]-exp([1]*pow(x,[2])+[3])", 0.2, 10);
  }
  if (fTRDvintage == 2013)
    SetupTRD2013(fTRDboundariesNeg, fTRDboundariesPos);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskNucleiYield::UserExec(Option_t *){
  /// The first check performed is on the particle type requested. If this type is AliPID::Unknown -
  /// the default one - the task does not process the information.
  if (fParticle == AliPID::kUnknown) {
    ::Error("AliAnalysisTaskNucleiYield::UserExec", "No particle type set");
    PostData(1, fList);
    if(fSaveTrees){
      PostData(2, fRTree);
      if(fIsMC){
        PostData(3, fSTree);
      }
    }
    return;
  }

  AliNanoAODHeader* nanoHeader = dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());
  
  AliVEvent *ev = InputEvent();

  fCentrality = -1.f;

  bool EventAccepted = true;
  if (!nanoHeader) {
    EventAccepted = fEventCut.AcceptEvent(ev);
    /// The centrality selection in PbPb uses the percentile determined with V0.
    fCentrality = fEventCut.GetCentrality(fEstimator);
  } else {
    if (fNanoPIDindexTPC == -1 || fNanoPIDindexTOF == -1) {
      AliNanoAODTrack::InitPIDIndex();
      fNanoPIDindexTPC  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, fParticle);
      fNanoPIDindexTOF  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, fParticle);
    }

    fCentrality = nanoHeader->GetCentralityV0M();
  }

  bool specialTrigger = true;
  if (fINT7intervals.size()) {
    unsigned int trigger = 0u;
    if (nanoHeader)
      trigger = nanoHeader->GetOfflineTrigger();
    else {
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
      trigger = handl->IsEventSelected() ;
    }
    bool kINT7trigger = (trigger & AliVEvent::kINT7) == AliVEvent::kINT7;

    for (int iInt = 0; iInt < fINT7intervals.size(); iInt +=2) {
      if (fCentrality >= fINT7intervals[iInt] && fCentrality < fINT7intervals[iInt+1]) {
        EventAccepted = kINT7trigger;
        specialTrigger = kINT7trigger;
        break;
      }
    }
  }
  
  if (!nanoHeader) {
    std::array <AliEventCuts::NormMask,4> norm_masks {
      AliEventCuts::kAnyEvent,
      AliEventCuts::kPassesNonVertexRelatedSelections,
      AliEventCuts::kHasReconstructedVertex,
      AliEventCuts::kPassesAllCuts
    };
    for (int iC = 0; iC < 4; ++iC) {
      if (fEventCut.CheckNormalisationMask(norm_masks[iC]) && (iC == 0 || specialTrigger)) {
        fNormalisationHist->Fill(fCentrality,iC);
      }
    }
    AliAODHeader* aodHeader = dynamic_cast<AliAODHeader*>(fInputEvent->GetHeader());
    fRefMult = aodHeader->GetRefMultiplicityComb08();
  }

  if (!EventAccepted) {
    PostData(1, fList);
    if(fSaveTrees){
      PostData(2, fRTree);
      if(fIsMC){
        PostData(3, fSTree);
      }
    }
    return;
  }

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  /// The magnetic field
  fMagField = ev->GetMagneticField();

  /// At the stage of event selection the Flattening is applied. This technique makes flat the
  /// centrality distribution using a pseudo-random selection based on prior probabilities.
  /// A complete description of this technique is present in the documentation of the Flatten
  /// function.

  if (Flatten(fCentrality) && fEnableFlattening) {
    PostData(1, fList);
    if(fSaveTrees){
      PostData(2, fRTree);
      if(fIsMC){
        PostData(3, fSTree);
      }
    }
    return;
  }

  TClonesArray *stack = nullptr;
  fRejectedParticles.clear();
  if (fIsMC || fPtShape) {
    // get branch "mcparticles"
    stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack)
      ::Fatal("AliAnalysisTaskNucleiYield::UserExec","MC analysis (or pt shape) requested on a sample without the MC particle array.");

    /// Making the list of the deuterons we want to measure
    for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
      AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
      const int pdg = std::abs(part->GetPdgCode());
      const int iC = part->Charge() > 0 ? 1 : 0;
      const int mult = -1 + 2 * iC;
      if (pdg != fPDG) continue;
      if (fPtShape) {
        if (part->IsPhysicalPrimary() && gRandom->Uniform(0, fPtShapeMaximum) > fPtShape->Eval(part->Pt())) {
          fRejectedParticles.push_back(iMC);
          continue;
        }
      }
      if (fSaveTrees) {
        SetSLightNucleus(part,fSimNucleus);
        fSTree->Fill();
      }
      if (fIsMC) fProduction->Fill(mult * part->P());
      if (part->Y() > fRequireYmax || part->Y() < fRequireYmin) continue;
      if (part->IsPhysicalPrimary() && fIsMC) fTotal[iC]->Fill(fCentrality,part->Pt());
    }
  }

  /// Checking how many deuterons in acceptance are reconstructed well
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliNanoAODTrack* nanoTrack = dynamic_cast<AliNanoAODTrack*>(ev->GetTrack(iT));
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));
    if (nanoHeader)
      TrackLoop(nanoTrack, true);
    else
      TrackLoop(aodTrack, false);
    
  } // End AOD track loop

  //  Post output data.
  PostData(1, fList);
  if (fSaveTrees) {
    PostData(2, fRTree);
    if (fIsMC) {
      PostData(3, fSTree);
    }
  }
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskNucleiYield::Terminate(Option_t *) {
  return;
}

/// This functions sets the centrality bins used in the analysis
///
/// \param nbins Number of centrality bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetCentBins(Int_t nbins, Float_t *bins) {
  fCentBins.Set(nbins + 1, bins);
}

/// This functions sets the \f$\mathrm{DCA}_{xy}\f$ bins used in the analysis
///
/// \param nbins Number of \f$\mathrm{DCA}_{xy}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetDCABins(Int_t nbins, Float_t min, Float_t max) {
  const float delta = (max - min) / nbins;
  fDCABins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB) {
    fDCABins[iB] = min + iB * delta;
  }
  fDCABins[nbins] = max;
}

/// This functions sets the \f$\mathrm{DCA}_{xy}\f$ bins used in the analysis
///
/// \param nbins Number of \f$\mathrm{DCA}_{xy}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetDCABins(Int_t nbins, Float_t *bins) {
  fDCABins.Set(nbins + 1, bins);
}

/// This functions sets the \f$p_{\mathrm{T}}\f$ bins used in the analysis
///
/// \param nbins Number of \f$p_{\mathrm{T}}\f$ bins
/// \param min Lower limit for the \f$p_{\mathrm{T}}\f$ axis
/// \param max Upper limit for the \f$p_{\mathrm{T}}\f$ axis
/// \return void
///
void AliAnalysisTaskNucleiYield::SetPtBins(Int_t nbins, Float_t min, Float_t max) {
  const float delta = (max - min) / nbins;
  fPtBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB) {
    fPtBins[iB] = min + iB * delta;
  }
  fPtBins[nbins] = max;
}

/// This functions sets the \f$p_{\mathrm{T}}\f$ bins used in the analysis
///
/// \param nbins Number of \f$p_{\mathrm{T}}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetPtBins(Int_t nbins, Float_t *bins) {
  fPtBins.Set(nbins + 1, bins);
}

/// This function allows to set a custom parametrisation for the TPC response function
///
/// \param par Array of 5 values corresponding to the Bethe Bloch parametrisation
/// \param sigma Sigma of the parametrisation
/// \return void
///
void AliAnalysisTaskNucleiYield::SetCustomTPCpid(Float_t *par, Float_t sigma) {
  if (par == 0x0 && sigma <= 0) {
    fCustomTPCpid.Set(1);
  } else {
    fCustomTPCpid.Set(6);
    for (int i = 0; i < 5; ++i)
      fCustomTPCpid.AddAt(par[i],i);
    fCustomTPCpid.AddAt(sigma, 5);
  }
}

float AliAnalysisTaskNucleiYield::GetTPCsigmas(AliVTrack* t) {
  if (fCustomTPCpid.GetSize() < 6 || fIsMC) {
    AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(t);
    return nanoT ? nanoT->GetVar(fNanoPIDindexTPC) : fPID->NumberOfSigmasTPC(t, fParticle);
  } else {
    const float p = t->GetTPCmomentum() / fPDGMass;
    const float r = AliExternalTrackParam::BetheBlochAleph(p, fCustomTPCpid[0], fCustomTPCpid[1],
        fCustomTPCpid[2], fCustomTPCpid[3], fCustomTPCpid[4]);
    return (t->GetTPCsignal() - r) / (fCustomTPCpid[5] * r); 
  }
}

float AliAnalysisTaskNucleiYield::GetTOFsigmas(AliVTrack* t) {
  AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(t);
  return nanoT ? nanoT->GetVar(fNanoPIDindexTOF) : fPID->NumberOfSigmasTOF(t, fParticle);
}

/// This function checks if the track passes the PID selection
///
/// \param t Track to be tested
/// \param sigmas Number of sigmas
/// \return Boolean value: true means that the track passes the PID selection
///
int AliAnalysisTaskNucleiYield::PassesPIDSelection(AliAODTrack *t) {
  bool tofPID = true, itsPID = true, tpcPID = true, electronRejection = true;

  if (fRequireITSpidSigmas > 0 && t->Pt() < fDisableITSatHighPt) {
    itsPID = TMath::Abs(fPID->NumberOfSigmasITS(t, fParticle)) < fRequireITSpidSigmas;
  }
  electronRejection = TMath::Abs(fPID->NumberOfSigmasITS(t, AliPID::kElectron)) > fITSelectronRejectionSigma;

  if (fRequireTOFpidSigmas > 0) {
    tofPID = TMath::Abs(fPID->NumberOfSigmasTOF(t, fParticle)) < fRequireTOFpidSigmas;
  }

  if (t->Pt() < fDisableTPCpidAtHighPt) {
    tpcPID = TMath::Abs(GetTPCsigmas(t)) < fRequireTPCpidSigmas;
  }

  return int(itsPID) | int(tpcPID) << 1 | int(tofPID) << 2| int(electronRejection) << 3;
}

int AliAnalysisTaskNucleiYield::PassesPIDSelection(AliNanoAODTrack *t) {
  bool tofPID = true, itsPID = true, tpcPID = true, electronRejection = true;

  if (fRequireITSpidSigmas > 0 && t->Pt() < fDisableITSatHighPt) {
    AliFatal("ITS PID not implemented for NanoAOD");
    // itsPID = TMath::Abs(fPID->NumberOfSigmasITS(t, fParticle)) < fRequireITSpidSigmas;
  }
  if (fITSelectronRejectionSigma > 0)
    AliFatal("Electron rejection not implemented for NanoAOD");
  // electronRejection = TMath::Abs(fPID->NumberOfSigmasITS(t, AliPID::kElectron)) > fITSelectronRejectionSigma;

  if (fRequireTOFpidSigmas > 0) {
    tofPID = TMath::Abs(t->GetVar(t->GetPIDIndex(AliNanoAODTrack::kSigmaTOF, fParticle))) < fRequireTOFpidSigmas;
  }

  if (t->Pt() < fDisableTPCpidAtHighPt) {
    tpcPID = TMath::Abs(GetTPCsigmas(t)) < fRequireTPCpidSigmas;
  }

  return int(itsPID) | int(tpcPID) << 1 | int(tofPID) << 2| int(electronRejection) << 3;
}

/// This function sets the number of TOF bins and the boundaries of the histograms
///
/// \param nbins Number of bins
/// \param min Lower boundary of the histogram
/// \param max Higher boundary of the histogram
/// \return void
///
void AliAnalysisTaskNucleiYield::SetTOFBins(Int_t nbins, Float_t min, Float_t max) {
  fTOFnBins = nbins;
  fTOFlowBoundary = min;
  fTOFhighBoundary = max;
}

/// This function sets the number of DCA\f$_{z}\f$ bins and the boundaries of the histogram
///
/// \param nbins Number of bins
/// \param limit Boundaries of the histogram (symmetrical with respect to zero)
/// \return void
///
void AliAnalysisTaskNucleiYield::SetDCAzBins(Int_t nbins, Float_t limit) {
  fDCAzNbins = nbins;
  fDCAzLimit = limit;
}

/// This function sets the number of n\f$_{sigma}\f$ bins and the boundaries of the histogram
///
/// \param nbins Number of bins
/// \param limit Boundaries of the histogram (symmetrical with respect to zero)
/// \return void
///
void AliAnalysisTaskNucleiYield::SetSigmaBins(Int_t nbins, Float_t limit) {
  fSigmaNbins = nbins;
  fSigmaLimit = limit;
}

/// This function sets the number of n\f$_{sigma_{TOF}}\f$ bins and the boundaries of the histogram
///
/// \param nbins Number of bins
/// \param limit Boundaries of the histogram (symmetrical with respect to zero)
/// \return void
///
void AliAnalysisTaskNucleiYield::SetTOFSigmaBins(Int_t nbins, Float_t limit) {
  fTOFSigmaNbins = nbins;
  fTOFSigmaLimit = limit;
}

/// This function sets the particle type to be analysed
///
/// \param part Particle type
/// \return void
///
void AliAnalysisTaskNucleiYield::SetParticleType(AliPID::EParticleType part) {
  fParticle = part;
  fPDGMass = AliPID::ParticleMass(part);
  fPDGMassOverZ = AliPID::ParticleMassZ(part);
  fCharge = TMath::Abs(AliPID::ParticleCharge(fParticle));
}

/// This function provides the flattening of the centrality distribution.
/// Please check the hardcoded values! It is better to provide those number by yourself: the
/// probability is computed as \f[\mathrm{Probability}=\frac{C_{i}}{C_{ref}} \f] where \f$C_{i}\f$
/// is the centrality in the bin _i_ and \f$C_{ref}\f$ is the centrality of the reference bin
/// (namely the value around you want the centrality to fluctuate).
///
/// \param cent Event centrality
/// \return Boolean value: true means that the event must be skipped
///
Bool_t AliAnalysisTaskNucleiYield::Flatten(float cent) {
  if (fFlatteningProbs.GetSize() <= 0) {
    Float_t prob[13] = {
      0.839266,0.822364,0.807522,0.804727,0.806675,
      0.828297,0.820842,0.834088,0.861455,1.,
      0.38112,0.661154,0.953928
    };
    fFlatteningProbs.Set(13, prob);
  }
  if (!fIsMC) {
    if (cent >= fFlatteningProbs.GetSize()) return false;
    else return gRandom->Rndm() > fFlatteningProbs[int(cent)];
  } else {
    // This flattening is required since a strange peak in VOM distribution is observed in MC
    if (fFixForLHC14a6) {
      if (cent < 1.5e-3f) {
        return true;
      } else if (cent < 0.05 && gRandom->Rndm() < 0.5) {
        return true;
      }
    }
  }
  return false;
}

/// This function provides the correction for wrongly calculated \f$p_{\mathrm{T}}\f$.
///
/// \param pt \f$p_{\mathrm{T}}\f$ of the track
/// \param positiveCharge True if the track has positive sign.
/// \return void
///
void AliAnalysisTaskNucleiYield::PtCorrection(float &pt, bool positiveCharge) {
  Float_t *par = (positiveCharge) ? fPtCorrectionM.GetArray() : fPtCorrectionA.GetArray();
  const Float_t correction = par[0] + par[1] * TMath::Exp(par[2] * pt * pt * pt);
  pt -= correction;
}

/// This function returns the number of clusters for each ITS subdetector
///
/// \param track
/// \param nSPD number of clusters in SPD
//  \param nSDD number of clusters in SDD
//  \param nSSD number of clusters in SSD
/// \return int number of clusters in ITS
///
int AliAnalysisTaskNucleiYield::GetNumberOfITSclustersPerLayer(AliVTrack *track, int &nSPD, int &nSDD, int &nSSD) {
  if (!track) return -1;
  nSPD = 0u;
  nSDD = 0u;
  nSSD = 0u;
  for (int i = 0; i < 6; ++i) {
    if (track->HasPointOnITSLayer(i)) {
      if (i < 2) nSPD++;
      else if (i < 4) nSDD++;
      else nSSD++;
    }
  }
  return nSPD + nSDD + nSSD;
}

void AliAnalysisTaskNucleiYield::SetSLightNucleus(AliAODMCParticle* part, SLightNucleus& snucl) {
  snucl.pt = part->Pt();
  snucl.eta = part->Eta();
  snucl.centrality = fCentrality;
  snucl.pdg = part->GetPdgCode();
  if (part->IsPhysicalPrimary())
    snucl.flag = SLightNucleus::kPrimary;
  else if (part->IsSecondaryFromWeakDecay())
    snucl.flag = SLightNucleus::kSecondaryWeakDecay;
  else
    snucl.flag = SLightNucleus::kSecondaryMaterial;
}


/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
/// \return \f$\beta\f$ of the particle, -1 means that there is no correct prolongation in TOF.
///
float AliAnalysisTaskNucleiYield::HasTOF(AliVTrack *track, AliPIDResponse *pid) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = hasTOFout && hasTOFtime && (len > 350.);


  if (!hasTOF) return -1.;
  const float tim = track->GetTOFsignal() - pid->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
  const float beta = len / (tim * LIGHT_SPEED);
  return beta;
}
/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
/// \return \f$\beta\f$ of the particle, -1 means that there is no correct prolongation in TOF.
///
float AliAnalysisTaskNucleiYield::HasTOF(AliNanoAODTrack *track, AliPIDResponse *pid) {
  const float len = track->GetIntegratedLength();
  bool hasTOF = track->HasTOFpid() && (len > 350.);

  if (!hasTOF) return -1.;
  const float tim = track->GetTOFsignal() - pid->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
  const float beta = len / (tim * LIGHT_SPEED);
  return beta;
}

/// This function checks whether a track pass TPC Geometrical cut
///
/// \param track Track that has to be checked
/// \return Boolean value: true means that track passed TPC Geometrical cut
///
Bool_t AliAnalysisTaskNucleiYield::IsSelectedTPCGeoCut(AliAODTrack *track) {
  Bool_t checkResult = kTRUE;
  AliESDtrack esdTrack(track);
  esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(track->GetTPCNclsF());

  float nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
  float lengthInActiveZoneTPC=esdTrack.GetLengthInActiveZone(0,fRequireDeadZoneWidth,220.,fMagField);
  double cutGeoNcrNclLength=fRequireCutGeoNcrNclLength-TMath::Power(TMath::Abs(esdTrack.GetSigned1Pt()),fRequireCutGeoNcrNclGeom1Pt);
  
  if (lengthInActiveZoneTPC < cutGeoNcrNclLength) checkResult = kFALSE;
  if (nCrossedRowsTPC<fCutGeoNcrNclFractionNcr*cutGeoNcrNclLength) checkResult=kFALSE;
  if (esdTrack.GetTPCncls()<fCutGeoNcrNclFractionNcl*cutGeoNcrNclLength) checkResult=kFALSE;
  
  return checkResult;
}
Bool_t AliAnalysisTaskNucleiYield::IsSelectedTPCGeoCut(AliNanoAODTrack *track) {
  static const Int_t tpcGeo_index = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCGeoLength");
  if(static_cast<AliNanoAODTrack*>(track)->GetVar(tpcGeo_index) > 0.5)
    return kTRUE;
  else
    return kFALSE;
}
