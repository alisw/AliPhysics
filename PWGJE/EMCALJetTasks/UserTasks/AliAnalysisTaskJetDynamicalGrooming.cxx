/**
 * Dynamical grooming jet substructure task. Adapted from AliAnalysisTaskJetDynamicalGrooming
 */

#include "AliAnalysisTaskJetDynamicalGrooming.h"

#include <cmath>

#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>

#include <AliAODEvent.h>
#include <AliAODMCHeader.h>
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisManager.h>
#include <AliGenPythiaEventHeader.h>
#include <AliLog.h>
#include <AliMCEvent.h>
#include <AliVCluster.h>
#include <AliVTrack.h>

#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliRhoParameter.h"

using std::cout;
using std::endl;

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming);

namespace PWGJE
{
namespace EMCALJetTasks
{
/**
 * Default constructor.
 */
AliAnalysisTaskJetDynamicalGrooming::AliAnalysisTaskJetDynamicalGrooming()
 : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetDynamicalGrooming", kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized(false),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fOneConstSelectOn(kFALSE),
  fTrackCheckPlots(kFALSE),
  fDoFillMCLund(kFALSE),
  fCheckResolution(kFALSE),
  fSubjetCutoff(0.1),
  fMinPtConst(1),
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fCutDoubleCounts(kTRUE),
  fDoAreaIterative(kTRUE),
  fPowerAlgo(1),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
  fDerivSubtrOrder(0),
  fStoreDetLevelJets(0),
  fPtJet(nullptr),
  fHLundIterative(nullptr),
  fHLundIterativeMC(nullptr),
  fHLundIterativeMCDet(nullptr),
  fHCheckResolutionSubjets(nullptr),
  fTreeSubstructure(0)
{
  for (Int_t i = 0; i < nSubstructureVariables; i++) {
    fShapesVar[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/**
 * Standard constructor.
 */
AliAnalysisTaskJetDynamicalGrooming::AliAnalysisTaskJetDynamicalGrooming(const char* name)
 : AliAnalysisTaskEmcalJet(name, kTRUE),
  fYAMLConfig(),
  fConfigurationInitialized(false),
  fContainer(0),
  fMinFractionShared(0),
  fJetShapeType(kData),
  fJetShapeSub(kNoSub),
  fJetSelection(kInclusive),
  fPtThreshold(-9999.),
  fRMatching(0.2),
  fCentSelectOn(kTRUE),
  fCentMin(0),
  fCentMax(10),
  fOneConstSelectOn(kFALSE),
  fTrackCheckPlots(kFALSE),
  fDoFillMCLund(kFALSE),
  fCheckResolution(kFALSE),
  fSubjetCutoff(0.1),
  fMinPtConst(1),
  fHardCutoff(0),
  fDoTwoTrack(kFALSE),
  fCutDoubleCounts(kTRUE),
  fDoAreaIterative(kTRUE),
  fPowerAlgo(1),
  fPhiCutValue(0.02),
  fEtaCutValue(0.02),
  fMagFieldPolarity(1),
  fDerivSubtrOrder(0),
  fStoreDetLevelJets(0),
  fPtJet(nullptr),
  fHLundIterative(nullptr),
  fHLundIterativeMC(nullptr),
  fHLundIterativeMCDet(nullptr),
  fHCheckResolutionSubjets(nullptr),
  fTreeSubstructure(0)
{
  // Standard constructor.
  for (Int_t i = 0; i < nSubstructureVariables; i++) {
    fShapesVar[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/**
 * Copy constructor
 * TODO: Update this!
 */
AliAnalysisTaskJetDynamicalGrooming::AliAnalysisTaskJetDynamicalGrooming(
 const AliAnalysisTaskJetDynamicalGrooming& other)
 : fYAMLConfig(other.fYAMLConfig), fConfigurationInitialized(other.fConfigurationInitialized)
{
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
AliAnalysisTaskJetDynamicalGrooming& AliAnalysisTaskJetDynamicalGrooming::operator=(
 AliAnalysisTaskJetDynamicalGrooming other)
{
  swap(*this, other);
  return *this;
}

/**
 * Initialize task.
 */
bool AliAnalysisTaskJetDynamicalGrooming::Initialize()
{
  fConfigurationInitialized = false;

  // Ensure that we have at least one configuration in the YAML config.
  /*if (fYAMLConfig.DoesConfigurationExist(0) == false) {
    // No configurations exist. Return immediately.
    return fConfigurationInitialized;
  }*/

  // Always initialize for streaming purposes
  fYAMLConfig.Initialize();

  // Setup task based on the properties defined in the YAML config
  /*AliDebugStream(2) << "Configuring task from the YAML configuration.\n";
  RetrieveAndSetTaskPropertiesFromYAMLConfig();
  SetupJetContainersFromYAMLConfig();
  SetupParticleContainersFromYAMLConfig();
  SetupClusterContainersFromYAMLConfig();
  AliDebugStream(2) << "Finished configuring via the YAML configuration.\n";*/

  // Print the results of the initialization
  // Print outside of the ALICE Log system to ensure that it is always available!
  std::cout << *this;

  fConfigurationInitialized = true;
  return fConfigurationInitialized;
}

/**
 *
 */
void AliAnalysisTaskJetDynamicalGrooming::SetupTree()
{
  const char* nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeSubstructure = new TTree(nameoutput, nameoutput);

  std::vector<std::string> branchNames = {
    "ptJet",
    "ktg",
    "kTLeading",
    "kTLeadingNSplittings",
    "kTDynamical",
    "kTDynamicalNSplittings",
    "ng",
    "zg",
    "rg",
    "ptJetMatch",
    "ktgMatch",
    "kTLeadingMatch",
    "kTLeadingMatchNSplittings",
    "kTDynamicalMatch",
    "kTDynamicalMatchNSplittings",
    "ngMatch",
    "zgMatch",
    "rgMatch",
    "LeadingTrackPt",
    "LeadingTrackPtMatch",
  };
  if (fStoreDetLevelJets) {
    std::vector<std::string> detectorLevelBranchNames = {
      "ptJetDet",
      "ktgDet",
      "kTLeadingDet",
      "kTLeadingDetNSplittings"
      "kTDynamicalDet",
      "kTDynamicalDetNSplittings"
      "ngDet",
      "zgDet",
      "rgDet",
      "LeadingTrackPtDet"
    };
    for (auto name : detectorLevelBranchNames) {
      branchNames.emplace_back(name);
    }
  }

  /*
  const Int_t nVar = nSubstructureVariables;
  TString* fShapesVarNames = new TString[nVar];

  fShapesVarNames[0] = "ptJet";
  fShapesVarNames[1] = "ktg";
  fShapesVarNames[2] = "ng";
  fShapesVarNames[3] = "zg";
  fShapesVarNames[4] = "rg";
  fShapesVarNames[5] = "ptJetMatch";
  fShapesVarNames[6] = "ktgMatch";
  fShapesVarNames[7] = "ngMatch";
  fShapesVarNames[8] = "zgMatch";
  fShapesVarNames[9] = "rgMatch";
  fShapesVarNames[10] = "LeadingTrackPt";
  fShapesVarNames[11] = "LeadingTrackPtMatch";
  if (fStoreDetLevelJets) {
    fShapesVarNames[12] = "ptJetDet";
    fShapesVarNames[13] = "ktgDet";
    fShapesVarNames[14] = "ngDet";
    fShapesVarNames[15] = "zgDet";
    fShapesVarNames[16] = "rgDet";
    fShapesVarNames[17] = "LeadingTrackPtDet";
  }
  */

  // for (Int_t ivar = 0; ivar < nVar; ivar++) {
  for (std::size_t i = 0; i < branchNames.size(); i++) {
    cout << "looping over variables" << endl;
    fTreeSubstructure->Branch(branchNames[i].c_str(), &fShapesVar[i],
                 TString::Format("%s/F", branchNames[i].c_str()));
  }

  // delete[] fShapesVarNames;
}

/**
 * Create output objects.
 */
void AliAnalysisTaskJetDynamicalGrooming::UserCreateOutputObjects()
{
  // First call the base class
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Check that the task was initialized
  if (!fConfigurationInitialized) {
    AliFatal("Task was not initialized. Please ensure that Initialize() was called!");
  }
  // Reinitialize the YAML configuration
  fYAMLConfig.Reinitialize();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fPtJet = new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec = 8;
  const Int_t nBinsSpec[8] = { 50, 100, 100, 20, 100, 50, 100, 2 };
  const Double_t lowBinSpec[8] = { 0., -5, 0, 0, 0, 0, 0, 0 };
  const Double_t hiBinSpec[8] = { 5., 10., 200, 20, 200, 50, 50, 2 };
  fHLundIterative = new THnSparseF("fHLundIterative", "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
                   dimSpec, nBinsSpec, lowBinSpec, hiBinSpec);
  fOutput->Add(fHLundIterative);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec2 = 7;
  const Int_t nBinsSpec2[7] = { 50, 100, 100, 20, 100, 50, 100 };
  const Double_t lowBinSpec2[7] = { 0., -5, 0, 0, 0, 0, 0 };
  const Double_t hiBinSpec2[7] = { 5., 10., 200, 20, 200, 100, 50 };
  fHLundIterativeMC =
   new THnSparseF("fHLundIterativeMC", "LundIterativePlotMC [log(1/theta),log(z*theta),pTjet,algo]", dimSpec2,
           nBinsSpec2, lowBinSpec2, hiBinSpec2);
  fOutput->Add(fHLundIterativeMC);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec3 = 7;
  const Int_t nBinsSpec3[7] = { 50, 100, 100, 20, 100, 50, 100 };
  const Double_t lowBinSpec3[7] = { 0., -5, 0, 0, 0, 0, 0 };
  const Double_t hiBinSpec3[7] = { 5., 10., 200, 20, 200, 100, 50 };
  fHLundIterativeMCDet =
   new THnSparseF("fHLundIterativeMCDet", "LundIterativePlotMCDet [log(1/theta),log(z*theta),pTjet,algo]", dimSpec3,
           nBinsSpec3, lowBinSpec3, hiBinSpec3);
  fOutput->Add(fHLundIterativeMCDet);

  ////
  const Int_t dimResol = 5;
  const Int_t nBinsResol[5] = { 10, 10, 80, 80, 80 };
  const Double_t lowBinResol[5] = { 0, 0, -1, -1, -1 };
  const Double_t hiBinResol[5] = { 200, 0.3, 1, 1, 1 };
  fHCheckResolutionSubjets = new THnSparseF("fHCheckResolutionSubjets", "Mom.Resolution of Subjets vs opening angle",
                       dimResol, nBinsResol, lowBinResol, hiBinResol);
  fOutput->Add(fHCheckResolutionSubjets);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i = 0; i < fOutput->GetEntries(); ++i) {
    TH1* h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1) {
      h1->Sumw2();
      continue;
    }
    THnSparse* hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if (hn)
      hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  // Setup the tree output.
  SetupTree();

  PostData(1, fOutput);
  PostData(2, fTreeSubstructure);
}

Bool_t AliAnalysisTaskJetDynamicalGrooming::Run()
{
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

double AliAnalysisTaskJetDynamicalGrooming::DynamicalGrooming(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R,
                               const double a) const
{
  double deltaR = subjet1.delta_R(subjet2);
  double z = subjet1.pt() / parent.pt();
  return z * (1 - z) * parent.pt() * std::pow(deltaR / R, a);
}

double AliAnalysisTaskJetDynamicalGrooming::CalculateZDrop(const fastjet::PseudoJet& subjet1,
                              const fastjet::PseudoJet& subjet2,
                              const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 0.1);
}

double AliAnalysisTaskJetDynamicalGrooming::CalculateKtDrop(const fastjet::PseudoJet& subjet1,
                              const fastjet::PseudoJet& subjet2,
                              const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 1);
}

double AliAnalysisTaskJetDynamicalGrooming::CalculateTimeDrop(const fastjet::PseudoJet& subjet1,
                               const fastjet::PseudoJet& subjet2,
                               const fastjet::PseudoJet& parent, const double R) const
{
  return DynamicalGrooming(subjet1, subjet2, parent, R, 2);
}

Bool_t AliAnalysisTaskJetDynamicalGrooming::FillHistograms()
{
  AliEmcalJet* jet1 = NULL;
  AliJetContainer* jetCont = GetJetContainer(0);
  // container zero is always the base container: the data container, the
  // embedded subtracted in the case of embedding or the detector level in case
  // of pythia

  if (fCentSelectOn)
    if ((fCent > fCentMax) || (fCent < fCentMin))
      return 0;

  Float_t rhoVal = 0, rhoMassVal = 0.;
  if (jetCont) {
    jetCont->ResetCurrentID();
    if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kDerivSub)) {
      // rho
      AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoSparseR020"));
      if (!rhoParam) {
        Printf(
         "%s: Could not retrieve rho %s (some histograms will be filled "
         "with zero)!",
         GetName(), jetCont->GetRhoName().Data());
      } else
        rhoVal = rhoParam->GetVal();
      // rhom
      AliRhoParameter* rhomParam =
       dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject("RhoMassSparseR020"));
      if (!rhomParam) {
        Printf(
         "%s: Could not retrieve rho_m %s (some histograms will be "
         "filled with zero)!",
         GetName(), jetCont->GetRhoMassName().Data());
      } else
        rhoMassVal = rhomParam->GetVal();
    }

    while ((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1)
        continue;
      AliEmcalJet* jet2 = nullptr;
      AliEmcalJet* jet3 = nullptr;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet* jetUS = NULL;
      Int_t ifound = 0, jfound = 0;
      Int_t ilab = -1, jlab = -1;

      // The embedding mode
      // the matching is done between unsubtracted embedded jets and detector
      // level jets unsubtracted and subtracted jets share the label. Once we
      // identify the corresponding unsubtracted jet, jetUS, then we fetch jet2,
      // which is the matched detector level jet In the case we are not
      // considering constituent subtraction, then the detector-level matched jet
      // is the one that was directly matched to the base jet1. Then, the
      // particle-level jet jet3 is obtained as the matched one to jet2 In short,
      // there are 2 consecutive matchinges, between particle-level (jet3) and
      // detector-level (jet2) pythia jets and between jet2 and the embedding
      // unsubtracted jet. Note that the matching obtained via ClosestJet is
      // purely geometrical. So below we require a fraction of the probe momentum
      // to be reconstructed in the embedded jet.
      if (fJetShapeType == kDetEmbPartPythia) {
        AliJetContainer* jetContUS = GetJetContainer(2);

        if (fJetShapeSub == kConstSub) {
          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1)
                ilab = i;
            }
          }
          if (ilab == -1)
            continue;
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();
        }

        if (fJetShapeSub == kEventSub) {
          jetUS = jet1->ClosestJet();
          if (!jetUS)
            continue;
          jet2 = jetUS->ClosestJet();
        }

        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub == kEventSub))
          jet2 = jet1->ClosestJet();

        if (!jet2) {
          Printf("jet2 does not exist, returning");
          continue;
        }

        // AliJetContainer *jetContPart=GetJetContainer(3);
        jet3 = jet2->ClosestJet();

        if (!jet3) {
          Printf("jet3 does not exist, returning");
          continue;
        }
        cout << "jet 3 exists" << jet3->Pt() << endl;

        AliJetContainer* jetContTrue = GetJetContainer(1);
        AliJetContainer* jetContPart = GetJetContainer(3);

        if (fCheckResolution)
          CheckSubjetResolution(jet2, jetContTrue, jet3, jetContPart);

        Double_t fraction = 0;
        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub == kEventSub))
          fraction = jetCont->GetFractionSharedPt(jet1);
        if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kEventSub))
          fraction = jetContUS->GetFractionSharedPt(jetUS);

        if (fraction < fMinFractionShared)
          continue;
      }

      // this is the mode to run over pythia to produce a det-part response
      // here we have also added the constituent-subtraction case, but we don't
      // use it normally in pp the matching is purely geometrical
      if (fJetShapeType == kPythiaDef) {
        AliJetContainer* jetContTrue = GetJetContainer(1);
        AliJetContainer* jetContUS = GetJetContainer(2);
        AliJetContainer* jetContPart = GetJetContainer(3);

        if (fJetShapeSub == kConstSub) {
          for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
            jetUS = jetContUS->GetJet(i);
            if (jetUS->GetLabel() == jet1->GetLabel()) {
              ifound++;
              if (ifound == 1)
                ilab = i;
            }
          }
          if (ilab == -1)
            continue;
          jetUS = jetContUS->GetJet(ilab);
          jet2 = jetUS->ClosestJet();

          if (!jet2) {
            Printf("jet2 does not exist, returning");
            continue;
          }

          for (Int_t j = 0; j < jetContPart->GetNJets(); j++) {
            jet3 = jetContPart->GetJet(j);
            if (!jet3)
              continue;
            if (jet3->GetLabel() == jet2->GetLabel()) {
              jfound++;
              if (jfound == 1)
                jlab = j;
            }
          }
          if (jlab == -1)
            continue;
          jet3 = jetContPart->GetJet(jlab);
          if (!jet3) {
            Printf("jet3 does not exist, returning");
            continue;
          }
        }
        if (!(fJetShapeSub == kConstSub))
          jet3 = jet1->ClosestJet();
        if (!jet3) {
          Printf("jet3 does not exist, returning");
          continue;
        }

        if (fCheckResolution)
          CheckSubjetResolution(jet1, jetCont, jet3, jetContTrue);
      }

      Double_t ptSubtracted = 0;
      if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub)
        ptSubtracted = jet1->Pt();

      else if (fJetShapeSub == kDerivSub) {
        ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
      }

      else if (fJetShapeSub == kNoSub) {
        if ((fJetShapeType == kData) || (fJetShapeType == kDetEmbPartPythia))
          ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
        else if ((fJetShapeType == kPythiaDef) || (fJetShapeType == kMCTrue) || (fJetShapeType == kGenOnTheFly))
          ptSubtracted = jet1->Pt();
      }

      if (ptSubtracted < fPtThreshold)
        continue;

      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1))
        continue;

      fShapesVar[0] = ptSubtracted;
      fShapesVar[18] = jet1->MaxTrackPt();

      if (fCutDoubleCounts == kTRUE && fJetShapeType == kDetEmbPartPythia)
        if (jet1->MaxTrackPt() > jet3->MaxTrackPt())
          continue;

      IterativeParents(jet1, jetCont);

      Float_t ptMatch = 0.;
      Float_t leadTrackMatch = 0.;
      Double_t ktgMatch = 0;
      Double_t kTLeadingMatch = 0;
      Double_t kTLeadingMatchNSplittings = 0;
      Double_t kTDynamicalMatch = 0;
      Double_t kTDynamicalMatchNSplittings = 0;

      Double_t nsdMatch = 0;
      Double_t zgMatch = 0;
      Double_t rgMatch = 0;
      Float_t ptDet = 0.;
      Float_t leadTrackDet = 0.;
      Double_t ktgDet = 0;
      Double_t kTLeadingDet = 0;
      Double_t kTLeadingDetNSplittings = 0;
      Double_t kTDynamicalDet = 0;
      Double_t kTDynamicalDetNSplittings = 0;

      Double_t nsdDet = 0;
      Double_t zgDet = 0;
      Double_t rgDet = 0;
      Double_t aver1 = 0;
      Double_t aver2 = 0;
      Double_t aver3 = 0;
      Double_t aver4 = 0;
      Double_t aver5 = 0;
      Double_t aver6 = 0;
      Double_t aver7 = 0;
      Double_t aver8 = 0;
      Int_t kMatched = 0;
      if (fJetShapeType == kPythiaDef) {
        kMatched = 1;
        if (fJetShapeSub == kConstSub)
          kMatched = 3;

        ptMatch = jet3->Pt();
        leadTrackMatch = jet3->MaxTrackPt();
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4, aver5, aver6, aver7, aver8);
        ktgMatch = aver1;
        kTLeadingMatch = aver2;
        kTLeadingMatchNSplittings = aver3;
        kTDynamicalMatch = aver4;
        kTDynamicalMatchNSplittings = aver5;
        nsdMatch = aver6;
        zgMatch = aver7;
        rgMatch = aver8;
      }

      if (fJetShapeType == kDetEmbPartPythia) {
        if (fJetShapeSub == kConstSub)
          kMatched = 3;
        if (fJetShapeSub == kDerivSub)
          kMatched = 2;
        ptMatch = jet3->Pt();
        leadTrackMatch = jet3->MaxTrackPt();
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4, aver5, aver6, aver7, aver8);
        ktgMatch = aver1;
        kTLeadingMatch = aver2;
        kTLeadingMatchNSplittings = aver3;
        kTDynamicalMatch = aver4;
        kTDynamicalMatchNSplittings = aver5;
        nsdMatch = aver6;
        zgMatch = aver7;
        rgMatch = aver8;
        if (fStoreDetLevelJets) {
          ptDet = jet2->Pt();
          leadTrackDet = jet2->MaxTrackPt();
          IterativeParentsMCAverage(jet2, 1, ktgDet, kTLeadingDet, kTLeadingDetNSplittings, kTDynamicalDet, kTDynamicalDetNSplittings, nsdDet, zgDet, rgDet);
        }
      }

      if (fJetShapeType == kMCTrue || fJetShapeType == kData || fJetShapeType == kGenOnTheFly) {
        ptMatch = 0.;
        leadTrackMatch = 0.;
        ktgMatch = 0.;
        kTLeadingMatch = 0;
        kTLeadingMatchNSplittings = 0;
        kTDynamicalMatch = 0;
        kTDynamicalMatchNSplittings = 0;
        nsdMatch = 0.;
        zgMatch = 0;
        rgMatch = 0;
      }

      /*{ 0: "ptJet",
       1: "ktg",
       2: "kTLeading",
       3: "kTLeadingNSplittings"
       4: "kTDynamical",
       5: "kTDynamicalNSplittings"
       6: "ng",
       7: "zg",
       8: "rg",
       9: "ptJetMatch",
       10: "ktgMatch",
       11: "kTLeadingMatch",
       12: "kTLeadingMatchNSplittings",
       13: "kTDynamicalMatch",
       14: "kTDynamicalMatchNSplittings",
       15: "ngMatch",
       16: "zgMatch",
       17: "rgMatch",
       18: "LeadingTrackPt",
       19: "LeadingTrackPtMatch",
       20: "ptJetDet",
       21: "ktgDet",
       22: "kTLeadingDet",
       23: "kTLeadingDetNSplittings",
       24: "kTDynamicalDet",
       25: "kTDynamicalDetNSplittings",
       26: "ngDet",
       27: "zgDet",
       28: "rgDet",
       29: "LeadingTrackPtDet"
      };*/

      fShapesVar[9] = ptMatch;
      fShapesVar[10] = ktgMatch;
      fShapesVar[11] = kTLeadingMatch;
      fShapesVar[12] = kTLeadingMatchNSplittings;
      fShapesVar[13] = kTDynamicalMatch;
      fShapesVar[14] = kTDynamicalMatchNSplittings;
      fShapesVar[15] = nsdMatch;
      fShapesVar[16] = zgMatch;
      fShapesVar[17] = rgMatch;
      fShapesVar[19] = leadTrackMatch;
      if (fStoreDetLevelJets) {
        fShapesVar[20] = ptDet;
        fShapesVar[21] = ktgDet;
        fShapesVar[22] = kTLeadingDet;
        fShapesVar[23] = kTLeadingDetNSplittings;
        fShapesVar[24] = kTDynamicalDet;
        fShapesVar[25] = kTDynamicalDetNSplittings;
        fShapesVar[26] = nsdDet;
        fShapesVar[27] = zgDet;
        fShapesVar[28] = rgDet;
        fShapesVar[29] = leadTrackDet;
      }

      fTreeSubstructure->Fill();
    }
  }

  return kTRUE;
}

Float_t AliAnalysisTaskJetDynamicalGrooming::GetJetMass(AliEmcalJet* jet, Int_t jetContNb = 0)
{
  // calc subtracted jet mass
  if ((fJetShapeSub == kDerivSub) && (jetContNb == 0))
    if (fDerivSubtrOrder == 1)
      return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
      return jet->GetShapeProperties()->GetSecondOrderSubtracted();
  else
    return jet->M();
}

Float_t AliAnalysisTaskJetDynamicalGrooming::Angularity(AliEmcalJet* jet, Int_t jetContNb = 0)
{
  AliJetContainer* jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den = 0.;
  Double_t num = 0.;
  AliVParticle* vp1 = nullptr;
  for (UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));

    if (!vp1) {
      Printf("AliVParticle associated to constituent not found");
      continue;
    }

    Double_t dphi = RelativePhi(vp1->Phi(), jet->Phi());
    Double_t dr2 = (vp1->Eta() - jet->Eta()) * (vp1->Eta() - jet->Eta()) + dphi * dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num = num + vp1->Pt() * dr;
    den = den + vp1->Pt();
  }
  return num / den;
}

Float_t AliAnalysisTaskJetDynamicalGrooming::GetJetAngularity(AliEmcalJet* jet, Int_t jetContNb = 0)
{
  if ((fJetShapeSub == kDerivSub) && (jetContNb == 0))
    if (fDerivSubtrOrder == 1)
      return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
    else
      return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
}

Double_t AliAnalysisTaskJetDynamicalGrooming::RelativePhi(Double_t mphi, Double_t vphi)
{
  if (vphi < -1 * TMath::Pi())
    vphi += (2 * TMath::Pi());
  else if (vphi > TMath::Pi())
    vphi -= (2 * TMath::Pi());
  if (mphi < -1 * TMath::Pi())
    mphi += (2 * TMath::Pi());
  else if (mphi > TMath::Pi())
    mphi -= (2 * TMath::Pi());
  double dphi = mphi - vphi;
  if (dphi < -1 * TMath::Pi())
    dphi += (2 * TMath::Pi());
  else if (dphi > TMath::Pi())
    dphi -= (2 * TMath::Pi());
  return dphi; // dphi in [-Pi, Pi]
}

/*
void AliAnalysisTaskJetDynamicalGrooming::IterativeParentsAreaBased(AliEmcalJet* fJet, AliJetContainer* fJetCont)
{
  // to still change and implement the 4 vector bkg subtraction to the subjets
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  AliParticleContainer* fTrackCont = fJetCont->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle* fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk)
        continue;
      if (fDoTwoTrack == kTRUE && CheckClosePartner(i, fJet, fTrk, fTrackCont))
        continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
      fInputVectors.push_back(PseudoTracks);
    }
  fastjet::JetAlgorithm jetalgo(fastjet::genkt_algorithm);
  fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);

  fastjet::JetDefinition fJetDef(jetalgo, 1., fPowerAlgo, static_cast<fastjet::RecombinationScheme>(0),
                  fastjet::BestFJ30);
  fastjet::AreaDefinition fAreaDef(fastjet::passive_area, ghost_spec);
  try {
    fastjet::ClusterSequenceArea fClustSeqSA(fInputVectors, fJetDef, fAreaDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets = fClustSeqSA.inclusive_jets(0);

    fastjet::PseudoJet jj;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    jj = fOutputJets[0];

    double nall = 0;
    double nsd = 0;
    int flagSubjet = 0;
    int flagSubjetkT = 0;
    double Rg = 0;
    double zg = 0;
    double xktg = 0;
    double z = 0;
    double cumtf = 0;
    fastjet::PseudoJet area1, area2;

    while (jj.has_parents(j1, j2) && z < fHardCutoff) {
      nall = nall + 1;

      flagSubjet = 0;
      flagSubjetkT = 0;
      area1 = j1.area_4vector();
      area2 = j2.area_4vector();
      fastjet::PseudoJet jet_sub1 = j1 - GetRhoVal(0) * area1;
      fastjet::PseudoJet jet_sub2 = j2 - GetRhoVal(0) * area2;

      if (jet_sub1.perp() < jet_sub2.perp())
        swap(jet_sub1, jet_sub2);
      if (jet_sub1.perp() < 0 && jet_sub2.perp() < 0)
        break;

      if (jet_sub2.perp() > 0) {
        double delta_R = jet_sub1.delta_R(jet_sub2);
        double xkt = jet_sub2.perp() * sin(delta_R);
        double lnpt_rel = log(xkt);
        double y = log(1. / delta_R);
        double form = 2 * 0.197 * jet_sub2.e() / (xkt * xkt);
        double rad = jet_sub2.e();

        z = jet_sub2.perp() / (jet_sub1.perp() + jet_sub2.perp());

        if (z > fHardCutoff)
          nsd = nsd + 1;
        if (z > fHardCutoff && flagSubjet == 0) {
          zg = z;
          xktg = xkt;
          Rg = delta_R;
          flagSubjet = 1;
        }
        if (lnpt_rel > 0) {
          cumtf = cumtf + form;
          if ((nsd == 0) && (flagSubjetkT == 0)) {
            xktg = xkt;
            flagSubjetkT = 1;
          }
        }
        Double_t LundEntries[7] = { y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, cumtf };
        fHLundIterative->Fill(LundEntries);
      }

      jj = jet_sub1;
    }

    fShapesVar[1] = xktg;
    fShapesVar[2] = kTLeading;
    fShapesVar[3] = kTDynamical;
    fShapesVar[4] = nsd;
    fShapesVar[5] = zg;
    fShapesVar[6] = Rg;

  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}*/

void AliAnalysisTaskJetDynamicalGrooming::IterativeParents(AliEmcalJet* fJet, AliJetContainer* fJetCont)
{
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  AliParticleContainer* fTrackCont = fJetCont->GetParticleContainer();
  cout << "is it really here?" << endl;
  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle* fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk)
        continue;
      if (fDoTwoTrack == kTRUE && CheckClosePartner(i, fJet, fTrk, fTrackCont))
        continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
      fInputVectors.push_back(PseudoTracks);
    }
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);

  fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);
  // fastjet::JetAlgorithm jetalgo(fastjet::genkt_algorithm);
  // fastjet::JetDefinition fJetDef(jetalgo, 1., fPowerAlgo,
  //                              static_cast<fastjet::RecombinationScheme>(0),
  //                             fastjet::BestFJ30);
  fastjet::AreaDefinition fAreaDef(fastjet::passive_area, ghost_spec);
  try {
    fastjet::ClusterSequenceArea fClustSeqSA(fInputVectors, fJetDef, fAreaDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets = fClustSeqSA.inclusive_jets(0);

    fastjet::PseudoJet jj;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    jj = fOutputJets[0];

    double nall = 0;
    double nsd = 0;
    int flagSubjet = 0;
    int flagSubjetkT = 0;
    double flagConst = 0;
    double Rg = 0;
    double zg = 0;
    double xktg = 0;
    double kTLeading = 0;
    double kTDynamical = 0;
    double cumtf = 0;

    int kTLeadingNSplittings = 0;
    int kTDynamicalNSplittings = 0;
    while (jj.has_parents(j1, j2)) {
      // Keep track of how many splits we have checked.
      nall = nall + 1;

      if (j1.perp() < j2.perp())
        swap(j1, j2);
      flagConst = 0;
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      double dynamicalKtGrooming = CalculateKtDrop(j1, j2, jj, fJetCont->GetJetRadius());
      double lnpt_rel = log(xkt);
      double y = log(1. / delta_R);
      double form = 2 * 0.197 * j2.e() / (xkt * xkt);
      double rad = j1.e() + j2.e();
      double z = j2.perp() / (j2.perp() + j1.perp());
      vector<fastjet::PseudoJet> constitj1 = sorted_by_pt(j1.constituents());
      if (constitj1[0].perp() > fMinPtConst)
        flagConst = 1;

      // Don't apply the fHardCutoff to dynamical grooming. It takes care of this itself.
      if (dynamicalKtGrooming > kTDynamical) {
        kTDynamical = dynamicalKtGrooming;
        kTDynamicalNSplittings = nall;
      }
      if (z > fHardCutoff) {
        nsd = nsd + 1;
        if (xkt > kTLeading) {
          kTLeading = xkt;
          // Used nsd instead of the splitting index because we want to know how mnay
          // passed the hard cutoff cnodition.
          kTLeadingNSplittings = nsd;
        }
      }
      if (z > fHardCutoff && flagSubjet == 0) {
        zg = z;
        xktg = xkt;
        Rg = delta_R;
        flagSubjet = 1;
      }
      if (lnpt_rel > 0) {
        cumtf = cumtf + form;
        if ((nsd == 0) && (flagSubjetkT == 0)) {
          xktg = xkt;
          flagSubjetkT = 1;
        }
      }

      Double_t LundEntries[8] = { y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, cumtf, flagConst };
      fHLundIterative->Fill(LundEntries);

      jj = j1;
    }

    fShapesVar[1] = xktg;
    fShapesVar[2] = kTLeading;
    fShapesVar[3] = static_cast<double>(kTLeadingNSplittings);
    fShapesVar[4] = kTDynamical;
    fShapesVar[5] = static_cast<double>(kTDynamicalNSplittings);
    fShapesVar[6] = nsd;
    fShapesVar[7] = zg;
    fShapesVar[8] = Rg;

  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

void AliAnalysisTaskJetDynamicalGrooming::IterativeParentsMCAverage(AliEmcalJet* fJet, Int_t km, Double_t& average1,
                                  Double_t& average2, Double_t& average3,
                                  Double_t& average4, Double_t& average5,
                                  Double_t& average6, Double_t& average7,
                                  Double_t& average8)
{
  AliJetContainer* jetCont = GetJetContainer(km);
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  AliParticleContainer* fTrackCont = jetCont->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle* fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk)
        continue;

      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
      fInputVectors.push_back(PseudoTracks);
    }
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);

  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);

  try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets = fClustSeqSA.inclusive_jets(0);

    fastjet::PseudoJet jj;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    jj = fOutputJets[0];
    int flagSubjet = 0;
    int flagSubjetkT = 0;
    double nall = 0;
    double nsd = 0;

    double zg = 0;
    double xktg = 0;
    double kTLeading = 0;
    double kTDynamical = 0;
    double Rg = 0;

    int kTLeadingNSplittings = 0;
    int kTDynamicalNSplittings = 0;

    double cumtf = 0;
    while (jj.has_parents(j1, j2)) {
      // Keep track of how many splits we have checked.
      nall = nall + 1;

      if (j1.perp() < j2.perp())
        swap(j1, j2);
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      double dynamicalKtGrooming = CalculateKtDrop(j1, j2, jj, jetCont->GetJetRadius());
      double lnpt_rel = log(xkt);
      double y = log(1. / delta_R);
      double form = 2 * 0.197 * j2.e() / (xkt * xkt);
      double rad = j1.e() + j2.e();
      double z = j2.perp() / (j2.perp() + j1.perp());

      // Don't apply the fHardCutoff to dynamical grooming. It takes care of this itself.
      if (dynamicalKtGrooming > kTDynamical) {
        kTDynamical = dynamicalKtGrooming;
        kTDynamicalNSplittings = nall;
      }
      if (z > fHardCutoff) {
        nsd = nsd + 1;
        if (xkt > kTLeading) {
          kTLeading = xkt;
          // Used nsd instead of the splitting index because we want to know how mnay
          // passed the hard cutoff cnodition.
          kTLeadingNSplittings = nsd;
        }
      }
      if (z > fHardCutoff && flagSubjet == 0) {
        zg = z;
        xktg = xkt;
        Rg = delta_R;
        flagSubjet = 1;
      }
      if (lnpt_rel > 0) {
        cumtf = cumtf + form;
        if ((nsd == 0) && (flagSubjetkT == 0)) {
          xktg = xkt;
          flagSubjetkT = 1;
        }
      }
      if (fDoFillMCLund == kTRUE) {
        Double_t LundEntries[7] = { y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, cumtf };
        fHLundIterativeMC->Fill(LundEntries);
        if (fStoreDetLevelJets) {
          fHLundIterativeMCDet->Fill(LundEntries);
        }
      }

      jj = j1;
    }

    average1 = xktg;
    average2 = kTLeading;
    average3 = static_cast<double>(kTLeadingNSplittings);
    average4 = kTDynamical;
    average5 = static_cast<double>(kTDynamicalNSplittings);
    average6 = nsd;
    average7 = zg;
    average8 = Rg;

  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

void AliAnalysisTaskJetDynamicalGrooming::CheckSubjetResolution(AliEmcalJet* fJet, AliJetContainer* fJetCont,
                                AliEmcalJet* fJetM, AliJetContainer* fJetContM)
{
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  std::vector<fastjet::PseudoJet> fInputVectorsM;
  fInputVectorsM.clear();
  fastjet::PseudoJet PseudoTracksM;

  AliParticleContainer* fTrackCont = fJetCont->GetParticleContainer();
  AliParticleContainer* fTrackContM = fJetContM->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle* fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk)
        continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
      fInputVectors.push_back(PseudoTracks);
    }
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);

  if (fTrackContM)
    for (Int_t i = 0; i < fJetM->GetNumberOfTracks(); i++) {
      AliVParticle* fTrk = fJetM->TrackAt(i, fTrackContM->GetArray());
      if (!fTrk)
        continue;
      PseudoTracksM.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracksM.set_user_index(fJetM->TrackAt(i) + 100);
      fInputVectorsM.push_back(PseudoTracksM);
    }
  fastjet::JetAlgorithm jetalgoM(fastjet::cambridge_algorithm);
  fastjet::JetDefinition fJetDefM(jetalgoM, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);

  try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets = fClustSeqSA.inclusive_jets(0);

    fastjet::ClusterSequence fClustSeqSAM(fInputVectorsM, fJetDefM);
    std::vector<fastjet::PseudoJet> fOutputJetsM;
    fOutputJetsM.clear();
    fOutputJetsM = fClustSeqSAM.inclusive_jets(0);

    fastjet::PseudoJet jj, jjM;
    fastjet::PseudoJet j1, j1M;
    fastjet::PseudoJet j2, j2M;
    jj = fOutputJets[0];
    jjM = fOutputJetsM[0];

    double z1 = 0;
    double z2 = 0;
    double zcut = 0.1;
    while ((jj.has_parents(j1, j2)) && (z1 < zcut)) {
      if (j1.perp() < j2.perp())
        swap(j1, j2);

      z1 = j2.perp() / (j1.perp() + j2.perp());
      jj = j1;
    }
    if (z1 < zcut)
      return;

    while ((jjM.has_parents(j1M, j2M)) && (z2 < zcut)) {
      if (j1M.perp() < j2M.perp())
        swap(j1M, j2M);

      z2 = j2M.perp() / (j1M.perp() + j2M.perp());
      jjM = j1M;
    }
    if (z2 < zcut)
      return;

    double delta_R1 = j1.delta_R(j1M);
    double delta_R2 = j2.delta_R(j2M);
    double delta_R = j1.delta_R(j2);
    double residz = (z1 - z2) / z2;
    double resid1 = (j1.perp() - j1M.perp()) / j1M.perp();
    double resid2 = (j2.perp() - j2M.perp()) / j2M.perp();

    if ((delta_R1 < fSubjetCutoff) && (delta_R2 < fSubjetCutoff)) {
      Double_t ResolEntries[5] = { fOutputJets[0].perp(), delta_R, resid1, resid2, residz };
      fHCheckResolutionSubjets->Fill(ResolEntries);
    }

  } catch (const fastjet::Error&) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

Bool_t AliAnalysisTaskJetDynamicalGrooming::CheckClosePartner(Int_t index, AliEmcalJet* fJet, AliVParticle* fTrk1,
                               AliParticleContainer* fTrackCont)
{
  // check if tracks are close//
  for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
    AliVParticle* fTrk2 = fJet->TrackAt(i, fTrackCont->GetArray());
    if (!fTrk2)
      continue;
    if (i == index)
      continue;
    Double_t phi1 = fTrk1->Phi();
    Double_t phi2 = fTrk2->Phi();
    Double_t chg1 = fTrk1->Charge();
    Double_t chg2 = fTrk2->Charge();
    Double_t ptv1 = fTrk1->Pt();
    Double_t ptv2 = fTrk2->Pt();
    Double_t deta = fTrk2->Eta() - fTrk1->Eta();
    const Float_t kLimit = fPhiCutValue * 3;

    if (TMath::Abs(fTrk1->Eta() - fTrk2->Eta()) < fEtaCutValue * 2.5 * 3) {
      Float_t initdpsinner = (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * 0.8 / ptv2) -
                  (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * 0.8 / ptv1)));

      Float_t initdpsouter = (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * 2.5 / ptv2) -
                  (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * 2.5 / ptv1)));

      initdpsinner = TVector2::Phi_mpi_pi(initdpsinner);
      initdpsouter = TVector2::Phi_mpi_pi(initdpsouter);

      if (TMath::Abs(initdpsinner) < kLimit || TMath::Abs(initdpsouter) < kLimit ||
        initdpsinner * initdpsouter < 0) {
        Double_t mindps = 1e5;

        for (Double_t rad = 0.8; rad < 2.51; rad += 0.01) {
          Double_t dps = (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * rad / ptv2) -
                  (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * rad / ptv1)));
          dps = TVector2::Phi_mpi_pi(dps);
          if (TMath::Abs(dps) < TMath::Abs(mindps))
            mindps = dps;
        }
        if (TMath::Abs(mindps) < fPhiCutValue && TMath::Abs(deta) < fEtaCutValue)
          return kTRUE;
      }
    }
  }
  return kFALSE;
}

Bool_t AliAnalysisTaskJetDynamicalGrooming::RetrieveEventObjects()
{
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

/**
 * Add task to an existing analysis manager.
 *
 * @param[in] suffix Suffix string to attach to the task name
 */
AliAnalysisTaskJetDynamicalGrooming* AliAnalysisTaskJetDynamicalGrooming::AddTaskJetDynamicalGrooming(
 const char* njetsBase, const char* njetsUS, const char* njetsTrue, const char* njetsPartLevel, const Double_t R,
 const char* nrhoBase, const char* ntracks, const char* ntracksUS, const char* ntracksPartLevel, const char* nclusters,
 const char* ntracksTrue, const char* type, const char* CentEst, Int_t pSel,
 AliAnalysisTaskJetDynamicalGrooming::JetShapeType_t jetShapeType,
 AliAnalysisTaskJetDynamicalGrooming::JetShapeSub_t jetShapeSub,
 AliAnalysisTaskJetDynamicalGrooming::JetSelectionType_t jetSelection,
 Float_t minpTHTrigger, Float_t maxpTHTrigger, Float_t acut,
 AliAnalysisTaskJetDynamicalGrooming::DerivSubtrOrder_t derivSubtrOrder,
 const std::string& suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliErrorClass("No analysis manager to connect to.");
    return nullptr;
  }

  // Setup task name
  std::string taskName = "AliAnalysisTaskJetDynamicalGrooming";
  std::string suffixName(suffix);
  if (suffixName != "") {
    taskName += "_";
    taskName += suffixName;
  }

  // Create task and configure as desired.
  AliAnalysisTaskJetDynamicalGrooming* task = new AliAnalysisTaskJetDynamicalGrooming(taskName.c_str());
  // Set a few general default.
  task->SetNCentBins(5);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);

  TString thename(njetsBase);
  // if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  // task->SetVzRange(-10.,10.);

  AliParticleContainer* trackCont; // = task->AddTrackContainer(ntracks);

  if ((jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub ||
     jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kEventSub) &&
    ((jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kData) ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia) ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef))) {
    trackCont = task->AddParticleContainer(ntracks);
  } else
    trackCont = task->AddTrackContainer(ntracks);

  // Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer* trackContUS = task->AddTrackContainer(ntracksUS);
  // Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer* trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  // Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia)
    trackContTrue->SetIsEmbedding(true);
  AliParticleContainer* trackContPartLevel = 0;

  if ((jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub) &&
    ((jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kMCTrue) ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef))) {
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  } else
    trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia)
    trackContPartLevel->SetIsEmbedding(true);
  // Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);

  AliClusterContainer* clusterCont = task->AddClusterContainer(nclusters);

  AliJetContainer* jetContBase = 0x0;
  AliJetContainer* jetContUS = 0x0;
  AliJetContainer* jetContTrue = 0x0;
  AliJetContainer* jetContPart = 0x0;
  TString strType(type);

  if ((jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kMCTrue ||
     (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kData) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub)
        jetContBase->SetAreaEmcCut(-2);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kDetEmbPartPythia) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);

      if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub)
        jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue, strType, R);
    if (jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    }

    if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub ||
      jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kEventSub) {
      jetContUS = task->AddJetContainer(njetsUS, strType, R);
      if (jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
      }
    }

    jetContPart = task->AddJetContainer(njetsPartLevel, strType, R);
    if (jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(acut);
    }
  }

  if (jetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef) {
    jetContBase = task->AddJetContainer(njetsBase, strType, R);
    if (jetContBase) {
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }

    jetContTrue = task->AddJetContainer(njetsTrue, strType, R);
    if (jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    }

    if (jetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub) {
      jetContUS = task->AddJetContainer(njetsUS, strType, R);
      if (jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
      }
    }

    jetContPart = task->AddJetContainer(njetsPartLevel, strType, R);
    if (jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(acut);
    }
  }

  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  // Configuration is via YAML, so nothing else to be done here.
  mgr->AddTask(task);

  // Create containers for input/output
  std::string listOutputName = task->DetermineOutputContainerName(taskName);
  std::string treeOutputName = task->DetermineOutputContainerName(taskName + "Tree");

  std::string outputFile = AliAnalysisManager::GetCommonFileName();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer* outputContainer = mgr->CreateContainer(
   listOutputName.c_str(), TList::Class(), AliAnalysisManager::kOutputContainer, outputFile.c_str());
  mgr->ConnectOutput(task, 1, outputContainer);

  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(treeOutputName.c_str(), TTree::Class(),
                               AliAnalysisManager::kOutputContainer, outputFile.c_str());
  mgr->ConnectOutput(task, 2, coutput2);

  return task;
}

std::string AliAnalysisTaskJetDynamicalGrooming::DetermineOutputContainerName(std::string containerName) const
{
  if (fJetShapeType == AliAnalysisTaskJetDynamicalGrooming::kMCTrue)
    containerName += "_MCTrue";
  if (fJetShapeType == AliAnalysisTaskJetDynamicalGrooming::kData)
    containerName += "_Data";
  if (fJetShapeType == AliAnalysisTaskJetDynamicalGrooming::kPythiaDef)
    containerName += "_PythiaDef";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kNoSub)
    containerName += "_NoSub";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kConstSub)
    containerName += "_ConstSub";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kEventSub)
    containerName += "_EventSub";
  if (fJetShapeSub == AliAnalysisTaskJetDynamicalGrooming::kDerivSub)
    containerName += "_DerivSub";
  if (fJetSelection == AliAnalysisTaskJetDynamicalGrooming::kInclusive)
    containerName += "_Incl";

  return containerName;
}

/**
 * Prints information about the jet-hadron performance task.
 *
 * @return std::string containing information about the task.
 */
std::string AliAnalysisTaskJetDynamicalGrooming::toString() const
{
  std::stringstream tempSS;
  tempSS << std::boolalpha;
  /*tempSS << "Particle collections:\n";
  TIter nextParticleCont(&fParticleCollArray);
  AliParticleContainer * particleCont;
  while ((particleCont = static_cast<AliParticleContainer *>(nextParticleCont()))) {
    tempSS << "\t" << particleCont->GetName() << ": " << particleCont->GetArrayName() << "\n";
  }
  tempSS << "Cluster collections:\n";
  TIter nextClusterCont(&fClusterCollArray);
  AliClusterContainer * clusterCont;
  while ((clusterCont = static_cast<AliClusterContainer *>(nextClusterCont()))) {
    tempSS << "\t" << clusterCont->GetName() << ": " << clusterCont->GetArrayName() << "\n";
  }
  tempSS << "Jet collections:\n";
  TIter nextJetCont(&fJetCollArray);
  AliJetContainer * jetCont;
  while ((jetCont = static_cast<AliJetContainer *>(nextJetCont()))) {
    tempSS << "\t" << jetCont->GetName() << ": " << jetCont->GetArrayName() << "\n";
  }
  // AliEventCuts
  tempSS << "AliEventCuts\n";
  tempSS << "\tEnabled: " << !fUseBuiltinEventSelection << "\n";
  // Efficiency
  tempSS << "Efficiency\n";
  tempSS << "\tSingle track efficiency identifier: " << fEfficiencyPeriodIdentifier << "\n";
  // QA
  tempSS << "QA Hists:\n";
  tempSS << "\tEnabled: " << fCreateQAHists << "\n";
  tempSS << "\tCreate cell QA hists: " << fCreateCellQAHists << "\n";
  // Use whether the pointer as null as a proxy. It's not ideal because it's not fully initialized
  // until after UserCreateOutputObjects(). But it's good enough for these purposes.
  tempSS << "\tCalculate event plane resolution (proxy of whether it's enabled - it may not be accurate): " <<
  (fFlowQnVectorManager != nullptr) << "\n";
  // Jet matching
  tempSS << "Jet matching:\n";
  tempSS << "\tEnabled: " << fPerformJetMatching << "\n";
  tempSS << "\tMax matching distance: " << fMaxJetMatchingDistance << "\n";
  // Response matrix
  tempSS << "Response matrix:\n";
  tempSS << "\tEnabled: " << fCreateResponseMatrix << "\n";
  tempSS << "\tConstruct response from 3 jet collections: " << fResponseFromThreeJetCollections << "\n";
  tempSS << "\tMin fraction shared pt: " << fMinFractionShared << "\n";
  tempSS << "\tJet leading hadron bias type: " << fLeadingHadronBiasType << "\n";
  tempSS << "\tResponse matrix fill map: \n";
  for (auto el : fResponseMatrixFillMap) {
    tempSS << "\t\tProperty " << el.first << " applied to jet " << el.second.first << "\n";
  }*/

  return tempSS.str();
}

/**
 * Print jet-hadron performance task information on an output stream using the string representation provided by
 * AliAnalysisTaskJetDynamicalGrooming::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream& AliAnalysisTaskJetDynamicalGrooming::Print(std::ostream& in) const
{
  in << toString();
  return in;
}

/**
 * Print basic jet-hadron performance task information using the string representation provided by
 * AliAnalysisTaskJetDynamicalGrooming::toString
 *
 * @param[in] opt Unused
 */
void AliAnalysisTaskJetDynamicalGrooming::Print(Option_t* opt) const { Printf("%s", toString().c_str()); }

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

/**
 * Implementation of the output stream operator for AliAnalysisTaskJetDynamicalGrooming. Printing
 * basic jet-hadron performance task information provided by function toString
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream& operator<<(std::ostream& in, const PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& myTask)
{
  std::ostream& result = myTask.Print(in);
  return result;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550.
 *
 * NOTE: We don't swap the base class values because the base class doesn't implement swap.
 */
void swap(PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& first,
     PWGJE::EMCALJetTasks::AliAnalysisTaskJetDynamicalGrooming& second)
{
  using std::swap;

  // Same ordering as in the constructors (for consistency)
  swap(first.fYAMLConfig, second.fYAMLConfig);
  swap(first.fConfigurationInitialized, second.fConfigurationInitialized);
}
