//
// Few substructure observables to try with the new PbPb data
//
//
//
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisDataSlot.h>
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
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>

#include "AliAODEvent.h"
#include "AliAnalysisTaskNewJetSubstructure.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNewJetSubstructure)

    //________________________________________________________________________
AliAnalysisTaskNewJetSubstructure::AliAnalysisTaskNewJetSubstructure()
: AliAnalysisTaskEmcalJet("AliAnalysisTaskNewJetSubstructure", kTRUE),
  fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
  fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.),
  fRMatching(0.2), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
  fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE),
  fDoFillMCLund(kFALSE), fCheckResolution(kFALSE), fSubjetCutoff(0.1),
  fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
  fDoAreaIterative(kTRUE), fPowerAlgo(1), fPhiCutValue(0.02),
  fEtaCutValue(0.02), fMagFieldPolarity(1), fDerivSubtrOrder(0),
  fPtJet(0x0), fHLundIterative(0x0), fHLundIterativeMC(0x0),
  fHLundIterativeMCDet(0x0), fHCheckResolutionSubjets(0x0),
  fStoreDetLevelJets(0), fTreeSubstructure(0), fDoSubJet(0)

{
  for (Int_t i = 0; i < 18; i++) {
    fShapesVar[i] = 0;
  }
  fShapesVar[18]  =  -1;
  fShapesVar[19]  =  -1;
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskNewJetSubstructure::AliAnalysisTaskNewJetSubstructure(
								     const char *name)
  : AliAnalysisTaskEmcalJet(name, kTRUE), fContainer(0),
    fMinFractionShared(0), fJetShapeType(kData), fJetShapeSub(kNoSub),
    fJetSelection(kInclusive), fPtThreshold(-9999.), fRMatching(0.2),
    fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
    fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE),
    fDoFillMCLund(kFALSE), fCheckResolution(kFALSE), fSubjetCutoff(0.1),
    fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
    fDoAreaIterative(kTRUE), fPowerAlgo(1), fPhiCutValue(0.02),
    fEtaCutValue(0.02), fMagFieldPolarity(1), fDerivSubtrOrder(0),
    fPtJet(0x0), fHLundIterative(0x0), fHLundIterativeMC(0x0),
    fHLundIterativeMCDet(0x0), fHCheckResolutionSubjets(0x0),
    fStoreDetLevelJets(0), fTreeSubstructure(0), fDoSubJet(0)
    
{
  // Standard constructor.
  for (Int_t i = 0; i < 18; i++) {
    fShapesVar[i] = -1;
  }
  fShapesVar[18]  =  -1;
  fShapesVar[19]  =  -1;
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskNewJetSubstructure::~AliAnalysisTaskNewJetSubstructure() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskNewJetSubstructure::UserCreateOutputObjects() {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fPtJet = new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec = 8;
  const Int_t nBinsSpec[8] = {50, 100, 100, 20, 100, 50, 100, 2};
  const Double_t lowBinSpec[8] = {0., -5, 0, 0, 0, 0, 0, 0};
  const Double_t hiBinSpec[8] = {5., 10., 200, 20, 200, 50, 50, 2};
  fHLundIterative =
      new THnSparseF("fHLundIterative",
                     "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
                     dimSpec, nBinsSpec, lowBinSpec, hiBinSpec);
  fOutput->Add(fHLundIterative);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec2 = 7;
  const Int_t nBinsSpec2[7] = {50, 100, 100, 20, 100, 50, 100};
  const Double_t lowBinSpec2[7] = {0., -5, 0, 0, 0, 0, 0};
  const Double_t hiBinSpec2[7] = {5., 10., 200, 20, 200, 100, 50};
  fHLundIterativeMC = new THnSparseF(
      "fHLundIterativeMC",
      "LundIterativePlotMC [log(1/theta),log(z*theta),pTjet,algo]", dimSpec2,
      nBinsSpec2, lowBinSpec2, hiBinSpec2);
  fOutput->Add(fHLundIterativeMC);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec3 = 7;
  const Int_t nBinsSpec3[7] = {50, 100, 100, 20, 100, 50, 100};
  const Double_t lowBinSpec3[7] = {0., -5, 0, 0, 0, 0, 0};
  const Double_t hiBinSpec3[7] = {5., 10., 200, 20, 200, 100, 50};
  fHLundIterativeMCDet = new THnSparseF(
      "fHLundIterativeMCDet",
      "LundIterativePlotMCDet [log(1/theta),log(z*theta),pTjet,algo]", dimSpec3,
      nBinsSpec3, lowBinSpec3, hiBinSpec3);
  fOutput->Add(fHLundIterativeMCDet);

  ////
  const Int_t dimResol = 5;
  const Int_t nBinsResol[5] = {10, 10, 80, 80, 80};
  const Double_t lowBinResol[5] = {0, 0, -1, -1, -1};
  const Double_t hiBinResol[5] = {200, 0.3, 1, 1, 1};
  fHCheckResolutionSubjets = new THnSparseF(
      "fHCheckResolutionSubjets", "Mom.Resolution of Subjets vs opening angle",
      dimResol, nBinsResol, lowBinResol, hiBinResol);
  fOutput->Add(fHCheckResolutionSubjets);

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i = 0; i < fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1 *>(fOutput->At(i));
    if (h1) {
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse *>(fOutput->At(i));
    if (hn)
      hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);
  const Int_t nVar = 20;
  const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeSubstructure = new TTree(nameoutput, nameoutput);
  TString *fShapesVarNames = new TString[nVar];

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
  //  if (fDoSubJet) fStoreDetLevelJets = true;
  if (fStoreDetLevelJets) {
    fShapesVarNames[12] = "ptJetDet";
    fShapesVarNames[13] = "ktgDet";
    fShapesVarNames[14] = "ngDet";
    fShapesVarNames[15] = "zgDet";
    fShapesVarNames[16] = "rgDet";
    fShapesVarNames[17] = "LeadingTrackPtDet";
  }
  if (fDoSubJet){
    if (fJetShapeType == kDetEmbPartPythia)
      {
	fShapesVarNames[18] = "subjet1";
	fShapesVarNames[19] = "subjet2";
      }
    if (fJetShapeType == kPythiaDef)
      {
	fShapesVarNames[12] = "subjet1";
        fShapesVarNames[13] = "subjet2";
      }
  }

  for (Int_t ivar = 0; ivar < nVar; ivar++) {
    fTreeSubstructure->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar],
                              Form("%s/F", fShapesVarNames[ivar].Data()));
  }

  PostData(1, fOutput);
  PostData(2, fTreeSubstructure);

  delete[] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskNewJetSubstructure::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskNewJetSubstructure::FillHistograms() {

  AliEmcalJet *jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(0);
  // container zero is always the base containe: the data container, the
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
      AliRhoParameter *rhoParam = dynamic_cast<AliRhoParameter *>(
								  InputEvent()->FindListObject("RhoSparseR020"));
      if (!rhoParam) {
        Printf("%s: Could not retrieve rho %s (some histograms will be filled "
               "with zero)!",
               GetName(), jetCont->GetRhoName().Data());
      } else
        rhoVal = rhoParam->GetVal();
      // rhom
      AliRhoParameter *rhomParam = dynamic_cast<AliRhoParameter *>(
								   InputEvent()->FindListObject("RhoMassSparseR020"));
      if (!rhomParam) {
        Printf("%s: Could not retrieve rho_m %s (some histograms will be "
               "filled with zero)!",
               GetName(), jetCont->GetRhoMassName().Data());
      } else
        rhoMassVal = rhomParam->GetVal();
    }

    while ((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1)
        continue;
      AliEmcalJet *jet2 = 0x0;
      AliEmcalJet *jet3 = 0x0;
      fPtJet->Fill(jet1->Pt());
      AliEmcalJet *jetUS = NULL;
      Int_t ifound = 0, jfound = 0;
      Int_t ilab = -1, jlab = -1;

      int sub1 = -1;
      int sub2 = -1;

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

        AliJetContainer *jetContUS = GetJetContainer(2);

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

	if(fJetShapeSub==kEventSub){
	  
	  jetUS = jet1->ClosestJet();
	  if (!jetUS) continue;
	  jet2 = jetUS->ClosestJet();
	}
	
        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub == kEventSub)) jet2 = jet1->ClosestJet();
        
	
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

	AliJetContainer *jetContTrue = GetJetContainer(1);
	AliJetContainer *jetContPart = GetJetContainer(3);
	
	if (fCheckResolution)
          CheckSubjetResolution(jet2, jetContTrue, jet3, jetContPart);

        Double_t fraction = 0;
        if (!(fJetShapeSub == kConstSub) && !(fJetShapeSub==kEventSub))
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

        AliJetContainer *jetContTrue = GetJetContainer(1);
        AliJetContainer *jetContUS = GetJetContainer(2);
        AliJetContainer *jetContPart = GetJetContainer(3);
	
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
        else if ((fJetShapeType == kPythiaDef) || (fJetShapeType == kMCTrue) ||
                 (fJetShapeType == kGenOnTheFly))
          ptSubtracted = jet1->Pt();
      }

      if (ptSubtracted < fPtThreshold)
        continue;
      
      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1))
        continue;

      fShapesVar[0] = ptSubtracted;
      fShapesVar[10] = jet1->MaxTrackPt();
      
      if(fCutDoubleCounts==kTRUE && fJetShapeType==kDetEmbPartPythia) if(jet1->MaxTrackPt()>jet2->MaxTrackPt()) continue;

      fastjet::PseudoJet *sub1Hyb=new fastjet::PseudoJet();
      fastjet::PseudoJet *sub2Hyb=new fastjet::PseudoJet();
      std::vector<fastjet::PseudoJet>* const1Hyb = new std::vector<fastjet::PseudoJet>();
      std::vector<fastjet::PseudoJet>* const2Hyb = new std::vector<fastjet::PseudoJet>();
      IterativeParents(jet1, jetCont, sub1Hyb, sub2Hyb, const1Hyb, const2Hyb);
    
      Float_t ptMatch = 0.;
      Float_t leadTrackMatch = 0.;
      Double_t ktgMatch = 0;
      ;
      Double_t nsdMatch = 0;
      Double_t zgMatch = 0;
      Double_t rgMatch = 0;
      Float_t ptDet = 0.;
      Float_t leadTrackDet = 0.;
      Double_t ktgDet = 0;
      ;
      Double_t nsdDet = 0;
      Double_t zgDet = 0;
      Double_t rgDet = 0;
      Double_t aver1 = 0;
      Double_t aver2 = 0;
      Double_t aver3 = 0;
      Double_t aver4 = 0;
      Int_t kMatched = 0;
      fastjet::PseudoJet *sub1Part = new fastjet::PseudoJet();
      fastjet::PseudoJet *sub2Part =new fastjet::PseudoJet();
      std::vector<fastjet::PseudoJet>* const1Part = new std::vector<fastjet::PseudoJet>();
      std::vector<fastjet::PseudoJet>* const2Part = new std::vector<fastjet::PseudoJet>();
      fastjet::PseudoJet *sub1Det = new fastjet::PseudoJet();
      fastjet::PseudoJet *sub2Det = new fastjet::PseudoJet();
      std::vector<fastjet::PseudoJet>* const1Det = new std::vector<fastjet::PseudoJet>();
      std::vector<fastjet::PseudoJet>* const2Det = new std::vector<fastjet::PseudoJet>();

      if (fJetShapeType == kPythiaDef) {
        kMatched = 1;
        if (fJetShapeSub == kConstSub)
          kMatched = 3;

        ptMatch = jet3->Pt();
	leadTrackMatch = jet3->MaxTrackPt();
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4, sub1Det, sub2Det, const1Det, const2Det);
        ktgMatch = aver1;
        nsdMatch = aver2;
        zgMatch = aver3;
        rgMatch = aver4;
      }

      if (fJetShapeType == kDetEmbPartPythia) {
        if (fJetShapeSub == kConstSub)
          kMatched = 3;
        if (fJetShapeSub == kDerivSub)
          kMatched = 2;
        ptMatch = jet3->Pt();
        leadTrackMatch = jet3->MaxTrackPt();
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4, sub1Part, sub2Part, const1Part, const2Part);
        ktgMatch = aver1;
        nsdMatch = aver2;
        zgMatch = aver3;
        rgMatch = aver4;
        if (fStoreDetLevelJets) {
          ptDet = jet2->Pt();
          leadTrackDet = jet2->MaxTrackPt();
          IterativeParentsMCAverage(jet2, 1, ktgDet, nsdDet, zgDet, rgDet, sub1Det, sub2Det, const1Det, const2Det);
	}
      }

      if (fJetShapeType == kMCTrue || fJetShapeType == kData ||
          fJetShapeType == kGenOnTheFly) {

        ptMatch = 0.;
        leadTrackMatch = 0.;
        ktgMatch = 0.;
        nsdMatch = 0.;
        zgMatch = 0;
        rgMatch = 0;
      }

      fShapesVar[5] = ptMatch;
      fShapesVar[6] = ktgMatch;
      fShapesVar[7] = nsdMatch;
      fShapesVar[8] = zgMatch;
      fShapesVar[9] = rgMatch;
      fShapesVar[11] = leadTrackMatch;
      if (fStoreDetLevelJets) {
        fShapesVar[12] = ptDet;
        fShapesVar[13] = ktgDet;
        fShapesVar[14] = nsdDet;
        fShapesVar[15] = zgDet;
        fShapesVar[16] = rgDet;
        fShapesVar[17] = leadTrackDet;
      }

      if (fDoSubJet)
	{
	  int nSDsub = -1;
	  if (fJetShapeType == kDetEmbPartPythia) nSDsub = nsdDet;	  	  
	  else if(fJetShapeType == kPythiaDef) nSDsub = nsdMatch;
	  if ((nSDsub != 0) && (fShapesVar[2] != 0))
	    {
	      if (CompareSubjets(sub1Det, sub1Hyb, const1Det, const1Hyb)) sub1 = 1;
	      else if (CompareSubjets(sub1Det, sub2Hyb, const1Det, const2Hyb)) sub1 = 2;
	      else sub1 = 3;
	      if (CompareSubjets(sub2Det, sub2Hyb, const2Det, const2Hyb)) sub2 = 2;
	      else if (CompareSubjets(sub2Det, sub1Hyb, const2Det, const1Hyb)) sub2 = 1;
	      else sub2 = 3;
	    }
	  else if (fShapesVar[2]  != 0) {sub1 = 0; sub2 = 0;}	      
	  if (fJetShapeType == kDetEmbPartPythia) 
	    {
	      fShapesVar[18] = sub1;
	      fShapesVar[19] = sub2;
	    }
	  else if(fJetShapeType == kPythiaDef)
	    {
	      fShapesVar[12] = sub1;
              fShapesVar[13] = sub2;
	    }   
	}


      fTreeSubstructure->Fill();
      delete sub1Det;
      delete sub2Det;
      delete sub1Hyb;
      delete sub2Hyb;
      delete const1Det;
      delete const2Det;
      delete const1Hyb;
      delete const2Hyb;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskNewJetSubstructure::GetJetMass(AliEmcalJet *jet,
                                                      Int_t jetContNb = 0) {
  // calc subtracted jet mass
  if ((fJetShapeSub == kDerivSub) && (jetContNb == 0))
    if (fDerivSubtrOrder == 1)
      return jet->GetShapeProperties()->GetFirstOrderSubtracted();
    else
      return jet->GetShapeProperties()->GetSecondOrderSubtracted();
  else
    return jet->M();
}

//________________________________________________________________________
Float_t AliAnalysisTaskNewJetSubstructure::Angularity(AliEmcalJet *jet,
                                                      Int_t jetContNb = 0) {

  AliJetContainer *jetCont = GetJetContainer(jetContNb);
  if (!jet->GetNumberOfTracks())
    return 0;
  Double_t den = 0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  for (UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle *>(
        jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));

    if (!vp1) {
      Printf("AliVParticle associated to constituent not found");
      continue;
    }

    Double_t dphi = RelativePhi(vp1->Phi(), jet->Phi());
    Double_t dr2 =
        (vp1->Eta() - jet->Eta()) * (vp1->Eta() - jet->Eta()) + dphi * dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num = num + vp1->Pt() * dr;
    den = den + vp1->Pt();
  }
  return num / den;
}

//________________________________________________________________________
Float_t
AliAnalysisTaskNewJetSubstructure::GetJetAngularity(AliEmcalJet *jet,
                                                    Int_t jetContNb = 0) {

  if ((fJetShapeSub == kDerivSub) && (jetContNb == 0))
    if (fDerivSubtrOrder == 1)
      return jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
    else
      return jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  else
    return Angularity(jet, jetContNb);
}

//__________________________________________________________________________________
Double_t AliAnalysisTaskNewJetSubstructure::RelativePhi(Double_t mphi,
                                                        Double_t vphi) {

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

//_________________________________________________________________________                                                                               
Bool_t AliAnalysisTaskNewJetSubstructure::CompareSubjets(fastjet::PseudoJet *subDet,  fastjet::PseudoJet *subHyb, std::vector<fastjet::PseudoJet> *constDet, std::vector<fastjet::PseudoJet>* constHyb)
{
  double pT_det = subDet->pt();
  double sumpT = 0;
  double delta =  0.01;

  for (int i = 0; i < constDet->size(); i++)
      {
	double eta_det = constDet->at(i).eta();
	double phi_det = constDet->at(i).phi();
	for (int j  = 0; j < constHyb->size(); j++)
	  {
	    double eta_hyb = constHyb->at(j).eta();
	    double phi_hyb = constHyb->at(j).phi();
	    double deta = eta_hyb - eta_det;
	    deta = std::sqrt(deta*deta);
	    if (deta > delta) continue;
	    double dphi = phi_hyb - phi_det;
	    dphi = std::sqrt(dphi*dphi);
	    if (dphi > delta) continue;
	    sumpT+=constDet->at(i).pt();
	  }
      }
  if (sumpT/pT_det > 0.5) return true;
  else return false;
}

//_________________________________________________________________________
void AliAnalysisTaskNewJetSubstructure::IterativeParentsAreaBased(
    AliEmcalJet *fJet, AliJetContainer *fJetCont) {
  // to still change and implement the 4 vector bkg subtraction to the subjets
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk)
        continue;
      //      if (fDoTwoTrack == kTRUE && CheckClosePartner(i, fJet, fTrk, fTrackCont))
      //        continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
      fInputVectors.push_back(PseudoTracks);
    }
  fastjet::JetAlgorithm jetalgo(fastjet::genkt_algorithm);
  fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);

  fastjet::JetDefinition fJetDef(jetalgo, 1., fPowerAlgo,
                                 static_cast<fastjet::RecombinationScheme>(0),
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
        Double_t LundEntries[7] = {
            y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, cumtf};
        fHLundIterative->Fill(LundEntries);
      }

      jj = jet_sub1;
    }

    fShapesVar[1] = xktg;
    fShapesVar[2] = nsd;
    fShapesVar[3] = zg;
    fShapesVar[4] = Rg;

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}
//_________________________________________________________________________
void AliAnalysisTaskNewJetSubstructure::IterativeParents(
							 AliEmcalJet *fJet, AliJetContainer *fJetCont,fastjet::PseudoJet *sub1,  fastjet::PseudoJet *sub2, std::vector < fastjet::PseudoJet > *const1, std::vector < fastjet::PseudoJet > *const2) {

  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;
  unsigned int constituentIndex = 0;
  for (auto part: fJet->GetParticleConstituents()) {
    if (fDoTwoTrack == kTRUE && CheckClosePartner(fJet, part))
      continue;
    PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E());
    PseudoTracks.set_user_index(constituentIndex);
    fInputVectors.push_back(PseudoTracks);
    constituentIndex++;
  }
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
  fastjet::JetDefinition fJetDef(jetalgo, 1.,
                                 static_cast<fastjet::RecombinationScheme>(0),
                                 fastjet::BestFJ30);

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
    fastjet::PseudoJet j1first;
    fastjet::PseudoJet j2first;
    jj = fOutputJets[0];

    double nall = 0;
    double nsd = 0;
    int flagSubjet = 0;
    int flagSubjetkT = 0;
    double flagConst=0;
    double Rg = 0;
    double zg = 0;
    double xktg = 0;
    double cumtf = 0;
    while (jj.has_parents(j1, j2)) {
      nall = nall + 1;
      if (j1.perp() < j2.perp())
        swap(j1, j2);
      flagConst=0;
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      double lnpt_rel = log(xkt);
      double y = log(1. / delta_R);
      double form = 2 * 0.197 * j2.e() / (xkt * xkt);
      double rad = j1.e()+j2.e();
      double z = j2.perp() / (j2.perp() + j1.perp());
      vector < fastjet::PseudoJet > constitj1 = sorted_by_pt(j1.constituents());
      if(constitj1[0].perp()>fMinPtConst) flagConst=1; 
      
      if (z > fHardCutoff)
        nsd = nsd + 1;
      if (z > fHardCutoff && flagSubjet == 0) {
        zg = z;
        xktg = xkt;
        Rg = delta_R;
	j1first =j1;
	*sub1 = j1first;
	j2first =j2;
	*sub2 = j2first;
        flagSubjet = 1;
      }
      if (lnpt_rel > 0) {
	cumtf = cumtf + form;
	if ((nsd == 0) && (flagSubjetkT == 0)) {
	  xktg = xkt;
	  flagSubjetkT = 1;
	}
      }
      
      Double_t LundEntries[8] = {
	y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, cumtf,flagConst};
      fHLundIterative->Fill(LundEntries);
      
      jj = j1;
    }
    if (sub1->has_constituents()) *const1 = sub1->constituents();
    if (sub2->has_constituents()) *const2 = sub2->constituents();

    fShapesVar[1] = xktg;
    fShapesVar[2] = nsd;
    fShapesVar[3] = zg;
    fShapesVar[4] = Rg;

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}
//_________________________________________________________________________
void AliAnalysisTaskNewJetSubstructure::IterativeParentsMCAverage(
    AliEmcalJet *fJet, Int_t km, Double_t &average1, Double_t &average2,
    Double_t &average3, Double_t &average4, fastjet::PseudoJet *sub1,  fastjet::PseudoJet *sub2, std::vector < fastjet::PseudoJet > *const1, std::vector < fastjet::PseudoJet > *const2) {
  AliJetContainer *jetCont = GetJetContainer(km);
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;
  unsigned int constituentIndex = 0;
  for (auto part: fJet->GetParticleConstituents()) {
    PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E());
    PseudoTracks.set_user_index(constituentIndex);
    fInputVectors.push_back(PseudoTracks);
  }
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);

  fastjet::JetDefinition fJetDef(jetalgo, 1.,
                                 static_cast<fastjet::RecombinationScheme>(0),
                                 fastjet::BestFJ30);


  try {
    fastjet::ClusterSequence fClustSeqSA(fInputVectors, fJetDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets = fClustSeqSA.inclusive_jets(0);


    fastjet::PseudoJet jj;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet j1first;
    fastjet::PseudoJet j2first;
    jj = fOutputJets[0];
    int flagSubjet = 0;
    int flagSubjetkT = 0;
    double nall = 0;
    double nsd = 0;

    double zg = 0;
    double xktg = 0;
    double Rg = 0;

    double cumtf = 0;
    while (jj.has_parents(j1, j2)) {
      nall = nall + 1;
      if (j1.perp() < j2.perp())
        swap(j1, j2);
      std::vector<fastjet::PseudoJet> v1 = j1.constituents();
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      double lnpt_rel = log(xkt);
      double y = log(1. / delta_R);
      double form = 2 * 0.197 * j2.e() / (xkt * xkt);
      double rad = j1.e()+j2.e();
      double z = j2.perp() / (j2.perp() + j1.perp());
      if (z > fHardCutoff)
        nsd = nsd + 1;
      if (z > fHardCutoff && flagSubjet == 0) {
	zg = z;
        xktg = xkt;
        Rg = delta_R;
	j1first = j1;
	*sub1 = j1first;
	j2first = j2;
	*sub2 = j2first;
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
        Double_t LundEntries[7] = {
            y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, cumtf};
        fHLundIterativeMC->Fill(LundEntries);
        if (fStoreDetLevelJets) {
          fHLundIterativeMCDet->Fill(LundEntries);
        }
      }

      jj = j1;
    }

    average1 = xktg;
    average2 = nsd;
    average3 = zg;
    average4 = Rg;
    if (sub1->has_constituents()) *const1 = sub1->constituents();
    if (sub2->has_constituents()) *const2 = sub2->constituents();

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

//_________________________________________________________________________
void AliAnalysisTaskNewJetSubstructure::CheckSubjetResolution(
    AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJetM,
    AliJetContainer *fJetContM) {

  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  std::vector<fastjet::PseudoJet> fInputVectorsM;
  fInputVectorsM.clear();
  fastjet::PseudoJet PseudoTracksM;

  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();
  AliParticleContainer *fTrackContM = fJetContM->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk)
        continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
      fInputVectors.push_back(PseudoTracks);
    }
  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
  fastjet::JetDefinition fJetDef(jetalgo, 1.,
                                 static_cast<fastjet::RecombinationScheme>(0),
                                 fastjet::BestFJ30);

  if (fTrackContM)
    for (Int_t i = 0; i < fJetM->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJetM->TrackAt(i, fTrackContM->GetArray());
      if (!fTrk)
        continue;
      PseudoTracksM.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracksM.set_user_index(fJetM->TrackAt(i) + 100);
      fInputVectorsM.push_back(PseudoTracksM);
    }
  fastjet::JetAlgorithm jetalgoM(fastjet::cambridge_algorithm);
  fastjet::JetDefinition fJetDefM(jetalgoM, 1.,
                                  static_cast<fastjet::RecombinationScheme>(0),
                                  fastjet::BestFJ30);

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
      Double_t ResolEntries[5] = {fOutputJets[0].perp(), delta_R, resid1,
                                  resid2, residz};
      fHCheckResolutionSubjets->Fill(ResolEntries);
    }

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

bool AliAnalysisTaskNewJetSubstructure::CheckClosePartner(AliEmcalJet* jet, PWG::JETFW::AliEmcalParticleJetConstituent & part1)
{
  auto jetConstituents = jet->GetParticleConstituents();
  for (auto part2 : jetConstituents) {
    if (part2 == part1) {
      continue;
    }
    Double_t phi1 = part1.Phi();
    Double_t phi2 = part2.Phi();
    Double_t chg1 = part1.GetParticle()->Charge();
    Double_t chg2 = part2.GetParticle()->Charge();
    Double_t ptv1 = part1.Pt();
    Double_t ptv2 = part2.Pt();
    Double_t deta = part2.Eta() - part1.Eta();
    const Float_t kLimit = fPhiCutValue * 3;

    if (TMath::Abs(part1.Eta() - part2.Eta()) < fEtaCutValue * 2.5 * 3) {
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

/*Bool_t AliAnalysisTaskNewJetSubstructure::CheckClosePartner(
    Int_t index, AliEmcalJet *fJet, AliVParticle *fTrk1,
    AliParticleContainer *fTrackCont) {
  // check if tracks are close//
  for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
    AliVParticle *fTrk2 = fJet->TrackAt(i, fTrackCont->GetArray());
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
      Float_t initdpsinner =
          (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * 0.8 / ptv2) -
           (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * 0.8 / ptv1)));

      Float_t initdpsouter =
          (phi2 - TMath::ASin(0.075 * chg2 * fMagFieldPolarity * 2.5 / ptv2) -
           (phi1 - TMath::ASin(0.075 * chg1 * fMagFieldPolarity * 2.5 / ptv1)));

      initdpsinner = TVector2::Phi_mpi_pi(initdpsinner);
      initdpsouter = TVector2::Phi_mpi_pi(initdpsouter);

      if (TMath::Abs(initdpsinner) < kLimit ||
          TMath::Abs(initdpsouter) < kLimit ||
          initdpsinner * initdpsouter < 0) {
        Double_t mindps = 1e5;

        for (Double_t rad = 0.8; rad < 2.51; rad += 0.01) {
          Double_t dps =
              (phi2 -
               TMath::ASin(0.075 * chg2 * fMagFieldPolarity * rad / ptv2) -
               (phi1 -
                TMath::ASin(0.075 * chg1 * fMagFieldPolarity * rad / ptv1)));
          dps = TVector2::Phi_mpi_pi(dps);
          if (TMath::Abs(dps) < TMath::Abs(mindps))
            mindps = dps;
        }
        if (TMath::Abs(mindps) < fPhiCutValue &&
            TMath::Abs(deta) < fEtaCutValue)
          return kTRUE;
      }
    }
  }
  return kFALSE;
  }*/

//________________________________________________________________________
Bool_t AliAnalysisTaskNewJetSubstructure::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskNewJetSubstructure::Terminate(Option_t *) {
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}
