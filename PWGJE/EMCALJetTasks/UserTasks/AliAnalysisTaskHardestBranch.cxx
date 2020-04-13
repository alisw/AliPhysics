//
// Checking the hardest branch
//  Leticia Cunqueiro
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
#include <vector>
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
#include <iostream>     // std::cout
#include <algorithm> 

#include "AliAODEvent.h"
#include "AliAnalysisTaskHardestBranch.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHardestBranch)

    //________________________________________________________________________
    AliAnalysisTaskHardestBranch::AliAnalysisTaskHardestBranch()
    : AliAnalysisTaskEmcalJet("AliAnalysisTaskHardestBranch", kTRUE),
      fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
      fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.),
      fRMatching(0.2), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
      fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE),
      fDoFillMCLund(kFALSE), fCheckResolution(kFALSE), fSubjetCutoff(0.1),
      fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE),
      fDoAreaIterative(kTRUE), fPowerAlgo(1), fPhiCutValue(0.02),
      fEtaCutValue(0.02), fMagFieldPolarity(1), fDerivSubtrOrder(0),
      fPtJet(0x0),fTreeSubstructure(0)

{
  for (Int_t i = 0; i < 10; i++) {
    fShapesVar[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHardestBranch::AliAnalysisTaskHardestBranch(
    const char *name)
    : AliAnalysisTaskEmcalJet(name, kTRUE), fContainer(0),
      fMinFractionShared(0), fJetShapeType(kData), fJetShapeSub(kNoSub),
      fJetSelection(kInclusive), fPtThreshold(-9999.), fRMatching(0.2),
      fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
      fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE),
      fDoFillMCLund(kFALSE), fCheckResolution(kFALSE), fSubjetCutoff(0.1),
      fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE),
      fDoAreaIterative(kTRUE), fPowerAlgo(1), fPhiCutValue(0.02),
      fEtaCutValue(0.02), fMagFieldPolarity(1), fDerivSubtrOrder(0),
      fPtJet(0x0),fTreeSubstructure(0)

{
  // Standard constructor.
  for (Int_t i = 0; i < 10; i++) {
    fShapesVar[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHardestBranch::~AliAnalysisTaskHardestBranch() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskHardestBranch::UserCreateOutputObjects() {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fPtJet = new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);

  

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
  const Int_t nVar = 13;
  const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeSubstructure = new TTree(nameoutput, nameoutput);
  TString *fShapesVarNames = new TString[nVar];

  fShapesVarNames[0] = "ptJet";
  fShapesVarNames[1] = "ktg";
  fShapesVarNames[2] = "tfg";
  fShapesVarNames[3] = "zg";
  fShapesVarNames[4] = "rg";
  fShapesVarNames[5] = "ng";
  fShapesVarNames[6] = "ptJetMatch";
  fShapesVarNames[7] = "ktgMatch";
  fShapesVarNames[8] = "tfgMatch";
  fShapesVarNames[9] = "zgMatch";
  fShapesVarNames[10] = "rgMatch";
  fShapesVarNames[11] = "ngMatch";
  fShapesVarNames[12] = "LeadingTrackPtMatch";
 

  for (Int_t ivar = 0; ivar < nVar; ivar++) {
    cout << "looping over variables" << endl;
    fTreeSubstructure->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar],
                              Form("%s/F", fShapesVarNames[ivar].Data()));
  }

  PostData(1, fOutput);
  PostData(2, fTreeSubstructure);

  delete[] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHardestBranch::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHardestBranch::FillHistograms() {

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

        if (!(fJetShapeSub == kConstSub))
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

        Double_t fraction = 0;
        if (!(fJetShapeSub == kConstSub))
          fraction = jetCont->GetFractionSharedPt(jet1);
        if (fJetShapeSub == kConstSub)
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
      if (fJetShapeSub == kConstSub)
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
      IterativeParents(jet1, jetCont);

      Float_t ptMatch = 0.;
      Float_t leadTrackMatch = 0.;
      Double_t ktgMatch = 0;
      Double_t tfgMatch =0;
      Double_t ngMatch = 0;
      Double_t zgMatch = 0;
      Double_t rgMatch = 0;
      Float_t ptDet = 0.;
     
   
      Double_t aver1 = 0;
      Double_t aver2 = 0;
      Double_t aver3 = 0;
      Double_t aver4 = 0;
      Double_t aver5 = 0;
      Int_t kMatched = 0;
      if (fJetShapeType == kPythiaDef) {
        kMatched = 1;
        if (fJetShapeSub == kConstSub)
          kMatched = 3;

        ptMatch = jet3->Pt();
        leadTrackMatch = jet3->MaxTrackPt();
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4, aver5);
        ktgMatch = aver1;
        tfgMatch = aver2;
        zgMatch = aver3;
        rgMatch = aver4;
	ngMatch = aver5;
      }

      if (fJetShapeType == kDetEmbPartPythia) {
        if (fJetShapeSub == kConstSub)
          kMatched = 3;
        if (fJetShapeSub == kDerivSub)
          kMatched = 2;
        ptMatch = jet3->Pt();
        leadTrackMatch = jet3->MaxTrackPt();
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4, aver5);
        ktgMatch = aver1;
        tfgMatch = aver2;
        zgMatch = aver3;
        rgMatch = aver4;
        ngMatch = aver5;
      }

      if (fJetShapeType == kMCTrue || fJetShapeType == kData ||
          fJetShapeType == kGenOnTheFly) {

        ptMatch = 0.;
        leadTrackMatch = 0.;
        ktgMatch = 0.;
        tfgMatch = 0.;
        zgMatch = 0;
        rgMatch = 0;
        ngMatch = 0;
      }

      fShapesVar[6] = ptMatch;
      fShapesVar[7] = ktgMatch;
      fShapesVar[8] = tfgMatch;
      fShapesVar[9] = zgMatch;
      fShapesVar[10] = rgMatch;
      fShapesVar[11] = ngMatch;
      fShapesVar[12] = leadTrackMatch;
      

      fTreeSubstructure->Fill();
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskHardestBranch::GetJetMass(AliEmcalJet *jet,
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
Float_t AliAnalysisTaskHardestBranch::Angularity(AliEmcalJet *jet,
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
AliAnalysisTaskHardestBranch::GetJetAngularity(AliEmcalJet *jet,
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
Double_t AliAnalysisTaskHardestBranch::RelativePhi(Double_t mphi,
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
void AliAnalysisTaskHardestBranch::IterativeParents(AliEmcalJet *fJet, AliJetContainer *fJetCont) {

  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
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
    std::vector<Double_t>  zvec;
    std::vector<Double_t>  ktvec;
    std::vector<Double_t> thetavec;
    std::vector<Double_t> tformvec;
     std::vector<Double_t> nvec;
    int indmax;
    while (jj.has_parents(j1, j2)) {
      nall = nall + 1;

      if (j1.perp() < j2.perp())
        swap(j1, j2);

      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
   
      double z = j2.perp() / (j2.perp() + j1.perp());
      double form = 2 * 0.197 * j2.e() / ((1.-z)*xkt * xkt);
      
      if(z>0.1) nsd=nsd+1;
      nvec.push_back(nsd);
      zvec.push_back(z);
      ktvec.push_back(xkt);
      thetavec.push_back(delta_R);
      tformvec.push_back(form);
            jj = j1;
    }
    if(nall>0){

     auto result = std::max_element(ktvec.begin(), ktvec.end());
     indmax=std::distance(ktvec.begin(), result);
   
      
    
    fShapesVar[1] = ktvec[indmax];
    fShapesVar[2] = tformvec[indmax];
    fShapesVar[3] = zvec[indmax];
    fShapesVar[4] = thetavec[indmax];
    fShapesVar[5] = nvec[indmax];}

     if(nall==0){
       
      fShapesVar[1] = -1;
      fShapesVar[2] = -1;;
      fShapesVar[3] = -1;;
      fShapesVar[4] = -1;;
      fShapesVar[5] = -1;}

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;

}
//_________________________________________________________________________
void AliAnalysisTaskHardestBranch::IterativeParentsMCAverage(AliEmcalJet *fJet, Int_t km, Double_t &average1, Double_t &average2, Double_t &average3, Double_t &average4, Double_t &average5) {
  AliJetContainer *jetCont = GetJetContainer(km);
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

  AliParticleContainer *fTrackCont = jetCont->GetParticleContainer();

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
    double nall = 0;
    double nsd = 0;
     
    std::vector<Double_t>  zvec;
    std::vector<Double_t> ktvec;
    std::vector<Double_t> thetavec;
    std::vector<Double_t>  tformvec;
    std::vector<Double_t> nvec;
    int indmax;

    
    while (jj.has_parents(j1, j2)) {
      nall = nall + 1;

      if (j1.perp() < j2.perp())
        swap(j1, j2);
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
   
       double z = j2.perp() / (j2.perp() + j1.perp());
       double form = 2 * 0.197 * j2.e() / ((1.-z)*xkt * xkt);
      double rad = j2.e();
     
      if(z>0.1) nsd=nsd+1;
       
      zvec.push_back(z);
      ktvec.push_back(xkt);
      thetavec.push_back(delta_R);
      tformvec.push_back(form);
      nvec.push_back(nsd);
      jj = j1;
    }

    if(nall>0){    
     auto result = std::max_element(ktvec.begin(), ktvec.end());
     indmax=std::distance(ktvec.begin(), result);


    average1 = ktvec[indmax];
    average2 = tformvec[indmax];
    average3 = zvec[indmax];
    average4 = thetavec[indmax];
    average5 = nvec[indmax];}

    if(nall==0){
       average1 = -1;
    average2 = -1;
    average3 = -1;
    average4 = -1;
    average5 = -1;}

    

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskHardestBranch::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskHardestBranch::Terminate(Option_t *) {
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}
