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
#include "AliAnalysisTaskLundPlane.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskLundPlane)

    //________________________________________________________________________
AliAnalysisTaskLundPlane::AliAnalysisTaskLundPlane()
: AliAnalysisTaskEmcalJet("AliAnalysisTaskLundPlane", kTRUE),
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
  fStoreDetLevelJets(0), fTreeSubstructure(0), fDoSubJet(0),fShapesVar_Splittings_angle(0),
  fShapesVar_Splittings_kt(0),fShapesVar_Splittings_z(0),fShapesVar_Splittings_energy(0),fShapesVar_Splittings_eta1(0),
  fShapesVar_Splittings_eta2(0),fShapesVar_Splittings_phi1(0),fShapesVar_Splittings_phi2(0)
{
  
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskLundPlane::AliAnalysisTaskLundPlane(
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
    fStoreDetLevelJets(0), fTreeSubstructure(0), fDoSubJet(0),fShapesVar_Splittings_angle(0),
    fShapesVar_Splittings_kt(0),fShapesVar_Splittings_z(0),fShapesVar_Splittings_energy(0),fShapesVar_Splittings_eta1(0),
    fShapesVar_Splittings_eta2(0),fShapesVar_Splittings_phi1(0),fShapesVar_Splittings_phi2(0)
    
{
 
  SetMakeGeneralHistograms(kTRUE);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskLundPlane::~AliAnalysisTaskLundPlane() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskLundPlane::UserCreateOutputObjects() {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);

 
  const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  
  fTreeSplittings = new TTree(nameoutput_Splittings, nameoutput_Splittings);
  TString *fShapesVarNames_Splittings=new TString[18];

  fShapesVarNames_Splittings[0] = "angle";
  fShapesVarNames_Splittings[1] = "kt";
  fShapesVarNames_Splittings[2] = "z";
  fShapesVarNames_Splittings[3] = "energy";
  fShapesVarNames_Splittings[4] = "eta1";
  fShapesVarNames_Splittings[5] = "phi1";
  fShapesVarNames_Splittings[6] = "eta2";
  fShapesVarNames_Splittings[7] = "phi2";
  fShapesVarNames_Splittings[8] = "angle_part";
  fShapesVarNames_Splittings[9] = "kt_part";
  fShapesVarNames_Splittings[10] = "z_part";
  fShapesVarNames_Splittings[11] = "energy_part";
  fShapesVarNames_Splittings[12] = "eta1_part";
  fShapesVarNames_Splittings[13] = "phi1_part";
  fShapesVarNames_Splittings[14] = "eta2_part";
  fShapesVarNames_Splittings[15] = "phi2_part";
  fShapesVarNames_Splittings[16] = "ptjet";
  fShapesVarNames_Splittings[17] = "ptjet_part"; 
  
  fTreeSplittings->Branch(fShapesVarNames_Splittings[0].Data(), &fShapesVar_Splittings_angle, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[1].Data(), &fShapesVar_Splittings_kt, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[2].Data(), &fShapesVar_Splittings_z, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[3].Data(), &fShapesVar_Splittings_energy, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[4].Data(), &fShapesVar_Splittings_eta1, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[5].Data(), &fShapesVar_Splittings_phi1, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[6].Data(), &fShapesVar_Splittings_eta2, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[7].Data(), &fShapesVar_Splittings_phi2, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[8].Data(), &fShapesVar_Splittings_angle_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[9].Data(), &fShapesVar_Splittings_kt_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[10].Data(), &fShapesVar_Splittings_z_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[11].Data(), &fShapesVar_Splittings_energy_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[12].Data(), &fShapesVar_Splittings_eta1_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[13].Data(), &fShapesVar_Splittings_phi1_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[14].Data(), &fShapesVar_Splittings_eta2_part, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[15].Data(), &fShapesVar_Splittings_phi2_part, 0,1);
   fTreeSplittings->Branch(fShapesVarNames_Splittings[16].Data(), &fShapesVar_Splittings_ptjet, 0,1);
  fTreeSplittings->Branch(fShapesVarNames_Splittings[17].Data(), &fShapesVar_Splittings_ptjet_part, 0,1);

  PostData(1, fOutput);
  PostData(2, fTreeSplittings);

  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskLundPlane::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskLundPlane::FillHistograms() {

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
      }
       

      Double_t ptSubtracted = 0;
      if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub)
        ptSubtracted = jet1->Pt();

      else if (fJetShapeSub == kDerivSub) {
        ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
      }

      else if (fJetShapeSub == kNoSub) ptSubtracted = jet1->Pt();
      

      if (ptSubtracted < fPtThreshold)
        continue;
      
      if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1))
        continue;

      fShapesVar[0] = ptSubtracted;
     
      IterativeDeclustering(jet1, jetCont);
    
    
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
    
      if (fJetShapeType == kPythiaDef) {
        kMatched = 1;
        if (fJetShapeSub == kConstSub)
          kMatched = 3;

        ptMatch = jet3->Pt();
	leadTrackMatch = jet3->MaxTrackPt();
        IterativeDeclustering(jet3, kMatched);
       
      }

   
      fTreeSplittings->Fill();

      fShapesVar_Splittings_angle.clear();
      fShapesVar_Splittings_kt.clear(); 
      fShapesVar_Splittings_z.clear();
      fShapesVar_Splittings_energy.clear();
      fShapesVar_Splittings_eta1.clear();
      fShapesVar_Splittings_phi1.clear();
      fShapesVar_Splittings_eta2.clear();
      fShapesVar_Splittings_phi2.clear();
     
      fShapesVar_Splittings_angle_part.clear();
      fShapesVar_Splittings_kt_part.clear(); 
      fShapesVar_Splittings_z_part.clear();
      fShapesVar_Splittings_energy_part.clear();
      fShapesVar_Splittings_eta1_part.clear();
      fShapesVar_Splittings_phi1_part.clear();
      fShapesVar_Splittings_eta2_part.clear();
      fShapesVar_Splittings_phi2_part.clear();
   


    
    }
  }

  return kTRUE;
}


//__________________________________________________________________________________
Double_t AliAnalysisTaskLundPlane::RelativePhi(Double_t mphi,
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

                                                            
//_____________________________
void AliAnalysisTaskLundPlane::IterativeDeclustering(AliEmcalJet *fJet, AliJetContainer *fJetCont) {

  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;
  unsigned int constituentIndex = 0;
  for (auto part: fJet->GetParticleConstituents()) {
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
 
  fastjet::AreaDefinition fAreaDef(fastjet::passive_area, ghost_spec);
  try {
    fastjet::ClusterSequenceArea fClustSeqSA(fInputVectors, fJetDef, fAreaDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets = fClustSeqSA.inclusive_jets(0);
  
    fastjet::PseudoJet jj;
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;

     std::vector<Double_t> delta_R_vec;
     std::vector<Double_t> xkt_vec;
     std::vector<Double_t> z_vec;
     std::vector<Double_t> rad_vec;
     std::vector<Double_t> eta1_vec;
     std::vector<Double_t> phi1_vec;
     std::vector<Double_t> eta2_vec;
     std::vector<Double_t> phi2_vec;
 
    jj = fOutputJets[0];
    while (jj.has_parents(j1, j2)) {
    
      if (j1.perp() < j2.perp())
        swap(j1, j2);
      
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      double rad = j1.e()+j2.e();
      double z = j2.perp() /jj.perp();
      double eta1=j1.eta();
      double eta2=j2.eta();
      double phi1=j1.phi();
      double phi2=j2.phi();
    
      delta_R_vec.push_back(delta_R);
      xkt_vec.push_back(xkt);
      z_vec.push_back(z);
      rad_vec.push_back(rad);
      eta1_vec.push_back(eta1);
      phi1_vec.push_back(phi1);
      eta2_vec.push_back(eta2);
      phi2_vec.push_back(phi2);
      
      jj = j1;
    }
    
          fShapesVar_Splittings_angle.push_back(delta_R_vec);
	  fShapesVar_Splittings_kt.push_back(xkt_vec); 
	  fShapesVar_Splittings_z.push_back(z_vec);
	  fShapesVar_Splittings_energy.push_back(rad_vec);
	  fShapesVar_Splittings_eta1.push_back(eta1_vec);
	  fShapesVar_Splittings_phi1.push_back(phi1_vec);
          fShapesVar_Splittings_eta2.push_back(eta2_vec);
	  fShapesVar_Splittings_phi2.push_back(phi2_vec);

	  delta_R_vec.clear();
	   xkt_vec.clear();
	    z_vec.clear();

            rad_vec.clear();
	    eta1_vec.clear();
	    eta2.clear();
            phi1_vec.clear();
	    phi2.clear();
	    
	    
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}
//_________________________________________________________________________
void AliAnalysisTaskLundPlane::IterativeParentsMCAverage(
    AliEmcalJet *fJet, Int_t km) {
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

     std::vector<Double_t> delta_R_vec;
     std::vector<Double_t> xkt_vec;
     std::vector<Double_t> z_vec;
     std::vector<Double_t> rad_vec;
     std::vector<Double_t> eta1_vec;
     std::vector<Double_t> phi1_vec;
     std::vector<Double_t> eta2_vec;
     std::vector<Double_t> phi2_vec;
 
    jj = fOutputJets[0];
    while (jj.has_parents(j1, j2)) {
    
      if (j1.perp() < j2.perp())
        swap(j1, j2);
      
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      double rad = j1.e()+j2.e();
      double z = j2.perp() /jj.perp();
      double eta1=j1.eta();
      double eta2=j2.eta();
      double phi1=j1.phi();
      double phi2=j2.phi();
    
      delta_R_vec.push_back(delta_R);
      xkt_vec.push_back(xkt);
      z_vec.push_back(z);
      rad_vec.push_back(rad);
      eta1_vec.push_back(eta1);
      phi1_vec.push_back(phi1);
      eta2_vec.push_back(eta2);
      phi2_vec.push_back(phi2);
      
      jj = j1;
    }
    
          fShapesVar_Splittings_angle_part.push_back(delta_R_vec);
	  fShapesVar_Splittings_kt_part.push_back(xkt_vec); 
	  fShapesVar_Splittings_z_part.push_back(z_vec);
	  fShapesVar_Splittings_energy_part.push_back(rad_vec);
	  fShapesVar_Splittings_eta1_part.push_back(eta1_vec);
	  fShapesVar_Splittings_phi1_part.push_back(phi1_vec);
          fShapesVar_Splittings_eta2_part.push_back(eta2_vec);
	  fShapesVar_Splittings_phi2_part.push_back(phi2_vec);

	  delta_R_vec.clear();
	   xkt_vec.clear();
	    z_vec.clear();

            rad_vec.clear();
	    eta1_vec.clear();
	    eta2.clear();
            phi1_vec.clear();
	    phi2.clear();

 

      



  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}




//________________________________________________________________________
Bool_t AliAnalysisTaskLundPlane::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskLundPlane::Terminate(Option_t *) {
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}
