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
#include "AliAnalysisTaskJetPull.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetPull)

    //________________________________________________________________________
AliAnalysisTaskJetPull::AliAnalysisTaskJetPull()
: AliAnalysisTaskEmcalJet("AliAnalysisTaskJetPull", kTRUE),
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
AliAnalysisTaskJetPull::AliAnalysisTaskJetPull(
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
AliAnalysisTaskJetPull::~AliAnalysisTaskJetPull() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetPull::UserCreateOutputObjects() {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fPtJet = new TH1F("fPtJet", "fPtJet", 100, 0, 200);
  fOutput->Add(fPtJet);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec = 7;
  const Int_t nBinsSpec[7] = {50, 100, 200, 20, 100, 50, 2};
  const Double_t lowBinSpec[7] = {0., -3, 0, 0, 0, 0, 0};
  const Double_t hiBinSpec[7] = {5., 2., 200, 20, 200, 50,2};
  fHLundIterative =
      new THnSparseF("fHLundIterative",
                     "LundIterativePlot [log(1/theta),log(z*theta),pTjet,algo]",
                     dimSpec, nBinsSpec, lowBinSpec, hiBinSpec);
  fOutput->Add(fHLundIterative);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec2 = 7;
  const Int_t nBinsSpec2[7] = {50, 100, 200, 20, 100, 50,2};
  const Double_t lowBinSpec2[7] = {0., -3, 0, 0, 0, 0, 0};
  const Double_t hiBinSpec2[7] = {5., 2., 200, 20, 200, 50, 2};
  fHLundIterativeMC = new THnSparseF(
      "fHLundIterativeMC",
      "LundIterativePlotMC [log(1/theta),log(z*theta),pTjet,algo]", dimSpec2,
      nBinsSpec2, lowBinSpec2, hiBinSpec2);
  fOutput->Add(fHLundIterativeMC);

  // log(1/theta),log(kt),jetpT,depth, tf, omega//
  const Int_t dimSpec3 = 7;
  const Int_t nBinsSpec3[7] = {50, 100, 200, 20, 100, 50, 2};
  const Double_t lowBinSpec3[7] = {0., -3, 0, 0, 0, 0,0};
  const Double_t hiBinSpec3[7] = {5., 2., 200, 20, 200, 50,2};
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
  const Int_t nVar = 10;
  const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeSubstructure = new TTree(nameoutput, nameoutput);
  TString *fShapesVarNames = new TString[nVar];

  fShapesVarNames[0] = "ptJet";
  fShapesVarNames[1] = "ktg";
  fShapesVarNames[2] = "pull";
  fShapesVarNames[3] = "zg";
  fShapesVarNames[4] = "rg";
  fShapesVarNames[5] = "ptJetMatch";
  fShapesVarNames[6] = "ktgMatch";
  fShapesVarNames[7] = "pullMatch";
  fShapesVarNames[8] = "zgMatch";
  fShapesVarNames[9] = "rgMatch";
 



  for (Int_t ivar = 0; ivar < nVar; ivar++) {
    fTreeSubstructure->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar],
                              Form("%s/F", fShapesVarNames[ivar].Data()));
  }

  PostData(1, fOutput);
  PostData(2, fTreeSubstructure);

  delete[] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPull::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPull::FillHistograms() {

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
     
      
      if(fCutDoubleCounts==kTRUE && fJetShapeType==kDetEmbPartPythia) if(jet1->MaxTrackPt()>jet2->MaxTrackPt()) continue;

    
      if ((fJetShapeType == kData && fCentSelectOn == false) || (fJetShapeType == kMCTrue) || (fJetShapeType == kPythiaDef)) IterativeParentsPP(jet1, jetCont);
      else IterativeParents(jet1, jetCont);
    
      Float_t ptMatch = 0.;
     
      Double_t ktgMatch = 0;
      Double_t pullMatch = 0;
      Double_t zgMatch = 0;
      Double_t rgMatch = 0;
     
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

        IterativeParentsMCAveragePP(jet3, kMatched, aver1, aver2, aver3, aver4);
        ktgMatch = aver1;
        pullMatch = aver2;
        zgMatch = aver3;
        rgMatch = aver4;
      }

      if (fJetShapeType == kDetEmbPartPythia) {
        if (fJetShapeSub == kConstSub)
          kMatched = 3;
        if (fJetShapeSub == kDerivSub)
          kMatched = 2;
        ptMatch = jet3->Pt();
      
        IterativeParentsMCAverage(jet3, kMatched, aver1, aver2, aver3, aver4);
        ktgMatch = aver1;
        pullMatch = aver2;
        zgMatch = aver3;
        rgMatch = aver4;
       
      }

      if (fJetShapeType == kMCTrue || fJetShapeType == kData ||
          fJetShapeType == kGenOnTheFly) {

        ptMatch = 0.;
     
        ktgMatch = 0.;
        pullMatch = 0.;
        zgMatch = 0;
        rgMatch = 0;
      }

      fShapesVar[5] = ptMatch;
      fShapesVar[6] = ktgMatch;
      fShapesVar[7] = pullMatch;
      fShapesVar[8] = zgMatch;
      fShapesVar[9] = rgMatch;
      fShapesVar[11] = leadTrackMatch;
          


      fTreeSubstructure->Fill();
   
    }
  }

  return kTRUE;
}




//__________________________________________________________________________________
Double_t AliAnalysisTaskJetPull::RelativePhi(Double_t mphi,
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
void AliAnalysisTaskJetPull::IterativeParents(
							 AliEmcalJet *fJet, AliJetContainer *fJetCont) {

  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;
 

   AliParticleContainer *fTrackCont = fJetCont->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;

      //if (fDoTwoTrack == kTRUE && CheckClosePartner(fJet, part)) continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(), fTrk->E());
      PseudoTracks.set_user_index(fJet->TrackAt(i) + 100);
    
   
    fInputVectors.push_back(PseudoTracks);
   
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

	j2first =j2;

        flagSubjet = 1;
      }
    
      
      Double_t LundEntries[7] = {
	y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, flagConst};
      fHLundIterative->Fill(LundEntries);
      
      jj = j1;
    }

     
    double disyjet=j2first.rap();
	 double disphijet=j2first.phi();
       double disty=0;
       double distphi=0;
     vector < fastjet::PseudoJet > constitfirst = sorted_by_pt(j1first.constituents());

     for(Int_t i=0;i<constitfirst.size();i++){
      double intery=constitfirst[i].rap()-j1first.rap();
      double interphi=constitfirst[i].phi()-j1first.phi();
      double deltar=TMath::Sqrt(intery*intery+interphi*interphi);
      
      disty=disty+intery*deltar*constitfirst[i].perp()/j1first.perp();
      distphi=distphi+interphi*deltar*constitfirst[i].perp()/j1first.perp();
       
       }
     double mod1=TMath::Sqrt(disty*disty+distphi*distphi);
      double mod2=TMath::Sqrt(disyjet*disyjet+disphijet*disphijet);
     double cospull= (disty*disyjet+distphi*disphijet)/(mod1*mod2);
     double pull=TMath::ACos(cospull);
     
    
    fShapesVar[1] = xktg;
    fShapesVar[2] = pull;
    fShapesVar[3] = zg;
    fShapesVar[4] = Rg;

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

void AliAnalysisTaskJetPull::IterativeParentsPP(
							 AliEmcalJet *fJet, AliJetContainer *fJetCont) {

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

	j2first =j2;

        flagSubjet = 1;
      }
     
      
      Double_t LundEntries[7] = {
	y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, flagConst};
      fHLundIterative->Fill(LundEntries);
      
      jj = j1;
    }
    
   
    double disyjet=j2first.rap();
	 double disphijet=j2first.phi();
       double disty=0;
       double distphi=0;
     vector < fastjet::PseudoJet > constitfirst = sorted_by_pt(j1first.constituents());
     for(Int_t i=0;i<constitfirst.size();i++){
      double intery=constitfirst[i].rap()-j1first.rap();
      double interphi=constitfirst[i].phi()-j1first.phi();
      double deltar=TMath::Sqrt(intery*intery+interphi*interphi);
      
      disty=disty+intery*deltar*constitfirst[i].perp()/j1first.perp();
      distphi=distphi+interphi*deltar*constitfirst[i].perp()/j1first.perp();
       
       }
     double mod1=TMath::Sqrt(disty*disty+distphi*distphi);
      double mod2=TMath::Sqrt(disyjet*disyjet+disphijet*disphijet);
     double cospull= (disty*disyjet+distphi*disphijet)/(mod1*mod2);
     double pull=TMath::ACos(cospull);
    fShapesVar[1] = xktg;
    fShapesVar[2] = pull;
    fShapesVar[3] = zg;
    fShapesVar[4] = Rg;

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

//_________________________________________________________________________
void AliAnalysisTaskJetPull::IterativeParentsMCAverage(
    AliEmcalJet *fJet, Int_t km, Double_t &average1, Double_t &average2,
    Double_t &average3, Double_t &average4) {
  AliJetContainer *jetCont = GetJetContainer(km);
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;

     AliParticleContainer *fTrackCont = jetCont->GetParticleContainer();

  if (fTrackCont)
    for (Int_t i = 0; i < fJet->GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = fJet->TrackAt(i, fTrackCont->GetArray());
      if (!fTrk) continue;

    
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
    fastjet::PseudoJet j1first;
    fastjet::PseudoJet j2first;
    jj = fOutputJets[0];
    int flagSubjet = 0;
  
     double flagConst=0;
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
	j1first = j1;

	j2first = j2;

        flagSubjet = 1;
      }
     
      if (fDoFillMCLund == kTRUE) {
        Double_t LundEntries[7] = {
            y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, flagConst};
        fHLundIterativeMC->Fill(LundEntries);
        if (fStoreDetLevelJets) {
          fHLundIterativeMCDet->Fill(LundEntries);
        }
      }

      jj = j1;
    }


    double disyjet=j2first.rap();
	 double disphijet=j2first.phi();
       double disty=0;
       double distphi=0;
     vector < fastjet::PseudoJet > constitfirst = sorted_by_pt(j1first.constituents());
     for(Int_t i=0;i<constitfirst.size();i++){
      double intery=constitfirst[i].rap()-j1first.rap();
      double interphi=constitfirst[i].phi()-j1first.phi();
      double deltar=TMath::Sqrt(intery*intery+interphi*interphi);
      
      disty=disty+intery*deltar*constitfirst[i].perp()/j1first.perp();
      distphi=distphi+interphi*deltar*constitfirst[i].perp()/j1first.perp();
       
       }
     double mod1=TMath::Sqrt(disty*disty+distphi*distphi);
      double mod2=TMath::Sqrt(disyjet*disyjet+disphijet*disphijet);
     double cospull= (disty*disyjet+distphi*disphijet)/(mod1*mod2);
     double pull=TMath::ACos(cospull);

    
 
    average1 = xktg;
    average2 = pull;
    average3 = zg;
    average4 = Rg;
  

  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

//_________________________________________________________________________
void AliAnalysisTaskJetPull::IterativeParentsMCAveragePP(
    AliEmcalJet *fJet, Int_t km, Double_t &average1, Double_t &average2,
    Double_t &average3, Double_t &average4) {
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
    double flagConst=0;
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
	j1first = j1;

	j2first = j2;

        flagSubjet = 1;
      }
      
      if (fDoFillMCLund == kTRUE) {
        Double_t LundEntries[7] = {
            y, lnpt_rel, fOutputJets[0].perp(), nall, form, rad, flagConst};
        fHLundIterativeMC->Fill(LundEntries);
        if (fStoreDetLevelJets) {
          fHLundIterativeMCDet->Fill(LundEntries);
        }
      }

      jj = j1;
    }


     
    double disyjet=j2first.rap();
	 double disphijet=j2first.phi();
       double disty=0;
       double distphi=0;
     vector < fastjet::PseudoJet > constitfirst = sorted_by_pt(j1first.constituents());
     for(Int_t i=0;i<constitfirst.size();i++){
      double intery=constitfirst[i].rap()-j1first.rap();
      double interphi=constitfirst[i].phi()-j1first.phi();
      double deltar=TMath::Sqrt(intery*intery+interphi*interphi);
      
      disty=disty+intery*deltar*constitfirst[i].perp()/j1first.perp();
      distphi=distphi+interphi*deltar*constitfirst[i].perp()/j1first.perp();
       
       }
     double mod1=TMath::Sqrt(disty*disty+distphi*distphi);
      double mod2=TMath::Sqrt(disyjet*disyjet+disphijet*disphijet);
     double cospull= (disty*disyjet+distphi*disphijet)/(mod1*mod2);
     double pull=TMath::ACos(cospull);




    
         
    average1 = xktg;
    average2 = pull;
    average3 = zg;
    average4 = Rg;
   
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetPull::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskJetPull::Terminate(Option_t *) {
  // Called once at the end of the analysis.

  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}
