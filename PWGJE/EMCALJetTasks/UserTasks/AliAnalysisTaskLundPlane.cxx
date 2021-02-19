//
// Task for the Lund plane in pp collisions
//authors: Leticia Cunqueiro, Laura Havener
//
//
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetContainer.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
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
  fStoreDetLevelJets(0), fDoSubJet(0),fTreeSplittings(0), fShapesVar_Splittings_angle(0),
    fShapesVar_Splittings_kt(0),fShapesVar_Splittings_z(0),fShapesVar_Splittings_eta1(0),
    fShapesVar_Splittings_eta2(0),fShapesVar_Splittings_phi1(0),fShapesVar_Splittings_phi2(0),
    fShapesVar_Splittings_angle_part(0),
    fShapesVar_Splittings_kt_part(0),fShapesVar_Splittings_z_part(0),fShapesVar_Splittings_eta1_part(0),
  fShapesVar_Splittings_eta2_part(0),fShapesVar_Splittings_phi1_part(0),fShapesVar_Splittings_phi2_part(0),fShapesVar_Splittings_ptjet(0),fShapesVar_Splittings_ptjet_part(0), fShapesVar_Splittings_mytrig(0), fShapesVar_Splittings_weightmb(0), fShapesVar_Splittings_weightej1(0), fShapesVar_Splittings_weightej2(0), fMatch(kFALSE), fTreeMatching(0), fHtrueMatch(0x0), fHtrueAll(0x0), fHrecoMatch(0x0), fHrecoAll(0x0), fHtrueMatch1D(0x0), fHtrueAll1D(0x0), fMatchR(0.1), fMomFrac(0.5), fStoreTrig(kFALSE)
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
    fStoreDetLevelJets(0),fDoSubJet(0),fTreeSplittings(0), fShapesVar_Splittings_angle(0),
    fShapesVar_Splittings_kt(0),fShapesVar_Splittings_z(0),fShapesVar_Splittings_eta1(0),
    fShapesVar_Splittings_eta2(0),fShapesVar_Splittings_phi1(0),fShapesVar_Splittings_phi2(0),
    fShapesVar_Splittings_angle_part(0),
    fShapesVar_Splittings_kt_part(0),fShapesVar_Splittings_z_part(0),fShapesVar_Splittings_eta1_part(0),
    fShapesVar_Splittings_eta2_part(0),fShapesVar_Splittings_phi1_part(0),fShapesVar_Splittings_phi2_part(0),fShapesVar_Splittings_ptjet(0),fShapesVar_Splittings_ptjet_part(0), fShapesVar_Splittings_mytrig(0), fShapesVar_Splittings_weightmb(0), fShapesVar_Splittings_weightej1(0), fShapesVar_Splittings_weightej2(0),fMatch(kFALSE), fTreeMatching(0), fHtrueMatch(0x0), fHtrueAll(0x0), fHrecoMatch(0x0), fHrecoAll(0x0), fHtrueMatch1D(0x0), fHtrueAll1D(0x0), fMatchR(0.1), fMomFrac(0.5), fStoreTrig(kFALSE)
    
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
  
  
  if (!fMatch)
    {
      fTreeSplittings = new TTree(nameoutput, nameoutput);
      TString *fShapesVarNames_Splittings=new TString[20];

      fShapesVarNames_Splittings[0] = "angle";
      fShapesVarNames_Splittings[1] = "kt";
      fShapesVarNames_Splittings[2] = "z";
    
      fShapesVarNames_Splittings[3] = "eta1";
      fShapesVarNames_Splittings[4] = "phi1";
      fShapesVarNames_Splittings[5] = "eta2";
      fShapesVarNames_Splittings[6] = "phi2";
      fShapesVarNames_Splittings[7] = "angle_part";
      fShapesVarNames_Splittings[8] = "kt_part";
      fShapesVarNames_Splittings[9] = "z_part";
   
      fShapesVarNames_Splittings[10] = "eta1_part";
      fShapesVarNames_Splittings[11] = "phi1_part";
      fShapesVarNames_Splittings[12] = "eta2_part";
      fShapesVarNames_Splittings[13] = "phi2_part";
      fShapesVarNames_Splittings[14] = "ptjet";
      fShapesVarNames_Splittings[15] = "ptjet_part"; 
      fShapesVarNames_Splittings[16] = "mytrig";
      fShapesVarNames_Splittings[17] = "weightmb";
       fShapesVarNames_Splittings[18] = "weightej1";
        fShapesVarNames_Splittings[19] = "weightej2";
	
      fTreeSplittings->Branch(fShapesVarNames_Splittings[0].Data(), &fShapesVar_Splittings_angle, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[1].Data(), &fShapesVar_Splittings_kt, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[2].Data(), &fShapesVar_Splittings_z, 0,1);
      
      fTreeSplittings->Branch(fShapesVarNames_Splittings[3].Data(), &fShapesVar_Splittings_eta1, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[4].Data(), &fShapesVar_Splittings_phi1, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[5].Data(), &fShapesVar_Splittings_eta2, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[6].Data(), &fShapesVar_Splittings_phi2, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[7].Data(), &fShapesVar_Splittings_angle_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[8].Data(), &fShapesVar_Splittings_kt_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[9].Data(), &fShapesVar_Splittings_z_part, 0,1);
      
      fTreeSplittings->Branch(fShapesVarNames_Splittings[10].Data(), &fShapesVar_Splittings_eta1_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[11].Data(), &fShapesVar_Splittings_phi1_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[12].Data(), &fShapesVar_Splittings_eta2_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[13].Data(), &fShapesVar_Splittings_phi2_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[14].Data(), &fShapesVar_Splittings_ptjet, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[15].Data(), &fShapesVar_Splittings_ptjet_part, 0,1);
      fTreeSplittings->Branch(fShapesVarNames_Splittings[16].Data(), &fShapesVar_Splittings_mytrig, 0,1);
       fTreeSplittings->Branch(fShapesVarNames_Splittings[17].Data(), &fShapesVar_Splittings_weightmb, 0,1);
        fTreeSplittings->Branch(fShapesVarNames_Splittings[18].Data(), &fShapesVar_Splittings_weightej1, 0,1);
	 fTreeSplittings->Branch(fShapesVarNames_Splittings[19].Data(), &fShapesVar_Splittings_weightej2, 0,1);
    }
  else {
    fTreeMatching = new TTree(nameoutput, nameoutput);
    int N = 0;
    if (!fDoSubJet) N = 7;
    else N = 9;
    TString *fShapesVarNames_Matching = new TString[N];
    fShapesVarNames_Matching[0] = "ptjet";
    fShapesVarNames_Matching[1] = "lnkt";
    fShapesVarNames_Matching[2] = "lnR";
    fShapesVarNames_Matching[3] = "ptjet_part";
    fShapesVarNames_Matching[4] = "lnkt_part";
    fShapesVarNames_Matching[5] = "lnR_part";
    fShapesVarNames_Matching[6] = "mytrig";
    fTreeMatching->Branch(fShapesVarNames_Matching[0].Data(), &fShapesVar_Matching_ptjet, 0,1);
    fTreeMatching->Branch(fShapesVarNames_Matching[1].Data(), &fShapesVar_Matching_lnkt, 0,1);
    fTreeMatching->Branch(fShapesVarNames_Matching[2].Data(), &fShapesVar_Matching_lnR, 0,1);
    fTreeMatching->Branch(fShapesVarNames_Matching[3].Data(), &fShapesVar_Matching_ptjet_part, 0,1);
    fTreeMatching->Branch(fShapesVarNames_Matching[4].Data(), &fShapesVar_Matching_lnkt_part, 0,1);
    fTreeMatching->Branch(fShapesVarNames_Matching[5].Data(), &fShapesVar_Matching_lnR_part, 0,1);
     fTreeMatching->Branch(fShapesVarNames_Matching[6].Data(), &fShapesVar_Matching_mytrig, 0,1);

    if (fDoSubJet) {
      fShapesVarNames_Matching[7] = "sub1";
      fShapesVarNames_Matching[8] = "sub2";
      fTreeMatching->Branch(fShapesVarNames_Matching[7].Data(), &fShapesVar_Matching_sub1, 0,1);
      fTreeMatching->Branch(fShapesVarNames_Matching[8].Data(), &fShapesVar_Matching_sub2, 0,1);
    }

    const Double_t ptbins_true[13] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 340};
    const Double_t ptbins_reco[11] = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220};
    const Double_t Rbins_true[9] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.5};
    const Double_t Rbins_reco[8] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
    const Double_t ktbins_true[15] = {-3, -1, -0.75, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3};
    const Double_t ktbins_reco[13] = {-1, -0.75, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2};    

    fHtrueMatch1D = new TH1D("fHtrueMatch1D", "matched truth jets", 12, ptbins_true);
    fOutput->Add(fHtrueMatch1D);
    fHtrueAll1D = new TH1D("fHtrueAll1D", "all truth jets", 12, ptbins_true);
    fOutput->Add(fHtrueAll1D);
    fHtrueMatch = new TH3D("fHtrueMatch", "matched truth splittings", 8, Rbins_true, 14, ktbins_true, 12, ptbins_true);
    fOutput->Add(fHtrueMatch);
    fHtrueAll = new TH3D("fHtrueAll", "all truth splitting", 8, Rbins_true, 14, ktbins_true, 12, ptbins_true);
    fOutput->Add(fHtrueAll);
    fHrecoMatch = new TH3D("fHrecoMatch", "matched reco splittings", 7, Rbins_reco, 12, ktbins_reco, 10, ptbins_reco);
    fOutput->Add(fHrecoMatch);
    fHrecoAll = new TH3D("fHrecoAll", "allreco splittings", 7, Rbins_reco, 12, ktbins_reco, 10, ptbins_reco);
    fOutput->Add(fHrecoAll);
  }

  PostData(1, fOutput);
  if (!fMatch) PostData(2, fTreeSplittings);
  else PostData(2, fTreeMatching);  
  
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


   
     Int_t mytrig=-1; 
      Bool_t mytrigmb=kFALSE;
      Bool_t mytrigej1=kFALSE;
      Bool_t mytrigej2=kFALSE;
      Double_t weightmb=-1;
      Double_t weightej1=-1;
      Double_t weightej2=-1;
      if(fStoreTrig==kTRUE){
        
        if(fInputHandler->IsEventSelected() & AliVEvent::kINT7 && fInputEvent->GetFiredTriggerClasses().Contains("INT7")) mytrigmb=kTRUE;
        if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE && fInputEvent->GetFiredTriggerClasses().Contains("EJ1")) mytrigej1=kTRUE;
	if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE && fInputEvent->GetFiredTriggerClasses().Contains("EJ2")) mytrigej2=kTRUE;

	if(mytrigmb==kTRUE && mytrigej1==kTRUE && mytrigej2==kTRUE) mytrig=3;
	if(mytrigmb==kTRUE && mytrigej1==kFALSE && mytrigej2==kFALSE) mytrig=0;
		if(mytrigmb==kFALSE && mytrigej1==kTRUE && mytrigej2==kFALSE) mytrig=1;
			if(mytrigmb==kFALSE && mytrigej1==kFALSE && mytrigej2==kTRUE) mytrig=2;
			if(mytrigmb==kTRUE && mytrigej1==kTRUE && mytrigej2==kFALSE) mytrig=4;
                 	if(mytrigmb==kTRUE && mytrigej1==kFALSE && mytrigej2==kTRUE) mytrig=5;
				if(mytrigmb==kFALSE && mytrigej1==kTRUE && mytrigej2==kTRUE) mytrig=6;
       if(mytrig==-1) return 0;

      }

     


    UInt_t newrun=InputEvent()->GetRunNumber();  
  if(fStoreTrig==kTRUE) {
    
    RunChanged(newrun); 

    if(mytrig==3){
      weightmb = 1./GetDownscaleWeight("INT7");
        weightej1 = 1./GetDownscaleWeight("EJ1");
	  weightej2 = 1./GetDownscaleWeight("EJ2");}

 if(mytrig==0){
   weightmb = 1./GetDownscaleWeight("INT7");}
    
  if(mytrig==1){
          weightej1 = 1./GetDownscaleWeight("EJ1");}

 if(mytrig==2){
          weightej2 = 1./GetDownscaleWeight("EJ2");}

if(mytrig==4){
       weightmb = 1./GetDownscaleWeight("INT7");
      weightej1 = 1./GetDownscaleWeight("EJ1");}         


if(mytrig==5){
       weightmb = 1./GetDownscaleWeight("INT7");
      weightej2 = 1./GetDownscaleWeight("EJ2");}  



 if(mytrig==6){
       weightej1 = 1./GetDownscaleWeight("EJ1");
      weightej2 = 1./GetDownscaleWeight("EJ2");}  

  }
      
      
    while ((jet1 = jetCont->GetNextAcceptJet())) {
      if (!jet1)
        continue;
      AliEmcalJet *jet2 = 0x0;
      AliEmcalJet *jet3 = 0x0;
      
      AliEmcalJet *jetUS = NULL;
      Int_t ifound = 0, jfound = 0;
      Int_t ilab = -1, jlab = -1;

     

      

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

     
	
      std::vector<std::vector<fastjet::PseudoJet>* >* constPart = new std::vector<std::vector<fastjet::PseudoJet>* >();
      std::vector<std::vector<fastjet::PseudoJet>* >* constDet = new std::vector<std::vector<fastjet::PseudoJet>* >();
      std::vector<fastjet::PseudoJet>* const1Part = new std::vector<fastjet::PseudoJet>();
      std::vector<fastjet::PseudoJet>* const1Det = new std::vector<fastjet::PseudoJet>();


	    
      IterativeDeclustering(jet1, jetCont, const1Det, constDet);
    
    
      
      
      Double_t ptMatch=0;
      Int_t kMatched = 0;
    
      if (fJetShapeType == kPythiaDef) {
        kMatched = 1;
        if (fJetShapeSub == kConstSub)
          kMatched = 3;

        ptMatch = jet3->Pt();

        IterativeDeclusteringMC(jet3, kMatched, const1Part, constPart);
       
      }
      
      fShapesVar_Splittings_ptjet=ptSubtracted;
      fShapesVar_Splittings_ptjet_part=ptMatch;
      fShapesVar_Splittings_mytrig=mytrig;
      fShapesVar_Splittings_weightmb=weightmb;
       fShapesVar_Splittings_weightej1=weightej1;
        fShapesVar_Splittings_weightej2=weightej2;
      
      if (fMatch) {
	Bool_t matched = SubjetMatching(const1Part,  constPart,  const1Det,  constDet);
      }
      else {	
	fTreeSplittings->Fill();
      }

      fShapesVar_Splittings_angle.clear();
      fShapesVar_Splittings_kt.clear(); 
      fShapesVar_Splittings_z.clear();
    
      fShapesVar_Splittings_eta1.clear();
      fShapesVar_Splittings_phi1.clear();
      fShapesVar_Splittings_eta2.clear();
      fShapesVar_Splittings_phi2.clear();
     
      fShapesVar_Splittings_angle_part.clear();
      fShapesVar_Splittings_kt_part.clear(); 
      fShapesVar_Splittings_z_part.clear();
    
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
void AliAnalysisTaskLundPlane::IterativeDeclustering(AliEmcalJet *fJet, AliJetContainer *fJetCont, std::vector < fastjet::PseudoJet > *const1, std::vector<std::vector < fastjet::PseudoJet > *> *constit) {

  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;
  unsigned int constituentIndex = 0;
  for (auto part: fJet->GetParticleConstituents()) {
    PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E());
    const AliVParticle* part2 = part.GetParticle();
    PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet));
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
    
     std::vector<Double_t> eta1_vec;
     std::vector<Double_t> phi1_vec;
     std::vector<Double_t> eta2_vec;
     std::vector<Double_t> phi2_vec;
 
    jj = fOutputJets[0];
    int index = 0;
    while (jj.has_parents(j1, j2)) {
    
      if (j1.perp() < j2.perp())
        swap(j1, j2);
      
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      //double rad = j1.e()+j2.e();
      double z = j2.perp() /jj.perp();
      double eta1=j1.eta();
      double eta2=j2.eta();
      double phi1=j1.phi();
      double phi2=j2.phi();
    
      delta_R_vec.push_back(delta_R);
      xkt_vec.push_back(xkt);
      z_vec.push_back(z);
   
      eta1_vec.push_back(eta1);
      phi1_vec.push_back(phi1);
      eta2_vec.push_back(eta2);
      phi2_vec.push_back(phi2);

      jj = j1;
      std::vector<fastjet::PseudoJet>* const2 = new std::vector<fastjet::PseudoJet>();
      if (j2.has_constituents()) *const2 = j2.constituents();
      constit->push_back(const2);
      if (index == 0) {
	if (j1.has_constituents()) *const1 = j1.constituents();
      }
      index++;
    }
    
          fShapesVar_Splittings_angle.push_back(delta_R_vec);
	  fShapesVar_Splittings_kt.push_back(xkt_vec); 
	  fShapesVar_Splittings_z.push_back(z_vec);

	  fShapesVar_Splittings_eta1.push_back(eta1_vec);
	  fShapesVar_Splittings_phi1.push_back(phi1_vec);
          fShapesVar_Splittings_eta2.push_back(eta2_vec);
	  fShapesVar_Splittings_phi2.push_back(phi2_vec);
	  

	  delta_R_vec.clear();
	   xkt_vec.clear();
	    z_vec.clear();

          
	    eta1_vec.clear();
	    eta2_vec.clear();
            phi1_vec.clear();
	    phi2_vec.clear();
	    
	    
  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}
//_________________________________________________________________________
void AliAnalysisTaskLundPlane::IterativeDeclusteringMC(
						       AliEmcalJet *fJet, Int_t km, std::vector < fastjet::PseudoJet > *const1, std::vector<std::vector < fastjet::PseudoJet > *> *constit) {
  AliJetContainer *jetCont = GetJetContainer(km);
  std::vector<fastjet::PseudoJet> fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet PseudoTracks;
  unsigned int constituentIndex = 0;
  for (auto  part: fJet->GetParticleConstituents()) {
    PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E());
    const AliVParticle* part2 = part.GetParticle();
    if (part2->Charge() == 0) continue;
    PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet));
    fInputVectors.push_back(PseudoTracks);
    constituentIndex++;
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
    
     std::vector<Double_t> eta1_vec;
     std::vector<Double_t> phi1_vec;
     std::vector<Double_t> eta2_vec;
     std::vector<Double_t> phi2_vec;
 
    jj = fOutputJets[0];
    int index = 0;
    while (jj.has_parents(j1, j2)) {
    
      if (j1.perp() < j2.perp())
        swap(j1, j2);
      
      double delta_R = j1.delta_R(j2);
      double xkt = j2.perp() * sin(delta_R);
      //double rad = j1.e()+j2.e();
      double z = j2.perp() /jj.perp();
      double eta1=j1.eta();
      double eta2=j2.eta();
      double phi1=j1.phi();
      double phi2=j2.phi();
    
      delta_R_vec.push_back(delta_R);
      xkt_vec.push_back(xkt);
      z_vec.push_back(z);
    
      eta1_vec.push_back(eta1);
      phi1_vec.push_back(phi1);
      eta2_vec.push_back(eta2);
      phi2_vec.push_back(phi2);
      
      jj = j1;

      std::vector<fastjet::PseudoJet>* const2 = new std::vector<fastjet::PseudoJet>();
      if (j2.has_constituents()) *const2 = j2.constituents();
      constit->push_back(const2);
      if (index == 0) {
        if (j1.has_constituents()) *const1 = j1.constituents();
      }
      index++;
    }
    
          fShapesVar_Splittings_angle_part.push_back(delta_R_vec);
	  fShapesVar_Splittings_kt_part.push_back(xkt_vec); 
	  fShapesVar_Splittings_z_part.push_back(z_vec);
	 
	  fShapesVar_Splittings_eta1_part.push_back(eta1_vec);
	  fShapesVar_Splittings_phi1_part.push_back(phi1_vec);
          fShapesVar_Splittings_eta2_part.push_back(eta2_vec);
	  fShapesVar_Splittings_phi2_part.push_back(phi2_vec);

	  delta_R_vec.clear();
	   xkt_vec.clear();
	    z_vec.clear();

           
	    eta1_vec.clear();
	    eta2_vec.clear();
            phi1_vec.clear();
	    phi2_vec.clear();

 

      



  } catch (fastjet::Error) {
    AliError(" [w] FJ Exception caught.");
    // return -1;
  }

  return;
}

//________________________________________________________________________                                                                                          
Bool_t AliAnalysisTaskLundPlane::SubjetMatching(std::vector < fastjet::PseudoJet > *constPart1, std::vector<std::vector < fastjet::PseudoJet > *> *constPart, std::vector < fastjet::PseudoJet > *constDet1, std::vector<std::vector < fastjet::PseudoJet > *> *constDet)
{
  fHtrueAll1D->Fill(fShapesVar_Splittings_ptjet_part);
  if ((fShapesVar_Splittings_ptjet_part < 0) || (fShapesVar_Splittings_ptjet_part > 340.)) return kFALSE;
  
  std::vector<int> reco_matches;
  float ptsub1_det = 0;

  for (int i = 0; i < fShapesVar_Splittings_kt_part.at(0).size(); i++)
    {
      float lnkt_part = std::log(fShapesVar_Splittings_kt_part.at(0).at(i));
      float lnr_part = std::log(0.4/fShapesVar_Splittings_angle_part.at(0).at(i));
      if (lnkt_part < -3. || lnkt_part > 3) continue;
      if (lnr_part < 0. || lnr_part > 2.5) continue;

      fHtrueAll->Fill(lnr_part, lnkt_part, fShapesVar_Splittings_ptjet_part);

      float dR_max = fMatchR;
      int ind_true = -1;
      int ind_reco = -1;
      
      for (int j = 0; j < fShapesVar_Splittings_kt.at(0).size(); j++)
	{
	  float deta = fShapesVar_Splittings_eta2.at(0).at(j) - fShapesVar_Splittings_eta2_part.at(0).at(i);
	  float dphi = fShapesVar_Splittings_phi2.at(0).at(j) - fShapesVar_Splittings_phi2_part.at(0).at(i);
	  if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
	  float dR = std::sqrt(dphi*dphi + deta*deta);
	  if (dR < dR_max) {
	    dR_max = dR;
	    ind_true = i;
	    ind_reco = j;
	  }
	}
      if (ind_reco == -1) continue;

      int ind_true_det = -1;
      float dR_max_det = fMatchR;
      for (int l = 0; l < fShapesVar_Splittings_kt_part.at(0).size(); l++)
	{
	  float deta = fShapesVar_Splittings_eta2.at(0).at(ind_reco) - fShapesVar_Splittings_eta2_part.at(0).at(l);
          float dphi = fShapesVar_Splittings_phi2.at(0).at(ind_reco) - fShapesVar_Splittings_phi2_part.at(0).at(l);
          if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
          float dR = std::sqrt(dphi*dphi + deta*deta);
          if (dR < dR_max_det) {
            dR_max_det = dR;
            ind_true_det = l;
          }
        }
      if (ind_true_det == -1) continue;
      if (ind_true!=ind_true_det) continue;
      if ((fShapesVar_Splittings_ptjet > 220) || (fShapesVar_Splittings_ptjet < 20.)) continue;
      float lnkt_det = std::log(fShapesVar_Splittings_kt.at(0).at(ind_reco));
      float lnr_det = std::log(0.4/fShapesVar_Splittings_angle.at(0).at(ind_reco));
      if (lnkt_det < -1. || lnkt_det > 2) continue;
      if (lnr_det < 0. || lnr_det > 1.4) continue;
      reco_matches.push_back(ind_reco);            
      fHtrueMatch->Fill(lnr_part, lnkt_part, fShapesVar_Splittings_ptjet_part);

      fShapesVar_Matching_ptjet = fShapesVar_Splittings_ptjet;
      fShapesVar_Matching_mytrig = fShapesVar_Splittings_mytrig;
      fShapesVar_Matching_lnR = lnr_det;
      fShapesVar_Matching_lnkt =	lnkt_det;
      fShapesVar_Matching_ptjet_part	= fShapesVar_Splittings_ptjet_part;
      fShapesVar_Matching_lnR_part =	lnr_part;
      fShapesVar_Matching_lnkt_part =        lnkt_part;
      if (i == 0) ptsub1_det = (fShapesVar_Splittings_kt.at(0).at(ind_reco)/sin(fShapesVar_Splittings_angle.at(0).at(ind_reco)))*((1/fShapesVar_Splittings_z.at(0).at(ind_reco)) - 1);
      float ptsub2_det = fShapesVar_Splittings_kt.at(0).at(ind_reco)/sin(fShapesVar_Splittings_angle.at(0).at(ind_reco));
      if (fDoSubJet) {
	fShapesVar_Matching_sub1 = CompareSubjets(ptsub1_det, constDet1, constPart1, true);
	fShapesVar_Matching_sub2 = CompareSubjets(ptsub2_det, constDet->at(ind_reco), constPart->at(i), true);
      }
      
      fTreeMatching->Fill();
    }

  if ((fShapesVar_Splittings_ptjet > 220) || (fShapesVar_Splittings_ptjet < 20.)) return kFALSE;
  fHtrueMatch1D->Fill(fShapesVar_Splittings_ptjet_part);

   for (int i = 0; i < fShapesVar_Splittings_kt.at(0).size(); i++)
     {
       float lnkt_det = std::log(fShapesVar_Splittings_kt.at(0).at(i));
      float lnr_det = std::log(0.4/fShapesVar_Splittings_angle.at(0).at(i));
      if (lnkt_det < -1. || lnkt_det > 2) continue;
      if (lnr_det < 0. || lnr_det > 1.4) continue;

      fHrecoAll->Fill(lnr_det, lnkt_det, fShapesVar_Splittings_ptjet);

      bool match = false;
      for (int j = 0; j < reco_matches.size(); j++)    
	{ 
	  if (i == reco_matches.at(j)) match = true;	
	}                                                                             
      if (!match) continue;

      fHrecoMatch->Fill(lnr_det, lnkt_det, fShapesVar_Splittings_ptjet);
     }

   return kTRUE;
}

//________________________________________________________________________                                                                                          
Bool_t AliAnalysisTaskLundPlane::CompareSubjets(float pT_det, std::vector<fastjet::PseudoJet> *constDet, std::vector<fastjet::PseudoJet>* constHyb, bool matchTag)
{
  //  double pT_det = subDet->pt();
  double sumpT = 0;
  double delta =  0.01;

  for (int i = 0; i < constDet->size(); i++)
      {
        double eta_det = constDet->at(i).eta();
        double phi_det = constDet->at(i).phi();
	int ind_det = constDet->at(i).user_index();
        for (int j  = 0; j < constHyb->size(); j++)
          {
            double eta_hyb = constHyb->at(j).eta();
            double phi_hyb = constHyb->at(j).phi();
	    int ind_hyb = constHyb->at(j).user_index();
            double deta = eta_hyb - eta_det;
            deta = std::sqrt(deta*deta);
	    double dphi = phi_hyb - phi_det;
            dphi = std::sqrt(dphi*dphi);
	    if (!matchTag) {
	      if (deta > delta) continue;
	      if (dphi > delta) continue;
	    }
	    else {
	      if (ind_det != ind_hyb) continue;
	    }
            sumpT+=constDet->at(i).pt();
          }
      }
  if (sumpT/pT_det > fMomFrac) return true;
  else return false;
}

//________________________________________________________________________                                                                                          
int AliAnalysisTaskLundPlane::GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet)
{
  // NOTE: Usually, we would use the global offset defined for the general subtracter extraction task. But we don't want to
  //       depend on that task, so we just define it here locally.
  int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + 2000000);
  return id;
}
///_______________________________________________________________________
Double_t AliAnalysisTaskLundPlane::GetDownscaleWeight(string trigString)
{ 
  Double_t weight = 1.;
  TString triggerclass;
  if(trigString == "INT7") triggerclass = "CINT7-B-NOPF-CENT";
  else if(trigString == "EJ1") triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if(trigString == "EJ2") triggerclass = "CEMC7EJ2-B-NOPF-CENT";
  if(triggerclass.Length()) weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}
//////////_________________________________________________________________
void AliAnalysisTaskLundPlane::RunChanged(Int_t newrun){
  if(fStoreTrig) {
    auto downscalehandler = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    if(downscalehandler->GetCurrentRun() != newrun){
      downscalehandler->SetRun(newrun);
    }
  }
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
