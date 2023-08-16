//Task for EEC in pp collisions
//Code inherited and modified from AliAnalysisTaskJetsEEC.cxx authored by Leticia Cunqueiro and Laura Havener
//Authors: Ananya Rai, Laura Havener
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
#include <TMath.h>
#include "AliAODEvent.h"
#include "AliAnalysisTaskJetsEEC.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetsEEC)

//________________________________________________________________________
AliAnalysisTaskJetsEEC::AliAnalysisTaskJetsEEC(): AliAnalysisTaskEmcalJet("AliAnalysisTaskJetsEEC", kTRUE),
fContainer(0), fMinFractionShared(0), fJetShapeType(kData),
fJetShapeSub(kNoSub), fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.02), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fStoreTrig(kFALSE), fpTcorr(0), fpaircut(0), jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), EEC_pt_hist_log(0), E3C_hist(0), E3C_pt_hist(0), E3C_pt_hist_log(0), EEC_det_pt_hist_3d(0), EEC_tru_pt_hist_3d(0), EEC_det_pt_hist_log_3d(0), EEC_tru_pt_hist_log_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), E3C_det_pt_hist_log_3d(0), E3C_tru_pt_hist_log_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), N2_det_pt_hist_log_3d(0), N2_tru_pt_hist_log_3d(0), N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), EEC_det_match_pt_det_log(0), EEC_tru_match_pt_tru_log(0), E3C_det_match_pt_det_log(0), E3C_tru_match_pt_tru_log(0), pt_tru(0), pt_tru_match(0), pt_det(0), pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JES_scaled(0), JER(0), pair_det_EEC(0), pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0)
{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEEC::AliAnalysisTaskJetsEEC(const char *name): AliAnalysisTaskEmcalJet(name, kTRUE), fContainer(0),
fMinFractionShared(0), fJetShapeType(kData), fJetShapeSub(kNoSub),
fJetSelection(kInclusive), fPtThreshold(-9999.), fMinENCtrackPt(1.0), fCentSelectOn(kTRUE), fCentMin(0), fCentMax(10),
fOneConstSelectOn(kFALSE), fTrackCheckPlots(kFALSE), fCheckResolution(kFALSE),
fMinPtConst(1), fHardCutoff(0), fDoTwoTrack(kFALSE), fCutDoubleCounts(kTRUE),
fPowerAlgo(1), fPhiCutValue(0.02),
fEtaCutValue(0.02), fDerivSubtrOrder(0),
fStoreDetLevelJets(0), fDoFillEncMC(kFALSE), fStoreTrig(kFALSE), fpTcorr(0), fpaircut(0), jet_pt_hist(0), EEC_hist(0), EEC_pt_hist(0), EEC_pt_hist_log(0), E3C_hist(0), E3C_pt_hist(0), E3C_pt_hist_log(0), EEC_det_pt_hist_3d(0), EEC_tru_pt_hist_3d(0), EEC_det_pt_hist_log_3d(0), EEC_tru_pt_hist_log_3d(0), E3C_det_pt_hist_3d(0), E3C_tru_pt_hist_3d(0), E3C_det_pt_hist_log_3d(0),E3C_tru_pt_hist_log_3d(0), N2_det_pt_hist_3d(0), N2_tru_pt_hist_3d(0), N2_det_pt_hist_log_3d(0), N2_tru_pt_hist_log_3d(0), N3_det_pt_hist_3d(0), N3_tru_pt_hist_3d(0), EEC_det_match_pt_det(0), EEC_tru_match_pt_tru(0), E3C_det_match_pt_det(0), E3C_tru_match_pt_tru(0), EEC_det_match_pt_det_log(0), EEC_tru_match_pt_tru_log(0), E3C_det_match_pt_det_log(0), E3C_tru_match_pt_tru_log(0), pt_tru(0), pt_tru_match(0), pt_det(0), pt_det_match(0), test_hist(0), R_matrix(0), JES(0), JER(0), pair_det_EEC(0), pair_tru_EEC(0), pair_det_E3C(0), pair_tru_E3C(0)

{
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskJetsEEC::~AliAnalysisTaskJetsEEC() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetsEEC::UserCreateOutputObjects() {

// Create user output.
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    TH1::AddDirectory(oldStatus);
    
    const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
    
    Double_t from = -4;
    Double_t to = 0;
    Int_t bins = 100;
    Double_t width = (to-from)/bins;
    Double_t new_bins[101] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
    }

    Double_t from_const = 15;
    Double_t to_const = 120;
    Int_t bins_const = 21;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t new_bins_const[22] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
    new_bins_const[i] = (from_const + i * width_const);
    }

    
    jet_pt_hist = new TH1D("jet_pt_hist", "Jet Pt", 21, 15, 120);
    fOutput->Add(jet_pt_hist);
    
    //EEC (data or det level)
    EEC_hist = new TH1D("EEC_hist","EEC", 100, new_bins);
    fOutput->Add(EEC_hist);
    
    EEC_pt_hist = new TH2D("EEC_pt_hist", "EEC and jet_pt 2D", 100, new_bins, 21, 15, 120);
    fOutput->Add(EEC_pt_hist);
    
    EEC_pt_hist_log = new TH2D("EEC_pt_hist_log", "EEC and jet_pt 2D", 100, -10,10,21, 15, 120);
    fOutput->Add(EEC_pt_hist_log);
    
    //EEEC histograms (data or det level)
    E3C_hist = new TH1D("E3C_hist","E3C", 100, new_bins);
    fOutput->Add(E3C_hist);
    
    E3C_pt_hist = new TH2D("E3C_pt_hist", "EEEC and jet_pt 2D", 100, new_bins, 21, 15, 120);
    fOutput->Add(E3C_pt_hist);
    
    E3C_pt_hist_log = new TH2D("E3C_pt_hist_log", "EEEC and jet_pt 2D", 100, -10,10,21, 15, 120);
    fOutput->Add(E3C_pt_hist_log);
    
    //MC true and det level histograms 3d (det level fills through the data loop)
    EEC_det_pt_hist_3d = new TH3D("EEC_det_pt_hist", "EEC det and jet_pt 3D", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(EEC_det_pt_hist_3d);
    
    EEC_tru_pt_hist_3d = new TH3D("EEC_tru_pt_hist", "EEC tru and jet_pt 3d", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(EEC_tru_pt_hist_3d);
    
    EEC_det_pt_hist_log_3d = new TH3D("EEC_det_pt_hist_log", "EEC det and jet_pt 3d",100, -10,10, 21, 15, 120,21, 15, 120);
    fOutput->Add(EEC_det_pt_hist_log_3d);
    
    EEC_tru_pt_hist_log_3d = new TH3D("EEC_tru_pt_hist_log", "EEC tru and jet_pt 3d",100, -10,10, 21, 15, 120,21, 15, 120);
    fOutput->Add(EEC_tru_pt_hist_log_3d);
    
    E3C_det_pt_hist_3d = new TH3D("E3C_det_pt_hist", "E3C det and jet_pt 3D", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(E3C_det_pt_hist_3d);
    
    E3C_tru_pt_hist_3d = new TH3D("E3C_tru_pt_hist", "E3C tru and jet_pt 3d)", 100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(E3C_tru_pt_hist_3d);
    
    E3C_det_pt_hist_log_3d = new TH3D("E3C_det_pt_hist_log", "E3C det and jet_pt 3d",100, -10,10,21, 15, 120,21, 15, 120);
    fOutput->Add(E3C_det_pt_hist_log_3d);
    
    E3C_tru_pt_hist_log_3d = new TH3D("E3C_tru_pt_hist_log", "E3C tru and jet_pt 3d",100, -10,10,21, 15, 120,21, 15, 120);
    fOutput->Add(E3C_tru_pt_hist_log_3d);
    
    //Num of pairs at det and true level with det pt and tru pt
    N2_det_pt_hist_3d = new TH3D("N2_det_pt_hist", "Num pairs det and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(N2_det_pt_hist_3d);
    
    N2_tru_pt_hist_3d = new TH3D("N2_tru_pt_hist", "Num pairs tru and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(N2_tru_pt_hist_3d);
    
    N3_det_pt_hist_3d = new TH3D("N3_det_pt_hist", "Num pairs det and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(N3_det_pt_hist_3d);
    
    N3_tru_pt_hist_3d = new TH3D("N3_tru_pt_hist", "Num pairs tru and jet_pt 3d",100, new_bins, 21, new_bins_const, 21, new_bins_const);
    fOutput->Add(N3_tru_pt_hist_3d);
    
    //MC true level histograms (detector level goes through the data loop so don't need det level histograms)
    EEC_tru_match_pt_tru = new TH2D("EEC_pt_tru_match_hist", "EEC and pT for tru matched jets", 100, new_bins, 21, 15, 120);
    fOutput->Add(EEC_tru_match_pt_tru);
    
    E3C_tru_match_pt_tru = new TH2D("E3C_pt_tru_match_hist", "E3C and pT for tru matched jets", 100, new_bins, 21, 15, 120);
    fOutput->Add(E3C_tru_match_pt_tru);

    EEC_tru_match_pt_tru_log = new TH2D("EEC_pt_tru_match_hist_log", "EEC and pT for tru matched jets",100, -10,10, 21, 15, 120);
    fOutput->Add(EEC_tru_match_pt_tru_log);
    
    E3C_tru_match_pt_tru_log = new TH2D("E3C_pt_tru_match_hist_log", "E3C and pT for tru matched jets",100, -10,10,21, 15, 120);
    fOutput->Add(E3C_tru_match_pt_tru_log);
    
    pt_tru = new TH1D("jet_pt_tru_hist", "Jet Pt", 21, 15, 120);
    fOutput->Add(pt_tru);
    
    test_hist = new TH1D("test_hist", "test R EEC", 100, 0, 1);
    fOutput->Add(test_hist);
    
    R_matrix = new TH2D("Response matrix", "Jet pT response matrix", 21, 15, 120, 21, 15, 120);
    fOutput->Add(R_matrix);
    
    JES = new TH2D("JES","Jet energy scale",21,15,120,200,-2,8);
    fOutput->Add(JES);
    
    JES_scaled = new TH2D("JES scaled","Jet energy scale for scaled det pT",21,15,120,200,-2,8);
    fOutput->Add(JES_scaled);
    
  PostData(1, fOutput);


}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEEC::Run() {
  // Run analysis code here, if needed. It will be executed before
  // FillHistograms().
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEEC::FillHistograms()
{
    //fill histogram goes through event loops
    AliEmcalJet *jet1 = NULL;
    AliJetContainer *jetCont = GetJetContainer(0);  //Get jet container
    // container zero is always the base container: the data container, the
    // embedded subtracted in the case of embedding or the detector level in case
    // of pythia.
    
    
    if (fCentSelectOn)
        if ((fCent > fCentMax) || (fCent < fCentMin))
            return 0;
    
    Float_t rhoVal = 0, rhoMassVal = 0.;
    if (jetCont)
    {
        jetCont->ResetCurrentID();
        if ((fJetShapeSub == kConstSub) || (fJetShapeSub == kDerivSub))
        {
            // rho
            AliRhoParameter *rhoParam = dynamic_cast<AliRhoParameter *>(
                                                                        InputEvent()->FindListObject("RhoSparseR020"));
            if (!rhoParam)
            {
                Printf("%s: Could not retrieve rho %s (some histograms will be filled "
                       "with zero)!",
                       GetName(), jetCont->GetRhoName().Data());
            }
            else
                rhoVal = rhoParam->GetVal();
            // rhom
            AliRhoParameter *rhomParam = dynamic_cast<AliRhoParameter *>(
                                                                         InputEvent()->FindListObject("RhoMassSparseR020"));
            if (!rhomParam)
            {
                Printf("%s: Could not retrieve rho_m %s (some histograms will be "
                       "filled with zero)!",
                       GetName(), jetCont->GetRhoMassName().Data());
            }
            else
                rhoMassVal = rhomParam->GetVal();
        }
        
        //Keeping triggers for full jet analysis (if desired)
        //     Int_t mytrig=-1;
        //     Bool_t mytrigmb=kFALSE;
        //     Bool_t mytrigej1=kFALSE;
        //     Bool_t mytrigej2=kFALSE;
        //     Double_t weightmb=-1;
        //     Double_t weightej1=-1;
        //     Double_t weightej2=-1;
        //
        //    if(fStoreTrig==kTRUE)
        //    {
        //    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7 && fInputEvent->GetFiredTriggerClasses().Contains("INT7")) mytrigmb=kTRUE;
        //    if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE && fInputEvent->GetFiredTriggerClasses().Contains("EJ1")) mytrigej1=kTRUE;
        //    if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE && fInputEvent->GetFiredTriggerClasses().Contains("EJ2")) mytrigej2=kTRUE;
        //    if(mytrigmb==kTRUE && mytrigej1==kTRUE && mytrigej2==kTRUE) mytrig=3;
        //    if(mytrigmb==kTRUE && mytrigej1==kFALSE && mytrigej2==kFALSE) mytrig=0;
        //    if(mytrigmb==kFALSE && mytrigej1==kTRUE && mytrigej2==kFALSE) mytrig=1;
        //    if(mytrigmb==kFALSE && mytrigej1==kFALSE && mytrigej2==kTRUE) mytrig=2;
        //    if(mytrigmb==kTRUE && mytrigej1==kTRUE && mytrigej2==kFALSE) mytrig=4;
        //    if(mytrigmb==kTRUE && mytrigej1==kFALSE && mytrigej2==kTRUE) mytrig=5;
        //    if(mytrigmb==kFALSE && mytrigej1==kTRUE && mytrigej2==kTRUE) mytrig=6;
        //    if(mytrig==-1) return 0;
        //    }
        //
        //  UInt_t newrun=InputEvent()->GetRunNumber();
        //
        //if(fStoreTrig==kTRUE)
        //  {
        //    RunChanged(newrun);
        //      if(mytrig==3)
        //    {
        //      weightmb = 1./GetDownscaleWeight("INT7");
        //        weightej1 = 1./GetDownscaleWeight("EJ1");
        //      weightej2 = 1./GetDownscaleWeight("EJ2");
        //    }
        //      if(mytrig==0)
        //    {
        //   weightmb = 1./GetDownscaleWeight("INT7");
        //    }
        //      if(mytrig==1)
        //      {
        //          weightej1 = 1./GetDownscaleWeight("EJ1");
        //      }
        //      if(mytrig==2)
        //      {
        //          weightej2 = 1./GetDownscaleWeight("EJ2");
        //      }
        //      if(mytrig==4)
        //      {
        //      weightmb = 1./GetDownscaleWeight("INT7");
        //      weightej1 = 1./GetDownscaleWeight("EJ1");
        //      }
        //      if(mytrig==5)
        //      {
        //    weightmb = 1./GetDownscaleWeight("INT7");
        //    weightej2 = 1./GetDownscaleWeight("EJ2");
        //      }
        //
        //      if(mytrig==6)
        //      {
        //      weightej1 = 1./GetDownscaleWeight("EJ1");
        //      weightej2 = 1./GetDownscaleWeight("EJ2");
        //      }
        //  }
        
        //Jet Loop
        while ((jet1 = jetCont->GetNextAcceptJet()))
        {
            if (!jet1) continue;
            if (jet1->Pt() < 15) {continue;} //Cuts on jet_pt
            if (fJetShapeType == kData)
            {
//            cout <<"yeah"<<endl;
            ComputeEEC(jet1, jetCont);//Computing the eec on the jet object
            }
            
            AliEmcalJet *jet2 = 0x0;
            AliEmcalJet *jet3 = 0x0;
            
            AliEmcalJet *jetUS = NULL;
            Int_t ifound = 0, jfound = 0;
            Int_t ilab = -1, jlab = -1;
            
            // All for MC. This is the mode to run over pythia to produce a det-part response. Here we have also added the constituent-subtraction case, but we don't use it normally in pp the matching is purely geometrical
            if (fJetShapeType == kPythiaDef)
            {
                AliJetContainer *jetContTrue = GetJetContainer(1); //we never use this in pp, usually an empty one
                AliJetContainer *jetContUS = GetJetContainer(2); //unsubtracted one, don't need for pp
                AliJetContainer *jetContPart = GetJetContainer(3); //this is pythia particle/true level
                if (fJetShapeSub == kConstSub)
                {
                    for (Int_t i = 0; i < jetContUS->GetNJets(); i++)
                    {
                        jetUS = jetContUS->GetJet(i);
                        if (jetUS->GetLabel() == jet1->GetLabel())
                        {
                            ifound++;
                            if (ifound == 1) ilab = i;
                        }
                        
                    }
                    if (ilab == -1) continue;
                    jetUS = jetContUS->GetJet(ilab);
                    jet2 = jetUS->ClosestJet();
                    if (!jet2)
                    {
                        Printf("jet2 does not exist, returning");
                        continue;
                        
                    }
                    for (Int_t j = 0; j < jetContPart->GetNJets(); j++)
                    {
                        jet3 = jetContPart->GetJet(j);
                        if (!jet3) continue;
                        if (jet3->GetLabel() == jet2->GetLabel())
                        {
                            jfound++;
                            if (jfound == 1) jlab = j;
                        }
                        
                    }
                    if (jlab == -1) continue;
                    jet3 = jetContPart->GetJet(jlab);
                    if (!jet3)
                    {
                        Printf("jet3 does not exist, returning");
                        continue;
                    }
                    
                }
                if (!(fJetShapeSub == kConstSub))
                    jet3 = jet1->ClosestJet();
                if (!jet3)
                {
                    //   Printf("jet3 does not exist, returning");
                    continue;
                }
                
            }
            
            //These are options for pb analysis
            Double_t ptSubtracted = 0;
            if (fJetShapeSub == kConstSub || fJetShapeSub == kEventSub) ptSubtracted = jet1->Pt();
            else if (fJetShapeSub == kDerivSub)
            {
                ptSubtracted = jet1->Pt() - GetRhoVal(0) * jet1->Area();
            }
            else if (fJetShapeSub == kNoSub) ptSubtracted = jet1->Pt();
            if (ptSubtracted < fPtThreshold) continue;
            if ((fCentSelectOn == kFALSE) && (jet1->GetNumberOfTracks() <= 1)) continue;
            
            Double_t ptMatch=0;
            Int_t kMatched = 0;
            
            
            if (fJetShapeType == kPythiaDef)
            {
                kMatched = 1;
                if (fJetShapeSub == kConstSub)
                    kMatched = 3;
                ptMatch = jet3->Pt();
//                cout<<"the matched jet "<<jet3->Pt()<<" "<<kMatched<<endl;
                ComputeEncMC(jet1, jetCont, jet3, kMatched);

            }
            
        } //close while loop
    } //close the jet cont loop
    
    return kTRUE;
}



//________________________________________________________________________
int AliAnalysisTaskJetsEEC::GetConstituentID(int constituentIndex, const AliVParticle* part, AliEmcalJet * jet)
{
  // NOTE: Usually, we would use the global offset defined for the general subtracter extraction task. But we don't want to
  //       depend on that task, so we just define it here locally.
  int id = part->GetLabel() != -1 ? part->GetLabel() : (jet->TrackAt(constituentIndex) + 2000000);
  return id;
}


//_______________________________________________________________________
Double_t AliAnalysisTaskJetsEEC::GetDownscaleWeight(string trigString)
{
  Double_t weight = 1.;
  TString triggerclass;
  if(trigString == "INT7") triggerclass = "CINT7-B-NOPF-CENT";
  else if(trigString == "EJ1") triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if(trigString == "EJ2") triggerclass = "CEMC7EJ2-B-NOPF-CENT";
  if(triggerclass.Length()) weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}


//______________________________________________________________________
void AliAnalysisTaskJetsEEC::ComputeEncMC(AliEmcalJet *fJet, AliJetContainer *fJetCont, AliEmcalJet *fJet_tru, Int_t km)
    //(jet, jet container, vector of jets, vector of constituents within each jet?)
{
    //Det level
    std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
    fConstituents.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex = 0;
    //The line below gets constituent particles within fjet. C++ syntax[ for (auto elem : container)    // capture elements by value ]
    for (auto part: fJet->GetParticleConstituents())
    {
        PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet)); //leave this in for the same reason as above
        if (PseudoTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents.push_back(PseudoTracks);
        constituentIndex++;
    }
    
    //Truth level
    AliJetContainer *jetCont = GetJetContainer(km);
    std::vector<fastjet::PseudoJet> fConstituents_tru; //Is a pseudojet object with constituents of the jet
    fConstituents_tru.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks_tru; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex_tru = 0;
    for (auto part_tru: fJet_tru->GetParticleConstituents())
    {
        PseudoTracks_tru.reset(part_tru.Px(), part_tru.Py(), part_tru.Pz(), part_tru.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part_tru2 = part_tru.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks_tru.set_user_index(GetConstituentID(constituentIndex_tru, part_tru2, fJet_tru)); //leave this in for the same reason as above
        if (PseudoTracks_tru.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents_tru.push_back(PseudoTracks_tru);
        constituentIndex_tru++;
    }
    
   
        double jet_pt = fJet->Pt();
        
        double jet_pt_tru = fJet_tru->Pt();
        pt_tru->Fill(jet_pt_tru); //filling histogram with momentum of jets
        
        
         if(fpTcorr == 1)
         {   jet_pt = jet_pt/(0.85); //applying JES correction to jet pT spectra to study spectra shape dependence of ENC
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
             
             double diff_scaled = (jet_pt - jet_pt_tru)/jet_pt_tru;
             JES_scaled->Fill(jet_pt_tru,diff_scaled);
             
             R_matrix->Fill(jet_pt_tru,jet_pt); //Filling the response matrix}
         }
         
         else
         {
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
             R_matrix->Fill(jet_pt_tru,jet_pt); //Filling the response matrix
             
             double diff = (jet_pt - jet_pt_tru)/jet_pt_tru;
             JES->Fill(jet_pt_tru,diff);
         }
        
        //Det level
        //Initializing objects
        std::vector<Double_t> delta_Rvec;
        std::vector<Double_t> energy_pairs_vec; //the weighting vector with EE
        std::vector<Double_t> energy_pairs_tri; //the weighting vector with EEE
        std::vector<Double_t> R_dist;
        std::vector<Double_t> logR_dist;
        std::vector<Double_t> max_R_distvec;
        std::vector<Double_t> max_logR_distvec;
        
        //Truth level
        //Initializing objects
        std::vector<Double_t> delta_Rvec_part;
        std::vector<Double_t> energy_pairs_vec_part; //the weighting vector with EE
        std::vector<Double_t> energy_pairs_tri_part; //the weighting vector with EEE
        std::vector<Double_t> R_dist_part;
        std::vector<Double_t> logR_dist_part;
        std::vector<Double_t> max_R_distvec_part;
        std::vector<Double_t> max_logR_distvec_part;
        
        //Looping over the det jet
        //For jets with 2 constituents
        if(int(fConstituents.size()) == 2)
        {
            for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
            {
                //For 3 point correlator
                for(int s=0; s<int(fConstituents.size()) ; s++)
                {
                    if(s==j) continue; //if s=j this would be 0
                    double eee_jss_2 =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                    double deltaR_jss_2 = fConstituents[j].delta_R(fConstituents[s]);
                    double delta_logR_jss_2 = log(deltaR_jss_2);
                    
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_2_det = fConstituents[j].eta();
                        double s_eta_3pt_2_det = fConstituents[s].eta();
                        double del_js_eta_3pt_2_det = abs(j_eta_3pt_2_det-s_eta_3pt_2_det);
                        if (del_js_eta_3pt_2_det < 0.008) continue;
                        else
                        {
                            energy_pairs_tri.push_back(eee_jss_2);
                            max_R_distvec.push_back(deltaR_jss_2);
                            max_logR_distvec.push_back(delta_logR_jss_2);
                            
                            //            EEEC_hist->Fill(deltaR_jss_2,eee_jss_2);
                            E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                            E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                            E3C_pt_hist_log->Fill(delta_logR_jss_2,jet_pt,eee_jss_2);
                            
                            N3_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru);
                            E3C_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru, eee_jss_2);
                            E3C_det_pt_hist_log_3d->Fill(delta_logR_jss_2, jet_pt, jet_pt_tru, eee_jss_2);
                        }
                    }
                    else
                    {
                        energy_pairs_tri.push_back(eee_jss_2);
                        max_R_distvec.push_back(deltaR_jss_2);
                        max_logR_distvec.push_back(delta_logR_jss_2);
                        
                        //            EEEC_hist->Fill(deltaR_jss_2,eee_jss_2);
                        E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                        E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                        E3C_pt_hist_log->Fill(delta_logR_jss_2,jet_pt,eee_jss_2);
                        
                        N3_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru);
                        E3C_det_pt_hist_3d->Fill(deltaR_jss_2, jet_pt, jet_pt_tru, eee_jss_2);
                        E3C_det_pt_hist_log_3d->Fill(delta_logR_jss_2, jet_pt, jet_pt_tru, eee_jss_2);
                    }
                    
                }//close s loop for the 3 point correlator
                //For 2 point correlator
                for(int s=0; s<j ; s++)
                {
                    double delta_R_js_2 = fConstituents[j].delta_R(fConstituents[s]);
                    double log_delta_R_js_2 = log(delta_R_js_2);
                    double ee_js_2 = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt_2_det = fConstituents[j].eta();
                        double s_eta_2pt_2_det = fConstituents[s].eta();
                        double del_js_eta_2pt_2_det = abs(j_eta_2pt_2_det-s_eta_2pt_2_det);
                        if (del_js_eta_2pt_2_det < 0.008) continue;
                        else
                        {
                            //Filling the vectors
                            delta_Rvec.push_back(delta_R_js_2);
                            energy_pairs_vec.push_back(ee_js_2);
                            
                            EEC_hist->Fill(delta_R_js_2,ee_js_2);
                            EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                            EEC_pt_hist_log->Fill(log_delta_R_js_2,jet_pt,ee_js_2);
                            
                            N2_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru);
                            EEC_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru, ee_js_2);
                            EEC_det_pt_hist_log_3d->Fill(log_delta_R_js_2, jet_pt, jet_pt_tru, ee_js_2);
                        }
                    }
                    else
                    {
                        //Filling the vectors
                        delta_Rvec.push_back(delta_R_js_2);
                        energy_pairs_vec.push_back(ee_js_2);
                        
                        EEC_hist->Fill(delta_R_js_2,ee_js_2);
                        EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                        EEC_pt_hist_log->Fill(log_delta_R_js_2,jet_pt,ee_js_2);
                        
                        N2_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru);
                        EEC_det_pt_hist_3d->Fill(delta_R_js_2, jet_pt, jet_pt_tru, ee_js_2);
                        EEC_det_pt_hist_log_3d->Fill(log_delta_R_js_2, jet_pt, jet_pt_tru, ee_js_2);
                    }
                }//close s loop for eec
            }//close j loop
        }//close if loop
        
        //For jets with more than 2 constituents
        else
        {
            for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
            {
                
                for(int s=0; s<int(fConstituents.size()) ; s++)
                {
                    if(s==j) continue; //This ensures I don't get stuff like (000) for (jss)
                    
                    double eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                    double deltaR_jss = fConstituents[j].delta_R(fConstituents[s]);
                    double delta_logR_jss = log(deltaR_jss);
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_det = fConstituents[j].eta();
                        double s_eta_3pt_det = fConstituents[s].eta();
                        double del_js_eta_3pt_det = abs(j_eta_3pt_det-s_eta_3pt_det);
                        if (del_js_eta_3pt_det < 0.008) continue;
                        else
                        {
                            energy_pairs_tri.push_back(eee_jss);
                            max_R_distvec.push_back(deltaR_jss);
                            max_logR_distvec.push_back(delta_logR_jss);
                            
                            E3C_hist->Fill(deltaR_jss,eee_jss);
                            E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                            E3C_pt_hist_log->Fill(delta_logR_jss,jet_pt,eee_jss);
                            
                            N3_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru);
                            E3C_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru, eee_jss);
                            E3C_det_pt_hist_log_3d->Fill(delta_logR_jss, jet_pt, jet_pt_tru, eee_jss);
                        }
                    }
                    else
                    {
                        energy_pairs_tri.push_back(eee_jss);
                        max_R_distvec.push_back(deltaR_jss);
                        max_logR_distvec.push_back(delta_logR_jss);
                        
                        E3C_hist->Fill(deltaR_jss,eee_jss);
                        E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                        E3C_pt_hist_log->Fill(delta_logR_jss,jet_pt,eee_jss);
                        
                        N3_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru);
                        E3C_det_pt_hist_3d->Fill(deltaR_jss, jet_pt, jet_pt_tru, eee_jss);
                        E3C_det_pt_hist_log_3d->Fill(delta_logR_jss, jet_pt, jet_pt_tru, eee_jss);
                    }
                    //For 3 point correlator
                    for( int m=0; m!=j && m!=s; m++)
                    {
                        if(s>j) continue;
                        
                        double eee_jsm = ((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt())/(pow(jet_pt,3)));
                        double deltaR_js = fConstituents[j].delta_R(fConstituents[s]);
                        double delta_logR_js = log(deltaR_js);
                        
                        double deltaR_jm = fConstituents[j].delta_R(fConstituents[m]);
                        double delta_logR_jm = log(deltaR_jm);
                        
                        double deltaR_sm = fConstituents[s].delta_R(fConstituents[m]);
                        double delta_logR_sm = log(deltaR_sm);
                        
                          //Pair cut
                        if(fpaircut == 1)
                        {
                            double j_eta_3pt_det = fConstituents[j].eta();
                            double m_eta_3pt_det = fConstituents[m].eta();
                            double s_eta_3pt_det = fConstituents[s].eta();
                            double del_jm_eta_3pt_det = abs(j_eta_3pt_det-m_eta_3pt_det);
                            double del_sm_eta_3pt_det = abs(s_eta_3pt_det-m_eta_3pt_det);
                            double del_js_eta_3pt_det = abs(j_eta_3pt_det-s_eta_3pt_det);
                            if (del_jm_eta_3pt_det < 0.008 || del_sm_eta_3pt_det < 0.008 || del_js_eta_3pt_det < 0.008) continue;
                            else
                            {
                                energy_pairs_tri.push_back(eee_jsm);
                                //                jetE.push_back(jet_pt);
                                
                                R_dist.push_back(deltaR_js);
                                R_dist.push_back(deltaR_jm);
                                R_dist.push_back(deltaR_sm);
                                
                                logR_dist.push_back(delta_logR_js);
                                logR_dist.push_back(delta_logR_jm);
                                logR_dist.push_back(delta_logR_sm);
                                
                                int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                                
                                max_R_distvec.push_back(R_dist[max_R]);
                                max_logR_distvec.push_back(logR_dist[max_R]);
                                
                                E3C_hist->Fill(R_dist[max_R],eee_jsm);
                                E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
                                E3C_pt_hist_log->Fill(logR_dist[max_R],jet_pt,eee_jsm);
                                
                                N3_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru);
                                E3C_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru, eee_jsm);
                                E3C_det_pt_hist_log_3d->Fill(logR_dist[max_R], jet_pt, jet_pt_tru, eee_jsm);
                                
                                
                                R_dist.clear();
                                logR_dist.clear();
                            }
                        }
                        else
                        {
                            energy_pairs_tri.push_back(eee_jsm);
                            //                jetE.push_back(jet_pt);
                            
                            R_dist.push_back(deltaR_js);
                            R_dist.push_back(deltaR_jm);
                            R_dist.push_back(deltaR_sm);
                            
                            logR_dist.push_back(delta_logR_js);
                            logR_dist.push_back(delta_logR_jm);
                            logR_dist.push_back(delta_logR_sm);
                            
                            int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                            
                            max_R_distvec.push_back(R_dist[max_R]);
                            max_logR_distvec.push_back(logR_dist[max_R]);
                            
                            E3C_hist->Fill(R_dist[max_R],eee_jsm);
                            E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
                            E3C_pt_hist_log->Fill(logR_dist[max_R],jet_pt,eee_jsm);
                            
                            N3_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru);
                            E3C_det_pt_hist_3d->Fill(R_dist[max_R], jet_pt, jet_pt_tru, eee_jsm);
                            E3C_det_pt_hist_log_3d->Fill(logR_dist[max_R], jet_pt, jet_pt_tru, eee_jsm);
                            
                            
                            R_dist.clear();
                            logR_dist.clear();
                        }
                    }//close m loop
                }//close s loop for the 3 point correlator
                //For loop for EEC
                for(int s=0; s<j ; s++)
                {
                
                    double delta_R_js = fConstituents[j].delta_R(fConstituents[s]);
                    double log_delta_R_js = log(delta_R_js);
                    double ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt_det = fConstituents[j].eta();
                        double s_eta_2pt_det = fConstituents[s].eta();
                        double del_js_eta_2pt_det = abs(j_eta_2pt_det-s_eta_2pt_det);
                        if (del_js_eta_2pt_det < 0.008) continue;
                        else
                        {
                            //Filling the vectors
                            delta_Rvec.push_back(delta_R_js);
                            energy_pairs_vec.push_back(ee_js);
                            
                            EEC_hist->Fill(delta_R_js,ee_js);
                            EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                            EEC_pt_hist_log->Fill(log_delta_R_js,jet_pt,ee_js);
                            
                            N2_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru);
                            EEC_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru, ee_js);
                            EEC_det_pt_hist_log_3d->Fill(log_delta_R_js, jet_pt, jet_pt_tru, ee_js);
                        }
                    }
                    else
                    {
                        //Filling the vectors
                        delta_Rvec.push_back(delta_R_js);
                        energy_pairs_vec.push_back(ee_js);
                        
                        EEC_hist->Fill(delta_R_js,ee_js);
                        EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                        EEC_pt_hist_log->Fill(log_delta_R_js,jet_pt,ee_js);
                        
                        N2_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru);
                        EEC_det_pt_hist_3d->Fill(delta_R_js, jet_pt, jet_pt_tru, ee_js);
                        EEC_det_pt_hist_log_3d->Fill(log_delta_R_js, jet_pt, jet_pt_tru, ee_js);
                    }
                    
                }//close s loop for the 2 point correlator
            } //close j loop
        }//close else loop
        
        
        //Looping over truth level jet
        //For jets with 2 constituents
        if(int(fConstituents_tru.size()) == 2)
        {
            for(int j=0; j<int(fConstituents_tru.size()); j++)  //looping over constituents of the fConstituents_tru object
            {
                //For 3 point correlator
                for(int s=0; s<int(fConstituents_tru.size()); s++)
                {
                    
                    if(s==j) continue; //if s=j this would be 0
                    
                    double eee_jss_2_tru =((3*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[s].pt())/(pow(jet_pt_tru,3)));
                    double deltaR_jss_2_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                    double delta_logR_jss_2_tru = log(deltaR_jss_2_tru);
                    
                      //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_2_tru = fConstituents[j].eta();
                        double s_eta_3pt_2_tru = fConstituents[s].eta();
                        double del_js_eta_3pt_2_tru = abs(j_eta_3pt_2_tru-s_eta_3pt_2_tru);
                        if (del_js_eta_3pt_2_tru < 0.008) continue;
                        else
                        {
                            energy_pairs_tri_part.push_back(eee_jss_2_tru);
                            max_R_distvec_part.push_back(deltaR_jss_2_tru);
                            max_logR_distvec_part.push_back(delta_logR_jss_2_tru);
                            
                            E3C_tru_match_pt_tru->Fill(deltaR_jss_2_tru, jet_pt_tru, eee_jss_2_tru);
                            E3C_tru_match_pt_tru_log->Fill(delta_logR_jss_2_tru, jet_pt_tru, eee_jss_2_tru);
                            
                            N3_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt);
                            E3C_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt, eee_jss_2_tru);
                            E3C_tru_pt_hist_log_3d->Fill(delta_logR_jss_2_tru, jet_pt_tru, jet_pt, eee_jss_2_tru);
                        }
                    }
                    else
                    {
                        energy_pairs_tri_part.push_back(eee_jss_2_tru);
                        max_R_distvec_part.push_back(deltaR_jss_2_tru);
                        max_logR_distvec_part.push_back(delta_logR_jss_2_tru);
                        
                        E3C_tru_match_pt_tru->Fill(deltaR_jss_2_tru, jet_pt_tru, eee_jss_2_tru);
                        E3C_tru_match_pt_tru_log->Fill(delta_logR_jss_2_tru, jet_pt_tru, eee_jss_2_tru);
                        
                        N3_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt);
                        E3C_tru_pt_hist_3d->Fill(deltaR_jss_2_tru, jet_pt_tru, jet_pt, eee_jss_2_tru);
                        E3C_tru_pt_hist_log_3d->Fill(delta_logR_jss_2_tru, jet_pt_tru, jet_pt, eee_jss_2_tru);
                    }
                }//close s loop for the 3 point correlator
                //For 2 point correlator
                for(int s=0; s<j ; s++)
                {
                    
                    double delta_R_js_2_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                    double log_delta_R_js_2_tru = log(delta_R_js_2_tru);
                    double ee_js_2_tru = (2*fConstituents_tru[j].pt()*fConstituents_tru[s].pt())/(pow((jet_pt_tru),2));
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt_2_tru = fConstituents[j].eta();
                        double s_eta_2pt_2_tru = fConstituents[s].eta();
                        double del_js_eta_2pt_2_tru = abs(j_eta_2pt_2_tru-s_eta_2pt_2_tru);
                        if (del_js_eta_2pt_2_tru < 0.008) continue;
                        else
                        {
                            //Filling the vectors
                            delta_Rvec_part.push_back(delta_R_js_2_tru);
                            energy_pairs_vec_part.push_back(ee_js_2_tru);
                            
                            EEC_tru_match_pt_tru->Fill(delta_R_js_2_tru, jet_pt_tru, ee_js_2_tru);
                            EEC_tru_match_pt_tru_log->Fill(log_delta_R_js_2_tru, jet_pt_tru, ee_js_2_tru);
                            
                            N2_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt);
                            EEC_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt, ee_js_2_tru);
                            EEC_tru_pt_hist_log_3d->Fill(log_delta_R_js_2_tru, jet_pt_tru, jet_pt, ee_js_2_tru);
                        }
                    }
                    else
                    {
                        //Filling the vectors
                        delta_Rvec_part.push_back(delta_R_js_2_tru);
                        energy_pairs_vec_part.push_back(ee_js_2_tru);
                        
                        EEC_tru_match_pt_tru->Fill(delta_R_js_2_tru, jet_pt_tru, ee_js_2_tru);
                        EEC_tru_match_pt_tru_log->Fill(log_delta_R_js_2_tru, jet_pt_tru, ee_js_2_tru);
                        
                        N2_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt);
                        EEC_tru_pt_hist_3d->Fill(delta_R_js_2_tru, jet_pt_tru, jet_pt, ee_js_2_tru);
                        EEC_tru_pt_hist_log_3d->Fill(log_delta_R_js_2_tru, jet_pt_tru, jet_pt, ee_js_2_tru);
                    }
                }//close s loop for eec
            }//close j loop
        }//close if loop
        
        //For jets with more than 2 constituents
        else
        {
            for(int j=0; j<int(fConstituents_tru.size()); j++)  //looping over constituents of the fConstituents_tru object
            {
                for(int s=0; s<int(fConstituents_tru.size()) ; s++)
                {
                    if(s==j) continue; //This ensures I don't get stuff like (000) for (jss)
                    
                    double eee_jss_tru =((3*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[s].pt())/(pow(jet_pt_tru,3)));
                    double deltaR_jss_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                    double delta_logR_jss_tru = log(deltaR_jss_tru);
                    
                     //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_tru = fConstituents[j].eta();
                        double s_eta_3pt_tru = fConstituents[s].eta();
                        double del_js_eta_3pt_tru = abs(j_eta_3pt_tru-s_eta_3pt_tru);
                        if (del_js_eta_3pt_tru < 0.008) continue;
                        else
                        {
                            energy_pairs_tri_part.push_back(eee_jss_tru);
                            max_R_distvec_part.push_back(deltaR_jss_tru);
                            max_logR_distvec_part.push_back(delta_logR_jss_tru);
                            
                            E3C_tru_match_pt_tru->Fill(deltaR_jss_tru,jet_pt_tru, eee_jss_tru);
                            E3C_tru_match_pt_tru_log->Fill(delta_logR_jss_tru,jet_pt_tru ,eee_jss_tru);
                            
                            N3_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt);
                            E3C_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt, eee_jss_tru);
                            E3C_tru_pt_hist_log_3d->Fill(delta_logR_jss_tru, jet_pt_tru, jet_pt, eee_jss_tru);
                        }
                    }
                    else
                    {
                        energy_pairs_tri_part.push_back(eee_jss_tru);
                        max_R_distvec_part.push_back(deltaR_jss_tru);
                        max_logR_distvec_part.push_back(delta_logR_jss_tru);
                        
                        E3C_tru_match_pt_tru->Fill(deltaR_jss_tru,jet_pt_tru, eee_jss_tru);
                        E3C_tru_match_pt_tru_log->Fill(delta_logR_jss_tru,jet_pt_tru ,eee_jss_tru);
                        
                        N3_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt);
                        E3C_tru_pt_hist_3d->Fill(deltaR_jss_tru, jet_pt_tru, jet_pt, eee_jss_tru);
                        E3C_tru_pt_hist_log_3d->Fill(delta_logR_jss_tru, jet_pt_tru, jet_pt, eee_jss_tru);
                    }
                    //For 3 point correlator
                    for( int m=0; m!=j && m!=s; m++)
                    {
                        if(s>j) continue;
                        
                        double eee_jsm_tru = ((6*fConstituents_tru[j].pt()*fConstituents_tru[s].pt()*fConstituents_tru[m].pt())/(pow(jet_pt_tru,3)));
                        double deltaR_js_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                        double delta_logR_js_tru = log(deltaR_js_tru);
                        
                        double deltaR_jm_tru = fConstituents_tru[j].delta_R(fConstituents_tru[m]);
                        double delta_logR_jm_tru = log(deltaR_jm_tru);
                        
                        double deltaR_sm_tru = fConstituents_tru[s].delta_R(fConstituents_tru[m]);
                        double delta_logR_sm_tru = log(deltaR_sm_tru);
                        
                        energy_pairs_tri_part.push_back(eee_jsm_tru);
                        //                jetE.push_back(jet_pt_tru);
                        
                        R_dist_part.push_back(deltaR_js_tru);
                        R_dist_part.push_back(deltaR_jm_tru);
                        R_dist_part.push_back(deltaR_sm_tru);
                        
                        logR_dist_part.push_back(delta_logR_js_tru);
                        logR_dist_part.push_back(delta_logR_jm_tru);
                        logR_dist_part.push_back(delta_logR_sm_tru);
                        
                        int max_R_tru = distance(R_dist_part.begin(), max_element(R_dist_part.begin(), R_dist_part.end()));//pick the longest side to compute the correlators with
                        
                        //Pair cut
                        if(fpaircut == 1)
                        {
                            double j_eta_3pt_tru = fConstituents[j].eta();
                            double m_eta_3pt_tru = fConstituents[m].eta();
                            double s_eta_3pt_tru = fConstituents[s].eta();
                            double del_jm_eta_3pt_tru = abs(j_eta_3pt_tru-m_eta_3pt_tru);
                            double del_sm_eta_3pt_tru = abs(s_eta_3pt_tru-m_eta_3pt_tru);
                            double del_js_eta_3pt_tru = abs(j_eta_3pt_tru-s_eta_3pt_tru);
                            if (del_jm_eta_3pt_tru < 0.008 || del_sm_eta_3pt_tru < 0.008 || del_js_eta_3pt_tru < 0.008) continue;
                            else
                            {
                                max_R_distvec_part.push_back(R_dist_part[max_R_tru]);
                                max_logR_distvec_part.push_back(logR_dist_part[max_R_tru]);
                                
                                E3C_tru_match_pt_tru->Fill(R_dist_part[max_R_tru], jet_pt_tru, eee_jsm_tru);
                                E3C_tru_match_pt_tru_log->Fill(logR_dist_part[max_R_tru], jet_pt_tru, eee_jsm_tru);
                                
                                N3_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt);
                                E3C_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt, eee_jsm_tru);
                                E3C_tru_pt_hist_log_3d->Fill(logR_dist_part[max_R_tru], jet_pt_tru, jet_pt, eee_jsm_tru);
                                
                                R_dist_part.clear();
                                logR_dist_part.clear();
                            }
                        }
                        else
                        {
                            max_R_distvec_part.push_back(R_dist_part[max_R_tru]);
                            max_logR_distvec_part.push_back(logR_dist_part[max_R_tru]);
                            
                            E3C_tru_match_pt_tru->Fill(R_dist_part[max_R_tru], jet_pt_tru, eee_jsm_tru);
                            E3C_tru_match_pt_tru_log->Fill(logR_dist_part[max_R_tru], jet_pt_tru, eee_jsm_tru);
                            
                            N3_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt);
                            E3C_tru_pt_hist_3d->Fill(R_dist_part[max_R_tru], jet_pt_tru, jet_pt, eee_jsm_tru);
                            E3C_tru_pt_hist_log_3d->Fill(logR_dist_part[max_R_tru], jet_pt_tru, jet_pt, eee_jsm_tru);
                            
                            R_dist_part.clear();
                            logR_dist_part.clear();
                        }
                    }//close m loop
                }//close s loop for the 3 point correlator
                //For loop for EEC
                for(int s=0; s<j ; s++)
                {
                    double delta_R_js_tru = fConstituents_tru[j].delta_R(fConstituents_tru[s]);
                    double log_delta_R_js_tru = log(delta_R_js_tru);
                    double ee_js_tru = (2*fConstituents_tru[j].pt()*fConstituents_tru[s].pt())/(pow((jet_pt_tru),2));
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt_tru = fConstituents[j].eta();
                        double s_eta_2pt_tru = fConstituents[s].eta();
                        double del_js_eta_2pt_tru = abs(j_eta_2pt_tru-s_eta_2pt_tru);
                        if (del_js_eta_2pt_tru < 0.008) continue;
                        else{
                            //Filling the vectors
                            delta_Rvec_part.push_back(delta_R_js_tru);
                            energy_pairs_vec_part.push_back(ee_js_tru);
                            
                            //            EEC_hist->Fill(delta_R_js,ee_js);
                            EEC_tru_match_pt_tru->Fill(delta_R_js_tru,jet_pt_tru,ee_js_tru);
                            EEC_tru_match_pt_tru_log->Fill(log_delta_R_js_tru,jet_pt_tru,ee_js_tru);
                            
                            N2_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt);
                            EEC_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt, ee_js_tru);
                            EEC_tru_pt_hist_log_3d->Fill(log_delta_R_js_tru, jet_pt_tru, jet_pt, ee_js_tru);
                        }
                    }
                    else
                    {
                        //Filling the vectors
                        delta_Rvec_part.push_back(delta_R_js_tru);
                        energy_pairs_vec_part.push_back(ee_js_tru);
                        
                        //            EEC_hist->Fill(delta_R_js,ee_js);
                        EEC_tru_match_pt_tru->Fill(delta_R_js_tru,jet_pt_tru,ee_js_tru);
                        EEC_tru_match_pt_tru_log->Fill(log_delta_R_js_tru,jet_pt_tru,ee_js_tru);
                        
                        N2_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt);
                        EEC_tru_pt_hist_3d->Fill(delta_R_js_tru, jet_pt_tru, jet_pt, ee_js_tru);
                        EEC_tru_pt_hist_log_3d->Fill(log_delta_R_js_tru, jet_pt_tru, jet_pt, ee_js_tru);
                    }
                }//close s loop for the 2 point correlator
            } //close j loop
        }//close else loop

    
    //catch (fastjet::Error)
    // {
    //    AliError(" [w] FJ Exception caught.");
    //    // return -1;
    // } //end error message
    return;
}

    
//EEC computation-------------------------------------------------------
void AliAnalysisTaskJetsEEC::ComputeEEC(AliEmcalJet *fJet, AliJetContainer *fJetCont)
    //(jet, jet container, vector of jets, vector of constituents within each jet?)
{
    //General EEC computation: Need to loop over jets, then loop within a single jet to compute EECs.
    //NOTE: Already in the jet loop defined in the Fill Histogram function. This puts you in the jet loop. The event loop stuff is being taken care of by AliAnalysis manager
    //fjet is my single jet. AliEmCalJet is pointing to the fJet object.
    
    std::vector<fastjet::PseudoJet> fConstituents; //Is a pseudojet object with constituents of the jet
    fConstituents.clear();
    //This snippet of code is getting particles within a single jet (fjet) and turning them into pseudojet objects so that fastjet capabilities can be used
    fastjet::PseudoJet PseudoTracks; //Creating a pseudojet object called PseduoTracks
    unsigned int constituentIndex = 0;
    //The line below gets constituent particles within fjet. C++ syntax[ for (auto elem : container)    // capture elements by value ]
    for (auto part: fJet->GetParticleConstituents())
    {
        PseudoTracks.reset(part.Px(), part.Py(), part.Pz(), part.E()); //part is the constituent at that point in the loop, part keeps getting redefined in each step.
        const AliVParticle* part2 = part.GetParticle(); //"hack", leave this in , to get the index of the jet from AliPhysics
        PseudoTracks.set_user_index(GetConstituentID(constituentIndex, part2, fJet)); //leave this in for the same reason as above
        if (PseudoTracks.pt() < fMinENCtrackPt) continue; //remove tracks below cut for ENCs
        fConstituents.push_back(PseudoTracks);
        constituentIndex++;
    }
    
        //Initializing objects for det level
        std::vector<Double_t> delta_Rvec;
        std::vector<Double_t> energy_pairs_vec; //the weighting vector with EE
        std::vector<Double_t> energy_pairs_tri; //the weighting vector with EEE
        std::vector<Double_t> R_dist;
        std::vector<Double_t> logR_dist;
        std::vector<Double_t> max_R_distvec;
        std::vector<Double_t> max_logR_distvec;
        
        //Looping over the jet
        double jet_pt = fJet->Pt();
        
        if(fpTcorr == 1)
         {   jet_pt = jet_pt/(0.85); //applying JES correction to jet pT spectra to study spectra shape dependence of ENC
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
             
         }
         
         else
         {
             jet_pt_hist->Fill(jet_pt); //filling histogram with momentum of jets
         
         }
        
        //For jets with 2 constituents
        if(int(fConstituents.size()) == 2)
        {
            for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
            {
                //For 3 point correlator
                for(int s=0; s<j ; s++)
                {
                    if(s==j) continue; //this ensure i dont get 000
                    
                    double eee_jss_2 =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                    double deltaR_jss_2 = fConstituents[j].delta_R(fConstituents[s]);
                    double delta_logR_jss_2 = log(deltaR_jss_2);
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt_2 = fConstituents[j].eta();
                        double s_eta_3pt_2 = fConstituents[s].eta();
                        double del_js_eta_3pt_2 = abs(j_eta_3pt_2-s_eta_3pt_2);
                        if (del_js_eta_3pt_2 < 0.008) continue;
                        else
                        {
                            energy_pairs_tri.push_back(eee_jss_2);
                            max_R_distvec.push_back(deltaR_jss_2);
                            max_logR_distvec.push_back(delta_logR_jss_2);
                            
                            E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                            E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                            E3C_pt_hist_log->Fill(delta_logR_jss_2,jet_pt,eee_jss_2);
                        }
                    }
                    else
                    {
                        energy_pairs_tri.push_back(eee_jss_2);
                        max_R_distvec.push_back(deltaR_jss_2);
                        max_logR_distvec.push_back(delta_logR_jss_2);
                        
                        E3C_hist->Fill(deltaR_jss_2,eee_jss_2);
                        E3C_pt_hist->Fill(deltaR_jss_2,jet_pt,eee_jss_2);
                        E3C_pt_hist_log->Fill(delta_logR_jss_2,jet_pt,eee_jss_2);
                    }
                }//close s loop for the 3 point correlator
                //For 2 point correlator
                for(int s=0; s<j ; s++)
                {
                    
                    double delta_R_js_2 = fConstituents[j].delta_R(fConstituents[s]);
                    double log_delta_R_js_2 = log(delta_R_js_2);
                    double ee_js_2 = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt_2 = fConstituents[j].eta();
                        double s_eta_2pt_2 = fConstituents[s].eta();
                        double del_js_eta_2pt_2 = abs(j_eta_2pt_2-s_eta_2pt_2);
                        if (del_js_eta_2pt_2 < 0.008) continue;
                        else
                        {
                            //Filling the vectors
                            delta_Rvec.push_back(delta_R_js_2);
                            energy_pairs_vec.push_back(ee_js_2);
                            
                            EEC_hist->Fill(delta_R_js_2,ee_js_2);
                            EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                            EEC_pt_hist_log->Fill(log_delta_R_js_2,jet_pt,ee_js_2);
                        }
                    }
                    else
                    {
                        //Filling the vectors
                        delta_Rvec.push_back(delta_R_js_2);
                        energy_pairs_vec.push_back(ee_js_2);
                        
                        EEC_hist->Fill(delta_R_js_2,ee_js_2);
                        EEC_pt_hist->Fill(delta_R_js_2,jet_pt,ee_js_2);
                        EEC_pt_hist_log->Fill(log_delta_R_js_2,jet_pt,ee_js_2);
                    }
                    
                }//close s loop for eec
            }//close j loop
        }//close if loop
        
        //For jets with more than 2 constituents
        else
        {
            for(int j=0; j<int(fConstituents.size()); j++)  //looping over constituents of the fConstituents object
            {
                
                for(int s=0; s<int(fConstituents.size()) ; s++)
                {
                    if(s==j) continue; //This ensures I don't get stuff like (000) for (jss)
                    
                    double eee_jss =((3*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[s].pt())/(pow(jet_pt,3)));
                    double deltaR_jss = fConstituents[j].delta_R(fConstituents[s]);
                    double delta_logR_jss = log(deltaR_jss);
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_3pt = fConstituents[j].eta();
                        double s_eta_3pt = fConstituents[s].eta();
                        double del_js_eta_3pt = abs(j_eta_3pt-s_eta_3pt);
                        if (del_js_eta_3pt < 0.008) continue;
                        else
                        {
                            energy_pairs_tri.push_back(eee_jss);
                            max_R_distvec.push_back(deltaR_jss);
                            max_logR_distvec.push_back(delta_logR_jss);
                            
                            E3C_hist->Fill(deltaR_jss,eee_jss);
                            E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                            E3C_pt_hist_log->Fill(delta_logR_jss,jet_pt,eee_jss);
                        }
                    }
                    else
                    {
                        energy_pairs_tri.push_back(eee_jss);
                        max_R_distvec.push_back(deltaR_jss);
                        max_logR_distvec.push_back(delta_logR_jss);
                        
                        E3C_hist->Fill(deltaR_jss,eee_jss);
                        E3C_pt_hist->Fill(deltaR_jss,jet_pt,eee_jss);
                        E3C_pt_hist_log->Fill(delta_logR_jss,jet_pt,eee_jss);
                    }
                    //For 3 point correlator
                    for( int m=0; m!=j && m!=s; m++)
                    {
                        if(s>j) continue;
                        
                            double eee_jsm = ((6*fConstituents[j].pt()*fConstituents[s].pt()*fConstituents[m].pt())/(pow(jet_pt,3)));
                            double deltaR_js = fConstituents[j].delta_R(fConstituents[s]);
                            double delta_logR_js = log(deltaR_js);
                            
                            double deltaR_jm = fConstituents[j].delta_R(fConstituents[m]);
                            double delta_logR_jm = log(deltaR_jm);
                            
                            double deltaR_sm = fConstituents[s].delta_R(fConstituents[m]);
                            double delta_logR_sm = log(deltaR_sm);
                            
                            energy_pairs_tri.push_back(eee_jsm);
                            //                jetE.push_back(jet_pt);
                            
                            R_dist.push_back(deltaR_js);
                            R_dist.push_back(deltaR_jm);
                            R_dist.push_back(deltaR_sm);
                            
                            logR_dist.push_back(delta_logR_js);
                            logR_dist.push_back(delta_logR_jm);
                            logR_dist.push_back(delta_logR_sm);
                            
                            int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                            
                            if(fpaircut==1)
                            {
                                //Pair cut
                                double j_eta_3pt = fConstituents[j].eta();
                                double m_eta_3pt = fConstituents[m].eta();
                                double s_eta_3pt = fConstituents[s].eta();
                                double del_jm_eta_3pt = abs(j_eta_3pt-m_eta_3pt);
                                double del_sm_eta_3pt = abs(s_eta_3pt-m_eta_3pt);
                                double del_js_eta_3pt = abs(j_eta_3pt-s_eta_3pt);
                                if (del_jm_eta_3pt < 0.008 || del_sm_eta_3pt < 0.008 || del_js_eta_3pt < 0.008) continue;
                                else
                                {
                                    max_R_distvec.push_back(R_dist[max_R]);
                                    max_logR_distvec.push_back(logR_dist[max_R]);
                                    
                                    E3C_hist->Fill(R_dist[max_R],eee_jsm);
                                    E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
                                    E3C_pt_hist_log->Fill(logR_dist[max_R],jet_pt,eee_jsm);
                                    
                                    R_dist.clear();
                                    logR_dist.clear();
                                }
                            }
                            else
                            {
                                max_R_distvec.push_back(R_dist[max_R]);
                                max_logR_distvec.push_back(logR_dist[max_R]);
                                
                                E3C_hist->Fill(R_dist[max_R],eee_jsm);
                                E3C_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm);
                                E3C_pt_hist_log->Fill(logR_dist[max_R],jet_pt,eee_jsm);
                                
                                R_dist.clear();
                                logR_dist.clear();
                            }

                    }//close m loop
                }//close s loop for the 3 point correlator
                //For loop for EEC
                for(int s=0; s<j ; s++)
                {
                    
                    double delta_R_js = fConstituents[j].delta_R(fConstituents[s]);
                    double log_delta_R_js = log(delta_R_js);
                    double ee_js = (2*fConstituents[j].pt()*fConstituents[s].pt())/(pow((jet_pt),2));
                    
                    //Pair cut
                    if(fpaircut == 1)
                    {
                        double j_eta_2pt = fConstituents[j].eta();
                        double s_eta_2pt = fConstituents[s].eta();
                        double del_js_eta_2pt = abs(j_eta_2pt-s_eta_2pt);
                        if (del_js_eta_2pt < 0.008) continue;
                        else
                        {
                            delta_Rvec.push_back(delta_R_js);
                            energy_pairs_vec.push_back(ee_js);
                            
                            EEC_hist->Fill(delta_R_js,ee_js);
                            EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                            EEC_pt_hist_log->Fill(log_delta_R_js,jet_pt,ee_js);
                        }
                    }
                    else
                    {
                        //Filling the vectors
                        delta_Rvec.push_back(delta_R_js);
                        energy_pairs_vec.push_back(ee_js);
                        
                        EEC_hist->Fill(delta_R_js,ee_js);
                        EEC_pt_hist->Fill(delta_R_js,jet_pt,ee_js);
                        EEC_pt_hist_log->Fill(log_delta_R_js,jet_pt,ee_js);
                    }
                }//close s loop for the 2 point correlator
            } //close j loop
        }//close else loop

    
    //catch (fastjet::Error)
    // {
    //    AliError(" [w] FJ Exception caught.");
    //    // return -1;
    // } //end error message
    return;
}


//_________________________________________________________________
void AliAnalysisTaskJetsEEC::RunChanged(Int_t newrun)
{
  if(fStoreTrig)
  {
    auto downscalehandler = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    if(downscalehandler->GetCurrentRun() != newrun)
    {
      downscalehandler->SetRun(newrun);
    }
  }
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetsEEC::RetrieveEventObjects()
{
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  return kTRUE;
}


//_______________________________________________________________________
void AliAnalysisTaskJetsEEC::Terminate(Option_t *)
{
  // Called once at the end of the analysis.
  // fTreeObservableTagging = dynamic_cast<TTree*>(GetOutputData(1));
  // if (!fTreeObservableTagging){
  //   Printf("ERROR: fTreeObservableTagging not available");
  //   return;
  // }
}

