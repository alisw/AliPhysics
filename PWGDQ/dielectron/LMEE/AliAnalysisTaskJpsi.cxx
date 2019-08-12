/*************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//    Class to produce jpsi cocktail                  (analysis task)   //
//                                        (description in .h file)       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// c++ includes
#include <iostream>
// ROOT includes
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include <TMath.h>
#include "TParticle.h"
#include "TLorentzVector.h"
// AliRoot includes
#include "AliAnalysisTaskJpsi.h" // this task
#include "AliAODMCParticle.h"


ClassImp(AliAnalysisTaskJpsi)

//______________________________| Default Constructor
AliAnalysisTaskJpsi::AliAnalysisTaskJpsi():
AliCocktailSmearing(),
  fTree(0x0),  
  fMotherOld(0x0),  
  fDaughter1Old(0x0),
  fDaughter2Old(0x0),
  fMother(0x0),  
  fDaughter1(0x0),
  fDaughter2(0x0),
  fOldTree(kTRUE),
  fYmin(-1.),
  fYmax(1.),
  fScaleFunction(0x0),
  fScaleFunctionLow(0x0),
  fScaleFunctionHigh(0x0),
  fJPsiPt(0x0), 
  fMee_Ptee_Jpsi_rec_die(0x0), 
  fMee_Ptee_Jpsi_gen_die(0x0),
  fMee_Ptee_Jpsi_rec_die_low(0x0),
  fMee_Ptee_Jpsi_gen_die_low(0x0),
  fMee_Ptee_Jpsi_rec_die_high(0x0),
  fMee_Ptee_Jpsi_gen_die_high(0x0),
  fMee_Ptee_Jpsi_rec_rad(0x0),
  fMee_Ptee_Jpsi_gen_rad(0x0),
  fMee_Ptee_Jpsi_rec_rad_low(0x0),
  fMee_Ptee_Jpsi_gen_rad_low(0x0),
  fMee_Ptee_Jpsi_rec_rad_high(0x0),
  fMee_Ptee_Jpsi_gen_rad_high(0x0),
  fOutputList(0x0),
  fOutputList_Low(0x0),
  fOutputList_High(0x0)
{
  //
  // Default constructor
  //

}
//______________________________ Destructor _________________
AliAnalysisTaskJpsi::~AliAnalysisTaskJpsi()
{
  // Destructor
  Info("~AliAnalysisTaskJpsi","Calling Destructor");
  if (fOutputList) delete fOutputList;
  if (fOutputList_Low) delete fOutputList_Low;
  if (fOutputList_High) delete fOutputList_High;
}

//______________________________ User Output ______________
void AliAnalysisTaskJpsi::Init()
{
  //
  // Create Histos
  //

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fOutputList_Low = new TList();
  fOutputList_Low->SetOwner(kTRUE);
  
  fOutputList_High = new TList();
  fOutputList_High->SetOwner(kTRUE);

  printf("\n\n========================================\n  Configuration of task: \n========================================\n\n");
  printf("  rapidity range for Jpsi normalization:         %f < y < %f \n", fYmin,fYmax);
  printf("  pt range single electrons:         %f < p_{T} < %f \n", fPtCutRange[0],fPtCutRange[1]);
  printf("  eta range single electrons:         %f < eta < %f \n", fEtaCutRange[0],fEtaCutRange[1]);
  printf("\n\n========================================\n  Configuration of smearing: \n========================================\n\n");
  AliCocktailSmearing::Print();
  printf("\n\n========================================\n  Configuration of param: \n========================================\n\n");
  if(!fScaleFunction) printf("Issue: no central parametrization given !!!!!\n");
  else printf("Central param is: %s\n",fScaleFunction->GetExpFormula("p").Data());
  if(!fScaleFunctionLow) printf("Issue: no low parametrization given !!!!!\n");
  else printf("Central param is: %s\n",fScaleFunctionLow->GetExpFormula("p").Data());
  if(!fScaleFunctionHigh) printf("Issue: no high parametrization given !!!!!\n");
  else printf("Central param is: %s\n",fScaleFunctionHigh->GetExpFormula("p").Data());
  
  CreateHistos();
  
}
//______________________________ User Exec ____________________________
void AliAnalysisTaskJpsi::ConnectTree(TTree *tree)
{

  //
  // Connect the tree
  //

  fTree = tree;

  if(!fOldTree) {
    fTree->SetBranchAddress("mother",&fMother);
    fTree->SetBranchAddress("daughter1",&fDaughter1);
    fTree->SetBranchAddress("daughter2",&fDaughter2);
  }
  else {
    fTree->SetBranchAddress("mother",&fMotherOld);
    fTree->SetBranchAddress("daughter1",&fDaughter1Old);
    fTree->SetBranchAddress("daughter2",&fDaughter2Old);
  }
}
//_______________________________________________________________
void AliAnalysisTaskJpsi::DetermineWeights()
{

  //
  // Compute MC pt spectrum of J/Psi for the given tree
  //
  
  if(!fTree) return;
  
  for (Int_t i = 0; i < fTree->GetEntries(); i++) {

    if(i % 10000 == 0) Printf("%12d of %12d",i,(Int_t) fTree->GetEntries());
    fTree->GetEntry(i);

    if(!fOldTree){
      if(!(fMother&&fDaughter1&&fDaughter2)){ Printf("particle corrupted"); continue; }
      if(fMother->Y() < fYmin || fMother->Y() > fYmax) continue;
      fJPsiPt->Fill(fMother->Pt());
    }
    else {
      if(!(fMotherOld&&fDaughter1Old&&fDaughter2Old)){ Printf("particle corrupted"); continue; }
      if(fMotherOld->Y() < fYmin || fMotherOld->Y() > fYmax) continue;
      fJPsiPt->Fill(fMotherOld->Pt());
    }
      
  }
  fJPsiPt->Scale(1./(fYmax-fYmin),"width"); // dN/dydpt

}
//__________________________________________________________________
void AliAnalysisTaskJpsi::Fill()
{
  //
  // Fill the mee,ptee histo
  //

  Double_t BR_jpsiee = 0.05971; // Branching ratio jpsi->ee
  Double_t BR_jpsirad =  0.0088; // Branching ratio jpsi->eegamma
  Double_t BR_jpsisum = BR_jpsiee + BR_jpsirad; // Total Branching ratio

  if(!fTree) return;
  if(!fScaleFunction) return;
  if(!fScaleFunctionLow) return;
  if(!fScaleFunctionHigh) return;

  //
  // Loop over the particles
  //
  
  for (Int_t i = 0; i < fTree->GetEntries(); i++) {

    if(i % 10000 == 0) printf("%12d of %12d\n",i,(Int_t) fTree->GetEntries());
    fTree->GetEntry(i);
    Double_t mpt = -1.;
    if(!fOldTree) {
      if(!(fMother&&fDaughter1&&fDaughter2)){ printf("particle corrupted\n"); continue; }
      if(fMother->Y() < fYmin || fMother->Y() > fYmax) continue;
      mpt = fMother->Pt();
    }
    else {
      if(!(fMotherOld&&fDaughter1Old&&fDaughter2Old)){ printf("particle corrupted\n"); continue; }
      if(fMotherOld->Y() < fYmin || fMotherOld->Y() > fYmax) continue;
      mpt = fMotherOld->Pt();
    }
    
    
    //
    // Determine the weight to fit the parametrization assuming forced decay into e+e- and e+e-gamma keeping the right propertion between the two
    //
    Int_t binpt = fJPsiPt->GetXaxis()->FindBin(mpt);
    Double_t lowedge = fJPsiPt->GetXaxis()->GetBinLowEdge(binpt);
    Double_t upedge = fJPsiPt->GetXaxis()->GetBinUpEdge(binpt);
    if(mpt < lowedge) binpt = binpt-1;
    if(mpt > upedge) binpt = binpt+1;
    Double_t bc = fJPsiPt->GetBinContent(binpt)/BR_jpsisum;
    Double_t sf(0.),sfLow(0.),sfHigh(0.);
    if(bc > 0.){ 
      Double_t bincenter = fJPsiPt->GetXaxis()->GetBinCenter(binpt);
      sf     = fScaleFunction     ->Eval(bincenter) / bc;
      sfLow  = fScaleFunctionLow  ->Eval(bincenter) / bc;
      sfHigh = fScaleFunctionHigh ->Eval(bincenter) / bc;
    }

    //
    // Variables without bremsstrahlung+momentum resolution
    //
    TLorentzVector ppo_1,ppo_2;
    if(!fOldTree){
      fDaughter1->Momentum(ppo_1);
      fDaughter2->Momentum(ppo_2);
    } else {
      fDaughter1Old->Momentum(ppo_1);
      fDaughter2Old->Momentum(ppo_2);
    }    
    Double_t mass_o     = (ppo_1 + ppo_2).M();
    Double_t pt_pair_o  = (ppo_1 + ppo_2).Pt();
    //Double_t opAngle_o  = ppo_1.Angle(ppo_2.Vect());

    // Charge
    int pid_1 = -999;
    int pid_2 = -999;
    if(!fOldTree){
      pid_1 = fDaughter1->PdgCode();
      pid_2 = fDaughter2->PdgCode();
    }
    else {
      pid_1 = fDaughter1Old->GetPdgCode();
      pid_2 = fDaughter2Old->GetPdgCode();
    }
    Short_t ch_1 = 1;
    if(pid_1>0) ch_1 = -1; // 11 is electron, -11 is positron
    Short_t ch_2 = 1;
    if(pid_2>0) ch_2 = -1; // 11 is electron, -11 is positron

    //
    // Apply smearing
    //
    TLorentzVector pp_1,pp_2;
    if(!GetSmearingOton()) {
      printf("Old smearing not supported !!\n");
      pp_1 = ppo_1;
      pp_2 = ppo_2;
    } else {
      pp_1 = ApplySmearingOton(ppo_1,ch_1);
      pp_2 = ApplySmearingOton(ppo_2,ch_2);
    }

    // Smeared values
    Double_t eta2  = pp_2.Eta();
    Double_t pt2   = pp_2.Pt();    
    Double_t eta1  = pp_1.Eta();
    Double_t pt1   = pp_1.Pt();
    Double_t mass     = (pp_1 + pp_2).M();
    Double_t pt_pair  = (pp_1 + pp_2).Pt();
    //Double_t opAngle  = pp_1.Angle(pp_2.Vect());
   
        
    //
    // Fiducial cuts
    //
    if(eta1 < fEtaCutRange[0] || eta1 > fEtaCutRange[1]) continue;
    if(pt1  < fPtCutRange[0]  || pt1  > fPtCutRange[1])  continue;
    if(eta2 < fEtaCutRange[0] || eta2 > fEtaCutRange[1]) continue;
    if(pt2  < fPtCutRange[0]  || pt2  > fPtCutRange[1])  continue;

    //
    // Type of decay: radiative or not
    //
    Int_t typeIdx(-1);
    if(mass_o > 3.09) typeIdx = 0;
    else                   typeIdx = 1;

    //
    // Fill Histos
    //
    if(typeIdx==0) {
      fMee_Ptee_Jpsi_rec_die->Fill(mass,pt_pair,sf);
      fMee_Ptee_Jpsi_gen_die->Fill(mass_o,pt_pair_o,sf);
      fMee_Ptee_Jpsi_rec_die_low->Fill(mass,pt_pair,sfLow);
      fMee_Ptee_Jpsi_gen_die_low->Fill(mass_o,pt_pair_o,sfLow);
      fMee_Ptee_Jpsi_rec_die_high->Fill(mass,pt_pair,sfHigh);
      fMee_Ptee_Jpsi_gen_die_high->Fill(mass_o,pt_pair_o,sfHigh);
    } else {
      fMee_Ptee_Jpsi_rec_rad->Fill(mass,pt_pair,sf);
      fMee_Ptee_Jpsi_gen_rad->Fill(mass_o,pt_pair_o,sf);
      fMee_Ptee_Jpsi_rec_rad_low->Fill(mass,pt_pair,sfLow);
      fMee_Ptee_Jpsi_gen_rad_low->Fill(mass_o,pt_pair_o,sfLow);
      fMee_Ptee_Jpsi_rec_rad_high->Fill(mass,pt_pair,sfHigh);
      fMee_Ptee_Jpsi_gen_rad_high->Fill(mass_o,pt_pair_o,sfHigh);
    }    
    
  } // Loop particles


}
//______________________________| Definitions of histo etc.
void AliAnalysisTaskJpsi::CreateHistos(){
  //
  // Init the histos
  //
  
  // binning 
  Int_t    nBinMass   = 800;
  Double_t MassMin    = 0.;
  Double_t MassMax    = 8.;
  Int_t    nBinPairPt = 400;
  Double_t PairPtMin  = 0.;
  Double_t PairPtMax  = 10.;

  // Jpsi histos
  fJPsiPt = new TH1D("hJpsi_pt","",nBinPairPt,PairPtMin,PairPtMax);
 
  // Electron histos
  fMee_Ptee_Jpsi_rec_die = new TH2D("MeePteeJPsi2ee",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_rec_rad = new TH2D("MeePteeJPsi2eeGamma",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_rec_die_low = new TH2D("MeePteeJPsi2ee_low",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_rec_rad_low = new TH2D("MeePteeJPsi2eeGamma_low",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_rec_die_high = new TH2D("MeePteeJPsi2ee_high",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_rec_rad_high = new TH2D("MeePteeJPsi2eeGamma_high",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);

  fMee_Ptee_Jpsi_gen_die = new TH2D("GenMeePteeJPsi2ee",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_gen_rad = new TH2D("GenMeePteeJPsi2eeGamma",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_gen_die_low = new TH2D("GenMeePteeJPsi2ee_low",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_gen_rad_low = new TH2D("GenMeePteeJPsi2eeGamma_low",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_gen_die_high = new TH2D("GenMeePteeJPsi2ee_high",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);
  fMee_Ptee_Jpsi_gen_rad_high = new TH2D("GenMeePteeJPsi2eeGamma_high",";m_{ee};p_{T,ee}",nBinMass,MassMin,MassMax,nBinPairPt,PairPtMin,PairPtMax);


  // Sumw2

  fJPsiPt->Sumw2();
  fMee_Ptee_Jpsi_rec_die->Sumw2();
  fMee_Ptee_Jpsi_rec_rad->Sumw2();
  fMee_Ptee_Jpsi_rec_die_low->Sumw2();
  fMee_Ptee_Jpsi_rec_rad_low->Sumw2();
  fMee_Ptee_Jpsi_rec_die_high->Sumw2();
  fMee_Ptee_Jpsi_rec_rad_high->Sumw2();
  fMee_Ptee_Jpsi_gen_die->Sumw2();
  fMee_Ptee_Jpsi_gen_rad->Sumw2();
  fMee_Ptee_Jpsi_gen_die_low->Sumw2();
  fMee_Ptee_Jpsi_gen_rad_low->Sumw2();
  fMee_Ptee_Jpsi_gen_die_high->Sumw2();
  fMee_Ptee_Jpsi_gen_rad_high->Sumw2();

  // Add to the list

  fOutputList->Add(fJPsiPt);
  fOutputList->Add(fMee_Ptee_Jpsi_rec_die);
  fOutputList->Add(fMee_Ptee_Jpsi_rec_rad);
  fOutputList_Low->Add(fMee_Ptee_Jpsi_rec_die_low);
  fOutputList_Low->Add(fMee_Ptee_Jpsi_rec_rad_low);
  fOutputList_High->Add(fMee_Ptee_Jpsi_rec_die_high);
  fOutputList_High->Add(fMee_Ptee_Jpsi_rec_rad_high);
  fOutputList->Add(fMee_Ptee_Jpsi_gen_die);
  fOutputList->Add(fMee_Ptee_Jpsi_gen_rad);
  fOutputList_Low->Add(fMee_Ptee_Jpsi_gen_die_low);
  fOutputList_Low->Add(fMee_Ptee_Jpsi_gen_rad_low);
  fOutputList_High->Add(fMee_Ptee_Jpsi_gen_die_high);
  fOutputList_High->Add(fMee_Ptee_Jpsi_gen_rad_high);
  
}

