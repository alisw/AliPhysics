/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//  * 20/04/2010 *
// Class for optimising and applying V0 cuts to obtain clean V0 samples
// Compatible with ESDs only
//
// Authors:
//    Matus Kalisky <matus.kalisky@cern.ch>
//

#include "TDatabasePDG.h"

#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"

#include "AliHFEV0cuts.h"

ClassImp(AliHFEV0cuts)

//________________________________________________________________
AliHFEV0cuts::AliHFEV0cuts():
  fQA(NULL)
  , fQAmc(NULL)
  , fMCEvent(NULL)
  , fInputEvent(NULL)
  , fPrimaryVertex(NULL)
  , fCurrentV0id(0)
  , fPdaughterPDG(0)
  , fNdaughterPDG(0)
  , fDestBits(0)
{
 
  //
  // Default constructor
  //
  

}
//________________________________________________________________
AliHFEV0cuts::~AliHFEV0cuts()
{
  //
  // destructor
  //
  
  if (fQA && TESTBIT(fDestBits, kBitQA)) delete fQA;
  if (fQAmc && TESTBIT(fDestBits, kBitQAmc)) delete fQAmc;
}

//________________________________________________________________
AliHFEV0cuts::AliHFEV0cuts(const AliHFEV0cuts &ref):
  TObject(ref)
  , fQA(NULL)
  , fQAmc(NULL)
  , fMCEvent(NULL)
  , fInputEvent(NULL)
  , fPrimaryVertex(NULL)
  , fCurrentV0id(0)
  , fPdaughterPDG(0)
  , fNdaughterPDG(0)
  , fDestBits(0)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);  
}
//________________________________________________________________
AliHFEV0cuts &AliHFEV0cuts::operator=(const AliHFEV0cuts &ref){
  //
  // Assignment operator
  //
  if(this != &ref)
    ref.Copy(*this);
  return *this;  
}
//________________________________________________________________
void AliHFEV0cuts::Copy(TObject &ref) const{
  //
  // Copy function
  // 
  AliHFEV0cuts &target = dynamic_cast<AliHFEV0cuts &>(ref);

  if(fQA) target.fQA = dynamic_cast<AliHFEcollection *>(fQA->Clone());  

  if(fQAmc) target.fQAmc = dynamic_cast<AliHFEcollection *>(fQAmc->Clone());  

  if(target.fMCEvent) delete target.fMCEvent;
  target.fMCEvent = new AliMCEvent;
  
  if(target.fPrimaryVertex) delete target.fPrimaryVertex;
  target.fPrimaryVertex = new AliKFVertex;

  TObject::Copy(ref);
  
}
//___________________________________________________________________
void AliHFEV0cuts::Init(const char* name){
  //
  // initialize the output objects and create histograms
  //

  //
  // all the "h_cut_XXX" histograms hare cut value distributions:
  // [0] for all candidates
  // [1] jus before the cut on given variable was applied, but after all the previous cuts
  //

  memset(&fDestBits, 0, sizeof(UInt_t));
  // now set the first two bits to 1
  SETBIT(fDestBits, kBitQA);
  SETBIT(fDestBits, kBitQAmc);
  

  fQA = new AliHFEcollection("fQA", name);

  fQAmc = new AliHFEcollection("fQAmc", name);

  // common for all V0s
  fQA->CreateTH2Fvector1(2, "h_all_AP", "armenteros plot for all V0 candidates", 200, -1, 1, 200, 0, 0.25);

  // gammas
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_CosPoint", "Gamma Cosine pointing angle; cos point. angle; counts", 100, 0, 0.1);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_DCA", "DCA between the gamma daughters; dca (cm); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_VtxR_old", "*old* Radius of the gamma conversion vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_VtxR", "Radius of the gamma conversion vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_PP", "gamma psi pair angle; psi pairangle (rad); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_Chi2", "gamma Chi2/NDF; Chi2/NDF; counts", 100, 0, 50);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_Sep", "gamma separation dist at TPC inned wall", 100, 0, 10);
  fQA->CreateTH1Fvector1(7, "h_Gamma_Mass", "Invariant mass of gammas; mass (GeV/c^{2}); counts", 100, 0, 0.2);
  
 
  // kaons
  fQA->CreateTH1Fvector1(2, "h_cut_K0_CosPoint", "K0 Cosine pointing angle; cos point. angle; counts", 100, 0, 0.1);
  fQA->CreateTH1Fvector1(2, "h_cut_K0_DCA", "DCA between the K0 daughters; dca (cm); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_K0_VtxR", "Radius of the K0 decay vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_K0_Chi2", "K0 Chi2/NDF; Chi2/NDF; counts", 100, 0, 50);
  fQA->CreateTH1Fvector1(5, "h_K0_Mass", "Invariant mass of K0; mass (GeV/c^{2}); counts", 125, 0.45, 0.55);

  // lambda
  fQA->CreateTH1Fvector1(2, "h_cut_L_CosPoint", "L Cosine pointing angle; cos point. angle; counts", 100, 0, 0.1);
  fQA->CreateTH1Fvector1(2, "h_cut_L_DCA", "DCA between the L daughters; dca (cm); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_L_VtxR", "Radius of the L decay vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_L_Chi2", "L Chi2/NDF; Chi2/NDF; counts", 100, 0, 50);
  fQA->CreateTH1Fvector1(5, "h_L_Mass", "Invariant mass of L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);
  fQA->CreateTH1Fvector1(5, "h_AL_Mass", "Invariant mass of anti L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);

  fQA->CreateTH2F("h_L_checks", "Lambda candidate check[0] -v- check[1]; check[0]; check[1]", 5, -0.75, 1.75, 6, -0.75, 1.75 );
  
  // electrons
  fQA->CreateTH1Fvector1(7, "h_Electron_P", "Momenta of conversion electrons -cuts-; P (GeV/c); counts", 50, 0.1, 20, 0);

  // K0 pions
  fQA->CreateTH1Fvector1(8, "h_PionK0_P", "Momenta of K0 pions -cuts-; P (GeV/c) counts;", 50, 0.1, 20, 0);
  
  // L pions
  fQA->CreateTH1Fvector1(9, "h_PionL_P", "Momenta of L pions -cuts-; P (GeV/c) counts;", 50, 0.1, 20, 0);
  
  // L protons
  fQA->CreateTH1Fvector1(9, "h_ProtonL_P", "Momenta of L protons -cuts-; P (GeV/c) counts;", 50, 0.1, 20, 0);    

  // single track cuts
  fQA->CreateTH1F("h_ST_NclsTPC", "Number of TPC clusters", 161, -1, 160);
  fQA->CreateTH1F("h_ST_TPCrefit", "TPC refit", 2, -0.5, 1.5);
  fQA->CreateTH1F("h_ST_chi2TPCcls", "chi2 per TPC cluster", 100, 0, 10);
  fQA->CreateTH1F("h_ST_TPCclsR", "TPC cluster ratio", 120, -0.1, 1.1);
  fQA->CreateTH1F("h_ST_kinks", "kinks", 2, -0.5, 1.5);
  fQA->CreateTH1F("h_ST_pt", "track pt", 100, 0.1, 20, 0);
  fQA->CreateTH1F("h_ST_eta", "track eta", 100, -1.5, 1.5);

  //
  // possibly new cuts
  //
 fQA->CreateTH2Fvector1(2, "h_cut_L_rdp_v_mp", "relative L daughter mom -v- mother mom; L mom (GeV/c); relative daughter mom p2/p1", 100, 0.1, 10, 100, 0, 1);

  // THnSparse histograms
  
  // THnSparse for the K0 mass
  // to be looked at after merging run by run
  // axes: mass, pt, theta, phi
  {
    Int_t nBin[4] = {100, 10, 10, 18};
    Double_t nMin[4] = {0.45, 0.1, 0., 0.};
    Double_t nMax[4] = {0.55, 10., TMath::Pi(), 2*TMath::Pi()};
    TString htitle = "K0 sparse; mass (GeV/c^{2}); p_{T} (GeV/c); theta (rad); phi(rad)";
    fQA->CreateTHnSparse("hK0", htitle, 4, nBin, nMin, nMax);
    fQA->BinLogAxis("hK0", 1);
  }


  // 
  // MC plots for checking and tuning the V0 cuts
  //
 
  const char *v0[4] = {"G", "K", "L"}; // to keep the names short
  // number of V0s left after each cut step - for signal and background - within given mass window
  for(Int_t i=0; i<3; ++i){
    fQAmc->CreateTH1F(Form("h_%s_cuts_S", v0[i]), Form("h_%s_cuts_S", v0[i]), 10, -0.5, 9.5);
    fQAmc->CreateTH1F(Form("h_%s_cuts_B", v0[i]), Form("h_%s_cuts_B", v0[i]), 10, -0.5, 9.5);
  }

  //
  // cut distributions for signal and background
  //

  const Float_t pMin = 0.1;
  const Float_t pMax = 10.;
  const Int_t pN = 12;


  // gamma signal
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_CosPoint_S", "S - Gamma Cosine pointing angle; mom (GeV/c); cos point. angle",  pN, pMin, pMax, 50, 0, 0.1, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_DCA_S", "S - DCA between the gamma daughters; mom (GeV/c); dca (cm)", pN, pMin, pMax, 50, 0, 2, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_VtxR_S", "S - Radius of the gamma conversion vertex; mom (GeV/c); r (cm)", pN, pMin, pMax, 100, 0, 100, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_PP_S", "S - gamma psi pair angle; mom (GeV/c); psi pairangle (rad)", pN, pMin, pMax, 50, 0, 0.5, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_Chi2_S", "S - gamma Chi2/NDF; mom (GeV/c); Chi2/NDF", pN, pMin, pMax, 50, 0, 25, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_Sep_S", "S - gamma separation TPC-inner; mom (GeV/c); tracks separatin (cm)", pN, pMin, pMax, 100, 0, 10, 0);
  // as a function of radius, not momentum
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_SepR_S", "S - gamma separation TPC-inner; radius cm; tracks separatin (cm)", 20, 0, 100, 100, 0, 20);

  fQAmc->CreateTH1Fvector1(9, "h_Gamma_Mass_S", "S - Invariant mass of gammas; mass (GeV/c^{2}); counts", 100, 0, 0.2);
  // gamma background
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_CosPoint_B", "B - Gamma Cosine pointing angle; mom (GeV/c); cos point. angle", pN, pMin, pMax, 50, 0, 0.1, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_DCA_B", "B - DCA between the gamma daughters; mom (GeV/c); dca (cm)", pN, pMin, pMax, 50, 0, 2, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_VtxR_B", "B - Radius of the gamma conversion vertex; mom (GeV/c); r (cm)", pN, pMin, pMax, 100, 0, 100, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_PP_B", "B - gamma psi pair angle; mom (GeV/c); psi pairangle (rad)", pN, pMin, pMax, 50, 0, 0.5, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_Chi2_B", "B - gamma Chi2/NDF; mom (GeV/c); Chi2/NDF", pN, pMin, pMax, 50, 0, 25, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_Sep_B", "B - gamma separation TPC-inner; mom (GeV/c); tracks separatin (cm)", pN, pMin, pMax, 100, 0, 50, 0);
  //
  fQAmc->CreateTH2Fvector1(2, "h_cut_Gamma_SepR_B", "S - gamma separation TPC-inner; radius cm; tracks separatin (cm)", 20, 0, 100, 100, 0, 20);
  fQAmc->CreateTH1Fvector1(9, "h_Gamma_Mass_B", "B - Invariant mass of gammas; mass (GeV/c^{2}); counts", 100, 0, 0.2);  
 
  // kaons signal
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_CosPoint_S", "S - K0 Cosine pointing angle; mom (GeV/c); cos point. angle", pN, pMin, pMax, 50, 0, 0.1, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_DCA_S", "S - DCA between the K0 daughters; mom (GeV/c); dca (cm)", pN, pMin, pMax, 50, 0, 2, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_VtxR_S", "S - Radius of the K0 decay vertex; mom (GeV/c); r (cm)", pN, pMin, pMax, 50, 0, 100, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_Chi2_S", "S - K0 Chi2/NDF; mom (GeV/c); Chi2/NDF", pN, pMin, pMax, 50, 0, 25, 0);

  fQAmc->CreateTH1Fvector1(5, "h_K0_Mass_S", "S - Invariant mass of K0; mass (GeV/c^{2}); counts", 125, 0.45, 0.55);
  // kaons background
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_CosPoint_B", "B - K0 Cosine pointing angle; mom (GeV/c); cos point. angle", pN, pMin, pMax, 50, 0, 0.1, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_DCA_B", "B - DCA between the K0 daughters; mom (GeV/c); dca (cm)", pN, pMin, pMax, 50, 0, 2, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_VtxR_B", "B - Radius of the K0 decay vertex; mom (GeV/c); r (cm)", pN, pMin, pMax, 50, 0, 100, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_K0_Chi2_B", "B - K0 Chi2/NDF; mom (GeV/c); Chi2/NDF", pN, pMin, pMax, 50, 0, 50, 0);

  fQAmc->CreateTH1Fvector1(5, "h_K0_Mass_B", "B - Invariant mass of K0; mass (GeV/c^{2}); counts", 125, 0.45, 0.55);

  // lambda signal
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_CosPoint_S", "S - L Cosine pointing angle; mom (GeV/c); cos point. angle", pN, pMin, pMax, 50, 0, 0.1, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_DCA_S", "S - DCA between the L daughters; mom (GeV/c); dca (cm)", pN, pMin, pMax, 50, 0, 2, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_VtxR_S", "S - Radius of the L decay vertex; mom (GeV/c); r (cm)", pN, pMin, pMax, 50, 0, 100, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_Chi2_S", "S - L Chi2/NDF; mom (GeV/c); Chi2/NDF", pN, pMin, pMax, 50, 0, 50, 0);

  fQAmc->CreateTH1Fvector1(5, "h_L_Mass_S", "S - Invariant mass of L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);
  fQAmc->CreateTH1Fvector1(5, "h_AL_Mass_S", "S - Invariant mass of anti L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);
  // lambda background
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_CosPoint_B", "B - L Cosine pointing angle; mom (GeV/c); cos point. angle", pN, pMin, pMax, 50, 0, 0.1, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_DCA_B", "B - DCA between the L daughters; mom (GeV/c); dca (cm)", pN, pMin, pMax, 50, 0, 2, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_VtxR_B", "B - Radius of the L decay vertex; mom (GeV/c); r (cm)", pN, pMin, pMax, 50, 0, 100, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_Chi2_B", "B - L Chi2/NDF; mom (GeV/c); Chi2/NDF", pN, pMin, pMax, 50, 0, 50, 0);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_rdp_v_mp_S", "S - relative L daughter mom -v- mother mom; L mom (GeV/c); relative daughter mom p2/p1", 100, 0.1, 10, 100, 0, 1);
  fQAmc->CreateTH2Fvector1(2, "h_cut_L_rdp_v_mp_B", "B - relative L daughter mom -v- mother mom; L mom (GeV/c); relative daughter mom p2/p1", 100, 0.1, 10, 100, 0, 1);
  fQAmc->CreateTH1Fvector1(5, "h_LAL_Mass_B", "B - Invariant mass of anti L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);


  // MC tagged daughter track momentum distribution after each cut step
//   fQAmc->CreateTH1Fvector1(10, "h_electron_p_S", "h_electron_p_S", 20, 0.1, 20, 0);
//   fQAmc->CreateTH1Fvector1(10, "h_K0pion_p_S", "h_K0pion_p_S", 20, 0.1, 20, 0);
//   fQAmc->CreateTH1Fvector1(10, "h_Lpion_p_S", "h_Lpion_p_S", 20, 0.1, 20, 0);
//   fQAmc->CreateTH1Fvector1(10, "h_proton_p_S", "h_proton_p_S", 20, 0.1, 20, 0);

  // V0 momnetum distribution of MC tagged signal and backglound after all cuts
  fQAmc->CreateTH1F("h_gamma_p_S", "true gammas after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_gamma_p_B", "true gamma BG after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_K0_p_S", "true K0s after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_K0_p_B", "true K0 BG after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_lambda_p_S", "MC true lambdas after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_lambda_p_B", "MC true lambda BG after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_alambda_p_S", "MC true anti-lambdas after all cuts", 20, 0.1, 10, 0);
  fQAmc->CreateTH1F("h_alambda_p_B", "MC true anti-lambda BG after all cuts", 20, 0.1, 10, 0);  

  // invariant mass ditributions for the V0 for different hypoteses (gamma, K0, L, AL)
  fQAmc->CreateTH1F("h_Mass_gamma_as_K0","h_Mass_gamma_as_K0", 200, 0, 2);
  fQAmc->CreateTH1F("h_Mass_gamma_as_L","h_Mass_gamma_as_L", 200, 0, 2);
  fQAmc->CreateTH1F("h_Mass_K0_as_G", "h_Mass_K0_as_gamma", 200, 0, 2);
  fQAmc->CreateTH1F("h_Mass_K0_as_L", "h_Mass_K0_as_Lambda", 200, 0, 2);
  fQAmc->CreateTH1F("h_Mass_L_as_G", "h_Mass_L_as_gamma", 200, 0, 2);
  fQAmc->CreateTH1F("h_Mass_L_as_K0", "h_Mass_L_as_K0", 200, 0, 2);

  // Invariant mass distribution of MC tagged signal for diffrent momenta
  fQAmc->CreateTH2F("h_gamma_MvP_S", "mc tagged gammas - signal; p (GeV/c); m (GeV/c^{2})", 12, 0.1, 20, 100, 0., 0.1, 0);
  fQAmc->CreateTH2F("h_K0_MvP_S", "mc tagged K0s - signal; p (GeV/c); m (GeV/c^{2})", 12, 0.1, 20, 100, 0.45, 0.55, 0);
  fQAmc->CreateTH2F("h_lambda_MvP_S", "mc tagged Lambdas - signal; p (GeV/c); m (GeV/c^{2})", 12, 0.1, 20, 100, 1.08, 1.14, 0);
    
  // electrons
  fQAmc->CreateTH1Fvector1(8, "h_Electron_P_S", "MC-S momenta of conversion electrons -cuts-; P (GeV/c); counts", 20, 0.1, 20, 0);
  fQAmc->CreateTH1Fvector1(8, "h_Electron_P_B", "MC-B momenta of conversion electrons -cuts-; P (GeV/c); counts", 20, 0.1, 20, 0);

  // K0 pions
  fQAmc->CreateTH1Fvector1(7, "h_PionK0_P_S", "MC-S momenta of K0 pions -cuts-; P (GeV/c) counts;", 20, 0.1, 20, 0);
  fQAmc->CreateTH1Fvector1(7, "h_PionK0_P_B", "MC-B momenta of K0 pions -cuts-; P (GeV/c) counts;", 20, 0.1, 20, 0);
  
  // L pions
  fQAmc->CreateTH1Fvector1(8, "h_PionL_P_S", "MC-S momenta of L pions -cuts-; P (GeV/c) counts;", 20, 0.1, 50, 0);
  fQAmc->CreateTH1Fvector1(8, "h_PionL_P_B", "MC-B momenta of L pions -cuts-; P (GeV/c) counts;", 20, 0.1, 50, 0);
  
  // L protons
  fQAmc->CreateTH1Fvector1(8, "h_ProtonL_P_S", "MC-S momenta of L protons -cuts-; P (GeV/c) counts;", 20, 0.1, 20, 0);    
  fQAmc->CreateTH1Fvector1(8, "h_ProtonL_P_B", "MC-B momenta of L protons -cuts-; P (GeV/c) counts;", 20, 0.1, 20, 0);    



  // cut efficiencies 
}
//________________________________________________________________
Bool_t AliHFEV0cuts::TrackCutsCommon(AliESDtrack* track){
  //
  // singe track cuts commom for all particle candidates
  //

  if(!track) return kFALSE;
 
  
  // status word
  ULong_t status = track->GetStatus();


  // No. of TPC clusters
  fQA->Fill("h_ST_NclsTPC", track->GetTPCNcls());
  if(track->GetTPCNcls() < 1) return kFALSE;   //

  // TPC refit
  if((status & AliESDtrack::kTPCrefit)){
    fQA->Fill("h_ST_TPCrefit", 1);
  }
  if(!(status & AliESDtrack::kTPCrefit)){
    fQA->Fill("h_ST_TPCrefit", 0);
    return kFALSE;
  }

  // Chi2 per TPC cluster
  Int_t nTPCclusters = track->GetTPCclusters(0);
  Float_t chi2perTPCcluster = track->GetTPCchi2()/Float_t(nTPCclusters);
  fQA->Fill("h_ST_chi2TPCcls", chi2perTPCcluster);
  if(chi2perTPCcluster > 4.0) return kFALSE;   // 4.0

  // TPC cluster ratio
  Float_t cRatioTPC = track->GetTPCNclsF() > 0. ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t> (track->GetTPCNclsF()) : 1.;
  fQA->Fill("h_ST_TPCclsR", cRatioTPC);
  if(cRatioTPC < 0.6) return kFALSE;

  // kinks
  fQA->Fill("h_ST_kinks", track->GetKinkIndex(0));
  if(track->GetKinkIndex(0) != 0) return kFALSE;

  // pt
  fQA->Fill("h_ST_pt",track->Pt());
  //if(track->Pt() < 0.1 || track->Pt() > 100) return kFALSE; //

  // eta
  fQA->Fill("h_ST_eta", track->Eta());
  //if(TMath::Abs(track->Eta()) > 0.9) return kFALSE;  

  return kTRUE;
}
//________________________________________________________________
Bool_t AliHFEV0cuts::V0CutsCommon(AliESDv0 *v0){
  //
  // V0 cuts common to all V0s
  //

  AliESDtrack* dN, *dP; 
 
  dP = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(v0->GetPindex()));
  dN = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(v0->GetNindex())); 
  
  if(!dN || !dP) return kFALSE;

  Int_t qP = dP->Charge();
  Int_t qN = dN->Charge();

  if((qP*qN) != -1) return kFALSE;

  return kTRUE;
}
//________________________________________________________________
Bool_t AliHFEV0cuts::GammaCuts(AliESDv0 *v0){
  //
  // gamma cuts 
  //
  
  if(!v0) return kFALSE;

  if(fMCEvent){
    if(1 == fCurrentV0id){
      fQAmc->Fill("h_Mass_gamma_as_K0",  v0->GetEffMass(2, 2));
      fQAmc->Fill("h_Mass_gamma_as_L",  v0->GetEffMass(2, 4));
      fQAmc->Fill("h_Mass_gamma_as_L",  v0->GetEffMass(4, 2));
    }
  }

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  if(CheckSigns(v0)){
    pIndex = v0->GetPindex();
    nIndex = v0->GetNindex();
  }
  else{
    pIndex = v0->GetNindex();
    nIndex = v0->GetPindex();    
  }
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  AliKFParticle *kfMother = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kElectron), TMath::Abs(kElectron));
  if(!kfMother) return kFALSE;
  
  // production vertex is set in the 'CreateMotherParticle' function
  //kfMother->SetMassConstraint(0, 0.001);

  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(pIndex));
  d[1] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(nIndex));

  Float_t iMass = v0->GetEffMass(0, 0);
  Float_t iP = v0->P();
  Float_t p[2] = {d[0]->GetP(), d[1]->GetP()};

  // Cut values
  const Double_t cutChi2NDF = 1.5;              // ORG [7.]  
  const Double_t cutCosPoint[2] = {0., 0.007};  // ORG [0., 0.02]
  const Double_t cutDCA[2] = {0., 0.25};       // ORG [0., 0.25]
  const Double_t cutProdVtxR[2] = {6., 90.};   // ORG [3., 90]
  const Double_t cutPsiPair[2] = {0., 0.05};   // ORG [0. 0.05]
  // mass cut
  const Double_t cutMass = 0.05;               // ORG [0.05]

  //
  // possible new cuts
  //
  // separation cut at the entrance to the TPC
  const Double_t cutSeparation = 1.;          // ORG 3.0 cm



  // Values
   
  // cos pointing angle
  Double_t cosPoint = v0->GetV0CosineOfPointingAngle();
  cosPoint = TMath::ACos(cosPoint);

  // DCA between daughters
  Double_t dca = v0->GetDcaV0Daughters();

  // Production vertex
  Double_t x, y, z; 
  v0->GetXYZ(x,y,z);
  Double_t r = TMath::Sqrt(x*x + y*y);

  Double_t xy[2];
  Double_t r2 = -1.;
  if ( GetConvPosXY(d[0], d[1], xy) ){
    r2 = TMath::Sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
  }

  // psi pair 
  Double_t psiPair = PsiPair(v0);
  
  // V0 chi2/ndf
  Double_t chi2ndf = kfMother->GetChi2()/kfMother->GetNDF();
  if(kfMother) delete kfMother; 

  // Separation
  AliExternalTrackParam const *param[2];
  param[0] = d[0]->GetInnerParam();
  param[1] = d[1]->GetInnerParam();
  Double_t sep = 999.;
  if(param[0] && param[1]){
    TVector3 xyz[3];
    xyz[0].SetXYZ(param[0]->GetX(), param[0]->GetY(), param[0]->GetZ());
    xyz[1].SetXYZ(param[1]->GetX(), param[1]->GetY(), param[1]->GetZ());
    xyz[2] = xyz[0] - xyz[1];
    sep = xyz[2].Mag();
  }


 
  //
  // Apply the cuts, produce QA plots (with mass cut)
  //
  fQA->Fill("h_Gamma_Mass", 0, iMass);
  
  // MC
  if(fMCEvent){
    if(1 == fCurrentV0id){
      fQAmc->Fill("h_Gamma_Mass_S", 0, iMass);
      fQAmc->Fill("h_gamma_MvP_S", iP, iMass);
    }
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_B", 0, iMass);
  }
  // cut distributions
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 0, p[0]);
    fQA->Fill("h_Electron_P", 0, p[1]);
    fQA->Fill("h_cut_Gamma_CosPoint", 0, cosPoint);
    fQA->Fill("h_cut_Gamma_DCA", 0, dca);
    fQA->Fill("h_cut_Gamma_VtxR_old", 0, r);
    fQA->Fill("h_cut_Gamma_VtxR", 0, r2);
    fQA->Fill("h_cut_Gamma_PP", 0, psiPair);
    fQA->Fill("h_cut_Gamma_Chi2", 0, chi2ndf);
    fQA->Fill("h_cut_Gamma_Chi2", 1, chi2ndf, iP);
    fQA->Fill("h_cut_Gamma_Sep", 0, iP, sep);
	
   
    if(fMCEvent){
      // MC signal
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_CosPoint_S", 0, iP, cosPoint);
	fQAmc->Fill("h_cut_Gamma_DCA_S", 0, iP, dca);
	fQAmc->Fill("h_cut_Gamma_VtxR_S", 0, iP, r2);
	fQAmc->Fill("h_cut_Gamma_PP_S", 0, iP, psiPair);
	fQAmc->Fill("h_cut_Gamma_Chi2_S", 0, iP, chi2ndf);
	fQAmc->Fill("h_cut_Gamma_Chi2_S", 1, iP, chi2ndf);
	fQAmc->Fill("h_cut_Gamma_Sep_S", 0, iP, sep);
	fQAmc->Fill("h_cut_Gamma_SepR_S", 0,r2, sep);
	fQAmc->Fill("h_Electron_P_S", 0, p[0]);
	fQAmc->Fill("h_Electron_P_S", 0, p[1]);
      }
      // MC background
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_CosPoint_B", 0, iP, cosPoint);
	fQAmc->Fill("h_cut_Gamma_DCA_B", 0, iP, dca);
	fQAmc->Fill("h_cut_Gamma_VtxR_B", 0, iP, r2);
	fQAmc->Fill("h_cut_Gamma_PP_B", 0, iP, psiPair);
	fQAmc->Fill("h_cut_Gamma_Chi2_B", 0, iP, chi2ndf);
	fQAmc->Fill("h_cut_Gamma_Chi2_B", 1, iP, chi2ndf);
	fQAmc->Fill("h_cut_Gamma_Sep_B", 0, iP, sep);
	fQAmc->Fill("h_cut_Gamma_SepR_B", 0,r2, sep);
	fQAmc->Fill("h_Electron_P_B", 0, p[0]);
	fQAmc->Fill("h_Electron_P_B", 0, p[1]);	
      }
    }
  }


  //
  // Chi2/NDF cut
  //
  if(chi2ndf > cutChi2NDF) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 1, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_cut_Gamma_CosPoint", 1, cosPoint);
    fQA->Fill("h_Electron_P", 1, p[0]);
    fQA->Fill("h_Electron_P", 1, p[1]);
  }
  if(fMCEvent){
    if(1 == fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_S", 1, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_B", 1, iMass);
    if(iMass < cutMass){
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_CosPoint_S", 1, iP, cosPoint);
	fQAmc->Fill("h_Electron_P_S", 1, p[0]);
	fQAmc->Fill("h_Electron_P_S", 1, p[1]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_CosPoint_B", 1, iP, cosPoint);
	fQAmc->Fill("h_Electron_P_B", 1, p[0]);
	fQAmc->Fill("h_Electron_P_B", 1, p[1]);
      }
    }
  }

  //
  // Cos point cut
  //
  if(cosPoint < cutCosPoint[0] || cosPoint > cutCosPoint[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 2, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 2, p[0]);
    fQA->Fill("h_Electron_P", 2, p[1]);
    fQA->Fill("h_cut_Gamma_DCA", 1, dca);
  }
  if(fMCEvent){
    if(1 == fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_S", 2, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_B", 2, iMass);
    if(iMass < cutMass){
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_DCA_S", 1, iP, dca);
	fQAmc->Fill("h_Electron_P_S", 2, p[0]);
	fQAmc->Fill("h_Electron_P_S", 2, p[1]);

      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_DCA_B", 1, iP, dca);
	fQAmc->Fill("h_Electron_P_B", 2, p[0]);
	fQAmc->Fill("h_Electron_P_B", 2, p[1]);

      }
    }
  }
  
  //
  // DCA cut
  //
  if(dca < cutDCA[0] || dca > cutDCA[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 3, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 3, p[0]);
    fQA->Fill("h_Electron_P", 3, p[1]);
    fQA->Fill("h_cut_Gamma_VtxR_old", 1, r);
    fQA->Fill("h_cut_Gamma_VtxR", 1, r2);
  }
  if(fMCEvent){
    if(1 == fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_S", 3, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_B", 3, iMass);
    if(iMass < cutMass){
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_VtxR_S", 1, iP, r2);
	fQAmc->Fill("h_Electron_P_S", 3, p[0]);
	fQAmc->Fill("h_Electron_P_S", 3, p[1]);

      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_VtxR_B", 1, iP, r2);
	fQAmc->Fill("h_Electron_P_B", 3, p[0]);
	fQAmc->Fill("h_Electron_P_B", 3, p[1]);	
      }
    }
  }

  //
  // Vertex radius cut
  //
  if(r < cutProdVtxR[0] || r > cutProdVtxR[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 4, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_cut_Gamma_PP", 1, psiPair);
    fQA->Fill("h_Electron_P", 4, p[0]);
    fQA->Fill("h_Electron_P", 4, p[1]);
  }
  if(fMCEvent){
    if(1 == fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_S", 4, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_B", 4, iMass);
    if(iMass < cutMass){
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_PP_S", 1, iP, psiPair);
	fQAmc->Fill("h_Electron_P_S", 4, p[0]);
	fQAmc->Fill("h_Electron_P_S", 4, p[1]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_PP_B", 1, iP, psiPair);
	fQAmc->Fill("h_Electron_P_B", 4, p[0]);
	fQAmc->Fill("h_Electron_P_B", 4, p[1]);
      }
    }
  }


  //
  // PsiPair cut
  //
  if(psiPair < cutPsiPair[0] || psiPair > cutPsiPair[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 5, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_cut_Gamma_Sep", 1, iP, sep);
    fQA->Fill("h_Electron_P", 5, p[0]);
    fQA->Fill("h_Electron_P", 5, p[1]);
  }
  if(fMCEvent){
    if(1 == fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_S", 5, iMass);
    else if(-2 != fCurrentV0id)fQAmc->Fill("h_Gamma_Mass_B", 5, iMass);
    
    if(iMass < cutMass){
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_Sep_S", 1, iP, sep);
	fQAmc->Fill("h_cut_Gamma_SepR_S", 1,r2, sep);
	fQAmc->Fill("h_Electron_P_S", 5, p[0]);
	fQAmc->Fill("h_Electron_P_S", 5, p[1]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_Gamma_Sep_B", 1, iP, sep);
	fQAmc->Fill("h_cut_Gamma_SepR_B", 1,r2, sep);
	fQAmc->Fill("h_Electron_P_B", 5, p[0]);
	fQAmc->Fill("h_Electron_P_B", 5, p[1]);
      }
    }
  }


  // TESTING NEW CUT
  //
  // distance of the tracks at the entrance of the TPC
  //
  if(sep < cutSeparation) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 6, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 6, p[0]);
    fQA->Fill("h_Electron_P", 6, p[1]);
  }
  if(fMCEvent){
    if(1 == fCurrentV0id) fQAmc->Fill("h_Gamma_Mass_S", 6, iMass);
    else if(-2 != fCurrentV0id)fQAmc->Fill("h_Gamma_Mass_B", 6, iMass);
    
    if(iMass < cutMass){
      if(1 == fCurrentV0id){
	fQAmc->Fill("h_Electron_P_S", 6, p[0]);
	fQAmc->Fill("h_Electron_P_S", 6, p[1]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_Electron_P_B", 6, p[0]);
	fQAmc->Fill("h_Electron_P_B", 6, p[1]);
      }
    }
  }
  
  // .. test
 

  if(iMass > cutMass) return kFALSE;

  // all cuts passed

  
  // some MC stuff
  //printf("**D: gamma V0id: %i, P: %i, N: %i \n", fCurrentV0id, fPdaughterPDG, fNdaughterPDG);
  if(1 == fCurrentV0id){
    fQAmc->Fill("h_gamma_p_S", iP);
    fQAmc->Fill("h_Electron_P_S", 7, p[0]);
    fQAmc->Fill("h_Electron_P_S", 7, p[1]);
  }
  else if (-2 != fCurrentV0id){
    fQAmc->Fill("h_gamma_p_B", iP);
    fQAmc->Fill("h_Electron_P_B", 7, p[0]);
    fQAmc->Fill("h_Electron_P_B", 7, p[1]);
  }


  return kTRUE;
}
//________________________________________________________________
Bool_t AliHFEV0cuts::K0Cuts(AliESDv0 *v0){
  //
  // K0 cuts
  //

  if(!v0) return kFALSE;

  if(fMCEvent){
    if(2 == fCurrentV0id){
      fQAmc->Fill("h_Mass_K0_as_G",  v0->GetEffMass(0, 0));
      fQAmc->Fill("h_Mass_K0_as_L",  v0->GetEffMass(2, 4));
      fQAmc->Fill("h_Mass_K0_as_L",  v0->GetEffMass(4, 2));
    }
  }

  //const Double_t cK0mass=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();  // PDG K0s mass
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  if(CheckSigns(v0)){
    pIndex = v0->GetPindex();
    nIndex = v0->GetNindex();
  }
  else{
    pIndex = v0->GetNindex();
    nIndex = v0->GetPindex();    
  }
 
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  AliKFParticle *kfMother = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kPiPlus));
  if(!kfMother) return kFALSE;
  // production vertex is set in the 'CreateMotherParticle' function
  //kfMother->SetMassConstraint(cK0mass, 0.);
  
  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(pIndex));
  d[1] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(nIndex));

  Float_t iMass = v0->GetEffMass(2, 2);
  Float_t iP = v0->P();
  Float_t p[2] = {d[0]->GetP(), d[1]->GetP()};
  Double_t theta = v0->Theta();
  Double_t phi = v0->Phi();
  Double_t pt = v0->Pt();
  Double_t data[4] = {0., 0., 0., 0.};

  // Cut values
  const Double_t cutChi2NDF = 2.;              // ORG [7.]
  const Double_t cutCosPoint[2] = {0., 0.02};  // ORG [0., 0.03]
  const Double_t cutDCA[2] = {0., 0.2};        // ORG [0., 0.1]
  const Double_t cutProdVtxR[2] = {2.0, 30.};   // ORG [0., 8.1]
  const Double_t cutMass[2] = {0.486, 0.508};   // ORG [0.485, 0.51]
  // Values

  // cos pointing angle
  Double_t cosPoint = v0->GetV0CosineOfPointingAngle();
  cosPoint = TMath::ACos(cosPoint);

  // DCA between daughters
  Double_t dca = v0->GetDcaV0Daughters();

  // Production vertex
  Double_t x, y, z; 
  v0->GetXYZ(x,y,z);

  Double_t r = TMath::Sqrt(x*x + y*y);  

  // V0 chi2/ndf
  Double_t chi2ndf = kfMother->GetChi2()/kfMother->GetNDF();
  
  if(kfMother) delete kfMother; 

  //
  // Apply the cuts, produce QA plots (with mass cut)
  //

  fQA->Fill("h_K0_Mass", 0, iMass);
  // MC
  if(fMCEvent){
    if(2 == fCurrentV0id){
      fQAmc->Fill("h_K0_Mass_S", 0, iMass);
      fQAmc->Fill("h_K0_MvP_S", iP, iMass);
    }
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_K0_Mass_B", 0, iMass);
  }

  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 0, p[0]);
    fQA->Fill("h_PionK0_P", 0, p[1]);
    fQA->Fill("h_cut_K0_CosPoint", 0, cosPoint);
    fQA->Fill("h_cut_K0_DCA", 0, dca);
    fQA->Fill("h_cut_K0_VtxR", 0, r);
    fQA->Fill("h_cut_K0_Chi2", 0, chi2ndf);
    fQA->Fill("h_cut_K0_Chi2", 1, chi2ndf);
  }
  
  // MC
  if(fMCEvent){
    if(iMass > cutMass[0] && iMass < cutMass[1]){   
      if(2 == fCurrentV0id){
	fQAmc->Fill("h_cut_K0_CosPoint_S", 0, iP, cosPoint);
	fQAmc->Fill("h_cut_K0_DCA_S", 0, iP, dca);
	fQAmc->Fill("h_cut_K0_VtxR_S", 0, iP, r);
	fQAmc->Fill("h_cut_K0_Chi2_S", 0, iP, chi2ndf);
	fQAmc->Fill("h_cut_K0_Chi2_S", 1, iP, chi2ndf);
	fQAmc->Fill("h_PionK0_P_S", 0, p[0]);
	fQAmc->Fill("h_PionK0_P_S", 0, p[1]);
       }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_K0_CosPoint_B", 0, iP, cosPoint);
	fQAmc->Fill("h_cut_K0_DCA_B", 0, iP, dca);
	fQAmc->Fill("h_cut_K0_VtxR_B", 0, iP, r);
	fQAmc->Fill("h_cut_K0_Chi2_B", 0, iP, chi2ndf);
	fQAmc->Fill("h_cut_K0_Chi2_B", 1, iP, chi2ndf);
	fQAmc->Fill("h_PionK0_P_B", 0, p[0]);
	fQAmc->Fill("h_PionK0_P_B", 0, p[1]);
      }  
    }
  }

  //
  // Chi2/NDF cut
  //
  if(chi2ndf > cutChi2NDF) return kFALSE;
  fQA->Fill("h_K0_Mass", 1, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_cut_K0_CosPoint", 1, cosPoint);
    fQA->Fill("h_PionK0_P", 1, p[0]);
    fQA->Fill("h_PionK0_P", 1, p[1]);
  }  
  if(fMCEvent){
    if(2 == fCurrentV0id) fQAmc->Fill("h_K0_Mass_S", 1, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_K0_Mass_B", 1, iMass);
     if(iMass > cutMass[0] && iMass < cutMass[1]){
       if(2 == fCurrentV0id){
	 fQAmc->Fill("h_cut_K0_CosPoint_S", 1, iP, cosPoint);
	 fQAmc->Fill("h_PionK0_P_S", 1, p[0]);
	 fQAmc->Fill("h_PionK0_P_S", 1, p[1]);
       }
       else if(-2 != fCurrentV0id){
	 fQAmc->Fill("h_cut_K0_CosPoint_B", 1, iP, cosPoint);
	 fQAmc->Fill("h_PionK0_P_B", 1, p[0]);
	 fQAmc->Fill("h_PionK0_P_B", 1, p[1]);
       }
     }
  }

  //
  // Cos point cut
  //
  if(cosPoint < cutCosPoint[0] || cosPoint > cutCosPoint[1]) return kFALSE;
  fQA->Fill("h_K0_Mass", 2, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 2, p[0]);
    fQA->Fill("h_PionK0_P", 2, p[1]);
    fQA->Fill("h_cut_K0_DCA", 1, dca);
  }
  if(fMCEvent){
    if(2 == fCurrentV0id) fQAmc->Fill("h_K0_Mass_S", 2, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_K0_Mass_B", 2, iMass);
     if(iMass > cutMass[0] && iMass < cutMass[1]){
       if(2 == fCurrentV0id){
	 fQAmc->Fill("h_cut_K0_DCA_S", 1, iP, dca);
	 fQAmc->Fill("h_PionK0_P_S", 2, p[0]);
	 fQAmc->Fill("h_PionK0_P_S", 2, p[1]);
       }
       else if(-2 != fCurrentV0id){
	 fQAmc->Fill("h_cut_K0_DCA_B", 1, iP, dca);
	 fQAmc->Fill("h_PionK0_P_B", 2, p[0]);
	 fQAmc->Fill("h_PionK0_P_B", 2, p[1]);
       }
     }
  }
  

  //
  // DCA cut
  //
  if(dca < cutDCA[0] || dca > cutDCA[1]) return kFALSE;
  fQA->Fill("h_K0_Mass", 3, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 3, p[0]);
    fQA->Fill("h_PionK0_P", 3, p[1]);
    fQA->Fill("h_cut_K0_VtxR", 1, r);
  }
  if(fMCEvent){
    if(2 == fCurrentV0id) fQAmc->Fill("h_K0_Mass_S", 3, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_K0_Mass_B", 3, iMass);
     if(iMass > cutMass[0] && iMass < cutMass[1]){
       if(2 == fCurrentV0id){
	 fQAmc->Fill("h_cut_K0_VtxR_S", 1, iP, r);
	 fQAmc->Fill("h_PionK0_P_S", 3, p[0]);
	 fQAmc->Fill("h_PionK0_P_S", 3, p[1]);
       }
       else if(-2 != fCurrentV0id){
	 fQAmc->Fill("h_cut_K0_VtxR_B", 1, iP, r);
	 fQAmc->Fill("h_PionK0_P_B", 3, p[0]);
	 fQAmc->Fill("h_PionK0_P_B", 3, p[1]);
       }
     }
  }

  
  //
  // Vertex R cut
  //
  if(r < cutProdVtxR[0] || r > cutProdVtxR[1]) return kFALSE;
  fQA->Fill("h_K0_Mass", 4, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 4, p[0]);
    fQA->Fill("h_PionK0_P", 4, p[1]);
  }
  if(fMCEvent){
    if(2 == fCurrentV0id) fQAmc->Fill("h_K0_Mass_S", 4, iMass);
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_K0_Mass_B", 4, iMass);
     if(iMass > cutMass[0] && iMass < cutMass[1]){
       if(2 == fCurrentV0id){
	 fQAmc->Fill("h_PionK0_P_S", 4, p[0]);
	 fQAmc->Fill("h_PionK0_P_S", 4, p[1]);
       }
       else if(-2 != fCurrentV0id){
	 fQAmc->Fill("h_PionK0_P_B", 4, p[0]);
	 fQAmc->Fill("h_PionK0_P_B", 4, p[1]);
       }
     }
  }

  data[0] = iMass;
  data[1] = pt;
  data[2] = theta;
  data[3] = phi;
  //printf("-D: m: %f, pT: %f, theta: %f, phi: %f\n", invMass, mPt, theta, phi);
  fQA->Fill("hK0", data);
  

  if(iMass < cutMass[0] || iMass > cutMass[1]) return kFALSE;

  // all cuts passed
  
  // some MC stuff
  if(2 == fCurrentV0id){
    fQAmc->Fill("h_K0_p_S", iP);
    fQAmc->Fill("h_PionK0_P_S", 5, p[0]);
    fQAmc->Fill("h_PionK0_P_S", 5, p[1]);
  }
  else if (-2 != fCurrentV0id){
    fQAmc->Fill("h_K0_p_B", iP);
    fQAmc->Fill("h_PionK0_P_B", 5, p[0]);
    fQAmc->Fill("h_PionK0_P_B", 5, p[1]);
  }

  return kTRUE;
}
//________________________________________________________________
Bool_t AliHFEV0cuts::LambdaCuts(AliESDv0 *v0, Bool_t &isLambda ){
  //
  // Lambda cuts - decision on Lambda - AntiLambda is taken too
  //
  // discrimination between lambda and antilambda - correlation of the following variables necessary:
  // - momentum of the proton AND momentum of the pion (proton momentum is allways larger)
  // - mass of the mother particle

  if(!v0) return kFALSE;

  if(fMCEvent){
    if(4 == fCurrentV0id){
      fQAmc->Fill("h_Mass_L_as_G",  v0->GetEffMass(0, 0));
      fQAmc->Fill("h_Mass_L_as_K0",  v0->GetEffMass(2, 0));
    }
  }

  const Double_t cL0mass=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();  // PDG lambda mass

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Float_t mMass[2] = {-1., -1.};
  if(CheckSigns(v0)){
    pIndex = v0->GetPindex();
    nIndex = v0->GetNindex();
    mMass[0] = v0->GetEffMass(4, 2);
    mMass[1] = v0->GetEffMass(2, 4);
  }
  else{
    pIndex = v0->GetNindex();
    nIndex = v0->GetPindex();    
    mMass[0] = v0->GetEffMass(2, 4);
    mMass[1] = v0->GetEffMass(4, 2);
  }
 
  daughter[0] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  AliKFParticle *kfMother[2] = {0x0, 0x0};
  // Lambda
  kfMother[0] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kProton), TMath::Abs(kPiPlus));
  if(!kfMother[0]) return kFALSE;
  
  // production vertex is set in the 'CreateMotherParticle' function
  //kfMother[0]->SetMassConstraint(cL0mass, 0.);

  // Anti Lambda
  kfMother[1] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kProton));
  if(!kfMother[1]) return kFALSE;
  // production vertex is set in the 'CreateMotherParticle' function
  //kfMother[1]->SetMassConstraint(cL0mass, 0.);

  Float_t dMass[2] = {TMath::Abs(mMass[0] - cL0mass), TMath::Abs(mMass[1] - cL0mass)};
  
  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(pIndex));
  d[1] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(nIndex));
  if(!d[0] || !d[1])    return kFALSE;
  
  Float_t p[2] = {d[0]->GetP(), d[1]->GetP()}; 

  // check the 3 lambda - antilambda variables
  Int_t check[2] = {-1, -1};   // 0 : lambda, 1 : antilambda
  // 1) momentum of the daughter particles - proton is expected to have higher momentum than pion
  check[0] = (p[0] > p[1]) ? 0 : 1;
  // 2) mass of the mother particle
  check[1] = (dMass[0] < dMass[1]) ? 0 : 1;
  fQA->Fill("h_L_checks", check[0]*1.0, check[1]*1.0);
 
  // if the two check do not agree
  if(check[0] != check[1]){
    if(kfMother[0]) delete kfMother[0]; 
    if(kfMother[1]) delete kfMother[1]; 
    return kFALSE;
  }

  // now that the check[0] == check[1]
  const Int_t type = check[0];

  Float_t iMass =0.;
  if(CheckSigns(v0)){
    iMass = (type == 0) ? v0->GetEffMass(4, 2) : v0->GetEffMass(2, 4);
  }
  else{
    iMass = (type == 0) ? v0->GetEffMass(2, 4) : v0->GetEffMass(4, 2);
  }
  Float_t iP = v0->P();

   // Cuts
  const Double_t cutChi2NDF = 2.;              // ORG [5.]
  const Double_t cutCosPoint[2] = {0., 0.01};  // ORG [0., 0.02]
  const Double_t cutDCA[2] = {0., 0.2};        // ORG [0., 0.2]
  const Double_t cutProdVtxR[2] = {2., 40.};   // ORG [0., 24.]
  const Double_t cutMass[2] = {1.11, 1.12};   // ORG [1.11, 1.12]
  // cundidate cuts
  // relative daughter momentum versusu mother momentum

  // compute the cut values
  
  // cos pointing angle
  Double_t cosPoint = v0->GetV0CosineOfPointingAngle();
  cosPoint = TMath::ACos(cosPoint);

  // DCA between daughters
  Double_t dca = v0->GetDcaV0Daughters();
  
  // Production vertex
  Double_t x, y, z; 
  v0->GetXYZ(x,y,z);
  Double_t r = TMath::Sqrt(x*x + y*y);

  // proton - pion indices
  Int_t ix[2] = {0, 1};
  if(1 == type){
    ix[0] = 1;
    ix[1] = 0;
  }

  // proton - pion indices - based on MC truth
  // for background use the reconstructed indices
  Int_t ixMC[2] = {-1, -1}; // {proton, pion}
  if(fMCEvent){
    if(4 == fCurrentV0id){
      ixMC[0] = 0;
      ixMC[1] = 1;
    }
    else if(-4 == fCurrentV0id){
      ixMC[0] = 1;
      ixMC[1] = 0;
    }
    else{
      ixMC[0] = ix[0];
      ixMC[1] = ix[1];
    }
  }

  // V0 chi2/ndf
  Double_t chi2ndf = kfMother[type]->GetChi2()/kfMother[type]->GetNDF();

  if(kfMother[0]) delete kfMother[0]; 
  if(kfMother[1]) delete kfMother[1]; 

  // Relative daughter momentum
  Double_t rP = (0 == check[0]) ? p[1]/p[0] : p[0]/p[1];
  

  //
  // Apply the cuts, produce QA plots (with mass cut)
  //
  
  (type == 0) ?   fQA->Fill("h_L_Mass", 0, iMass) :  fQA->Fill("h_AL_Mass", 0, iMass);

 

  // MC
  if(fMCEvent){
    if(4 == fCurrentV0id){
      fQAmc->Fill("h_L_Mass_S", 0, iMass);
      fQAmc->Fill("h_lambda_MvP_S", iP, iMass);
    }
    else if(-4 == fCurrentV0id){
      fQAmc->Fill("h_AL_Mass_S", 0, iMass);
      fQAmc->Fill("h_lambda_MvP_S", iP, iMass);
    }
    else if(-2 != fCurrentV0id) fQAmc->Fill("h_LAL_Mass_B", 0, iMass);
  }


  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 0, p[ix[0]]);
    fQA->Fill("h_PionL_P", 0, p[ix[1]]);
    fQA->Fill("h_cut_L_Chi2", 0, chi2ndf);
    fQA->Fill("h_cut_L_Chi2", 1, chi2ndf);
    fQA->Fill("h_cut_L_CosPoint", 0, cosPoint);
    fQA->Fill("h_cut_L_DCA", 0, dca);
    fQA->Fill("h_cut_L_VtxR", 0, r);
    fQA->Fill("h_cut_L_rdp_v_mp", 0, iP, rP);
  }
  if(fMCEvent){
    if(iMass > cutMass[0] && iMass < cutMass[1]){
      if(4 == TMath::Abs(fCurrentV0id)){
	fQAmc->Fill("h_cut_L_Chi2_S", 0, iP, chi2ndf);
	fQAmc->Fill("h_cut_L_Chi2_S", 1, iP, chi2ndf);
	fQAmc->Fill("h_cut_L_CosPoint_S", 0, iP, cosPoint);
	fQAmc->Fill("h_cut_L_DCA_S", 0, iP, dca);
	fQAmc->Fill("h_cut_L_VtxR_S", 0, iP, r);
	fQAmc->Fill("h_cut_L_rdp_v_mp_S", 0, iP, rP);	
	fQAmc->Fill("h_ProtonL_P_S", 0, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_S", 0, p[ixMC[1]]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_L_Chi2_B", 0, iP, chi2ndf);
	fQAmc->Fill("h_cut_L_Chi2_B", 1, iP, chi2ndf);
	fQAmc->Fill("h_cut_L_CosPoint_B", 0, iP, cosPoint);
	fQAmc->Fill("h_cut_L_DCA_B", 0, iP, dca);
	fQAmc->Fill("h_cut_L_VtxR_B", 0, iP, r);
	fQAmc->Fill("h_cut_L_rdp_v_mp_B", 0, iP, rP);	
	fQAmc->Fill("h_ProtonL_P_B", 0, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_B", 0, p[ixMC[1]]);
      }
    }
  }
  //
  // Chi2/NDF cut
  //
  if(chi2ndf > cutChi2NDF) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 1, iMass) :  fQA->Fill("h_AL_Mass", 1, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_cut_L_CosPoint", 1, cosPoint);
    fQA->Fill("h_ProtonL_P", 1, p[ix[0]]);
    fQA->Fill("h_PionL_P", 1, p[ix[1]]);
  }
  if(fMCEvent){
    if(4 == fCurrentV0id) fQAmc->Fill("h_L_Mass_S", 1, iMass);
    else if(-4 == fCurrentV0id)  fQAmc->Fill("h_AL_Mass_S", 1, iMass);
    else if(-2 != fCurrentV0id)  fQAmc->Fill("h_LAL_Mass_B", 1, iMass);
    if(iMass > cutMass[0] && iMass < cutMass[1]){
      if(4 == TMath::Abs(fCurrentV0id)){
	fQAmc->Fill("h_cut_L_CosPoint_S", 1, iP, cosPoint);
	fQAmc->Fill("h_ProtonL_P_S", 1, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_S", 1, p[ixMC[1]]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_L_CosPoint_B", 1, iP, cosPoint);
	fQAmc->Fill("h_ProtonL_P_B", 1, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_B", 1, p[ixMC[1]]);
      }
    }
  }

  //
  // Cos point cut
  //
  if(cosPoint < cutCosPoint[0] || cosPoint > cutCosPoint[1]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 2, iMass) :  fQA->Fill("h_AL_Mass", 2, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 2, p[ix[0]]);
    fQA->Fill("h_PionL_P", 2, p[ix[1]]);
    fQA->Fill("h_cut_L_DCA", 1, dca);
  }
  if(fMCEvent){
    if(4 == fCurrentV0id) fQAmc->Fill("h_L_Mass_S", 2, iMass);
    else if(-4 == fCurrentV0id)  fQAmc->Fill("h_AL_Mass_S", 2, iMass);
    else if(-2 != fCurrentV0id)  fQAmc->Fill("h_LAL_Mass_B", 2, iMass);
    if(iMass > cutMass[0] && iMass < cutMass[1]){
      if(4 == TMath::Abs(fCurrentV0id)){
	fQAmc->Fill("h_cut_L_DCA_S", 1, iP, dca);
	fQAmc->Fill("h_ProtonL_P_S", 2, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_S", 2, p[ixMC[1]]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_L_DCA_B", 1, iP, dca);
	fQAmc->Fill("h_ProtonL_P_B", 2, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_B", 2, p[ixMC[1]]);
      }
    }
  }

  //
  // DCA cut
  //  
  if(dca < cutDCA[0] || dca > cutDCA[1]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 3, iMass) :  fQA->Fill("h_AL_Mass", 3, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 3, p[ix[0]]);
    fQA->Fill("h_PionL_P", 3, p[ix[1]]);
    fQA->Fill("h_cut_L_VtxR", 1, r);
  }
  if(fMCEvent){
    if(4 == fCurrentV0id) fQAmc->Fill("h_L_Mass_S", 3, iMass);
    else if(-4 == fCurrentV0id)  fQAmc->Fill("h_AL_Mass_S", 3, iMass);
    else if(-2 != fCurrentV0id)  fQAmc->Fill("h_LAL_Mass_B", 3, iMass);
    if(iMass > cutMass[0] && iMass < cutMass[1]){
      if(4 == TMath::Abs(fCurrentV0id)){
	fQAmc->Fill("h_cut_L_VtxR_S", 1, iP, r);
	fQAmc->Fill("h_ProtonL_P_S", 3, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_S", 3, p[ixMC[1]]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_cut_L_VtxR_B", 1, iP, r);
	fQAmc->Fill("h_ProtonL_P_B", 3, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_B", 3, p[ixMC[1]]);
      }
    }
  }

  //
  // Vertex radius cut
  //
  if(r < cutProdVtxR[0] || r > cutProdVtxR[1]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 4, iMass) :  fQA->Fill("h_AL_Mass", 4, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 4, p[ix[0]]);
    fQA->Fill("h_PionL_P", 4, p[ix[1]]);
  }
  if(fMCEvent){
    if(4 == fCurrentV0id) fQAmc->Fill("h_L_Mass_S", 4, iMass);
    else if(-4 == fCurrentV0id)  fQAmc->Fill("h_AL_Mass_S", 4, iMass);
    else if(-2 != fCurrentV0id)  fQAmc->Fill("h_LAL_Mass_B", 4, iMass);
    if(iMass > cutMass[0] && iMass < cutMass[1]){
      if(4 == TMath::Abs(fCurrentV0id)){
	fQAmc->Fill("h_ProtonL_P_S", 4, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_S", 4, p[ixMC[1]]);
      }
      else if(-2 != fCurrentV0id){
	fQAmc->Fill("h_ProtonL_P_B", 4, p[ixMC[0]]);
	fQAmc->Fill("h_PionL_P_B", 4, p[ixMC[1]]);
      }
    }
  }

  if(iMass < cutMass[0] || iMass > cutMass[1]) {
    return kFALSE;
  }

  // all cuts passed

  // assign the lambda type value: Lambda: kTRUE, Anti-Lambda: kFALSE
  isLambda = (0 == type) ? kTRUE : kFALSE;

  // some MC stuff
  if(4 == fCurrentV0id){
    fQAmc->Fill("h_lambda_p_S", iP);
  }
  else if(-4 == fCurrentV0id){
    fQAmc->Fill("h_alambda_p_S", iP);
  }
  else if (-2 != fCurrentV0id && 0 == type){
    fQAmc->Fill("h_lambda_p_B", iP);
  }
  else if(-2 != fCurrentV0id && 0 != type ){
    fQAmc->Fill("h_alambda_p_B", iP);
  }
  //
  if(fMCEvent){
    if(4 == TMath::Abs(fCurrentV0id)){
      fQAmc->Fill("h_ProtonL_P_S", 5, p[ixMC[0]]);
      fQAmc->Fill("h_PionL_P_S", 5, p[ixMC[1]]);
    }
    else if(-2 != fCurrentV0id){
      fQAmc->Fill("h_ProtonL_P_B", 5, p[ixMC[0]]);
      fQAmc->Fill("h_PionL_P_B", 5, p[ixMC[1]]);
    }
  }
  return kTRUE;
}
//________________________________________________________________
Double_t AliHFEV0cuts::OpenAngle(AliESDv0 const *v0){
  //
  // Opening angle between two daughter tracks
  //
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
    
  
  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

  
  Double_t openAngle = TMath::ACos((mp[0]*mn[0] + mp[1]*mn[1] + mp[2]*mn[2])/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1] + mp[2]*mp[2])*TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1] + mn[2]*mn[2])));
  
  return TMath::Abs(openAngle);
}
//________________________________________________________________
Double_t AliHFEV0cuts::PsiPair(const AliESDv0 *v0) {
  //
  // Angle between daughter momentum plane and plane 
  // 

  if(!fInputEvent) return -1.;

  Float_t magField = fInputEvent->GetMagneticField();

  Int_t pIndex = -1;
  Int_t nIndex = -1;
  if(CheckSigns(v0)){
    pIndex = v0->GetPindex();
    nIndex = v0->GetNindex();
  }
  else{
    pIndex = v0->GetNindex();
    nIndex = v0->GetPindex();    
  }
 

  AliESDtrack* daughter[2];

  daughter[0] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(nIndex));

  Double_t x, y, z;
  v0->GetXYZ(x,y,z);//Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  

  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter; 


  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) -  TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis

  Double_t radiussum = TMath::Sqrt(x*x + y*y) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3];
  Double_t momNegProp[3];
    
  AliExternalTrackParam pt(*daughter[0]), nt(*daughter[1]);
    
  Double_t psiPair = 4.;

  if(nt.PropagateTo(radiussum,magField) == 0)//propagate tracks to the outside
    psiPair =  -5.;
  if(pt.PropagateTo(radiussum,magField) == 0)
    psiPair = -5.;
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);
  
  Double_t pEle =
    TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos =
    TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
    
  Double_t scalarproduct =
    momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
    
  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));  

  return psiPair; 
}
//________________________________________________________________
AliKFParticle *AliHFEV0cuts::CreateMotherParticle(AliVTrack const *pdaughter, AliVTrack const *ndaughter, Int_t pspec, Int_t nspec){
  //
  // Creates a mother particle
  //
  AliKFParticle pkfdaughter(*pdaughter, pspec);
  AliKFParticle nkfdaughter(*ndaughter, nspec);
  
  // - check if the daughter particles are coming from the primary vertex
  // - check the number of tracks that take part in the creaton of primary vertex.
  //   important: after removeal of candidate tracks there must be at least 2 tracks left
  //   otherwise the primary vertex will be corrupted  
 
  
  // Create the mother particle 
  AliKFParticle *m = new AliKFParticle(pkfdaughter, nkfdaughter);
  // important !!!
  m->SetField(fInputEvent->GetMagneticField());
  if(TMath::Abs(kElectron) == pspec && TMath::Abs(kElectron) == nspec) m->SetMassConstraint(0, 0.001);
  else if(TMath::Abs(kPiPlus) == pspec && TMath::Abs(kPiPlus) == nspec) m->SetMassConstraint(TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(), 0.);
  else if(TMath::Abs(kProton) == pspec && TMath::Abs(kPiPlus) == nspec) m->SetMassConstraint(TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass(), 0.);
  else if(TMath::Abs(kPiPlus) == pspec && TMath::Abs(kProton) == nspec) m->SetMassConstraint(TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass(), 0.);
  else{
    AliError("Wrong daughter ID - mass constraint can not be set");
  }

  AliKFVertex improvedVertex = *fPrimaryVertex;
  improvedVertex += *m;
  m->SetProductionVertex(improvedVertex);
  
  // update 15/06/2010
  // mother particle will not be added to primary vertex but only to its copy 
  // as this confilcts with calling
  // m->SetPrimaryVertex() function and
  // subsequently removing the mother particle afterwards
  // Sourse: Sergey Gorbunov

  return m;
}
//___________________________________________________________________
void  AliHFEV0cuts::Armenteros(const AliESDv0 *v0, Float_t val[2]){
  //
  // computes the Armenteros variables for given V0
  // fills the histogram
  // returns the values via "val"
  //
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};  
  Double_t mm[3] = {0,0,0};  

  if(CheckSigns(v0)){
    v0->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
    v0->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
  }
  else{
    v0->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
    v0->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
  }
  v0->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother

  TVector3 vecN(mn[0],mn[1],mn[2]);
  TVector3 vecP(mp[0],mp[1],mp[2]);
  TVector3 vecM(mm[0],mm[1],mm[2]);
  
  Double_t thetaP = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
  Double_t thetaN = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
  
  Double_t alfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/
    ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) ;
  Double_t qt = vecP.Mag()*sin(thetaP);

  val[0] = alfa;
  val[1] = qt;

}
//___________________________________________________________________
Bool_t AliHFEV0cuts::CheckSigns(AliESDv0 const *v0){
  //
  // check wheter the sign was correctly applied to 
  // V0 daughter tracks
  //
  
  Bool_t correct = kFALSE;

  Int_t pIndex = 0, nIndex = 0;
  pIndex = v0->GetPindex();
  nIndex = v0->GetNindex();
  
  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(pIndex));
  d[1] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(nIndex));

  Int_t sign[2];
  sign[0] = (int)d[0]->GetSign();
  sign[1] = (int)d[1]->GetSign();
  
  if(-1 == sign[0] && 1 == sign[1]){
    correct = kFALSE;
    //v0->SetIndex(0, pIndex);  // set the index of the negative v0 track
    //v0->SetIndex(1, nIndex);  // set the index of the positive v0 track
  }
  else{
    correct = kTRUE;
  }
  
  //pIndex = v0->GetPindex();
  //nIndex = v0->GetNindex();
  //printf("-D2: P: %i, N: %i\n", pIndex, nIndex);

  return correct;
}
//___________________________________________________________________
Bool_t AliHFEV0cuts::GetConvPosXY(AliESDtrack * const ptrack, AliESDtrack * const ntrack, Double_t convpos[2]){
  //
  // recalculate the gamma conversion XY postition
  //

  const Double_t b = fInputEvent->GetMagneticField();

  Double_t helixcenterpos[2];
  GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

  Double_t helixcenterneg[2];
  GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

  Double_t  poshelix[6];
  ptrack->GetHelixParameters(poshelix,b);
  Double_t posradius = TMath::Abs(1./poshelix[4]);

  Double_t  neghelix[6];
  ntrack->GetHelixParameters(neghelix,b);
  Double_t negradius = TMath::Abs(1./neghelix[4]);

  Double_t xpos = helixcenterpos[0];
  Double_t ypos = helixcenterpos[1];
  Double_t xneg = helixcenterneg[0];
  Double_t yneg = helixcenterneg[1];

  convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
  convpos[1] = (ypos*negradius+  yneg*posradius)/(negradius+posradius);

  return 1;
}
//___________________________________________________________________
Bool_t AliHFEV0cuts::GetHelixCenter(AliESDtrack * const track, Double_t b,Int_t charge, Double_t center[2]){
  // see header file for documentation
  
  Double_t pi = TMath::Pi();
  
  Double_t  helix[6];
  track->GetHelixParameters(helix,b);
  
  Double_t xpos =  helix[5];
  Double_t ypos =  helix[0];
  Double_t radius = TMath::Abs(1./helix[4]);
  Double_t phi = helix[2];

  if(phi < 0){
    phi = phi + 2*pi;
  }

  phi -= pi/2.;
  Double_t xpoint =  radius * TMath::Cos(phi);
  Double_t ypoint =  radius * TMath::Sin(phi);

  if(b<0){
    if(charge > 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
  }
  if(b>0){
    if(charge < 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
  }
  center[0] =  xpos + xpoint;
  center[1] =  ypos + ypoint;

  return 1;
}
