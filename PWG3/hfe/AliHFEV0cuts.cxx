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

#include "AliHFEcollection.h"

#include "AliHFEV0cuts.h"

ClassImp(AliHFEV0cuts)

//________________________________________________________________
AliHFEV0cuts::AliHFEV0cuts():
  fQA(NULL)
  , fMCEvent(NULL)
  , fInputEvent(NULL)
  , fPrimaryVertex(NULL)
  , fIsMC(kFALSE)
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
  if (fQA) delete fQA;
}

//________________________________________________________________
AliHFEV0cuts::AliHFEV0cuts(const AliHFEV0cuts &ref):
  TObject(ref)
  , fQA(NULL)
  , fMCEvent(NULL)
  , fInputEvent(NULL)
  , fPrimaryVertex(NULL)
  , fIsMC(kFALSE)
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

  fQA = new AliHFEcollection("fQA", name);


  // common for all V0s
  fQA->CreateTH2Fvector1(2, "h_all_AP", "armenteros plot for all V0 candidates", 200, -1, 1, 200, 0, 0.25);

  // gammas
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_CosPoint", "Gamma Cosine pointing angle; cos point. angle; counts", 100, 0, 0.1);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_DCA", "DCA between the gamma daughters; dca (cm); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_VtxR", "Radius of the gamma conversion vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_OA", "opening angle of the gamma products; opening angle (rad); counts", 100, 0, 1);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_PP", "gamma psi pair angle; psi pairangle (rad); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_Gamma_Chi2", "gamma Chi2/NDF; Chi2/NDF; counts", 100, 0, 25);
  fQA->CreateTH1Fvector1(9, "h_Gamma_Mass", "Invariant mass of gammas; mass (GeV/c^{2}); counts", 100, 0, 0.2);

 
  // kaons
  fQA->CreateTH1Fvector1(2, "h_cut_K0_CosPoint", "K0 Cosine pointing angle; cos point. angle; counts", 100, 0, 0.1);
  fQA->CreateTH1Fvector1(2, "h_cut_K0_DCA", "DCA between the K0 daughters; dca (cm); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_K0_VtxR", "Radius of the K0 decay vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_K0_Chi2", "K0 Chi2/NDF; Chi2/NDF; counts", 100, 0, 25);
  fQA->CreateTH2Fvector1(2, "h_cut_K0_AP", "Armenteros plot for K0 candidates; #alpha; q_{T} (GeV/c)", 100, -1, 1, 100, 0, 0.3);
  fQA->CreateTH1Fvector1(8, "h_K0_Mass", "Invariant mass of K0; mass (GeV/c^{2}); counts", 125, 0.45, 0.55);

  // lambda
  fQA->CreateTH1Fvector1(2, "h_cut_L_CosPoint", "L Cosine pointing angle; cos point. angle; counts", 100, 0, 0.1);
  fQA->CreateTH1Fvector1(2, "h_cut_L_DCA", "DCA between the L daughters; dca (cm); counts", 100, 0, 2);
  fQA->CreateTH1Fvector1(2, "h_cut_L_VtxR", "Radius of the L decay vertex; r (cm); counts", 1000, 0, 100);
  fQA->CreateTH1Fvector1(2, "h_cut_L_Chi2", "L Chi2/NDF; Chi2/NDF; counts", 100, 0, 25);
  fQA->CreateTH2Fvector1(2, "h_cut_L_AP", "Armenteros plot for Lambda candidates; #alpha; q_{T} (GeV/c)", 100, -1, 1, 100, 0, 0.3);
  fQA->CreateTH2Fvector1(2, "h_cut_AL_AP", "Armenteros plot for anti Lambda candidates; #alpha; q_{T} (GeV/c)", 100, -1, 1, 100, 0, 0.3);
  fQA->CreateTH2Fvector1(2, "h_cut_Gamma_AP", "Armenteros plot for gamma candidates; #alpha; q_{T} (GeV/c)", 100, -1, 1, 100, 0, 0.3);
  fQA->CreateTH1Fvector1(9, "h_L_Mass", "Invariant mass of L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);
  fQA->CreateTH1Fvector1(9, "h_AL_Mass", "Invariant mass of anti L; mass (GeV/c^{2}); counts", 60, 1.1, 1.13);

  fQA->CreateTH2F("h_L_checks", "Lambda candidate check[0] -v- check[1]; check[0]; check[1]", 5, -0.75, 1.75, 6, -0.75, 1.75 );
  
  // electrons
  fQA->CreateTH1Fvector1(9, "h_Electron_P", "Momenta of conversion electrons -cuts-; P (GeV/c); counts", 50, 0.1, 20, 0);

  // K0 pions
  fQA->CreateTH1Fvector1(8, "h_PionK0_P", "Momenta of K0 pions -cuts-; P (GeV/c) counts;", 50, 0.1, 20, 0);
  
  // L pions
  fQA->CreateTH1Fvector1(9, "h_PionL_P", "Momenta of L pions -cuts-; P (GeV/c) counts;", 50, 0.1, 50, 0);
  
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
  
  // Gamma
  fQA->CreateTH2Fvector1(2, "h_cut_Gamma_OAvP", "open. ang. of the Gamma daughters versus Gamma mom; Gamma p (GeV/c); opening angle (pions) (rad)", 100, 0.1, 10, 200, 0., 0.2);
  // K0
  fQA->CreateTH2Fvector1(2, "h_cut_K0_OAvP", "open. ang. of the K0 daughters versus K0 momentum; K0 p (GeV/c); opening angle (pions) (rad)", 100, 0.1, 10, 100, 0, 3.5);
  // Lambda
  fQA->CreateTH2Fvector1(2, "h_cut_L_OAvP", "open. ang. of the L daughters versus L momentum; Lambda p (GeV/c); openeing angle pion-proton (rad)", 100, 0.1, 10, 100, 0, 3.5);
  fQA->CreateTH2Fvector1(2, "h_cut_L_rdp_v_mp", "relative L daughter mom -v- mother mom; L mom (GeV/c); relative daughter mom p2/p1", 100, 0.1, 10, 100, 0, 1);
  fQA->CreateTH2Fvector1(2, "h_cut_L_qT_v_mp", "A-P q_{T} -v- Lambda momentum; mom (GeV/c); q_{T} GeV/c", 100, 0.1, 10, 50, 0., 0.12);


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
  if(track->GetTPCNcls() < 80) return kFALSE;

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
  if(chi2perTPCcluster > 3.5) return kFALSE;

  // TPC cluster ratio
  Float_t cRatioTPC = track->GetTPCNclsF() > 0. ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t> (track->GetTPCNclsF()) : 1.;
  fQA->Fill("h_ST_TPCclsR", cRatioTPC);
  if(cRatioTPC < 0.6) return kFALSE;

  // kinks
  fQA->Fill("h_ST_kinks", track->GetKinkIndex(0));
  if(track->GetKinkIndex(0) != 0) return kFALSE;

  // pt
  fQA->Fill("h_ST_pt",track->Pt());
  if(track->Pt() < 0.1 || track->Pt() > 100) return kFALSE;

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

  // loose cuts first
  if(LooseRejectK0(v0) || LooseRejectLambda(v0)) return kFALSE;

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
  kfMother->SetMassConstraint(0, 0.001);

  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(pIndex));
  d[1] = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(nIndex));

  Float_t iMass = v0->GetEffMass(0, 0);
  Float_t iP = v0->P();
  Float_t p[2] = {d[0]->GetP(), d[1]->GetP()};

  // Armenteros
  Float_t ap[2];
  Armenteros(v0, ap);

  // Cut values
  const Double_t cutCosPoint[2] = {0., 0.03};  // ORG [0., 0.03]
  const Double_t cutDCA[2] = {0., 0.25};       // ORG [0., 0.25]
  const Double_t cutProdVtxR[2] = {6., 50.};   // ORG [6., 9999]
  const Double_t cutOAngle[2] = {0, 0.1};      // ORG [0., 0.1]
  const Double_t cutPsiPair[2] = {0., 0.05};   // ORG [0. 0.05]
  const Double_t cutMass = 0.05;               // ORG [0.05]
  const Double_t cutChi2NDF = 7.;              // ORG [7.]
  // armenteros cuts
  const Double_t cutAlpha[2] = {0.35, 0.45};   // [0.35, 0.45]
  const Double_t cutQT = 0.015;
  
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

  // Opening angle
  Double_t oAngle = OpenAngle(v0);

  // psi pair 
  Double_t psiPair = PsiPair(v0);
  
  // V0 chi2/ndf
  Double_t chi2ndf = kfMother->GetChi2()/kfMother->GetNDF();

  if(kfMother) delete kfMother; 
 
  //
  // Apply the cuts, produce QA plots (with mass cut)
  //
  fQA->Fill("h_Gamma_Mass", 0, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_all_AP", 0, ap[0], ap[1]);
    fQA->Fill("h_Electron_P", 0, p[0]);
    fQA->Fill("h_Electron_P", 0, p[1]);
    fQA->Fill("h_cut_Gamma_CosPoint", 0, cosPoint);
    fQA->Fill("h_cut_Gamma_CosPoint", 1, cosPoint);
    fQA->Fill("h_cut_Gamma_DCA", 0, dca);
    fQA->Fill("h_cut_Gamma_VtxR", 0, r);
    fQA->Fill("h_cut_Gamma_OA", 0, oAngle);
    fQA->Fill("h_cut_Gamma_PP", 0, psiPair);
    fQA->Fill("h_cut_Gamma_Chi2", 0, chi2ndf);
    fQA->Fill("h_cut_Gamma_OAvP", 0, iP, oAngle);	
    fQA->Fill("h_cut_Gamma_AP", 0, ap[0], ap[1]);
 }

  if(cosPoint < cutCosPoint[0] || cosPoint > cutCosPoint[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 1, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 1, p[0]);
    fQA->Fill("h_Electron_P", 1, p[1]);
    fQA->Fill("h_cut_Gamma_DCA", 1, dca);
  }
  
  if(dca < cutDCA[0] || dca > cutDCA[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 2, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 2, p[0]);
    fQA->Fill("h_Electron_P", 2, p[1]);
    fQA->Fill("h_cut_Gamma_VtxR", 1, r);
  }

  if(r < cutProdVtxR[0] || r > cutProdVtxR[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 3, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 3, p[0]);
    fQA->Fill("h_Electron_P", 3, p[1]);
    fQA->Fill("h_cut_Gamma_OA", 1, oAngle);
  }

  if(oAngle < cutOAngle[0] || oAngle > cutOAngle[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 4, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 4, p[0]);
    fQA->Fill("h_Electron_P", 4, p[1]);
    fQA->Fill("h_cut_Gamma_PP", 1, psiPair);
  }

  if(psiPair < cutPsiPair[0] || psiPair > cutPsiPair[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 5, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 5, p[0]);
    fQA->Fill("h_Electron_P", 5, p[1]);
    fQA->Fill("h_cut_Gamma_Chi2", 1, chi2ndf);
  }

  if(chi2ndf > cutChi2NDF) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 6, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 6, p[0]);
    fQA->Fill("h_Electron_P", 6, p[1]);
    fQA->Fill("h_cut_Gamma_OAvP", 1, iP, oAngle);	
    fQA->Fill("h_all_AP", 1, ap[0], ap[1]);
    fQA->Fill("h_cut_Gamma_AP", 1, ap[0], ap[1]);
  }

  if(TMath::Abs(ap[0]) > cutAlpha[0] && TMath::Abs(ap[0]) < cutAlpha[1]) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 7, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 7, p[0]);
    fQA->Fill("h_Electron_P", 7, p[1]);
  }

  if(ap[1] > cutQT) return kFALSE;
  fQA->Fill("h_Gamma_Mass", 8, iMass);
  if(iMass < cutMass){
    fQA->Fill("h_Electron_P", 8, p[0]);
    fQA->Fill("h_Electron_P", 8, p[1]);
  }
  

  if(iMass > cutMass) return kFALSE;

  // all cuts passed
    
  return kTRUE;
}
//________________________________________________________________
Bool_t AliHFEV0cuts::K0Cuts(AliESDv0 *v0){
  //
  // K0 cuts
  //

  if(!v0) return kFALSE;

  const Double_t cK0mass=TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();  // PDG K0s mass
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
  kfMother->SetMassConstraint(cK0mass, 0.);
  
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
  const Double_t cutCosPoint[2] = {0., 0.03};  // ORG [0., 0.03]
  const Double_t cutDCA[2] = {0., 0.2};        // ORG [0., 0.1]
  const Double_t cutProdVtxR[2] = {2.0, 30.};   // ORG [0., 8.1]
  const Double_t cutMass[2] = {0.49, 0.51};   // ORG [0.485, 0.51]
  const Double_t cutChi2NDF = 5.;              // ORG [7.]
  const Double_t cutOAngleP = (1.0/(iP + 0.3) - 0.1); // momentum dependent min. OAngle ~ 1/x
  // cundidate cuts
  // armenteros plot
  const Double_t cutQT = 0.1075;
  // elipse cut - see bellow

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

  // Opening angle
  Double_t oAngle = OpenAngle(v0);
  
  // Armenteros
  Float_t ap[2];
  Armenteros(v0, ap);
  const Double_t cutAP = 0.22 * TMath::Sqrt( TMath::Abs( (1-ap[0]*ap[0]/(0.92*0.92)) ) );
  

  //
  // Apply the cuts, produce QA plots (with mass cut)
  //

  fQA->Fill("h_K0_Mass", 0, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_all_AP", 0, ap[0], ap[1]);
    fQA->Fill("h_PionK0_P", 0, p[0]);
    fQA->Fill("h_PionK0_P", 0, p[1]);
    fQA->Fill("h_cut_K0_CosPoint", 0, cosPoint);
    fQA->Fill("h_cut_K0_CosPoint", 1, cosPoint);
    fQA->Fill("h_cut_K0_DCA", 0, dca);
    fQA->Fill("h_cut_K0_VtxR", 0, r);
    fQA->Fill("h_cut_K0_Chi2", 0, chi2ndf);
    fQA->Fill("h_cut_K0_OAvP", 0, iP, oAngle);
    fQA->Fill("h_cut_K0_AP", 0, ap[0], ap[1]);
  }

  if(cosPoint < cutCosPoint[0] || cosPoint > cutCosPoint[1]) return kFALSE;
  fQA->Fill("h_K0_Mass", 1, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 1, p[0]);
    fQA->Fill("h_PionK0_P", 1, p[1]);
    fQA->Fill("h_cut_K0_DCA", 1, dca);
  }
  
  if(dca < cutDCA[0] || dca > cutDCA[1]) return kFALSE;
  fQA->Fill("h_K0_Mass", 2, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 2, p[0]);
    fQA->Fill("h_PionK0_P", 2, p[1]);
    fQA->Fill("h_cut_K0_VtxR", 1, r);
  }

  
  if(r < cutProdVtxR[0] || r > cutProdVtxR[1]) return kFALSE;
  fQA->Fill("h_K0_Mass", 3, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 3, p[0]);
    fQA->Fill("h_PionK0_P", 3, p[1]);
    fQA->Fill("h_cut_K0_Chi2", 1, chi2ndf);
  }

  if(chi2ndf > cutChi2NDF) return kFALSE;
  fQA->Fill("h_K0_Mass", 4, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 4, p[0]);
    fQA->Fill("h_PionK0_P", 4, p[1]);
    fQA->Fill("h_cut_K0_OAvP", 1, iP, oAngle);
  }  

  if(oAngle < cutOAngleP) return kFALSE;
  fQA->Fill("h_K0_Mass", 5, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 5, p[0]);
    fQA->Fill("h_PionK0_P", 5, p[1]);
    fQA->Fill("h_cut_K0_AP", 1, ap[0], ap[1]);
    fQA->Fill("h_all_AP", 1, ap[0], ap[1]);
  }

  if(ap[1] < cutQT) return kFALSE;
  fQA->Fill("h_K0_Mass", 6, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 6, p[0]);
    fQA->Fill("h_PionK0_P", 6, p[1]);
  }
    
  if(ap[1] > cutAP) return kFALSE;
  fQA->Fill("h_K0_Mass", 7, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_PionK0_P", 7, p[0]);
    fQA->Fill("h_PionK0_P", 7, p[1]);
  }

  data[0] = iMass;
  data[1] = pt;
  data[2] = theta;
  data[3] = phi;
  //printf("-D: m: %f, pT: %f, theta: %f, phi: %f\n", invMass, mPt, theta, phi);
  fQA->Fill("hK0", data);
  

  if(iMass < cutMass[0] || iMass > cutMass[1]) return kFALSE;

  // all cuts passed
  
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

  // loose cuts first
  if(LooseRejectK0(v0) || LooseRejectGamma(v0)) return kFALSE;
  
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
  kfMother[0]->SetMassConstraint(cL0mass, 0.);

  // Anti Lambda
  kfMother[1] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kProton));
  if(!kfMother[1]) return kFALSE;
  // production vertex is set in the 'CreateMotherParticle' function
  kfMother[1]->SetMassConstraint(cL0mass, 0.);

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
  const Double_t cutCosPoint[2] = {0., 0.03};  // ORG [0., 0.03]
  const Double_t cutDCA[2] = {0., 0.2};        // ORG [0., 0.2]
  const Double_t cutProdVtxR[2] = {2., 30.};   // ORG [0., 24.]
  const Double_t cutMass[2] = {1.11, 1.12};   // ORG [1.11, 1.12]
  const Double_t cutChi2NDF = 5.;              // ORG [5.]
  // cundidate cuts
  // opening angle as a function of L momentum
  const Double_t cutOAngleP = 0.3 - 0.2*iP; // momentum dependent min. OAngle linear cut
  // relative daughter momentum versusu mother momentum
  // armenteros plot cuts
  const Double_t cutQT = 0.03;
  const Double_t cutAlpha = 0.7;  // VERY strong - should supress the overlap with K0
  // next cut see below

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

  Int_t ix[2] = {0, 1};
  if(1 == type){
    ix[0] = 1;
    ix[1] = 0;
  }
  
  // V0 chi2/ndf
  Double_t chi2ndf = kfMother[type]->GetChi2()/kfMother[type]->GetNDF();

  if(kfMother[0]) delete kfMother[0]; 
  if(kfMother[1]) delete kfMother[1]; 

  // Opening angle
  Double_t oAngle = OpenAngle(v0);

  // Relative daughter momentum
  Double_t rP = (0 == check[0]) ? p[1]/p[0] : p[0]/p[1];
  
  // Armenteros
  Float_t ap[2];
  Armenteros(v0, ap);
  
  Double_t cutAP[2]; // a bit of woodoo :-)
  cutAP[0] = 1.0 - (ap[0]-0.7 * ap[0]-0.7)*1.1 - 0.87;
  cutAP[1] = 1.0 - (ap[0]+0.7 * ap[0]+0.7)*1.1 - 0.87;


  //
  // Apply the cuts, produce QA plots (with mass cut)
  //
  
  (type == 0) ?   fQA->Fill("h_L_Mass", 0, iMass) :  fQA->Fill("h_AL_Mass", 0, iMass);

  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_all_AP", 0, ap[0], ap[1]);
    fQA->Fill("h_ProtonL_P", 0, p[ix[0]]);
    fQA->Fill("h_PionL_P", 0, p[ix[1]]);
    fQA->Fill("h_cut_L_Chi2", 0, chi2ndf);
    fQA->Fill("h_cut_L_CosPoint", 0, cosPoint);
    fQA->Fill("h_cut_L_CosPoint", 1, cosPoint);
    fQA->Fill("h_cut_L_DCA", 0, dca);
    fQA->Fill("h_cut_L_VtxR", 0, r);
    fQA->Fill("h_cut_L_OAvP", 0, iP, oAngle);
    fQA->Fill("h_cut_L_rdp_v_mp", 0, iP, rP);
    if(0 ==type)  fQA->Fill("h_cut_L_AP", 0, ap[0], ap[1]);
    else fQA->Fill("h_cut_AL_AP", 0, ap[0], ap[1]);
    fQA->Fill("h_cut_L_qT_v_mp", 0, iP, ap[1]);
  }

  if(cosPoint < cutCosPoint[0] || cosPoint > cutCosPoint[1]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 1, iMass) :  fQA->Fill("h_AL_Mass", 1, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 1, p[ix[0]]);
    fQA->Fill("h_PionL_P", 1, p[ix[1]]);
    fQA->Fill("h_cut_L_DCA", 1, dca);
  }
  
  if(dca < cutDCA[0] || dca > cutDCA[1]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 2, iMass) :  fQA->Fill("h_AL_Mass", 2, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 2, p[ix[0]]);
    fQA->Fill("h_PionL_P", 2, p[ix[1]]);
     fQA->Fill("h_cut_L_VtxR", 1, r);
  }

  if(r < cutProdVtxR[0] || r > cutProdVtxR[1]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 3, iMass) :  fQA->Fill("h_AL_Mass", 3, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 3, p[ix[0]]);
    fQA->Fill("h_PionL_P", 3, p[ix[1]]);
    fQA->Fill("h_cut_L_Chi2", 1, chi2ndf);
  }


  if(chi2ndf > cutChi2NDF) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 4, iMass) :  fQA->Fill("h_AL_Mass", 4, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 4, p[ix[0]]);
    fQA->Fill("h_PionL_P", 4, p[ix[1]]);
    fQA->Fill("h_cut_L_OAvP", 1, iP, oAngle);
  }

  if(oAngle < cutOAngleP) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 5, iMass) :  fQA->Fill("h_AL_Mass", 5, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 5, p[ix[0]]);
    fQA->Fill("h_PionL_P", 5, p[ix[1]]);
    fQA->Fill("h_cut_L_rdp_v_mp", 1, iP, rP);
    if(0 == type)  fQA->Fill("h_cut_L_AP", 1, ap[0], ap[1]);
    else fQA->Fill("h_cut_AL_AP", 1, ap[0], ap[1]);
    fQA->Fill("h_cut_L_qT_v_mp", 1, iP, ap[1]);
    fQA->Fill("h_all_AP", 1, ap[0], ap[1]);
  }
  
  if(ap[1] < cutQT) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 6, iMass) :  fQA->Fill("h_AL_Mass", 6, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 6, p[ix[0]]);
    fQA->Fill("h_PionL_P", 6, p[ix[1]]);
  }

  if(ap[1] > cutAP[type]) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 7, iMass) :  fQA->Fill("h_AL_Mass", 7, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 7, p[ix[0]]);
    fQA->Fill("h_PionL_P", 7, p[ix[1]]);
  }
  if(TMath::Abs(ap[0]) > cutAlpha) return kFALSE;
  (type == 0) ?   fQA->Fill("h_L_Mass", 8, iMass) :  fQA->Fill("h_AL_Mass", 8, iMass);
  if(iMass > cutMass[0] && iMass < cutMass[1]){
    fQA->Fill("h_ProtonL_P", 8, p[ix[0]]);
    fQA->Fill("h_PionL_P", 8, p[ix[1]]);
  }

  if(iMass < cutMass[0] || iMass > cutMass[1]) {
    return kFALSE;
  }

  // all cuts passed

  // assign the lambda type value: Lambda: kTRUE, Anti-Lambda: kFALSE
  isLambda = (0 == type) ? kTRUE : kFALSE;


  return kTRUE;
}
//________________________________________________________________
Double_t AliHFEV0cuts::OpenAngle(AliESDv0 *v0) const {
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
Double_t AliHFEV0cuts::PsiPair(AliESDv0 *v0) {
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
AliKFParticle *AliHFEV0cuts::CreateMotherParticle(AliVTrack* const pdaughter, AliVTrack* const ndaughter, Int_t pspec, Int_t nspec){
  //
  // Creates a mother particle
  //
  AliKFParticle pkfdaughter(*pdaughter, pspec);
  AliKFParticle nkfdaughter(*ndaughter, nspec);
  
  // - check if the daughter particles are coming from the primary vertex
  // - check the number of tracks that take part in the creaton of primary vertex.
  //   important: after removeal of candidate tracks there must be at least 2 tracks left
  //   otherwise the primary vertex will be corrupted  
  
  // ESD Analyis
  //const AliESDVertex *esdvertex = dynamic_cast<const AliESDVertex *>(fInputEvent->GetPrimaryVertex());
  //if(!esdvertex) return NULL;
  //UShort_t *contrib = esdvertex->GetIndices();
  
  //
  // not using the removal of the daughter track now
  //
  //   Int_t nTracks = esdvertex->GetNIndices();
  //   printf(" -D: N Vertex tracks: %i\n", nTracks);
  //   printf(" -D: N Contributors: %i\n", fPrimaryVertex->GetNContributors()); 
  //   Int_t nfound = 0;
  //   for(Int_t id = 0; id < esdvertex->GetNIndices(); id++){
  //     if(contrib[id] == pdaughter->GetID()){
  //       if( (nTracks - nfound) <= 2 ) return NULL;
  //       *fPrimaryVertex -= pkfdaughter;
  //       removed[0] = kTRUE;
  //       nfound++;
  //     }
  //     if(contrib[id] == ndaughter->GetID()){
  //       if( (nTracks - nfound) <=2 ) return NULL;
  //       *fPrimaryVertex -= nkfdaughter;
  //       removed[1] = kTRUE;
  //       nfound++;
  //     }
  //     if(nfound == 2) break;
  //   }
  
  //  printf(" -D: n removed: %i\n", nfound);
  
  // Create the mother particle 
  AliKFParticle *m = new AliKFParticle(pkfdaughter, nkfdaughter);

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
//_________________________________________________
Bool_t AliHFEV0cuts::LooseRejectK0(AliESDv0 * const v0) const {
  //
  // Reject K0 based on loose cuts
  //
  Double_t mass = v0->GetEffMass(AliPID::kPion, AliPID::kPion);
  if(mass > 0.494 && mass < 0.501) return kTRUE;
  return kFALSE;
}

//_________________________________________________
Bool_t AliHFEV0cuts::LooseRejectLambda(AliESDv0 * const v0) const {
  //
  // Reject Lambda based on loose cuts
  //
  Double_t mass1 = v0->GetEffMass(AliPID::kPion, AliPID::kProton);
  Double_t mass2 = v0->GetEffMass(AliPID::kProton, AliPID::kPion);
  
  if(mass1 > 1.1 && mass1 < 1.12) return kTRUE;
  if(mass2 > 1.1 && mass2 < 1.12) return kTRUE;
  return kFALSE;
}

//_________________________________________________
Bool_t AliHFEV0cuts::LooseRejectGamma(AliESDv0 * const v0) const {
  //
  // Reject Lambda based on loose cuts
  //
 
  Double_t mass = v0->GetEffMass(AliPID::kElectron, AliPID::kElectron);
  
  if(mass < 0.02) return kTRUE;
  return kFALSE;
}
//___________________________________________________________________
void  AliHFEV0cuts::Armenteros(AliESDv0 *v0, Float_t val[2]){
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
Bool_t AliHFEV0cuts::CheckSigns(AliESDv0* const v0){
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
