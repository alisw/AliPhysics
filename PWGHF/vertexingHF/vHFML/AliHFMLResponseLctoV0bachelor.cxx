// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseLctoV0bachelor
// \brief helper class to handle application of ML models for LctopK0s analyses trained
// with python libraries
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseLctoV0bachelor.h"
#include "AliAODRecoCascadeHF.h"
#include "AliVertexingHFUtils.h"
#include "AliPIDCombined.h"
#ifndef HomogeneousField
#define HomogeneousField 
#endif

/// \cond CLASSIMP
ClassImp(AliHFMLResponseLctoV0bachelor);
/// \endcond

//--------------------------------------------------------------------------
AliHFMLResponseLctoV0bachelor::AliHFMLResponseLctoV0bachelor() : AliHFMLResponse()
{
  //
  // Default constructor
  //
}

//--------------------------------------------------------------------------
AliHFMLResponseLctoV0bachelor::AliHFMLResponseLctoV0bachelor(const Char_t *name, const Char_t *title,
                                                             const std::string configfilepath) : AliHFMLResponse(name, title, configfilepath),
                                                                                                 fIsMixedEvent(false)
{
  //
  // Standard constructor
  //
}

//--------------------------------------------------------------------------
AliHFMLResponseLctoV0bachelor::~AliHFMLResponseLctoV0bachelor()
{
  //
  // Destructor
  //
}

//--------------------------------------------------------------------------
AliHFMLResponseLctoV0bachelor::AliHFMLResponseLctoV0bachelor(const AliHFMLResponseLctoV0bachelor &source) : AliHFMLResponse(source),
                                                                                                            fIsMixedEvent(source.fIsMixedEvent)
{
  //
  // Copy constructor
  //
}

//--------------------------------------------------------------------------
AliHFMLResponseLctoV0bachelor &AliHFMLResponseLctoV0bachelor::operator=(const AliHFMLResponseLctoV0bachelor &source)
{
  //
  // assignment operator
  //
  if (&source == this)
    return *this;
  
  AliHFMLResponse::operator=(source);
  
  return *this;
}

//--------------------------------------------------------------------------
void AliHFMLResponseLctoV0bachelor::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/)
{
  fVars["pt_cand"] = cand->Pt();
  fVars["d_len"] = cand->DecayLength();
  fVars["d_len_xy"] = cand->DecayLengthXY();
  fVars["norm_dl"] = cand->NormalizedDecayLength();
  fVars["norm_dl_xy"] = cand->NormalizedDecayLengthXY();
  fVars["cos_p"] = cand->CosPointingAngle();
  fVars["cos_p_xy"] = cand->CosPointingAngleXY();
  fVars["imp_par_xy"] = cand->ImpParXY();
  fVars["dca"] = cand->GetDCA();
  fVars["inv_mass"] = dynamic_cast<AliAODRecoCascadeHF *>(cand)->InvMassLctoK0sP();
  
  AliAODv0 * v0part = ((AliAODRecoCascadeHF*)cand)->Getv0();
  fVars["imp_par_prong0"] = cand->Getd0Prong(0);
  fVars["imp_par_prong1"] = v0part->Getd0Prong(0);
  fVars["imp_par_prong2"] = v0part->Getd0Prong(1);

  fVars["inv_mass_K0s"] = v0part->MassK0Short();
  fVars["dca_K0s"] = v0part->GetDCA();
  fVars["imp_par_K0s"] = cand->Getd0Prong(1);
  fVars["d_len_K0s"] = ((AliAODRecoCascadeHF*)cand)->DecayLengthV0();
  fVars["DecayLxy_Ks0"] = ((AliAODRecoCascadeHF*)cand)->DecayLengthXYV0();
  fVars["ctau_K0s"] = ((AliAODRecoCascadeHF*)cand)->DecayLengthV0()*0.497/(v0part->P());
  fVars["cos_p_K0s"] = ((AliAODRecoCascadeHF*)cand)->CosV0PointingAngle();
  fVars["pt_K0s"] = v0part->Pt();
  
  // Cosine of proton emission angle (theta*) in the rest frame of the mother particle
  // (from AliRDHFCutsLctoV0)
  TLorentzVector vpr, vk0s,vlc;
  Double_t massK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();    // mass K0S PDG
  Double_t massPrPDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();    // mass Proton PDG
  vpr.SetXYZM(cand->PxProng(0), cand->PyProng(0), cand->PzProng(0), massPrPDG);
  vk0s.SetXYZM(cand->PxProng(1), cand->PyProng(1), cand->PzProng(1), massK0SPDG);
  vlc = vpr + vk0s;
  TVector3 vboost = vlc.BoostVector();
  vpr.Boost(-vboost);
  fVars["cos_t_star"] = TMath::Cos(vpr.Angle(vlc.Vect()));

  AliAODTrack *bachelor = (AliAODTrack*)((AliAODRecoCascadeHF*)cand)->GetBachelor();
  fVars["pt_Pr"] = bachelor->Pt();
  AliAODVertex *primvert = dynamic_cast<AliAODVertex*>(cand->GetPrimaryVtx());
  if(!fIsMixedEvent){
    // Sign of d0 proton (different from regular d0)
    // (from AliRDHFCutsLctoV0)
    Double_t d0z0bach[2], covd0z0bach[3];
    bachelor->PropagateToDCA(primvert, bfield, kVeryBig, d0z0bach, covd0z0bach); // HOW DO WE SET THE B FIELD?; kVeryBig should come from AliExternalTrackParam
    Double_t tx[3];
    bachelor->GetXYZ(tx);
    tx[0] -= primvert->GetX();
    tx[1] -= primvert->GetY();
    tx[2] -= primvert->GetZ();
    Double_t innerpro = tx[0]*cand->Px()+tx[1]*cand->Py();
    Double_t signd0 = 1.;
    if(innerpro<0.) signd0 = -1.;
    fVars["signd0"] = signd0*TMath::Abs(d0z0bach[0]);
  }
  //Armenteros qT/|alpha|
  fVars["armenteros_K0s"] = v0part->PtArmV0()/TMath::Abs(v0part->AlphaV0());



  if(!fIsMixedEvent){
    for(int iProng = 0; iProng < 3; iProng++){
      AliAODTrack *dautrack;
      if(iProng == 0) dautrack = dynamic_cast<AliAODTrack *>(cand->GetDaughter(iProng));
      else            dautrack = dynamic_cast<AliAODTrack *>(v0part->GetDaughter(iProng-1));
      
      double nsigmaTPCpi = -999., nsigmaTPCK = -999., nsigmaTPCp = -999., nsigmaTOFpi = -999., nsigmaTOFK = -999., nsigmaTOFp = -999.;
      pidHF->GetnSigmaTPC(dautrack, 2, nsigmaTPCpi);
      pidHF->GetnSigmaTPC(dautrack, 3, nsigmaTPCK);
      pidHF->GetnSigmaTPC(dautrack, 4, nsigmaTPCp);
      pidHF->GetnSigmaTOF(dautrack, 2, nsigmaTOFpi);
      pidHF->GetnSigmaTOF(dautrack, 3, nsigmaTOFK);
      pidHF->GetnSigmaTOF(dautrack, 4, nsigmaTOFp);
      
      fVars[Form("nsigTPC_Pi_%d", iProng)] = nsigmaTPCpi;
      fVars[Form("nsigTPC_K_%d", iProng)]  = nsigmaTPCK;
      fVars[Form("nsigTPC_Pr_%d", iProng)]  = nsigmaTPCp;
      fVars[Form("nsigTOF_Pi_%d", iProng)] = nsigmaTOFpi;
      fVars[Form("nsigTOF_K_%d", iProng)]  = nsigmaTOFK;
      fVars[Form("nsigTOF_Pr_%d", iProng)]  = nsigmaTOFp;
      
      fVars[Form("nsigComb_Pi_%d", iProng)] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCpi, nsigmaTOFpi);
      fVars[Form("nsigComb_K_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCK, nsigmaTOFK);
      fVars[Form("nsigComb_Pr_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCp, nsigmaTOFp);
    }
  }
    
//  variables used by KFparticle

      ///! KFParticle kfpLc
  KFPVertex pVertex;
  Double_t pos[3],cov[6];
   
  ///fpVtx  === primvert
  primvert->GetXYZ(pos);
  primvert->GetCovarianceMatrix(cov);
  pVertex.SetXYZ((Float_t)pos[0],(Float_t)pos[1], (Float_t)pos[2]);
  Float_t covF[6];
  for (Int_t i=0;i<6;i++) {covF[i] = (Float_t)cov[i]; }
  pVertex.SetCovarianceMatrix(covF);
  pVertex.SetChi2(primvert->GetChi2());
  pVertex.SetNDF(primvert->GetNDF());
  pVertex.SetNContributors(primvert->GetNContributors());
  KFParticle PV(pVertex);
  // define missing TOF response as a value to be ignored by XGBoost
  Double_t nSigmaTOFPr_NaN = fVars["nsigTOF_Pr_0"];
  if (nSigmaTOFPr_NaN == -999.) nSigmaTOFPr_NaN = TMath::QuietNaN();
  fVars["nSigmaTOFPr_NaN"] = nSigmaTOFPr_NaN;
  AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(((AliAODRecoCascadeHF*)cand)->Getv0PositiveTrack());
  AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(((AliAODRecoCascadeHF*)cand)->Getv0NegativeTrack());
  
   // check charge of the first daughter, if negative, define it as the second one
  if (v0Pos->Charge()<0) {
    v0Pos = (AliAODTrack*) (((AliAODRecoCascadeHF*)cand)->Getv0NegativeTrack());
    v0Neg = (AliAODTrack*) (((AliAODRecoCascadeHF*)cand)->Getv0PositiveTrack());
  }     
  
  KFParticle kfpPionPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos,211);
  KFParticle kfpPionMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg,-211);
  
  KFParticle kfpKs0;
  const KFParticle *Ks0Daughters[2] = {&kfpPionPlus, &kfpPionMinus};
  kfpKs0.Construct(Ks0Daughters,2); 
   
  fVars["ldl_Ks0"] = AliVertexingHFUtils::ldlFromKF(kfpKs0,PV);
  KFParticle kfpKs0_massConstraint = kfpKs0;
  const Float_t massKs0_PDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  kfpKs0_massConstraint.SetNonlinearMassConstraint(massKs0_PDG);
  
  //! NEED kfpProton then go back to 882
  KFParticle kfpProton;
  AliAODTrack *bachPart = dynamic_cast<AliAODTrack*>(((AliAODRecoCascadeHF*)cand)->GetBachelor());
  if (bachPart->Charge() > 0) kfpProton = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart,2212);
  if (bachPart->Charge() < 0) kfpProton = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart,-2212);      
  
  KFParticle kfpLc;
  const KFParticle *LcDaughters[2] = {&kfpProton, &kfpKs0_massConstraint};
  kfpLc.Construct(LcDaughters,2);
  
  Double_t cosPA_Ks0 = AliVertexingHFUtils::CosPointingAngleFromKF(kfpKs0_massConstraint, kfpLc);
  Double_t cosPA_Lc  = AliVertexingHFUtils::CosPointingAngleFromKF(kfpLc, PV);
  fVars["PA_Ks0"] = TMath::ACos(cosPA_Ks0);
  fVars["PA_Lc"] = TMath::ACos(cosPA_Lc);
  
  KFParticle kfpLc_PV = kfpLc;
  kfpLc_PV.SetProductionVertex(PV);
  
  fVars["chi2topo_Lc"] = kfpLc_PV.GetChi2()/kfpLc_PV.GetNDF();
  
  if(!fIsMixedEvent){
    // Bayesian PID probability for proton
    AliPIDCombined *PIDComb  = new AliPIDCombined();
    PIDComb->SetDefaultTPCPriors();
    PIDComb->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    Double_t probTPCTOF[AliPID::kSPECIES] = {-1.};
    AliPIDResponse *pidresp = (AliPIDResponse*)pidHF->GetPidResponse();
    UInt_t detUsed = PIDComb->ComputeProbabilities(bachelor, pidresp, probTPCTOF);
    Double_t probProton = -1.;
    if (detUsed == (UInt_t)PIDComb->GetDetectorMask()) { // TPC+TOF combined
      probProton = probTPCTOF[AliPID::kProton];
    }
    else {   /// if TOF not available, try with TPC-only
      PIDComb->SetDetectorMask(AliPIDResponse::kDetTPC);
      detUsed= PIDComb->ComputeProbabilities(bachelor, pidresp, probTPCTOF);
      if (detUsed == (UInt_t)PIDComb->GetDetectorMask()) { //TPC-only worked. Else, probability returns as -1
        probProton = probTPCTOF[AliPID::kProton];
      }
    }
    fVars["CombinedPIDProb_Pr"] = probProton;
    
    delete PIDComb;
  } else {
    fIsMixedEvent = false;
  }

}

//--------------------------------------------------------------------------
void AliHFMLResponseLctoV0bachelor::SetMapOfProtonMixedEventVariables(float nsigmaTPCp, float nsigmaTOFp, float ncombp)
{
  //
  // Keep the possibility to access proton PID variables in case of mixed event background
  //

  for(int iProng = 0; iProng < 3; iProng++){

    fVars[Form("nsigTPC_Pi_%d", iProng)] = -999;
    fVars[Form("nsigTPC_K_%d", iProng)] = -999;
    fVars[Form("nsigTOF_Pi_%d", iProng)] = -999;
    fVars[Form("nsigTOF_K_%d", iProng)] = -999;

    fVars[Form("nsigComb_Pi_%d", iProng)] = -999;
    fVars[Form("nsigComb_K_%d", iProng)] = -999;
    
    if(iProng == 0){
      fVars[Form("nsigTPC_Pr_%d", iProng)] = nsigmaTPCp;
      fVars[Form("nsigTOF_Pr_%d", iProng)] = nsigmaTOFp;
      fVars[Form("nsigComb_Pr_%d", iProng)] = ncombp;
    } else {
      fVars[Form("nsigTPC_Pr_%d", iProng)] = -999;
      fVars[Form("nsigTOF_Pr_%d", iProng)] = -999;
      fVars[Form("nsigComb_Pr_%d", iProng)] = -999;
    }
  }
  fVars["CombinedPIDProb_Pr"] = -999;
  fVars["signd0"] = -9999;
  
  fIsMixedEvent = true;
}
