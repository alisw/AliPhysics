// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseLambdactopK0s
// \brief helper class to handle application of ML models for Lc->pK0s analyses trained
// with python libraries
// \authors:
// A. Palasciano, antonio.palasciano@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseLambdactopK0s.h"
#include "AliVertexingHFUtils.h"
#include "AliAODRecoCascadeHF.h"

// add includes for KFParticle
#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponseLambdactopK0s);
/// \endcond

//________________________________________________________________
AliHFMLResponseLambdactopK0s::AliHFMLResponseLambdactopK0s() : AliHFMLResponse()
{
    //
    // Default constructor
    //
   // fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLResponseLambdactopK0s::AliHFMLResponseLambdactopK0s(const Char_t *name, const Char_t *title, 
                                                         const std::string configfilepath) : AliHFMLResponse(name, title, configfilepath)
{
    //
    // Standard constructor
    //
   // fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLResponseLambdactopK0s::~AliHFMLResponseLambdactopK0s()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponseLambdactopK0s::AliHFMLResponseLambdactopK0s(const AliHFMLResponseLambdactopK0s &source) : AliHFMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliHFMLResponseLambdactopK0s &AliHFMLResponseLambdactopK0s::operator=(const AliHFMLResponseLambdactopK0s &source)
{
    //
    // assignment operator
    //
    if (&source == this)
        return *this;

    AliHFMLResponse::operator=(source);

    return *this;
}

//________________________________________________________________
void AliHFMLResponseLambdactopK0s::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int)
{
    /* // massHypo == 0 --> pK0s, massHypo == 1 --> pL
    if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpK0s)
        fCandType |= kLctopK0s;
    else
        fCandType &= ~kLctopK0s;

    if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpiL)
        fCandType |= kLctopiL;
    else
        fCandType &= ~kLctopiL;  */

    // topological variables
    AliAODRecoCascadeHF* candCasc = (AliAODRecoCascadeHF*)cand;
    AliAODv0 *v0part = candCasc->Getv0();
    AliAODTrack *bachPart = dynamic_cast<AliAODTrack*>(candCasc->GetBachelor());
    AliAODVertex *primVert = dynamic_cast<AliAODVertex*>(candCasc->GetPrimaryVtx());
    float massK0s = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    float massL = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    //Use KF Particle
        // convert primary vertex
        KFPVertex pVertex;
        double pos[3], cov[6];
        primVert->GetXYZ(pos);
        primVert->GetCovarianceMatrix(cov);
        float posF[3], covF[6];
        for(int iEl = 0; iEl < 3; iEl++)
            posF[iEl] = (float)pos[iEl];
        for(int iEl = 0; iEl < 6; iEl++)
            covF[iEl] = (float)cov[iEl];

        pVertex.SetXYZ((float)pos[0], (float)pos[1], (float)pos[2]);
        pVertex.SetCovarianceMatrix(covF);
        pVertex.SetChi2(primVert->GetChi2());
        pVertex.SetNDF(primVert->GetNDF());
        pVertex.SetNContributors(primVert->GetNContributors());
        KFParticle primVertKF(pVertex);

        double covP[21], covN[21], covB[21];

        int pdgBach = -1;
        int pdgV0dau[2] = {-1, -1};
        /* if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpK0s)
        { */
            pdgBach = 2212;
            pdgV0dau[0] = 211;
            pdgV0dau[1] = 211;
 /*        }
        else if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpiL)
        {
            pdgBach = 211;
            pdgV0dau[0] = 2212;
            pdgV0dau[1] = 211;
        } */

        AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(candCasc->Getv0PositiveTrack());
        AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(candCasc->Getv0NegativeTrack());
        // check charge of the first daughter, if negative, define it as the second one
        if (v0Pos->Charge()<0) {
            v0Pos = dynamic_cast<AliAODTrack*>(candCasc->Getv0NegativeTrack());
            v0Neg = dynamic_cast<AliAODTrack*>(candCasc->Getv0PositiveTrack());
        }

      //This have been performed separatly in CheckCovarianceMatrix inside the correlation task
        /* if ( !bachPart->GetCovarianceXYZPxPyPz(covB) || !v0Pos->GetCovarianceXYZPxPyPz(covP) || !v0Neg->GetCovarianceXYZPxPyPz(covN) )
            return false;
        if ( !AliVertexingHFUtils::CheckAODtrackCov(bachPart) || !AliVertexingHFUtils::CheckAODtrackCov(v0Pos) || !AliVertexingHFUtils::CheckAODtrackCov(v0Neg) )
            return false; */

        KFParticle KFBachelor;
        if(bachPart->Charge() > 0) KFBachelor = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, pdgBach);
        if(bachPart->Charge() < 0) KFBachelor = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, -pdgBach);

        KFParticle KFV0DauPlus, KFV0DauMinus;
        /* if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpK0s)
        { */
            KFV0DauPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, pdgV0dau[0]);
            KFV0DauMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -pdgV0dau[1]);
        /* }
        else if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpiL)
        {
            if (bachPart->Charge() > 0)
            {
                KFV0DauPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, pdgV0dau[0]);
                KFV0DauMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -pdgV0dau[1]);
            }
            else
            {
                KFV0DauPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, pdgV0dau[1]);
                KFV0DauMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -pdgV0dau[0]);
            }
        } */

        KFParticle KFV0;
        const KFParticle *Ks0Daughters[2] = {&KFV0DauPlus, &KFV0DauMinus};
        KFV0.Construct(Ks0Daughters, 2);
        float massV0Reco = 0., massV0RecoUnc = 0.;
        KFV0.GetMass(massV0Reco, massV0RecoUnc);

        // check V0 covariance matrix
        /* if(!AliVertexingHFUtils::CheckKFParticleCov(KFV0))
            return false; */

        // set V0 variables
        fDecayLengthV0 = AliVertexingHFUtils::DecayLengthFromKF(KFV0, primVertKF);
        fImpParV0 = cand->Getd0Prong(1); // no need to re reconstruct with KF for this
        /* if(massHypo == kpK0s)
        { */
            fMassK0s = massV0Reco;
            fMassL = -1.;
        /* }
        else if(massHypo == kpiL)
        {
            fMassK0s = -1.;
            fMassL = massV0Reco;
        } */
        fPtV0 = KFV0.GetPt();
        fcTauK0s = fDecayLengthV0 * massK0s / fPtV0;
        fcTauL = fDecayLengthV0 * massL / fPtV0;
        
        float dcaPointV0[8], dcaPointV0Cov[36];
        KFV0.GetParametersAtPoint(posF, covF, dcaPointV0, dcaPointV0Cov);
        float dcaV02 = 0;
        for(int i = 0; i < 3; i++){
            dcaV02 += (dcaPointV0[i] - pos[i]) * (dcaPointV0[i] - pos[i]);
        }
        fDCAV0 = TMath::Sqrt(dcaV02);
        fImpParV0 = KFV0.GetDistanceFromVertexXY(primVertKF);
        
        fCosPV0 = AliVertexingHFUtils::CosPointingAngleFromKF(KFV0, primVertKF);

        // Armenteros qT/|alpha|
        fArmqTOverAlphaV0 = v0part->PtArmV0()/TMath::Abs(v0part->AlphaV0()); // not done with KF in any case

        // set mass constrain before reconstructing the Lc
/*         if(massHypo == AliHFMLVarHandlerNonPromptLc2V0bachelor::kpK0s)
 */            KFV0.SetNonlinearMassConstraint(massK0s);
        /* else
            KFV0.SetNonlinearMassConstraint(massL); */

        // reconstruct Lc
        KFParticle KFLc;
        const KFParticle *LcDaughters[2] = {&KFBachelor, &KFV0};
        KFLc.Construct(LcDaughters, 2);

        // check Lc covariance matrix
        /* if (!AliVertexingHFUtils::CheckKFParticleCov(KFLc))
            return false; */

        // Lc -> K0sp variables
        fPt = KFLc.GetPt();
        fCosP = AliVertexingHFUtils::CosPointingAngleFromKF(KFLc, primVertKF);
        fCosPXY = AliVertexingHFUtils::CosPointingAngleXYFromKF(KFLc, primVertKF);
        fImpParXY = KFLc.GetDistanceFromVertexXY(primVertKF);
        fDecayLength = AliVertexingHFUtils::DecayLengthFromKF(KFLc, primVertKF);
        fDecayLengthXY = AliVertexingHFUtils::DecayLengthXYFromKF(KFLc, primVertKF);
        fNormDecayLengthXY = AliVertexingHFUtils::ldlXYFromKF(KFLc, primVertKF);
        fImpParV0 = cand->Getd0Prong(1);
        fDCAV0 = TMath::Sqrt(dcaV02);
        fImpParV0 = KFV0.GetDistanceFromVertexXY(primVertKF);
        fCosPV0 = AliVertexingHFUtils::CosPointingAngleFromKF(KFV0, primVertKF);
        fPtV0 = KFV0.GetPt();
        fcTauK0s = fDecayLengthV0 * massK0s / fPtV0;
        fcTauL = fDecayLengthV0 * massL / fPtV0;

        float dcaPoint[8], dcaPointCov[36];
        KFLc.GetParametersAtPoint(posF, covF, dcaPoint, dcaPointCov);
        float dca2 = 0;
        for(int i = 0; i < 3; i++){
            dca2 += (dcaPoint[i] - pos[i]) * (dcaPoint[i] - pos[i]);
        }
        fDCA = TMath::Sqrt(dca2);
        
        float massLcReco = 0., massLcRecoUnc = 0.;
        KFLc.GetMass(massLcReco, massLcRecoUnc);
        fInvMass = massLcReco;

        // Sign of d0 proton (different from regular d0)
        // only Lc momentum from KF (bachelor in any case not reconstructed with KF)
        double d0z0bach[2], covd0z0bach[3];
        bachPart->PropagateToDCA(primVert, bfield, kVeryBig, d0z0bach, covd0z0bach);
        double tx[3];
        bachPart->GetXYZ(tx);
        tx[0] -= primVert->GetX();
        tx[1] -= primVert->GetY();
        tx[2] -= primVert->GetZ();
        double innerpro = tx[0]*KFLc.Px()+tx[1]*KFLc.Py();
        double signd0 = 1.;
        if(innerpro<0.) signd0 = -1.;
        signd0 = signd0*TMath::Abs(d0z0bach[0]);
        fsignd0 = signd0;

        // chi2 after topological constrain to primary vertex
        KFLc.SetProductionVertex(primVertKF);
        fKFTopoChi2 = KFLc.GetChi2()/KFLc.GetNDF();

    //End use KFParticles 
    //DeltaMass is Used
    //if(fMassV0Opt == kDeltaMassV0)
    //{
        if(fMassK0s > 0)
            fMassK0s = TMath::Abs(fMassK0s-massK0s);
        if(fMassL > 0)
            fMassL = TMath::Abs(fMassL-massL);

    // single-track variables   --> GetdoProng() from here!
    AliAODTrack* prongtracks[3];
    for(unsigned int iProng = 0; iProng < fNProngs; iProng++){
        if(iProng==0)
        {
            fImpParProng[iProng] = cand->Getd0Prong(iProng);
            prongtracks[iProng] = (AliAODTrack*)candCasc->GetDaughter(iProng);
        }
        else
        {
            fImpParProng[iProng] = v0part->Getd0Prong(iProng-1);
            prongtracks[iProng] = (AliAODTrack*)v0part->GetDaughter(iProng-1);
        }
    }

    // variables used in model definition: ["d_len", "norm_dl_xy", "cos_p", "cos_p_xy", "imp_par_xy", "dca", "nsigComb_Pi_0", "nsigComb_Pr_0", "signd0", "delta_mass_K0s", "imp_par_V0", "d_len_V0", "ctau_K0s", "cos_p_V0", "KF_chi2_topo", "imp_par_prong2"]
    //fVars["pt_cand"] = cand->Pt();
    fVars["d_len"] = fDecayLength;
    //fVars["d_len_xy"] = cand->DecayLengthXY();
    //fVars["norm_dl"] = cand->NormalizedDecayLength();
    fVars["norm_dl_xy"] = fNormDecayLengthXY;
    fVars["cos_p"] = fCosP;
    fVars["cos_p_xy"] = fCosPXY;
    fVars["imp_par_xy"] = fImpParXY;
    fVars["dca"] = fDCA;

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
      
    /*fVars[Form("nsigTPC_Pi_%d", iProng)] = nsigmaTPCpi;
      fVars[Form("nsigTPC_K_%d", iProng)]  = nsigmaTPCK;
      fVars[Form("nsigTPC_Pr_%d", iProng)]  = nsigmaTPCp;
      fVars[Form("nsigTOF_Pi_%d", iProng)] = nsigmaTOFpi;
      fVars[Form("nsigTOF_K_%d", iProng)]  = nsigmaTOFK;
      fVars[Form("nsigTOF_Pr_%d", iProng)]  = nsigmaTOFp; */
      if(iProng!=0) continue;
      fVars[Form("nsigComb_Pi_%d", iProng)] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCpi, nsigmaTOFpi);
      //fVars[Form("nsigComb_K_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCK, nsigmaTOFK);
      fVars[Form("nsigComb_Pr_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCp, nsigmaTOFp);
  }
  double fKFLcInvMass = KFLc.GetMass();
  double fKFLcPt = KFLc.GetPt();
  double fKFLcEta = KFLc.GetEta();  
  double fKFLcPhi = KFLc.GetPhi();  
  double fKFLcRapidity = KFLc.GetRapidity();
  //printf("Model variables\n");
  //printf("fKFLcPhi=%2f, fKFLcPt=%.2f, fKFLcEta=%.2f, fKFLcRapidity=%.2f, fKFLcInvMass=%.2f\n",fKFLcPhi, fKFLcPt, fKFLcEta, fKFLcRapidity, fKFLcInvMass);

  fVars["signd0"] = fsignd0;
  fVars["delta_mass_K0s"] = fMassK0s;
  fVars["imp_par_V0"] = fImpParV0;
  fVars["d_len_V0"] = fDecayLengthV0;
  fVars["ctau_K0s"] = fcTauK0s;
  fVars["cos_p_V0"] = fCosPV0;
  fVars["KF_chi2_topo"] = fKFTopoChi2;
  fVars["imp_par_prong2"] = fImpParProng[2];
  

    std::cout <<"d_len =============================="<<fVars["d_len"]<<std::endl;           
    std::cout <<"norm_dl_xy =============================="<<fVars["norm_dl_xy"]<<std::endl;
    std::cout <<"cos_p =============================="<<fVars["cos_p"]<<std::endl;
    std::cout <<"cos_p_xy ==========================="<<fVars["cos_p_xy"]<<std::endl;
    std::cout <<"imp_par_xy =============================="<<fVars["imp_par_xy"]<<std::endl;
    std::cout <<"dca ================================"<<fVars["dca"]<<std::endl;    
    std::cout <<"nsigComb_Pi_0 =========================="<<fVars["nsigComb_Pi_0"]<<std::endl;    
    std::cout <<"nsigComb_Pr_0 =========================="<<fVars["nsigComb_Pr_0"]<<std::endl;    
    std::cout <<"imp_par_prong2 =========================="<<fVars["imp_par_prong2"]<<std::endl;

    std::cout <<"signd0 =========================="<<fVars["signd0"]<<std::endl;
    std::cout <<"delta_mass_K0s =========================="<<fVars["delta_mass_K0s"]<<std::endl;
    std::cout <<"imp_par_V0 =========================="<<fVars["imp_par_V0"]<<std::endl;
    std::cout <<"d_len_V0 =========================="<<fVars["d_len_V0"]<<std::endl;
    std::cout <<"ctau_K0s =========================="<<fVars["ctau_K0s"]<<std::endl;
    std::cout <<"cos_p_V0 =========================="<<fVars["cos_p_V0"]<<std::endl;
    std::cout <<"KF_chi2_topo =========================="<<fVars["KF_chi2_topo"]<<std::endl;

  //return true;
  
}
