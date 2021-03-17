// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandlerNonPromptLc2V0bachelor
// \brief helper class to handle a tree and variables for non-prompt Lc->V0bachelor ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////


#include <TDatabasePDG.h>

#include "AliHFMLVarHandlerNonPromptLc2V0bachelor.h"
#include "AliAODRecoCascadeHF.h"
#include "AliVertexingHFUtils.h"

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
ClassImp(AliHFMLVarHandlerNonPromptLc2V0bachelor);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerNonPromptLc2V0bachelor::AliHFMLVarHandlerNonPromptLc2V0bachelor() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerNonPromptLc2V0bachelor::AliHFMLVarHandlerNonPromptLc2V0bachelor(int PIDopt, int massOpt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
    fNProngs = 3; // --> cannot be changed
    SetMassV0Option(massOpt);
}

//________________________________________________________________
AliHFMLVarHandlerNonPromptLc2V0bachelor::~AliHFMLVarHandlerNonPromptLc2V0bachelor()
{
    //
    // Default Destructor
    //
}

//________________________________________________________________
TTree* AliHFMLVarHandlerNonPromptLc2V0bachelor::BuildTree(TString name, TString title)
{
    if(fTreeVar) {
        delete fTreeVar;
        fTreeVar = nullptr;
    }
    fTreeVar = new TTree(name.Data(), title.Data());

    //set common variables
    AddCommonDmesonVarBranches();
    //set global event variables
    AddGlobalEventVarBranches();
    //set single-track variables
    if(fAddSingleTrackVar)
        AddSingleTrackBranches();
    //sed pid variables
    if(fPidOpt != kNoPID)
        AddPidBranches(true, false, true, true, true);

    //set Lc->V0bachelor specific variables
    TString massK0sname = "", massLname = "";
    if(fMassV0Opt == kMassV0)
    {
        massLname = "mass_L";
        massK0sname = "mass_K0s";
    }
    else if(fMassV0Opt == kDeltaMassV0)
    {
        massLname = "delta_mass_L";
        massK0sname = "delta_mass_K0s";
    }

    fTreeVar->Branch("signd0", &fsignd0);
    fTreeVar->Branch(massK0sname.Data(), &fMassK0s);
    fTreeVar->Branch(massLname.Data(), &fMassL);
    fTreeVar->Branch("dca_V0", &fDCAV0);
    fTreeVar->Branch("imp_par_V0", &fImpParV0);
    fTreeVar->Branch("imp_par_V0_xy", &fImpParV0XY);
    fTreeVar->Branch("d_len_V0", &fDecayLengthV0);
    fTreeVar->Branch("armenteros_V0", &fArmqTOverAlphaV0);
    fTreeVar->Branch("ctau_K0s", &fcTauK0s);
    fTreeVar->Branch("ctau_L", &fcTauL);
    fTreeVar->Branch("cos_p_V0", &fCosPV0);
    fTreeVar->Branch("pt_V0", &fPtV0);
    if(fUseKFParticle)
        fTreeVar->Branch("KF_chi2_topo", &fKFTopoChi2);

    if(fAddSingleTrackVar) {
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
    }

    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerNonPromptLc2V0bachelor::SetVariables(AliAODRecoDecayHF* cand, float bfield, int massHypo, AliAODPidHF *pidrespo)
{
    if(!cand)
        return false;
    if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
        if(!(fCandType&kSignal || fCandType&kRefl))
            return true;
    }

    // massHypo == 0 --> pK0s, massHypo == 1 --> pL
    if(massHypo == kpK0s)
        fCandType |= kLctopK0s;
    else
        fCandType &= ~kLctopK0s;

    if(massHypo == kpiL)
        fCandType |= kLctopiL;
    else
        fCandType &= ~kLctopiL;

    // topological variables
    AliAODRecoCascadeHF* candCasc = (AliAODRecoCascadeHF*)cand;
    AliAODv0 *v0part = candCasc->Getv0();
    AliAODTrack *bachPart = dynamic_cast<AliAODTrack*>(candCasc->GetBachelor());
    AliAODVertex *primVert = dynamic_cast<AliAODVertex*>(candCasc->GetPrimaryVtx());
    float massK0s = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    float massL = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    if(!fUseKFParticle)
    {
        // common
        fPt = candCasc->Pt();
        fDecayLength = candCasc->DecayLength();
        fDecayLengthXY = candCasc->DecayLengthXY();
        fNormDecayLengthXY = candCasc->NormalizedDecayLengthXY();
        fCosP = candCasc->CosPointingAngle();
        fCosPXY = candCasc->CosPointingAngleXY();
        fImpParXY = candCasc->ImpParXY();

        double d0[2], cov[3];
        cand->PropagateToDCA(cand->GetPrimaryVtx(), 0., 3., d0, cov); //propagate as a straight line
        fDCA = TMath::Sqrt(d0[0]*d0[0]+d0[1]*d0[1]);

        // Lc -> K0sp variables
        if(massHypo == kpK0s)
            fInvMass = candCasc->InvMassLctoK0sP();
        else if(massHypo == kpiL)
            fInvMass = candCasc->InvMassLctoLambdaPi();

        fImpParV0 = cand->Getd0Prong(1);
        fDecayLengthV0 = candCasc->DecayLengthV0();
        fMassK0s = v0part->MassK0Short();
        fMassL = v0part->MassLambda();
        fDCAV0 = v0part->GetDCA();
        fPtV0 = v0part->Pt();

        fcTauK0s = candCasc->DecayLengthV0() * massK0s / v0part->P();
        fcTauL = candCasc->DecayLengthV0() * massL / v0part->P();
        fCosPV0 = candCasc->CosV0PointingAngle();

        // Sign of d0 proton (different from regular d0)
        // (from AliRDHFCutsLctoV0)
        double d0z0bach[2], covd0z0bach[3];
        bachPart->PropagateToDCA(primVert, bfield, kVeryBig, d0z0bach, covd0z0bach);
        double tx[3];
        bachPart->GetXYZ(tx);
        tx[0] -= primVert->GetX();
        tx[1] -= primVert->GetY();
        tx[2] -= primVert->GetZ();
        double innerpro = tx[0]*cand->Px()+tx[1]*cand->Py();
        double signd0 = 1.;
        if(innerpro<0.) signd0 = -1.;
        signd0 = signd0*TMath::Abs(d0z0bach[0]);
        fsignd0 = signd0;

        // Armenteros qT/|alpha|
        fArmqTOverAlphaV0 = v0part->PtArmV0()/TMath::Abs(v0part->AlphaV0());
    }
    else
    {
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
        if(massHypo == kpK0s)
        {
            pdgBach = 2212;
            pdgV0dau[0] = 211;
            pdgV0dau[1] = 211;
        }
        else if(massHypo == kpiL)
        {
            pdgBach = 211;
            pdgV0dau[0] = 2212;
            pdgV0dau[1] = 211;
        }

        AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(candCasc->Getv0PositiveTrack());
        AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(candCasc->Getv0NegativeTrack());
        // check charge of the first daughter, if negative, define it as the second one
        if (v0Pos->Charge()<0) {
            v0Pos = dynamic_cast<AliAODTrack*>(candCasc->Getv0NegativeTrack());
            v0Neg = dynamic_cast<AliAODTrack*>(candCasc->Getv0PositiveTrack());
        }

        if ( !bachPart->GetCovarianceXYZPxPyPz(covB) || !v0Pos->GetCovarianceXYZPxPyPz(covP) || !v0Neg->GetCovarianceXYZPxPyPz(covN) )
            return false;
        if ( !AliVertexingHFUtils::CheckAODtrackCov(bachPart) || !AliVertexingHFUtils::CheckAODtrackCov(v0Pos) || !AliVertexingHFUtils::CheckAODtrackCov(v0Neg) )
            return false;

        KFParticle KFBachelor;
        if(bachPart->Charge() > 0) KFBachelor = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, pdgBach);
        if(bachPart->Charge() < 0) KFBachelor = AliVertexingHFUtils::CreateKFParticleFromAODtrack(bachPart, -pdgBach);

        KFParticle KFV0DauPlus, KFV0DauMinus;
        if(massHypo == kpK0s)
        {
            KFV0DauPlus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Pos, pdgV0dau[0]);
            KFV0DauMinus = AliVertexingHFUtils::CreateKFParticleFromAODtrack(v0Neg, -pdgV0dau[1]);
        }
        else if(massHypo == kpiL)
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
        }

        KFParticle KFV0;
        const KFParticle *Ks0Daughters[2] = {&KFV0DauPlus, &KFV0DauMinus};
        KFV0.Construct(Ks0Daughters, 2);
        float massV0Reco = 0., massV0RecoUnc = 0.;
        KFV0.GetMass(massV0Reco, massV0RecoUnc);

        // check V0 covariance matrix
        if(!AliVertexingHFUtils::CheckKFParticleCov(KFV0))
            return false;

        // set V0 variables
        fDecayLengthV0 = AliVertexingHFUtils::DecayLengthFromKF(KFV0, primVertKF);
        fImpParV0 = cand->Getd0Prong(1); // no need to re reconstruct with KF for this
        if(massHypo == kpK0s)
        {
            fMassK0s = massV0Reco;
            fMassL = -1.;
        }
        else if(massHypo == kpiL)
        {
            fMassK0s = -1.;
            fMassL = massV0Reco;
        }
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
        fImpParV0XY = KFV0.GetDistanceFromVertexXY(primVertKF);
        
        fCosPV0 = AliVertexingHFUtils::CosPointingAngleFromKF(KFV0, primVertKF);

        // Armenteros qT/|alpha|
        fArmqTOverAlphaV0 = v0part->PtArmV0()/TMath::Abs(v0part->AlphaV0()); // not done with KF in any case

        // set mass constrain before reconstructing the Lc
        if(massHypo == kpK0s)
            KFV0.SetNonlinearMassConstraint(massK0s);
        else
            KFV0.SetNonlinearMassConstraint(massL);

        // reconstruct Lc
        KFParticle KFLc;
        const KFParticle *LcDaughters[2] = {&KFBachelor, &KFV0};
        KFLc.Construct(LcDaughters, 2);

        // check Lc covariance matrix
        if (!AliVertexingHFUtils::CheckKFParticleCov(KFLc))
            return false;

        // Lc -> K0sp variables
        fPt = KFLc.GetPt();
        fCosP = AliVertexingHFUtils::CosPointingAngleFromKF(KFLc, primVertKF);
        fCosPXY = AliVertexingHFUtils::CosPointingAngleXYFromKF(KFLc, primVertKF);
        fImpParXY = KFLc.GetDistanceFromVertexXY(primVertKF);
        fDecayLength = AliVertexingHFUtils::DecayLengthFromKF(KFLc, primVertKF);
        fDecayLengthXY = AliVertexingHFUtils::DecayLengthXYFromKF(KFLc, primVertKF);
        fNormDecayLengthXY = AliVertexingHFUtils::ldlXYFromKF(KFLc, primVertKF);

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
    }

    if(fMassV0Opt == kDeltaMassV0)
    {
        if(fMassK0s > 0)
            fMassK0s = TMath::Abs(fMassK0s-massK0s);
        if(fMassL > 0)
            fMassL = TMath::Abs(fMassL-massL);
    }

    // single-track variables
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

    if(fAddSingleTrackVar) {
        bool setsingletrack = SetSingleTrackVars(prongtracks);
        if(!setsingletrack)
            return false;
    }

    // pid variables
    if(fPidOpt != kNoPID) {
        bool setpid = SetPidVars(prongtracks, pidrespo, true, false, true, true, true);
        if(!setpid)
            return false;
    }

    return true;
}
