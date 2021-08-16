// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandlerD0toKpi
// \brief helper class to handle a tree and variables for D+ ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandlerD0toKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandlerD0toKpi);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerD0toKpi::AliHFMLVarHandlerD0toKpi() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs = 2; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerD0toKpi::AliHFMLVarHandlerD0toKpi(int PIDopt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
   fNProngs = 2; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerD0toKpi::~AliHFMLVarHandlerD0toKpi()
{
    //
    // Default Destructor
    //
}

//________________________________________________________________
TTree* AliHFMLVarHandlerD0toKpi::BuildTree(TString name, TString title) 
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
    //set PID variables
    if(fPidOpt != kNoPID) 
        AddPidBranches(true, true, false, true, true);

    //set D+ variables
    fTreeVar->Branch("imp_par_prod", &fImpParProd);
    fTreeVar->Branch("max_norm_d0d0exp", &fNormd0MeasMinusExp);
    fTreeVar->Branch("cos_t_star", &fCosThetaStar);

    if(fAddSingleTrackVar) {
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
    }

    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerD0toKpi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF *pidrespo) 
{
    if(!cand) 
        return false;
    if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
        if(!(fCandType & kSignal)) 
            return true;
    }

    //topological variables
    //common
    fPt = cand->Pt();
    fDecayLength = cand->DecayLength();
    fDecayLengthXY = cand->DecayLengthXY();
    fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    fCosP = cand->CosPointingAngle();
    fCosPXY = cand->CosPointingAngleXY();
    fImpParXY = cand->ImpParXY();
    fNormd0MeasMinusExp = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(cand, bfield);

    double d0[2], cov[3];
    cand->PropagateToDCA(cand->GetPrimaryVtx(), 0., 3., d0, cov); //propagate as a straight line
    fDCA = TMath::Sqrt(d0[0]*d0[0]+d0[1]*d0[1]);

    //D0 -> Kpi variables
    AliAODRecoDecayHF2Prong* cand2p = (AliAODRecoDecayHF2Prong*)cand;
    fImpParProd = cand2p->Prodd0d0();
    if(masshypo == kD0) {
        fInvMass = cand2p->InvMassD0();
        fCosThetaStar = cand2p->CosThetaStarD0();
    }
    else {
        fInvMass = cand2p->InvMassD0bar();
        fCosThetaStar = cand2p->CosThetaStarD0bar();
    }

    //single track variables
    AliAODTrack* prongtracks[3];
    for(unsigned int iProng = 0; iProng < fNProngs; iProng++) 
        prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);

    if(fAddSingleTrackVar) {
        bool setsingletrack = SetSingleTrackVars(prongtracks);  
        if(!setsingletrack) 
            return false;
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fImpParProng[iProng] = cand->Getd0Prong(iProng);
    }

    //pid variables
    if(fPidOpt != kNoPID) {
        bool setpid = SetPidVars(prongtracks, pidrespo, true, true, false, true, true);
        if(!setpid) 
            return false;
    }

    return true;
}
