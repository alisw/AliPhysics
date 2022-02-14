// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandlerLctopKpi
// \brief helper class to handle a tree and variables for Lc->pKpi ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////


#include <TDatabasePDG.h>

#include "AliHFMLVarHandlerLctopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandlerLctopKpi);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerLctopKpi::AliHFMLVarHandlerLctopKpi() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerLctopKpi::AliHFMLVarHandlerLctopKpi(int PIDopt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
    fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerLctopKpi::~AliHFMLVarHandlerLctopKpi()
{
    //
    // Default Destructor
    //
}

//________________________________________________________________
TTree* AliHFMLVarHandlerLctopKpi::BuildTree(TString name, TString title) 
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
        AddPidBranches(true, true, true, true, true);

    //set Lc->pKpi variables
    fTreeVar->Branch("sig_vert", &fSigmaVertex);
    fTreeVar->Branch("max_norm_d0d0exp", &fNormd0MeasMinusExp);

    if(fAddSingleTrackVar) {
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
    }

    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerLctopKpi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF *pidrespo) 
{
    if(!cand) 
        return false;
    if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
        if(!(fCandType&kSignal || fCandType&kRefl)) 
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

    //Lc -> pKpi variables
    AliAODRecoDecayHF3Prong* cand3p = (AliAODRecoDecayHF3Prong*)cand;
    fSigmaVertex = cand3p->GetSigmaVert();

    if(masshypo==kpKpi)
        fInvMass = cand3p->InvMassLcpKpi();
    else if(masshypo==kpiKp)
        fInvMass = cand3p->InvMassLcpiKp();
        
    //single-track variables
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
        bool setpid = SetPidVars(prongtracks, pidrespo, true, true, true, true, true);
        if(!setpid) 
            return false;
    }

    return true;
}
//----------------------------------------------------------------
void AliHFMLVarHandlerLctopKpi::SetIsLcpKpiRes(int restype){
    if(restype==1)
        fCandType  |= kLcNonRes;
    else
        fCandType  &= ~kLcNonRes;
    if(restype==2)
        fCandType  |= kLcKStar;
    else
        fCandType  &= ~kLcKStar;
    if(restype==3)
        fCandType  |= kLcDelta;
    else
        fCandType  &= ~kLcDelta;
    if(restype==4)
        fCandType  |= kLcLambda1520;
    else
        fCandType  &= ~kLcLambda1520;
}
