// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandlerDplustoKKpi
// \brief helper class to handle a tree and variables for D+ ML analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandlerDplustoKpipi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandlerDplustoKpipi);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerDplustoKpipi::AliHFMLVarHandlerDplustoKpipi() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerDplustoKpipi::AliHFMLVarHandlerDplustoKpipi(int PIDopt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
   fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerDplustoKpipi::~AliHFMLVarHandlerDplustoKpipi()
{
    //
    // Default Destructor
    //
}

//________________________________________________________________
TTree* AliHFMLVarHandlerDplustoKpipi::BuildTree(TString name, TString title) 
{
    if(fTreeVar) {
        delete fTreeVar;
        fTreeVar = nullptr;
    }
    fTreeVar = new TTree(name.Data(), title.Data());

    //set common variables
    AddCommonDmesonVarBranches();
    //set single-track variables
    if(fAddSingleTrackVar)
        AddSingleTrackBranches();
    //set PID variables
    if(fPidOpt != kNoPID) 
        AddPidBranches(true, true, false, true, true);

    //set D+ variables
    fTreeVar->Branch("sig_vert", &fSigmaVertex);
    fTreeVar->Branch("max_norm_d0d0exp", &fNormd0MeasMinusExp);

    if(fAddSingleTrackVar) {
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
    }

    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerDplustoKpipi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliAODPidHF *pidrespo) 
{
    if(!cand) 
        return false;
    if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
        if(!(fCandType&kSignal)) 
            return true;
    }
    
    fCandType &= ~kRefl; //protection --> D+->Kpipi cannot be reflected

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

    //D+ -> Kpipi variables
    AliAODRecoDecayHF3Prong* cand3p = (AliAODRecoDecayHF3Prong*)cand;
    fInvMass = cand3p->InvMassDplus();
    fSigmaVertex = cand3p->GetSigmaVert();
        
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
