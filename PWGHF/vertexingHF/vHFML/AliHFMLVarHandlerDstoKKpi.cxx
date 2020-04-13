// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandlerDstoKKpi
// \brief helper class to handle a tree and variables for Ds ML analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////


#include <TDatabasePDG.h>

#include "AliHFMLVarHandlerDstoKKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandlerDstoKKpi);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerDstoKKpi::AliHFMLVarHandlerDstoKKpi() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerDstoKKpi::AliHFMLVarHandlerDstoKKpi(int PIDopt, int massopt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
    fNProngs = 3; // --> cannot be changed
    SetMassKKOption(massopt);
}

//________________________________________________________________
AliHFMLVarHandlerDstoKKpi::~AliHFMLVarHandlerDstoKKpi()
{
    //
    // Default Destructor
    //
}

//________________________________________________________________
TTree* AliHFMLVarHandlerDstoKKpi::BuildTree(TString name, TString title) 
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
    //sed pid variables
    if(fPidOpt != kNoPID) 
        AddPidBranches(true, true, false, true, true);

    //set Ds variables
    TString massKKname = "";
    if(fMassKKOpt == kMassKK) 
        massKKname = "mass_KK";
    else if(fMassKKOpt == kDeltaMassKKPhi) 
        massKKname = "delta_mass_KK";

    fTreeVar->Branch("sig_vert", &fSigmaVertex);
    fTreeVar->Branch(massKKname.Data(), &fMassKK);
    fTreeVar->Branch("cos_PiDs", &fCosPiDs);
    fTreeVar->Branch("cos_PiKPhi_3", &fCosPiKPhi);
    fTreeVar->Branch("max_norm_d0d0exp", &fNormd0MeasMinusExp);

    if(fAddSingleTrackVar) {
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
    }

    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerDstoKKpi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF *pidrespo) 
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

    //Ds+ -> KKpi variables
    AliAODRecoDecayHF3Prong* cand3p = (AliAODRecoDecayHF3Prong*)cand;
    fSigmaVertex = cand3p->GetSigmaVert();
    float massPhi = 0;
    float cospikphi =-2;
    if(fMassKKOpt == kDeltaMassKKPhi)
        massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    if(masshypo == kKKpi){ //phiKKpi
        fInvMass = cand3p->InvMassDsKKpi();
        fMassKK = TMath::Abs(cand3p->InvMass2Prongs(0, 1, 321, 321) - massPhi);
        fCosPiDs = cand3p->CosPiDsLabFrameKKpi();
        cospikphi = cand3p->CosPiKPhiRFrameKKpi();
    }
    else if(masshypo == kpiKK){ //phipiKK
        fInvMass = cand3p->InvMassDspiKK();
        fMassKK = TMath::Abs(cand3p->InvMass2Prongs(1, 2, 321, 321) - massPhi);
        fCosPiDs = cand3p->CosPiDsLabFramepiKK();
        cospikphi = cand3p->CosPiKPhiRFramepiKK();
    }
    fCosPiKPhi = cospikphi * cospikphi * cospikphi;
        
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
        bool setpid = SetPidVars(prongtracks, pidrespo, true, true, false, true, true);
        if(!setpid) 
            return false;
    }

    return true;
}
