// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandlerDstartoD0pi
// \brief helper class to handle a tree and variables for D+ ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandlerDstartoD0pi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandlerDstartoD0pi);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerDstartoD0pi::AliHFMLVarHandlerDstartoD0pi() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerDstartoD0pi::AliHFMLVarHandlerDstartoD0pi(int PIDopt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
   fNProngs = 3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerDstartoD0pi::~AliHFMLVarHandlerDstartoD0pi()
{
    //
    // Default Destructor
    //
}

//________________________________________________________________
TTree* AliHFMLVarHandlerDstartoD0pi::BuildTree(TString name, TString title) 
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

    //set D*+ variables
    fTreeVar->Branch("cos_t_star", &fCosThetaStar);
    fTreeVar->Branch("imp_par_prod", &fImpParProd);
    fTreeVar->Branch("max_norm_d0d0exp", &fNormd0MeasMinusExp);
    fTreeVar->Branch("angle_D0dkpPisoft", &fAngleD0dkpPisoft);
    fTreeVar->Branch("delta_mass_D0", &fDeltaInvMassD0);

    if(fAddSingleTrackVar) {
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
            fTreeVar->Branch(Form("imp_par_prong%d", iProng), &fImpParProng[iProng]);
    }

    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerDstartoD0pi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int /*masshypo*/, AliAODPidHF *pidrespo) 
{
    if(!cand) 
        return false;
    if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
        if(!(fCandType & kSignal)) 
            return true;
    }
    
    fCandType &= ~kRefl; //protection --> D*+->Kpipi cannot be reflected

    AliAODRecoDecayHF2Prong *dzero = ((AliAODRecoCascadeHF*)cand)->Get2Prong();

    //topological variables
    //common
    fPt = cand->Pt();
    fDecayLength = dzero->DecayLength();
    fDecayLengthXY = dzero->DecayLengthXY();
    fNormDecayLengthXY = dzero->NormalizedDecayLengthXY();
    fCosP = dzero->CosPointingAngle();
    fCosPXY = dzero->CosPointingAngleXY();
    fImpParXY = dzero->ImpParXY();
    fNormd0MeasMinusExp = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(dzero, bfield);

    double d0[2], cov[3];
    dzero->PropagateToDCA(dzero->GetPrimaryVtx(), 0., 3., d0, cov); //propagate as a straight line
    fDCA = TMath::Sqrt(d0[0]*d0[0]+d0[1]*d0[1]);

    //D*+ -> Kpipi variables
    double massD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    AliAODTrack* prongtracks[3];
    prongtracks[0] = (AliAODTrack*)((AliAODRecoCascadeHF*)cand)->GetBachelor();

    fInvMass = ((AliAODRecoCascadeHF*)cand)->DeltaInvMass();
    fAngleD0dkpPisoft=((AliAODRecoCascadeHF*)cand)->AngleD0dkpPisoft();
    fImpParProd = dzero->Prodd0d0();
    if( (((AliAODRecoCascadeHF*)cand)->Charge()) > 0) {
        fCosThetaStar = dzero->CosThetaStarD0();
        fDeltaInvMassD0 = TMath::Abs(dzero->InvMassD0() - massD0);
        prongtracks[1] = (AliAODTrack*)dzero->GetDaughter(0);
        prongtracks[2] = (AliAODTrack*)dzero->GetDaughter(1);
    }
    else {
        fCosThetaStar = dzero->CosThetaStarD0bar();
        fDeltaInvMassD0 = TMath::Abs(dzero->InvMassD0bar() - massD0);
        prongtracks[1] = (AliAODTrack*)dzero->GetDaughter(1);
        prongtracks[2] = (AliAODTrack*)dzero->GetDaughter(0);
    }

    if(fAddSingleTrackVar) {
        bool setsingletrack = SetSingleTrackVars(prongtracks);  
        if(!setsingletrack) 
            return false;
        for(unsigned int iProng = 0; iProng < fNProngs; iProng++){
            if(iProng == 0)
                fImpParProng[iProng] = cand->Getd0Prong(iProng);
            else
                fImpParProng[iProng] = dzero->Getd0Prong(iProng-1);
        }
    }

    //pid variables
    if(fPidOpt != kNoPID) {
        bool setpid = SetPidVars(prongtracks, pidrespo, true, true, false, true, true);
        if(!setpid) 
            return false;
    }

    return true;
}
