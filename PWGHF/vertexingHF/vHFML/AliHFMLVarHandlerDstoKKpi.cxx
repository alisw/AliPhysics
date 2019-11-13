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

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandlerDstoKKpi);
/// \endcond

//________________________________________________________________
AliHFMLVarHandlerDstoKKpi::AliHFMLVarHandlerDstoKKpi() : AliHFMLVarHandler()
{
    //
    // Default constructor
    //
    fNProngs=3; // --> cannot be changed
}

//________________________________________________________________
AliHFMLVarHandlerDstoKKpi::AliHFMLVarHandlerDstoKKpi(int PIDopt) : AliHFMLVarHandler(PIDopt)
{
    //
    // Standard constructor
    //
    fNProngs=3; // --> cannot be changed
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
        fTreeVar=nullptr;
    }
    fTreeVar = new TTree(name.Data(),title.Data());

    //set common variables
    AddCommonDmesonVarBranches();
    //set single-track variables
    AddSingleTrackBranches();
    //sed pid variables
    if(fPidOpt!=kNoPID) 
        AddPidBranches(true,true,false,true,true);

    //set Ds variables
    TString massKKname="";
    if(fMassKKOpt==kMassKK) 
        massKKname = "mass_KK";
    else if(fMassKKOpt==kDeltaMassKKPhi) 
        massKKname = "delta_mass_KK";

    fTreeVar->Branch("sig_vert",&fSigmaVertex);
    fTreeVar->Branch(massKKname.Data(),&fMassKK);
    fTreeVar->Branch("cos_PiDs",&fCosPiDs);
    fTreeVar->Branch("cos_PiKPhi_3",&fCosPiKPhi);
    fTreeVar->Branch("max_norm_d0d0exp",&fNormd0MeasMinusExp);
    for(unsigned int iProng=0; iProng<fNProngs; iProng++)
        fTreeVar->Branch(Form("imp_par_prong%d",iProng),&fImpParProng[iProng]);
        
    return fTreeVar;
}

//________________________________________________________________
bool AliHFMLVarHandlerDstoKKpi::SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo, AliAODPidHF *pidrespo) 
{
    if(!cand) return false;
    if(fFillOnlySignal) { //if fill only signal and not signal candidate, do not store
        if(!(fCandType&kSignal || fCandType&kRefl)) return true;
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
    fDCA  =cand->GetDCA();
    fNormd0MeasMinusExp = ComputeMaxd0MeasMinusExp(cand,bfield);

    //Ds+ -> KKpi variables
    fSigmaVertex = ((AliAODRecoDecayHF3Prong*)cand)->GetSigmaVert();
    float massPhi = 0;
    float cospikphi =-2;
    if(fMassKKOpt==kDeltaMassKKPhi) 
        massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    if(masshypo==0){ //phiKKpi
        fInvMass = ((AliAODRecoDecayHF3Prong*)cand)->InvMassDsKKpi();
        fMassKK = TMath::Abs(((AliAODRecoDecayHF3Prong*)cand)->InvMass2Prongs(0,1,321,321)-massPhi);
        fCosPiDs = ((AliAODRecoDecayHF3Prong*)cand)->CosPiDsLabFrameKKpi();
        cospikphi = ((AliAODRecoDecayHF3Prong*)cand)->CosPiKPhiRFrameKKpi();
    }
    else if(masshypo==1){ //phipiKK
        fInvMass = ((AliAODRecoDecayHF3Prong*)cand)->InvMassDspiKK();
        fMassKK = TMath::Abs(((AliAODRecoDecayHF3Prong*)cand)->InvMass2Prongs(1,2,321,321)-massPhi);
        fCosPiDs = ((AliAODRecoDecayHF3Prong*)cand)->CosPiDsLabFramepiKK();
        cospikphi = ((AliAODRecoDecayHF3Prong*)cand)->CosPiKPhiRFramepiKK();
    }
    fCosPiKPhi = cospikphi*cospikphi*cospikphi;
    for(unsigned int iProng=0; iProng<fNProngs; iProng++)
        fImpParProng[iProng] = cand->Getd0Prong(iProng);
        
    //single-track variables
    AliAODTrack* prongtracks[3];
    for(unsigned int iProng=0; iProng<fNProngs; iProng++) prongtracks[iProng] = (AliAODTrack*)cand->GetDaughter(iProng);
    bool setsingletrack = SetSingleTrackVars(prongtracks);  
    if(!setsingletrack) return false;

    //pid variables
    if(fPidOpt==kNoPID) return true;

    bool setpid = SetPidVars(prongtracks,pidrespo,true,true,false,true,true);
    if(!setpid) return false;

    return true;
}
