// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseDstartoD0pi
// \brief helper class to handle application of ML models for D*+ analyses trained
// with python libraries
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseDstartoD0pi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponseDstartoD0pi);
/// \endcond

//________________________________________________________________
AliHFMLResponseDstartoD0pi::AliHFMLResponseDstartoD0pi() : AliHFMLResponse()
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponseDstartoD0pi::AliHFMLResponseDstartoD0pi(const Char_t *name, const Char_t *title, 
                                                       const std::string configfilepath) : AliHFMLResponse(name, title, configfilepath)
{
    //
    // Standard constructor
    //
}

//________________________________________________________________
AliHFMLResponseDstartoD0pi::~AliHFMLResponseDstartoD0pi()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponseDstartoD0pi::AliHFMLResponseDstartoD0pi(const AliHFMLResponseDstartoD0pi &source) : AliHFMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliHFMLResponseDstartoD0pi &AliHFMLResponseDstartoD0pi::operator=(const AliHFMLResponseDstartoD0pi &source)
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
void AliHFMLResponseDstartoD0pi::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/)
{

    int charge = ((AliAODRecoCascadeHF*)cand)->Charge();
    AliAODRecoDecayHF2Prong *dzero = ((AliAODRecoCascadeHF*)cand)->Get2Prong();
    double massD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    fVars["pt_cand"] = cand->Pt();
    fVars["d_len"] = dzero->DecayLength();
    fVars["d_len_xy"] = dzero->DecayLengthXY();
    fVars["norm_dl"] = dzero->NormalizedDecayLength();
    fVars["norm_dl_xy"] = dzero->NormalizedDecayLengthXY();
    fVars["cos_p"] = dzero->CosPointingAngle();
    fVars["cos_p_xy"] = dzero->CosPointingAngleXY();
    fVars["imp_par_xy"] = dzero->ImpParXY();
    fVars["imp_par_prod"] = dzero->Prodd0d0();
    fVars["max_norm_d0d0exp"] = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(dzero, bfield);
    fVars["angle_D0dkpPisoft"] = ((AliAODRecoCascadeHF*)cand)->AngleD0dkpPisoft();

    AliAODTrack *dautrack[3];
    dautrack[0] = dynamic_cast<AliAODTrack *>(((AliAODRecoCascadeHF*)cand)->GetBachelor());
    if(charge > 0) {
        fVars["cos_t_star"] = dzero->CosThetaStarD0();
        fVars["delta_mass_D0"] = TMath::Abs(dzero->InvMassD0() - massD0);
        dautrack[1] = (AliAODTrack*)dzero->GetDaughter(0);
        dautrack[2] = (AliAODTrack*)dzero->GetDaughter(1);
    }
    else {
        fVars["cos_t_star"] = dzero->CosThetaStarD0bar();
        fVars["delta_mass_D0"] = TMath::Abs(dzero->InvMassD0bar() - massD0);
        dautrack[1] = (AliAODTrack*)dzero->GetDaughter(1);
        dautrack[2] = (AliAODTrack*)dzero->GetDaughter(0);
    }

    double d0[2], cov[3];
    dzero->PropagateToDCA(dzero->GetPrimaryVtx(), 0., 3., d0, cov); //propagate as a straight line
    fVars["dca"] = TMath::Sqrt(d0[0]*d0[0]+d0[1]*d0[1]);

    for (int iProng = 0; iProng < 3; iProng++)
    {
        double nsigmaTPCpi = -999., nsigmaTPCK = -999., nsigmaTOFpi = -999., nsigmaTOFK = -999.;
        pidHF->GetnSigmaTPC(dautrack[iProng], 2, nsigmaTPCpi);
        pidHF->GetnSigmaTPC(dautrack[iProng], 3, nsigmaTPCK);
        pidHF->GetnSigmaTOF(dautrack[iProng], 2, nsigmaTOFpi);
        pidHF->GetnSigmaTOF(dautrack[iProng], 3, nsigmaTOFK);

        fVars[Form("nsigTPC_Pi_%d", iProng)] = nsigmaTPCpi;
        fVars[Form("nsigTPC_K_%d", iProng)]  = nsigmaTPCK;
        fVars[Form("nsigTOF_Pi_%d", iProng)] = nsigmaTOFpi;
        fVars[Form("nsigTOF_K_%d", iProng)]  = nsigmaTOFK;

        fVars[Form("nsigComb_Pi_%d", iProng)] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCpi, nsigmaTOFpi);
        fVars[Form("nsigComb_K_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCK, nsigmaTOFK);
    }
}
