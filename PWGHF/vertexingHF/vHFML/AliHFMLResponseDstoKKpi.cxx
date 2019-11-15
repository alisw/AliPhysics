// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseDstoKKpi
// \brief helper class to handle application of ML models for Ds analyses trained
// with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseDstoKKpi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponseDstoKKpi);
/// \endcond

//________________________________________________________________
AliHFMLResponseDstoKKpi::AliHFMLResponseDstoKKpi() : AliHFMLResponse()
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponseDstoKKpi::AliHFMLResponseDstoKKpi(string configfilename) : AliHFMLResponse(configfilename)
{
    //
    // Standard constructor
    //

    if (configfilename != "")
        SetConfigFile(configfilename);
}

//________________________________________________________________
AliHFMLResponseDstoKKpi::~AliHFMLResponseDstoKKpi()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponseDstoKKpi::AliHFMLResponseDstoKKpi(const AliHFMLResponseDstoKKpi &source) : AliHFMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliHFMLResponseDstoKKpi &AliHFMLResponseDstoKKpi::operator=(const AliHFMLResponseDstoKKpi &source)
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
void AliHFMLResponseDstoKKpi::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo)
{
    fVars["pt_cand"] = cand->Pt();
    fVars["d_len"] = cand->DecayLength();
    fVars["d_len_xy"] = cand->DecayLengthXY();
    fVars["norm_dl"] = cand->NormalizedDecayLength();
    fVars["norm_dl_xy"] = cand->NormalizedDecayLengthXY();
    fVars["cos_p"] = cand->CosPointingAngle();
    fVars["cos_p_xy"] = cand->CosPointingAngleXY();
    fVars["imp_par_xy"] = cand->ImpParXY();
    fVars["sig_vert"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->GetSigmaVert();
    fVars["max_norm_d0d0exp"] = ComputeMaxd0MeasMinusExp(cand, bfield);

    double massPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    if (masshypo == 0)
    {
        fVars["inv_mass"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->InvMassDsKKpi();
        fVars["delta_mass_KK"] = TMath::Abs(dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->InvMass2Prongs(0, 1, 321, 321) - massPhi);
        fVars["cos_PiDs"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->CosPiDsLabFrameKKpi();
        fVars["cos_PiKPhi_3"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->CosPiKPhiRFrameKKpi();
    }
    else if (masshypo == 1)
    {
        fVars["inv_mass"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->InvMassDspiKK();
        fVars["delta_mass_KK"] = TMath::Abs(dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->InvMass2Prongs(1, 2, 321, 321) - massPhi);
        fVars["cos_PiDs"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->CosPiDsLabFramepiKK();
        fVars["cos_PiKPhi_3"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->CosPiKPhiRFramepiKK();
    }
    fVars["cos_PiKPhi_3"] = fVars["cos_PiKPhi_3"] * fVars["cos_PiKPhi_3"] * fVars["cos_PiKPhi_3"];

    for (int iProng = 0; iProng < 3; iProng++)
    {
        AliAODTrack *dautrack = dynamic_cast<AliAODTrack *>(cand->GetDaughter(iProng));

        double nsigma = -999.;
        pidHF->GetnSigmaTPC(dautrack, 2, nsigma);
        fVars[Form("nsigTPC_Pi_%d", iProng)] = nsigma;
        pidHF->GetnSigmaTPC(dautrack, 3, nsigma);
        fVars[Form("nsigTPC_K_%d", iProng)] = nsigma;
        pidHF->GetnSigmaTOF(dautrack, 2, nsigma);
        fVars[Form("nsigTOF_Pi_%d", iProng)] = nsigma;
        pidHF->GetnSigmaTOF(dautrack, 3, nsigma);
        fVars[Form("nsigTOF_K_%d", iProng)] = nsigma;

        fVars[Form("nsigComb_Pi_%d", iProng)] = CombineNsigmaTPCTOF(fVars[Form("nsigTPC_Pi_%d", iProng)], fVars[Form("nsigTOF_Pi_%d", iProng)]);
        fVars[Form("nsigComb_K_%d", iProng)] = CombineNsigmaTPCTOF(fVars[Form("nsigTPC_K_%d", iProng)], fVars[Form("nsigTOF_K_%d", iProng)]);
    }
}
