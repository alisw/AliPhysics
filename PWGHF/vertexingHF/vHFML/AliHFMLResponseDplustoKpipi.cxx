// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseDplustoKpipi
// \brief helper class to handle application of ML models for D+ analyses trained
// with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseDplustoKpipi.h"
#include "AliAODRecoDecayHF3Prong.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponseDplustoKpipi);
/// \endcond

//________________________________________________________________
AliHFMLResponseDplustoKpipi::AliHFMLResponseDplustoKpipi() : AliHFMLResponse()
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponseDplustoKpipi::AliHFMLResponseDplustoKpipi(string configfilename) : AliHFMLResponse(configfilename)
{
    //
    // Standard constructor
    //

    if (configfilename != "")
        SetConfigFile(configfilename);
}

//________________________________________________________________
AliHFMLResponseDplustoKpipi::~AliHFMLResponseDplustoKpipi()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponseDplustoKpipi::AliHFMLResponseDplustoKpipi(const AliHFMLResponseDplustoKpipi &source) : AliHFMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliHFMLResponseDplustoKpipi &AliHFMLResponseDplustoKpipi::operator=(const AliHFMLResponseDplustoKpipi &source)
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
void AliHFMLResponseDplustoKpipi::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/)
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

    for (int iProng = 0; iProng < 3; iProng++)
    {
        AliAODTrack *dautrack = dynamic_cast<AliAODTrack *>(cand->GetDaughter(iProng));

        double nsigma = -999.
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
