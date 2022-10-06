// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseD0toKpi
// \brief helper class to handle application of ML models for D0 analyses trained
// with python libraries
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseD0toKpi.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliHFMLVarHandlerD0toKpi.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponseD0toKpi);
/// \endcond

//________________________________________________________________
AliHFMLResponseD0toKpi::AliHFMLResponseD0toKpi() : AliHFMLResponse()
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponseD0toKpi::AliHFMLResponseD0toKpi(const Char_t *name, const Char_t *title, 
                                               const std::string configfilepath) : AliHFMLResponse(name, title, configfilepath)
{
    //
    // Standard constructor
    //
}

//________________________________________________________________
AliHFMLResponseD0toKpi::~AliHFMLResponseD0toKpi()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponseD0toKpi::AliHFMLResponseD0toKpi(const AliHFMLResponseD0toKpi &source) : AliHFMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliHFMLResponseD0toKpi &AliHFMLResponseD0toKpi::operator=(const AliHFMLResponseD0toKpi &source)
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
void AliHFMLResponseD0toKpi::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo)
{
    fVars["pt_cand"] = cand->Pt();
    fVars["d_len"] = cand->DecayLength();
    fVars["d_len_xy"] = cand->DecayLengthXY();
    fVars["norm_dl"] = cand->NormalizedDecayLength();
    fVars["norm_dl_xy"] = cand->NormalizedDecayLengthXY();
    fVars["cos_p"] = cand->CosPointingAngle();
    fVars["cos_p_xy"] = cand->CosPointingAngleXY();
    fVars["imp_par_xy"] = cand->ImpParXY();
    fVars["max_norm_d0d0exp"] = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(cand, bfield);

    AliAODRecoDecayHF2Prong* cand2p = (AliAODRecoDecayHF2Prong*)cand;
    fVars["imp_par_prod"] = cand2p->Prodd0d0();
    if(masshypo == AliHFMLVarHandlerD0toKpi::kD0)
        fVars["cos_t_star"] = cand2p->CosThetaStarD0();
    else
        fVars["cos_t_star"] = cand2p->CosThetaStarD0bar();

    for (int iProng = 0; iProng < 2; iProng++)
    {
        AliAODTrack *dautrack = dynamic_cast<AliAODTrack *>(cand->GetDaughter(iProng));

        double nsigmaTPCpi = -999., nsigmaTPCK = -999., nsigmaTOFpi = -999., nsigmaTOFK = -999.;
        pidHF->GetnSigmaTPC(dautrack, 2, nsigmaTPCpi);
        pidHF->GetnSigmaTPC(dautrack, 3, nsigmaTPCK);
        pidHF->GetnSigmaTOF(dautrack, 2, nsigmaTOFpi);
        pidHF->GetnSigmaTOF(dautrack, 3, nsigmaTOFK);

        fVars[Form("nsigTPC_Pi_%d", iProng)] = nsigmaTPCpi;
        fVars[Form("nsigTPC_K_%d", iProng)]  = nsigmaTPCK;
        fVars[Form("nsigTOF_Pi_%d", iProng)] = nsigmaTOFpi;
        fVars[Form("nsigTOF_K_%d", iProng)]  = nsigmaTOFK;

        fVars[Form("nsigComb_Pi_%d", iProng)] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCpi, nsigmaTOFpi);
        fVars[Form("nsigComb_K_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCK, nsigmaTOFK);
    }
}
