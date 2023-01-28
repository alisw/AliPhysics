// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseLambdactopKpi
// \brief helper class to handle application of ML models for Lc analyses trained
// with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliHFMLResponseLambdactopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponseLambdactopKpi);
/// \endcond

//________________________________________________________________
AliHFMLResponseLambdactopKpi::AliHFMLResponseLambdactopKpi() : AliHFMLResponse()
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponseLambdactopKpi::AliHFMLResponseLambdactopKpi(const Char_t *name, const Char_t *title, 
                                                         const std::string configfilepath) : AliHFMLResponse(name, title, configfilepath)
{

    //
    // Standard constructor
    //
}

//________________________________________________________________
AliHFMLResponseLambdactopKpi::~AliHFMLResponseLambdactopKpi()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponseLambdactopKpi::AliHFMLResponseLambdactopKpi(const AliHFMLResponseLambdactopKpi &source) : AliHFMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliHFMLResponseLambdactopKpi &AliHFMLResponseLambdactopKpi::operator=(const AliHFMLResponseLambdactopKpi &source)
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
void AliHFMLResponseLambdactopKpi::SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/)
{

   //fVars["pt_cand"] = cand->Pt();
    fVars["d_len"] = cand->DecayLength();
    fVars["d_len_xy"] = cand->DecayLengthXY();
    //fVars["norm_dl"] = cand->NormalizedDecayLength();
    fVars["norm_dl_xy"] = cand->NormalizedDecayLengthXY();
    fVars["cos_p"] = cand->CosPointingAngle();
    fVars["cos_p_xy"] = cand->CosPointingAngleXY();
    fVars["imp_par_xy"] = cand->ImpParXY();
    double d0[2], cov[3];
    cand->PropagateToDCA(cand->GetPrimaryVtx(), 0., 3., d0, cov); //propagate as a straight line
    fVars["dca"] = TMath::Sqrt(d0[0]*d0[0]+d0[1]*d0[1]);
    
   
    
    for (int iProng = 0; iProng < 3; iProng++)
    {
    AliAODTrack *dautrack = dynamic_cast<AliAODTrack *>(cand->GetDaughter(iProng));
	//fVars[Form("pt_prong%d", iProng)] = dautrack->Pt();
	
        double nsigmaTPCpi = -999., nsigmaTPCK = -999., nsigmaTPCp = -999., nsigmaTOFpi = -999., nsigmaTOFK = -999., nsigmaTOFp = -999.;
        pidHF->GetnSigmaTPC(dautrack, 2, nsigmaTPCpi);
        pidHF->GetnSigmaTPC(dautrack, 3, nsigmaTPCK);
        pidHF->GetnSigmaTPC(dautrack, 4, nsigmaTPCp);
        pidHF->GetnSigmaTOF(dautrack, 2, nsigmaTOFpi);
        pidHF->GetnSigmaTOF(dautrack, 3, nsigmaTOFK);
        pidHF->GetnSigmaTOF(dautrack, 4, nsigmaTOFp);
  
     
        fVars[Form("nsigComb_Pi_%d", iProng)] = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCpi, nsigmaTOFpi);
        fVars[Form("nsigComb_K_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCK, nsigmaTOFK);
        fVars[Form("nsigComb_Pr_%d", iProng)]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCp, nsigmaTOFp);
        
        //fVars["nsigComb_Pr_0"]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCp, nsigmaTOFp);
        //fVars["nsigComb_K_1"]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCK, nsigmaTOFK);
        //fVars["nsigComb_Pr_2"]  = AliVertexingHFUtils::CombineNsigmaTPCTOF(nsigmaTPCp, nsigmaTOFp);
    }
    
    fVars["sig_vert"] = dynamic_cast<AliAODRecoDecayHF3Prong *>(cand)->GetSigmaVert();
    fVars["max_norm_d0d0exp"] = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(cand, bfield);
    
    /* for(unsigned int iProng1 = 0; iProng1 < 3; iProng1++){
            fVars[Form("imp_par_prong%d", iProng1)] = cand->Getd0Prong(iProng1);
    } */
    std::cout <<"d_len =============================="<<fVars["d_len"]<<std::endl;           
    std::cout <<"d_len_xy =============================="<<fVars["d_len_xy"]<<std::endl;       
    std::cout <<"norm_dl_xy =============================="<<fVars["norm_dl_xy"]<<std::endl;
    std::cout <<"cos_p =============================="<<fVars["cos_p"]<<std::endl;
    std::cout <<"cos_p_xy ==========================="<<fVars["cos_p_xy"]<<std::endl;
    std::cout <<"imp_par_xy =============================="<<fVars["imp_par_xy"]<<std::endl;
    std::cout <<"dca ================================"<<fVars["dca"]<<std::endl;    
    /*std::cout <<"pt_prong0 =========================="<<fVars["pt_prong0"]<<std::endl;    
    std::cout <<"pt_prong1 =========================="<<fVars["pt_prong1"]<<std::endl;    
    std::cout <<"pt_prong2 =========================="<<fVars["pt_prong2"]<<std::endl;    */ 
    std::cout <<"nsigComb_Pi_0 =========================="<<fVars["nsigComb_Pi_0"]<<std::endl;    
    std::cout <<"nsigComb_K_0 =========================="<<fVars["nsigComb_K_0"]<<std::endl;    
    std::cout <<"nsigComb_Pr_0 =========================="<<fVars["nsigComb_Pr_0"]<<std::endl;    
    std::cout <<"nsigComb_Pi_1 =========================="<<fVars["nsigComb_Pi_1"]<<std::endl;    
    std::cout <<"nsigComb_K_1 =========================="<<fVars["nsigComb_K_1"]<<std::endl;    
    std::cout <<"nsigComb_Pr_1 =========================="<<fVars["nsigComb_Pr_1"]<<std::endl;    
    std::cout <<"nsigComb_Pi_2 =========================="<<fVars["nsigComb_Pi_2"]<<std::endl;    
    std::cout <<"nsigComb_K_2 =========================="<<fVars["nsigComb_K_2"]<<std::endl;    
    std::cout <<"nsigComb_Pr_2 =========================="<<fVars["nsigComb_Pr_2"]<<std::endl;    
    std::cout <<"sig_vert =========================="<<fVars["sig_vert"]<<std::endl;
    std::cout <<"max_norm_d0d0exp =========================="<<fVars["max_norm_d0d0exp"]<<std::endl;        
    /* std::cout <<"imp_par_prong0 =========================="<<fVars["imp_par_prong0"]<<std::endl;    
    std::cout <<"imp_par_prong1 =========================="<<fVars["imp_par_prong1"]<<std::endl;    
    std::cout <<"imp_par_prong2 =========================="<<fVars["imp_par_prong2"]<<std::endl; */
    
}
