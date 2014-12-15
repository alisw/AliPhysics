//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtModel.cc
//
// Description:
//
// Modification history:
//
//    RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include "EvtGenBase/EvtModel.hh"

#include "EvtGenModels/EvtTauVectornu.hh"
#include "EvtGenModels/EvtVVP.hh"
#include "EvtGenModels/EvtSLN.hh"
#include "EvtGenModels/EvtISGW2.hh"
#include "EvtGenModels/EvtMelikhov.hh"
#include "EvtGenModels/EvtSLPole.hh"
#include "EvtGenModels/EvtPropSLPole.hh"
#include "EvtGenModels/EvtSLBKPole.hh"
#include "EvtGenModels/EvtISGW.hh"
#include "EvtGenModels/EvtBHadronic.hh"
#include "EvtGenModels/EvtVSS.hh"
#include "EvtGenModels/EvtVSSMix.hh"
#include "EvtGenModels/EvtVSSBMixCPT.hh"
#include "EvtGenModels/EvtVSPPwave.hh"
#include "EvtGenModels/EvtGoityRoberts.hh"
#include "EvtGenModels/EvtSVS.hh"
#include "EvtGenModels/EvtTSS.hh"
#include "EvtGenModels/EvtTVSPwave.hh"
#include "EvtGenModels/EvtSVVHelAmp.hh"
#include "EvtGenModels/EvtSVPHelAmp.hh"
#include "EvtGenModels/EvtSVPCP.hh"
#include "EvtGenModels/EvtVVSPwave.hh"
#include "EvtGenModels/EvtDDalitz.hh"
#include "EvtGenModels/EvtOmegaDalitz.hh"
#include "EvtGenModels/EvtPi0Dalitz.hh"
#include "EvtGenModels/EvtEtaDalitz.hh"
#include "EvtGenModels/EvtPhsp.hh"
#include "EvtGenModels/EvtBtoXsgamma.hh"
#include "EvtGenModels/EvtBtoXsll.hh"
#include "EvtGenModels/EvtBtoXsEtap.hh"
#include "EvtGenModels/EvtSSSCP.hh"
#include "EvtGenModels/EvtSSSCPpng.hh"
#include "EvtGenModels/EvtSTSCP.hh"
#include "EvtGenModels/EvtSTS.hh"
#include "EvtGenModels/EvtSSSCPT.hh"
#include "EvtGenModels/EvtSVSCP.hh"
#include "EvtGenModels/EvtSSDCP.hh"
#include "EvtGenModels/EvtSVSNONCPEIGEN.hh"
#include "EvtGenModels/EvtSVVNONCPEIGEN.hh"
#include "EvtGenModels/EvtSVVCP.hh"
#include "EvtGenModels/EvtSVVCPLH.hh"
#include "EvtGenModels/EvtSVSCPLH.hh"
#include "EvtGenModels/EvtSll.hh"
#include "EvtGenModels/EvtVll.hh"
#include "EvtGenModels/EvtTaulnunu.hh"
#include "EvtGenModels/EvtTauHadnu.hh"
#include "EvtGenModels/EvtTauScalarnu.hh"
#include "EvtGenModels/EvtKstarnunu.hh"
#include "EvtGenModels/EvtbTosllBall.hh"
#include "EvtGenModels/EvtSingleParticle.hh"
#include "EvtGenModels/EvtVectorIsr.hh"
#include "EvtGenModels/EvtBToPlnuBK.hh"
#include "EvtGenModels/EvtBToVlnuBall.hh"

#include "EvtGenModels/EvtSVVHelCPMix.hh"
#include "EvtGenModels/EvtLb2Lll.hh"
#include "EvtGenModels/EvtHypNonLepton.hh"

#include "EvtGenModels/EvtbTosllAli.hh"
#include "EvtGenModels/EvtBToDDalitzCPK.hh"
#include "EvtGenModels/EvtPVVCPLH.hh"
#include "EvtGenModels/EvtLambdaB2LambdaV.hh"
#include "EvtGenModels/EvtSSD_DirectCP.hh"

#include "EvtGenModels/EvtHQET.hh"
#include "EvtGenModels/EvtHQET2.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtHelAmp.hh"
#include "EvtGenModels/EvtPartWave.hh"

#include "EvtGenModels/EvtBto2piCPiso.hh"
#include "EvtGenModels/EvtBtoKpiCPiso.hh"
#include "EvtGenModels/EvtSVSCPiso.hh"
#include "EvtGenModels/EvtVVpipi.hh"
#include "EvtGenModels/EvtY3SToY1SpipiMoxhay.hh"
#include "EvtGenModels/EvtYmSToYnSpipiCLEO.hh"
#include "EvtGenModels/EvtVVPIPI_WEIGHTED.hh"
#include "EvtGenModels/EvtVPHOtoVISRHi.hh"

#include "EvtGenModels/EvtBTo4piCP.hh"
#include "EvtGenModels/EvtBTo3piCP.hh"
#include "EvtGenModels/EvtCBTo3piP00.hh"
#include "EvtGenModels/EvtCBTo3piMPP.hh"
#include "EvtGenModels/EvtBToKpipiCP.hh"
#include "EvtGenModels/EvtBsquark.hh"
#include "EvtGenModels/EvtPhiDalitz.hh"

#include "EvtGenModels/EvtLNuGamma.hh"

#include "EvtGenModels/EvtVub.hh"
#include "EvtGenModels/EvtVubHybrid.hh"
#include "EvtGenModels/EvtVubNLO.hh"
#include "EvtGenModels/EvtVubBLNP.hh"
#include "EvtGenModels/EvtVubBLNPHybrid.hh"

#include "EvtGenModels/EvtPto3P.hh"
#include "EvtGenModels/EvtBtoKD3P.hh"
#include "EvtGenModels/EvtKstarstargamma.hh"
#include "EvtGenModels/EvtFlatQ2.hh"
#include "EvtGenModels/EvtLambdaP_BarGamma.hh"
#include "EvtGenModels/EvtBBScalar.hh"
#include "EvtGenModels/EvtKKLambdaC.hh"
#include "EvtGenModels/EvtMultibody.hh"

#include "EvtGenModels/EvtBaryonPCR.hh"

#include "EvtGenModels/EvtDMix.hh"
#include "EvtGenModels/EvtD0mixDalitz.hh"
#include "EvtGenModels/EvtD0gammaDalitz.hh"
#include "EvtGenModels/EvtEta2MuMuGamma.hh"

#include "EvtGenModels/EvtBcToNPi.hh"
#include "EvtGenModels/EvtBcPsiNPi.hh"
#include "EvtGenModels/EvtBcBsNPi.hh"
#include "EvtGenModels/EvtBcBsStarNPi.hh"

#include "EvtGenModels/EvtBcVMuNu.hh"
#include "EvtGenModels/EvtBcVNpi.hh"
#include "EvtGenModels/EvtSVP.hh"
#include "EvtGenModels/EvtTVP.hh"
#include "EvtGenModels/EvtXPsiGamma.hh"

#include "EvtGenModels/EvtGenericDalitz.hh"

#include "EvtGenModels/EvtModelReg.hh"
using std::fstream;
using std::cout;
using std::endl;

EvtModelReg::EvtModelReg(const std::list<EvtDecayBase*>* extraModels) 
{

    EvtModel &modelist=EvtModel::instance();

    if(extraModels){
      for(std::list<EvtDecayBase*>::const_iterator it = extraModels->begin(); 
	  it != extraModels->end(); ++it){
	modelist.registerModel(*it);
      }
    }

    modelist.registerModel(new EvtBBScalar);
    modelist.registerModel(new EvtLambdaP_BarGamma);
    modelist.registerModel(new EvtFlatQ2);
    modelist.registerModel(new EvtTauHadnu);
    modelist.registerModel(new EvtTauVectornu);
    modelist.registerModel(new EvtVVP);
    modelist.registerModel(new EvtSLN);
    modelist.registerModel(new EvtISGW2);
    modelist.registerModel(new EvtMelikhov);
    modelist.registerModel(new EvtSLPole);
    modelist.registerModel(new EvtPropSLPole);
    modelist.registerModel(new EvtSLBKPole);
    modelist.registerModel(new EvtHQET);
    modelist.registerModel(new EvtHQET2);
    modelist.registerModel(new EvtISGW);
    modelist.registerModel(new EvtBHadronic);
    modelist.registerModel(new EvtVSS);
    modelist.registerModel(new EvtVSSMix);
    modelist.registerModel(new EvtVSSBMixCPT);
    modelist.registerModel(new EvtVSPPwave);
    modelist.registerModel(new EvtGoityRoberts);
    modelist.registerModel(new EvtSVS);
    modelist.registerModel(new EvtTSS);
    modelist.registerModel(new EvtTVSPwave);
    modelist.registerModel(new EvtSVVHelAmp);
    modelist.registerModel(new EvtSVPHelAmp);
    modelist.registerModel(new EvtSVPCP);
    modelist.registerModel(new EvtVVSPwave);
    modelist.registerModel(new EvtDDalitz);
    modelist.registerModel(new EvtOmegaDalitz);
    modelist.registerModel(new EvtEtaDalitz);
    modelist.registerModel(new EvtPhsp);
    modelist.registerModel(new EvtBtoXsgamma);
    modelist.registerModel(new EvtBtoXsll);
    modelist.registerModel(new EvtBtoXsEtap);
    modelist.registerModel(new EvtSSSCP);
    modelist.registerModel(new EvtSSSCPpng);
    modelist.registerModel(new EvtSTSCP);
    modelist.registerModel(new EvtSTS);
    modelist.registerModel(new EvtSSSCPT);
    modelist.registerModel(new EvtSVSCP);
    modelist.registerModel(new EvtSSDCP);
    modelist.registerModel(new EvtSVSNONCPEIGEN);
    modelist.registerModel(new EvtSVVNONCPEIGEN);
    modelist.registerModel(new EvtSVVCP);
    modelist.registerModel(new EvtSVVCPLH);
    modelist.registerModel(new EvtSVSCPLH);
    modelist.registerModel(new EvtSll);
    modelist.registerModel(new EvtVll);
    modelist.registerModel(new EvtTaulnunu);
    modelist.registerModel(new EvtTauScalarnu);
    modelist.registerModel(new EvtKstarnunu);
    modelist.registerModel(new EvtbTosllBall);
    modelist.registerModel(new EvtBto2piCPiso);
    modelist.registerModel(new EvtBtoKpiCPiso);
    modelist.registerModel(new EvtSVSCPiso);
    modelist.registerModel(new EvtSingleParticle);
    modelist.registerModel(new EvtVectorIsr);
    modelist.registerModel(new EvtPi0Dalitz);
    modelist.registerModel(new EvtHelAmp);
    modelist.registerModel(new EvtPartWave);
    modelist.registerModel(new EvtVVpipi);
    modelist.registerModel(new EvtY3SToY1SpipiMoxhay);
    modelist.registerModel(new EvtYmSToYnSpipiCLEO);
    modelist.registerModel(new EvtBsquark);
    modelist.registerModel(new EvtPhiDalitz); 
    modelist.registerModel(new EvtBToPlnuBK); 
    modelist.registerModel(new EvtBToVlnuBall); 
    modelist.registerModel(new EvtVVPIPI_WEIGHTED);
    modelist.registerModel(new EvtVPHOtoVISRHi);

    modelist.registerModel(new EvtBTo4piCP);
    modelist.registerModel(new EvtBTo3piCP);
    modelist.registerModel(new EvtCBTo3piP00);
    modelist.registerModel(new EvtCBTo3piMPP);
    modelist.registerModel(new EvtBToKpipiCP);

    modelist.registerModel(new EvtLb2Lll);
    modelist.registerModel(new EvtHypNonLepton);
    modelist.registerModel(new EvtSVVHelCPMix);

    modelist.registerModel(new EvtLNuGamma);
    modelist.registerModel(new EvtKstarstargamma);

    modelist.registerModel(new EvtVub); 

    modelist.registerModel(new EvtVubHybrid);
    modelist.registerModel(new EvtVubNLO);
    modelist.registerModel(new EvtVubBLNP);
    modelist.registerModel(new EvtVubBLNPHybrid);

    modelist.registerModel(new EvtPto3P);
    modelist.registerModel(new EvtBtoKD3P);
    modelist.registerModel(new EvtKKLambdaC);
    modelist.registerModel(new EvtMultibody);
    modelist.registerModel(new EvtDMix);
    modelist.registerModel(new EvtD0mixDalitz);
    modelist.registerModel(new EvtD0gammaDalitz);
    modelist.registerModel(new EvtEta2MuMuGamma);

    modelist.registerModel(new EvtbTosllAli);
    modelist.registerModel(new EvtBaryonPCR);

    modelist.registerModel(new EvtBToDDalitzCPK);
    modelist.registerModel(new EvtLambdaB2LambdaV);
    modelist.registerModel(new EvtLambda2PPiForLambdaB2LambdaV);
    modelist.registerModel(new EvtV2VpVmForLambdaB2LambdaV);
    modelist.registerModel(new EvtPVVCPLH);
    modelist.registerModel(new EvtSSD_DirectCP);

    modelist.registerModel(new EvtBcToNPi(true)); // true = print author info
    modelist.registerModel(new EvtBcPsiNPi);
    modelist.registerModel(new EvtBcBsNPi);
    modelist.registerModel(new EvtBcBsStarNPi);

    modelist.registerModel(new EvtBcVMuNu);
    modelist.registerModel(new EvtBcVNpi);
    modelist.registerModel(new EvtSVP);
    modelist.registerModel(new EvtTVP);
    modelist.registerModel(new EvtXPsiGamma);

    modelist.registerModel(new EvtGenericDalitz());

}
