#ifndef ALIHFMLVARHANDLERD0TOKPI_H
#define ALIHFMLVARHANDLERD0TOKPI_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandlerD0toKpi
// \brief helper class to handle a tree and variables for D+ ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"

class AliHFMLVarHandlerD0toKpi : public AliHFMLVarHandler
{
    public:

        enum massHypo {
            kD0,
            kD0bar
        };

        AliHFMLVarHandlerD0toKpi();
        AliHFMLVarHandlerD0toKpi(int PIDopt);
        virtual ~AliHFMLVarHandlerD0toKpi();

        AliHFMLVarHandlerD0toKpi(const AliHFMLVarHandlerD0toKpi &source) = delete;
        AliHFMLVarHandlerD0toKpi& operator=(const AliHFMLVarHandlerD0toKpi &source) = delete;

        virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
        virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliAODPidHF *pidrespo = nullptr);

    private:
        float fImpParProng[knMaxProngs] = {-999., -999., -999., -999.};  /// prong impact parameter
        float fNormd0MeasMinusExp = -999.;                               /// candidate topomatic variable
        float fImpParProd = -999.;                                       /// d0K*d0pi product
        float fCosThetaStar = -999.;                                     /// cosine of theta star

        /// \cond CLASSIMP
        ClassDef(AliHFMLVarHandlerD0toKpi, 1); /// 
        /// \endcond
};
#endif
