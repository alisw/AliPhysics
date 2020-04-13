#ifndef ALIHFMLVARHANDLERDPLUSTOKPIPI_H
#define ALIHFMLVARHANDLERDPLUSTOKPIPI_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandlerDplustoKKpi
// \brief helper class to handle a tree and variables for D+ ML analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"

class AliHFMLVarHandlerDplustoKpipi : public AliHFMLVarHandler
{
    public:
        AliHFMLVarHandlerDplustoKpipi();
        AliHFMLVarHandlerDplustoKpipi(int PIDopt);
        virtual ~AliHFMLVarHandlerDplustoKpipi();

        AliHFMLVarHandlerDplustoKpipi(const AliHFMLVarHandlerDplustoKpipi &source) = delete;
        AliHFMLVarHandlerDplustoKpipi& operator=(const AliHFMLVarHandlerDplustoKpipi &source) = delete;

        virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
        virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliAODPidHF *pidrespo = nullptr);

    private:
        float fImpParProng[knMaxProngs] = {-999., -999., -999., -999.};  ///prong impact parameter
        float fSigmaVertex = -999.;                                      /// candidate sigma vertex
        float fNormd0MeasMinusExp = -999.;                               ///candidate topomatic variable

        /// \cond CLASSIMP
        ClassDef(AliHFMLVarHandlerDplustoKpipi, 1); /// 
        /// \endcond
};
#endif
