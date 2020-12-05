#ifndef ALIHFMLVARHANDLERLCTOPKPI_H
#define ALIHFMLVARHANDLERLCTOPKPI_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandlerLctopKpi
// \brief helper class to handle a tree and variables for Ds ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"

class AliHFMLVarHandlerLctopKpi : public AliHFMLVarHandler
{
    public:

        enum massHypo {
            kpKpi,
            kpiKp
        };

        enum restype {
            kLcNonRes     = BIT(9),
            kLcLambda1520 = BIT(10),
            kLcKStar      = BIT(11),
            kLcDelta      = BIT(12)
        };



        AliHFMLVarHandlerLctopKpi();
        AliHFMLVarHandlerLctopKpi(int PIDopt);
        virtual ~AliHFMLVarHandlerLctopKpi();

        AliHFMLVarHandlerLctopKpi(const AliHFMLVarHandlerLctopKpi &source) = delete;
        AliHFMLVarHandlerLctopKpi& operator=(const AliHFMLVarHandlerLctopKpi &source) = delete;

        virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
        virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliAODPidHF *pidrespo = nullptr);
        void SetIsLcpKpiRes(int restype);

    private:
        float fImpParProng[knMaxProngs] = {-999., -999., -999., -999.};  /// prong impact parameter
        float fSigmaVertex = -999.;                                      /// candidate sigma vertex
        float fNormd0MeasMinusExp = -999.;                               /// candidate topomatic variable

        /// \cond CLASSIMP
        ClassDef(AliHFMLVarHandlerLctopKpi, 1); /// 
        /// \endcond
};
#endif
