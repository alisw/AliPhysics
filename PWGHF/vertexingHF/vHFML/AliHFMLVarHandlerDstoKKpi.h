#ifndef ALIHFMLVARHANDLERDSTOKKPI_H
#define ALIHFMLVARHANDLERDSTOKKPI_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandlerDstoKKpi
// \brief helper class to handle a tree and variables for Ds ML analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"

class AliHFMLVarHandlerDstoKKpi : public AliHFMLVarHandler
{
    public:

        enum massKKopt {kDeltaMassKKPhi,kMassKK};
        static const int kDplustoKKpi = BIT(9);

        AliHFMLVarHandlerDstoKKpi();
        AliHFMLVarHandlerDstoKKpi(int PIDopt);
        virtual ~AliHFMLVarHandlerDstoKKpi();

        AliHFMLVarHandlerDstoKKpi(const AliHFMLVarHandlerDstoKKpi &source) = delete;
        AliHFMLVarHandlerDstoKKpi& operator=(const AliHFMLVarHandlerDstoKKpi &source) = delete;

        virtual TTree* BuildTree(TString name="tree", TString title="tree");
        virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliAODPidHF *pidrespo=nullptr);

        void SetMassKKOption(int opt) {fMassKKOpt=opt;}
        void SetIsDplustoKKpi(bool isDplus) {
            if(isDplus) 
                fCandType |= kDplustoKKpi;
            else 
                fCandType &= ~kDplustoKKpi;
        }

    private:
        float fImpParProng[knMaxProngs]; ///prong impact parameter
        float fSigmaVertex; /// candidate sigma vertex
        float fMassKK; /// candidate massKK
        float fCosPiDs; /// candidate cos3piDs
        float fCosPiKPhi; /// candidate cospiKphi
        float fNormd0MeasMinusExp; ///candidate topomatic variable
        int fMassKKOpt; /// option for massKK variable (mass or delta mass wrt phi)

        /// \cond CLASSIMP
        ClassDef(AliHFMLVarHandlerDstoKKpi,1); /// 
        /// \endcond
};
#endif
