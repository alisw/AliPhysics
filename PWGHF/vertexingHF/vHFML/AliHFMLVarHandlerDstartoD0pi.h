#ifndef ALIHFMLVARHANDLERDSTARTOD0PI_H
#define ALIHFMLVARHANDLERDSTARTOD0PI_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliHFMLVarHandlerDstartoD0pi
// \brief helper class to handle a tree and variables for D*+ ML analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"

class AliHFMLVarHandlerDstartoD0pi : public AliHFMLVarHandler
{
    public:
        AliHFMLVarHandlerDstartoD0pi();
        AliHFMLVarHandlerDstartoD0pi(int PIDopt);
        virtual ~AliHFMLVarHandlerDstartoD0pi();

        AliHFMLVarHandlerDstartoD0pi(const AliHFMLVarHandlerDstartoD0pi &source) = delete;
        AliHFMLVarHandlerDstartoD0pi& operator=(const AliHFMLVarHandlerDstartoD0pi &source) = delete;

        virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
        virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliAODPidHF *pidrespo = nullptr);

    private:
        float fImpParProng[knMaxProngs] = {-999., -999., -999., -999.};  /// prong impact parameter
        float fCosThetaStar; ///candidate cos theta star
        float fImpParProd; ///D0 product of impact parameter
        float fNormd0MeasMinusExp; ///candidate topomatic variable
        float fAngleD0dkpPisoft; ///angle between D0 decay plane and soft pion
        float fDeltaInvMassD0; /// delta mass of D0

        /// \cond CLASSIMP
        ClassDef(AliHFMLVarHandlerDstartoD0pi, 1); /// 
        /// \endcond
};
#endif
