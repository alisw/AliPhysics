#ifndef ALIHFTREEHANDLERD0TOKPI_H
#define ALIHFTREEHANDLERD0TOKPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerD0toKpi
// \brief helper class to handle a tree for D0 cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

class AliHFTreeHandlerD0toKpi : public AliHFTreeHandler
{
  public:
    
    enum isDzeroDzeroBar {
        kDzeroTopo         = BIT(11),
        kDzeroBarTopo      = BIT(12),
        kDzeroPID          = BIT(13),
        kDzeroBarPID       = BIT(14),
        kDzeroComb         = BIT(15),
        kDzeroBarComb      = BIT(16),
        kDzeroTopoFilt     = BIT(17),
        kDzeroBarTopoFilt  = BIT(18),
        kDzeroPIDFilt      = BIT(19),
        kDzeroBarPIDFilt   = BIT(20),
        kDzeroCombFilt     = BIT(21),
        kDzeroBarCombFilt  = BIT(22)
    };
    
    AliHFTreeHandlerD0toKpi();
    AliHFTreeHandlerD0toKpi(int PIDopt);

    virtual ~AliHFTreeHandlerD0toKpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfiled, int masshypo=0, AliPIDResponse *pidrespo=nullptr);
    
    void SetIsDzeroDzeroBar(int isSel, int isSelTopo, int isSelPID, int isSelFilt, int isSelTopoFilt, int isSelPIDFilt);

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; /// candidate costhetastar
    float fImpParProd; /// daughter impact-parameter product
    float fNormd0MeasMinusExp; ///candidate topomatic variable
    float fImpParErrProng[knMaxProngs]; ///error on prongs rphi impact param [cm]
    float fNormDecayLength; ///candidate normalised decay length

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerD0toKpi,3); /// 
    /// \endcond
};
#endif
