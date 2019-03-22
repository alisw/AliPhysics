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

using std::vector;

class AliHFTreeHandlerD0toKpi : public AliHFTreeHandler
{
  public:
    
    enum isDzeroDzeroBar {
        kDzeroTopo     = BIT(11),
        kDzeroBarTopo  = BIT(12),
        kDzeroPID      = BIT(13),
        kDzeroBarPID   = BIT(14),
        kDzeroComb     = BIT(15),
        kDzeroBarComb  = BIT(16),
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
    virtual bool SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfiled, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    virtual void FillTree(); //to be called for each event, not each candidate!
    
    void SetIsDzeroDzeroBar(int isSel, int isSelTopo, int isSelPID, int isSelFilt, int isSelTopoFilt, int isSelPIDFilt);

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fCosThetaStar; /// vector of candidate costhetastar
    vector<float> fImpParProd; /// vector of daughter impact-parameter product
    vector<float> fNormd0MeasMinusExp; ///vector of candidate topomatic variable
    vector<float> fImpParErrProng[knMaxProngs]; ///vector of error on prongs rphi impact param [cm]
    vector<float> fNormDecayLength; ///vector of candidate normalised decay length

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerD0toKpi,2); /// 
    /// \endcond
};
#endif
