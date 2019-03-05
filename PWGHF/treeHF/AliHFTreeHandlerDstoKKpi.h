#ifndef ALIHFTREEHANDLERDSTOKKPI_H
#define ALIHFTREEHANDLERDSTOKKPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDstoKKpi
// \brief helper class to handle a tree for Ds cut optimisation and MVA analyses
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

class AliHFTreeHandlerDstoKKpi : public AliHFTreeHandler
{
  public:

    enum massKKopt {kMassKK,kDeltaMassKKPhi};

    static const int kDplustoKKpi = BIT(11);

    AliHFTreeHandlerDstoKKpi();
    AliHFTreeHandlerDstoKKpi(int PIDopt);

    virtual ~AliHFTreeHandlerDstoKKpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    virtual void FillTree();

    void SetMassKKOption(int opt) {fMassKKOpt=opt;}
    void SetIsDplustoKKpi(bool isDplus) {
      if(isDplus) fCandTypeMap |= kDplustoKKpi;
      else fCandTypeMap &= ~kDplustoKKpi;
    }
  
    static bool IsDplustoKKpi(int candtype) {
      if(candtype>>11&1) return true;
      return false;
    }

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fSigmaVertex; /// vector of candidate sigma vertex
    vector<float> fMassKK; /// vector of candidate massKK
    vector<float> fCosPiDs; /// vector of candidate cos3piDs
    vector<float> fCosPiKPhi; /// vector of candidate cospiKphi
    vector<float> fNormd0MeasMinusExp; ///vector of candidate topomatic variable
    int fMassKKOpt; /// option for massKK variable (mass or delta mass wrt phi)

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerDstoKKpi,1); /// 
    /// \endcond
};
#endif
