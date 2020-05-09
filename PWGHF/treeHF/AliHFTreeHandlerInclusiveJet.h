#ifndef ALIHFTREEHANDLERINCLUSIVE_H
#define ALIHFTREEHANDLERINCLUSIVE_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerInclusiveJet
// \brief helper class to handle a tree for inclusive analyses
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

class AliHFTreeHandlerInclusiveJet : public AliHFTreeHandler
{
  public:
    
    
    AliHFTreeHandlerInclusiveJet();

    virtual ~AliHFTreeHandlerInclusiveJet();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo, AliPIDResponse *pidrespo);
    virtual bool SetMCGenVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long);


  private:


    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerInclusiveJet,1); /// 
    /// \endcond
};
#endif
