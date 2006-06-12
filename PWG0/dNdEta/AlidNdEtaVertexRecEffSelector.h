/* $Id$ */

#ifndef ALIDNDETAVERTEXRECEFFSELECTOR_H
#define ALIDNDETAVERTEXRECEFFSELECTOR_H

// This class plots the vertex reconstruction efficiency

#include "AliSelectorRL.h"

class TH1F;

class AlidNdEtaVertexRecEffSelector : public AliSelectorRL {
  public:
    AlidNdEtaVertexRecEffSelector();
    virtual ~AlidNdEtaVertexRecEffSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    static const Float_t fkEtaRange;

    Bool_t CheckVertex();

    TH1F* fdNGen;  //! generated multiplicity
    TH1F* fdNRec;  //! generated multiplicity of events with reconstructed vertex

    TH1F* fVtxGen;  //! generated vertex z 
    TH1F* fVtxRec;  //! generated vertex z of events with reconstructed vertex

 private:

  ClassDef(AlidNdEtaVertexRecEffSelector, 0);
};

#endif
