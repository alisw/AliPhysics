/* $Id$ */

#ifndef AliVertexSelector_H
#define AliVertexSelector_H

#include "AliSelectorRL.h"

class TH1F;
class TH2F;
class AliESDtrackCuts;

class AliVertexSelector : public AliSelectorRL {
  public:
    AliVertexSelector();
    virtual ~AliVertexSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 private:
    AliVertexSelector(const AliVertexSelector&);
    AliVertexSelector& operator=(const AliVertexSelector&);

    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

    TH1F* fVertexMC;   // MC vertex distribution
    TH1F* fVertexESD;  // ESD vertex distribution
    TH2F* fVertexCorr; // vertex correlation (esd vtx - mc vtx vs. mc vtx)
    TH2F* fVertexCorr2; // vertex correlation (esd vtx - mc vtx vs. mc vtx)

    TH2F* fFakes; // counts the fakes

  ClassDef(AliVertexSelector, 0);
};

#endif
