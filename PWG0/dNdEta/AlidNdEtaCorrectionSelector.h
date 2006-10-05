/* $Id$ */

#ifndef ALIDNDETACORRECTIONSELECTOR_H
#define ALIDNDETACORRECTIONSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class AlidNdEtaCorrection;
class TH1F;
class TParticlePDG;

class AlidNdEtaCorrectionSelector : public AliSelectorRL {
  public:
    AlidNdEtaCorrectionSelector();
    virtual ~AlidNdEtaCorrectionSelector();

    void ReadUserObjects(TTree* tree);
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    Bool_t SignOK(TParticlePDG* particle);

    AliESDtrackCuts*  fEsdTrackCuts;          // Object containing the parameters of the esd track cuts

    AlidNdEtaCorrection* fdNdEtaCorrection;      // contains the intermediate histograms (on each slave)

    TH1F* fPIDParticles; // pid of primary particles
    TH1F* fPIDTracks; // pid of reconstructed tracks

    TH1F* fClustersITSPos; //
    TH1F* fClustersTPCPos; //

    TH1F* fClustersITSNeg; //
    TH1F* fClustersTPCNeg; //

    Int_t fSignMode;  // if 0 process all particles, if +-1 process only particles with that sign

 private:
    AlidNdEtaCorrectionSelector(const AlidNdEtaCorrectionSelector&);
    AlidNdEtaCorrectionSelector& operator=(const AlidNdEtaCorrectionSelector&);

  ClassDef(AlidNdEtaCorrectionSelector, 0);
};

#endif
