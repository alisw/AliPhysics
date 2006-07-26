/* $Id$ */

#ifndef ALIDNDETASYSTEMATICSSELECTOR_H
#define ALIDNDETASYSTEMATICSSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class AlidNdEtaCorrection;
class TH3F;
class TH1F;

class AlidNdEtaSystematicsSelector : public AliSelectorRL {
  public:
    AlidNdEtaSystematicsSelector();
    virtual ~AlidNdEtaSystematicsSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    void FillCorrectionMaps(TObjArray* listOfTracks);
    void FillSecondaries(TObjArray* listOfTracks);
    void FillSigmaVertex();

    TH3F* fSecondaries; // (accepted tracks) vs (tracks from sec)/(n * tracks from sec) vs pT
    AlidNdEtaCorrection* fdNdEtaCorrection[4];      // correction for different particle species: here pi, K, p, others
    TH1F* fSigmaVertex; // (accepted tracks) vs (n of sigma to vertex cut)

    AliESDtrackCuts* fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

 private:

  ClassDef(AlidNdEtaSystematicsSelector, 0);
};

#endif
