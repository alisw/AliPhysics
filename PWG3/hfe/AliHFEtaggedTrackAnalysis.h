/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class AliHFEtaggedTrackAnalysis
// Analyses tracks with an apriori PID information (i.e. using the daugther
// tracks from well-identified decays of neutral charged particles).
// More information can be found inside the implementation file.
//
#ifndef ALIHFETAGGEDTRACKANALYSIS_H
#define ALIHFETAGGEDTRACKANALYSIS_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliHFEvarManager;

class AliHFEtaggedTrackAnalysis : public TObject{
  public:
    AliHFEtaggedTrackAnalysis();
    AliHFEtaggedTrackAnalysis(const AliHFEtaggedTrackAnalysis &ref);
    AliHFEtaggedTrackAnalysis &operator=(const AliHFEtaggedTrackAnalysis &ref);
    ~AliHFEtaggedTrackAnalysis();
    
    void InitContainer();
    void ProcessTrack(AliVParticle *track, Int_t abinitioPID);

    AliHFEcontainer *GetContainer() const { return fContainer; }
    AliHFEpidQAmanager *GetPIDqa() const { return fPIDqa; }
    TList * GetPIDQA() const;
    TList * GetCutQA() const;
    Bool_t  GetClean() const { return fClean; }; 

    void SetCuts(AliHFEcuts *cuts);
    void SetPID(AliHFEpid *pid);
    void SetClean(Bool_t clean) { fClean = clean; };

  private:
    enum{
      kIsOwner = BIT(14),
      kIsOwnerCuts = BIT(15)
    };
    AliHFEvarManager    *fVarManager;   // Variable Manager
    AliHFEcontainer     *fContainer;    // Output container
    AliHFEpid           *fPID;          // PID selection
    AliHFEpidQAmanager  *fPIDqa;        // PID monitoring
    AliHFEcuts          *fCuts;         // Single track cuts
    AliCFManager        *fCFM;          // CF Manager used for the track filtering
    Bool_t               fClean;        // Clean
    
  ClassDef(AliHFEtaggedTrackAnalysis, 0)
};
#endif
