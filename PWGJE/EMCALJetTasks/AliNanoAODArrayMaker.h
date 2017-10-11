/// \class AliNanoAODArrayMaker
/// \brief creates input arrays for the analysis from NanoAODs
///
/// The class converts the nanoAOD information into the standard arrays which are used in the PWGJE EMCAL analysis
/// this includes embedded tracks in measured events
///
/// \author M.Zimmermann

#ifndef ALIANALYSISTASKARRAYMAKER_H
#define ALIANALYSISTASKARRAYMAKER_H

class AliNanoAODTrack;
class TClonesArray;
class TString;

#include "AliAnalysisTaskSE.h"

class AliNanoAODArrayMaker : public AliAnalysisTaskSE {
 public:
  AliNanoAODArrayMaker() : AliAnalysisTaskSE(), fOutputArrayName(), fOutputArray(0), fOutputList(0x0) {}
  AliNanoAODArrayMaker(const char *name);
  virtual ~AliNanoAODArrayMaker() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetOutputArrayName(const char* name) {fOutputArrayName = name;}
  void SetOutputArrayPythiaName(const char* name) {fOutputArrayPythiaName = name;}
  void SetOutputArrayDataName(const char* name) {fOutputArrayDataName = name;}

  AliAODTrack*        GetAODTrack(AliNanoAODTrack* track, Int_t index);
  AliAODTrack*        GetAODTrack(AliNanoAODTrack* track);

 private:
  TString             fOutputArrayName; /// name of the output array with all particles
  TClonesArray*       fOutputArray; //!<! array with all particles

  TString             fOutputArrayPythiaName; /// name of the output array with pythia particles
  TClonesArray*       fPythiaArray; //!<! output array with pythia particles

  TString             fOutputArrayDataName; /// name of the output array with data particles
  TClonesArray*       fDataArray; //!<! output array with data particles
    
  TList* fOutputList;

  AliNanoAODArrayMaker(const AliNanoAODArrayMaker&); // not implemented
  AliNanoAODArrayMaker& operator=(const AliNanoAODArrayMaker&); // not implemented   

  
  ClassDef(AliNanoAODArrayMaker, 1); // example of analysis
};

#endif
