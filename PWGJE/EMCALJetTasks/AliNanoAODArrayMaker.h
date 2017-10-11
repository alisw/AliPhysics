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
  TString             fOutputArrayName; ///
  TClonesArray*       fOutputArray; //!<!

  TString             fOutputArrayPythiaName; ///
  TClonesArray*       fPythiaArray; //!<!

  TString             fOutputArrayDataName; ///
  TClonesArray*       fDataArray; //!<!
    
  TList* fOutputList;

  AliNanoAODArrayMaker(const AliNanoAODArrayMaker&); // not implemented
  AliNanoAODArrayMaker& operator=(const AliNanoAODArrayMaker&); // not implemented   

  
  ClassDef(AliNanoAODArrayMaker, 1); // example of analysis
};

#endif
