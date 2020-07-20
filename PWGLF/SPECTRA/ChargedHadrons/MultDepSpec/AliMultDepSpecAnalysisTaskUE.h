/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis in the underlying event


#ifndef AliMultDepSpecAnalysisTaskUE_cxx
#define AliMultDepSpecAnalysisTaskUE_cxx
#include "AliMultDepSpecAnalysisTask.h"

class AliMultDepSpecAnalysisTaskUE : public AliMultDepSpecAnalysisTask
{
public:
  AliMultDepSpecAnalysisTaskUE();
  AliMultDepSpecAnalysisTaskUE(const char *name);
  virtual ~AliMultDepSpecAnalysisTaskUE();
  
  //virtual void   UserCreateOutputObjects() {AliMultDepSpecAnalysisTask::UserCreateOutputObjects();}
  //virtual void   UserExec(Option_t*);
  //virtual void   Terminate(Option_t*);

  AliMultDepSpecAnalysisTaskUE* AddTaskMultDepSpecUE(const std::string& dataSet, int cutModeLow, int cutModeHigh, TString options, bool isMC);
  
protected:

  virtual void DefineDefaultAxes(int maxMult = 100); // called in AddTask
  virtual void BookHistograms();    // called in UserCreateOutputObjects
  virtual bool InitEvent();
  virtual bool InitTrack(AliVTrack* track);
  virtual bool InitParticle(AliMCParticle* particle);     // called for ESDs
  virtual bool InitParticle(AliAODMCParticle* particle);  // called for AODs

  virtual bool SelectTrack();
  virtual bool SelectParticle();
    
  virtual void FillEventHistos();
  virtual void FillMeasTrackHistos();
  virtual void FillMeasParticleHistos();
  virtual void FillTrueParticleHistos();
  
private:
  
  AliMultDepSpecAnalysisTaskUE(const AliMultDepSpecAnalysisTaskUE&); // not implemented
  AliMultDepSpecAnalysisTaskUE& operator=(const AliMultDepSpecAnalysisTaskUE&); // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTaskUE, 1); // example of analysis
  /// \endcond
};

#endif
