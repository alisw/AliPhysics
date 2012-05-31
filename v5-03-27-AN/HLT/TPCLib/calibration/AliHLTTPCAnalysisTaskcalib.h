//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCANALYSISTASKCALIB_H
#define ALIHLTTPCANALYSISTASKCALIB_H

#include "AliTPCAnalysisTaskcalib.h"
#include "TObjArray.h"
#include "AliTPCcalibBase.h"
#include "AliExternalTrackParam.h"

class AliESDEvent;
class AliESDtrack;
class AliESDfriend;
class AliTPCseed;

class AliHLTTPCAnalysisTaskcalib : public AliTPCAnalysisTaskcalib
{
public:
  /** constructor */
  AliHLTTPCAnalysisTaskcalib();
  AliHLTTPCAnalysisTaskcalib(const char *name);
 
  /** destructor */
  virtual ~AliHLTTPCAnalysisTaskcalib();
 
  void AddJob(AliTPCcalibBase *job) {fCalibJobs->Add(job);}
  TObjArray* GetJobs() {return fCalibJobs;}
  

  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void FinishTaskOutput();
  void         SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}

//protected:
  virtual void     Process(AliESDEvent *event);
  virtual void     Process(AliTPCseed *track);
  virtual void     Process(AliESDtrack *track, Int_t run);
  virtual Long64_t Merge(TCollection *li);
  virtual void     Analyze();
  void             RegisterDebugOutput();

private:
  TObjArray    *fCalibJobs;   // array of calibration objects - WE ARE NOT OWNER?
  AliESDEvent  *fESD;         //! current esd
  AliESDfriend *fESDfriend;   //! current esd friend
  TString       fDebugOutputPath; // debug output path   
  
  /** copy constructor prohibited */
  AliHLTTPCAnalysisTaskcalib(const AliHLTTPCAnalysisTaskcalib&); 
  /** assignment operator prohibited */
  AliHLTTPCAnalysisTaskcalib& operator=(const AliHLTTPCAnalysisTaskcalib&);
  
  ClassDef(AliHLTTPCAnalysisTaskcalib,0)
};

#endif
