#ifndef ALIMESPP13_H
#define ALIMESPP13_H

////////////////////////////////////////////////////////////////////////////
//  Task for multiplicity and event shape studies in pp @ 13 TeV          //
//                                                                        //
//  Authors:                                                              //
//    Amelia Lindner <amelia.lindner@cern.ch>                             //
//    Alex Bercuci (a.bercuci@gsi.de)                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIANALYSISTASKSE_H
#include <AliAnalysisTaskSE.h>
#endif

class AliAnalysisFilter;
class AliPIDCombined;
// class AliESDevent;
class AliMCEvent;
class TObjArray;
class TClonesArray;
class AliMESeventInfo;
class AliAnalysisUtils;
class TTree;
class TTreeSRedirector;

#include "AliEventCuts.h"
#include "TTreeStream.h"

class AliMESpp13 : public AliAnalysisTaskSE
{
public:
  //_______________________________________________________________________________________
  enum AliMESpp13Steering
  {
    kMCdata = BIT(18) // MC presence bit
    ,
    kPP = BIT(19) // pp/pA data
    ,
    kPostProcess = BIT(20) // run pos processing of QA histos
  };
  enum EMESpp13QA
  {
    kEfficiency,
    kEvInfo,
    kTrkInfo,
    kMCevInfo,
    kMCtrkInfo,
    kNqa
  };
  enum MESpp13Containers
  {
    kQA = 1,
    kTree = 1,
    kMCGenTree,
    // kMCMissTree,
    kNcontainers
  };
  AliMESpp13();
  AliMESpp13(const char *name);
  virtual ~AliMESpp13();

  Bool_t HasMCdata() const { return TestBit(kMCdata); };
  virtual void SetMCdata(Bool_t mc = kTRUE);
  virtual void FinishTaskOutput();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);

protected:
  Bool_t BuildQAHistos();

private:
  AliMESpp13(const AliMESpp13 &);
  AliMESpp13 &operator=(const AliMESpp13 &);

  AliAnalysisFilter *fTrackFilter;    // working track filter
  TObjArray *fTracks;                 //!
  AliMESeventInfo *fEvInfo;           //!
  TObjArray *fMCtracks;               //!
  AliMESeventInfo *fMCevInfo;         //!
  TTreeSRedirector *fTreeSRedirector; //! temp tree to dump output
  TClonesArray *fTracksIO;            //!
  TClonesArray *fMCtracksIO;          //!
  TClonesArray *fMCGenTracksIO;       //!
  TClonesArray *fMCtracksMissIO;      //!
  TTree *fTree;                       //!
  TTree *fMCGenTree;                  //!
  TTree *fMCMissTree;                 //!

  AliPPVsMultUtils *fUtils;
  AliEventCuts fEventCutsQA; //!

  ClassDef(AliMESpp13, 3) // MESpp task
};

#endif
