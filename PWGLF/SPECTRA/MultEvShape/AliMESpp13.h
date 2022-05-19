#ifndef ALIMESPP13_H
#define ALIMESPP13_H

////////////////////////////////////////////////////////////////////////////
//  Tender for Multiplicity and Event Shape group                         //
//  Tender configuration                                                  //
//  Authors:                                                              //
//    Amelia Lindner <amelia.lindner@cern.ch>                             //
//    Alex Bercuci (a.bercuci@gsi.de)                                     //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
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
    kConfig = 0,
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
    kEventInfo = 1,
    kTracks,
    kEventTree,
    kTracksTree,
    kMCeventInfo,
    kMCtracks,
    kMCeventTree,
    kMCtracksTree,
    kMCGenTracksTree,
    kMCMissedTracksTree,
    kNcontainers
  };
  AliMESpp13();
  AliMESpp13(const char *name);
  virtual ~AliMESpp13();

  // static Int_t    MakeMultiplicityESD(AliESDEvent* const, const char *opt);
  static Int_t MakeMultiplicityMC(AliMCEvent *const);
  static Int_t MakeMultiplicity0408MC(AliMCEvent *const);
  static Int_t MakeMultiplicityV0MMC(AliMCEvent *const);
  static Double_t ComputeDeltaPhi(Double_t, Double_t);

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

  AliAnalysisFilter *fTrackFilter; // working track filter
  TObjArray *fTracks;                 //!
  AliMESeventInfo *fEvInfo;           //!
  TObjArray *fMCtracks;               //!
  AliMESeventInfo *fMCevInfo;         //!
  TTreeSRedirector *fTreeSRedirector; //! temp tree to dump output
  TTree *fEventTree;                  //!
  TTree *fTracksTree;                 //!
  TTree *fMCeventTree;                //!
  TTree *fMCtracksTree;               //!
  TTree *fMCGenTracksTree;            //!
  TTree *fMCMissedTracksTree;         //!
  AliPPVsMultUtils *fUtils;
  AliEventCuts fEventCutsQA; //!

  ClassDef(AliMESpp13, 3) // MESpp task
};

#endif