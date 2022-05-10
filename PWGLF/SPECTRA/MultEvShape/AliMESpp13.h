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
  class AliMESconfigTender : public TObject
  {
  public:
    friend class AliMESpp13;
    enum EMESconfigEventCuts
    {
      kNoEC = 0 // no event cuts
      ,
      k7TeV // vertex and trigger cuts for 7TeV
      ,
      k13TeV // vertex and trigger cuts for 13TeV
    };
    enum EMESconfigTrackCuts
    {
      kNoTC = 0 // no track cuts
      ,
      kStandardITSTPCTrackCuts2010 // 2010 standard cuts
      ,
      kStandardITSTPCTrackCuts2011 // 2011 standard cuts
    };
    enum EMESconfigPIDpriors
    {
      kNoPP = 0 // no priors
      ,
      kTPC // TPC priors
      ,
      kIterative // iterative priors
    };

    AliMESconfigTender();
    void Print(Option_t *o = "") const; // *MENU*

  protected:
    UChar_t fTrackCuts; // track cuts selector
    UChar_t fEventCuts; // event cuts selector
    UChar_t fPIDpriors; // PID prior selector

    ClassDef(AliMESconfigTender, 1)
  };

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

  virtual Bool_t ConfigTask(AliMESconfigTender::EMESconfigEventCuts ec,
                            AliMESconfigTender::EMESconfigTrackCuts tc,
                            AliMESconfigTender::EMESconfigPIDpriors pp);
  Bool_t HasMCdata() const { return TestBit(kMCdata); };
  virtual void SetMCdata(Bool_t mc = kTRUE);
  virtual void SetPriors();

  virtual Bool_t PostProcess();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void DeleteStreamer();

protected:
  Bool_t BuildQAHistos();

private:
  AliMESpp13(const AliMESpp13 &);
  AliMESpp13 &operator=(const AliMESpp13 &);

  AliMESconfigTender fConfig; // currrent configuration of task

  AliAnalysisFilter *fTrackFilter; // working track filter
  AliPIDCombined *fPIDcomb;        // working PID combined service

  TObjArray *fTracks;                 //!
  AliMESeventInfo *fEvInfo;           //!
  TObjArray *fMCtracks;               //!
  AliMESeventInfo *fMCevInfo;         //!
  TTreeSRedirector *fTreeSRedirector; //! temp tree to dump output
  TTree *fEventTree;                  //!
  TTree *fTracksTree;                 //!
  TTree *fMCeventTree;                //!
  TTree *fMCtracksTree;               //!
  TTree *fMCGenTracksTree;    //!
  TTree *fMCMissedTracksTree; //!

  AliPPVsMultUtils *fUtils;
  AliEventCuts fEventCutsQA; //!

  ClassDef(AliMESpp13, 5) // Tender task for the Multi Event Shape
};

#endif