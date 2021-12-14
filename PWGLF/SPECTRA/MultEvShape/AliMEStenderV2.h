#ifndef ALIMESTENDERV2_H
#define ALIMESTENDERV2_H

////////////////////////////////////////////////////////////////////////////
//  Tender for Multiplicity and Event Shape group                         //
//  Tender configuration                                                  //
//  Authors:                                                              //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
//    Andrei Herghelegiu <aherghe@niham.nipne.ro>                         //
//    Madalina Tarzila <mtarzila@niham.nipne.ro>                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIANALYSISTASKSE_H
#include <AliAnalysisTaskSE.h>
#endif

class AliAnalysisFilter;
class AliPIDCombined;
//class AliESDevent;
class AliMCEvent;
class TObjArray;
class TClonesArray;
class AliMESeventInfo;
class AliAnalysisUtils;
#include "AliEventCuts.h"
class AliMEStenderV2 : public AliAnalysisTaskSE
{
public:
  class AliMESconfigTender : public TObject
  {
  public:
    friend class AliMEStenderV2;
    enum EMESconfigEventCuts{
      kNoEC = 0     // no event cuts
			,k7TeV    // vertex and trigger cuts for 7TeV
      ,k13TeV       // vertex and trigger cuts for 13TeV
    };
    enum EMESconfigTrackCuts{
      kNoTC = 0                      // no track cuts
			,kStandardITSTPCTrackCuts2010  // 2010 standard cuts
      ,kStandardITSTPCTrackCuts2011  // 2011 standard cuts
    };
    enum EMESconfigPIDpriors{
      kNoPP = 0                    // no priors
      ,kTPC                        // TPC priors
      ,kIterative                  // iterative priors
    };

    AliMESconfigTender();
    void    Print(Option_t *o="") const;   // *MENU*

  protected:
    UChar_t fTrackCuts;    // track cuts selector
    UChar_t fEventCuts;    // event cuts selector
    UChar_t fPIDpriors;    // PID prior selector

    ClassDef(AliMESconfigTender, 1)
  };

  //_______________________________________________________________________________________
  enum AliMEStenderV2Steering{
     kMCdata      = BIT(18)     // MC presence bit
    ,kPP          = BIT(19)     // pp/pA data
    ,kPostProcess = BIT(20)     // run pos processing of QA histos
  };
  enum EMEStenderV2QA{
     kConfig = 0
    ,kEfficiency
    ,kEvInfo
    ,kTrkInfo
    ,kMCevInfo
    ,kMCtrkInfo
    ,kNqa
  };
  AliMEStenderV2();
  AliMEStenderV2(const char *name);
  virtual ~AliMEStenderV2();

  //static Int_t    MakeMultiplicityESD(AliESDEvent* const, const char *opt);
  static Int_t    MakeMultiplicityMC(AliMCEvent* const);
  static Int_t    MakeMultiplicity0408MC(AliMCEvent* const);
  static Int_t    MakeMultiplicityV0MMC(AliMCEvent* const);

  virtual Bool_t  ConfigTask(AliMESconfigTender::EMESconfigEventCuts ec,
                             AliMESconfigTender::EMESconfigTrackCuts tc,
                             AliMESconfigTender::EMESconfigPIDpriors pp);
  Bool_t          HasMCdata() const       { return TestBit(kMCdata);};
  virtual void    SetDebugLevel(Int_t level);
  virtual void    SetMCdata(Bool_t mc = kTRUE);
  virtual void    SetPriors();

  virtual Bool_t  PostProcess();
  virtual void    UserCreateOutputObjects();
  virtual void    UserExec(Option_t *opt);

protected:
  Bool_t          BuildQAHistos();

private:
  AliMEStenderV2(const AliMEStenderV2&);
  AliMEStenderV2& operator=(const AliMEStenderV2&);

  AliMESconfigTender  fConfig;       // currrent configuration of task

  AliAnalysisFilter  *fTrackFilter;  // working track filter
  AliPIDCombined     *fPIDcomb;      // working PID combined service

  TObjArray *fTracks;       //!
  TClonesArray *fTracksIO;       //!
  AliMESeventInfo *fEvInfo; //!
  TObjArray *fMCtracks;  //!
  TClonesArray *fMCtracksIO; //!
  AliMESeventInfo *fMCevInfo;  //!
  TClonesArray *fMCGenTracksIO;   //!
  TClonesArray *fMCtracksMissIO; //!

  AliPPVsMultUtils *fUtils;
  AliEventCuts fEventCutsQA; //!

  ClassDef(AliMEStenderV2, 5)          // Tender task for the Multi Event Shape
};

#endif
