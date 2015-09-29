#ifndef ALIMESTENDER_H
#define ALIMESTENDER_H

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
class AliMESeventInfo;
class AliAnalysisUtils;
class AliMEStender : public AliAnalysisTaskSE
{
public:
  class AliMESconfigTender : public TObject
  {
  public:
    friend class AliMEStender;
    enum EMESconfigEventCuts{
        kNoEC = 0           // no event cuts
      ,kStandard           // vertex and trigger cuts
    };
    enum EMESconfigTrackCuts{
        kNoTC = 0                     // no track cuts
      ,kStandardITSTPCTrackCuts2010  // 2010 standard cuts
    };
    enum EMESconfigPIDpriors{
        kNoPP = 0                     // no priors
      ,kTPC                          // TPC priors
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
  enum AliMEStenderSteering{
     kMCdata      = BIT(18)     // MC presence bit
    ,kPP          = BIT(19)     // pp/pA data
    ,kPostProcess = BIT(20)     // run pos processing of QA histos
  };
  enum EMEStenderQA{
     kConfig = 0
    ,kEfficiency
    ,kEvInfo
    ,kTrkInfo
    ,kMCevInfo
    ,kMCtrkInfo
    ,kNqa
  };
  AliMEStender();
  AliMEStender(const char *name);
  virtual ~AliMEStender();

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
  AliMEStender(const AliMEStender&);
  AliMEStender& operator=(const AliMEStender&);

  AliMESconfigTender  fConfig;       // currrent configuration of task

  AliAnalysisFilter  *fTrackFilter;  // working track filter
  AliPIDCombined     *fPIDcomb;      // working PID combined service

  TObjArray *fTracks;  //!
  AliMESeventInfo *fEvInfo;  //!
  TObjArray *fMCtracks;  //!
  AliMESeventInfo *fMCevInfo;  //!

  AliPPVsMultUtils *fUtils;

  ClassDef(AliMEStender, 5)          // Tender task for the Multi Event Shape
};

#endif
