#ifndef ALIMESBASETASK_H
#define ALIMESBASETASK_H

////////////////////////////////////////////////////////////////////////////
//  Base task for Multiplicity and Event Shape group                      //
//  Authors:                                                              //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
//    Andrei Herghelegiu <aherghe@niham.nipne.ro>                         //
//    Madalina Tarzila <mtarzila@niham.nipne.ro>                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIANALYSISTASKSE_H
#include <AliAnalysisTaskSE.h>
#endif

class AliMEStender;
class AliMESeventInfo;
class TObjArray;
class TTreeSRedirector;
class AliMESbaseTask : public AliAnalysisTaskSE
{
public:
  friend class AliMEStender;
  enum AliMESbaseSteering{
     kMCdata      = BIT(18)     // MC presence bit
    ,kPP          = BIT(19)     // pp/pA data 
    ,kPostProcess = BIT(20)     // run pos processing of QA histos
  };
  enum EMEScontainers{
     kQA        = 1
    ,kEventInfo = 1
    ,kTracks
    ,kMCeventInfo
    ,kMCtracks
    ,kNcontainers
  };

  AliMESbaseTask();
  virtual ~AliMESbaseTask();
  static Int_t   DebugUsers()            { return fgDebugUsers; }
  Bool_t         IsPP() const            { return TestBit(kPP);};
  Bool_t         HasMCdata() const       { return TestBit(kMCdata);};
  Bool_t         HasPostProcess() const  { return TestBit(kPostProcess);};

  virtual void   SetDebugLevel(Int_t level);
  virtual void   SetMCdata(Bool_t mc = kTRUE);

  virtual void   SetPostProcess(Bool_t pp = kTRUE) { SetBit(kPostProcess, pp);}
  virtual void   SetPP(Bool_t pp = kTRUE)          { SetBit(kPP, pp);}

protected:
  AliMESbaseTask(const char *name);
  virtual Bool_t BuildQAHistos();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual Bool_t PostProcess() = 0;
  static TTreeSRedirector* DebugStream() { return fgDebugStream;}  
  static void              CloseDebugStream();
  static void              OpenDebugStream();
  static void              AddDebugUser(const char *name);
  
  AliMESeventInfo    *fEvInfo;       // event info for reconstructed event
  TObjArray          *fTracks;       // array of reconstructed tracks/event
  AliMESeventInfo    *fMCevInfo;     // event info for MC event
  TObjArray          *fMCtracks;     // array of MC tracks/event
  
private:
  AliMESbaseTask(const AliMESbaseTask&);
  AliMESbaseTask& operator=(const AliMESbaseTask&);
  static Int_t             fgDebugUsers;   //! no. of user tasks for the debug stream
  static TTreeSRedirector *fgDebugStream;  //! Debug stream
  
  ClassDef(AliMESbaseTask, 1)        // Basic Analisys task for the Multi Event Shape
};

#endif

