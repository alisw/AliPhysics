#ifndef AliRsnAnalysisAT_h
#define AliRsnAnalysisAT_h

#include "AliRsnBaseAT.h"

class AliRsnEvent;

/**
  @author Martin Vala <Martin.Vala@cern.ch>
*/
class AliRsnAnalysisAT : public AliRsnBaseAT
{
  public:
    AliRsnAnalysisAT ( const char *name = "AliRsnAnalysisAT");

    ~AliRsnAnalysisAT();

    virtual void    InitIOVars();
    virtual void    LocalInit();
    virtual Bool_t  Notify();
    virtual void    CreateOutputObjects();
    virtual void    Exec ( Option_t *option );
    virtual void    Terminate ( Option_t * );
    virtual void    Cleanup ();

  private:
    AliRsnPairMgr *fPairMgr;

    TList           *fOutList;              // List of output
    TH1F            *fHist[100][100];       // output histograms
    
    AliRsnEventBuffer *fRsnMVEventBuffer;
    
    void            ProcessEventAnalysis(AliRsnEvent *curEvent);
    void            PostEventProcess(const Short_t &index=0);
    AliRsnEvent*  GetRsnMVEventFromInputType(const Short_t &index=0);
    
    AliRsnEvent*  GetRsnMVFromAOD(const Short_t &index=0);
    AliRsnEvent*  GetRsnMVFromESD(const Short_t &index=0);
    AliRsnEvent*  GetRsnMVFromESDMC(const Short_t &index=0);
    AliRsnEvent*  GetRsnMVFromRSN(const Short_t &index=0);
    
    ClassDef ( AliRsnAnalysisAT, 1 );
};

#endif
