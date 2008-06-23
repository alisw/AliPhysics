#ifndef AliRsnBaseAT_cxx
#define AliRsnBaseAT_cxx

#include <TChain.h>

#include "AliAnalysisTask.h"

class AliAnalysisManager;

class AliESDEvent;
class AliAODEvent;
class AliRsnEvent;


class AliRsnBaseAT : public AliAnalysisTask
{
  public:
    AliRsnBaseAT ( const char *name = "AliRsnBaseAT" );
    virtual ~AliRsnBaseAT() {}

    enum EInputType
    {
      kAOD = 0,
      kESD,
      kESDMC,
      kMC,
      kRSN,
      kLastIndex
    };

    virtual void   InitIOVars();
    virtual void   LocalInit() {;}
    virtual Bool_t Notify();
    virtual void   ConnectInputData ( Option_t * );
    virtual void   CreateOutputObjects() {;}
    virtual void   Exec ( Option_t *option ) {;}
    virtual void   Terminate ( Option_t * ) {;}

    void SetInputType (EInputType& theValue, Short_t inputIndex=0) { fInputType[inputIndex] = theValue; }
    EInputType GetInputType ( Short_t inputIndex=0 ) { return fInputType[inputIndex]; }

    TChain* GetChain ( const Int_t& index = 0) const { return fChain[index]; }

    AliRsnEvent *GetRSNEvent ( Int_t index=0 ) { return fRSN[index]; }

    void SetAnalysisMgr ( AliAnalysisManager* theValue ) { fAnalysisMgr = theValue; }
    AliAnalysisManager* GetAnalysisMgr() const { return fAnalysisMgr; }

  protected:
    
    Long64_t      fNumOfEvents;

    TChain        *fChain[2];         // input chain
    EInputType    fInputType[2];      // input type

    AliRsnEvent *fRSN[2];           // RsnMV event
    AliESDEvent   *fESD[2];           // ESD event
    AliAODEvent   *fAOD[2];           // AOD event

    AliAnalysisManager *fAnalysisMgr; // pointer to current AnalysisMgr

    virtual void  ConnectInputDataByInputType ( EInputType type ,Short_t inputIndex=0 );
    virtual void  ConnectRSN ( Short_t inputIndex );
    virtual void  ConnectESD ( Short_t inputIndex );
    virtual void  ConnectESDMC ( Short_t inputIndex );
    virtual void  ConnectAOD ( Short_t inputIndex );

    ClassDef ( AliRsnBaseAT, 1 );
};

#endif
