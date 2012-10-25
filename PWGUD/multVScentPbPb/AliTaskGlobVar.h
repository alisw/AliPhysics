#ifndef ALITASKGLOBVAR_H
#define ALITASKGLOBVAR_H

///////////////////////////////////////////////////////////////////////////
// Analysis task to extract global variables to the tree                 //
///////////////////////////////////////////////////////////////////////////

class AliESDEvent;
class TList;
class AliESDtrackCuts;
class AliTriggerAnalysis;
#include "AliAnalysisTaskSE.h"


typedef struct {
  enum {kTDCNA=0x1,kTDCPA=0x1<<1,kTDCNC=0x1<<2,kTDCPC=0x1<<3,kSPDVTXOK=0x1<<4};
  Int_t   runID;
  UInt_t  timeStamp;
  Float_t zdcNA;
  Float_t zdcPA;
  Float_t zdcNC;
  Float_t zdcPC;
  Float_t zdcNAC;//common PM
  Float_t zdcNCC;//common PM
  Float_t zem1;
  Float_t zem2;
  //
  Float_t zvSPD;
  Float_t zvTPC;
  Short_t chunk;
  Short_t flags;
  Short_t spd1;
  Short_t spd2;
  Short_t ncontSPDV;
  Short_t ncontTPCV;
  Short_t nTrTPC;
  Short_t nTrTPCITS;
  Short_t nTracklets;
  Short_t v0A;
  Short_t v0C;
  Short_t v0Corr;
  //  Short_t v0CorrResc;
  Float_t mcZV;
  Short_t mcdNdEta;
  Short_t mcNPart;
  Short_t mcNBColl;
} GloVars_t;


class AliTaskGlobVar : public AliAnalysisTaskSE {
  //
 public:
  AliTaskGlobVar(const char *name = "AliTaskGlobVar");
  virtual ~AliTaskGlobVar(); 
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  //
  void    SetUseMC(Bool_t mc=kTRUE)                                   {fUseMC = mc;}
  Float_t GetCorrV0(const AliESDEvent* esd, float &v0CorrResc) const;
  Bool_t  ZDCTimeTrigger(const AliESDEvent *aEsd) const;
  AliESDtrackCuts* CreatedNdPtTrackCuts(Int_t cutMode, Bool_t fieldOn=kTRUE);
  //
 protected:
  Bool_t       fUseMC;                    // do we use MC info
  TList*       fOutput;                   // output list send on output slot 1 
  //
  TTree*       fOutTree;                   // output tree
  AliESDtrackCuts* fTrackCuts;             //! optional track cuts
  AliESDtrackCuts* fTrackCuts1;            //! optional track cuts
  //
  GloVars_t fGlobVars;                     // data container
 private:    
  AliTaskGlobVar(const AliTaskGlobVar&); // not implemented
  AliTaskGlobVar& operator=(const AliTaskGlobVar&); // not implemented 
  
  ClassDef(AliTaskGlobVar, 1);  
};


#endif
