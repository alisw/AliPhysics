#ifndef ALIANALYSISTASKNOISETPC_H
#define ALIANALYSISTASKNOISETPC_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis flags rare noise events in the TPC.                        //
//                                                                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2F;
class TH3F;
class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliESDpid;


#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisNoiseTPC : public AliAnalysisTaskSE {
 public:
  AliAnalysisNoiseTPC(const char *name,  UInt_t StartTime, UInt_t EndTime, Int_t deltaTime);
  AliAnalysisNoiseTPC();
  virtual ~AliAnalysisNoiseTPC() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  //
  //
  
 private:
  //
  //
  AliESDEvent *fESD;                  //! ESD object
  TList       *fListHist;             //! list for histograms
  //
  AliESDtrackCuts * fESDtrackCuts;    // basic cut variables
  //
  //
  //
  THnSparseS * fHistNoiseTracks;      //! histogram with all necessary information for real tracks

  //
  Int_t    fTimeBins;                   //Bins time
  Double_t fTimeStart;                  //Start time
  Double_t fTimeEnd;                    //End time

  //
  AliAnalysisNoiseTPC(const AliAnalysisNoiseTPC&); 
  AliAnalysisNoiseTPC& operator=(const AliAnalysisNoiseTPC&); 

  ClassDef(AliAnalysisNoiseTPC, 1); 
};

#endif
