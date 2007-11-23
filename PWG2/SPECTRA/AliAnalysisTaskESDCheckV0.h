#ifndef AliAnalysisTaskESDCheckV0_cxx
#define AliAnalysisTaskESDCheckV0_cxx

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskESDCheckV0 class
//            This task is for QAing the V0s from the ESD
//              Origin: B.H. Nov2007, hippolyt@in2p3.fr
//-----------------------------------------------------------------

class TList;
class TH1F;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskESDCheckV0 : public AliAnalysisTask {
 public:
  AliAnalysisTaskESDCheckV0(const char *name = "AliAnalysisTaskESDCheckV0");
  virtual ~AliAnalysisTaskESDCheckV0() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;                            //! ESD object
  TList       *fListHist;                       //! List of histograms
  TH1F        *fHistTrackMultiplicity;          //! Track multiplicity distribution
  TH1F        *fHistV0Multiplicity;             //! V0 multiplicity distribution
  TH1F        *fHistV0OnFlyStatus;              //! V0 on fly status distribution

              // V0 offline distributions
  TH1F        *fHistV0MultiplicityOff;          //! V0 multiplicity distribution offline
  TH1F        *fHistV0Chi2Off;                  //! V0 chi2 distribution
  TH1F        *fHistDcaV0DaughtersOff;          //! Dca between V0 daughters
  TH1F        *fHistV0CosineOfPointingAngleOff; //! Cosine of V0 pointing angle
  TH1F        *fHistV0RadiusOff;                //! V0 radial distance distribution
  TH1F        *fHistDcaV0ToPrimVertexOff;       //! Dca of V0 to primary vertex
  TH1F        *fHistDcaPosToPrimVertexOff;      //! Dca of V0 positive daughter to primary vertex
  TH1F        *fHistDcaNegToPrimVertexOff;      //! Dca of V0 negative daughter to primary vertex

  TH1F        *fHistMassK0Off;                  //! Invariant Mass of K0s
  TH1F        *fHistMassLambdaOff;              //! Invariant Mass of Lambda
  TH1F        *fHistMassAntiLambdaOff;          //! Invariant Mass of Anti-Lambda

              // V0 on-the-fly distributions
  TH1F        *fHistV0MultiplicityOn;           //! V0 multiplicity distribution offline
  TH1F        *fHistV0Chi2On;                   //! V0 chi2 distribution
  TH1F        *fHistDcaV0DaughtersOn;           //! Dca between V0 daughters
  TH1F        *fHistV0CosineOfPointingAngleOn;  //! Cosine of V0 pointing angle
  TH1F        *fHistV0RadiusOn;                 //! V0 radial distance distribution
  TH1F        *fHistDcaV0ToPrimVertexOn;        //! Dca of V0 to primary vertex
  TH1F        *fHistDcaPosToPrimVertexOn;       //! Dca of V0 positive daughter to primary vertex
  TH1F        *fHistDcaNegToPrimVertexOn;       //! Dca of V0 negative daughter to primary vertex

  TH1F        *fHistMassK0On;                   //! Invariant Mass of K0s
  TH1F        *fHistMassLambdaOn;               //! Invariant Mass of Lambda
  TH1F        *fHistMassAntiLambdaOn;           //! Invariant Mass of Anti-Lambda
   
  AliAnalysisTaskESDCheckV0(const AliAnalysisTaskESDCheckV0&);            // not implemented
  AliAnalysisTaskESDCheckV0& operator=(const AliAnalysisTaskESDCheckV0&); // not implemented
  
  ClassDef(AliAnalysisTaskESDCheckV0, 1);
};

#endif
