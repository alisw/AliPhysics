#ifndef ALIANALYSISTASKITSTPCALIGNMENT_H
#define ALIANALYSISTASKITSTPCALIGNMENT_H

///////////////////////////////////////////////////////////////////////////
//  Class AliAnalysisTaskITSTPCalignment
//  runs ITS-TPC alignment with TPC vdrift calibration
//
//    Origin: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
///////////////////////////////////////////////////////////////////////////

#include<AliAnalysisTask.h>
class TTree;
class AliESDEvent;
class AliRelAlignerKalman;
class AliRelAlignerKalmanArray;
class TH2F;

class AliAnalysisTaskITSTPCalignment : public AliAnalysisTask
{
public:
  AliAnalysisTaskITSTPCalignment();
  AliAnalysisTaskITSTPCalignment(const char *name);
  virtual ~AliAnalysisTaskITSTPCalignment() {}

  void SetSaveInterval( const UInt_t si ) {fSaveInterval = si;}
  void SetTimeMatchingTolerance( const UInt_t tol ) {fTimeMatchingTolerance = tol; }
  void SetDoQA( const Bool_t qa ) {fDoQA=qa;}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

private:
  AliESDEvent* fESD;    //ESD object
  AliRelAlignerKalmanArray* fArray; //array of aligners
  TH2F* fYZResidualsHist; //2D histogram with the yz residuals
  TH2F* fPLResidualsHist; //2D histogram with the phi lambda residuals
  TList* fListOfHistos;   //list with QA histograms
  UInt_t fSaveInterval;   //save interveal
  UInt_t fTimeMatchingTolerance; //time matching tolerance

  Bool_t fDoQA; //whether to fill QA histograms

  AliAnalysisTaskITSTPCalignment(const AliAnalysisTaskITSTPCalignment&); // not implemented
  AliAnalysisTaskITSTPCalignment& operator=(const AliAnalysisTaskITSTPCalignment&); // not implemented

  ClassDef(AliAnalysisTaskITSTPCalignment, 1);
};

#endif
