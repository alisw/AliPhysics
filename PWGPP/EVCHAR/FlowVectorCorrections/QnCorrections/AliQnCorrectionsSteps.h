#ifndef ALIQNCORRECTIONS_STEPS_H
#define ALIQNCORRECTIONS_STEPS_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2015                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/
 
 

#include <TClonesArray.h>

class AliQnCorrectionsHistograms;
class AliQnCorrectionsConfiguration;
class AliQnCorrectionsQnVector;


//_________________________________________________________________________
class AliQnCorrectionsSteps  {
  
  

 public: 

 AliQnCorrectionsSteps();
 virtual ~AliQnCorrectionsSteps();


  //enum CorrectionSteps {
  //  kNothing=-1,
  //  kPass0,
  //  kDataVectorEqualization,
  //  kRecentering,
  //  kAlignment,
  //  kTwist,
  //  kRescaling,
  //  kNcorrectionSteps,
  //  kUndefined,
  //  kCopy
  //};


 //static void BuildQnVectors(AliQnCorrectionsQnVector* QvectorOut, TClonesArray* dataVectorArray, Int_t QnConfIndex, Int_t minHar, Int_t maxHar, Int_t EqualizationMethod=-1);
 static void CalibrateDataVector(TClonesArray* dataVectorArray, AliQnCorrectionsConfiguration* QnConf,  AliQnCorrectionsHistograms* inputHistos, Double_t* fillValues) ;
 static void RecenterQvec(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorOut, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t useStep, Int_t minHar, Int_t maxHar) ;
 static void TwistAndRescale2nQn(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorTwist, AliQnCorrectionsQnVector* QvectorRescale, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Bool_t doTwist, Bool_t doRescaling);
 static void TwistAndRescale3DetectorCorrelation(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorTwist, AliQnCorrectionsQnVector* QvectorRescale, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Bool_t doTwist, Bool_t doRescaling, Int_t eventClassParameter);



  //void CalibrateChannels(Bool_t useFriendChain=kFALSE);
  //void RecenterQvec(Bool_t useFriendChain=kFALSE);
  //void RotateQvec();
  //void CorrelationTwistAndRescalingQvec(Float_t* values, Int_t corpar);
  //void TwoDetectorCorrelationTwistQvec(Float_t* values, Int_t corpar);
  //void CorrelationRescalingQvec(Float_t* values, Int_t corpar);
  //void ThreeDetectorCorrelationTPCTwistQvec(Float_t* values, Int_t corpar);
  //void ThreeDetectorCorrelationTPCRescalingQvec(Float_t* values, Int_t corpar);
  //void ThreeDetectorCorrelationTPCTwistAndRescalingQvec(Float_t* values, Int_t corpar);
  //void 2nTwistQvec(Float_t* values, Int_t corpar);
  //void 2nRescalingQvec(Float_t* values, Int_t corpar);
  //void U2nTwistAndRescalingQvec(Float_t* values, Int_t corpar);


 private:

  ClassDef(AliQnCorrectionsSteps, 1);
};



#endif
