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


 static void CalibrateDataVector(TClonesArray* dataVectorArray, AliQnCorrectionsConfiguration* QnConf,  AliQnCorrectionsHistograms* inputHistos, Double_t* fillValues) ;
 static void RecenterQvec(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorOut, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t useStep, Int_t minHar, Int_t maxHar) ;
 static void TwistAndRescale2nQn(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorTwist, AliQnCorrectionsQnVector* QvectorRescale, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Bool_t doTwist, Bool_t doRescaling);
 static void TwistAndRescale3DetectorCorrelation(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorTwist, AliQnCorrectionsQnVector* QvectorRescale, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Bool_t doTwist, Bool_t doRescaling, Int_t eventClassParameter);
 static void RotateQvec(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorOut, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Int_t alignmentHarmonic);

 private:

  ClassDef(AliQnCorrectionsSteps, 1);
};



#endif
