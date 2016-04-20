/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2015                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
 
 
  
/*
***********************************************************
    Event plane class that contains correction setting for specific event plane
***********************************************************
*/


#include "AliQnCorrectionsSteps.h"
#include "AliQnCorrectionsConfiguration.h"
#include "AliQnCorrectionsAxes.h"
#include "AliQnCorrectionsCuts.h"


ClassImp(AliQnCorrectionsConfiguration)


//_______________________________________________________________________________
AliQnCorrectionsConfiguration::AliQnCorrectionsConfiguration() :
  TObject(),
  fQnConfigurationCorrelationIndices(),
  fCuts(0x0),
  fQnNormalization(-1),
  fEqualizationMethod(-1),
  fTwistAndRescalingMethod(-1),
  fLocalIndex(-1),
  fGlobalIndex(-1),
  fCalibrationStep(-1),
  fMinimumHarmonic(2),
  fMaximumHarmonic(2),
  fAlignmentHarmonic(-1),
  fDetectorType(-1),
  fChannelList(0x0),
  fChannelGroups(0x0),
  fCalibrationDetectorNames(""),
  fEqualizationDetectorNames(""),
  fQnConfigurationName(""),
  fQnConfigurationCorrelationNames(),
  fAlignmentReferenceDetector(""),
  fCommonCorrectionAxes(0x0),
  fDataVectorEqualizationAxes(0x0),
  fRecenteringAxes(0x0),
  fAlignmentAxes(0x0),
  fTwistAndRescalingAxes(0x0),
  fEqualizationHistPath(""),
  fRecenteringHistPath(""),
  fCorrelationHistPath(""),
  fChannelEqualization(kFALSE),
  fRecenterQvec(kFALSE),
  fRotateQvec(kFALSE),
  fTwistQvec(kFALSE),
  fScaleQvec(kFALSE),
  fFillChannelEqualization(kFALSE),
  fFillRecenterQvec(kFALSE),
  fFillRotateQvec(kFALSE),
  fFillTwistQvec(kFALSE),
  fFillScaleQvec(kFALSE),
  fIsTracking(kFALSE),
  fUseLabel(kTRUE),
  fRequestedCorrectionMap(),
  fApplyCorrectionMap(),
  fRequestedHistogramMap(),
  fFillHistogramMap()
  //fCorrectionPasses()
{   
  //
  // Constructor
  //


    for(Int_t ic=0; ic<2; ++ic){
      fQnConfigurationCorrelationNames[ic]="";
      fQnConfigurationCorrelationIndices[ic]=-1;
    }


  for(Int_t is=0; is<AliQnCorrectionsConstants::nCorrectionSteps; is++) {fRequestedCorrectionMap[is]=kFALSE;fRequestedHistogramMap[is]=kFALSE;fApplyCorrectionMap[is]=kFALSE;fFillHistogramMap[is]=kFALSE;}


}



//_______________________________________________________________________________
AliQnCorrectionsConfiguration::~AliQnCorrectionsConfiguration()
{
  //
  // De-Constructor
  //

}


//_______________________________________________________________________________
AliQnCorrectionsConfiguration::AliQnCorrectionsConfiguration(const AliQnCorrectionsConfiguration &c) :
  TObject(),
  fCuts(c.fCuts),
  fQnNormalization(c.fQnNormalization),
  fEqualizationMethod(c.fEqualizationMethod),
  fTwistAndRescalingMethod(c.fTwistAndRescalingMethod),
  fLocalIndex(c.fLocalIndex),
  fGlobalIndex(c.fGlobalIndex),
  fCalibrationStep(c.fCalibrationStep),
  fMinimumHarmonic(c.fMinimumHarmonic),
  fMaximumHarmonic(c.fMaximumHarmonic),
  fAlignmentHarmonic(c.fAlignmentHarmonic),
  fDetectorType(c.fDetectorType),
  fChannelList(c.fChannelList),
  fChannelGroups(c.fChannelGroups),
  fCalibrationDetectorNames(c.fCalibrationDetectorNames),
  fEqualizationDetectorNames(c.fEqualizationDetectorNames),
  fQnConfigurationName(c.fQnConfigurationName),
  fAlignmentReferenceDetector(c.fAlignmentReferenceDetector),
  fCommonCorrectionAxes(c.fCommonCorrectionAxes),
  fDataVectorEqualizationAxes(c.fDataVectorEqualizationAxes),
  fRecenteringAxes(c.fRecenteringAxes),
  fAlignmentAxes(c.fAlignmentAxes),
  fTwistAndRescalingAxes(c.fTwistAndRescalingAxes),
  fEqualizationHistPath(c.fEqualizationHistPath),
  fRecenteringHistPath(c.fRecenteringHistPath),
  fCorrelationHistPath(c.fCorrelationHistPath),
  fChannelEqualization(c.fChannelEqualization),
  fRecenterQvec(c.fRecenterQvec),
  fRotateQvec(c.fRotateQvec),
  fTwistQvec(c.fTwistQvec),
  fScaleQvec(c.fScaleQvec),
  fFillChannelEqualization(c.fFillChannelEqualization),
  fFillRecenterQvec(c.fFillRecenterQvec),
  fFillRotateQvec(c.fFillRotateQvec),
  fFillTwistQvec(c.fFillTwistQvec),
  fFillScaleQvec(c.fFillScaleQvec),
  fIsTracking(c.fIsTracking),
  fUseLabel(c.fUseLabel)
{   
  //
  // Constructor
  //
  fQnConfigurationCorrelationNames[0] = c.QnConfigurationCorrelationName(0);
  fQnConfigurationCorrelationNames[1] = c.QnConfigurationCorrelationName(1);
  fQnConfigurationCorrelationIndices[0] = c.QnConfigurationCorrelationIndex(0);
  fQnConfigurationCorrelationIndices[1] = c.QnConfigurationCorrelationIndex(1);

  for(Int_t is=0; is<AliQnCorrectionsConstants::nCorrectionSteps; is++) {
    fRequestedCorrectionMap[is] = c.fRequestedCorrectionMap[is];
    fRequestedHistogramMap[is]  = c.fRequestedHistogramMap[is] ;
    fApplyCorrectionMap[is]     = c.fApplyCorrectionMap[is]    ;
    fFillHistogramMap[is]       = c.fFillHistogramMap[is]      ;
  }

}


//_______________________________________________________________________________
AliQnCorrectionsConfiguration & AliQnCorrectionsConfiguration::operator=(const AliQnCorrectionsConfiguration &c) {
  if (this == &c) return *this;
  else {
    fCuts=c.fCuts;
    fQnNormalization=c.fQnNormalization;
    fEqualizationMethod=c.fEqualizationMethod;
    fTwistAndRescalingMethod=c.fTwistAndRescalingMethod;
    fAlignmentHarmonic=c.fAlignmentHarmonic;
    fLocalIndex=c.fLocalIndex;
    fGlobalIndex=c.fGlobalIndex;
    fCalibrationStep=c.fCalibrationStep;
    fMinimumHarmonic=c.fMinimumHarmonic;
    fMaximumHarmonic=c.fMaximumHarmonic;
    fDetectorType=c.fDetectorType;
    fChannelList=c.fChannelList;
    fChannelGroups=c.fChannelGroups;
    fCalibrationDetectorNames=c.fCalibrationDetectorNames;
    fEqualizationDetectorNames=c.fEqualizationDetectorNames;
    fAlignmentReferenceDetector=c.fAlignmentReferenceDetector;
    fQnConfigurationName=c.fQnConfigurationName;
    fCommonCorrectionAxes=c.fCommonCorrectionAxes;
    fDataVectorEqualizationAxes=c.fDataVectorEqualizationAxes;
    fRecenteringAxes=c.fRecenteringAxes;
    fAlignmentAxes=c.fAlignmentAxes;
    fTwistAndRescalingAxes=c.fTwistAndRescalingAxes;
    fChannelEqualization=c.fChannelEqualization;
    fRecenterQvec=c.fRecenterQvec;
    fRotateQvec=c.fRotateQvec;
    fTwistQvec=c.fTwistQvec;
    fScaleQvec=c.fScaleQvec;
    fIsTracking=c.fIsTracking;
    fUseLabel=c.fUseLabel;

    fQnConfigurationCorrelationNames[0] = c.QnConfigurationCorrelationName(0);
    fQnConfigurationCorrelationNames[1] = c.QnConfigurationCorrelationName(1);
    fQnConfigurationCorrelationIndices[0] = c.QnConfigurationCorrelationIndex(0);
    fQnConfigurationCorrelationIndices[1] = c.QnConfigurationCorrelationIndex(1);

    for(Int_t is=0; is<AliQnCorrectionsConstants::nCorrectionSteps; is++) {
      fRequestedCorrectionMap[is] = c.fRequestedCorrectionMap[is];
      fRequestedHistogramMap[is]  = c.fRequestedHistogramMap[is] ;
      fApplyCorrectionMap[is]     = c.fApplyCorrectionMap[is]    ;
      fFillHistogramMap[is]       = c.fFillHistogramMap[is]      ;
    }

    return *this;
  }
    
}

//_______________________________________________________________________________
Bool_t AliQnCorrectionsConfiguration::PassCuts(Float_t* values) {
  return fCuts->IsSelected(values);
}


