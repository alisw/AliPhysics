#ifndef ALIQNCORRECTIONS_CONFIGURATION_H
#define ALIQNCORRECTIONS_CONFIGURATION_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2015                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/
 
 
 

#include <Rtypes.h>
#include <TArrayS.h>
#include "AliQnCorrectionsAxes.h"
#include "AliQnCorrectionsSteps.h"

class AliQnCorrectionsAxes;
class AliQnCorrectionsCuts;
class TObject;

//_____________________________________________________________________
class AliQnCorrectionsConfiguration : public TObject {


 public:
  AliQnCorrectionsConfiguration();
  AliQnCorrectionsConfiguration(const AliQnCorrectionsConfiguration &c);
  AliQnCorrectionsConfiguration & operator=(const AliQnCorrectionsConfiguration &c);
  ~AliQnCorrectionsConfiguration();
    

  // setters for users
  //void SetCorrectionSteps(AliQnCorrectionsSteps::CorrectionSteps flag1, AliQnCorrectionsSteps::CorrectionSteps flag2=AliQnCorrectionsConstants::kCopy, AliQnCorrectionsSteps::CorrectionSteps flag3=AliQnCorrectionsConstants::kCopy, AliQnCorrectionsSteps::CorrectionSteps flag4=AliQnCorrectionsConstants::kCopy, AliQnCorrectionsSteps::CorrectionSteps flag5=AliQnCorrectionsConstants::kCopy, AliQnCorrectionsSteps::CorrectionSteps flag6=AliQnCorrectionsConstants::kCopy); // which correction has to be applied at which step
  void SetQnCorrectionDataVectorEqualization(Bool_t set=kTRUE) {SetMaps(AliQnCorrectionsConstants::kDataVectorEqualization,set);}
  void SetQnCorrectionRecentering(Bool_t set=kTRUE) {SetMaps(AliQnCorrectionsConstants::kRecentering,set);}
  void SetQnCorrectionTwist(Bool_t set=kTRUE) {SetMaps(AliQnCorrectionsConstants::kTwist,set);}
  void SetQnCorrectionRescaling(Bool_t set=kTRUE) {SetMaps(AliQnCorrectionsConstants::kRescaling,set);}
  void SetQnCorrectionAlignment(Bool_t set=kTRUE) {SetMaps(AliQnCorrectionsConstants::kAlignment,set);}

  void SetQnNormalization(UShort_t QnNormalization)      {  fQnNormalization     = QnNormalization;}
  void SetDataVectorEqualizationMethod(Float_t equalizationMethod)    {  fEqualizationMethod    = equalizationMethod;}
  void SetTwistAndRescalingMethod(Float_t method)                     {  fTwistAndRescalingMethod = method;}
  void SetHarmonicForAlignment(Int_t har)                                {  fAlignmentHarmonic=har;}
  
//TODO Add possibility to enable an individual harmonic (need to check feasibility for CPU etc)
  void SetQnHarmonicsRange(Int_t minHar, Int_t maxHar)             { fMinimumHarmonic = minHar;fMaximumHarmonic = maxHar;}
  void SetDataVectorIdList(TArrayS * list)                       { fChannelList     = list;} 
  void SetDataVectorEqualizationGroups(TArrayS * list)                     { fChannelGroups   = list;} 
  void SetQnConfigurationName(TString name)                 { fQnConfigurationName=name;}
  void SetReferenceQnForTwistAndRescaling(TString detA, TString detB)   {fQnConfigurationCorrelationNames[0]  =detA;fQnConfigurationCorrelationNames[1]  =detB;}
  void SetReferenceQnForAlignment(TString det)                         {fQnConfigurationCorrelationNames[2]  =det;} 
  
  void SetQnCorrectionsCommonAxes(AliQnCorrectionsAxes* axes)       { if(!fDataVectorEqualizationAxes) fDataVectorEqualizationAxes=axes; if(!fRecenteringAxes) fRecenteringAxes=axes; if(!fAlignmentAxes) fAlignmentAxes=axes; if(!fTwistAndRescalingAxes) fTwistAndRescalingAxes=axes;}
  void SetDataVectorEqualizationAxes(AliQnCorrectionsAxes* axes)    { fDataVectorEqualizationAxes= axes ;}
  void SetRecenteringAxes(AliQnCorrectionsAxes* axes)               { fRecenteringAxes= axes ;}
  void SetAlignmentAxes(AliQnCorrectionsAxes* axes)                 { fAlignmentAxes= axes ;}
  void SetTwistAndRescalingAxes(AliQnCorrectionsAxes* axes)         { fTwistAndRescalingAxes= axes ;}
//TODO Add switch to SetCorrelationAxes for Twist/Rescaling/Rotation
  void SetDataVectorCuts(AliQnCorrectionsCuts* trackCuts)                { fCuts       = trackCuts;}

//TODO Extend QA of the DataContainer variables before and after cuts
  void SetTracking(Bool_t set=kTRUE)                        {fIsTracking         =set;}  // obsolete?



  // INTERNAL FRAMEWORK FUNCTIONS
  void SetLocalIndex(Int_t index)                           { fLocalIndex      = index;}
  void SetGlobalIndex(Int_t index)                          { fGlobalIndex     = index;}
  //void SetQnConfigurationCorrelationIndices(Int_t detA, Int_t detB)     {fQnConfigurationCorrelationIndices[0]=detA;fQnConfigurationCorrelationIndices[1]=detB;}
  void SetQnConfigurationCorrelationIndex(Int_t det, Int_t cor)         {fQnConfigurationCorrelationIndices[det]=cor;}
  void SetDetectorId(UShort_t detectorType)                 { fDetectorType    = detectorType;} // SetDetectorId-> SetDetectorId
  void SetFillHistogram(Int_t step, Bool_t b=kTRUE)   {fFillHistogramMap[step]=b;}
  void SetApplyCorrection(Int_t step, Bool_t b=kTRUE) {fApplyCorrectionMap[step]=b;}
  void SetMinimumHarmonic(Int_t har)                  {fMinimumHarmonic=har;}





  //void SetPassNumbers();
  //void SetPassCorrections(Int_t pass, AliQnCorrectionsSteps::CorrectionSteps flag) { fCorrectionPasses[(Int_t) flag]=pass;};
  //void DisableCorrection(AliQnCorrectionsSteps::CorrectionSteps flag) { fCorrectionPasses[(Int_t) flag]=-1;};
  //void EnableCorrection(AliQnCorrectionsSteps::CorrectionSteps flag) { fCorrectionPasses[(Int_t) flag]=-1;};

  //void SetChannelEqualization(Bool_t set) {fChannelEqualization=set;}
  //void SetRecentering(Bool_t set)         {fRecenterQvec       =set;}
  //void SetRotation(Bool_t set)            {fRotateQvec         =set;}
  //void SetTwist(Bool_t set)               {fTwistQvec          =set;}
  //void SetScaling(Bool_t set)             {fScaleQvec          =set;}

  
  void UseCalibrationDirectoryNameAllEvents(Bool_t b=kFALSE)      {fUseLabel=(b ? kFALSE : kTRUE);}

  //void SetCorrectionFlag(AliQnCorrectionsSteps::CorrectionSteps flag);  // which correction has to be applied at some point

  //void EnableChannelEqualization(Bool_t b=kTRUE)    {fChannelEqualization=b;}
  //void EnableRecentering(Bool_t b=kTRUE)            {fRecenterQvec=b;}
  //void EnableRotation(Bool_t b=kTRUE)               {fRotateQvec=b;}
  //void EnableTwist(Bool_t b=kTRUE)                  {fTwistQvec=b;}
  //void EnableScaling(Bool_t b=kTRUE)                {fScaleQvec=b;}

  //void EnableFillChannelEqualization(Bool_t b=kTRUE)   {fFillChannelEqualization=b;}
  //void EnableFillRecentering(Bool_t b=kTRUE)           {fFillRecenterQvec=b;}
  //void EnableFillRotation(Bool_t b=kTRUE)              {fFillRotateQvec=b;}
  //void EnableFillTwist(Bool_t b=kTRUE)                 {fFillTwistQvec=b;}
  //void EnableFillScaling(Bool_t b=kTRUE)               {fFillScaleQvec=b;}

  // OBSOLETE
  void SetCalibrationDetectorName(TString name)     { fCalibrationDetectorNames =name;}  // obsolete
  void SetEqualizationDetectorName(TString name)    { fEqualizationDetectorNames=name;}  // obsolete

  // getters
  UShort_t GetQnNormalizationMethod()       const {return fQnNormalization ;}
  Short_t GetDataVectorEqualizationMethod()       const {return fEqualizationMethod ;}
  Short_t GetTwistAndRescalingMethod()    const {return fTwistAndRescalingMethod ;}
  
  Int_t LocalIndex()                 const {return fLocalIndex ;}
  Int_t GlobalIndex()                const {return fGlobalIndex ;}
  Int_t CalibrationStep()            const {return fCalibrationStep ;}
  UShort_t DetectorType()            const {return fDetectorType ;}
  UShort_t UseChannel(Int_t ch)      const {return fChannelList->At(ch);}
  TArrayS* ChannelList()             const {return fChannelList;} 
  TArrayS* ChannelGroups()           const {return fChannelGroups;}
  Int_t ChannelGroup(Int_t i)        const {return fChannelGroups->At(i);}
  AliQnCorrectionsCuts* Cuts()          const {return fCuts;}
  Bool_t PassCuts(Float_t* values);
  TString CalibrationDetectorName()  const {return fCalibrationDetectorNames ;}
  TString EqualizationDetectorName() const {return fEqualizationDetectorNames ;}
  TString QnConfigurationName()      const {return fQnConfigurationName;}
  TString QnConfigurationCorrelationName(Int_t detCor) const {return fQnConfigurationCorrelationNames[detCor];}
  TString GetReferenceQnForAlignment()  const  {return fQnConfigurationCorrelationNames[2];} 
  Int_t QnConfigurationCorrelationIndex(Int_t detCor)  const {return fQnConfigurationCorrelationIndices[detCor];}
  Int_t MinimumHarmonic()            const {return fMinimumHarmonic;}
  Int_t MaximumHarmonic()            const {return fMaximumHarmonic;}
  Int_t AlignmentHarmonic()          const {return fAlignmentHarmonic;}

  AliQnCorrectionsAxes* EqualizationBinning()  const {return fDataVectorEqualizationAxes;}
  AliQnCorrectionsAxes* CalibrationBinning()   const {return fCommonCorrectionAxes; }
  AliQnCorrectionsAxes* GetRecenteringAxes()           const {return fRecenteringAxes; }
  AliQnCorrectionsAxes* GetAlignmentAxes()           const {return fAlignmentAxes; }
  AliQnCorrectionsAxes* GetTwistAndRescalingAxes()   const {return fTwistAndRescalingAxes; }


  Bool_t IsEnabledChannelEqualization()   const {return fChannelEqualization;}
  Bool_t IsEnabledRecentering()           const {return fRecenterQvec;}
  Bool_t IsEnabledRotation()              const {return fRotateQvec;}
  Bool_t IsEnabledTwist()                 const {return fTwistQvec;}
  Bool_t IsEnabledScaling()               const {return fScaleQvec;}
  Bool_t IsEnabledFillChannelEqualization()   const {return fFillChannelEqualization;}
  Bool_t IsEnabledFillRecentering()           const {return fFillRecenterQvec;}
  Bool_t IsEnabledFillRotation()              const {return fFillRotateQvec;}
  Bool_t IsEnabledFillTwist()                 const {return fFillTwistQvec;}
  Bool_t IsEnabledFillScaling()               const {return fFillScaleQvec;}
  Bool_t IsTracking()              const {return fIsTracking;}

  Bool_t CorrectWithEventLabel()   const {return fUseLabel;}

  Bool_t IsRequestedCorrection(Int_t step)      const {return fRequestedCorrectionMap[step];}
  Bool_t IsRequestedFillHistogram(Int_t step)   const {return fRequestedHistogramMap[step];}
  Bool_t IsApplyCorrection(Int_t step)          const {return fApplyCorrectionMap[step];}
  Bool_t IsFillHistogram(Int_t step)            const {return fFillHistogramMap[step];}
  //Int_t GetCorrectionPass(Int_t step) const {return fCorrectionPasses[step];}


  
 private:

  // setters for framework
  void SetCalibrationStep(Int_t index)                      { fCalibrationStep = index;}
  
  void SetMaps(AliQnCorrectionsConstants::CorrectionSteps step, Bool_t b=kTRUE) {fRequestedCorrectionMap[(Int_t) step]=b;fRequestedHistogramMap[(Int_t) step]=b;fApplyCorrectionMap[(Int_t) step]=b;fFillHistogramMap[(Int_t) step]=b;};

  Int_t  fQnConfigurationCorrelationIndices[3];

  AliQnCorrectionsCuts* fCuts;

  UShort_t fQnNormalization;        
  Short_t fEqualizationMethod;        
  Short_t fTwistAndRescalingMethod;        
  Int_t fLocalIndex;
  Int_t fGlobalIndex;
  Int_t fCalibrationStep;
  Int_t fMinimumHarmonic;
  Int_t fMaximumHarmonic;
  Int_t fAlignmentHarmonic;
  UShort_t fDetectorType;  
  TArrayS * fChannelList;
  TArrayS * fChannelGroups;
  TString  fCalibrationDetectorNames;
  TString  fEqualizationDetectorNames;
  TString  fQnConfigurationName;
  TString  fQnConfigurationCorrelationNames[3];
  TString fAlignmentReferenceDetector;

  AliQnCorrectionsAxes* fCommonCorrectionAxes;
  AliQnCorrectionsAxes* fDataVectorEqualizationAxes;
  AliQnCorrectionsAxes* fRecenteringAxes;
  AliQnCorrectionsAxes* fAlignmentAxes;
  AliQnCorrectionsAxes* fTwistAndRescalingAxes;

  TString fEqualizationHistPath;
  TString fRecenteringHistPath;
  TString fCorrelationHistPath;
  Bool_t fChannelEqualization;
  Bool_t fRecenterQvec;
  Bool_t fRotateQvec;
  Bool_t fTwistQvec;
  Bool_t fScaleQvec;
  Bool_t fFillChannelEqualization;
  Bool_t fFillRecenterQvec;
  Bool_t fFillRotateQvec;
  Bool_t fFillTwistQvec;
  Bool_t fFillScaleQvec;
  Bool_t fIsTracking;
  Bool_t fUseLabel;

  Bool_t fRequestedCorrectionMap[AliQnCorrectionsConstants::nCorrectionSteps];  // which correction procedures are applied
  Bool_t fApplyCorrectionMap[AliQnCorrectionsConstants::nCorrectionSteps];  // which correction procedures are applied
  Bool_t fRequestedHistogramMap[AliQnCorrectionsConstants::nCorrectionSteps];  // which correction procedures are applied
  Bool_t fFillHistogramMap[AliQnCorrectionsConstants::nCorrectionSteps];  // which correction procedures are applied
  //Int_t fCorrectionPasses[AliQnCorrectionsConstants::nCorrectionSteps];  // which correction has to be applied at which pass over data

  ClassDef(AliQnCorrectionsConfiguration, 1);
};


#endif
