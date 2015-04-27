/*
***********************************************************
    Event plane class that contains correction setting for specific event plane
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
***********************************************************
*/

#ifndef ALIEVENTPLANECONFIGURATION_H
#define ALIEVENTPLANECONFIGURATION_H

#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>
#include <TList.h>
#include <iostream>
#include <Rtypes.h>
#include <TArrayS.h>
#include <THn.h>
#include <TVector.h>
#include <TArrayI.h>
#include <TProfile.h>
//#include "AliReducedEvent.h"
#include "AliEventPlaneBinning.h"
#include "AliEventPlaneQvector.h"
class AliEventPlaneBinning;
//class AliEventPlaneQvector;
class AliEventPlaneCuts;
class AliEventPlaneVarManager;




//_____________________________________________________________________
class AliEventPlaneConfiguration : public TObject {


 public:
  AliEventPlaneConfiguration();
  ~AliEventPlaneConfiguration();
    

  // setters
  void SetCalibrationMethod(Float_t calibrationMethod) {  fCalibrationMethod = calibrationMethod;}
  void SetEqualizationMethod(Float_t equalizationMethod) {  fEqualizationMethod = equalizationMethod; fChannelEqualization=kTRUE;}
  void SetTwistAndScalingMethod(Float_t equalizationMethod) {  fTwistAndScalingMethod = equalizationMethod; fTwistQvec=kTRUE; fScaleQvec=kTRUE;}
  //void SetMinimumSignal(UShort_t minimumSignal) {fMinimumSignal = minimumSignal;}
  void SetCalibrationStep(Int_t calibrationStep)  { fCalibrationStep = calibrationStep;}
  void SetLocalIndex(Int_t index)  { fLocalIndex = index;}
  void SetGlobalIndex(Int_t index)  { fGlobalIndex = index;}
  //void SetEvenPlaneDetectorId(UShort_t eventPlaneDetectorStep) {fEvenPlaneDetectorId = eventPlaneDetectorStep;}
  void SetDetectorType(UShort_t detectorType) { fDetectorType = detectorType;}
  void SetChannelList(TArrayS * list)   { fChannelList = list;} 
  void SetTrackCuts(AliEventPlaneCuts* trackCuts) { fTrackCuts = trackCuts;}
  void SetCalibrationDetectorName(TString name)  {  fCalibrationDetectorNames=name;}
  void SetEqualizationDetectorName(TString name)  {  fEqualizationDetectorNames=name;}
  void SetEventPlaneDetectorName(TString name)  {  fEventPlaneDetectorNames=name;}
  void SetCorrelationDetectorNames(TString detA, TString detB)   {fCorrelationDetectorNames[0]=detA;fCorrelationDetectorNames[1]=detB;}
  void SetCorrelationDetectorIndices(Int_t detA, Int_t detB)   {fCorrelationDetectorIndices[0]=detA;fCorrelationDetectorIndices[1]=detB;}
  void SetEqualizationHistogramM(THnF* eqHistQ)  {  fInputEqualizationHistogramM=eqHistQ;}
  void SetEqualizationHistogramE(THnF* eqHistE)  {  fInputEqualizationHistogramE=eqHistE;}
  void SetCalibrationHistogramQ(Int_t ih, Int_t is, Int_t ic, THnF* eqHistQ)  { if (ih>=1 && ih<=fgkEPMaxHarmonics && (ic==0||ic==1) )fInputCalibrationHistogramsQ[is][ih-fgkEPMinHarmonics][ic]=eqHistQ;}
  void SetCalibrationHistogramE(Int_t is, THnF* eqHistE)  {  fInputCalibrationHistogramsE[is]=eqHistE;}

  void SetU2nHistogram(Int_t ih, Int_t ic, THnF* eqHistQ)  { if (ih>=fgkEPMinHarmonics && ih<=fgkEPMaxHarmonics && (ic==0||ic==1) )fInputU2nHistograms[ih-fgkEPMinHarmonics][ic]=eqHistQ;}
  void SetU2nHistogramE(THnF* eqHistE)  {  fInputU2nHistogramsE=eqHistE;}

  void SetCorrelationProfile(Int_t combination, Int_t har, Int_t component, TProfile* corProf)  {  fCorrelationProfiles[combination][har-fgkEPMinHarmonics][component]=corProf;}
  void CreateCalibrationHistograms() {CreateMultiplicityHistograms();CreateQvectorHistograms();CreateCorrelationHistograms();}
  void CreateQvectorHistograms();
  void CreateCorrelationHistograms();
  void CreateMultiplicityHistograms();
  void SetU2nProfile(Int_t har, Int_t component, TProfile* corProf)  {  fU2nProfiles[har-1][component]=corProf;}
  //void SetEqualizationBinning(AliEventPlaneBinning* eqBins) {  fEqualizationBinning= new AliEventPlaneBinning(eqBins );}
  //void SetCalibrationBinning(AliEventPlaneBinning* calBins) {  fCalibrationBinning = new AliEventPlaneBinning(calBins);}
  void SetEqualizationBinning(AliEventPlaneBinning* eqBins) {  fEqualizationBinning= eqBins ;}
  void SetCalibrationBinning(AliEventPlaneBinning* calBins) {  fCalibrationBinning = calBins;}
  void SetEqualizationHistPath(const Char_t* path)  {  fEqualizationHistPath=Form("%s",path);}
  void SetRecenteringHistPath(const Char_t* path)  {  fRecenteringHistPath=Form("%s",path); fRecenterQvec=kTRUE;}
  void SetCorrelationHistPath(const Char_t* path)  {  fCorrelationHistPath=Form("%s",path);}
  void SetChannelEqualization(Bool_t set) {fChannelEqualization=set;}
  void SetRecentering(Bool_t set) {fRecenterQvec=set;}
  void SetRotation(Bool_t set) {fRotateQvec=set;}
  void SetTwist(Bool_t set) {fTwistQvec=set;}
  void SetScaling(Bool_t set) {fScaleQvec=set;}
  void ConnectInputMultiplicityHistograms(TString file);
  void ConnectInputCalibrationHistograms(TString file);

  // getters
  UShort_t CalibrationMethod()   const  {return  fCalibrationMethod ;}
  Short_t EqualizationMethod()    const {return   fEqualizationMethod ;}
  Short_t TwistAndScalingMethod()    const {return   fTwistAndScalingMethod ;}
  //Float_t MinimumSignal()        {return fMinimumSignal;}
  Int_t CalibrationStep()    const   {return  fCalibrationStep ;}
  Int_t LocalIndex()    const   {return  fLocalIndex ;}
  Int_t GlobalIndex()    const   {return  fGlobalIndex ;}
  //UShort_t EvenPlaneDetectorId(()      {return fEvenPlaneDetectorId;}
  UShort_t DetectorType()     const     {return  fDetectorType ;}
  //Double_t * ChannelList()        {return fChannelList->At(ch);}
  UShort_t UseChannel(Int_t ch)    const     {return  fChannelList->At(ch);}
  TArrayS* ChannelList()  const {return fChannelList;} 
  AliEventPlaneCuts* TrackCuts() const {return fTrackCuts;}
  Bool_t IsTrackSelected(Float_t* values);
  TString CalibrationDetectorName()  const {return  fCalibrationDetectorNames ;}
  TString EqualizationDetectorName()  const {return  fEqualizationDetectorNames ;}
  TString EventPlaneDetectorName()  const {return  fEventPlaneDetectorNames;}
  TString CorrelationDetectorName(Int_t detCor) const {return fCorrelationDetectorNames[detCor];}
  Int_t CorrelationDetectorIndex(Int_t detCor) const {return fCorrelationDetectorIndices[detCor];}
  THnF* EqualizationHistM()  const {return fEqualizationHistM;}
  THnF* EqualizationHistE()  const {return fEqualizationHistE;}
  THnF* EqualizationHistogramM(Int_t step)  const { return (step < 3 ?  fEqualizationHistogramsM[step]:  0x0);}
  THnF* EqualizationHistogramE(Int_t step)  const { return (step < 3 ?  fEqualizationHistogramsE[step]:  0x0);}
  THnF* InputEqualizationHistogramM()  const { return fInputEqualizationHistogramM;}
  THnF* InputEqualizationHistogramE()  const { return fInputEqualizationHistogramE;}
  THnF* CalibrationHistQ(Int_t ih, Int_t ic)  const { if( ih>=1 && ih<=fgkEPMaxHarmonics && (ic==0||ic==1)) return fCalibrationHistQ[ih-1][ic]; else return 0x0;}
  THnF* CalibrationHistogramQ(Int_t step, Int_t ih, Int_t ic)  const { if(step>=0 && step <6 && ih>=1 && ih<=fgkEPMaxHarmonics && (ic==0||ic==1)) return fCalibrationHistogramsQ[step][ih-fgkEPMinHarmonics][ic]; else return 0x0;}
  THnF* CalibrationHistE()  const {return fCalibrationHistE;}
  THnF* CalibrationHistogramE(Int_t step)  const {return fCalibrationHistogramsE[step];}


  THnF* InputCalibrationHistogramQ(Int_t step, Int_t ih, Int_t ic)  const { if(step>=0 && step <6 && ih>=1 && ih<=fgkEPMaxHarmonics && (ic==0||ic==1)) return fInputCalibrationHistogramsQ[step][ih-fgkEPMinHarmonics][ic]; else return 0x0;}
  THnF* InputCalibrationHistogramE(Int_t step)  const {return fInputCalibrationHistogramsE[step];}

  TProfile* CorrelationProfile(Int_t combination, Int_t har, Int_t component)  const {return fCorrelationProfiles[combination][har-1][component];}
  TProfile* CorrelationProf(Int_t stage, Int_t det, Int_t har, Int_t component)  const {return fCorrelationProfs[stage][det][har-fgkEPMinHarmonics][component];}
  TProfile* CorrelationEpProf(Int_t stage, Int_t det, Int_t har, Int_t component)  const {return fCorrelationEpProfs[stage][det][har-fgkEPMinHarmonics][component];}
  TProfile* U2nProfile(Int_t har, Int_t component)  const {return fU2nProfiles[har-1][component];}
  THnF* U2nHistogram(Int_t har, Int_t component)  const {return fU2nHistograms[har-fgkEPMinHarmonics][component];}
  THnF* U2nHistogramE()  const {return fU2nHistogramsE;}
  THnF* InputU2nHistogram(Int_t har, Int_t component)  const {return fInputU2nHistograms[har-fgkEPMinHarmonics][component];}
  THnF* InputU2nHistogramE()  const {return fInputU2nHistogramsE;}
  AliEventPlaneBinning* EqualizationBinning()  const {return fEqualizationBinning;}
  AliEventPlaneBinning* CalibrationBinning()  const {return fCalibrationBinning;}
  void Copy(AliEventPlaneConfiguration* epConf, Int_t localIndex, Int_t globalIndex);
  TString EqualizationHistPath()  const {return fEqualizationHistPath;}
  TString RecenteringHistPath()  const {return fRecenteringHistPath;}
  TString CorrelationHistPath()  const {return fCorrelationHistPath;}
  TClonesArray* Qvectors()  const {return fQvectors;}

  Bool_t doChannelEqualization() const {return fChannelEqualization;}
  Bool_t doRecentering() const {return fRecenterQvec;}
  Bool_t doRotation() const {return fRotateQvec;}
  Bool_t doTwist() const {return fTwistQvec;}
  Bool_t doScaling() const {return fScaleQvec;}
  
 private:

  Int_t  fCorrelationDetectorIndices[2];
  AliEventPlaneCuts* fTrackCuts;
  THnF* fEqualizationHistM;
  THnF* fEqualizationHistE;
  THnF* fCalibrationHistQ[fgkNHarmonics][2];
  THnF* fCalibrationHistE;
  TProfile* fCorrelationProfiles[3][fgkNHarmonics][4];
  TProfile* fU2nProfiles[fgkNHarmonics*2][2];


  UShort_t fCalibrationMethod;        
  Short_t fEqualizationMethod;        
  Short_t fTwistAndScalingMethod;        
  Int_t fCalibrationStep;        
  Int_t fLocalIndex;
  Int_t fGlobalIndex;
  UShort_t fDetectorType;  
  TArrayS * fChannelList;
  TString  fCalibrationDetectorNames;
  TString  fEqualizationDetectorNames;
  TString  fEventPlaneDetectorNames;
  TString  fCorrelationDetectorNames[2];
  THnF* fEqualizationHistogramsM[3];
  THnF* fEqualizationHistogramsE[3];
  THnF* fCalibrationHistogramsQ[6][fgkNHarmonics][2];
  THnF* fCalibrationHistogramsE[6];
  TProfile* fCorrelationProfs[6][3][fgkNHarmonics][4];
  TProfile* fCorrelationEpProfs[6][3][fgkNHarmonics][4];
  THnF* fU2nHistograms[fgkNHarmonics][2];
  THnF* fU2nHistogramsE;
  THnF* fInputEqualizationHistogramM;
  THnF* fInputEqualizationHistogramE;
  THnF* fInputCalibrationHistogramsQ[6][fgkNHarmonics][2];
  THnF* fInputCalibrationHistogramsE[6];
  THnF* fInputU2nHistograms[fgkNHarmonics][2];
  THnF* fInputU2nHistogramsE;
  AliEventPlaneBinning* fEqualizationBinning;
  AliEventPlaneBinning* fCalibrationBinning;
  //UShort_t fEventPlaneDetectorId[50];
  TString fEqualizationHistPath;
  TString fRecenteringHistPath;
  TString fCorrelationHistPath;
  Bool_t fChannelEqualization;
  Bool_t fRecenterQvec;
  Bool_t fRotateQvec;
  Bool_t fTwistQvec;
  Bool_t fScaleQvec;
  TString fStages[6];
  TClonesArray* fQvectors;





  //AliESDtrackCuts* fESDTrackCuts;
  //AliAODtrackCuts* fAODTrackCuts;

  ClassDef(AliEventPlaneConfiguration, 1);
};


#endif
