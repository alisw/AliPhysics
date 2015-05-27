/*
***********************************************************
  Manager for event plane corrections framework
  Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2014/12/10
  *********************************************************
*/

#ifndef ALIEVENTPLANEMANAGER_H
#define ALIEVENTPLANEMANAGER_H

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
//#include <AliESDEvent.h>
//#include <AliESDtrackCuts.h>
//#include "AliReducedEvent.h"
#include "AliEventPlaneBinning.h"
#include "AliEventPlaneQvector.h"
#include "AliEventPlaneConfiguration.h"
#include "AliEventPlaneDetector.h"

//const Int_t fgkEPMaxHarmonics = 6;
//const Int_t fgkEPMaxDetectors = 20;




//_________________________________________________________________________
class AliEventPlaneManager : public TObject {
  
  

 public: 
  enum EventPlaneStatus {
    kRaw=0,
    kUndefined,
    kEqualized,
    kRecentered,
    kDiagonalized,
    kRescaled,
    kNMaxFlowFlags
  };
  
    
    enum Detector {
    kVZERO=0,
    kTPC,    // 1
    kZDC,    // 2
    kTZERO,   // 3
    kFMD,    // 4
    kNdetectors
  };

  AliEventPlaneManager();
  ~AliEventPlaneManager();

  void SetEventStatus(Bool_t b) {fUseEvent=b;}
  
  void SetInputFriendFileName(TString name) {fInputFriendFileName = name;}
  void SetOutputFriendFileName(TString name) {fOutputFriendFileName = name;}
  void SetRunLightWeight(Bool_t light)   {fRunLightWeight = light;}
  void SetCorrelationCombinations(Int_t combi, Int_t det1, Int_t det2) {fCorrelationHists[combi][0]=det1;fCorrelationHists[combi][1]=det2;}
  void SetEqualizationIndex(Int_t eq, Int_t det1) {fEqualizationHists[eq]=det1;}

  void GetQvector(Bool_t useFriendChain=kFALSE, Bool_t useEqualizedWeights=kFALSE, Float_t* values=0x0);
  void CalibrateChannels(Float_t * values, Bool_t useFriendChain=kFALSE);
  void RecenterQvec(Float_t * values, Bool_t useFriendChain=kFALSE);
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
  void U2nRescalingQvec(Float_t* values, Int_t corpar);
  void U2nTwistQvec(Float_t* values, Int_t corpar);
  void FillCorrelationHistograms(Float_t value, Int_t step);


  //void FillTPC(const AliESDEvent& esd, Float_t* values);
  //void FillTPC(AliReducedEvent* event, Float_t* values);
  //void FillVZERO(TObject* event);
  //void FillZDC(AliReducedEvent* event);
  //void FillTZERO(TObject* event);
  //void FillFMD(AliReducedEvent* event);

  void ClearEvent();
  void SetEmptyQvectors();
  void ClearUnusedDetectors()  {for(Int_t idet=0; idet<kNdetectors; idet++) if(fNEventPlaneDetectors[idet]==0) fDetector[idet]=0x0;};
  void AddEventPlaneConfiguration(AliEventPlaneConfiguration* epconf, Detector type)  {
          epconf->SetDetectorType( (UShort_t) type);
          TClonesArray& eparr = *(fEventPlaneConfigurations[(Int_t) type]);
          AliEventPlaneConfiguration *epConf=new(eparr[fNEventPlaneDetectors[(Int_t) type]]) AliEventPlaneConfiguration();
          epConf->Copy(epconf, fNEventPlaneDetectors[(Int_t) type], NEventPlaneDetectors());
          fIndex[NEventPlaneDetectors()][0] = (Int_t) type;
          fIndex[NEventPlaneDetectors()][1] = fNEventPlaneDetectors[(Int_t) type];
          fNEventPlaneDetectors[(Int_t) type]++;
          fNdetectors++;
  }
//  void AddEventPlaneMasterDetector(AliEventPlaneMasterDetector* epconf)  {
//          TClonesArray& eparr = *(fEventPlaneMasters);
//          AliEventPlaneMasterDetector *epConf=new(eparr[fNMasterDetectors]) AliEventPlaneMasterDetector();
//          epConf->Copy(epconf);
//          fNMasterDetectors++;
//  }
  void AddQvectors(Int_t q1, Int_t q2) ;
  

  AliEventPlaneConfiguration* EventPlaneConfiguration(Int_t globalIndex) ;
  AliEventPlaneConfiguration* EventPlaneConfiguration(TString name) ;
  AliEventPlaneConfiguration* EventPlaneConfiguration(Int_t det, Int_t localIndex)  {return  static_cast<AliEventPlaneConfiguration*>(fEventPlaneConfigurations[det]->At(localIndex));}
  //AliEventPlaneMasterDetector* EventPlaneMaster(Int_t det)  {return static_cast<AliEventPlaneMasterDetector*>(fEventPlaneMasters->At(det));}
  TClonesArray* EventPlaneQvector(Int_t det)  {return static_cast<TClonesArray*>(fQvectors->At(det));}
  //TClonesArray* GetEventPlaneMasterDetectors()  {return fEventPlaneMasters;}
  Int_t Ndetectors()  const                    {return fNdetectors;}
  Int_t NEventPlaneDetectors()  const          {Int_t total=0; for(Int_t idet=0; idet<kNdetectors; idet++) total+=fNEventPlaneDetectors[idet]; return total; }
  Int_t NEventPlaneDetectors(Int_t det)  const          {return fNEventPlaneDetectors[det];}
  TClonesArray* GetReducedDetector(Int_t det)  {return (det>-1&&det<kNdetectors ? fDetector[det] : 0x0);}
  TClonesArray* GetEventPlaneConfigurations(Int_t det)  {return fEventPlaneConfigurations[det];}
  TClonesArray* GetQvectors(Int_t det)         {return (TClonesArray*)fQvectors->At(det) ;}
  TList* GetQvectors()                         {return fQvectors;}

  Int_t CorrelationDetectorIndex(Int_t combi, Int_t det) {return fCorrelationHists[combi][det];}
  Int_t EqualizationIndex(Int_t eq) {return fEqualizationHists[eq];}

  TString InputFriendFileName() {return fInputFriendFileName;}
  TString OutputFriendFileName() {return fOutputFriendFileName;}
  Bool_t RunLightWeight()   {return fRunLightWeight;}   

  //static AliEventPlaneManager* Instance() { if(!fgAliEventPlaneManager) fgAliEventPlaneManager = new AliEventPlaneManager(); return fgAliEventPlaneManager;}

  //static AliEventPlaneManager* fgAliEventPlaneManager;

 private:
  //AliEventPlaneConfiguration* fEventPlaneConfiguration;
  AliEventPlaneManager(const AliEventPlaneManager &c);
  AliEventPlaneManager& operator= (const AliEventPlaneManager &c);

  Int_t fNdetectors;
  Int_t fNEventPlaneDetectors[kNdetectors];
  //Int_t fNMasterDetectors;
  Int_t fIndex[fgkEPMaxDetectors][2];

  Int_t fCorrelationHists[10000][2];
  Int_t fEqualizationHists[1000];

  Bool_t fRunLightWeight;
  Bool_t fUseEvent;       // use event for calibration histograms
  TString fInputFriendFileName;
  TString fOutputFriendFileName;

  TClonesArray* fDetector[kNdetectors];            //->   array containing global tracks
  static TClonesArray* fgDetector[kNdetectors];    //       global tracks

  TList* fQvectors;             //->   array containing Qvec
  static TList* fgQvectors;    //       global qvectors

  TClonesArray* fEventPlaneConfigurations[kNdetectors];             //->   array containing EventPlane configurations
  static TClonesArray* fgEventPlaneConfigurations[kNdetectors];    //       global EventPlane configurations

  //TClonesArray* fEventPlaneMasters;             //->   array containing EventPlane configurations
  //static TClonesArray* fgEventPlaneMasters;    //       global EventPlane configurations

  Double_t fVZEROminMult;
  Double_t fTZEROminMult;
  Double_t fZDCminMult;
  Double_t fFMDminMult;


  ClassDef(AliEventPlaneManager, 1);
};




////_______________________________________________________________________________
inline AliEventPlaneConfiguration* AliEventPlaneManager::EventPlaneConfiguration(Int_t globalIndex) {

  if(globalIndex>=(NEventPlaneDetectors())) return 0x0;
  return (AliEventPlaneConfiguration*) fEventPlaneConfigurations[fIndex[globalIndex][0]]->At(fIndex[globalIndex][1]);

 // if(globalIndex>(NEventPlaneDetectors()-1)) return 0x0;
 // Bool_t ok=kFALSE;
 // Int_t idet=0;
 // while(globalIndex > fNEventPlaneDetectors[idet]-1){
 //   globalIndex=globalIndex-fNEventPlaneDetectors[idet];
 //   idet++;
 // }
 // Int_t localIndex = globalIndex;

 // AliEventPlaneConfiguration* epConf = static_cast<AliEventPlaneConfiguration*>(fEventPlaneConfigurations[idet]->At(localIndex));
 // return epConf;

}
   

////_______________________________________________________________________________
inline AliEventPlaneConfiguration* AliEventPlaneManager::EventPlaneConfiguration(TString name) {

 AliEventPlaneConfiguration* EPconf = 0x0;
 for(Int_t idet=0; idet<kNdetectors; idet++){
   TClonesArray* epConfList=GetEventPlaneConfigurations(idet);
   TIter nextEPconf(epConfList);
   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
    if(!EPconf) continue;
    
    if(name.EqualTo(EPconf->EventPlaneDetectorName())) return EPconf;
   }
 }
 return 0x0;
}
   
#endif
