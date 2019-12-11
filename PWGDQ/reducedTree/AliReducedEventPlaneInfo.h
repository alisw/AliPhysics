// Class for reduced event plane information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDEVENTPLANEINFO_H
#define ALIREDUCEDEVENTPLANEINFO_H

#include <TObject.h>
#include <TMath.h>

//_________________________________________________________________________
class AliReducedEventPlaneInfo : public TObject {
  
  friend class AliAnalysisTaskReducedTreeMaker;    // friend analysis task which fills the object
  
 public: 

  enum EventPlaneStatus {
    kRaw=0,                  // Not-calibrated 
    kCalibrated,             // Calibrated (if its the case)
    kRecentered,             // Recentered
    kShifted,                // Shifted  
    kUnset,                  // Q-vector not computed
    kNMaxFlowFlags
  };

  enum EventPlaneDetector {
    kTPC=0,                  // All TPC tracks
    kTPCptWeights,           // All TPC tracks, using pt weights  
    kTPCpos,                 // Positive TPC tracks 
    kTPCneg,                 // Negative TPC tracks 
    kTPCsideA,               // TPC tracks from A-side 
    kTPCsideC,               // TPC tracks from C-side
    kVZEROA,                 // VZERO A-side channels
    kVZEROC,                 // VZERO C-side channels
    kFMD,                    // FMD 
    kZDCA,                   // ZDC A-side channels
    kZDCC,                   // ZDC C-side channels
    kNdetectors
  };
  
  static const Int_t fgkNMaxHarmonics = 10;
  
  AliReducedEventPlaneInfo();
  virtual ~AliReducedEventPlaneInfo();
  
  Double_t Qx(Int_t det, Int_t harmonic)  const {return (det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics ? fQvector[det][harmonic-1][0] : -999.);}
  Double_t Qy(Int_t det, Int_t harmonic)  const {return (det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics ? fQvector[det][harmonic-1][1] : -999.);}
  Double_t EventPlane(Int_t det, Int_t h) const;
  UChar_t GetEventPlaneStatus(Int_t det, Int_t h) const {return (det>=0 && det<kNdetectors && h>0 && h<=fgkNMaxHarmonics ? fEventPlaneStatus[det][h-1] : kUnset);} 
  Bool_t  CheckEventPlaneStatus(Int_t det, Int_t h, EventPlaneStatus flag) const;
  void    CopyEvent(const AliReducedEventPlaneInfo* event);

  void SetQx(Int_t det, Int_t harmonic, Float_t qx) { if(det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics) fQvector[det][harmonic-1][0]=qx;}
  void SetQy(Int_t det, Int_t harmonic, Float_t qy) { if(det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics) fQvector[det][harmonic-1][1]=qy;}
  void SetEventPlaneStatus(Int_t det, Int_t harmonic, EventPlaneStatus status) { 
    if(det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics) 
      fEventPlaneStatus[det][harmonic-1] |= (1<<status);
  }
  
 protected:
  // Q-vectors for the first 10 harmonics from TPC, VZERO, FMD and ZDC detectors
  Double_t fQvector[kNdetectors][fgkNMaxHarmonics][2];     // Q vector components for all detectors and 6 harmonics
  UChar_t fEventPlaneStatus[kNdetectors][fgkNMaxHarmonics];  // Bit maps for the event plane status (1 char per detector and per harmonic)
   
  void ClearEvent();
  AliReducedEventPlaneInfo(const AliReducedEventPlaneInfo &c);
  AliReducedEventPlaneInfo& operator= (const AliReducedEventPlaneInfo &c);

  ClassDef(AliReducedEventPlaneInfo, 1);
};

//_______________________________________________________________________________
inline Double_t AliReducedEventPlaneInfo::EventPlane(Int_t det, Int_t harmonic) const
{
  //
  // Event plane from detector "det" and harmonic "harmonic"
  //
  if(det<0 || det>=kNdetectors || harmonic<1 || harmonic>fgkNMaxHarmonics) return -999.;
  return TMath::ATan2(fQvector[det][harmonic-1][1], fQvector[det][harmonic-1][0])/Double_t(harmonic);
}

//_______________________________________________________________________________
inline Bool_t AliReducedEventPlaneInfo::CheckEventPlaneStatus(Int_t det, Int_t h, EventPlaneStatus flag) const {
  //
  // Check the status of the event plane for a given detector and harmonic
  //
  if(det<0 || det>=kNdetectors || h<1 || h>fgkNMaxHarmonics) return kFALSE;
  return (flag<kNMaxFlowFlags ? (fEventPlaneStatus[det][h-1]&(1<<flag)) : kFALSE);
}

#endif
