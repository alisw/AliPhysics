#ifndef ALIESDTRDTRACK_H
#define ALIESDTRDTRACK_H

//
//  Tracks from the TRD Global Tracking Unit (GTU, trigger)
//  

#include "TObject.h"

class AliESDTrdTrack : public TObject {

 public:

  AliESDTrdTrack();
  virtual ~AliESDTrdTrack(){};
  AliESDTrdTrack(const AliESDTrdTrack& track);
  AliESDTrdTrack& operator=(const AliESDTrdTrack& track);

  Float_t GetYproj()     { return fYproj; };
  Float_t GetZproj()     { return fZproj; };
  Float_t GetSlope()     { return fSlope; };
  Int_t   GetDetector()  { return fDetector; };
  Int_t   GetTracklets() { return fNtracklets; };
  Int_t   GetPlanes()    { return fNplanes; };
  Int_t   GetClusters()  { return fNclusters; };
  Float_t GetPt()        { return fPt; };
  Float_t GetPhi()       { return fPhi; };
  Float_t GetEta()       { return fEta; };
  Int_t   GetLabel()     { return fLabel; };
  Float_t GetPID()       { return fPID; };
  Bool_t  IsElectron()   { return fIsElectron; }

  void SetYproj(Float_t val)     { fYproj = val; };
  void SetZproj(Float_t val)     { fZproj = val; };
  void SetSlope(Float_t val)     { fSlope = val; };
  void SetDetector(Int_t det)    { fDetector = det; };
  void SetTracklets(Int_t val)   { fNtracklets = val; };
  void SetPlanes(Int_t val)      { fNplanes = val; };
  void SetClusters(Int_t val)    { fNclusters = val; };
  void SetPt(Float_t val)        { fPt = val; };
  void SetPhi(Float_t val)       { fPhi = val; };
  void SetEta(Float_t val)       { fEta = val; };
  void SetLabel(Int_t val)       { fLabel = val; };
  void SetPID(Float_t val)       { fPID = val; };
  void SetIsElectron(Bool_t val) { fIsElectron = val; };

 protected:

  Float_t fYproj;                                   // Average values calculated
  Float_t fZproj;                                   // from the tracklets 
  Float_t fSlope;                                   //
  Int_t   fDetector;                                // First detector in the module
  Int_t   fNtracklets;                              // Number of tracklets
  Int_t   fNplanes;                                 // Number of TRD planes
  Int_t   fNclusters;                               // Total number of clusters
  Float_t fPt;                                      // Transverse momentum
  Float_t fPhi;                                     // Phi angle at the vertex
  Float_t fEta;                                     // Eta at the vertex
  Int_t   fLabel;                                   // Track label
  Float_t fPID;                                     // PID electron likelihood
  Bool_t  fIsElectron;                              // Electron flag

  ClassDef(AliESDTrdTrack,1)

};

#endif
