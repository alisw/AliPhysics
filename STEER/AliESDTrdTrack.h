#ifndef ALIESDTRDTRACK_H
#define ALIESDTRDTRACK_H

//  Tracks from the TRD Global Tracking Unit (GTU, trigger)
//  Stored in the ESD object
//  Author: B.Vulpescu

#include "TObject.h"

class AliESDTrdTrack : public TObject {

 public:

  AliESDTrdTrack();
  virtual ~AliESDTrdTrack(){};
  AliESDTrdTrack(const AliESDTrdTrack& track);
  AliESDTrdTrack& operator=(const AliESDTrdTrack& track);

  Double_t GetYproj()     const { return fYproj; };
  Double_t GetZproj()     const { return fZproj; };
  Double_t GetSlope()     const { return fSlope; };
  Char_t   GetDetector()  const { return fDetector; };
  Short_t   GetTracklets() const { return fNtracklets; };
  Char_t   GetPlanes()    const { return fNplanes; };
  Short_t   GetClusters()  const { return fNclusters; };
  Double_t GetPt()        const { return fPt; };
  Double_t GetPhi()       const { return fPhi; };
  Double_t GetEta()       const { return fEta; };
  Int_t   GetLabel()     const { return fLabel; };
  Double_t GetPID()       const { return fPID; };

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

 protected:

  Double32_t fYproj;                                   // Average values calculated
  Double32_t fZproj;                                   // from the tracklets 
  Double32_t fSlope;                                   // slope of the tracklet
  Double32_t fPt;                                      // Transverse momentum
  Double32_t fPhi;                                     //  Phi angle at the vertex
  Double32_t fEta;                                     // Eta at the vertex
  Double32_t fPID;                                     //[0.,1.,8] PID electron likelihood

  Int_t fLabel;                                        // Track label
  Short_t   fNtracklets;                               // Number of tracklets
  Short_t   fNclusters;                                // Total number of clusters
  Char_t   fNplanes;                                  // Number of TRD planes
  Char_t   fDetector;                                 // First detector in the module


  ClassDef(AliESDTrdTrack,2)

};

#endif
