#ifndef ALITRDGTUTRACK_H
#define ALITRDGTUTRACK_H

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD module global track (GTU)                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>

class AliTRDltuTracklet;

class AliTRDgtuTrack : public TObject {

 public:

  enum { kNmaxTrk = 12, kNplan = 6 };

  AliTRDgtuTrack();
  AliTRDgtuTrack(const AliTRDgtuTrack& track);
  virtual ~AliTRDgtuTrack();
  AliTRDgtuTrack &operator=(const AliTRDgtuTrack &t);
  virtual void Copy(TObject &t) const;

  Bool_t        IsSortable() const { return kTRUE; }
  virtual Int_t Compare(const TObject *o) const;

  void Reset();
  void ResetTracklets() { if(fTracklets) fTracklets->Delete(); };
  void AddTracklet(AliTRDltuTracklet *trk);
  AliTRDltuTracklet *GetTracklet(Int_t pos) const;
  TObjArray *Tracklets() { 
    if(!fTracklets) fTracklets = new TObjArray(400); return fTracklets; 
  };
  Int_t GetNtracklets() const {
    if (fTracklets) return fTracklets->GetEntriesFast();
    return 0;
  };

  Float_t GetYproj()     const { return fYproj; };
  Float_t GetZproj()     const { return fZproj; };
  Float_t GetSlope()     const { return fSlope; };
  Int_t   GetTracklets() const { return fNtracklets; };
  Int_t   GetPlanes()    const { return fNplanes; };
  Int_t   GetClusters()  const { return fNclusters; };
  Float_t GetPt()        const { return fPt; };
  Float_t GetPhi()       const { return fPhi; };
  Float_t GetEta()       const { return fEta; };
  Int_t   GetLabel()     const { return fLabel; };
  Int_t   GetDetector()  const { return fDetector; };
  Float_t GetPID()       const { return fPID; };
  Bool_t  IsElectron()   const { return fIsElectron; };

  void Track(Float_t xpl, Float_t field);
  void CookLabel();
  void SetDetector(Int_t det) { fDetector = det; };
  void MakePID();

 protected:

  TObjArray  *fTracklets;                           //! Array of LTU tracklets

  Float_t fYproj;                                   // Average y-projection
  Float_t fZproj;                                   // Average z-projection
  Float_t fSlope;                                   // Average slope

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

  ClassDef(AliTRDgtuTrack,2)

};

#endif
