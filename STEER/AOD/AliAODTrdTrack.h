#ifndef ALIAODTRDTRACK_H
#define ALIAODTRDTRACK_H

/// \class AliAODTrdTrack
/// \brief format for the TRD tracks calculated in the
///
/// Global Tracking Unit, used for the TRD L1 trigger
/// \author Jochen Klein <jochen.klein@cern.ch>

#include "TRef.h"
#include "TClonesArray.h"
#include "AliVTrack.h"
#include "AliVTrdTrack.h"
#include "AliAODTrdTracklet.h"

class AliAODTrdTrack : public AliVTrdTrack {
 public:

  AliAODTrdTrack();
  AliAODTrdTrack(const AliVTrdTrack &rhs);
  virtual ~AliAODTrdTrack() {};
  AliAODTrdTrack(const AliAODTrdTrack& track);
  AliAODTrdTrack& operator=(const AliAODTrdTrack& track);
  virtual void Copy(TObject &obj) const;

  virtual Int_t GetA()         const { return fA; }
  virtual Int_t GetLayerMask() const { return fLayerMask; }
  virtual Int_t GetPID()       const { return fPID; }
  virtual Int_t GetPt()        const;
  virtual Int_t GetStack()     const { return fGlobalStack%5; }
  virtual Int_t GetSector()    const { return fGlobalStack/5; }

  virtual Bool_t GetTrackInTime() const { return (fFlagsTiming & 0x1); }
  virtual UChar_t GetFlagsTiming() const { return fFlagsTiming; }

  virtual Int_t GetLabel()     const { return fLabel; }

  virtual Double_t Pt()        const { return GetPt() / 128.; }

  Int_t GetNTracklets() const {
    Int_t count = 0;
    for (Int_t iLayer = 0; iLayer < 6; ++iLayer)
      count += (fLayerMask >> iLayer) & 1;
    return count;
  }
  virtual AliAODTrdTracklet* GetTracklet(Int_t idx) const { return (AliAODTrdTracklet*) fTracklets[idx]; }

  virtual AliVTrack* GetTrackMatch() const { return (AliVTrack*) fTrackMatch.GetObject(); }

  virtual void SetA(Int_t a) { fA = a; }
  virtual void SetLayerMask(Int_t mask) { fLayerMask = mask; }
  virtual void SetPID(Int_t pid) { fPID = pid; }
  virtual void SetLabel(Int_t label) { fLabel = label; }
  virtual void SetSector(Int_t sector) { fGlobalStack = 5*sector + (fGlobalStack%5); }
  virtual void SetStack(Int_t stack) { fGlobalStack = 5*(fGlobalStack%5) + stack; }

  void AddTracklet(const AliVTrdTracklet& trkl, Int_t layer) { new (fTracklets[layer]) AliAODTrdTracklet(trkl); }
  void SetTrackMatchReference(AliVTrack *trk) { fTrackMatch = trk; }

  virtual Bool_t IsSortable() const { return kFALSE; }
  virtual Int_t Compare(const TObject* /* obj */) const { return 0; }

 protected:

  Char_t   fGlobalStack;		  ///< stack (0-89) in which the track was found
					  // (unique because of stack-wise tracking)
  UChar_t  fPID;			  ///< electron PID for this track
  UChar_t  fLayerMask;			  ///< mask of contributing tracklets
  Int_t    fA;				  ///< transverse offset from nominal primary vertex
  UChar_t  fFlagsTiming;                  ///< timing flags

  TClonesArray fTracklets;                ///< array of contributing tracklets
  TRef fTrackMatch;                       ///< reference to matched global track

  Int_t fLabel;				  ///< Track label

  ClassDef(AliAODTrdTrack,1)
};

#endif
