#ifndef AliMFTCAHit_H
#define AliMFTCAHit_H

#include "TObject.h"

//_________________________________________________________________________________
class AliMFTCAHit : public TObject {
  
public:
  
  AliMFTCAHit();
  ~AliMFTCAHit() {};
  
  AliMFTCAHit (const AliMFTCAHit &hit);
  AliMFTCAHit &operator=(const AliMFTCAHit&);
  
  void SetPos(Double_t x, Double_t y, Double_t z) {
    fPos[0] = x; fPos[1] = y; fPos[2] = z; }
  const Double_t *GetPos() { return fPos; }
  void SetTrackGID(Int_t gid, Int_t la, Int_t id, Int_t detid) {
    fTrackGID = gid;
    fLayer = la;
    fID = id;
    fDetElemID = detid;
    fIsFace = fLayer%2;
  }
  const Int_t GetTrackGID()  { return fTrackGID;  }
  const Int_t GetLayer()     { return fLayer;     }
  const Int_t GetID()        { return fID;        }
  const Int_t GetDetElemID() { return fDetElemID; }
  const Int_t IsFace()       { return fIsFace;    }
  
  void SetUsed() { fIsUsed = kTRUE; }
  const Bool_t IsUsed() { return fIsUsed; }
  
  void SetInRoad(Int_t i) { fInRoad[fNRoads++] = i; }
  const Int_t GetNRoads() { return fNRoads; }
  const Int_t GetInRoad(Int_t i) { return fInRoad[i]; }

  void SetMFTClsId(Int_t id) { fMFTClsId = id; }
  const Int_t GetMFTClsId() { return fMFTClsId; }
  
private:
  
  Int_t    fTrackGID;     // From MC track with global identifier
  Int_t    fLayer;        // Layer number
  Int_t    fID;           // Index of the hit in the layer
  Double_t fPos[3];       // X,Y,Z position
  Bool_t   fIsUsed;       // for TrackFinder
  Int_t    fDetElemID;    // ladder ID
  Int_t    fNRoads;       // the number of found roads which contain thhis hit
  Int_t    fInRoad[100];  // index of the roads
  Int_t    fIsFace;       // "0" if on the disk side towards the IP, or "1"
  Int_t    fNInL;         // number of hit in layer
  Int_t    fMFTClsId;     // ID of MFT cluster, to combine with IsFace()
  
  ClassDef(AliMFTCAHit,2);
  
};
#endif
