#ifndef AliToyMCEvent_H
#define AliToyMCEvent_H


#include <TClonesArray.h>
#include "AliToyMCTrack.h"

class AliToyMCEvent : public TObject {
 public:
  AliToyMCEvent();
  AliToyMCEvent(const AliToyMCEvent &event);
  virtual ~AliToyMCEvent() {}
  AliToyMCEvent& operator = (const AliToyMCEvent &event);

  enum EEventType {
    kPhysics=0,
    kLaser
  };
  
  static void SetInitialEventCounter(Int_t iev) { fgEvCounter=iev; }

  AliToyMCTrack* AddTrack(const AliToyMCTrack &track);
  AliToyMCTrack* AddTrack(Double_t xyz[3],Double_t pxpypz[3],
                          Double_t cv[21],Short_t sign);
  
  Int_t GetNumberOfTracks() const { return fTracks.GetEntriesFast(); }
  const AliToyMCTrack* GetTrack(Int_t track) const { return static_cast<const AliToyMCTrack*>(fTracks.At(track)); }
    
  void SetT0 (Float_t time)              { fT0   = time;          }
  void SetX(Float_t var)                 { fX = var;              }
  void SetY(Float_t var)                 { fY = var;              }
  void SetZ(Float_t var)                 { fZ = var;              }

  UInt_t GetEventNumber()    const {return fEventNumber;  }
  Float_t GetT0()            const {return fT0;           }
  Float_t GetX()             const {return fX;            }
  Float_t GetY()             const {return fY;            }
  Float_t GetZ()             const {return fZ;            }

  void SetEventType(EEventType type) { fEventType=type; }
  EEventType GetEventType() const    { return fEventType; }
  
 private:
  static Int_t fgEvCounter;
    
  UInt_t fEventNumber;

  EEventType fEventType;
    
  Float_t fT0;
  Float_t fX;
  Float_t fY;
  Float_t fZ;

  TClonesArray fTracks;

  ClassDef(AliToyMCEvent, 1);

};

#endif
