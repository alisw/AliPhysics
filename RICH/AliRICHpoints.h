#ifndef AliRICHpoints_H
#define AliRICHpoints_H

#include "TPolyMarker3D.h"
#include "AliRICH.h"
#include "AliPoints.h"

class AliRICHpoints : public AliPoints {
protected:
   Int_t            fHitIndex;         // Link to hit number 
   Int_t            fTrackIndex;       // Link to track number 
   Int_t            fDigitIndex;       // Link to digit 

public:
  AliRICHpoints();
  AliRICHpoints(Int_t npoints);
  virtual ~AliRICHpoints();

  Int_t                 GetHitIndex() {return fHitIndex;}
  Int_t                 GetTrackIndex(); // *MENU*
  Int_t                 GetDigitIndex() {return fDigitIndex;}
  AliRICHhit           *GetHit() const;
  AliRICHdigit         *GetDigit() const;
  virtual void          InspectHit(); // *MENU*
  virtual void          DumpHit(); // *MENU*
  virtual void          InspectDigit(); // *MENU*
  virtual void          DumpDigit(); // *MENU*
  virtual void          GetCenterOfGravity(); // *MENU*
  virtual void          SetHitIndex(Int_t hitindex) {fHitIndex = hitindex;}
  virtual void          SetTrackIndex(Int_t trackindex) {fTrackIndex = trackindex;}
  virtual void          SetDigitIndex(Int_t digitindex) {fDigitIndex = digitindex;}
  
  ClassDef(AliRICHpoints,1) //Class to draw detector clusters (is PolyMarker3D)
};
#endif
