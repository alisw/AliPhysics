#ifndef ALIESDTRACKEXT_H
#define ALIESDTRACKEXT_H

// $Id$

#include "AliESDtrack.h"

class AliEsdTrackExt : public AliESDtrack {
 public: 
  AliEsdTrackExt();
  AliEsdTrackExt(const AliESDtrack &t);

  void        DeleteParams();
  Double_t    GetEmcEta() const { return fEmcEta; }
  Double_t    GetEmcPhi() const { return fEmcPhi; }
  void        MakeMiniTrack(Bool_t dall=0, Bool_t dcon=1, Bool_t dtrp=1, Bool_t dmap=1, 
                            Bool_t dits=1, Bool_t dtpc=1, Bool_t dtrd=1, Bool_t dtof=1, 
                            Bool_t dhmp=1);
  void        Setup();
 protected: 
   Double32_t fEmcEta;                      //[0,0,16]  
   Double32_t fEmcPhi;                      //[0,0,16]  
   Double32_t fNCrossedRows;                //[0,0,16]  
   Double32_t fChi2TPCConstrainedVsGlobal;  //[0,0,16]  

  ClassDef(AliEsdTrackExt,1) // Extended ESD track class
};
#endif
