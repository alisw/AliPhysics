#ifndef ALIRSNCUTPHI_H
#define ALIRSNCUTPHI_H

/***************************************************************************
 Class for single track/daughter phi cut.
 Author: Francesca Bellini (fbellini@cern.ch)

Options:
- "OuterTPC" apply cut to phi estimated at the TPC outer radius 
- "InTRD" select phi (at outer TPC radius) region corresponding to the TRD modules
   (2010 config) where a 5° region is excluded to avoid border effects 
- "OutTRD" select phi (at outer TPC radius) region without TRD modules 
   (2010 config) where a 5° region is excluded to avoid border effects 
****************************************************************************/
#include <TMath.h>
#include <TClonesArray.h>

#include "AliESDtrack.h"
#include "AliRsnCut.h"
//#include "AliRsnValueDaughter.h";
#include "AliVTrack.h"

class AliRsnCutPhi : public AliRsnCut {
 public:
  
  AliRsnCutPhi();
  AliRsnCutPhi(const char *name, TString opt);
  AliRsnCutPhi(const AliRsnCutPhi &copy);
  AliRsnCutPhi &operator=(const AliRsnCutPhi &copy);
  virtual ~AliRsnCutPhi() { }
  
  Bool_t   IsSelected(TObject *object);
  void     SetPhiRange(Double_t a, Double_t b) {fPhiRange[0]=a; fPhiRange[1]=b; return;};
  
 protected:
  Bool_t   IsInsideTRD(AliVTrack *vtrack);
  Bool_t   IsOutsideTRD(AliVTrack *vtrack);
  Double_t   GetTrackPhi(AliVTrack * vtrack, Double_t radius);

 private:
  TString  fOption;
  Double_t fPhiRange[2];
  //AliRsnValueDaughter *fPhi;
  //  Double_t fPhi;
  ClassDef(AliRsnCutPhi, 1)
    
};

#endif
