// $Id$

#ifndef AliESDJet_H
#define AliESDJet_H

class AliESDJet: public AliVTrack { 
 public:
  AliESDJet();
  ~AliESDJet() {;}
 protected:
  Double32_t       fPt;       //[0,0,12]   pt at vertex
  Double32_t       fEta;      //[-1,1,12]  eta at vertex
  Double32_t       fPhi;      //[0,6.3,12] phi at vertex

  ClassDef(AliESDJet, 1) // ESD jet class
};
#endif
