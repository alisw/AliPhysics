#ifndef AliFTrack_H
#define AliFTrack_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast track class                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class AliFTrack : public TObject {

private:
   Int_t          fIdTrack;       //Track code
   Double_t       fChTrack;       //Track charge
   Double_t       fPT;            //Track transverse momenta
   Double_t       fEta;           //Track eta
   Double_t       fPhi;           //Track phi
   Double_t       fV11;           //
   Double_t       fV22;           //
   Double_t       fV33;           //
   Double_t       fV12;           //
   Double_t       fV13;           //
   Double_t       fV23;           //
   Int_t          fIFlag;         //

public:
                  AliFTrack() {;}
                  AliFTrack(Int_t code, Double_t charge, Double_t pT, Double_t eta, Double_t phi,
                            Double_t v11, Double_t v22, Double_t v33,
                            Double_t v12, Double_t v13, Double_t v23, Int_t iFlag);
   virtual       ~AliFTrack() {;}
   virtual void   Draw(Option_t *option="");
   virtual void   Paint(Option_t *option="");

   //getters
   Int_t          IdTrack() {return fIdTrack;}
   Double_t       ChTrack() {return fChTrack;}
   Double_t       PT()  {return fPT;}
   Double_t       Eta() {return fEta;}
   Double_t       Phi() {return fPhi;}
   Double_t       V11() {return fV11;}
   Double_t       V22() {return fV22;}
   Double_t       V33() {return fV33;}
   Double_t       V12() {return fV12;}
   Double_t       V13() {return fV13;}
   Double_t       V23() {return fV23;}
   Int_t          IFlag() {return fIFlag;}


   ClassDef(AliFTrack, 1)   //AliFast track class
};

#endif
