#ifndef AliFTrack_H
#define AliFTrack_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast track class                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliFTrack : public TObject {

public:
                  AliFTrack() {;}
                  AliFTrack(Int_t code, Double_t charge, Double_t pT, Double_t eta, Double_t phi,
                            Double_t v11, Double_t v22, Double_t v33,
                            Double_t v12, Double_t v13, Double_t v23, Int_t iFlag);
   virtual       ~AliFTrack() {;}
   virtual void   Draw(Option_t *option="");
   virtual void   Paint(Option_t *option="");

   //getters
   Int_t          IdTrack() const {return fIdTrack;}
   Double_t       ChTrack() const {return fChTrack;}
   Double_t       PT()  const {return fPT;}
   Double_t       Eta() const {return fEta;}
   Double_t       Phi() const {return fPhi;}
   Double_t       V11() const {return fV11;}
   Double_t       V22() const {return fV22;}
   Double_t       V33() const {return fV33;}
   Double_t       V12() const {return fV12;}
   Double_t       V13() const {return fV13;}
   Double_t       V23() const {return fV23;}
   Int_t          IFlag() const {return fIFlag;}

private:
   Int_t          fIdTrack;       //Track code
   Double_t       fChTrack;       //Track charge
   Double_t       fPT;            //Track transverse momenta
   Double_t       fEta;           //Track eta
   Double_t       fPhi;           //Track phi
   Double_t       fV11;           //Element of the covariance matrix
   Double_t       fV22;           //Element of the covariance matrix
   Double_t       fV33;           //Element of the covariance matrix
   Double_t       fV12;           //Element of the covariance matrix
   Double_t       fV13;           //Element of the covariance matrix
   Double_t       fV23;           //Element of the covariance matrix
   Int_t          fIFlag;         //Status flag


   ClassDef(AliFTrack, 1)   //AliFast track class
};

#endif
