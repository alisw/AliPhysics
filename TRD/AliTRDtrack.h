#ifndef ALITRDTRACK_H
#define ALITRDTRACK_H  

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  TRD track object                                                       //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <TObject.h> 

class AliTRDcluster;

class AliTRDtrack : public TObject {

// Represents reconstructed TRD track

public:

   AliTRDtrack() { fN=0;}
   AliTRDtrack(const AliTRDcluster *c, UInt_t index, const Double_t xx[5],
               const Double_t cc[15], Double_t xr, Double_t alpha);  
   AliTRDtrack(const AliTRDtrack& t);    

   Int_t    Compare(const TObject *o) const;
   void     CookdEdx(Double_t low=0.05, Double_t up=0.70);   

   Double_t GetAlpha() const {return fAlpha;}
   Double_t GetC()     const {return fC;}
   Int_t    GetClusterIndex(Int_t i) const {return fIndex[i];}    
   Float_t  GetClusterdQdl(Int_t i) const {return fdQdl[i];}    
   Double_t GetChi2()  const {return fChi2;}
   Double_t GetZchi2()  const {return fZchi2;}
   void     GetCovariance(Double_t cov[15]) const;  
   Double_t GetdEdx()  const {return fdEdx;}
   Double_t GetEta()   const {return fE;}
   Int_t    GetLabel() const {return fLab;}
   Int_t    GetNclusters() const {return fN;}
   Double_t GetP()     const {  
     return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
   }
   Double_t GetPredictedChi2(const AliTRDcluster *c) const ;
   Double_t GetPt()    const {return 0.299792458*0.4/GetC()/100;}
   void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const ;
   Double_t GetSigmaC2()   const {return fCcc;}
   Double_t GetSigmaTgl2() const {return fCtt;}
   Double_t GetSigmaY2()   const {return fCyy;}
   Double_t GetSigmaZ2()   const {return fCzz;}
   Double_t GetTgl() const {return fT;}
   Double_t GetX()   const {return fX;}
   Double_t GetY()   const {return fY;} // returns running Y
   Double_t GetZ()   const {return fZ;}
 
   Float_t  GetLikelihoodElectron() const { return fLhElectron; };
 
   Bool_t   IsSortable() const {return kTRUE;}

   Int_t    PropagateTo(Double_t xr,
                   Double_t x0=8.72,Double_t rho=5.86e-3,Double_t pm=0.139);
   Int_t    Rotate(Double_t angle);

   void     SetLabel(Int_t lab=0) {fLab=lab;}
   void     SetdEdx(Float_t dedx) {fdEdx=dedx;}  

   void     Update(const AliTRDcluster* c, Double_t chi2, UInt_t i);

   void     SetLikelihoodElectron(Float_t l) { fLhElectron = l; };   

protected:

   Int_t    fLab;         // track label  
   Double_t fChi2;        // total chi2 value for R*phi measurements
   Double_t fZchi2;       // total chi2 value for Z measurements  
   Float_t  fdEdx;        // dE/dx 

   Double_t fAlpha;       // rotation angle
   Double_t fX;           // running local X-coordinate of the track (time bin)

   Double_t fY;           // Y-coordinate of the track
   Double_t fZ;           // Z-coordinate of the track
   Double_t fC;           // track curvature
   Double_t fE;           // C*x0
   Double_t fT;           // tangent of the track dip angle   

   Double_t fCyy;                         // covariance
   Double_t fCzy, fCzz;                   // matrix
   Double_t fCcy, fCcz, fCcc;             // of the
   Double_t fCey, fCez, fCec, fCee;       // track
   Double_t fCty, fCtz, fCtc, fCte, fCtt; // parameters   

   Short_t fN;             // number of clusters associated with the track
   UInt_t  fIndex[200];    // global indexes of these clusters  
   Float_t fdQdl[200];     // cluster amplitudes corrected for track angles    
			   
   Float_t fLhElectron;    // Likelihood to be an electron    

   ClassDef(AliTRDtrack,2) // TRD reconstructed tracks

};                     


#endif   
