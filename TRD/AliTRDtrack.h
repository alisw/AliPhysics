#ifndef ALITRDTRACK_H
#define ALITRDTRACK_H  

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

#include <AliKalmanTrack.h>
#include <TMath.h>   

class AliTRDcluster;
class AliTPCtrack; 

const unsigned kMAX_CLUSTERS_PER_TRACK=210; 

class AliTRDtrack : public AliKalmanTrack {

// Represents reconstructed TRD track

public:

   AliTRDtrack():AliKalmanTrack(){}
   AliTRDtrack(const AliTRDcluster *c, UInt_t index, const Double_t xx[5],
               const Double_t cc[15], Double_t xr, Double_t alpha);  
   AliTRDtrack(const AliTRDtrack& t);    
   AliTRDtrack(const AliKalmanTrack& t, Double_t alpha); 

   Int_t    Compare(const TObject *o) const;
   void     CookdEdx(Double_t low=0.05, Double_t up=0.70);   

   Double_t GetAlpha() const {return fAlpha;}
   Double_t GetC()     const {return fC;}
   Int_t    GetClusterIndex(Int_t i) const {return fIndex[i];}    
   Float_t  GetClusterdQdl(Int_t i) const {return fdQdl[i];}    

   void     GetCovariance(Double_t cov[15]) const;  
   Float_t  GetdEdx()  const {return fdEdx;}
   Double_t GetEta()   const {return fE;}

   void     GetExternalCovariance(Double_t cov[15]) const ;   
   void     GetExternalParameters(Double_t& xr, Double_t x[5]) const ;

   Double_t GetLikelihoodElectron() const { return fLhElectron; };

   Double_t Get1Pt()   const {return TMath::Abs(fC*GetConvConst());} 
   Double_t GetP()     const {  
     return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
   }
   Double_t GetPredictedChi2(const AliTRDcluster*) const ;
   Double_t GetPt()    const {return 1./Get1Pt();}   
   void     GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const ;
   void     GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z) const ;
   Int_t    GetSeedLabel() const { return fSeedLab; }
   Double_t GetSigmaC2()   const {return fCcc;}
   Double_t GetSigmaTgl2() const {return fCtt;}
   Double_t GetSigmaY2()   const {return fCyy;}
   Double_t GetSigmaZ2()   const {return fCzz;}
   Double_t GetSnp()  const {return fX*fC - fE;}
   Double_t GetTgl()  const {return fT;}
   Double_t GetX()    const {return fX;}
   Double_t GetY()    const {return fY;}
   Double_t GetZ()    const {return fZ;}

   Int_t    PropagateTo(Double_t xr,
                   Double_t x0=8.72,Double_t rho=5.86e-3,Double_t pm=0.139);
   void     ResetCovariance();   
   Int_t    Rotate(Double_t angle);

   void     SetdEdx(Float_t dedx) {fdEdx=dedx;}  
   void     SetLikelihoodElectron(Float_t l) { fLhElectron = l; };  

   void     SetSampledEdx(Float_t q, Int_t i) {
              Double_t s=GetSnp(), t=GetTgl();
              q*= TMath::Sqrt((1-s*s)/(1+t*t));
              fdQdl[i]=q;
            }     

   void     SetSeedLabel(Int_t lab) { fSeedLab=lab; }

   Int_t    Update(const AliTRDcluster* c, Double_t chi2, UInt_t i);

protected:

   Int_t    fSeedLab;     // track label taken from seeding  
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

   UInt_t  fIndex[kMAX_CLUSTERS_PER_TRACK];  // global indexes of clusters  
   Float_t fdQdl[kMAX_CLUSTERS_PER_TRACK];   // cluster amplitudes corrected 
                                             // for track angles    
			   
   Float_t fLhElectron;    // Likelihood to be an electron    

   ClassDef(AliTRDtrack,2) // TRD reconstructed tracks

};                     


#endif   
