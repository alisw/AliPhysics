#ifndef AliFTrackMaker_H
#define AliFTrackMaker_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFast TrackMaker class.                                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TH1D; 

#include "AliFDet.h"
#include "AliFMaker.h"

class AliFTrack;

class AliFTrackMaker : public AliFMaker {
public:
                  AliFTrackMaker();
                  AliFTrackMaker(const char *name, const char *title);
		  AliFTrackMaker(const AliFTrackMaker &aftmk);

   AliFTrackMaker& operator = (const AliFTrackMaker &aftmk) 
     {aftmk.Copy(*this); return (*this);}

   virtual       ~AliFTrackMaker();
   virtual void   Clear(Option_t *option="");
   virtual void   Draw(Option_t *option="");
   virtual void   Finish();
   virtual void   Init();
   virtual void   Make();
   virtual void   MakeTest(Int_t n);
   virtual void   PrintInfo() {}
   AliFTrack     *AddTrack(Int_t code, Double_t charge, Double_t pT, Double_t eta, Double_t phi,
                  Double_t v11, Double_t v22, Double_t v33, 
                  Double_t v12, Double_t v13, Double_t v23, Int_t iFlag);
           void   LogLikelyhood(Int_t idTrack, Double_t pInvers,  Double_t lambda);
           void   TPCResolution(Double_t ptransv, Double_t radiPad, Double_t lambda); 
   Double_t       ParticleMass(Int_t idTrack) const;  
   Double_t       Rapidity(Double_t pT, Double_t pZ);  
   Double_t       Angle(Double_t pX, Double_t pY);  
   Int_t          Charge(Int_t kf);
   Int_t          Compress(Int_t kf);
           void   ErrorMatrix(Int_t idTrack, Double_t pT,  Double_t eta,
                              Double_t &v11, Double_t &v22, Double_t &v33, 
                              Double_t &v12, Double_t &v13, Double_t &v23, Int_t &iFlag);

//    Getters
   Double_t       HH(Int_t id1, Int_t id2) const {return fHH[id1][id2];}
   Double_t       HHI(Int_t id1, Int_t id2) const {return fHHI[id1][id2];}
   Double_t       SigmaRPhiSQ() const {return fSigmaRPhiSQ;}
   Double_t       SigmaZSQ() const {return fSigmaZSQ;}
   Int_t          NTracks() const {return fNTracks;}


//     Getters Tracks histograms



//    Setters for tracks
   //masses
   void           SetPionMass(Double_t val=0.1395679e0) {fPionMass=val;}
   void           SetKaonMass(Double_t val=0.493646e0) {fKaonMass=val;}
   void           SetElectronMass(Double_t val=0.51099906e-3) {fElectronMass=val;}
   void           SetProtonMass(Double_t val=0.93827231e0) {fProtonMass=val;}
   void           SetHH(Int_t id1, Int_t id2, Double_t hh) {fHH[id1][id2]=hh;}
   void           SetHHI(Int_t id1, Int_t id2, Double_t hhi) {fHHI[id1][id2]=hhi;}
   void           SetSigmaRPhiSQ(Double_t val){fSigmaRPhiSQ=val;}
   void           SetSigmaZSQ(Double_t val){fSigmaZSQ=val;}
   void           SetRecTrack(Int_t val=100){fRecTrack=val;}

protected:
   void Copy(TObject &aftmk) const;

   Int_t           fNTracks;          // Number of tracks
   Int_t           fRecTrack;         // Tracks reconstruction  on/off
   //masses
   Double_t        fPionMass;         // Mass of pion
   Double_t        fKaonMass;         // Mass of kaon
   Double_t        fElectronMass;     // Mass of electron
   Double_t        fProtonMass;       // Mass of proton
   //matrices
   Double_t        fHH[kNMaxDet2][kNMaxDet2];   // Matrix 
   Double_t        fHHI[kNMaxDet2][kNMaxDet2];  // Matrix
   //TPC resolution
   Double_t        fSigmaRPhiSQ;                // Sigma R-PHI
   Double_t        fSigmaZSQ;                   // Sigma Z^2
   //     Tracks histograms (control)
   TH1D          *fResID11;          //histogram ID11 Elec: delta(1/pTot)/pTot
   TH1D          *fResID12;          //histogram ID12 Elec: delta(lambda)/lambda
   TH1D          *fResID13;          //histogram ID13 Elec: delta(phi)/phi
   TH1D          *fResID21;          //histogram ID21 Pion: delta(1/pTot)/pTot
   TH1D          *fResID22;          //histogram ID22 Pion: delta(lambda)/lambda
   TH1D          *fResID23;          //histogram ID23 Pion: delta(phi)/phi
   TH1D          *fResID31;          //histogram ID31 Kaon: delta(1/pTot)/pTot
   TH1D          *fResID32;          //histogram ID32 Kaon: delta(lambda)/lambda
   TH1D          *fResID33;          //histogram ID33 Kaon: delta(phi)/phi
   TH1D          *fResID41;          //histogram ID41 Proton: delta(1/pTot)/pTot
   TH1D          *fResID42;          //histogram ID42 Proton: delta(lambda)/lambda
   TH1D          *fResID43;          //histogram ID43 Proton: delta(phi)/phi
   //     Tracks histograms (Test job)
   TH1D          *fResID1Test;       //histogram ID1 in res.f 
   TH1D          *fResID2Test;       //histogram ID2 in res.f 
   TH1D          *fResID3Test;       //histogram ID3 in res.f 
   TH1D          *fResID4Test;       //histogram ID4 in res.f 
   TH1D          *fResID5Test;       //histogram ID5 in res.f 

   // 
   ClassDef(AliFTrackMaker, 1)   //AliFast TrackMaker
};

#endif








