#ifndef ALITRDMCTRACK_H
#define ALITRDMCTRACK_H  

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

#include <TObject.h> 

class AliTRDgeometry;

class AliTRDmcTrack : public TObject {

// Represents TRD related info about generated track

public:

  AliTRDmcTrack();
  AliTRDmcTrack(Int_t label, Bool_t primary, Float_t mass, Int_t charge, Int_t pdg); 

  void SetPin(Int_t plane, Double_t px, Double_t py, Double_t pz)
              { Pin[plane][0]  = px; Pin[plane][1]  = py; Pin[plane][2]  = pz; }

  void SetPout(Int_t plane, Double_t px, Double_t py, Double_t pz)
              { Pout[plane][0] = px; Pout[plane][1] = py; Pout[plane][2] = pz; }

  void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz, Int_t opt = 0) const;
  void GetPlanePxPyPz(Double_t &px, Double_t &py, Double_t &pz
		     ,Int_t plane, Int_t opt = 0) const;

  void Update(Int_t index) { if (fN < 199) fIndex[fN++] = index; }

  Int_t   GetTrackIndex()          const { return fLab;      }
  Bool_t  IsPrimary()              const { return fPrimary;  }
  Float_t GetMass()                const { return fMass;     }
  Int_t   GetCharge()              const { return fCharge;   }
  Int_t   GetPdgCode()             const { return fPDG;      }

  Int_t   GetNumberOfClusters()    const { return fN;        }
  Int_t   GetClusterIndex(Int_t i) const { return fIndex[i]; }  

protected:

   Int_t    fLab;             // Track index  
   Bool_t   fPrimary;         // TRUE if it's a primary particle
   Float_t  fMass;            // Mass of the MC track
   Int_t    fCharge;          // Charge of the MC track
   Int_t    fPDG;             // PDG code of the MC track

   Int_t    fN;               // Number of TRD clusters associated with the track
   Int_t    fIndex[200];      // Indices of these clusters  
			   
   Double_t Pin[6][3];        // Px,Py,Pz at the entrance of each TRD plane   
   Double_t Pout[6][3];       // Px,Py,Pz at the exit of each TRD plane

   ClassDef(AliTRDmcTrack,1)  // TRD MC track

};                   

#endif   
