#ifndef ALITRDMCTRACK_H
#define ALITRDMCTRACK_H  

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD MC track                                                          //
//  Used for efficiency estimates and matching of reconstructed tracks    //
//  to MC particles                                                       //                    
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h> 

class AliTRDgeometry;

const Int_t kMAX_TB = 30;  

class AliTRDmcTrack : public TObject {

 public:

  AliTRDmcTrack();
  AliTRDmcTrack(Int_t label, Int_t seedLabel, Bool_t primary
              , Float_t mass, Int_t charge, Int_t pdg); 

          void     SetSeedLabel(Int_t l)                    { fSeedLab           = l;    }
          void     SetNumberOfClusters(Int_t n)             { fN                 = n;    }
          void     SetPin(Int_t plane, Double_t px, Double_t py, Double_t pz)
                                                            { Pin[plane][0]      = px; 
                                                              Pin[plane][1]      = py; 
                                                              Pin[plane][2]      = pz;   }
          void     SetPout(Int_t plane, Double_t px, Double_t py, Double_t pz)
                                                            { Pout[plane][0]     = px; 
                                                              Pout[plane][1]     = py; 
                                                              Pout[plane][2]     = pz;   }
          void     SetXYZin(Int_t plane, Double_t x, Double_t y, Double_t z)
                                                            { XYZin[plane][0]    = x; 
                                                              XYZin[plane][1]    = y; 
                                                              XYZin[plane][2]    = z;    }
          void     SetXYZout(Int_t plane, Double_t x, Double_t y, Double_t z)
                                                            { XYZout[plane][0]   = x; 
                                                              XYZout[plane][1]   = y; 
                                                              XYZout[plane][2]   = z;    }

          Int_t    GetTrackIndex() const                    { return fLab;               }
          Int_t    GetSeedLabel() const                     { return fSeedLab;           }
          Float_t  GetMass() const                          { return fMass;              }
          Int_t    GetCharge() const                        { return fCharge;            }
          Int_t    GetPdgCode() const                       { return fPDG;               }
          Int_t    GetNumberOfClusters() const              { return fN;                 }
          Int_t    GetClusterIndex(Int_t ltb, Int_t p, Int_t n) const 
                                                            { return fIndex[ltb][p][n];  }
          void     GetPxPyPzXYZ(Double_t &px, Double_t &py, Double_t &pz
                              , Double_t &x,  Double_t &y,  Double_t &z
		              , Int_t opt = 0) const;
          void     GetPlanePxPyPz(Double_t &px, Double_t &py, Double_t &pz
		                , Int_t plane, Int_t opt = 0) const;
          void     GetXYZin(Int_t plane, Double_t &x, Double_t &y, Double_t &z) const 
                                                            { x = XYZin[plane][0]; 
                                                              y = XYZin[plane][1]; 
                                                              z = XYZin[plane][2]; 
                                                              return;                    }
          void     GetXYZout(Int_t plane, Double_t &x, Double_t &y, Double_t &z) const
                                                            { x = XYZout[plane][0]; 
                                                              y = XYZout[plane][1]; 
                                                              z = XYZout[plane][2]; 
                                                              return;                    }

          Bool_t   IsPrimary() const                        { return fPrimary;           }
          void     Update(Int_t ltb, Int_t p, Int_t n, Int_t index) 
                                                            { fIndex[ltb][p][n] = index; }

 protected:

          Int_t    fLab;                  //  Track index  
          Int_t    fSeedLab;              //  Seed track index  
          Bool_t   fPrimary;              //  TRUE if it's a primary particle
          Float_t  fMass;                 //  Mass of the MC track
          Int_t    fCharge;               //  Charge of the MC track
          Int_t    fPDG;                  //  PDG code of the MC track

          Int_t    fN;                    //  Number of TRD clusters associated with the track
          Int_t    fIndex[kMAX_TB][6][2]; //  Indices of these clusters  
			   
          Double_t Pin[6][3];             //  Px,Py,Pz at the entrance of each TRD plane   
          Double_t Pout[6][3];            //  Px,Py,Pz at the exit of each TRD plane

          Double_t XYZin[6][3];           //  X,Y,Z at the entrance of the TRD  
          Double_t XYZout[6][3];          //  X,Y,Z at the exit of the TRD    

  ClassDef(AliTRDmcTrack,1)               //  TRD MC track

};                   

#endif   
