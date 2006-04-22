#ifndef ALITRDTRACK_H
#define ALITRDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliKalmanTrack.h>
#include <AliTRDtracklet.h>

class AliTRDcluster;
class AliTPCtrack;
class AliESDtrack;
class AliTrackReference;

const unsigned kMAXCLUSTERSPERTRACK = 210; 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class AliTRDtrack : public AliKalmanTrack {
   
  enum { kNdet = 540, kNstacks = 90, kNplane = 6, kNcham = 5, kNsect = 18 };

  friend class AliTRDtracker;

 public:

  AliTRDtrack();
   AliTRDtrack(const AliTRDcluster *c, UInt_t index, const Double_t xx[5],
               const Double_t cc[15], Double_t xr, Double_t alpha);  
   AliTRDtrack(const AliTRDtrack& t);    
   AliTRDtrack(const AliKalmanTrack& t, Double_t alpha); 
   AliTRDtrack(const AliESDtrack& t);
   ~AliTRDtrack();

   AliTRDtrack     &operator=(const AliTRDtrack &t);

   //static AliTRDtrack* MakeTrack(const AliTrackReference *ref, Double_t mass);
   Int_t           Compare(const TObject *o) const;
   void            CookdEdx(Double_t low=0.05, Double_t up=0.7);   
   Float_t         StatusForTOF();
   Double_t        GetAlpha() const                     { return fAlpha;         }
   Int_t           GetSector() const;
   Double_t        GetC() const                         { return fC;             }
   Int_t           GetClusterIndex(Int_t i) const       { return fIndex[i];      }    
   Float_t         GetClusterdQdl(Int_t i) const        { return fdQdl[i];       }    
   void            GetCovariance(Double_t cov[15]) const;  
   Double_t        GetdEdx() const                      { return fdEdx;          }
   Double_t        GetPIDsignal() const                 { return GetdEdx();      }
   Float_t         GetPIDsignals(Int_t i) const         { return fdEdxPlane[i];  }
   Int_t           GetPIDTimBin(Int_t i) const          { return fTimBinPlane[i];}
   Double_t        GetEta() const                       { return fE;             }

   void            GetExternalCovariance(Double_t cov[15]) const;   
   void            GetExternalParameters(Double_t& xr, Double_t x[5]) const;

   Double_t        GetLikelihoodElectron() const        { return fLhElectron;    }
   Double_t        Get1Pt() const;
   Double_t        GetP() const;
   Double_t        GetPredictedChi2(const AliTRDcluster* c, Double_t h01) const;
   Double_t        GetPt() const                        { return 1.0 / Get1Pt(); }   
   void            GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
   void            GetGlobalXYZ(Double_t &x, Double_t &y, Double_t &z) const;
   Int_t           GetSeedLabel() const                 { return fSeedLab;       }
   Double_t        GetSigmaC2() const                   { return fCcc;           }
   Double_t        GetSigmaTgl2() const                 { return fCtt;           }
   Double_t        GetSigmaY2() const                   { return fCyy;           }
   Double_t        GetSigmaZ2() const                   { return fCzz;           }
   Double_t        GetSnp() const                       { return fX*fC - fE;     }
   Double_t        GetTgl() const                       { return fT;             }
   Double_t        GetX() const                         { return fX;             }
   Double_t        GetY() const                         { return fY;             }
   Double_t        GetZ() const                         { return fZ;             }
   Int_t*          GetBackupIndexes()                   { return fIndexBackup;   }
   Int_t*          GetIndexes()                         { return fIndex;         }

   Double_t        GetYat(Double_t xk) const; 
   Int_t           GetProlongation(Double_t xk, Double_t &y, Double_t &z);

   void            SetStop(Bool_t stop)                 { fStopped        = stop;   }
   Bool_t          GetStop() const                      { return fStopped;          }

   Int_t           PropagateTo(Double_t xr, Double_t x0=8.72, Double_t rho=5.86e-3);
   Int_t           PropagateToX(Double_t xr, Double_t step);
   Int_t           PropagateToR(Double_t xr, Double_t step);
   void            ResetCovariance();   
   void            ResetCovariance(Float_t mult);   
   void            ResetClusters()                      { SetChi2(0.); SetNumberOfClusters(0); }
   Int_t           Rotate(Double_t angle, Bool_t absolute=kFALSE);

   void            SetdEdx(Float_t dedx)                { fdEdx           = dedx;   }  
   void            SetPIDsignals(Float_t dedx, Int_t i) { fdEdxPlane[i]   = dedx;   }
   void            SetPIDTimBin(Int_t timbin, Int_t i)  { fTimBinPlane[i] = timbin; }
   void            SetLikelihoodElectron(Float_t l)     { fLhElectron     = l;      }

   void            SetSampledEdx(Float_t q, Int_t i);
   void            SetSampledEdx(Float_t q);
   void            SetSeedLabel(Int_t lab)              { fSeedLab        = lab;    }

   Int_t           Update(const AliTRDcluster* c, Double_t chi2, UInt_t i, Double_t h01);
   Int_t           UpdateMI(const AliTRDcluster* c, Double_t chi2, UInt_t i, Double_t h01, Int_t plane);
   Int_t           UpdateMI(const AliTRDtracklet & tracklet);

   void            AddNWrong()                          { fNWrong++;                }
  
   Int_t           GetNWrong() const                    { return fNWrong;           }
   Int_t           GetNRotate() const                   { return fNRotate;          }
   Int_t           GetNCross() const                    { return fNCross;           }
   void            IncCross()                           { fNCross++; if (fBackupTrack) fBackupTrack->IncCross(); }
   AliTRDtrack*    GetBackupTrack()                     { return fBackupTrack;      }
   void            MakeBackupTrack();

 protected:

   inline void     GetXYZ(Float_t r[3]) const;

   Double_t        GetPredictedChi2(const AliCluster*/*c*/) const                  { return 0.0; }
   Int_t           Update(const AliCluster*/*c*/, Double_t /*chi2*/, UInt_t /*i*/) { return 0;   }

   Int_t           fSeedLab;                               // track label taken from seeding  
   Float_t         fdEdx;                                  // dE/dx 
   Float_t         fdEdxT;                                 // dE/dx  - truncated mean
   Float_t         fDE;                                    // integrated delta energy
   Float_t         fdEdxPlane[kNplane];                    // dE/dx from all 6 planes
   Int_t           fTimBinPlane[kNplane];                  // time bin of Max cluster from all 6 planes

   Double_t        fAlpha;                                 // rotation angle
   Double_t        fX;                                     // running local X-coordinate of the track (time bin)
   Bool_t          fStopped;                               // track stop indication

   Double_t        fY;                                     // Y-coordinate of the track
   Double_t        fZ;                                     // Z-coordinate of the track
   Double_t        fE;                                     // C*x0
   Double_t        fT;                                     // tangent of the track momentum dip angle
   Double_t        fC;                                     // track curvature

   Double_t        fCyy;                                   // covariance
   Double_t        fCzy, fCzz;                             // matrix
   Double_t        fCey, fCez, fCee;                       // of the
   Double_t        fCty, fCtz, fCte, fCtt;                 // track
   Double_t        fCcy, fCcz, fCce, fCct, fCcc;           // parameters   
   
   Int_t           fIndex[kMAXCLUSTERSPERTRACK];           // global indexes of clusters  
   Int_t           fIndexBackup[kMAXCLUSTERSPERTRACK];     // backup indexes of clusters - used in iterations
   Float_t         fdQdl[kMAXCLUSTERSPERTRACK];            // cluster amplitudes corrected for track angles    
                           
   Float_t         fLhElectron;                            // Likelihood to be an electron    
   Int_t           fNWrong;                                // number of wrong clusters
   Int_t           fNRotate;                               // number of rotation
   Int_t           fNCross;                                // number of the cross materials
   Int_t           fNExpected;                             // expected number of cluster
   Int_t           fNLast;                                 // number of clusters in last 2 layers
   Int_t           fNExpectedLast;                         // number of expected clusters on last 2 layers
   Int_t           fNdedx;                                 // number of clusters for dEdx measurment
   Float_t         fChi2Last;                              // chi2 in the  last 2 layers
   AliTRDtracklet  fTracklets[6];                          // tracklets
   Float_t         fBudget[3];                             // integrated material budget
   AliTRDtrack    *fBackupTrack;                           //! backup track

   ClassDef(AliTRDtrack,4)                                 // TRD reconstructed tracks

};                     

#endif   
