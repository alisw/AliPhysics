#ifndef ALITRDTRACK_H
#define ALITRDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliESDtrack.h"
#include "AliKalmanTrack.h"

#include "AliTRDtracklet.h"

#ifndef ALITRDSEEDV1_H
#include "AliTRDseedV1.h"
#endif

class AliTrackReference;
class AliTRDcluster;
class AliTRDtrack : public AliKalmanTrack {

 public:

  enum { kMAXCLUSTERSPERTRACK = 210 };
   
  enum { kNdet      = 540
       , kNstacks   =  90
       , kNplane    =   AliESDtrack::kTRDnPlanes
       , kNcham     =   5
       , kNsect     =  18
       , kNslice    =   3
       , kNMLPslice =   8 };
  
  enum AliTRDPIDMethod {
         kNN = 0
       , kLQ = 1 };

         AliTRDtrack();
         AliTRDtrack(/*const*/ AliTRDcluster *c, Int_t index
                   , const Double_t xx[5], const Double_t cc[15]
                   , Double_t xr, Double_t alpha);
         AliTRDtrack(const AliTRDtrack &t/*, const Bool_t owner = kTRUE*/);    
         AliTRDtrack(const AliESDtrack &t);
 virtual ~AliTRDtrack();
         AliTRDtrack(const AliKalmanTrack &t, Double_t alpha); 

         void            ResetClusters()                              { SetChi2(0.0); 
                                                                        SetNumberOfClusters(0);              }
         Int_t           Compare(const TObject *o) const;
         void            AddNWrong()                                  { fNWrong++;                           }  
         void            IncCross()                                   { fNCross++; 
                                                                        if (fBackupTrack) 
                                                                          fBackupTrack->IncCross();          }

         Int_t           GetSector() const;
         Float_t         GetClusterdQdl(Int_t i) const                { return fdQdl[i];                     }    
         Double_t        GetdEdx() const                              { return fdEdx;                        }
         Int_t           GetNdedx() const                             { return fNdedx;                       }
         Double_t        GetPIDsignal() const                         { return GetdEdx();                    }
         Int_t           GetClusterIndex(Int_t i) const               { return fIndex[i];                    }
	 using AliExternalTrackParam::GetC;
         Double_t        GetC() const                                 { return AliExternalTrackParam
                                                                             ::GetC(GetBz());                }
         Double_t        GetPredictedChi2(const AliTRDcluster *c
                                        , Double_t h01) const;
         AliTRDPIDMethod GetPIDMethod() const                         { return fPIDmethod;                   }
         Float_t         GetPIDsignals(Int_t iPlane
                                     , Int_t iSlice) const            { return fdEdxPlane[iPlane][iSlice];   }
         Int_t           GetPIDTimBin(Int_t i) const                  { return fTimBinPlane[i];              }
         Double_t        GetLikelihoodElectron() const                { return fLhElectron;                  }
         Int_t           GetSeedLabel() const                         { return fSeedLab;                     }
         Int_t          *GetBackupIndexes()                           { return fIndexBackup;                 }
         Int_t          *GetIndexes()                                 { return fIndex;                       }
         Int_t           GetProlongation(Double_t xk
                                       , Double_t &y
                                       , Double_t &z);
         Bool_t          GetStop() const                              { return fStopped;                     }
         Int_t           GetNRotate() const                           { return fNRotate;                     }
         Int_t           GetNWrong() const                            { return fNWrong;                      }
         Int_t           GetNCross() const                            { return fNCross;                      }
         Int_t           GetNExpected() const                         { return fNExpected;                   }
         Int_t           GetNLast() const                             { return fNLast;                       }
         Int_t           GetNExpectedLast() const                     { return fNExpectedLast;               }
         AliTRDtracklet  GetTracklets(Int_t i) const                  { return fTracklets[i];                }
         Float_t         GetBudget(Int_t i) const                     { return fBudget[i];                   }
         Float_t         GetChi2Last() const                          { return fChi2Last;                    }
         AliTRDtrack    *GetBackupTrack()                             { return fBackupTrack;                 }
         // dummy to bridge the function in AliTRDtrackV1
	 //Int_t          GetNumberOfClusters() const                   { printf("AliTRDtrack::GetNumberOfClusters()\n"); 
         //                                                               return AliKalmanTrack::GetNumberOfClusters();   }
 inline virtual Int_t    GetNumberOfTracklets() const;
 virtual Int_t           GetTrackletIndex(Int_t plane) const          { return plane>=0 && plane<6 
                                                                             ? fTrackletIndex[plane] : -1;   } 


         void            SetdEdx(Double_t dedx)                       { fdEdx                      = dedx;   }
         void            SetStop(Bool_t stop)                         { fStopped                   = stop;   }
         void            SetPIDsignals(Float_t dedx
                                     , Int_t iPlane
                                     , Int_t iSlice)                  { fdEdxPlane[iPlane][iSlice] = dedx;   }
         void            SetPIDTimBin(Int_t timbin, Int_t i)          { fTimBinPlane[i]            = timbin; }
         void            SetLikelihoodElectron(Float_t l)             { fLhElectron                = l;      }
         void            SetSampledEdx(Float_t q, Int_t i);
         void            SetSampledEdx(Float_t q);
         void            SetSeedLabel(Int_t lab)                      { fSeedLab                   = lab;    }
         void            SetNWrong(Int_t nwrong)                      { fNWrong                    = nwrong; }
         void            SetNCross(Int_t ncross)                      { fNCross                    = ncross; }
         void            SetNExpected(Int_t nexp)                     { fNExpected                 = nexp;   }
         void            SetNLast(Int_t nlast)                        { fNLast                     = nlast;  }
         void            SetNExpectedLast(Int_t nexp)                 { fNExpectedLast             = nexp;   }
         void            SetChi2Last(Float_t chi2)                    { fChi2Last                  = chi2;   }
         void            SetTracklets(Int_t i, AliTRDtracklet t)      { fTracklets[i]              = t;      }
	 void            SetBudget(Int_t i, Float_t budget)           { fBudget[i]                 = budget; }
         void            SetPIDMethod(AliTRDPIDMethod method)         { fPIDmethod                 = method; }

         void	         SetTrackSegmentDirMom(const Int_t plane);
         void            CookdEdx(Double_t low = 0.05, Double_t up = 0.7);
         void            CookdEdxTimBin(const Int_t tid);
         Bool_t          CookPID(Int_t &pidQuality); 
         void            SetCluster(AliTRDcluster* cl
                                  , Int_t index = -1)                 { fClusters[(index == -1) 
                                                                                  ? GetNumberOfClusters()-1 
                                                                                  : index]         = cl;     }
         AliTRDcluster*  GetCluster(Int_t layer) const                { return (layer >= 0 &&
                                                                                layer <  GetNumberOfClusters()) 
                                                                             ? fClusters[layer] 
                                                                             : 0x0;                          }
         Float_t         GetMomentumPlane(Int_t plane) const          { return (plane >= 0 && 
                                                                                plane <  kNplane) 
                                                                             ? fMom[plane] 
                                                                             : 0.0;                          }
   const Double_t*       GetPID() const                               { return fPID;                         }
         Float_t         GetSnpPlane(Int_t plane) const               { return (plane >= 0 && 
                                                                                plane <  kNplane) 
                                                                             ? fSnp[plane] 
                                                                             : 0.0;                          }
         Float_t         GetTglPlane(Int_t plane) const               { return (plane >= 0 && 
                                                                                plane <  kNplane) 
                                                                             ? fTgl[plane]  
                                                                             : 0.0;                          }
         Float_t         GetTrackLengthPlane(Int_t plane) const;

         void            MakeBackupTrack();
         Bool_t          PropagateTo(Double_t xr, Double_t x0 = 8.72, Double_t rho = 5.86e-3);
         Int_t           PropagateToR(Double_t xr, Double_t step);
         Int_t           PropagateToX(Double_t xr, Double_t step);
         Bool_t          Rotate(Double_t angle, Bool_t absolute = kFALSE);
         Float_t         StatusForTOF();
         Bool_t          Update(const AliTRDcluster *c, Double_t chi2, Int_t index, Double_t h01);
         Bool_t          Update(const AliTRDtracklet &tracklet, Double_t chi2, Int_t index);
         Int_t           UpdateMI(/*const */AliTRDcluster *c, Double_t chi2
                                , Int_t index, Double_t h01
                                , Int_t plane, Int_t tid = 0);

 protected:

  AliTRDtrack           &operator=(const AliTRDtrack &t);

	 void	         CookdEdxNN(Float_t *dedx);
         Double_t        GetBz() const;
         Bool_t          Update(const AliCluster */*c*/, Double_t /*chi2*/, Int_t /*idx*/) { return 0;   }
         Double_t        GetPredictedChi2(const AliCluster* /*c*/) const                   { return 0.0; }

 protected:

         Int_t    fSeedLab;                           //  Track label taken from seeding  
         Float_t  fdEdx;                              //  dE/dx (truncated mean)
         Float_t  fDE;                                //  Integrated delta energy
         Float_t  fdEdxPlane[kNplane][kNslice];       //  dE/dx from all 6 planes in 3 slices each
         Int_t    fTimBinPlane[kNplane];              //  Time bin of Max cluster from all 6 planes
         UChar_t  fPIDquality;                        //  No of planes used for PID calculation	
         Double_t fPID[AliPID::kSPECIES];             //  PID probabilities
	 Float_t  fMom[kNplane];                      //  Track momentum at chamber entrance
	 Float_t  fSnp[kNplane];                      //  Track direction
	 Float_t  fTgl[kNplane];                      //  Track direction
  AliTRDcluster  *fClusters[kMAXCLUSTERSPERTRACK];    //  List of assigned clusters
         Bool_t   fClusterOwner;                      //  Indicates the track is owner of cluster
  AliTRDPIDMethod fPIDmethod;                         //  Switch between different PID methods
         Bool_t   fStopped;                           //  Track stop indication
         Int_t    fIndex[kMAXCLUSTERSPERTRACK];       //  Global indexes of clusters  
         Int_t    fIndexBackup[kMAXCLUSTERSPERTRACK]; //  Backup indexes of clusters - used in iterations
         Float_t  fdQdl[kMAXCLUSTERSPERTRACK];        //  Cluster amplitudes corrected for track angles    
           
	 Float_t  fLhElectron;                        //  Likelihood to be an electron
         Int_t    fNWrong;                            //  Number of wrong clusters
         Int_t    fNRotate;                           //  Number of rotation
         Int_t    fNCross;                            //  Number of the cross materials
         Int_t    fNExpected;                         //  Expected number of cluster
         Int_t    fNLast;                             //  Number of clusters in last 2 layers
         Int_t    fNExpectedLast;                     //  Number of expected clusters on last 2 layers
         Int_t    fNdedx;                             //  Number of clusters for dEdx measurment
         Float_t  fChi2Last;                          //  Chi2 in the  last 2 layers
  AliTRDtracklet  fTracklets[6];                      //  Tracklets
         Float_t  fBudget[3];                         //  Integrated material budget
  AliTRDtrack    *fBackupTrack;                       //! Backup track
	
	 Int_t    fTrackletIndex[6];                  //  Tracklets index in the tracker list
  AliTRDseedV1    fTracklet[6];                       //  Tracklets array defining the track

  ClassDef(AliTRDtrack,9)                             //  TRD reconstructed tracks

};                     

//___________________________________________________________
inline Int_t AliTRDtrack::GetNumberOfTracklets() const
{
	Int_t ntrklt = 0;
	for(int ip=0; ip<6; ip++) if(fTrackletIndex[ip] >= 0) ntrklt++;
	return ntrklt;
}


#endif   

