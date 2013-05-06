#ifndef ALICALOCLUSTERINFO_H
#define ALICALOCLUSTERINFO_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliCaloClusterInfo
// class used to extract and store reco info of calo cluster
// Author: H-S. Zhu, hongsheng.zhu@cern.ch
//                   hszhu@iopp.ccnu.edu.cn
//***********************************************************

#include <TObject.h>
#include <TString.h>
#include <TLorentzVector.h>

class AliESDEvent;
class AliVCluster;
class AliPHOSGeoUtils;

class AliCaloClusterInfo : public TObject{
 public:
 
  AliCaloClusterInfo();
  AliCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, AliPHOSGeoUtils* const phosGeo, Double_t vtx[3]);
  AliCaloClusterInfo(const AliCaloClusterInfo &src);
  AliCaloClusterInfo& operator=(const AliCaloClusterInfo &src);
  virtual ~AliCaloClusterInfo();

  TLorentzVector LorentzVector() const{ return fLorentzVector; }

  Int_t    GetModule()         const { return fModule;          }
  Int_t    GetNCells()         const { return fNCells;          }
  Int_t    GetTRUNumber()      const { return fTRUNumber;       }
  Int_t    GetNTracksMatched() const { return fNTracksMatched;  }
  Short_t  GetTrackCharge()    const { return fTrackCharge;     }
  UInt_t   GetPIDBit()         const { return fPIDBit;          }
  Double_t GetDistToBad()      const { return fDistToBad;       }
  Double_t GetEmcCpvDistance() const { return fEmcCpvDistance;  }
  Double_t GetM02()            const { return fM02;             }
  Double_t GetM20()            const { return fM20;             }
  Double_t GetTOF()            const { return fTOF;             }
  Double_t GetTrackDz()        const { return fTrackDz;         }
  Double_t GetTrackDx()        const { return fTrackDx;         }
  Double_t GetTrackPt()        const { return fTrackPt;         }

  void SetPIDBit(UInt_t bit)         { fPIDBit    |= bit;       }

  Bool_t TestCPV(Double_t mf);

  static void SetLogWeight(Float_t logWeight=0.05) { fgLogWeight = logWeight; }

 private:

  void FillCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, AliPHOSGeoUtils* const phosGeo, Double_t vtx[3]);
  void Reclusterize(AliVCluster *clust, AliPHOSGeoUtils* const phosGeo);

  Bool_t TestDisp();
  Bool_t AreNeighbors(Int_t id1,Int_t id2, AliPHOSGeoUtils* const phosGeo);
  Bool_t IsInFiducialRegion(Int_t cellX, Int_t cellZ);
  Int_t  GetTRUNumber(Int_t cellX, Int_t cellZ);

  static Float_t fgLogWeight; 

  TLorentzVector fLorentzVector;

  Int_t    fModule;
  Int_t    fNCells;                      // Number of cells in cluster
  Int_t    fTRUNumber;                   // TRU Number
  Int_t    fNTracksMatched;
  Short_t  fTrackCharge;
  UInt_t   fPIDBit;                      // PID Bit
  Double_t fDistToBad ;                  // Distance to nearest bad channel
  Double_t fEmcCpvDistance;              // Distance from PHOS EMC rec.point to the closest CPV rec.point
  Double_t fM02;                         // lambda0
  Double_t fM20;                         // lambda1
  Double_t fTOF;
  Double_t fTrackDz;
  Double_t fTrackDx;
  Double_t fTrackPt;

  ClassDef(AliCaloClusterInfo,1);
};

#endif // #ifdef ALICALOCLUSTERINFO_H
