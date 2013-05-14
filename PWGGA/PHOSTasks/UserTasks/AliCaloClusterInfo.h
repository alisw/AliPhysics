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

class AliCaloClusterInfo : public TObject{
 public:
 
  AliCaloClusterInfo();
  AliCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, Int_t relID[4], Double_t mf);
  AliCaloClusterInfo(const AliCaloClusterInfo &src);
  AliCaloClusterInfo& operator=(const AliCaloClusterInfo &src);
  virtual ~AliCaloClusterInfo();

  TLorentzVector LorentzVector() const { return fLorentzVector; }

  Int_t    GetModule()         const { return fModule;          }
  Int_t    GetTRUNumber()      const { return fTRUNumber;       }
  Int_t    GetNCells()         const { return fNCells;          }
  UInt_t   GetPIDBit()         const { return fPIDBit;          }
  Double_t GetDistToBad()      const { return fDistToBad;       }
  Double_t GetEmcCpvDistance() const { return fEmcCpvDistance;  }
  Double_t GetM02()            const { return fM02;             }
  Double_t GetM20()            const { return fM20;             }
  Double_t GetTOF()            const { return fTOF;             }

  void SetLorentzVector(TLorentzVector momentum) { fLorentzVector = momentum; }
  void SetPIDBit(UInt_t bit)                     { fPIDBit       |= bit;      }

  Bool_t   IsInFiducialRegion(Int_t cellX, Int_t cellZ);
  Double_t TestDisp();

 private:

  void FillCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, Int_t relID[4], Double_t mf);

  Int_t    GetTRUNumber(Int_t cellX, Int_t cellZ);
  Double_t TestCpv(Double_t trkPt, Short_t trkCharge, Double_t trkDz, Double_t trkDx, Double_t mf);

  TLorentzVector fLorentzVector;

  Int_t    fModule;
  Int_t    fTRUNumber;          // TRU Number
  Int_t    fNCells;             // Number of cells in cluster
  UInt_t   fPIDBit;             // PID Bit
  Double_t fDistToBad ;         // Distance to nearest bad channel
  Double_t fEmcCpvDistance;     // Distance from PHOS EMC rec.point to the closest CPV rec.point
  Double_t fM02;                // lambda0
  Double_t fM20;                // lambda1
  Double_t fTOF;

  ClassDef(AliCaloClusterInfo,1);
};

#endif // #ifdef ALICALOCLUSTERINFO_H
