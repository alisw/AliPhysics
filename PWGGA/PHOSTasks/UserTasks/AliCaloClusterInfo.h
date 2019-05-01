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
#include <TArrayI.h>
#include <TLorentzVector.h>
#include <TParticle.h>

class AliStack;
class AliESDEvent;
class AliVCluster;

class AliCaloClusterInfo : public TObject{
 public:

  AliCaloClusterInfo();
  AliCaloClusterInfo(AliVCluster* const clust, Int_t relID[4]);
  AliCaloClusterInfo(const AliCaloClusterInfo &src);
  AliCaloClusterInfo& operator=(const AliCaloClusterInfo &src);
  virtual ~AliCaloClusterInfo();

  TLorentzVector LorentzVector() const { return fLorentzVector; }

  Int_t    GetModule()          const { return fModule;          }
  Int_t    GetTRUNumber()       const { return fTRUNumber;       }
  Int_t    GetNCells()          const { return fNCells;          }
  TArrayI* GetLabelsArray()     const { return fLabels;          }
  Int_t    GetLabel()           const { if (fLabels && fLabels->GetSize() >0)        return fLabels->At(0);     else return -1;   }
  Int_t    GetLabelAt(UInt_t i) const { if (fLabels && i<(UInt_t)fLabels->GetSize()) return fLabels->At(i);     else return -999; }
  UInt_t   GetNLabels()         const { if (fLabels)                                 return fLabels->GetSize(); else return (0);  }

  UInt_t   GetPIDBit()          const { return fPIDBit;          }
  Double_t GetDistToBad()       const { return fDistToBad;       }
  Double_t GetM02()             const { return fM02;             }
  Double_t GetM20()             const { return fM20;             }
  Double_t GetTOF()             const { return fTOF;             }

  void SetLorentzVector(TLorentzVector momentum) { fLorentzVector = momentum; }
  void SetLabels(UInt_t size, Int_t* labels)     { if (fLabels) delete fLabels; fLabels = new TArrayI(size, labels); }

  Bool_t IsInFiducialRegion(Int_t cellX, Int_t cellZ);
  Bool_t IsMergedClusterFromPi0(AliStack* const stack, Int_t &pi0Indx);
  Bool_t IsClusterFromCvtedPi0(AliStack* const stack, Bool_t &isConverted, Int_t &pi0Indx);

 private:

  void FillCaloClusterInfo(AliVCluster* const clust, Int_t relID[4]);
  void SetPIDBit(UInt_t bit)   { fPIDBit |= bit; }
  Double_t TestDisp();

  Bool_t IsClusterFromPi0(AliStack* const stack, Int_t label, Int_t &pi0Indx);
  Bool_t IsClusterFromPi0Pure(AliStack* const stack, Int_t label, Int_t &pi0Indx);
  Bool_t IsClusterFromPi0Converted(AliStack* const stack, Int_t label, Int_t &pi0Indx);
  Int_t  GetTRUNumber(Int_t cellX, Int_t cellZ);

  TLorentzVector fLorentzVector;

  Int_t    fModule;
  Int_t    fTRUNumber;          // TRU Number
  Int_t    fNCells;             // Number of cells in cluster
  UInt_t   fPIDBit;             // PID Bit
  TArrayI* fLabels;             // list of primaries that generated the cluster, ordered in deposited energy
  Double_t fDistToBad;          // Distance to nearest bad channel
  Double_t fM02;                // lambda0
  Double_t fM20;                // lambda1
  Double_t fTOF;

  ClassDef(AliCaloClusterInfo, 3);
};

#endif // #ifdef ALICALOCLUSTERINFO_H
