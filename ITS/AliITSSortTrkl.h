#ifndef ALIITSSORTTRKL_H 
#define ALIITSSORTTRKL_H 

/* Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Helper class for finding multiple primary vertices       //
//           To be used by AliITSVertexer3D                           //
//           Origin M. Masera masera@to.infn.it                       //
////////////////////////////////////////////////////////////////////////


#include<TBits.h>
#include "AliLog.h"
#include "AliITSTracklPairs.h"

class TClonesArray;

class AliITSSortTrkl : public TObject {

 public:

  AliITSSortTrkl();
  AliITSSortTrkl(Int_t n, Double_t cut = 0.05);
  AliITSSortTrkl(TClonesArray &tclo, Int_t n, Double_t cut, Double_t rcut);
  virtual ~AliITSSortTrkl();
  Int_t AddPairs(Int_t t1, Int_t t2, Double_t dca, Double_t *coo);
  Int_t GetIndex() const {return fIndex;}
  Int_t FindClusters();
  void SetCut(Double_t cut){fCut = cut;}
  Double_t GetCut() const {return fCut; }
  Int_t* GetClusters(Int_t index) const {if(index>=0 && index<fNoClus){return fClusters[index];} else {return NULL;}}
  Int_t GetNumberOfClusters() const {return fNoClus;}
  Int_t GetSizeOfCluster(Int_t index) const {if(index>=0 && index<fNoClus){return fSize[index];} else {return -1;}}
  static void SortAndClean(Int_t numb, Int_t *arr, Int_t& numb2);
  Int_t* GetTrackletsLab(Int_t index, Int_t& dim) const;

  // FOR DEBUGGING PURPOSES
  Int_t* GetClustersTmp(Int_t index){return fClustersTmp[index];}
  AliITSTracklPairs* GetPairsAt(Int_t i) const {if(!(i>=0 && i<fIndex)){AliError(Form("Index %d out of bounds",i)); return NULL;} else{ return fPairs[i];} }


 protected:

  AliITSSortTrkl(const AliITSSortTrkl& pa);
  AliITSSortTrkl& operator=(const AliITSSortTrkl& /* pa */);
  void Cleanup();
  void DeleteClustersTmp();
  void PrepareClustersTmp();
  void Clustering(Int_t i, Int_t *v);
  Int_t* AliITSSortTrkl::FindLabels(Int_t *v, Int_t dimmax, Int_t& dim) const;

  const Int_t fkSize;         // Maximum number of tracklet pairs
  Int_t fIndex;               // Total number of tracklet pairs (<=fkSize)
  AliITSTracklPairs **fPairs;  // array of tracklet pairs (pointers to)
  Int_t **fClustersTmp;      // Temporary list of clusters of tracklet pairs
  Int_t **fClusters;      // List of clusters of tracklet pairs after cleanup
  Int_t fNoClus;         // Number of clusters of tracklet pairs
  Int_t *fSize;          // Number of pairs for each cluster
  TBits fMark;           // Used to mask used pairs
  Double_t fCut;         // cut on distance of DCAs of pairs for association
  Double_t fCoarseMaxRCut;  // cut on distance from beam axis

 ClassDef(AliITSSortTrkl,0);
};

#endif
