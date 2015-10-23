#ifndef ALIITSUCLUSTERLINES_H
#define ALIITSUCLUSTERLINES_H
/* Copyright(c) 2009-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/////////////////////////////////////////////////////////////////////////////
// Class intended to gather and compute the parameters of vertex candidate //
/////////////////////////////////////////////////////////////////////////////
 
#include <TObject.h>
#include <vector>

class AliStrLine;

using std::vector;


class AliITSUClusterLines : public TObject {
 public:
  AliITSUClusterLines();
  AliITSUClusterLines(UInt_t first, AliStrLine *firstL, UInt_t second, AliStrLine *secondL,Bool_t=kFALSE);
  virtual ~AliITSUClusterLines();
  
  void Add(UInt_t label, AliStrLine *line, Bool_t weight=kFALSE);
  void ComputeClusterCentroid();
  inline UInt_t GetSize() const { return fLabels.size(); }
  inline Int_t* GetLabels(UInt_t &n) {n=fLabels.size(); return &fLabels[0]; }
  inline void GetA(Float_t a[3]) { for(Int_t i=0; i<3; ++i) a[i]=fA[i]; }
  inline void GetA(Double_t a[3]) { for(Int_t i=0; i<3; ++i) a[i]=fA[i]; }
  inline void GetB(Float_t b[3]) { for(Int_t i=0; i<3; ++i) b[i]=fB[i]; }
  inline void GetB(Double_t b[3]) { for(Int_t i=0; i<3; ++i) b[i]=fB[i]; }
  void GetCovMatrix(Float_t cov[6]);
  void GetCovMatrix(Double_t cov[6]);
  inline void GetVertex(Float_t p[3]) { for(Int_t i=0; i<3; ++i) p[i]=fV[i]; }
  inline void GetVertex(Double_t p[3]) { for(Int_t i=0; i<3; ++i) p[i]=fV[i]; }
  inline void GetWeight(Float_t w[9]) { for(Int_t i=0; i<9; ++i) w[i]=fW[i]; }
  inline void GetWeight(Double_t w[9]) { for(Int_t i=0; i<9; ++i) w[i]=fW[i]; }
  //
  virtual Bool_t IsSortable()                 const {return kTRUE;}
  virtual Bool_t IsEqual(const TObject* obj)  const;
  virtual Int_t	 Compare(const TObject* obj)  const;

 protected:
  Double_t fA[6];         // AX=B 
  Double_t fB[3];         // AX=B
  vector<Int_t> fLabels;  // labels
  Double_t fV[3];         // vertex candidate
  Double_t fW[9];         // weight matrix
  ClassDef(AliITSUClusterLines,1);
};

#endif
