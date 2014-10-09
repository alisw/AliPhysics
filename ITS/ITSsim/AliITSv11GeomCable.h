#ifndef ALIITSV11GEOMCABLE_H
#define ALIITSV11GEOMCABLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


class TGeoVolume;
class TGeoNode;

#include <TObjArray.h>

//*************************************************************************
//   Base class of cable classes
//
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************


class AliITSv11GeomCable : public TNamed {

 public:
  AliITSv11GeomCable();
  AliITSv11GeomCable(const char* name);

  virtual ~AliITSv11GeomCable();
  void SetDebug(Int_t debug = 1) {fDebug = debug;};

  void  SetInitialNode(TGeoVolume *vol);
  void  ResetInitialNode();

  void  AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt, Double_t *coord);
  virtual Int_t       GetNCheckPoints() const;
  virtual Int_t       GetPoint(Int_t iCheckPt, Double_t *coord) const;
  virtual Int_t       GetVect(Int_t iCheckPt, Double_t *coord) const;
  virtual TGeoVolume* GetVolume( Int_t iCheckPt ) const;

  virtual Int_t       GetCheckPoint( Int_t iCheckPt, Int_t nOccur,
				     Int_t motherLevel, Double_t *coord);
  virtual Int_t       GetCheckVect( Int_t iCheckPt, Int_t nOccur,
				    Int_t motherLevel, Double_t *coord);
  virtual Int_t       GetCheckVect( const Double_t *localCoord,
				    TGeoVolume *vol, Int_t nOccur,
				    Int_t motherLevel, Double_t *coord);
  void ResetPoints();

 protected:
  AliITSv11GeomCable(const AliITSv11GeomCable &source);
  AliITSv11GeomCable& operator=(const AliITSv11GeomCable &source);
  bool     CheckDaughter(const TGeoNode* node, Int_t i = 0);
  void     ResetCheckDaughter();
  void     CopyFrom(Double_t *c, const Double_t *o) const;
  Double_t ScalProd(const Double_t *a, const Double_t *b) const;

  static const Int_t fgkCableMaxNodeLevel = 50; // max. number of levels
  static const Int_t fgkCableMaxLayer = 15;     // max. number of layers

  Int_t fDebug;                         // debug flag
  Int_t fNodeInd[fgkCableMaxNodeLevel]; // index of nodes in the node tree
  TObjArray fPointArray;                // array of points
  TObjArray fVolumeArray;               // volumes containing the points
  TGeoVolume *fCurrentVol;              // volume to search in the node tree
  TGeoNode *fInitialNode;               // initial node to start searching

  ClassDef(AliITSv11GeomCable,1)
};

inline Int_t AliITSv11GeomCable::GetNCheckPoints() const{
  return fVolumeArray.GetEntriesFast(); }

inline void AliITSv11GeomCable::ResetCheckDaughter() {
  for (Int_t i=0; i<fgkCableMaxNodeLevel; i++) fNodeInd[i] = -1; }

inline void AliITSv11GeomCable::CopyFrom(Double_t *c, const Double_t *o)
const { *(c++)=*(o++); *(c++)=*(o++); *c=*o; }

inline Double_t AliITSv11GeomCable::ScalProd(const Double_t *a,
						const Double_t *b) const {
  Double_t s = *(a++)*(*(b++)); s+=*(a++)*(*(b++)); s+=*a*(*b);
  return s;
}





#endif
