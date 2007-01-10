#ifndef ALIMUONDETELEMENT_H
#define ALIMUONDETELEMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONDetElement
/// \brief Detection element object containing information for combined 
/// cluster / track finder in the MUON arm 
///
/// \author Alexander Zinchenko, JINR Dubna
 
#include <TObject.h>
class TObjArray;
class TClonesArray;
class AliMUONDigit;
class AliMUONHitMapA1;
class AliMUONGeometrySegmentation;
class AliMUONData;
class AliMUONRawCluster;
class AliMUONClusterFinderAZ;

class AliMUONDetElement : public TObject 
{
 public:

  AliMUONDetElement();
  AliMUONDetElement(Int_t idDE, AliMUONDigit *dig, AliMUONClusterFinderAZ *recModel); // constructor
  virtual ~AliMUONDetElement(); // Destructor

  Int_t IdDE(void) const { return fidDE; } // det. elem. ID
  Int_t Chamber(void) const { return fChamber; } // chamber No
  Double_t Z(void) const { return fZ; } // Z-coordinate
  Int_t Left(Int_t cath) const { return fLeft[cath]; } // number of unused digits
  //Int_t GetMapElem(AliMUONDigit *digit) const; // get map element
  TObjArray *Digits(Int_t cath) const { return fDigits[cath]; } // array of digits
  TObjArray *RawClusters() const { return fRawClus; } // array of raw clusters
  TObjArray *HitsForRec() const { return fHitsForRec; } // hits for rec.
  Int_t NHitsForRec() const { return fNHitsForRec; } // No. of hits for rec.
  //Bool_t Inside(Double_t x, Double_t y, Double_t z) const; // check if point inside DE
  Bool_t Inside(Double_t x, Double_t y, Double_t z, Double_t dx, Double_t dy) const; // check if point inside DE
  
  void SetID(Int_t idDE) { fidDE = idDE; } // set det. elem. ID
  void SetIndex(Int_t index) { fIndex = index; } // set position index
  void SetZ(Double_t z) { fZ = z; } // set Z-coord.
  void SetLeft(Int_t cath, Int_t left) { fLeft[cath] = left; } // set No. of digits
  //void SetMapElem(const AliMUONDigit *digit, Int_t flag); // set map element 
  void AddDigit(AliMUONDigit *dig); // add digit
  void Fill(AliMUONData *data); // fill hit maps 
  void ClusterReco(Double_t x, Double_t y, Double_t dx, Double_t dy); // run cluster reco around (x,y)
  void AddHitForRec(AliMUONRawCluster *clus); // make HitForRec
  // What is necessary for sorting TObjArray's
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* detElem) const; // "Compare" function for sorting

 private:
 
  Int_t fidDE; ///< det. elem. ID
  Int_t fIndex; ///< det. elem. position index in container
  Int_t fChamber; ///< chamber No
  Double_t fZ; ///< det. elem. Z-coordinate
  Int_t fLeft[2]; ///< numbers of digits not used for clustering
  Int_t fNHitsForRec; ///< number of hits for rec.
  AliMUONGeometrySegmentation* fSeg[2]; ///< segmentation
  AliMUONHitMapA1 *fHitMap[2]; ///< map of digits
  TObjArray *fDigits[2]; ///< container of digits from this det. elem.
  TObjArray *fRawClus; ///< raw clusters
  TObjArray *fHitsForRec; ///< HitForRec's
  AliMUONClusterFinderAZ *fRecModel; ///< cluster finder

  // Functions
  AliMUONDetElement(const AliMUONDetElement & rhs); // copy constructor
  AliMUONDetElement& operator = (const AliMUONDetElement& rhs); // assignment operator

  ClassDef(AliMUONDetElement,0) // detection element object
    };
#endif
