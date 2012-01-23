/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice      
 */

//                                                                      
// Class for AliHFEcollection                                           
// Serves as a data container - currently based on internal TList      
//                                                                      
// Authors:                                                             
//   Markus Fasel  <M.Fasel@gsi.de>                                    
//   Matus Kalisky <matus.kalisky@cern.ch>  (contact)                   
//


//
// Provides an option for storing and creating histograms outside the
// analysis class
// the performance will be improved once the TMap is used insted of TTree
//

/*
 * vesion: 1.0.1
 */


#ifndef ALIHFECOLLECTION_H
#define ALIHFECOLLECTION_H

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#ifndef ROOT_THashList
#include "THashList.h"
#endif

class TCollection;
class TBrowser;

class AliHFEcollection : public TNamed{

 public:
  AliHFEcollection();
  AliHFEcollection(const char* name, const char* title);
  AliHFEcollection(const AliHFEcollection &c);
  AliHFEcollection &operator=(const AliHFEcollection &c);
  virtual ~AliHFEcollection();
  

  virtual void Browse(TBrowser *b);
  virtual Bool_t IsFolder() const { return kTRUE; }

  // Set & Create functions
  Bool_t CreateTH1F(const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax, Int_t logAxis = -1);
  Bool_t CreateTH1Farray(const char* name, const char* title, Int_t nBin, const Double_t* xbins);

  Bool_t CreateTH2F(const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY, Int_t logAxis = -1);
  Bool_t CreateTH3F(const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY, Int_t nBinZ, Float_t minZ, Float_t maxZ, Int_t logAxis = -1);

  Bool_t CreateTH1Fvector1(Int_t X, const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax, Int_t logAxis = -1);
  Bool_t CreateTH1Fvector2(Int_t X, Int_t Y, const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax, Int_t logAxis = -1);
  Bool_t CreateTH2Fvector1(Int_t X, const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY, Int_t logAxis = -1);
  Bool_t CreateProfile(const char* name, const char* title, Int_t nbins, Double_t xmin, Double_t xmax);
  Bool_t CreateTHnSparse(const char* name, const char* title, Int_t dim, const Int_t* nbins, const Double_t* xmin, const Double_t* xmax);

  Bool_t BinLogAxis(const char* name, Int_t dim);
  Bool_t Sumw2(const char*name);
    

  Long64_t Merge(const TCollection *list);
  virtual void Print(Option_t *) const;

  // Get functions
  TList* GetList() const { return fList; }
  TObject* Get(const char* name); 

  // Fill functions
  Bool_t Fill(const char* name, Double_t v);
  Bool_t Fill(const char* name, Int_t v);
  Bool_t Fill(const char* name, Int_t X, Double_t v);
  Bool_t Fill(const char* name, Int_t X, Int_t Y, Double_t v);
  Bool_t Fill(const char* name, Double_t v1, Double_t v2);
  Bool_t Fill(const char* name, Int_t X, Double_t v1, Double_t v2);
  Bool_t Fill(const char* name, Double_t v1, Double_t v2, Double_t v3);
  Bool_t Fill(const char* name, Double_t* entry, Double_t weight = 1);
 private:
  Bool_t CheckObject(const char* name);
   void Copy(TObject &ref) const;

 private:
  THashList*                           fList;      // Object container

  ClassDef(AliHFEcollection, 1)

};

#endif
