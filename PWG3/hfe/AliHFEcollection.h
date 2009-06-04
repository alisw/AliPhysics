/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice      
 */

/************************************************************************
 *                                                                      *
 * Class for AliHFEcollection                                           *
 * Serves as a data container - currently based on internal TList       *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel  <M.Fasel@gsi.de>                                     *
 *   Matus Kalisky <matus.kalisky@cern.ch>  (contact)                   *
 ************************************************************************/

/*
 * Provides an option for storing and creating histograms outside the
 * analysis class
 * the performance will be improved once the TMap is used insted of TTree
 */

/*
 * vesion: 1.0.1
 */


#ifndef __ALIHFECOLLECTION_H__
#define __ALIHFECOLLECTION_H__

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TList;
class TCollection;
class TBrowser;

class AliHFEcollection : public TNamed{

 public:
  AliHFEcollection();
  AliHFEcollection(char* name, char* title);
  AliHFEcollection(const AliHFEcollection &c);
  AliHFEcollection &operator=(const AliHFEcollection &c);
  virtual ~AliHFEcollection();
  

  virtual void Browse(TBrowser*);

  // Set & Create functions
  Bool_t CreateTH1F(const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax);

  Bool_t CreateTH2F(const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY);

  Bool_t CreateTH1Fvector1(Int_t _X, const char* _name, const char* _title, Int_t _nBin, Float_t _nMin, Float_t _nMax);
  Bool_t CreateTH2Fvector1(Int_t _X, const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY);

  Bool_t CreateTH1Fvector2(Int_t _X, Int_t _Y, const char* _name, const char* _title, Int_t _nBin, Float_t _nMin, Float_t _nMax);
  

  Long64_t Merge(TCollection *list);

  // Get functions
  TList* GetList()  const  { return fListE; }
  TObject* Get(const char* name); 
  TObject* Get(const char* name, Int_t _X);
  TObject* Get(const char* name, Int_t _X, Int_t _Y);

 private:
  Bool_t CheckObject(const char*);
   void Copy(TObject &ref) const;

 private:
  TList*                           fListE;

  ClassDef(AliHFEcollection, 1)

};

#endif
