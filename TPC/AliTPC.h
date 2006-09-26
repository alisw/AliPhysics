#ifndef ALITPC_H
#define ALITPC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for TPC                     //
////////////////////////////////////////////////

class TFile;
class TTree;
#include <TMatrixFfwd.h>
#include <TVector.h>

class AliTPCDigitsArray;
class AliTPCLoader;
class AliTPCParam;
class AliTPCTrackHits; // M.I.  -MI4 old hits - to be removed later
class AliTPCTrackHitsV2; // M.I.
class AliRawReader;
#include "AliDetector.h"
#include "AliDigit.h" 
#include "AliHit.h" 


class AliTPC : public AliDetector {

public:
  AliTPC(); 
  AliTPC(const char *name, const char *title);
  AliTPC(const AliTPC& t);
  
  virtual AliLoader* MakeLoader(const char* topfoldername);
  
  virtual      ~AliTPC();
  virtual void  AddHit(Int_t a1, Int_t *a2, Float_t *a3);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials();
  virtual void  AddAlignableVolumes() const {}
  
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual void  SDigits2Digits(){;} //MI change -cycling to production
  virtual void  SDigits2Digits2(Int_t eventnumber=0);

  virtual void  Hits2SDigits(); // MI - cycling around
  virtual void  Hits2SDigits2(Int_t eventnumber=0);

  virtual void  Hits2Digits();
  virtual void  Hits2Digits(Int_t eventnumber);   //MI change
  virtual void  Hits2DigitsSector(Int_t isec);  //MI change
  virtual void  Init();
  virtual Int_t IsVersion() const =0;
  virtual void  Digits2Raw();
  virtual Bool_t Raw2SDigits(AliRawReader* rawReader);
  Int_t         GetNsectors() const  {return fNsectors;}
  virtual void  MakeBranch(Option_t *opt=" ");
  virtual void  ResetDigits();
  virtual void  SetSens(Int_t sens);
  virtual void  SetSide(Float_t side);


  virtual void  StepManager()=0;
  virtual void  DrawDetector() {}
  AliTPCDigitsArray*  GetDigitsArray() {return fDigitsArray;} //MI change
  AliTPCParam *GetParam(){return fTPCParam;} // M.K, M.I changes
  void SetParam(AliTPCParam *param){fTPCParam=param;} // M.K, M.I changes
  void SetDigitsArray(AliTPCDigitsArray* param) {fDigitsArray=param;}  //MI change

// additional function neccesary for the new hits 
   virtual void MakeBranch2(Option_t *opt=" ", const char *file=0);  //
   virtual void SetTreeAddress();
   virtual void SetTreeAddress2();
   virtual void AddHit2(Int_t a1,  Int_t *a2, Float_t *a3);  //
   virtual void ResetHits();
   virtual void ResetHits2();     
   virtual AliHit* FirstHit(Int_t track);
   virtual AliHit* NextHit();
   virtual AliHit* FirstHit2(Int_t track);
   virtual AliHit* NextHit2();
   virtual void LoadPoints(Int_t dummy);
   virtual void LoadPoints2(Int_t dummy);
   virtual void LoadPoints3(Int_t dumy);
   virtual void FinishPrimary();
   virtual void RemapTrackHitIDs(Int_t *map);
   void SetHitType(Int_t type){fHitType =type;} //set type of hit container
   void SetDigitsSwitch(Int_t sw){fDigitsSwitch = sw;}
   void SetDefSwitch(Int_t def){fDefaults = def;}
   Float_t GetNoise();  //get Current noise  
   void    GenerNoise(Int_t tablasize);  // make noise table
   Bool_t  IsSectorActive(Int_t sec) const;    // check if the sector is active
   void    SetActiveSectors(Int_t * sectors, Int_t n);  //set active sectors
   Int_t GetHitType() const {return fHitType;}
   void    SetActiveSectors(Int_t flag=0); //loop over al hits and set active only hitted sectors
   Bool_t  TrackInVolume(Int_t id,Int_t track);  //return true if current track is in volume

// static functions
   static AliTPCParam* LoadTPCParam(TFile *file); 
protected:
   Int_t          fDefaults; // defaults switch
  Int_t          fSens;             // ISENS
  Int_t          fNsectors;         // Number of sectors in TPC
  //MI changes
  AliTPCDigitsArray * fDigitsArray;              //! detector digit object  
  AliTPCParam *fTPCParam;           // pointer to TPC parameters 
  AliTPCTrackHitsV2 *fTrackHits;      //! hits for given track M.I.
  //  AliTPCTrackHits *fTrackHitsOld;      //! hits for given track M.I. MIold -

  Int_t  fHitType; // if fNewHit = 1 old data structure if 2 new hits  if 4  old MI stucture
  //  3 both types 
  Int_t fDigitsSwitch; // digits type, 0->normal, 1->summable

  //MK changes

  Float_t        fSide;  // selects left(-1), right(+1), or both(0) sides of the TPC
   
 private:
  //
  void SetDefaults();
  void DigitizeRow(Int_t irow,Int_t isec,TObjArray **rowTriplet);
  Float_t GetSignal(TObjArray *p1, Int_t ntr, TMatrixF *m1, 
                   TMatrixF *m2,Int_t *IndexRange);
  void GetList (Float_t label,Int_t np,TMatrixF *m,Int_t *IndexRange,
                Float_t **pList);
  void MakeSector(Int_t isec,Int_t nrows,TTree *TH,Stat_t ntracks,TObjArray **row);
  void TransportElectron(Float_t *xyz, Int_t *index);
  Int_t fCurrentIndex[4];// index[0] indicates coordinate system, 
                         // index[1] sector number, 
                         // index[2] pad row number  
                         // index[3] pad row number for which signal is calculated
  Int_t      fNoiseDepth;  //! noise table
  Float_t *  fNoiseTable;  //! table with noise
  Int_t      fCurrentNoise; //! index of the noise in  the noise table 
  Bool_t*    fActiveSectors; //! bool indicating which sectors are active

  ClassDef(AliTPC,12)  // Time Projection Chamber class
};


//_____________________________________________________________________________
 
class AliTPChit : public AliHit {
public:
   Int_t     fSector;     //sector number
   Int_t     fPadRow;     //Pad Row number
   Float_t   fQ ;         //charge
   Float_t   fTime;       //hit time
 
public:
   AliTPChit() {fTime = 0.;}
   AliTPChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
   virtual ~AliTPChit() {}
   void SetX(Float_t x){fX = x;}
   void SetY(Float_t y){fY = y;}
   void SetZ(Float_t z){fZ = z;}

   Float_t Time() const {return fTime;}
 
   ClassDef(AliTPChit,2)  // Time Projection Chamber hits
};


#endif









