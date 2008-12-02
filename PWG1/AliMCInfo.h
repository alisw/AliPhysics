#ifndef ALIMCINFO_H
#define ALIMCINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliGenInfo                               //
//   collect together MC info for comparison purposes - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliTPCdigitRow
//
////////////////////////////////////////////////////////////////////////

#include <TParticle.h>
#include "AliTrackReference.h"

class TFile;
class AliRunLoader;
class AliStack;


class AliTPCdigitRow: public TObject {
public:
  AliTPCdigitRow();
  virtual ~AliTPCdigitRow(){;}
  void SetRow(Int_t row);
  Bool_t TestRow(Int_t row) const ;
  AliTPCdigitRow & operator=(const AliTPCdigitRow &digOld);
  Int_t RowsOn(Int_t upto=8*32) const;
  Int_t Last() const;
  Int_t First() const ;
  void Reset();

private:
  UChar_t fDig[32];   // bitmask of the digits presence
  ClassDef(AliTPCdigitRow,1)  // container for digit pattern
};


////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliMCInfo
//
////////////////////////////////////////////////////////////////////////

class AliMCInfo: public TObject {
  friend class  AliGenInfoMaker;
  friend class  AliRecInfoMaker;
  friend class  AliESDRecInfo;
public:
  AliMCInfo();
  ~AliMCInfo();   
  AliMCInfo(const AliMCInfo& info);
  AliMCInfo& operator=(const AliMCInfo& info);
  //
  void Update();
  Int_t     GetEventNr() const   {return fEventNr;}
  const AliTrackReference&  GetTrackRef() const {return fTrackRef;}
  const AliTrackReference&  GetTrackRefOut() const {return fTrackRefOut;}
  const AliTrackReference&  GetTRdecay() const {return fTRdecay;} 
  TParticle& GetParticle()   {return fParticle;}
  //
  Int_t     GetPrimPart() const  {return fPrimPart;}
  Float_t   GetMass()   const    {return fMass;}                  
  Float_t   GetCharge() const    {return fCharge;}
  Int_t     GetLabel()  const    {return fLabel;}

  Int_t     GetMCtracks() const  {return fMCtracks;}
  Int_t     GetPdg()      const  {return fPdg;}
  const Float_t*   GetDecayCoord() const {return fDecayCoord;}
  const Double_t*  GetVDist()      const {return fVDist;}

  Bool_t   IsTPCdecay() const   {return fTPCdecay;}

  Int_t    GetRowsWithDigitsInn() const {return fRowsWithDigitsInn;}
  Int_t    GetRowsWithDigits() const  {return fRowsWithDigits;}
  Int_t    GetRowsTrackLength() const {return fRowsTrackLength;}
  Float_t GetPrim() const { return fPrim;}
  
  AliTPCdigitRow & GetTPCRow() {return fTPCRow;}
  Int_t GetNTPCRef() const {return fNTPCRef;}      
  Int_t GetNITSRef() const {return fNITSRef;}
  Int_t GetNTRDRef() const {return fNTRDRef;}
  Int_t GetNTOFRef() const {return fNTOFRef;}
  const TClonesArray *GetTPCReferences() const { return fTPCReferences;}  
  const TClonesArray * GetTRDReferences() const { return fTRDReferences;}  
  const TClonesArray * GetITSReferences() const { return fITSReferences;}  
  const TClonesArray * GetTOFReferences() const { return fTOFReferences;}  
  void CalcTPCrows(TClonesArray *arrayTR);
private:
  AliTrackReference  fTrackRef;      // track reference saved in the output tree
  AliTrackReference  fTrackRefOut;   // decay track reference saved in the output tree
  AliTrackReference  fTRdecay;       // track reference at decay point
  //
  Int_t     fPrimPart;               // index of primary particle in TreeH
  TParticle fParticle;               // generated particle 
  Float_t   fMass;                   // mass of the particle
  Float_t   fCharge;                 // charge of the particle
  Int_t     fLabel;                  // track label
  Int_t     fEventNr;                // event number
  Int_t     fMCtracks;               // indication of how many times the track is retuturned back
  Int_t fPdg;                        //pdg code
  Float_t fDecayCoord[3];            // position of particle decay
  Double_t fVDist[4];                //distance of the particle vertex from primary vertex
  Bool_t fTPCdecay;                  //indicates decay in TPC
  //
  // TPC row information using digits
  Int_t fRowsWithDigitsInn;          // number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;             // number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;            // last - first row with digit
  //
  // TPC track refernce information
  Float_t fTPCtrackLength;           // distance between first and last track reference
  //
  Float_t fPrim;                     // theoretical dedx in tpc according particle momenta and mass
  AliTPCdigitRow fTPCRow;                  // information about digits row pattern
  Int_t fNTPCRef;                    // tpc references counter
  Int_t fNITSRef;                    // ITS references counter
  Int_t fNTRDRef;                    // TRD references counter
  Int_t fNTOFRef;                    // TOF references counter
  TClonesArray * fTPCReferences;     //containner with all track references -in the TPC
  TClonesArray * fITSReferences;     //container with ITS references
  TClonesArray * fTRDReferences;     //container with TRD references  
  TClonesArray * fTOFReferences;     //container with TRD references  
  //
  ClassDef(AliMCInfo,1);  // container for 
};



#endif
