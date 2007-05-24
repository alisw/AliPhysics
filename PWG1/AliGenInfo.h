#ifndef ALIGENINFO_H
#define ALIGENINFO_H
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
class AliTPCParam;

const Int_t kgRowBytes = 32;

class AliTPCdigitRow: public TObject {
public:
  AliTPCdigitRow();
  virtual ~AliTPCdigitRow(){;}
  void SetRow(Int_t row);
  Bool_t TestRow(Int_t row) const ;
  AliTPCdigitRow & operator=(const AliTPCdigitRow &digOld);
  Int_t RowsOn(Int_t upto=8*kgRowBytes) const;
  Int_t Last() const;
  Int_t First() const ;
  void Reset();

private:
  UChar_t fDig[kgRowBytes];
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
  void Update();
  Int_t     GetEventNr() const   {return fEventNr;}
  const AliTrackReference&  GetTrackRef() const {return fTrackRef;}
  const AliTrackReference&  GetTrackRefOut() const {return fTrackRefOut;}
  const AliTrackReference&  GetTRdecay() const {return fTRdecay;} 
  TParticle& GetParticle()   {return fParticle;}
  Float_t TPCBetheBloch(Float_t bg);
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
private:
  AliTrackReference  fTrackRef;      // track reference saved in the output tree
  AliTrackReference  fTrackRefOut;   // decay track reference saved in the output tree
  AliTrackReference  fTRdecay;       // track reference at decay point
  //
  Int_t     fPrimPart;               // index of primary particle in TreeH
  TParticle fParticle;               // generated particle 
  Float_t   fMass;                   // mass of the particle
  Float_t   fCharge;                 //
  Int_t     fLabel;                  // track label
  Int_t     fEventNr;                // event number
  Int_t     fMCtracks;               // indication of how many times the track is retuturned back
  Int_t fPdg;                        //pdg code
  Float_t fDecayCoord[3];            // position of particle decay
  Double_t fVDist[4];                //distance of the particle vertex from primary vertex
  Bool_t fTPCdecay;                  //indicates decay in TPC
  Int_t fRowsWithDigitsInn;          // number of rows with digits in the inner sectors
  Int_t fRowsWithDigits;             // number of rows with digits in the outer sectors
  Int_t fRowsTrackLength;            // last - first row with digit
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



class AliGenV0Info: public TObject {
public:
  AliGenV0Info();       //
  void Update(Float_t vertex[3]);       
  AliMCInfo &  GetPlus()      {return fMCd;}
  AliMCInfo &  GetMinus()     {return fMCm;}
  TParticle &  GetMopther()   {return fMotherP;}
  Double_t    GetMCDist1() const { return fMCDist1;}
  Double_t    GetMCDist2() const {return fMCDist2;}  
  const Double_t*  GetMCPdr() const {return fMCPdr;}
  const Double_t*  GetMCPd()  const {return fMCPd;}
  const Double_t*  GetMCX()  const {return fMCX;}
  //  const Double_t    fMCXr;
  //
//   Double_t     fMCPm[3];    
//   Double_t     fMCAngle[3]; 
//   Double_t     fMCRr;       
//   Double_t     fMCR;       
//   Int_t        fPdg[2];   
//   Int_t        fLab[2];   
//   //
//   Double_t       fInvMass;  
//   Float_t        fPointAngleFi;
//   Float_t        fPointAngleTh;
//   Float_t        fPointAngle;  

  void SetInfoP(AliMCInfo &plus) {fMCd=plus;}
  void SetInfoM(AliMCInfo &minus){fMCm=minus;}
  void SetMother(TParticle&mother){fMotherP=mother;}
private:
  AliMCInfo   fMCd;       //info about daughter particle - second particle for V0
  AliMCInfo   fMCm;       //info about mother particle   - first particle for V0
  TParticle   fMotherP;   //particle info about mother particle
  Double_t    fMCDist1;    //info about closest distance according closest MC - linear DCA
  Double_t    fMCDist2;    //info about closest distance parabolic DCA
  //
  Double_t    fMCPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t    fMCPd[4];     //exact momentum from MC info
  Double_t    fMCX[3];      //exact position of the vertex
  Double_t    fMCXr[3];     //rec. position according helix
  //
  Double_t     fMCPm[3];    //momentum at the vertex mother
  Double_t     fMCAngle[3]; //three angels
  Double_t     fMCRr;       // rec position of the vertex 
  Double_t     fMCR;        //exact r position of the vertex
  Int_t        fPdg[2];   //pdg code of mother and daugter particles
  Int_t        fLab[2];   //MC label of the partecle  
  //
  Double_t       fInvMass;  //reconstructed invariant mass -
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  //
  ClassDef(AliGenV0Info,1)  // container for  
};



class AliGenKinkInfo: public TObject {
public:
  AliGenKinkInfo();          //default cosntructor
  void    Update();          // put some derived info to special field 
  Float_t GetQt();           //
  AliMCInfo &  GetPlus()      {return fMCd;}
  AliMCInfo &  GetMinus()     {return fMCm;}
  void SetInfoDaughter(AliMCInfo &daughter) {fMCd=daughter;}
  void SetInfoMother(AliMCInfo &mother){fMCm=mother;}
private:
  AliMCInfo   fMCd;          //info about daughter particle - second particle for V0
  AliMCInfo   fMCm;          //info about mother particle   - first particle for V0
  Double_t    fMCDist1;      //info about closest distance according closest MC - linear DCA
  Double_t    fMCDist2;      //info about closest distance parabolic DCA
  //
  Double_t     fMCPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t     fMCPd[4];     //exact momentum from MC info
  Double_t     fMCX[3];      //exact position of the vertex
  Double_t     fMCXr[3];     //rec. position according helix
  //
  Double_t     fMCPm[3];     //momentum at the vertex mother
  Double_t     fMCAngle[3];  //three angels
  Double_t     fMCRr;        // rec position of the vertex 
  Double_t     fMCR;         //exact r position of the vertex
  Int_t        fPdg[2];      //pdg code of mother and daugter particles
  Int_t        fLab[2];      //MC label of the partecle
  ClassDef(AliGenKinkInfo,1) // container for  
};

#endif
