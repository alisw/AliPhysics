#ifndef THIGZ_H 
#define THIGZ_H 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////// 
//  Emulation of HIGZ for Root
//////////////////////////////////////////////// 
 
#include <TCanvas.h> 

class THIGZ : public TCanvas { 

protected:
   Int_t     fFAIS;   //Fill Area Interior Style (0,1,2,3)
   Int_t     fFASI;   //Fill Area Style Index
   Int_t     fLTYP;   //Line TYPe
   Float_t   fBASL;   //BAsic Segment Length
   Float_t   fLWID;   //Line WIDth
   Int_t     fMTYP;   //Marker TYPe
   Float_t   fMSCF;   //Marker SCale Factor
   Int_t     fPLCI;   //PolyLine Color Index
   Int_t     fPMCI;   //PolyMarker Color Index
   Int_t     fFACI;   //Fill Area Color Index
   Int_t     fTXCI;   //TeXt Color Index
   Int_t     fTXAL;   //10*(alignment horizontal) + (alignment vertical)
   Float_t   fCHHE;   //CHaracter HEight)
   Float_t   fTANG;   //Text ANGle
   Int_t     fTXFP;   //10*(TeXt Font) + (TeXt Precision)
   Int_t     fBORD;   //Border for IGBOX, IGFBOX and IGARC (0=No , 1=Yes)
   Int_t     fNCOL;   //Number of entry in the COLor map.
   Int_t     fDRMD;   //Drawing mode: 1.=copy 2.=xor
   Int_t     fSYNC;   //Synchronise the graphics in X11 1.=yes 0.=no
   Int_t     fCLIP;   //Clipping mode: 1.=on 0.=off
   Int_t     f2BUF;   //10*(WKID)+(double buffer mode: 1.=on 0.=off)
   Int_t     fPID;    //integer identifier f current primitive
   TString   fPname;  //Name of current primitive ID
      
public: 
   THIGZ(); 
   THIGZ(Int_t size); 
   virtual ~THIGZ(); 
   Float_t        Get(const char *name);
   virtual  void  Reset(Option_t *option="");
   virtual  void  Set(const char *name, Float_t rval);
//
   virtual  void  Gdopt(const char *name,const char *value); // *MENU*
   virtual  void  Gdraw(const char *name,Float_t theta=30, Float_t phi=30, Float_t psi=0,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01); // *MENU*
   virtual  void  Gdrawc(const char *name,Int_t axis=1, Float_t cut=0,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01); // *MENU*
   virtual  void  Gdspec(const char *name); // *MENU*
   virtual  void  Gdtree(const char *name,Int_t levmax=15,Int_t ispec=0); // *MENU*
   virtual  void  Gsatt(const char *name, const char *att, Int_t val); // *MENU*
   virtual  void  SetBOMB(Float_t bomb=1); // *MENU*

   //dummies
   virtual void   Divide(Int_t nx=1, Int_t ny=1, Float_t xmargin=0.01, Float_t ymargin=0.01, Int_t color=0);
   virtual void   SetGrid(Int_t valuex = 1, Int_t valuey = 1);
   virtual void   SetGridx(Int_t value = 1);
   virtual void   SetGridy(Int_t value = 1);
   virtual void   SetLogx(Int_t value = 1);
   virtual void   SetLogy(Int_t value = 1);
   virtual void   SetLogz(Int_t value = 1);
   virtual void   SetTickx(Int_t value = 1);
   virtual void   SetTicky(Int_t value = 1);
   virtual void   X3d(Option_t *option="");
        
  virtual  Int_t     FAIS() const {return fFAIS;}
  virtual  Int_t     FASI() const {return fFASI;}
  virtual  Int_t     LTYP() const {return fLTYP;}
  virtual  Float_t   BASL() const {return fBASL;}
  virtual  Float_t   LWID() const {return fLWID;}
  virtual  Int_t     MTYP() const {return fMTYP;}
  virtual  Float_t   MSCF() const {return fMSCF;}
  virtual  Int_t     PLCI() const {return fPLCI;}
  virtual  Int_t     PMCI() const {return fPMCI;}
  virtual  Int_t     FACI() const {return fFACI;}
  virtual  Int_t     TXCI() const {return fTXCI;}
  virtual  Int_t     TXAL() const {return fTXAL;}
  virtual  Float_t   CHHE() const {return fCHHE;}
  virtual  Float_t   TANG() const {return fTANG;}
  virtual  Int_t     TXFP() const {return fTXFP;}
  virtual  Int_t     BORD() const {return fBORD;}
  virtual  Int_t     NCOL() const {return fNCOL;}
  virtual  Int_t     DRMD() const {return fDRMD;}
  virtual  Int_t     SYNC() const {return fSYNC;}
  virtual  Int_t     CLIP() const {return fCLIP;}
  virtual  Int_t     I2BUF() const {return f2BUF;}
  virtual  Int_t     PID() const {return fPID;}
  virtual  const char* Pname() const {return fPname.Data();}
  virtual  void      SetPID(Int_t pid) {fPID=pid;}
  virtual  void      SetFACI(Int_t faci) {fFACI=faci;}
  virtual  void      SetFAIS(Int_t fais) {fFAIS=fais;}
  virtual  void      SetLTYP(Int_t lin) {fLTYP=lin;}
  virtual  void      SetMTYP(Int_t mk) {fMTYP=mk;}
  virtual  void      SetLWID(Int_t lw) {fLWID=lw;}
  virtual  void      SetPLCI(Int_t lcol) {fPLCI=lcol;}
  virtual  void      SetPMCI(Int_t mcol) {fPMCI=mcol;}
  virtual  void      SetTXCI(Int_t tcol) {fTXCI=tcol;}
  virtual  void      SetPname(const char *name) {fPname=name;}

   ClassDef(THIGZ,1)  //Emulation of HIGZ for Root
}; 

   R__EXTERN THIGZ *gHigz;

inline void THIGZ::Divide(Int_t, Int_t, Float_t, Float_t, Int_t) { }
inline void THIGZ::SetGrid(Int_t, Int_t) { }
inline void THIGZ::SetGridx(Int_t) { }
inline void THIGZ::SetGridy(Int_t) { }
inline void THIGZ::SetLogx(Int_t) { }
inline void THIGZ::SetLogy(Int_t) { }
inline void THIGZ::SetLogz(Int_t) { }
inline void THIGZ::SetTickx(Int_t) { }
inline void THIGZ::SetTicky(Int_t) { }
inline void THIGZ::X3d(Option_t *) { }
    
#endif 
