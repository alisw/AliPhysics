#ifndef THIGZ_H 
#define THIGZ_H 
//////////////////////////////////////////////// 
//  Emulation of HIGZ for Root
//////////////////////////////////////////////// 
 
#include <TCanvas.h> 

class THIGZ : public TCanvas { 

public:
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
   virtual void   x3d(Option_t *option="");
        
   ClassDef(THIGZ,1)  //Emulation of HIGZ for Root
}; 

   R__EXTERN THIGZ *higz;

inline void THIGZ::Divide(Int_t, Int_t, Float_t, Float_t, Int_t) { }
inline void THIGZ::SetGrid(Int_t, Int_t) { }
inline void THIGZ::SetGridx(Int_t) { }
inline void THIGZ::SetGridy(Int_t) { }
inline void THIGZ::SetLogx(Int_t) { }
inline void THIGZ::SetLogy(Int_t) { }
inline void THIGZ::SetLogz(Int_t) { }
inline void THIGZ::SetTickx(Int_t) { }
inline void THIGZ::SetTicky(Int_t) { }
inline void THIGZ::x3d(Option_t *) { }
    
#endif 
