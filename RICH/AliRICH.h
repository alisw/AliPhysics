#ifndef RICH_H
#define RICH_H
////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"

class AliRICH : public AliDetector {
  
protected:
  Int_t               fNmips;            //Number of mips in RICH
  Int_t               fNckovs;           //Number of cerenkovs in RICH
  Int_t               fNpadhits;         //Number of pad hits in RICH
  
  TClonesArray        *fMips;              //List of mips
  TClonesArray        *fCkovs;             //List of cerenkovs
  TClonesArray        *fPadhits;           //List of Padhits
  
  Float_t             fChslope;            //Charge slope
  Float_t             fAlphaFeed;          //Feed-back coefficient
  Float_t             fSxcharge;           //Charge slope along x
  Int_t               fIritri;             //Trigger flag

public:
  AliRICH();
  AliRICH(const char *name, const char *title);
  virtual              ~AliRICH();
  virtual void         AddHit(Int_t, Int_t*, Float_t*);
  virtual void         AddMipHit(Int_t, Int_t*, Float_t*);
  virtual void         AddCkovHit(Int_t, Int_t*, Float_t*);
  virtual void         AddPadHit(Int_t, Int_t*, Float_t*);
  virtual void         BuildGeometry();
  virtual void         CreateGeometry() {}
  virtual void         CreateMaterials() {}
  Int_t                DistancetoPrimitive(Int_t px, Int_t py);
  inline virtual int   GetNmips() {return fNmips;}
  inline virtual int   GetNckovs() {return fNckovs;}
  inline virtual int   GetNpadhits() {return fNpadhits;}
  virtual Int_t        IsVersion() const =0;
  virtual void         Init();
  inline TClonesArray  *Mips()    {return fMips;}
  inline TClonesArray  *Ckovs()    {return fCkovs;}
  inline TClonesArray  *Padhits()    {return fPadhits;}
  void                 FinishEvent(void){;} 
  virtual void         MakeBranch(Option_t *); 
  void                 SetTreeAddress(void);
  virtual void         StepManager();
  virtual void         PreTrack();
  
  virtual void         SetSP(Float_t chslope){ fChslope=chslope;}
  virtual void         SetFEED(Float_t alphafeed){fAlphaFeed=alphafeed;}
  virtual void         SetSIGM(Float_t sxcharge){fSxcharge=sxcharge;}
  virtual void         SetTRIG(Int_t iritri) {fIritri=iritri;}
  virtual void         ResetHits();
  virtual void         UpdateMipHit(Float_t*);
  virtual void         RichIntegration();
  virtual void         AnodicWires(Float_t &);
  virtual void         GetChargeMip(Float_t &);
  virtual void         GetCharge(Float_t &);
  virtual void         FeedBack(Float_t *, Float_t );
  virtual Float_t      FMathieson(Float_t , Float_t );
  
  ClassDef(AliRICH,1)  // Base class for RICH
};

class AliRICHv1 : public AliRICH {
  
public:
  AliRICHv1();
  AliRICHv1(const char *name, const char *title);
  virtual              ~AliRICHv1();
  virtual void          CreateGeometry();
  virtual void          CreateMaterials();
  virtual Int_t         IsVersion() const {return 1;}
  virtual void          DrawDetector();
  
  
  ClassDef(AliRICHv1,1)  // RICH version 1
};

//_____________________________________________________________________________
class AliRICHhit: public AliHit {
public:
  Int_t     fVolume[2];  //array of volumes
  
  //Pad informations
  Int_t     fFirstpad;   //First index in padhits
  Int_t     fLastpad;    //Last index in padhits
  
  //Hit information
  Int_t     fModule;     //Module number
  Float_t   fTheta;      //Theta of the particle generating the hit
  
  Float_t   fArrivaltime;// Time of hit.
  Int_t     fPart;       //Particle type
  
  // we don't know what is this for :
  Int_t     fFeed;       //Type of feedback (origin of charge deposition)
  
public:
  AliRICHhit() {}
  AliRICHhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits,
	     Int_t fNpadhits);
  virtual ~AliRICHhit(){}
  
  void SetLastpad(Int_t lastpad){fLastpad=lastpad;}
  
  ClassDef(AliRICHhit,1)  // Hits for set RICH
};

//_____________________________________________________________________________
class AliRICHmip: public AliRICHhit 
{
public:
  // Hit information keep
  Float_t   fPhi;        //Phi of the particle generating the hit
  Float_t   fPs;         //Momentum of the particle generating the hit
  Float_t   fQ;          //Charge of the particle
  
  // Generated cerenkov information (Z of generation stored in fZ of AliHit)
  Int_t     fFirstCkov;  //Index in the ckov TcloneArray of the first generated 
  Int_t     fLastCkov;   //Here the last.
  
public:
  AliRICHmip() {}
  AliRICHmip(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits,
	     Int_t fNckovs, Int_t fNpadhits);
  virtual ~AliRICHmip() {}
  
  Float_t GetZ()                 { return fZ;}
  void    SetX(Float_t x)        { fX        = x;     }
  void    SetY(Float_t y)        { fY        = y;     }
  void    SetZ(Float_t z)        { fZ        = z;     }
  void    SetLastCkov(Int_t last){ fLastCkov = last;  }
  void    SetModule(Int_t module){ fModule   = module;}
  void    SetTheta(Float_t theta){ fTheta    = theta; }
  void    SetPhi(Float_t phi)    { fPhi      = phi;   }
  
  ClassDef(AliRICHmip,1)  //Mip hits for RICH
};

//_____________________________________________________________________________
class AliRICHckov: public AliRICHhit 
{
public:
  // Hit information keep
  Float_t   fEnergy;     //Photon energy
  Int_t     fStop;       //Stop mechanism (cut, threshold, ...)
  
  //Parent info
  Int_t     fParent;     //Index in array of mips of parent which generatethis
public:
  AliRICHckov() {}
  AliRICHckov(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits,
	      Int_t fNmips, Int_t fNpadhits);
  virtual ~AliRICHckov() {}
  
  ClassDef(AliRICHckov,1) //Cerenkov hits for RICH
};

//_____________________________________________________________________________
class AliRICHpadhit: public AliHit
{
public:
  Int_t     fVolume[2];  //array of volumes
  
  // Hit information
  Int_t     fX;          //Integer x position in pad
  Int_t     fY;          //Integer y position in pad
  Int_t     fModule;     //Module number
  // Particle info
  Int_t     fParentMip;  //Parent particle
  Int_t     fParentCkov; //Parent CKOV
  // physics info
  Int_t     fProd;       //Production mechanism
  Float_t   fCharge;     //Charge deposited

public:
  AliRICHpadhit(){}
  AliRICHpadhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits,
		Int_t fNmips,Int_t fNckovs);
  virtual ~AliRICHpadhit() {}
  
  ClassDef(AliRICHpadhit,1) //Pad hits for RICH
};  

#endif



