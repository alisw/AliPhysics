#ifndef CPV_H
#define CPV_H
////////////////////////////////////////////////
//  Manager and hits classes for set:CPV      //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 17 September 1999          //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include <TArray.h> 
#include <TRandom.h> 
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>

// --- galice header files ---
#include "AliDetector.h"
#include "AliHit.h"
#include "AliRun.h"

//==============================================================================
//                              AliCPVExactHit
//==============================================================================

class AliCPVExactHit : public TObject {
  
public:
  TLorentzVector fMomentum;   // 4-momentum of the particle
  Float_t        fXYhit[2];   // Hit's X,Y-coordinates
  Int_t          fIpart;      // Hit's particle type
  
public:
  virtual ~AliCPVExactHit() {}
           AliCPVExactHit() {}
           AliCPVExactHit(TLorentzVector p, Float_t *xy, Int_t ipart);

  TLorentzVector GetMomentum()  {return  fMomentum; }
  Float_t        GetXY(Int_t i) {return  fXYhit[i]; }
  Int_t          GetIpart()     {return  fIpart;    }
  void           Print();

  ClassDef(AliCPVExactHit,1)  // Hits object for set:CPV
};
 
//==============================================================================
//                              AliCPVHit
//==============================================================================

class AliCPVHit : public TObject {
  
public:
  Float_t        fXYhit[2];   // Hit's X,Y-coordinates
  
public:
  virtual ~AliCPVHit() {}
           AliCPVHit() {}
           AliCPVHit(Float_t *xy);

  Float_t  GetXY(Int_t i) {return  fXYhit[i]; }
  void     Print();

  ClassDef(AliCPVHit,1)  // Hits object for set:CPV
};
 
//==============================================================================
//                              AliCPVCradle
//==============================================================================

class AliCPVCradle : public TObject {

public:

  virtual   ~AliCPVCradle(void);
             AliCPVCradle(void);
             AliCPVCradle(Int_t   Geometry  ,
                          Float_t PadZSize  ,
                          Float_t PadPhiSize,
                          Float_t Radius    ,
                          Float_t Thickness ,
                          Float_t Angle     ,
                          Int_t   Nz        ,
                          Int_t   Nphi     );

  Float_t  GetPadZSize    (void) const {return fPadZSize;  }
  Float_t  GetPadPhiSize  (void) const {return fPadPhiSize;}
  Float_t  GetRadius      (void) const {return fRadius;    }
  Float_t  GetThickness   (void) const {return fThickness; }
  Int_t    GetNz          (void) const {return fNz;        }
  Int_t    GetNphi        (void) const {return fNphi;      }
  Float_t  GetAngle       (void) const {return fAngle;     }

  void     AddHit(TLorentzVector p, Float_t *xy, Int_t ipart);
  void     Clear(Option_t *opt="");
  void     Print(Option_t *opt="");
  
  void     Reconstruction(Float_t min_distance, Float_t min_signal);
  
  TClonesArray *GetHitExact         (void) {return fHitExact;        }
  TClonesArray *GetHitReconstructed (void) {return fHitReconstructed;}

private:

  Float_t       fPadZSize;          // Pad size in beam direction in cm
  Float_t       fPadPhiSize;        // Pad size in phi direction in cm
  Float_t       fRadius;            // Distance from IP to CPV in cm
  Float_t       fThickness;         // CPV thickness in cm
  Float_t       fAngle;             // Position of CRADLE center in degrees
  Int_t         fNz;                // Cells amount in beam direction
  Int_t         fNphi;              // Cells amount in phi direction
  
  TClonesArray *fHitExact;          // List of exact hits in the cradle
  TClonesArray *fHitReconstructed;  // List of reconstructed hits inthe cradle
  Int_t         fNhitsExact;        // Number of exact hits
  Int_t         fNhitsReconstructed;// Number of reconstructed hits
  
  ClassDef(AliCPVCradle,1)          // CPV cradle
};

//==============================================================================
//                            AliCPV
//==============================================================================

class AliCPV : public AliDetector {

public:

                AliCPV();
                AliCPV(const char *name, const char *title);
  virtual      ~AliCPV();

  virtual void  Init           ();
  virtual void  SetGeometry    (Int_t ncradles, Int_t nz, Int_t nphi, Float_t angle);
  virtual void  BuildGeometry  ();
  virtual void  CreateGeometry () {}
  virtual void  CreateMaterials();
  void          FinishEvent    (void);
  void          FinishRun      (void);

  void          ResetDigits    (void);
  void          Print          (Option_t *opt="");

  AliCPVCradle &GetCradle(int n) {return *(AliCPVCradle*)fCradles->operator[](n);}
  void          Reconstruction(Float_t min_distance, Float_t min_signal);

  virtual void  StepManager()=0;
  virtual void  AddCPVCradles();

  virtual Int_t GetCPVIdtmed (void){Int_t *idtmed = fIdtmed->GetArray()-1999;
				    return idtmed[2000];}

  Float_t  GetPadZSize          (void) const {return fPadZSize;  }
  Float_t  GetPadPhiSize        (void) const {return fPadPhiSize;}
  Float_t  GetRadius            (void) const {return fRadius;    }
  Float_t  GetThickness         (void) const {return fThickness; }
  Int_t    GetNz                (void) const {return fNz;        }
  Int_t    GetNphi              (void) const {return fNphi;      }
  Int_t    GetNofCradles        (void) const {return fNCradles;  }
  Float_t  GetAngle             (void) const {return fAngle;     }

  Int_t                 fDebugLevel;

private:

  Float_t               fPadZSize;           // Pad size along beam       [cm]
  Float_t               fPadPhiSize;         // Pad size across beam      [cm]
  Float_t               fRadius;             // Distance of CPV from IP   [cm]
  Float_t               fThickness;          // CPV thickness             [cm]
  Float_t               fAngle;              // Angle between CPV cradles [deg]
  Int_t                 fNz;                 // Number of pads along beam
  Int_t                 fNphi;               // Number of pads across beam

  Int_t                 fNCradles;           // Number of cradles
  TClonesArray         *fCradles;            // Array of CPV cradles

  ClassDef(AliCPV,1)  // Detector CPV
};
 
#endif
