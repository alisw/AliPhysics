#ifndef PHOS_H
#define PHOS_H
////////////////////////////////////////////////
//  Manager and hits classes for set:PHOS     //
////////////////////////////////////////////////
 
// --- CLHEP ---
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>

// --- ROOT system ---
#include <TArray.h> 
#include <TRandom.h> 
#include <TH2.h>

// --- galice header files ---
#include "AliDetector.h"
#include "AliHit.h"
#include "AliRun.h"

class AliPHOSgamma : public TObject {

  public:
                        virtual ~AliPHOSgamma(void) {}
                        AliPHOSgamma(void) {}
                        AliPHOSgamma(const AliPHOSgamma &g) { *this=g; }
                        AliPHOSgamma(Float_t X, Float_t Xsigma,
                                     Float_t Y, Float_t Ysigma,
                                     Float_t E, Float_t Esigma,
                                     Float_t Px, Float_t Py, Float_t Pz) :
                          fX(X), fXsigma(Xsigma),
                          fY(Y), fYsigma(Ysigma),
                          fE(E), fEsigma(Esigma),
                          fPx(Px), fPy(Py), fPz(Pz)
                        {}

    Float_t             fX;             // cm. x-coordinate (in beam direction)
    Float_t             fXsigma;        // cm. x-coordinate error

    Float_t             fY;             // cm. y-coordinate (around beam)
    Float_t             fYsigma;        // cm. y-coordinate error

    Float_t             fE;             // GeV. energy
    Float_t             fEsigma;        // GeV. energy error

    Float_t             fPx;            // GeV. Gamma momentum Px
    Float_t             fPy;            // GeV. Gamma momentum Py
    Float_t             fPz;            // GeV. Gamma momentum Pz

    void                Print(Option_t *options=NULL);
    AliPHOSgamma       &operator=(const AliPHOSgamma &g);

 private:

    ClassDef(AliPHOSgamma,1)            // Gamma particle in PHOS cradle
};

//______________________________________________________________________________

class AliPHOShit : public AliHit {
   
public:
  Int_t     fVolume[5];  //array of volumes
  Float_t   fELOS;       //ELOS
 
public:
  AliPHOShit() {}
  AliPHOShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliPHOShit() {}
 
  ClassDef(AliPHOShit,1)  //Hits object for set:PHOS
};
 
//______________________________________________________________________________


class AliPHOSCradle : public TObject {

  public:

    virtual            ~AliPHOSCradle(void) {}
                        AliPHOSCradle(void);
                        AliPHOSCradle(int   Geometry           ,
                                      float CrystalSideSize    ,
                                      float CrystalLength      ,
                                      float WrapThickness      ,
                                      float AirThickness       ,
                                      float PIN_SideSize       ,
                                      float PIN_Length         ,
                                      float Radius             ,
                                      float CPV_Thickness      ,
                                      float CPV_PHOS_Distance  ,
                                      int   Nz                 ,
                                      int   Nphi               ,
                                      float Angle );

    void                AddCPVHit(float x, float y);

    Float_t             GetCrystalSideSize     (void) const {return fCrystalSideSize;}
    Float_t             GetCellSideSize        (void) const {return fCrystalSideSize+2*fWrapThickness+2*fAirThickness;}
    Float_t             GetCrystalLength       (void) const {return fCrystalLength;}
    Float_t             GetWrapThickness       (void) const {return fWrapThickness;}
    Float_t             GetAirThickness        (void) const {return fAirThickness;}
    Float_t             GetPIN_SideSize        (void) const {return fPIN_SideSize;}
    Float_t             GetPIN_Length          (void) const {return fPIN_Length;}
    Float_t             GetRadius              (void) const {return fRadius;}
    Float_t             GetCPV_PHOS_Distance   (void) const {return fCPV_PHOS_Distance;}
    Float_t             GetCPV_Thikness        (void) const {return fCPV_Thickness;}
    Int_t               GetNz                  (void) const {return fNz;}
    Int_t               GetNphi                (void) const {return fNphi;}
    Float_t             GetPhi                 (void) const {return fPhi;}

    void                Clear(Option_t *opt="");                            // Clear all data.
    void                Print(Option_t *opt="");
    void                Distortion(const TH2F *Noise=NULL, const TH2F *Stochastic=NULL, const TH2F *Calibration=NULL);
    TH2F               *CreateHistForDistortion(const char *name, const char *title, Int_t Nx, Int_t Ny, 
						 Float_t MU_mu, Float_t MU_sigma, Float_t SIGMA_mu, Float_t SIGMA_sigma);
    Float_t             GetDistortedValue(const TH2F *h, UInt_t n);

    void                Reconstruction(Float_t signal_step, UInt_t min_signal_reject);

    void                GetXY(const Hep3Vector &p,const Hep3Vector &v,float R,float &x,float &y,float &l) const;

    TObjArray          &GetGammasReconstructed (void)       {return fGammasReconstructed;}
    TObjArray          &GetParticles           (void)       {return fParticles;}

    TH2F                fCellEnergy;                            // GeV. Energy in cells
    TH2S                fChargedTracksInPIN;                    // amount. hits in PIN

    TArrayF             fCPV_hitsX;                             // cm. X-hits in CPV detector. (0 - cradle center)
    TArrayF             fCPV_hitsY;                             // cm. Y-hits in CPV detector. (0 - cradle center)

  private:

    Int_t               fGeometry;                              // Geometry type: 1 or 2
    Float_t             fCrystalSideSize;                       // cm.
    Float_t             fCrystalLength;                         // cm.
    Float_t             fWrapThickness;                         // cm.
    Float_t             fAirThickness;                          // cm.
    Float_t             fPIN_SideSize;                          // cm.
    Float_t             fPIN_Length;                            // cm.

    Float_t             fRadius;                                // cm. Distance to PHOS

    Float_t             fCPV_PHOS_Distance;                     // cm. Distance from CPV to PHOS.
    Float_t             fCPV_Thickness;                         // cm. CPV thikness.

    Int_t               fNz;                                    // Cells amount in beam direction
    Int_t               fNphi;                                  // Cells amount around beam

    Float_t             fPhi;                                   // degree. Position of CRADLE center

    TObjArray           fGammasReconstructed;                   // List of reconstructed gammas
    TObjArray           fParticles;                             // List of particles in the direction of this cradle

    TRandom             r;                                      //! Random number class, do not stream

//  friend class AliPHOS;

    ClassDef(AliPHOSCradle,1)   // PHOS cradle
};

class AliPHOS : public AliDetector {

 public:

    enum {CradlesBranch_Bit=1};

                        AliPHOS();
                        AliPHOS(const char *name, const char *title);
  virtual              ~AliPHOS();
  virtual void          AddHit(Int_t, Int_t*, Float_t*);
  virtual void          BuildGeometry();
  virtual void          CreateGeometry() {}
  virtual void          CreateMaterials();
  Int_t                 DistancetoPrimitive(Int_t px, Int_t py);
  void                  FinishEvent(void);

  virtual void          Init();
  virtual Int_t         IsVersion() const =0;
  void                  MakeBranch(Option_t *option);
  void                  SetTreeAddress(void);
  void                  FinishRun(void);
  void                  ResetDigits(void);
  void                  Print(Option_t *opt="");
  AliPHOSCradle        *GetCradleOfTheParticle(const Hep3Vector &p,const Hep3Vector &v) const;
  AliPHOSCradle        &GetCradle(int n) {return *(AliPHOSCradle*)fCradles->operator[](n);}
  //  AliPHOSCradle        &GetCradle(int n) {return *((AliPHOSCradle*) (*fCradles)[n]) ;}
  void                  Reconstruction(Float_t signal_step, UInt_t min_signal_reject);
  virtual void          SetFlags(Float_t p1,Float_t p2=0,Float_t p3=0,Float_t p4=0,
                               Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  virtual void          SetCell(Float_t p1,Float_t p2=0,Float_t p3=0,Float_t p4=0,
                               Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  virtual void          SetRadius(Float_t radius);
  virtual void          SetCradleSize(Int_t nz, Int_t nphi, Int_t ncradles);
  virtual void          SetCradleA(Float_t angle);
  virtual void          SetCPV(Float_t p1,Float_t p2=0,Float_t p3=0,Float_t p4=0,
                               Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  virtual void          SetExtra(Float_t p1,Float_t p2=0,Float_t p3=0,Float_t p4=0,
                               Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  virtual void          SetTextolitWall(Float_t dx, Float_t dy, Float_t dz);
  virtual void          SetInnerAir(Float_t dx, Float_t dy, Float_t dz);
  virtual void          SetFoam(Float_t dx, Float_t dy, Float_t dz, Float_t dr);
  virtual void          StepManager()=0;
  virtual void          DefPars();
  virtual void          AddPHOSCradles();



  virtual Int_t GetPHOS_IDTMED_PbWO4         (void){return gAlice->Idtmed()[700-1];}
  virtual Int_t GetPHOS_IDTMED_CPV           (void){return gAlice->Idtmed()[701-1];}
  virtual Int_t GetPHOS_IDTMED_Al            (void){return gAlice->Idtmed()[702-1];}
  virtual Int_t GetPHOS_IDTMED_Tyvek         (void){return gAlice->Idtmed()[703-1];}
  virtual Int_t GetPHOS_IDTMED_PIN           (void){return gAlice->Idtmed()[706-1];}
  virtual Int_t GetPHOS_IDTMED_AIR           (void){return gAlice->Idtmed()[799-1];}
  
  
  virtual Int_t   &GetPHOS_Ndiv_magic    (void) {return PHOS_Ndiv_magic;}
  virtual Float_t  GetCrystalSideSize    (void) const {return PHOScell[0]; }
  virtual Float_t  GetCrystalLength      (void) const {return PHOScell[1]; }
  virtual Float_t  GetWrapThickness      (void) const {return PHOScell[2]; }
  virtual Float_t  GetAirThickness       (void) const {return PHOScell[3]; }
  virtual Float_t  GetPIN_SideSize       (void) const {return PHOScell[4]; }
  virtual Float_t  GetPIN_Length         (void) const {return PHOScell[5]; }
  virtual Float_t  GetRadius             (void) const {return PHOSradius;  }
  virtual Float_t  GetCPV_Thickness      (void) const {return PHOSCPV[0];  }
  virtual Float_t  GetCPV_PHOS_Distance  (void) const {return PHOSCPV[1];  }
  virtual Int_t    GetNz                 (void) const {return PHOSsize[0]; }
  virtual Int_t    GetNphi               (void) const {return PHOSsize[1]; }
  virtual Int_t    GetCradlesAmount      (void) const {return PHOSsize[2]; }
  virtual Float_t  GetAngleBetweenCradles(void) const {return PHOScradlesA;}
  virtual Float_t  GetPHOS_flag          (Int_t n) const {return PHOSflags[n];}
  virtual Float_t  GetPHOSextra          (Int_t n) const {return PHOSextra[n];}
  virtual Float_t  GetPHOSFoam           (Int_t n) const {return PHOSFTI[n];}
  virtual Float_t  GetPHOStxwall         (Int_t n) const {return PHOSTXW[n];}
  virtual Float_t  GetPHOSAir            (Int_t n) const {return PHOSAIR[n];}
  virtual Float_t &GetCradleAngle        (Int_t n) {return PHOSangle[n];}


  TObjArray            *fCradles;               //!  Cradles in PHOS
  Int_t                 fDebugLevel;

  TTree                *fTreePHOS;              //! Pointer to PHOS tree.

private:

  TString               fBranchNameOfCradles;   // 
  TString               fTreeName;              // Name of PHOS tree: "PHOS"

#define MAXCRAD 100

  Float_t PHOSflags[9], PHOScell[9], PHOSradius, PHOSCPV[9];
  Int_t   PHOSsize[3];
  Float_t PHOScradlesA,PHOSTXW[3],PHOSAIR[3],PHOSFTI[4],PHOSextra[9],
        PHOSangle[MAXCRAD];
  Int_t   PHOS_Ndiv_magic;

 ClassDef(AliPHOS,1)  //Hits manager for set:PHOS
};
 
#endif

