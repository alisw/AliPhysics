#ifndef ALI_TPC_CORRECTION_H
#define ALI_TPC_CORRECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////////////
// AliTPCCorrection class                                                     //
////////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>
#include "TMatrixD.h"
#include "TMatrixF.h"
class TH2F;
class TTimeStamp;
class TCollection;
class TTreeSRedirector;
class AliExternalTrackParam;
class TTree;
class THnSparse;
class AliESDVertex;


class AliTPCCorrection : public TNamed {
public:
  enum CompositionType {kParallel,kQueue};

  AliTPCCorrection();
  AliTPCCorrection(const char *name,const char *title);
  virtual ~AliTPCCorrection();
  

  // functions to correct a space point
          void CorrectPoint (      Float_t x[],const Short_t roc);
          void CorrectPointLocal(Float_t x[],const Short_t roc);
          void CorrectPoint (const Float_t x[],const Short_t roc,Float_t xp[]);
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

  // functions to distort a space point
          void DistortPoint (      Float_t x[],const Short_t roc);
          void DistortPointLocal(Float_t x[],const Short_t roc);
          void DistortPoint (const Float_t x[],const Short_t roc,Float_t xp[]);
  virtual void GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]);

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);

  // convenience functions
  virtual void Print(Option_t* option="") const;
 
  TH2F* CreateHistoDRinXY   (Float_t z=10.,Int_t nx=100,Int_t ny=100);
  TH2F* CreateHistoDRPhiinXY(Float_t z=10.,Int_t nx=100,Int_t nphi=100);
  TH2F* CreateHistoDZinXY   (Float_t z=10.,Int_t nx=100,Int_t ny=100);

  TH2F* CreateHistoDRinZR   (Float_t phi=0.,Int_t nZ=100,Int_t nR=100);
  TH2F* CreateHistoDRPhiinZR(Float_t phi=0.,Int_t nZ=100,Int_t nR=100);
  TH2F* CreateHistoDZinZR   (Float_t phi=0.,Int_t nZ=100,Int_t nR=100);


  TTree* CreateDistortionTree(Double_t step=5);
  static void  MakeDistortionMap(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run,  Float_t refX, Int_t type, Int_t integ=1);
  static void  MakeDistortionMapCosmic(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run,  Float_t refX, Int_t type);
  static void  MakeDistortionMapSector(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run, Int_t type);
  // normally called directly in the correction classes which inherit from this class
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2);
  AliExternalTrackParam * FitDistortedTrack(AliExternalTrackParam & trackIn, Double_t refX, Int_t dir,TTreeSRedirector *pcstream);
  void StoreInOCDB(Int_t startRun, Int_t endRun, const char *comment=0);
  static void MakeTrackDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step=1, Int_t offset=0, Bool_t debug=0);
  static void MakeSectorDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step=1, Int_t offset=0, Bool_t debug=0);
  static void MakeLaserDistortionTreeOld(TTree* tree, TObjArray *corrArray, Int_t itype);
  static void MakeLaserDistortionTree(TTree* tree, TObjArray *corrArray, Int_t itype);

  void FastSimDistortedVertex(Double_t orgVertex[3], Int_t nTracks, AliESDVertex &aV, AliESDVertex &avOrg, AliESDVertex &cV, AliESDVertex &cvOrg, TTreeSRedirector * const pcstream, Double_t etaCuts);

  static void AddVisualCorrection(AliTPCCorrection* corr, Int_t position);
  static Double_t GetCorrSector(Double_t sector, Double_t r, Double_t kZ, Int_t axisType, Int_t corrType=0);
  static Double_t GetCorrXYZ(Double_t gx, Double_t gy, Double_t gz, Int_t axisType, Int_t corrType=0);

protected:
  TH2F* CreateTH2F(const char *name,const char *title,
		   const char *xlabel,const char *ylabel,const char *zlabel,
		   Int_t nbinsx,Double_t xlow,Double_t xup,
		   Int_t nbinsy,Double_t ylow,Double_t yup);
 
  static const Double_t fgkTPCZ0;       // nominal gating grid position 
  static const Double_t fgkIFCRadius;   // Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static const Double_t fgkOFCRadius;   // Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static const Double_t fgkZOffSet;     // Offset from CE: calculate all distortions closer to CE as if at this point
  static const Double_t fgkCathodeV;    // Cathode Voltage (volts)
  static const Double_t fgkGG;          // Gating Grid voltage (volts)
  static const Double_t fgkdvdE;        // [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)
  static const Double_t fgkEM;          // charge/mass in [C/kg]
  static const Double_t fgke0;          // vacuum permittivity [A·s/(V·m)]


  enum {kNR=   72  };     // Number of R points in the table for interpolating distortion data
  enum {kNPhi= 18*10+1};  // Number of Phi points in the table for interpolating distortion data ( plus one extra for 360 == 0 ) 
  enum {kNZ=   166};      // Number of Z points in the table for interpolating distortion data

  Double_t fgkRList[kNR]; // points in the radial direction (for the lookup table)
  Double_t fgkPhiList[kNPhi]; // points in the phi direction (for the lookup table)
  Double_t fgkZList[kNZ]; // points in the z direction (for the lookup table)

  // Simple Interpolation functions: e.g. with tricubic interpolation (not yet in TH3)
  Int_t fILow, fJLow, fKLow;          // variable to help in the interpolation 
  // Double_t versions
  void Interpolate2DEdistortion( const Int_t order, const Double_t r, const Double_t z, 
				 const Double_t er[kNZ][kNR], Double_t &erValue );
  void Interpolate3DEdistortion( const Int_t order, const Double_t r, const Float_t phi, const Double_t z, 
				 const Double_t er[kNZ][kNPhi][kNR], const Double_t ephi[kNZ][kNPhi][kNR], 
				 const Double_t ez[kNZ][kNPhi][kNR],
				 Double_t &erValue, Double_t &ephiValue, Double_t &ezValue);
  // TMatrixD versions (for e.g. Poisson relaxation)
  Double_t Interpolate2DTable( const Int_t order, const Double_t x, const Double_t y, 
			      const Int_t nx,  const Int_t ny, const Double_t xv[], const Double_t yv[], 
			      const TMatrixD &array );
  Double_t Interpolate3DTable( const Int_t order, const Double_t x,   const Double_t y,   const Double_t z,
			      const Int_t  nx,    const Int_t  ny,    const Int_t  nz,
			      const Double_t xv[], const Double_t yv[], const Double_t zv[],
			      TMatrixD **arrayofArrays );
  Double_t Interpolate( const Double_t xArray[], const Double_t yArray[], 
			const Int_t order, const Double_t x );
  void Search( const Int_t n, const Double_t xArray[], const Double_t x, Int_t &low );
 
  // TMatrixF versions (smaller size, e.g. for final look up table)
  Float_t Interpolate2DTable( const Int_t order, const Double_t x, const Double_t y, 
			      const Int_t nx,  const Int_t ny, const Double_t xv[], const Double_t yv[], 
			      const TMatrixF &array );
  Float_t Interpolate3DTable( const Int_t order, const Double_t x,   const Double_t y,   const Double_t z,
			      const Int_t  nx,    const Int_t  ny,    const Int_t  nz,
			      const Double_t xv[], const Double_t yv[], const Double_t zv[],
			       TMatrixF **arrayofArrays ); 
  Float_t Interpolate( const Double_t xArray[], const Float_t yArray[], 
			const Int_t order, const Double_t x );

  virtual Int_t IsPowerOfTwo ( Int_t i ) const  ;

  
 
  // Algorithms to solve the laplace or poisson equation 
  void PoissonRelaxation2D(TMatrixD &arrayV, TMatrixD &chargeDensity, 
			   TMatrixD &arrayErOverEz, TMatrixD &arrayDeltaEz,
			   const Int_t rows, const Int_t columns, const Int_t iterations,
			   const Bool_t rocDisplacement = kTRUE);

  void PoissonRelaxation3D( TMatrixD **arrayofArrayV, TMatrixD **arrayofChargeDensities, 
			    TMatrixD **arrayofEroverEz, TMatrixD **arrayofEPhioverEz, TMatrixD **arrayofEz,
			    const Int_t rows, const Int_t columns,  const Int_t phislices, 
			    const Float_t deltaphi, const Int_t iterations, const Int_t summetry,
			    const Bool_t rocDisplacement = kTRUE); 
    
protected:
  Double_t fT1;         // tensor term of wt - T1
  Double_t fT2;         // tensor term of wt - T2
  static TObjArray *fgVisualCorrection;  // array of orrection for visualization
private:
  AliTPCCorrection(const AliTPCCorrection &);               // not implemented
  AliTPCCorrection &operator=(const AliTPCCorrection &);    // not implemented

  void InitLookUpfulcrums();   // to initialize the grid of the look up table

  ClassDef(AliTPCCorrection,4);
};

#endif
