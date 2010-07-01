#ifndef ALI_TPC_CORRECTION_H
#define ALI_TPC_CORRECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliTPCCorrection class                                                     //
//                                                                            //
// This class provides a general framework to deal with space point           //
// distortions. An correction class which inherits from here is for example   //
// AliTPCExBBShape or AliTPCExBTwist                                          //
//                                                                            //
// General functions are (for example):                                       //
//   CorrectPoint(x,roc) where x is the vector of inital positions in         //
//   cartesian coordinates and roc represents the Read Out chamber number     //
//   according to the offline naming convention. The vector x is overwritten  //
//   with the corrected coordinates.                                          //
//                                                                            //
// An alternative usage would be CorrectPoint(x,roc,dx), which leaves the     //
//   vector x untouched, put returns the distortions via the vector dx        //
//                                                                            //
// The class allows "effective Omega Tau" corrections to be shifted to the    //
// single distortion classes.                                                 //
//                                                                            //
// Note: This class is normally used via the class AliTPCComposedCorrection   //
//                                                                            //
// date: 27/04/2010                                                           //
// Authors: Magnus Mager, Stefan Rossegger, Jim Thomas                        //
////////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>
#include "TMatrixD.h"
class TH2F;
class TTimeStamp;
class TCollection;
class TTreeSRedirector;
class AliExternalTrackParam;
class TTree;
class THnSparse;



class AliTPCCorrection : public TNamed {
public:
  enum CompositionType {kParallel,kQueue};

  AliTPCCorrection();
  AliTPCCorrection(const char *name,const char *title);
  virtual ~AliTPCCorrection();
  

  // functions to correct a space point
          void CorrectPoint (      Float_t x[],const Short_t roc);
          void CorrectPoint (const Float_t x[],const Short_t roc,Float_t xp[]);
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

  // functions to distort a space point
          void DistortPoint (      Float_t x[],const Short_t roc);
          void DistortPoint (const Float_t x[],const Short_t roc,Float_t xp[]);
  virtual void GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]);

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);

  // convenience functions
  virtual void Print(Option_t* option="") const;
 
  TH2F* CreateHistoDRinXY   (Float_t z=10.,Int_t nx=100,Int_t ny=100);
  TH2F* CreateHistoDRPhiinXY(Float_t z=10.,Int_t nx=100,Int_t nphi=100);
  TH2F* CreateHistoDRinZR   (Float_t phi=0.,Int_t nZ=100,Int_t nR=100);
  TH2F* CreateHistoDRPhiinZR(Float_t phi=0.,Int_t nZ=100,Int_t nR=100);
  TTree* CreateDistortionTree(Double_t step=5);
  static void  MakeDistortionMap(THnSparse * his0, TTreeSRedirector *pcstream, const char* hname, Int_t run);
  // normally called directly in the correction classes which inherit from this class
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2);
  AliExternalTrackParam * FitDistortedTrack(AliExternalTrackParam & trackIn, Double_t refX, Int_t dir,TTreeSRedirector *pcstream);
  void StoreInOCDB(Int_t startRun, Int_t endRun, const char *comment=0);
  static void MakeTrackDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step=1, Bool_t debug=0);
  static void MakeLaserDistortionTree(TTree* tree, TObjArray *corrArray, Int_t itype);
protected:
  TH2F* CreateTH2F(const char *name,const char *title,
		   const char *xlabel,const char *ylabel,const char *zlabel,
		   Int_t nbinsx,Double_t xlow,Double_t xup,
		   Int_t nbinsy,Double_t ylow,Double_t yup);
 
  static const Double_t fgkTPCZ0;      // nominal gating grid position 
  static const Double_t fgkIFCRadius;   // Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static const Double_t fgkOFCRadius;   // Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static const Double_t fgkZOffSet;     // Offset from CE: calculate all distortions closer to CE as if at this point
  static const Double_t fgkCathodeV;    // Cathode Voltage (volts)
  static const Double_t fgkGG;          // Gating Grid voltage (volts)

  enum {kNR=   92};              // Number of R points in the table for interpolating distortion data
  enum {kNZ=  270};              // Number of Z points in the table for interpolating distortion data
  static const Double_t fgkRList[kNR]; // points in the radial direction (for the lookup table)
  static const Double_t fgkZList[kNZ]; // points in the z direction (for the lookup table)

  // Simple Interpolation functions: e.g. with tricubic interpolation (not yet in TH3)
  Int_t fJLow;         // variable to help in the interpolation 
  Int_t fKLow;         // variable to help in the interpolation 
  void Interpolate2DEdistortion( const Int_t order, const Double_t r, const Double_t z, 
				 const Double_t er[kNZ][kNR], Double_t &erValue );
  Double_t Interpolate( const Double_t xArray[], const Double_t yArray[], 
			const Int_t order, const Double_t x );
  void Search( const Int_t n, const Double_t xArray[], const Double_t x, Int_t &low );
  virtual Int_t IsPowerOfTwo ( Int_t i ) const  ;
    
  // Algorithms to solve the laplace or possion equation 
  void PoissonRelaxation2D(TMatrixD &arrayV, const TMatrixD &chargeDensity, TMatrixD &arrayErOverEz, const Int_t rows, const Int_t columns, const Int_t iterations );


protected:
  Double_t fT1;         // tensor term of wt - T1
  Double_t fT2;         // tensor term of wt - T2
private:
  ClassDef(AliTPCCorrection,2);
};

#endif
