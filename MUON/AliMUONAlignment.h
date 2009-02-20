#ifndef ALIMUONALIGNMENT_H
#define ALIMUONALIGNMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONAlignment
/// \brief Class for alignment of muon spectrometer
//
// Authors: Bruce Becker, Javier Castillo

#include <TObject.h>

class TGeoCombiTrans;
class TClonesArray;
class AliMillepede;
class AliMUONGeometryTransformer;
class AliMUONTrack;
class AliMUONTrackParam;
class AliMUONVCluster;

class AliMUONAlignment:public TObject
{

public:
  AliMUONAlignment();
  virtual ~AliMUONAlignment();

  void ProcessTrack(AliMUONTrack *track);
  /// Set geometry transformer
  void SetGeometryTransformer(AliMUONGeometryTransformer * transformer) {
    fTransform = transformer;
  }

  /// Set flag for Magnetic field On/Off
  void SetBFieldOn(Bool_t bBFieldOn) {
    fBFieldOn =  bBFieldOn;
  }
  /// Define chambers to align
  void SetChOnOff(Bool_t *bChOnOff) {
    for(int iCh=0; iCh<10; iCh++)
      fChOnOff[iCh] =  bChOnOff[iCh];
  }
  /// Possibility to align only one side of the detector
  void SetSpecLROnOff(Bool_t *bSpecLROnOff) {
    fSpecLROnOff[0] =  bSpecLROnOff[0];    
    fSpecLROnOff[1] =  bSpecLROnOff[1];    
  }
  void FixStation(Int_t iSt);
  void FixChamber(Int_t iCh);
  void FixHalfSpectrometer(Bool_t *bChOnOff, Bool_t *bSpecLROnOff);
  void AllowVariations(Bool_t *bChOnOff);
  void SetNonLinear(Bool_t *bChOnOff, Bool_t *bVarXYT);
  void AddConstraints(Bool_t *bChOnOff, Bool_t *bVarXYT);
  void AddConstraints(Bool_t *bChOnOff, Bool_t *bVarXYT, Bool_t *bDetTLBR, Bool_t *bSpecLROnOff);
  void ResetConstraints();

  void FixParameter(Int_t param, Double_t value);
  void SetNonLinear(Int_t param);
  void AddConstraint(Double_t *factor, Double_t value );
  void InitGlobalParameters(Double_t *par);   
  /// Set array of local derivatives
  void SetLocalDerivative(Int_t index, Double_t value) {      
    fLocalDerivatives[index] = value;
  }
  /// Set array of global derivatives
  void SetGlobalDerivative(Int_t index, Double_t value) {
    fGlobalDerivatives[index] = value;
  }  
  void LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t lSingleFit);
  void GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls);
  void PrintGlobalParameters();
  Double_t GetParError(Int_t iPar);
  
  AliMUONGeometryTransformer* 
    ReAlign(const AliMUONGeometryTransformer * transformer, double *misAlignments, Bool_t verbose);

  void SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t chId, Double_t chResX, Double_t chResY, Double_t deResX, Double_t deResY);

 private:
  /// Not implemented
  AliMUONAlignment(const AliMUONAlignment& right);
  /// Not implemented
  AliMUONAlignment&  operator = (const AliMUONAlignment& right);

  void Init(Int_t nGlobal, Int_t nLocal, Int_t nStdDev);
  void ConstrainT(Int_t lDetElem, Int_t lCh, Double_t *lConstraintT, Int_t iVar, Double_t lWeight=1.0);
  void ConstrainL(Int_t lDetElem, Int_t lCh, Double_t *lConstraintL, Int_t iVar, Double_t lWeight=1.0);
  void ConstrainB(Int_t lDetElem, Int_t lCh, Double_t *lConstraintB, Int_t iVar, Double_t lWeight=1.0);
  void ConstrainR(Int_t lDetElem, Int_t lCh, Double_t *lConstraintR, Int_t iVar, Double_t lWeight=1.0);
  void FillDetElemData();
  void FillRecPointData();
  void FillTrackParamData();
  void ResetLocalEquation();
  void LocalEquationX();
  void LocalEquationY();

  TGeoCombiTrans ReAlign(const TGeoCombiTrans& transform, double *detElemMisAlignment) const;

  Bool_t fBFieldOn;        ///< Flag for Magnetic filed On/Off
  Bool_t fChOnOff[10];     ///< Flags for chamber On/Off
  Bool_t fSpecLROnOff[2];  ///< Flags for left right On/Off		       	 		       	       		       
  Bool_t fDoF[3];          ///< Flags degrees of freedom to align (x,y,phi)
  Double_t fAllowVar[3];   ///< "Encouraged" variation for degrees of freedom 
  Double_t fStartFac;      ///< Initial value for chi2 cut 
                           ///< if > 1 Iterations in AliMillepede are turned on
  Double_t fResCutInitial; ///< Cut on residual for first iteration
  Double_t fResCut;        ///< Cut on residual for other iterations 

  AliMillepede *fMillepede; ///< Detector independent alignment class
  
  TClonesArray *fTrackParamAtCluster; ///< Array of track parameters 
  AliMUONTrack *fTrack;               ///< AliMUONTrack 
  AliMUONVCluster *fCluster;          ///< AliMUONVCluster
  AliMUONTrackParam *fTrackParam;     ///< Track parameters 

  Int_t fNGlobal;  ///< Number of global parameters
  Int_t fNLocal;   ///< Number of local parameters
  Int_t fNStdDev;  ///< Number of standard deviations for chi2 cut
  Double_t fClustPos[3];    ///< Cluster position
  Double_t fClustPosLoc[3]; ///< Cluster position in local coordinates
  Double_t fTrackSlope0[2]; ///< Track slope at reference point
  Double_t fTrackSlope[2];  ///< Track slope at current point
  Double_t fTrackPos0[3];   ///< Track intersection at reference point
  Double_t fTrackPos[3];    ///< Track intersection at current point
  Double_t fTrackPosLoc[3]; ///< Track intersection at current point in local coordinates 
  Double_t fMeas[2];        ///< Current measurement (depend on B field On/Off)  
  Double_t fSigma[2];       ///< Estimated resolution on measurement

  Double_t fGlobalDerivatives[468]; ///< Array of global derivatives
  Double_t fLocalDerivatives[4];    ///< Array of local derivatives

  Double_t fConstraintX[468];   ///< Array for constraint equation all X
  Double_t fConstraintY[468];   ///< Array for constraint equation all Y
  Double_t fConstraintP[468];   ///< Array for constraint equation all P
  Double_t fConstraintXT[468];  ///< Array for constraint equation X Top half
  Double_t fConstraintYT[468];  ///< Array for constraint equation Y Top half
  Double_t fConstraintPT[468];  ///< Array for constraint equation P Top half
  Double_t fConstraintXZT[468];  ///< Array for constraint equation X vs Z Top half
  Double_t fConstraintYZT[468];  ///< Array for constraint equation Y vs Z Top half
  Double_t fConstraintPZT[468];  ///< Array for constraint equation P vs Z Top half
  Double_t fConstraintXYT[468];  ///< Array for constraint equation X vs Y Top half
  Double_t fConstraintYYT[468];  ///< Array for constraint equation Y vs Y Top half
  Double_t fConstraintPYT[468];  ///< Array for constraint equation P vs Y Top half
  Double_t fConstraintXB[468];  ///< Array for constraint equation X Bottom half
  Double_t fConstraintYB[468];  ///< Array for constraint equation Y Bottom half
  Double_t fConstraintPB[468];  ///< Array for constraint equation P Bottom half
  Double_t fConstraintXZB[468];  ///< Array for constraint equation X vs Z Bottom half
  Double_t fConstraintYZB[468];  ///< Array for constraint equation Y vs Z Bottom half
  Double_t fConstraintPZB[468];  ///< Array for constraint equation P vs Z Bottom half
  Double_t fConstraintXYB[468];  ///< Array for constraint equation X vs Y Bottom half
  Double_t fConstraintYYB[468];  ///< Array for constraint equation Y vs Y Bottom half
  Double_t fConstraintPYB[468];  ///< Array for constraint equation P vs Y Bottom half
  Double_t fConstraintXR[468];  ///< Array for constraint equation X Right half
  Double_t fConstraintYR[468];  ///< Array for constraint equation Y Right half
  Double_t fConstraintPR[468];  ///< Array for constraint equation P Right half
  Double_t fConstraintXZR[468];  ///< Array for constraint equation X vs Z Right half
  Double_t fConstraintYZR[468];  ///< Array for constraint equation Y vs Z Right half
  Double_t fConstraintPZR[468];  ///< Array for constraint equation P vs Z Right half
  Double_t fConstraintXYR[468];  ///< Array for constraint equation X vs Y Right half
  Double_t fConstraintYYR[468];  ///< Array for constraint equation Y vs Y Right half
  Double_t fConstraintPYR[468];  ///< Array for constraint equation P vs Y Right half
  Double_t fConstraintXL[468];  ///< Array for constraint equation X Left half
  Double_t fConstraintYL[468];  ///< Array for constraint equation Y Left half
  Double_t fConstraintPL[468];  ///< Array for constraint equation P Left half
  Double_t fConstraintXZL[468];  ///< Array for constraint equation X vs Z Left half
  Double_t fConstraintYZL[468];  ///< Array for constraint equation Y vs Z Left half
  Double_t fConstraintPZL[468];  ///< Array for constraint equation P vs Z Left half
  Double_t fConstraintXYL[468];  ///< Array for constraint equation X vs Y Left half
  Double_t fConstraintYYL[468];  ///< Array for constraint equation Y vs Y Left half
  Double_t fConstraintPYL[468];  ///< Array for constraint equation P vs Y Left half
  Double_t fConstraintX3[468];  ///< Array for constraint equation St3 X
  Double_t fConstraintY3[468];  ///< Array for constraint equation St3 Y
  Double_t fConstraintX4[468];  ///< Array for constraint equation St4 X
  Double_t fConstraintY4[468];  ///< Array for constraint equation St4 Y
  Double_t fConstraintP4[468];  ///< Array for constraint equation St4 P
  Double_t fConstraintX5[468];  ///< Array for constraint equation St5 X
  Double_t fConstraintY5[468];  ///< Array for constraint equation St5 Y

  Int_t fDetElemId;        ///< Detection element id
  Int_t fDetElemNumber;    ///< Detection element number
  Double_t fPhi;           ///< Azimuthal tilt of detection element 
  Double_t fCosPhi;        ///< Cosine of fPhi
  Double_t fSinPhi;        ///< Sine of fPhi
  Double_t fDetElemPos[3]; ///< Position of detection element

  AliMUONGeometryTransformer *fTransform; ///< Geometry transformation

  static Int_t fgNSt;            ///< Number tracking stations
  static Int_t fgNCh;            ///< Number tracking chambers
  static Int_t fgNTrkMod;        ///< Number of tracking modules (4 ch + 6*2 half-ch)
  static Int_t fgNParCh;         ///< Number of degrees of freedom per chamber
  static Int_t fgNDetElem;       ///< Total number of detection elements
  static Int_t fgNDetElemCh[10]; ///< Number of detection elements per chamber
  static Int_t fgSNDetElemCh[10];///< Sum of detection elements up to this chamber (inc)

ClassDef(AliMUONAlignment, 0) //Class for alignment of muon spectrometer
};

#endif
