#ifndef ALIITSALIGNMILLE_H
#define ALIITSALIGNMILLE_H
/* Copyright(c) 2007-2009 , ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


/// \ingroup rec
/// \class AliITSAlignMille
/// \brief Class for alignment of ITS
//
// Authors: Marcello Lunardon

#include <TString.h>
#include <TObject.h>
#include "AliTrackPointArray.h"

class AliMillepede;
class AliAlignObjParams;
class TGeoManager;
class TGeoHMatrix;

// number of used objects
#define ITSMILLE_NDETELEM    2198
#define ITSMILLE_NPARCH         6
#define ITSMILLE_NLOCAL         4
#define ITSMILLE_NSTDEV         3       

class AliITSAlignMille:public TObject
{

public:
  AliITSAlignMille(const Char_t *configFilename="AliITSAlignMille.conf", Bool_t initmille=kTRUE);
  virtual ~AliITSAlignMille();
  
  // geometry methods 
  Int_t     GetModuleIndex(const Char_t *symname);
  Int_t     GetModuleIndex(UShort_t voluid);
  UShort_t  GetModuleVolumeID(const Char_t *symname);
  UShort_t  GetModuleVolumeID(Int_t index);
  void      SetCurrentModule(Int_t index);

  // configuration methods
  void      SetGeometryFileName(const Char_t* filename="geometry.root") 
    { fGeometryFileName = filename; }
  const Char_t* GetGeometryFileName() {return fGeometryFileName.Data();}
  void      PrintCurrentModuleInfo();
  void      Print();
  
  // fitting methods
  void      SetMinNPtsPerTrack(Int_t pts=3) {fMinNPtsPerTrack=pts;}
  //Bool_t    CheckTrack(AliTrackPointArray *track);
  Int_t     ProcessTrack(AliTrackPointArray *track);
  void      InitTrackParams(int meth=1);
  Int_t     InitModuleParams();
  Int_t     CheckCurrentTrack();
  Bool_t    CheckVolumeID(UShort_t voluid) const ;
  Int_t     CalcIntersectionPoint(Double_t *lpar, Double_t *gpar);
  Int_t     CalcDerivatives(Int_t paridx, Bool_t islpar);
  Double_t* GetLocalIntersectionPoint() {return fPintLoc;}
  Double_t* GetGlobalIntersectionPoint() {return fPintGlo;}

  // millepede methods
  void      FixParameter(Int_t param, Double_t value);
  void      AddConstraint(Double_t *factor, Double_t value );
  void      InitGlobalParameters(Double_t *par);   
  void      SetLocalDerivative(Int_t index, Double_t value) 
    {fLocalDerivatives[index] = value;}
  void      SetGlobalDerivative(Int_t index, Double_t value) 
    {fGlobalDerivatives[index] = value;}  
  void      LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t lSingleFit);
  void      GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls);
  void      PrintGlobalParameters();
  Double_t  GetParError(Int_t iPar);
  Int_t     SetLocalEquations();
  
  // fitting stuffs
  AliTrackPointArray *GetCurrentTrack() {return fTrack;}
  AliTrackPoint      *GetCurrentCluster() {return &fCluster;}

  // geometry stuffs
  Int_t  GetNModules() const {return fNModules;}
  Int_t  GetCurrentModuleIndex() const {return fCurrentModuleIndex;}
  TGeoHMatrix *GetCurrentModuleHMatrix() {return fCurrentModuleHMatrix;}
  Double_t    *GetCurrentModuleTranslation() {return fCurrentModuleTranslation;}
  Int_t  GetCurrentModuleInternalIndex() const {return fCurrentModuleInternalIndex;}
  Int_t       *GetModuleIndexArray() {return fModuleIndex;}
  UShort_t    *GetModuleVolumeIDArray() {return fModuleVolumeID;}
  
 private:

  // configuration methods
  Int_t     LoadConfig(const Char_t *cfile="AliITSAlignMille.conf");
  void      ResetLocalEquation();
  void      InitGeometry();

  // millepede methods
  void      Init(Int_t nGlobal, Int_t nLocal, Int_t nStdDev);

  // millepede stuffs
  AliMillepede *fMillepede;   ///< Detector independent alignment class
  static Int_t  fgNParCh;      ///< Number of degrees of freedom per chamber
  static Int_t  fgNDetElem;    ///< Total number of detection elements
  Double_t      fStartFac;      ///< Initial value for chi2 cut 
                              ///< if > 1 Iterations in AliMil. are turned on
  Double_t      fResCutInitial; ///< Cut on residual for first iteration
  Double_t      fResCut;        ///< Cut on residual for other iterations 
  Int_t         fNGlobal;       ///< Number of global parameters
  Int_t         fNLocal;        ///< Number of local parameters
  Int_t         fNStdDev;       ///< Number of standard deviations for chi2 cut
  Bool_t        fIsMilleInit;  ///
  Double_t      fParSigTranslations; ///< init sigma for transl. params [cm]
  Double_t      fParSigRotations; ///< init sigma for rot. params [deg]

  // fitting stuffs
  AliTrackPointArray *fTrack;       ///< pointer to current track 
  AliTrackPoint fCluster;           ///< current cluster
  Double_t     *fGlobalDerivatives;   ///< Array of global derivatives
  Double_t      fLocalDerivatives[ITSMILLE_NLOCAL]; ///< Array of local deriv.
  Double_t      fLocalInitParam[ITSMILLE_NLOCAL];   ///< Array with inital values for local parameters for current track
  Double_t      fModuleInitParam[ITSMILLE_NPARCH];  ///< Array with inital values for current module parameters (init geometry)
  Double_t      fPintLoc[3]; ///
  Double_t      fPintLoc0[3]; ///
  Double_t      fPintGlo[3]; ///
  Double_t      fMeasLoc[3]; // current point local coordinates (the original ones)
  Double_t      fMeasGlo[3]; // current point glob. coord (AliTrackPoint)
  Double_t      fSigmaLoc[3]; // stdev current point
  TGeoHMatrix  *fTempHMat; ///
  AliAlignObjParams *fTempAlignObj; ///
  Double_t      fDerivativeXLoc; // localX deriv.
  Double_t      fDerivativeZLoc; // localZ deriv.
  Double_t      fDeltaPar; ///
  Int_t         fMinNPtsPerTrack; ///
  
  // geometry stuffs
  TString       fGeometryFileName;  ///
  TGeoManager  *fGeoManager;        ///
  Int_t         fCurrentModuleIndex;   ///
  Int_t         fCurrentModuleInternalIndex;  ///
  Double_t      fCurrentModuleTranslation[3]; ///
  Int_t         fNModules;  /// number of defined modules from config file
  Int_t         fModuleIndex[ITSMILLE_NDETELEM]; ///
  UShort_t      fModuleVolumeID[ITSMILLE_NDETELEM];  ///
  Bool_t        fFreeParam[ITSMILLE_NDETELEM][ITSMILLE_NPARCH];  ///
  Bool_t        fUseLocalShifts; /// 
  TGeoHMatrix  *fCurrentModuleHMatrix; /// 

  AliITSAlignMille(const AliITSAlignMille& rhs);
  AliITSAlignMille& operator=(const AliITSAlignMille& rhs);


ClassDef(AliITSAlignMille, 0)};

#endif
