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
#include <TArray.h>
#include "AliTrackPointArray.h"

class AliMillepede;
class AliAlignObjParams;
class TGeoManager;
class TGeoHMatrix;
class AliITSAlignMilleModule;
class AliTrackFitterRieman;

// number of used objects
#define ITSMILLE_NDETELEM    2198
#define ITSMILLE_NPARCH         6
#define ITSMILLE_NLOCAL         5
#define ITSMILLE_NSTDEV         3       


struct MilleData {
  /// structure to store data for 2 LocalEquations (X and Z)
  Double_t measX;
  Double_t sigmaX;
  Int_t    idxlocX[ITSMILLE_NLOCAL];
  Double_t derlocX[ITSMILLE_NLOCAL];
  Int_t    idxgloX[ITSMILLE_NPARCH];  
  Double_t dergloX[ITSMILLE_NPARCH];

  Double_t measZ;
  Double_t sigmaZ;
  Int_t    idxlocZ[ITSMILLE_NLOCAL];
  Double_t derlocZ[ITSMILLE_NLOCAL];
  Int_t    idxgloZ[ITSMILLE_NPARCH];  
  Double_t dergloZ[ITSMILLE_NPARCH];
};

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
  void      SetCurrentSensitiveModule(Int_t index); // set as current the SENSITIVE module with index 'index'

  // configuration methods
  void      SetGeometryFileName(const Char_t* filename="geometry.root") 
    { fGeometryFileName = filename; }
  const Char_t* GetGeometryFileName() {return fGeometryFileName.Data();}
  const Char_t* GetPreAlignmentFileName() {return fPreAlignmentFileName.Data();}
  void      PrintCurrentModuleInfo();
  void      Print(Option_t*) const;
  Bool_t    IsConfigured() const {return fIsConfigured;}
  void      SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts);
  
  // fitting methods
  void      SetMinNPtsPerTrack(Int_t pts=3) {fMinNPtsPerTrack=pts;}
  Int_t     ProcessTrack(AliTrackPointArray *track);
  AliTrackPointArray *PrepareTrack(AliTrackPointArray *track); // build a new AliTrackPointArray with selected conditions
  void      InitTrackParams(int meth=1);
  Bool_t    InitRiemanFit();
  AliTrackFitterRieman  *GetRiemanFitter() const {return fRieman;}
  Int_t     InitModuleParams();
  Int_t     CheckCurrentTrack();
  Bool_t    CheckVolumeID(UShort_t voluid) const; // checks voluid for sensitive volumes
  Int_t     IsDefined(UShort_t voluid) const;
  Int_t     IsContained(UShort_t voluid) const;
  Int_t     CalcIntersectionPoint(Double_t *lpar, Double_t *gpar);
  Int_t     CalcDerivatives(Int_t paridx, Bool_t islpar);
  Double_t* GetLocalIntersectionPoint() {return fPintLoc;}
  Double_t* GetGlobalIntersectionPoint() {return fPintGlo;}
  void      SetInitTrackParamsMeth(Int_t meth=1) {fInitTrackParamsMeth=meth;}
  AliTrackPointArray *SortTrack(AliTrackPointArray *atp);
  void      SetTemporaryExcludedModule(Int_t index) {fTempExcludedModule=index;}

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
  Int_t     AddLocalEquation(MilleData &m);
  void      SetLocalEquations(MilleData *m, Int_t neq);
  
  // fitting stuffs
  AliTrackPointArray *GetCurrentTrack() {return fTrack;}
  AliTrackPoint      *GetCurrentCluster() {return &fCluster;}
  void      SetCurrentTrack(AliTrackPointArray *atp) {fTrack=atp;}
  void      SetCurrentCluster(AliTrackPoint &atp) {fCluster=atp;}

  // geometry stuffs
  Int_t     GetNModules() const {return fNModules;}
  Int_t     GetCurrentModuleIndex() const {return fCurrentModuleIndex;}
  TGeoHMatrix *GetCurrentModuleHMatrix() {return fCurrentModuleHMatrix;}
  Double_t *GetCurrentModuleTranslation() {return fCurrentModuleTranslation;}
  Int_t     GetCurrentModuleInternalIndex() const {return fCurrentModuleInternalIndex;}
  Int_t    *GetModuleIndexArray() {return fModuleIndex;}
  Int_t    *GetProcessedPoints() {return fProcessedPoints;}
  Int_t     GetTotBadLocEqPoints() const {return fTotBadLocEqPoints;}
  AliITSAlignMilleModule  *GetMilleModule(UShort_t voluid); // get pointer to the defined supermodule
  AliITSAlignMilleModule  *GetCurrentModule();
  UShort_t *GetModuleVolumeIDArray() {return fModuleVolumeID;}

  // debug stuffs
  Double_t  *GetMeasLoc() { return fMeasLoc;}
  Double_t  *GetSigmaLoc() { return fSigmaLoc;}
  Double_t   GetBField() const {return fBField;}
  Double_t  *GetLocalInitParam() {return fLocalInitParam;}
  Double_t   GetLocalDX() const {return fDerivativeXLoc;}
  Double_t   GetLocalDZ() const {return fDerivativeZLoc;}
  Double_t   GetParSigTranslations() const {return fParSigTranslations;}
  Double_t   GetParSigRotations() const {return fParSigRotations;}
  Int_t      GetPreAlignmentQualityFactor(Int_t index); // if not prealign. return -1
  void       SetBug(Int_t bug) {fBug=bug;} // 1:SSD inversion sens.18-19

 private:

  // configuration methods
  Int_t     LoadConfig(const Char_t *cfile="AliITSAlignMille.conf");
  Int_t     LoadSuperModuleFile(const Char_t *cfile="ITSMilleSuperModules.root");
  void      ResetLocalEquation();
  void      InitGeometry();
  Int_t     ApplyToGeometry();

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
  Double_t      fSigmaXfactor; ///
  Double_t      fSigmaZfactor; ///
  AliAlignObjParams *fTempAlignObj; ///
  Double_t      fDerivativeXLoc; // localX deriv.
  Double_t      fDerivativeZLoc; // localZ deriv.
  Int_t         fMinNPtsPerTrack; ///
  Int_t         fInitTrackParamsMeth; ///
  Int_t        *fProcessedPoints; /// array of statistics of used points per module
  Int_t         fTotBadLocEqPoints; /// total number of reject points because of bad EqLoc
  AliTrackFitterRieman *fRieman; /// riemann fitter for helices
  Bool_t        fRequirePoints;  // required points in specific layers
  Int_t         fNReqLayUp[6];    /// number of points required in layer[n] with Y>0
  Int_t         fNReqLayDown[6];  /// number of points required in layer[n] with Y<0
  Int_t         fNReqLay[6];      /// number of points required in layer[n] 
  Int_t         fNReqDetUp[3];    /// number of points required in Detector[n] with Y>0
  Int_t         fNReqDetDown[3];  /// number of points required in Detector[n] with Y<0
  Int_t         fNReqDet[3];      /// number of points required in Detector[n]
  Int_t         fTempExcludedModule; /// single module temporary excluded from initial fit

  // geometry stuffs
  TString       fGeometryFileName;  ///
  TString       fPreAlignmentFileName;  ///
  TGeoManager  *fGeoManager;        ///
  Int_t         fCurrentModuleIndex;   /// SuperModule index
  Int_t         fCurrentModuleInternalIndex;  /// SuperModule internal index
  Int_t         fCurrentSensVolIndex;   /// Current point (sens. vol.) index
  Double_t      fCurrentModuleTranslation[3]; ///
  Int_t         fNModules;  /// number of defined modules from config file
  Int_t         fModuleIndex[ITSMILLE_NDETELEM*2]; ///
  UShort_t      fModuleVolumeID[ITSMILLE_NDETELEM*2];  ///
  Bool_t        fFreeParam[ITSMILLE_NDETELEM*2][ITSMILLE_NPARCH];  ///
  Bool_t        fUseLocalShifts; /// 
  Bool_t        fUseSuperModules; /// 
  Bool_t        fUsePreAlignment; /// 
  Bool_t        fUseSortedTracks; /// default is kTRUE 
  Bool_t        fBOn; /// magentic field ON
  Double_t      fBField; /// value of magnetic field
  Int_t         fNSuperModules; /// number of custom supermodules in SM file
  TGeoHMatrix  *fCurrentModuleHMatrix; /// SuperModule matrix
  Bool_t        fIsConfigured; ///
  Int_t         fPreAlignQF[ITSMILLE_NDETELEM*2]; ///
  Double_t      fSensVolSigmaXfactor[ITSMILLE_NDETELEM*2]; ///
  Double_t      fSensVolSigmaZfactor[ITSMILLE_NDETELEM*2]; ///
  Int_t         fBug; /// tag for temporary bug correction

  AliITSAlignMilleModule *fMilleModule[ITSMILLE_NDETELEM*2]; /// array of super modules to be aligned

  AliITSAlignMilleModule *fSuperModule[ITSMILLE_NDETELEM*2]; /// array of super modules defined in supermodule file

  AliITSAlignMille(const AliITSAlignMille& rhs);
  AliITSAlignMille& operator=(const AliITSAlignMille& rhs);


  ClassDef(AliITSAlignMille, 0)

};

#endif
