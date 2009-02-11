#ifndef ALIITSALIGNMILLE2_H
#define ALIITSALIGNMILLE2_H
/* Copyright(c) 2007-2009 , ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                                */

/* $Id$ */

/// \ingroup rec
/// \class AliITSAlignMille2
/// \brief Class for alignment of ITS
//
// Authors: Marcello Lunardon

#include <TString.h>
#include <TObject.h>
#include "AliTrackPointArray.h"
#include "AliITSAlignMille2Module.h"

class AliMillePede2;
class AliAlignObjParams;
class TGeoManager;
class TGeoHMatrix;
//class AliITSAlignMille2Module;
class AliTrackFitterRieman;
class TVirtualFitter;
// number of used objects

#define ITSMILLE2_NPARCH         6
#define ITSMILLE2_NLOCAL         5

struct Mille2Data {
  /// structure to store data for 2 LocalEquations (X and Z)
  //
  Int_t    moduleIDX[10];  // max 10 hierarchy levels!!!
  Int_t    levFilled;
  //
  Double_t measX;
  Double_t measZ;
  Double_t sigmaX;
  Double_t sigmaZ;
  //
  Double_t derlocX[ITSMILLE2_NLOCAL];
  Double_t derlocZ[ITSMILLE2_NLOCAL];
  //
  Double_t dergloX[ITSMILLE2_NPARCH*10];
  Double_t dergloZ[ITSMILLE2_NPARCH*10];
};

class AliITSAlignMille2: public TObject
{
 public:
  enum {kNLocal=ITSMILLE2_NLOCAL,kNParCh=ITSMILLE2_NPARCH,
	kMaxITSSensID=2197,kMaxITSSensVID=14300,kMinITSSupeModuleID=14336};
  enum {kDOFTX,kDOFTY,kDOFTZ,kDOFPH,kDOFTH,kDOFPS};
  //
 public:
  //
  AliITSAlignMille2(const Char_t *configFilename="AliITSAlignMille.conf", Bool_t initmille=kTRUE);
  virtual ~AliITSAlignMille2();
  //
  AliMillePede2* GetMillePede()                                  const {return fMillepede;}

  // geometry methods 
  Int_t     GetModuleIndex(const Char_t *symname);
  Int_t     GetModuleIndex(UShort_t voluid);
  UShort_t  GetModuleVolumeID(const Char_t *symname);
  UShort_t  GetModuleVolumeID(Int_t index);
  void      SetParSigTranslations(double v)                           {fParSigTranslations = v;}
  void      SetParSigRotations(double v)                              {fParSigRotations = v;}
  //
  // configuration methods
  void      SetGeometryFileName(const Char_t* filename="geometry.root") { fGeometryFileName = filename; }
  const Char_t* GetGeometryFileName()                                   {return fGeometryFileName.Data();}
  const Char_t* GetPreAlignmentFileName()                               {return fPreAlignmentFileName.Data();}
  TClonesArray* GetPreAlignmentDeltas()                           const {return fPrealignment;}
  void      SetCurrentModule(Int_t id)                                  {fCurrentModule = GetMilleModule(id);}
  void      PrintCurrentModuleInfo()                              const {if (fCurrentModule) fCurrentModule->Print();}
  void      Print(Option_t*)                                      const;
  Bool_t    IsConfigured()                                        const {return fIsConfigured;}
  void      SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts);
  void      SetUseGlobalDelta(Bool_t v=kTRUE)                           {fUseGlobalDelta = v;}
  Bool_t    GetUseGlobalDelta()                                   const {return fUseGlobalDelta;}
  //
  // fitting methods
  AliTrackFitterRieman *GetRiemanFitter()                         const {return fRieman;}
  AliTrackPointArray   *PrepareTrack(const AliTrackPointArray *track); 
  Int_t     ProcessTrack(const AliTrackPointArray *track);
  void      InitTrackParams(int meth=1);
  void      SetMinNPtsPerTrack(Int_t pts=3)                             {fMinNPtsPerTrack=pts;}
  void      SetInitTrackParamsMeth(Int_t meth=1)                        {fInitTrackParamsMeth=meth;}
  Bool_t    InitRiemanFit();
  Int_t     InitModuleParams();
  Int_t     CheckCurrentTrack();
  Bool_t    CheckVolumeID(UShort_t voluid)                        const;
  Int_t     IsDefined(UShort_t voluid)                            const;
  Int_t     IsContained(UShort_t voluid)                          const;
  Int_t     CalcIntersectionPoint(Double_t *lpar, Double_t *gpar);
  Int_t     CalcDerivatives(Int_t paridx, Bool_t islpar);
  Double_t* GetLocalIntersectionPoint()                           const {return (Double_t*)fPintLoc;}
  Double_t* GetGlobalIntersectionPoint()                          const {return (Double_t*)fPintGlo;}
  AliTrackPointArray *SortTrack(const AliTrackPointArray *atp);
  void      SetTemporaryExcludedModule(Int_t index)                     {fTempExcludedModule=index;}
  Int_t     GetTemporaryExcludedModule()                          const {return fTempExcludedModule;}
  //
  // Hierarchical contraints
  void      ConstrainModuleSubUnits(Int_t idm, Double_t val=0, UInt_t pattern=0xfffffff);
  void      ConstrainOrphans(Double_t val=0,UInt_t pattern=0xfffffff);
  void      ConstrainLinComb(const Int_t *modLst,const Float_t *wghLst, Int_t nmd, Double_t val=0, UInt_t pattern=0xfffffff);
  //
  void      PostConstrainModuleSubUnitsMedian(Int_t idm, Double_t val=0, UInt_t pattern=0xfffffff);
  void      PostConstrainOrphansMedian(Double_t val=0, UInt_t pattern=0xfffffff);
  //
  //  Double_t* GetModuleConstrDerivs(Int_t matRow,Int_t matCol, TGeoHMatrix &matLoc);
  //
  // millepede methods
  void      FixParameter(Int_t param, Double_t value);
  void      AddConstraint(Double_t *factor, Double_t value );
  void      InitGlobalParameters(Double_t *par);   
  void      SetLocalDerivative(Int_t index, Double_t value)             {fLocalDerivatives[index] = value;}
  void      SetGlobalDerivative(Int_t index, Double_t value)            {fGlobalDerivatives[index] = value;}  
  //
  Int_t     GlobalFit(Double_t *parameters=0,Double_t *errors=0,Double_t *pulls=0);
  void      PrintGlobalParameters();
  Double_t  GetParError(Int_t iPar);
  Int_t     AddLocalEquation(Mille2Data &m);
  void      SetLocalEquations(const Mille2Data *marr, Int_t neq);
  //
  // fitting stuffs
  AliTrackPointArray *GetCurrentTrack()                                 {return fTrack;}
  AliTrackPoint      *GetCurrentCluster()                               {return &fCluster;}
  void      SetCurrentTrack(AliTrackPointArray *atp)                    {fTrack=atp;}
  void      SetCurrentCluster(AliTrackPoint &atp)                       {fCluster=atp;}

  // geometry stuffs
  Int_t     GetNModules()                   const {return fNModules;}
  Int_t     GetCurrentModuleIndex()         const {return fCurrentModule ? fCurrentModule->GetIndex():-1;}
  TGeoHMatrix *GetCurrentModuleHMatrix()    const {return fCurrentModule ? fCurrentModule->GetMatrix():0;}
  Double_t *GetCurrentModuleTranslation()   const {return fCurrentModule ? fCurrentModule->GetMatrix()->GetTranslation():0;}
  Int_t     GetCurrentModuleInternalIndex() const {return fCurrentModule ? static_cast<Int_t>(fCurrentModule->GetUniqueID()):-1;}
  Int_t     GetTotBadLocEqPoints()          const {return fTotBadLocEqPoints;}
  //
  AliITSAlignMille2Module*  GetMilleModuleByVID(UShort_t voluid) const; // get pointer to the defined supermodule
  AliITSAlignMille2Module*  GetMilleModule(Int_t id)             const {return (AliITSAlignMille2Module*)fMilleModule[id];}
  AliITSAlignMille2Module*  GetCurrentModule()                   const {return fCurrentModule;}
  AliITSAlignMille2Module*  GetSuperModule(Int_t id)             const {return (AliITSAlignMille2Module*)fSuperModule[id];}
  //
  // debug stuffs
  Double_t  *GetMeasLoc()                                        const {return (Double_t*)fMeasLoc;}
  Double_t  *GetSigmaLoc()                                       const {return (Double_t*)fSigmaLoc;}
  Double_t   GetBField()                                         const {return fBField;}
  Double_t  *GetLocalInitParam()                                 const {return (Double_t*)fLocalInitParam;}
  Double_t  *GetLocalInitParEr()                                 const {return (Double_t*)fLocalInitParEr;}
  Double_t   GetLocalDX()                                        const {return fDerivativeLoc[0];}
  Double_t   GetLocalDY()                                        const {return fDerivativeLoc[1];}
  Double_t   GetLocalDZ()                                        const {return fDerivativeLoc[2];}
  Double_t   GetLocalDif(int i)                                  const {return fDerivativeLoc[i];}
  Double_t   GetParSigTranslations()                             const {return fParSigTranslations;}
  Double_t   GetParSigRotations()                                const {return fParSigRotations;}
  Int_t      GetPreAlignmentQualityFactor(Int_t index)           const;// if not prealign. return -1
  void       SetBug(Int_t bug) {fBug=bug;}                             // 1:SSD inversion sens.18-19
  static AliITSAlignMille2* GetInstance()                              {return fgInstance;}
  //
 private:
  //
  // configuration methods
  Int_t     LoadConfig(const Char_t *cfile="AliITSAlignMille.conf");
  Int_t     LoadSuperModuleFile(const Char_t *cfile="ITSMilleSuperModules.root");
  void      ResetLocalEquation();
  void      InitGeometry();
  Int_t     ApplyToGeometry();
  //
  // millepede methods
  void      Init(Int_t nGlobal, Int_t nLocal, Int_t nStdDev);
  //
  AliITSAlignMille2(const AliITSAlignMille2& rhs);
  AliITSAlignMille2& operator=(const AliITSAlignMille2& rhs);
  //
 protected:
  //
  // millepede stuffs
  AliMillePede2 *fMillepede;                    // Detector independent alignment class
  Double_t      fStartFac;                      // Initial value for chi2 cut 
  Double_t      fResCutInitial;                 // Cut on residual for first iteration
  Double_t      fResCut;                        // Cut on residual for other iterations 
  Int_t         fNGlobal;                       // Number of global parameters
  Int_t         fNLocal;                        // Number of local parameters
  Int_t         fNStdDev;                       // Number of standard deviations for chi2 cut
  Bool_t        fIsMilleInit;                   // Flag for initialization
  Bool_t        fSensorsIn;                     // Are sensors explicitly provieded by the conf file?
  Double_t      fParSigTranslations;            // init sigma for transl. params [cm]
  Double_t      fParSigRotations;               // init sigma for rot. params [deg]
  //
  // fitting stuffs
  AliITSAlignMille2Module *fCurrentModule;      // Current SuperModule index
  AliTrackPointArray *fTrack;                   // pointer to current track 
  AliTrackPoint fCluster;                       // current cluster
  Double_t     *fGlobalDerivatives;             // Array of global derivatives
  Double_t      fLocalDerivatives[kNLocal];     // Array of local deriv.
  Double_t      fLocalInitParam[kNLocal];       // Array with inital values for local parameters for current track
  Double_t      fLocalInitParEr[kNLocal][kNLocal];// Array with inital values for local parameters for current track
  Double_t      fModuleInitParam[kNParCh];      // Array with inital values for current module parameters (init geometry)
  Double_t      fPintLoc[3]; 
  Double_t      fPintLoc0[3];
  Double_t      fPintGlo[3]; 
  Double_t      fMeasLoc[3];                    // current point local coordinates (the original ones)
  Double_t      fMeasGlo[3];                    // current point glob. coord (AliTrackPoint)
  Double_t      fSigmaLoc[3];                   // stdev current point
  Double_t      fSigmaFactor[3];                // multiplicative factor for cluster sigmaX,Y,Z
  //
  Double_t      fDerivativeLoc[3];              // localXYZ deriv.
  Int_t         fMinNPtsPerTrack;               // min number of points per track to accept it
  Int_t         fInitTrackParamsMeth;           // method for track fit
  Int_t         fTotBadLocEqPoints;             // total number of reject points because of bad EqLoc
  AliTrackFitterRieman *fRieman;                // riemann fitter for helices
  //
  // >> new members
  Bool_t        fUseGlobalDelta;  // intetpret deltas as global 
  Bool_t        fRequirePoints;   // required points in specific layers
  Int_t         fNReqLayUp[6];    /// number of points required in layer[n] with Y>0
  Int_t         fNReqLayDown[6];  /// number of points required in layer[n] with Y<0
  Int_t         fNReqLay[6];      /// number of points required in layer[n] 
  Int_t         fNReqDetUp[3];    /// number of points required in Detector[n] with Y>0
  Int_t         fNReqDetDown[3];  /// number of points required in Detector[n] with Y<0
  Int_t         fNReqDet[3];      /// number of points required in Detector[n]
  Int_t         fTempExcludedModule; /// single module temporary excluded from initial fit
  // << new members
  //
  // geometry stuffs
  TString       fGeometryFileName;              // Geometry file name
  TString       fPreAlignmentFileName;          // file with prealigned objects
  TGeoManager  *fGeoManager;
  Bool_t        fIsConfigured;
  TArrayS       fPreAlignQF;
  //
  TClonesArray* fPrealignment; // array of prealignment global deltas
  TObjArray     fMilleModule; /// array of super modules to be aligned
  TObjArray     fSuperModule; /// array of super modules defined in supermodule file
  Int_t         fNModules;                      // number of defined modules from config file
  Int_t         fNSuperModules; /// number of custom supermodules in SM file
  Bool_t        fUsePreAlignment;               // start from prealigned setup 
  Bool_t        fUseSortedTracks;               // default is kTRUE 
  Bool_t        fBOn;                           // magentic field ON
  Double_t      fBField;                        // value of magnetic field
  Int_t         fBug;                           /// tag for temporary bug correction
  //
  static AliITSAlignMille2* fgInstance;         // global pointer on itself
  static Int_t              fgInstanceID;       // global counter of the instances
  ClassDef(AliITSAlignMille2, 0)

};

#endif
