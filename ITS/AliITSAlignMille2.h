#ifndef ALIITSALIGNMILLE2_H
#define ALIITSALIGNMILLE2_H

/* Copyright(c) 2007-2009 , ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                                */

/* $Id$ */
//-----------------------------------------------------------------------------
//
//  Interface to AliMillePede2 alignment class for the ALICE ITS detector
// 
//  ITS specific alignment class which interface to AliMillepede.   
//  For each track ProcessTrack calculates the local and global derivatives
//  at each hit and fill the corresponding local equations. Provide methods for
//  fixing or constraining detection elements for best results. 
// 
//  author M. Lunardon (thanks to J. Castillo), ruben.shahoyan@cern.ch
//-----------------------------------------------------------------------------

#include <TString.h>
#include <TObject.h>
#include <TGeoMatrix.h>
#include <TArrayS.h>
#include <TArrayD.h>
#include "AliTrackPointArray.h"
#include "AliITSAlignMille2Module.h"

class TSystem;
class TGeoManager;
class TVirtualFitter;
class AliMillePede2;
class AliAlignObjParams;
class AliTrackFitterRieman;
class AliITSAlignMille2Constraint;
class AliITSAlignMille2ConstrArray;
class AliITSresponseSDD;
class AliITSTPArrayFit;
class AliITSsegmentationSDD;
class AliITSDriftSpeedArraySDD;
class AliCDBEntry;

class AliITSAlignMille2: public TObject
{
 public:
  enum {kX,kY,kZ};
  enum {kCosmics, kCollision, kNDataType};
  enum {kNLocal=5,kMaxPoints=100,
	kNParChGeom = AliITSAlignMille2Module::kMaxParGeom,
	kNParCh     = AliITSAlignMille2Module::kMaxParTot,
	kMaxITSSensID=2197,kVtxSensID=kMaxITSSensID+1,kMaxITSSensVID=14300,kVtxSensVID=14371,
	kMinITSSupeModuleID=14336,
	kSDDoffsID=240,kNSDDmod=260};
  //
  enum {kSameInitDeltasBit=BIT(14),kSameInitSDDRespBit=BIT(15),kSameInitSDDVDriftBit=BIT(16),kSameDiamondBit=BIT(17)};
 public:
  //
  AliITSAlignMille2(const Char_t *configFilename="AliITSAlignMille.conf",TList* userInfo=0);
  virtual ~AliITSAlignMille2();
  //
  AliMillePede2* GetMillePede()                                   const {return fMillepede;}
  AliITSTPArrayFit* GetTPAFitter()                                const {return fTPAFitter;}
  //
  // configuration methods
  //
  Int_t     IsVIDDefined(UShort_t voluid)                         const;
  Int_t     IsVIDContained(UShort_t voluid)                       const;
  Int_t     IsSymDefined(const Char_t* name)                      const;
  Int_t     IsSymContained(const Char_t* name)                    const;
  Int_t     GetRequestedModID(UShort_t voluid)                    const;
  //
  Int_t     GetModuleIndex(const Char_t *symname);
  Int_t     GetModuleIndex(UShort_t voluid);
  UShort_t  GetModuleVolumeID(const Char_t *symname);
  UShort_t  GetModuleVolumeID(Int_t index);
  AliITSAlignMille2Module*  GetMilleModuleByVID(UShort_t voluid) const; // get pointer to the defined supermodule
  AliITSAlignMille2Module*  GetMilleModuleBySymName(const Char_t* symname) const; // get pointer to the defined supermodule
  AliITSAlignMille2Module*  GetMilleModuleIfContained(const Char_t* symname) const;
  AliITSAlignMille2Module*  GetMilleModule(Int_t id)             const {return (AliITSAlignMille2Module*)fMilleModule[id];}
  AliITSAlignMille2Module*  GetCurrentModule()                   const {return fCurrentModule;}
  AliITSAlignMille2Module*  GetSuperModule(Int_t id)             const {return (AliITSAlignMille2Module*)fSuperModule[id];}
  AliITSAlignMille2Module*  CreateVertexModule();
  //
  AliAlignObjParams*        GetPrealignedObject(const Char_t* symname) const;
  AliAlignObjParams*        GetConstrRefObject(const Char_t* symname) const;
  //
  void          ConvertParamsToGlobal();
  void          ConvertParamsToLocal();
  //
  const Char_t* GetGeometryPath()                                   {return fGeometryPath.Data();}
  const Char_t* GetPreAlignmentPath()                               {return fPreDeltaPath.Data();}
  TClonesArray* GetPreAlignmentDeltas()                           const {return fPrealignment;}
  AliITSresponseSDD* GetSDDPrecalResp()                           const {return fPreRespSDD;}
  AliITSresponseSDD* GetSDDInitResp()                             const {return fIniRespSDD;}
  TObjArray*         GetSDDInitVDrift()                           const {return fIniVDriftSDD;}
  void      PrintCurrentModuleInfo()                              const {if (fCurrentModule) fCurrentModule->Print();}
  void      Print(Option_t*)                                      const;
  Bool_t    IsConfigured()                                        const {return fIsConfigured;}
  Bool_t    GetUseGlobalDelta()                                   const {return fUseGlobalDelta;}
  Bool_t    IsConstraintWrtRef()                                  const {return fConstrRef!=0;}
  Bool_t    FixedOrphans()                                        const;
  Bool_t    IsLocalYError()                                       const {return fUseLocalYErr;}
  //
  // geometry stuffs
  Int_t     GetNModules()                   const {return fNModules;}
  Int_t     GetCurrentModuleIndex()         const {return fCurrentModule ? fCurrentModule->GetIndex():-1;}
  TGeoHMatrix *GetCurrentModuleHMatrix()    const {return fCurrentModule ? fCurrentModule->GetMatrix():0;}
  Double_t *GetCurrentModuleTranslation()   const {return fCurrentModule ? fCurrentModule->GetMatrix()->GetTranslation():0;}
  Int_t     GetCurrentModuleInternalIndex() const {return fCurrentModule ? Int_t(fCurrentModule->GetUniqueID()):-1;}
  Int_t     GetTotBadLocEqPoints()          const {return fTotBadLocEqPoints;}
  Int_t     GetNConstraints()               const {return fConstraints.GetLast()+1;}
  Int_t     InitModuleParams();
  //
  // fitting methods
  AliTrackFitterRieman *GetRiemanFitter()   const                       {return fRieman;}
  AliTrackPointArray   *PrepareTrack(const AliTrackPointArray *track); 
  AliTrackPointArray *GetCurrentTrack()     const                       {return (AliTrackPointArray*)fTrack;}
  AliTrackPoint      *GetCurrentCluster()   const                       {return (AliTrackPoint*)&fCluster;}
  void      ProcessSDDPointInfo(const AliTrackPoint* pnt,Int_t sID, Int_t pntID);
  void      SetCurrentTrack(const AliTrackPointArray *atp)              {fTrack = (AliTrackPointArray*)atp;}
  void      SetCurrentCluster(const AliTrackPoint &atp);
  void      InitTrackParams(int meth=1);
  Int_t     ProcessTrack(const AliTrackPointArray *track, Double_t wgh=1.0);
  Int_t     FitTrack();
  Int_t     CheckCurrentTrack();
  //
  Int_t     CalcIntersectionPoint(Double_t *lpar, Double_t *gpar);
  Int_t     CalcDerivatives(Int_t paridx, Bool_t islpar);
  void      JacobianPosGloLoc(int locid,double* jacobian);
  Double_t* GetLocalIntersectionPoint()                           const {return (Double_t*)fPintLoc;}
  Double_t* GetGlobalIntersectionPoint()                          const {return (Double_t*)fPintGlo;}
  AliTrackPointArray *SortTrack(const AliTrackPointArray *atp);
  void      SetTemporaryExcludedModule(Int_t index)                     {fTempExcludedModule=index;}
  Int_t     GetTemporaryExcludedModule()                          const {return fTempExcludedModule;}
  Double_t  GetMeasGlo(Int_t dim)                                 const {return fMeasGlo[dim];}
  Double_t  GetMeasLoc(Int_t dim)                                 const {return fMeasLoc[dim];}
  Int_t     GetCurrentLayer()                                     const;
  void      SetBField(Double_t b=0);
  void      SetTypeCosmics()                                            {fDataType = kCosmics;}
  void      SetTypeCollision()                                          {fDataType = kCollision;}
  void      SetDataType(Int_t tp=kCosmics)                              {fDataType = tp>=0&&tp< kNDataType ? tp:kCosmics;}
  void      SetUseLocalYErrors(Bool_t v=kTRUE)                          {fUseLocalYErr = v && fTPAFitter;}
  void      SetMinPointsPerSensor( Int_t n )                            {fMinPntPerSens = n>0 ? n:0;}
  Int_t     GetMinPointsPerSensor()                               const {return fMinPntPerSens;}
  void      ConstrainHelixFitPT(  Int_t q=0,Double_t pt=-1, Double_t err=-1);
  void      ConstrainHelixFitCurv(Int_t q=0,Double_t crv=-1,Double_t crverr=-1);
  void      RemoveHelixFitConstraint();
  Double_t  GetHelixContraintCharge()                             const {return fConstrCharge;}
  Double_t  GetHelixContraintPT()                                 const {return fConstrPT;}
  Double_t  GetHelixContraintPTErr()                              const {return fConstrPTErr;}
  Int_t     GetDataType()                                         const {return fDataType;}
  //
  TGeoHMatrix* GetSensorOrigMatrixSID(Int_t sid)                  const;
  TGeoHMatrix* GetSensorOrigMatrixVID(Int_t vid)                  const;
  //
  TGeoHMatrix* GetSensorCurrMatrixSID(Int_t sid)                  const;
  TGeoHMatrix* GetSensorCurrMatrixVID(Int_t vid)                  const;
  //
  AliCDBEntry* GetCDBEntry(const char* path);
  // Hierarchical contraints
  void      TieSDDVDriftsLR(AliITSAlignMille2Module* mod);
  Bool_t    PseudoParentsAllowed()                                const {return fAllowPseudoParents;}
  void      ConstrainModuleSubUnitsMean(Int_t idm, Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainModuleSubUnitsMedian(Int_t idm, Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainOrphansMean(Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainOrphansMedian(Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainLocal(const Char_t* name,Double_t *parcf,Int_t npar,Double_t val,Double_t err);
  //
  void      ApplyGaussianConstraint(const AliITSAlignMille2ConstrArray* cstr);
  void      ApplyPreConstraints();
  void      ApplyPostConstraints();
  //
  void      SetWeightPt(Double_t w=1)                                   {fWeightPt = w;}
  void      SetSDDVDCorrMult(Bool_t v=kTRUE)                            {fIsSDDVDriftMult=v;}   
  Double_t  GetWeightPt()                                         const {return fWeightPt;}
  Bool_t    IsSDDVDCorrMult()                                     const {return fIsSDDVDriftMult;}
  Bool_t    IsParModConstrained(const AliITSAlignMille2Module* mod,Int_t par, Bool_t &meanmed, Bool_t &gaussian) const;
  Bool_t    IsParModFamilyVaried(const AliITSAlignMille2Module* mod,Int_t par,Int_t depth=999)             const;
  Bool_t    IsParFamilyFree(const AliITSAlignMille2Module* mod,Int_t par,Int_t depth=999)                  const;
  //
  // millepede methods
  Int_t     GlobalFit();
  void      FixParameter(Int_t param, Double_t value);
  void      PrintGlobalParameters();
  //
  TClonesArray*      CreateDeltas();
  AliITSresponseSDD* CreateSDDResponse();
  // module specific 
  //
  Double_t  GetTDriftSDD()                  const;
  Double_t  GetVDriftSDD()                  const;
  Double_t  GetDriftSpeed(Int_t id)         const {return fDriftSpeed[id];}
  Double_t  GetDriftSpeed0(Int_t id)        const {return fDriftSpeed0[id];}
  Double_t  GetDriftTime0(Int_t id)         const {return fDriftTime0[id];}

  //
  AliITSAlignMille2Constraint* GetConstraint(Int_t i)            const {return (AliITSAlignMille2Constraint*)fConstraints.At(i);}
  AliITSAlignMille2Constraint* GetConstraint(const char* name)   const {return (AliITSAlignMille2Constraint*)fConstraints.FindObject(name);}
  //
  // debug stuffs
  void       FetchCluster(int ip)                                      {fTrack->GetPoint(fCluster,ip);fCluster.SetUniqueID(ip);} 
  void       SetLocalInitParams(const Double_t *par)                   {for (int i=kNLocal;i--;) fLocalInitParam[i]=par[i];}
  Bool_t     IsTypeCosmics()                                     const {return fDataType==kCosmics;}
  Bool_t     IsTypeCollision()                                   const {return fDataType==kCollision;}
  Double_t  *GetMeasLoc()                                        const {return (Double_t*)fMeasLoc;}
  Double_t  *GetSigmaLoc()                                       const {return (Double_t*)fSigmaLoc;}
  Double_t   GetBField()                                         const {return fBField;}
  Bool_t     IsFieldON()                                         const {return fBOn;}
  Bool_t     IsDiamondUsed()                                     const {return fUseDiamond;}
  Double_t  *GetLocalInitParam()                                 const {return (Double_t*)fLocalInitParam;}
  Double_t  *GetLocalInitParEr()                                 const {return (Double_t*)fLocalInitParEr;}
  Double_t   GetLocalDif(int par, int coor)                      const {return fDerivativeLoc[par][coor];}
  Double_t   GetGlobalDif(int par, int coor)                     const {return fDerivativeGlo[par][coor];}
  Int_t      GetPreAlignmentQualityFactor(Int_t index)           const;// if not prealign. return -1
  void       SetBug(Int_t bug) {fBug=bug;}                             // 1:SSD inversion sens.18-19
  static     AliITSAlignMille2* GetInstance()                          {return fgInstance;}

  // pepo270809
  Int_t      GetExtraClustersMode() const {return fExtraClustersMode;}
  void       SetExtraClustersMode(Int_t mode) {fExtraClustersMode=mode;}  
  // endpepo270809

  // pepo
  // flag for AliITSAlignMille compatibility
  Int_t      GetMilleVersion() const {return fMilleVersion;}
  void       SetMilleVersion(Int_t m1) {fMilleVersion=m1;}
  // modified existing methods
  void      SetCurrentModule(Int_t id);
  // old methods recovered
  Int_t     IsDefined(UShort_t voluid) const {return IsVIDDefined(voluid);}
  Int_t     IsContained(UShort_t voluid) const {return IsVIDContained(voluid);}
  // moved from private to public
  void      SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts,Int_t runtype=-1); 
  Bool_t    InitRiemanFit();
  void      SetMinNPtsPerTrack(Int_t pts=3)  {fMinNPtsPerTrack=pts;}
  //
  static Bool_t    IsZero(Double_t v,Double_t threshold = 1e-15)       { return TMath::Abs(v)<threshold; }
  static void      SetWordBit(UInt_t word,Int_t bitID)                 { word |= (1<<bitID);}
  static void      ResetWordBit(UInt_t word,Int_t bitID)               { word &= ~(1<<bitID);}
  static Bool_t    TestWordBit(UInt_t word,Int_t bitID)                { return (Bool_t)(word&(1<<bitID));}      
  //
 protected:
  //
  struct Mille2Data { // structure to store data for LocalEquations (X and Z, optionally Y)
    enum {kMaxLev = 7};
    Double_t fMeas[3];                // measured coordinates
    Double_t fSigma[3];               // measured errors
    Double_t fDerLoc[kNLocal][3];     // calculated local derivatives
    Int_t    fNModFilled, fNGlobFilled, fModuleID[kMaxLev]; // used module info
    Int_t    fParMilleID[AliITSAlignMille2Module::kMaxParTot*kMaxLev]; // param id's
    Double_t fDerGlo[AliITSAlignMille2Module::kMaxParTot*kMaxLev][3]; // global derivatives
  };
  //
  // configuration methods
  void      Init();
  Int_t     CacheMatricesOrig();
  Int_t     CacheMatricesCurr();
  Int_t     ProcessUserInfo(TList *userInfo=0);
  Int_t     LoadConfig(const Char_t *cfile="AliITSAlignMille.conf");
  TObjArray* GetConfigRecord(FILE* stream, TString& recTitle, TString& recOpt, Bool_t rew);
  Int_t     CheckConfigRecords(FILE* stream);
  Int_t     ReloadInitCalib(TList *userInfo);
  Int_t     ReloadInitCalib();
  //
  void      BuildHierarchy();
  Int_t     LoadSuperModuleFile(const Char_t *cfile="ITSMilleSuperModules.root");
  Int_t     LoadSDDResponse(TString& path, AliITSresponseSDD *&resp);
  Int_t     LoadSDDVDrift(TString& path, TObjArray *&arr);
  Int_t     LoadDeltas(TString& path, TClonesArray *&arr);
  Int_t     LoadDiamond(TString& path);
  void      ResetLocalEquation();
  Int_t     InitGeometry();
  Int_t     ApplyToGeometry();
  //
  void      ConstrainModuleSubUnits(Int_t idm, Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainOrphans(Double_t val=0,UInt_t pattern=0xff);
  void      PostConstrainModuleSubUnits(Int_t type,Int_t idm, Double_t val, UInt_t pattern);
  void      PostConstrainOrphans(Int_t type,Double_t val, UInt_t pattern);
  //
  void      SetGeometryPath(const Char_t* filename="geometry.root") { fGeometryPath = filename; }

  void      SetInitTrackParamsMeth(Int_t meth=1)                        {fIniTrackParamsMeth=meth;}
  //
  void      AddConstraint(Double_t *factor, Double_t value, Double_t sigma=0);
  void      InitGlobalParameters(Double_t *par);
  Bool_t    SetLocalDerivative(Int_t index, Double_t value)             {return IsZero(fLocalDerivatives[index]=value);}
  Bool_t    SetGlobalDerivative(Int_t index, Double_t value)            {return IsZero(fGlobalDerivatives[index]=value);}  
  //
  // millepede methods
  //
  Int_t     AddLocalEquation(Mille2Data &m);
  Int_t     AddLocalEquationTPA(Mille2Data &m);
  void      SetLocalEquations(const Mille2Data *marr, Int_t neq);
  void      SetUseGlobalDelta(Bool_t v=kTRUE)                           {fUseGlobalDelta = v;}
  void      SetAllowPseudoParents(Bool_t v=kTRUE)                       {fAllowPseudoParents = v;} 
  Int_t     SetConstraintWrtRef(const char* reffname);
  //
  AliITSAlignMille2(const AliITSAlignMille2& rhs);
  AliITSAlignMille2& operator=(const AliITSAlignMille2& rhs);
  //
 protected:
  //
  enum {
    kOCDBDefaultPath,
    kOCDBSpecificPath,
    kGeomFile,
    kSuperModileFile,
    kConstrRefFile,
    kPreDeltaFile,
    kPreCalSDDFile,
    kPreVDriftSDDFile,
    kInitCalSDDFile,
    kInitVDriftSDDFile,
    kInitDeltaFile,
    kGlobalDeltas,
    kConstrLocal,
    kModVolID,
    kModIndex,
    kPseudoParents,
    kTrackFitMethod,
    kMinPntTrack,
    kNStDev,
    kResCutInit,
    kResCutOther,
    kLocalSigFactor,
    kStartFactor,
    kFinalFactor,
    kBField,
    kSparseMatrix,
    kRequirePoint,
    kConstrOrphans,
    kConstrSubunits,
    kApplyConstr,
    kExtraClustersMode,
    kTPAFitter,
    kUseLocalYErr,
    kMinPointsSens,
    kSDDVDCorrMult,
    kWeightPt,
    kUseDiamond,
    kSameSDDT0,
    //
    kNKeyWords
  };                                            // id's of the keywirds for config file records

  // millepede stuffs
  AliMillePede2 *fMillepede;                    // Detector independent alignment class
  Double_t      fStartFac;                      // Initial factor for chi2 cut 
  Double_t      fFinalFac;                      // Final factor for chi2 cut 
  Double_t      fResCutInitial;                 // Cut on residual for first iteration
  Double_t      fResCut;                        // Cut on residual for other iterations 
  Int_t         fNGlobal;                       // Number of global parameters
  Int_t         fNLocal;                        // Number of local parameters
  Int_t         fNStdDev;                       // Number of standard deviations for chi2 cut
  Bool_t        fIsMilleInit;                   // Flag for initialization
  Bool_t        fAllowPseudoParents;            // For simple constraints don't involve parents into the fit
  //
  // fitting stuffs
  AliITSTPArrayFit        *fTPAFitter;          // TPArrayFitter       
  AliITSAlignMille2Module *fCurrentModule;      // Current SuperModule index
  AliTrackPointArray *fTrack;                   // pointer to current track 
  TObjArray     fTrackBuff;                     // buffer for tracks of min length
  AliTrackPoint fCluster;                       // current cluster
  Int_t         fCurrentSensID;                 // sensor index for current cluster
  TArrayD       fClusLoc;                       // local  coordinates of the clusters
  TArrayD       fClusGlo;                       // global coordinates of the clusters
  TArrayD       fClusSigLoc;                    // local cov matrix of the clusters
  Double_t     *fGlobalDerivatives;             // Array of global derivatives
  Double_t      fLocalDerivatives[kNLocal];     // Array of local deriv.
  Double_t      fLocalInitParam[kNLocal];       // Array with inital values for local parameters for current track
  Double_t      fLocalInitParEr[kNLocal][kNLocal];// Array with inital values for local parameters for current track
  Double_t      fModuleInitParam[kNParCh];      // Array with inital values for current module parameters (init geometry)
  Double_t      fPintLoc[3];                    // track/module intersection point in local coordinates
  Double_t      fPintLoc0[3];                   // track/module intersection point in local coordinates (before variation)
  Double_t      fPintGlo[3];                    // track/module intersection point in global coordinates
  Double_t     *fMeasLoc;                       // current point local coordinates (the original ones)
  Double_t     *fMeasGlo;                       // current point glob. coord (AliTrackPoint)
  Double_t     *fSigmaLoc;                      // stdev current point
  Double_t      fSigmaFactor[3];                // multiplicative factor for cluster sigmaX,Y,Z
  Double_t      fConstrPT;                      // optional PT constraint for helix (abs value)
  Double_t      fConstrPTErr;                   // error on this constraint (0 - exact)
  Int_t         fConstrCharge;                  // optional constraint on charge of Helix track (0 - no constraint)
  //
  Double_t      fDerivativeLoc[kNLocal][3];     // XYZ deriv. over local params
  Double_t      fDerivativeGlo[kNParCh][3];     // XYZ deriv. over global params
  Int_t         fMinNPtsPerTrack;               // min number of points per track to accept it
  Int_t         fIniTrackParamsMeth;            // method for track fit
  Int_t         fTotBadLocEqPoints;             // total number of reject points because of bad EqLoc
  AliTrackFitterRieman *fRieman;                // riemann fitter for helices
  //
  TObjArray     fConstraints;                   // list of constraints
  TObjArray     fCacheMatrixOrig;               // cach for original geom matrices
  TObjArray     fCacheMatrixCurr;               // cach for prealigned geom matrices
  // >> new members
  Bool_t        fUseGlobalDelta;                // intetpret deltas as global 
  Bool_t        fRequirePoints[kNDataType];     // required points in specific layers
  Int_t         fNReqLayUp[kNDataType][6];      // number of points required in layer[n] with Y>0
  Int_t         fNReqLayDown[kNDataType][6];    // number of points required in layer[n] with Y<0
  Int_t         fNReqLay[kNDataType][6];        // number of points required in layer[n] 
  Int_t         fNReqDetUp[kNDataType][3];      // number of points required in Detector[n] with Y>0
  Int_t         fNReqDetDown[kNDataType][3];    // number of points required in Detector[n] with Y<0
  Int_t         fNReqDet[kNDataType][3];        // number of points required in Detector[n]
  Int_t         fTempExcludedModule; /// single module temporary excluded from initial fit
  // << new members
  //
  // OCDB stuff
  TList        *fIniUserInfo;                   // initial user info (validity is not guaranteed after initialization)
  TString       fIniDeltaPath;                  // where to take the deltas used to produce the points
  TString       fIniSDDRespPath;                // where to take the initial SDD response used to produce the points
  TString       fPreCalSDDRespPath;             // precalibration SDD response file name
  TString       fIniSDDVDriftPath;              // initial SDD vdrift file name
  TString       fPreSDDVDriftPath;              // initial SDD vdrift file name
  // geometry stuffs
  TString       fGeometryPath;                  // Geometry file name
  TString       fPreDeltaPath;                  // file with prealigned objects
  TString       fConstrRefPath;                 // file with prealigned objects wrt which constraints are defined
  TString       fDiamondPath;                   // file with diamond constraint
  TGeoManager  *fGeoManager;                    // pointer to Alice geomanager
  Bool_t        fIsConfigured;                  // flag for loaded config file
  TArrayS       fPreAlignQF;                    // prealignment flags (not used?)
  //
  AliITSresponseSDD* fIniRespSDD;               // array of SDD t0/vdrift calib params used to create the track points
  AliITSresponseSDD* fPreRespSDD;               // array of SDD t0/vdrift calib params
  TObjArray*         fIniVDriftSDD;             // array of AliITSDriftSpeedArraySDD objects used for original reco
  TObjArray*         fPreVDriftSDD;             // array of AliITSDriftSpeedArraySDD objects to be used as a starting point instead of fIniVDriftSDD
  AliITSsegmentationSDD* fSegmentationSDD;      // extraction of SDD segmentation params
  TClonesArray* fPrealignment; // array of prealignment global deltas
  TClonesArray* fConstrRef;    // array of refererence deltas with respect to which the constraint are defined (survey?)
  TObjArray     fMilleModule; /// array of super modules to be aligned
  TObjArray     fSuperModule; /// array of super modules defined in supermodule file
  Int_t         fNModules;                      // number of defined modules from config file
  Int_t         fNSuperModules; /// number of custom supermodules in SM file
  Bool_t        fUsePreAlignment;               // start from prealigned setup 
  Bool_t        fUseLocalYErr;                  // use local Yerror due to the sensor thickness
  Bool_t        fBOn;                           // magentic field ON
  Double_t      fBField;                        // value of magnetic field
  Int_t         fDataType;                      // is this cosmics or collision processing?
  Int_t         fMinPntPerSens;                 // min number of points per module to vary it
  Int_t         fBug;                           /// tag for temporary bug correction
  // pepo
  Int_t         fMilleVersion; /// tag for backward compatibility
  // endpepo
  // pepo270809
  Int_t         fExtraClustersMode; /// 1=remove random / 2=remove internal / 10=use only tracks with xcl
  // endpepo270809
  //
  Double_t      fTrackWeight;                      //weight given by the user to current track
  Double_t      fWeightPt;                         //weight track equations by pT in this power
  Bool_t        fIsSDDVDriftMult;                  //use multiplicative correction for SDD vdrift
  Double_t      fDriftSpeed[50];                   //temporary array for corrected drift speed of SDD alitrackpoints
  Double_t      fDriftSpeed0[50];                  //temporary array for original  drift speed of SDD alitrackpoints
  Double_t      fDriftTime0[50];                   //temporary array for drift time 0's used for SDD alitrackpoints
  Double_t      fExtClusterPar[9];                 //array to store the parameters of the externally imposed cluster
  AliTrackPoint fDiamond;                          //optional constraint on the vertex
  AliTrackPoint fDiamondI;                         //constraint on the vertex with inverted error matrix
  Bool_t        fUseDiamond;                       //use diamond as a vertex constraint
  Int_t         fDiamondPointID;                   //ID of the diamond point in the track
  Int_t         fDiamondModID;                     //id of the fake diamond module
  //
  static AliITSAlignMille2* fgInstance;         // global pointer on itself
  static Int_t              fgInstanceID;       // global counter of the instances
  static const Char_t     * fgkRecKeys[];       // keywords for config file records
  static const Char_t       fgkXYZ[];           // XYZ labels
  //
  ClassDef(AliITSAlignMille2, 0)
};


//______________________________________________________________________________________
inline void AliITSAlignMille2::SetCurrentCluster(const AliTrackPoint &atp) 
{
  // set current cluster
  fCluster = atp; 
  fCurrentSensID = AliITSAlignMille2Module::GetIndexFromVolumeID(fCluster.GetVolumeID());
}

//______________________________________________________________________________________
inline TGeoHMatrix* AliITSAlignMille2::GetSensorOrigMatrixSID(Int_t sid) const 
{
  // get cached original matrix by sensor ID
  return sid<0 ? 0 : (TGeoHMatrix*) fCacheMatrixOrig[sid];
}

//______________________________________________________________________________________
inline TGeoHMatrix* AliITSAlignMille2::GetSensorOrigMatrixVID(Int_t vid) const
{
  // get cached original matrix by sensor volume ID
  return GetSensorOrigMatrixSID( AliITSAlignMille2Module::GetIndexFromVolumeID(vid) );
}

//______________________________________________________________________________________
inline TGeoHMatrix* AliITSAlignMille2::GetSensorCurrMatrixSID(Int_t sid) const 
{
  // get cached current matrix by sensor ID
  return sid<0 ? 0 : (TGeoHMatrix*) fCacheMatrixCurr[sid];
}

//______________________________________________________________________________________
inline TGeoHMatrix* AliITSAlignMille2::GetSensorCurrMatrixVID(Int_t vid) const
{
  // get cached current matrix by sensor volume ID
  return GetSensorCurrMatrixSID( AliITSAlignMille2Module::GetIndexFromVolumeID(vid) );
}

#endif

