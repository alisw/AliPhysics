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

#include <TArrayI.h>
#include <TArrayD.h>
#include <TString.h>
#include <TObject.h>
#include <TGeoMatrix.h>
#include "AliITSresponseSDD.h"
#include "AliTrackPointArray.h"
#include "AliITSAlignMille2Module.h"
#include "AliITSAlignMille2Constraint.h"
#include "AliITSAlignMille2ConstrArray.h"

class AliMillePede2;
class AliAlignObjParams;
class TGeoManager;
//class AliITSAlignMille2Module;
class AliTrackFitterRieman;
class TVirtualFitter;
// number of used objects



class AliITSAlignMille2: public TObject
{
 public:
 enum {kNLocal=5,kMaxPoints=100,
       kNParChGeom = AliITSAlignMille2Module::kMaxParGeom,
       kNParCh     = AliITSAlignMille2Module::kMaxParTot,
       kMaxITSSensID=2197,kMaxITSSensVID=14300,kMinITSSupeModuleID=14336,kSDDoffsID=240};
  //
 protected:
 struct Mille2Data { // structure to store data for 2 LocalEquations (X and Z)
   enum {kMaxLev = 7};
   Double_t measX, measZ, sigmaX, sigmaZ;
   Double_t derlocX[kNLocal], derlocZ[kNLocal];
   Int_t    nModFilled, nGlobFilled, moduleID[kMaxLev];
   Int_t    parMilleID[AliITSAlignMille2Module::kMaxParTot*kMaxLev];
   Double_t dergloX[AliITSAlignMille2Module::kMaxParTot*kMaxLev];
   Double_t dergloZ[AliITSAlignMille2Module::kMaxParTot*kMaxLev];
 };
 //
 public:
  //
  AliITSAlignMille2(const Char_t *configFilename="AliITSAlignMille.conf");
  virtual ~AliITSAlignMille2();
  //
  AliMillePede2* GetMillePede()                                   const {return fMillepede;}
  //
  // configuration methods
  //
  Int_t     IsVIDDefined(UShort_t voluid)                         const;
  Int_t     IsVIDContained(UShort_t voluid)                       const;
  Int_t     IsSymDefined(const Char_t* name)                      const;
  Int_t     IsSymContained(const Char_t* name)                    const;
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
  //
  AliAlignObjParams*        GetPrealignedObject(const Char_t* symname) const;
  AliAlignObjParams*        GetConstrRefObject(const Char_t* symname) const;
  //
  void          ConvertParamsToGlobal();
  void          ConvertParamsToLocal();
  //
  const Char_t* GetGeometryFileName()                                   {return fGeometryFileName.Data();}
  const Char_t* GetPreAlignmentFileName()                               {return fPreAlignmentFileName.Data();}
  TClonesArray* GetPreAlignmentDeltas()                           const {return fPrealignment;}
  AliITSresponseSDD* GetSDDPrecalibration()                       const {return fCorrectSDD;}
  AliITSresponseSDD* GetSDDInit()                                 const {return fInitialRecSDD;}
  void      PrintCurrentModuleInfo()                              const {if (fCurrentModule) fCurrentModule->Print();}
  void      Print(Option_t*)                                      const;
  Bool_t    IsConfigured()                                        const {return fIsConfigured;}
  Bool_t    GetUseGlobalDelta()                                   const {return fUseGlobalDelta;}
  Bool_t    IsConstraintWrtRef()                                  const {return fConstrRef!=0;}
  Bool_t    FixedOrphans()                                        const;
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
  AliTrackFitterRieman *GetRiemanFitter()                         const {return fRieman;}
  AliTrackPointArray   *PrepareTrack(const AliTrackPointArray *track); 
  AliTrackPointArray *GetCurrentTrack()                                 {return fTrack;}
  AliTrackPoint      *GetCurrentCluster()                               {return &fCluster;}
  void      SetCurrentTrack(AliTrackPointArray *atp)                    {fTrack=atp;}
  void      SetCurrentCluster(AliTrackPoint &atp)                       {fCluster=atp;}
  void      InitTrackParams(int meth=1);
  Int_t     ProcessTrack(const AliTrackPointArray *track);
  Int_t     CheckCurrentTrack();
  //
  Int_t     CalcIntersectionPoint(Double_t *lpar, Double_t *gpar);
  Int_t     CalcDerivatives(Int_t paridx, Bool_t islpar);
  Double_t* GetLocalIntersectionPoint()                           const {return (Double_t*)fPintLoc;}
  Double_t* GetGlobalIntersectionPoint()                          const {return (Double_t*)fPintGlo;}
  AliTrackPointArray *SortTrack(const AliTrackPointArray *atp);
  void      SetTemporaryExcludedModule(Int_t index)                     {fTempExcludedModule=index;}
  Int_t     GetTemporaryExcludedModule()                          const {return fTempExcludedModule;}
  Double_t  GetMeasGlo(Int_t dim)                                 const {return fMeasGlo[dim];}
  Double_t  GetMeasLoc(Int_t dim)                                 const {return fMeasLoc[dim];}
  Int_t     GetCurrentLayer()                                     const;
  //
  // Hierarchical contraints
  Bool_t    PseudoParentsAllowed()                                const {return fAllowPseudoParents;}
  void      ConstrainModuleSubUnitsMean(Int_t idm, Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainModuleSubUnitsMedian(Int_t idm, Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainOrphansMean(Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainOrphansMedian(Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainLocal(const Char_t* name,Double_t *parcf,Int_t npar,Double_t val,Double_t err);
  //
  void      ApplyGaussianConstraint(AliITSAlignMille2ConstrArray* cstr);
  void      ApplyPreConstraints();
  void      ApplyPostConstraints();
  //
  Bool_t    IsParModConstrained(AliITSAlignMille2Module* mod,Int_t par, Bool_t &meanmed, Bool_t &gaussian) const;
  Bool_t    IsParModFamilyVaried(AliITSAlignMille2Module* mod,Int_t par,Int_t depth=999)                   const;
  Bool_t    IsParFamilyFree(AliITSAlignMille2Module* mod,Int_t par,Int_t depth=999)                        const;
  //
  // millepede methods
  Int_t     GlobalFit();
  void      FixParameter(Int_t param, Double_t value);
  void      PrintGlobalParameters();
  Int_t     AddLocalEquation(Mille2Data &m);
  void      SetLocalEquations(const Mille2Data *marr, Int_t neq);
  //
  // module specific 
  //
  Double_t  GetTDriftSDD()                  const;
  Double_t  GetVDriftSDD()                  const;
  //
  AliITSAlignMille2Constraint* GetConstraint(Int_t i)            const {return (AliITSAlignMille2Constraint*)fConstraints.At(i);}
  AliITSAlignMille2Constraint* GetConstraint(const char* name)   const {return (AliITSAlignMille2Constraint*)fConstraints.FindObject(name);}
  //
  // debug stuffs
  void       FetchCluster(const AliTrackPointArray *trc,int ip)        {trc->GetPoint(fCluster,ip);}
  void       SetLocalInitParams(Double_t *par)                         {for (int i=kNLocal;i--;) fLocalInitParam[i]=par[i];}
  Double_t  *GetMeasLoc()                                        const {return (Double_t*)fMeasLoc;}
  Double_t  *GetSigmaLoc()                                       const {return (Double_t*)fSigmaLoc;}
  Double_t   GetBField()                                         const {return fBField;}
  Double_t  *GetLocalInitParam()                                 const {return (Double_t*)fLocalInitParam;}
  Double_t  *GetLocalInitParEr()                                 const {return (Double_t*)fLocalInitParEr;}
  Double_t   GetLocalDif(int par, int coor)                      const {return fDerivativeLoc[par][coor];}
  Double_t   GetGlobalDif(int par, int coor)                     const {return fDerivativeGlo[par][coor];}
  Int_t      GetPreAlignmentQualityFactor(Int_t index)           const;// if not prealign. return -1
  void       SetBug(Int_t bug) {fBug=bug;}                             // 1:SSD inversion sens.18-19
  static     AliITSAlignMille2* GetInstance()                          {return fgInstance;}

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
  void      SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts); 
  Bool_t    InitRiemanFit();
  void      SetMinNPtsPerTrack(Int_t pts=3)  {fMinNPtsPerTrack=pts;}
  //
 private:
  //
  // configuration methods
  void      Init();
  Int_t     LoadConfig(const Char_t *cfile="AliITSAlignMille.conf");
  TObjArray* GetConfigRecord(FILE* stream, TString& recTitle, TString& recOpt, Bool_t rew);
  //
  void      BuildHierarchy();
  Int_t     LoadSuperModuleFile(const Char_t *cfile="ITSMilleSuperModules.root");
  void      ResetLocalEquation();
  Int_t     InitGeometry();
  Int_t     ApplyToGeometry();
  //
  void      ConstrainModuleSubUnits(Int_t idm, Double_t val=0, UInt_t pattern=0xff);
  void      ConstrainOrphans(Double_t val=0,UInt_t pattern=0xff);
  void      PostConstrainModuleSubUnits(Int_t type,Int_t idm, Double_t val, UInt_t pattern);
  void      PostConstrainOrphans(Int_t type,Double_t val, UInt_t pattern);
  //
  void      SetGeometryFileName(const Char_t* filename="geometry.root") { fGeometryFileName = filename; }

  void      SetInitTrackParamsMeth(Int_t meth=1)                        {fInitTrackParamsMeth=meth;}
  //
  void      AddConstraint(Double_t *factor, Double_t value, Double_t sigma=0);
  void      InitGlobalParameters(Double_t *par);   
  void      SetLocalDerivative(Int_t index, Double_t value)             {fLocalDerivatives[index] = value;}
  void      SetGlobalDerivative(Int_t index, Double_t value)            {fGlobalDerivatives[index] = value;}  
  //
  // millepede methods
  //
  void      SetUseGlobalDelta(Bool_t v=kTRUE)                           {fUseGlobalDelta = v;}
  void      SetAllowPseudoParents(Bool_t v=kTRUE)                       {fAllowPseudoParents = v;} 
  Int_t     SetConstraintWrtRef(const char* reffname);
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
  Bool_t        fAllowPseudoParents;            // For simple constraints don't involve parents into the fit
  //
  // fitting stuffs
  AliITSAlignMille2Module *fCurrentModule;      // Current SuperModule index
  AliTrackPointArray *fTrack;                   // pointer to current track 
  TObjArray     fTrackBuff;                     // buffer for tracks of min length
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
  Double_t      fDerivativeLoc[kNLocal][3];    // XYZ deriv. over local params
  Double_t      fDerivativeGlo[kNParCh][3];     // XYZ deriv. over global params
  Int_t         fMinNPtsPerTrack;               // min number of points per track to accept it
  Int_t         fInitTrackParamsMeth;           // method for track fit
  Int_t         fTotBadLocEqPoints;             // total number of reject points because of bad EqLoc
  AliTrackFitterRieman *fRieman;                // riemann fitter for helices
  //
  TObjArray     fConstraints;                   // list of constraints
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
  TString       fConstrRefFileName;             // file with prealigned objects wrt which constraints are defined
  TGeoManager  *fGeoManager;
  Bool_t        fIsConfigured;
  TArrayS       fPreAlignQF;
  //
  AliITSresponseSDD* fCorrectSDD;   // array of SDD t0/vdrift calib params
  AliITSresponseSDD* fInitialRecSDD;   // array of SDD t0/vdrift calib params used to create the track points
  TClonesArray* fPrealignment; // array of prealignment global deltas
  TClonesArray* fConstrRef;    // array of refererence deltas with respect to which the constraint are defined (survey?)
  TObjArray     fMilleModule; /// array of super modules to be aligned
  TObjArray     fSuperModule; /// array of super modules defined in supermodule file
  Int_t         fNModules;                      // number of defined modules from config file
  Int_t         fNSuperModules; /// number of custom supermodules in SM file
  Bool_t        fUsePreAlignment;               // start from prealigned setup 
  Bool_t        fBOn;                           // magentic field ON
  Double_t      fBField;                        // value of magnetic field
  Int_t         fBug;                           /// tag for temporary bug correction
  // pepo
  Int_t         fMilleVersion; /// tag for backward compatibility
  // endpepo
  //
  Double_t      fDriftSpeed[50];                   //temporary array for drift times of SDD alitrackpoints
  Double_t      fDriftTime0[50];                   //temporary array for drift time 0's used for SDD alitrackpoints
  //
  static AliITSAlignMille2* fgInstance;         // global pointer on itself
  static Int_t              fgInstanceID;       // global counter of the instances
  //
  ClassDef(AliITSAlignMille2, 0)
};

#endif
