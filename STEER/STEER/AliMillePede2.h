#ifndef ALIMILLEPEDE2_H
#define ALIMILLEPEDE2_H

/**********************************************************************************************/
/* General class for alignment with large number of degrees of freedom                        */
/* Based on the original milliped2 by Volker Blobel                                           */
/* http://www.desy.de/~blobel/mptalks.html                                                    */
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include "AliMinResSolve.h"
#include "AliMillePedeRecord.h"
class TFile;
class AliMatrixSq;
class AliSymMatrix;
class AliRectMatrix;
class AliMatrixSparse;
class AliLog;
class TStopwatch;
class TArrayL;


class AliMillePede2: public TObject
{
 public:
  //
  enum {gkFailed,gkInvert,gkNoInversion};    // used global matrix solution methods
  //
  AliMillePede2();
  AliMillePede2(const AliMillePede2& src);
  virtual ~AliMillePede2();
  AliMillePede2& operator=(const AliMillePede2& )             {printf("Dummy\n"); return *this;}
  //
  Int_t                InitMille(int nGlo, int nLoc, Int_t lNStdDev=-1,double lResCut=-1., double lResCutInit=-1.);
  //
  Int_t                GetNGloPar()                     const {return fNGloPar;}
  Int_t                GetNLocPar()                     const {return fNLocPar;}
  Long_t               GetNLocalEquations()             const {return fNLocEquations;}
  Int_t                GetCurrentIteration()            const {return fIter;}
  Int_t                GetNMaxIterations()              const {return fMaxIter;}
  Int_t                GetNStdDev()                     const {return fNStdDev;} 
  Int_t                GetNGlobalConstraints()          const {return fNGloConstraints;}
  Int_t                GetNLagrangeConstraints()        const {return fNLagrangeConstraints;}
  Long_t               GetNLocalFits()                  const {return fNLocFits;}
  Long_t               GetNLocalFitsRejected()          const {return fNLocFitsRejected;}
  Int_t                GetNGlobalsFixed()               const {return fNGloFix;}
  Int_t                GetGlobalSolveStatus()           const {return fGloSolveStatus;}
  Float_t              GetChi2CutFactor()               const {return fChi2CutFactor;}
  Float_t              GetChi2CutRef()                  const {return fChi2CutRef;}
  Float_t              GetResCurInit()                  const {return fResCutInit;}
  Float_t              GetResCut()                      const {return fResCut;}
  Int_t                GetMinPntValid()                 const {return fMinPntValid;}
  Int_t                GetProcessedPoints(Int_t i)      const {return fProcPnt[i];}
  Int_t*               GetProcessedPoints()             const {return fProcPnt;}
  Int_t                GetParamGrID(Int_t i)            const {return fParamGrID[i];}
  //
  AliMatrixSq*         GetGlobalMatrix()                const {return fMatCGlo;}
  AliSymMatrix*        GetLocalMatrix()                 const {return fMatCLoc;}
  Double_t*            GetGlobals()                     const {return fVecBGlo;}
  Double_t*            GetDeltaPars()                   const {return fDeltaPar;}
  Double_t*            GetInitPars()                    const {return fInitPar;}
  Double_t*            GetSigmaPars()                   const {return fSigmaPar;}
  Bool_t*              GetIsLinear()                    const {return fIsLinear;}
  Double_t             GetFinalParam(int i)             const {return fDeltaPar[i]+fInitPar[i];}
  Double_t             GetFinalError(int i)             const {return GetParError(i);}
  //
  Double_t             GetGlobal(Int_t i)               const {return fVecBGlo[i];}
  Double_t             GetInitPar(Int_t i)              const {return fInitPar[i];}
  Double_t             GetSigmaPar(Int_t i)             const {return fSigmaPar[i];}
  Bool_t               GetIsLinear(Int_t i)             const {return fIsLinear[i];}
  static Bool_t        IsGlobalMatSparse()                    {return fgIsMatGloSparse;}
  static Bool_t        IsWeightSigma()                        {return fgWeightSigma;}
  //
  void                 SetParamGrID(Int_t grID,Int_t i)       {fParamGrID[i] = grID; if(fNGroupsSet<grID)fNGroupsSet=grID;}
  void                 SetNGloPar(Int_t n)                    {fNGloPar = n;}
  void                 SetNLocPar(Int_t n)                    {fNLocPar = n;}
  void                 SetNMaxIterations(Int_t n=10)          {fMaxIter = n;}
  void                 SetNStdDev(Int_t n)                    {fNStdDev = n;}
  void                 SetChi2CutFactor(Float_t v)            {fChi2CutFactor = v;}
  void                 SetChi2CutRef(Float_t v)               {fChi2CutRef = v;}
  void                 SetResCurInit(Float_t v)               {fResCutInit = v;}
  void                 SetResCut(Float_t v)                   {fResCut = v;}
  void                 SetMinPntValid(Int_t n)                {fMinPntValid = n>0 ? n:1;}
  static void          SetGlobalMatSparse(Bool_t v=kTRUE)     {fgIsMatGloSparse = v;}
  static void          SetWeightSigma(Bool_t v=kTRUE)         {fgWeightSigma = v;}
  //
  void                 SetInitPars(const Double_t* par)       {memcpy(fInitPar,par,fNGloPar*sizeof(Double_t));}
  void                 SetSigmaPars(const Double_t* par)      {memcpy(fSigmaPar,par,fNGloPar*sizeof(Double_t));}
  void                 SetInitPar(Int_t i,Double_t par)       {fInitPar[i] = par;;}
  void                 SetSigmaPar(Int_t i,Double_t par)      {fSigmaPar[i] = par;}
  //
  Int_t                GlobalFit(Double_t *par=0, Double_t *error=0, Double_t *pull=0);
  Int_t                GlobalFitIteration();
  Int_t                SolveGlobalMatEq();
  static void          SetInvChol(Bool_t v=kTRUE)             {fgInvChol = v;}
  static void          SetMinResPrecondType(Int_t tp=0)       {fgMinResCondType = tp;}
  static void          SetMinResTol(Double_t val=1e-12)       {fgMinResTol = val;}
  static void          SetMinResMaxIter(Int_t val=2000)       {fgMinResMaxIter = val;}
  static void          SetIterSolverType(Int_t val=AliMinResSolve::kSolMinRes) {fgIterSol = val;}
  static void          SetNKrylovV(Int_t val=60)              {fgNKrylovV = val;}
  //
  static Bool_t        GetInvChol()                           {return fgInvChol;}
  static Int_t         GetMinResPrecondType()                 {return fgMinResCondType;}
  static Double_t      GetMinResTol()                         {return fgMinResTol;}
  static Int_t         GetMinResMaxIter()                     {return fgMinResMaxIter;}
  static Int_t         GetIterSolverType()                    {return fgIterSol;}
  static Int_t         GetNKrylovV()                          {return fgNKrylovV;}
  //
  Double_t             GetParError(int iPar)           const;
  Int_t                PrintGlobalParameters()         const;
  void                 SetRejRunList(const UInt_t *runs, Int_t nruns);
  void                 SetAccRunList(const UInt_t *runs, Int_t nruns);
  Bool_t               IsRecordAcceptable() const;
  //
  //
  Int_t                SetIterations(double lChi2CutFac);

  //
  // constraints
  void                 SetGlobalConstraint(double *dergb, double val, double sigma=0);
  void                 SetGlobalConstraint(const int *indgb,double *dergb, int ngb, double val, double sigma=0);
  //
  // processing of the local measurement
  void                 SetRecordRun(Int_t run);
  void                 SetRecordWeight(double wgh);
  void                 SetLocalEquation(double *dergb, double *derlc, double lMeas, double lSigma);
  void                 SetLocalEquation(int *indgb, double *dergb, int ngb, int *indlc, 
					double *derlc,int nlc,double lMeas,double lSigma);
  //
  // manipilation with processed data and costraints records and its buffer
  void                 SetDataRecFName(const char* flname)   {fDataRecFName = flname;}
  const Char_t*        GetDataRecFName()               const {return fDataRecFName.Data();}
  void                 SetConsRecFName(const char* flname)   {fConstrRecFName = flname;}
  const Char_t*        GetConsRecFName()               const {return fConstrRecFName.Data();}

  Bool_t               InitDataRecStorage(Bool_t read=kFALSE);
  Bool_t               InitConsRecStorage(Bool_t read=kFALSE);
  Bool_t               ImposeDataRecFile(const char* fname);
  Bool_t               ImposeConsRecFile(const char* fname);
  void                 CloseDataRecStorage();
  void                 CloseConsRecStorage();
  void                 ReadRecordData(Long_t recID);
  void                 ReadRecordConstraint(Long_t recID);
  Bool_t               ReadNextRecordData();
  Bool_t               ReadNextRecordConstraint();
  void                 SaveRecordData();
  void                 SaveRecordConstraint();
  AliMillePedeRecord*  GetRecord()                      const {return fRecord;}
  Long_t               GetSelFirst()                    const {return fSelFirst;}
  Long_t               GetSelLast()                     const {return fSelLast;}
  void                 SetSelFirst(Long_t v)                  {fSelFirst = v;}
  void                 SetSelLast(Long_t v)                   {fSelLast = v;}
  //
  Float_t              Chi2DoFLim(int nSig, int nDoF)   const;
  //
  // aliases for compatibility with millipede1
  void                 SetParSigma(Int_t i,Double_t par)      {SetSigmaPar(i,par);}
  void                 SetGlobalParameters(Double_t *par)     {SetInitPars(par);}
  //
 protected:
  //
  Int_t                LocalFit(double *localParams=0);
  Bool_t               IsZero(Double_t v,Double_t eps=1e-16)   const {return TMath::Abs(v)<eps;}                  
  //
 protected:
  //
  Int_t                 fNLocPar;                        // number of local parameters
  Int_t                 fNGloPar;                        // number of global parameters
  Int_t                 fNGloSize;                       // final size of the global matrix (NGloPar+NConstraints)
  //
  Long_t                fNLocEquations;                  // Number of local equations 
  Int_t                 fIter;                           // Current iteration
  Int_t                 fMaxIter;                        // Maximum number of iterations
  Int_t                 fNStdDev;                        // Number of standard deviations for chi2 cut 
  Int_t                 fNGloConstraints;                // Number of constraint equations
  Int_t                 fNLagrangeConstraints;           // Number of constraint equations requiring Lagrange multiplier
  Long_t                fNLocFits;                       // Number of local fits
  Long_t                fNLocFitsRejected;               // Number of local fits rejected
  Int_t                 fNGloFix;                        // Number of globals fixed by user
  Int_t                 fGloSolveStatus;                 // Status of global solver at current step
  //
  Float_t               fChi2CutFactor;                  // Cut factor for chi2 cut to accept local fit 
  Float_t               fChi2CutRef;                     // Reference cut for chi2 cut to accept local fit 
  Float_t               fResCutInit;                     // Cut in residual for first iterartion
  Float_t               fResCut;                         // Cut in residual for other iterartiona
  Int_t                 fMinPntValid;                    // min number of points for global to vary
  //
  Int_t                 fNGroupsSet;                     // number of groups set
  Int_t                *fParamGrID;                      // group id for the every parameter
  Int_t                *fProcPnt;                        // N of processed points per global variable
  Double_t             *fVecBLoc;                        // Vector B local (parameters) 
  Double_t             *fDiagCGlo;                       // Initial diagonal elements of C global matrix
  Double_t             *fVecBGlo;                        // Vector B global (parameters)
  //
  Double_t             *fInitPar;                        // Initial global parameters
  Double_t             *fDeltaPar;                       // Variation of global parameters
  Double_t             *fSigmaPar;                       // Sigma of allowed variation of global parameter
  //
  Bool_t               *fIsLinear;                       // Flag for linear parameters
  Bool_t               *fConstrUsed;                     // Flag for used constraints
  //
  Int_t                *fGlo2CGlo;                       // global ID to compressed ID buffer
  Int_t                *fCGlo2Glo;                       // compressed ID to global ID buffer
  //
  // Matrices
  AliSymMatrix         *fMatCLoc;                        // Matrix C local
  AliMatrixSq          *fMatCGlo;                        // Matrix C global
  AliRectMatrix        *fMatCGloLoc;                     // Rectangular matrix C g*l 
  Int_t                *fFillIndex;                      // auxilary index array for fast matrix fill
  Double_t             *fFillValue;                      // auxilary value array for fast matrix fill
  //
  // processed data record bufferization   
  TString               fDataRecFName;                   // Name of File for data records               
  AliMillePedeRecord   *fRecord;                         // Buffer of measurements records
  TFile                *fDataRecFile;                    // File of processed measurements records
  TTree                *fTreeData;                       // Tree of processed measurements records
  Int_t                 fRecFileStatus;                  // state of the record file (0-no, 1-read, 2-rw)
  //
  TString               fConstrRecFName;                 // Name of File for constraints records               
  TTree                *fTreeConstr;                     // Tree of constraint records
  TFile                *fConsRecFile;                    // File of processed constraints records
  Long_t                fCurrRecDataID;                  // ID of the current data record
  Long_t                fCurrRecConstrID;                // ID of the current constraint record
  Bool_t                fLocFitAdd;                      // Add contribution of carrent track (and not eliminate it)
  Int_t                 fSelFirst;                       // event selection start
  Int_t                 fSelLast;                        // event selection end
  TArrayL*              fRejRunList;                     // list of runs to reject (if any)
  TArrayL*              fAccRunList;                     // list of runs to select (if any)
  //
  static Bool_t         fgInvChol;                       // Invert global matrix in Cholesky solver
  static Bool_t         fgWeightSigma;                   // weight parameter constraint by statistics
  static Bool_t         fgIsMatGloSparse;                // Type of the global matrix (sparse ...)
  static Int_t          fgMinResCondType;                // Type of the preconditioner for MinRes method 
  static Double_t       fgMinResTol;                     // Tolerance for MinRes solution
  static Int_t          fgMinResMaxIter;                 // Max number of iterations for the MinRes method
  static Int_t          fgIterSol;                       // type of iterative solution: MinRes or FGMRES
  static Int_t          fgNKrylovV;                      // size of Krylov vectors buffer in FGMRES
  //
  ClassDef(AliMillePede2,0)
};

//_____________________________________________________________________________________________
inline void AliMillePede2::ReadRecordData(Long_t recID)       {fTreeData->GetEntry(recID); fCurrRecDataID=recID;}

//_____________________________________________________________________________________________
inline void AliMillePede2::ReadRecordConstraint(Long_t recID) {fTreeConstr->GetEntry(recID); fCurrRecConstrID=recID;}

//_____________________________________________________________________________________________
inline void AliMillePede2::SaveRecordData()                   {fTreeData->Fill(); fRecord->Reset(); fCurrRecDataID++;}

//_____________________________________________________________________________________________
inline void AliMillePede2::SaveRecordConstraint()             {fTreeConstr->Fill(); fRecord->Reset();fCurrRecConstrID++;}

#endif
