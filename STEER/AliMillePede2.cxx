/**********************************************************************************************/
/* General class for alignment with large number of degrees of freedom                        */
/* Based on the original milliped2 by Volker Blobel                                           */
/* http://www.desy.de/~blobel/mptalks.html                                                    */
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

#include "AliMillePede2.h"
#include "AliLog.h"
#include <TStopwatch.h>
#include <TFile.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TArrayL.h>
#include "AliMatrixSq.h"
#include "AliSymMatrix.h"
#include "AliRectMatrix.h"
#include "AliMatrixSparse.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;


ClassImp(AliMillePede2)

Bool_t   AliMillePede2::fgInvChol        = kTRUE;     // Invert global matrix with Cholesky solver
Bool_t   AliMillePede2::fgWeightSigma    = kTRUE;     // weight local constraint by module statistics
Bool_t   AliMillePede2::fgIsMatGloSparse = kFALSE;    // use faster dense matrix by default
Int_t    AliMillePede2::fgMinResCondType = 1;         // Jacoby preconditioner by default
Double_t AliMillePede2::fgMinResTol      = 1.e-11;    // default tolerance
Int_t    AliMillePede2::fgMinResMaxIter  = 10000;     // default max number of iterations 
Int_t    AliMillePede2::fgIterSol        = AliMinResSolve::kSolMinRes; // default iterative solver
Int_t    AliMillePede2::fgNKrylovV       = 240;        // default number of Krylov vectors to keep

//_____________________________________________________________________________________________
AliMillePede2::AliMillePede2() 
: fNLocPar(0),
  fNGloPar(0),
  fNGloSize(0),
//
  fNLocEquations(0),
  fIter(0),
  fMaxIter(10),
  fNStdDev(3),
  fNGloConstraints(0),
  fNLagrangeConstraints(0),
  fNLocFits(0),
  fNLocFitsRejected(0),
  fNGloFix(0),
  fGloSolveStatus(gkFailed),
//
  fChi2CutFactor(1.),
  fChi2CutRef(1.),
  fResCutInit(100.),
  fResCut(100.),
  fMinPntValid(1),
//
  fNGroupsSet(0),
  fParamGrID(0),
  fProcPnt(0),
  fVecBLoc(0),
  fDiagCGlo(0),
  fVecBGlo(0),
  fInitPar(0),
  fDeltaPar(0),
  fSigmaPar(0),
  fIsLinear(0),
  fConstrUsed(0),
//
  fGlo2CGlo(0),
  fCGlo2Glo(0),
//
  fMatCLoc(0),
  fMatCGlo(0),
  fMatCGloLoc(0),
  //
  fFillIndex(0),
  fFillValue(0),
  //
  fDataRecFName("/tmp/mp2_data_records.root"),
  fRecord(0),
  fDataRecFile(0),  
  fTreeData(0),
  fRecFileStatus(0),
  //
  fConstrRecFName("/tmp/mp2_constraints_records.root"),
  fTreeConstr(0),
  fConsRecFile(0),
  fCurrRecDataID(0),
  fCurrRecConstrID(0),
  fLocFitAdd(kTRUE),
  fSelFirst(1),
  fSelLast(-1),
  fRejRunList(0),
  fAccRunList(0)
{}

//_____________________________________________________________________________________________
AliMillePede2::AliMillePede2(const AliMillePede2& src) : 
  TObject(src),fNLocPar(0),fNGloPar(0),fNGloSize(0),fNLocEquations(0),fIter(0),
  fMaxIter(10),fNStdDev(3),fNGloConstraints(0),fNLagrangeConstraints(0),
  fNLocFits(0),fNLocFitsRejected(0),
  fNGloFix(0),fGloSolveStatus(0),fChi2CutFactor(0),fChi2CutRef(0),fResCutInit(0),
  fResCut(0),fMinPntValid(1),fNGroupsSet(0),fParamGrID(0),fProcPnt(0),fVecBLoc(0),fDiagCGlo(0),fVecBGlo(0),
  fInitPar(0),fDeltaPar(0),fSigmaPar(0),fIsLinear(0),fConstrUsed(0),fGlo2CGlo(0),fCGlo2Glo(0),
  fMatCLoc(0),fMatCGlo(0),fMatCGloLoc(0),fFillIndex(0),fFillValue(0),
  fDataRecFName(0),fRecord(0),fDataRecFile(0),
  fTreeData(0),fRecFileStatus(0),fConstrRecFName(0),fTreeConstr(0),fConsRecFile(0),fCurrRecDataID(0),
  fCurrRecConstrID(0),fLocFitAdd(kTRUE),
  fSelFirst(1),
  fSelLast(-1),
  fRejRunList(0),
  fAccRunList(0)
{printf("Dummy\n");}

//_____________________________________________________________________________________________
AliMillePede2::~AliMillePede2() 
{
  CloseDataRecStorage();
  CloseConsRecStorage();
  //
  delete[] fParamGrID;
  delete[] fProcPnt;
  delete[] fVecBLoc;
  delete[] fDiagCGlo;
  delete[] fVecBGlo;
  delete[] fInitPar;
  delete[] fDeltaPar;
  delete[] fSigmaPar;
  delete[] fGlo2CGlo;
  delete[] fCGlo2Glo;
  delete[] fIsLinear;
  delete[] fConstrUsed;
  delete[] fFillIndex;
  delete[] fFillValue;
  //
  delete fRecord;
  delete fMatCLoc;
  delete fMatCGlo;
  delete fMatCGloLoc;
  delete fRejRunList;
  delete fAccRunList;
} 

//_____________________________________________________________________________________________
Int_t AliMillePede2::InitMille(int nGlo, int nLoc, int lNStdDev,double lResCut, double lResCutInit)
{
  //
  if (nLoc>0)        fNLocPar = nLoc;
  if (nGlo>0)        fNGloPar = nGlo;
  if (lResCutInit>0) fResCutInit = lResCutInit; 
  if (lResCut>0)     fResCut     = lResCut;
  if (lNStdDev>0)    fNStdDev    = lNStdDev;
  //
  fNGloSize = fNGloPar;
  //
  try {
    //
    if (fgIsMatGloSparse) {fMatCGlo = new AliMatrixSparse(fNGloPar); fMatCGlo->SetSymmetric(kTRUE);}
    else                   fMatCGlo = new AliSymMatrix(fNGloPar);
    //
    fFillIndex    = new Int_t[fNGloPar];
    fFillValue    = new Double_t[fNGloPar];
    //
    fMatCLoc      = new AliSymMatrix(fNLocPar);
    fMatCGloLoc   = new AliRectMatrix(fNGloPar,fNLocPar);
    //
    fParamGrID    = new Int_t[fNGloPar];
    fProcPnt      = new Int_t[fNGloPar];
    fVecBLoc      = new Double_t[fNLocPar];
    fDiagCGlo     = new Double_t[fNGloPar];
    //
    fInitPar      = new Double_t[fNGloPar];
    fDeltaPar     = new Double_t[fNGloPar];
    fSigmaPar     = new Double_t[fNGloPar];
    fIsLinear     = new Bool_t[fNGloPar];
    //
    fGlo2CGlo     = new Int_t[fNGloPar];
    fCGlo2Glo     = new Int_t[fNGloPar];
  }
  catch(bad_alloc&) {
    AliInfo(Form("Failed to allocate the memory for %d global and %d local parameters",fNGloPar,fNLocPar));
    return 0;
  }
  //
  memset(fVecBLoc   ,0,fNLocPar*sizeof(Double_t));
  memset(fDiagCGlo  ,0,fNGloPar*sizeof(Double_t));
  memset(fInitPar   ,0,fNGloPar*sizeof(Double_t));
  memset(fDeltaPar  ,0,fNGloPar*sizeof(Double_t));
  memset(fSigmaPar  ,0,fNGloPar*sizeof(Double_t));
  memset(fProcPnt   ,0,fNGloPar*sizeof(Int_t));
  //
  for (int i=fNGloPar;i--;) {
    fGlo2CGlo[i] = fCGlo2Glo[i] = -1;
    fIsLinear[i] = kTRUE;
    fParamGrID[i] = -1;
  }
  //
  return 1;
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::ImposeDataRecFile(const char* fname)
{
  CloseDataRecStorage();
  SetDataRecFName(fname);
  return InitDataRecStorage(kTRUE); // open in read mode
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::ImposeConsRecFile(const char* fname)
{
  CloseConsRecStorage();
  SetConsRecFName(fname);
  return InitConsRecStorage(kTRUE); // open in read mode
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::InitDataRecStorage(Bool_t read)
{
  // initialize the buffer for processed measurements records
  //
  if (fDataRecFile) {AliInfo("Data Records File is already initialized"); return kFALSE;} 
  //
  if (!fRecord) fRecord = new AliMillePedeRecord();
  //
  fDataRecFile = TFile::Open(GetDataRecFName(),read ? "":"recreate");
  if (!fDataRecFile) {AliFatal(Form("Failed to initialize data records file %s",GetDataRecFName())); return kFALSE;}
  //
  AliInfo(Form("File %s used for derivatives records",GetDataRecFName()));
  if (read) {
    fTreeData = (TTree*)fDataRecFile->Get("AliMillePedeRecords_Data");
    if (!fTreeData) {AliFatal(Form("Did not find data records tree in %s",GetDataRecFName())); return kFALSE;}
    fTreeData->SetBranchAddress("Record_Data",&fRecord);
    AliInfo(Form("Found %d derivatives records",fTreeData->GetEntries()));
  }
  else {
    fTreeData = new TTree("AliMillePedeRecords_Data","Data Records for AliMillePede2");
    fTreeData->Branch("Record_Data","AliMillePedeRecord",&fRecord,32000,99);
  }
  fCurrRecDataID = -1;
  fRecFileStatus = read ? 1:2;
  //
  return kTRUE;
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::InitConsRecStorage(Bool_t read)
{
  // initialize the buffer for processed measurements records
  //
  if (fConsRecFile) {AliInfo("Constraints Records File is already initialized"); return kFALSE;} 
  //
  if (!fRecord) fRecord = new AliMillePedeRecord();
  //
  fConsRecFile = TFile::Open(GetConsRecFName(),read ? "":"recreate");
  if (!fConsRecFile) {AliInfo(Form("Failed to initialize constraints records file %s",GetConsRecFName())); return kFALSE;}
  //
  AliInfo(Form("File %s used for constraints records",GetConsRecFName()));
  if (read) {
    fTreeConstr = (TTree*)fConsRecFile->Get("AliMillePedeRecords_Constraints");
    if (!fTreeConstr) {AliInfo(Form("Did not find constraints records tree in %s",GetConsRecFName())); return kFALSE;}
    fTreeConstr->SetBranchAddress("Record_Constraints",&fRecord);
    AliInfo(Form("Found %d constraints records",fTreeConstr->GetEntries()));
    //
  }
  else {
    //
    fTreeConstr = new TTree("AliMillePedeRecords_Constraints","Constraints Records for AliMillePede2");
    fTreeConstr->Branch("Record_Constraints","AliMillePedeRecord",&fRecord,32000,99);
  }
  fCurrRecConstrID = -1;
  //
  return kTRUE;
}

//_____________________________________________________________________________________________
void AliMillePede2::CloseDataRecStorage()
{
  if (fTreeData) {
    if (fDataRecFile->IsWritable()) {
      fDataRecFile->cd();
      fTreeData->Write();
    }
    delete fTreeData;  
    fTreeData = 0;
    fDataRecFile->Close();
    delete fDataRecFile;
    fDataRecFile = 0;
  }
  fRecFileStatus = 0;
  //
}

//_____________________________________________________________________________________________
void AliMillePede2::CloseConsRecStorage()
{
  if (fTreeConstr) {
    if (fConsRecFile->IsWritable()) {
      fConsRecFile->cd();
      fTreeConstr->Write();
    }
    delete fTreeConstr;  
    fTreeConstr = 0;
    fConsRecFile->Close();
    delete fConsRecFile;
    fConsRecFile = 0;
  }
  //
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::ReadNextRecordData()
{
  // read next data record (if any)
  if (!fTreeData || ++fCurrRecDataID >= fTreeData->GetEntries()) { fCurrRecDataID--; return kFALSE;}
  fTreeData->GetEntry(fCurrRecDataID);
  return kTRUE;
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::ReadNextRecordConstraint()
{
  // read next constraint record (if any)
  if (!fTreeConstr || ++fCurrRecConstrID >= fTreeConstr->GetEntries()) { fCurrRecConstrID--; return kFALSE;}
  fTreeConstr->GetEntry(fCurrRecConstrID);
  return kTRUE;
}

//_____________________________________________________________________________________________
void AliMillePede2::SetRecordWeight(double wgh)
{
  if (fRecFileStatus<2) InitDataRecStorage(); // create a buffer to store the data
  fRecord->SetWeight(wgh);
}

//_____________________________________________________________________________________________
void AliMillePede2::SetRecordRun(Int_t run)
{
  if (fRecFileStatus<2) InitDataRecStorage(); // create a buffer to store the data
  fRecord->SetRunID(run);
}

//_____________________________________________________________________________________________
void AliMillePede2::SetLocalEquation(double *dergb, double *derlc, double lMeas, double lSigma)
{
  if (fRecFileStatus<2) InitDataRecStorage(); // create a buffer to store the data
  //
  // write data of single measurement
  if (lSigma<=0.0) { // If parameter is fixed, then no equation
    for (int i=fNLocPar; i--;) derlc[i] = 0.0;
    for (int i=fNGloPar; i--;) dergb[i] = 0.0;
    return;
  }
  //
  fRecord->AddResidual(lMeas);
  //
  // Retrieve local param interesting indices
  for (int i=0;i<fNLocPar;i++) if (!IsZero(derlc[i])) {fRecord->AddIndexValue(i,derlc[i]); derlc[i] = 0.0;}
  //
  fRecord->AddWeight( 1.0/lSigma/lSigma );
  //
  // Idem for global parameters
  for (int i=0;i<fNGloPar;i++) if (!IsZero(dergb[i])) {
    fRecord->AddIndexValue(i,dergb[i]); dergb[i] = 0.0;
    fRecord->MarkGroup(fParamGrID[i]);
  }
  //
}

//_____________________________________________________________________________________________
void AliMillePede2::SetLocalEquation(int *indgb, double *dergb, int ngb, int *indlc,
				     double *derlc,int nlc,double lMeas,double lSigma)
{	
  // write data of single measurement
  if (lSigma<=0.0) { // If parameter is fixed, then no equation
    for (int i=nlc;i--;)  derlc[i] = 0.0;
    for (int i=ngb;i--;)  dergb[i] = 0.0;
    return;
  }
  //
  if (fRecFileStatus<2) InitDataRecStorage(); // create a buffer to store the data
  //
  fRecord->AddResidual(lMeas);
  //
  // Retrieve local param interesting indices
  for (int i=0;i<nlc;i++) if (!IsZero(derlc[i])) {fRecord->AddIndexValue(indlc[i],derlc[i]); derlc[i]=0.; indlc[i]=0;}
  //
  fRecord->AddWeight( 1./lSigma/lSigma );
  //
  // Idem for global parameters
  for (int i=0;i<ngb;i++) if (!IsZero(dergb[i])) {fRecord->AddIndexValue(indgb[i],dergb[i]); dergb[i]=0.; indgb[i]=0;} 
  //
}


//_____________________________________________________________________________________________
void AliMillePede2::SetGlobalConstraint(double *dergb, double val, double sigma)
{	
  // Define a constraint equation.
  if (!fConsRecFile || !fConsRecFile->IsWritable()) InitConsRecStorage(); // create a buffer to store the data
  //
  fRecord->Reset();
  fRecord->AddResidual(val);
  fRecord->AddWeight(sigma);
  for (int i=0; i<fNGloPar; i++) if (!IsZero(dergb[i]))  fRecord->AddIndexValue(i,dergb[i]);
  fNGloConstraints++;
  if (IsZero(sigma)) fNLagrangeConstraints++;
  //  printf("NewConstraint:\n"); fRecord->Print(); //RRR
  SaveRecordConstraint();
  //
}

//_____________________________________________________________________________________________
void AliMillePede2::SetGlobalConstraint(const int *indgb, double *dergb, int ngb, double val,double sigma)
{	
  // Define a constraint equation.
  if (!fConsRecFile || !fConsRecFile->IsWritable()) InitConsRecStorage(); // create a buffer to store the data
  fRecord->Reset();
  fRecord->AddResidual(val);
  fRecord->AddWeight(sigma);   // dummy
  for (int i=0; i<ngb; i++) if (!IsZero(dergb[i]))  fRecord->AddIndexValue(indgb[i],dergb[i]);
  fNGloConstraints++;
  if (IsZero(sigma)) fNLagrangeConstraints++;
  SaveRecordConstraint();
  //
}

//_____________________________________________________________________________________________
Int_t AliMillePede2::LocalFit(double *localParams)
{
  /*
    Perform local parameters fit once all the local equations have been set
    -----------------------------------------------------------
    localParams = (if !=0) will contain the fitted track parameters and
    related errors
  */
  static int     nrefSize = 0;
  //  static TArrayI refLoc,refGlo,nrefLoc,nrefGlo;
  static Int_t  *refLoc=0,*refGlo=0,*nrefLoc=0,*nrefGlo=0;
  int nPoints = 0;
  //
  AliSymMatrix    &matCLoc    = *fMatCLoc;
  AliMatrixSq     &matCGlo    = *fMatCGlo;
  AliRectMatrix   &matCGloLoc = *fMatCGloLoc;
  //
  memset(fVecBLoc,0,fNLocPar*sizeof(double));
  matCLoc.Reset();	
  //
  int cnt=0;
  int recSz = fRecord->GetSize();
  //
  while(cnt<recSz) {  // Transfer the measurement records to matrices
    //
    // extract addresses of residual, weight and pointers on local and global derivatives for each point
    if (nrefSize<=nPoints) {
      int *tmpA = 0;
      nrefSize = 2*(nPoints+1); 
      tmpA = refLoc;  refLoc  = new Int_t[nrefSize]; if (tmpA) memcpy(refLoc,tmpA,nPoints*sizeof(int));
      tmpA = refGlo;  refGlo  = new Int_t[nrefSize]; if (tmpA) memcpy(refGlo,tmpA,nPoints*sizeof(int));
      tmpA = nrefLoc; nrefLoc = new Int_t[nrefSize]; if (tmpA) memcpy(nrefLoc,tmpA,nPoints*sizeof(int));
      tmpA = nrefGlo; nrefGlo = new Int_t[nrefSize]; if (tmpA) memcpy(nrefGlo,tmpA,nPoints*sizeof(int));
    }
    //
    refLoc[nPoints]  = ++cnt;
    int nLoc = 0;
    while(!fRecord->IsWeight(cnt)) {nLoc++; cnt++;}
    nrefLoc[nPoints] = nLoc;
    //
    refGlo[nPoints]  = ++cnt;
    int nGlo = 0;
    while(!fRecord->IsResidual(cnt) && cnt<recSz) {nGlo++; cnt++;} 
    nrefGlo[nPoints] = nGlo;
    //
    nPoints++;
  }
  double vl;
  //
  double gloWgh = fRecord->GetWeight(); // global weight for this set
  Int_t maxLocUsed = 0;
  //
  for (int ip=nPoints;ip--;) {  // Transfer the measurement records to matrices
    double  resid  = fRecord->GetValue( refLoc[ip]-1 );
    double  weight = fRecord->GetValue( refGlo[ip]-1 )*gloWgh;    
    double *derLoc = fRecord->GetValue()+refLoc[ip];
    double *derGlo = fRecord->GetValue()+refGlo[ip];
    int    *indLoc = fRecord->GetIndex()+refLoc[ip];
    int    *indGlo = fRecord->GetIndex()+refGlo[ip];
    //
    for (int i=nrefGlo[ip];i--;) {      // suppress the global part (only relevant with iterations)
      int iID = indGlo[i];              // Global param indice
      if (fSigmaPar[iID]<=0.) continue;                                    // fixed parameter RRRCheck
      if (fIsLinear[iID]) resid -= derGlo[i]*(fInitPar[iID]+fDeltaPar[iID]);  // linear parameter
      else                resid -= derGlo[i]*fDeltaPar[iID];                  // nonlinear parameter
    }
    //
    // Symmetric matrix, don't bother j>i coeffs
    for (int i=nrefLoc[ip];i--;) {         // Fill local matrix and vector
      fVecBLoc[ indLoc[i] ] += weight*resid*derLoc[i];
      if (indLoc[i]>maxLocUsed) maxLocUsed = indLoc[i];  
      for (int j=i+1;j--;) matCLoc(indLoc[i] ,indLoc[j]) += weight*derLoc[i]*derLoc[j];
    }
    //
  } // end of the transfer of the measurement record to matrices
  //
  matCLoc.SetSizeUsed(++maxLocUsed);   // data with B=0 may use less than declared nLocals 
  //
  // first try to solve by faster Cholesky decomposition, then by Gaussian elimination
  if (!matCLoc.SolveChol(fVecBLoc,kTRUE)) {
    AliInfo("Failed to solve locals by Cholesky, trying Gaussian Elimination");
    if (!matCLoc.SolveSpmInv(fVecBLoc,kTRUE)) { 
      AliInfo("Failed to solve locals by Gaussian Elimination, skip...");
      matCLoc.Print("d");
      return 0; // failed to solve
    }
  }
  //
  // If requested, store the track params and errors
  //  printf("locfit: "); for (int i=0;i<fNLocPar;i++) printf("%+e |",fVecBLoc[i]); printf("\n");
  if (localParams) for (int i=maxLocUsed; i--;) {
      localParams[2*i]   = fVecBLoc[i];
      localParams[2*i+1] = TMath::Sqrt(TMath::Abs(matCLoc.QueryDiag(i)));
    }
  //
  float lChi2 = 0;
  int   nEq   = 0;
  //
  for (int ip=nPoints;ip--;) {  // Calculate residuals
    double  resid  = fRecord->GetValue( refLoc[ip]-1 );
    double  weight = fRecord->GetValue( refGlo[ip]-1 )*gloWgh;    
    double *derLoc = fRecord->GetValue()+refLoc[ip];
    double *derGlo = fRecord->GetValue()+refGlo[ip];
    int    *indLoc = fRecord->GetIndex()+refLoc[ip];
    int    *indGlo = fRecord->GetIndex()+refGlo[ip];
    //
    // Suppress local and global contribution in residuals;
    for (int i=nrefLoc[ip];i--;) resid -= derLoc[i]*fVecBLoc[ indLoc[i] ];     // local part 
    //
    for (int i=nrefGlo[ip];i--;) {                                             // global part
      int iID = indGlo[i];
      if ( fSigmaPar[iID] <= 0.) continue;                                     // fixed parameter RRRCheck
      if (fIsLinear[iID]) resid -= derGlo[i]*(fInitPar[iID]+fDeltaPar[iID]);   // linear parameter
      else                resid -= derGlo[i]*fDeltaPar[iID];                   // nonlinear parameter
    }
    //
    // reject the track if the residual is too large (outlier)
    double absres = TMath::Abs(resid);
    if ( (absres >= fResCutInit && fIter ==1 ) ||
	 (absres >= fResCut     && fIter > 1)) {
      if (fLocFitAdd) fNLocFitsRejected++;      
      //      printf("reject res %5ld %+e\n",fCurrRecDataID,resid);
      return 0;
    }
    // 
    lChi2 += weight*resid*resid ; // total chi^2
    nEq++;                        // number of equations			
  } // end of Calculate residuals
  //
  lChi2 /= gloWgh;  
  int nDoF = nEq-maxLocUsed;
  lChi2 = (nDoF>0) ? lChi2/nDoF : 0;  // Chi^2/dof  
  //
  if (fNStdDev != 0 && nDoF>0 && lChi2 > Chi2DoFLim(fNStdDev,nDoF)*fChi2CutFactor) { // check final chi2
    if (fLocFitAdd) fNLocFitsRejected++;      
    //    printf("reject chi2 %5ld: %+e\n",fCurrRecDataID, lChi2);
    return 0;
  }
  //
  if (fLocFitAdd) {
    fNLocFits++;
    fNLocEquations += nEq;
  }
  else {
    fNLocFits--;
    fNLocEquations -= nEq;
  }
  //
  //  local operations are finished, track is accepted 
  //  We now update the global parameters (other matrices)
  //
  int nGloInFit = 0;
  //
  for (int ip=nPoints;ip--;) {  // Update matrices
    double  resid  = fRecord->GetValue( refLoc[ip]-1 );
    double  weight = fRecord->GetValue( refGlo[ip]-1 )*gloWgh;    
    double *derLoc = fRecord->GetValue()+refLoc[ip];
    double *derGlo = fRecord->GetValue()+refGlo[ip];
    int    *indLoc = fRecord->GetIndex()+refLoc[ip];
    int    *indGlo = fRecord->GetIndex()+refGlo[ip];
    //
    for (int i=nrefGlo[ip];i--;) {    // suppress the global part
      int iID = indGlo[i];           // Global param indice
      if ( fSigmaPar[iID] <= 0.) continue;                                         // fixed parameter RRRCheck
      if (fIsLinear[iID])  resid -= derGlo[i]*(fInitPar[iID]+fDeltaPar[iID]);      // linear parameter
      else                 resid -= derGlo[i]*fDeltaPar[iID];                      // nonlinear parameter
    }
    //
    for (int ig=nrefGlo[ip];ig--;) {
      int iIDg = indGlo[ig];   // Global param indice (the matrix line)          
      if ( fSigmaPar[iIDg] <= 0.) continue;                                    // fixed parameter RRRCheck
      if (fLocFitAdd) fVecBGlo[ iIDg ] += weight*resid*derGlo[ig]; //!!!
      else            fVecBGlo[ iIDg ] -= weight*resid*derGlo[ig]; //!!!
      //
      // First of all, the global/global terms (exactly like local matrix)
      int nfill = 0;
      for (int jg=ig+1;jg--;) {	// matCGlo is symmetric by construction  
	int jIDg = indGlo[jg];
	if ( fSigmaPar[jIDg] <= 0.) continue;                                    // fixed parameter RRRCheck
	if ( !IsZero(vl = weight*derGlo[ig]*derGlo[jg]) ) {
	  fFillIndex[nfill]   = jIDg;
	  fFillValue[nfill++] = fLocFitAdd ? vl:-vl;
	}
      }
      if (nfill) matCGlo.AddToRow(iIDg,fFillValue,fFillIndex,nfill);
      //
      // Now we have also rectangular matrices containing global/local terms.
      int iCIDg = fGlo2CGlo[iIDg];  // compressed Index of index          
      if (iCIDg == -1) {
	Double_t *rowGL = matCGloLoc(nGloInFit);
	for (int k=maxLocUsed;k--;) rowGL[k] = 0.0;  // reset the row
	iCIDg = fGlo2CGlo[iIDg] = nGloInFit;
	fCGlo2Glo[nGloInFit++] = iIDg;
      }
      //
      Double_t *rowGLIDg = matCGloLoc(iCIDg);
      for (int il=nrefLoc[ip];il--;) rowGLIDg[ indLoc[il] ] += weight*derGlo[ig]*derLoc[il];
      fProcPnt[iIDg] += fLocFitAdd ? 1:-1;       // update counter
      //
    }
  } // end of Update matrices
  //
  // calculate fMatCGlo -= fMatCGloLoc * fMatCLoc * fMatCGloLoc^T
  // and       fVecBGlo -= fMatCGloLoc * fVecBLoc
  //
  //-------------------------------------------------------------- >>>
  double vll;
  for (int iCIDg=0; iCIDg<nGloInFit; iCIDg++) {
    int iIDg = fCGlo2Glo[iCIDg];
    //
    vl = 0;
    Double_t *rowGLIDg =  matCGloLoc(iCIDg);
    for (int kl=0;kl<maxLocUsed;kl++) if (rowGLIDg[kl]) vl += rowGLIDg[kl]*fVecBLoc[kl];
    if  (!IsZero(vl)) fVecBGlo[iIDg] -= fLocFitAdd ? vl : -vl;
    //
    int nfill = 0;
    for (int jCIDg=0;jCIDg<=iCIDg; jCIDg++) {  
      int jIDg = fCGlo2Glo[jCIDg];
      //
      vl = 0;
      Double_t *rowGLJDg =  matCGloLoc(jCIDg);
      for (int kl=0;kl<maxLocUsed;kl++) {
	// diag terms
	if ( (!IsZero(vll=rowGLIDg[kl]*rowGLJDg[kl])) ) vl += matCLoc.QueryDiag(kl)*vll;
	//
	// off-diag terms
	for (int ll=0;ll<kl;ll++) {
	  if ( !IsZero(vll=rowGLIDg[kl]*rowGLJDg[ll]) ) vl += matCLoc(kl,ll)*vll;
	  if ( !IsZero(vll=rowGLIDg[ll]*rowGLJDg[kl]) ) vl += matCLoc(kl,ll)*vll;
	}
      }
      if (!IsZero(vl)) {
	fFillIndex[nfill]   = jIDg;
	fFillValue[nfill++] = fLocFitAdd ? -vl : vl;
      }
    }
    if (nfill) 	matCGlo.AddToRow(iIDg,fFillValue,fFillIndex,nfill);
  }
  //
  // reset compressed index array
  //
  for (int i=nGloInFit;i--;) { 
    fGlo2CGlo[ fCGlo2Glo[i] ] = -1;
    fCGlo2Glo[i] = -1;
  }
  //
  //---------------------------------------------------- <<<
  return 1;
}

//_____________________________________________________________________________________________
Int_t AliMillePede2::GlobalFit(Double_t *par, Double_t *error, Double_t *pull)
{
  // performs a requested number of global iterations
  fIter = 1;
  //
  TStopwatch sw; sw.Start();
  //
  int res = 0;
  AliInfo("Starting Global fit.");
  while (fIter<=fMaxIter) {
    //
    res = GlobalFitIteration();
    if (!res) break;
    //
    if (!IsZero(fChi2CutFactor-fChi2CutRef)) {    
      fChi2CutFactor = TMath::Sqrt(fChi2CutFactor);
      if (fChi2CutFactor < 1.2*fChi2CutRef) {
	fChi2CutFactor = fChi2CutRef;
	//RRR	fIter = fMaxIter - 1;     // Last iteration
      }
    }
    fIter++;
  }
  //
  sw.Stop();
  AliInfo(Form("Global fit %s, CPU time: %.1f",res ? "Converged":"Failed",sw.CpuTime()));  
  if (!res) return 0;
  //
  if (par) for (int i=fNGloPar;i--;) par[i] = fInitPar[i]+fDeltaPar[i]; 
  //
  if (fGloSolveStatus==gkInvert) { // errors on params are available
    if (error) for (int i=fNGloPar;i--;) error[i] = fProcPnt[i]>0 ? TMath::Sqrt(TMath::Abs(fMatCGlo->QueryDiag(i))) : 0.;
    if (pull)  for (int i=fNGloPar;i--;) pull[i]  = fProcPnt[i]>0 && (fSigmaPar[i]*fSigmaPar[i]-fMatCGlo->QueryDiag(i))>0. && fSigmaPar[i]>0 
					   ? fDeltaPar[i]/TMath::Sqrt(fSigmaPar[i]*fSigmaPar[i]-fMatCGlo->QueryDiag(i)) : 0;
  }
  //
  return 1;
}

//_____________________________________________________________________________________________
Int_t AliMillePede2::GlobalFitIteration()
{
  // perform global parameters fit once all the local equations have been fitted
  //
  AliInfo(Form("Global Fit Iteration#%2d (Local Fit Chi^2 cut factor: %.2f)",fIter,fChi2CutFactor));
  //
  if (!fNGloPar || !fTreeData) {
    AliInfo("No data was stored, aborting iteration");
    return 0;
  }
  TStopwatch sw,sws; 
  sw.Start();
  sws.Stop();
  //
  if (!fConstrUsed) {
    fConstrUsed = new Bool_t[fNGloConstraints];
    memset(fConstrUsed,0,fNGloConstraints*sizeof(Bool_t));
  }
  // Reset all info specific for this step
  AliMatrixSq& matCGlo = *fMatCGlo;
  matCGlo.Reset();
  memset(fProcPnt,0,fNGloPar*sizeof(Int_t));
  //
  fNGloConstraints = fTreeConstr ? fTreeConstr->GetEntries() : 0;
  //
  // count number of Lagrange constraints: they need new row/cols to be added
  fNLagrangeConstraints = 0;
  for (int i=0; i<fNGloConstraints; i++) {
    ReadRecordConstraint(i);
    if ( IsZero(fRecord->GetValue(1)) ) fNLagrangeConstraints++; // exact constraint (no error) -> Lagrange multiplier 
  }
  //
  // if needed, readjust the size of the global vector (for matrices this is done automatically)
  if (!fVecBGlo || fNGloSize!=fNGloPar+fNLagrangeConstraints) {
    delete[] fVecBGlo;   // in case some constraint was added between the two manual iterations
    fNGloSize = fNGloPar+fNLagrangeConstraints;
    fVecBGlo = new Double_t[fNGloSize];
  }
  memset(fVecBGlo,0,fNGloSize*sizeof(double));
  //
  fNLocFits         = 0;
  fNLocFitsRejected = 0;
  fNLocEquations    = 0;
  //
  //  Process data records and build the matrices
  Long_t ndr = fTreeData->GetEntries();
  Long_t first = fSelFirst>0  ? fSelFirst : 0;
  Long_t last  = fSelLast<1   ? ndr : (fSelLast>=ndr ? ndr : fSelLast+Long_t(1));
  ndr = last - first;
  //
  AliInfo(Form("Building the Global matrix from data records %d : %d",first,last));
  if (ndr<1) return 0;
  //
  TStopwatch swt; swt.Start();
  fLocFitAdd = kTRUE;  // add contributions of matching tracks
  for (Long_t i=0;i<ndr;i++) {
    Long_t iev = i+first;
    ReadRecordData(iev);
    if (!IsRecordAcceptable()) continue;
    LocalFit();
    if ( (i%int(0.2*ndr)) == 0) printf("%.1f%% of local fits done\n", double(100.*i)/ndr);
  }
  swt.Stop();
  printf("%ld local fits done: ", ndr);
  swt.Print();
  sw.Start(kFALSE);
  //
  //
  // ---------------------- Reject parameters with low statistics ------------>>
  fNGloFix = 0;
  if (fMinPntValid>1 && fNGroupsSet) {
    //
    printf("Checking parameters with statistics < %d\n",fMinPntValid);
    TStopwatch swsup;
    swsup.Start();
    // 1) build the list of parameters to fix
    Int_t fixArrSize = 10;
    Int_t nFixedGroups = 0;
    TArrayI fixGroups(fixArrSize);
    //
    int grIDold = -2;
    int oldStart = -1;
    double oldMin = 1.e20;
    double oldMax =-1.e20;
    //
    for (int i=fNGloPar;i--;) { // // Reset row and column of fixed params and add 1/sig^2 to free ones
      int grID = fParamGrID[i];
      if (grID<0) continue; // not in the group
      //
      if (grID!=grIDold) { // starting new group
	if (grIDold>=0) { // decide if the group has enough statistics
	  if (oldMin<fMinPntValid && oldMax<2*fMinPntValid) { // suppress group
	    for (int iold=oldStart;iold>i;iold--) fProcPnt[iold] = 0;
	    Bool_t fnd = kFALSE;    // check if the group is already accounted
	    for (int j=nFixedGroups;j--;) if (fixGroups[j]==grIDold) {fnd=kTRUE; break;}
	    if (!fnd) {
	      if (nFixedGroups>=fixArrSize) {fixArrSize*=2; fixGroups.Set(fixArrSize);}
	      fixGroups[nFixedGroups++] = grIDold; // add group to fix
	    }
	  }	  
	}
	grIDold = grID; // mark the start of the new group
	oldStart = i;
	oldMin =  1.e20;
	oldMax = -1.e20;
      }
      if (oldMin>fProcPnt[i]) oldMin = fProcPnt[i];
      if (oldMax<fProcPnt[i]) oldMax = fProcPnt[i];
      //
    }
    // extra check for the last group
    if (grIDold>=0 && oldMin<fMinPntValid && oldMax<2*fMinPntValid) { // suppress group
      for (int iold=oldStart;iold--;) fProcPnt[iold] = 0;
      Bool_t fnd = kFALSE;    // check if the group is already accounted
      for (int j=nFixedGroups;j--;) if (fixGroups[j]==grIDold) {fnd=kTRUE; break;}
      if (!fnd) {
	if (nFixedGroups>=fixArrSize) {fixArrSize*=2; fixGroups.Set(fixArrSize);}
	fixGroups[nFixedGroups++] = grIDold; // add group to fix
      }
    }
    //
    // 2) loop over records and add contributions of fixed groups with negative sign
    fLocFitAdd = kFALSE;
    //
    for (Long_t i=0;i<ndr;i++) {
      Long_t iev = i+first;
      ReadRecordData(iev);
      if (!IsRecordAcceptable()) continue;
      Bool_t suppr = kFALSE;
      for (int ifx=nFixedGroups;ifx--;)if (fRecord->IsGroupPresent(fixGroups[ifx])) suppr = kTRUE;
      if (suppr) LocalFit();
    }
    fLocFitAdd = kTRUE;
    //
    if (nFixedGroups) {
      printf("Suppressed contributions of groups with NPoints<%d :\n",fMinPntValid);
      for (int i=0;i<nFixedGroups;i++) printf("%d ",fixGroups[i]); printf("\n");
    }
    swsup.Stop();
    swsup.Print();
  }
  // ---------------------- Reject parameters with low statistics ------------<<
  //
  // add large number to diagonal of fixed params  
  //
  for (int i=fNGloPar;i--;) { // // Reset row and column of fixed params and add 1/sig^2 to free ones
    //    printf("#%3d : Nproc : %5d   grp: %d\n",i,fProcPnt[i],fParamGrID[i]);
    if (fProcPnt[i]<1) {
      fNGloFix++; 
      fVecBGlo[i] = 0.;
      matCGlo.DiagElem(i) = 1.;//float(fNLocEquations*fNLocEquations);
      //      matCGlo.DiagElem(i) = float(fNLocEquations*fNLocEquations);
    }
    else matCGlo.DiagElem(i) += (fgWeightSigma ? fProcPnt[i] : 1.)/(fSigmaPar[i]*fSigmaPar[i]);
  }
  //
  for (int i=fNGloPar;i--;) fDiagCGlo[i] = matCGlo.QueryDiag(i); // save the diagonal elements  
  //
  // add constraint equations
  int nVar = fNGloPar;                    // Current size of global matrix	
  for (int i=0; i<fNGloConstraints; i++) {
    ReadRecordConstraint(i);
    double val   = fRecord->GetValue(0);  
    double sig   = fRecord->GetValue(1);  
    int    *indV = fRecord->GetIndex()+2;
    double *der  = fRecord->GetValue()+2;
    int    csize = fRecord->GetSize()-2;
    //
    // check if after suppression of fixed variables there are non-0 derivatives
    // and determine the max statistics of involved params
    int nSuppressed = 0;
    int maxStat = 1;
    for (int j=csize;j--;) {
      if (fProcPnt[indV[j]]<1) nSuppressed++; 
      else {
	maxStat = TMath::Max(maxStat,fProcPnt[indV[j]]);
      }
    }
    //
    if (nSuppressed==csize) {
      //      AliInfo(Form("Neglecting constraint %d of %d derivatives since no free parameters left",i,csize));
      //
      // was this constraint ever created ? 
      if ( sig==0 && fConstrUsed[i] ) { // this is needed only for constraints with Lagrange multiplier
	// to avoid empty row impose dummy constraint on "Lagrange multiplier"
	matCGlo.DiagElem(nVar) = 1.;
	fVecBGlo[nVar++] = 0;
      }
      continue;
    }
    //
    // account for already accumulated corrections
    for (int j=csize; j--;) val -= der[j]*(fInitPar[ indV[j] ]+fDeltaPar[ indV[j] ]);
    //
    if (sig > 0) {  // this is a gaussian constriant: no Lagrange multipliers are added
      //
      double sig2i = (fgWeightSigma ? TMath::Sqrt(maxStat) : 1.)/sig/sig;
      for (int ir=0;ir<csize;ir++) {
	int iID = indV[ir];
	for (int ic=0;ic<=ir;ic++) { // matrix is symmetric
	  int jID = indV[ic];
	  double vl = der[ir]*der[ic]*sig2i;
	  if (!IsZero(vl)) matCGlo(iID,jID) += vl;
	}
	fVecBGlo[iID] += val*der[ir]*sig2i;
      }
    }
    else {   // this is exact constriant:  Lagrange multipliers must be added
      for (int j=csize; j--;) {
	int jID = indV[j];
	if (fProcPnt[jID]<1) continue;                      // this parameter was fixed, don't put it into constraint
	matCGlo(nVar,jID) = float(fNLocEquations)*der[j];   // fMatCGlo is symmetric, only lower triangle is filled  
      }
      //
      if (matCGlo.QueryDiag(nVar)) matCGlo.DiagElem(nVar) = 0.0;
      fVecBGlo[nVar++] = float(fNLocEquations)*val; //RS ? should we use here fNLocFits ? 
      fConstrUsed[i] = kTRUE;
    }
  }
  //
  AliInfo(Form("Obtained %-7ld equations from %-7ld records (%-7ld rejected). Fixed %-4d globals",
	       fNLocEquations,fNLocFits,fNLocFitsRejected,fNGloFix));

  //
  sws.Start();
  fGloSolveStatus = SolveGlobalMatEq();                     // obtain solution for this step
  sws.Stop();
  printf("Solve %d |",fIter); sws.Print();
  //
  sw.Stop();
  AliInfo(Form("Iteration#%2d %s. CPU time: %.1f",fIter,fGloSolveStatus==gkFailed ? "Failed":"Converged",sw.CpuTime()));
  if (fGloSolveStatus==gkFailed) return 0;
  //
  for (int i=fNGloPar;i--;) fDeltaPar[i] += fVecBGlo[i];    // Update global parameters values (for iterations)
  //
  //  PrintGlobalParameters();
  return 1;
}

//_____________________________________________________________________________________________
Int_t AliMillePede2::SolveGlobalMatEq()
{
  //
  // solve global matrix equation MatCGlob*X=VecBGlo and store the result in the VecBGlo
  //
  /*
  printf("GlobalMatrix\n");
  fMatCGlo->Print();
  printf("RHS\n");
  for (int i=0;i<fNGloPar;i++) printf("%d %+e\n",i,fVecBGlo[i]);
  */
  //
  if (!fgIsMatGloSparse) {
    //
    if (fNLagrangeConstraints==0) { // pos-def systems are faster to solve by Cholesky
      if ( ((AliSymMatrix*)fMatCGlo)->SolveChol(fVecBGlo, fgInvChol) ) return fgInvChol ? gkInvert:gkNoInversion;
      else AliInfo("Solution of Global Dense System by Cholesky failed, trying Gaussian Elimiation");
    }
    //
    if (((AliSymMatrix*)fMatCGlo)->SolveSpmInv(fVecBGlo, kTRUE)) return gkInvert;
    else AliInfo("Solution of Global Dense System by Gaussian Elimination failed, trying iterative methods");
  }
  // try to solve by minres
  TVectorD sol(fNGloSize);
  //
  AliMinResSolve *slv = new AliMinResSolve(fMatCGlo,fVecBGlo);
  if (!slv) return gkFailed;
  //
  Bool_t res = kFALSE;
  if      (fgIterSol == AliMinResSolve::kSolMinRes) 
    res =  slv->SolveMinRes(sol,fgMinResCondType,fgMinResMaxIter,fgMinResTol);
  else if (fgIterSol == AliMinResSolve::kSolFGMRes) 
    res =  slv->SolveFGMRES(sol,fgMinResCondType,fgMinResMaxIter,fgMinResTol,fgNKrylovV);
  else 
    AliInfo(Form("Undefined Iteritive Solver ID=%d, only %d are defined",fgIterSol,AliMinResSolve::kNSolvers));
  //
  if (!res) {
    const char* faildump = "fgmr_failed.dat";
    int defout = dup(1);
    int slvDump = open(faildump, O_RDWR|O_CREAT, 0666);
    dup2(slvDump,1);
    //
    printf("#Failed to solve using solver %d with PreCond: %d MaxIter: %d Tol: %e NKrylov: %d\n",
	   fgIterSol,fgMinResCondType,fgMinResMaxIter,fgMinResTol,fgNKrylovV);
    printf("#Dump of matrix:\n");
    fMatCGlo->Print("10");
    printf("#Dump of RHS:\n");
    for (int i=0;i<fNGloSize;i++) printf("%d %+.10f\n",i,fVecBGlo[i]);
    //
    dup2(defout,1);
    close(slvDump);
    close(defout);
    printf("#Dumped failed matrix and RHS to %s\n",faildump);
    return gkFailed;
  }
  for (int i=fNGloSize;i--;) fVecBGlo[i] = sol[i];
  return gkNoInversion;
  //
}

//_____________________________________________________________________________________________
Float_t AliMillePede2::Chi2DoFLim(int nSig, int nDoF) const
{
  /// return the limit in chi^2/nd for n sigmas stdev authorized
  // Only n=1, 2, and 3 are expected in input
  int lNSig;
  float sn[3]        =	{0.47523, 1.690140, 2.782170};
  float table[3][30] = {{1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630,
			 1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321,
			 1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136,
			 1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040},
			{4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736,
			 1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658,
			 1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346,
			 1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742},
			{9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468,
			 2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635,
			 2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901,
			 1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681}};
  
  if (nDoF < 1) {
    return 0.0;
  }
  else {  
    lNSig = TMath::Max(1,TMath::Min(nSig,3));
    
    if (nDoF <= 30) {    
      return table[lNSig-1][nDoF-1];
    }
    else { // approximation
      return ((sn[lNSig-1]+TMath::Sqrt(float(2*nDoF-3)))*
	      (sn[lNSig-1]+TMath::Sqrt(float(2*nDoF-3))))/float(2*nDoF-2);
    }
  }
}

//_____________________________________________________________________________________________
Int_t AliMillePede2::SetIterations(double lChi2CutFac)
{
  // Number of iterations is calculated from lChi2CutFac 
  fChi2CutFactor = TMath::Max(1.0, lChi2CutFac);
  //
  AliInfo(Form("Initial cut factor is %f",fChi2CutFactor));
  fIter = 1; // Initializes the iteration process
  return 1;
}

//_____________________________________________________________________________________________
Double_t AliMillePede2::GetParError(int iPar) const
{
  // return error for parameter iPar
  if (fGloSolveStatus==gkInvert) {
    double res = fMatCGlo->QueryDiag(iPar);
    if (res>=0) return TMath::Sqrt(res);
  } 
  return -1.;
}


//_____________________________________________________________________________________________
Int_t AliMillePede2::PrintGlobalParameters() const
{
  ///  Print the final results into the logfile
  double lError = 0.;
  double lGlobalCor =0.;
	
  AliInfo("");
  AliInfo("   Result of fit for global parameters");
  AliInfo("   ===================================");
  AliInfo("    I       initial       final       differ        lastcor        error       gcor       Npnt");
  AliInfo("----------------------------------------------------------------------------------------------");
  //
  for (int i=0; i<fNGloPar; i++) {
    lError = GetParError(i);
    lGlobalCor = 0.0;
    //		
    double dg;
    if (fGloSolveStatus==gkInvert && TMath::Abs( (dg=fMatCGlo->QueryDiag(i)) *fDiagCGlo[i]) > 0) {    
      lGlobalCor = TMath::Sqrt(TMath::Abs(1.0-1.0/(dg*fDiagCGlo[i])));
      AliInfo(Form("%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %6d",
		   i,fInitPar[i],fInitPar[i]+fDeltaPar[i],fDeltaPar[i],fVecBGlo[i],lError,lGlobalCor,fProcPnt[i]));
    }
    else {    
      AliInfo(Form("%d\t %.6f\t %.6f\t %.6f\t %.6f\t OFF\t OFF\t %6d",i,fInitPar[i],fInitPar[i]+fDeltaPar[i],
		   fDeltaPar[i],fVecBGlo[i],fProcPnt[i]));
    }
  }
  return 1;
}

//_____________________________________________________________________________________________
Bool_t AliMillePede2::IsRecordAcceptable() const
{
  // validate record according run lists set by the user
  static Long_t prevRunID = kMaxInt;
  static Bool_t prevAns   = kTRUE;
  Long_t runID = fRecord->GetRunID();
  if (runID!=prevRunID) {
    int n = 0;
    prevRunID = runID;
    // is run to be rejected?
    if (fRejRunList && (n=fRejRunList->GetSize())) {
      prevAns = kTRUE;
      for (int i=n;i--;) if (runID == (*fRejRunList)[i]) {
	  prevAns = kFALSE; 
	  printf("New Run to reject: %ld -> %d\n",runID,prevAns);
	  break;
	}
    }
    else if (fAccRunList && (n=fAccRunList->GetSize())) {     // is run specifically selected
      prevAns = kFALSE;
      for (int i=n;i--;) if (runID == (*fAccRunList)[i]) {prevAns = kTRUE; break;}
    }
  }
  //
  return prevAns;
  //
}

//_____________________________________________________________________________________________
void AliMillePede2::SetRejRunList(const UInt_t *runs, Int_t nruns)
{
  // set the list of runs to be rejected
  if (fRejRunList) delete fRejRunList; 
  fRejRunList = 0;
  if (nruns<1 || !runs) return;
  fRejRunList = new TArrayL(nruns);
  for (int i=0;i<nruns;i++) (*fRejRunList)[i] = runs[i];
}

//_____________________________________________________________________________________________
void AliMillePede2::SetAccRunList(const UInt_t *runs, Int_t nruns)
{
  // set the list of runs to be selected
  if (fAccRunList) delete fAccRunList;
  fAccRunList = 0;
  if (nruns<1 || !runs) return;
  fAccRunList = new TArrayL(nruns);
  for (int i=0;i<nruns;i++) (*fAccRunList)[i] = runs[i];
}
