#include "AliMillePede2.h"
#include "AliLog.h"
#include <TStopwatch.h>

using namespace std;

ClassImp(AliMillePedeRecord)

//_____________________________________________________________________________________________
AliMillePedeRecord::AliMillePedeRecord() : 
fSize(0),fIndex(0),fValue(0) {SetBufferSize(0);}

//_____________________________________________________________________________________________
AliMillePedeRecord::AliMillePedeRecord(const AliMillePedeRecord& src) : 
  TObject(src),fSize(src.fSize),fIndex(0),fValue(0)
{
  fIndex = new int[GetBufferSize()];
  memcpy(fIndex,src.fIndex,fSize*sizeof(int));
  fValue = new double[GetBufferSize()];
  memcpy(fValue,src.fValue,fSize*sizeof(double));
}

//_____________________________________________________________________________________________
AliMillePedeRecord& AliMillePedeRecord::operator=(const AliMillePedeRecord& rhs)
{
  if (this!=&rhs) {
    Reset();
    for (int i=0;i<rhs.GetSize();i++) {
      double val;
      int    ind;
      rhs.GetIndexValue(i,ind,val);
      AddIndexValue(ind,val);
    }
  }
  return *this;
}

//_____________________________________________________________________________________________
AliMillePedeRecord::~AliMillePedeRecord() {delete[] fIndex; delete[] fValue;}

//_____________________________________________________________________________________________
void AliMillePedeRecord::Print(const Option_t *) const
{
  if (!fSize) {AliInfo("No data"); return;}
  int cnt=0,point=0;
  //  
  while(cnt<fSize) {
    //
    double resid = fValue[cnt++];
    double *derLoc = GetValue()+cnt;
    int    *indLoc = GetIndex()+cnt;
    int     nLoc = 0;
    while(!IsWeight(cnt)) {nLoc++;cnt++;}
    double weight = GetValue(cnt++);
    double *derGlo = GetValue()+cnt;
    int    *indGlo = GetIndex()+cnt;
    int     nGlo = 0;
    while(!IsResidual(cnt) && cnt<fSize) {nGlo++; cnt++;} 
    //
    printf("\n*** Point#%2d | Residual = %+.4e | Weight = %+.4e\n",point++,resid,weight);
    printf("Locals : "); 
    for (int i=0;i<nLoc;i++) printf("[%5d] %+.4e|",indLoc[i],derLoc[i]); printf("\n");
    printf("Globals: "); 
    for (int i=0;i<nGlo;i++) printf("[%5d] %+.4e|",indGlo[i],derGlo[i]); printf("\n");
    //
  }
  //
}

//_____________________________________________________________________________________________
void AliMillePedeRecord::Expand(int bfsize)
{
  // add extra space 
  bfsize = TMath::Max(bfsize,GetBufferSize());
  int *tmpI = new int[bfsize];
  memcpy(tmpI,fIndex, fSize*sizeof(int));
  delete fIndex;
  fIndex = tmpI;
  //
  double *tmpD = new double[bfsize];
  memcpy(tmpD,fValue, fSize*sizeof(double));
  delete fValue;
  fValue = tmpD;
  //
  SetBufferSize(bfsize);
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
ClassImp(AliMillePede2)

Bool_t   AliMillePede2::fgIsMatGloSparse = kFALSE;   // use faster dense matrix by default
Int_t    AliMillePede2::fgMinResCondType = kTRUE;        // No preconditioner by default
Double_t AliMillePede2::fgMinResTol      = 1.e-8;   // default tolerance
Int_t    AliMillePede2::fgMinResMaxIter  = 3000;     // default max number of iterations 
Int_t    AliMillePede2::fgIterSol        = AliMinResSolve::kSolMinRes; // default iterative solver
Int_t    AliMillePede2::fgNKrylovV       = 30;       // default number of Krylov vectors to keep

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
  fNLocFits(0),
  fNLocFitsRejected(0),
  fNGloFix(0),
  fGloSolveStatus(gkFailed),
//
  fChi2CutFactor(1.),
  fChi2CutRef(1.),
  fResCutInit(100.),
  fResCut(100.),
  fMinPntValid(0),
//
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
  fDataRecFName("/tmp/mp2_data_records.root"),
  fRecord(0),
  fDataRecFile(0),  
  fTreeData(0),
  //
  fConstrRecFName("/tmp/mp2_constraints_records.root"),
  fTreeConstr(0),
  fConsRecFile(0),
  fCurrRecDataID(0),
  fCurrRecConstrID(0)
{}

//_____________________________________________________________________________________________
AliMillePede2::AliMillePede2(const AliMillePede2& src) : 
  TObject(src),fNLocPar(0),fNGloPar(0),fNGloSize(0),fNLocEquations(0),fIter(0),
  fMaxIter(10),fNStdDev(3),fNGloConstraints(0),fNLocFits(0),fNLocFitsRejected(0),
  fNGloFix(0),fGloSolveStatus(0),fChi2CutFactor(0),fChi2CutRef(0),fResCutInit(0),
  fResCut(0),fMinPntValid(0),fProcPnt(0),fVecBLoc(0),fDiagCGlo(0),fVecBGlo(0),
  fInitPar(0),fDeltaPar(0),fSigmaPar(0),fIsLinear(0),fConstrUsed(0),fGlo2CGlo(0),fCGlo2Glo(0),
  fMatCLoc(0),fMatCGlo(0),fMatCGloLoc(0),fDataRecFName(0),fRecord(0),fDataRecFile(0),
  fTreeData(0),fConstrRecFName(0),fTreeConstr(0),fConsRecFile(0),fCurrRecDataID(0),fCurrRecConstrID(0)
{printf("Dummy\n");}

//_____________________________________________________________________________________________
AliMillePede2::~AliMillePede2() 
{
  CloseDataRecStorage();
  CloseConsRecStorage();
  //
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
  //
  delete fRecord;
  delete fMatCLoc;
  delete fMatCGlo;
  delete fMatCGloLoc;
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
    fMatCLoc      = new AliSymMatrix(fNLocPar);
    fMatCGloLoc   = new TMatrixD(fNGloPar,fNLocPar);
    //
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
    fGlo2CGlo[i]=fCGlo2Glo[i] = -1;
    fIsLinear[i] = kTRUE;
  }
  //
  return 1;
}

//_____________________________________________________________________________________________
void AliMillePede2::ResetConstraints()
{
  // reset constraints tree. Allows to redefine the constraints for preprocessed data 
  CloseConsRecStorage();
  InitConsRecStorage(kFALSE);
  fNGloConstraints = 0;
  AliInfo("Constraints are reset");
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
  if (!fDataRecFile) {AliInfo(Form("Failed to initialize data records file %s",GetDataRecFName())); return kFALSE;}
  //
  AliInfo(Form("File %s used for derivatives records",GetDataRecFName()));
  if (read) {
    fTreeData = (TTree*)fDataRecFile->Get("AliMillePedeRecords_Data");
    if (!fTreeData) {AliInfo(Form("Did not find data records tree in %s",GetDataRecFName())); return kFALSE;}
    fTreeData->SetBranchAddress("Record_Data",&fRecord);
    AliInfo(Form("Found %d derivatives records",fTreeData->GetEntries()));
  }
  else {
    fTreeData = new TTree("AliMillePedeRecords_Data","Data Records for AliMillePede2");
    fTreeData->Branch("Record_Data","AliMillePedeRecord",&fRecord,32000,99);
  }
  fCurrRecDataID = -1;
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
void AliMillePede2::SetLocalEquation(double *dergb, double *derlc, double lMeas, double lSigma)
{
  if (!fDataRecFile || !fDataRecFile->IsWritable()) InitDataRecStorage(); // create a buffer to store the data
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
  for (int i=0;i<fNLocPar;i++) if (derlc[i]!=0.0) {fRecord->AddIndexValue(i,derlc[i]); derlc[i] = 0.0;}
  //
  fRecord->AddWeight( 1.0/lSigma/lSigma );
  //
  // Idem for global parameters
  for (int i=0;i<fNGloPar;i++) if (dergb[i]!=0.0) {fRecord->AddIndexValue(i,dergb[i]); dergb[i] = 0.0;}
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
  if (!fDataRecFile || !fDataRecFile->IsWritable()) InitDataRecStorage(); // create a buffer to store the data
  //
  fRecord->AddResidual(lMeas);
  //
  // Retrieve local param interesting indices
  for (int i=0;i<nlc;i++) if (derlc[i]!=0.0) {fRecord->AddIndexValue(indlc[i],derlc[i]); derlc[i]=0.; indlc[i]=0;}
  //
  fRecord->AddWeight( 1./lSigma/lSigma );
  //
  // Idem for global parameters
  for (int i=0;i<ngb;i++) if (dergb[i]!=0.0) {fRecord->AddIndexValue(indgb[i],dergb[i]); dergb[i]=0.; indgb[i]=0;} 
  //
}

//_____________________________________________________________________________________________
void AliMillePede2::SetGlobalConstraint(double *dergb, double val)
{	
  // Define a constraint equation.
  if (!fConsRecFile || !fConsRecFile->IsWritable()) InitConsRecStorage(); // create a buffer to store the data
  //
  fRecord->Reset();
  fRecord->AddResidual(val);
  fRecord->AddWeight(val);   // dummy
  for (int i=0; i<fNGloPar; i++) if (dergb[i]!=0)  fRecord->AddIndexValue(i,dergb[i]);
  fNGloConstraints++;
  SaveRecordConstraint();
  //
}

//_____________________________________________________________________________________________
void AliMillePede2::SetGlobalConstraint(const int *indgb, double *dergb, int ngb, double val)
{	
  // Define a constraint equation.
  if (!fConsRecFile || !fConsRecFile->IsWritable()) InitConsRecStorage(); // create a buffer to store the data
  fRecord->Reset();
  fRecord->AddResidual(val);
  fRecord->AddWeight(val);   // dummy
  for (int i=0; i<ngb; i++) if (dergb[i]!=0)  fRecord->AddIndexValue(indgb[i],dergb[i]);
  fNGloConstraints++;
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
  static int     nRefSize = 0;
  static TArrayI RefLoc,RefGlo,nRefLoc,nRefGlo;
  int nPoints = 0;
  //
  AliSymMatrix    &MatCLoc    = *fMatCLoc;
  AliMatrixSq     &MatCGlo    = *fMatCGlo;
  TMatrixD        &MatCGloLoc = *fMatCGloLoc;
  //
  memset(fVecBLoc,0,fNLocPar*sizeof(double));
  MatCLoc.Reset();	
  //
  int cnt=0;
  int recSz = fRecord->GetSize();
  //
  while(cnt<recSz) {  // Transfer the measurement records to matrices
    //
    // extract addresses of residual, weight and pointers on local and global derivatives for each point
    if (nRefSize<=nPoints) {
      nRefSize = 2*(nPoints+1); 
      RefLoc.Set(nRefSize); 
      RefGlo.Set(nRefSize);
      nRefLoc.Set(nRefSize); 
      nRefGlo.Set(nRefSize);
    }
    //
    RefLoc[nPoints]  = ++cnt;
    int nLoc = 0;
    while(!fRecord->IsWeight(cnt)) {nLoc++; cnt++;}
    nRefLoc[nPoints] = nLoc;
    //
    RefGlo[nPoints]  = ++cnt;
    int nGlo = 0;
    while(!fRecord->IsResidual(cnt) && cnt<recSz) {nGlo++; cnt++;} 
    nRefGlo[nPoints] = nGlo;
    //
    nPoints++;
  }
  double vl;
  //
  for (int ip=nPoints;ip--;) {  // Transfer the measurement records to matrices
    double  resid  = fRecord->GetValue( RefLoc[ip]-1 );
    double  weight = fRecord->GetValue( RefGlo[ip]-1 );    
    double *derLoc = fRecord->GetValue()+RefLoc[ip];
    double *derGlo = fRecord->GetValue()+RefGlo[ip];
    int    *indLoc = fRecord->GetIndex()+RefLoc[ip];
    int    *indGlo = fRecord->GetIndex()+RefGlo[ip];
    //
    for (int i=nRefGlo[ip];i--;) {       // suppress the global part (only relevant with iterations)
      int iID = indGlo[i];              // Global param indice
      if ( fSigmaPar[iID] <= 0.) continue;                                    // fixed parameter RRRCheck
      if (fIsLinear[iID]) resid -= derGlo[i]*(fInitPar[iID]+fDeltaPar[iID]);  // linear parameter
      else                resid -= derGlo[i]*fDeltaPar[iID];                  // nonlinear parameter
    }
    //
    for (int i=nRefLoc[ip];i--;) {         // Fill local matrix and vector
      int iID = indLoc[i];              // Local param indice (the matrix line) 
      if ( (vl=weight*resid*derLoc[i])!=0 ) fVecBLoc[iID] += vl;
      //			
      for (int j=i+1;j--;) {           // Symmetric matrix, don't bother j>i coeffs
	int jID = indLoc[j];	
	if ( (vl=weight*derLoc[i]*derLoc[j])!=0 ) MatCLoc(iID,jID) += vl;
      }
    }   
    //
  } // end of the transfer of the measurement record to matrices
  //
  // first try to solve by faster Cholesky decomposition, then by Gaussian elimination
  if (!MatCLoc.SolveChol(fVecBLoc,kTRUE)) {
    AliInfo("Failed to solve locals by Cholesky, traying Gaussian Elimination");
    if (!MatCLoc.SolveSpmInv(fVecBLoc,kTRUE)) { 
      AliInfo("Failed to solve locals by Gaussian Elimination, skip...");
      MatCLoc.Print("d");
      return 0; // failed to solve
    }
  }
  //
  // If requested, store the track params and errors
  if (localParams) for (int i=fNLocPar; i--;) {
      localParams[2*i]   = fVecBLoc[i];
      localParams[2*i+1] = TMath::Sqrt(TMath::Abs(MatCLoc.QuerryDiag(i)));
    }
  //
  float lChi2 = 0;
  int   nEq   = 0;
  //
  for (int ip=nPoints;ip--;) {  // Calculate residuals
    double  resid  = fRecord->GetValue( RefLoc[ip]-1 );
    double  weight = fRecord->GetValue( RefGlo[ip]-1 );    
    double *derLoc = fRecord->GetValue()+RefLoc[ip];
    double *derGlo = fRecord->GetValue()+RefGlo[ip];
    int    *indLoc = fRecord->GetIndex()+RefLoc[ip];
    int    *indGlo = fRecord->GetIndex()+RefGlo[ip];
    //
    // Suppress local and global contribution in residuals;
    for (int i=nRefLoc[ip];i--;) resid -= derLoc[i]*fVecBLoc[ indLoc[i] ];     // local part 
    //
    for (int i=nRefGlo[ip];i--;) {                                             // global part
      int iID = indGlo[i];
      if ( fSigmaPar[iID] <= 0.) continue;                                     // fixed parameter RRRCheck
      if (fIsLinear[iID]) resid -= derGlo[i]*(fInitPar[iID]+fDeltaPar[iID]);   // linear parameter
      else                resid -= derGlo[i]*fDeltaPar[iID];                   // nonlinear parameter
    }
    //
    // reject the track if the residual is too large (outlier)
    if ( (TMath::Abs(resid) >= fResCutInit && fIter ==1 ) ||
	 (TMath::Abs(resid) >= fResCut     && fIter > 1)) {
      fNLocFitsRejected++;      
      return 0;
    }
    // 
    lChi2 += weight*resid*resid ; // total chi^2
    nEq++;                        // number of equations			
  } // end of Calculate residuals
  //
  //
  int nDoF = nEq-fNLocPar;
  lChi2 = (nDoF>0) ? lChi2/nDoF : 0;  // Chi^2/dof  
  //
  if (fNStdDev != 0 && nDoF>0 && lChi2 > Chi2DoFLim(fNStdDev,nDoF)*fChi2CutFactor) { // check final chi2
    fNLocFitsRejected++;      
    return 0;
  }
  //
  fNLocFits++;
  fNLocEquations += nEq;
  //
  //  local operations are finished, track is accepted 
  //  We now update the global parameters (other matrices)
  //
  int nGloInFit = 0;
  //
  for (int ip=nPoints;ip--;) {  // Update matrices
    double  resid  = fRecord->GetValue( RefLoc[ip]-1 );
    double  weight = fRecord->GetValue( RefGlo[ip]-1 );    
    double *derLoc = fRecord->GetValue()+RefLoc[ip];
    double *derGlo = fRecord->GetValue()+RefGlo[ip];
    int    *indLoc = fRecord->GetIndex()+RefLoc[ip];
    int    *indGlo = fRecord->GetIndex()+RefGlo[ip];
    //
    for (int i=nRefGlo[ip];i--;) {    // suppress the global part
      int iID = indGlo[i];           // Global param indice
      if ( fSigmaPar[iID] <= 0.) continue;                                         // fixed parameter RRRCheck
      if (fIsLinear[iID] == 0)  resid -= derGlo[i]*(fInitPar[iID]+fDeltaPar[iID]); // linear parameter
      else                      resid -= derGlo[i]*fDeltaPar[iID];                 // nonlinear parameter
    }
    //
    for (int ig=nRefGlo[ip];ig--;) {
      int iIDg = indGlo[ig];   // Global param indice (the matrix line)          
      if ( fSigmaPar[iIDg] <= 0.) continue;                                    // fixed parameter RRRCheck
      if ( (vl = weight*resid*derGlo[ig])!=0 ) fVecBGlo[ iIDg ] += vl;
      //
      // First of all, the global/global terms (exactly like local matrix)
      for (int jg=ig+1;jg--;) {	// MatCGlo is symmetric by construction  
	int jIDg = indGlo[jg];
	if ( fSigmaPar[jIDg] <= 0.) continue;                                    // fixed parameter RRRCheck
	if ( (vl = weight*derGlo[ig]*derGlo[jg])!=0 ) MatCGlo(iIDg,jIDg) += vl; 
      } 
      //
      // Now we have also rectangular matrices containing global/local terms.
      int iCIDg = fGlo2CGlo[iIDg];  // compressed Index of index          
      if (iCIDg == -1) {
	for (int k=fNLocPar;k--;) MatCGloLoc(nGloInFit,k) = 0.0;  // reset the row
	iCIDg = fGlo2CGlo[iIDg] = nGloInFit;
	fCGlo2Glo[nGloInFit++] = iIDg;
      }
      //
      for (int il=nRefLoc[ip];il--;) {
	int iIDl = indLoc[il];
	if ( (vl = weight*derGlo[ig]*derLoc[il])!=0) MatCGloLoc(iCIDg,iIDl) += vl;
      }
      // update counter
      fProcPnt[iIDg]++;
      //
    }
  } // end of Update matrices
  //
  // calculate fMatCGlo -= fMatCGloLoc * fMatCLoc * fMatCGloLoc^T
  // and       fVecBGlo -= fMatCGloLoc * fVecBLoc
  //
  double vll;
  for (int iCIDg=0; iCIDg<nGloInFit; iCIDg++) {
    int iIDg = fCGlo2Glo[iCIDg];
    //
    vl = 0;
    for (int kl=0;kl<fNLocPar;kl++) if ( (vll = MatCGloLoc(iCIDg,kl))!=0) vl += vll*fVecBLoc[kl];
    if  (vl!=0) fVecBGlo[iIDg] -= vl;
    //
    for (int jCIDg=0;jCIDg<=iCIDg; jCIDg++) {  
      int jIDg = fCGlo2Glo[jCIDg];
      //
      vl = 0;
      for (int kl=0;kl<fNLocPar;kl++) {
	// diag terms
	if ( (vll=MatCGloLoc(iCIDg,kl)*MatCGloLoc(jCIDg,kl))!=0 ) vl += MatCLoc.QuerryDiag(kl)*vll;
	//
	// off-diag terms
	for (int ll=0;ll<kl;ll++) {
	  if ( (vll=MatCGloLoc(iCIDg,kl)*MatCGloLoc(jCIDg,ll))!=0 ) vl += MatCLoc(kl,ll)*vll;
	  if ( (vll=MatCGloLoc(iCIDg,ll)*MatCGloLoc(jCIDg,kl))!=0 ) vl += MatCLoc(kl,ll)*vll;
	}
      }
      if (vl!=0) MatCGlo(iIDg,jIDg) -= vl; 
      //
    }
  }
  //
  // reset compressed index array
  //
  for (int i=nGloInFit;i--;) { fGlo2CGlo[ fCGlo2Glo[i] ] = -1; fCGlo2Glo[i] = -1;}
  //
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
    if (fChi2CutFactor != fChi2CutRef) {    
      fChi2CutFactor = TMath::Sqrt(fChi2CutFactor);
      if (fChi2CutFactor < 1.2*fChi2CutRef) {
	fChi2CutFactor = fChi2CutRef;
	fIter = fMaxIter - 1;     // Last iteration
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
    if (error) for (int i=fNGloPar;i--;) error[i] = fProcPnt[i]>0 ? TMath::Sqrt(TMath::Abs(fMatCGlo->QuerryDiag(i))) : 0.;
    if (pull)  for (int i=fNGloPar;i--;) pull[i]  = fProcPnt[i]>0 && (fSigmaPar[i]*fSigmaPar[i]-fMatCGlo->QuerryDiag(i))>0. && fSigmaPar[i]>0 
					   ? fDeltaPar[i]/TMath::Sqrt(fSigmaPar[i]*fSigmaPar[i]-fMatCGlo->QuerryDiag(i)) : 0;
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
  TStopwatch sw; sw.Start();
  //
  if (!fConstrUsed) {
    fConstrUsed = new Bool_t[fNGloConstraints];
    memset(fConstrUsed,0,fNGloConstraints*sizeof(Bool_t));
  }
  // Reset all info specific for this step
  AliMatrixSq& MatCGlo = *fMatCGlo;
  MatCGlo.Reset();
  memset(fProcPnt,0,fNGloPar*sizeof(Int_t));
  //
  fNGloConstraints = fTreeConstr ? fTreeConstr->GetEntries() : 0;
  //
  // if needed, readjust the size of the global vector (for matrices this is done automatically)
  if (!fVecBGlo || fNGloSize!=fNGloPar+fNGloConstraints) {
    delete[] fVecBGlo;   // in case some constraint was added between the two manual iterations
    fNGloSize = fNGloPar+fNGloConstraints;
    fVecBGlo = new Double_t[fNGloSize];
  }
  memset(fVecBGlo,0,fNGloSize*sizeof(double));
  //
  fNLocFits         = 0;
  fNLocFitsRejected = 0;
  fNLocEquations      = 0;
  //
  //  Process data records and build the matrices
  Long_t ndr = fTreeData->GetEntries();
  AliInfo(Form("Building the Global matrix from %d data records",ndr));
  for (Long_t i=0;i<ndr;i++) {
    ReadRecordData(i);
    LocalFit();
  }
  //
  fNGloFix = 0;
  for (int i=fNGloPar;i--;) fDiagCGlo[i] = MatCGlo.QuerryDiag(i); // save the diagonal elements  
  //
  //  Account for fixed parameters
  for (int i=fNGloPar;i--;) { // // Reset row and column of fixed params and add 1/sig^2 to free ones
    if ( fSigmaPar[i] <= 0. || fProcPnt[i]<fMinPntValid) {
      fNGloFix++; 
      fProcPnt[i] *= -1;
      fVecBGlo[i] = 0.;
      if (IsGlobalMatSparse()) {
	AliVectorSparse& row = *((AliMatrixSparse*)fMatCGlo)->GetRow(i);
	row.Clear();
	row(i) = float(fNLocEquations);
	for (int j=i+1;j<fNGloPar;j++) ((AliMatrixSparse*)fMatCGlo)->Zero(i,j);
      }
      else 
	for (int j=fNGloPar;j--;) if (MatCGlo.Querry(i,j)) MatCGlo(i,j)=0;
	MatCGlo.DiagElem(i) = float(fNLocEquations);
    }
    else MatCGlo.DiagElem(i) += 1.0/(fSigmaPar[i]*fSigmaPar[i]);
  }
  //
  // add constraint equations
  int nVar = fNGloPar;                    // Current size of global matrix	
  for (int i=0; i<fNGloConstraints; i++) {
    ReadRecordConstraint(i);
    double val   = fRecord->GetValue(0);  
    int    *indV = fRecord->GetIndex()+2;
    double *der  = fRecord->GetValue()+2;
    int    csize = fRecord->GetSize()-2;
    //
    // check if after suppression of fixed variables this there are non-0 derivatives
    int NSuppressed = 0;
    for (int j=csize;j--;)  if (fProcPnt[indV[j]]<1) NSuppressed++;
    //
    if (NSuppressed==csize) {
      AliInfo(Form("Neglecting constraint %d of %d derivatives since no free parameters left",i,csize));
      //
      // was this constraint ever created ? 
      if ( fConstrUsed[i] ) {
	// to avoid empty row impose dummy constraint on "Lagrange multiplier"
	MatCGlo.DiagElem(nVar) = 1.;
	fVecBGlo[nVar++] = 0;
      }
      continue;
    }
    //
    for (int j=csize; j--;) {
      int jID = indV[j];
      val -= der[j]*(fInitPar[jID]+fDeltaPar[jID]);
      if (fProcPnt[jID]<1) continue;                      // this parameter was fixed, don't put it into constraint
      MatCGlo(nVar,jID) = float(fNLocEquations)*der[j];   // fMatCGlo is symmetric, only lower triangle is filled  
    }
    //
    if (MatCGlo.QuerryDiag(nVar)) MatCGlo.DiagElem(nVar) = 0.0;
    fVecBGlo[nVar++] = float(fNLocEquations)*val; //RS ? should we use here fNLocFits ? 
    fConstrUsed[i] = kTRUE;
  }
  //
  AliInfo(Form("Obtained %-7ld equations from %-7ld records (%-7ld rejected). Fixed %-4d globals",
	       fNLocEquations,fNLocFits,fNLocFitsRejected,fNGloFix));

  //
  fGloSolveStatus = SolveGlobalMatEq();                     // obtain solution for this step
  //
  sw.Stop();
  AliInfo(Form("Iteration#%2d %s. CPU time: %.1f",fIter,fGloSolveStatus==gkFailed ? "Failed":"Converged",sw.CpuTime()));
  if (fGloSolveStatus==gkFailed) return 0;
  //
  for (int i=fNGloPar;i--;) fDeltaPar[i] += fVecBGlo[i];    // Update global parameters values (for iterations)
  //
  return 1;
}

//_____________________________________________________________________________________________
Int_t AliMillePede2::SolveGlobalMatEq()
{
  //
  // solve global matrix equation MatCGlob*X=VecBGlo and store the result in the VecBGlo
  //
  if (!fgIsMatGloSparse) {
    //
    if (fNGloConstraints==0) { // pos-def systems are faster to solve by Cholesky
      if ( ((AliSymMatrix*)fMatCGlo)->SolveChol(fVecBGlo, kTRUE) ) return gkInvert;
      else AliInfo("Solution of Global Dense System by Cholesky failed, trying Gaussian Elimiation");
    }
    //
    if (((AliSymMatrix*)fMatCGlo)->SolveSpmInv(fVecBGlo, kTRUE)) return gkInvert;
    else AliInfo("Solution of Global Dense System by Gaussian Elimiation failed, trying iterative methods");
  }
  // try to solve by minres
  TVectorD SOL(fNGloSize);
  //
  AliMinResSolve *slv = new AliMinResSolve(fMatCGlo,fVecBGlo);
  if (!slv) return gkFailed;
  //
  Bool_t res = kFALSE;
  if      (fgIterSol == AliMinResSolve::kSolMinRes) 
    res =  slv->SolveMinRes(SOL,fgMinResCondType,fgMinResMaxIter,fgMinResTol);
  else if (fgIterSol == AliMinResSolve::kSolFGMRes) 
    res =  slv->SolveFGMRES(SOL,fgMinResCondType,fgMinResMaxIter,fgMinResTol,fgNKrylovV);
  else 
    AliInfo(Form("Undefined Iteritive Solver ID=%d, only %d are defined",fgIterSol,AliMinResSolve::kNSolvers));
  //
  if (!res) return gkFailed;
  for (int i=fNGloSize;i--;) fVecBGlo[i] = SOL[i];
  return gkMinRes;
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
    double res = fMatCGlo->QuerryDiag(iPar);
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
    if (fGloSolveStatus==gkInvert && TMath::Abs( (dg=fMatCGlo->QuerryDiag(i)) *fDiagCGlo[i]) > 0) {    
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
