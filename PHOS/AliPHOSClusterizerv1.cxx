/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/* $Log:
   1 October 2000. Yuri Kharlov:
     AreNeighbours()
     PPSD upper layer is considered if number of layers>1

   18 October 2000. Yuri Kharlov:
     AliPHOSClusterizerv1()
     CPV clusterizing parameters added

     MakeClusters()
     After first PPSD digit remove EMC digits only once
*/
//*-- Author: Yves Schutz (SUBATECH)  & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//////////////////////////////////////////////////////////////////////////////
//  Clusterization class. Performs clusterization (collects neighbouring active cells) and 
//  unfolding of the clusters with several local maxima.  
//  results are stored in TreeR#, branches PHOSEmcRP (EMC recPoints),
//  PHOSCpvRP (CPV RecPoints) and AliPHOSClusterizer (Clusterizer with all 
//  parameters including input digits branch file name thresholds etc.)
//  This TTask normally called from Reconstructioner, but as well can be used it in 
//  standalone mode:
// root [0] AliPHOSClusterizerv1 * cl = new AliPHOSClusterizerv1("galice.root")  
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//               //reads gAlice from header file "..."                      
// root [1] cl->ExecuteTask()  
//               //finds RecPoints in all events stored in galice.root
// root [2] cl->SetDigitsBranch("digits2") 
//               //sets another input file
// root [3] cl->SetRecPointsBranch("recp2")  
//               //sets another aouput file
// root [4] cl->SetEmcLocalMaxCut(0.03)  
//               //set clusterization parameters
// root [5] cl->ExecuteTask("deb all time")  
//               //once more finds RecPoints options are 
//               // deb - print number of found rec points
//               // deb all - print number of found RecPoints and some their characteristics 
//               // time - print benchmarking results

// --- ROOT system ---

#include "TROOT.h" 
#include "TFile.h" 
#include "TFolder.h" 
#include "TMath.h" 
#include "TMinuit.h"
#include "TTree.h" 
#include "TSystem.h" 
#include "TBenchmark.h"

// --- Standard library ---

#include <iostream.h>
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliPHOSClusterizerv1.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSDigit.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOS.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliRun.h"

ClassImp(AliPHOSClusterizerv1)

//____________________________________________________________________________
  AliPHOSClusterizerv1::AliPHOSClusterizerv1():AliPHOSClusterizer()
{
  // default ctor (to be used)
  SetName("AliPHOSClusterizer");
  SetTitle("Version 1") ;

  fNumberOfCpvClusters     = 0 ; 
  fNumberOfEmcClusters     = 0 ; 
    
  fCpvClusteringThreshold  = 0.0;
  fEmcClusteringThreshold  = 0.2;   
  fPpsdClusteringThreshold = 0.0000002 ;
  
  fEmcLocMaxCut            = 0.03 ;
  fCpvLocMaxCut            = 0.03 ;
  
  fW0                      = 4.5 ;
  fW0CPV                   = 4.0 ;

  fGeom  = 0 ;
  
  fDigits = 0 ;
  fDigitizer = 0 ;
  fEmcRecPoints = 0 ;
  fCpvRecPoints = 0 ;

  fIsInitialized = kFALSE ;
  
}
//____________________________________________________________________________
  AliPHOSClusterizerv1::AliPHOSClusterizerv1(const char* headerFile,const char* digitsFile):AliPHOSClusterizer()
{
  SetName("AliPHOSClusterizer");
  SetTitle("Version 1") ;
  
  fNumberOfCpvClusters     = 0 ; 
  fNumberOfEmcClusters     = 0 ; 
    
  fCpvClusteringThreshold  = 0.0;
  fEmcClusteringThreshold  = 0.2;   
  fPpsdClusteringThreshold = 0.0000002 ;
  
  fEmcLocMaxCut            = 0.03 ;
  fCpvLocMaxCut            = 0.03 ;
  
  fW0                      = 4.5 ;
  fW0CPV                   = 4.0 ;
  
  fToUnfold = kTRUE ;
  
  fHeaderFileName = headerFile  ;
  fDigitsBranchTitle = digitsFile ;
  
  TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data() ) ;
  
  if(file == 0){
    file = new TFile(fHeaderFileName.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }
  
  AliPHOS * phos = (AliPHOS *) gAlice->GetDetector("PHOS") ;    
  fGeom  = AliPHOSGeometry::GetInstance(phos->GetGeometry()->GetName(),phos->GetGeometry()->GetTitle() );
  
  fDigits = new TClonesArray("AliPHOSDigit",10) ;
  fDigitizer = new AliPHOSDigitizer() ;
  fEmcRecPoints = new TObjArray(200) ;
  fCpvRecPoints = new TObjArray(200) ;
  
  if(!gMinuit) gMinuit = new TMinuit(100) ;
  
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 


  fIsInitialized = kTRUE ;

}
//____________________________________________________________________________
void AliPHOSClusterizerv1::Exec(Option_t * option){
  // Steerign function

  if(!fIsInitialized) Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSClusterizer");  
  
  Int_t nEvents = (Int_t) gAlice->TreeE()->GetEntries() ;
  
  for(fEvent = 0;fEvent< nEvents; fEvent++){
    if(!ReadDigits())  //reads digits for event fEvent
      return;
    MakeClusters() ;
    
    if(fToUnfold) MakeUnfolding() ;
    WriteRecPoints() ;
    if(strstr(option,"deb"))
      PrintRecPoints(option) ;
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSClusterizer");
    cout << "AliPHOSClusterizer:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSClusterizer") << " seconds for Clusterizing " 
	 <<  gBenchmark->GetCpuTime("PHOSClusterizer")/nEvents << " seconds per event " << endl ;
    cout << endl ;
  }
  
  
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::FindFit(AliPHOSEmcRecPoint * emcRP, int * maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters)
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 

  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliPHOSClusterizerv1::UnfoldingChiSquare) ;  
                                         // To set the address of the minimization function 

  TList * toMinuit = new TList();
  toMinuit->AddAt(emcRP,0) ;
  toMinuit->AddAt(fDigits,1) ;
  
  gMinuit->SetObjectFit(toMinuit) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliPHOSDigit * digit ;

  Int_t ierflg  = 0; 
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;


  for(iDigit = 0; iDigit < nDigits; iDigit++){
    digit = (AliPHOSDigit *) maxAt[iDigit]; 

    Int_t relid[4] ;
    Float_t x ;
    Float_t z ;
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    fGeom->RelPosInModule(relid, x, z) ;

    Float_t energy = maxAtEnergy[iDigit] ;

    gMinuit->mnparm(index, "x",  x, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){ 
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : x = " << x << endl ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "z",  z, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : z = " << z << endl ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;
    index++ ;   
    if(ierflg != 0){
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : energy = " << energy << endl ;      
      return kFALSE;
    }
  }

  Double_t p0 = 0.1 ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; The number of function call slightly
                      //  depends on it. 
  Double_t p1 = 1.0 ;
  Double_t p2 = 0.0 ;

  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls  
  gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient  
  gMinuit->SetMaxIterations(5);
  gMinuit->mnexcm("SET NOW", &p2 , 0, ierflg) ;  // No Warnings

  gMinuit->mnexcm("MIGRAD", &p0, 0, ierflg) ;    // minimize 

  if(ierflg == 4){  // Minimum not found   
    cout << "PHOS Unfolding>  Fit not converged, cluster abandoned "<< endl ;      
    return kFALSE ;
  }            
  for(index = 0; index < nPar; index++){
    Double_t err ;
    Double_t val ;
    gMinuit->GetParameter(index, val, err) ;    // Returns value and error of parameter index
    fitparameters[index] = val ;
   }

  delete toMinuit ;
  return kTRUE;

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::Init(){

  if(!fIsInitialized){
    if(fHeaderFileName.IsNull())
      fHeaderFileName = "galice.root" ;
    

    TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data() ) ;

    if(file == 0){
      file = new TFile(fHeaderFileName.Data(),"update") ;
      gAlice = (AliRun *) file->Get("gAlice") ;
    }

    AliPHOS * phos = (AliPHOS *) gAlice->GetDetector("PHOS") ;    
    fGeom  = AliPHOSGeometry::GetInstance(phos->GetGeometry()->GetName(),phos->GetGeometry()->GetTitle() );

    fDigits = new TClonesArray("AliPHOSDigit",10) ;
    fDigitizer = new AliPHOSDigitizer() ;
    fEmcRecPoints = new TObjArray(200) ;
    fCpvRecPoints = new TObjArray(200) ;
     
    if(!gMinuit) gMinuit = new TMinuit(100) ;

    // add Task to //root/Tasks folder
    TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
    roottasks->Add(this) ; 

    fIsInitialized = kTRUE ;
  }
}
//____________________________________________________________________________
Int_t AliPHOSClusterizerv1::AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)const
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 are not neighbour but do not continue searching
  // neighbours are defined as digits having at least common vertex 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  

  Int_t rv = 0 ; 

  Int_t relid1[4] ; 
  fGeom->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  fGeom->AbsToRelNumbering(d2->GetId(), relid2) ; 
 
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same PHOS module and the same PPSD Module 
    Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( coldiff <= 1 )  && ( rowdiff <= 1 )){
      rv = 1 ; 
    }
    else {
      if((relid2[2] > relid1[2]) && (relid2[3] > relid1[3]+1)) 
	rv = 2; //  Difference in row numbers is too large to look further 
    }

  } 
  else {
    
    if( (relid1[0] < relid2[0]) || (relid1[1] < relid2[1]) )  
      rv=2 ;

  }

  //Do NOT clusterize upper PPSD  
  if( IsInPpsd(d1) && IsInPpsd(d2) &&
     relid1[1] > 0                 &&
     relid1[1] < fGeom->GetNumberOfPadsPhi()*fGeom->GetNumberOfPadsPhi() ) rv = 2 ;

  return rv ; 
}


//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInEmc(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-EMC module
 
  Bool_t rv = kFALSE ; 

  Int_t relid[4] ; 
  fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] == 0  ) rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInPpsd(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-PPSD module
 
  Bool_t rv = kFALSE ; 

  Int_t relid[4] ; 
  fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] > 0 && relid[0] > fGeom->GetNCPVModules() ) rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInCpv(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-CPV module
 
  Bool_t rv = kFALSE ; 

  Int_t relid[4] ; 
  fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] > 0 && relid[0] <= fGeom->GetNCPVModules() ) rv = kTRUE; 

  return rv ; 
}
//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::ReadDigits(){

  fNumberOfEmcClusters  = 0 ;
  fNumberOfCpvClusters  = 0 ;

  // Get Digits Tree header from file
  char treeName[20]; 
  sprintf(treeName,"TreeD%d",fEvent);
  gAlice->GetEvent(fEvent) ;
  gAlice->SetEvent(fEvent) ;

  TTree * treeD = gAlice->TreeD()  ; // (TTree*)file->Get(treeName);

  if(treeD==0){
    cout << "Error in AliPHOSClusterizerv1 : no "<<treeName << endl  ;
    cout << "    Do nothing " << endl ;
    return kFALSE ;
  }

  TBranch * digitsBranch = 0;
  TBranch * digitizerBranch = 0;

  TObjArray * branches = treeD->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t phosNotFound = kTRUE ;
  Bool_t digitizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){

    if(phosNotFound){
      digitsBranch=(TBranch *) branches->At(ibranch) ;
      if( fDigitsBranchTitle.CompareTo(digitsBranch->GetTitle())==0 )
	if( strcmp(digitsBranch->GetName(),"PHOS") == 0)
	  phosNotFound = kFALSE ;
    }
    
    if(digitizerNotFound){
      digitizerBranch = (TBranch *) branches->At(ibranch) ;
      if( fDigitsBranchTitle.CompareTo(digitizerBranch->GetTitle()) == 0)
	if( strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) 
	  digitizerNotFound = kFALSE ;
    }
    
  }
  
  if(digitizerNotFound || phosNotFound){
    cout << "ERROR in AliPHOSClusterizerv1: " << endl ;
    cout << "      Can't find Branch with digits or Digitizer "<< endl ; ;
    cout << "      Do nothing" <<endl  ;
    return kFALSE ;
  }
  
  digitsBranch->SetAddress(&fDigits) ;
  digitizerBranch->SetAddress(&fDigitizer) ;
  
  treeD->GetEvent(0) ;
  
  fPedestal = fDigitizer->GetPedestal() ;
  fSlope    = fDigitizer->GetSlope() ;
  return kTRUE ;
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::WriteRecPoints(){
  
  Int_t index ;
  //Evaluate poisition, dispersion and other RecPoint properties...
  for(index = 0; index < fEmcRecPoints->GetEntries(); index++)
    ((AliPHOSEmcRecPoint *)fEmcRecPoints->At(index))->EvalAll(fW0,fDigits) ;

  fEmcRecPoints->Sort() ;

  for(index = 0; index < fEmcRecPoints->GetEntries(); index++)
    ((AliPHOSEmcRecPoint *)fEmcRecPoints->At(index))->SetIndexInList(index) ;

  //Now the same for CPV
  for(index = 0; index < fCpvRecPoints->GetEntries(); index++)
    ((AliPHOSRecPoint *)fCpvRecPoints->At(index))->EvalAll(fW0CPV,fDigits)  ;

  fCpvRecPoints->Sort() ;

  for(index = 0; index < fCpvRecPoints->GetEntries(); index++)
    ((AliPHOSRecPoint *)fCpvRecPoints->At(index))->SetIndexInList(index) ;

  if(gAlice->TreeR()==0)
    gAlice->MakeTree("R") ;
  

  //Check, if branches already exist
  TBranch * emcBranch = 0;
  TBranch * cpvBranch = 0;
  TBranch * clusterizerBranch = 0;
  
  TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t emcNotFound = kTRUE ;
  Bool_t cpvNotFound = kTRUE ;  
  Bool_t clusterizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    
    if(emcNotFound){
      emcBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsBranchTitle.CompareTo(emcBranch->GetTitle())==0 ){
	if( strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0) {
	  emcNotFound = kFALSE ;
	}
      }
    }
    if(cpvNotFound){
      cpvBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsBranchTitle.CompareTo(cpvBranch->GetTitle())==0 ){
	if( strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) {
	  cpvNotFound = kFALSE ;
	}
      }
    }

    if(clusterizerNotFound){
      clusterizerBranch = (TBranch *) branches->At(ibranch) ;
      if( fRecPointsBranchTitle.CompareTo(clusterizerBranch->GetTitle()) == 0){
	if( strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) {
	  clusterizerNotFound = kFALSE ;
	}
      }
    }
    
  }
  
  if(!(clusterizerNotFound && emcNotFound && cpvNotFound)){
    cout << "AliPHOSClusterizer error" << endl;
    cout << "       Branches PHOSEmcRP, PHOSCpvRP and AliPHOSClusterizer " << endl ;
    cout << "       with title '" << fRecPointsBranchTitle.Data() <<"' already exist" << endl ;
    cout << "       can not overwrite " << endl ;
    return ;
  }

  //Make branches in TreeR for RecPoints and Clusterizer
  char * filename = 0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")!=0){   //generating file name
    filename = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(filename,"%s/PHOS.Reco.root",gAlice->GetBaseFile()) ;
  }
  
  //Make new branches
  TDirectory *cwd = gDirectory;
  
  //First EMC
  Int_t bufferSize = 32000 ;    
  Int_t splitlevel = 0 ;
  emcBranch = gAlice->TreeR()->Branch("PHOSEmcRP","TObjArray",&fEmcRecPoints,bufferSize,splitlevel);
  emcBranch->SetTitle(fRecPointsBranchTitle.Data());
  if (filename) {
    emcBranch->SetFile(filename);
    TIter next( emcBranch->GetListOfBranches());
    while ((emcBranch=(TBranch*)next())) {
      emcBranch->SetFile(filename);
    }   
    cwd->cd();
  }
    
  //Now CPV branch
  cpvBranch = gAlice->TreeR()->Branch("PHOSCpvRP","TObjArray",&fCpvRecPoints,bufferSize,splitlevel);
  cpvBranch->SetTitle(fRecPointsBranchTitle.Data());
  if (filename) {
    cpvBranch->SetFile(filename);
    TIter next( cpvBranch->GetListOfBranches());
    while ((cpvBranch=(TBranch*)next())) {
      cpvBranch->SetFile(filename);
    }   
    cwd->cd();
  } 
    
  //And Finally  clusterizer branch
  AliPHOSClusterizerv1 * cl = this ;
  clusterizerBranch = gAlice->TreeR()->Branch("AliPHOSClusterizer","AliPHOSClusterizerv1",
					      &cl,bufferSize,splitlevel);
  clusterizerBranch->SetTitle(fRecPointsBranchTitle.Data());
  if (filename) {
    clusterizerBranch->SetFile(filename);
    TIter next( clusterizerBranch->GetListOfBranches());
    while ((clusterizerBranch=(TBranch*)next())) {
      clusterizerBranch->SetFile(filename);
    }   
    cwd->cd();
  }
  
  gAlice->TreeR()->Fill() ;
  
  gAlice->TreeR()->Write(0,kOverwrite) ;  
  
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
  fEmcRecPoints->Clear() ;
  fCpvRecPoints->Clear() ;
  

  // Clusterization starts  
  TClonesArray * digits =  (TClonesArray*)fDigits->Clone() ;

  TIter nextdigit(digits) ; 
  AliPHOSDigit * digit ; 
  Bool_t notremoved = kTRUE ;

  while ( (digit = (AliPHOSDigit *)nextdigit()) ) { // scan over the list of digits
    AliPHOSRecPoint * clu = 0 ; 

    TArrayI clusterdigitslist(1000) ;   
    Int_t index ;

    if (( IsInEmc (digit) && Calibrate(digit->GetAmp()) > fEmcClusteringThreshold  ) || 
        ( IsInPpsd(digit) && Calibrate(digit->GetAmp()) > fPpsdClusteringThreshold ) ||
        ( IsInCpv (digit) && Calibrate(digit->GetAmp()) > fCpvClusteringThreshold  ) ) {
      
      Int_t iDigitInCluster = 0 ; 

      if  ( IsInEmc(digit) ) {   
	// start a new EMC RecPoint
	if(fNumberOfEmcClusters >= fEmcRecPoints->GetSize()) fEmcRecPoints->Expand(2*fNumberOfEmcClusters+1) ;
	fEmcRecPoints->AddAt(new  AliPHOSEmcRecPoint(), fNumberOfEmcClusters) ;
	clu = (AliPHOSEmcRecPoint *) fEmcRecPoints->At(fNumberOfEmcClusters) ; 
	fNumberOfEmcClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp())) ; 
	clusterdigitslist[iDigitInCluster] = digit->GetIndexInList() ;	
	iDigitInCluster++ ; 
	digits->Remove(digit) ; 

      } else { 
	
	// start a new PPSD/CPV cluster
	if(fNumberOfCpvClusters >= fCpvRecPoints->GetSize()) fCpvRecPoints->Expand(2*fNumberOfCpvClusters+1);
	if(IsInPpsd(digit)) 
	  fCpvRecPoints->AddAt(new AliPHOSPpsdRecPoint(),fNumberOfCpvClusters) ;
	else
	  fCpvRecPoints->AddAt(new AliPHOSCpvRecPoint(), fNumberOfCpvClusters) ;
	clu =  (AliPHOSPpsdRecPoint *) fCpvRecPoints->At(fNumberOfCpvClusters)  ;  
	fNumberOfCpvClusters++ ; 

	clu->AddDigit(*digit, Calibrate(digit->GetAmp()) ) ;	
	clusterdigitslist[iDigitInCluster] = digit->GetIndexInList()  ;	
	iDigitInCluster++ ; 
	digits->Remove(digit) ; 
        nextdigit.Reset() ;
	
	// Here we remove resting EMC digits, which cannot make cluster
	
        if( notremoved ) { 
	  while( ( digit = (AliPHOSDigit *)nextdigit() ) ) {
            if( IsInEmc(digit) ) 
	      digits->Remove(digit) ;
            else 
	      break ;
	  }
	  notremoved = kFALSE ;
	}
	
      } // else        
      
      nextdigit.Reset() ;
      
      AliPHOSDigit * digitN ; 
      index = 0 ;
      while (index < iDigitInCluster){ // scan over digits already in cluster 
	digit =  (AliPHOSDigit*)fDigits->At(clusterdigitslist[index])  ;      
	index++ ; 
        while ( (digitN = (AliPHOSDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
          switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;
	  case 1 :   // are neighbours 
	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp() ) ) ;
	    clusterdigitslist[iDigitInCluster] = digitN->GetIndexInList() ; 
	    iDigitInCluster++ ; 
	    digits->Remove(digitN) ;
	    break ;
          case 2 :   // too far from each other
	    goto endofloop;   
	  } // switch
	  
	} // while digitN
	
      endofloop: ;
	nextdigit.Reset() ; 
	
      } // loop over cluster     

    } // energy theshold  

    
  } // while digit

  delete digits ;

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeUnfolding(){
  //Unfolds clusters using the shape of ElectroMagnetic shower
  // Performs unfolding of all EMC/CPV but NOT ppsd clusters

  //Unfold first EMC clusters 
  if(fNumberOfEmcClusters > 0){

    Int_t nModulesToUnfold = fGeom->GetNModules() ; 

    Int_t numberofNotUnfolded = fNumberOfEmcClusters ; 
    Int_t index ;   
    for(index = 0 ; index < numberofNotUnfolded ; index++){
      
      AliPHOSEmcRecPoint * emcRecPoint = (AliPHOSEmcRecPoint *) fEmcRecPoints->At(index) ;
      if(emcRecPoint->GetPHOSMod()> nModulesToUnfold)
	break ;
      
      Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
      Int_t * maxAt = new Int_t[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fEmcLocMaxCut,fDigits) ;
      
      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       
	UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;
	fEmcRecPoints->Remove(emcRecPoint); 
	fEmcRecPoints->Compress() ;
	index-- ;
	fNumberOfEmcClusters -- ;
	numberofNotUnfolded-- ;
      }
      
      delete[] maxAt ; 
      delete[] maxAtEnergy ; 
    }
  } 
  //Unfolding of EMC clusters finished


  //Unfold now CPV clusters
  if(fNumberOfCpvClusters > 0){
    
    Int_t nModulesToUnfold = fGeom->GetNCPVModules() ;

    Int_t numberofCpvNotUnfolded = fNumberOfCpvClusters ;     
    Int_t index ;   
    for(index = 0 ; index < numberofCpvNotUnfolded ; index++){
      
      AliPHOSRecPoint * recPoint = (AliPHOSRecPoint *) fCpvRecPoints->At(index) ;

      if(recPoint->GetPHOSMod()> nModulesToUnfold)
	break ;
      
      AliPHOSEmcRecPoint * emcRecPoint = (AliPHOSEmcRecPoint*) recPoint ; 
      
      Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
      Int_t * maxAt = new Int_t[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fCpvLocMaxCut,fDigits) ;
      
      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       
	UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;
	fCpvRecPoints->Remove(emcRecPoint); 
	fCpvRecPoints->Compress() ;
	index-- ;
	numberofCpvNotUnfolded-- ;
	fNumberOfCpvClusters-- ;
      }
      
      delete[] maxAt ; 
      delete[] maxAtEnergy ; 
    } 
  }
  //Unfolding of Cpv clusters finished
  
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::SetDigitsBranch(const char * title){
  
    fDigitsBranchTitle = title  ; 

}
//____________________________________________________________________________
void AliPHOSClusterizerv1::SetRecPointsBranch(const char * title){
  
    fRecPointsBranchTitle = title;

}

//____________________________________________________________________________
Double_t  AliPHOSClusterizerv1::ShowerShape(Double_t r)
{ 
  // Shape of the shower (see PHOS TDR)
  // If you change this function, change also the gradien evaluation  in ChiSquare()

  Double_t r4    = r*r*r*r ;
  Double_t r295  = TMath::Power(r, 2.95) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliPHOSClusterizerv1::UnfoldCluster(AliPHOSEmcRecPoint * iniEmc, 
						 Int_t nMax, 
						 int * maxAt, 
						 Float_t * maxAtEnergy)
{ 
  // Performs the unfolding of a cluster with nMax overlapping showers 

  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;

  Bool_t rv = FindFit(iniEmc, maxAt, maxAtEnergy, nPar, fitparameters) ;
  if( !rv ) {
    // Fit failed, return and remove cluster
    delete[] fitparameters ; 
    return ;
  }

  // create ufolded rec points and fill them with new energy lists
  // First calculate energy deposited in each sell in accordance with fit (without fluctuations): efit[]
  // and later correct this number in acordance with actual energy deposition

  Int_t nDigits = iniEmc->GetMultiplicity() ;  
  Float_t * efit = new Float_t[nDigits] ;
  Float_t xDigit,zDigit,distance ;
  Float_t xpar,zpar,epar  ;
  Int_t relid[4] ;
  AliPHOSDigit * digit ;
  Int_t * emcDigits = iniEmc->GetDigitsList() ;

  Int_t iparam ;
  Int_t iDigit ;
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
    digit = (AliPHOSDigit*) fDigits->At(emcDigits[iDigit] ) ;   
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    fGeom->RelPosInModule(relid, xDigit, zDigit) ;
    efit[iDigit] = 0;

    iparam = 0 ;    
    while(iparam < nPar ){
      xpar = fitparameters[iparam] ;
      zpar = fitparameters[iparam+1] ;
      epar = fitparameters[iparam+2] ;
      iparam += 3 ;
      distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      distance =  TMath::Sqrt(distance) ;
      efit[iDigit] += epar * ShowerShape(distance) ;
    }
  }
  

  // Now create new RecPoints and fill energy lists with efit corrected to fluctuations
  // so that energy deposited in each cell is distributed betwin new clusters proportionally
  // to its contribution to efit

  Float_t * emcEnergies = iniEmc->GetEnergiesList() ;
  Float_t ratio ;

  iparam = 0 ;
  while(iparam < nPar ){
    xpar = fitparameters[iparam] ;
    zpar = fitparameters[iparam+1] ;
    epar = fitparameters[iparam+2] ;
    iparam += 3 ;    
    
    AliPHOSEmcRecPoint * emcRP ;  

    if(iniEmc->IsEmc()){ //create new entries in fEmcRecPoints...
      
      if(fNumberOfEmcClusters >= fEmcRecPoints->GetSize())
	fEmcRecPoints->Expand(2*fNumberOfEmcClusters) ;
      
      (*fEmcRecPoints)[fNumberOfEmcClusters] = new AliPHOSEmcRecPoint() ;
      emcRP = (AliPHOSEmcRecPoint *) fEmcRecPoints->At(fNumberOfEmcClusters);
      fNumberOfEmcClusters++ ;
    }
    else{//create new entries in fCpvRecPoints
      if(fNumberOfCpvClusters >= fCpvRecPoints->GetSize())
	fCpvRecPoints->Expand(2*fNumberOfCpvClusters) ;
      
      (*fCpvRecPoints)[fNumberOfCpvClusters] = new AliPHOSCpvRecPoint() ;
      emcRP = (AliPHOSEmcRecPoint *) fCpvRecPoints->At(fNumberOfCpvClusters);
      fNumberOfCpvClusters++ ;
    }
    
    Float_t eDigit ;
    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
      digit = (AliPHOSDigit*) fDigits->At( emcDigits[iDigit] ) ; 
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      fGeom->RelPosInModule(relid, xDigit, zDigit) ;
      distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      distance =  TMath::Sqrt(distance) ;
      ratio = epar * ShowerShape(distance) / efit[iDigit] ; 
      eDigit = emcEnergies[iDigit] * ratio ;
      emcRP->AddDigit( *digit, eDigit ) ;
    }	
  }
 
  delete[] fitparameters ; 
  delete[] efit ; 
  
}

//_____________________________________________________________________________
void AliPHOSClusterizerv1::UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates th Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do


  TList * toMinuit = (TList*) gMinuit->GetObjectFit() ;

  AliPHOSEmcRecPoint * emcRP = (AliPHOSEmcRecPoint*) toMinuit->At(0)  ;
  TClonesArray * digits = (TClonesArray*)toMinuit->At(1)  ;


  
  //  AliPHOSEmcRecPoint * emcRP = (AliPHOSEmcRecPoint *) gMinuit->GetObjectFit() ; // EmcRecPoint to fit

  Int_t * emcDigits     = emcRP->GetDigitsList() ;

  Int_t nOfDigits = emcRP->GetDigitsMultiplicity() ; 

  Float_t * emcEnergies = emcRP->GetEnergiesList() ;

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  fret = 0. ;     
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < nPar ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t efit ;    

  AliPHOSDigit * digit ;
  Int_t iDigit ;

  for( iDigit = 0 ; iDigit < nOfDigits ; iDigit++) {

    digit = (AliPHOSDigit*) digits->At( emcDigits[iDigit] ) ; 

    Int_t relid[4] ;
    Float_t xDigit ;
    Float_t zDigit ;

    geom->AbsToRelNumbering(digit->GetId(), relid) ;

    geom->RelPosInModule(relid, xDigit, zDigit) ;

     if(iflag == 2){  // calculate gradient
       Int_t iParam = 0 ;
       efit = 0 ;
       while(iParam < nPar ){
	 Double_t distance = (xDigit - x[iParam]) * (xDigit - x[iParam]) ;
	 iParam++ ; 
	 distance += (zDigit - x[iParam]) * (zDigit - x[iParam]) ; 
	 distance = TMath::Sqrt( distance ) ; 
	 iParam++ ; 	 
	 efit += x[iParam] * ShowerShape(distance) ;
	 iParam++ ;
       }
       Double_t sum = 2. * (efit - emcEnergies[iDigit]) / emcEnergies[iDigit] ; // Here we assume, that sigma = sqrt(E) 
       iParam = 0 ;
       while(iParam < nPar ){
	 Double_t xpar = x[iParam] ;
	 Double_t zpar = x[iParam+1] ;
	 Double_t epar = x[iParam+2] ;
	 Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );
	 Double_t shape = sum * ShowerShape(dr) ;
	 Double_t r4 = dr*dr*dr*dr ;
	 Double_t r295 = TMath::Power(dr,2.95) ;
	 Double_t deriv =-4. * dr*dr * ( 2.32 / ( (2.32 + 0.26 * r4) * (2.32 + 0.26 * r4) ) +
					 0.0316 * (1. + 0.0171 * r295) / ( ( 1. + 0.0652 * r295) * (1. + 0.0652 * r295) ) ) ;
	 
	 Grad[iParam] += epar * shape * deriv * (xpar - xDigit) ;  // Derivative over x    
	 iParam++ ; 
	 Grad[iParam] += epar * shape * deriv * (zpar - zDigit) ;  // Derivative over z         
	 iParam++ ; 
	 Grad[iParam] += shape ;                                  // Derivative over energy     	
	 iParam++ ; 
       }
     }
     efit = 0;
     iparam = 0 ;

     while(iparam < nPar ){
       Double_t xpar = x[iparam] ;
       Double_t zpar = x[iparam+1] ;
       Double_t epar = x[iparam+2] ;
       iparam += 3 ;
       Double_t distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
       distance =  TMath::Sqrt(distance) ;
       efit += epar * ShowerShape(distance) ;
     }

     fret += (efit-emcEnergies[iDigit])*(efit-emcEnergies[iDigit])/emcEnergies[iDigit] ; 
     // Here we assume, that sigma = sqrt(E)
  }

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::Print(Option_t * option)const
{
  if(fIsInitialized){
    
    // Print parameters
    
    cout << "---------------"<< GetName() << " " << GetTitle()<< "-----------" << endl 
	 << "Clusterizing digits from the file: " << fHeaderFileName.Data() << endl 
	 << "                           Branch: " << fDigitsBranchTitle.Data() << endl 
	 << endl 
	 << "                       EMC Clustering threshold = " << fEmcClusteringThreshold << endl
	 << "                       EMC Local Maximum cut    = " << fEmcLocMaxCut << endl
	 << "                       EMC Logarothmic weight   = " << fW0 << endl
	 << endl
	 << "                       CPV Clustering threshold = " << fCpvClusteringThreshold << endl
	 << "                       CPV Local Maximum cut    = " << fCpvLocMaxCut << endl
       << "                       CPV Logarothmic weight   = " << fW0CPV << endl
	 << endl
	 << "                      PPSD  Clustering threshold = " << fPpsdClusteringThreshold << endl;
    if(fToUnfold)
      cout << " Unfolding on " << endl ;
    else
      cout << " Unfolding off " << endl ;
    
    cout << "------------------------------------------------------------------" <<endl ;
  }
  else
    cout << " AliPHOSClusterizerv1 not initialized " << endl ;
}
//____________________________________________________________________________
void AliPHOSClusterizerv1::PrintRecPoints(Option_t * option){
  //Prints list of RecPoints produced at the current pass of AliPHOSClusterizer

  cout << "AliPHOSClusterizerv1: " << endl ;
  cout << "       Found "<< fEmcRecPoints->GetEntriesFast() << " EMC Rec Points and " 
	   << fCpvRecPoints->GetEntriesFast() << " CPV RecPoints" << endl ;

  if(strstr(option,"all")) {
    cout << "EMC clusters " << endl ;
    cout << "  Index    " 
	 << "  Ene(MeV) " 
	 << "  Multi    "
 	 << "  Module   "  
	 << "    X      "
	 << "    Y      "
	 << "    Z      "
	 << " Lambda 1  "
	 << " Lambda 2  "
	 << " MaxEnergy "
	 << " # of prim "
	 << " Primaries list "      <<  endl;      
    
    Int_t index ;
    for (index = 0 ; index < fEmcRecPoints->GetEntries() ; index++) {
      AliPHOSEmcRecPoint * rp = (AliPHOSEmcRecPoint * )fEmcRecPoints->At(index) ; 

      cout << setw(6) << rp->GetIndexInList() << "     ";
      cout << setw(6) << rp->GetEnergy()      << "     ";
      cout << setw(6) << rp->GetMultiplicity()<< "     ";
      cout << setw(6) << rp->GetPHOSMod()     << "     ";

      TVector3  locpos;  
      rp->GetLocalPosition(locpos);
      cout << setw(8) <<  locpos.X()          << "     ";
      cout << setw(8) <<  locpos.Y()          << "     ";
      cout << setw(8) <<  locpos.Z()          << "     ";

      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      cout << setw(10)<<  lambda[0]           << "     ";
      cout << setw(10)<<  lambda[1]           << "     ";
      
      
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      cout << setw(8) <<    primaries         << "     ";

      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++)
	cout << setw(4)  <<  primaries[iprimary] << " ";
      cout << endl;  	 
    }
    cout << endl ;

    //Now plot PCV/PPSD recPoints
    cout << "EMC clusters " << endl ;
    cout << "  Index    " 
	 << "  Multi    "
 	 << "  Module   "  
 	 << "  Layer    "  
	 << "    X      "
	 << "    Y      "
	 << "    Z      "
	 << " # of prim "
	 << " Primaries list "      <<  endl;      
    
    for (index = 0 ; index < fEmcRecPoints->GetEntries() ; index++) {
      AliPHOSRecPoint * rp = (AliPHOSRecPoint * )fCpvRecPoints->At(index) ; 
      cout << setw(6) << rp->GetIndexInList() << "     ";
      cout << setw(6) << rp->GetPHOSMod()     << "     ";

      if( (strcmp(rp->ClassName() , "AliPHOSPpsdRecPoint" )) == 0){
	AliPHOSPpsdRecPoint * ppsd = (AliPHOSPpsdRecPoint*) rp ;
	if(ppsd->GetUp())
	  cout <<"        CPV     ";
	else
	  cout <<"       PPSD     ";
      }
      else
	cout <<"        CPV     ";
      
      TVector3  locpos;  
      rp->GetLocalPosition(locpos);
      cout << setw(8) <<  locpos.X()          << "     ";
      cout << setw(8) <<  locpos.Y()          << "     ";
      cout << setw(8) <<  locpos.Z()          << "     ";
      
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      cout << setw(8) <<    primaries         << "     ";

      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++)
	cout << setw(4)  <<  primaries[iprimary] << " ";
      cout << endl;  	 
    }


    cout << "-------------------------------------------------"<<endl ;
  }
}

