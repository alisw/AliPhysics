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
     AliEMCALClusterizerv1()
     CPV clusterizing parameters added
     MakeClusters()
     After first PPSD digit remove EMC digits only once
*/
//*-- Author: Yves Schutz (SUBATECH)  & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//////////////////////////////////////////////////////////////////////////////
//  Clusterization class. Performs clusterization (collects neighbouring active cells) and 
//  unfolds the clusters having several local maxima.  
//  Results are stored in TreeR#, branches EMCALTowerRP (EMC recPoints),
//  EMCALPreShoRP (CPV RecPoints) and AliEMCALClusterizer (Clusterizer with all 
//  parameters including input digits branch title, thresholds etc.)
//  This TTask is normally called from Reconstructioner, but can as well be used in 
//  standalone mode.
// Use Case:
//  root [0] AliEMCALClusterizerv1 * cl = new AliEMCALClusterizerv1("galice.root")  
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//               //reads gAlice from header file "..."                      
//  root [1] cl->ExecuteTask()  
//               //finds RecPoints in all events stored in galice.root
//  root [2] cl->SetDigitsBranch("digits2") 
//               //sets another title for Digitis (input) branch
//  root [3] cl->SetRecPointsBranch("recp2")  
//               //sets another title four output branches
//  root [4] cl->SetTowerLocalMaxCut(0.03)  
//               //set clusterization parameters
//  root [5] cl->ExecuteTask("deb all time")  
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


// --- AliRoot header files ---

#include "AliEMCALClusterizerv1.h"
#include "AliEMCALDigit.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALTowerRecPoint.h"
#include "AliEMCAL.h"
#include "AliEMCALGetter.h"
#include "AliEMCALGeometry.h"
#include "AliRun.h"

ClassImp(AliEMCALClusterizerv1)
  
//____________________________________________________________________________
  AliEMCALClusterizerv1::AliEMCALClusterizerv1() : AliEMCALClusterizer()
{
  // default ctor (to be used mainly by Streamer)
  
  InitParameters() ; 
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(const char* headerFile, const char* name, const Bool_t toSplit)
:AliEMCALClusterizer(headerFile, name, toSplit)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 

}

//____________________________________________________________________________
  AliEMCALClusterizerv1::~AliEMCALClusterizerv1()
{
  // dtor
  fSplitFile = 0 ; 
  
}

//____________________________________________________________________________
const TString AliEMCALClusterizerv1::BranchName() const 
{ 
  TString branchName(GetName() ) ;
  branchName.Remove(branchName.Index(Version())-1) ;
  return branchName ;
}

//____________________________________________________________________________
Float_t  AliEMCALClusterizerv1::Calibrate(Int_t amp, Int_t where) const
{
  //To be replased later by the method, reading individual parameters from the database
  // where = 0 == PRE ; where = 1 == ECAL ; where = 2 == HCAL
  if ( where == 0 ) // calibrate as PRE section
    return -fADCpedestalPRE + amp * fADCchannelPRE ; 
  else if (where == 1) //calibrate as EC section 
    return -fADCpedestalEC + amp * fADCchannelEC ;
  else if (where == 2) //calibrate as HC section
    return -fADCpedestalHC + amp * fADCchannelHC ;
  else 
    Fatal("Calibrate", "Something went wrong!") ;
  return -9999999. ; 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Exec(Option_t * option)
{
  // Steering method

  if( strcmp(GetName(), "")== 0 ) 
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALClusterizer"); 
  
  if(strstr(option,"print"))
    Print("") ; 

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ;
  if(gime->BranchExists("RecPoints"))
    return ;
  Int_t nevents = gime->MaxEvent() ;
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){

    gime->Event(ievent,"D") ;

    if(ievent == 0)
      GetCalibrationParameters() ;

    fNumberOfPREClusters = fNumberOfECClusters = fNumberOfHCClusters = 0 ;
           
    MakeClusters() ;

    if(fToUnfold)
      MakeUnfolding() ;

    WriteRecPoints(ievent) ;

    if(strstr(option,"deb"))  
      PrintRecPoints(option) ;

    //increment the total number of digits per run 
    fRecPointsInRun += gime->PRERecPoints()->GetEntriesFast() ;  
    fRecPointsInRun += gime->ECALRecPoints()->GetEntriesFast() ;  
    fRecPointsInRun += gime->HCALRecPoints()->GetEntriesFast() ;  
 }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALClusterizer");
    Info("Exec", "took %f seconds for Clusterizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALClusterizer"), gBenchmark->GetCpuTime("EMCALClusterizer")/nevents ) ;
  }
  
}

//____________________________________________________________________________
Bool_t AliEMCALClusterizerv1::FindFit(AliEMCALTowerRecPoint * emcRP, AliEMCALDigit ** maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters) const
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 
  // The initial values for fitting procedure are set equal to the positions of local maxima.
  // Cluster will be fitted as a superposition of nPar/3 electromagnetic showers

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TClonesArray * digits = gime->Digits() ; 
  

  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliEMCALClusterizerv1::UnfoldingChiSquare) ;  
                                         // To set the address of the minimization function 
  TList * toMinuit = new TList();
  toMinuit->AddAt(emcRP,0) ;
  toMinuit->AddAt(digits,1) ;
  
  gMinuit->SetObjectFit(toMinuit) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliEMCALDigit * digit ;

  Int_t ierflg  = 0; 
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;

  AliEMCALGeometry * geom = gime->EMCALGeometry() ; 

  for(iDigit = 0; iDigit < nDigits; iDigit++){
    digit = maxAt[iDigit]; 

    Int_t relid[4] ;
    Float_t x = 0.;
    Float_t z = 0.;
    geom->AbsToRelNumbering(digit->GetId(), relid) ;
    geom->PosInAlice(relid, x, z) ;

    Float_t energy = maxAtEnergy[iDigit] ;

    gMinuit->mnparm(index, "x",  x, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){ 
      Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : x = %f",  x ) ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "z",  z, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){
       Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : z = %f", z) ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;
    index++ ;   
    if(ierflg != 0){
     Error("FindFit", "EMCAL Unfolding  Unable to set initial value for fit procedure : energy = %f", energy) ;      
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
    Error("FindFit", "EMCAL Unfolding  Fit not converged, cluster abandoned " ) ;      
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
void AliEMCALClusterizerv1::GetCalibrationParameters() 
{
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ;
  const AliEMCALDigitizer * dig = gime->Digitizer(BranchName()) ;

  fADCchannelPRE  = dig->GetPREchannel() ;
  fADCpedestalPRE = dig->GetPREpedestal() ; 

  fADCchannelEC   = dig->GetECchannel() ;
  fADCpedestalEC  = dig->GetECpedestal();

  fADCchannelHC   = dig->GetHCchannel() ;
  fADCpedestalHC  = dig->GetHCpedestal();
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of EMCAL tasks
  
  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;

  TString branchname = GetName() ;
  branchname.Remove(branchname.Index(Version())-1) ;

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance(GetTitle(), branchname.Data(), fToSplit ) ; 
  if ( gime == 0 ) {
    Error("Init", "Could not obtain the Getter object !" ) ; 
    return ;
  } 

  fSplitFile = 0 ;
  if(fToSplit){
    // construct the name of the file as /path/EMCAL.SDigits.root
    //First - extract full path if necessary
    TString fileName(GetTitle()) ;
    Ssiz_t islash = fileName.Last('/') ;
    if(islash<fileName.Length())
      fileName.Remove(islash+1,fileName.Length()) ;
    else
      fileName="" ;
    // Next - append the file name 
    fileName+="EMCAL.RecData." ;
    if((strcmp(branchname.Data(),"Default")!=0)&&(strcmp(branchname.Data(),"")!=0)){
      fileName+=branchname ;
      fileName+="." ;
    }
    fileName+="root" ;
    // Finally - check if the file already opened or open the file
    fSplitFile = static_cast<TFile*>(gROOT->GetFile(fileName.Data()));   
    if(!fSplitFile)
      fSplitFile =  TFile::Open(fileName.Data(),"update") ;
  }

  const AliEMCALGeometry * geom = gime->EMCALGeometry() ;
  fNTowers = geom->GetNZ() *  geom->GetNPhi() ;
  if(!gMinuit) 
    gMinuit = new TMinuit(100) ;

  gime->PostClusterizer(this) ;
  gime->PostRecPoints(branchname ) ;
 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::InitParameters()
{
  fNumberOfPREClusters = fNumberOfECClusters = fNumberOfHCClusters = 0 ;   
  fPREClusteringThreshold  = 0.0001; // must be adjusted according to the noise leve set by digitizer
  fECClusteringThreshold   = 0.0045;  // must be adjusted according to the noise leve set by digitizer
  fHCClusteringThreshold   = 0.001;  // must be adjusted according to the noise leve set by digitizer  
  fPRELocMaxCut = 0.03 ;
  fECLocMaxCut  = 0.03 ;
  fHCLocMaxCut  = 0.03 ;
  
  fPREW0    = 4.0 ;
  fECW0     = 4.5 ;
  fHCW0     = 4.5 ;

  fTimeGate = 1.e-8 ; 
  
  fToUnfold = kFALSE ;
  
  TString clusterizerName( GetName()) ; 
  if (clusterizerName.IsNull() ) 
    clusterizerName = "Default" ; 
  clusterizerName.Append(":") ; 
  clusterizerName.Append(Version()) ; 
  SetName(clusterizerName) ;
  fRecPointsInRun          = 0 ; 

}

//____________________________________________________________________________
Int_t AliEMCALClusterizerv1::AreNeighbours(AliEMCALDigit * d1, AliEMCALDigit * d2)const
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 are not neighbour but do not continue searching
  // neighbours are defined as digits having at least a common vertex 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  

   AliEMCALGeometry * geom = AliEMCALGetter::GetInstance()->EMCALGeometry() ;

  Int_t rv = 0 ; 

  Int_t relid1[4] ; 
  geom->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  geom->AbsToRelNumbering(d2->GetId(), relid2) ; 
  
  if ( (relid1[0] == relid2[0]) && // inside the same EMCAL Arm  
       (relid1[1]==relid2[1]) ) {  // and same tower section
    Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( coldiff <= 1 )  && ( rowdiff <= 1 )){
      if((relid1[1] != 0) || (TMath::Abs(d1->GetTime() - d2->GetTime() ) < fTimeGate))
      rv = 1 ; 
    }
    else {
      if((relid2[2] > relid1[2]) && (relid2[3] > relid1[3]+1)) 
	rv = 2; //  Difference in row numbers is too large to look further 
    }

  } 
  else {
    
    if( (relid1[0] < relid2[0]) || (relid1[1] != relid2[1]) )  
      rv=0 ;
  }

  if (gDebug == 2 ) 
    Info("AreNeighbours", "neighbours=%d, id1=%d, relid1=%d,%d,%d,%d \n id2=%d, relid2=%d,%d,%d,%d", 
	 rv, d1->GetId(), relid1[0], relid1[1], relid1[2], relid1[3], 
	 d2->GetId(), relid2[0], relid2[1], relid2[2], relid2[3]) ;   
  
  return rv ; 
}

 
//____________________________________________________________________________
void AliEMCALClusterizerv1::WriteRecPoints(Int_t event)
{

  // Creates new branches with given title
  // fills and writes into TreeR.

  AliEMCALGetter *gime = AliEMCALGetter::GetInstance() ; 

  TObjArray * aPRERecPoints = gime->PRERecPoints() ; 
  TObjArray * aECRecPoints = gime->ECALRecPoints() ; 
  TObjArray * aHCRecPoints = gime->HCALRecPoints() ; 

  TClonesArray * digits = gime->Digits() ; 
  TTree * treeR ; 
  
  if(fToSplit){
    if(!fSplitFile)
      return ;
    fSplitFile->cd() ;
    TString name("TreeR") ;
    name += event ; 
    treeR = dynamic_cast<TTree*>(fSplitFile->Get(name)); 
  }
  else{
    treeR = gAlice->TreeR();
  }

  if(!treeR){
    gAlice->MakeTree("R", fSplitFile);
    treeR = gAlice->TreeR() ;
  }
 
  Int_t index ;

  //Evaluate position, dispersion and other RecPoint properties for PRE section
  for(index = 0; index < aPRERecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALRecPoint *>(aPRERecPoints->At(index)))->EvalAll(fPREW0,digits)  ;
  aPRERecPoints->Sort() ;
  
  for(index = 0; index < aPRERecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALRecPoint *>(aPRERecPoints->At(index)))->SetIndexInList(index) ;
  
  aPRERecPoints->Expand(aPRERecPoints->GetEntriesFast()) ;
  
  //Evaluate position, dispersion and other RecPoint properties for EC section
  for(index = 0; index < aECRecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALTowerRecPoint *>(aECRecPoints->At(index)))->EvalAll(fECW0,digits) ;
  
  aECRecPoints->Sort() ;

  for(index = 0; index < aECRecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALTowerRecPoint *>(aECRecPoints->At(index)))->SetIndexInList(index) ;

  aECRecPoints->Expand(aECRecPoints->GetEntriesFast()) ; 
  
  //Evaluate position, dispersion and other RecPoint properties for HC section
  for(index = 0; index < aHCRecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALTowerRecPoint *>(aHCRecPoints->At(index)))->EvalAll(fHCW0,digits) ;
  
  aHCRecPoints->Sort() ;

  for(index = 0; index < aHCRecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALTowerRecPoint *>(aHCRecPoints->At(index)))->SetIndexInList(index) ;

  aHCRecPoints->Expand(aHCRecPoints->GetEntriesFast()) ; 
 
  Int_t bufferSize = 32000 ;    
  Int_t splitlevel = 0 ; 
  
  //PRE section branch 
  TBranch * branchPRE = treeR->Branch("EMCALPRERP","TObjArray",&aPRERecPoints,bufferSize,splitlevel);
  branchPRE->SetTitle(BranchName());

  //EC section branch
  TBranch * branchEC = treeR->Branch("EMCALECRP","TObjArray",&aECRecPoints,bufferSize,splitlevel);
  branchEC->SetTitle(BranchName());

  //HC section branch
  TBranch * branchHC = treeR->Branch("EMCALHCRP","TObjArray",&aHCRecPoints,bufferSize,splitlevel);
  branchHC->SetTitle(BranchName());
    
  //And Finally  clusterizer branch
  AliEMCALClusterizerv1 * cl = (AliEMCALClusterizerv1*)gime->Clusterizer(BranchName()) ;
  TBranch * clusterizerBranch = treeR->Branch("AliEMCALClusterizer","AliEMCALClusterizerv1",
					      &cl,bufferSize,splitlevel);
  clusterizerBranch->SetTitle(BranchName());

  branchPRE->Fill() ;
  branchEC->Fill() ;
  branchHC->Fill() ;
  clusterizerBranch->Fill() ;

  treeR->AutoSave() ; //Write(0,kOverwrite) ;  
  if(gAlice->TreeR()!=treeR)
    treeR->Delete(); 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
    
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 

  AliEMCALGeometry * geom = gime->EMCALGeometry() ; 


  TObjArray * aPRERecPoints = gime->PRERecPoints(BranchName()) ; 
  TObjArray * aECRecPoints  = gime->ECALRecPoints(BranchName()) ; 
  TObjArray * aHCRecPoints  = gime->HCALRecPoints(BranchName()) ; 

  aPRERecPoints->Delete() ;
  aECRecPoints->Delete() ;
  aHCRecPoints->Delete() ;

  TClonesArray * digits = gime->Digits() ; 
  if ( !digits ) {
    Fatal("MakeClusters -> Digits with name %s not found", GetName() ) ; 
  } 

  TIter next(digits) ; 
  AliEMCALDigit * digit ; 
  Int_t ndigEC=0, ndigPRE=0, ndigHC=0 ; 

  // count the number of digits in EC section
  while ( (digit = dynamic_cast<AliEMCALDigit *>(next())) ) { // scan over the list of digits 
    if (geom->IsInECAL(digit->GetId())) 
      ndigEC++ ; 
    else if (geom->IsInPRE(digit->GetId()))
      ndigPRE++; 
    else if (geom->IsInHCAL(digit->GetId()))
      ndigHC++;
    else {
      Error("MakeClusters", "id = %d is a wrong ID!", digit->GetId()) ; 
      abort() ;
    }
  }

  // add amplitude of PRE and EC sections
  Int_t digEC ; 
  for (digEC = 0 ; digEC < ndigEC ; digEC++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(digEC)) ;
    Int_t digPRE ;
    for (digPRE = ndigEC ; digPRE < ndigEC+ndigPRE ; digPRE++) {
      AliEMCALDigit *  digitPRE = dynamic_cast<AliEMCALDigit *>(digits->At(digPRE)) ;
      if ( geom->AreInSameTower(digit->GetId(), digitPRE->GetId()) ){
	Float_t  amp = static_cast<Float_t>(digit->GetAmp()) + geom->GetSummationFraction() * static_cast<Float_t>(digitPRE->GetAmp()) + 0.5 ; 
	digit->SetAmp(static_cast<Int_t>(amp)) ; 
	break ; 
      }
    }
    if (gDebug) 
      Info("MakeClusters","id = %d amp = %d", digit->GetId(), digit->GetAmp()) ; 
  }  

  TClonesArray * digitsC =  dynamic_cast<TClonesArray*>(digits->Clone()) ;
  
  
  // Clusterization starts  
  
  TIter nextdigit(digitsC) ; 
  Bool_t notremovedEC = kTRUE, notremovedPRE = kTRUE ;
  
  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigit())) ) { // scan over the list of digitsC
    AliEMCALRecPoint * clu = 0 ; 
    
    TArrayI clusterPREdigitslist(50), clusterECdigitslist(50), clusterHCdigitslist(50);   
 
    Bool_t inPRE = kFALSE, inECAL = kFALSE, inHCAL = kFALSE ;
    if( geom->IsInPRE(digit->GetId()) ) {
      inPRE = kTRUE ; 
    }
    else if( geom->IsInECAL(digit->GetId()) ) {
      inECAL = kTRUE ;
    }
    else if( geom->IsInHCAL(digit->GetId()) ) {
      inHCAL = kTRUE ;
    }
    
    if (gDebug == 2) { 
      if (inPRE)
	Info("MakeClusters","id = %d, ene = %f , thre = %f ", 
	     digit->GetId(),Calibrate(digit->GetAmp(), 0), fPREClusteringThreshold) ;  
      if (inECAL)
	Info("MakeClusters","id = %d, ene = %f , thre = %f", 
	     digit->GetId(),Calibrate(digit->GetAmp(), 1), fECClusteringThreshold) ;  
      if (inHCAL)
	Info("MakeClusters","id = %d, ene = %f , thre = %f", 
	     digit->GetId(),Calibrate(digit->GetAmp(), 2), fHCClusteringThreshold ) ;  
    }
    
    if ( (inPRE  && (Calibrate(digit->GetAmp(), 0) > fPREClusteringThreshold  )) || 
	 (inECAL && (Calibrate(digit->GetAmp(), 1) > fECClusteringThreshold  ))  || 
	 (inHCAL && (Calibrate(digit->GetAmp(), 2) > fHCClusteringThreshold  )) ) {
      
      Int_t  iDigitInPRECluster = 0, iDigitInECCluster = 0, iDigitInHCCluster = 0; 
      Int_t where ; // PRE = 0, ECAl = 1, HCAL = 2

      // Find the seed in each of the section ECAL/PRE/HCAL

      if( geom->IsInECAL(digit->GetId()) ) {   
	where = 1 ; // to tell we are in ECAL
	// start a new Tower RecPoint
	if(fNumberOfECClusters >= aECRecPoints->GetSize()) 
	  aECRecPoints->Expand(2*fNumberOfECClusters+1) ;
	AliEMCALTowerRecPoint * rp = new  AliEMCALTowerRecPoint("") ; 
	rp->SetECAL() ; 
	aECRecPoints->AddAt(rp, fNumberOfECClusters) ;
	clu = dynamic_cast<AliEMCALTowerRecPoint *>(aECRecPoints->At(fNumberOfECClusters)) ; 
  	fNumberOfECClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp(), where)) ; 
	clusterECdigitslist[iDigitInECCluster] = digit->GetIndexInList() ;	
	iDigitInECCluster++ ; 
	digitsC->Remove(digit) ; 
	if (gDebug == 2 ) 
	  Info("MakeClusters","OK id = %d, ene = %f , thre = %f ", digit->GetId(),Calibrate(digit->GetAmp(), 1), fECClusteringThreshold) ;  
	
      } 
      else if( geom->IsInPRE(digit->GetId()) ) { 
	where = 0 ; // to tell we are in PRE
	// start a new Pre Shower cluster
	if(fNumberOfPREClusters >= aPRERecPoints->GetSize()) 
	  aPRERecPoints->Expand(2*fNumberOfPREClusters+1);
	AliEMCALTowerRecPoint * rp = new AliEMCALTowerRecPoint("") ;	
	rp->SetPRE() ; 
	aPRERecPoints->AddAt(rp, fNumberOfPREClusters) ;
	clu =  dynamic_cast<AliEMCALTowerRecPoint *>(aPRERecPoints->At(fNumberOfPREClusters))  ;  
	fNumberOfPREClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp(), where));	
	clusterPREdigitslist[iDigitInPRECluster] = digit->GetIndexInList()  ;	
	iDigitInPRECluster++ ; 
	digitsC->Remove(digit) ;
	if (gDebug == 2 ) 
	  Info("MakeClusters","OK id = %d, ene = %f , thre = %f ", digit->GetId(),Calibrate(digit->GetAmp(), 0), fPREClusteringThreshold) ;  
	
	nextdigit.Reset() ;

	// Here we remove remaining EC digits, which cannot make a cluster
	
	if( notremovedEC ) { 
	  while( ( digit = dynamic_cast<AliEMCALDigit *>(nextdigit()) ) ) {
	    if( geom->IsInECAL(digit->GetId()) )
	      digitsC->Remove(digit) ;
	    else 
	      break ; 
	  }
	  notremovedEC = kFALSE ;
	}

      } 
      else if( geom->IsInHCAL(digit->GetId()) ) { 
	where = 2 ; // to tell we are in HCAL
	// start a new HCAL cluster
	if(fNumberOfHCClusters >= aHCRecPoints->GetSize()) 
	  aHCRecPoints->Expand(2*fNumberOfHCClusters+1);
	AliEMCALTowerRecPoint * rp = new AliEMCALTowerRecPoint("") ;	
	rp->SetHCAL() ; 
	aHCRecPoints->AddAt(rp, fNumberOfHCClusters) ;
	clu =  dynamic_cast<AliEMCALTowerRecPoint *>(aHCRecPoints->At(fNumberOfHCClusters))  ;  
	fNumberOfHCClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp(), where));	
	clusterHCdigitslist[iDigitInHCCluster] = digit->GetIndexInList()  ;	
	iDigitInHCCluster++ ; 
	digitsC->Remove(digit) ;
	if (gDebug == 2 ) 
	  Info("MakeClusters","OK id = %d, ene = %f , thre = %f ", digit->GetId(),Calibrate(digit->GetAmp(), 2), fHCClusteringThreshold) ;  
 
	nextdigit.Reset() ;
   
	// Here we remove remaining PRE digits, which cannot make a cluster
	
	if( notremovedPRE ) { 
	  while( ( digit = dynamic_cast<AliEMCALDigit *>(nextdigit()) ) ) {
	    if( geom->IsInPRE(digit->GetId()) )
	      digitsC->Remove(digit) ;
	    else 
	      break ; 
	  }
	  notremovedPRE = kFALSE ;
	}	
      }        
      
      nextdigit.Reset() ;
      
      AliEMCALDigit * digitN ; 
      Int_t index = 0 ;

      // Do the Clustering in each of the three section ECAL/PRE/HCAL

      while (index < iDigitInECCluster){ // scan over digits already in cluster 
	digit =  (AliEMCALDigit*)digits->At(clusterECdigitslist[index])  ;      
	index++ ; 
        while ( (digitN = (AliEMCALDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
	  // Info("MakeClusters","id1 = %d, id2 = %d , neighbours = %d", digit->GetId(), digitN->GetId(), ineb) ;  
         switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;
	  case 1 :   // are neighbours 
	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp(), 1) ) ;
	    clusterECdigitslist[iDigitInECCluster] = digitN->GetIndexInList() ; 
	    iDigitInECCluster++ ; 
	    digitsC->Remove(digitN) ;
	    break ;
          case 2 :   // too far from each other
	    goto endofloop1;   
	  } // switch
	  
	} // while digitN
	
      endofloop1: ;
	nextdigit.Reset() ; 
      } // loop over EC cluster
      
      index = 0 ; 
      while (index < iDigitInPRECluster){ // scan over digits already in cluster 
	digit =  (AliEMCALDigit*)digits->At(clusterPREdigitslist[index])  ;      
	index++ ; 
        while ( (digitN = (AliEMCALDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
	  //	  Info("MakeClusters","id1 = %d, id2 = %d , neighbours = %d", digit->GetId(), digitN->GetId(), ineb) ;  
	  switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;
	  case 1 :   // are neighbours 
	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp(), 0) ) ;
	    clusterPREdigitslist[iDigitInPRECluster] = digitN->GetIndexInList() ; 
	    iDigitInPRECluster++ ; 
	    digitsC->Remove(digitN) ;
	    break ;
          case 2 :   // too far from each other
	    goto endofloop2;   
	  } // switch
	  
	} // while digitN
	
      endofloop2: ;
       	nextdigit.Reset() ; 
      } // loop over PRE cluster
    
      index = 0 ; 
      while (index < iDigitInHCCluster){ // scan over digits already in cluster 
	digit =  (AliEMCALDigit*)digits->At(clusterHCdigitslist[index])  ;      
	index++ ; 
        while ( (digitN = (AliEMCALDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
	  //Info("MakeClusters","id1 = %d, id2 = %d , neighbours = %d", digit->GetId(), digitN->GetId(), ineb) ;  
	  switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;
	  case 1 :   // are neighbours 
	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp(), 2) ) ;
	    clusterHCdigitslist[iDigitInHCCluster] = digitN->GetIndexInList() ; 
	    iDigitInHCCluster++ ; 
	    digitsC->Remove(digitN) ;
	    break ;
          case 2 :   // too far from each other
	    goto endofloop3;   
	  } // switch  
	} // while digitN
	
      endofloop3: ;
	nextdigit.Reset() ; 
      } // loop over HC cluster

    } // energy theshold     
  } // while digit  
  delete digitsC ;
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeUnfolding()
{
  Fatal("AliEMCALClusterizerv1::MakeUnfolding", "--> Unfolding not implemented") ;
 
}

//____________________________________________________________________________
Double_t  AliEMCALClusterizerv1::ShowerShape(Double_t r)
{ 
  // Shape of the shower (see EMCAL TDR)
  // If you change this function, change also the gradient evaluation in ChiSquare()

  Double_t r4    = r*r*r*r ;
  Double_t r295  = TMath::Power(r, 2.95) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliEMCALClusterizerv1::UnfoldCluster(AliEMCALTowerRecPoint * iniTower, 
						 Int_t nMax, 
						 AliEMCALDigit ** maxAt, 
						 Float_t * maxAtEnergy)
{
  // Performs the unfolding of a cluster with nMax overlapping showers 
  
  Fatal("UnfoldCluster", "--> Unfolding not implemented") ;

}

//_____________________________________________________________________________
void AliEMCALClusterizerv1::UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do
  
  ::Fatal("UnfoldingChiSquare","Unfolding not implemented") ;
}
//____________________________________________________________________________
void AliEMCALClusterizerv1::Print(Option_t * option)const
{
  // Print clusterizer parameters

  TString message("\n") ; 
  
  if( strcmp(GetName(), "") !=0 ){
    
    // Print parameters
 
    TString taskName(GetName()) ; 
    taskName.ReplaceAll(Version(), "") ;
    
    message += "--------------- " ; 
    message += taskName.Data() ; 
    message += " " ; 
    message += GetTitle() ; 
    message += "-----------\n" ;  
    message += "Clusterizing digits from the file: " ; 
    message += taskName.Data() ;  
    message += "\n                           Branch: " ; 
    message += GetName() ;  
    message += "\n                       Pre Shower Clustering threshold = " ; 
    message += fPREClusteringThreshold ;
    message += "\n                       Pre Shower  Local Maximum cut    = " ;
    message += fPRELocMaxCut ;
    message += "\n                       Pre Shower Logarothmic weight   = " ; 
    message += fPREW0 ;
    message += "\n                       EC Clustering threshold = " ; 
    message += fECClusteringThreshold ; 
    message += "\n                       EC Local Maximum cut    = " ;
    message += fECLocMaxCut ; 
    message += "\n                       EC Logarothmic weight   = " ;
    message += fECW0 ;
    message += "\n                       Pre Shower Clustering threshold = " ; 
    message += fHCClusteringThreshold ; 
    message += "\n                       HC Local Maximum cut    = " ;
    message += fHCLocMaxCut ; 
    message += "\n                       HC Logarothmic weight   = " ;
    message += fHCW0 ;
    if(fToUnfold)
      message +="\nUnfolding on\n" ;
    else
      message += "\nUnfolding off\n";
    
    message += "------------------------------------------------------------------" ; 
  }
  else
    message += "AliEMCALClusterizerv1 not initialized " ;
  
  Info("Print", message.Data() ) ; 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::PrintRecPoints(Option_t * option)
{
  // Prints list of RecPoints produced at the current pass of AliEMCALClusterizer

  TObjArray * aPRERecPoints = AliEMCALGetter::GetInstance()->PRERecPoints() ; 
  TObjArray * aECRecPoints = AliEMCALGetter::GetInstance()->ECALRecPoints() ; 
  TObjArray * aHCRecPoints = AliEMCALGetter::GetInstance()->HCALRecPoints() ; 

  Info("PrintRecPoints", "Clusterization result:") ; 
  
  printf("event # %d\n", gAlice->GetEvNumber() ) ;
  printf("           Found %d PRE SHOWER RecPoints, %d EC Rec Points and %d HC Rec Points\n ", 
	 aPRERecPoints->GetEntriesFast(), aECRecPoints->GetEntriesFast(), aHCRecPoints->GetEntriesFast() ) ; 

  fRecPointsInRun +=  aPRERecPoints->GetEntriesFast() ; 
  fRecPointsInRun +=  aECRecPoints->GetEntriesFast() ; 
  fRecPointsInRun +=  aHCRecPoints->GetEntriesFast() ; 
  
  if(strstr(option,"all")) {

    //Pre shower recPoints

    printf("-----------------------------------------------------------------------\n") ;
    printf("Clusters in PRE section\n") ;
    printf("Index    Ene(GeV) Multi Module     phi     r   theta    X    Y      Z   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;      

    Int_t index ;
    
    for (index = 0 ; index < aPRERecPoints->GetEntries() ; index++) {
      AliEMCALTowerRecPoint * rp = dynamic_cast<AliEMCALTowerRecPoint *>(aPRERecPoints->At(index)) ; 
      TVector3  globalpos;  
      rp->GetGlobalPosition(globalpos);
      TVector3  localpos;  
      rp->GetLocalPosition(localpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      Int_t * primaries;
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      printf("\n%6d  %8.4f  %3d     %2d     %4.1f    %4.1f %4.1f  %4.1f %4.1f %4.1f    %4.1f   %4f  %4f    %2d     : ", 
	     rp->GetIndexInList(), rp->GetEnergy(), rp->GetMultiplicity(), rp->GetEMCALArm(), 
	     globalpos.X(), globalpos.Y(), globalpos.Z(), localpos.X(), localpos.Y(), localpos.Z(), 
	     rp->GetDispersion(), lambda[0], lambda[1], nprimaries) ; 
      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	printf("%d ", primaries[iprimary] ) ; 
      }  	 
    }
    
    printf("\n-----------------------------------------------------------------------\n") ;
    printf("Clusters in ECAL section\n") ;
    printf("Index    Ene(GeV) Multi Module     phi     r   theta    X    Y      Z   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;      
    
    for (index = 0 ; index < aECRecPoints->GetEntries() ; index++) {
      AliEMCALTowerRecPoint * rp = dynamic_cast<AliEMCALTowerRecPoint * >(aECRecPoints->At(index)) ; 
      TVector3  globalpos;  
      rp->GetGlobalPosition(globalpos);
      TVector3  localpos;  
      rp->GetLocalPosition(localpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      printf("\n%6d  %8.4f  %3d     %2d     %4.1f    %4.1f %4.1f  %4.1f %4.1f %4.1f    %4.1f   %4f  %4f    %2d     : ", 
	     rp->GetIndexInList(), rp->GetEnergy(), rp->GetMultiplicity(), rp->GetEMCALArm(), 
	     globalpos.X(), globalpos.Y(), globalpos.Z(), localpos.X(), localpos.Y(), localpos.Z(), 
	     rp->GetDispersion(), lambda[0], lambda[1], nprimaries) ; 
      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	printf("%d ", primaries[iprimary] ) ; 
      } 
    }

    printf("\n-----------------------------------------------------------------------\n") ;
    printf("Clusters in HCAL section\n") ;
    printf("Index    Ene(GeV) Multi Module     phi     r   theta    X    Y      Z   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;      
    
    for (index = 0 ; index < aHCRecPoints->GetEntries() ; index++) {
      AliEMCALTowerRecPoint * rp = dynamic_cast<AliEMCALTowerRecPoint * >(aHCRecPoints->At(index)) ; 
      TVector3  globalpos;  
      rp->GetGlobalPosition(globalpos);
      TVector3  localpos;  
      rp->GetLocalPosition(localpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      printf("\n%6d  %8.4f  %3d     %2d     %4.1f    %4.1f %4.1f  %4.1f %4.1f %4.1f    %4.1f   %4f  %4f    %2d     : ", 
	     rp->GetIndexInList(), rp->GetEnergy(), rp->GetMultiplicity(), rp->GetEMCALArm(), 
	     globalpos.X(), globalpos.Y(), globalpos.Z(), localpos.X(), localpos.Y(), localpos.Z(), 
	     rp->GetDispersion(), lambda[0], lambda[1], nprimaries) ;      
      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	printf("%d ", primaries[iprimary] ) ; 
      } 
    }

    printf("\n-----------------------------------------------------------------------\n");
  }
}
