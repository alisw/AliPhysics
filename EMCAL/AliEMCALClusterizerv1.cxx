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



#include <Riostream.h>



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

Float_t  AliEMCALClusterizerv1::Calibrate(Int_t amp, Bool_t inpresho) const

{//To be replased later by the method, reading individual parameters from the database



  if ( inpresho ) // calibrate as pre shower

     return -fADCpedestalPreSho + amp * fADCchannelPreSho ; 



  else //calibrate as tower 

    return -fADCpedestalTower + amp * fADCchannelTower ;                

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



    fNumberOfTowerClusters = fNumberOfPreShoClusters = 0 ;

           

    MakeClusters() ;

    

    if(fToUnfold)

      MakeUnfolding() ;



    WriteRecPoints(ievent) ;



    if(strstr(option,"deb"))  

      PrintRecPoints(option) ;



    //increment the total number of digits per run 

    fRecPointsInRun += gime->TowerRecPoints()->GetEntriesFast() ;  

    fRecPointsInRun += gime->PreShowerRecPoints()->GetEntriesFast() ;  

 }

  

  if(strstr(option,"tim")){

    gBenchmark->Stop("EMCALClusterizer");

    cout << "AliEMCALClusterizer:" << endl ;

    cout << "  took " << gBenchmark->GetCpuTime("EMCALClusterizer") << " seconds for Clusterizing " 

	 <<  gBenchmark->GetCpuTime("EMCALClusterizer")/nevents << " seconds per event " << endl ;

    cout << endl ;

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

      cout << "EMCAL Unfolding>  Unable to set initial value for fit procedure : x = " << x << endl ;

      return kFALSE;

    }

    gMinuit->mnparm(index, "z",  z, 0.1, 0, 0, ierflg) ;

    index++ ;   

    if(ierflg != 0){

      cout << "EMCAL Unfolding>  Unable to set initial value for fit procedure : z = " << z << endl ;

      return kFALSE;

    }

    gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;

    index++ ;   

    if(ierflg != 0){

      cout << "EMCAL Unfolding>  Unable to set initial value for fit procedure : energy = " << energy << endl ;      

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

    cout << "EMCAL Unfolding>  Fit not converged, cluster abandoned "<< endl ;      

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



  fADCchannelTower   = dig->GetTowerchannel() ;

  fADCpedestalTower  = dig->GetTowerpedestal();



  fADCchannelPreSho  = dig->GetPreShochannel() ;

  fADCpedestalPreSho = dig->GetPreShopedestal() ; 



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

    cerr << "ERROR: AliEMCALClusterizerv1::Init -> Could not obtain the Getter object !" << endl ; 

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

  fNumberOfPreShoClusters = fNumberOfTowerClusters = 0 ; 



 

  

  fPreShoClusteringThreshold  = 0.0001;

  fTowerClusteringThreshold   = 0.2;   

  

  fTowerLocMaxCut  = 0.03 ;

  fPreShoLocMaxCut = 0.03 ;

  

  fW0     = 4.5 ;

  fW0CPV  = 4.0 ;



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

 

  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same EMCAL Arm 

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

      rv=2 ;



  }



  return rv ; 

}





//____________________________________________________________________________

Bool_t AliEMCALClusterizerv1::IsInTower(AliEMCALDigit * digit) const

{

  // Tells if (true) or not (false) the digit is in a EMCAL-Tower 

 

  Bool_t rv = kFALSE ; 

  if (!digit->IsInPreShower()) 

    rv = kTRUE; 

  return rv ; 

}



//____________________________________________________________________________

Bool_t AliEMCALClusterizerv1::IsInPreShower(AliEMCALDigit * digit) const

{

  // Tells if (true) or not (false) the digit is in a EMCAL-PreShower

 

  Bool_t rv = kFALSE ; 

  if (digit->IsInPreShower()) 

    rv = kTRUE; 

  return rv ; 

}



//____________________________________________________________________________

void AliEMCALClusterizerv1::WriteRecPoints(Int_t event)

{



  // Creates new branches with given title

  // fills and writes into TreeR.



  AliEMCALGetter *gime = AliEMCALGetter::GetInstance() ; 

  TObjArray * towerRecPoints = gime->TowerRecPoints() ; 

  TObjArray * preshoRecPoints = gime->PreShowerRecPoints() ; 

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

  //Evaluate position, dispersion and other RecPoint properties...

  for(index = 0; index < towerRecPoints->GetEntries(); index++)

    (dynamic_cast<AliEMCALTowerRecPoint *>(towerRecPoints->At(index)))->EvalAll(fW0,digits) ;

  

  towerRecPoints->Sort() ;



  for(index = 0; index < towerRecPoints->GetEntries(); index++)

    (dynamic_cast<AliEMCALTowerRecPoint *>(towerRecPoints->At(index)))->SetIndexInList(index) ;



  towerRecPoints->Expand(towerRecPoints->GetEntriesFast()) ; 



  //Now the same for pre shower

  for(index = 0; index < preshoRecPoints->GetEntries(); index++)

    (dynamic_cast<AliEMCALRecPoint *>(preshoRecPoints->At(index)))->EvalAll(fW0CPV,digits)  ;



  preshoRecPoints->Sort() ;



  for(index = 0; index < preshoRecPoints->GetEntries(); index++)

    (dynamic_cast<AliEMCALRecPoint *>(preshoRecPoints->At(index)))->SetIndexInList(index) ;



  preshoRecPoints->Expand(preshoRecPoints->GetEntriesFast()) ;

  

  Int_t bufferSize = 32000 ;    

  Int_t splitlevel = 0 ;



  //First Tower branch

  TBranch * towerBranch = treeR->Branch("EMCALTowerRP","TObjArray",&towerRecPoints,bufferSize,splitlevel);

  towerBranch->SetTitle(BranchName());

  

  //Now Pre Shower branch 

  TBranch * preshoBranch = treeR->Branch("EMCALPreShoRP","TObjArray",&preshoRecPoints,bufferSize,splitlevel);

  preshoBranch->SetTitle(BranchName());

    

  //And Finally  clusterizer branch

  AliEMCALClusterizerv1 * cl = (AliEMCALClusterizerv1*)gime->Clusterizer(BranchName()) ;

  TBranch * clusterizerBranch = treeR->Branch("AliEMCALClusterizer","AliEMCALClusterizerv1",

					      &cl,bufferSize,splitlevel);

  clusterizerBranch->SetTitle(BranchName());



  towerBranch        ->Fill() ;

  preshoBranch        ->Fill() ;

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

  

  TObjArray * towerRecPoints  = gime->TowerRecPoints(BranchName()) ; 

  TObjArray * preshoRecPoints = gime->PreShowerRecPoints(BranchName()) ; 

  towerRecPoints->Delete() ;

  preshoRecPoints->Delete() ;

  

  TClonesArray * digits = gime->Digits() ; 

  if ( !digits ) {

    cerr << "ERROR:  AliEMCALClusterizerv1::MakeClusters -> Digits with name " 

	 << GetName() << " not found ! " << endl ; 

    abort() ; 

  } 

  TClonesArray * digitsC =  dynamic_cast<TClonesArray*>(digits->Clone()) ;

  

  

  // Clusterization starts  

  

  TIter nextdigit(digitsC) ; 

  AliEMCALDigit * digit ; 

  Bool_t notremoved = kTRUE ;

  

  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigit())) ) { // scan over the list of digitsC

    AliEMCALRecPoint * clu = 0 ; 

    

    TArrayI clusterdigitslist(1500) ;   

    Int_t index ;

 

    if (( IsInTower (digit)  && Calibrate(digit->GetAmp(),digit->IsInPreShower()) > fTowerClusteringThreshold  ) || 

        ( IsInPreShower (digit) && Calibrate(digit->GetAmp(),digit->IsInPreShower()) > fPreShoClusteringThreshold  ) ) {

      

      Int_t iDigitInCluster = 0 ; 

      

      if  ( IsInTower(digit) ) {   

	// start a new Tower RecPoint

	if(fNumberOfTowerClusters >= towerRecPoints->GetSize()) 

	  towerRecPoints->Expand(2*fNumberOfTowerClusters+1) ;

  	

	towerRecPoints->AddAt(new  AliEMCALTowerRecPoint(""), fNumberOfTowerClusters) ;

	clu = dynamic_cast<AliEMCALTowerRecPoint *>(towerRecPoints->At(fNumberOfTowerClusters)) ; 

  	fNumberOfTowerClusters++ ; 

	clu->AddDigit(*digit, Calibrate(digit->GetAmp(),digit->IsInPreShower())) ; 

	clusterdigitslist[iDigitInCluster] = digit->GetIndexInList() ;	

	iDigitInCluster++ ; 

	digitsC->Remove(digit) ; 

	

      } else { 

	

	// start a new Pre Shower cluster

	if(fNumberOfPreShoClusters >= preshoRecPoints->GetSize()) 

	  preshoRecPoints->Expand(2*fNumberOfPreShoClusters+1);

	

	preshoRecPoints->AddAt(new AliEMCALTowerRecPoint(""), fNumberOfPreShoClusters) ;

	

	clu =  dynamic_cast<AliEMCALTowerRecPoint *>(preshoRecPoints->At(fNumberOfPreShoClusters))  ;  

	fNumberOfPreShoClusters++ ; 

	clu->AddDigit(*digit, Calibrate(digit->GetAmp(),digit->IsInPreShower() ) );	

	clusterdigitslist[iDigitInCluster] = digit->GetIndexInList()  ;	

	iDigitInCluster++ ; 

	digitsC->Remove(digit) ; 

	nextdigit.Reset() ;

	

	// Here we remove remaining Tower digits, which cannot make a cluster

	

	if( notremoved ) { 

	  while( ( digit = dynamic_cast<AliEMCALDigit *>(nextdigit()) ) ) {

	    if( IsInTower(digit) )

	      digitsC->Remove(digit) ;

	    else 

	      break ; 

	  }

	  notremoved = kFALSE ;

	}

	

      } // else        

      

      nextdigit.Reset() ;

      

      AliEMCALDigit * digitN ; 

      index = 0 ;

      while (index < iDigitInCluster){ // scan over digits already in cluster 

	digit =  (AliEMCALDigit*)digits->At(clusterdigitslist[index])  ;      

	index++ ; 

        while ( (digitN = (AliEMCALDigit *)nextdigit()) ) { // scan over the reduced list of digits 

	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!

         switch (ineb ) {

          case 0 :   // not a neighbour

	    break ;

	  case 1 :   // are neighbours 

	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp(), digitN->IsInPreShower() ) ) ;

	    clusterdigitslist[iDigitInCluster] = digitN->GetIndexInList() ; 

	    iDigitInCluster++ ; 

	    digitsC->Remove(digitN) ;

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

  

  delete digitsC ;

  

}



//____________________________________________________________________________

void AliEMCALClusterizerv1::MakeUnfolding()

{

  Fatal("AliEMCALClusterizerv1::MakeUnfolding", "--> Unfolding not implemented") ;

  

//   // Unfolds clusters using the shape of an ElectroMagnetic shower

//   // Performs unfolding of all EMC/CPV clusters



//   AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 

  

//   const AliEMCALGeometry * geom = gime->EMCALGeometry() ;

//   TObjArray * emcRecPoints = gime->TowerRecPoints() ; 

//   TObjArray * cpvRecPoints = gime->PreShoRecPoints() ; 

//   TClonesArray * digits = gime->Digits() ; 

  

//   // Unfold first EMC clusters 

//   if(fNumberOfTowerClusters > 0){



//     Int_t nModulesToUnfold = geom->GetNModules() ; 



//     Int_t numberofNotUnfolded = fNumberOfTowerClusters ; 

//     Int_t index ;   

//     for(index = 0 ; index < numberofNotUnfolded ; index++){

      

//       AliEMCALTowerRecPoint * emcRecPoint = (AliEMCALTowerRecPoint *) emcRecPoints->At(index) ;

//       if(emcRecPoint->GetEMCALMod()> nModulesToUnfold)

// 	break ;

      

//       Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 

//       Int_t * maxAt = new Int_t[nMultipl] ;

//       Float_t * maxAtEnergy = new Float_t[nMultipl] ;

//       Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fTowerLocMaxCut,digits) ;

      

//       if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       

// 	UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;

// 	emcRecPoints->Remove(emcRecPoint); 

// 	emcRecPoints->Compress() ;

// 	index-- ;

// 	fNumberOfTowerClusters -- ;

// 	numberofNotUnfolded-- ;

//       }

      

//       delete[] maxAt ; 

//       delete[] maxAtEnergy ; 

//     }

//   } 

//   // Unfolding of EMC clusters finished





//   // Unfold now CPV clusters

//   if(fNumberOfPreShoClusters > 0){

    

//     Int_t nModulesToUnfold = geom->GetNModules() ;



//     Int_t numberofPreShoNotUnfolded = fNumberOfPreShoClusters ;     

//     Int_t index ;   

//     for(index = 0 ; index < numberofPreShoNotUnfolded ; index++){

      

//       AliEMCALRecPoint * recPoint = (AliEMCALRecPoint *) cpvRecPoints->At(index) ;



//       if(recPoint->GetEMCALMod()> nModulesToUnfold)

// 	break ;

      

//       AliEMCALTowerRecPoint * emcRecPoint = (AliEMCALTowerRecPoint*) recPoint ; 

      

//       Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 

//       Int_t * maxAt = new Int_t[nMultipl] ;

//       Float_t * maxAtEnergy = new Float_t[nMultipl] ;

//       Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fPreShoLocMaxCut,digits) ;

      

//       if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       

// 	UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;

// 	cpvRecPoints->Remove(emcRecPoint); 

// 	cpvRecPoints->Compress() ;

// 	index-- ;

// 	numberofPreShoNotUnfolded-- ;

// 	fNumberOfPreShoClusters-- ;

//       }

      

//       delete[] maxAt ; 

//       delete[] maxAtEnergy ; 

//     } 

//   }

//   //Unfolding of PreSho clusters finished

  

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

  

  Fatal("AliEMCALClusterizerv1::UnfoldCluster", "--> Unfolding not implemented") ;



 //  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 

//   const AliEMCALGeometry * geom = gime->EMCALGeometry() ;

//   const TClonesArray * digits = gime->Digits() ; 

//   TObjArray * emcRecPoints = gime->TowerRecPoints() ; 

//   TObjArray * cpvRecPoints = gime->PreShoRecPoints() ; 



//   Int_t nPar = 3 * nMax ;

//   Float_t * fitparameters = new Float_t[nPar] ;



//   Bool_t rv = FindFit(iniTower, maxAt, maxAtEnergy, nPar, fitparameters) ;

//   if( !rv ) {

//     // Fit failed, return and remove cluster

//     delete[] fitparameters ; 

//     return ;

//   }



//   // create ufolded rec points and fill them with new energy lists

//   // First calculate energy deposited in each sell in accordance with fit (without fluctuations): efit[]

//   // and later correct this number in acordance with actual energy deposition



//   Int_t nDigits = iniTower->GetMultiplicity() ;  

//   Float_t * efit = new Float_t[nDigits] ;

//   Float_t xDigit=0.,zDigit=0.,distance=0. ;

//   Float_t xpar=0.,zpar=0.,epar=0.  ;

//   Int_t relid[4] ;

//   AliEMCALDigit * digit = 0 ;

//   Int_t * emcDigits = iniTower->GetDigitsList() ;



//   Int_t iparam ;

//   Int_t iDigit ;

//   for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){

//     digit = (AliEMCALDigit*) digits->At(emcDigits[iDigit] ) ;   

//     geom->AbsToRelNumbering(digit->GetId(), relid) ;

//     geom->RelPosInModule(relid, xDigit, zDigit) ;

//     efit[iDigit] = 0;



//     iparam = 0 ;    

//     while(iparam < nPar ){

//       xpar = fitparameters[iparam] ;

//       zpar = fitparameters[iparam+1] ;

//       epar = fitparameters[iparam+2] ;

//       iparam += 3 ;

//       distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;

//       distance =  TMath::Sqrt(distance) ;

//       efit[iDigit] += epar * ShowerShape(distance) ;

//     }

//   }

  



//   // Now create new RecPoints and fill energy lists with efit corrected to fluctuations

//   // so that energy deposited in each cell is distributed betwin new clusters proportionally

//   // to its contribution to efit



//   Float_t * emcEnergies = iniTower->GetEnergiesList() ;

//   Float_t ratio ;



//   iparam = 0 ;

//   while(iparam < nPar ){

//     xpar = fitparameters[iparam] ;

//     zpar = fitparameters[iparam+1] ;

//     epar = fitparameters[iparam+2] ;

//     iparam += 3 ;    

    

//     AliEMCALTowerRecPoint * emcRP = 0 ;  



//     if(iniTower->IsTower()){ //create new entries in fTowerRecPoints...

      

//       if(fNumberOfTowerClusters >= emcRecPoints->GetSize())

// 	emcRecPoints->Expand(2*fNumberOfTowerClusters) ;

      

//       (*emcRecPoints)[fNumberOfTowerClusters] = new AliEMCALTowerRecPoint("") ;

//       emcRP = (AliEMCALTowerRecPoint *) emcRecPoints->At(fNumberOfTowerClusters);

//       fNumberOfTowerClusters++ ;

//     }

//     else{//create new entries in fPreShoRecPoints

//       if(fNumberOfPreShoClusters >= cpvRecPoints->GetSize())

// 	cpvRecPoints->Expand(2*fNumberOfPreShoClusters) ;

      

//       (*cpvRecPoints)[fNumberOfPreShoClusters] = new AliEMCALPreShoRecPoint("") ;

//       emcRP = (AliEMCALTowerRecPoint *) cpvRecPoints->At(fNumberOfPreShoClusters);

//       fNumberOfPreShoClusters++ ;

//     }

    

//     Float_t eDigit ;

//     for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){

//       digit = (AliEMCALDigit*) digits->At( emcDigits[iDigit] ) ; 

//       geom->AbsToRelNumbering(digit->GetId(), relid) ;

//       geom->RelPosInModule(relid, xDigit, zDigit) ;

//       distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;

//       distance =  TMath::Sqrt(distance) ;

//       ratio = epar * ShowerShape(distance) / efit[iDigit] ; 

//       eDigit = emcEnergies[iDigit] * ratio ;

//       emcRP->AddDigit( *digit, eDigit ) ;

//     }	

//   }

 

//   delete[] fitparameters ; 

//   delete[] efit ; 

  

}



//_____________________________________________________________________________

void AliEMCALClusterizerv1::UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)

{

  // Calculates the Chi square for the cluster unfolding minimization

  // Number of parameters, Gradient, Chi squared, parameters, what to do



  abort() ; 

 //  Fatal("AliEMCALClusterizerv1::UnfoldingChiSquare","-->Unfolding not implemented") ;



//   TList * toMinuit = (TList*) gMinuit->GetObjectFit() ;



//   AliEMCALTowerRecPoint * emcRP = (AliEMCALTowerRecPoint*) toMinuit->At(0)  ;

//   TClonesArray * digits = (TClonesArray*)toMinuit->At(1)  ;





  

//   //  AliEMCALTowerRecPoint * emcRP = (AliEMCALTowerRecPoint *) gMinuit->GetObjectFit() ; // TowerRecPoint to fit



//   Int_t * emcDigits     = emcRP->GetDigitsList() ;



//   Int_t nOdigits = emcRP->GetDigitsMultiplicity() ; 



//   Float_t * emcEnergies = emcRP->GetEnergiesList() ;



//   const AliEMCALGeometry * geom = AliEMCALGetter::GetInstance()->EMCALGeometry() ; 

//   fret = 0. ;     

//   Int_t iparam ;



//   if(iflag == 2)

//     for(iparam = 0 ; iparam < nPar ; iparam++)    

//       Grad[iparam] = 0 ; // Will evaluate gradient

  

//   Double_t efit ;    



//   AliEMCALDigit * digit ;

//   Int_t iDigit ;



//   for( iDigit = 0 ; iDigit < nOdigits ; iDigit++) {



//     digit = (AliEMCALDigit*) digits->At( emcDigits[iDigit] ) ; 



//     Int_t relid[4] ;

//     Float_t xDigit ;

//     Float_t zDigit ;



//     geom->AbsToRelNumbering(digit->GetId(), relid) ;



//     geom->RelPosInModule(relid, xDigit, zDigit) ;



//      if(iflag == 2){  // calculate gradient

//        Int_t iParam = 0 ;

//        efit = 0 ;

//        while(iParam < nPar ){

// 	 Double_t distance = (xDigit - x[iParam]) * (xDigit - x[iParam]) ;

// 	 iParam++ ; 

// 	 distance += (zDigit - x[iParam]) * (zDigit - x[iParam]) ; 

// 	 distance = TMath::Sqrt( distance ) ; 

// 	 iParam++ ; 	 

// 	 efit += x[iParam] * ShowerShape(distance) ;

// 	 iParam++ ;

//        }

//        Double_t sum = 2. * (efit - emcEnergies[iDigit]) / emcEnergies[iDigit] ; // Here we assume, that sigma = sqrt(E) 

//        iParam = 0 ;

//        while(iParam < nPar ){

// 	 Double_t xpar = x[iParam] ;

// 	 Double_t zpar = x[iParam+1] ;

// 	 Double_t epar = x[iParam+2] ;

// 	 Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );

// 	 Double_t shape = sum * ShowerShape(dr) ;

// 	 Double_t r4 = dr*dr*dr*dr ;

// 	 Double_t r295 = TMath::Power(dr,2.95) ;

// 	 Double_t deriv =-4. * dr*dr * ( 2.32 / ( (2.32 + 0.26 * r4) * (2.32 + 0.26 * r4) ) +

// 					 0.0316 * (1. + 0.0171 * r295) / ( ( 1. + 0.0652 * r295) * (1. + 0.0652 * r295) ) ) ;

	 

// 	 Grad[iParam] += epar * shape * deriv * (xpar - xDigit) ;  // Derivative over x    

// 	 iParam++ ; 

// 	 Grad[iParam] += epar * shape * deriv * (zpar - zDigit) ;  // Derivative over z         

// 	 iParam++ ; 

// 	 Grad[iParam] += shape ;                                  // Derivative over energy     	

// 	 iParam++ ; 

//        }

//      }

//      efit = 0;

//      iparam = 0 ;



//      while(iparam < nPar ){

//        Double_t xpar = x[iparam] ;

//        Double_t zpar = x[iparam+1] ;

//        Double_t epar = x[iparam+2] ;

//        iparam += 3 ;

//        Double_t distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;

//        distance =  TMath::Sqrt(distance) ;

//        efit += epar * ShowerShape(distance) ;

//      }



//      fret += (efit-emcEnergies[iDigit])*(efit-emcEnergies[iDigit])/emcEnergies[iDigit] ; 

//      // Here we assume, that sigma = sqrt(E)

//   }



}



//____________________________________________________________________________

void AliEMCALClusterizerv1::Print(Option_t * option)const

{

  // Print clusterizer parameters



  if( strcmp(GetName(), "") !=0 ){

    

    // Print parameters

 

    TString taskName(GetName()) ; 

    taskName.ReplaceAll(Version(), "") ;



    cout << "---------------"<< taskName.Data() << " " << GetTitle()<< "-----------" << endl 

	 << "Clusterizing digits from the file: " << taskName.Data() << endl 

	 << "                           Branch: " << GetName() << endl 

	 << endl 

	 << "                       EMC Clustering threshold = " << fTowerClusteringThreshold << endl

	 << "                       EMC Local Maximum cut    = " << fTowerLocMaxCut << endl

	 << "                       EMC Logarothmic weight   = " << fW0 << endl

	 << endl

	 << "                       CPV Clustering threshold = " << fPreShoClusteringThreshold << endl

	 << "                       CPV Local Maximum cut    = " << fPreShoLocMaxCut << endl

       << "                       CPV Logarothmic weight   = " << fW0CPV << endl

	 << endl ;

    if(fToUnfold)

      cout << " Unfolding on " << endl ;

    else

      cout << " Unfolding off " << endl ;

    

    cout << "------------------------------------------------------------------" <<endl ;

  }

  else

    cout << " AliEMCALClusterizerv1 not initialized " << endl ;

}

//____________________________________________________________________________

void AliEMCALClusterizerv1::PrintRecPoints(Option_t * option)

{

  // Prints list of RecPoints produced at the current pass of AliEMCALClusterizer



  TObjArray * towerRecPoints = AliEMCALGetter::GetInstance()->TowerRecPoints() ; 

  TObjArray * preshoRecPoints = AliEMCALGetter::GetInstance()->PreShowerRecPoints() ; 



  cout << "AliEMCALClusterizerv1: : event "<<gAlice->GetEvNumber() << endl ;

  cout << "       Found "<< towerRecPoints->GetEntriesFast() << " TOWER Rec Points and " 

	   << preshoRecPoints->GetEntriesFast() << " PRE SHOWER RecPoints" << endl ;



  fRecPointsInRun +=  towerRecPoints->GetEntriesFast() ; 

  fRecPointsInRun +=  preshoRecPoints->GetEntriesFast() ; 



  if(strstr(option,"all")) {



    cout << "Tower clusters " << endl ;

    cout << " Index  Ene(MeV)   Multi  Module     phi     r  theta    Lambda 1   Lambda 2  # of prim  Primaries list "      <<  endl;      

    

    Int_t index ;

    for (index = 0 ; index < towerRecPoints->GetEntries() ; index++) {

      AliEMCALTowerRecPoint * rp = dynamic_cast<AliEMCALTowerRecPoint * >(towerRecPoints->At(index)) ; 

      TVector3  globalpos;  

      rp->GetGlobalPosition(globalpos);

      Float_t lambda[2]; 

      rp->GetElipsAxis(lambda);

      Int_t * primaries; 

      Int_t nprimaries;

      primaries = rp->GetPrimaries(nprimaries);



      cout << setw(4) << rp->GetIndexInList() << "   " 

	   << setw(7) << setprecision(3) << rp->GetEnergy() << "      "

	   << setw(3) << rp->GetMultiplicity() << "      " 

	   << setw(1) << rp->GetEMCALArm() << "     " 

	   << setw(5) << setprecision(4) << globalpos.X() << "  " 

	   << setw(5) << setprecision(4) << globalpos.Y() << "  " 

	   << setw(5) << setprecision(4) << globalpos.Z() << "     "

	   << setw(4) << setprecision(2) << lambda[0]  << "  "

	   << setw(4) << setprecision(2) << lambda[1]  << "  "

	   << setw(2) << nprimaries << "  " ;

     

      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++)

	cout << setw(4) <<   primaries[iprimary] << "  "  ;

      cout << endl ;  	 

    }



    //Now plot Pre shower recPoints



    cout << "-----------------------------------------------------------------------"<<endl ;



    cout << "PreShower clusters " << endl ;

    cout << " Index  Ene(MeV)   Multi  Module     phi     r  theta    Lambda 1   Lambda 2  # of prim  Primaries list "      <<  endl;      

    

    for (index = 0 ; index < preshoRecPoints->GetEntries() ; index++) {

      AliEMCALTowerRecPoint * rp = dynamic_cast<AliEMCALTowerRecPoint *>(preshoRecPoints->At(index)) ; 

      TVector3  globalpos;  

      rp->GetGlobalPosition(globalpos);

      Float_t lambda[2]; 

      rp->GetElipsAxis(lambda);

      Int_t * primaries; 

      Int_t nprimaries;

      primaries = rp->GetPrimaries(nprimaries);



      cout << setw(4) << rp->GetIndexInList() << "   " 

	   << setw(7) << setprecision(3) << rp->GetEnergy() << "      "

	   << setw(3) << rp->GetMultiplicity() << "      " 

	   << setw(1) << rp->GetEMCALArm() << "     " 

	   << setw(5) << setprecision(4) << globalpos.X() << "  " 

	   << setw(5) << setprecision(4) << globalpos.Y() << "  " 

	   << setw(5) << setprecision(4) << globalpos.Z() << "     "

	   << setw(4) << setprecision(2) << lambda[0]  << "  "

	   << setw(4) << setprecision(2) << lambda[1]  << "  "

	   << setw(2) << nprimaries << "  " ;

     

      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++)

	cout << setw(4) <<   primaries[iprimary] << "  "  ;

      cout << endl ;  	 

    }



    cout << "-----------------------------------------------------------------------"<<endl ;

  }

}



