
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
#include "AliEMCALGetter.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCAL.h"
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALClusterizerv1)

Int_t addOn[20][60][60]; 

//____________________________________________________________________________
  AliEMCALClusterizerv1::AliEMCALClusterizerv1() : AliEMCALClusterizer()
{
  // default ctor (to be used mainly by Streamer)
  
  InitParameters() ; 
  fDefaultInit = kTRUE ;
  for(Int_t is=0;is<20;is++){ 
    for(Int_t i=0;i<60;i++){ 
      for(Int_t j=0;j<60;j++){ 
	addOn[is][i][j]=0;
      }
    }
  }
//PH   cout<<"file to read 1"<<endl;
  ReadFile();
//PH   cout<<"file read 1"<<endl;
}

//____________________________________________________________________________
AliEMCALClusterizerv1::AliEMCALClusterizerv1(const TString alirunFileName, const TString eventFolderName)
:AliEMCALClusterizer(alirunFileName, eventFolderName)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
  for(Int_t is=0;is<20;is++){ 
    for(Int_t i=0;i<60;i++){ 
      for(Int_t j=0;j<60;j++){ 
	addOn[is][i][j]=0;
      }
    }
  }
//PH   cout<<"file to read 2"<<endl;
  ReadFile();
//PH   cout<<"file read 2"<<endl;

}

//____________________________________________________________________________
  AliEMCALClusterizerv1::~AliEMCALClusterizerv1()
{
  // dtor
  
}

//____________________________________________________________________________
const TString AliEMCALClusterizerv1::BranchName() const 
{ 
   return GetName();

}

//____________________________________________________________________________
Float_t  AliEMCALClusterizerv1::Calibrate(Int_t amp) const
{
  //To be replased later by the method, reading individual parameters from the database 
  return -fADCpedestalECA + amp * fADCchannelECA ; 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Exec(Option_t * option)
{
  // Steering method to perform clusterization for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1 (by default), then process events until the end.

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALClusterizer"); 
  
  if(strstr(option,"print"))
    Print("") ; 

  AliEMCALGetter * gime = AliEMCALGetter::Instance() ;

  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1;
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ;
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    gime->Event(ievent,"D") ;
    GetCalibrationParameters() ;

    fNumberOfECAClusters = 0;
           
    MakeClusters() ;

    if(fToUnfold)
      MakeUnfolding() ;

    WriteRecPoints() ;

    if(strstr(option,"deb"))  
      PrintRecPoints(option) ;

    //increment the total number of recpoints per run   
    fRecPointsInRun += gime->ECARecPoints()->GetEntriesFast() ;  
  }
  
  Unload();

  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALClusterizer");
    printf("Exec took %f seconds for Clusterizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALClusterizer"), gBenchmark->GetCpuTime("EMCALClusterizer")/nEvents ) ;
  }

  SaveHists();

}

//____________________________________________________________________________
Bool_t AliEMCALClusterizerv1::FindFit(AliEMCALRecPoint * emcRP, AliEMCALDigit ** maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters) const
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 
  // The initial values for fitting procedure are set equal to the positions of local maxima.
  // Cluster will be fitted as a superposition of nPar/3 electromagnetic showers

  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
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

    Int_t relid[2] ;
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
  // Gets the parameters for the calibration from the digitizer
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ;

  if ( !gime->Digitizer() ) 
    gime->LoadDigitizer();
  AliEMCALDigitizer * dig = gime->Digitizer();  

  fADCchannelECA   = dig->GetECAchannel() ;
  fADCpedestalECA  = dig->GetECApedestal();
//PH  cout<<"ChannelECA, peds "<<fADCchannelECA<<" "<<fADCpedestalECA<<endl;
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of EMCAL tasks
  
  AliEMCALGetter * gime = AliEMCALGetter::Instance();
  if(!gime)
    gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName.Data());

  AliEMCALGeometry * geom = gime->EMCALGeometry() ;
//PH   cout<<"gime,geom "<<gime<<" "<<geom<<endl;

//Sub  fNTowers = geom->GetNZ() *  geom->GetNPhi() ;
  fNTowers =400;
  if(!gMinuit) 
    gMinuit = new TMinuit(100) ;
 //Sub if ( !gime->Clusterizer() ) 
 //Sub   gime->PostClusterizer(this); 
 BookHists();
//PH   cout<<"hists booked "<<endl;
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::InitParameters()
{ 
  // Initializes the parameters for the Clusterizer
  fNumberOfECAClusters = 0;

  fECAClusteringThreshold   = 0.5;  // value obtained from Alexei
  fECALocMaxCut = 0.03;

  fECAW0     = 4.5 ;
  fTimeGate = 1.e-8 ; 
  fToUnfold = kFALSE ;
  fRecPointsInRun  = 0 ;
  fMinECut = 0.3;
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

   AliEMCALGeometry * geom = AliEMCALGetter::Instance()->EMCALGeometry() ;

  Int_t rv = 0 ; 

  Int_t relid1[2] ; 
 //Sub geom->AbsToRelNumbering(d1->GetId(), relid1) ; 
    Int_t nSupMod=0;
    Int_t nTower=0;
    Int_t nIphi=0;
    Int_t nIeta=0;
    Int_t iphi=0;
    Int_t ieta=0;
   geom->GetCellIndex(d1->GetId(), nSupMod,nTower,nIphi,nIeta);
   geom->GetCellPhiEtaIndexInSModule(nTower,nIphi,nIeta, iphi,ieta);
   relid1[0]=ieta;
   relid1[1]=iphi;


    Int_t nSupMod1=0;
    Int_t nTower1=0;
    Int_t nIphi1=0;
    Int_t nIeta1=0;
    Int_t iphi1=0;
    Int_t ieta1=0;
   geom->GetCellIndex(d2->GetId(), nSupMod1,nTower1,nIphi1,nIeta1);
   geom->GetCellPhiEtaIndexInSModule(nTower1,nIphi1,nIeta1, iphi1,ieta1);
   Int_t relid2[2] ; 
   relid2[0]=ieta1;
   relid2[1]=iphi1;

  //Sub  geom->AbsToRelNumbering(d2->GetId(), relid2) ; 
  

  Int_t rowdiff = TMath::Abs( relid1[0] - relid2[0] ) ;  
  Int_t coldiff = TMath::Abs( relid1[1] - relid2[1] ) ;  
  
  if (( coldiff <= 1 )  && ( rowdiff <= 1 )){
      rv = 1 ; 
  }
  else {
    if((relid2[0] > relid1[0]) && (relid2[1] > relid1[1]+1)) 
      rv = 2; //  Difference in row numbers is too large to look further 
  }
 
  if (gDebug == 2 ) 
if(rv==1)printf("AreNeighbours: neighbours=%d, id1=%d, relid1=%d,%d \n id2=%d, relid2=%d,%d \n", 
	 rv, d1->GetId(), relid1[0], relid1[1],
	 d2->GetId(), relid2[0], relid2[1]) ;   
  
  return rv ; 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::Unload() 
{
  // Unloads the Digits and RecPoints
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  gime->EmcalLoader()->UnloadDigits() ; 
  gime->EmcalLoader()->UnloadRecPoints() ; 
}
 
//____________________________________________________________________________
void AliEMCALClusterizerv1::WriteRecPoints()
{

  // Creates new branches with given title
  // fills and writes into TreeR.

  AliEMCALGetter *gime = AliEMCALGetter::Instance() ; 

  TObjArray * aECARecPoints = gime->ECARecPoints() ; 

  TClonesArray * digits = gime->Digits() ; 
  TTree * treeR = gime->TreeR(); ; 
  
  Int_t index ;

  //Evaluate position, dispersion and other RecPoint properties for EC section
  for(index = 0; index < aECARecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(index)))->EvalAll(fECAW0,digits) ;
  
  aECARecPoints->Sort() ;

  for(index = 0; index < aECARecPoints->GetEntries(); index++)
    (dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(index)))->SetIndexInList(index) ;

  aECARecPoints->Expand(aECARecPoints->GetEntriesFast()) ; 
  
  Int_t bufferSize = 32000 ;    
  Int_t splitlevel = 0 ; 

  //EC section branch
  TBranch * branchECA = treeR->Branch("EMCALECARP","TObjArray",&aECARecPoints,bufferSize,splitlevel);
  branchECA->SetTitle(BranchName());

  branchECA->Fill() ;

  gime->WriteRecPoints("OVERWRITE");
  gime->WriteClusterizer("OVERWRITE");
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
    
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 

  AliEMCALGeometry * geom = gime->EMCALGeometry() ; 

  TObjArray * aECARecPoints  = gime->ECARecPoints() ; 

  aECARecPoints->Delete() ;

  TClonesArray * digits = gime->Digits() ; 
  TClonesArray * digitsC =  dynamic_cast<TClonesArray*>(digits->Clone()) ;

  // Clusterization starts    
  TIter nextdigit(digitsC) ;
  AliEMCALDigit * digit;

//just for hist
/*
  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigit())) ) { // scan over the list of digitsC
	  digitAmp->Fill(digit->GetAmp());
  }
  */
/////////// 

  while ( (digit = dynamic_cast<AliEMCALDigit *>(nextdigit())) ) { // scan over the list of digitsC
    AliEMCALRecPoint * clu = 0 ; 
    
    TArrayI clusterECAdigitslist(50000);   
   
//Sub    if (gDebug == 2) { 
   //  printf("MakeClusters: id = %d, ene = %f , thre = %f \n", digit->GetId(),Calibrate(digit->GetAmp()),
//         fECAClusteringThreshold) ;  
//    }
////////////////////////temp solution, embedding///////////////////
   int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
   int iphi=0, ieta=0;
          geom->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
          geom->GetCellPhiEtaIndexInSModule(nTower,nIphi,nIeta, iphi,ieta);

         /////////////////////////////////// 

//cout<<ieta<<" "<<iphi<<endl;

//    if ( geom->IsInECA(digit->GetId()) && (Calibrate(digit->GetAmp()) > fECAClusteringThreshold  ) ){
    if (geom->CheckAbsCellId(digit->GetId()) && (Calibrate(digit->GetAmp()+addOn[nSupMod-1][ieta-1][iphi-1]) > fECAClusteringThreshold  ) ){
      //if(addOn[nSupMod-1][ieta-1][iphi-1]>0)cout<<"11 digit, add "<<ieta<<" "<<iphi<<" "<<addOn[nSupMod-1][ieta-1][iphi-1]<<" "<<digit->GetAmp()<<endl;
//	    cout<<"crossed the threshold "<<endl;
      Int_t iDigitInECACluster = 0;
      // start a new Tower RecPoint
      if(fNumberOfECAClusters >= aECARecPoints->GetSize()) 
	  aECARecPoints->Expand(2*fNumberOfECAClusters+1) ;
      AliEMCALRecPoint * rp = new  AliEMCALRecPoint("") ; 
      aECARecPoints->AddAt(rp, fNumberOfECAClusters) ;
      clu = dynamic_cast<AliEMCALRecPoint *>(aECARecPoints->At(fNumberOfECAClusters)) ; 
      fNumberOfECAClusters++ ; 
      clu->AddDigit(*digit, Calibrate(digit->GetAmp()+addOn[nSupMod-1][ieta-1][iphi-1])) ; 
      clusterECAdigitslist[iDigitInECACluster] = digit->GetIndexInList() ;	
      iDigitInECACluster++ ; 
//    cout<<" 1st setp:cluno, digNo "<<fNumberOfECAClusters<<" "<<iDigitInECACluster<<endl;
      digitsC->Remove(digit) ; 
      if (gDebug == 2 ) 
	printf("MakeClusters: OK id = %d, ene = %f , thre = %f \n", digit->GetId(),Calibrate(digit->GetAmp()), fECAClusteringThreshold) ;  
      nextdigit.Reset() ;
      
      AliEMCALDigit * digitN ; 
      Int_t index = 0 ;

      // Find the neighbours
      while (index < iDigitInECACluster){ // scan over digits already in cluster 
	digit =  (AliEMCALDigit*)digits->At(clusterECAdigitslist[index])  ;      
	index++ ; 
        while ( (digitN = (AliEMCALDigit *)nextdigit())) { // scan over the reduced list of digits 
//	cout<<"we have new digit "<<endl;
          // check that the digit is above the min E Cut
         ////////////////////
          geom->GetCellIndex(digitN->GetId(), nSupMod,nTower,nIphi,nIeta);
          geom->GetCellPhiEtaIndexInSModule(nTower,nIphi,nIeta, iphi,ieta);
	  ////////////////
          if(Calibrate(digitN->GetAmp()+addOn[nSupMod-1][ieta-1][iphi-1]) < fMinECut  )  digitsC->Remove(digitN);
//	cout<<" new digit above ECut "<<endl;
	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!! 
//	cout<<" new digit neighbour?? "<<ineb<<endl;
         switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;
	  case 1 :   // are neighbours 
	    //if(addOn[nSupMod-1][ieta-1][iphi-1]>0)cout<<"22 digit, add "<<nSupMod<<" "<<ieta<<" "<<iphi<<" "<<addOn[nSupMod-1][ieta-1][iphi-1]<<" "<<digit->GetAmp()<<endl;
	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp()+addOn[nSupMod-1][ieta-1][iphi-1]) ) ;
	    clusterECAdigitslist[iDigitInECACluster] = digitN->GetIndexInList() ; 
	    iDigitInECACluster++ ; 
//    cout<<"when neighbour: cluno, digNo "<<digit->GetId()<<" "<<fNumberOfECAClusters<<" "<<iDigitInECACluster<<endl;
	    digitsC->Remove(digitN) ;
//	    break ;
//          case 2 :   // too far from each other
//Subh	    goto endofloop1;   
//	    cout<<"earlier go to end of loop"<<endl;   
	  } // switch
    //cout<<"in nextDigit loop "<<fNumberOfECAClusters<<" "<<iDigitInECACluster<<endl;
	} // while digitN
	
//Sub      endofloop1: ;
	nextdigit.Reset() ; 
      } // loop over ECA cluster
    } // energy threshold
    else if(Calibrate(digit->GetAmp()+addOn[nSupMod-1][ieta-1][iphi-1]) < fMinECut  ){
      //if(addOn[nSupMod-1][ieta-1][iphi-1]>0)cout<<"33 digit, add  "<<ieta<<" "<<iphi<<" "<<addOn[nSupMod-1][ieta-1][iphi-1]<<" "<<digit->GetAmp()<<endl;
      digitsC->Remove(digit);
    }
    //cout<<"after endofloop: cluno, digNo "<<fNumberOfECAClusters<<endl;
  } // while digit  
  delete digitsC ;
cout<<"total no of clusters "<<fNumberOfECAClusters<<" from "<<digits->GetEntriesFast()<<" digits"<<endl; 
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::MakeUnfolding() const
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
void  AliEMCALClusterizerv1::UnfoldCluster(AliEMCALRecPoint * /*iniTower*/, 
					   Int_t /*nMax*/, 
					   AliEMCALDigit ** /*maxAt*/, 
					   Float_t * /*maxAtEnergy*/) const
{
  // Performs the unfolding of a cluster with nMax overlapping showers 
  
  Fatal("UnfoldCluster", "--> Unfolding not implemented") ;

}

//_____________________________________________________________________________
void AliEMCALClusterizerv1::UnfoldingChiSquare(Int_t & /*nPar*/, Double_t * /*Grad*/,
					       Double_t & /*fret*/,
					       Double_t * /*x*/, Int_t /*iflag*/)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do
  
  ::Fatal("UnfoldingChiSquare","Unfolding not implemented") ;
}
//____________________________________________________________________________
void AliEMCALClusterizerv1::Print(Option_t * /*option*/)const
{
  // Print clusterizer parameters

  TString message("\n") ; 
  
  if( strcmp(GetName(), "") !=0 ){
    
    // Print parameters
 
    TString taskName(GetName()) ;
    taskName.ReplaceAll(Version(), "") ;
    
    printf("--------------- "); 
    printf(taskName.Data()) ; 
    printf(" "); 
    printf(GetTitle()) ; 
    printf("-----------\n");  
    printf("Clusterizing digits from the file: "); 
    printf(taskName.Data());  
    printf("\n                           Branch: "); 
    printf(GetName()); 
    printf("\n                       ECA Local Maximum cut    = %f", fECALocMaxCut); 
    printf("\n                       ECA Logarithmic weight   = %f", fECAW0); 
    if(fToUnfold)
      printf("\nUnfolding on\n");
    else
      printf("\nUnfolding off\n");
    
    printf("------------------------------------------------------------------"); 
  }
  else
    printf("AliEMCALClusterizerv1 not initialized ") ;
}

//____________________________________________________________________________
void AliEMCALClusterizerv1::PrintRecPoints(Option_t * option)
{
  // Prints list of RecPoints produced at the current pass of AliEMCALClusterizer

  TObjArray * aECARecPoints = AliEMCALGetter::Instance()->ECARecPoints() ; 
  printf("PrintRecPoints: Clusterization result:") ; 
  
  printf("event # %d\n", gAlice->GetEvNumber() ) ;
  printf("           Found %d ECA Rec Points\n ", 
	 aECARecPoints->GetEntriesFast()) ; 

  fRecPointsInRun +=  aECARecPoints->GetEntriesFast() ; 
  
  if(strstr(option,"all")) {
    Int_t index =0;
    printf("\n-----------------------------------------------------------------------\n") ;
    printf("Clusters in ECAL section\n") ;
    printf("Index    Ene(GeV) Multi Module     phi     r   theta    X    Y      Z   Dispersion Lambda 1   Lambda 2  # of prim  Primaries list\n") ;      
   Float_t maxE=0; 
   Float_t maxL1=0; 
   Float_t maxL2=0; 
   Float_t maxDis=0; 

    for (index = 0 ; index < aECARecPoints->GetEntries() ; index++) {
      AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint * >(aECARecPoints->At(index)) ; 
      TVector3  globalpos;  
      //rp->GetGlobalPosition(globalpos);
      TVector3  localpos;  
      rp->GetLocalPosition(localpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      printf("\n%6d  %8.4f  %3d     %4.1f    %4.1f %4.1f  %4.1f %4.1f %4.1f    %4.1f   %4f  %4f    %2d     : ", 
	     rp->GetIndexInList(), rp->GetEnergy(), rp->GetMultiplicity(),
	     globalpos.X(), globalpos.Y(), globalpos.Z(), localpos.X(), localpos.Y(), localpos.Z(), 
	     rp->GetDispersion(), lambda[0], lambda[1], nprimaries) ; 
  /////////////
      if(rp->GetEnergy()>maxE){
	      maxE=rp->GetEnergy();
	      maxL1=lambda[0];
	      maxL2=lambda[1];
	      maxDis=rp->GetDispersion();
      }
      pointE->Fill(rp->GetEnergy());
      pointL1->Fill(lambda[0]);
      pointL2->Fill(lambda[1]);
      pointDis->Fill(rp->GetDispersion());
      pointMult->Fill(rp->GetMultiplicity());
      ///////////// 
      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	printf("%d ", primaries[iprimary] ) ; 
      } 
    }

      MaxE->Fill(maxE);
      MaxL1->Fill(maxL1);
      MaxL2->Fill(maxL2);
      MaxDis->Fill(maxDis);


    printf("\n-----------------------------------------------------------------------\n");
  }
}
  void AliEMCALClusterizerv1::BookHists()
{
	pointE = new TH1F("pointE","point energy", 2000, 0.0, 150.);
	pointL1 = new TH1F("pointL1","point L1", 1000, 0.0, 3.);
	pointL2 = new TH1F("pointL2","point L2", 1000, 0.0, 3.);
	pointDis = new TH1F("pointDis","point Dis", 1000, 0.0, 3.);
	pointMult = new TH1F("pointMult","point Mult", 100, 0.0, 100.);
	digitAmp = new TH1F("digitAmp","Digit Amplitude", 2000, 0.0, 5000.);
	MaxE = new TH1F("maxE","Max point energy", 2000, 0.0, 150.);
	MaxL1 = new TH1F("maxL1","Max point L1", 1000, 0.0, 3.);
	MaxL2 = new TH1F("maxL2","Max point L2", 1000, 0.0, 3.);
	MaxDis = new TH1F("maxDis","Max point Dis", 1000, 0.0, 3.);
}
void AliEMCALClusterizerv1::SaveHists()
{
 recofile=new TFile("reco.root","RECREATE"); 
	pointE->Write();
	pointL1->Write();
	pointL2->Write();
	pointDis->Write();
	pointMult->Write();
	digitAmp->Write();
      MaxE->Write();
      MaxL1->Write();
      MaxL2->Write();
      MaxDis->Write();
recofile->Close();
}

void AliEMCALClusterizerv1::ReadFile()
{
  return; // 3-jan-05
	FILE *fp = fopen("hijing1.dat","r");
        for(Int_t line=0;line<9113;line++){
         Int_t eg,l1,l2,sm;
         Int_t ncols0;
         ncols0 = fscanf(fp,"%d %d %d %d",&sm,&l1,&l2,&eg);
	// cout<<eg<<" "<<l1<<" "<<l2<<endl;
	 addOn[sm-1][l1-1][l2-1]=eg;
	 //addOn[sm-1][l1-1][l2-1]=0;
	}
}

