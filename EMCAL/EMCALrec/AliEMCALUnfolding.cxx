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

//_________________________________________________________________________
//  Base class for the cluster unfolding algorithm 
//*-- Author: Adam Matyja (SUBATECH)
//  Based on unfolding in clusterizerv1 done by Cynthia Hadjidakis
//-- Unfolding for eta~0: Cynthia Hadjidakis - still in AliEMCALCLusterizerv1
//-- Unfolding extension for whole EMCAL: Adam Matyja (SUBATECH & INP PAN)
//
//  unfolds the clusters having several local maxima. 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TClonesArray.h"
#include <TMath.h> 
#include <TMinuit.h>

// --- Standard library ---
#include <cassert>

// --- AliRoot header files ---
#include "AliEMCALUnfolding.h"
#include "AliEMCALGeometry.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliEMCAL.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALReconstructor.h"

#include "AliLog.h"
#include "AliCDBManager.h"
class AliCDBStorage;
#include "AliCDBEntry.h"

Double_t AliEMCALUnfolding::fgSSPars[8]={0.9262,3.365,1.548,0.1625,-0.4195,0.,0.,2.332};
Double_t AliEMCALUnfolding::fgPar5[3]={12.31,-0.007381,-0.06936};
Double_t AliEMCALUnfolding::fgPar6[3]={0.05452,0.0001228,0.001361};

ClassImp(AliEMCALUnfolding)
  
//____________________________________________________________________________
AliEMCALUnfolding::AliEMCALUnfolding():
  fNumberOfECAClusters(0),
  fECALocMaxCut(0),
  fThreshold(0.01),//10 MeV
  fRejectBelowThreshold(0),//split
  fGeom(NULL),
  fRecPoints(NULL),
  fDigitsArr(NULL)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  Init() ;
}

//____________________________________________________________________________
AliEMCALUnfolding::AliEMCALUnfolding(AliEMCALGeometry* geometry):
  fNumberOfECAClusters(0),
  fECALocMaxCut(0),
  fThreshold(0.01),//10 MeV
  fRejectBelowThreshold(0),//split
  fGeom(geometry),
  fRecPoints(NULL),
  fDigitsArr(NULL)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  // use this contructor to avoid usage of Init() which uses runloader
  // change needed by HLT - MP
  if (!fGeom)
  {
    AliFatal("AliEMCALUnfolding: Geometry not initialized.");
  }

}

//____________________________________________________________________________
AliEMCALUnfolding::AliEMCALUnfolding(AliEMCALGeometry* geometry,Float_t ECALocMaxCut,Double_t *SSPars,Double_t *Par5,Double_t *Par6):
  fNumberOfECAClusters(0),
  fECALocMaxCut(ECALocMaxCut),
  fThreshold(0.01),//10 MeV
  fRejectBelowThreshold(0),//split
  fGeom(geometry),
  fRecPoints(NULL),
  fDigitsArr(NULL)
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  // use this contructor to avoid usage of Init() which uses runloader
  // change needed by HLT - MP
  if (!fGeom)
  {
    AliFatal("AliEMCALUnfolding: Geometry not initialized.");
  }
  Int_t i=0;
  for (i = 0; i < 8; i++) fgSSPars[i] = SSPars[i];
  for (i = 0; i < 3; i++) {
    fgPar5[i] = Par5[i];
    fgPar6[i] = Par6[i];
  }

}

//____________________________________________________________________________
void AliEMCALUnfolding::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of EMCAL tasks

  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl && rl->GetAliRun()){
    AliEMCAL* emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if(emcal)fGeom = emcal->GetGeometry();
  }
  
  if(!fGeom)
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  
  AliDebug(1,Form("geom %p",fGeom));
  
  if(!gMinuit) 
    //    gMinuit = new TMinuit(100) ;//the same is in FindFitV2
    gMinuit = new TMinuit(30) ;//the same is in FindFitV2
  
}

//____________________________________________________________________________
  AliEMCALUnfolding::~AliEMCALUnfolding()
{
  // dtor
}

//____________________________________________________________________________
void AliEMCALUnfolding::SetInput(Int_t numberOfECAClusters,TObjArray *recPoints,TClonesArray *digitsArr)
{
  //
  //Set input for unfolding purposes
  //
  SetNumberOfECAClusters(numberOfECAClusters);
  SetRecPoints(recPoints);
  SetDigitsArr(digitsArr);
}

//____________________________________________________________________________
void AliEMCALUnfolding::MakeUnfolding()
{
  // Unfolds clusters using the shape of an ElectroMagnetic shower
  // Performs unfolding of all clusters
  
  AliDebug(4,Form(" V1: total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast()));
  if(fNumberOfECAClusters > 0){
    if (fGeom==0)
      AliFatal("Did not get geometry from EMCALLoader") ;
    //Int_t nModulesToUnfold = fGeom->GetNCells();
    
    Int_t numberOfClustersToUnfold=fNumberOfECAClusters;
    //we unfold only clusters present in the array untill now
    //fNumberOfECAClusters may change due to unfilded clusters
    //so 0 to numberOfClustersToUnfold-1: clusters before unfolding
    //numberOfClustersToUnfold to the end: new clusters from unfolding
    //of course numberOfClustersToUnfold also is decreased but we don't loop over clusters added in UF 
    Int_t index ;
    for(index = 0 ; index < numberOfClustersToUnfold ; index++){
      AliEMCALRecPoint * recPoint = dynamic_cast<AliEMCALRecPoint *>( fRecPoints->At(index) ) ;
      if(recPoint){
        Int_t nMultipl = recPoint->GetMultiplicity() ;
        AliEMCALDigit ** maxAt = new AliEMCALDigit*[nMultipl] ;
        Float_t * maxAtEnergy = new Float_t[nMultipl] ;
        Int_t nMax = recPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fECALocMaxCut,fDigitsArr) ;
        if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0
	  AliDebug(4,Form("  *** V1+UNFOLD *** Cluster index before UF %d",fNumberOfECAClusters));
          if(UnfoldClusterV2(recPoint, nMax, maxAt, maxAtEnergy) ){
	    //if unfolding correct remove old recPoint
            fRecPoints->Remove(recPoint);
            fRecPoints->Compress() ;//is it really needed
            index-- ;
            fNumberOfECAClusters-- ;
	    numberOfClustersToUnfold--;
          }
	  AliDebug(4,Form("  Cluster index after UF %d",fNumberOfECAClusters));
	} else{
          recPoint->SetNExMax(1) ; //Only one local maximum
        }
        
        delete[] maxAt ;
        delete[] maxAtEnergy ;
      } else {
	//AliError("RecPoint NULL"); //end of check if recPoint exist
	Error("MakeUnfolding", "RecPoint NULL, index = %d, fNumberOfECAClusters = %d, numberOfClustersToUnfold = %d",index,fNumberOfECAClusters,numberOfClustersToUnfold) ;
      }
    } // rec point loop
  }//end of check fNumberOfECAClusters
  // End of Unfolding of clusters

  AliDebug(4,Form(" V1+UNFOLD: total no of clusters %d from %d digits",fNumberOfECAClusters,fDigitsArr->GetEntriesFast()));
//  for(Int_t i=0;i<fNumberOfECAClusters;i++){
//    AliEMCALRecPoint * recPoint = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(i));
//    Int_t nMultipl = recPoint->GetMultiplicity() ;
//    Double_t energy=recPoint->GetEnergy();
//    Int_t absIdMaxDigit=recPoint->GetAbsIdMaxDigit();
//    Int_t sm=recPoint->GetSuperModuleNumber();
//    Double_t pointEne=recPoint->GetPointEnergy();
//    Float_t maxEne=recPoint->GetMaximalEnergy();
//    Int_t maxEneInd=recPoint->GetMaximalEnergyIndex();
//    printf("  cluster %d,ncells %d,ene %f,absIdMaxCell %d,sm %d,pointEne %f,maxEne %f,maxEneInd %d\n",i,nMultipl,energy,absIdMaxDigit,sm,pointEne,maxEne,maxEneInd);
//  }

}

//____________________________________________________________________________
Int_t AliEMCALUnfolding::UnfoldOneCluster(AliEMCALRecPoint * iniTower, 
					  Int_t nMax, 
					  AliEMCALDigit ** maxAt, 
					  Float_t * maxAtEnergy,
					  TObjArray *list)
{
  // Input one cluster
  // Output list of clusters
  // returns number of clusters
  // if fit failed or unfolding is not applicable returns 0 and empty list
  
  //**************************** part 1 *******************************************
  // Performs the unfolding of a cluster with nMax overlapping showers 
  
  //cout<<"unfolding check here part 1"<<endl;
  AliDebug(5,Form("  Original cluster E %f, nMax = %d",iniTower->GetEnergy(),nMax ));

  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;
  
  if (fGeom==0)
    AliFatal("Did not get geometry from EMCALLoader") ;
  
  Bool_t rv = FindFitV2(iniTower, maxAt, maxAtEnergy, nPar, fitparameters) ;
  if( !rv ) 
  {
    // Fit failed, return (and remove cluster? - why? I leave the cluster)
    iniTower->SetNExMax(-1) ;
    delete[] fitparameters ;
    return 0;//changed here
  }
  
  //speed up solution for clusters with 2 maxima where one maximum is below threshold fThreshold
  if(nMax==2){
    if(fitparameters[2]<fThreshold || fitparameters[5]<fThreshold){
      AliDebug(1,"One of fitted energy below threshold");
      iniTower->SetNExMax(1) ;
      delete[] fitparameters ;
      return 0;//changed here
    }
  }

  //**************************** part 2 *******************************************
  // create unfolded rec points and fill them with new energy lists
  // First calculate energy deposited in each sell in accordance with
  // fit (without fluctuations): efit[]
  // and later correct this number in acordance with actual energy
  // deposition
  
  //  cout<<"unfolding check here part 2"<<endl;
  Int_t nDigits = iniTower->GetMultiplicity() ;
  Float_t * efit = new Float_t[nDigits] ;//new fitted energy in cells
  Float_t xpar=0.,zpar=0.,epar=0.  ;//center of gravity in cell units
  
  AliEMCALDigit * digit = 0 ;
  Int_t * digitsList = iniTower->GetDigitsList() ;
  
  Int_t iSupMod =  0 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;//x direction
  Int_t ieta    =  0 ;//z direstion
  
  Int_t iparam = 0 ;
  Int_t iDigit = 0 ;
  
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++)
  {
    digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At(digitsList[iDigit] ) ) ;
    if(digit)
    {
      fGeom->GetCellIndex(digit->GetId(),iSupMod,iTower,iIphi,iIeta); 
      fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                         iIphi, iIeta,iphi,ieta);
      EvalParsPhiDependence(digit->GetId(),fGeom);
      
      efit[iDigit] = 0.;
      iparam = 0;
      while(iparam < nPar )
      {
        xpar = fitparameters[iparam] ;
        zpar = fitparameters[iparam+1] ;
        epar = fitparameters[iparam+2] ;

        efit[iDigit] += epar * ShowerShapeV2((Float_t)iphi - xpar,(Float_t)ieta - zpar) ;
        iparam += 3 ;
      }

    } else AliDebug(1,"Digit NULL part 2!");
    
  }//digit loop
  
  //**************************** part 3 *******************************************
  // Now create new RecPoints and fill energy lists with efit corrected to fluctuations
  // so that energy deposited in each cell is distributed between new clusters proportionally
  // to its contribution to efit
  
  Float_t * energiesList = iniTower->GetEnergiesList() ;
  Float_t ratio = 0. ;
  Float_t eDigit = 0. ;
  Int_t nSplittedClusters=(Int_t)nPar/3;
  
  Float_t * correctedEnergyList = new Float_t[nDigits*nSplittedClusters];
  //above - temporary table with energies after unfolding.
  //the order is following: 
  //first cluster <first cell - last cell>, 
  //second cluster <first cell - last cell>, etc.
  
  //**************************** sub-part 3.1 *************************************
  //If not the energy from a given cell in the cluster is divided in correct proportions 
  //in accordance to the other clusters and added to them and set to 0.
  
  //  cout<<"unfolding check here part 3.1"<<endl;

  iparam = 0 ;
  while(iparam < nPar )
  {
    xpar = fitparameters[iparam] ;
    zpar = fitparameters[iparam+1] ;
    epar = fitparameters[iparam+2] ;
    
    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++)
    {
      digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At( digitsList[iDigit] ) ) ;
      if(digit)
      {
        fGeom->GetCellIndex(digit->GetId(),iSupMod,iTower,iIphi,iIeta); 
        fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                           iIphi, iIeta,iphi,ieta);
        
        EvalParsPhiDependence(digit->GetId(),fGeom);
       
        if(efit[iDigit]==0) 
        {//just for sure
          correctedEnergyList[iparam/3*nDigits+iDigit] = 0.;//correction here
          continue;
        }
        
        ratio = epar * ShowerShapeV2((Float_t)iphi - xpar,(Float_t)ieta - zpar) / efit[iDigit] ;
        eDigit = energiesList[iDigit] * ratio ;
        
        //add energy to temporary matrix
        correctedEnergyList[iparam/3*nDigits+iDigit] = eDigit;
        
      } else AliDebug(1,"NULL digit part 3");
    }//digit loop 
    iparam += 3 ;
  }//while
  
  //**************************** sub-part 3.2 *************************************
  //here we check if energy of the cell in the cluster after unfolding is above threshold. 
  //here we correct energy for each cell and cluster
  //  cout<<"unfolding check here part 3.2"<<endl;


  //here we have 3 possibilities
  //when after UF cell energy in cluster is below threshold:
  //1 - keep it associated to cluster - equivalent of threshold=0
  //2 - default - split (or add) energy of that cell into that cell in the other cluster(s)
  //3 - reject that cell from cluster - fraction energy in cell=0 - breaks energy conservation
  //Bool_t rejectBelowThreshold=kTRUE;//default option = 2 - split = kFALSE

  if(fThreshold > 0){//option 2 or 3
    if(fRejectBelowThreshold){//option 3
      for(iDigit = 0 ; iDigit < nDigits ; iDigit++){//digit loop
	for(iparam = 0 ; iparam < nPar ; iparam+=3){//param0 loop = energy loop
	  if(correctedEnergyList[iparam/3*nDigits+iDigit] < fThreshold ) correctedEnergyList[iparam/3*nDigits+iDigit]=0.;
	}
      }
    }else{//option 2
      Float_t maximumEne=0.;
      Int_t maximumIndex=0;
      Bool_t isAnyBelowThreshold=kFALSE;
      //  Float_t Threshold=0.01;
      Float_t * energyFraction = new Float_t[nSplittedClusters];
      Int_t iparam2 = 0 ;
      for(iDigit = 0 ; iDigit < nDigits ; iDigit++){
	isAnyBelowThreshold=kFALSE;
	maximumEne=0.;
	for(iparam = 0 ; iparam < nPar ; iparam+=3){
	  if(correctedEnergyList[iparam/3*nDigits+iDigit] < fThreshold ) isAnyBelowThreshold = kTRUE;
	  if(correctedEnergyList[iparam/3*nDigits+iDigit] > maximumEne) 
	    {
	      maximumEne = correctedEnergyList[iparam/3*nDigits+iDigit];
	      maximumIndex = iparam;
	    }
	}//end of loop over clusters after unfolding
	
	if(!isAnyBelowThreshold) continue; //no cluster-cell below threshold 
	
	if(maximumEne < fThreshold) 
	  {//add all cluster cells and put energy into max index, other set to 0
	    maximumEne=0.;
	    for(iparam = 0 ; iparam < nPar ; iparam+=3)
	      {
		maximumEne+=correctedEnergyList[iparam/3*nDigits+iDigit];
		correctedEnergyList[iparam/3*nDigits+iDigit]=0;
	      }
	    correctedEnergyList[maximumIndex/3*nDigits+iDigit]=maximumEne;
	    continue;
	  }//end if
	
	//divide energy of cell below threshold in the correct proportion and add to other cells
	maximumEne=0.;//not used any more so use it for the energy sum 
	for(iparam = 0 ; iparam < nPar ; iparam+=3)
	  {//calculate energy sum
	    if(correctedEnergyList[iparam/3*nDigits+iDigit] < fThreshold) energyFraction[iparam/3]=0;
	    else 
	      {
		energyFraction[iparam/3]=1.;
		maximumEne+=correctedEnergyList[iparam/3*nDigits+iDigit];
	      }
	  }//end of loop over clusters after unfolding
	if(maximumEne>0.) {
	  for(iparam = 0 ; iparam < nPar ; iparam+=3){//calculate fraction
	    energyFraction[iparam/3] = energyFraction[iparam/3] * correctedEnergyList[iparam/3*nDigits+iDigit] / maximumEne;
	  }
	  
	  for(iparam = 0 ; iparam < nPar ; iparam+=3)
	    {//add energy from cells below threshold to others
	      if(energyFraction[iparam/3]>0.) continue;
	      else
		{
		  for(iparam2 = 0 ; iparam2 < nPar ; iparam2+=3)
		    {
		      correctedEnergyList[iparam2/3*nDigits+iDigit] += (energyFraction[iparam2/3] * 
									correctedEnergyList[iparam/3*nDigits+iDigit]) ;
		    }//inner loop
		  correctedEnergyList[iparam/3*nDigits+iDigit] = 0.;
		}
	    }
	} else {
	  //digit energy to be set to 0
	  for(iparam = 0 ; iparam < nPar ; iparam+=3)
	    {
	      correctedEnergyList[iparam/3*nDigits+iDigit] = 0.;
	    }
	}//correction for: is energy>0
	
      }//end of loop over digits
      delete[] energyFraction;
      
    }//end of option 2 or 3 
  } else {//option 1
    //do nothing
  }

  
  //**************************** sub-part 3.3 *************************************
  //here we add digits to recpoints with corrected energy
  //  cout<<"unfolding check here part 3.3"<<endl;

  Int_t newClusterIndex=0;
  iparam = 0 ;
  while(iparam < nPar )
  {
    AliEMCALRecPoint * recPoint = 0 ;
    
    if(nSplittedClusters >= list->GetSize())
      list->Expand(nSplittedClusters);
    
    //add recpoint
    (*list)[newClusterIndex] = new AliEMCALRecPoint("") ;
    recPoint = dynamic_cast<AliEMCALRecPoint *>( list->At(newClusterIndex) ) ;
    
    if(recPoint){//recPoint present -> good
      recPoint->SetNExMax(nSplittedClusters) ;//can be wrong number, to be corrected in outer method
      
      for(iDigit = 0 ; iDigit < nDigits ; iDigit ++) {
        digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At( digitsList[iDigit] ) ) ;
        if(digit && correctedEnergyList[iparam/3*nDigits+iDigit]>0. ){
          //if(correctedEnergyList[iparam/3*nDigits+iDigit]<fThreshold) printf("Final E cell %f < %f\n",correctedEnergyList[iparam/3*nDigits+iDigit],fThreshold);
          recPoint->AddDigit( *digit, correctedEnergyList[iparam/3*nDigits+iDigit], kFALSE ) ; //FIXME, need to study the shared case
        } else {
          AliDebug(1,Form("NULL digit part3.3 or NULL energy=%f",correctedEnergyList[iparam/3*nDigits+iDigit]));
        }
      }//digit loop

      if(recPoint->GetMultiplicity()==0){//recpoint exists but no digits associated -> remove from list
	delete (*list)[newClusterIndex];
	list->RemoveAt(newClusterIndex);
	nSplittedClusters--;
	newClusterIndex--;//decrease cluster number
      }else {//recPoint exists and has digits associated -> very good increase number of clusters 
	AliDebug(5,Form("cluster %d, digit no %d, energy %f",iparam/3,(recPoint->GetDigitsList())[0],(recPoint->GetEnergiesList())[0]));
      }
      
    } else {//recPoint empty -> remove from list
      AliError("NULL RecPoint");
      //protection from recpoint with no digits
      delete (*list)[newClusterIndex];
      list->RemoveAt(newClusterIndex);
      nSplittedClusters--;
      newClusterIndex--;//decrease cluster number
    }

    iparam += 3 ;
    newClusterIndex++;
  }//while
  
  delete[] fitparameters ;
  delete[] efit ;
  delete[] correctedEnergyList ;

//  print 
  AliDebug(5,Form("  nSplittedClusters %d, fNumberOfECAClusters %d, newClusterIndex %d,list->Entries() %d\n",nSplittedClusters,fNumberOfECAClusters,newClusterIndex,list->GetEntriesFast() ));
  
  //  cout<<"end of unfolding check part 3.3"<<endl;
  return nSplittedClusters;
}

//____________________________________________________________________________
Bool_t AliEMCALUnfolding::UnfoldClusterV2(AliEMCALRecPoint * iniTower, 
					  Int_t nMax, 
					  AliEMCALDigit ** maxAt, 
					  Float_t * maxAtEnergy)
{
  // Extended to whole EMCAL 
  // Performs the unfolding of a cluster with nMax overlapping showers 
  // Returns true if success (1->several clusters), otherwise false (fit failed)

  TObjArray *list =new TObjArray(2);//temporary object
  Int_t nUnfoldedClusters=UnfoldOneCluster(iniTower,nMax,maxAt,maxAtEnergy,list);

  // here we write new clusters from list to fRecPoints 
  AliDebug(5,Form("Number of clusters after unfolding %d",list->GetEntriesFast()));
  Int_t iDigit=0;
  AliEMCALDigit * digit = 0 ;
  for(Int_t i=0;i<list->GetEntriesFast();i++) {
    AliEMCALRecPoint * recPoint = 0 ;
    
    if(fNumberOfECAClusters >= fRecPoints->GetSize())
      fRecPoints->Expand(2*fNumberOfECAClusters) ;
    
    //add recpoint
    (*fRecPoints)[fNumberOfECAClusters] = new AliEMCALRecPoint("") ;//fNumberOfECAClusters-1 is old cluster before unfolding
    recPoint = dynamic_cast<AliEMCALRecPoint *>( fRecPoints->At(fNumberOfECAClusters) ) ;
		AliEMCALRecPoint * rpUFOne = dynamic_cast<AliEMCALRecPoint *>(list->At(i)) ;
		
    if( recPoint && rpUFOne ){//recPoint present -> good
      
			recPoint->SetNExMax(list->GetEntriesFast()) ;
			
      Int_t   *digitsList = rpUFOne->GetDigitsList();
      Float_t *energyList = rpUFOne->GetEnergiesList();

      if(!digitsList || ! energyList)
      {
        AliDebug(-1,"No digits index or energy available");
				delete (*fRecPoints)[fNumberOfECAClusters];
				fRecPoints->RemoveAt(fNumberOfECAClusters);
        continue;
      }
      
      AliDebug(5,Form("cluster %d, digit no %d, energy %f\n",i,digitsList[0],energyList[0]));

      for(iDigit = 0 ; iDigit < rpUFOne->GetMultiplicity(); iDigit ++) {
        digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At( digitsList[iDigit] ) ) ;
        if(digit) recPoint->AddDigit( *digit, energyList[iDigit], kFALSE ) ; //FIXME, need to study the shared case
      }//digit loop
      fNumberOfECAClusters++ ; 
    } else {//recPoint empty -> remove from list
      AliError("NULL RecPoint");
      delete (*fRecPoints)[fNumberOfECAClusters];
      fRecPoints->RemoveAt(fNumberOfECAClusters);
    }

  }//loop over unfolded clusters
  
  //print energy of new unfolded clusters
  AliDebug(5,Form("  nUnfoldedClusters %d, fNumberOfECAClusters %d",nUnfoldedClusters,fNumberOfECAClusters ));
  for(Int_t inewclus=0; inewclus<nUnfoldedClusters;inewclus++){
		AliEMCALRecPoint * rp = dynamic_cast<AliEMCALRecPoint *>(fRecPoints->At(fNumberOfECAClusters-1-inewclus));
    if(rp) AliDebug(5,Form("  Unfolded cluster %d E %f",inewclus, rp->GetEnergy() ));
  }

  //clear tables  
  list->SetOwner(kTRUE);
  list->Delete();
  delete list;
  if(nUnfoldedClusters>1) return kTRUE;
  return kFALSE;
}



//____________________________________________________________________________
Bool_t AliEMCALUnfolding::UnfoldClusterV2old(AliEMCALRecPoint * iniTower, 
					  Int_t nMax, 
					  AliEMCALDigit ** maxAt, 
					  Float_t * maxAtEnergy)
{
  // Extended to whole EMCAL 
  // Performs the unfolding of a cluster with nMax overlapping showers 
  
  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;
  
  if (fGeom==0)
    AliFatal("Did not get geometry from EMCALLoader") ;
  
  Bool_t rv = FindFitV2(iniTower, maxAt, maxAtEnergy, nPar, fitparameters) ;
  if( !rv ) {
    // Fit failed, return (and remove cluster? - why? I leave the cluster)
    iniTower->SetNExMax(-1) ;
    delete[] fitparameters ;
    return kFALSE;
  }
  
  // create unfolded rec points and fill them with new energy lists
  // First calculate energy deposited in each sell in accordance with
  // fit (without fluctuations): efit[]
  // and later correct this number in acordance with actual energy
  // deposition
  
  Int_t nDigits = iniTower->GetMultiplicity() ;
  Float_t * efit = new Float_t[nDigits] ;//new fitted energy in cells
  Float_t xpar=0.,zpar=0.,epar=0.  ;//center of gravity in cell units
  
  AliEMCALDigit * digit = 0 ;
  Int_t * digitsList = iniTower->GetDigitsList() ;
  
  Int_t iSupMod =  0 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;//x direction
  Int_t ieta    =  0 ;//z direstion
  
  Int_t iparam = 0 ;
  Int_t iDigit = 0 ;
  
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
    digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At(digitsList[iDigit] ) ) ;
    if(digit){
      fGeom->GetCellIndex(digit->GetId(),iSupMod,iTower,iIphi,iIeta); 
      fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                         iIphi, iIeta,iphi,ieta);
      EvalParsPhiDependence(digit->GetId(),fGeom);
      
      efit[iDigit] = 0.;
      iparam = 0;
      while(iparam < nPar ){
        xpar = fitparameters[iparam] ;
        zpar = fitparameters[iparam+1] ;
        epar = fitparameters[iparam+2] ;
        iparam += 3 ;
        
        efit[iDigit] += epar * ShowerShapeV2((Float_t)iphi - xpar,(Float_t)ieta - zpar) ;
      }
    } else AliError("Digit NULL!");
    
  }//digit loop
  
  // Now create new RecPoints and fill energy lists with efit corrected to fluctuations
  // so that energy deposited in each cell is distributed between new clusters proportionally
  // to its contribution to efit
  
  Float_t * energiesList = iniTower->GetEnergiesList() ;
  Float_t ratio = 0 ;
  
  iparam = 0 ;
  while(iparam < nPar ){
    xpar = fitparameters[iparam] ;
    zpar = fitparameters[iparam+1] ;
    epar = fitparameters[iparam+2] ;
    iparam += 3 ;
    
    AliEMCALRecPoint * recPoint = 0 ;
    
    if(fNumberOfECAClusters >= fRecPoints->GetSize())
      fRecPoints->Expand(2*fNumberOfECAClusters) ;
    
    //add recpoint
    (*fRecPoints)[fNumberOfECAClusters] = new AliEMCALRecPoint("") ;
    recPoint = dynamic_cast<AliEMCALRecPoint *>( fRecPoints->At(fNumberOfECAClusters) ) ;
    
    if(recPoint){
      
      fNumberOfECAClusters++ ;
      recPoint->SetNExMax((Int_t)nPar/3) ;
      
      Float_t eDigit = 0. ;
      for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
        digit = dynamic_cast<AliEMCALDigit*>( fDigitsArr->At( digitsList[iDigit] ) ) ;
        if(digit){
          fGeom->GetCellIndex(digit->GetId(),iSupMod,iTower,iIphi,iIeta); 
          fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                             iIphi, iIeta,iphi,ieta);
          EvalParsPhiDependence(digit->GetId(),fGeom);
          if(efit[iDigit]==0) continue;//just for sure
          ratio = epar * ShowerShapeV2((Float_t)iphi - xpar,(Float_t)ieta - zpar) / efit[iDigit] ;
          eDigit = energiesList[iDigit] * ratio ;
          recPoint->AddDigit( *digit, eDigit, kFALSE ) ; //FIXME, need to study the shared case
        } else AliError("NULL digit");
      }//digit loop 
    } else AliError("NULL RecPoint");
  }//while
  
  delete[] fitparameters ;
  delete[] efit ;
  
  return kTRUE;
}


//____________________________________________________________________________
Bool_t AliEMCALUnfolding::FindFitV2(AliEMCALRecPoint * recPoint, AliEMCALDigit ** maxAt, 
					const Float_t* maxAtEnergy,
					Int_t nPar, Float_t * fitparameters) const
{
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima
  // The initial values for fitting procedure are set equal to the
  // positions of local maxima.       
  // Cluster will be fitted as a superposition of nPar/3
  // electromagnetic showers

  if (fGeom==0) AliFatal("Did not get geometry from EMCALLoader");
	
  if(!gMinuit){
    //    gMinuit = new TMinuit(100) ;//max 100 parameters
    if(nPar<30) gMinuit = new TMinuit(30);
    else gMinuit = new TMinuit(nPar) ;//max nPar parameters
    //
  } else {
    if(gMinuit->fMaxpar < nPar) {
      delete gMinuit;
      gMinuit = new TMinuit(nPar);
    }
  }

  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliEMCALUnfolding::UnfoldingChiSquareV2) ;
  // To set the address of the minimization function
  TList * toMinuit = new TList();
  toMinuit->AddAt(recPoint,0) ;
  toMinuit->AddAt(fDigitsArr,1) ;
  toMinuit->AddAt(fGeom,2) ;

  gMinuit->SetObjectFit(toMinuit) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliEMCALDigit * digit ;

  Int_t ierflg  = 0;
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;

  Int_t iSupMod =  0 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;//x direction
  Int_t ieta    =  0 ;//z direstion

  for(iDigit = 0; iDigit < nDigits; iDigit++)
  {
    digit = maxAt[iDigit];
    if(!digit)
    {
      AliError("Null digit pointer");
      continue;
    }
    
    fGeom->GetCellIndex(digit->GetId(),iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
				       iIphi, iIeta,iphi,ieta);

    Float_t energy = maxAtEnergy[iDigit] ;

    //gMinuit->mnparm(index, "x",  iphi, 0.1, 0, 0, ierflg) ;//original
    gMinuit->mnparm(index, "x",  iphi, 0.05, 0, 0, ierflg) ;
    index++ ;
    if(ierflg != 0){
      Error("FindFit", "EMCAL Unfolding unable to set initial value for fit procedure: x=%d, param.id=%d, nMaxima=%d",iphi,index-1,nPar/3 ) ;
      toMinuit->Clear();
      delete toMinuit ;
      return kFALSE;
    }
    //gMinuit->mnparm(index, "z",  ieta, 0.1, 0, 0, ierflg) ;//original
    gMinuit->mnparm(index, "z",  ieta, 0.05, 0, 0, ierflg) ;
    index++ ;
    if(ierflg != 0){
      Error("FindFit", "EMCAL Unfolding unable to set initial value for fit procedure: z=%d, param.id=%d, nMaxima=%d", ieta, index-1,nPar/3) ;
      toMinuit->Clear();
      delete toMinuit ;
      return kFALSE;
    }
    //gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;//original
    gMinuit->mnparm(index, "Energy",  energy , 0.001*energy, 0., 5.*energy, ierflg) ;//was 0.05
    index++ ;
    if(ierflg != 0){
      Error("FindFit", "EMCAL Unfolding unable to set initial value for fit procedure: energy = %f, param.id=%d, nMaxima=%d", energy, index-1, nPar/3) ;
      toMinuit->Clear();
      delete toMinuit ;
      return kFALSE;
    }
  }

  Double_t p0 = 0.1 ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; 
                      // The number of function call slightly depends on it.
  //  Double_t p1 = 1.0 ;// par to gradient 
  Double_t p2 = 0.0 ;
  //  Double_t p3 = 3.0 ;
  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls
  //  gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient
  gMinuit->SetMaxIterations(5);//was 5
  gMinuit->mnexcm("SET NOW", &p2 , 0, ierflg) ;  // No Warnings
  //gMinuit->mnexcm("SET PRI", &p3 , 3, ierflg) ;  // printouts

  gMinuit->mnexcm("MIGRAD", &p0, 0, ierflg) ;    // minimize
  //gMinuit->mnexcm("MINI", &p0, 0, ierflg) ;    // minimize
  if(ierflg == 4){  // Minimum not found
    AliDebug(1,"EMCAL Unfolding  Fit not converged, cluster abandoned " ) ;
    toMinuit->Clear();
    delete toMinuit ;
    return kFALSE ;
  }
  for(index = 0; index < nPar; index++){
    Double_t err = 0. ;
    Double_t val = 0. ;
    gMinuit->GetParameter(index, val, err) ;    // Returns value and error ofOA parameter index
    fitparameters[index] = val ;
  }

  toMinuit->Clear();
  delete toMinuit ;

  if(gMinuit->fMaxpar>30) delete gMinuit;

  return kTRUE;

}

//____________________________________________________________________________
Double_t  AliEMCALUnfolding::ShowerShapeV2(Double_t x, Double_t y)
{ 
  // extended to whole EMCAL 
  // Shape of the shower
  // If you change this function, change also the gradient evaluation in ChiSquare()

  Double_t r = fgSSPars[7]*TMath::Sqrt(x*x+y*y);
  Double_t rp1  = TMath::Power(r, fgSSPars[1]) ;
  Double_t rp5  = TMath::Power(r, fgSSPars[5]) ;
  Double_t shape = fgSSPars[0]*TMath::Exp( -rp1 * (1. / (fgSSPars[2] + fgSSPars[3] * rp1) + fgSSPars[4] / (1 + fgSSPars[6] * rp5) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void AliEMCALUnfolding::UnfoldingChiSquareV2(Int_t & nPar, Double_t * Grad,
					     Double_t & fret,
					     Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do
  
  TList * toMinuit = dynamic_cast<TList*>( gMinuit->GetObjectFit() ) ;
  if(toMinuit){
    AliEMCALRecPoint * recPoint = dynamic_cast<AliEMCALRecPoint*>( toMinuit->At(0) )  ;
    TClonesArray * digits = dynamic_cast<TClonesArray*>( toMinuit->At(1) )  ;
    // A bit buggy way to get an access to the geometry
    // To be revised!
    AliEMCALGeometry *geom = dynamic_cast<AliEMCALGeometry *>(toMinuit->At(2));
    
    if(recPoint && digits && geom){
      
      Int_t * digitsList     = recPoint->GetDigitsList() ;
      
      Int_t nOdigits = recPoint->GetDigitsMultiplicity() ;
      
      Float_t * energiesList = recPoint->GetEnergiesList() ;
      
      fret = 0. ;
      Int_t iparam = 0 ;
      
      if(iflag == 2)
        for(iparam = 0 ; iparam < nPar ; iparam++)
          Grad[iparam] = 0 ; // Will evaluate gradient
      
      Double_t efit = 0. ;
      
      AliEMCALDigit * digit ;
      Int_t iDigit ;
      
      Int_t iSupMod =  0 ;
      Int_t iTower  =  0 ;
      Int_t iIphi   =  0 ;
      Int_t iIeta   =  0 ;
      Int_t iphi    =  0 ;//x direction
      Int_t ieta    =  0 ;//z direstion
      
      
      for( iDigit = 0 ; iDigit < nOdigits ; iDigit++) {
        if(energiesList[iDigit]==0) continue;
        
        digit = dynamic_cast<AliEMCALDigit*>( digits->At( digitsList[iDigit] ) );
        
        if(digit){
	  geom->GetCellIndex(digit->GetId(),iSupMod,iTower,iIphi,iIeta); 
	  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					    iIphi, iIeta,iphi,ieta);
	  EvalParsPhiDependence(digit->GetId(),geom);
	  
	  if(iflag == 2){  // calculate gradient
	    Int_t iParam = 0 ;
	    efit = 0. ;
	    while(iParam < nPar ){
	      Double_t dx = ((Float_t)iphi - x[iParam]) ;
	      iParam++ ;
	      Double_t dz = ((Float_t)ieta - x[iParam]) ;
	      iParam++ ;
	      efit += x[iParam] * ShowerShapeV2(dx,dz) ;
	      iParam++ ;
	    }
	    
	    Double_t sum = 2. * (efit - energiesList[iDigit]) / energiesList[iDigit] ; // Here we assume, that sigma = sqrt(E)
	    iParam = 0 ;
	    while(iParam < nPar ){
	      Double_t xpar = x[iParam] ;
	      Double_t zpar = x[iParam+1] ;
	      Double_t epar = x[iParam+2] ;
	      
	      Double_t dr = fgSSPars[7]*TMath::Sqrt( ((Float_t)iphi - xpar) * ((Float_t)iphi - xpar) + ((Float_t)ieta - zpar) * ((Float_t)ieta - zpar) );
	      Double_t shape = sum * ShowerShapeV2((Float_t)iphi - xpar,(Float_t)ieta - zpar) ;
	      Double_t rp1  = TMath::Power(dr, fgSSPars[1]) ;
	      Double_t rp5  = TMath::Power(dr, fgSSPars[5]) ;
	      
	      Double_t deriv = -2 * TMath::Power(dr,fgSSPars[1]-2.) * fgSSPars[7] * fgSSPars[7] * 
		(fgSSPars[1] * ( 1/(fgSSPars[2]+fgSSPars[3]*rp1) + fgSSPars[4]/(1+fgSSPars[6]*rp5) ) - 
		 (fgSSPars[1]*fgSSPars[3]*rp1/( (fgSSPars[2]+fgSSPars[3]*rp1)*(fgSSPars[2]+fgSSPars[3]*rp1) ) + 
		  fgSSPars[4]*fgSSPars[5]*fgSSPars[6]*rp5/( (1+fgSSPars[6]*rp5)*(1+fgSSPars[6]*rp5) ) ) );
	      
	      //Double_t deriv =-1.33 * TMath::Power(dr,0.33)*dr * ( 1.57 / ( (1.57 + 0.0860 * r133) * (1.57 + 0.0860 * r133) )
	      //                                                   - 0.55 / (1 + 0.000563 * r669) / ( (1 + 0.000563 * r669) * (1 + 0.000563 * r669) ) ) ;
	      
	      Grad[iParam] += epar * shape * deriv * ((Float_t)iphi - xpar) ;  // Derivative over x
	      iParam++ ;
	      Grad[iParam] += epar * shape * deriv * ((Float_t)ieta - zpar) ;  // Derivative over z
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
	    efit += epar * ShowerShapeV2((Float_t)iphi - xpar,(Float_t)ieta - zpar) ;
	  }
	  
	  fret += (efit-energiesList[iDigit])*(efit-energiesList[iDigit])/energiesList[iDigit] ;
	  // Here we assume, that sigma = sqrt(E) 
        } else printf("AliEMCALUnfoding::UnfoldingChiSquareV2 - NULL digit!, nPar %d \n", nPar); // put nPar here to cheat coverity and rule checker
      } // digit loop
    } // recpoint, digits and geom not NULL
  }// List is not NULL
  
}


//____________________________________________________________________________
void AliEMCALUnfolding::SetShowerShapeParams(Double_t *pars){
  for(UInt_t i=0;i<7;++i)
    fgSSPars[i]=pars[i];
  if(pars[2]==0. && pars[3]==0.) fgSSPars[2]=1.;//to avoid dividing by 0
}

//____________________________________________________________________________
void AliEMCALUnfolding::SetPar5(Double_t *pars){
  for(UInt_t i=0;i<3;++i)
    fgPar5[i]=pars[i];
}

//____________________________________________________________________________
void AliEMCALUnfolding::SetPar6(Double_t *pars){
  for(UInt_t i=0;i<3;++i)
    fgPar6[i]=pars[i];
}

//____________________________________________________________________________
void AliEMCALUnfolding::EvalPar5(Double_t phi){
  //
  //Evaluate the 5th parameter of the shower shape function
  //phi in degrees range (-10,10)
  //
  //fSSPars[5] = 12.31 - phi*0.007381 - phi*phi*0.06936;
  fgSSPars[5] = fgPar5[0] + phi * fgPar5[1] + phi*phi * fgPar5[2];
}

//____________________________________________________________________________
void AliEMCALUnfolding::EvalPar6(Double_t phi){
  //
  //Evaluate the 6th parameter of the shower shape function
  //phi in degrees range (-10,10)
  //
  //fSSPars[6] = 0.05452 + phi*0.0001228 + phi*phi*0.001361;
  fgSSPars[6] = fgPar6[0] + phi * fgPar6[1] + phi*phi * fgPar6[2];
}

//____________________________________________________________________________
void AliEMCALUnfolding::EvalParsPhiDependence(Int_t absId, const AliEMCALGeometry *geom){
  //
  // calculate params p5 and p6 depending on the phi angle in global coordinate
  // for the cell with given absId index
  // output phiglob should be in range -10 to 10 degree
  //
  Double_t etaGlob = 0.;//eta in global c.s. - unused
  Double_t phiGlob = 0.;//phi in global c.s. in radians
  geom->EtaPhiFromIndex(absId, etaGlob, phiGlob);
  if(phiGlob<0)phiGlob+=TMath::TwoPi();
  phiGlob*=180./TMath::Pi();
  phiGlob-=20.*(Int_t)(phiGlob/20.);
  Int_t superModule=geom->GetSuperModuleNumber(absId);
  if(superModule==10 || superModule==11 || superModule==18 || superModule==19) //sm 10,11,18,19 shift by 3.5 deg other 10 deg
    phiGlob-=3.5;
  else
    phiGlob-=10.;

  EvalPar5(phiGlob);
  EvalPar6(phiGlob);
}

