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

/* History of cvs commits:
 *
 * $Log: AliPHOSClusterizerv1.cxx,v $
 * Revision 1.118  2007/12/11 21:23:26  kharlov
 * Added possibility to swith off unfolding
 *
 * Revision 1.117  2007/10/18 08:42:05  kharlov
 * Bad channels cleaned before clusterization
 *
 * Revision 1.116  2007/10/01 20:24:08  kharlov
 * Memory leaks fixed
 *
 * Revision 1.115  2007/09/26 14:22:17  cvetan
 * Important changes to the reconstructor classes. Complete elimination of the run-loaders, which are now steered only from AliReconstruction. Removal of the corresponding Reconstruct() and FillESD() methods.
 *
 * Revision 1.114  2007/09/06 16:06:44  kharlov
 * Absence of sorting results in loose of all unfolded clusters
 *
 * Revision 1.113  2007/08/28 12:55:07  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.112  2007/08/22 09:20:50  hristov
 * Updated QA classes (Yves)
 *
 * Revision 1.111  2007/08/08 12:11:28  kharlov
 * Protection against uninitialized fQADM
 *
 * Revision 1.110  2007/08/07 14:16:00  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.109  2007/07/24 17:20:35  policheh
 * Usage of RecoParam objects instead of hardcoded parameters in reconstruction.
 * (See $ALICE_ROOT/PHOS/macros/BeamTest2006/RawReconstruction.C).
 *
 * Revision 1.108  2007/06/18 07:00:51  kharlov
 * Bug fix for attempt to use AliPHOSEmcRecPoint after its deletion
 *
 * Revision 1.107  2007/05/25 14:12:26  policheh
 * Local to tracking CS transformation added for CPV rec. points
 *
 * Revision 1.106  2007/05/24 13:01:22  policheh
 * Local to tracking CS transformation invoked for each EMC rec.point
 *
 * Revision 1.105  2007/05/02 13:41:22  kharlov
 * Mode protection against absence of calib.data from AliPHOSCalibData to AliPHOSClusterizerv1::GetCalibrationParameters()
 *
 * Revision 1.104  2007/04/27 16:55:53  kharlov
 * Calibration stops if PHOS CDB objects do not exist
 *
 * Revision 1.103  2007/04/11 11:55:45  policheh
 * SetDistancesToBadChannels() added.
 *
 * Revision 1.102  2007/03/28 19:18:15  kharlov
 * RecPoints recalculation in TSM removed
 *
 * Revision 1.101  2007/03/06 06:51:27  kharlov
 * Calculation of cluster properties dep. on vertex posponed to TrackSegmentMaker
 *
 * Revision 1.100  2007/01/10 11:07:26  kharlov
 * Raw digits writing to file (B.Polichtchouk)
 *
 * Revision 1.99  2006/11/07 16:49:51  kharlov
 * Corrections for next event switch in case of raw data (B.Polichtchouk)
 *
 * Revision 1.98  2006/10/27 17:14:27  kharlov
 * Introduce AliDebug and AliLog (B.Polichtchouk)
 *
 * Revision 1.97  2006/08/29 11:41:19  kharlov
 * Missing implementation of ctors and = operator are added
 *
 * Revision 1.96  2006/08/25 16:56:30  kharlov
 * Compliance with Effective C++
 *
 * Revision 1.95  2006/08/11 12:36:26  cvetan
 * Update of the PHOS code needed in order to read and reconstruct the beam test raw data (i.e. without an existing galice.root)
 *
 * Revision 1.94  2006/08/07 12:27:49  hristov
 * Removing obsolete code which affected the event numbering scheme
 *
 * Revision 1.93  2006/08/01 12:20:17  cvetan
 * 1. Adding a possibility to read and reconstruct an old rcu formatted raw data. This is controlled by an option of AliReconstruction and AliPHOSReconstructor. 2. In case of raw data processing (without galice.root) create the default AliPHOSGeometry object. Most likely this should be moved to the CDB
 *
 * Revision 1.92  2006/04/29 20:26:46  hristov
 * Separate EMC and CPV calibration (Yu.Kharlov)
 *
 * Revision 1.91  2006/04/22 10:30:17  hristov
 * Add fEnergy to AliPHOSDigit and operate with EMC amplitude in energy units (Yu.Kharlov)
 *
 * Revision 1.90  2006/04/11 15:22:59  hristov
 * run number in query set to -1: forces AliCDBManager to use its run number (A.Colla)
 *
 * Revision 1.89  2006/03/13 14:05:42  kharlov
 * Calibration objects for EMC and CPV
 *
 * Revision 1.88  2006/01/11 08:54:52  hristov
 * Additional protection in case no calibration entry was found
 *
 * Revision 1.87  2005/11/22 08:46:43  kharlov
 * Updated to new CDB (Boris Polichtchouk)
 *
 * Revision 1.86  2005/11/14 21:52:43  hristov
 * Coding conventions
 *
 * Revision 1.85  2005/09/27 16:08:08  hristov
 * New version of CDB storage framework (A.Colla)
 *
 * Revision 1.84  2005/09/21 10:02:47  kharlov
 * Reading calibration from CDB (Boris Polichtchouk)
 *
 * Revision 1.82  2005/09/02 15:43:13  kharlov
 * Add comments in GetCalibrationParameters and Calibrate
 *
 * Revision 1.81  2005/09/02 14:32:07  kharlov
 * Calibration of raw data
 *
 * Revision 1.80  2005/08/24 15:31:36  kharlov
 * Setting raw digits flag
 *
 * Revision 1.79  2005/07/25 15:53:53  kharlov
 * Read raw data
 *
 * Revision 1.78  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//*-- Author: Yves Schutz (SUBATECH)  & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//////////////////////////////////////////////////////////////////////////////
//  Clusterization class. Performs clusterization (collects neighbouring active cells) and 
//  unfolds the clusters having several local maxima.  
//  Results are stored in TreeR#, branches PHOSEmcRP (EMC recPoints),
//  PHOSCpvRP (CPV RecPoints) and AliPHOSClusterizer (Clusterizer with all 
//  parameters including input digits branch title, thresholds etc.)
//  This TTask is normally called from Reconstructioner, but can as well be used in 
//  standalone mode.
// Use Case:
//  root [0] AliPHOSClusterizerv1 * cl = new AliPHOSClusterizerv1(<pointer_to_phos_geometry_onject>)  
//  root [1] cl->Digits2Clusters(digitsTree,clusterTree)
//               //finds RecPoints in the current event
//  root [2] cl->SetDigitsBranch("digits2") 
//               //sets another title for Digitis (input) branch
//  root [3] cl->SetRecPointsBranch("recp2")  
//               //sets another title four output branches
//  root [4] cl->SetEmcLocalMaxCut(0.03)  
//               //set clusterization parameters

// --- ROOT system ---

#include "TMath.h" 
#include "TMinuit.h"
#include "TTree.h" 
#include "TBenchmark.h"
#include "TClonesArray.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliConfig.h"
#include "AliPHOSGeometry.h" 
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSDigit.h"
#include "AliPHOSDigitizer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSReconstructor.h"
#include "AliPHOSCalibData.h"

ClassImp(AliPHOSClusterizerv1)
  
//____________________________________________________________________________
AliPHOSClusterizerv1::AliPHOSClusterizerv1() :
  AliPHOSClusterizer(),
  fDefaultInit(0),            fEmcCrystals(0),          fToUnfold(0),
  fWrite(0),                  
  fNumberOfEmcClusters(0),    fNumberOfCpvClusters(0),
  fEmcClusteringThreshold(0), fCpvClusteringThreshold(0), 
  fEmcLocMaxCut(0),           fW0(0),                   fCpvLocMaxCut(0),
  fW0CPV(0),                 
  fTimeGateLowAmp(0.),        fTimeGateLow(0.),         fTimeGateHigh(0.),  
  fEcoreRadius(0)
{
  // default ctor (to be used mainly by Streamer)
  
  fDefaultInit = kTRUE ;
  
  for(Int_t i=0; i<53760; i++){
    fDigitsUsed[i]=0 ;
  }
}

//____________________________________________________________________________
AliPHOSClusterizerv1::AliPHOSClusterizerv1(AliPHOSGeometry *geom) :
  AliPHOSClusterizer(geom),
  fDefaultInit(0),            fEmcCrystals(0),          fToUnfold(0),
  fWrite(0),                
  fNumberOfEmcClusters(0),    fNumberOfCpvClusters(0),
  fEmcClusteringThreshold(0), fCpvClusteringThreshold(0), 
  fEmcLocMaxCut(0),           fW0(0),                   fCpvLocMaxCut(0),
  fW0CPV(0),                  
  fTimeGateLowAmp(0.),        fTimeGateLow(0.),         fTimeGateHigh(0.),  
  fEcoreRadius(0) 
{
  // ctor with the indication of the file where header Tree and digits Tree are stored
  
  for(Int_t i=0; i<53760; i++){
    fDigitsUsed[i]=0 ;
  }
  
  Init() ;
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________
  AliPHOSClusterizerv1::~AliPHOSClusterizerv1()
{
  // dtor

}
//____________________________________________________________________________
void AliPHOSClusterizerv1::Digits2Clusters(Option_t *option)
{
  // Steering method to perform clusterization for one event
  // The input is the tree with digits.
  // The output is the tree with clusters.

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSClusterizer"); 
  
  if(strstr(option,"print")) {
    Print() ; 
    return ;
  }

  MakeClusters() ;
    
  AliDebug(2,Form(" ---- Printing clusters (%d)\n",
		  fEMCRecPoints->GetEntries()));
  if(AliLog::GetGlobalDebugLevel()>1)
    fEMCRecPoints->Print();

  if(fToUnfold)             
    MakeUnfolding();

  WriteRecPoints();

  if(strstr(option,"deb"))  
    PrintRecPoints(option) ;

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSClusterizer");
    AliInfo(Form("took %f seconds for Clusterizing\n",
		 gBenchmark->GetCpuTime("PHOSClusterizer"))); 
  }
  fEMCRecPoints->Delete();
  fCPVRecPoints->Delete();
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::FindFit(AliPHOSEmcRecPoint * emcRP, AliPHOSDigit ** maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters) const
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 
  // The initial values for fitting procedure are set equal to the positions of local maxima.
  // Cluster will be fitted as a superposition of nPar/3 electromagnetic showers

  
  if(!gMinuit) //it was deleted by someone else
    gMinuit = new TMinuit(100) ;
  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliPHOSClusterizerv1::UnfoldingChiSquare) ;  
                                         // To set the address of the minimization function 

  TList * toMinuit = new TList();
  toMinuit->AddAt(emcRP,0) ;
  toMinuit->AddAt(fDigitsArr,1) ;
  toMinuit->AddAt(fGeom,2) ;
  
  gMinuit->SetObjectFit(toMinuit) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliPHOSDigit * digit ;

  Int_t ierflg  = 0; 
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;

  for(iDigit = 0; iDigit < nDigits; iDigit++){
    digit = maxAt[iDigit]; 

    Int_t relid[4] ;
    Float_t x = 0.;
    Float_t z = 0.;
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    fGeom->RelPosInModule(relid, x, z) ;

    Float_t energy = maxAtEnergy[iDigit] ;

    gMinuit->mnparm(index, "x",  x, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){ 
      Warning("FindFit", "PHOS Unfolding unable to set initial value for fit procedure : x = %f\n", x ) ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "z",  z, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){
       Warning("FindFit", "PHOS Unfolding unable to set initial value for fit procedure : z =%f\n", z ) ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;
    index++ ;   
    if(ierflg != 0){
      Warning("FindFit", "PHOS Unfolding unable to set initial value for fit procedure : energy = %f\n", energy ) ;      
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
    Warning("FindFit", "PHOS Unfolding fit not converged, cluster abandoned\n" );      
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
void AliPHOSClusterizerv1::Init()
{
  // Make all memory allocations which can not be done in default constructor.
  // Attach the Clusterizer task to the list of PHOS tasks
 
  fEmcCrystals = fGeom->GetNModules() *  fGeom->GetNCristalsInModule() ;

  if(!gMinuit) 
    gMinuit = new TMinuit(100);

  if (!fgCalibData) 
    fgCalibData = new AliPHOSCalibData(-1); //use AliCDBManager's run number
  if (fgCalibData->GetCalibDataEmc() == 0)
    AliFatal("Calibration parameters for PHOS EMC not found. Stop reconstruction.\n");
  if (fgCalibData->GetCalibDataCpv() == 0)   
    AliFatal("Calibration parameters for PHOS CPV not found. Stop reconstruction.\n");   

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::InitParameters()
{

  fNumberOfCpvClusters     = 0 ; 
  fNumberOfEmcClusters     = 0 ; 

  const AliPHOSRecoParam* recoParam = AliPHOSReconstructor::GetRecoParam();
  if(!recoParam) AliFatal("Reconstruction parameters are not set!");

  recoParam->Print();

  fEmcClusteringThreshold  = recoParam->GetEMCClusteringThreshold();
  fCpvClusteringThreshold  = recoParam->GetCPVClusteringThreshold();
  
  fEmcLocMaxCut            = recoParam->GetEMCLocalMaxCut();
  fCpvLocMaxCut            = recoParam->GetCPVLocalMaxCut();

  fW0                      = recoParam->GetEMCLogWeight();
  fW0CPV                   = recoParam->GetCPVLogWeight();

  fTimeGateLowAmp          = recoParam->GetTimeGateAmpThresh() ;
  fTimeGateLow             = recoParam->GetTimeGateLow() ;
  fTimeGateHigh            = recoParam->GetTimeGateHigh() ;

  fEcoreRadius             = recoParam->GetEMCEcoreRadius();
  
  fToUnfold                = recoParam->EMCToUnfold() ;
    
  fWrite                   = kTRUE ;
}

//____________________________________________________________________________
Int_t AliPHOSClusterizerv1::AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)const
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 are not neighbour but do not continue searching
  //                                       =-1 are not neighbour, continue searching, but do not look before d2 next time
  // neighbours are defined as digits having at least a common vertex 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  

  Int_t relid1[4] ; 
  fGeom->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  fGeom->AbsToRelNumbering(d2->GetId(), relid2) ; 
 
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same PHOS module 
    Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( coldiff <= 1 )  && ( rowdiff <= 1 )){   //At least common vertex
      //    if (( relid1[2]==relid2[2] && coldiff <= 1 )  || ( relid1[3]==relid2[3] &&  rowdiff <= 1 )){ //common side
      if((relid1[1] != 0) || CheckTimeGate(d1->GetTime(),d1->GetEnergy(),d2->GetTime(),d2->GetEnergy())) 
      return 1 ; 
    }
    else {
      if((relid2[2] > relid1[2]) && (relid2[3] > relid1[3]+1)) 
        return 2; //  Difference in row numbers is too large to look further 
    }
    return 0 ;

  } 
  else {
    if(relid1[0] > relid2[0] && relid1[1]==relid2[1] ) //we switched to the next module
      return -1 ;
    if(relid1[1] < relid2[1]) //we switched from EMC(0) to CPV(-1)
      return -1 ;
    
    return 2 ;

  }

  return 0 ; 
}
//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::CheckTimeGate(Float_t t1, Float_t amp1, Float_t t2, Float_t amp2)const{
  //Check if two cells have reasonable time difference
  //Note that at low amplitude time is defined up to 1 tick == 100 ns.
  if(amp1<fTimeGateLowAmp || amp2<fTimeGateLowAmp){
   return (TMath::Abs(t1 - t2 ) < fTimeGateLow) ;
  }
  else{ //Time should be measured with good accuracy
   return (TMath::Abs(t1 - t2 ) < fTimeGateHigh) ;
  }

}
//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInEmc(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-EMC module
 
  Bool_t rv = kFALSE ; 

  Int_t nEMC = fGeom->GetNModules()*fGeom->GetNPhi()*fGeom->GetNZ();  

  if(digit->GetId() <= nEMC )   rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInCpv(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-CPV module
 
  Bool_t rv = kFALSE ; 
  
  Int_t nEMC = fGeom->GetNModules()*fGeom->GetNPhi()*fGeom->GetNZ();

  if(digit->GetId() > nEMC )   rv = kTRUE;

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::WriteRecPoints()
{

  // Creates new branches with given title
  // fills and writes into TreeR.

  Int_t index ;
  //Evaluate position, dispersion and other RecPoint properties..
  Int_t nEmc = fEMCRecPoints->GetEntriesFast();
  Float_t emcMinE= AliPHOSReconstructor::GetRecoParam()->GetEMCMinE(); //Minimal digit energy
  TVector3 fakeVtx(0.,0.,0.) ;
  for(index = 0; index < nEmc; index++){
    AliPHOSEmcRecPoint * rp =
      static_cast<AliPHOSEmcRecPoint *>( fEMCRecPoints->At(index) );
    rp->Purify(emcMinE) ;
    if(rp->GetMultiplicity()==0){
      fEMCRecPoints->RemoveAt(index) ;
      delete rp ;
      continue;
    }

    // No vertex is available now, calculate corrections in PID
    rp->EvalAll(fDigitsArr) ;
    rp->EvalCoreEnergy(fW0,fEcoreRadius,fDigitsArr) ;
    rp->EvalAll(fW0,fakeVtx,fDigitsArr) ;
    rp->EvalLocal2TrackingCSTransform();
  }
  fEMCRecPoints->Compress() ;
  fEMCRecPoints->Sort() ; 
  //  fEMCRecPoints->Expand(fEMCRecPoints->GetEntriesFast()) ;
  for(index = 0; index < fEMCRecPoints->GetEntries(); index++){
    static_cast<AliPHOSEmcRecPoint *>( fEMCRecPoints->At(index) )->SetIndexInList(index) ;
  }
  
  //For each rec.point set the distance to the nearest bad crystal (BVP)
  SetDistancesToBadChannels();

  //Now the same for CPV
  for(index = 0; index < fCPVRecPoints->GetEntries(); index++){
    AliPHOSCpvRecPoint * rp = static_cast<AliPHOSCpvRecPoint *>( fCPVRecPoints->At(index) );
    rp->EvalAll(fDigitsArr) ;
    rp->EvalAll(fW0CPV,fakeVtx,fDigitsArr) ;
    rp->EvalLocal2TrackingCSTransform();
  }
  fCPVRecPoints->Sort() ;
  
  for(index = 0; index < fCPVRecPoints->GetEntries(); index++)
    static_cast<AliPHOSCpvRecPoint *>( fCPVRecPoints->At(index) )->SetIndexInList(index) ;
  
  fCPVRecPoints->Expand(fCPVRecPoints->GetEntriesFast()) ;
  
  if(fWrite){ //We write TreeR
    fTreeR->Fill();
  }
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeClusters()
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits

  fNumberOfCpvClusters     = 0 ;
  fNumberOfEmcClusters     = 0 ;

  //Mark all digits as unused yet
  const Int_t maxNDigits = 3584; // There is no clusters larger than PHOS module ;)
  Int_t nDigits=fDigitsArr->GetEntriesFast() ;

  for(Int_t i=0; i<nDigits; i++){
    fDigitsUsed[i]=0 ;
  }
  Int_t iFirst = 0 ; //first index of digit which potentially can be a part of cluster
                     //e.g. first digit in this module, first CPV digit etc.
  AliPHOSDigit * digit ; 
  TArrayI clusterdigitslist(maxNDigits) ;   
  AliPHOSRecPoint * clu = 0 ; 
  for(Int_t i=0; i<nDigits; i++){
    if(fDigitsUsed[i])
      continue ;

    digit=static_cast<AliPHOSDigit*>(fDigitsArr->At(i)) ;

    clu=0 ;

    Int_t index ;

    //is this digit so energetic that start cluster?
    if (( IsInEmc(digit) &&  Calibrate(digit->GetEnergy(),digit->GetId()) > fEmcClusteringThreshold ) || 
        ( IsInCpv(digit) &&  Calibrate(digit->GetEnergy(),digit->GetId()) > fCpvClusteringThreshold ) ) {
      Int_t iDigitInCluster = 0 ; 
      if  ( IsInEmc(digit) ) {   
        // start a new EMC RecPoint
        if(fNumberOfEmcClusters >= fEMCRecPoints->GetSize()) 
          fEMCRecPoints->Expand(2*fNumberOfEmcClusters+1) ;
          
        fEMCRecPoints->AddAt(new  AliPHOSEmcRecPoint(""), fNumberOfEmcClusters) ;
        clu = static_cast<AliPHOSEmcRecPoint *>( fEMCRecPoints->At(fNumberOfEmcClusters) ) ; 
	fNumberOfEmcClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetEnergy(),digit->GetId()),CalibrateT(digit->GetTime(),digit->GetId())) ;
        clusterdigitslist[iDigitInCluster] = digit->GetIndexInList() ;
        iDigitInCluster++ ;
        fDigitsUsed[i]=kTRUE ; 
      } else { 
        // start a new CPV cluster
        if(fNumberOfCpvClusters >= fCPVRecPoints->GetSize()) 
          fCPVRecPoints->Expand(2*fNumberOfCpvClusters+1);

        fCPVRecPoints->AddAt(new AliPHOSCpvRecPoint(""), fNumberOfCpvClusters) ;
        clu =  static_cast<AliPHOSCpvRecPoint *>( fCPVRecPoints->At(fNumberOfCpvClusters) ) ;  
        fNumberOfCpvClusters++ ; 
        clu->AddDigit(*digit,  Calibrate(digit->GetEnergy(),digit->GetId()),0.) ; // no timing information in CPV
        clusterdigitslist[iDigitInCluster] = digit->GetIndexInList()  ;        
        iDigitInCluster++ ; 
        fDigitsUsed[i]=kTRUE ;
      } // else        

      //Now scan remaining digits in list to find neigbours of our seed
 
      AliPHOSDigit * digitN ; 
      index = 0 ;
      while (index < iDigitInCluster){ // scan over digits already in cluster 
        digit =  static_cast<AliPHOSDigit*>( fDigitsArr->At(clusterdigitslist[index]) )  ;      
        index++ ; 
        for(Int_t j=iFirst; j<nDigits; j++){
	  if (iDigitInCluster >= maxNDigits) {
	    AliError(Form("The number of digits in cluster is more than %d, skip the rest of event",
			  maxNDigits));
	    return;
	  }
          if(fDigitsUsed[j]) 
            continue ;        //look through remaining digits
          digitN = static_cast<AliPHOSDigit*>( fDigitsArr->At(j) ) ;
          Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
          switch (ineb ) {
          case -1:   //too early (e.g. previous module), do not look before j at subsequent passes
            iFirst=j ;
            break ;
          case 0 :   // not a neighbour
            break ;
          case 1 :   // are neighbours 
	    clu->AddDigit(*digitN, Calibrate(digitN->GetEnergy(),digitN->GetId()),CalibrateT(digitN->GetTime(),digitN->GetId())) ;
            clusterdigitslist[iDigitInCluster] = j ; 
            iDigitInCluster++ ; 
            fDigitsUsed[j]=kTRUE ;
            break ;
          case 2 :   // too far from each other
            goto endOfLoop;   
          } // switch
          
        }
        
        endOfLoop: ; //scanned all possible neighbours for this digit
        
      } // loop over cluster     
    } // energy theshold  
  }

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeUnfolding()
{
  // Unfolds clusters using the shape of an ElectroMagnetic shower
  // Performs unfolding of all EMC/CPV clusters

  // Unfold first EMC clusters 
  if(fNumberOfEmcClusters > 0){

    Int_t nModulesToUnfold = fGeom->GetNModules() ; 

    Int_t numberofNotUnfolded = fNumberOfEmcClusters ; 
    Int_t index ;   
    for(index = 0 ; index < numberofNotUnfolded ; index++){
      
      AliPHOSEmcRecPoint * emcRecPoint = static_cast<AliPHOSEmcRecPoint *>( fEMCRecPoints->At(index) ) ;
      if(emcRecPoint->GetPHOSMod()> nModulesToUnfold)
        break ;
      
      Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
      AliPHOSDigit ** maxAt = new AliPHOSDigit*[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fEmcLocMaxCut,fDigitsArr) ;
      
      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       
        UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;

        fEMCRecPoints->Remove(emcRecPoint); 
        fEMCRecPoints->Compress() ;
        index-- ;
        fNumberOfEmcClusters -- ;
        numberofNotUnfolded-- ;
      }
      else{
        emcRecPoint->SetNExMax(1) ; //Only one local maximum
      }
      
      delete[] maxAt ; 
      delete[] maxAtEnergy ; 
    }
  } 
  // Unfolding of EMC clusters finished


  // Unfold now CPV clusters
  if(fNumberOfCpvClusters > 0){
    
    Int_t nModulesToUnfold = fGeom->GetNModules() ;

    Int_t numberofCpvNotUnfolded = fNumberOfCpvClusters ;     
    Int_t index ;   
    for(index = 0 ; index < numberofCpvNotUnfolded ; index++){
      
      AliPHOSRecPoint * recPoint = static_cast<AliPHOSRecPoint *>( fCPVRecPoints->At(index) ) ;

      if(recPoint->GetPHOSMod()> nModulesToUnfold)
        break ;
      
      AliPHOSEmcRecPoint * emcRecPoint = static_cast<AliPHOSEmcRecPoint*>(recPoint) ; 
      
      Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
      AliPHOSDigit ** maxAt = new AliPHOSDigit*[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fCpvLocMaxCut,fDigitsArr) ;
      
      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       
        UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;
        fCPVRecPoints->Remove(emcRecPoint); 
        fCPVRecPoints->Compress() ;
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
Double_t  AliPHOSClusterizerv1::ShowerShape(Double_t x, Double_t z)
{ 
  // Shape of the shower (see PHOS TDR)
  // If you change this function, change also the gradient evaluation in ChiSquare()

  //for the moment we neglect dependence on the incident angle.  

  Double_t r2    = x*x + z*z ;
  Double_t r4    = r2*r2 ;
  Double_t r295  = TMath::Power(r2, 2.95/2.) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliPHOSClusterizerv1::UnfoldCluster(AliPHOSEmcRecPoint * iniEmc, 
                                          Int_t nMax, 
                                          AliPHOSDigit ** maxAt, 
                                          Float_t * maxAtEnergy)
{ 
  // Performs the unfolding of a cluster with nMax overlapping showers 

  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;

  Bool_t rv = FindFit(iniEmc, maxAt, maxAtEnergy, nPar, fitparameters) ;

  if( !rv ) {
    // Fit failed, return and remove cluster
    iniEmc->SetNExMax(-1) ;
    delete[] fitparameters ; 
    return ;
  }

  // create ufolded rec points and fill them with new energy lists
  // First calculate energy deposited in each sell in accordance with fit (without fluctuations): efit[]
  // and later correct this number in acordance with actual energy deposition

  Int_t nDigits = iniEmc->GetMultiplicity() ;  
  Float_t * efit = new Float_t[nDigits] ;
  Float_t xDigit=0.,zDigit=0. ;
  Float_t xpar=0.,zpar=0.,epar=0.  ;
  Int_t relid[4] ;
  AliPHOSDigit * digit = 0 ;
  Int_t * emcDigits = iniEmc->GetDigitsList() ;

  TVector3 vIncid ;

  Int_t iparam ;
  Int_t iDigit ;
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
    digit = static_cast<AliPHOSDigit*>( fDigitsArr->At(emcDigits[iDigit] ) ) ;   
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    fGeom->RelPosInModule(relid, xDigit, zDigit) ;
    efit[iDigit] = 0;

    iparam = 0 ;    
    while(iparam < nPar ){
      xpar = fitparameters[iparam] ;
      zpar = fitparameters[iparam+1] ;
      epar = fitparameters[iparam+2] ;
      iparam += 3 ;
//      fGeom->GetIncidentVector(fVtx,relid[0],xpar,zpar,vIncid) ;
//      efit[iDigit] += epar * ShowerShape(xDigit - xpar,zDigit - zpar,vIncid) ;
      efit[iDigit] += epar * ShowerShape(xDigit - xpar,zDigit - zpar) ;
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
//    fGeom->GetIncidentVector(fVtx,iniEmc->GetPHOSMod(),xpar,zpar,vIncid) ;

    AliPHOSEmcRecPoint * emcRP = 0 ;  

    if(iniEmc->IsEmc()){ //create new entries in fEMCRecPoints...
      
      if(fNumberOfEmcClusters >= fEMCRecPoints->GetSize())
        fEMCRecPoints->Expand(2*fNumberOfEmcClusters) ;
      
      (*fEMCRecPoints)[fNumberOfEmcClusters] = new AliPHOSEmcRecPoint("") ;
      emcRP = static_cast<AliPHOSEmcRecPoint *>( fEMCRecPoints->At(fNumberOfEmcClusters) ) ;
      fNumberOfEmcClusters++ ;
      emcRP->SetNExMax((Int_t)nPar/3) ;
    }
    else{//create new entries in fCPVRecPoints
      if(fNumberOfCpvClusters >= fCPVRecPoints->GetSize())
        fCPVRecPoints->Expand(2*fNumberOfCpvClusters) ;
      
      (*fCPVRecPoints)[fNumberOfCpvClusters] = new AliPHOSCpvRecPoint("") ;
      emcRP = static_cast<AliPHOSEmcRecPoint *>( fCPVRecPoints->At(fNumberOfCpvClusters) ) ;
      fNumberOfCpvClusters++ ;
    }
    
    Float_t eDigit ;
    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
      digit = static_cast<AliPHOSDigit*>( fDigitsArr->At( emcDigits[iDigit] ) ) ; 
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      fGeom->RelPosInModule(relid, xDigit, zDigit) ;
//      ratio = epar * ShowerShape(xDigit - xpar,zDigit - zpar,vIncid) / efit[iDigit] ; 
      ratio = epar * ShowerShape(xDigit - xpar,zDigit - zpar) / efit[iDigit] ; 
      eDigit = emcEnergies[iDigit] * ratio ;
      emcRP->AddDigit( *digit, eDigit,CalibrateT(digit->GetTime(),digit->GetId()) ) ;
    }        
  }
 
  delete[] fitparameters ; 
  delete[] efit ; 
  
}

//_____________________________________________________________________________
void AliPHOSClusterizerv1::UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  TList * toMinuit = static_cast<TList*>( gMinuit->GetObjectFit() ) ;

  AliPHOSEmcRecPoint * emcRP = static_cast<AliPHOSEmcRecPoint*>( toMinuit->At(0) )  ;
  TClonesArray * digits = static_cast<TClonesArray*>( toMinuit->At(1) )  ;
  // A bit buggy way to get an access to the geometry
  // To be revised!
  AliPHOSGeometry *geom = static_cast<AliPHOSGeometry *>(toMinuit->At(2));

//  TVector3 * vtx = static_cast<TVector3*>(toMinuit->At(3)) ;  //Vertex position
  
  //  AliPHOSEmcRecPoint * emcRP = static_cast<AliPHOSEmcRecPoint *>( gMinuit->GetObjectFit() ) ; // EmcRecPoint to fit

  Int_t * emcDigits     = emcRP->GetDigitsList() ;

  Int_t nOdigits = emcRP->GetDigitsMultiplicity() ; 

  Float_t * emcEnergies = emcRP->GetEnergiesList() ;
  
//  TVector3 vInc ;
  fret = 0. ;     
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < nPar ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t efit ;    

  AliPHOSDigit * digit ;
  Int_t iDigit ;

  for( iDigit = 0 ; iDigit < nOdigits ; iDigit++) {

    digit = static_cast<AliPHOSDigit*>( digits->At( emcDigits[iDigit] ) ); 

    Int_t relid[4] ;
    Float_t xDigit ;
    Float_t zDigit ;

    geom->AbsToRelNumbering(digit->GetId(), relid) ;

    geom->RelPosInModule(relid, xDigit, zDigit) ;

     if(iflag == 2){  // calculate gradient
       Int_t iParam = 0 ;
       efit = 0 ;
       while(iParam < nPar ){
         Double_t dx = (xDigit - x[iParam]) ;
         iParam++ ; 
         Double_t dz = (zDigit - x[iParam]) ; 
         iParam++ ;          
//         fGeom->GetIncidentVector(*vtx,emcRP->GetPHOSMod(),x[iParam-2],x[iParam-1],vInc) ;
//         efit += x[iParam] * ShowerShape(dx,dz,vInc) ;
         efit += x[iParam] * ShowerShape(dx,dz) ;
         iParam++ ;
       }
       Double_t sum = 2. * (efit - emcEnergies[iDigit]) / emcEnergies[iDigit] ; // Here we assume, that sigma = sqrt(E) 
       iParam = 0 ;
       while(iParam < nPar ){
         Double_t xpar = x[iParam] ;
         Double_t zpar = x[iParam+1] ;
         Double_t epar = x[iParam+2] ;
//         fGeom->GetIncidentVector(*vtx,emcRP->GetPHOSMod(),xpar,zpar,vInc) ;
         Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );
//         Double_t shape = sum * ShowerShape(xDigit - xpar,zDigit - zpar,vInc) ;
         Double_t shape = sum * ShowerShape(xDigit - xpar,zDigit - zpar) ;
//DP: No incident angle dependence in gradient yet!!!!!!
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
//       fGeom->GetIncidentVector(*vtx,emcRP->GetPHOSMod(),xpar,zpar,vInc) ;
//       efit += epar * ShowerShape(xDigit - xpar,zDigit - zpar,vInc) ;
       efit += epar * ShowerShape(xDigit - xpar,zDigit - zpar) ;
     }

     fret += (efit-emcEnergies[iDigit])*(efit-emcEnergies[iDigit])/emcEnergies[iDigit] ; 
     // Here we assume, that sigma = sqrt(E)
  }

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::Print(const Option_t *)const
{
  // Print clusterizer parameters

  TString message ; 
  TString taskName(GetName()) ; 
  taskName.ReplaceAll(Version(), "") ;
  
  if( strcmp(GetName(), "") !=0 ) {  
    // Print parameters
    message  = "\n--------------- %s %s -----------\n" ; 
    message += "Clusterizing digits from the file: %s\n" ;
    message += "                           Branch: %s\n" ; 
    message += "                       EMC Clustering threshold = %f\n" ; 
    message += "                       EMC Local Maximum cut    = %f\n" ; 
    message += "                       EMC Logarothmic weight   = %f\n" ;
    message += "                       CPV Clustering threshold = %f\n" ;
    message += "                       CPV Local Maximum cut    = %f\n" ;
    message += "                       CPV Logarothmic weight   = %f\n" ;
    if(fToUnfold)
      message += " Unfolding on\n" ;
    else
      message += " Unfolding off\n" ;
    
    message += "------------------------------------------------------------------" ;
  }
  else 
    message = " AliPHOSClusterizerv1 not initialized " ;
  
  AliInfo(Form("%s, %s %s %s %s %f %f %f %f %f %f", message.Data(),  
       taskName.Data(), 
       GetTitle(),
       taskName.Data(), 
       GetName(), 
       fEmcClusteringThreshold, 
       fEmcLocMaxCut, 
       fW0, 
       fCpvClusteringThreshold, 
       fCpvLocMaxCut, 
       fW0CPV )) ; 
}
//____________________________________________________________________________
void AliPHOSClusterizerv1::PrintRecPoints(Option_t * option)
{
  // Prints list of RecPoints produced at the current pass of AliPHOSClusterizer

  AliInfo(Form("\nFound %d EMC RecPoints and %d CPV RecPoints", 
	       fEMCRecPoints->GetEntriesFast(),
	       fCPVRecPoints->GetEntriesFast() ))  ;
 
  if(strstr(option,"all")) {
    printf("\n EMC clusters \n") ;
    printf("Index    Ene(MeV) Multi Module     X    Y   Z    Lambdas_1  Lambda_2  # of prim  Primaries list\n") ;      
    Int_t index ;
    for (index = 0 ; index < fEMCRecPoints->GetEntries() ; index++) {
      AliPHOSEmcRecPoint * rp = (AliPHOSEmcRecPoint * )fEMCRecPoints->At(index) ; 
      TVector3  locpos;  
      rp->GetLocalPosition(locpos);
      Float_t lambda[2]; 
      rp->GetElipsAxis(lambda);
      Int_t * primaries; 
      Int_t nprimaries;
      primaries = rp->GetPrimaries(nprimaries);
      printf("\n%6d  %8.2f  %3d     %2d     %4.1f   %4.1f   %4.1f   %4f  %4f    %2d     : ", 
	      rp->GetIndexInList(), rp->GetEnergy(), rp->GetMultiplicity(), rp->GetPHOSMod(), 
	      locpos.X(), locpos.Y(), locpos.Z(), lambda[0], lambda[1], nprimaries) ; 
      
      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++) {
	printf("%d ", primaries[iprimary] ) ; 
      }
      printf("\n") ;
    }

    //Now plot CPV recPoints
    printf("\n CPV clusters \n") ;
    printf("Index    Ene(MeV) Module     X     Y    Z  \n") ;      
    for (index = 0 ; index < fCPVRecPoints->GetEntries() ; index++) {
      AliPHOSCpvRecPoint * rp = (AliPHOSCpvRecPoint * )fCPVRecPoints->At(index) ; 
      
      TVector3  locpos;  
      rp->GetLocalPosition(locpos);
      
      printf("\n%6d  %8.2f  %2d     %4.1f    %4.1f %4.1f \n", 
	     rp->GetIndexInList(), rp->GetEnergy(), rp->GetPHOSMod(), 
	     locpos.X(), locpos.Y(), locpos.Z()) ; 
    }
  }
}


//____________________________________________________________________________
void AliPHOSClusterizerv1::SetDistancesToBadChannels()
{
  //For each EMC rec. point set the distance to the nearest bad crystal.
  //Author: Boris Polichtchouk 

  if(!fgCalibData->GetNumOfEmcBadChannels()) return;

  Int_t badIds[8000];
  memset(badIds,0,8000*sizeof(Int_t));
  fgCalibData->EmcBadChannelIds(badIds);

  AliPHOSEmcRecPoint* rp;

  TMatrixF gmat;
  TVector3 gposRecPoint; // global (in ALICE frame) position of rec. point
  TVector3 gposBadChannel; // global position of bad crystal
  TVector3 dR;

  Float_t dist,minDist;
  Int_t relid[4]={0,0,0,0} ;
  TVector3 lpos ;
  for(Int_t iRP=0; iRP<fEMCRecPoints->GetEntries(); iRP++){
    rp = (AliPHOSEmcRecPoint*)fEMCRecPoints->At(iRP);
    //evaluate distance to border
    relid[0]=rp->GetPHOSMod() ;
    relid[2]=1 ;
    relid[3]=1 ;
    Float_t xcorner,zcorner;
    fGeom->RelPosInModule(relid,xcorner,zcorner) ; //coordinate of the corner cell
    rp->GetLocalPosition(lpos) ;
    minDist = 2.2+TMath::Min(-xcorner-TMath::Abs(lpos.X()),-zcorner-TMath::Abs(lpos.Z())); //2.2 - crystal size
    for(Int_t iBad=0; iBad<fgCalibData->GetNumOfEmcBadChannels(); iBad++) {
      fGeom->AbsToRelNumbering(badIds[iBad],relid)  ;
      if(relid[0]!=rp->GetPHOSMod()) //We can not evaluate global position directly since 
        continue ;                   //bad channels can be in the module which does not exist in simulations.
      rp->GetGlobalPosition(gposRecPoint,gmat);
      fGeom->RelPosInAlice(badIds[iBad],gposBadChannel);
      AliDebug(2,Form("BC position:[%.3f,%.3f,%.3f], RP position:[%.3f,%.3f,%.3f]. E=%.3f\n",
		      gposBadChannel.X(),gposBadChannel.Y(),gposBadChannel.Z(),
		      gposRecPoint.X(),gposRecPoint.Y(),gposRecPoint.Z(),rp->GetEnergy()));
      dR = gposBadChannel-gposRecPoint;
      dist = dR.Mag();
      if(dist<minDist) minDist = dist;
    }

    rp->SetDistanceToBadCrystal(minDist); 
  }

}
//==================================================================================
Float_t AliPHOSClusterizerv1::Calibrate(Float_t amp, Int_t absId) const{
  // Calibrate EMC digit, i.e. multiply its Amp by a factor read from CDB

  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  if(relId[1]){ //CPV
    Float_t calibration = fgCalibData->GetADCchannelCpv(module,column,row);
    return amp*calibration ;
  }   
  else{ //EMC
    Float_t calibration = fgCalibData->GetADCchannelEmc(module,column,row);
    return amp*calibration ;
  }
}
//==================================================================================
Float_t AliPHOSClusterizerv1::CalibrateT(Float_t time, Int_t absId)const{
  // Calibrate time in EMC digit

  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  if(relId[1]){ //CPV
    return 0. ;
  }
  else{ //EMC
    time += fgCalibData->GetTimeShiftEmc(module,column,row);
    return time ;
  }
}

