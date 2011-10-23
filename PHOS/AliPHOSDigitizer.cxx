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
 * $Log: AliPHOSDigitizer.cxx,v $
 * Revision 1.104  2007/12/18 09:08:18  hristov
 * Splitting of the QA maker into simulation and reconstruction dependent parts (Yves)
 *
 * Revision 1.103  2007/11/07 11:25:06  schutz
 * Comment out the QA checking before starting digitization
 *
 * Revision 1.102  2007/10/19 18:04:29  schutz
 * The standalone QA data maker is called from AliSimulation and AliReconstruction outside the event loop; i.e. re-reading the data. The QA data making in the event loop has been commented out.
 *
 * Revision 1.101  2007/10/14 21:08:10  schutz
 * Introduced the checking of QA results from previous step before entering the event loop
 *
 * Revision 1.100  2007/10/10 09:05:10  schutz
 * Changing name QualAss to QA
 *
 * Revision 1.99  2007/09/30 17:08:20  schutz
 * Introducing the notion of QA data acquisition cycle (needed by online)
 *
 * Revision 1.98  2007/09/26 14:22:17  cvetan
 * Important changes to the reconstructor classes. Complete elimination of the run-loaders, which are now steered only from AliReconstruction. Removal of the corresponding Reconstruct() and FillESD() methods.
 *
 * Revision 1.97  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.96  2007/04/28 10:43:36  policheh
 * Dead channels simulation: digit energy sets to 0.
 *
 * Revision 1.95  2007/04/10 07:20:52  kharlov
 * Decalibration should use the same CDB as calibration in AliPHOSClusterizerv1
 *
 * Revision 1.94  2007/02/01 10:34:47  hristov
 * Removing warnings on Solaris x86
 *
 * Revision 1.93  2006/10/17 13:17:01  kharlov
 * Replace AliInfo by AliDebug
 *
 * Revision 1.92  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
 * Revision 1.91  2006/04/29 20:25:30  hristov
 * Decalibration is implemented (Yu.Kharlov)
 *
 * Revision 1.90  2006/04/22 10:30:17  hristov
 * Add fEnergy to AliPHOSDigit and operate with EMC amplitude in energy units (Yu.Kharlov)
 *
 * Revision 1.89  2006/04/11 15:22:59  hristov
 * run number in query set to -1: forces AliCDBManager to use its run number (A.Colla)
 *
 * Revision 1.88  2006/03/13 14:05:43  kharlov
 * Calibration objects for EMC and CPV
 *
 * Revision 1.87  2005/08/24 15:33:49  kharlov
 * Calibration data for raw digits
 *
 * Revision 1.86  2005/07/12 20:07:35  hristov
 * Changes needed to run simulation and reconstrruction in the same AliRoot session
 *
 * Revision 1.85  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//*-- Author :  Dmitri Peressounko (SUBATECH & Kurchatov Institute) 
//////////////////////////////////////////////////////////////////////////////
// This class performs digitization of Summable digits (in the PHOS case it is just
// the sum of contributions from all primary particles into a given cell). 
// In addition it performs mixing of summable digits from different events.
// The name of the class is also the title of the branch that will contain 
// the created SDigits
// The title of the class is the name of the file that contains the hits from
// which the SDigits are created
//
// For each event two branches are created in TreeD:
//   "PHOS" - list of digits
//   "AliPHOSDigitizer" - AliPHOSDigitizer with all parameters used in digitization
//
// Note, that one can set a title for new digits branch, and repeat digitization with
// another set of parameters.
//
// Use case:
// root[0] AliPHOSDigitizer * d = new AliPHOSDigitizer() ;
// root[1] d->Digitize()             
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       //Digitizes SDigitis in all events found in file galice.root 
//
// root[2] AliPHOSDigitizer * d1 = new AliPHOSDigitizer("galice1.root") ;  
//                       // Will read sdigits from galice1.root
// root[3] d1->MixWith("galice2.root")       
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       // Reads another set of sdigits from galice2.root
// root[3] d1->MixWith("galice3.root")       
//                       // Reads another set of sdigits from galice3.root
// root[4] d->Digitize("deb timing")    
//                       // Reads SDigits from files galice1.root, galice2.root ....
//                       // mixes them and stores produced Digits in file galice1.root          
//                       // deb - prints number of produced digits
//                       // deb all - prints list of produced digits
//                       // timing  - prints time used for digitization
//

// --- ROOT system ---
#include "TTree.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TMath.h"

// --- Standard library ---

// --- AliRoot header files ---
#include <TGeoManager.h>                                                                                                                   
#include "AliLog.h"
#include "AliDigitizationInput.h"
#include "AliPHOSDigit.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSTick.h"
#include "AliPHOSSimParam.h"
#include "AliPHOSCalibData.h"
#include "AliRunLoader.h"
#include "AliPHOSLoader.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSDigitizer)


//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer() :
  AliDigitizer("",""),
  fDefaultInit(kTRUE),
  fDigitsInRun(0),
  fInit(kFALSE),
  fInput(0),
  fInputFileNames(0x0),
  fEventNames(0x0),
  fEmcCrystals(0),
  fEventFolderName(""),
  fFirstEvent(0),
  fLastEvent(0), 
  fcdb(0x0),
  fEventCounter(0),
  fPulse(0),
  fADCValuesLG(0),
  fADCValuesHG(0)
{
  // ctor
  InitParameters() ; 
  fDigInput = 0 ;                     // We work in the standalong mode
}

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(TString alirunFileName, 
				   TString eventFolderName):
  AliDigitizer("PHOSDigitizer", alirunFileName), 
  fDefaultInit(kFALSE),
  fDigitsInRun(0),
  fInit(kFALSE),
  fInput(0),
  fInputFileNames(0x0),
  fEventNames(0x0),
  fEmcCrystals(0),
  fEventFolderName(eventFolderName),
  fFirstEvent(0),
  fLastEvent(0), 
  fcdb(0x0),
  fEventCounter(0),
  fPulse(0),
  fADCValuesLG(0),
  fADCValuesHG(0)
{
  // ctor
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
  fDigInput = 0 ;                     // We work in the standalone mode
  fcdb = new AliPHOSCalibData(-1);
}

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(const AliPHOSDigitizer & d) : 
  AliDigitizer(d),
  fDefaultInit(d.fDefaultInit),
  fDigitsInRun(d.fDigitsInRun),
  fInit(d.fInit),
  fInput(d.fInput),
  fInputFileNames(0x0),//?
  fEventNames(0x0),//?
  fEmcCrystals(d.fEmcCrystals),
  fEventFolderName(d.fEventFolderName),
  fFirstEvent(d.fFirstEvent),
  fLastEvent(d.fLastEvent), 
  fcdb (0x0), 
  fEventCounter(0),
  fPulse(0),
  fADCValuesLG(0),
  fADCValuesHG(0)
{
  // copyy ctor 
  SetName(d.GetName()) ; 
  SetTitle(d.GetTitle()) ; 
  fcdb = new AliPHOSCalibData(-1);
}

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(AliDigitizationInput * rd) :
  AliDigitizer(rd,"PHOSDigitizer"),
  fDefaultInit(kFALSE),
  fDigitsInRun(0),
  fInit(kFALSE),
  fInput(0),
  fInputFileNames(0x0),
  fEventNames(0x0),
  fEmcCrystals(0),
  fEventFolderName(fDigInput->GetInputFolderName(0)),
  fFirstEvent(0),
  fLastEvent(0), 
  fcdb (0x0), 
  fEventCounter(0),
  fPulse(0),
  fADCValuesLG(0),
  fADCValuesHG(0)

{
  // ctor Init() is called by RunDigitizer
  fDigInput = rd ; 
  SetTitle(static_cast<AliStream*>(fDigInput->GetInputStream(0))->GetFileName(0));
  InitParameters() ; 
  fDefaultInit = kFALSE ; 
  fcdb = new AliPHOSCalibData(-1);
}

//____________________________________________________________________________ 
  AliPHOSDigitizer::~AliPHOSDigitizer()
{
  // dtor  
  delete [] fInputFileNames ; 
  delete [] fEventNames ; 

  delete fPulse;
  delete [] fADCValuesLG;
  delete [] fADCValuesHG;

  if(fcdb){ delete fcdb ; fcdb=0;} 

}

//____________________________________________________________________________
void AliPHOSDigitizer::Digitize(Int_t event) 
{ 
  
  // Makes the digitization of the collected summable digits.
  //  It first creates the array of all PHOS modules
  //  filled with noise (different for EMC, and CPV) and
  //  then adds contributions from SDigits. 
  // This design avoids scanning over the list of digits to add 
  // contribution to new SDigits only.

  Bool_t toMakeNoise = kTRUE ; //Do not create noisy digits if merge with real data

  //First stream 
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));

  Int_t readEvent = event ; 
  if (fDigInput) 
    readEvent = static_cast<AliStream*>(fDigInput->GetInputStream(0))->GetCurrentEventNumber() ; 
  AliDebug(1,Form("Adding event %d from input stream 0 %s %s", 
		  readEvent, GetTitle(), fEventFolderName.Data())) ; 
  rl->GetEvent(readEvent) ;
  phosLoader->CleanSDigits() ; 
  phosLoader->LoadSDigits("READ") ;

  //Prepare Output
  TClonesArray * digits = phosLoader->Digits() ;
  if( !digits ) {
    phosLoader->MakeDigitsArray() ; 
    digits = phosLoader->Digits() ;
  }
 
  digits->Clear() ;

  //
  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;
  //Making digits with noise, first EMC
  //Check which PHOS modules are present
  Bool_t isPresent[5] ;
  TString volpath ;
  Int_t nmod=0 ;
  for(Int_t i=0; i<5; i++){
    volpath = "/ALIC_1/PHOS_";
    volpath += i+1;
    if (gGeoManager->CheckPath(volpath.Data())) {
      isPresent[i]=1 ;
      nmod++ ;
    }
    else{
      isPresent[i]=0 ;
    }
  }

  Int_t nEMC = nmod*geom->GetNPhi()*geom->GetNZ();
  
  Int_t nCPV ;
  Int_t absID ;
  
  //check if CPV exists
  Bool_t isCPVpresent=0 ;
  for(Int_t i=1; i<=5 && !isCPVpresent; i++){
    volpath = "/ALIC_1/PHOS_";
    volpath += i;
    volpath += "/PCPV_1";
    if (gGeoManager->CheckPath(volpath.Data())) 
      isCPVpresent=1 ;
  } 
  
  if(isCPVpresent){
    nCPV = nEMC + geom->GetNumberOfCPVPadsZ() * geom->GetNumberOfCPVPadsPhi() * nmod ;
  }
  else{
     nCPV = nEMC ;
  }  

  digits->Expand(nCPV) ;

  //take all the inputs to add together and load the SDigits
  TObjArray * sdigArray = new TObjArray(fInput) ;
  sdigArray->AddAt(phosLoader->SDigits(), 0) ;
 
  for(Int_t i = 1 ; i < fInput ; i++){
    TString tempo(fEventNames[i]) ; 
    tempo += i ;
    AliRunLoader* rl2 = AliRunLoader::GetRunLoader(tempo) ;
    if(!rl2){ 
      rl2 = AliRunLoader::Open(fInputFileNames[i], tempo) ;
      if (!rl2) {
        Fatal("AliPHOSDigitizer", "Could not find the Run Loader for %s - %s",fInputFileNames[i].Data(), tempo.Data()) ;
        return ;
      }
      rl2->LoadHeader();
    }
    AliPHOSLoader * phosLoader2 = static_cast<AliPHOSLoader*>(rl2->GetLoader("PHOSLoader"));
 
    if(fDigInput){ 
      readEvent = static_cast<AliStream*>(fDigInput->GetInputStream(i))->GetCurrentEventNumber() ; 
    }
    TClonesArray * digs ;
    if(AliPHOSSimParam::GetInstance()->IsStreamDigits(i)){ //This is Digits Stream
      AliInfo(Form("Adding event %d from input stream %d %s %s", 
		 readEvent, i, fInputFileNames[i].Data(), tempo.Data())) ; 
       rl2->GetEvent(readEvent) ;
       phosLoader2->CleanDigits() ;
       phosLoader2->LoadDigits("READ") ;
       digs = phosLoader2->Digits() ;
       toMakeNoise=0 ; //Do not add noise, it is already in stream
    }
    else{
      AliInfo(Form("Adding event %d (SDigits) from input stream %d %s %s",
                 readEvent, i, fInputFileNames[i].Data(), tempo.Data())) ;
       rl2->GetEvent(readEvent) ;
       phosLoader2->CleanSDigits() ; 
       phosLoader2->LoadSDigits("READ") ;
       digs = phosLoader2->SDigits() ;
    } 
    sdigArray->AddAt(digs, i) ;
  }

  //Find the first crystal with signal
  Int_t nextSig = 200000 ; 
  TClonesArray * sdigits ;  
  for(Int_t i = 0 ; i < fInput ; i++){
    sdigits = static_cast<TClonesArray *>(sdigArray->At(i)) ;
    if ( !sdigits->GetEntriesFast() )
      continue ; 
    Int_t curNext = static_cast<AliPHOSDigit *>(sdigits->At(0))->GetId() ;
    if(curNext < nextSig) 
      nextSig = curNext ;
  }
  
  TArrayI index(fInput) ;
  index.Reset() ;  //Set all indexes to zero
  
  AliPHOSDigit * digit ;
  AliPHOSDigit * curSDigit ;
  
//  TClonesArray * ticks = new TClonesArray("AliPHOSTick",1000) ;
  
  //Put Noise contribution
  Double_t apdNoise = 0. ;
  if(toMakeNoise)
     apdNoise = AliPHOSSimParam::GetInstance()->GetAPDNoise() ; 

  Int_t emcpermod=geom->GetNPhi()*geom->GetNZ();
  Int_t idigit= 0;
  for(Int_t imod=0; imod<5; imod++){
    if(!isPresent[imod])
      continue ;
    Int_t firstAbsId=imod*emcpermod+1 ;
    Int_t lastAbsId =(imod+1)*emcpermod ; 
    for(absID = firstAbsId ; absID <= lastAbsId ; absID++){
      Float_t noise = gRandom->Gaus(0.,apdNoise) ; 
      new((*digits)[idigit]) AliPHOSDigit( -1, absID, noise, TimeOfNoise() ) ;
      //look if we have to add signal?
      digit = static_cast<AliPHOSDigit *>(digits->At(idigit)) ;
      idigit++ ;
    
      if(absID==nextSig){
      //Add SDigits from all inputs 
//      ticks->Clear() ;
//      Int_t contrib = 0 ;

//New Timing model is necessary
//      Float_t a = digit->GetEnergy() ;
//      Float_t b = TMath::Abs( a / fTimeSignalLength) ;
//      //Mark the beginning of the signal
//      new((*ticks)[contrib++]) AliPHOSTick(digit->GetTime(),0, b);  
//      //Mark the end of the signal     
//      new((*ticks)[contrib++]) AliPHOSTick(digit->GetTime()+fTimeSignalLength, -a, -b); 

// Calculate time as time of the largest digit
        Float_t time = digit->GetTime() ;
        Float_t eTime= digit->GetEnergy() ;
      
        //loop over inputs
        for(Int_t i = 0 ; i < fInput ; i++){
  	  if( static_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] ){
  	    curSDigit = static_cast<AliPHOSDigit*>(static_cast<TClonesArray *>(sdigArray->At(i))->At(index[i])) ; 	
            if(AliPHOSSimParam::GetInstance()->IsStreamDigits(i)){ //This is Digits Stream
              curSDigit->SetEnergy(Calibrate(curSDigit->GetEnergy(),curSDigit->GetId())) ;
              curSDigit->SetTime(CalibrateT(curSDigit->GetTime(),curSDigit->GetId())) ;
            }
          }
	  else
	    curSDigit = 0 ;
	  //May be several digits will contribute from the same input
	  while(curSDigit && curSDigit->GetId() == absID){	   
	    //Shift primary to separate primaries belonging different inputs
	    Int_t primaryoffset ;
	    if(fDigInput)
	      primaryoffset = fDigInput->GetMask(i) ; 
	    else
	      primaryoffset = 10000000*i ;
	    curSDigit->ShiftPrimary(primaryoffset) ;
	  
//New Timing model is necessary
//	  a = curSDigit->GetEnergy() ;
//	  b = a /fTimeSignalLength ;
//	  new((*ticks)[contrib++]) AliPHOSTick(curSDigit->GetTime(),0, b);  
//	  new((*ticks)[contrib++]) AliPHOSTick(curSDigit->GetTime()+fTimeSignalLength, -a, -b); 
          if(curSDigit->GetEnergy()>eTime){
            eTime=curSDigit->GetEnergy() ;
            time=curSDigit->GetTime() ;
          }
	    *digit += *curSDigit ;  //add energies

	    index[i]++ ;
	    if( static_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] )
	      curSDigit = static_cast<AliPHOSDigit*>(static_cast<TClonesArray *>(sdigArray->At(i))->At(index[i])) ; 	
	    else
	      curSDigit = 0 ;
	  }
        }
      
        //calculate and set time
//New Timing model is necessary
//      Float_t time = FrontEdgeTime(ticks) ;
        digit->SetTime(time) ;
      
        //Find next signal module
        nextSig = 200000 ;
        for(Int_t i = 0 ; i < fInput ; i++){
	  sdigits = static_cast<TClonesArray *>(sdigArray->At(i)) ;
	  Int_t curNext = nextSig ;
	  if(sdigits->GetEntriesFast() > index[i] ){
	    curNext = static_cast<AliPHOSDigit *>(sdigits->At(index[i]))->GetId() ;
	  }
	  if(curNext < nextSig) nextSig = curNext ;
        }
      }
    }
  }


  //Apply non-linearity
  if(AliPHOSSimParam::GetInstance()->IsCellNonlinearityOn()){ //Apply non-lineairyt on cell level
    const Double_t aNL = AliPHOSSimParam::GetInstance()->GetCellNonLineairyA() ;
    const Double_t bNL = AliPHOSSimParam::GetInstance()->GetCellNonLineairyB() ;
    const Double_t cNL = AliPHOSSimParam::GetInstance()->GetCellNonLineairyC() ;
    for(Int_t i = 0 ; i < nEMC ; i++){
      digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ;
      Double_t e= digit->GetEnergy() ;
      // version(1)      digit->SetEnergy(e*(1+a*TMath::Exp(-e/b))) ;
      digit->SetEnergy(e*cNL*(1.+aNL*TMath::Exp(-e*e/2./bNL/bNL))) ; //Better agreement with data...
    }  
  }


  //distretize energy if necessary
  if(AliPHOSSimParam::GetInstance()->IsEDigitizationOn()){
    Float_t adcW=AliPHOSSimParam::GetInstance()->GetADCchannelW() ;
    for(Int_t i = 0 ; i < nEMC ; i++){                                                                                                       
      digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ;
      digit->SetEnergy(adcW*ceil(digit->GetEnergy()/adcW)) ;
    } 
  }
 
  //Apply decalibration if necessary
  for(Int_t i = 0 ; i < nEMC ; i++){
    digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ;
    Decalibrate(digit) ;
  }
  
//  ticks->Delete() ;
//  delete ticks ;
  
  //Now CPV digits (different noise and no timing)
  Int_t cpvpermod = geom->GetNumberOfCPVPadsZ() * geom->GetNumberOfCPVPadsPhi() ;
  Int_t nEMCtotal=emcpermod*5 ;
  Float_t cpvNoise = AliPHOSSimParam::GetInstance()->GetCPVNoise() ;
  if(isCPVpresent){  //CPV is present in current geometry
    for(Int_t imod=0; imod<5; imod++){ //module is present in current geometry
      if(!isPresent[imod])
        continue ;
      Int_t firstAbsId=nEMCtotal+imod*cpvpermod+1 ;
      Int_t lastAbsId =nEMCtotal+(imod+1)*cpvpermod ;
      for(absID = firstAbsId; absID <= lastAbsId; absID++){
        Float_t noise = gRandom->Gaus(0., cpvNoise) ; 
        new((*digits)[idigit]) AliPHOSDigit( -1,absID,noise, TimeOfNoise() ) ;
        idigit++ ;
        //look if we have to add signal?
        if(absID==nextSig){
          digit = static_cast<AliPHOSDigit *>(digits->At(idigit-1)) ;
          //Add SDigits from all inputs
          for(Int_t i = 0 ; i < fInput ; i++){
	     if( static_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] )
	       curSDigit = static_cast<AliPHOSDigit*>( static_cast<TClonesArray *>(sdigArray->At(i))->At(index[i])) ; 	
	     else
	       curSDigit = 0 ;

	     //May be several digits will contribute from the same input
	     while(curSDigit && curSDigit->GetId() == absID){	   
	       //Shift primary to separate primaries belonging different inputs
	       Int_t primaryoffset ;
	       if(fDigInput)
	         primaryoffset = fDigInput->GetMask(i) ; 
	       else
	         primaryoffset = 10000000*i ;
	       curSDigit->ShiftPrimary(primaryoffset) ;

	       //add energies
	       *digit += *curSDigit ;  
	       index[i]++ ;
	       if( static_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] )
	         curSDigit = static_cast<AliPHOSDigit*>( static_cast<TClonesArray *>(sdigArray->At(i))->At(index[i]) ) ; 	
	       else
	         curSDigit = 0 ;
	     }
          }

          //Find next signal module
          nextSig = 200000 ;
          for(Int_t i = 0 ; i < fInput ; i++){
	    sdigits = static_cast<TClonesArray *>(sdigArray->At(i)) ;
	    Int_t curNext = nextSig ;
	    if(sdigits->GetEntriesFast() > index[i] )
	      curNext = static_cast<AliPHOSDigit *>( sdigits->At(index[i]) )->GetId() ;
	    if(curNext < nextSig) nextSig = curNext ;
          }
      
        }
      }
    }
  }

  delete sdigArray ; //We should not delete its contents 
  
  Int_t relId[4];

  //set amplitudes in bad channels to zero

  for(Int_t i = 0 ; i <digits->GetEntriesFast(); i++){
    digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ;
    geom->AbsToRelNumbering(digit->GetId(),relId);
    if(relId[1] == 0) // Emc
      if(fcdb->IsBadChannelEmc(relId[0],relId[3],relId[2])) digit->SetEnergy(0.); 
  }

  //remove digits below thresholds
  Float_t emcThreshold = AliPHOSSimParam::GetInstance()->GetEmcDigitsThreshold() ;
  for(Int_t i = 0 ; i < nEMC ; i++){
    digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ;

    if(digit->GetEnergy() < emcThreshold){
      digits->RemoveAt(i) ;
      continue ;
    }

    geom->AbsToRelNumbering(digit->GetId(),relId);

//    digit->SetEnergy(TMath::Ceil(digit->GetEnergy())-0.9999) ;

    Float_t tres = TimeResolution(digit->GetEnergy()) ; 
    digit->SetTime(gRandom->Gaus(digit->GetTime(), tres) ) ;

    fPulse->Reset();
    fPulse->SetAmplitude(digit->GetEnergy()/
			 fcdb->GetADCchannelEmc(relId[0],relId[3],relId[2]));
    fPulse->SetTZero(digit->GetTimeR());
    fPulse->MakeSamples();
    fPulse->GetSamples(fADCValuesHG, fADCValuesLG) ; 
    Int_t nSamples = fPulse->GetRawFormatTimeBins();
    digit->SetALTROSamplesHG(nSamples,fADCValuesHG);
    digit->SetALTROSamplesLG(nSamples,fADCValuesLG);
  }

  Float_t cpvDigitThreshold = AliPHOSSimParam::GetInstance()->GetCpvDigitsThreshold() ;
  for(Int_t i = nEMC; i < nCPV ; i++){
    if( static_cast<AliPHOSDigit*>(digits->At(i))->GetEnergy() < cpvDigitThreshold )
      digits->RemoveAt(i) ;
  } 
    
  digits->Compress() ;  
  Int_t ndigits = digits->GetEntriesFast() ;
  
  //Set indexes in list of digits and make true digitization of the energy
  for (Int_t i = 0 ; i < ndigits ; i++) { 
    digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ; 
    digit->SetIndexInList(i) ;     
    if(digit->GetId() > fEmcCrystals){ //digitize CPV only
      digit->SetAmp(DigitizeCPV(digit->GetEnergy(),digit->GetId()) ) ;
    }
  }

}
//____________________________________________________________________________
Float_t AliPHOSDigitizer::Calibrate(Float_t amp,Int_t absId){
  //Apply calibration
  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  if(relId[1]==0){ //This Is EMC
    Float_t calibration = fcdb->GetADCchannelEmc(module,column,row);
    return amp*calibration ;
  }
  return 0 ;
}
//____________________________________________________________________________
void AliPHOSDigitizer::Decalibrate(AliPHOSDigit *digit)
{
  // Decalibrate EMC digit, i.e. change its energy by a factor read from CDB

  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(digit->GetId(),relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  if(relId[1]==0){ //This Is EMC
    Float_t decalib     = fcdb->GetADCchannelEmcDecalib(module,column,row); // O(1)
    Float_t calibration = fcdb->GetADCchannelEmc(module,column,row)*decalib;
    Float_t energy = digit->GetEnergy()/calibration;
    digit->SetEnergy(energy); //Now digit measures E in ADC counts
    Float_t time = digit->GetTime() ;
    time-=fcdb->GetTimeShiftEmc(module,column,row);
    digit->SetTime(time) ;
  }
}
//____________________________________________________________________________
Float_t AliPHOSDigitizer::CalibrateT(Float_t time,Int_t absId){
  //Apply time calibration
  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absolute ID
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  time += fcdb->GetTimeShiftEmc(module,column,row);
  return time ;
}
//____________________________________________________________________________
Int_t AliPHOSDigitizer::DigitizeCPV(Float_t charge, Int_t absId)
{
  // Returns digitized value of the CPV charge in a pad absId

  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //Determine rel.position of the cell absId
  Int_t relId[4];
  geom->AbsToRelNumbering(absId,relId);
  Int_t module=relId[0];
  Int_t row   =relId[2];
  Int_t column=relId[3];
  
  Int_t channel = 0;
  
  if(absId > fEmcCrystals){ //digitize CPV only

    //reading calibration data for cell absId.
    Float_t adcPedestalCpv = fcdb->GetADCpedestalCpv(module,column,row);
    Float_t adcChanelCpv   = fcdb->GetADCchannelCpv( module,column,row);

    channel = (Int_t) TMath::Ceil((charge - adcPedestalCpv)/adcChanelCpv) ;       
    Int_t nMax = AliPHOSSimParam::GetInstance()->GetNADCcpv() ;
    if(channel > nMax ) channel = nMax ;
  }
  return channel ;
}

//____________________________________________________________________________
void AliPHOSDigitizer::Digitize(Option_t *option) 
{ 
  // Steering method to process digitization for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1, then process events until the end.
  // by default fLastEvent = fFirstEvent (process only one event)

  if (!fInit) { // to prevent overwrite existing file
    AliError(Form("Give a version name different from %s", 
		  fEventFolderName.Data() )) ;
    return ;
  }   

  if (strstr(option,"print")) {
    Print();
    return ; 
  }
  
 if(strstr(option,"tim"))
     gBenchmark->Start("PHOSDigitizer");
  
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
  
  if (fLastEvent == -1) 
    fLastEvent = rl->GetNumberOfEvents() - 1 ;
  else if (fDigInput) 
    fLastEvent = fFirstEvent ; 
 
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;
  
  Int_t ievent ;

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    fEventCounter++ ; 

    Digitize(ievent) ; //Add prepared SDigits to digits and add the noise

    WriteDigits() ;

    if(strstr(option,"deb"))
      PrintDigits(option);
    
    //increment the total number of Digits per run 
    fDigitsInRun += phosLoader->Digits()->GetEntriesFast() ;  
 }
 
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSDigitizer");
    TString message ; 
    message = "  took %f seconds for Digitizing %f seconds per event\n" ; 
    AliInfo(Form( message.Data(), 
	 gBenchmark->GetCpuTime("PHOSDigitizer"), 
	 gBenchmark->GetCpuTime("PHOSDigitizer")/nEvents )); 
  } 
}
//____________________________________________________________________________ 
Float_t AliPHOSDigitizer::TimeResolution(Float_t e){
  //calculate TOF resolution using beam-test resutls
  Float_t a=AliPHOSSimParam::GetInstance()->GetTOFa() ;
  Float_t b=AliPHOSSimParam::GetInstance()->GetTOFb() ;
  return TMath::Sqrt(a*a+b*b/e) ;
}

////____________________________________________________________________________ 
//Float_t AliPHOSDigitizer::FrontEdgeTime(TClonesArray * ticks) const
//{
//  // Returns the shortest time among all time ticks
//
//  ticks->Sort() ; //Sort in accordance with times of ticks
//  TIter it(ticks) ;
//  AliPHOSTick * ctick = (AliPHOSTick *) it.Next() ;
//  Float_t time = ctick->CrossingTime(fTimeThreshold) ;    
//
//  AliPHOSTick * t ;  
//  while((t=(AliPHOSTick*) it.Next())){
//    if(t->GetTime() < time)  //This tick starts before crossing
//      *ctick+=*t ;
//    else
//      return time ;
//
//    time = ctick->CrossingTime(fTimeThreshold) ;    
//  }
//  return time ;
//}

//____________________________________________________________________________ 
Bool_t AliPHOSDigitizer::Init()
{
  // Makes all memory allocations
  fInit = kTRUE ; 

  AliPHOSGeometry *geom;
  if (!(geom = AliPHOSGeometry::GetInstance())) 
        geom = AliPHOSGeometry::GetInstance("IHEP","");
//   const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  fEmcCrystals = geom->GetNModules() * geom->GetNCristalsInModule() ;
  
  fFirstEvent = 0 ; 
  fLastEvent = fFirstEvent ; 
  if (fDigInput) 
    fInput = fDigInput->GetNinputs() ; 
  else 
    fInput           = 1 ;
  
  fInputFileNames  = new TString[fInput] ; 
  fEventNames      = new TString[fInput] ; 
  fInputFileNames[0] = GetTitle() ; 
  fEventNames[0]     = fEventFolderName.Data() ; 
  Int_t index ; 
  for (index = 1 ; index < fInput ; index++) {
    fInputFileNames[index] = static_cast<AliStream*>(fDigInput->GetInputStream(index))->GetFileName(0); 
    TString tempo = fDigInput->GetInputFolderName(index) ;
    fEventNames[index] = tempo.Remove(tempo.Length()-1) ; // strip of the stream number added by fDigInput
  }

  //to prevent cleaning of this object while GetEvent is called
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  if(!rl){
    rl = AliRunLoader::Open(GetTitle(), fEventFolderName) ; 
  }
  return fInit ; 
}

//____________________________________________________________________________ 
void AliPHOSDigitizer::InitParameters()
{
  // Set initial parameters Digitizer

  fDigitsInRun  = 0 ; 
  SetEventRange(0,-1) ;
  fPulse = new AliPHOSPulseGenerator();
  fADCValuesLG = new Int_t[fPulse->GetRawFormatTimeBins()];
  fADCValuesHG = new Int_t[fPulse->GetRawFormatTimeBins()];
    
}

//__________________________________________________________________
void AliPHOSDigitizer::Print(const Option_t *)const 
{
  // Print Digitizer's parameters
  AliInfo(Form("\n------------------- %s -------------", GetName() )) ; 
  if( strcmp(fEventFolderName.Data(), "") != 0 ){
    printf(" Writing Digits to branch with title  %s\n", fEventFolderName.Data()) ;
    
    Int_t nStreams ; 
    if (fDigInput) 
      nStreams =  GetNInputStreams() ;
    else 
      nStreams = fInput ; 
    
    Int_t index = 0 ;  
    for (index = 0 ; index < nStreams ; index++) {  
      TString tempo(fEventNames[index]) ; 
      tempo += index ;
      printf ("Adding SDigits from %s \n", fInputFileNames[index].Data()) ; 
    }
    
 //   printf("\nWith following parameters:\n") ;
 //   printf("  Electronics noise in EMC (fPinNoise)    = %f GeV\n", fPinNoise ) ; 
 //   printf("  Threshold  in EMC  (fEMCDigitThreshold) = %f GeV\n", fEMCDigitThreshold ) ;  
 //   printf("  Noise in CPV (fCPVNoise)                = %f aux units\n", fCPVNoise ) ; 
 //   printf("  Threshold in CPV (fCPVDigitThreshold)   = %f aux units\n",fCPVDigitThreshold ) ;  
    printf(" ---------------------------------------------------\n") ;   
  }
  else
    AliInfo(Form("AliPHOSDigitizer not initialized" )) ;
  
}

//__________________________________________________________________
 void AliPHOSDigitizer::PrintDigits(Option_t * option)
{
  // Print a table of digits
  

  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
  TClonesArray * digits = phosLoader->Digits() ; 
  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;
  
  AliInfo(Form("%d", digits->GetEntriesFast())) ; 
  printf("\nevent %d", gAlice->GetEvNumber()) ;
  printf("\n       Number of entries in Digits list %d", digits->GetEntriesFast() )  ;  

 
  if(strstr(option,"all")||strstr(option,"EMC")){  
    //loop over digits
    AliPHOSDigit * digit;
    printf("\nEMC digits (with primaries):\n")  ;
    printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
    Int_t maxEmc = geom->GetNModules()*geom->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; (index < digits->GetEntriesFast()) && 
	   (static_cast<AliPHOSDigit *>(digits->At(index))->GetId() <= maxEmc) ; index++) {
      digit = (AliPHOSDigit * )  digits->At(index) ;
      if(digit->GetNprimary() == 0) 
	continue;
//       printf("%6d  %8d    %6.5e %4d      %2d :",
// 	      digit->GetId(), digit->GetAmp(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  // YVK
      printf("%6d  %.4f    %6.5e %4d      %2d :",
	      digit->GetId(), digit->GetEnergy(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	printf("%d ",digit->GetPrimary(iprimary+1) ) ; 
      }    
      printf("\n") ;
    }
  }
  
  if(strstr(option,"all")||strstr(option,"CPV")){
    
    //loop over CPV digits
    AliPHOSDigit * digit;
    printf("\nCPV digits (with primaries):\n")  ;
    printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
    Int_t maxEmc = geom->GetNModules()*geom->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; index < digits->GetEntriesFast(); index++) {
      digit = (AliPHOSDigit * )  digits->At(index) ;
      if(digit->GetId() > maxEmc){
	printf("%6d  %8d    %4d      %2d :",
		digit->GetId(), digit->GetAmp(), digit->GetIndexInList(), digit->GetNprimary()) ;  
	Int_t iprimary;
	for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	  printf("%d ",digit->GetPrimary(iprimary+1) ) ; 
	}    
	printf("\n") ;
      }
    }
  }

}

//__________________________________________________________________
Float_t AliPHOSDigitizer::TimeOfNoise(void) const
{  // Calculates the time signal generated by noise
  //PH  Info("TimeOfNoise", "Change me") ; 
  return gRandom->Rndm() * 1.28E-5;
}

//__________________________________________________________________
void AliPHOSDigitizer::Unload() 
{  
  
  Int_t i ; 
  for(i = 1 ; i < fInput ; i++){
    TString tempo(fEventNames[i]) ; 
    tempo += i ;
    AliRunLoader* rl = AliRunLoader::GetRunLoader(tempo) ; 
    AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
    phosLoader->UnloadSDigits() ; 
  }
  
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
  phosLoader->UnloadDigits() ; 
}

//____________________________________________________________________________
void AliPHOSDigitizer::WriteDigits()
{

  // Makes TreeD in the output file. 
  // Check if branch already exists: 
  //   if yes, exit without writing: ROOT TTree does not support overwriting/updating of 
  //      already existing branches. 
  //   else creates branch with Digits, named "PHOS", title "...",
  //      and branch "AliPHOSDigitizer", with the same title to keep all the parameters
  //      and names of files, from which digits are made.

  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
 
  const TClonesArray * digits = phosLoader->Digits() ; 
  TTree * treeD = phosLoader->TreeD();
  if(!treeD){
    phosLoader->MakeTree("D");
    treeD = phosLoader->TreeD();
  }
  
  // -- create Digits branch
  Int_t bufferSize = 32000 ;    
  TBranch * digitsBranch = treeD->Branch("PHOS","TClonesArray",&digits,bufferSize);
  digitsBranch->SetTitle(fEventFolderName);
  digitsBranch->Fill() ;
  
  phosLoader->WriteDigits("OVERWRITE");

  Unload() ; 

}
