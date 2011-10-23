/*************************************************************************
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

//_________________________________________________________________________
// This is a Class that makes SDigits out of Hits
// A Summable Digits is the sum of all hits originating 
// from one in one tower of the EMCAL 
// A threshold for assignment of the primary to SDigit is applied 
//
// JLK 26-Jun-2008 Added explanation:
// SDigits need to hold the energy sum of the hits, but AliEMCALDigit
// can (should) only store amplitude.  Therefore, the SDigit energy is
// "digitized" before being stored and must be "calibrated" back to an
// energy before SDigits are summed to form true Digits
//
// SDigits are written to TreeS, branch "EMCAL"
// AliEMCALSDigitizer with all current parameters is written 
// to TreeS branch "AliEMCALSDigitizer".
// Both branches have the same title. If necessary one can produce 
// another set of SDigits with different parameters. Two versions
// can be distunguished using titles of the branches.
// User case:
// root [0] AliEMCALSDigitizer * s = new AliEMCALSDigitizer("galice.root")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root [1] s->Digitize()
//             // Makes SDigitis for all events stored in galice.root
// root [2] s->SetPedestalParameter(0.001)
//             // One can change parameters of digitization
// root [3] s->SetSDigitsBranch("Redestal 0.001")
//             // and write them into the new branch
// root [4] s->Digitize("deb all tim")
//             // available parameters:
//             deb - print # of produced SDigitis
//             deb all  - print # and list of produced SDigits
//             tim - print benchmarking information
//
//*-- Author : Sahal Yacoob (LBL)
// based on  : AliPHOSSDigitzer 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TBenchmark.h>
#include <TBrowser.h>
//#include <Riostream.h>
#include <TMath.h>
#include <TROOT.h>

// --- Standard library ---
#include "stdlib.h"

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliEMCALDigit.h"
#include "AliEMCALLoader.h"
#include "AliEMCALHit.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALSimParam.h"

ClassImp(AliEMCALSDigitizer)
           
//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer()
  : TNamed("",""),
    fA(0.),fB(0.),fECPrimThreshold(0.),
    fDefaultInit(kTRUE),
    fEventFolderName(0),
    fInit(0),
    fSDigitsInRun(0),
    fFirstEvent(0),
    fLastEvent(0),
    fSampling(0.),
	fHits(0)
{
  // ctor
  InitParameters();
}

//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const char * alirunFileName, 
				       const char * eventFolderName)
  : TNamed("EMCALSDigitizer", alirunFileName),
    fA(0.),fB(0.),fECPrimThreshold(0.),
    fDefaultInit(kFALSE),
    fEventFolderName(eventFolderName),
    fInit(0),
    fSDigitsInRun(0),
    fFirstEvent(0),
    fLastEvent(0),
    fSampling(0.),
    fHits(0)
{
  // ctor
  Init();
  InitParameters() ; 

}


//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) 
  : TNamed(sd.GetName(),sd.GetTitle()),
    fA(sd.fA),
    fB(sd.fB),
    fECPrimThreshold(sd.fECPrimThreshold),
    fDefaultInit(sd.fDefaultInit),
    fEventFolderName(sd.fEventFolderName),
    fInit(sd.fInit),
    fSDigitsInRun(sd.fSDigitsInRun),
    fFirstEvent(sd.fFirstEvent),
    fLastEvent(sd.fLastEvent),
    fSampling(sd.fSampling),
    fHits(0)
{
  //cpy ctor 
}


//____________________________________________________________________________ 
AliEMCALSDigitizer::~AliEMCALSDigitizer() {
  //dtor
	
  if(fHits){
	  fHits->Clear();
	  delete fHits;
  }
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Init(){
  // Initialization: open root-file, allocate arrays for hits and sdigits
  // 
  // Initialization can not be done in the default constructor
  //============================================================= YS
  //  The initialisation is now done by the getter

  fInit = kTRUE; 
   
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));

  if ( emcalLoader == 0 ) {
    Fatal("Init", "Could not obtain the AliEMCALLoader");
    return ;
  } 
  
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::InitParameters()
{
  //initialize parameters for sdigitization

  const AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  if (geom->GetSampling() == 0.) {
    Fatal("InitParameters", "Sampling factor not set !") ; 
  }

  // Get the parameters from the OCDB via the loader
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  AliEMCALSimParam * simParam = 0x0;
  if(emcalLoader) simParam = emcalLoader->SimulationParameters();
	
  if(!simParam){
	simParam = AliEMCALSimParam::GetInstance();
	AliWarning("Simulation Parameters not available in OCDB?");
  }
	
  //
  //JLK 26-Jun-2008 THIS SHOULD HAVE BEEN EXPLAINED AGES AGO:
  //
  //In order to be able to store SDigit Energy info into
  //AliEMCALDigit, we need to convert it temporarily to an ADC amplitude
  //and later when summing SDigits to form digits, convert it back to
  //energy.  These fA and fB parameters accomplish this through the
  //Digitize() and Calibrate() methods
  //
  // Initializes parameters
  fA         = simParam->GetA(); //0;
  fB         = simParam->GetB(); //1.e+6;  // Changed 24 Apr 2007. Dynamic range now 2 TeV
  fSampling  = geom->GetSampling();

  // threshold for deposit energy of hit
  fECPrimThreshold  = simParam->GetECPrimaryThreshold();//0.05;// GeV // 22-may-07 was 0// 24-nov-04 - was 1.e-6;
  
  AliDebug(2,Form("Print: \n------------------- %s -------------\n",GetName()));
  AliDebug(2,Form("   fInit                                 %i\n", int(fInit)));
  AliDebug(2,Form("   fFirstEvent                           %i\n", fFirstEvent));
  AliDebug(2,Form("   fLastEvent                            %i\n", fLastEvent));
  AliDebug(2,Form("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()));
  AliDebug(2,Form("   with digitization parameters       A = %f\n", fA));
  AliDebug(2,Form("                                      B = %f\n", fB));
  AliDebug(2,Form("   Threshold for EC Primary assignment  = %f\n", fECPrimThreshold));
  AliDebug(2,Form("   Sampling                             = %f\n", fSampling));
  AliDebug(2,Form("---------------------------------------------------\n"));

}

//____________________________________________________________________________
void AliEMCALSDigitizer::Digitize(Option_t *option) 
{ 
	// Collects all hit of the same tower into digits
	TString o(option); o.ToUpper();
	if (strstr(option, "print") ) {
		
		AliDebug(2,Form("Print: \n------------------- %s -------------\n",GetName()));
		AliDebug(2,Form("   fInit                                 %i\n", int(fInit)));
		AliDebug(2,Form("   fFirstEvent                           %i\n", fFirstEvent));
		AliDebug(2,Form("   fLastEvent                            %i\n", fLastEvent));
		AliDebug(2,Form("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()));
		AliDebug(2,Form("   with digitization parameters       A = %f\n", fA));
		AliDebug(2,Form("                                      B = %f\n", fB));
		AliDebug(2,Form("   Threshold for EC Primary assignment  = %f\n", fECPrimThreshold));
		AliDebug(2,Form("   Sampling                             = %f\n", fSampling));
		AliDebug(2,Form("---------------------------------------------------\n"));
		
		return ; 
	}
	
	
	if(strstr(option,"tim"))
		gBenchmark->Start("EMCALSDigitizer");
	
	AliRunLoader *rl = AliRunLoader::Instance();
	AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
	
	if (!fInit) { // to prevent overwrite existing file
	  AliError( Form("Give a version name different from %s", fEventFolderName.Data()) ) ;
	  return ;
	}
	
	if (fLastEvent == -1) 
		fLastEvent = rl->GetNumberOfEvents() - 1 ;
	else {
		fLastEvent = TMath::Min(fLastEvent, rl->GetNumberOfEvents()-1);
	}
	Int_t nEvents   = fLastEvent - fFirstEvent + 1;
	
	Float_t energy=0.; // de * fSampling - 23-nov-04
	rl->LoadKinematics();
	rl->LoadHits("EMCAL");
  
  Int_t ievent;
	for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
		rl->GetEvent(ievent);
		TTree * treeS = emcalLoader->TreeS();
		if ( !treeS ) { 
			emcalLoader->MakeSDigitsContainer();
			treeS = emcalLoader->TreeS();
		}
    
		TClonesArray * sdigits = emcalLoader->SDigits() ;
		sdigits->Clear("C");
		
		Int_t nSdigits = 0 ;
		Int_t iHit, iTrack, iSDigit;
		
		AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(); 
		
		TTree *treeH = emcalLoader->TreeH();	
		if (treeH) {
			Int_t nTrack = treeH->GetEntries();  // TreeH has array of hits for every primary
			TBranch * branchH = treeH->GetBranch("EMCAL");
			//if(fHits)fHits->Clear();
			branchH->SetAddress(&fHits);
			for (iTrack = 0; iTrack < nTrack; iTrack++) {
				branchH->GetEntry(iTrack); 
        
        if(fHits){
          
          Int_t nHit = fHits->GetEntriesFast();
          for(iHit = 0; iHit< nHit;iHit++){
            
            AliEMCALHit * hit = dynamic_cast<AliEMCALHit*>(fHits->At(iHit)) ;
            AliEMCALDigit * curSDigit = 0 ;
            AliEMCALDigit * sdigit = 0 ;
            Bool_t newsdigit = kTRUE; 
            
            // hit->GetId() - Absolute Id number EMCAL segment
            if(hit){
              if(geom->CheckAbsCellId(hit->GetId())) { // was IsInECA(hit->GetId())
                energy = hit->GetEnergy() * fSampling; // 23-nov-04
                if(energy >  fECPrimThreshold )
                  // Assign primary number only if deposited energy is significant
                  curSDigit =  new AliEMCALDigit(hit->GetPrimary(),
                                                 hit->GetIparent(), hit->GetId(), 
                                                 Digitize(energy), hit->GetTime(),kFALSE,
                                                 -1, 0,0,energy ) ;
                else
                  curSDigit =  new AliEMCALDigit(-1, 
                                                 -1,
                                                 hit->GetId(), 
                                                 Digitize(energy), hit->GetTime(),kFALSE,
                                                 -1, 0,0,energy ) ;
              } else {
                Warning("Digitize"," abs id %i is bad \n", hit->GetId());
                newsdigit = kFALSE;
                curSDigit = 0;
              }
              
              if(curSDigit != 0){
                for(Int_t check= 0; check < nSdigits ; check++) {
                  sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(check)) ;
                  if(sdigit){
                    if( sdigit->GetId() == curSDigit->GetId()) { // Are we in the same ECAL tower ?              
                      *sdigit = *sdigit + *curSDigit;
                      newsdigit = kFALSE;
                    }
                  }// sdigit exists
                  else {
                    AliWarning("Sdigit do not exist");
                    newsdigit = kFALSE;  
                  }// sdigit does not exist
                }//sdigit loop
              }// currsdigit exists
              
              if (newsdigit) {
                new((*sdigits)[nSdigits])  AliEMCALDigit(*curSDigit);
                nSdigits++ ;  
              }
              delete curSDigit ;
              
            }// hit exists
            else AliFatal("Hit is NULL!");
            
          }  // loop over all hits (hit = deposited energy/entering particle)
          
        }//fHits is not NULL
        else AliFatal("fHit is NULL!");
        
				sdigits->Sort() ;
				
				nSdigits = sdigits->GetEntriesFast() ;
				fSDigitsInRun += nSdigits ;  
				
				for (iSDigit = 0 ; iSDigit < sdigits->GetEntriesFast() ; iSDigit++) { 
					AliEMCALDigit * sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(iSDigit)) ;
					if(sdigit)sdigit->SetIndexInList(iSDigit) ;
          else AliFatal("sdigit is NULL!");
				}	
				if(fHits)fHits->Clear();
			}//track loop
		}// tree exists
		
		// Now write SDigits    
		
		Int_t bufferSize = 32000 ;    
		TBranch * sdigitsBranch = treeS->GetBranch("EMCAL");
		if (sdigitsBranch)
			sdigitsBranch->SetAddress(&sdigits);
		else
			treeS->Branch("EMCAL",&sdigits,bufferSize);
		
		treeS->Fill();
		
		emcalLoader->WriteSDigits("OVERWRITE");
		
		//NEXT - SDigitizer
		//emcalLoader->WriteSDigitizer("OVERWRITE");  // why in event cycle ?
		
		if(strstr(option,"deb"))
			PrintSDigits(option) ;  
	}
	
	Unload();
		
	if(strstr(option,"tim")){
		gBenchmark->Stop("EMCALSDigitizer"); 
		printf("\n Digitize: took %f seconds for SDigitizing %f seconds per event\n", 
           gBenchmark->GetCpuTime("EMCALSDigitizer"), gBenchmark->GetCpuTime("EMCALSDigitizer")/nEvents ) ; 
	}
}

//__________________________________________________________________
Float_t AliEMCALSDigitizer::Digitize(Float_t energy)const {
  // Digitize the energy
  //
  //JLK 26-Jun-2008 EXPLANATION LONG OVERDUE:
  //
  //We have to digitize the SDigit energy so that it can be stored in
  //AliEMCALDigit, which has only an ADC amplitude member and
  //(rightly) no energy member.  This method converts the energy to an
  //integer which can be re-converted back to an energy with the
  //Calibrate(energy) method when it is time to create Digits from SDigits 
  //
  Double_t aSignal = fA + energy*fB;
  if (TMath::Abs(aSignal)>2147483647.0) { 
    //PH 2147483647 is the max. integer
    //PH This apparently is a problem which needs investigation
    AliWarning(Form("Too big or too small energy %f",aSignal));
    aSignal = TMath::Sign((Double_t)2147483647,aSignal);
  }

  return (Float_t ) aSignal;
}

//__________________________________________________________________
Float_t AliEMCALSDigitizer::Calibrate(Float_t amp)const {
  //
  // Convert the amplitude back to energy in GeV
  //
  //JLK 26-Jun-2008 EXPLANATION LONG OVERDUE:
  //
  //We have to digitize the SDigit energy with the method Digitize() 
  //so that it can be stored in AliEMCALDigit, which has only an ADC 
  //amplitude member and (rightly) no energy member.  This method is
  //just the reverse of Digitize(): it converts the stored amplitude 
  //back to an energy value in GeV so that the SDigit energies can be 
  //summed before adding noise and creating digits out of them
  //
  return (Float_t)(amp - fA)/fB;

}
 

//__________________________________________________________________
void AliEMCALSDigitizer::Print1(Option_t * option)
{
  Print(); 
  PrintSDigits(option);
}

//__________________________________________________________________
void AliEMCALSDigitizer::Print(Option_t *option) const
{ 
  // Prints parameters of SDigitizer
  printf("Print: \n------------------- %s ------------- option %s\n", GetName() , option) ; 
  printf("   fInit                                 %i\n", int(fInit));
  printf("   fFirstEvent                           %i\n", fFirstEvent);
  printf("   fLastEvent                            %i\n", fLastEvent);
  printf("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()) ;
  printf("   with digitization parameters       A = %f\n", fA) ; 
  printf("                                      B = %f\n", fB) ;
  printf("   Threshold for EC Primary assignment  = %f\n", fECPrimThreshold)  ;
  printf("   Sampling                             = %f\n", fSampling);
  printf("---------------------------------------------------\n") ;
}

//__________________________________________________________________
Bool_t AliEMCALSDigitizer::operator==( AliEMCALSDigitizer const &sd )const
{
  // Equal operator.
  // SDititizers are equal if their pedestal, slope and threshold are equal
  if( (fA==sd.fA)&&(fB==sd.fB)&&
      (fECPrimThreshold==sd.fECPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}

//__________________________________________________________________
void AliEMCALSDigitizer::PrintSDigits(Option_t * option)
{
  //Prints list of digits produced at the current pass of AliEMCALDigitizer
  
  AliEMCALLoader *rl = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  if(rl){
    const TClonesArray * sdigits = rl->SDigits() ; 
    
    printf("\n") ;  
    printf("event %i", rl->GetRunLoader()->GetEventNumber());
    printf(" Number of entries in SDigits list %i", sdigits->GetEntriesFast()); 
    if(strstr(option,"all")||strstr(option,"EMC")){
      
      //loop over digits
      AliEMCALDigit * digit;
      printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
      Int_t   index = 0;
      Float_t isum  = 0.;
      const Int_t bufferSize = 8192;
      char * tempo = new char[bufferSize]; 
      for (index = 0 ; index < sdigits->GetEntries() ; index++) {
        digit = dynamic_cast<AliEMCALDigit *>( sdigits->At(index) ) ;
        if(digit){
          snprintf(tempo, bufferSize,"\n%6d  %8f    %6.5e %4d      %2d :",
                  digit->GetId(), digit->GetAmplitude(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
          printf("%s",tempo);
          isum += digit->GetAmplitude();
          
          Int_t iprimary;
          for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
            snprintf(tempo,bufferSize, "%d ",digit->GetPrimary(iprimary+1) ) ; 
            printf("%s",tempo); 
          }  
        } //sdigit exists
        else AliFatal("SDigit is NULL!");
      }//loop
      delete [] tempo ;
      printf("\n** Sum %2.3f : %10.3f GeV/c **\n ", isum, Calibrate(isum));
    } else printf("\n");
  }
  else AliFatal("EMCALLoader is NULL!");
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Unload() const
{
  // Unload Hits and SDigits from the folder
  AliEMCALLoader *rl = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  if(rl){
    rl->UnloadHits() ; 
    rl->UnloadSDigits() ;
  }
  else AliFatal("EMCALLoader is NULL!");
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Browse(TBrowser* b)
{
  TNamed::Browse(b);
}
