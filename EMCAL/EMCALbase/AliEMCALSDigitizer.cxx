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

// --- ROOT system ---
#include <TBenchmark.h>
#include <TBrowser.h>
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
#include "AliSort.h"

/// \cond CLASSIMP
ClassImp(AliEMCALSDigitizer) ;
/// \endcond
    
///
/// Default Constructor
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
  InitParameters();
}

///
/// Constructor
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
  Init();
  InitParameters() ; 
}

///
/// Copy constructor
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
{ }

///
/// Assignment operator; use copy constructor
//_____________________________________________________________________
AliEMCALSDigitizer& AliEMCALSDigitizer::operator = (const AliEMCALSDigitizer& source)
{ 
  if (&source == this) return *this;

  new (this) AliEMCALSDigitizer(source);
  return *this;
}

///
/// Destructor
//____________________________________________________________________________ 
AliEMCALSDigitizer::~AliEMCALSDigitizer() 
{	
  if(fHits)
  {
	  fHits->Clear();
	  delete fHits;
  }
}

///
/// Initialization of loader
/// 
/// Initialization can not be done in the default constructor
//____________________________________________________________________________ 
void AliEMCALSDigitizer::Init()
{
  fInit = kTRUE; 
   
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));

  if ( emcalLoader == 0 ) 
  {
    AliFatal("Could not obtain the AliEMCALLoader");
    return ;
  } 
}

///
/// Initialize parameters for sdigitization
//____________________________________________________________________________ 
void AliEMCALSDigitizer::InitParameters()
{
  const AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  if (geom->GetSampling() == 0.) 
    AliFatal("Sampling factor not set!") ; 
  
  
  // Get the parameters from the OCDB via the loader
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  AliEMCALSimParam * simParam = 0x0;
  if(emcalLoader) simParam = emcalLoader->SimulationParameters();
  
  if(!simParam)
  {
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

///
/// Collect all hit of the same tower into a summable digit
//____________________________________________________________________________
void AliEMCALSDigitizer::Digitize(Option_t *option) 
{ 
  TString o(option); o.ToUpper();
  if (strstr(option, "print") ) 
  {
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
  
  if (!fInit) 
  { 
    // to prevent overwrite existing file
    AliError( Form("Give a version name different from %s", fEventFolderName.Data()) ) ;
    return ;
  }
  
  if (fLastEvent == -1) 
    fLastEvent = rl->GetNumberOfEvents() - 1 ;
  else 
    fLastEvent = TMath::Min(fLastEvent, rl->GetNumberOfEvents()-1);
  
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;
  
  Float_t energy=0.; // de * fSampling - 23-nov-04
  rl->LoadKinematics();
  rl->LoadHits("EMCAL");
  
  Int_t ievent;
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) 
  {
    rl->GetEvent(ievent);
    
    TTree * treeS = emcalLoader->TreeS();
    if ( !treeS ) 
    { 
      emcalLoader->MakeSDigitsContainer();
      treeS = emcalLoader->TreeS();
    }
    
    TClonesArray * sdigits = emcalLoader->SDigits() ;
    sdigits->Clear("C");
    
    Int_t nSdigits = 0 ;
    Int_t iHit, iTrack, iSDigit;
    
    AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(); 
    
    TTree *treeH = emcalLoader->TreeH();	
    if (treeH) 
    {
      Int_t nTrack = treeH->GetEntries();  // TreeH has array of hits for every primary
      TBranch * branchH = treeH->GetBranch("EMCAL");
      //if(fHits)fHits->Clear();
      branchH->SetAddress(&fHits);
      for (iTrack = 0; iTrack < nTrack; iTrack++) 
      {
        branchH->GetEntry(iTrack); 
        
        if(fHits)
        {
          Int_t nHit = fHits->GetEntriesFast();
          for(iHit = 0; iHit< nHit;iHit++)
          {
            AliEMCALHit * hit = static_cast<AliEMCALHit*>(fHits->UncheckedAt(iHit)) ;
            AliEMCALDigit * curSDigit = 0 ;
            AliEMCALDigit * sdigit = 0 ;
            Bool_t newsdigit = kTRUE; 
            
            // hit->GetId() - Absolute Id number EMCAL segment
            if(hit)
            {
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
              } 
              else
              {
                Warning("Digitize"," abs id %i is bad \n", hit->GetId());
                newsdigit = kFALSE;
                curSDigit = 0;
              }
              
              if(curSDigit != 0)
              {
                for(Int_t check= 0; check < nSdigits ; check++) 
                {
                  sdigit = static_cast<AliEMCALDigit *>(sdigits->UncheckedAt(check)) ;
                  
                  if(sdigit)
                  {
                    if( sdigit->GetId() == curSDigit->GetId()) 
                    { // Are we in the same ECAL tower ?              
                      *sdigit = *sdigit + *curSDigit;
                      newsdigit = kFALSE;
                    }
                  }// sdigit exists
                  else 
                  {
                    AliWarning("Sdigit do not exist");
                    newsdigit = kFALSE;  
                  }// sdigit does not exist
                }//sdigit loop
              }// currsdigit exists
              
              if (newsdigit) 
              {
                new((*sdigits)[nSdigits])  AliEMCALDigit(*curSDigit);
                nSdigits++ ;  
              }
              
              delete curSDigit ;
              
            }// hit exists
            else AliFatal("Hit is NULL!");
            
          }  // loop over all hits (hit = deposited energy/entering particle)
          
        }//fHits is not NULL
        else AliFatal("fHit is NULL!");
        
        // sdigits->Sort() ;
        AliSort::TClonesArraySort<AliEMCALDigit>(sdigits);
        
        nSdigits = sdigits->GetEntriesFast() ;
        fSDigitsInRun += nSdigits ;  
        
        for (iSDigit = 0 ; iSDigit < sdigits->GetEntriesFast() ; iSDigit++) 
        { 
          AliEMCALDigit * sdigit = static_cast<AliEMCALDigit *>(sdigits->UncheckedAt(iSDigit)) ;
          
          if  ( sdigit ) sdigit->SetIndexInList(iSDigit) ;
          else           AliFatal("sdigit is NULL!");
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
		
  if(strstr(option,"tim"))
  {
    gBenchmark->Stop("EMCALSDigitizer"); 
    printf("\n Digitize: took %f seconds for SDigitizing %f seconds per event\n", 
           gBenchmark->GetCpuTime("EMCALSDigitizer"), gBenchmark->GetCpuTime("EMCALSDigitizer")/nEvents ) ; 
  }
}

///
/// Digitize the energy
//
///JLK 26-Jun-2008 EXPLANATION LONG OVERDUE:
//
/// We have to digitize the SDigit energy so that it can be stored in
/// AliEMCALDigit, which has only an ADC amplitude member and
/// (rightly) no energy member.  This method converts the energy to an
/// integer which can be re-converted back to an energy with the
/// Calibrate(energy) method when it is time to create Digits from SDigits 
///
/// \param energy: summed energy in GeV
//__________________________________________________________________
Float_t AliEMCALSDigitizer::Digitize(Float_t energy) const 
{
  Double_t aSignal = fA + energy*fB;
  if (TMath::Abs(aSignal)>2147483647.0) 
  { 
    //PH 2147483647 is the max. integer
    //PH This apparently is a problem which needs investigation
    AliWarning(Form("Too big or too small energy %f",aSignal));
    
    aSignal = TMath::Sign((Double_t)2147483647,aSignal);
  }
  
  return (Float_t ) aSignal;
}

///
/// Convert the amplitude back to energy in GeV
//
///JLK 26-Jun-2008 EXPLANATION LONG OVERDUE:
//
/// We have to digitize the SDigit energy with the method Digitize() 
/// so that it can be stored in AliEMCALDigit, which has only an ADC 
/// amplitude member and (rightly) no energy member.  This method is
/// just the reverse of Digitize(): it converts the stored amplitude 
/// back to an energy value in GeV so that the SDigit energies can be 
/// summed before adding noise and creating digits out of them
///
/// \param amp: energy in ADC of the digit
//__________________________________________________________________
Float_t AliEMCALSDigitizer::Calibrate(Float_t amp) const 
{
  return (Float_t)(amp - fA)/fB;
}

///
/// Print info, call all prints
//__________________________________________________________________
void AliEMCALSDigitizer::Print1(Option_t * option)
{
  Print(); 
  PrintSDigits(option);
}

///
/// Prints parameters of SDigitizer
//__________________________________________________________________
void AliEMCALSDigitizer::Print(Option_t *option) const
{ 
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

/// Equal operator.
/// SDititizers are equal if their pedestal, slope and threshold are equal
//__________________________________________________________________
Bool_t AliEMCALSDigitizer::operator==( AliEMCALSDigitizer const &sd )const
{
  if( (fA==sd.fA)&&(fB==sd.fB)&&
     (fECPrimThreshold==sd.fECPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}

///
/// Prints list of digits produced at the current pass of AliEMCALDigitizer
//__________________________________________________________________
void AliEMCALSDigitizer::PrintSDigits(Option_t * option)
{  
  AliEMCALLoader *rl = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  
  if(rl)
  {
    const TClonesArray * sdigits = rl->SDigits() ; 
    
    printf("\n") ;  
    printf("event %i", rl->GetRunLoader()->GetEventNumber());
    printf(" Number of entries in SDigits list %i", sdigits->GetEntriesFast()); 
    
    if(strstr(option,"all")||strstr(option,"EMC"))
    {
      // loop over digits
      AliEMCALDigit * digit;
      printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
      Int_t   index = 0;
      Float_t isum  = 0.;
      const Int_t bufferSize = 8192;
      char * tempo = new char[bufferSize]; 
      
      for (index = 0 ; index < sdigits->GetEntries() ; index++) 
      {
        digit = dynamic_cast<AliEMCALDigit *>( sdigits->At(index) ) ;
        if(digit)
        {
          snprintf(tempo, bufferSize,"\n%6d  %8f    %6.5e %4d      %2d :",
                   digit->GetId(), digit->GetAmplitude(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
          printf("%s",tempo);
          isum += digit->GetAmplitude();
          
          Int_t iprimary;
          for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) 
          {
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

///
/// Unload Hits and SDigits from the folder
//____________________________________________________________________________ 
void AliEMCALSDigitizer::Unload() const
{
  AliEMCALLoader *rl = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  
  if(rl)
  {
    rl->UnloadHits() ; 
    rl->UnloadSDigits() ;
  }
  else AliFatal("EMCALLoader is NULL!");
}

///
/// Browse (obsolete?)
//____________________________________________________________________________ 
void AliEMCALSDigitizer::Browse(TBrowser* b)
{
  TNamed::Browse(b);
}
