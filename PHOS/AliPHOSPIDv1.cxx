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
 * $Log$
 * Revision 1.103  2006/09/07 18:31:08  kharlov
 * Effective c++ corrections (T.Pocheptsov)
 *
 * Revision 1.102  2006/01/23 17:51:48  hristov
 * Using the recommended way of forward declarations for TVector and TMatrix (see v5-08-00 release notes). Additional clean-up
 *
 * Revision 1.101  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version v1 of the PHOS particle identifier 
// Particle identification based on the 
//     - RCPV: distance from CPV recpoint to EMCA recpoint.
//     - TOF 
//     - PCA: Principal Components Analysis..
// The identified particle has an identification number corresponding 
// to a 9 bits number:
//     -Bit 0 to 2: bit set if RCPV > CpvEmcDistance (each bit corresponds
//      to a different efficiency-purity point of the photon identification) 
//     -Bit 3 to 5: bit set if TOF  < TimeGate (each bit corresponds
//      to a different efficiency-purity point of the photon identification) 
//     -Bit 6 to 9: bit set if Principal Components are 
//      inside an ellipse defined by the parameters a, b, c, x0 and y0.
//      (each bit corresponds to a different efficiency-purity point of the 
//      photon identification)
//      The PCA (Principal components analysis) needs a file that contains
//      a previous analysis of the correlations between the particles. This 
//      file is $ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root. Analysis done for 
//      energies between 0.5 and 100 GeV.
//      A calibrated energy is calculated. The energy of the reconstructed
//      cluster is corrected with the formula A + B * E  + C * E^2, whose 
//      parameters where obtained through the study of the reconstructed 
//      energy distribution of monoenergetic photons. 
//
//      All the parameters (RCPV(2 rows-3 columns),TOF(1r-3c),PCA(5r-4c) 
//      and calibration(1r-3c))are stored in a file called 
//      $ALICE_ROOT/PHOS/Parameters.dat. Each time that AliPHOSPIDv1 is 
//      initialized, this parameters are copied to a Matrix (9,4), a 
//      TMatrixD object.  
//
// use case:
//  root [0] AliPHOSPIDv1 * p = new AliPHOSPIDv1("galice1.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//          // reading headers from file galice1.root and create  RecParticles 
            // TrackSegments and RecPoints are used 
//          // set file name for the branch RecParticles
//  root [1] p->ExecuteTask("deb all time")
//          // available options
//          // "deb" - prints # of reconstructed particles
//          // "deb all" -  prints # and list of RecParticles
//          // "time" - prints benchmarking results
//                  
//  root [2] AliPHOSPIDv1 * p2 = new AliPHOSPIDv1("galice1.root","v1",kTRUE)
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                //Split mode.  
//  root [3] p2->ExecuteTask()
//


//*-- Author: Yves Schutz (SUBATECH)  & Gines Martinez (SUBATECH) & 
//            Gustavo Conesa April 2002
//            PCA redesigned by Gustavo Conesa October 2002:
//            The way of using the PCA has changed. Instead of 2
//            files with the PCA, each one with different energy ranges 
//            of application, we use the wide one (0.5-100 GeV), and instead
//            of fixing 3 ellipses for different ranges of energy, it has been
//            studied the dependency of the ellipses parameters with the 
//            energy, and they are implemented in the code as a funtion 
//            of the energy. 
//
//
//
// --- ROOT system ---


// --- Standard library ---
#include <TMatrixF.h>
#include "TFormula.h"
#include "TBenchmark.h"
#include "TPrincipal.h"
#include "TFile.h" 
#include "TSystem.h"

// --- AliRoot header files ---
	      //#include "AliLog.h"
#include "AliGenerator.h"
#include "AliPHOS.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSGetter.h"

ClassImp( AliPHOSPIDv1) 

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1() :
  fBayesian(kFALSE),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNEvent(0),
  fFileNamePrincipalPhoton(),
  fFileNamePrincipalPi0(),
  fFileNameParameters(),
  fPrincipalPhoton(0),
  fPrincipalPi0(0),
  fX(0),
  fPPhoton(0),
  fPPi0(0),
  fRecParticlesInRun(0),
  fParameters(0),
  fTFphoton(0),
  fTFpiong(0),
  fTFkaong(0),
  fTFkaonl(0),
  fTFhhadrong(0),
  fTFhhadronl(0),
  fDFmuon(0),
  fERecWeight(0),
  fChargedNeutralThreshold(0.),
  fTOFEnThreshold(0),
  fDispEnThreshold(0),
  fDispMultThreshold(0)
{ 
  // default ctor
 
  InitParameters() ; 
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const AliPHOSPIDv1 & pid ) : 
  AliPHOSPID(pid),
  fBayesian(kFALSE),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNEvent(0),
  fFileNamePrincipalPhoton(),
  fFileNamePrincipalPi0(),
  fFileNameParameters(),
  fPrincipalPhoton(0),
  fPrincipalPi0(0),
  fX(0),
  fPPhoton(0),
  fPPi0(0),
  fRecParticlesInRun(0),
  fParameters(0),
  fTFphoton(0),
  fTFpiong(0),
  fTFkaong(0),
  fTFkaonl(0),
  fTFhhadrong(0),
  fTFhhadronl(0),
  fDFmuon(0),
  fERecWeight(0),
  fChargedNeutralThreshold(0.),
  fTOFEnThreshold(0),
  fDispEnThreshold(0),
  fDispMultThreshold(0)

{ 
  // ctor
  InitParameters() ; 
  Init() ;

}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const TString alirunFileName, const TString eventFolderName) :
  AliPHOSPID(alirunFileName, eventFolderName),
  fBayesian(kFALSE),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNEvent(0),
  fFileNamePrincipalPhoton(),
  fFileNamePrincipalPi0(),
  fFileNameParameters(),
  fPrincipalPhoton(0),
  fPrincipalPi0(0),
  fX(0),
  fPPhoton(0),
  fPPi0(0),
  fRecParticlesInRun(0),
  fParameters(0),
  fTFphoton(0),
  fTFpiong(0),
  fTFkaong(0),
  fTFkaonl(0),
  fTFhhadrong(0),
  fTFhhadronl(0),
  fDFmuon(0),
  fERecWeight(0),
  fChargedNeutralThreshold(0.),
  fTOFEnThreshold(0),
  fDispEnThreshold(0),
  fDispMultThreshold(0)

{ 
  //ctor with the indication on where to look for the track segments
 
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________
AliPHOSPIDv1::~AliPHOSPIDv1()
{ 
  // dtor
  fPrincipalPhoton = 0;
  fPrincipalPi0 = 0;

  delete [] fX ;       // Principal input 
  delete [] fPPhoton ; // Photon Principal components
  delete [] fPPi0 ;    // Pi0 Principal components

  delete fParameters;
  delete fTFphoton;
  delete fTFpiong;
  delete fTFkaong;
  delete fTFkaonl;
  delete fTFhhadrong;
  delete fTFhhadronl;
  delete fDFmuon;
}
//____________________________________________________________________________
const TString AliPHOSPIDv1::BranchName() const 
{  

  return GetName() ;
}
 
//____________________________________________________________________________
void AliPHOSPIDv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of PHOS tasks

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  if(!gime)
    gime = AliPHOSGetter::Instance(GetTitle(), fEventFolderName.Data()) ; 

  if ( !gime->PID() ) 
    gime->PostPID(this) ;
}

//____________________________________________________________________________
void AliPHOSPIDv1::InitParameters()
{
  // Initialize PID parameters
  fWrite                   = kTRUE ;
  fRecParticlesInRun = 0 ; 
  fNEvent            = 0 ;            
  fRecParticlesInRun = 0 ;
  fBayesian          = kTRUE ;
  SetParameters() ; // fill the parameters matrix from parameters file
  SetEventRange(0,-1) ;

  // initialisation of response function parameters
  // Tof

//   // Photons
//   fTphoton[0] = 0.218    ;
//   fTphoton[1] = 1.55E-8  ; 
//   fTphoton[2] = 5.05E-10 ;
//   fTFphoton = new TFormula("ToF response to photons" , "gaus") ; 
//   fTFphoton->SetParameters( fTphoton[0], fTphoton[1], fTphoton[2]) ; 

//   // Pions
//   //Gaus (0 to max probability)
//   fTpiong[0] = 0.0971    ; 
//   fTpiong[1] = 1.58E-8  ; 
//   fTpiong[2] = 5.69E-10 ;
//   fTFpiong = new TFormula("ToF response to pions" , "gaus") ; 
//   fTFpiong->SetParameters( fTpiong[0], fTpiong[1], fTpiong[2]) ; 

//   // Kaons
//   //Gaus (0 to max probability)
//   fTkaong[0] = 0.0542  ; 
//   fTkaong[1] = 1.64E-8 ; 
//   fTkaong[2] = 6.07E-10 ;
//   fTFkaong = new TFormula("ToF response to kaon" , "gaus") ; 
//   fTFkaong->SetParameters( fTkaong[0], fTkaong[1], fTkaong[2]) ; 
//   //Landau (max probability to inf) 
//   fTkaonl[0] = 0.264   ;
//   fTkaonl[1] = 1.68E-8  ; 
//   fTkaonl[2] = 4.10E-10 ;
//   fTFkaonl = new TFormula("ToF response to kaon" , "landau") ; 
//   fTFkaonl->SetParameters( fTkaonl[0], fTkaonl[1], fTkaonl[2]) ; 

//   //Heavy Hadrons
//   //Gaus (0 to max probability)
//   fThhadrong[0] = 0.0302   ;  
//   fThhadrong[1] = 1.73E-8  ; 
//   fThhadrong[2] = 9.52E-10 ;
//   fTFhhadrong = new TFormula("ToF response to heavy hadrons" , "gaus") ; 
//   fTFhhadrong->SetParameters( fThhadrong[0], fThhadrong[1], fThhadrong[2]) ; 
//   //Landau (max probability to inf) 
//   fThhadronl[0] = 0.139    ;  
//   fThhadronl[1] = 1.745E-8  ; 
//   fThhadronl[2] = 1.00E-9  ;
//   fTFhhadronl = new TFormula("ToF response to heavy hadrons" , "landau") ; 
//   fTFhhadronl->SetParameters( fThhadronl[0], fThhadronl[1], fThhadronl[2]) ; 

  // Photons
  fTphoton[0] = 7.83E8   ;
  fTphoton[1] = 1.55E-8  ; 
  fTphoton[2] = 5.09E-10 ;
  fTFphoton = new TFormula("ToF response to photons" , "gaus") ; 
  fTFphoton->SetParameters( fTphoton[0], fTphoton[1], fTphoton[2]) ; 

  // Pions
  //Gaus (0 to max probability)
  fTpiong[0] = 6.73E8    ; 
  fTpiong[1] = 1.58E-8  ; 
  fTpiong[2] = 5.87E-10 ;
  fTFpiong = new TFormula("ToF response to pions" , "gaus") ; 
  fTFpiong->SetParameters( fTpiong[0], fTpiong[1], fTpiong[2]) ; 

  // Kaons
  //Gaus (0 to max probability)
  fTkaong[0] = 3.93E8  ; 
  fTkaong[1] = 1.64E-8 ; 
  fTkaong[2] = 6.07E-10 ;
  fTFkaong = new TFormula("ToF response to kaon" , "gaus") ; 
  fTFkaong->SetParameters( fTkaong[0], fTkaong[1], fTkaong[2]) ; 
  //Landau (max probability to inf) 
  fTkaonl[0] = 2.0E9    ;
  fTkaonl[1] = 1.68E-8  ; 
  fTkaonl[2] = 4.10E-10 ;
  fTFkaonl = new TFormula("ToF response to kaon" , "landau") ; 
  fTFkaonl->SetParameters( fTkaonl[0], fTkaonl[1], fTkaonl[2]) ; 

  //Heavy Hadrons
  //Gaus (0 to max probability)
  fThhadrong[0] = 2.02E8   ;  
  fThhadrong[1] = 1.73E-8  ; 
  fThhadrong[2] = 9.52E-10 ;
  fTFhhadrong = new TFormula("ToF response to heavy hadrons" , "gaus") ; 
  fTFhhadrong->SetParameters( fThhadrong[0], fThhadrong[1], fThhadrong[2]) ; 
  //Landau (max probability to inf) 
  fThhadronl[0] = 1.10E9    ;  
  fThhadronl[1] = 1.74E-8   ; 
  fThhadronl[2] = 1.00E-9   ;
  fTFhhadronl = new TFormula("ToF response to heavy hadrons" , "landau") ; 
  fTFhhadronl->SetParameters( fThhadronl[0], fThhadronl[1], fThhadronl[2]) ; 



  // Shower shape: dispersion gaussian parameters
  // Photons
  
//   fDphoton[0] = 4.62e-2;  fDphoton[1] = 1.39e-2 ; fDphoton[2] = -3.80e-2;//constant
//   fDphoton[3] = 1.53   ;  fDphoton[4] =-6.62e-2 ; fDphoton[5] = 0.339   ;//mean
//   fDphoton[6] = 6.89e-2;  fDphoton[7] =-6.59e-2 ; fDphoton[8] = 0.194   ;//sigma
  
//   fDpi0[0] = 0.0586  ;  fDpi0[1] = 1.06E-3 ; fDpi0[2] = 0.      ;//constant
//   fDpi0[3] = 2.67    ;  fDpi0[4] =-2.00E-2 ; fDpi0[5] = 9.37E-5 ;//mean
//   fDpi0[6] = 0.153   ;  fDpi0[7] = 9.34E-4 ; fDpi0[8] =-1.49E-5 ;//sigma
  
//   fDhadron[0] = 1.61E-2 ;  fDhadron[1] = 3.03E-3 ; fDhadron[2] = 1.01E-2 ;//constant
//   fDhadron[3] = 3.81    ;  fDhadron[4] = 0.232   ; fDhadron[5] =-1.25    ;//mean
//   fDhadron[6] = 0.897   ;  fDhadron[7] = 0.0987  ; fDhadron[8] =-0.534   ;//sigma
  
  fDphoton[0] = 1.5    ;  fDphoton[1] = 0.49    ; fDphoton[2] =-1.7E-2 ;//constant
  fDphoton[3] = 1.5    ;  fDphoton[4] = 4.0E-2  ; fDphoton[5] = 0.21   ;//mean
  fDphoton[6] = 4.8E-2 ;  fDphoton[7] =-0.12    ; fDphoton[8] = 0.27   ;//sigma
  fDphoton[9] = 16.; //for E>  fDphoton[9] parameters calculated at  fDphoton[9]

  fDpi0[0] = 0.25      ;  fDpi0[1] = 3.3E-2     ; fDpi0[2] =-1.0e-5    ;//constant
  fDpi0[3] = 1.50      ;  fDpi0[4] = 398.       ; fDpi0[5] = 12.       ;//mean
  fDpi0[6] =-7.0E-2    ;  fDpi0[7] =-524.       ; fDpi0[8] = 22.       ;//sigma
  fDpi0[9] = 110.; //for E>  fDpi0[9] parameters calculated at  fDpi0[9]

  fDhadron[0] = 6.5    ;  fDhadron[1] =-5.3     ; fDhadron[2] = 1.5    ;//constant
  fDhadron[3] = 3.8    ;  fDhadron[4] = 0.23    ; fDhadron[5] =-1.2    ;//mean
  fDhadron[6] = 0.88   ;  fDhadron[7] = 9.3E-2  ; fDhadron[8] =-0.51   ;//sigma
  fDhadron[9] = 2.; //for E>  fDhadron[9] parameters calculated at  fDhadron[9]

  fDmuon[0] = 0.0631 ;
  fDmuon[1] = 1.4    ; 
  fDmuon[2] = 0.0557 ;
  fDFmuon = new TFormula("Shower shape response to muons" , "landau") ; 
  fDFmuon->SetParameters( fDmuon[0], fDmuon[1], fDmuon[2]) ; 


  // x(CPV-EMC) distance gaussian parameters
  
//   fXelectron[0] = 8.06e-2 ;  fXelectron[1] = 1.00e-2; fXelectron[2] =-5.14e-2;//constant
//   fXelectron[3] = 0.202   ;  fXelectron[4] = 8.15e-3; fXelectron[5] = 4.55   ;//mean
//   fXelectron[6] = 0.334   ;  fXelectron[7] = 0.186  ; fXelectron[8] = 4.32e-2;//sigma
  
//   //charged hadrons gaus
//   fXcharged[0] = 6.43e-3 ;  fXcharged[1] =-4.19e-5; fXcharged[2] = 1.42e-3;//constant
//   fXcharged[3] = 2.75    ;  fXcharged[4] =-0.40   ; fXcharged[5] = 1.68   ;//mean
//   fXcharged[6] = 3.135   ;  fXcharged[7] =-9.41e-2; fXcharged[8] = 1.31e-2;//sigma
  
//   // z(CPV-EMC) distance gaussian parameters
  
//   fZelectron[0] = 8.22e-2 ;  fZelectron[1] = 5.11e-3; fZelectron[2] =-3.05e-2;//constant
//   fZelectron[3] = 3.09e-2 ;  fZelectron[4] = 5.87e-2; fZelectron[5] =-9.49e-2;//mean
//   fZelectron[6] = 0.263   ;  fZelectron[7] =-9.02e-3; fZelectron[8] = 0.151 ;//sigma
  
//   //charged hadrons gaus
  
//   fZcharged[0] = 1.00e-2 ;  fZcharged[1] = 2.82E-4 ; fZcharged[2] = 2.87E-3 ;//constant
//   fZcharged[3] =-4.68e-2 ;  fZcharged[4] =-9.21e-3 ; fZcharged[5] = 4.91e-2 ;//mean
//   fZcharged[6] = 1.425   ;  fZcharged[7] =-5.90e-2 ; fZcharged[8] = 5.07e-2 ;//sigma


  fXelectron[0] =-1.6E-2 ;  fXelectron[1] = 0.77  ; fXelectron[2] =-0.15 ;//constant
  fXelectron[3] = 0.35   ;  fXelectron[4] = 0.25  ; fXelectron[5] = 4.12 ;//mean
  fXelectron[6] = 0.30   ;  fXelectron[7] = 0.11  ; fXelectron[8] = 0.16 ;//sigma
  fXelectron[9] = 3.; //for E>  fXelectron[9] parameters calculated at  fXelectron[9]

  //charged hadrons gaus
  fXcharged[0] = 0.14    ;  fXcharged[1] =-3.0E-2 ; fXcharged[2] = 0     ;//constant
  fXcharged[3] = 1.4     ;  fXcharged[4] =-9.3E-2 ; fXcharged[5] = 1.4   ;//mean
  fXcharged[6] = 5.7     ;  fXcharged[7] = 0.27   ; fXcharged[8] =-1.8   ;//sigma
  fXcharged[9] = 1.2; //for E>  fXcharged[9] parameters calculated at  fXcharged[9]

  // z(CPV-EMC) distance gaussian parameters
  
  fZelectron[0] = 0.49   ;  fZelectron[1] = 0.53   ; fZelectron[2] =-9.8E-2 ;//constant
  fZelectron[3] = 2.8E-2 ;  fZelectron[4] = 5.0E-2 ; fZelectron[5] =-8.2E-2 ;//mean
  fZelectron[6] = 0.25   ;  fZelectron[7] =-1.7E-2 ; fZelectron[8] = 0.17   ;//sigma
  fZelectron[9] = 3.; //for E>  fZelectron[9] parameters calculated at  fZelectron[9]

  //charged hadrons gaus
  
  fZcharged[0] = 0.46    ;  fZcharged[1] =-0.65    ; fZcharged[2] = 0.52    ;//constant
  fZcharged[3] = 1.1E-2  ;  fZcharged[4] = 0.      ; fZcharged[5] = 0.      ;//mean
  fZcharged[6] = 0.60    ;  fZcharged[7] =-8.2E-2  ; fZcharged[8] = 0.45    ;//sigma
  fZcharged[9] = 1.2; //for E>  fXcharged[9] parameters calculated at  fXcharged[9]

  //Threshold to differentiate between charged and neutral
  fChargedNeutralThreshold = 1e-5;
  fTOFEnThreshold          = 2;          //Maximum energy to use TOF
  fDispEnThreshold         = 0.5;       //Minimum energy to use shower shape
  fDispMultThreshold       = 3;       //Minimum multiplicity to use shower shape

  //Weight to hadrons recontructed energy

  fERecWeightPar[0] = 0.32 ; 
  fERecWeightPar[1] = 3.8  ;
  fERecWeightPar[2] = 5.4E-3 ; 
  fERecWeightPar[3] = 5.6E-2 ;
  fERecWeight = new TFormula("Weight for hadrons" , "[0]*exp(-x*[1])+[2]*exp(-x*[3])") ; 
  fERecWeight ->SetParameters(fERecWeightPar[0],fERecWeightPar[1] ,fERecWeightPar[2] ,fERecWeightPar[3]) ; 


  for (Int_t i =0; i<  AliPID::kSPECIESN ; i++)
    fInitPID[i] = 1.;
 
}

//________________________________________________________________________
void  AliPHOSPIDv1::Exec(Option_t *option)
{
  // Steering method to perform particle reconstruction and identification
  // for the event range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1 (by default), then process events until the end.
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSPID");
  
  if(strstr(option,"print")) {
    Print() ; 
    return ; 
  }


  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
 
  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1 ;
  else 
    fLastEvent = TMath::Min(fLastEvent,gime->MaxEvent());
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ; 
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    gime->Event(ievent,"TR") ;
    if(gime->TrackSegments() && //Skip events, where no track segments made
       gime->TrackSegments()->GetEntriesFast()) {

      MakeRecParticles() ;

      if(fBayesian)
	MakePID() ; 
      
      WriteRecParticles();
      if(strstr(option,"deb"))
	PrintRecParticles(option) ;
      //increment the total number of rec particles per run 
      fRecParticlesInRun += gime->RecParticles()->GetEntriesFast() ; 
    }
  }
  if(strstr(option,"deb"))
      PrintRecParticles(option);
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSPID");
    AliInfo(Form("took %f seconds for PID %f seconds per event", 
	 gBenchmark->GetCpuTime("PHOSPID"),  
	 gBenchmark->GetCpuTime("PHOSPID")/nEvents)) ;
  }
  if(fWrite)
    Unload();
}

//________________________________________________________________________
Double_t  AliPHOSPIDv1::GausF(Double_t  x, Double_t  y, Double_t * par)
{
  //Given the energy x and the parameter y (tof, shower dispersion or cpv-emc distance), 
  //this method returns a density probability of this parameter, given by a gaussian 
  //function whose parameters depend with the energy  with a function: a/(x*x)+b/x+b
  //Float_t xorg = x;
  if (x > par[9]) x = par[9];
  
  //Double_t cnt    = par[1] / (x*x) + par[2] / x + par[0] ;
  Double_t cnt    = par[0] + par[1] * x + par[2] * x * x ;
  Double_t mean   = par[4] / (x*x) + par[5] / x + par[3] ;
  Double_t sigma  = par[7] / (x*x) + par[8] / x + par[6] ;
 
//   if(xorg > 30)
//     cout<<"En_in = "<<xorg<<"; En_out = "<<x<<"; cnt = "<<cnt
// 	<<"; mean = "<<mean<<"; sigma = "<<sigma<<endl;
      
  //  Double_t arg    = - (y-mean) * (y-mean) / (2*sigma*sigma) ;
  //  return cnt * TMath::Exp(arg) ;
  if(TMath::Abs(sigma) > 1.e-10){
    return cnt*TMath::Gaus(y,mean,sigma);
  }
  else
    return 0.;
 
}
//________________________________________________________________________
Double_t  AliPHOSPIDv1::GausPol2(Double_t  x, Double_t y, Double_t * par)
{
  //Given the energy x and the parameter y (tof, shower dispersion or cpv-emc distance), 
  //this method returns a density probability of this parameter, given by a gaussian 
  //function whose parameters depend with the energy like second order polinomial

  Double_t cnt    = par[0] + par[1] * x + par[2] * x * x ;
  Double_t mean   = par[3] + par[4] * x + par[5] * x * x ;
  Double_t sigma  = par[6] + par[7] * x + par[8] * x * x ;

  if(TMath::Abs(sigma) > 1.e-10){
    return cnt*TMath::Gaus(y,mean,sigma);
  }
  else
    return 0.;
 


}

//____________________________________________________________________________
const TString AliPHOSPIDv1::GetFileNamePrincipal(TString particle) const
{
  //Get file name that contains the PCA for a particle ("photon or pi0")
  particle.ToLower();
  TString name;
  if      (particle=="photon") 
    name = fFileNamePrincipalPhoton ;
  else if (particle=="pi0"   ) 
    name = fFileNamePrincipalPi0    ;
  else    
    AliError(Form("Wrong particle name: %s (choose from pi0/photon)\n",
		  particle.Data()));
  return name;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterCalibration(Int_t i) const 
{
  // Get the i-th parameter "Calibration"
  Float_t param = 0.;
  if (i>2 || i<0) { 
    AliError(Form("Invalid parameter number: %d",i));
  } else
    param = (*fParameters)(0,i);
  return param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetCalibratedEnergy(Float_t e) const
{
//      It calibrates Energy depending on the recpoint energy.
//      The energy of the reconstructed cluster is corrected with 
//      the formula A + B* E  + C* E^2, whose parameters where obtained 
//      through the study of the reconstructed energy distribution of 
//      monoenergetic photons.
 
  Float_t p[]={0.,0.,0.};
  for (Int_t i=0; i<3; i++) p[i] = GetParameterCalibration(i);
  Float_t enerec = p[0] +  p[1]*e + p[2]*e*e;
  return enerec ;

}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterCpv2Emc(Int_t i, TString axis) const 
{
  // Get the i-th parameter "CPV-EMC distance" for the specified axis
  Float_t param = 0.;
  if(i>2 || i<0) {
    AliError(Form("Invalid parameter number: %d",i));
  } else {
    axis.ToLower();
    if      (axis == "x") 
      param = (*fParameters)(1,i);
    else if (axis == "z") 
      param = (*fParameters)(2,i);
    else { 
      AliError(Form("Invalid axis name: %s",axis.Data()));
    }
  }
  return  param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetCpv2EmcDistanceCut(TString axis, Float_t e) const
{
  // Get CpvtoEmcDistance Cut depending on the cluster energy, axis and 
  // Purity-Efficiency point 

  axis.ToLower();
  Float_t p[]={0.,0.,0.};
  for (Int_t i=0; i<3; i++) p[i] = GetParameterCpv2Emc(i,axis);
  Float_t sig = p[0] + TMath::Exp(p[1] - p[2]*e);
  return sig;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetEllipseParameter(TString particle, TString param, Float_t e) const 
{
  // Calculates the parameter param of the ellipse

  particle.ToLower();
  param.   ToLower();
  Float_t p[4]={0.,0.,0.,0.};
  Float_t value = 0.0;
  for (Int_t i=0; i<4; i++) p[i] = GetParameterToCalculateEllipse(particle,param,i);
  if (particle == "photon") {
    if      (param.Contains("a"))  e = TMath::Min((Double_t)e,70.);
    else if (param.Contains("b"))  e = TMath::Min((Double_t)e,70.);
    else if (param.Contains("x0")) e = TMath::Max((Double_t)e,1.1);
  }

 if (particle == "photon")
    value = p[0]/TMath::Sqrt(e) + p[1]*e + p[2]*e*e + p[3];
  else if (particle == "pi0")
    value = p[0] + p[1]*e + p[2]*e*e;

  return value;
}

//_____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterPhotonBoundary (Int_t i) const
{ 
  // Get the parameter "i" to calculate the boundary on the moment M2x
  // for photons at high p_T
  Float_t param = 0;
  if (i>3 || i<0) {
    AliError(Form("Wrong parameter number: %d\n",i));
  } else
    param = (*fParameters)(14,i) ;
  return param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterPi0Boundary (Int_t i) const
{ 
  // Get the parameter "i" to calculate the boundary on the moment M2x
  // for pi0 at high p_T
  Float_t param = 0;
  if (i>2 || i<0) {
    AliError(Form("Wrong parameter number: %d\n",i));
  } else
    param = (*fParameters)(15,i) ;
  return param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterTimeGate(Int_t i) const
{
  // Get TimeGate parameter depending on Purity-Efficiency i:
  // i=0 - Low purity, i=1 - Medium purity, i=2 - High purity
  Float_t param = 0.;
  if(i>2 || i<0) {
    AliError(Form("Invalid Efficiency-Purity choice %d",i));
  } else
    param = (*fParameters)(3,i) ; 
  return param;
}

//_____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterToCalculateEllipse(TString particle, TString param, Int_t i) const
{ 
  // Get the parameter "i" that is needed to calculate the ellipse 
  // parameter "param" for the particle "particle" ("photon" or "pi0")

  particle.ToLower();
  param.   ToLower();
  Int_t offset = -1;
  if      (particle == "photon") 
    offset=0;
  else if (particle == "pi0")    
    offset=5;
  else
    AliError(Form("Wrong particle name: %s (choose from pi0/photon)\n",
		  particle.Data()));

  Int_t p= -1;
  Float_t par = 0;

  if     (param.Contains("a")) p=4+offset; 
  else if(param.Contains("b")) p=5+offset; 
  else if(param.Contains("c")) p=6+offset; 
  else if(param.Contains("x0"))p=7+offset; 
  else if(param.Contains("y0"))p=8+offset;

  if      (i>4 || i<0) {
    AliError(Form("No parameter with index %d", i)) ; 
  } else if (p==-1) {
    AliError(Form("No parameter with name %s", param.Data() )) ; 
  } else
    par = (*fParameters)(p,i) ;
  
  return par;
}


//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSCpvRecPoint * cpv, Option_t *  axis)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  
  const AliPHOSGeometry * geom = AliPHOSGetter::Instance()->PHOSGeometry() ; 
  TVector3 vecEmc ;
  TVector3 vecCpv ;
  if(cpv){
    emc->GetLocalPosition(vecEmc) ;
    cpv->GetLocalPosition(vecCpv) ; 
    
    if(emc->GetPHOSMod() == cpv->GetPHOSMod()){      
      // Correct to difference in CPV and EMC position due to different distance to center.
      // we assume, that particle moves from center
      Float_t dCPV = geom->GetIPtoOuterCoverDistance();
      Float_t dEMC = geom->GetIPtoCrystalSurface() ;
      dEMC         = dEMC / dCPV ;
      vecCpv = dEMC * vecCpv  - vecEmc ; 
      if (axis == "X") return vecCpv.X();
      if (axis == "Y") return vecCpv.Y();
      if (axis == "Z") return vecCpv.Z();
      if (axis == "R") return vecCpv.Mag();
    }
    return 100000000 ;
  }
  return 100000000 ;
}
//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetCPVBit(AliPHOSEmcRecPoint * emc,AliPHOSCpvRecPoint * cpv, Int_t effPur, Float_t e) const
{
  //Calculates the pid bit for the CPV selection per each purity.
  if(effPur>2 || effPur<0)
    AliError(Form("Invalid Efficiency-Purity choice %d",effPur));
  
  Float_t sigX = GetCpv2EmcDistanceCut("X",e);
  Float_t sigZ = GetCpv2EmcDistanceCut("Z",e);
  
  Float_t deltaX = TMath::Abs(GetDistance(emc, cpv,  "X"));
  Float_t deltaZ = TMath::Abs(GetDistance(emc, cpv,  "Z"));
  //Info("GetCPVBit"," xdist %f, sigx %f, zdist %f, sigz %f",deltaX, sigX, deltaZ,sigZ) ;
 
  //if(deltaX>sigX*(effPur+1))
  //if((deltaX>sigX*(effPur+1)) || (deltaZ>sigZ*(effPur+1)))
  if((deltaX>sigX*(effPur+1)) && (deltaZ>sigZ*(effPur+1)))
    return 1;//Neutral
  else
    return 0;//Charged
}

//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetPrincipalBit(TString particle, const Double_t* p, Int_t effPur, Float_t e)const
{
  //Is the particle inside de PCA ellipse?
  
  particle.ToLower();
  Int_t    prinbit  = 0 ;
  Float_t a  = GetEllipseParameter(particle,"a" , e); 
  Float_t b  = GetEllipseParameter(particle,"b" , e);
  Float_t c  = GetEllipseParameter(particle,"c" , e);
  Float_t x0 = GetEllipseParameter(particle,"x0", e); 
  Float_t y0 = GetEllipseParameter(particle,"y0", e);
  
  Float_t r = TMath::Power((p[0] - x0)/a,2) + 
              TMath::Power((p[1] - y0)/b,2) +
            c*(p[0] -  x0)*(p[1] - y0)/(a*b) ;
  //3 different ellipses defined
  if((effPur==2) && (r<1./2.)) prinbit= 1;
  if((effPur==1) && (r<2.   )) prinbit= 1;
  if((effPur==0) && (r<9./2.)) prinbit= 1;

  if(r<0)
    AliError("Negative square?") ;

  return prinbit;

}
//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetHardPhotonBit(AliPHOSEmcRecPoint * emc) const
{
  // Set bit for identified hard photons (E > 30 GeV)
  // if the second moment M2x is below the boundary

  Float_t e   = emc->GetEnergy();
  if (e < 30.0) return 0;
  Float_t m2x = emc->GetM2x();
  Float_t m2xBoundary = GetParameterPhotonBoundary(0) *
    TMath::Exp(-TMath::Power(e-GetParameterPhotonBoundary(1),2)/2.0/
	        TMath::Power(GetParameterPhotonBoundary(2),2)) +
    GetParameterPhotonBoundary(3);
  AliDebug(1, Form("GetHardPhotonBit","E=%f, m2x=%f, boundary=%f",
		       e,m2x,m2xBoundary));
  if (m2x < m2xBoundary)
    return 1;// A hard photon
  else
    return 0;// Not a hard photon
}

//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetHardPi0Bit(AliPHOSEmcRecPoint * emc) const
{
  // Set bit for identified hard pi0  (E > 30 GeV)
  // if the second moment M2x is above the boundary

  Float_t e   = emc->GetEnergy();
  if (e < 30.0) return 0;
  Float_t m2x = emc->GetM2x();
  Float_t m2xBoundary = GetParameterPi0Boundary(0) +
                    e * GetParameterPi0Boundary(1);
  AliDebug(1,Form("E=%f, m2x=%f, boundary=%f",e,m2x,m2xBoundary));
  if (m2x > m2xBoundary)
    return 1;// A hard pi0
  else
    return 0;// Not a hard pi0
}

//____________________________________________________________________________
TVector3 AliPHOSPIDv1::GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSCpvRecPoint * )const 
{ 
  // Calculates the momentum direction:
  //   1. if only a EMC RecPoint, direction is given by IP and this RecPoint
  //   2. if a EMC RecPoint and CPV RecPoint, direction is given by the line through the 2 recpoints 
  //  However because of the poor position resolution of PPSD the direction is always taken as if we were 
  //  in case 1.

  TVector3 dir(0,0,0) ; 
  TMatrixF  dummy ;
  
  emc->GetGlobalPosition(dir, dummy) ;

  //account correction to the position of IP
  Float_t xo,yo,zo ; //Coordinates of the origin
  if(gAlice && gAlice->GetMCApp() && gAlice->Generator()){
    gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  }
  else{
    xo=yo=zo=0.;
  }
  TVector3 origin(xo,yo,zo);
  dir = dir - origin ;
  dir.SetMag(1.) ;

  return dir ;  
}

//________________________________________________________________________
Double_t  AliPHOSPIDv1::LandauF(Double_t  x, Double_t y, Double_t * par)
{
  //Given the energy x and the parameter y (tof, shower dispersion or cpv-emc distance), 
  //this method returns a density probability of this parameter, given by a landau 
  //function whose parameters depend with the energy  with a function: a/(x*x)+b/x+b

  if (x > par[9]) x = par[9];

  //Double_t cnt    = par[1] / (x*x) + par[2] / x + par[0] ;
  Double_t cnt    = par[0] + par[1] * x + par[2] * x * x ;
  Double_t mean   = par[4] / (x*x) + par[5] / x + par[3] ;
  Double_t sigma  = par[7] / (x*x) + par[8] / x + par[6] ;

  if(TMath::Abs(sigma) > 1.e-10){
    return cnt*TMath::Landau(y,mean,sigma);
  }
  else
    return 0.;

}
//________________________________________________________________________
Double_t  AliPHOSPIDv1::LandauPol2(Double_t  x, Double_t y, Double_t * par)
{

  //Given the energy x and the parameter y (tof, shower dispersion or cpv-emc distance), 
  //this method returns a density probability of this parameter, given by a landau 
  //function whose parameters depend with the energy like second order polinomial

  Double_t cnt    = par[2] * (x*x) + par[1] * x + par[0] ;
  Double_t mean   = par[5] * (x*x) + par[4] * x + par[3] ;
  Double_t sigma  = par[8] * (x*x) + par[7] * x + par[6] ;

   if(TMath::Abs(sigma) > 1.e-10){
    return cnt*TMath::Landau(y,mean,sigma);
  }
  else
    return 0.;


}
// //________________________________________________________________________
// Double_t  AliPHOSPIDv1::ChargedHadronDistProb(Double_t  x, Double_t y, Double_t * parg, Double_t * parl)
// {
//   Double_t cnt   = 0.0 ;
//   Double_t mean  = 0.0 ;
//   Double_t sigma = 0.0 ;
//   Double_t arg   = 0.0 ;
//   if (y < parl[4] / (x*x) + parl[5] / x + parl[3]){
//     cnt    = parg[1] / (x*x) + parg[2] / x + parg[0] ;
//     mean   = parg[4] / (x*x) + parg[5] / x + parg[3] ;
//     sigma  = parg[7] / (x*x) + parg[8] / x + parg[6] ;
//     TF1 * f = new TF1("gaus","gaus",0.,100.);
//     f->SetParameters(cnt,mean,sigma);
//     arg  = f->Eval(y) ;
//   }
//   else{
//     cnt    = parl[1] / (x*x) + parl[2] / x + parl[0] ;
//     mean   = parl[4] / (x*x) + parl[5] / x + parl[3] ;
//     sigma  = parl[7] / (x*x) + parl[8] / x + parl[6] ;
//     TF1 * f = new TF1("landau","landau",0.,100.);
//     f->SetParameters(cnt,mean,sigma);
//     arg  = f->Eval(y) ;
//   }
//   //  Double_t mean   = par[3] + par[4] * x + par[5] * x * x ;
//   //   Double_t sigma  = par[6] + par[7] * x + par[8] * x * x ;
  
//   //Double_t arg    = -(y-mean)*(y-mean)/(2*sigma*sigma) ;
//   //return cnt * TMath::Exp(arg) ;
  
//   return arg;
  
// }
//____________________________________________________________________________
void  AliPHOSPIDv1::MakePID()
{
  // construct the PID weight from a Bayesian Method
  
  const Int_t kSPECIES = AliPID::kSPECIESN ;
 
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 

  Int_t nparticles = gime->RecParticles()->GetEntriesFast() ;

  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ;
  if ( !emcRecPoints || !cpvRecPoints || !trackSegments ) {
    AliFatal("RecPoints or TrackSegments not found !") ;  
  }
  TIter next(trackSegments) ; 
  AliPHOSTrackSegment * ts ; 
  Int_t index = 0 ; 

  Double_t * stof[kSPECIES] ;
  Double_t * sdp [kSPECIES]  ;
  Double_t * scpv[kSPECIES] ;
  Double_t * sw  [kSPECIES] ;
  //Info("MakePID","Begin MakePID"); 
  
  for (Int_t i =0; i< kSPECIES; i++){
    stof[i] = new Double_t[nparticles] ;
    sdp [i] = new Double_t[nparticles] ;
    scpv[i] = new Double_t[nparticles] ;
    sw  [i] = new Double_t[nparticles] ;
  }
  

  while ( (ts = (AliPHOSTrackSegment *)next()) ) {
    
    //cout<<">>>>>> Bayesian Index "<<index<<endl;

    AliPHOSEmcRecPoint * emc = 0 ;
    if(ts->GetEmcIndex()>=0)
      emc = (AliPHOSEmcRecPoint *) emcRecPoints->At(ts->GetEmcIndex()) ;
    
    AliPHOSCpvRecPoint * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSCpvRecPoint *) cpvRecPoints->At(ts->GetCpvIndex()) ;
    
//     Int_t track = 0 ; 
//     track = ts->GetTrackIndex() ; //TPC tracks ?
    
    if (!emc) {
      AliFatal(Form("-> emc(%d) = %d", ts->GetEmcIndex(), emc )) ;
    }


    // ############Tof#############################

    //    Info("MakePID", "TOF");
    Float_t  en   = emc->GetEnergy();    
    Double_t time = emc->GetTime() ;
    //    cout<<">>>>>>>Energy "<<en<<"Time "<<time<<endl;
   
    // now get the signals probability
    // s(pid) in the Bayesian formulation
    
    stof[AliPID::kPhoton][index]   = 1.; 
    stof[AliPID::kElectron][index] = 1.;
    stof[AliPID::kEleCon][index]   = 1.;
    //We assing the same prob to charged hadrons, sum is 1
    stof[AliPID::kPion][index]     = 1./3.; 
    stof[AliPID::kKaon][index]     = 1./3.; 
    stof[AliPID::kProton][index]   = 1./3.;
    //We assing the same prob to neutral hadrons, sum is 1
    stof[AliPID::kNeutron][index]  = 1./2.;
    stof[AliPID::kKaon0][index]    = 1./2.;
    stof[AliPID::kMuon][index]     = 1.; 
 
    if(en <  fTOFEnThreshold) {

      Double_t pTofPion = fTFpiong ->Eval(time) ; //gaus distribution
      Double_t pTofKaon = 0;

      if(time < fTkaonl[1])
	pTofKaon = fTFkaong  ->Eval(time) ; //gaus distribution
      else 
	pTofKaon = fTFkaonl  ->Eval(time) ; //landau distribution

      Double_t pTofNucleon = 0;

      if(time < fThhadronl[1])
	pTofNucleon = fTFhhadrong   ->Eval(time) ; //gaus distribution
      else
	pTofNucleon = fTFhhadronl   ->Eval(time) ; //landau distribution
      //We assing the same prob to neutral hadrons, sum is the average prob
      Double_t pTofNeHadron =  (pTofKaon + pTofNucleon)/2. ;
      //We assing the same prob to charged hadrons, sum is the average prob
      Double_t pTofChHadron =  (pTofPion + pTofKaon + pTofNucleon)/3. ;

      stof[AliPID::kPhoton][index]   = fTFphoton     ->Eval(time) ; 
      //gaus distribution
      stof[AliPID::kEleCon][index]   = stof[AliPID::kPhoton][index] ; 
      //a conversion electron has the photon ToF
      stof[AliPID::kMuon][index]     = stof[AliPID::kPhoton][index] ;
 
      stof[AliPID::kElectron][index] = pTofPion  ;                             

      stof[AliPID::kPion][index]     =  pTofChHadron ; 
      stof[AliPID::kKaon][index]     =  pTofChHadron ;
      stof[AliPID::kProton][index]   =  pTofChHadron ;

      stof[AliPID::kKaon0][index]    =  pTofNeHadron ;     
      stof[AliPID::kNeutron][index]  =  pTofNeHadron ;            
    } 
    
    //    Info("MakePID", "Dispersion");
    
    // ###########Shower shape: Dispersion####################
    Float_t dispersion = emc->GetDispersion();
    //dispersion is not well defined if the cluster is only in few crystals
    
    sdp[AliPID::kPhoton][index]   = 1. ;
    sdp[AliPID::kElectron][index] = 1. ;
    sdp[AliPID::kPion][index]     = 1. ; 
    sdp[AliPID::kKaon][index]     = 1. ; 
    sdp[AliPID::kProton][index]   = 1. ;
    sdp[AliPID::kNeutron][index]  = 1. ;
    sdp[AliPID::kEleCon][index]   = 1. ; 
    sdp[AliPID::kKaon0][index]    = 1. ; 
    sdp[AliPID::kMuon][index]     = 1. ; 
    
    if(en > fDispEnThreshold && emc->GetMultiplicity() >  fDispMultThreshold){
      sdp[AliPID::kPhoton][index]   = GausF(en , dispersion, fDphoton) ;
      sdp[AliPID::kElectron][index] = sdp[AliPID::kPhoton][index] ;
      sdp[AliPID::kPion][index]     = LandauF(en , dispersion, fDhadron ) ; 
      sdp[AliPID::kKaon][index]     = sdp[AliPID::kPion][index]  ; 
      sdp[AliPID::kProton][index]   = sdp[AliPID::kPion][index]  ;
      sdp[AliPID::kNeutron][index]  = sdp[AliPID::kPion][index]  ;
      sdp[AliPID::kEleCon][index]   = sdp[AliPID::kPhoton][index]; 
      sdp[AliPID::kKaon0][index]    = sdp[AliPID::kPion][index]  ; 
      sdp[AliPID::kMuon][index]     = fDFmuon ->Eval(dispersion) ; 
      //landau distribution
    }
    
//      Info("MakePID","multiplicity %d, dispersion %f", emc->GetMultiplicity(), dispersion);
//      Info("MakePID","ss: photon %f, hadron %f ",  sdp[AliPID::kPhoton][index],  sdp[AliPID::kPion][index]);
//       cout<<">>>>>multiplicity "<<emc->GetMultiplicity()<<", dispersion "<< dispersion<<endl ;
//       cout<<"<<<<<ss: photon   "<<sdp[AliPID::kPhoton][index]<<", hadron    "<<sdp[AliPID::kPion][index]<<endl;

    //########## CPV-EMC  Distance#######################
    //     Info("MakePID", "Distance");

    Float_t x             = TMath::Abs(GetDistance(emc, cpv,  "X")) ;
    Float_t z             = GetDistance(emc, cpv,  "Z") ;
   
    Double_t pcpv         = 0 ;
    Double_t pcpvneutral  = 0. ;
   
    Double_t elprobx      = GausF(en , x, fXelectron) ;
    Double_t elprobz      = GausF(en , z, fZelectron) ;
    Double_t chprobx      = GausF(en , x, fXcharged)  ;
    Double_t chprobz      = GausF(en , z, fZcharged)  ;
    Double_t pcpvelectron = elprobx * elprobz;
    Double_t pcpvcharged  = chprobx * chprobz;
  
//     cout<<">>>>energy "<<en<<endl;
//     cout<<">>>>electron : x "<<x<<" xprob "<<elprobx<<" z "<<z<<" zprob "<<elprobz<<endl;
//     cout<<">>>>hadron   : x "<<x<<" xprob "<<chprobx<<" z "<<z<<" zprob "<<chprobz<<endl;
//     cout<<">>>>electron : px*pz "<<pcpvelectron <<" hadron: px*pz "<<pcpvcharged<<endl;  

    // Is neutral or charged?
    if(pcpvelectron >= pcpvcharged)  
      pcpv = pcpvelectron ;
    else
      pcpv = pcpvcharged ;
    
    if(pcpv < fChargedNeutralThreshold)
      {
	pcpvneutral  = 1. ;
	pcpvcharged  = 0. ;
	pcpvelectron = 0. ;
      }
    //    else
    //      cout<<">>>>>>>>>>>CHARGED>>>>>>>>>>>"<<endl;
    
    scpv[AliPID::kPion][index]     =  pcpvcharged  ; 
    scpv[AliPID::kKaon][index]     =  pcpvcharged  ; 
    scpv[AliPID::kProton][index]   =  pcpvcharged  ;

    scpv[AliPID::kMuon][index]     =  pcpvelectron ; 
    scpv[AliPID::kElectron][index] =  pcpvelectron ;
    scpv[AliPID::kEleCon][index]   =  pcpvelectron ; 

    scpv[AliPID::kPhoton][index]   =  pcpvneutral  ;
    scpv[AliPID::kNeutron][index]  =  pcpvneutral  ; 
    scpv[AliPID::kKaon0][index]    =  pcpvneutral  ; 

    
    //   Info("MakePID", "CPV passed");

    //############## Pi0 #############################
    stof[AliPID::kPi0][index]      = 0. ;  
    scpv[AliPID::kPi0][index]      = 0. ;
    sdp [AliPID::kPi0][index]      = 0. ;

    if(en > 30.){
      // pi0 are detected via decay photon
      stof[AliPID::kPi0][index]  =   stof[AliPID::kPhoton][index];
      scpv[AliPID::kPi0][index]  = pcpvneutral  ;
      if(emc->GetMultiplicity() >  fDispMultThreshold)
	sdp [AliPID::kPi0][index]  = GausF(en , dispersion, fDpi0) ;
	//sdp [AliPID::kPi0][index]  = GausPol2(en , dispersion, fDpi0) ;
//       cout<<"E = "<<en<<" GeV; disp = "<<dispersion<<"; mult = "
// 	  <<emc->GetMultiplicity()<<endl;
//       cout<<"PDF: photon = "<<sdp [AliPID::kPhoton][index]<<"; pi0 = "
// 	  <<sdp [AliPID::kPi0][index]<<endl;
    }
    
  

    
    //############## muon #############################

    if(en > 0.5){
      //Muons deposit few energy
      scpv[AliPID::kMuon][index]     =  0 ;
      stof[AliPID::kMuon][index]     =  0 ;
      sdp [AliPID::kMuon][index]     =  0 ;
    }

    //Weight to apply to hadrons due to energy reconstruction

    Float_t weight = fERecWeight ->Eval(en) ;
 
    sw[AliPID::kPhoton][index]   = 1. ;
    sw[AliPID::kElectron][index] = 1. ;
    sw[AliPID::kPion][index]     = weight ; 
    sw[AliPID::kKaon][index]     = weight ; 
    sw[AliPID::kProton][index]   = weight ;
    sw[AliPID::kNeutron][index]  = weight ;
    sw[AliPID::kEleCon][index]   = 1. ; 
    sw[AliPID::kKaon0][index]    = weight ; 
    sw[AliPID::kMuon][index]     = weight ; 
    sw[AliPID::kPi0][index]      = 1. ;

//     if(en > 0.5){
//       cout<<"######################################################"<<endl;
//       //cout<<"MakePID: energy "<<en<<", tof "<<time<<", distance "<<distance<<", dispersion "<<dispersion<<endl ;
//       cout<<"MakePID: energy "<<en<<", tof "<<time<<", dispersion "<<dispersion<<", x "<<x<<", z "<<z<<endl ;
//       cout<<">>>>>multiplicity "<<emc->GetMultiplicity()<<endl;
//       cout<<">>>>electron : xprob "<<elprobx<<" zprob "<<elprobz<<endl;
//       cout<<">>>>hadron   : xprob "<<chprobx<<" zprob "<<chprobz<<endl;
//       cout<<">>>>electron : px*pz "<<pcpvelectron <<" hadron: px*pz "<<pcpvcharged<<endl;  
      
//        cout<<"Photon   , pid "<< fInitPID[AliPID::kPhoton]<<" tof "<<stof[AliPID::kPhoton][index]
//  	  <<", cpv "<<scpv[AliPID::kPhoton][index]<<", ss "<<sdp[AliPID::kPhoton][index]<<endl;
//       cout<<"EleCon   , pid "<< fInitPID[AliPID::kEleCon]<<", tof "<<stof[AliPID::kEleCon][index]
// 	  <<", cpv "<<scpv[AliPID::kEleCon][index]<<" ss "<<sdp[AliPID::kEleCon][index]<<endl;
//       cout<<"Electron , pid "<< fInitPID[AliPID::kElectron]<<", tof "<<stof[AliPID::kElectron][index]
// 	  <<", cpv "<<scpv[AliPID::kElectron][index]<<" ss "<<sdp[AliPID::kElectron][index]<<endl;
//       cout<<"Muon     , pid "<< fInitPID[AliPID::kMuon]<<", tof "<<stof[AliPID::kMuon][index]
// 	  <<", cpv "<<scpv[AliPID::kMuon][index]<<" ss "<<sdp[AliPID::kMuon][index]<<endl;
//        cout<<"Pi0      , pid "<< fInitPID[AliPID::kPi0]<<", tof "<<stof[AliPID::kPi0][index]
//  	  <<", cpv "<<scpv[AliPID::kPi0][index]<<" ss "<<sdp[AliPID::kPi0][index]<<endl;
//       cout<<"Pion     , pid "<< fInitPID[AliPID::kPion]<<", tof "<<stof[AliPID::kPion][index]
// 	  <<", cpv "<<scpv[AliPID::kPion][index]<<" ss "<<sdp[AliPID::kPion][index]<<endl;
//       cout<<"Kaon0    , pid "<< fInitPID[AliPID::kKaon0]<<", tof "<<stof[AliPID::kKaon0][index]
// 	  <<", cpv "<<scpv[AliPID::kKaon0][index]<<" ss "<<sdp[AliPID::kKaon0][index]<<endl;
//       cout<<"Kaon     , pid "<< fInitPID[AliPID::kKaon]<<", tof "<<stof[AliPID::kKaon][index]
// 	  <<", cpv "<<scpv[AliPID::kKaon][index]<<" ss "<<sdp[AliPID::kKaon][index]<<endl;
//       cout<<"Neutron  , pid "<< fInitPID[AliPID::kNeutron]<<", tof "<<stof[AliPID::kNeutron][index]
// 	  <<", cpv "<<scpv[AliPID::kNeutron][index]<<" ss "<<sdp[AliPID::kNeutron][index]<<endl;
//       cout<<"Proton   , pid "<< fInitPID[AliPID::kProton]<<", tof "<<stof[AliPID::kProton][index]
// 	  <<", cpv "<<scpv[AliPID::kProton][index]<<" ss "<<sdp[AliPID::kProton][index]<<endl;
//       cout<<"######################################################"<<endl;
//     }
      index++;
  }
  
  //for (index = 0 ; index < kSPECIES ; index++) 
  // pid[index] /= nparticles ; 
  

  //  Info("MakePID", "Total Probability calculation");
  
  for(index = 0 ; index < nparticles ; index ++) {
    
    AliPHOSRecParticle * recpar = gime->RecParticle(index) ;  
    
    //Conversion electron?
    
    if(recpar->IsEleCon()){
      fInitPID[AliPID::kEleCon]   = 1. ;
      fInitPID[AliPID::kPhoton]   = 0. ;
      fInitPID[AliPID::kElectron] = 0. ;
    }
    else{
      fInitPID[AliPID::kEleCon]   = 0. ;
      fInitPID[AliPID::kPhoton]   = 1. ;
      fInitPID[AliPID::kElectron] = 1. ;
    }
    //	fInitPID[AliPID::kEleCon]   = 0. ;
    
    
    // calculates the Bayesian weight
    
    Int_t jndex ;
    Double_t wn = 0.0 ; 
    for (jndex = 0 ; jndex < kSPECIES ; jndex++) 
      wn += stof[jndex][index] * sdp[jndex][index]  * scpv[jndex][index] * 
	sw[jndex][index] * fInitPID[jndex] ;
    
    //    cout<<"*************wn "<<wn<<endl;
    if (TMath::Abs(wn)>0)
      for (jndex = 0 ; jndex < kSPECIES ; jndex++) {
	//cout<<"jndex "<<jndex<<" wn "<<wn<<" SetPID * wn"
	//<<stof[jndex][index] * sdp[jndex][index] * pid[jndex]  << endl;
	//cout<<" tof "<<stof[jndex][index] << " disp " <<sdp[jndex][index] << " pid "<< fInitPID[jndex] << endl;
	// 	if(jndex ==  AliPID::kPi0 || jndex ==  AliPID::kPhoton){
	// 	  cout<<"Particle "<<jndex<<"  final prob * wn   "
	// 	      <<stof[jndex][index] * sdp[jndex][index] * scpv[jndex][index] * 
	// 	    fInitPID[jndex] <<"  wn  "<< wn<<endl;
	// 	  cout<<"pid "<< fInitPID[jndex]<<", tof "<<stof[jndex][index]
	// 	      <<", cpv "<<scpv[jndex][index]<<" ss "<<sdp[jndex][index]<<endl;
	// 	}
	recpar->SetPID(jndex, stof[jndex][index] * sdp[jndex][index] * 
		       sw[jndex][index] * scpv[jndex][index] * 
		       fInitPID[jndex] / wn) ; 
      }
  }
  //  Info("MakePID", "Delete");
  
  for (Int_t i =0; i< kSPECIES; i++){
    delete [] stof[i];
    delete [] sdp [i];
    delete [] scpv[i];
    delete [] sw  [i];
  }
  //  Info("MakePID","End MakePID"); 
}

//____________________________________________________________________________
void  AliPHOSPIDv1::MakeRecParticles()
{
  // Makes a RecParticle out of a TrackSegment
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 
  if ( !emcRecPoints || !cpvRecPoints || !trackSegments ) {
    AliFatal("RecPoints or TrackSegments not found !") ;  
  }
  TClonesArray * recParticles  = gime->RecParticles() ; 
  recParticles->Clear();

  TIter next(trackSegments) ; 
  AliPHOSTrackSegment * ts ; 
  Int_t index = 0 ; 
  AliPHOSRecParticle * rp ; 
  while ( (ts = (AliPHOSTrackSegment *)next()) ) {
    //  cout<<">>>>>>>>>>>>>>>PCA Index "<<index<<endl;
    new( (*recParticles)[index] ) AliPHOSRecParticle() ;
    rp = (AliPHOSRecParticle *)recParticles->At(index) ; 
    rp->SetTrackSegment(index) ;
    rp->SetIndexInList(index) ;
    	
    AliPHOSEmcRecPoint * emc = 0 ;
    if(ts->GetEmcIndex()>=0)
      emc = (AliPHOSEmcRecPoint *) emcRecPoints->At(ts->GetEmcIndex()) ;
    
    AliPHOSCpvRecPoint * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSCpvRecPoint *) cpvRecPoints->At(ts->GetCpvIndex()) ;
    
    Int_t track = 0 ; 
    track = ts->GetTrackIndex() ; 
      
    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    
    if (!emc) {
      AliFatal(Form("-> emc(%d) = %d", ts->GetEmcIndex(), emc )) ;
    }

    Float_t e = emc->GetEnergy() ;   
    
    Float_t  lambda[2] ;
    emc->GetElipsAxis(lambda) ;
    
    if((lambda[0]>0.01) && (lambda[1]>0.01)){
      // Looking PCA. Define and calculate the data (X),
      // introduce in the function X2P that gives the components (P).  

      Float_t  spher = 0. ;
      Float_t  emaxdtotal = 0. ; 
      
      if((lambda[0]+lambda[1])!=0) 
	spher=TMath::Abs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
      
      emaxdtotal=emc->GetMaximalEnergy()/emc->GetEnergy(); 
      
      fX[0] = lambda[0] ;  
      fX[1] = lambda[1] ; 
      fX[2] = emc->GetDispersion() ; 
      fX[3] = spher ; 
      fX[4] = emc->GetMultiplicity() ;  
      fX[5] = emaxdtotal ;  
      fX[6] = emc->GetCoreEnergy() ;  
      
      fPrincipalPhoton->X2P(fX,fPPhoton);
      fPrincipalPi0   ->X2P(fX,fPPi0);

    }
    else{
      fPPhoton[0]=-100.0;  //We do not accept clusters with 
      fPPhoton[1]=-100.0;  //one cell as a photon-like
      fPPi0[0]   =-100.0;
      fPPi0[1]   =-100.0;
    }
    
    Float_t time = emc->GetTime() ;
    rp->SetTof(time) ; 
    
    // Loop of Efficiency-Purity (the 3 points of purity or efficiency 
    // are taken into account to set the particle identification)
    for(Int_t effPur = 0; effPur < 3 ; effPur++){
      
      // Looking at the CPV detector. If RCPV greater than CpvEmcDistance, 
      // 1st,2nd or 3rd bit (depending on the efficiency-purity point )
      // is set to 1
      if(GetCPVBit(emc, cpv, effPur,e) == 1 ){  
	rp->SetPIDBit(effPur) ;
	//cout<<"CPV bit "<<effPur<<endl;
      }
      // Looking the TOF. If TOF smaller than gate,  4th, 5th or 6th 
      // bit (depending on the efficiency-purity point )is set to 1             
      if(time< (*fParameters)(3,effPur)) 
	rp->SetPIDBit(effPur+3) ;		    
  
      //Photon PCA
      //If we are inside the ellipse, 7th, 8th or 9th 
      // bit (depending on the efficiency-purity point )is set to 1 
      if(GetPrincipalBit("photon",fPPhoton,effPur,e) == 1) 
	rp->SetPIDBit(effPur+6) ;

      //Pi0 PCA
      //If we are inside the ellipse, 10th, 11th or 12th 
      // bit (depending on the efficiency-purity point )is set to 1 
      if(GetPrincipalBit("pi0"   ,fPPi0   ,effPur,e) == 1) 
	rp->SetPIDBit(effPur+9) ;
    }
    if(GetHardPhotonBit(emc))
      rp->SetPIDBit(12) ;
    if(GetHardPi0Bit   (emc))
      rp->SetPIDBit(13) ;
    
    if(track >= 0) 
      rp->SetPIDBit(14) ; 

    //Set momentum, energy and other parameters 
    Float_t  encal = GetCalibratedEnergy(e);
    TVector3 dir   = GetMomentumDirection(emc,cpv) ; 
    dir.SetMag(encal) ;
    rp->SetMomentum(dir.X(),dir.Y(),dir.Z(),encal) ;
    rp->SetCalcMass(0);
    rp->Name(); //If photon sets the particle pdg name to gamma
    rp->SetProductionVertex(0,0,0,0);
    rp->SetFirstMother(-1);
    rp->SetLastMother(-1);
    rp->SetFirstDaughter(-1);
    rp->SetLastDaughter(-1);
    rp->SetPolarisation(0,0,0);
    //Set the position in global coordinate system from the RecPoint
    AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
    AliPHOSTrackSegment * ts  = gime->TrackSegment(rp->GetPHOSTSIndex()) ; 
    AliPHOSEmcRecPoint  * erp = gime->EmcRecPoint(ts->GetEmcIndex()) ; 
    TVector3 pos ; 
    geom->GetGlobal(erp, pos) ; 
    rp->SetPos(pos);
    index++ ; 
  }
}
  
//____________________________________________________________________________
void  AliPHOSPIDv1::Print(const Option_t *) const
{
  // Print the parameters used for the particle type identification

    AliInfo("=============== AliPHOSPIDv1 ================") ;
    printf("Making PID\n") ;
    printf("    Pricipal analysis file from 0.5 to 100 %s\n", fFileNamePrincipalPhoton.Data() )   ; 
    printf("    Name of parameters file     %s\n", fFileNameParameters.Data() )  ;
    printf("    Matrix of Parameters: 14x4\n") ;
    printf("        Energy Calibration  1x3 [3 parametres to calibrate energy: A + B* E + C * E^2]\n") ;
    printf("        RCPV 2x3 rows x and z, columns function cut parameters\n") ;
    printf("        TOF  1x3 [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]\n") ;
    printf("        PCA  5x4 [5 ellipse parametres and 4 parametres to calculate them: A/Sqrt(E) + B* E + C * E^2 + D]\n") ;
    Printf("    Pi0 PCA  5x3 [5 ellipse parametres and 3 parametres to calculate them: A + B* E + C * E^2]\n") ;
    fParameters->Print() ;
}



//____________________________________________________________________________
void AliPHOSPIDv1::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliPHOSGetter *gime = AliPHOSGetter::Instance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 

  TString message ; 
  message  = "\nevent " ;
  message += gAlice->GetEvNumber() ; 
  message += "       found " ; 
  message += recParticles->GetEntriesFast(); 
  message += " RecParticles\n" ; 

  if(strstr(option,"all")) {  // printing found TS
    message += "\n  PARTICLE         Index    \n" ; 
    
    Int_t index ;
    for (index = 0 ; index < recParticles->GetEntries() ; index++) {
      AliPHOSRecParticle * rp = (AliPHOSRecParticle * ) recParticles->At(index) ;       
      message += "\n" ;
      message += rp->Name().Data() ;  
      message += " " ;
      message += rp->GetIndexInList() ;  
      message += " " ;
      message += rp->GetType()  ;
    }
  }
  AliInfo(message.Data() ) ; 
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameters() 
{
  // PCA : To do the Principal Components Analysis it is necessary 
  // the Principal file, which is opened here
  fX       = new double[7]; // Data for the PCA 
  fPPhoton = new double[7]; // Eigenvalues of the PCA
  fPPi0    = new double[7]; // Eigenvalues of the Pi0 PCA

  // Read photon principals from the photon file
  
  fFileNamePrincipalPhoton = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root" ; 
  TFile f( fFileNamePrincipalPhoton.Data(), "read" ) ;
  fPrincipalPhoton = dynamic_cast<TPrincipal*> (f.Get("principal")) ; 
  f.Close() ; 

  // Read pi0 principals from the pi0 file

  fFileNamePrincipalPi0    = "$ALICE_ROOT/PHOS/PCA_pi0_40-120.root" ;
  TFile fPi0( fFileNamePrincipalPi0.Data(), "read" ) ;
  fPrincipalPi0    = dynamic_cast<TPrincipal*> (fPi0.Get("principal")) ; 
  fPi0.Close() ;

  // Open parameters file and initialization of the Parameters matrix. 
  // In the File Parameters.dat are all the parameters. These are introduced 
  // in a matrix of 16x4  
  // 
  // All the parameters defined in this file are, in order of row: 
  // line   0   : calibration 
  // lines  1,2 : CPV rectangular cat for X and Z
  // line   3   : TOF cut
  // lines  4-8 : parameters to calculate photon PCA ellipse
  // lines  9-13: parameters to calculate pi0 PCA ellipse
  // lines 14-15: parameters to calculate border for high-pt photons and pi0

  fFileNameParameters = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters.dat");
  fParameters = new TMatrixF(16,4) ;
  const Int_t kMaxLeng=255;
  char string[kMaxLeng];

  // Open a text file with PID parameters
  FILE *fd = fopen(fFileNameParameters.Data(),"r");
  if (!fd)
    AliFatal(Form("File %s with a PID parameters cannot be opened\n",
	  fFileNameParameters.Data()));

  Int_t i=0;
  // Read parameter file line-by-line and skip empty line and comments
  while (fgets(string,kMaxLeng,fd) != NULL) {
    if (string[0] == '\n' ) continue;
    if (string[0] == '!'  ) continue;
    sscanf(string, "%f %f %f %f",
	   &(*fParameters)(i,0), &(*fParameters)(i,1), 
	   &(*fParameters)(i,2), &(*fParameters)(i,3));
    i++;
    AliDebug(1, Form("SetParameters", "line %d: %s",i,string));
  }
  fclose(fd);
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterCalibration(Int_t i,Float_t param) 
{
  // Set parameter "Calibration" i to a value param
  if(i>2 || i<0) {
    AliError(Form("Invalid parameter number: %d",i));
  } else
    (*fParameters)(0,i) = param ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterCpv2Emc(Int_t i, TString axis, Float_t cut) 
{
  // Set the parameters to calculate Cpv-to-Emc Distance Cut depending on 
  // Purity-Efficiency point i

  if(i>2 || i<0) {
    AliError(Form("Invalid parameter number: %d",i));
  } else {
    axis.ToLower();
    if      (axis == "x") (*fParameters)(1,i) = cut;
    else if (axis == "z") (*fParameters)(2,i) = cut;
    else { 
      AliError(Form("Invalid axis name: %s",axis.Data()));
    }
  }
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterPhotonBoundary(Int_t i,Float_t param) 
{
  // Set parameter "Hard photon boundary" i to a value param
  if(i>4 || i<0) {
    AliError(Form("Invalid parameter number: %d",i));
  } else
    (*fParameters)(14,i) = param ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterPi0Boundary(Int_t i,Float_t param) 
{
  // Set parameter "Hard pi0 boundary" i to a value param
  if(i>1 || i<0) {
    AliError(Form("Invalid parameter number: %d",i));
  } else
    (*fParameters)(15,i) = param ;
}

//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterTimeGate(Int_t i, Float_t gate) 
{
  // Set the parameter TimeGate depending on Purity-Efficiency point i 
  if (i>2 || i<0) {
    AliError(Form("Invalid Efficiency-Purity choice %d",i));
  } else
    (*fParameters)(3,i)= gate ; 
} 

//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterToCalculateEllipse(TString particle, TString param, Int_t i, Float_t par) 
{  
  // Set the parameter "i" that is needed to calculate the ellipse 
  // parameter "param" for a particle "particle"
  
  particle.ToLower();
  param.   ToLower();
  Int_t p= -1;
  Int_t offset=0;

  if      (particle == "photon") offset=0;
  else if (particle == "pi0")    offset=5;
  else
    AliError(Form("Wrong particle name: %s (choose from pi0/photon)\n",
		  particle.Data()));

  if     (param.Contains("a")) p=4+offset; 
  else if(param.Contains("b")) p=5+offset; 
  else if(param.Contains("c")) p=6+offset; 
  else if(param.Contains("x0"))p=7+offset; 
  else if(param.Contains("y0"))p=8+offset;
  if((i>4)||(i<0)) {
    AliError(Form("No parameter with index %d", i)) ; 
  } else if(p==-1) {
    AliError(Form("No parameter with name %s", param.Data() )) ; 
  } else
    (*fParameters)(p,i) = par ;
} 

//____________________________________________________________________________
void AliPHOSPIDv1::Unload() 
{
  //Unloads RecPoints, Tracks and RecParticles
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;  
  gime->PhosLoader()->UnloadRecPoints() ;
  gime->PhosLoader()->UnloadTracks() ;
  gime->PhosLoader()->UnloadRecParticles() ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles()
{
  //It writes reconstructed particles and pid to file

  AliPHOSGetter *gime = AliPHOSGetter::Instance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 
  recParticles->Expand(recParticles->GetEntriesFast() ) ;
  if(fWrite){
    TTree * treeP =  gime->TreeP();
    
    //First rp
    Int_t bufferSize = 32000 ;
    TBranch * rpBranch = treeP->Branch("PHOSRP",&recParticles,bufferSize);
    rpBranch->SetTitle(BranchName());
    
    rpBranch->Fill() ;
    
    gime->WriteRecParticles("OVERWRITE");
    gime->WritePID("OVERWRITE");
  }
}


//_______________________________________________________________________
void AliPHOSPIDv1::SetInitPID(const Double_t *p) {
  // Sets values for the initial population of each particle type 
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) fInitPID[i] = p[i];
}
//_______________________________________________________________________
void AliPHOSPIDv1::GetInitPID(Double_t *p) const {
  // Gets values for the initial population of each particle type 
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) p[i] = fInitPID[i];
}
