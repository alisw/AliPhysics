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
 * Revision 1.116  2007/10/10 09:05:10  schutz
 * Changing name QualAss to QA
 *
 * Revision 1.115  2007/08/22 09:20:50  hristov
 * Updated QA classes (Yves)
 *
 * Revision 1.114  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.113  2007/07/18 16:29:54  policheh
 * Raw Sdigits energy converted to GeV.
 *
 * Revision 1.112  2007/02/25 22:59:13  policheh
 * Digits2Raw(): ALTRO buffer and mapping created per each DDL.
 *
 * Revision 1.111  2007/02/18 15:21:47  kharlov
 * Corrections for memory leak in Digits2Raw due to AliAltroMapping
 *
 * Revision 1.110  2007/02/13 10:52:08  policheh
 * Raw2SDigits() implemented
 *
 * Revision 1.109  2007/02/05 10:43:25  hristov
 * Changes for correct initialization of Geant4 (Mihaela)
 *
 * Revision 1.108  2007/02/01 10:34:47  hristov
 * Removing warnings on Solaris x86
 *
 * Revision 1.107  2007/01/29 16:29:37  kharlov
 * Digits2Raw(): special workaround for digits with time out of range
 *
 * Revision 1.106  2007/01/17 17:28:56  kharlov
 * Extract ALTRO sample generation to a separate class AliPHOSPulseGenerator
 *
 * Revision 1.105  2007/01/12 21:44:29  kharlov
 * Simulate and reconstruct two gains simulaneouslsy
 */

//_________________________________________________________________________
// Base Class for PHOS description:
//   PHOS consists of a PbWO4 calorimeter (EMCA) and a gazeous charged 
//    particles detector (CPV or PPSD).
//   The only provided method here is CreateMaterials, 
//    which defines the materials common to all PHOS versions.   
// 
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
class TFile;
#include <TF1.h> 
#include <TFolder.h> 
#include <TGeoGlobalMagField.h>
#include <TH1F.h> 
#include <TRandom.h> 
#include <TTree.h>
#include <TVirtualMC.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliMagF.h"
#include "AliPHOS.h"
#include "AliPHOSLoader.h"
#include "AliRun.h"
#include "AliRawReader.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSDigit.h"
#include "AliAltroBuffer.h"
#include "AliAltroMapping.h"
#include "AliCaloAltroMapping.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSPulseGenerator.h"
#include "AliDAQ.h"
#include "AliPHOSRawFitterv0.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSQAChecker.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSSimParam.h"

ClassImp(AliPHOS)

//____________________________________________________________________________
  AliPHOS:: AliPHOS() : AliDetector(),fgCalibData(0)
{
  // Default ctor
  fName   = "PHOS" ;

}

//____________________________________________________________________________
AliPHOS::AliPHOS(const char* name, const char* title): AliDetector(name, title),
fgCalibData(0)
{
  //   ctor : title is used to identify the layout
}

//____________________________________________________________________________
AliPHOS::~AliPHOS() 
{  
  if(fgCalibData) delete fgCalibData ;
}

//____________________________________________________________________________
AliDigitizer* AliPHOS::CreateDigitizer(AliDigitizationInput* digInput) const
{
  return new AliPHOSDigitizer(digInput);
}

//____________________________________________________________________________
void AliPHOS::CreateMaterials()
{
  // Definitions of materials to build PHOS and associated tracking media.
  // media number in idtmed are 699 to 798.

  // --- The PbWO4 crystals ---
  Float_t aX[3] = {207.19, 183.85, 16.0} ;
  Float_t zX[3] = {82.0, 74.0, 8.0} ;
  Float_t wX[3] = {1.0, 1.0, 4.0} ;
  Float_t dX = 8.28 ;

  AliMixture(0, "PbWO4$", aX, zX, dX, -3, wX) ;


  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;

  AliMixture(1, "Polystyrene$", aP, zP, dP, -2, wP) ;

  // --- Aluminium ---
  AliMaterial(2, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0) ;
  // ---         Absorption length is ignored ^

 // --- Tyvek (CnH2n) ---
  Float_t aT[2] = {12.011, 1.00794} ;
  Float_t zT[2] = {6.0, 1.0} ;
  Float_t wT[2] = {1.0, 2.0} ;
  Float_t dT = 0.331 ;

  AliMixture(3, "Tyvek$", aT, zT, dT, -2, wT) ;

  // --- Polystyrene foam ---
  Float_t aF[2] = {12.011, 1.00794} ;
  Float_t zF[2] = {6.0, 1.0} ;
  Float_t wF[2] = {1.0, 1.0} ;
  Float_t dF = 0.12 ;

  AliMixture(4, "Foam$", aF, zF, dF, -2, wF) ;

 // --- Titanium ---
  Float_t aTIT[3] = {47.88, 26.98, 54.94} ;
  Float_t zTIT[3] = {22.0, 13.0, 25.0} ;
  Float_t wTIT[3] = {69.0, 6.0, 1.0} ;
  Float_t dTIT = 4.5 ;

  AliMixture(5, "Titanium$", aTIT, zTIT, dTIT, -3, wTIT);

 // --- Silicon ---
  AliMaterial(6, "Si$", 28.0855, 14., 2.33, 9.36, 42.3, 0, 0) ;



  // --- Foam thermo insulation ---
  Float_t aTI[2] = {12.011, 1.00794} ;
  Float_t zTI[2] = {6.0, 1.0} ;
  Float_t wTI[2] = {1.0, 1.0} ;
  Float_t dTI = 0.04 ;

  AliMixture(7, "Thermo Insul.$", aTI, zTI, dTI, -2, wTI) ;

  // --- Textolith ---
  Float_t aTX[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zTX[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wTX[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dTX    = 1.75 ;

  AliMixture(8, "Textolit$", aTX, zTX, dTX, -4, wTX) ;

  //--- FR4  ---
  Float_t aFR[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zFR[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wFR[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dFR = 1.8 ; 

  AliMixture(9, "FR4$", aFR, zFR, dFR, -4, wFR) ;

  // --- The Composite Material for  micromegas (so far polyetylene) ---                                       
  Float_t aCM[2] = {12.01, 1.} ; 
  Float_t zCM[2] = {6., 1.} ; 
  Float_t wCM[2] = {1., 2.} ; 
  Float_t dCM = 0.935 ; 

  AliMixture(10, "Compo Mat$", aCM, zCM, dCM, -2, wCM) ;

  // --- Copper ---                                                                    
  AliMaterial(11, "Cu$", 63.546, 29, 8.96, 1.43, 14.8, 0, 0) ;
 
  // --- G10 : Printed Circuit material ---                                                  
  Float_t aG10[4] = { 12., 1., 16., 28.} ;
  Float_t zG10[4] = { 6., 1., 8., 14.} ;
  Float_t wG10[4] = { .259, .288, .248, .205} ;
  Float_t dG10  = 1.7 ; 

  
  AliMixture(12, "G10$", aG10, zG10, dG10, -4, wG10);

  // --- Lead ---                                                                     
  AliMaterial(13, "Pb$", 207.2, 82, 11.35, 0.56, 0., 0, 0) ;

 // --- The gas mixture ---                                                                
 // Co2
  Float_t aCO[2] = {12.0, 16.0} ; 
  Float_t zCO[2] = {6.0, 8.0} ; 
  Float_t wCO[2] = {1.0, 2.0} ; 
  Float_t dCO = 0.001977 ; 

  AliMixture(14, "CO2$", aCO, zCO, dCO, -2, wCO);

 // Ar
  Float_t dAr = 0.001782 ; 
  AliMaterial(15, "Ar$", 39.948, 18.0, dAr, 14.0, 0., 0, 0) ;   
 
  // Ar+CO2 Mixture (80% / 20%)
  Float_t arContent = 0.80 ;  // Ar-content of the ArCO2-mixture
  Float_t aArCO[3]  = {39.948, 12.0, 16.0} ;
  Float_t zArCO[3]  = {18.0  ,  6.0,  8.0} ;
  Float_t wArCO[3];
  wArCO[0] = arContent;
  wArCO[1] = (1-arContent)*1;
  wArCO[2] = (1-arContent)*2;
  Float_t dArCO = arContent*dAr + (1-arContent)*dCO ;
  AliMixture(16, "ArCO2$", aArCO, zArCO, dArCO,  -3, wArCO) ; 


  // --- Stainless steel (let it be pure iron) ---
  AliMaterial(17, "Steel$", 55.845, 26, 7.87, 1.76, 0., 0, 0) ;


  // --- Fiberglass ---
  Float_t aFG[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zFG[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wFG[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dFG    = 1.9 ;

  AliMixture(18, "Fibergla$", aFG, zFG, dFG, -4, wFG) ;

  // --- Cables in Air box  ---
  // SERVICES

  Float_t aCA[4] = { 1.,12.,55.8,63.5 };
  Float_t zCA[4] = { 1.,6.,26.,29. }; 
  Float_t wCA[4] = { .014,.086,.42,.48 };
  Float_t dCA    = 0.8 ;  //this density is raw estimation, if you know better - correct

  AliMixture(19, "Cables  $", aCA, zCA, dCA, -4, wCA) ;


  // --- Air ---
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
 
  AliMixture(99, "Air$", aAir, zAir, dAir, 4, wAir) ;

  // DEFINITION OF THE TRACKING MEDIA

  // for PHOS: idtmed[699->798] equivalent to fIdtmed[0->100]
  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ() ;
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max() ;

  // The scintillator of the calorimeter made of PBW04                              -> idtmed[699]
  AliMedium(0, "PHOS Xtal    $", 0, 1,
	    isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[700]
  AliMedium(1, "CPV scint.   $", 1, 1,
	    isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Various Aluminium parts made of Al                                             -> idtmed[701]
  AliMedium(2, "Al parts     $", 2, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // The Tywek which wraps the calorimeter crystals                                 -> idtmed[702]
  AliMedium(3, "Tyvek wrapper$", 3, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // The Polystyrene foam around the calorimeter module                             -> idtmed[703]
  AliMedium(4, "Polyst. foam $", 4, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The Titanium around the calorimeter crystal                                    -> idtmed[704]
  AliMedium(5, "Titan. cover $", 5, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.0001, 0.0001, 0, 0) ;

  // The Silicon of the pin diode to read out the calorimeter crystal               -> idtmed[705] 
 AliMedium(6, "Si PIN       $", 6, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.01, 0.01, 0, 0) ;

 // The thermo insulating material of the box which contains the calorimeter module -> idtmed[706]
  AliMedium(7, "Thermo Insul.$", 7, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The Textolit which makes up the box which contains the calorimeter module      -> idtmed[707]
  AliMedium(8, "Textolit     $", 8, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // FR4: The Plastic which makes up the frame of micromegas                        -> idtmed[708]
  AliMedium(9, "FR4 $", 9, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ; 


  // The Composite Material for  micromegas                                         -> idtmed[709]
  AliMedium(10, "CompoMat   $", 10, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Copper                                                                         -> idtmed[710]
  AliMedium(11, "Copper     $", 11, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ;

  // G10: Printed Circuit material                                                  -> idtmed[711]
 
  AliMedium(12, "G10        $", 12, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.01, 0, 0) ;

  // The Lead                                                                       -> idtmed[712]
 
  AliMedium(13, "Lead      $", 13, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The gas mixture: ArCo2                                                         -> idtmed[715]
 
  AliMedium(16, "ArCo2      $", 16, 1,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.01, 0, 0) ;
 
  // Stainless steel                                                                -> idtmed[716]
  AliMedium(17, "Steel     $", 17, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ;

  // Fibergalss                                                                     -> idtmed[717]
  AliMedium(18, "Fiberglass$", 18, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Cables in air                                                                  -> idtmed[718]
  AliMedium(19, "Cables    $", 19, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Air                                                                            -> idtmed[798] 
  AliMedium(99, "Air          $", 99, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;
}

//_____________________________________________________________________________
void AliPHOS::Init()
{
}

//____________________________________________________________________________
void AliPHOS::Digits2Raw()
{
// convert digits of the current event to raw data
  
  if(AliPHOSSimParam::GetInstance()->IsEDigitizationOn()){
    AliError("Energy digitization should be OFF if use Digits2Raw") ;
  }

  AliPHOSLoader * loader = static_cast<AliPHOSLoader*>(fLoader) ; 

  // get the digits
  loader->LoadDigits();
  TClonesArray* digits = loader->Digits() ;

  if (!digits) {
    AliError(Form("No digits found !"));
    return;
  }

  // get the geometry
  AliPHOSGeometry* geom = GetGeometry();
  if (!geom) {
    AliError(Form("No geometry found !"));
    return;
  }

  // get mapping from OCDB
  const TObjArray* maps = AliPHOSRecoParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  // some digitization constants
  const Float_t    kThreshold = 1.; // skip digits below 1 ADC channel
  const Int_t      kAdcThreshold = 1;  // Lower ADC threshold to write to raw data

  Int_t prevDDL = -1;

  if(fgCalibData==0)
    fgCalibData= new AliPHOSCalibData(-1) ;

  // Create a shaper pulse object
  AliPHOSPulseGenerator *pulse = new AliPHOSPulseGenerator();

  //Set Time step of sample
  pulse->SetTimeStep(TMath::Abs(fgCalibData->GetSampleTimeStep())) ;
  
  Int_t *adcValuesLow = new Int_t[pulse->GetRawFormatTimeBins()];
  Int_t *adcValuesHigh= new Int_t[pulse->GetRawFormatTimeBins()];
  
  const Int_t maxDDL = 20;
  AliAltroBuffer  *buffer[maxDDL];
  AliAltroMapping *mapping[maxDDL];

  for(Int_t jDDL=0; jDDL<maxDDL; jDDL++) {
    buffer[jDDL]=0;
    mapping[jDDL]=0;
  }

  //!!!!for debug!!!
  Int_t modMax=-111;
  Int_t colMax=-111;
  Int_t rowMax=-111;
  Float_t eMax=-333;
  //!!!for debug!!!

  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliPHOSDigit* digit = static_cast<AliPHOSDigit *>(digits->At(iDigit)) ;

    // Skip small energy below treshold
    if (digit->GetEnergy() < kThreshold) 
      continue;

    Int_t relId[4];
    geom->AbsToRelNumbering(digit->GetId(), relId);
    Int_t module = 5-relId[0];
 
    // Begin FIXME 
    if (relId[1] != 0) 
      continue;    // ignore digits from CPV
    // End FIXME 

    Int_t row = relId[2]-1;
    Int_t col = relId[3]-1;
    
    Int_t iRCU = -111;
    
    //RCU0
    if(0<=row&&row<16 && 0<=col&&col<56) iRCU=0;
    
    //RCU1
    if(16<=row&&row<32 && 0<=col&&col<56) iRCU=1;

    //RCU2
    if(32<=row&&row<48 && 0<=col&&col<56) iRCU=2;

    //RCU3
    if(48<=row&&row<64 && 0<=col&&col<56) iRCU=3;
    
    
    // PHOS EMCA has 4 DDL per module. Splitting is based on the (row,column) numbers.
    // here module already in PHOS online convention: 0<module<4 and inverse order.
    Int_t iDDL = 4 * module  + iRCU;

    // new DDL
    if (iDDL != prevDDL) {
      if (buffer[iDDL] == 0) {
	// open new file and write dummy header
	TString fileName = AliDAQ::DdlFileName("PHOS",iDDL);

	mapping[iDDL] = (AliAltroMapping*)maps->At(iDDL);
	buffer[iDDL]  = new AliAltroBuffer(fileName.Data(),mapping[iDDL]);
	buffer[iDDL]->WriteDataHeader(kTRUE, kFALSE);  //Dummy;
      }
      prevDDL = iDDL;
    }

    AliDebug(2,Form("digit E=%.4f GeV, t=%g s, (mod,col,row)=(%d,%d,%d)\n",
		    digit->GetEnergy(),digit->GetTimeR(),
		    relId[0]-1,relId[3]-1,relId[2]-1));
    // if a signal is out of time range, write only trailer
    if (digit->GetTimeR() > pulse->GetRawFormatTimeMax()*0.5 ) {
      AliDebug(2,"Signal is out of time range.\n");
      buffer[iDDL]->FillBuffer(0);
      buffer[iDDL]->FillBuffer(pulse->GetRawFormatTimeBins() );  // time bin
      buffer[iDDL]->FillBuffer(3);                               // bunch length
      buffer[iDDL]->WriteTrailer(3, relId[3]-1, relId[2]-1, 0);  // trailer
      
    // calculate the time response function
    } else {
      Double_t energy = 0 ;
      if (digit->GetId() <= geom->GetNModules() * geom->GetNCristalsInModule()) {
	energy=digit->GetEnergy();
	if(energy>eMax) {eMax=energy; modMax=relId[0]; colMax=col; rowMax=row;}
      }
      else {
 	energy = 0; // CPV raw data format is now know yet
      }
      pulse->SetAmplitude(energy);
      pulse->SetTZero(digit->GetTimeR());
      Double_t r =fgCalibData->GetHighLowRatioEmc(relId[0],relId[3],relId[2]) ;
      pulse->SetHG2LGRatio(r) ;
      pulse->MakeSamples();
      pulse->GetSamples(adcValuesHigh, adcValuesLow) ; 

      buffer[iDDL]->WriteChannel(relId[3]-1, relId[2]-1, 0, 
			   pulse->GetRawFormatTimeBins(), adcValuesLow , kAdcThreshold);
      buffer[iDDL]->WriteChannel(relId[3]-1, relId[2]-1, 1, 
			   pulse->GetRawFormatTimeBins(), adcValuesHigh, kAdcThreshold);
    }
  }
  delete [] adcValuesLow;
  delete [] adcValuesHigh;

  // write real header and close last file
  for (Int_t iDDL=0; iDDL<maxDDL; iDDL++) {
    if (buffer[iDDL]) {
      buffer[iDDL]->Flush();
      buffer[iDDL]->WriteDataHeader(kFALSE, kFALSE);
      delete buffer[iDDL];
      //if (mapping[iDDL]) delete mapping[iDDL];
    }
  }
  
  AliDebug(1,Form("Digit with max. energy:  modMax %d colMax %d rowMax %d  eMax %f\n",
	 modMax,colMax,rowMax,eMax));

  delete pulse;
  loader->UnloadDigits();
}

//____________________________________________________________________________
void AliPHOS::Hits2SDigits()  
{ 
// create summable digits

  AliPHOSSDigitizer phosDigitizer(fLoader->GetRunLoader()->GetFileName().Data()) ;
  phosDigitizer.SetEventRange(0, -1) ; // do all the events
 
  phosDigitizer.Digitize("all") ; 
}


//____________________________________________________________________________
AliLoader* AliPHOS::MakeLoader(const char* topfoldername)
{
//different behaviour than standard (singleton getter)
// --> to be discussed and made eventually coherent
 fLoader = new AliPHOSLoader(GetName(),topfoldername);
 return fLoader;
}

//____________________________________________________________________________
void AliPHOS::SetTreeAddress()
{ 
  // Links Hits in the Tree to Hits array
  TBranch *branch;
  char branchname[20];
  snprintf(branchname,20,"%s",GetName());
  // Branch address for hit tree
    TTree *treeH = fLoader->TreeH();
  if (treeH) {
    branch = treeH->GetBranch(branchname);
    if (branch) 
     { 
       if (fHits == 0x0) fHits= new TClonesArray("AliPHOSHit",1000);
       //AliInfo(Form("<%s> Setting Hits Address",GetName()));
       branch->SetAddress(&fHits);
     }
  }
}

//____________________________________________________________________________ 	 
Bool_t AliPHOS::Raw2SDigits(AliRawReader* rawReader) 	 
{ 	 
	  	 
  AliPHOSLoader * loader = static_cast<AliPHOSLoader*>(fLoader) ; 	 
	  	 
  TTree * tree = 0 ; 	 
  tree = loader->TreeS() ; 	 
  if ( !tree ) { 	 
    loader->MakeTree("S"); 	 
    tree = loader->TreeS() ; 	 
  } 	 
	  	 
  TClonesArray * sdigits = loader->SDigits() ; 	 
  if(!sdigits) { 	 
    loader->MakeSDigitsArray(); 	 
    sdigits = loader->SDigits(); 	 
  } 	 
  sdigits->Clear(); 	 
	  	 
  rawReader->Reset() ;

  const TObjArray* maps = AliPHOSRecoParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  AliAltroMapping *mapping[20];
  for(Int_t i = 0; i < 20; i++) {
    mapping[i] = (AliAltroMapping*)maps->At(i);
  }

  AliPHOSRawFitterv0 fitter;

  fitter.SubtractPedestals(AliPHOSSimParam::GetInstance()->EMCSubtractPedestals());
  fitter.SetAmpOffset(AliPHOSSimParam::GetInstance()->GetGlobalAltroOffset());
  fitter.SetAmpThreshold(AliPHOSSimParam::GetInstance()->GetGlobalAltroThreshold());

  AliPHOSRawDigiProducer pr(rawReader,mapping);
  pr.SetSampleQualityCut(AliPHOSSimParam::GetInstance()->GetEMCSampleQualityCut()); 	 
  pr.MakeDigits(sdigits,&fitter);
	  	 
  Int_t bufferSize = 32000 ; 	 
  // TBranch * sdigitsBranch = tree->Branch("PHOS",&sdigits,bufferSize); 	 
  tree->Branch("PHOS",&sdigits,bufferSize); 	 
  tree->Fill(); 	 
	  	 
  fLoader->WriteSDigits("OVERWRITE"); 	 
  return kTRUE; 	 
  
}
