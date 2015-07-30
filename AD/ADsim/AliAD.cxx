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

/* $Id: AliAD.cxx  $ */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                  AD (ALICE Diffractive)  Detector                     //
//                                                                       //
//  This class contains the base procedures for the AD  detector         //
//  All comments should be sent to :                                     //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


// --- Standard libraries ---
#include <Riostream.h>
#include <stdlib.h>

// --- ROOT libraries ---
#include <TNamed.h>
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "TStopwatch.h"
#include "TParameter.h"
#include "TF1.h"

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliMC.h"
#include "AliAD.h"
#include "AliADhit.h"
#include "AliADLoader.h"
#include "AliADDigitizer.h"
#include "AliDigitizationInput.h"
#include "AliADdigit.h"
#include "AliADSDigit.h"
#include "AliADBuffer.h"
#include "AliDAQ.h"
#include "AliRawReader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliADReconstructor.h"
#include "AliADCalibData.h"
#include "AliADConst.h"

ClassImp(AliAD)
 //__________________________________________________________________
AliAD::AliAD()
   : AliDetector(),
     fCalibData(NULL),
     fSetADATwoInstalled(0),
     fSetADCTwoInstalled(0)
{
	/// Default Constructor
    

}

//_____________________________________________________________________________
AliAD::AliAD(const char *name, const char *title)
   : AliDetector(name,title),
     fCalibData(NULL),
     fSetADATwoInstalled(kTRUE),
     fSetADCTwoInstalled(kTRUE)

{
  
   // Standard constructor for AD Detector
  

}

//_____________________________________________________________________________
AliAD::~AliAD()
{
   //
   // Default destructor for AD Detector
   //
  
}

//_____________________________________________________________________________
void AliAD::CreateMaterials()
{
  //
  // MATERIALS FOR ADC AND ADA
  //
  
  // Parameters for simulation scope for ADA and ADC (stolen from AliVZEROv7::CreateMaterials )
  Int_t    fieldType       = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();     // Field type 
  Double_t maxField        = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();       // Field max.
  // Int_t   isxfld1 = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  // Float_t sxmgmx  = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  Int_t    isxfld2         = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->PrecInteg();
  Double_t maxBending      = 10;    // Max Angle
  Double_t maxStepSize     = 0.01;  // Max step size 
  Double_t maxEnergyLoss   = 1;     // Max Delta E
  Double_t precision       = 0.003; // Precision
  Double_t minStepSize     = 0.003; // Minimum step size 
  Float_t  density,  as[11], zs[11], ws[11];
  Double_t radLength, absLength, a_ad, z_ad;
  Int_t    id;
  
  //
  // Air
  //
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;
  // Steel  
  Float_t asteel[4] = { 55.847,51.9961, 58.6934, 28.0855 };
  Float_t zsteel[4] = {    26.,    24.,     28., 14.     };
  Float_t wsteel[4] = {   .715,    .18,      .1, .005    };
  //
  // Cast iron
  //
  Float_t acasti[4] = {55.847,12.011,28.085,54.938};
  Float_t zcasti[4] = {26.,6.,14.,25.};
  Float_t wcasti[4] = {0.929,0.035,0.031,0.005};
  // --- Define the various materials for GEANT --- 
  //     Aluminum 
  AliMaterial(9,  "ALU1      ", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALU2      ", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALU3      ", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "IRON1     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON2     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON3     ", 55.85, 26., 7.87, 1.76, 17.1);

  //
  //     Copper
  AliMaterial(11, "COPPER1   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(31, "COPPER2   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(51, "COPPER3   ", 63.55, 29., 8.96, 1.43, 15.1);
  
  //     Vacuum 
  AliMixture(16, "VACUUM1 ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(36, "VACUUM2 ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(56, "VACUUM3 ", aAir, zAir, dAir1, 4, wAir);
  
  //     Stainless Steel 
  AliMixture(19, "STAINLESS STEEL1", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(39, "STAINLESS STEEL2", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(59, "STAINLESS STEEL3", asteel, zsteel, 7.88, 4, wsteel);
  
  //     Cast iron
  AliMixture(18, "CAST IRON1", acasti, zcasti, 7.2, 4, wcasti);
  AliMixture(38, "CAST IRON2", acasti, zcasti, 7.2, 4, wcasti);
  AliMixture(58, "CAST IRON3", acasti, zcasti, 7.2, 4, wcasti);
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  
  //    Aluminum 
  AliMedium(9,  "ALU_C0          ", 9, 0,  fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1          ", 29, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(49, "ALU_C2          ", 49, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Iron 
  AliMedium(10, "FE_C0           ", 10, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1           ", 30, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2           ", 50, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);

  //    Copper 
  AliMedium(11, "Cu_C0           ", 11, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(31, "Cu_C1           ", 31, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(51, "Cu_C2           ", 51, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  AliMedium(16, "VA_C0           ", 16, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "VA_C1           ", 36, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(56, "VA_C2           ", 56, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  AliMedium(19, "ST_C0           ", 19, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(39, "ST_C1           ", 39, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(59, "ST_C3           ", 59, 0, fieldType, maxField, tmaxfd, stemax, deemax, epsil, stmin);

   //
   // Parameters  for AD scintillator: NE-102 (BC400)
   //
   // NE-102, has the following properties : (from internet, need better reference)
   //    Density : ca. 1.03 g/cm3
   //    Electrons/cm3: 3.39 x 10^23
   //    H atoms/cm3: 5.28 x 10^22
   //    C atoms/cm3: 4.78 x 10^22
   //    Ratio of H to C : 1.104 .
   //    wavelength of emission : ~4.23 nm.
   //    Decay time : ~2.4 ns.
   //    Luminescent efficiency : typically 18% of NaI(Tl)
   //    Photons/MeV: 2.5 x 10^4 
   //
   // H                // C 
   // as[0] = 1.00794;    as[1] = 12.011;
   // zs[0] = 1.;         zs[1] = 6.;
   // ws[0] = 5.23;       ws[1] = 4.74;
   // density = 1.032;
   // id      = 1;
   // AliMixture( id, "NE102", as, zs, density, -2, ws );
   // AliMedium( id, "NE102", id, 1, fieldType, maxField, maxBending, maxStepSize,
   //            maxEnergyLoss, precision, minStepSize );

   // ecalvovi@cern.ch
   // Parameters  for AD scintillator: BC404
   //
   // NE-102, has the following properties : (from internet, need better reference)
   //    Density : ca. 1.032 g/cm3
   //    Electrons/cm3: 3.37 x 10^23
   //    H atoms/cm3: 5.21 x 10^22
   //    C atoms/cm3: 4.74 x 10^22
   //    Ratio of H to C : 1.100 
   //    wavelength of emission : 408 nm.
   //    Decay time : 1.8 ns.
   //    Luminescent efficiency : typically 18% of NaI(Tl)
   //    Photons/MeV: ??
   //
   // H                // C 
   as[0] = 1.00794;    as[1] = 12.011;
   zs[0] = 1.;         zs[1] = 6.;
   ws[0] = 5.21;       ws[1] = 4.74;
   density = 1.032;
   id      = 1;
   AliMixture( id, "BC404", as, zs, density, -2, ws );
   AliMedium ( id, "BC404", id, 1, fieldType, maxField, maxBending, maxStepSize,
              maxEnergyLoss, precision, minStepSize );
   // parameters AliMedium: numed  name   nmat   isvol  ifield fieldm tmaxfd stemax deemax epsil  stmin  
   // ... 
   // isvol       sensitive volume if isvol!=0
   // ifield      magnetic field flag (see below)
   // fieldm      maximum magnetic field
   // ...
   // ifield =  0       no magnetic field
   //        = -1       user decision in guswim
   //        =  1       tracking performed with Runge Kutta
   //        =  2       tracking performed with helix
   //        =  3       constant magnetic field along z


   //
   // Parameters for lightGuide:  
   //     TODO check material 
   // Should be Poly(methyl methacrylate) (PMMA) acrylic 
   // (C5O2H8)n 
   // Density  1.18 g/cm3
   // Mixture PMMA    Aeff=12.3994 Zeff=6.23653 rho=1.18 radlen=34.0677 intlen=63.3073
   // Element #0 : C  Z=  6.00 A= 12.01 w= 0.600 natoms=5
   // Element #1 : H  Z=  1.00 A=  1.01 w= 0.081 natoms=8
   // Element #2 : O  Z=  8.00 A= 16.00 w= 0.320 natoms=2

   // Carbon          Hydrogen          Oxygen
   as[0] = 12.0107;   as[1] = 1.00794;  as[2] = 15.9994;
   zs[0] = 6.;        zs[1] = 1.;       zs[2] = 8.;
   ws[0] = 0.60;      ws[1] = 0.081;    ws[2] = 0.32;
   density = 1.18;
   id      = 2;
   AliMixture( id, "PMMA", as, zs, density, 3, ws );
   AliMedium( id,"PMMA", id, 1, fieldType, maxField, maxBending, maxStepSize,
              maxEnergyLoss, precision, minStepSize );

   
   // mu-metal
   // Niquel          Iron              Molybdenum        Manganese
   as[0] = 58.6934;   as[1] = 55.845;   as[2] = 95.94;    as[3] = 54.9380;  
   zs[0] = 28.;       zs[1] = 26.;      zs[2] = 42.;      zs[3] = 25.;   
   ws[0] = 0.802;     ws[1] = 0.14079;  ws[2] = 0.0485;   ws[3] = 0.005;
   // Silicon         Chromium          Cobalt            Aluminium
   as[4] = 28.0855;   as[5] = 51.9961;  as[6] = 58.9332;  as[7] = 26.981539;   
   zs[4] = 14.;       zs[5] = 24.;      zs[6] = 27.;      zs[7] = 13.;   
   ws[4] = 0.003;     ws[5] = 0.0002;   ws[6] = 0.0002;   ws[7] = 0.0001;
   // Carbon          Phosphorus        Sulfur
   as[8] = 12.0107;   as[9] = 30.97376; as[10] = 32.066;
   zs[8] = 6.;        zs[9] = 15.;      zs[10] = 16.;
   ws[8] = 0.00015;   ws[9] = 0.00005;  ws[10] = 0.00001;
   density = 8.25;
   id      = 3;
   AliMixture( id, "MuMetal", as, zs, density, 11, ws );
   AliMedium( id,"MuMetal", id, 1, fieldType, maxField, maxBending, maxStepSize,
              maxEnergyLoss, precision, minStepSize );

   // Parameters for ADCPMA: Aluminium
   a_ad = 26.98; 
   z_ad = 13.00;
   density   = 2.7;
   radLength = 8.9;
   absLength = 37.2;
   id = 4;
   AliMaterial (id, "Alum",  a_ad, z_ad, density, radLength, absLength, 0, 0 );
   AliMedium( id, "Alum", id, 1, fieldType, maxField, maxBending, maxStepSize,
              maxEnergyLoss, precision, minStepSize );

   // Parameters for ADCPMG: Glass for the simulation Aluminium 
   // TODO fix material
   a_ad = 26.98; 
   z_ad = 13.00;
   density   = 2.7;
   radLength = 8.9;
   absLength = 37.2;
   id = 5;
   AliMaterial( id, "Glass",  a_ad, z_ad, density, radLength, absLength, 0, 0 );
   AliMedium( id, "Glass", id, 1, fieldType, maxField, maxBending, maxStepSize,
              maxEnergyLoss, precision, minStepSize );


}
//_____________________________________________________________________________
void AliAD::SetTreeAddress()
{
   //
   // Sets tree address for hits.
   //

	TBranch *branch;
  	char branchname[20];
  	snprintf(branchname,19,"%s",GetName());
  	// Branch address for hit tree
  	TTree *treeH = fLoader->TreeH();
  	if (treeH ) 
	{
    		branch = treeH->GetBranch(branchname);
    		if (branch) branch->SetAddress(&fHits);
  	}
}


//_____________________________________________________________________________
void AliAD::MakeBranch(Option_t* opt)
{
	const char* oH = strstr(opt,"H");
	if (fLoader->TreeH() && oH && (fHits==0x0))
	{
		fHits = new TClonesArray("AliADhit",1000);
		fNhits = 0;
	}
	AliDetector::MakeBranch(opt);
}
//_____________________________________________________________________________
AliLoader* AliAD::MakeLoader(const char* topfoldername)
{ 
 
   AliDebug(1,Form("Creating AliADLoader, Top folder is %s ",topfoldername));
   fLoader = new AliADLoader(GetName(),topfoldername);
   return fLoader;
}

//_____________________________________________________________________________
AliDigitizer* AliAD::CreateDigitizer(AliDigitizationInput* digInput) const
{
   //
   // Creates a digitizer for AD
   //
   return new AliADDigitizer(digInput);
}
//_____________________________________________________________________________
void AliAD::Hits2Digits(){
  //
  // Converts hits to digits
  //
  // Creates the AD digitizer 
  AliADDigitizer* dig = new AliADDigitizer(this,AliADDigitizer::kHits2Digits);

  // Creates the digits
  dig->Digitize("");

  // deletes the digitizer
  delete dig;
}

//_____________________________________________________________________________
void AliAD::Hits2SDigits(){
  //
  // Converts hits to summable digits
  //
  // Creates the AD digitizer 
  AliADDigitizer* dig = new AliADDigitizer(this,AliADDigitizer::kHits2SDigits);

  // Creates the sdigits
  dig->Digitize("");

  // deletes the digitizer
  delete dig;
}

//_____________________________________________________________________________

void AliAD::Digits2Raw()
{
   //
   //  Converts digits of the current event to raw data
   //
   GetCalibData();
   
   AliAD *fAD = (AliAD*)gAlice->GetDetector("AD");
   fLoader->LoadDigits();
   TTree* digits = fLoader->TreeD();
   if (!digits) {
      Error("Digits2Raw", "no digits tree");
      return;
   }
   TClonesArray * ADdigits = new TClonesArray("AliADdigit",1000);
   fAD->SetTreeAddress();  		
   digits->GetBranch("ADDigit")->SetAddress(&ADdigits); 
  
   const char *fileName    = AliDAQ::DdlFileName("AD",0);
   AliADBuffer* buffer  = new AliADBuffer(fileName);
   
   // Get Trigger information first
   // Read trigger inputs from trigger-detector object
   AliDataLoader * dataLoader = fLoader->GetDigitsDataLoader();
   if( !dataLoader->IsFileOpen() ) 
        dataLoader->OpenFile("READ");
   AliTriggerDetector* trgdet = (AliTriggerDetector*)dataLoader->GetDirectory()->Get("Trigger");
   UInt_t triggerInfo = 0;
   if(trgdet) {
      triggerInfo = trgdet->GetMask() & 0xffff;
   }
   else {
      AliError(Form("There is no trigger object for %s",fLoader->GetName()));
   }

   buffer->WriteTriggerInfo((UInt_t)triggerInfo); 
   buffer->WriteTriggerScalers(); 
   buffer->WriteBunchNumbers(); 
  
   // Now retrieve the channel information: charge smaples + time and 
   // dump it into ADC and Time arrays
   Int_t nEntries = Int_t(digits->GetEntries());
   Short_t aADC[16][kADNClocks];
   Short_t aTime[16];
   Short_t aWidth[16];
   Bool_t  aIntegrator[16];
   Bool_t  aBBflag[16];
   Bool_t  aBGflag[16];
  
   for (Int_t i = 0; i < nEntries; i++) {
     fAD->ResetDigits();
     digits->GetEvent(i);
     Int_t ndig = ADdigits->GetEntriesFast(); 
   
     if(ndig == 0) continue;
     for(Int_t k=0; k<ndig; k++){
         AliADdigit* fADDigit = (AliADdigit*) ADdigits->At(k);

	 Int_t iChannel       = kOnlineChannel[fADDigit->PMNumber()];
	 
	 for(Int_t iClock = 0; iClock < kADNClocks; ++iClock) aADC[iChannel][iClock] = fADDigit->ChargeADC(20-iClock);
	 Int_t board = AliADCalibData::GetBoardNumber(iChannel);
	 aTime[iChannel]      = TMath::Nint(fADDigit->Time() / fCalibData->GetTimeResolution(board));
	 aWidth[iChannel]     = TMath::Nint(fADDigit->Width() / fCalibData->GetWidthResolution(board));
	 aIntegrator[iChannel]= fADDigit->Integrator();
	 aBBflag[iChannel]    = fADDigit->GetBBflag();
	 aBGflag[iChannel]    = fADDigit->GetBGflag();
	 
         //AliDebug(1,Form("DDL: %s\tdigit number: %d\tPM number: %d\tADC: %d\tTime: %f",
			 //fileName,k,fADDigit->PMNumber(),aADC[iChannel][AliADdigit::kADNClocks/2],aTime[iChannel])); 
     }        
   }

   // Now fill raw data	
   Int_t iCIU=0;
   for (Int_t  iV0CIU = 0; iV0CIU < 8; iV0CIU++) {
   
      if(iV0CIU != 2 && iV0CIU != 5) {
       	buffer->WriteEmptyCIU();
        continue;
       }
      // decoding of one Channel Interface Unit numbered iCIU - there are 8 channels per CIU (and 2 CIUs) : 
      for(Int_t iChannel_Offset = iCIU*8; iChannel_Offset < (iCIU*8)+8; iChannel_Offset=iChannel_Offset+4) { 
         for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {
             buffer->WriteChannel(iChannel, aADC[iChannel], aIntegrator[iChannel]);       
         }
         buffer->WriteBeamFlags(&aBBflag[iChannel_Offset],&aBGflag[iChannel_Offset]); 
         buffer->WriteMBInfo(); 
         buffer->WriteMBFlags();   
         buffer->WriteBeamScalers(); 
      } 

      for(Int_t iChannel = iCIU*8 + 7; iChannel >= iCIU*8; iChannel--) {
          buffer->WriteTiming(aTime[iChannel], aWidth[iChannel]); 
      } // End of decoding of one CIU card
      iCIU++;    
  } // end of decoding the eight CIUs

     
  delete buffer;
  fLoader->UnloadDigits();  
}

//_____________________________________________________________________________

Bool_t AliAD::Raw2SDigits(AliRawReader* rawReader)
{
	// reads raw data to produce sdigits
	// for AD not implemented yet 
	return kTRUE;
}

//_____________________________________________________________________________
void AliAD::GetCalibData()
{
  // Gets calibration object for AD set
  // Do nothing in case it is already loaded
  if (fCalibData) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("AD/Calib/Data");
  if (entry) fCalibData = (AliADCalibData*) entry->GetObject();
  if (!fCalibData)  AliFatal("No calibration data from calibration database !");
  
  return;
}
