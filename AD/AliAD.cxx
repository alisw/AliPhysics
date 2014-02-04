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
#include "AliDAQ.h"
#include "AliRawReader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliADReconstructor.h"

ClassImp(AliAD)
 //__________________________________________________________________
AliAD::AliAD()
   : AliDetector(),
     fSetADAToInstalled(0),
     fSetADCToInstalled(0)
{
	/// Default Constructor
    

}

//_____________________________________________________________________________
AliAD::AliAD(const char *name, const char *title)
   : AliDetector(name,title),
     fSetADAToInstalled(kTRUE),
     fSetADCToInstalled(kTRUE)

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
   Double_t maxBending      = 10;    // Max Angle
   Double_t maxStepSize     = 0.01;  // Max step size 
   Double_t maxEnergyLoss   = 1;     // Max Delta E
   Double_t precision       = 0.003; // Precision
   Double_t minStepSize     = 0.003; // Minimum step size 
   Float_t  density,  as[11], zs[11], ws[11];
   Double_t radLength, absLength, a_ad, z_ad;
   Int_t    id;
   
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
   as[0] = 1.00794;    as[1] = 12.011;
   zs[0] = 1.;         zs[1] = 6.;
   ws[0] = 5.23;       ws[1] = 4.74;
   density = 1.032;
   id      = 1;
   AliMixture( id, "NE102", as, zs, density, -2, ws );
   AliMedium( id, "NE102", id, 1, fieldType, maxField, maxBending, maxStepSize,
              maxEnergyLoss, precision, minStepSize );

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

void AliAD::Digits2Raw()
{
	// produces raw data from digits
	// for AD not implemented yet (needs detailed hardware info)
}

//_____________________________________________________________________________

Bool_t AliAD::Raw2SDigits(AliRawReader* rawReader)
{
	// reads raw data to produce digits
	// for AD not implemented yet (needs detailed hardware info)
	return kTRUE;
}
