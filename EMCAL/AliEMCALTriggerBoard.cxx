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

/*
 

EMCal trigger board super class
run the sliding window algorithm
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerBoard.h"
#include "AliEMCALTriggerPatch.h"

#include <TClonesArray.h>

#include <iostream>
#include <cstdlib>

using namespace std;

ClassImp(AliEMCALTriggerBoard)

//_______________
AliEMCALTriggerBoard::AliEMCALTriggerBoard() : TObject(),
fRegion(0x0),
fMap(0x0),
fRegionSize(0x0),
fSubRegionSize(0x0),
fPatchSize(0x0),
fPatches( new TClonesArray("AliEMCALTriggerPatch",10) )
{
	
}	

//_______________
AliEMCALTriggerBoard::AliEMCALTriggerBoard(const TVector2& RS) : TObject(),
fRegion(0x0),
fMap(0x0),
fRegionSize(    new TVector2( RS ) ),
fSubRegionSize( new TVector2() ),
fPatchSize(     new TVector2() ),
fPatches( new TClonesArray("AliEMCALTriggerPatch",10) )
{
   fRegion = (int**)malloc( (int)fRegionSize->X() * sizeof( int* ) );  

   if (!fRegion) printf("Error: malloc could not allocate %d bytes for fRegion\n",
                        int(fRegionSize->X() * sizeof( int* )));

      fMap = (int**)malloc( (int)fRegionSize->X() * sizeof( int* ) );

      if (!fMap) printf("Error: malloc could not allocate %d bytes for fMap\n",
                        int(fRegionSize->X() * sizeof( int* )));

   for (Int_t i=0;i<fRegionSize->X();i++)
   {
      fRegion[i] = (int*)malloc( (int)fRegionSize->Y() * sizeof( int ) );
      
      if (!fRegion[i]) printf("Error: malloc could not allocate %d bytes for fRegion[%d]\n",
                              i,int(fRegionSize->Y() * sizeof( int )));

         fMap[i] = (int*)malloc( (int)fRegionSize->Y() * sizeof( int ) );

            if (!fMap[i]) printf("Error: malloc could not allocate %d bytes for fMap[%d]\n",
                              i,int(fRegionSize->Y() * sizeof( int )));
   }

	// Initialize region matrix
	ZeroRegion();
	
	for (int i=0; i<fRegionSize->X(); ++i)
		for (int j=0; j<fRegionSize->Y(); ++j) fMap[i][j] = 0;
}

//_______________
AliEMCALTriggerBoard::~AliEMCALTriggerBoard()
{
   for (Int_t i=0;i<fRegionSize->X();i++) 
   {
      if (fRegion[i]) {delete fRegion[i]; fRegion[i] = 0;}
      if (   fMap[i]) {delete    fMap[i];    fMap[i] = 0;}
   }
   
   delete [] fRegion; fRegion = 0x0;
   delete []    fMap;    fMap = 0x0;
   
   if(fPatches)fPatches->Delete();
   delete fPatches;
}

//_______________
void AliEMCALTriggerBoard::ZeroRegion()
{
	//
	for (Int_t i=0;i<int(fRegionSize->X());i++) for (Int_t j=0;j<int(fRegionSize->Y());j++) fRegion[i][j] = 0;
}

//_______________
void AliEMCALTriggerBoard::SlidingWindow( L1TriggerType_t /*type*/, Int_t thres )
{
	//
	Int_t ipatch = 0;
	
	for (Int_t i=0; i<=int(fRegionSize->X()-fPatchSize->X()*fSubRegionSize->X()); i+=int(fSubRegionSize->X()))
	{
		for (Int_t j=0; j<=int(fRegionSize->Y()-fPatchSize->Y()*fSubRegionSize->Y()); j+=int(fSubRegionSize->Y()))
		{
			ipatch++;
			
			Int_t sum = 0;
			
			for (Int_t k=0; k<int(fPatchSize->X()*fSubRegionSize->X()); k++)
			{
				for (Int_t l=0; l<int(fPatchSize->Y()*fSubRegionSize->Y()); l++)
				{
					sum += fRegion[i+k][j+l];
				}
			}

			if ( sum > thres ) 
			{
				//if ( type == kJet ) sum /= 4; // truncate patch sum for jet case
				
				new((*fPatches)[fPatches->GetLast()+1]) 
						AliEMCALTriggerPatch( int(i/fSubRegionSize->X()), int(j/fSubRegionSize->Y()), int(sum) );
			}
		}
	}
}

//__________
void AliEMCALTriggerBoard::Scan()
{
	//
	cout << "     ";
	for (Int_t i=0; i<int(fRegionSize->X()); i++) printf("%8d ",i);
	cout << "\n";
	for (Int_t i=0; i<int(fRegionSize->X())-5; i++) printf("-------");
	cout << "\n";
	
	for (Int_t i=0; i<int(fRegionSize->Y()); i++)
	{
		if (i && !(i%12))
		{
			for (Int_t j=0; j<int(fRegionSize->X())-5; j++) printf("-------");
			cout << endl;
		}
		
		printf("%3d |",i);
		for (Int_t j=0; j<int(fRegionSize->X()); j++) 
		{
			printf("%2d/%5d ", fMap[j][i], fRegion[j][i]);
		}
		cout << endl;
	}
}

//__________
void AliEMCALTriggerBoard::Reset()
{
	//
	fPatches->Delete();
	ZeroRegion();
}

