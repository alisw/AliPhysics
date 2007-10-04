//
// TOFDigitsInfo
//

#include <Reve/TTreeTools.h>

#include "TOFDigitsInfo.h"
#include <AliTOFdigit.h>
#include <AliTOFGeometry.h>
//#include <AliTOFDigitMap.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

//_________________________________________________________

ClassImp(TOFDigitsInfo)

  TOFDigitsInfo::TOFDigitsInfo(): 
    TObject(),
    ReferenceCount(),
    fTree (0),
    fNewTree (0),
    fGeom (0),
    fTOFdigitMap(new AliTOFDigitMap())
{}
/* ******************************************************* */

TOFDigitsInfo:: ~TOFDigitsInfo() 
{

  delete fGeom;
  delete fTree;
  delete fNewTree;
  delete fTOFdigitMap;

}
/* ******************************************************* */

void TOFDigitsInfo::SetTree(TTree* tree)
{
  static const Exc_t eH("TOFDigitsInfo::SetTree ");
  
  if(fGeom == 0) {
    fGeom = new AliTOFGeometry();
  }
  
  fTree = tree;
  /*
  DecRefCount();
  IncRefCount();
  */
}
/* ******************************************************* */

void TOFDigitsInfo::LoadDigits()
{

  TClonesArray *digitsTOF = 0x0;
  AliTOFdigit *digs;

  fTree->SetBranchAddress("TOF",&digitsTOF);
  fTree->GetEntry(0);

  Int_t vol[5] = {-1,-1,-1,-1,-1};

  for (Int_t digitNumber=0; digitNumber<digitsTOF->GetEntries(); digitNumber++) {

    //if (digitNumber==digitsTOF->GetEntries()-1) printf(" Hello  4 -> %3i digit of %i \n", digitNumber+1, digitsTOF->GetEntries());
  
    digs = (AliTOFdigit*)digitsTOF->UncheckedAt(digitNumber);

    vol[0] = digs->GetSector(); // Sector Number (0-17)
    vol[1] = digs->GetPlate();  // Plate Number (0-4)
    vol[2] = digs->GetStrip();  // Strip Number (0-14/18)
    vol[3] = digs->GetPadx();   // Pad Number in x direction (0-47)
    vol[4] = digs->GetPadz();   // Pad Number in z direction (0-1)

    fTOFdigitMap->AddDigit(vol, digitNumber);
    //if (digitNumber==digitsTOF->GetEntries()-1) printf(" I am inside LoadDigits %3i \n", digitNumber);

  }

}

/* ******************************************************* */

void TOFDigitsInfo::GetDigits(Int_t nSector, Int_t nPlate,
			      Int_t nStrip, Int_t nPadZ, Int_t nPadX,
			      Int_t indexDigit[3])
{

  Int_t vol[5] = {nSector,nPlate,nStrip,nPadX,nPadZ};

  fTOFdigitMap->GetDigitIndex(vol, indexDigit);

}
/* ******************************************************* */

TClonesArray* TOFDigitsInfo::GetDigits(Int_t nSector, Int_t nPlate,
				       Int_t nStrip)
{

  Int_t newCounter = 0;
  Int_t nDigitsInVolume[3] = {-1, -1, -1};
  Int_t dummy[3] = {-1, -1, -1};
  Int_t informations[4] = {-1, -1, -1, -1};

  TClonesArray* digitsTOFnew = new TClonesArray("AliTOFdigit",  300);
  TClonesArray &ldigits = *digitsTOFnew;

  AliTOFdigit *digs;

  TClonesArray *digitsTOF = 0x0;
  fTree->SetBranchAddress("TOF",&digitsTOF);
  fTree->GetEntry(0);


  Int_t vol[5] = {nSector,nPlate,nStrip,-1,-1};

  for(Int_t iPadZ=0; iPadZ<fGeom->NpadZ(); iPadZ++){
    vol[4] = iPadZ;
    for(Int_t iPadX=0; iPadX<fGeom->NpadX(); iPadX++) {
      vol[3] = iPadX;

      //GetDigits(vol[0], vol[1], vol[2], vol[3], vol[4], nDigitsInVolume)

      fTOFdigitMap->GetDigitIndex(vol, nDigitsInVolume);

      for (Int_t ii=0; ii<3; ii++) {

	if (nDigitsInVolume[ii]>=0 ) {
	  //printf("  nDigitsInVolume[%2i]  = %3i\n ", ii, nDigitsInVolume[ii]);
	  digs = (AliTOFdigit*)digitsTOF->UncheckedAt(nDigitsInVolume[ii]);
	  informations[0] = digs->GetTdc();
	  informations[1] = digs->GetAdc();
	  informations[2] = digs->GetToT();
	  informations[3] = digs->GetTdcND();
	  new (ldigits[newCounter++]) AliTOFdigit(dummy, vol, informations);
	}

      }

      for (Int_t ii=0; ii<4; ii++) informations[ii]=-1;
      for (Int_t ii=0; ii<3; ii++) dummy[ii]=-1;
      for (Int_t ii=0; ii<3; ii++) nDigitsInVolume[ii]=-1;

    }
  }

  /*
  if (digitsTOFnew)
    printf("Sector %2i   Plate %1i  Strip %2i  -> number of digits %3i \n",
	   nSector, nPlate, nStrip, digitsTOFnew->GetEntries());
  */
  return digitsTOFnew;

}
/* ******************************************************* */

TClonesArray* TOFDigitsInfo::GetDigits(Int_t nSector)
{

  Int_t newCounter = 0;
  Int_t nDigitsInVolume[3] = {-1, -1, -1};
  Int_t dummy[3] = {-1, -1, -1};
  Int_t informations[4] = {-1, -1, -1, -1};

  Int_t nStrips=19;

  TClonesArray* digitsTOFnew = new TClonesArray("AliTOFdigit",  300);
  TClonesArray &ldigits = *digitsTOFnew;

  AliTOFdigit *digs;

  TClonesArray *digitsTOF = 0x0;
  fTree->SetBranchAddress("TOF",&digitsTOF);
  fTree->GetEntry(0);

  //Int_t nSector = 1;
  Int_t vol[5] = {nSector,-1,-1,-1,-1};
 
  for(Int_t iPlate=0; iPlate<fGeom->NPlates(); iPlate++){
    vol[1] = iPlate;
    if(iPlate==2) nStrips=15;
    else nStrips=19;
      
    for(Int_t iStrip=0; iStrip<nStrips; iStrip++){
      vol[2] = iStrip;
	
      for(Int_t iPadZ=0; iPadZ<fGeom->NpadZ(); iPadZ++){
	vol[4] = iPadZ;

	for(Int_t iPadX=0; iPadX<fGeom->NpadX(); iPadX++) {
	  vol[3] = iPadX;

	  //GetDigits(vol[0], vol[1], vol[2], vol[3], vol[4], nDigitsInVolume)

	  fTOFdigitMap->GetDigitIndex(vol, nDigitsInVolume);

	  for (Int_t ii=0; ii<3; ii++) {

	    if (nDigitsInVolume[ii]>=0 ) {
	      //printf("  nDigitsInVolume[%2i]  = %3i\n ", ii, nDigitsInVolume[ii]);
	      digs = (AliTOFdigit*)digitsTOF->UncheckedAt(nDigitsInVolume[ii]);
	      informations[0] = digs->GetTdc();
	      informations[1] = digs->GetAdc();
	      informations[2] = digs->GetToT();
	      informations[3] = digs->GetTdcND();
	      new (ldigits[newCounter++]) AliTOFdigit(dummy, vol, informations);
	    }

	  }

	  for (Int_t ii=0; ii<4; ii++) informations[ii]=-1;
	  for (Int_t ii=0; ii<3; ii++) dummy[ii]=-1;
	  for (Int_t ii=0; ii<3; ii++) nDigitsInVolume[ii]=-1;
	    
	}
      }
    }
  }

  /*
  if (digitsTOFnew)
    printf("Sector %2i   Plate %1i  Strip %2i  -> number of digits %3i \n",
	   nSector, nPlate, nStrip, digitsTOFnew->GetEntries());
  */
  return digitsTOFnew;

}
/* ******************************************************* */

void TOFDigitsInfo::GetDigits()
{

  for (Int_t iSector=0; iSector<fGeom->NSectors(); iSector++) {

    fNewTree = new TTree();




  }

}
