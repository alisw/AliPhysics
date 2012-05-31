/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

//
// Class to map TOF digit/raw data information
//
// Author: A. De Caro (email: decaro@sa.infn.it)
//

#include <TClonesArray.h>
#include <TTree.h>

//#include <TEveTreeTools.h>

#include <AliDAQ.h>
#include <AliLog.h>
#include <AliRawReader.h>

#include <AliTOFCableLengthMap.h>
#include <AliTOFdigit.h>
#include <AliTOFGeometry.h>
#include <AliTOFrawData.h>
#include <AliTOFRawStream.h>
#include <AliTOFDigitMap.h>

#include "AliEveTOFDigitsInfo.h"

//_________________________________________________________

ClassImp(AliEveTOFDigitsInfo)

  AliEveTOFDigitsInfo::AliEveTOFDigitsInfo(): 
    TObject(),
    TEveRefCnt(),
    fTree (0),
    fNewTree (0),
    fGeom (new AliTOFGeometry()),
    fTOFdigitMap(new AliTOFDigitMap())
{}
/* ******************************************************* */

AliEveTOFDigitsInfo:: ~AliEveTOFDigitsInfo() 
{
  //dtr

  delete fGeom;
  delete fTree;
  delete fNewTree;
  delete fTOFdigitMap;

}
/* ******************************************************* */

void AliEveTOFDigitsInfo::SetTree(TTree * const tree)
{
  //
  // Set fTree global variable
  //

  static const TEveException kEH("AliEveTOFDigitsInfo::SetTree ");
  
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
void AliEveTOFDigitsInfo::ReadRaw(AliRawReader* rawReader, Int_t newDecoder)
{
  //
  // Read raw-data. AliTOFdigit is used to
  // store raw-adata for all sub-detectors.
  //

  //AliTOFCableLengthMap *cableLength = new AliTOFCableLengthMap();

  //ofstream ftxt;
  //Char_t fileName[100];
  //sprintf(fileName,"TOFrawDataReadingFromEVE%d.txt",nEvent);

  //ftxt.open(fileName,ios::app);
  //ftxt << endl;
  //ftxt << "  " << nEvent << endl;

  //if (nEvent<0) printf("%3i\n", nEvent); // only to use nEvent variable

  const Int_t kDDL = AliDAQ::NumberOfDdls("TOF");

  TClonesArray *tofDigits = new TClonesArray("AliTOFdigit",10000);
  fTree = new TTree();
  fTree->Branch("TOF", &tofDigits, 32000);
  fTree->GetEntry(0);

  TClonesArray * clonesRawData = 0x0;

  Int_t detectorIndex[5];
  Int_t digit[4];

  AliTOFRawStream stream(rawReader);

  for (Int_t indexDDL = 0; indexDDL < kDDL; indexDDL++) {

    rawReader->Reset();
    if (newDecoder==0) stream.LoadRawData(indexDDL);
    else if (newDecoder==1) stream.LoadRawDataBuffers(indexDDL);
    else if (newDecoder==2) stream.LoadRawDataBuffersV2(indexDDL);

    clonesRawData = (TClonesArray*)stream.GetRawData();

    if (clonesRawData->GetEntriesFast()) AliDebug(2, Form(" Number of TOF digits in the sector number %2i: %5i", indexDDL, clonesRawData->GetEntriesFast()));

    for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

      AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

      if (tofRawDatum->GetTOF()==-1) continue;

      //Int_t cLenInt = Int_t(cableLength->GetCableTimeShift(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),tofRawDatum->GetTDC())*1000./AliTOFGeometry::TdcBinWidth());
      digit[0] = tofRawDatum->GetTOF();// - cLenInt;
      digit[1] = tofRawDatum->GetTOT();
      digit[2] = tofRawDatum->GetTOT();
      digit[3] = -1;

      /*
      if (indexDDL<10) ftxt << "  " << indexDDL;
      else             ftxt << " " << indexDDL;
      if (tofRawDatum->GetTRM()<10) ftxt << "  " << tofRawDatum->GetTRM();
      else                          ftxt << " " << tofRawDatum->GetTRM();
      ftxt << "  " << tofRawDatum->GetTRMchain();
      if (tofRawDatum->GetTDC()<10) ftxt << "  " << tofRawDatum->GetTDC();
      else                          ftxt << " " << tofRawDatum->GetTDC();
      ftxt << "  " << tofRawDatum->GetTDCchannel();
      */

      stream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
				  tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);

      /* check valid index */
      if (detectorIndex[0]==-1||detectorIndex[1]==-1||detectorIndex[2]==-1||detectorIndex[3]==-1||detectorIndex[4]==-1) continue;

      // Do not reconstruct anything in the holes
      if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
	if (detectorIndex[1]==2) { // plate with holes
	  continue;
	}
      }

      /*
      if (detectorIndex[0]<10) ftxt  << "  ->  " << detectorIndex[0];
      else                     ftxt  << "  -> " << detectorIndex[0];
      ftxt << "  " << detectorIndex[1];
      if (detectorIndex[2]<10) ftxt << "  " << detectorIndex[2];
      else                     ftxt << " " << detectorIndex[2];
      ftxt << "  " << detectorIndex[3];
      if (detectorIndex[4]<10) ftxt << "  " << detectorIndex[4];
      else                     ftxt << " " << detectorIndex[4];

      if (tofRawDatum->GetTOT()<10)                                            ftxt << "        " << tofRawDatum->GetTOT();
      else if (tofRawDatum->GetTOT()>=10 && tofRawDatum->GetTOT()<100)         ftxt << "       " << tofRawDatum->GetTOT();
      else if (tofRawDatum->GetTOT()>=100 && tofRawDatum->GetTOT()<1000)       ftxt << "      " << tofRawDatum->GetTOT();
      else if (tofRawDatum->GetTOT()>=1000 && tofRawDatum->GetTOT()<10000)     ftxt << "     " << tofRawDatum->GetTOT();
      else if (tofRawDatum->GetTOT()>=10000 && tofRawDatum->GetTOT()<100000)   ftxt << "    " << tofRawDatum->GetTOT();
      else if (tofRawDatum->GetTOT()>=100000 && tofRawDatum->GetTOT()<1000000) ftxt << "   " << tofRawDatum->GetTOT();
      else                                                                     ftxt << "  " << tofRawDatum->GetTOT();
      if (tofRawDatum->GetTOF()<10)                                            ftxt << "        " << tofRawDatum->GetTOF() << endl;
      else if (tofRawDatum->GetTOF()>=10 && tofRawDatum->GetTOF()<100)         ftxt << "       " << tofRawDatum->GetTOF() << endl;
      else if (tofRawDatum->GetTOF()>=100 && tofRawDatum->GetTOF()<1000)       ftxt << "      " << tofRawDatum->GetTOF() << endl;
      else if (tofRawDatum->GetTOF()>=1000 && tofRawDatum->GetTOF()<10000)     ftxt << "     " << tofRawDatum->GetTOF() << endl;
      else if (tofRawDatum->GetTOF()>=10000 && tofRawDatum->GetTOF()<100000)   ftxt << "    " << tofRawDatum->GetTOF() << endl;
      else if (tofRawDatum->GetTOF()>=100000 && tofRawDatum->GetTOF()<1000000) ftxt << "   " << tofRawDatum->GetTOF() << endl;
      else                                                                     ftxt << "  " << tofRawDatum->GetTOF() << endl;
      */

      TClonesArray &aDigits = *tofDigits;
      Int_t last = tofDigits->GetEntriesFast();

      fTOFdigitMap->AddDigit(detectorIndex, last);

      AliDebug(2,Form(" %3i -> %2i %2i %2i %2i %2i   %i  %i\n",
		      last, detectorIndex[0], detectorIndex[1],
		      detectorIndex[2], detectorIndex[4], detectorIndex[3],
		      digit[1], digit[0]));

      Int_t tracknum[3]={-1,-1,-1};
      new (aDigits[last]) AliTOFdigit(tracknum, detectorIndex, digit);

    } // while loop

    clonesRawData->Clear();

  } // DDL Loop

  fTree->Fill();

  //ftxt.close();

  //delete cableLength;

}


/* ******************************************************* */
void AliEveTOFDigitsInfo::LoadDigits()
{
  //
  // Load TOF digits
  //

  TClonesArray *digitsTOF = 0x0;
  AliTOFdigit *digs;

  fTree->SetBranchAddress("TOF",&digitsTOF);
  fTree->GetEntry(0);

  Int_t vol[5] = {-1,-1,-1,-1,-1};

  for (Int_t digitNumber=0; digitNumber<digitsTOF->GetEntries(); digitNumber++) {

    if (digitNumber==digitsTOF->GetEntries()-1)
      AliDebug(2,Form(" Hello  4 -> %3i digit of %i \n", digitNumber+1, digitsTOF->GetEntries()));

    digs = (AliTOFdigit*)digitsTOF->UncheckedAt(digitNumber);

    vol[0] = digs->GetSector(); // Sector Number (0-17)
    vol[1] = digs->GetPlate();  // Plate Number (0-4)
    vol[2] = digs->GetStrip();  // Strip Number (0-14/18)
    vol[3] = digs->GetPadx();   // Pad Number in x direction (0-47)
    vol[4] = digs->GetPadz();   // Pad Number in z direction (0-1)

    fTOFdigitMap->AddDigit(vol, digitNumber);
    if (digitNumber==digitsTOF->GetEntries()-1)
      AliDebug(2,Form(" I am inside LoadDigits %3i \n", digitNumber));

  }

}

/* ******************************************************* */

void AliEveTOFDigitsInfo::GetDigits(Int_t nSector, Int_t nPlate,
				    Int_t nStrip, Int_t nPadZ, Int_t nPadX,
				    Int_t indexDigit[3])
{
  //
  // Get TOF digit indices in the TOF volume
  // (nSector, nPlate,nStrip,nPadZ,nPadX)
  //

  Int_t vol[5] = {nSector,nPlate,nStrip,nPadX,nPadZ};

  fTOFdigitMap->GetDigitIndex(vol, indexDigit);
  //for (Int_t ii=1; ii<3; ii++) indexDigit[ii]=-1;

}
/* ******************************************************* */

TClonesArray* AliEveTOFDigitsInfo::GetDigits(Int_t nSector, Int_t nPlate,
					     Int_t nStrip)
{
  //
  // Get TOF digits in the TOF volume
  // (nSector, nPlate,nStrip)
  //

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

      fTOFdigitMap->GetDigitIndex(vol, nDigitsInVolume);

      for (Int_t ii=0; ii<3; ii++) {
	//if (ii!=0) continue;
	if (nDigitsInVolume[ii]>=0 ) {
	  AliDebug(2,Form("  nDigitsInVolume[%2i]  = %3i\n ", ii, nDigitsInVolume[ii]));
	  digs = (AliTOFdigit*)digitsTOF->UncheckedAt(nDigitsInVolume[ii]);
	  informations[0] = digs->GetTdc();
	  informations[1] = digs->GetAdc();
	  informations[2] = digs->GetToT();
	  informations[3] = digs->GetTdcND();
	  for(Int_t kk=0; kk<3; kk++) dummy[kk] = digs->GetTrack(kk);
	  new (ldigits[newCounter++]) AliTOFdigit(dummy, vol, informations);
	}

      }

      for (Int_t ii=0; ii<4; ii++) informations[ii]=-1;
      for (Int_t ii=0; ii<3; ii++) dummy[ii]=-1;
      for (Int_t ii=0; ii<3; ii++) nDigitsInVolume[ii]=-1;

    }
  }

  if (digitsTOFnew)
    AliDebug(2, Form("Sector %2i   Plate %1i  Strip %2i  -> number of digits %3i \n",
		     nSector, nPlate, nStrip, digitsTOFnew->GetEntries()));

  return digitsTOFnew;

}
/* ******************************************************* */

TClonesArray* AliEveTOFDigitsInfo::GetDigits(Int_t nSector)
{
  //
  // Get TOF digits in the TOF SM nSector
  //

  const Int_t kND = AliTOFDigitMap::kMaxDigitsPerPad;

  Int_t newCounter = 0;
  Int_t nDigitsInVolume[kND];
  Int_t dummy[3];
  Int_t informations[4];

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

	  for (Int_t ii=0; ii<4; ii++) informations[ii]=-1;
	  for (Int_t ii=0; ii<3; ii++) dummy[ii]=-1;
	  for (Int_t ii=0; ii<kND; ii++) nDigitsInVolume[ii]=-1;

	  vol[3] = iPadX;

	  fTOFdigitMap->GetDigitIndex(vol, nDigitsInVolume);

	  for (Int_t ii=0; ii<kND; ii++) {
	    //if (ii!=0) continue;
	    if (nDigitsInVolume[ii]>=0 ) {

	      digs = (AliTOFdigit*)digitsTOF->UncheckedAt(nDigitsInVolume[ii]);
	      informations[0] = digs->GetTdc();
	      informations[1] = digs->GetAdc();
	      informations[2] = digs->GetToT();
	      informations[3] = digs->GetTdcND();
	      for(Int_t kk=0; kk<3; kk++) dummy[kk] = digs->GetTrack(kk);
	      new (ldigits[newCounter++]) AliTOFdigit(dummy, vol, informations);

	      AliDebug(2,Form(" %2i -> %2i %2i %2i %2i %2i %7i %7i\n",
			      nDigitsInVolume[ii],
			      vol[0], vol[1], vol[2], vol[4], vol[3],
			      informations[1], informations[0]));

	    }

	  }

	}
      }
    }
  }

  if (digitsTOFnew)
    AliDebug(2,Form("Sector %2i  -> number of digits %3i \n",
		    nSector, digitsTOFnew->GetEntries()));

  return digitsTOFnew;

}
/* ******************************************************* */

Int_t AliEveTOFDigitsInfo::GetTOFInfos() const
{
  //
  // Return number of TOF digits
  //

  return fTOFdigitMap->GetFilledCellNumber();

}

/* ******************************************************* */
Int_t AliEveTOFDigitsInfo::IsStripFilled(Int_t iSector, Int_t iPlate, Int_t iStrip)
{
  //
  // Return number of TOF digits
  // in volume (iSector,iPlate,iStrip)
  //

  Int_t vol[5] = {iSector, iPlate, iStrip, -1, -1};

  Int_t index = 0;

  for (Int_t iPadZ=0; iPadZ<fGeom->NpadZ(); iPadZ++) 
    for (Int_t iPadX=0; iPadX<fGeom->NpadX(); iPadX++) 
      {
	vol[3] = iPadX;
	vol[4] = iPadZ;
	if (fTOFdigitMap->GetDigitIndex(vol,0)>=0) index++;
      }

  return index;

}

/* ******************************************************* */
/*
void AliEveTOFDigitsInfo::GetDigits()
{

  for (Int_t iSector=0; iSector<fGeom->NSectors(); iSector++) {

    fNewTree = new TTree();




  }

}
*/
