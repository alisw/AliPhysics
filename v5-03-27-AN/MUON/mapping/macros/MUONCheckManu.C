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

// $Id$

// By Ch. Finck, Subatech

// 17/03/09  Does not compile

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

#include "TExMap.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "AliCDBManager.h"
#include "AliMpHelper.h"
#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpIntPair.h"
#include "AliMpFiles.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpDDLStore.h"
#include "AliMpCDB.h"

#endif

const Int_t kStrangeBin = 26;
const Int_t kUndefBin   = 36;
const Int_t kBadCalBin  = 46;
const Int_t kBadBin       = 5;

const Int_t kNotScanSerial  = 40000;

void MUONCheckManu(Int_t iCh = 10, Bool_t rootFile = kFALSE)
{

    // Macro to check the bin number for manu in the CDB for each chamber
    // chamber = 0 means all chambers 
    // The normal bin number ranges from 0-4.
    // the abnormal bin number are as follow:
    // 25 -> strange, 
    // 35 ->undefined to be re-tested, 
    // 45 -> bad calib0,
    // 5-10 -> reading error on manu.
    // cath: 0 -> bending, 1 -> non-bending
    // author: Ch. Finck

    static Int_t nBad    = 0;
    static Int_t nBadS  = 0;
    static Int_t nBadU  = 0;
    static Int_t nDB     = 0;
    static Int_t nManu = 0;
    static Int_t badDE = 0;


    // Load mapping
    if ( !AliMpCDB::LoadMpSegmentation() ) {
       printf("Could not access mapping from OCDB !\n");
    }
    
    // Load DDL store
    if ( !AliMpCDB::LoadDDLStore() ) {
      printf("Could not access DDL Store from OCDB !\n");
      exit(0);
    }
    AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
			  
    TString fileNameIn;
    TExMap serialNbtoBinId;
    TExMap serialNbtoGain;
    TExMap serial;

    char fileNameOut[255];
    sprintf(fileNameOut,"serial_number_%d.dat",iCh);

    char fileNameOut1[255];
    sprintf(fileNameOut1,"serial_number_bad_%d.dat",iCh);


    char histoFile[255];
    TFile* fileOut = 0x0;

    if (rootFile) {
      sprintf(histoFile,"GainDE_%d.root",iCh);
      fileOut = new TFile(histoFile,"RECREATE","DE manu gain");
    }


    fileNameIn  = AliMpFiles::SerialToBinFilePath();
  
    ifstream in(fileNameIn, ios::in);
    if (!in) printf("file %s not found\n", fileNameIn.Data());

    ofstream out(fileNameOut, ios::out);
  
    ofstream out1(fileNameOut1, ios::out);

    char line[255];
    Int_t binId          = 0;
    Int_t manuSerial = 0;
    Float_t gain        = 0;

    // reading file
    while ( in.getline(line,255) ) {

      if ( line[0] == '#' ) continue;
      if ( line[0] == '\n' ) continue;
      if ( line[0] == '\0' ) continue;

      TString tmp(AliMpHelper::Normalize(line));

      TObjArray* stringList = tmp.Tokenize(TString(" "));

      TString sSM   = ((TObjString*)stringList->At(2))->GetString();

      TString sBin  = ((TObjString*)stringList->At(4))->GetString();

      TString sGain  = ((TObjString*)stringList->At(6))->GetString();

      manuSerial =  atoi(sSM.Data());
      binId  =  atoi(sBin.Data()) + 1;

      if (sBin.CompareTo("S") == 0) {
	binId = kStrangeBin;
      }

      if (sBin.CompareTo("I") == 0) {
	binId = kUndefBin;
      }

      if (sBin.CompareTo("C") == 0) {
	binId = kBadCalBin;
      }

      gain  =  atof(sGain.Data());
      //printf("%s %d   %s %d  %s %f\n", sSM.Data(), manuSerial, sBin.Data(),  binId, sGain.Data(), gain);

      fflush(stdout);

  
      // to avoid that a manu that is at the same time in a normal bin and in a bin S, stays in bin S
      if (serialNbtoBinId.GetValue(manuSerial)  != 0 && binId == kStrangeBin)  
	  continue;
     

      if (serialNbtoBinId.GetValue(manuSerial) != 0) {
	serialNbtoBinId.Remove((Long_t)manuSerial);
	serialNbtoGain.Remove((Long_t)manuSerial);
      }


      serialNbtoBinId.Add((Long_t)manuSerial, (Long_t)binId); 
      serialNbtoGain.Add((Long_t)manuSerial, (Long_t)(gain*100000)); 
      delete stringList;

    }
    in.close();

//  return;

// ///////////
// // iterates over DE
  
    char title[256];
    char name[256];

    AliMpDEIterator it;

    out1 << "Bin:" << endl;
    out1 << "  25 -> Strange, 35 -> re-tested, 45 -> bad calib0" << endl;
    out1 << "  >= 5 error id, 10 -> 1" << endl;
    out1 << "cath: 0 -> bending, 1 -> non-bending" << endl;


    TH1F* hslat_manu_bp[1100];
    TH1F* hslat_manu_nbp[1100];

    TH1F* hCh_manu_bp[11];
    TH1F* hCh_manu_nbp[11];

    static TH1F* hAll_manu_bp;
    static TH1F* hAll_manu_nbp;
    static TH1F* hAll_manu;
    static TH1F* hBinS;

    TH2F* hBin_DE = 0x0;

    if (rootFile) {

	sprintf(name,"hAll_manu_bp");
	sprintf(title,"MANU gain total for bending");
	hAll_manu_bp = new TH1F(name,title,300,1.5,5.0);
	hAll_manu_bp->SetDirectory(fileOut);

	sprintf(name,"hAll_manu_npb");
	sprintf(title,"MANU gain total for non-bending");
	hAll_manu_nbp = new TH1F(name,title,300,1.5,5.0);
	hAll_manu_nbp->SetDirectory(fileOut);

	sprintf(name,"hAll_manu");
	sprintf(title,"MANU gain total");
	hAll_manu = new TH1F(name,title,300,1.5,5.0);
	hAll_manu->SetDirectory(fileOut);

	sprintf(name,"hBin_DE");
	sprintf(title,"DE versus Bin");
	hBin_DE = new TH2F(name,title,5,0,5, 1100, 1, 1100);
	hBin_DE->SetDirectory(fileOut);

	sprintf(name,"hBin_S");
	sprintf(title,"Bin_S");
	hBinS = new TH1F(name,title, 300, 600, 1100);
	hBinS->SetDirectory(fileOut);
    }

    Int_t begCh;
    Int_t endCh;

    if (iCh == 0) {
	begCh = 1;
	endCh = 11;
    } else {
	begCh = iCh;
	endCh = iCh+1;
    }

    for (iCh = begCh; iCh < endCh; ++iCh) {

      if (rootFile) {

	sprintf(name,"hCh_manu_bp_CH%d", iCh);
	sprintf(title,"MANU gain total for CH%d bending", iCh);
	hCh_manu_bp[iCh] = new TH1F(name,title,300,1.5,5.0);
	hCh_manu_bp[iCh]->SetDirectory(fileOut);

	sprintf(name,"hCh_manu_nbp_CH%d", iCh);
	sprintf(title,"MANU gain total for CH%d non-bending", iCh);
	hCh_manu_nbp[iCh] = new TH1F(name,title,300,1.5,5.0);
	hCh_manu_nbp[iCh]->SetDirectory(fileOut);
      }

      for ( it.First(iCh-1); ! it.IsDone(); it.Next() ) {

	Int_t flag = 0;
	Int_t iDE = it.CurrentDE()->GetId();

	AliMpDetElement* detElem =  AliMpDEManager::GetDetElement(iDE);
	TString nameDE = detElem->GetDEName();

	out << endl;
	out << " DE: " << iDE << endl;
	out << "manuId  serial  binId" << endl;

	out1 << endl;
	out1 << " DE: " << iDE << " name: " << nameDE.Data() << endl;
	out1 << "manuId  serial  binId cath" << endl;

	printf("\nDE:%d\n", iDE);

	if (rootFile) {

	  sprintf(name,"hslat_manu_bp%d",iDE);
	  sprintf(title,"MANU gain for DE %d bending",iDE);
	  hslat_manu_bp[iDE] = new TH1F(name,title,300,1.5,4.5);
	  hslat_manu_bp[iDE]->SetDirectory(fileOut);

	  sprintf(name,"hslat_manu_nbp%d",iDE);
	  sprintf(title,"MANU gain for DE %d non-bending",iDE);
	  hslat_manu_nbp[iDE] = new TH1F(name,title,300,1.5,4.5);
	  hslat_manu_nbp[iDE]->SetDirectory(fileOut);
	}

	TList manuList;

	for ( Int_t cath = 0; cath <=1 ; ++cath ) {
	  const AliMpVSegmentation* seg 
	      = AliMpSegmentation::Instance()->GetMpSegmentation(iDE,AliMp::GetCathodType(cath));
        
	  TArrayI manus;

	  seg->GetAllElectronicCardIDs(manus);
          
	  // filling
	  for ( Int_t im = 0; im < manus.GetSize(); ++im ) {
	    AliMpIntPair* manu = 0x0;
	    if (manus[im] > AliMpConstants::ManuMask(AliMp::kNonBendingPlane))
		manu = new AliMpIntPair(manus[im], 1, kTRUE);
	    else
		manu = new AliMpIntPair(manus[im], 0, kTRUE);

	    manuList.Add(manu);
	  }

	}
	manuList.Sort();
   

	// check manu serial
	for (Int_t iEntry = 0; iEntry < manuList.GetEntries(); ++iEntry) {
	  AliMpIntPair* manuPtr = (AliMpIntPair*)manuList.At(iEntry);

	  AliMpDetElement* detElem =  ddlStore->GetDetElement(iDE);
	  manuSerial = detElem->GetManuSerialFromId(manuPtr->GetFirst());

	  binId =  (Int_t)serialNbtoBinId.GetValue(manuSerial);
	  if (manuSerial == 0) {
	    printf("manu %d not found in mapping\n",manuPtr->GetFirst()); 
	    continue;
	  }
	  if (manuSerial != 0)
	      nManu++;

	  if (!binId) {
	    printf("Bin for manu %d with serial %d not available\n",manuPtr->GetFirst(), manuSerial); 
	  }
	  if (binId > kBadBin && binId < kStrangeBin)
	      printf("Bin for manu %d with serial %d is bad\n",manuPtr->GetFirst(), manuSerial); 

	  if (binId == kStrangeBin)
	      printf("Bin for manu %d with serial %d is strange\n",manuPtr->GetFirst(), manuSerial); 

	  if (binId == kUndefBin)
	      printf("Bin for manu %d with serial %d is unidentified\n",manuPtr->GetFirst(), manuSerial); 

	  gain  =  (Float_t)(serialNbtoGain.GetValue(manuSerial)/100000.);

	  out  << setw(4) << ( manuPtr->GetFirst() & 0x3FF) << " " << setw(5) << manuSerial 
	       << " " << setw(2) << binId-1 << " " <<   manuPtr->GetSecond() 
	       << " " << setw(7) << gain << endl;

	  // bad manu
	  if ((binId == 0 || binId == kUndefBin) && manuSerial < kNotScanSerial) nBadU++;
	  if (binId == 0 && manuSerial >= kNotScanSerial) nDB++;

	  // check whether serial number appears twice
	  if (serial.GetValue(manuSerial) == 0) 
	    serial.Add((Long_t)manuSerial, (Long_t)binId);
	  else {
	    out1 << setw(4) << ( manuPtr->GetFirst() & 0x3FF) << " " << setw(5) << manuSerial 
		   << " " << setw(2) << binId-1 << " " <<   manuPtr->GetSecond() 
		   << " " << setw(7) << gain << " Serial number Twice" << endl;  
	  }

	  // count bad/strange manus
	  if (!binId || binId > kBadBin) {
	    if (binId == kStrangeBin) nBadS++;
	    if (binId > kBadBin) {
	      nBad++;
	      flag = 1;
	    }

	    out1 << setw(4) << ( manuPtr->GetFirst() & 0x3FF) << " " << setw(5) << manuSerial 
		 << " " << setw(2) << binId-1 << " " <<   manuPtr->GetSecond() 
		 << " " << setw(7) << gain << endl;
	  }

	  if (rootFile && binId == kStrangeBin) 
	      hBinS->Fill(iDE);

	  if (rootFile && binId != kStrangeBin) {
	    hAll_manu->Fill(gain);
	    hBin_DE->Fill(binId-1, iDE);

	    if(manuPtr->GetSecond()) {
		hslat_manu_nbp[iDE]->Fill(gain);
		hAll_manu_nbp->Fill(gain);
		hCh_manu_nbp[iCh]->Fill(gain);
	    } else {
		hslat_manu_bp[iDE]->Fill(gain);
		hAll_manu_bp->Fill(gain);
		hCh_manu_bp[iCh]->Fill(gain);
	    }
	
	  }
	} // manu

	if (flag) badDE++;

	manuList.Delete();
      }// DE

    } // ich

    serialNbtoBinId.Delete();
    serialNbtoGain.Delete();

    if (rootFile) {
      fileOut->Write();
      fileOut->Close();
    }
    
    out.close();
    out1.close();

    printf("\n");
    printf("Number of bad manus %d and strange manus %d total %d\n", nBad, nBadS, nManu);
    printf("Number of unidentified manus %d\n", nBadU);
    printf("Number of manus not scanned %d\n", nDB);
    printf("Number of bad DE %d\n", badDE);
}

