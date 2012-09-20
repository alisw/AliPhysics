#if !defined(__CINT__) || defined(__MAKECINT__)


#include <iostream>

// root
#include <TStyle.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>



// aliroot
#include "AliRun.h"
#include "AliTRDgeometry.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"
#include "AliTRDfeeParam.h"
#include "AliTRDdigitsParam.h"



#endif



void AliTRDanadigits()
{
  gStyle->SetNdivisions(005,"X");

  Int_t iev1 = 1;
  Int_t iev2 = 20;
  //
  // Analyzes the digits
  //
  const Int_t kNdet = AliTRDgeometry::Ndet();

  //open run loader and load gAlice, kinematics and header
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (!rl)   return;
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  if (!gAlice) return;

  // Import the Trees for the event nEvent in the file
  rl->LoadKinematics();
  rl->GetEvent(0);
  rl->GetHeader();

  AliLoader* loader = rl->GetLoader("TRDLoader");
  if (!loader) {
    cout << "<AliTRDanalyzeHits> No TRDLoader found" << endl;
    return;
  }

  //TH2D
  Int_t detA = 240;//AliTRDgeometry::GetDetector(0,4,1); //layer,stac,sm
  Int_t detB = 250;//AliTRDgeometry::GetDetector(5,4,1);
  printf("detA %d \t detB %d \n",detA,detB);
  TH2D *chA = new TH2D("chA", "chA", kNdet, 0, kNdet, 200, -50, 50);
  TH2D *chB = new TH2D("chB", "chB", kNdet, 0, kNdet, 200, -50, 50);

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();

  for (Int_t ievent = iev1; ievent < iev2; ievent++) {
    printf("Process event %d\n",ievent);
    rl->GetEvent(ievent);
    if (!loader->TreeD()){
      printf("loader Loading Digits ... \n");
      loader->LoadDigits();
    }

    // Read the digits from the file
    if (!(digitsManager->ReadDigits(loader->TreeD()))) {
      cout << "<AliTRDanalyzeDigits> Cannot read the digits" << endl;
      return;
    }
    else {
      printf("digitsManager Read Digits Done\n");
    }

    // loop over selected chambers
    for(Int_t det = detA; det < detB; det++) {

      printf("try in %d :\t",det);
      Int_t        plan = AliTRDgeometry::GetLayer( det );   // Plane
      Int_t        cham = AliTRDgeometry::GetStack( det ); // Chamber
      Int_t        sect = AliTRDgeometry::GetSector( det );  // Sector (=iDDL)
      Int_t        nRow = AliTRDgeometry::GetRowMax( plan, cham, sect );
      Int_t        nCol = AliTRDgeometry::GetColMax( plan );
      printf(" sm %d \t lay %d \t stac %d \n",sect,plan,cham);


      AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(det); //mod
      digits->Expand();
      if (digits->HasData()) { //array

	printf("data\n");

	Int_t nbtimebin  = digitsManager->GetDigitsParam()->GetNTimeBins(det);

	AliTRDfeeParam *feeParam = AliTRDfeeParam::Instance();  
	// Signal for regular pads
	for (Int_t irow = 0; irow < nRow; irow++ ) {
	  for (Int_t icol = 0; icol < nCol; icol++ ) {

	    Int_t rob = feeParam->GetROBfromPad(irow,icol);

	    for(Int_t k = 0; k < nbtimebin; k++){
	      //	      printf(" row %d \t col %d \t timebin %d \t --> signal %d \n",irow,icol,k,signal);
	      Short_t signal = 999;
	      signal = (Short_t) digits->GetData(irow,icol,k);

	      if(signal==999) continue;
	      if(rob%2==0) { chB->Fill(det,signal); }
	      else {chA->Fill(det,signal); }
	    } // loop: ntimebin
	  } // loop: nCol
	} //loop: nRow
      } // if: hasdata
    } // loop: det
    loader->UnloadDigits();  
  }// event

  // Display the detector matrix
  TCanvas *c = new TCanvas("c","",50,50,900,600);
  c->Divide(2,1);

  chA->SetAxisRange(detA,detB-1,"X");
  chA->SetTitle("HalfChamberSide A");
  chA->SetXTitle("Detector number");
  chA->SetYTitle("charge deposit [a.u]");
  chB->SetAxisRange(detA,detB-1,"X");
  chB->SetTitle("HalfChamberSide B");
  chB->SetXTitle("Detector number");
  chB->SetYTitle("charge deposit [a.u]");

  c->cd(1);
  chA->Draw("COLZ");

  c->cd(2);
  chB->Draw("COLZ");

  //c->SaveAs("ChargeDeposit.pdf");
}

