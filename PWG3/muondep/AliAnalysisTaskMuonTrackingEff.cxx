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

//Class to calculate the intrinsic efficiency of the detection elements of the
//MUON tracking chambers in function of the position in the detection element.
//Author:  Nicolas LE BRIS - SUBATECH Nantes


//ROOT includes
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TH2F.h>

//ANALYSIS includes
#include "AliAnalysisManager.h"

//STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliAnalysisManager.h"
#include "AliGeomManager.h"


//PWG3/muondep includes
#include "AliAnalysisTaskMuonTrackingEff.h"
#include "AliCheckMuonDetEltResponse.h"

//MUON includes
#include "AliMUONTrackExtrap.h"
#include "AliMUONGeometryTransformer.h"

ClassImp(AliAnalysisTaskMuonTrackingEff)

const Int_t AliAnalysisTaskMuonTrackingEff::fgkTotNbrOfDetectionElt  = 156;
const Int_t AliAnalysisTaskMuonTrackingEff::fgkTotNbrOfChamber  = 10;

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff()
  :
  AliAnalysisTask(),
  fESD(0x0),
  fTransformer(0x0),
  fDetEltEffHistList(0x0),
  fDetEltTDHistList(0x0),
  fDetEltTTHistList(0x0),
  fChamberEffHistList(0x0),
  fChamberTDHistList(0x0),
  fChamberTTHistList(0x0),
  fChamberEff(0x0),
  fIsCosmicData(kFALSE)
{
/// Default constructor
}
//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& src)
  :
  AliAnalysisTask(src),
  fESD(0x0),
  fTransformer(0x0),
  fDetEltEffHistList(0x0),
  fDetEltTDHistList(0x0),
  fDetEltTTHistList(0x0),
  fChamberEffHistList(0x0),
  fChamberTDHistList(0x0),
  fChamberTTHistList(0x0),
  fChamberEff(0x0),
  fIsCosmicData(kFALSE)
{
  /// copy ctor
  src.Copy(*this);
}
//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff& AliAnalysisTaskMuonTrackingEff::operator=(const AliAnalysisTaskMuonTrackingEff& src)
{
  /// assignement operator
  if ( this != &src ) 
  {
    src.Copy(*this);
  }
  return *this;
}

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff(const char* name,
							       Bool_t isCosmic)
  :
  AliAnalysisTask(name, "AnalysisTaskESD"),
  fESD(0x0),
  fTransformer(0x0),
  fDetEltEffHistList(0x0),
  fDetEltTDHistList(0x0),
  fDetEltTTHistList(0x0),
  fChamberEffHistList(0x0),
  fChamberTDHistList(0x0),
  fChamberTTHistList(0x0),
  fChamberEff(0x0),
  fIsCosmicData(kFALSE)
{
//Constructor
//-----------
  
  fIsCosmicData = isCosmic;  

//Load the geometry
  if (!AliGeomManager::GetGeometry())
    AliGeomManager::LoadGeometry();
  fTransformer = new AliMUONGeometryTransformer();
  fTransformer->LoadGeometryData();  


//Define detection element efficiency histograms
//----------------------------------------------

    fDetEltEffHistList = new TClonesArray("TH2F",fgkTotNbrOfDetectionElt + 3); //!<+3 for: 1.the total efficiency chamber by chamber histogram.
                                                                //!<        2.the total number of tracks detected by the tracking system.
                                                                //!<        3.the number of track used for the efficiency calculation chamber by chamber
    fDetEltTDHistList  = new TClonesArray("TH2F",fgkTotNbrOfDetectionElt + 1);
    fDetEltTTHistList  = new TClonesArray("TH2F",fgkTotNbrOfDetectionElt + 1);


    fChamberEffHistList = new TClonesArray("TH1F", fgkTotNbrOfChamber + 3);
    fChamberTDHistList = new TClonesArray("TH1F", fgkTotNbrOfChamber + 1);
    fChamberTTHistList = new TClonesArray("TH1F", fgkTotNbrOfChamber + 1);
 


    for (Int_t i = 0; i<fgkTotNbrOfDetectionElt; ++i)
    {
      Int_t iDetElt = 0;
      if      (i<16) iDetElt = 100*(i/4+1)       + i -  int(i/4)*4;
      else if (i<52) iDetElt = 100*((i-16)/18+5) + i - (int((i-16)/18)*18 + 16);
      if      (i>51) iDetElt = 100*((i-52)/26+7) + i - (int((i-52)/26)*26 + 52); 

      Char_t histName[255]; 
      Char_t histTitle[255];

      sprintf(histName,"Eff_detEltNbr%d",iDetElt);
      sprintf(histTitle,"detEltNbr %d",iDetElt);
      if(i<16) new((*fDetEltEffHistList)[i]) TH2F(histName,histTitle,12,-10,110,12,-10,110);//!<Stations 1 & 2.
      else     new((*fDetEltEffHistList)[i]) TH2F(histName,histTitle,28,-140,140,8,-40,40); //!<Stations 3 -> 5.
      
      sprintf(histName,"TD_detEltNbr%d",iDetElt);
      if(i<16) new((*fDetEltTDHistList)[i]) TH2F(histName,histTitle,12,-10,110,12,-10,110);//!<Stations 1 & 2.
      else     new((*fDetEltTDHistList)[i]) TH2F(histName,histTitle,28,-140,140,8,-40,40); //!<Stations 3 -> 5.

      sprintf(histName,"TT_detEltNbr%d",iDetElt);
      if(i<16) new((*fDetEltTTHistList)[i]) TH2F(histName,histTitle,12,-10,110,12,-10,110); //!<Stations 1 & 2.
      else     new((*fDetEltTTHistList)[i]) TH2F(histName,histTitle,28,-140,140,8,-40,40);  //!<Stations 3 -> 5.
      
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> GetXaxis() -> SetTitle("X (cm)");
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> GetYaxis() -> SetTitle("Y (cm)");
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> GetZaxis() -> SetTitle("Efficiency (%)");  
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> GetXaxis() -> SetTitleOffset(1.8);
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> GetYaxis() -> SetTitleOffset(1.8); 
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> GetZaxis() -> SetTitleOffset(1.2); 
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i)) -> SetOption("LEGO");
    }

    for (Int_t j = 0; j < fgkTotNbrOfChamber; j++)
      {
	Char_t histName[255]; 
	Char_t histTitle[255];


	sprintf(histName,"Eff_ChamberNbr%d",j+1);
	sprintf(histTitle,"ChamberNbr %d",j+1);
	if (j<4) new ((*fChamberEffHistList)[j]) TH1F(histName, histTitle, 4, 0.0, 4.0);
	else if (j<6) new ((*fChamberEffHistList)[j]) TH1F(histName, histTitle, 18, 0.0, 18.0);
	if (j>=6) new ((*fChamberEffHistList)[j]) TH1F(histName, histTitle, 26, 0.0, 26.0);

	sprintf(histName,"TD_ChamberNbr%d",j+1);
	if (j<4) new ((*fChamberTDHistList)[j]) TH1F(histName, histTitle, 4, 0.0, 4.0);
	else if (j<6) new ((*fChamberTDHistList)[j]) TH1F(histName, histTitle, 18, 0.0, 18.0);
	if (j>=6) new ((*fChamberTDHistList)[j]) TH1F(histName, histTitle, 26, 0.0, 26.0);

	sprintf(histName,"TT_ChamberNbr%d",j+1);
	if (j<4) new ((*fChamberTTHistList)[j]) TH1F(histName, histTitle, 4, 0.0, 4.0);
	else if (j<6) new ((*fChamberTTHistList)[j]) TH1F(histName, histTitle, 18, 0.0, 18.0);
	if (j>=6) new ((*fChamberTTHistList)[j]) TH1F(histName, histTitle, 26, 0.0, 26.0);

	((TH1F*) fChamberEffHistList->UncheckedAt(j)) -> GetXaxis() -> SetTitle("DetElement");
	((TH1F*) fChamberEffHistList->UncheckedAt(j)) -> GetYaxis() -> SetTitle("Efficiency (%)");
	((TH1F*) fChamberEffHistList->UncheckedAt(j)) -> GetXaxis() -> SetTitleOffset(1.8);
	((TH1F*) fChamberEffHistList->UncheckedAt(j)) -> GetYaxis() -> SetTitleOffset(1.8); 
	((TH1F*) fChamberEffHistList->UncheckedAt(j)) -> Sumw2();	
      }

    new((*fDetEltTDHistList )[fgkTotNbrOfDetectionElt]) TH2F("TD_Chamber" ,"TD_Chamber" ,10,0,10,1,0,1); //!<Detected tracks.
    new((*fDetEltTTHistList )[fgkTotNbrOfDetectionElt]) TH2F("TT_Chamber" ,"TT_Chamber" ,10,0,10,1,0,1); //!<Tracks total number.
    new((*fDetEltEffHistList)[fgkTotNbrOfDetectionElt]) TH2F("fChamberEff","fChamberEff",10,0,10,1,0,1); //!<Chamber efficiency.

    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> GetXaxis() -> SetTitle("Chamber number");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> GetYaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> GetZaxis() -> SetTitle("Efficiency (%)");  
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> GetZaxis() -> SetTitleOffset(1.2); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt)) -> SetOption("LEGO");
  
    new((*fDetEltEffHistList)[157]) TH2F("TT_Chamber" ,"TT_Chamber" ,10,0,10,1,0,1); //!<Tracks total number by chamber.

    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> GetXaxis() -> SetTitle("Chamber number");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> GetYaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> GetZaxis() -> SetTitle("Number of tracks");  
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> GetZaxis() -> SetTitleOffset(1.2); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 1)) -> SetOption("LEGO");

    new((*fDetEltEffHistList)[158]) TH2F("Total_Number_of_Tracks" ,"Total_Number_of_Tracks" ,1,0,1,1,0,1); //!<Tracks total number.

    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> GetXaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> GetYaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> GetZaxis() -> SetTitle("Number of tracks");  
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> GetZaxis() -> SetTitleOffset(1.2); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fgkTotNbrOfDetectionElt + 2)) -> SetOption("LEGO");    


    new((*fChamberTDHistList )[fgkTotNbrOfChamber]) TH1F("TD_Chamber_2" ,"TD_Chamber_2" ,10,0,10); //!<Detected tracks.
    new((*fChamberTTHistList )[fgkTotNbrOfChamber]) TH1F("TT_Chamber_2" ,"TT_Chamber_2" ,10,0,10); //!<Tracks total number.
    new((*fChamberEffHistList)[fgkTotNbrOfChamber]) TH1F("fChamberEff_2","fChamberEff_2",10,0,10); //!<Chamber efficiency.

    ((TH1F*) fChamberTDHistList->UncheckedAt(fgkTotNbrOfChamber)) -> Sumw2();
    ((TH1F*) fChamberTTHistList->UncheckedAt(fgkTotNbrOfChamber)) -> Sumw2();

    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber)) -> GetXaxis() -> SetTitle("Chamber number");
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber)) -> GetYaxis() -> SetTitle("");
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber)) -> Sumw2();
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber)) -> SetOption("");
  
    new((*fChamberEffHistList)[11]) TH1F("TT_Chamber_2" ,"TT_Chamber_2" ,10,0,10); //!<Tracks total number by chamber.

    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 1)) -> GetXaxis() -> SetTitle("Chamber number");
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 1)) -> GetYaxis() -> SetTitle("");
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 1)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 1)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 1)) -> Sumw2();
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 1)) -> SetOption("");

    new((*fChamberEffHistList)[12]) TH1F("Total_Number_of_Tracks_2" ,"Total_Number_of_Tracks_2" ,1,0,1); //!<Tracks total number.

    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 2)) -> GetXaxis() -> SetTitle("");
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 2)) -> GetYaxis() -> SetTitle("");
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 2)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 2)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 2)) -> Sumw2();    
    ((TH1F*) fChamberEffHistList->UncheckedAt(fgkTotNbrOfChamber + 2)) -> SetOption("");    




 
//Define input & output
//---------------------

// -Input slot 0 works with a TChain:
    DefineInput(0, TChain::Class());

// -Output slots 0 to 5 writes into a TClonesArray:
    DefineOutput(0, TClonesArray::Class());
    DefineOutput(1, TClonesArray::Class());
    DefineOutput(2, TClonesArray::Class());
    DefineOutput(3, TClonesArray::Class());
    DefineOutput(4, TClonesArray::Class());
    DefineOutput(5, TClonesArray::Class());
}



//______________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::~AliAnalysisTaskMuonTrackingEff()
{
// Destructor.
  delete fTransformer;
  delete fDetEltEffHistList;
  delete fDetEltTDHistList;
  delete fDetEltTTHistList;
  delete fChamberEffHistList;
  delete fChamberTDHistList;
  delete fChamberTTHistList;
  delete fChamberEff;
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::CreateOutputObjects()
{
    OpenFile(0);
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::ConnectInputData(Option_t */*option*/)

{
  //Set Field Map for track extrapolation
  AliMUONTrackExtrap::SetField();

  //Connect input
    AliESDInputHandler* esdHandler = (AliESDInputHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if (!esdHandler)
    {
      Printf("ERROR: Could not get ESDInputHandler");
    } 
    else fESD = esdHandler->GetEvent();
}


//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Exec(Option_t */*option*/)
{
//Execute analysis for current event
  if (fChamberEff == 0x0)
    fChamberEff = new AliCheckMuonDetEltResponse(fTransformer, fESD, fDetEltTDHistList, fDetEltTTHistList, fChamberTDHistList, fChamberTTHistList);
  fChamberEff->CheckDetEltResponse();

    for( Int_t i = 0; i<156; ++i)
    {
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i))-> Divide( (TH2F*) fDetEltTDHistList->UncheckedAt(i),
							     (TH2F*) fDetEltTTHistList->UncheckedAt(i), 100., 1.); 
    }

    for (Int_t j = 0; j < 10; j++)
      {
	((TH1F*) fChamberEffHistList->UncheckedAt(j))-> Divide( (TH1F*) fChamberTDHistList->UncheckedAt(j),
								(TH1F*) fChamberTTHistList->UncheckedAt(j), 100., 1.); 	
      }


    ((TH2F*) fDetEltEffHistList->UncheckedAt(156))-> Divide( (TH2F*) fDetEltTDHistList->UncheckedAt(156),
							     (TH2F*) fDetEltTTHistList->UncheckedAt(156), 100., 1.); 


    ((TH2F*) fDetEltEffHistList->UncheckedAt(157))-> Add   ( (TH2F*) fDetEltTTHistList ->UncheckedAt(156),
							     (TH2F*) fDetEltEffHistList->UncheckedAt(157), 1., 0.); 


    ((TH2F*) fDetEltEffHistList->UncheckedAt(158))-> Fill(0., 0., ((Double_t)fESD -> GetNumberOfMuonTracks()));


    ((TH1F*) fChamberEffHistList->UncheckedAt(10))-> Divide( (TH1F*) fChamberTDHistList->UncheckedAt(10),
							     (TH1F*) fChamberTTHistList->UncheckedAt(10), 100., 1.); 


    ((TH1F*) fChamberEffHistList->UncheckedAt(11))-> Add   ( (TH1F*) fChamberTTHistList ->UncheckedAt(10),
							     (TH1F*) fChamberEffHistList->UncheckedAt(11), 1., 0.); 


    ((TH1F*) fChamberEffHistList->UncheckedAt(12))-> Fill(0., ((Double_t)fChamberEff -> GetNbrUsableTracks()));
    fChamberEff->SetNbrUsableTracks(0);
    //((TH1F*) fChamberEffHistList->UncheckedAt(12))-> Fill(0., ((Double_t)fESD -> GetNumberOfMuonTracks()));

    ComputeErrors();

//Post the output data:
    PostData(0, fDetEltTDHistList);  
    PostData(1, fDetEltTTHistList);  
    PostData(2, fDetEltEffHistList);  
    PostData(3, fChamberTDHistList);
    PostData(4, fChamberTTHistList);
    PostData(5, fChamberEffHistList);
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Terminate(Option_t */*option*/)
{
//Terminate analysis

}


//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::ComputeErrors()
{
  // Compute error on the efficiency
  // eff = Ntd/Ntt
  // error = max {1/Ntt, sqrt(eff*(1-eff)/Ntt)}

  for (Int_t ii = 0; ii <= 10; ii++)
    {
      Int_t NumberOfBins = ((TH1F*) fChamberEffHistList->UncheckedAt(ii))->GetNbinsX();
      for (Int_t jj = 1; jj <= NumberOfBins; jj++)
	{
	  Double_t Ntd = ((TH1F*) fChamberTDHistList->UncheckedAt(ii))->GetBinContent(jj);
	  Double_t Ntt = ((TH1F*) fChamberTTHistList->UncheckedAt(ii))->GetBinContent(jj);

	  if (Ntt > 0.0 && Ntd > 0.0)
	    {
	      Double_t eff = ((TH1F*) fChamberEffHistList->UncheckedAt(ii))->GetBinContent(jj)/100.0;
	      Double_t err1 = 1.0/Ntt;
	      Double_t err2 = TMath::Sqrt(eff*(1.0 - eff)/Ntt);
	      Double_t error = TMath::Max(err1, err2);

	      ((TH1F*) fChamberEffHistList->UncheckedAt(ii))->SetBinError(jj, error*100.0);
	    }
	}
    }
}
