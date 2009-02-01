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
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TH2F.h>

//ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"

//STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliAnalysisManager.h"

//PWG3/muon includes
#include "AliAnalysisTaskMuonTrackingEff.h"
#include "AliCheckMuonDetEltResponse.h"

//MUON includes
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrackExtrap.h"

ClassImp(AliAnalysisTaskMuonTrackingEff)

const Int_t AliAnalysisTaskMuonTrackingEff::fTotNbrOfDetectionElt  = 156;

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff()
    :
    AliAnalysisTask(),
    fTransformer(0x0),
    fESD(0x0),
    fDetEltEffHistList(0x0),
    fDetEltTDHistList(0x0),
    fDetEltTTHistList(0x0)
{
/// Default constructor
}
//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& src)
    :
    AliAnalysisTask(src),
    fTransformer(0x0),
    fESD(0x0),
    fDetEltEffHistList(0x0),
    fDetEltTDHistList(0x0),
    fDetEltTTHistList(0x0)
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
				       const AliMUONGeometryTransformer* transformer)
    :
    AliAnalysisTask(name, "AnalysisTaskESD"),
    fTransformer(transformer),
    fESD(0x0),
    fDetEltEffHistList(0x0),
    fDetEltTDHistList(0x0),
    fDetEltTTHistList(0x0)
{
//Constructor
//-----------
    
//Define detection element efficiency histograms
//----------------------------------------------

    fDetEltEffHistList = new TClonesArray("TH2F",fTotNbrOfDetectionElt + 3); //!<+3 for: 1.the total efficiency chamber by chamber histogram.
                                                                //!<        2.the total number of tracks detected by the tracking system.
                                                                //!<        3.the number of track used for the efficiency calculation chamber by chamber
    fDetEltTDHistList  = new TClonesArray("TH2F",fTotNbrOfDetectionElt + 1);
    fDetEltTTHistList  = new TClonesArray("TH2F",fTotNbrOfDetectionElt + 1);

    for (Int_t i = 0; i<fTotNbrOfDetectionElt; ++i)
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

    new((*fDetEltTDHistList )[fTotNbrOfDetectionElt]) TH2F("TD_Chamber" ,"TD_Chamber" ,10,0,10,1,0,1); //!<Detected tracks.
    new((*fDetEltTTHistList )[fTotNbrOfDetectionElt]) TH2F("TT_Chamber" ,"TT_Chamber" ,10,0,10,1,0,1); //!<Tracks total number.
    new((*fDetEltEffHistList)[fTotNbrOfDetectionElt]) TH2F("fChamberEff","fChamberEff",10,0,10,1,0,1); //!<Chamber efficiency.

    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> GetXaxis() -> SetTitle("Chamber number");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> GetYaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> GetZaxis() -> SetTitle("Efficiency (%)");  
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> GetZaxis() -> SetTitleOffset(1.2); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt)) -> SetOption("LEGO");
  
    new((*fDetEltEffHistList)[157]) TH2F("TT_Chamber" ,"TT_Chamber" ,10,0,10,1,0,1); //!<Tracks total number by chamber.

    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> GetXaxis() -> SetTitle("Chamber number");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> GetYaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> GetZaxis() -> SetTitle("Number of tracks");  
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> GetZaxis() -> SetTitleOffset(1.2); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 1)) -> SetOption("LEGO");

    new((*fDetEltEffHistList)[158]) TH2F("Total_Number_of_Tracks" ,"Total_Number_of_Tracks" ,1,0,1,1,0,1); //!<Tracks total number.

    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> GetXaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> GetYaxis() -> SetTitle("");
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> GetZaxis() -> SetTitle("Number of tracks");  
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> GetXaxis() -> SetTitleOffset(1.8);
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> GetYaxis() -> SetTitleOffset(1.8); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> GetZaxis() -> SetTitleOffset(1.2); 
    ((TH2F*) fDetEltEffHistList->UncheckedAt(fTotNbrOfDetectionElt + 2)) -> SetOption("LEGO");    

//Define input & output
//---------------------

// -Input slot 0 works with a TChain:
    DefineInput(0, TChain::Class());

// -Output slot 0 writes into a TClonesArray:
    DefineOutput(0, TClonesArray::Class());
}



//______________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::~AliAnalysisTaskMuonTrackingEff()
{
// Destructor.
    delete fDetEltEffHistList;
    delete fDetEltTDHistList;
    delete fDetEltTTHistList;
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
    
   
    AliCheckMuonDetEltResponse* chamberEff;
    chamberEff = new AliCheckMuonDetEltResponse(fTransformer, fESD, fDetEltTDHistList, fDetEltTTHistList);
    chamberEff->CheckDetEltResponse();


    for( Int_t i = 0; i<156; ++i)
    {
      ((TH2F*) fDetEltEffHistList->UncheckedAt(i))-> Divide( (TH2F*) fDetEltTDHistList->UncheckedAt(i),
							     (TH2F*) fDetEltTTHistList->UncheckedAt(i), 100., 1.); 
    }

    ((TH2F*) fDetEltEffHistList->UncheckedAt(156))-> Divide( (TH2F*) fDetEltTDHistList->UncheckedAt(156),
							     (TH2F*) fDetEltTTHistList->UncheckedAt(156), 100., 1.); 


    ((TH2F*) fDetEltEffHistList->UncheckedAt(157))-> Add   ( (TH2F*) fDetEltTTHistList ->UncheckedAt(156),
							     (TH2F*) fDetEltEffHistList->UncheckedAt(157), 1., 0.); 


    ((TH2F*) fDetEltEffHistList->UncheckedAt(158))-> Fill(0., 0., ((Double_t)fESD -> GetNumberOfMuonTracks()));


//Post the output data:
    PostData(0, fDetEltEffHistList);  
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Terminate(Option_t */*option*/)
{
//Terminate analysis

}
