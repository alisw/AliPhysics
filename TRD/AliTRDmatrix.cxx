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
$Log$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Contains the pixel information for one TRD chamber                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDmatrix.h"

ClassImp(AliTRDmatrix)

//_____________________________________________________________________________
AliTRDmatrix::AliTRDmatrix():TObject()
{
  //
  // Create a TRD detector matrix
  // 
 
  fRow        = 0;
  fCol        = 0;
  fTime       = 0;
  fPixel      = 0;
  fSector     = 0;
  fChamber    = 0;
  fPlane      = 0;
  fPixelArray = NULL;

}

//_____________________________________________________________________________
AliTRDmatrix::AliTRDmatrix(Int_t nRow, Int_t nCol, Int_t nTime
                         , Int_t iSec, Int_t iCha, Int_t iPla)
             :TObject()
{
  //
  // Create a TRD detector matrix with a given size
  // 

  fRow        = nRow;
  fCol        = nCol;
  fTime       = nTime;
  fPixel      = nRow * nCol * nTime;
  fSector     = iSec;
  fChamber    = iCha;
  fPlane      = iPla;
  fPixelArray = new TObjArray(fPixel);
  for (Int_t iPixel = 0; iPixel < fPixel; iPixel++) {
    AliTRDpixel *pixel = new AliTRDpixel();
    fPixelArray->Add(pixel);
  }

}

//_____________________________________________________________________________
AliTRDmatrix::~AliTRDmatrix()
{

  if (fPixelArray) {
    fPixelArray->Delete();
    delete fPixelArray;
  }

}

//_____________________________________________________________________________
void AliTRDmatrix::AddSignal(Int_t iRow, Int_t iCol, Int_t iTime, Float_t signal)
{
  //
  // Add a value to the amplitude of the signal for one specific pixel
  //

  AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
  if (pixel) {
    signal += pixel->GetSignal();
    pixel->SetSignal(signal);
  }

}

//_____________________________________________________________________________
void AliTRDmatrix::Draw()
{
  //
  // Draws a 3D view of the detector matrix
  //

  Char_t ctitle[50];
  sprintf(ctitle,"Matrix (Sector:%d Chamber:%d Plane:%d)"
                ,fSector,fChamber,fPlane);
  TH3F *hMatrix = new TH3F("hMatrix",ctitle,fRow ,-0.5,fRow +0.5
                                           ,fCol ,-0.5,fCol +0.5
                                           ,fTime,-0.5,fTime+0.5);

  for (Int_t iRow  = 0; iRow  < fRow;  iRow++ ) {
    for (Int_t iCol  = 0; iCol  < fCol;  iCol++ ) {
      for (Int_t iTime = 0; iTime < fTime; iTime++) {
        AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
        if (pixel) hMatrix->Fill3(iRow,iCol,iTime,pixel->GetSignal());
      }
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *cMatrix = new TCanvas("cMatrix","Detector matrix 3D-view"
                                          ,50,50,600,400);
  cMatrix->ToggleEventStatus();
  hMatrix->SetXTitle("Pad-row (z)");
  hMatrix->SetYTitle("Pad-column (rphi)");
  hMatrix->SetZTitle("Timebucket");
  hMatrix->Draw("BOX");

}

//_____________________________________________________________________________
void AliTRDmatrix::DrawRow(Int_t iRow)
{
  //
  // Draws a 2D slice of the detector matrix along one row
  //

  if ((iRow < 0) || (iRow >= fRow)) {
    printf("Index out of bounds (%d/%d)\n",iRow,fRow);
    return;
  }

  Char_t ctitle[50];
  sprintf(ctitle,"Pad-row %d (Sector:%d Chamber:%d Plane:%d)"
                ,iRow,fSector,fChamber,fPlane);
  TH2F *hSliceRow = new TH2F("hSliceRow",ctitle,fCol ,-0.5,fCol +0.5
                                               ,fTime,-0.5,fTime+0.5);

  for (Int_t iCol  = 0; iCol  < fCol;  iCol++ ) {
    for (Int_t iTime = 0; iTime < fTime; iTime++) {
      AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
      if (pixel) hSliceRow->Fill(iCol,iTime,pixel->GetSignal());
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *cSliceRow = new TCanvas("cSliceRow","Detector matrix 2D-slice"
                                              ,50,50,600,400);
  cSliceRow->ToggleEventStatus();
  hSliceRow->SetXTitle("Pad-column (rphi)");
  hSliceRow->SetYTitle("Timebucket");
  hSliceRow->Draw("COLZ");

}

//_____________________________________________________________________________
void AliTRDmatrix::DrawCol(Int_t iCol)
{
  //
  // Draws a 2D slice of the detector matrix along one column
  //

  if ((iCol < 0) || (iCol >= fCol)) {
    printf("Index out of bounds (%d/%d)\n",iCol,fCol);
    return;
  }

  Char_t ctitle[50];
  sprintf(ctitle,"Pad-column %d (Sector:%d Chamber:%d Plane:%d)"
                ,iCol,fSector,fChamber,fPlane);
  TH2F *hSliceCol = new TH2F("hSliceCol",ctitle,fRow ,-0.5,fRow +0.5
                                               ,fTime,-0.5,fTime+0.5);

  for (Int_t iRow  = 0; iRow  < fRow;  iRow++ ) {
    for (Int_t iTime = 0; iTime < fTime; iTime++) {
      AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
      if (pixel) hSliceCol->Fill(iRow,iTime,pixel->GetSignal());
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *cSliceCol = new TCanvas("cSliceCol","Detector matrix 2D-slice"
                                              ,50,50,600,400);
  cSliceCol->ToggleEventStatus();
  hSliceCol->SetXTitle("Pad-row (z)");
  hSliceCol->SetYTitle("Timebucket");
  hSliceCol->Draw("COLZ");

}

//_____________________________________________________________________________
void AliTRDmatrix::DrawTime(Int_t iTime)
{
  //
  // Draws a 2D slice of the detector matrix along one time slice
  //

  if ((iTime < 0) || (iTime >= fTime)) {
    printf("Index out of bounds (%d/%d)\n",iTime,fTime);
    return;
  }

  Char_t ctitle[50];
  sprintf(ctitle,"Time-slice %d (Sector:%d Chamber:%d Plane:%d)"
                ,iTime,fSector,fChamber,fPlane);
  TH2F *hSliceTime = new TH2F("hSliceTime",ctitle,fRow,-0.5,fRow+0.5
                                                 ,fCol,-0.5,fCol+0.5);

  for (Int_t iRow = 0; iRow < fRow; iRow++) {
    for (Int_t iCol = 0; iCol < fCol; iCol++) {
      AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
      if (pixel) hSliceTime->Fill(iRow,iCol,pixel->GetSignal());
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *cSliceTime = new TCanvas("cSliceTime","Detector matrix 2D-slice"
                                                ,50,50,600,400);
  cSliceTime->ToggleEventStatus();
  hSliceTime->SetXTitle("Pad-row (z)");
  hSliceTime->SetYTitle("Pad-column (rphi)");
  hSliceTime->Draw("COLZ");

}

//_____________________________________________________________________________
void AliTRDmatrix::SetSignal(Int_t iRow, Int_t iCol, Int_t iTime, Float_t signal)
{
  //
  // Set the amplitude of the signal for one specific pixel
  //

  AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
  if (pixel) {
    pixel->SetSignal(signal);
  }

}

//_____________________________________________________________________________
Bool_t AliTRDmatrix::AddTrack(Int_t iRow, Int_t iCol, Int_t iTime, Int_t track)
{
  // 
  // Add this track number to the stored tracks passing through this pixel. 
  // If there are already three stored the return status is FALSE.  
  //

  AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
  if (!(pixel)) return kTRUE;

  Bool_t trackSet = kFALSE;
  for (Int_t i = 0; i < kTrackPixel; i++) {
    if (pixel->GetTrack(i) == track) {
      trackSet = kTRUE;
      break;
    }
    if (pixel->GetTrack(i) ==     0) {
      pixel->SetTrack(i,track);
      trackSet = kTRUE;
      break;
    }
  }    

  return trackSet;

}

//_____________________________________________________________________________
void AliTRDmatrix::SetTrack(Int_t iRow, Int_t iCol, Int_t iTime
                          , Int_t iTrack, Int_t track)
{
  //
  // Store the number of a track which is passing through this pixel
  //

  AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
  if (pixel) {
    pixel->SetTrack(iTrack,track);
  }

}

//_____________________________________________________________________________
Float_t AliTRDmatrix::GetSignal(Int_t iRow, Int_t iCol, Int_t iTime)
{
  //
  // Returns the amplitude of the signal for one specific pixel
  //

  AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
  if (pixel) {
    return (pixel->GetSignal());
  }
  else {
    return 0;
  }

}

//_____________________________________________________________________________
Int_t AliTRDmatrix::GetTrack(Int_t iRow, Int_t iCol, Int_t iTime, Int_t iTrack)
{
  //
  // Returns the numbers of the tracks passing through one specific pixel
  //

  if ((iTrack < 0) || (iTrack >= kTrackPixel)) {
    printf("GetTrack: Index out of bounds (%d)\n",iTrack);
    return 0;
  }

  AliTRDpixel *pixel = GetPixel(iRow,iCol,iTime);
  if (pixel) {
    return (pixel->GetTrack(iTrack));
  }
  else {
    return 0;
  }

}

//_____________________________________________________________________________
Int_t AliTRDmatrix::GetIndex(Int_t iRow, Int_t iCol, Int_t iTime)
{

  if ((iRow  >= 0) && (iRow  < fRow ) &&
      (iCol  >= 0) && (iCol  < fCol ) &&
      (iTime >= 0) && (iTime < fTime)) {
    return (iTime + iCol * fTime + iRow * fTime * fCol);
  }
  else {
    return -1;
  }

}

//_____________________________________________________________________________
AliTRDpixel *AliTRDmatrix::GetPixel(Int_t iRow, Int_t iCol, Int_t iTime)
{

  Int_t iPixel = GetIndex(iRow,iCol,iTime);
  if (iPixel < 0) {
    return NULL;
  }
  else {
    return ((AliTRDpixel *) fPixelArray->At(iPixel));
  }

}
