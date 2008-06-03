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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//  TRD 6-chambers stack                                                     //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

#include "AliRun.h"
#include "AliLog.h"

#include "AliTRDgeometry.h"
#include "AliTRDmodule.h"
#include "AliTRDltuTracklet.h"
#include "AliTRDgtuTrack.h"
#include "AliTRDtrigParam.h"
#include "AliTRDzmaps.h"

ClassImp(AliTRDmodule)

//_____________________________________________________________________________
AliTRDmodule::AliTRDmodule()
  :TObject()
  ,fXprojPlane(0)
  ,fField(0)
  ,fTracklets(new TObjArray(400))
  ,fTracks(new TObjArray(400))
  ,fDeltaY(0)
  ,fDeltaS(0)
  ,fLTUtrk(0)
  ,fGTUtrk(0)
{
  //
  // AliTRDmodule default constructor
  //

  fXprojPlane = AliTRDtrigParam::Instance()->GetXprojPlane();
  fDeltaY     = AliTRDtrigParam::Instance()->GetDeltaY();
  fDeltaS     = AliTRDtrigParam::Instance()->GetDeltaS();

  // The magnetic field strength
  Double_t x[3] = { 0.0, 0.0, 0.0 };
  Double_t b[3];
  gAlice->Field(x,b);  // b[] is in kilo Gauss
  fField = b[2] * 0.1; // Tesla

  for (Int_t iLayer = 0; iLayer < AliTRDgeometry::Nlayer(); iLayer++) {
    for (Int_t i = 0; i < kNsubZchan; i++) {
      fZnchan[iLayer][i] = 0;
      for (Int_t j = 0; j < kNmaxZchan; j++) {
        fZtrkid[iLayer][j][i] = -1;
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDmodule::AliTRDmodule(const AliTRDmodule &m)
  :TObject(m)
  ,fXprojPlane(m.fXprojPlane)
  ,fField(m.fField)
  ,fTracklets(NULL)
  ,fTracks(NULL)
  ,fDeltaY(m.fDeltaY)
  ,fDeltaS(m.fDeltaS)
  ,fLTUtrk(NULL)
  ,fGTUtrk(NULL)
{
  //
  // AliTRDmodule copy constructor
  //

  for (Int_t iLayer = 0; iLayer < AliTRDgeometry::Nlayer(); iLayer++) {
    for (Int_t i = 0; i < kNsubZchan; i++) {
      ((AliTRDmodule &) m).fZnchan[iLayer][i] = 0;
      for (Int_t j = 0; j < kNmaxZchan; j++) {
        ((AliTRDmodule &) m).fZtrkid[iLayer][j][i] = -1;
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDmodule::~AliTRDmodule()
{
  //
  // Destructor
  //

}

//_____________________________________________________________________________
AliTRDmodule &AliTRDmodule::operator=(const AliTRDmodule &m)
{
  //
  // Assignment operator
  //

  if (this != &m) ((AliTRDmodule &) m).Copy(*this); 
  return *this;

}

//_____________________________________________________________________________
void AliTRDmodule::Copy(TObject &m) const
{
  //
  // copy function
  //

  ((AliTRDmodule &) m).fXprojPlane = fXprojPlane;
  ((AliTRDmodule &) m).fField      = fField;
  ((AliTRDmodule &) m).fTracklets  = NULL;
  ((AliTRDmodule &) m).fTracks     = NULL;
  ((AliTRDmodule &) m).fDeltaY     = fDeltaY;
  ((AliTRDmodule &) m).fDeltaS     = fDeltaS;
  ((AliTRDmodule &) m).fLTUtrk     = NULL;
  ((AliTRDmodule &) m).fGTUtrk     = NULL;

  for (Int_t iLayer = 0; iLayer < AliTRDgeometry::Nlayer(); iLayer++) {
    for (Int_t i = 0; i < kNsubZchan; i++) {
      ((AliTRDmodule &) m).fZnchan[iLayer][i] = 0;
      for (Int_t j = 0; j < kNmaxZchan; j++) {
        ((AliTRDmodule &) m).fZtrkid[iLayer][j][i] = -1;
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDmodule::Reset() 
{
  //
  // Reset the tracks and tracklets in the module
  //

  ResetTracklets();
  ResetTracks();

  fLTUtrk    = 0;
  fGTUtrk    = 0;
  fTracklets = new TObjArray(400);
  fTracks    = new TObjArray(400);

}

//_____________________________________________________________________________
void AliTRDmodule::ResetTracks() 
{
  // 
  // Reset the tracks in the module
  //

  if (fTracks) {

    AliTRDgtuTrack *trk;
    for (Int_t i = 0; i < GetNtracks(); i++) {
      
      trk = GetTrack(i);
      trk->Reset();
      
    }

    fTracks->Delete();

  }

}

//_____________________________________________________________________________
AliTRDgtuTrack *AliTRDmodule::GetTrack(Int_t pos) const
{
  //
  // Return track at position "pos"
  //

  if (fTracks == 0) {
    return 0;
  }

  void *trk = fTracks->UncheckedAt(pos);
  if (trk == 0) {
    return 0;
  }

  return (AliTRDgtuTrack *) trk;

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveTrack(Int_t pos)
{
  //
  // Remove the track at position "pos"
  //

  if (fTracks == 0) {
    return;
  }

  fTracks->RemoveAt(pos);
  fTracks->Compress();

}

//_____________________________________________________________________________
void AliTRDmodule::AddTracklet(Int_t det, Int_t row, Float_t rowz, Float_t slope 
			     , Float_t offset, Float_t time, Int_t ncl
                             , Int_t label, Float_t q) 
{
  // 
  // Add a tracklet to this track
  //
  
  fLTUtrk = new AliTRDltuTracklet(det,row,rowz,slope,offset,time,ncl,label,q);
  Tracklets()->Add(fLTUtrk);

}

//_____________________________________________________________________________
AliTRDltuTracklet *AliTRDmodule::GetTracklet(Int_t pos) const
{
  //
  // Get the tracklet at position "pos"
  //

  if (fTracklets == 0) {
    return 0;
  }

  void *trk = fTracklets->UncheckedAt(pos);
  if (trk == 0) {
    return 0;
  }

  return (AliTRDltuTracklet *) trk;

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveTracklet(Int_t pos)
{
  //
  // Remove the tracklet at position "pos"
  //

  if (fTracklets == 0) {
    return;
  }

  fTracklets->RemoveAt(pos);
  fTracklets->Compress();

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveMultipleTracklets()
{
  //
  // Remove multiple found tracklets
  //

  Float_t offDiffMin = 0.5;  // [cm]

  AliTRDltuTracklet *trk;
  Int_t   det1, det2, row1, row2, ncl1, ncl2, label1, label2;
  Float_t off1, off2;
  Int_t   itrk = 0;
  while (itrk < (GetNtracklets() - 1)) {

    trk    = GetTracklet(itrk);
    det1   = trk->GetDetector();
    row1   = trk->GetRow();
    off1   = trk->GetOffset();
    ncl1   = trk->GetNclusters();
    label1 = trk->GetLabel();

    trk    = GetTracklet(itrk+1);
    det2   = trk->GetDetector();
    row2   = trk->GetRow();
    off2   = trk->GetOffset();
    ncl2   = trk->GetNclusters();
    label2 = trk->GetLabel();

    if ((det1 == det2) && (row1 == row2)) {
      if ((off2 - off1) < offDiffMin) {
	if (ncl1 < ncl2) {
	  RemoveTracklet(itrk  );
	}    
        else {
	  RemoveTracklet(itrk+1);
	}
      }
    }

    itrk++;

  }

}

//_____________________________________________________________________________
void AliTRDmodule::SortZ(Int_t cha)
{
  //
  // Match tracklets in the x-z plane (pad row sorting)
  //

  InitZLUT();

  AliTRDltuTracklet *trk;
  Int_t row, pla, det;

  for (Int_t iTrk = 0; iTrk < GetNtracklets(); iTrk++) {

    trk = GetTracklet(iTrk);
    row = trk->GetRow();
    det = trk->GetDetector();
    pla = trk->GetPlane(det);

    for (Int_t iZchan = 0; iZchan < kNsubZchan; iZchan++) {
      if (fZChannelMap[cha][iZchan][pla][row] == 1) {
	fZtrkid[pla][fZnchan[pla][iZchan]][iZchan] = iTrk;
	fZnchan[pla][iZchan]++;
      }
    }

  }

}

//_____________________________________________________________________________
void AliTRDmodule::InitZLUT()
{
  //
  // Initialize the pad row sorting look-up-table
  //

  for (Int_t iLayer = 0; iLayer < AliTRDgeometry::Nlayer(); iLayer++) {
    for (Int_t i = 0; i < kNsubZchan; i++) {
      fZnchan[iLayer][i] = 0;
      for (Int_t j = 0; j < kNmaxZchan; j++) {
        fZtrkid[iLayer][j][i] = -1;
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDmodule::FindTracks() 
{
  //
  // Find tracks from tracklets
  //

  for (Int_t iZchan = 0; iZchan < kNsubZchan; iZchan++) {
    FindTracksCombi(iZchan);
  }

}

//_____________________________________________________________________________
void AliTRDmodule::FindTracksCombi(Int_t zchan) 
{
  //
  // Find tracks by pure combinatorics...
  //
  
  static Int_t trkTrack[12];
  
  Int_t   nTracklets;
  Int_t   nPlanes;
  Int_t   ntrk1;
  Int_t   trkId1;
  Int_t   ntrk2;
  Int_t   trkId2;

  Float_t y1;
  Float_t y1min;
  Float_t y1max;
  Float_t s1;
  Float_t z1;
  Float_t s1min;
  Float_t s1max;
  Float_t y2;
  Float_t s2;
  Float_t z2;

  AliTRDltuTracklet *trk1;
  AliTRDltuTracklet *trk2;
  AliTRDltuTracklet *trk ;

  Bool_t isPlane[kNplan];

  for (Int_t iPlan1 = 0; iPlan1 < kNplan; iPlan1++) {

    ntrk1 = fZnchan[iPlan1][zchan];

    for (Int_t iTrk1 = 0; iTrk1 < ntrk1; iTrk1++) {

      for (Int_t iPlan = 0; iPlan < kNplan; iPlan++) {
        isPlane[iPlan] = kFALSE;
      }

      trkId1 = fZtrkid[iPlan1][iTrk1][zchan];

      nTracklets = 0;
      for (Int_t iList = 0; iList < kNmaxTrk; iList++) {
	trkTrack[iList] = -1;
      }
      trkTrack[nTracklets++] = trkId1;

      isPlane[iPlan1] = kTRUE;

      trk1  = GetTracklet(trkId1);
      y1    = trk1->GetYproj(fXprojPlane);
      y1min = y1 - fDeltaY;
      y1max = y1 + fDeltaY;
      s1    = trk1->GetSlope();
      s1min = s1 - fDeltaS;
      s1max = s1 + fDeltaS;
      z1    = trk1->GetZproj(fXprojPlane);      

      for (Int_t iPlan2 = 0; iPlan2 < kNplan; iPlan2++) {

	if (iPlan2 == iPlan1) continue;

	ntrk2 = fZnchan[iPlan2][zchan];

	for (Int_t iTrk2 = 0; iTrk2 < ntrk2; iTrk2++) {

	  trkId2 = fZtrkid[iPlan2][iTrk2][zchan];

	  if (trkId2 == trkId1) continue;

	  trk2 = GetTracklet(trkId2);
	  y2   = trk2->GetYproj(fXprojPlane);
	  s2   = trk2->GetSlope();
	  z2   = trk2->GetZproj(fXprojPlane);

	  if ((y1min < y2 && y2 < y1max) && 
	      (s1min < s2 && s2 < s1max)) {

	    if (nTracklets >= kNmaxTrk) {
	      AliWarning("Too many tracklets for this track.");
	    }    
            else {
	      trkTrack[nTracklets++] = trkId2;
	      isPlane[iPlan2] = kTRUE;
	    }

	  }

	}  // end trk 2

      }  // end plan 2

      nPlanes = 0;
      for (Int_t iPlan = 0; iPlan < kNplan; iPlan++) {
	nPlanes += (Int_t) isPlane[iPlan];
      }
      
      if (nPlanes >= 4) {

	Int_t cha1, cha2, npoints1, npoints2;
	for (Int_t iList = 0; iList < (nTracklets - 1); iList++) {

	  if (trkTrack[iList] == -1 || trkTrack[iList+1] == -1) continue;
	  trk1 = GetTracklet(trkTrack[iList  ]);
	  trk2 = GetTracklet(trkTrack[iList+1]);

	  cha1 = trk1->GetDetector();
	  cha2 = trk2->GetDetector();
	  if (cha1 != cha2) continue;

	  npoints1 = trk1->GetNclusters();
	  npoints2 = trk2->GetNclusters();

	  if (npoints1 == npoints2) {
	    trkTrack[iList] = -1;
	  } 
          else {
	    if (npoints1 > npoints2) trkTrack[iList+1] = -1;
	    if (npoints1 < npoints2) trkTrack[iList  ] = -1;
	  }

	}

	fGTUtrk = new AliTRDgtuTrack();
	for (Int_t iList = 0; iList < nTracklets; iList++) {
	  if (trkTrack[iList] == -1) continue;
	  trk = GetTracklet(trkTrack[iList]);
	  fGTUtrk->AddTracklet(trk);
	}
	fGTUtrk->Track(fXprojPlane,fField);
	AddTrack();

      }
           
    }  // end trk 1

  }  // end plan 1

}

//_____________________________________________________________________________
void AliTRDmodule::AddTrack() 
{
  //
  // Add a found track to the module
  //
  
  Tracks()->Add(fGTUtrk);

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveMultipleTracks()
{
  //
  // Remove multiple found tracks
  //

  AliTRDgtuTrack *trk1;
  AliTRDgtuTrack *trk2;

  Float_t yproj1;
  Float_t yproj2;
  Float_t alpha1;
  Float_t alpha2;
  Int_t   ntrk1;
  Int_t   ntrk2;
  Int_t   iTrack = 0;

  while (iTrack < (GetNtracks()-1)) {

    trk1   = GetTrack(iTrack  );
    trk2   = GetTrack(iTrack+1);

    ntrk1  = trk1->GetNtracklets();
    yproj1 = trk1->GetYproj();
    alpha1 = trk1->GetSlope();
    ntrk2  = trk2->GetNtracklets();
    yproj2 = trk2->GetYproj();
    alpha2 = trk2->GetSlope();

    if ((TMath::Abs(yproj1-yproj2) < fDeltaY) && 
        (TMath::Abs(alpha1-alpha2) < fDeltaS)) {
      if (ntrk1 < ntrk2) {
	RemoveTrack(iTrack  );
      } 
      else {
	RemoveTrack(iTrack+1);
      }
    } 
    else {
      iTrack++;
    }
    
  }

}

//_____________________________________________________________________________
TObjArray *AliTRDmodule::Tracklets() 
{ 
  //
  // Returns the list of tracklets
  //

  if (!fTracklets) {
    fTracklets = new TObjArray(400); 
  }

  return fTracklets; 

}

//_____________________________________________________________________________
void AliTRDmodule::ResetTracklets() 
{
  //
  // Resets the list of tracklets
  //
 
  if (fTracklets) {
    fTracklets->Delete();
  } 

}

//_____________________________________________________________________________
void AliTRDmodule::SortTracklets()  
{
  //
  // Sorts the list of tracklets
  //

  if (fTracklets) {
    fTracklets->Sort();
  }
 
}

//_____________________________________________________________________________
Int_t AliTRDmodule::GetNtracklets() const 
{
  //
  // Returns the number of tracklets
  //

  if (fTracklets) {
    return fTracklets->GetEntriesFast();
  }

  return 0;

}

//_____________________________________________________________________________
TObjArray *AliTRDmodule::Tracks() 
{
  //
  // Returns the list of tracks
  //
 
  if (!fTracks) {
    fTracks = new TObjArray(400);
  }
 
  return fTracks; 

}

//_____________________________________________________________________________
void AliTRDmodule::SortTracks()  
{ 
  //
  // Sort the list of tracks
  //

  if (fTracks) {
    fTracks->Sort();
  } 

}

//_____________________________________________________________________________
Int_t AliTRDmodule::GetNtracks() const 
{
  //
  // Returns the number of tracks
  //

  if (fTracks) {
    return fTracks->GetEntriesFast();
  }

  return 0;

}
