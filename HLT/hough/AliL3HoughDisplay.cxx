// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TFile.h>

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3HoughTrack.h"
#include "AliL3TrackArray.h"
#include "AliL3MemHandler.h"
#include "AliL3HoughDisplay.h"

#if __GNUC__ == 3
using namespace std;
#endif

//_____________________________________________________________
// Display class for Hough transform code

ClassImp(AliL3HoughDisplay)


AliL3HoughDisplay::AliL3HoughDisplay()
{
  //default ctor
  fTracks = 0;
  fDigitRowData = 0;
  fNDigitRowData = 0;
  fShowSlice = -1;
  fPatch = -1;
}

AliL3HoughDisplay::~AliL3HoughDisplay()
{
  //dtor
  if(fTracks)
    delete fTracks;
}

void AliL3HoughDisplay::Init(Char_t *trackfile, Char_t *gfile)
{
  //Init hough display
  TFile *file = TFile::Open(gfile);
  if(!file->IsOpen())
    cerr<<"AliL3HoughDisplay::AliL3HoughDisplay : Geometry file " << gfile << " does not exist"<<endl;
  fGeom = (TGeometry*)file->Get("AliceGeom");
  file->Close();
  
  fTracks = new AliL3TrackArray();
  AliL3MemHandler *tfile = new AliL3MemHandler();
  tfile->SetBinaryInput(trackfile);
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
}

void AliL3HoughDisplay::GenerateHits(AliL3Track *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n)
{
  //Generate hits according to the track parameters
  n=0;
  Float_t xyz[3];
  for(Int_t i=AliL3Transform::GetFirstRow(0); i<AliL3Transform::GetLastRow(5); i++)
    {
      if(track->GetCrossingPoint(i,xyz))
	{
	  AliL3Transform::Local2Global(xyz,0);
	  x[n] = xyz[0];
	  y[n] = xyz[1];
	  z[n] = xyz[2];
	  n++;
	}
      else
	break;
    }
}

TPolyMarker3D *AliL3HoughDisplay::LoadDigits()
{
  //Load digits  
  AliL3DigitRowData *tempPt = fDigitRowData;
  if(!tempPt)
    {
      cerr<<"AliL3HoughDisplay::LoadDigits : No data"<<endl;
      return 0;
    }
  
  UInt_t nrows = AliL3Transform::GetNRows(fPatch);
  Int_t count=0;
  for(UInt_t i=0; i<nrows; i++)
    {
      count += tempPt->fNDigit;
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  tempPt = fDigitRowData;
  TPolyMarker3D *pm = new TPolyMarker3D(count);
  Float_t xyz[3];
  Int_t sector,row;
  count=0;
  for(UInt_t i=0; i<nrows; i++)
    {
      AliL3DigitData *digPt = tempPt->fDigitData;
      Int_t padrow = (Int_t)tempPt->fRow;
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  AliL3Transform::Slice2Sector(fShowSlice,padrow,sector,row);
	  AliL3Transform::Raw2Global(xyz,sector,row,(Int_t)digPt->fPad,(Int_t)digPt->fTime);
	  pm->SetPoint(count,xyz[0],xyz[1],xyz[2]);
	  count++;
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }

  cout<<"Displaying "<<count<<" digits"<<endl;
  return pm;
}

void AliL3HoughDisplay::DisplayEvent()
{
  //Display the found tracks.
  
  if(!fTracks)
    {
      cerr<<"AliL3HoughDisplay::DisplayTracks() : No tracks"<<endl;
      return;
    }
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  c1->cd();
  
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);

  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.);
  c1->SetPhi(0.);
    
  Int_t ntracks = fTracks->GetNTracks();
  TPolyLine3D *line = new TPolyLine3D[ntracks];
  
  Int_t n;
  Float_t x[176],y[176],z[176];
  for(Int_t j=0; j<ntracks; j++)
    {
      AliL3Track *track = fTracks->GetCheckedTrack(j); 
      if(!track) continue;        
      track->CalculateHelix();
      GenerateHits(track,x,y,z,n);
      TPolyMarker3D *pm = new TPolyMarker3D(n);
      
      for(Int_t h=0; h<n; h++)
	pm->SetPoint(h,x[h],y[h],z[h]);
      
      pm->SetMarkerColor(2);
      pm->Draw();
      TPolyLine3D *currentline = &(line[j]);
      currentline = new TPolyLine3D(n,x,y,z,"");
      currentline->SetLineColor(4);
      currentline->Draw("same");
      
    }
  
  if(fShowSlice>=0)
    {
      TPolyMarker3D *pm = LoadDigits();
      pm->SetMarkerColor(2);
      pm->Draw("same");
    }
  
  cout<<"Displaying...."<<endl;
  fGeom->Draw("same");
  c1->x3d();
}
  

