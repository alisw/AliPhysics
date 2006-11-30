// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group


#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TFile.h>

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTDigitData.h"
#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTTrackArray.h"
#include "AliHLTMemHandler.h"
#include "AliHLTHoughDisplay.h"
#include "AliHLTDigitData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// Display class for Hough transform code

ClassImp(AliHLTHoughDisplay)


AliHLTHoughDisplay::AliHLTHoughDisplay()
{
  //default ctor
  fTracks = 0;
  fDigitRowData = 0;
  fNDigitRowData = 0;
  fShowSlice = -1;
  fPatch = -1;
}

AliHLTHoughDisplay::~AliHLTHoughDisplay()
{
  //dtor
  if(fTracks)
    delete fTracks;
}

void AliHLTHoughDisplay::Init(Char_t *trackfile, Char_t *gfile)
{
  //Init hough display
  TFile *file = TFile::Open(gfile);
  if(!file->IsOpen())
    cerr<<"AliHLTHoughDisplay::AliHLTHoughDisplay : Geometry file " << gfile << " does not exist"<<endl;
  fGeom = (TGeometry*)file->Get("AliceGeom");
  file->Close();
  
  fTracks = new AliHLTTrackArray();
  AliHLTMemHandler *tfile = new AliHLTMemHandler();
  tfile->SetBinaryInput(trackfile);
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
}

void AliHLTHoughDisplay::GenerateHits(AliHLTTrack *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n)
{
  //Generate hits according to the track parameters
  n=0;
  Float_t xyz[3];
  for(Int_t i=AliHLTTransform::GetFirstRow(0); i<AliHLTTransform::GetLastRow(5); i++)
    {
      if(track->GetCrossingPoint(i,xyz))
	{
	  AliHLTTransform::Local2Global(xyz,0);
	  x[n] = xyz[0];
	  y[n] = xyz[1];
	  z[n] = xyz[2];
	  n++;
	}
      else
	break;
    }
}

TPolyMarker3D *AliHLTHoughDisplay::LoadDigits()
{
  //Load digits  
  AliHLTDigitRowData *tempPt = fDigitRowData;
  if(!tempPt)
    {
      cerr<<"AliHLTHoughDisplay::LoadDigits : No data"<<endl;
      return 0;
    }
  
  UInt_t nrows = AliHLTTransform::GetNRows(fPatch);
  Int_t count=0;
  for(UInt_t i=0; i<nrows; i++)
    {
      count += tempPt->fNDigit;
      AliHLTMemHandler::UpdateRowPointer(tempPt);
    }
  tempPt = fDigitRowData;
  TPolyMarker3D *pm = new TPolyMarker3D(count);
  Float_t xyz[3];
  Int_t sector,row;
  count=0;
  for(UInt_t i=0; i<nrows; i++)
    {
      AliHLTDigitData *digPt = tempPt->fDigitData;
      Int_t padrow = (Int_t)tempPt->fRow;
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  AliHLTTransform::Slice2Sector(fShowSlice,padrow,sector,row);
	  AliHLTTransform::Raw2Global(xyz,sector,row,(Int_t)digPt->fPad,(Int_t)digPt->fTime);
	  pm->SetPoint(count,xyz[0],xyz[1],xyz[2]);
	  count++;
	}
      AliHLTMemHandler::UpdateRowPointer(tempPt);
    }

  cout<<"Displaying "<<count<<" digits"<<endl;
  return pm;
}

void AliHLTHoughDisplay::DisplayEvent()
{
  //Display the found tracks.
  
  if(!fTracks)
    {
      cerr<<"AliHLTHoughDisplay::DisplayTracks() : No tracks"<<endl;
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
      AliHLTTrack *track = fTracks->GetCheckedTrack(j); 
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
  

