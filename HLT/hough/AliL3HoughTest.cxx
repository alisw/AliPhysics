//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 

#include "AliL3StandardIncludes.h"
#include "AliL3HoughTest.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"
#include "AliL3Histogram.h"
#include "TRandom.h"
#include "TMath.h"
#include "TH2.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTest


ClassImp(AliL3HoughTest)

AliL3HoughTest::AliL3HoughTest()
{
  fData=0;
}


AliL3HoughTest::~AliL3HoughTest()
{
  if(fData)
    delete [] fData;
}

void AliL3HoughTest::GenerateTrackData(Double_t pt,Double_t psi,Int_t patch)
{
  fCurrentPatch=patch;
  if(fData)
    delete fData;
  fData = new SimData[AliL3Transform::GetNRows(patch)];
  memset(fData,0,AliL3Transform::GetNRows(patch)*sizeof(SimData));
  
  AliL3ModelTrack *track = new AliL3ModelTrack();
  track->Init(0,patch);
  track->SetPt(pt);
  track->SetPsi(psi);
  track->SetCharge(1);
  track->CalculateHelix();
  

  Int_t temp[200];
  Int_t entries=100;
  Int_t clustercharge=100;
  for(Int_t i=AliL3Transform::GetFirstRow(patch); i<=AliL3Transform::GetLastRow(patch); i++)
    {
      Float_t xyz[3];
      
      if(!track->GetCrossingPoint(i,xyz))
	continue;

      Int_t sector,row;
      AliL3Transform::Slice2Sector(0,i,sector,row);
      AliL3Transform::Local2Raw(xyz,sector,row);

      if(xyz[1] < 0 || xyz[1] >= AliL3Transform::GetNPads(i))
	continue;
      track->SetPadHit(row,xyz[1]);

      memset(temp,0,200*sizeof(Int_t));
      Double_t xysigma = sqrt(track->GetParSigmaY2(i));
      Int_t minpad=200,j;
      for(j=0; j<entries; j++)
	{
	  Int_t pad = TMath::Nint(gRandom->Gaus(xyz[1],xysigma));
	  if(pad < 0 || pad >= AliL3Transform::GetNPads(i))
	    continue;
	  temp[pad]++;
	  if(pad < minpad)
	    minpad=pad;
	}
      for(j=0; j<200; j++)
	{
	  if(temp[j]==0) continue;
	  Int_t index = j - minpad;
	  if(index < 0 || index >= 10)
	    {
	      cerr<<"AliL3HoughTest::GenerateTrackData : Wrong index "<<index<<endl;
	      exit(5);
	    }
	  Int_t charge = clustercharge*temp[j]/entries;
	  if(charge < 3 ) 
	    continue;
	  if(charge > 1023)
	    charge=1023;
	  fData[i].pads[index]=charge;
	  fData[i].npads++;
	}
      fData[i].minpad=minpad;
    }
  delete track;
}

void AliL3HoughTest::Transform(AliL3Histogram *hist)
{
  if(!fData)
    {
      cerr<<"AliL3HoughTest::Transform : No data"<<endl;
      return;
    }
  Float_t R,phi,phi0,kappa,xyz[3];
  Int_t pad,charge,sector,row;
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      AliL3Transform::Slice2Sector(0,i,sector,row);
      for(Int_t j=0; j<fData[i].npads; j++)
	{
	  pad = j + fData[i].minpad;
	  charge = fData[i].pads[j];
	  
	  AliL3Transform::Raw2Local(xyz,sector,row,pad,0);

	  R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);

	  phi = AliL3Transform::GetPhi(xyz);
	  
	  for(Int_t k=hist->GetFirstYbin(); k<=hist->GetLastYbin(); k++)
	    {
	      phi0 = hist->GetBinCenterY(k);
	      kappa = 2*sin(phi-phi0)/R;
	      hist->Fill(kappa,phi0,charge);
	    }
	}
    }
}

void AliL3HoughTest::TransformC(AliL3Histogram *hist)
{
  if(!fData)
    {
      cerr<<"AliL3HoughTest::TransformC : No data"<<endl;
      return;
    }
  Int_t pad1,pad2,charge1,charge2,sector,row;
  Float_t r1,r2,phi1,phi2,phi_0,kappa,hit[3],hit2[3];
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      for(Int_t d1=0; d1<fData[i].npads; d1++)
	{
	  pad1 = d1 + fData[i].minpad;
	  charge1 = fData[i].pads[d1];
	  AliL3Transform::Slice2Sector(0,i,sector,row);
	  AliL3Transform::Raw2Local(hit,sector,row,pad1,0);
	  r1 = sqrt(hit[0]*hit[0]+hit[1]*hit[1]);
	  phi1 = atan2(hit[1],hit[0]);
	  
	  for(Int_t j=i+1; j<=AliL3Transform::GetLastRow(fCurrentPatch); j++)
	    {
	      for(Int_t d2=0; d2<fData[j].npads; d2++)
		{
		  pad2 = d2 + fData[j].minpad;
		  charge2 = fData[j].pads[d2];
		  AliL3Transform::Slice2Sector(0,j,sector,row);
		  AliL3Transform::Raw2Local(hit2,sector,row,pad2,0);
		  r2 = sqrt(hit2[0]*hit2[0]+hit2[1]*hit2[1]);
		  phi2 = atan2(hit2[1],hit2[0]);
		  phi_0 = atan( (r2*sin(phi1) - r1*sin(phi2)) / (r2*cos(phi1) - r1*cos(phi2)) );
		  kappa = 2*sin(phi1-phi_0) / r1;
		  hist->Fill(kappa,phi_0,charge1+charge2);
		}
	    }
	}
    }
}

void AliL3HoughTest::FillImage(TH2 *hist)
{
  if(!fData)
    {
      cerr<<"AliL3HoughTest::FillImage : No data to fill"<<endl;
      return;
    }
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      for(Int_t j=0; j<fData[i].npads; j++)
	{
	  Int_t pad = j + fData[i].minpad;
	  Int_t charge = fData[i].pads[j];
	  Float_t xyz[3];
	  Int_t sector,row;
	  AliL3Transform::Slice2Sector(0,i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,pad,0);
	  hist->Fill(xyz[0],xyz[1],charge);
	}
    }
}
