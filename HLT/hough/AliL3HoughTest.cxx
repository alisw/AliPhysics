// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"
#include "AliL3HoughTest.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"
#include "AliL3Histogram.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"

#include <TRandom.h>
#include <TMath.h>
#include <TH2.h>
#include <TH3.h>
#include <TAxis.h>

#if __GNUC__ == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTest


ClassImp(AliL3HoughTest)

AliL3HoughTest::AliL3HoughTest()
{
  //ctor
  fData=0;
}


AliL3HoughTest::~AliL3HoughTest()
{
  //dtor
  if(fData)
    delete [] fData;
}

Bool_t AliL3HoughTest::GenerateTrackData(Double_t pt,Double_t psi,Double_t tgl,Int_t sign,Int_t patch,Int_t minhits)
{
  //Generate digits according to given track parameters?
  fCurrentPatch=patch;
  if(fData)
    delete fData;
  fData = new AliL3SimData[AliL3Transform::GetNRows(patch)];
  memset(fData,0,AliL3Transform::GetNRows(patch)*sizeof(AliL3SimData));
  
  AliL3ModelTrack *track = new AliL3ModelTrack();
  track->Init(0,patch);
  track->SetPt(pt);
  track->SetPsi(psi);
  track->SetTgl(tgl);
  track->SetCharge(sign);
  track->SetFirstPoint(0,0,0);
  track->CalculateHelix();

  Int_t temp[200];
  //  Int_t temp2[AliL3Transform::GetNTimeBins()];
  Int_t * temp2 = new Int_t[AliL3Transform::GetNTimeBins()];
  Int_t entries=100;
  Int_t clustercharge=100;
  Int_t hitcounter=0;
  for(Int_t i=AliL3Transform::GetFirstRow(patch); i<=AliL3Transform::GetLastRow(patch); i++)
    {
      
      Float_t xyz[3];
      if(!track->GetCrossingPoint(i,xyz))
	continue;

      Int_t rowindex = i - AliL3Transform::GetFirstRow(patch);
      Int_t sector,row;
      AliL3Transform::Slice2Sector(0,i,sector,row);
      AliL3Transform::Local2Raw(xyz,sector,row);

      if(xyz[1] < 0 || xyz[1] >= AliL3Transform::GetNPads(i) || xyz[2] < 0 || xyz[2] >= AliL3Transform::GetNTimeBins())
	continue;
      hitcounter++;
      track->SetPadHit(i,xyz[1]);
      track->SetTimeHit(i,xyz[2]);

      Double_t beta = track->GetCrossingAngle(i);
      track->SetCrossingAngleLUT(i,beta);
      track->CalculateClusterWidths(i);
      
      memset(temp,0,200*sizeof(Int_t));
      memset(temp2,0,AliL3Transform::GetNTimeBins()*sizeof(Int_t));
      Double_t xysigma = sqrt(track->GetParSigmaY2(i));
      Double_t zsigma = sqrt(track->GetParSigmaZ2(i));
      
      Int_t minpad=200,j;
      Int_t mintime = 1000;
      for(j=0; j<entries; j++)
	{
	  Int_t pad = TMath::Nint(gRandom->Gaus(xyz[1],xysigma));
	  Int_t time = TMath::Nint(gRandom->Gaus(xyz[2],zsigma));
	  if(pad < 0 || pad >= AliL3Transform::GetNPads(i) || time < 0 || time >= AliL3Transform::GetNTimeBins())
	    continue;
	  temp[pad]++;
	  temp2[time]++;
	  if(pad < minpad)
	    minpad=pad;
	  if(time < mintime)
	    mintime=time;
	}
      Int_t npads=0;
      for(j=0; j<200; j++)
	{
	  if(temp[j]==0) continue;
	  
	  Int_t index = j - minpad;
	  
	  if(index < 0 || index >= 10)
	    {
	      cerr<<"AliL3HoughTest::GenerateTrackData : Wrong index "<<index<<endl;
	      exit(5);
	    }
	  npads++;
	  Int_t seqcharge = clustercharge*temp[j]/entries;
	  Int_t ntimes=0;
	  for(Int_t k=0; k<AliL3Transform::GetNTimeBins(); k++)
	    {
	      if(temp2[k]==0) continue;
	      Int_t tindex = k - mintime;
	      if(tindex < 0 || tindex >= 10)
		{
		  cerr<<"AliL3HoughTest::GenerateTrackData : Wrong timeindex "<<tindex<<" "<<k<<" "<<mintime<<endl;
		  exit(5);
		}
	      Int_t charge = seqcharge*temp2[k]/entries;
	      if(charge < 3 ) 
		continue;
	      if(charge > 1023)
		charge=1023;
	      //cout<<"row "<<i<<" pad "<<j<<" time "<<k<<" charge "<<charge<<endl;
	      ntimes++;
	      fData[rowindex].fPads[index][tindex]=charge;
	    }
	}
      fData[rowindex].fMinpad=minpad;
      fData[rowindex].fMintime=mintime;
      fData[rowindex].fNPads=npads;
    }
  delete track;
  delete [] temp2;
  if(hitcounter < minhits)
    return kFALSE;
  return kTRUE;
}

void AliL3HoughTest::Transform2Circle(AliL3Histogram *hist)
{
  //Hough trasnform
  if(!fData)
    {
      cerr<<"AliL3HoughTest::Transform : No data"<<endl;
      return;
    }
  Float_t r,phi,phi0,kappa,xyz[3];
  Int_t pad,time,charge,sector,row;
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      Int_t rowindex = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      AliL3Transform::Slice2Sector(0,i,sector,row);
      for(Int_t j=0; j<fData[rowindex].fNPads; j++)
	{
	  pad = j + fData[rowindex].fMinpad;
	  for(Int_t k=0; k<10; k++)
	    {
	      time = k + fData[rowindex].fMintime;
	      charge = fData[rowindex].fPads[j][k];
	      if(charge == 0) continue;
	      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	      
	      r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	      
	      phi = AliL3Transform::GetPhi(xyz);
	      
	      for(Int_t k=hist->GetFirstYbin(); k<=hist->GetLastYbin(); k++)
		{
		  phi0 = hist->GetBinCenterY(k);
		  kappa = 2*sin(phi-phi0)/r;
		  hist->Fill(kappa,phi0,charge);
		}
	    }
	}
    }
}

void AliL3HoughTest::Transform2CircleC(AliL3Histogram *hist)
{
  //Hough transform
  if(!fData)
    {
      cerr<<"AliL3HoughTest::TransformC : No data"<<endl;
      return;
    }
  Int_t pad1,pad2,time1,time2,charge1,charge2,sector,row;
  Float_t r1,r2,phi1,phi2,phi0,kappa,hit[3],hit2[3];
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      Int_t rowindex1 = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      for(Int_t d1=0; d1<fData[rowindex1].fNPads; d1++)
	{
	  pad1 = d1 + fData[rowindex1].fMinpad;
	  for(Int_t j=0; j<10; j++)
	    {
	      time1 = j + fData[rowindex1].fMintime;
	      charge1 = fData[rowindex1].fPads[d1][j];
	      if(charge1==0) continue;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(hit,sector,row,pad1,time1);
	      r1 = sqrt(hit[0]*hit[0]+hit[1]*hit[1]);
	      phi1 = atan2(hit[1],hit[0]);
	      
	      for(Int_t j=i+1; j<=AliL3Transform::GetLastRow(fCurrentPatch); j++)
		{
		  Int_t rowindex2 = j - AliL3Transform::GetFirstRow(fCurrentPatch);
		  for(Int_t d2=0; d2<fData[rowindex2].fNPads; d2++)
		    {
		      pad2 = d2 + fData[rowindex2].fMinpad;
		      for(Int_t k=0; k<10; k++)
			{
			  time2 = k + fData[rowindex2].fMintime;
			  charge2 = fData[rowindex2].fPads[d2][k];
			  if(charge2==0) continue;
			  AliL3Transform::Slice2Sector(0,j,sector,row);
			  AliL3Transform::Raw2Local(hit2,sector,row,pad2,time2);
			  r2 = sqrt(hit2[0]*hit2[0]+hit2[1]*hit2[1]);
			  phi2 = atan2(hit2[1],hit2[0]);
			  phi0 = atan( (r2*sin(phi1) - r1*sin(phi2)) / (r2*cos(phi1) - r1*cos(phi2)) );
			  
			  kappa = 2*sin(phi1-phi0) / r1;
			  hist->Fill(kappa,phi0,charge1+charge2);
			}
		    }
		}
	      return;
	    }
	}
    }
}

void AliL3HoughTest::Transform2CircleF(AliL3Histogram *hist)
{
  //Fix one point in the middle of the tpc
  
  if(!fData)
    {
      cerr<<"AliL3HoughTest::TransformF : No data"<<endl;
      return;
    }
  Int_t pad1,pad2,time1,time2,charge1,charge2,sector,row;
  Float_t r1,r2,phi1,phi2,phi0,kappa,hit[3],hit2[3];
  
  Int_t rowindex1 = 80 - AliL3Transform::GetFirstRow(fCurrentPatch);
  for(Int_t d1=1; d1<fData[rowindex1].fNPads; d1++)
    {
      pad1 = d1 + fData[rowindex1].fMinpad;
      
      AliL3Transform::Slice2Sector(0,80,sector,row);
      AliL3Transform::Raw2Local(hit,sector,row,pad1,0);
      r1 = sqrt(hit[0]*hit[0]+hit[1]*hit[1]);
      phi1 = atan2(hit[1],hit[0]);
      
      for(Int_t j=0; j<10; j++)
	{
	  time1 = j + fData[rowindex1].fMintime;
	  charge1 = fData[rowindex1].fPads[d1][j];
	  if(charge1==0) continue;
	  
	  for(Int_t j=0; j<=AliL3Transform::GetLastRow(fCurrentPatch); j++)
	    {
	      if(j==80) continue;
	      Int_t rowindex2 = j - AliL3Transform::GetFirstRow(fCurrentPatch);
	      for(Int_t d2=0; d2<fData[rowindex2].fNPads; d2++)
		{
		  pad2 = d2 + fData[rowindex2].fMinpad;
		  for(Int_t k=0; k<10; k++)
		    {
		      time2 = k + fData[rowindex2].fMintime;
		      charge2 = fData[rowindex2].fPads[d2][k];
		      if(charge2==0) continue;
		      AliL3Transform::Slice2Sector(0,j,sector,row);
		      AliL3Transform::Raw2Local(hit2,sector,row,pad2,time2);
		      r2 = sqrt(hit2[0]*hit2[0]+hit2[1]*hit2[1]);
		      phi2 = atan2(hit2[1],hit2[0]);
		      phi0 = atan( (r2*sin(phi1) - r1*sin(phi2)) / (r2*cos(phi1) - r1*cos(phi2)) );
		      
		      kappa = 2*sin(phi1-phi0) / r1;
		      hist->Fill(kappa,phi0,charge1+charge2);
		      //cout<<"Filling "<<kappa<<" psi "<<phi0<<" charge "<<charge1<<endl;
		    }
		}
	    }
	}
      return;
    }
}

void AliL3HoughTest::Transform2Line(AliL3Histogram *hist,Int_t *rowrange)
{
  //Hough transform
  if(!fData)
    {
      cerr<<"AliL3HoughTest::Transform2Line : No data"<<endl;
      return;
    }
  
  Int_t pad,time,charge,sector,row;
  Float_t hit[3],theta,rho;
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      if(i<rowrange[0])
	continue;
      if(i>rowrange[1])
	break;
      Int_t rowindex = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      for(Int_t d=0; d<fData[rowindex].fNPads; d++)
	{
	  pad = d + fData[rowindex].fMinpad;
	  for(Int_t j=0; j<10; j++)
	    {
	      time = j + fData[rowindex].fMintime;
	      charge = fData[rowindex].fPads[d][j];
	      if(charge==0) continue;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(hit,sector,row,pad,time);
	      
	      hit[0] = hit[0] - AliL3Transform::Row2X(rowrange[0]);
	      
	      for(Int_t xbin=hist->GetFirstXbin(); xbin<hist->GetLastXbin(); xbin++)
		{
		  theta = hist->GetBinCenterX(xbin);
		  rho = hit[0]*cos(theta) + hit[1]*sin(theta);
		  hist->Fill(theta,rho,charge);
		}
	    }
	}
    }
}

void AliL3HoughTest::Transform2LineC(AliL3Histogram *hist,Int_t *rowrange)
{
  //Hough transform?
  if(!fData)
    {
      cerr<<"AliL3HoughTest::Transform2Line : No data"<<endl;
      return;
    }
  
  Int_t pad1,pad2,time1,time2,charge1,charge2,sector,row;
  Float_t theta,rho,hit[3],hit2[3];
  for(Int_t i=rowrange[0]; i<=rowrange[1]; i++)
    {
      Int_t rowindex1 = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      for(Int_t d1=0; d1<fData[rowindex1].fNPads; d1++)
	{
	  pad1 = d1 + fData[rowindex1].fMinpad;
	  for(Int_t j=0; j<10; j++)
	    {
	      time1 = j + fData[rowindex1].fMintime;
	      charge1 = fData[rowindex1].fPads[d1][j];
	      if(charge1==0) continue;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(hit,sector,row,pad1,time1);
	      
	      hit[0] = hit[0] - AliL3Transform::Row2X(rowrange[0]);
	      
	      for(Int_t i2=i+1; i2<=rowrange[1]; i2++)
		{
		  Int_t rowindex2 = i2 - AliL3Transform::GetFirstRow(fCurrentPatch);
		  for(Int_t d2=0; d2<fData[rowindex2].fNPads; d2++)
		    {
		      pad2 = d2 + fData[rowindex2].fMinpad;
		      for(Int_t k=0; k<10; k++)
			{
			  time2 = k + fData[rowindex2].fMintime;
			  charge2 = fData[rowindex2].fPads[d2][k];
			  if(charge2==0) continue;
			  AliL3Transform::Slice2Sector(0,i2,sector,row);
			  AliL3Transform::Raw2Local(hit2,sector,row,pad2,time2);
			  
			  hit2[0] = hit2[0] - AliL3Transform::Row2X(rowrange[0]);
			  
			  theta = atan2(hit2[0]-hit[0],hit[1]-hit2[1]);
			  rho = hit[0]*cos(theta)+hit[1]*sin(theta);
			  hist->Fill(theta,rho,charge1+charge2);
			}
		    }
		}
	    }
	}
    }
}

void AliL3HoughTest::FillImage(TH2 *hist,Int_t displayrow)
{
  //Draw digits data
  if(!fData)
    {
      cerr<<"AliL3HoughTest::FillImage : No data to fill"<<endl;
      return;
    }
  
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      Int_t rowindex = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      if(displayrow >=0)
	if(i != displayrow) continue;

      //cout<<"row "<<i<<" fNPads "<<fData[rowindex].fNPads<<endl;
      for(Int_t j=0; j<fData[rowindex].fNPads; j++)
	{
	  Int_t pad = j + fData[rowindex].fMinpad;
	  for(Int_t k=0; k<10; k++)
	    {
	      Int_t time = k + fData[rowindex].fMintime;
	      Int_t charge = fData[rowindex].fPads[j][k];
	      if(charge==0) continue;
	      //cout<<i<<" "<<pad<<" "<<time<<" "<<charge<<endl;
	      Float_t xyz[3];
	      Int_t sector,row;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	      if(displayrow>=0)
		hist->Fill(pad,time,charge);
	      else
		hist->Fill(xyz[0],xyz[1],charge);
	    }
	}
      if(displayrow>=0)
	break;
    }
}

void AliL3HoughTest::Transform2Line3D(TH3 *hist,Int_t *rowrange,Float_t *phirange)
{
  //HT in 3D?
  if(!fData)
    {
      cerr<<"AliL3HoughTest::Transform2Line : No data"<<endl;
      return;
    }
  
  Int_t pad,time,charge,sector,row;
  Float_t hit[3],theta,rho,r,delta;
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      
      if(i<rowrange[0])
	continue;
      if(i>rowrange[1])
	break;
      
      Int_t rowindex = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      for(Int_t d=0; d<fData[rowindex].fNPads; d++)
	{
	  pad = d + fData[rowindex].fMinpad;
	  for(Int_t j=0; j<10; j++)
	    {
	      time = j + fData[rowindex].fMintime;
	      charge = fData[rowindex].fPads[d][j];
	      if(charge==0) continue;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(hit,sector,row,pad,time);
	      
	      Float_t phi = AliL3Transform::GetPhi(hit);
	      if(phi < phirange[0] || phi > phirange[1])
		continue;
	      
	      hit[0] = hit[0] - AliL3Transform::Row2X(rowrange[0]);
	      Float_t x = hit[0] + AliL3Transform::Row2X(rowrange[0]);
	      r = sqrt(x*x + hit[1]*hit[1]);
	      delta = atan(hit[2]/r);
	      
	      for(Int_t xbin=hist->GetXaxis()->GetFirst(); xbin<=hist->GetXaxis()->GetLast(); xbin++)
		{
		  theta = hist->GetXaxis()->GetBinCenter(xbin);
		  rho = hit[0]*cos(theta) + hit[1]*sin(theta);
		  hist->Fill(theta,rho,delta,charge);
		}
	    }
	}
    }
}

void AliL3HoughTest::Transform2LineC3D(TH3 *hist,Int_t *rowrange)
{
  //HT in 3D?
  if(!fData)
    {
      cerr<<"AliL3HoughTest::Transform2Line : No data"<<endl;
      return;
    }
  
  Int_t pad1,pad2,time1,time2,charge1,charge2,sector,row;
  Float_t theta,rho,hit[3],hit2[3],r1,r2,delta,delta1,delta2;
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      if(i<rowrange[0])
	continue;
      if(i>rowrange[1])
	break;
      Int_t rowindex1 = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      for(Int_t d1=0; d1<fData[rowindex1].fNPads; d1++)
	{
	  pad1 = d1 + fData[rowindex1].fMinpad;
	  for(Int_t j=0; j<10; j++)
	    {
	      time1 = j + fData[rowindex1].fMintime;
	      charge1 = fData[rowindex1].fPads[d1][j];
	      if(charge1==0) continue;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(hit,sector,row,pad1,time1);
	      r1 = sqrt(hit[0]*hit[0]+hit[1]*hit[1]);
	      delta1 = atan(hit[2]/r1);
	      hit[0] = hit[0] - AliL3Transform::Row2X(rowrange[0]);
	      
	      for(Int_t i2=i+1; i2<=rowrange[1]; i2++)
		{
		  Int_t rowindex2 = i2 - AliL3Transform::GetFirstRow(fCurrentPatch);
		  for(Int_t d2=0; d2<fData[rowindex2].fNPads; d2++)
		    {
		      pad2 = d2 + fData[rowindex2].fMinpad;
		      for(Int_t k=0; k<10; k++)
			{
			  time2 = k + fData[rowindex2].fMintime;
			  charge2 = fData[rowindex2].fPads[d2][k];
			  if(charge2==0) continue;
			  AliL3Transform::Slice2Sector(0,i2,sector,row);
			  AliL3Transform::Raw2Local(hit2,sector,row,pad2,time2);
			  r2 = sqrt(hit2[0]*hit2[0]+hit2[1]*hit2[1]);
			  delta2 = atan(hit2[2]/r2);
			  delta = (charge1*delta1 + charge2*delta2)/(charge1+charge2);
			  hit2[0] = hit2[0] - AliL3Transform::Row2X(rowrange[0]);
			  
			  theta = atan2(hit2[0]-hit[0],hit[1]-hit2[1]);
			  rho = hit[0]*cos(theta)+hit[1]*sin(theta);
			  hist->Fill(theta,rho,delta,charge1+charge2);
			}
		    }
		}
	    }
	}
    }  
}

void AliL3HoughTest::TransformLines2Circle(TH3 *hist,AliL3TrackArray *tracks)
{
  //Another HT? 
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *tr = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!tr) continue;
      Int_t middlerow = (Int_t)(tr->GetFirstRow()+(tr->GetLastRow()-tr->GetFirstRow())/2);
      Float_t hit[3];
      tr->GetLineCrossingPoint(middlerow,hit);
      hit[0] += AliL3Transform::Row2X(tr->GetFirstRow());
      Float_t r = sqrt(hit[0]*hit[0] + hit[1]*hit[1]);
      hit[2] = r*tr->GetTgl();
      Float_t phi = atan2(hit[1],hit[0]);
      Float_t theta = tr->GetPsiLine() - AliL3Transform::Pi()/2;
      Float_t psi = 2*phi - theta;
      Float_t kappa = 2/r*sin(phi-psi);
      Float_t delta = atan(hit[2]/r);
      hist->Fill(kappa,psi,delta,tr->GetWeight());
      
    }
}

void AliL3HoughTest::Transform2Center(AliL3Histogram *hist)
{
  //Choose parameter space to be center of curvature.
  
  if(!fData)
    {
      cerr<<"AliL3HoughTest::TransformC : No data"<<endl;
      return;
    }
  Int_t pad1,pad2,time1,time2,charge1,charge2,sector,row;
  Float_t r1,r2,phi1,phi2,hit[3],hit2[3];
  //Float_t phi_0,kappa;
  for(Int_t i=AliL3Transform::GetFirstRow(fCurrentPatch); i<=AliL3Transform::GetLastRow(fCurrentPatch); i++)
    {
      Int_t rowindex1 = i - AliL3Transform::GetFirstRow(fCurrentPatch);
      for(Int_t d1=0; d1<fData[rowindex1].fNPads; d1++)
	{
	  pad1 = d1 + fData[rowindex1].fMinpad;
	  for(Int_t j=0; j<10; j++)
	    {
	      time1 = j + fData[rowindex1].fMintime;
	      charge1 = fData[rowindex1].fPads[d1][j];
	      if(charge1==0) continue;
	      AliL3Transform::Slice2Sector(0,i,sector,row);
	      AliL3Transform::Raw2Local(hit,sector,row,pad1,time1);
	      r1 = sqrt(hit[0]*hit[0]+hit[1]*hit[1])/2;
	      phi1 = atan2(hit[1],hit[0]);
	      
	      for(Int_t j=i+1; j<=AliL3Transform::GetLastRow(fCurrentPatch); j++)
		{
		  Int_t rowindex2 = j - AliL3Transform::GetFirstRow(fCurrentPatch);
		  for(Int_t d2=0; d2<fData[rowindex2].fNPads; d2++)
		    {
		      pad2 = d2 + fData[rowindex2].fMinpad;
		      for(Int_t k=0; k<10; k++)
			{
			  time2 = k + fData[rowindex2].fMintime;
			  charge2 = fData[rowindex2].fPads[d2][k];
			  if(charge2==0) continue;
			  AliL3Transform::Slice2Sector(0,j,sector,row);
			  AliL3Transform::Raw2Local(hit2,sector,row,pad2,time2);
			  r2 = sqrt(hit2[0]*hit2[0]+hit2[1]*hit2[1])/2;
			  phi2 = atan2(hit2[1],hit2[0]);
			  Float_t yc = (r2-(r1/cos(phi1))*cos(phi2))/(sin(phi2)-tan(phi1)*cos(phi2));
			  Float_t xc = r1/cos(phi1) - (r2*tan(phi1)-r1*sin(phi1)*cos(phi2))/(sin(phi2)-tan(phi1)*cos(phi2));
			  hist->Fill(xc,yc,charge1+charge2);
			}
		    }
		}
	    }
	}
    }
}


void AliL3HoughTest::FindAbsMaxima(TH3 *hist,Int_t zsearch,Float_t &maxx,Float_t &maxy,Float_t &maxz,Int_t &maxvalue) const
{
  //Find peaks in the Hough space  
  TH1 *h1 = hist->Project3D("z");
      
  TAxis *z = hist->GetZaxis();
  Float_t zpeak[50];
  Int_t n=0,i,zbin[50];
  for(i=z->GetFirst()+1; i<z->GetLast()-1; i++)
    {
      int bin1 = h1->GetBin(i-1);
      int bin2 = h1->GetBin(i);
      int bin3 = h1->GetBin(i+1);
      if(h1->GetBinContent(bin2) > h1->GetBinContent(bin1) && h1->GetBinContent(bin2) > h1->GetBinContent(bin3))
	{
	  zbin[n]=bin2;
	  zpeak[n++]=h1->GetBinCenter(bin2);
	}
    }
  
  Int_t zrange[2] = {z->GetFirst(),z->GetLast()};
  z->SetRange(zbin[0]-zsearch,zbin[0]+zsearch);

  TH1 *h2 = hist->Project3D("yx");
  z->SetRange(zrange[0],zrange[1]);

  TAxis *x = h2->GetXaxis();
  TAxis *y = h2->GetYaxis();
  maxvalue=0;
  for(i=0; i<n; i++)
    {
      for(Int_t xbin=x->GetFirst(); xbin<=x->GetLast(); xbin++)
	{
	  Float_t xvalue = x->GetBinCenter(xbin);
	  for(Int_t ybin=y->GetFirst(); ybin<=y->GetLast(); ybin++)
	    {
	      Float_t yvalue = y->GetBinCenter(ybin);
	      Int_t value = (Int_t)h2->GetBinContent(xbin,ybin);
	      if(value > maxvalue)
		{
		  maxvalue=value;
		  maxx = xvalue;
		  maxy = yvalue;
		  maxz = zpeak[i];
		}
	    }
	}
    }
  
  delete h1;
  delete h2;
}

