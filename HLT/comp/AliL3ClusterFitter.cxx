//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV

#include "AliL3StandardIncludes.h"
#include "AliL3ClusterFitter.h"
#include "AliL3FitUtilities.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3ModelTrack.h"
#include "AliL3TrackArray.h"
#include "AliL3MemHandler.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliL3ClusterFitter
//


ClassImp(AliL3ClusterFitter)

AliL3ClusterFitter::AliL3ClusterFitter()
{
  fPadFitRange=2;
  fTimeFitRange=3;
  plane=0;
  fNmaxOverlaps = 3;
  fChiSqMax=12;
}

AliL3ClusterFitter::~AliL3ClusterFitter()
{

}

void AliL3ClusterFitter::FindClusters()
{
  if(!fTracks)
    {
      cerr<<"AliL3ClusterFitter::Process : No tracks"<<endl;
      return;
    }
  if(!fRowData)
    {
      cerr<<"AliL3ClusterFitter::Process : No data "<<endl;
      return;
    }
  
  AliL3DigitRowData *rowPt = fRowData;
  AliL3DigitData *digPt=0;

  Int_t pad,time;
  Short_t charge;

  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      fCurrentPadRow = i;
      memset((void*)fRow,0,(AliL3Transform::GetNTimeBins()+1)*(AliL3Transform::GetNPads(i)+1)*sizeof(Digit));
      digPt = (AliL3DigitData*)rowPt->fDigitData;
      //cout<<"Loading row "<<i<<" with "<<(Int_t)rowPt->fNDigit<<" digits"<<endl;
      for(UInt_t j=0; j<rowPt->fNDigit; j++)
	{
	  pad = digPt[j].fPad;
	  time = digPt[j].fTime;
	  charge = digPt[j].fCharge;
	  fRow[(AliL3Transform::GetNTimeBins()+1)*pad+time].fCharge = charge;
	  fRow[(AliL3Transform::GetNTimeBins()+1)*pad+time].fUsed = kFALSE;
	  //cout<<"Row "<<i<<" pad "<<pad<<" time "<<time<<" charge "<<charge<<endl;
	}
      
      for(Int_t k=0; k<fTracks->GetNTracks(); k++)
	{
	  AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  if(track->GetPadHit(i) < 0 || track->GetTimeHit(i) < 0)
	    {
	      track->SetCluster(i,0,0,0,0,0,0);
	      continue;
	    }
	  FitClusters(track);
	}
      FillZeros(rowPt);
      AliL3MemHandler::UpdateRowPointer(rowPt);
    }
}

void AliL3ClusterFitter::FitCluster(AliL3ModelTrack *track)
{
  cout<<"Track had no overlaps"<<endl;
  Int_t minpad = (Int_t)rint(track->GetPadHit(fCurrentPadRow)) - fPadFitRange;
  Int_t maxpad = (Int_t)rint(track->GetPadHit(fCurrentPadRow)) + fPadFitRange;
  Int_t mintime = (Int_t)rint(track->GetTimeHit(fCurrentPadRow)) - fTimeFitRange;
  Int_t maxtime = (Int_t)rint(track->GetTimeHit(fCurrentPadRow)) + fTimeFitRange;
  
  Int_t size = FIT_PTS;

  if(minpad <= 0)
    minpad=0;
  if(mintime<=0)
    mintime=0;
  if(maxpad>=AliL3Transform::GetNPads(fCurrentPadRow)-1)
    maxpad=AliL3Transform::GetNPads(fCurrentPadRow)-1;
  if(maxtime>=AliL3Transform::GetNTimeBins()-1)
    maxtime=AliL3Transform::GetNTimeBins()-1;

  if(plane)
    delete [] plane;
  plane = new DPOINT[FIT_PTS];
  memset(plane,0,FIT_PTS*sizeof(DPOINT));

  Double_t x[FIT_PTS],y[FIT_PTS],s[FIT_PTS];
  
  //Fill the fit parameters:
  Double_t a[FIT_MAXPAR];
  a[2] = track->GetPadHit(fCurrentPadRow);
  a[4] = track->GetTimeHit(fCurrentPadRow);
  a[1] = fRow[(AliL3Transform::GetNTimeBins()+1)*((Int_t)rint(a[2])) + (Int_t)rint(a[4])].fCharge;
  a[3] = sqrt(track->GetParSigmaY2(fCurrentPadRow))*sqrt(2);
  a[5] = sqrt(track->GetParSigmaZ2(fCurrentPadRow))*sqrt(2);
  a[6] = sqrt(track->GetParSigmaZ2(fCurrentPadRow))*sqrt(2);
  
  if(!a[1])
    {
      cout<<"No charge"<<endl;
      if(track->GetNClusters() == fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
	track->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
      return;
    }
  
  Int_t pad_num=0;
  Int_t time_num_max=0;
  Int_t ndata=0;
  Int_t charge,tot_charge=0;

  //Fill the proper arrays:
  for(Int_t i=minpad; i<=maxpad; i++)
    {
      Int_t max_charge = 0;
      Int_t time_num=0;
      for(Int_t j=mintime; j<=maxtime; j++)
	{
	  charge = fRow[(AliL3Transform::GetNTimeBins()+1)*i + j].fCharge;
	  if(charge > max_charge)
	    {
	      max_charge = charge;
	      time_num++;
	     }
	  if(charge <= 0) continue;
	  //cout<<"Filling padrow "<<fCurrentPadRow<<" pad "<<i<<" time "<<j<<" charge "<<charge<<endl;
	  tot_charge += charge;
	  ndata++;
	  if(ndata >= size)
	    cerr<<"Too many points"<<endl;
	  plane[ndata].u = (Double_t)i;
	  plane[ndata].v = (Double_t)j;
	  x[ndata]=ndata;
	  y[ndata]=charge;
	  s[ndata]= 1 + sqrt((Double_t)charge);
	}
      if(max_charge) //there was charge on this pad
	pad_num++;
      if(time_num_max < time_num)
	time_num_max = time_num;
    }

  if(pad_num <= 1 || time_num_max <=1) //too few to do fit
    {
      //cout<<"To few pads for fit on row "<<fCurrentPadRow<<endl;
      if(track->GetNClusters() == fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
	track->SetCluster(fCurrentPadRow,0,0,0,0,0,0);

      return;
    }
  
  Int_t npars = NUM_PARS;
  Int_t lista[FIT_MAXPAR];
  Double_t dev[FIT_MAXPAR],chisq_f;
  lista[1] = 1;
  lista[2] = 1;
  lista[3] = 0;
  lista[4] = 1;
  lista[5] = 0;
  lista[6] = 0;
  
  //Declare a pointer to the C fitting function
  void (*funcs) ( double, double *, double *,double *,int );
  funcs = f2gauss5;
  
  //cout<<"Doing fit with parameters "<<a[2]<<" "<<a[4]<<" "<<a[3]<<" "<<a[5]<<" "<<a[1]<<endl;
  
  Int_t ret = lev_marq_fit( x, y, s, ndata, a, lista, dev, npars, &chisq_f,f2gauss5);
  if(ret)
    cerr<<"Fit error"<<endl;
  
  tot_charge = (Int_t)(a[1] * a[3] * a[5]);
  chisq_f /= (ndata - 3);
  //  cout<<"Chisq per degree of freedom : "<<chisq_f<<endl;
  //cout<<"pad "<<a[2]<<" time "<<a[4]<<" adc "<<a[1]<<endl;
  //cout<<"Setting cluster on row "<<fCurrentPadRow<<" pad "<<a[2]<<" time "<<a[4]<<" charge "<<tot_charge<<endl;
  
  //Make sure that the cluster has not already been set before setting it:
  if(track->GetNClusters() == fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
    track->SetCluster(fCurrentPadRow,a[2],a[4],tot_charge,0,0,pad_num);
}

void AliL3ClusterFitter::FitClusters(AliL3ModelTrack *track)
{
  //Handle single and overlapping clusters
    
  if(!track->IsPresent(fCurrentPadRow) && track->GetNClusters() != fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
    {    
      if(!track->IsPresent(fCurrentPadRow) && 
	 track->GetNClusters() != fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch) + 1)//debug
	{
	  cerr<<"AliL3ClusterFitter::FitClusters() : Mismatching clustercount "
	      <<track->GetNClusters()<<" "<<fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch) + 1<<endl;
	  exit(5);
	}
      return; //This cluster has been set before.
    }
  
  
  Int_t minpad,maxpad,mintime,maxtime;
  
  Int_t size = FIT_PTS;
  Int_t max_tracks = FIT_MAXPAR/NUM_PARS;
  if(track->GetNOverlaps(fCurrentPadRow) > max_tracks)
    {
      cerr<<"AliL3ClusterFitter::FitOverlappingClusters : Too many overlapping tracks"<<endl;
      return;
    }
  Int_t *overlaps = track->GetOverlaps(fCurrentPadRow);
  
  //Check if at least one cluster is not already fitted
  Bool_t all_fitted=kTRUE;
  Int_t k=-1;
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliL3ModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliL3ModelTrack*)fTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      if(!tr->IsPresent(fCurrentPadRow))
	{
	  all_fitted = kFALSE;
	  break;
	}
    }
  if(all_fitted)
    {
      cout<<"But all the clusters were already fitted on row "<<fCurrentPadRow<<endl;
      return;
    }
  
  plane = new DPOINT[FIT_PTS];
  memset(plane,0,FIT_PTS*sizeof(DPOINT));

  Double_t x[FIT_PTS],y[FIT_PTS],s[FIT_PTS];
  
  //Fill the fit parameters:
  Double_t a[FIT_MAXPAR];
  Int_t lista[FIT_MAXPAR];
  Double_t dev[FIT_MAXPAR],chisq_f;
  
  minpad = mintime = 999;
  maxpad = maxtime = 0;
  Int_t fit_pars=0;
  
  Int_t ntracks = track->GetNOverlaps(fCurrentPadRow)+1;
  Bool_t *do_fit = new Bool_t[ntracks];
  for(Int_t i=0; i<ntracks; i++)
    do_fit[i]=kTRUE;
  
  Int_t n_overlaps=0;
  k=-1;
  //Fill the overlapping tracks:
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliL3ModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliL3ModelTrack*)fTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      
      Int_t hitpad = (Int_t)rint(tr->GetPadHit(fCurrentPadRow));
      Int_t hittime = (Int_t)rint(tr->GetTimeHit(fCurrentPadRow));
      Int_t charge = fRow[(AliL3Transform::GetNTimeBins()+1)*hitpad + hittime].fCharge;
      if(!charge) //There is not charge here, so the cluster is non-existing -- remove it.
	{
	  if(tr->GetNClusters() == fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
	    tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);		  
	  do_fit[k] = kFALSE;
	  continue;
	}
      
      if(k==0)
	LocateCluster(tr,minpad,maxpad,mintime,maxtime);

      if(tr->GetPadHit(fCurrentPadRow) < (Double_t)minpad || tr->GetPadHit(fCurrentPadRow) > (Double_t)maxpad ||
	 tr->GetTimeHit(fCurrentPadRow) < (Double_t)mintime || tr->GetTimeHit(fCurrentPadRow) > (Double_t)maxtime)
	{
	  do_fit[k] = kFALSE;//This cluster is outside the region already specified, so it will not be included in this fit.

	  //If this is the first: remove it, because it will not be checked again.
	  if(k==0)
	    if(tr->GetNClusters() == fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
	      tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);		  
	  continue;
	}
      
      
      cout<<"Fitting track cluster, pad "<<tr->GetPadHit(fCurrentPadRow)<<" time "
	  <<tr->GetTimeHit(fCurrentPadRow)<<" charge "<<charge<<" xywidth "<<sqrt(tr->GetParSigmaY2(fCurrentPadRow))
	  <<" zwidth "<<sqrt(tr->GetParSigmaZ2(fCurrentPadRow))<<endl;
      
      a[n_overlaps*NUM_PARS+2] = tr->GetPadHit(fCurrentPadRow);
      a[n_overlaps*NUM_PARS+4] = tr->GetTimeHit(fCurrentPadRow);
      
      if(!tr->IsPresent(fCurrentPadRow)) //Cluster is not fitted before
	{
	  a[n_overlaps*NUM_PARS+1] = charge;
	  a[n_overlaps*NUM_PARS+3] = sqrt(tr->GetParSigmaY2(fCurrentPadRow))*sqrt(2);
	  a[n_overlaps*NUM_PARS+5] = sqrt(tr->GetParSigmaZ2(fCurrentPadRow))*sqrt(2);
	  a[n_overlaps*NUM_PARS+6] = sqrt(tr->GetParSigmaZ2(fCurrentPadRow))*sqrt(2);
	  lista[n_overlaps*NUM_PARS + 1] = 1;
	  lista[n_overlaps*NUM_PARS + 2] = 1;
	  lista[n_overlaps*NUM_PARS + 3] = 0;
	  lista[n_overlaps*NUM_PARS + 4] = 1;
	  lista[n_overlaps*NUM_PARS + 5] = 0;
	  lista[n_overlaps*NUM_PARS + 6] = 0;
	  fit_pars             += 3;
	}
      else
	{
	  Int_t charge;
	  Float_t xywidth,zwidth,pad,time;
	  tr->GetPad(fCurrentPadRow,pad);
	  tr->GetTime(fCurrentPadRow,time);
	  tr->GetClusterCharge(fCurrentPadRow,charge);
	  xywidth = sqrt(tr->GetParSigmaY2(fCurrentPadRow));
	  zwidth = sqrt(tr->GetParSigmaZ2(fCurrentPadRow));
	  cout<<"Cluster had been fitted before, pad "<<pad<<" time "<<time<<" charge "<<charge<<" width "<<xywidth<<" "<<zwidth<<endl;
	  
	  a[n_overlaps*NUM_PARS+2] = pad;
	  a[n_overlaps*NUM_PARS+4] = time;
	  a[n_overlaps*NUM_PARS+1] = charge;
	  a[n_overlaps*NUM_PARS+3] = sqrt(xywidth)*sqrt(2);
	  a[n_overlaps*NUM_PARS+5] = sqrt(zwidth)*sqrt(2);
	  a[n_overlaps*NUM_PARS+6] = sqrt(zwidth)*sqrt(2);

	  lista[n_overlaps*NUM_PARS + 1] = 1;
	  lista[n_overlaps*NUM_PARS + 2] = 0;
	  lista[n_overlaps*NUM_PARS + 3] = 0;
	  lista[n_overlaps*NUM_PARS + 4] = 0;
	  lista[n_overlaps*NUM_PARS + 5] = 0;
	  lista[n_overlaps*NUM_PARS + 6] = 0;
	  fit_pars             += 1;
	}
      n_overlaps++;
    }
  
  if(n_overlaps==0) //No clusters here
    return;
  cout<<"Setting init searchrange; pad "<<minpad<<" "<<maxpad<<" time "<<mintime<<" "<<maxtime<<endl;

  if(minpad <= 0)
    minpad=0;
  if(mintime<=0)
    mintime=0;
  if(maxpad>=AliL3Transform::GetNPads(fCurrentPadRow)-1)
    maxpad=AliL3Transform::GetNPads(fCurrentPadRow)-1;
  if(maxtime>=AliL3Transform::GetNTimeBins()-1)
    maxtime=AliL3Transform::GetNTimeBins()-1;
  
  Int_t pad_num=0;
  Int_t time_num_max=0;
  Int_t ndata=0;
  Int_t tot_charge=0;

  for(Int_t i=minpad; i<=maxpad; i++)
    {
      Int_t max_charge = 0;
      Int_t time_num=0;
      for(Int_t j=mintime; j<=maxtime; j++)
	{
	  Int_t charge = fRow[(AliL3Transform::GetNTimeBins()+1)*i + j].fCharge;
	  
	  if(charge <= 0) continue;

	  if(charge > max_charge)
	    {
	      max_charge = charge;
	      time_num++;
	    }
	  cout<<"Filling padrow "<<fCurrentPadRow<<" pad "<<i<<" time "<<j<<" charge "<<charge<<endl;
	  tot_charge += charge;
	  ndata++;
	  if(ndata >= size)
	    cerr<<"Too many points"<<endl;
	  
	  //This digit will most likely be used:
	  fRow[(AliL3Transform::GetNTimeBins()+1)*i + j].fUsed = kTRUE;
	  plane[ndata].u = (Double_t)i;
	  plane[ndata].v = (Double_t)j;
	  x[ndata]=ndata;
	  y[ndata]=charge;
	  s[ndata]= 1 + sqrt((Double_t)charge);
	}
      if(max_charge) //there was charge on this pad
	pad_num++;
      if(time_num_max < time_num)
	time_num_max = time_num;
    }
  
  k=-1;
  if(pad_num <= 1 || time_num_max <=1 || n_overlaps > fNmaxOverlaps) //too few to do fit
    {
      while(k < track->GetNOverlaps(fCurrentPadRow))
	{
	  AliL3ModelTrack *tr=0;
	  if(k==-1)
	    tr = track;
	  else
	    tr = (AliL3ModelTrack*)fTracks->GetCheckedTrack(overlaps[k]);
	  k++;
	  if(!tr) continue;
	  if(do_fit[k] == kFALSE) continue;
	  
	  if(tr->GetNClusters() == fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch))
	    tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
	}
      cout<<"Too few digits or too many overlaps: "<<pad_num<<" "<<time_num_max<<" "<<n_overlaps<<endl;
      return;
    }
  
  Int_t npars = n_overlaps * NUM_PARS;
  cout<<"Number of overlapping clusters "<<n_overlaps<<endl;
  Int_t ret = lev_marq_fit( x, y, s, ndata, a, lista, dev, npars, &chisq_f, f2gauss5 );
  if(ret)
    {
      cerr<<"Fit error"<<endl;
      exit(5);
    }

  chisq_f /= (ndata-fit_pars);
  cout<<"Chisq "<<chisq_f<<endl;

  k=-1;
  n_overlaps=0;
  while(k < track->GetNOverlaps(fCurrentPadRow))
    {
      AliL3ModelTrack *tr=0;
      if(k==-1)
	tr = track;
      else
	tr = (AliL3ModelTrack*)fTracks->GetCheckedTrack(overlaps[k]);
      k++;
      if(!tr) continue;
      if(do_fit[k] == kFALSE) continue;
      
      if(!tr->IsPresent(fCurrentPadRow))
	{
	  if(tr->GetNClusters() != fCurrentPadRow - AliL3Transform::GetFirstRow(fPatch)) continue;//This cluster has been set before
	  if(chisq_f < fChiSqMax)//cluster fit is good enough
	    {
	      tot_charge = (Int_t)(a[n_overlaps*NUM_PARS+1] * a[n_overlaps*NUM_PARS+3] * a[n_overlaps*NUM_PARS+5]);
	      tr->SetCluster(fCurrentPadRow,a[n_overlaps*NUM_PARS+2],a[n_overlaps*NUM_PARS+4],tot_charge,0,0,pad_num);
	      cout<<"Setting cluster in pad "<<a[n_overlaps*NUM_PARS+2]<<" time "<<a[n_overlaps*NUM_PARS+4]<<" charge "<<tot_charge<<endl;
	    }
	  else //fit was too bad
	    tr->SetCluster(fCurrentPadRow,0,0,0,0,0,0);
	}
      n_overlaps++;
    }
  
  delete [] plane;
  delete [] do_fit;
}

