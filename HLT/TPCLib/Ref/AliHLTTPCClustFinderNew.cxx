// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCClustFinderNew.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCMemHandler.h"

#if GCCVERSION == 3
using namespace std;
#endif

/** \class AliHLTTPCClustFinderNew
<pre>
//_____________________________________________________________
// AliHLTTPCClustFinderNew
//
// The current cluster finder for HLT
// (Based on STAR L3)
// 
// The cluster finder is initialized with the Init function, 
// providing the slice and patch information to work on. 
// The input is a AliHLTTPCDigitRowData structure using the 
// Read function. The resulting space points will be in the
// array given by the SetOutputArray function.
// 
// There are several setters which control the behaviour:
//
// - SetXYError(Float_t):   set fixed error in XY direction
// - SetZError(Float_t):    set fixed error in Z  direction
//                            (used if errors are not calculated) 
// - SetDeconv(Bool_t):     switch on/off deconvolution
// - SetThreshold(UInt_t):  set charge threshold for cluster
// - SetMatchWidth(UInt_t): set the match distance in 
//                            time for sequences to be merged 
// - SetSTDOutput(Bool_t):  switch on/off output about found clusters   
// - SetCalcErr(Bool_t):    switch on/off calculation of 
//                          space point errors (or widths in raw system)
// - SetRawSP(Bool_t):      switch on/off convertion to raw system
//
//
// Example Usage:
//
// AliHLTTPCFileHandler *file = new AliHLTTPCFileHandler();
// file->SetAliInput(digitfile); //give some input file
// for(int slice=0; slice<=35; slice++){
//   for(int patch=0; pat<6; pat++){
//     file->Init(slice,patch);
//     UInt_t ndigits=0;
//     UInt_t maxclusters=100000;
//     UInt_t pointsize = maxclusters*sizeof(AliHLTTPCSpacePointData);
//     AliHLTTPCSpacePointData *points = (AliHLTTPCSpacePointData*)memory->Allocate(pointsize);
//     AliHLTTPCDigitRowData *digits = (AliHLTTPCDigitRowData*)file->AliAltroDigits2Memory(ndigits,event);
//     AliHLTTPCClustFinderNew *cf = new AliHLTTPCClustFinderNew();
//     cf->SetMatchWidth(2);
//     cf->InitSlice(slice,patch,maxclusters);
//     cf->SetSTDOutput(kTRUE);    //Some output to standard IO
//     cf->SetRawSP(kFALSE);       //Convert space points to local system
//     cf->SetThreshold(5);        //Threshold of cluster charge
//     cf->SetDeconv(kTRUE);       //Deconv in pad and time direction
//     cf->SetCalcErr(kTRUE);      //Calculate the errors of the spacepoints
//     cf->SetOutputArray(points); //Move the spacepoints to the array
//     cf->Read(ndigits,digits);   //give the data to the cf
//     cf->ProcessDigits();        //process the rows given by init
//     Int_t npoints = cf->GetNumberOfClusters();
//     AliHLTTPCMemHandler *out= new AliHLTTPCMemHandler();
//     out->SetBinaryOutput(fname);
//     out->Memory2Binary(npoints,points); //store the spacepoints
//     out->CloseBinaryOutput();
//     delete out;
//     file->free();
//     delete cf;
//   }
// }
</pre> 
*/

ClassImp(AliHLTTPCClustFinderNew)

AliHLTTPCClustFinderNew::AliHLTTPCClustFinderNew()
{
  fMatch = 1;
  fThreshold = 10;
  fXYErr = 0.2;
  fZErr = 0.3;
  fDeconvPad = kTRUE;
  fDeconvTime = kTRUE;
  fStdout = kFALSE;
  fCalcerr = kTRUE;
  fRawSP = kFALSE;
  fFirstRow=0;
  fLastRow=0;
}

AliHLTTPCClustFinderNew::~AliHLTTPCClustFinderNew()
{
}

void AliHLTTPCClustFinderNew::InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t nmaxpoints)
{
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow = firstrow;
  fLastRow = lastrow;
}

void AliHLTTPCClustFinderNew::InitSlice(Int_t slice,Int_t patch,Int_t nmaxpoints)
{
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow=AliHLTTPCTransform::GetFirstRow(patch);
  fLastRow=AliHLTTPCTransform::GetLastRow(patch);
}

void AliHLTTPCClustFinderNew::SetOutputArray(AliHLTTPCSpacePointData *pt)
{
  fSpacePointData = pt;
}

void AliHLTTPCClustFinderNew::Read(UInt_t ndigits,AliHLTTPCDigitRowData *ptr)
{
  fNDigitRowData = ndigits;
  fDigitRowData = ptr;
}

void AliHLTTPCClustFinderNew::ProcessDigits()
{
  //Loop over rows, and call processrow
  
  AliHLTTPCDigitRowData *tempPt = (AliHLTTPCDigitRowData*)fDigitRowData;
  
  for(Int_t i=fFirstRow; i<=fLastRow; i++)
    {
      fCurrentRow = i;
      if((Int_t)tempPt->fRow!=fCurrentRow){
	LOG(AliHLTTPCLog::kWarning,"AliHLTTPCClustFinderNew::ProcessDigits","Digits")
	  <<"Row number should match! "<<tempPt->fRow<<" "<<fCurrentRow<<ENDLOG;
	continue;
      }
#if 0
      LOG(AliHLTTPCLog::kDebug,"AliHLTTPCClustFinderNew::ProcessDigits","Digits")
	  << "Row " << AliHLTTPCLog::kDec << tempPt->fRow << " digits: "
	  << tempPt->fNDigit << " (Offset: " << ((unsigned long)tempPt) - ((unsigned long)fDigitRowData)
	  << ")." << ENDLOG;
#endif
      ProcessRow(tempPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliHLTTPCDigitRowData) + tempPt->fNDigit*sizeof(AliHLTTPCDigitData);
      tmp += size;
      tempPt = (AliHLTTPCDigitRowData*)tmp;
    }
  LOG(AliHLTTPCLog::kDebug,"AliHLTTPCClustFinderNew::ProcessDigits","Space points")
    <<"Cluster finder found "<<fNClusters<<" clusters in slice "<<fCurrentSlice
    <<" patch "<<fCurrentPatch<<ENDLOG;
}

void AliHLTTPCClustFinderNew::ProcessRow(AliHLTTPCDigitRowData *tempPt)
{

  UInt_t last_pad = 123456789;

  ClusterData *pad1[5000]; //2 lists for internal memory=2pads
  ClusterData *pad2[5000]; //2 lists for internal memory=2pads
  ClusterData clusterlist[10000]; //Clusterlist

  ClusterData **currentPt;  //List of pointers to the current pad
  ClusterData **previousPt; //List of pointers to the previous pad
  currentPt = pad2;
  previousPt = pad1;
  UInt_t n_previous=0,n_current=0,n_total=0;

  //Loop over sequences of this row:
  for(UInt_t bin=0; bin<tempPt->fNDigit; bin++)
    {
      AliHLTTPCDigitData *data = tempPt->fDigitData;
#if 0
      LOG( AliHLTTPCLog::kDebug, "AliHLTTPCClustFinderNew::ProcessRow", "Signal" )
	  << "Pad " << AliHLTTPCLog::kDec << data[bin].fPad 
	  << " (" << bin << ")"
	  << " time " << data[bin].fTime << ": " << data[bin].fCharge
	  << " (Offset: " << ((unsigned long)(data+bin)) - ((unsigned long)fDigitRowData)
	  << ")." << ENDLOG;
#ifdef do_mc
      LOG( AliHLTTPCLog::kWarning, "AliHLTTPCClustFinderNew::ProcessRow", "Signal" )
	  << "Compiled with do_mc" << ENDLOG;
#endif
#endif
      if(data[bin].fPad != last_pad)
	{
	  //This is a new pad
	  
	  //Switch the lists:
	  if(currentPt == pad2)
	    {
	      currentPt = pad1;
	      previousPt = pad2;
	    }
	  else 
	    {
	      currentPt = pad2;
	      previousPt = pad1;
	    }
	  n_previous = n_current;
	  n_current = 0;
	  if(bin[data].fPad != last_pad+1)
	    {
	      //this happens if there is a pad with no signal.
	      n_previous = n_current = 0;
	    }
	  last_pad = data[bin].fPad;
	}

      Bool_t new_cluster = kTRUE;
      UInt_t seq_charge=0,seq_average=0,seq_error=0;
      UInt_t last_charge=0,last_was_falling=0;
      Int_t new_bin=-1;

      if(fDeconvTime)
	{
	redo: //This is a goto.
	  if(new_bin > -1)
	    {
	      bin = new_bin;
	      new_bin = -1;
	    }
	  
	  last_charge=0;
	  last_was_falling = 0;
	}
      
      while(1) //Loop over current sequence
	{
	  if(data[bin].fTime >= AliHLTTPCTransform::GetNTimeBins())
	    {
	      LOG(AliHLTTPCLog::kFatal,"AliHLTTPCClustFinderNew::ProcessRow","Digits")
		<<"Timebin out of range "<<(Int_t)data[bin].fTime<<ENDLOG;
	      break;
	    }

	  //Get the current ADC-value
	  UInt_t charge = data[bin].fCharge;
	  
	  if(fDeconvTime)
	    {
	      //Check if the last pixel in the sequence is smaller than this
	      if(charge > last_charge)
		{
		  if(last_was_falling)
		    {
		      new_bin = bin;
		      break;
		    }
		}
	      else last_was_falling = 1; //last pixel was larger than this
	      last_charge = charge;
	    }
	  
	  //Sum the total charge of this sequence
	  seq_charge += charge;
	  seq_average += data[bin].fTime*charge;
	  seq_error += data[bin].fTime*data[bin].fTime*charge;

	  //Check where to stop:
	  if(bin >= tempPt->fNDigit - 1) //out of range
	    break;
	  if(data[bin+1].fPad != data[bin].fPad) //new pad
	    break;
	  if(data[bin+1].fTime != data[bin].fTime+1) //end of sequence
	    break;
	  
	  bin++;
	}//end loop over sequence
      
      //Calculate mean of sequence:
      Int_t seq_mean=0;
      if(seq_charge)
	seq_mean = seq_average/seq_charge;
      else
	{
	  LOG(AliHLTTPCLog::kFatal,"AliHLTTPCClustFinderNew::ProcessRow","Data")
	    <<"Error in data given to the cluster finder"<<ENDLOG;
	  seq_mean = 1;
	  seq_charge = 1;
	}
      
      //Calculate mean in pad direction:
      Int_t pad_mean = seq_charge*data[bin].fPad;
      Int_t pad_error = data[bin].fPad*pad_mean;

      //Compare with results on previous pad:
      for(UInt_t p=0; p<n_previous; p++)
	{
	  //dont merge sequences on the same pad twice
	  if(previousPt[p]->fLastMergedPad==data[bin].fPad) continue;

	  Int_t difference = seq_mean - previousPt[p]->fMean;
	  if(difference < -fMatch) break;

	  if(difference <= fMatch) //There is a match here!!
	    {
	      ClusterData *local = previousPt[p];

	      if(fDeconvPad)
		{
		  if(seq_charge > local->fLastCharge)
		    {
		      if(local->fChargeFalling) //The previous pad was falling
			{			
			  break; //create a new cluster
			}		    
		    }
		  else
		    local->fChargeFalling = 1;
		  local->fLastCharge = seq_charge;
		}
	      
	      //Don't create a new cluster, because we found a match
	      new_cluster = kFALSE;
	      
	      //Update cluster on current pad with the matching one:
	      local->fTotalCharge += seq_charge;
	      local->fPad += pad_mean;
	      local->fPad2 += pad_error;
	      local->fTime += seq_average;
	      local->fTime2 += seq_error;
	      local->fMean = seq_mean;
	      local->fFlags++; //means we have more than one pad 
	      local->fLastMergedPad = data[bin].fPad;

	      currentPt[n_current] = local;
	      n_current++;
	      
	      break;
	    } //Checking for match at previous pad
	} //Loop over results on previous pad.
      
      if(new_cluster)
	{
	  //Start a new cluster. Add it to the clusterlist, and update
	  //the list of pointers to clusters in current pad.
	  //current pad will be previous pad on next pad.

	  //Add to the clusterlist:
	  ClusterData *tmp = &clusterlist[n_total];
	  tmp->fTotalCharge = seq_charge;
	  tmp->fPad = pad_mean;
	  tmp->fPad2 = pad_error;
	  tmp->fTime = seq_average;
	  tmp->fTime2 = seq_error;
	  tmp->fMean = seq_mean;
	  tmp->fFlags = 0;  //flags for single pad clusters
	  tmp->fLastMergedPad = data[bin].fPad;

	  if(fDeconvPad)
	    {
	      tmp->fChargeFalling = 0;
	      tmp->fLastCharge = seq_charge;
	    }

	  //Update list of pointers to previous pad:
	  currentPt[n_current] = &clusterlist[n_total];
	  n_total++;
	  n_current++;
	}

      if(fDeconvTime)
	if(new_bin >= 0) goto redo;
    }//Loop over digits on this padrow
  
  WriteClusters(n_total,clusterlist);
}

void AliHLTTPCClustFinderNew::WriteClusters(Int_t n_clusters,ClusterData *list)
{
  Int_t thisrow,thissector;
  UInt_t counter = fNClusters;
  
  for(int j=0; j<n_clusters; j++)
    {
      if(!list[j].fFlags) continue; //discard single pad clusters
      if(list[j].fTotalCharge < fThreshold) continue; //noise cluster

      Float_t xyz[3];      
      Float_t fpad =(Float_t)list[j].fPad /(Float_t)list[j].fTotalCharge;
      Float_t fpad2=fXYErr*fXYErr; //fixed given error
      Float_t ftime =(Float_t)list[j].fTime /(Float_t)list[j].fTotalCharge;
      Float_t ftime2=fZErr*fZErr;  //fixed given error

      if(fCalcerr) { //calc the errors, otherwice take the fixed error 
	Int_t patch = AliHLTTPCTransform::GetPatch(fCurrentRow);

	Float_t sy2=(Float_t)list[j].fPad2/(Float_t)list[j].fTotalCharge - fpad*fpad;
	if(sy2 < 0) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinderNew::WriteClusters","Cluster width")
	      <<"SigmaY2 negative "<<sy2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	    continue;
	} else {
	  if(!fRawSP){
	    fpad2 = (sy2 + 1./12)*AliHLTTPCTransform::GetPadPitchWidth(patch)*AliHLTTPCTransform::GetPadPitchWidth(patch);
	    if(sy2 != 0){
	      fpad2*=0.108; //constants are from offline studies
	      if(patch<2)
		fpad2*=2.07;
	    }
	  } else fpad2=sy2; //take the width not the error
	}
	Float_t sz2=(Float_t)list[j].fTime2/(Float_t)list[j].fTotalCharge - ftime*ftime;
	if(sz2 < 0){
	  LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinderNew::WriteClusters","Cluster width")
	    <<"SigmaZ2 negative "<<sz2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	  continue;
	} else {
	  if(!fRawSP){
	    ftime2 = (sz2 + 1./12)*AliHLTTPCTransform::GetZWidth()*AliHLTTPCTransform::GetZWidth();
	    if(sz2 != 0) {
	      ftime2 *= 0.169; //constants are from offline studies
	      if(patch<2)
		ftime2 *= 1.77;
	    }
	  } else ftime2=sz2; //take the width, not the error
	}
      }
      if(fStdout==kTRUE)
	  {
#if 0
	cout<<"WriteCluster: padrow "<<fCurrentRow<<" pad "<<fpad << " +- "<<fpad2<<" time "<<ftime<<" +- "<<ftime2<<" charge "<<list[j].fTotalCharge<<endl;
#endif
	  }
      
      if(!fRawSP){
	AliHLTTPCTransform::Slice2Sector(fCurrentSlice,fCurrentRow,thissector,thisrow);
	AliHLTTPCTransform::Raw2Local(xyz,thissector,thisrow,fpad,ftime);
	
	if(xyz[0]==0) LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinder","Cluster Finder")
	  <<AliHLTTPCLog::kDec<<"Zero cluster"<<ENDLOG;
	if(fNClusters >= fMaxNClusters)
	  {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinder::WriteClusters","Cluster Finder")
	      <<AliHLTTPCLog::kDec<<"Too many clusters "<<fNClusters<<ENDLOG;
	    return;
	  }  
	
	fSpacePointData[counter].fX = xyz[0];
	fSpacePointData[counter].fY = xyz[1];
	fSpacePointData[counter].fZ = xyz[2];
	
      } else {
	fSpacePointData[counter].fX = fCurrentRow;
	fSpacePointData[counter].fY = fpad;
	fSpacePointData[counter].fZ = ftime;
      }
      
      fSpacePointData[counter].fCharge = list[j].fTotalCharge;
      fSpacePointData[counter].fPadRow = fCurrentRow;
      fSpacePointData[counter].fSigmaY2 = fpad2;
      fSpacePointData[counter].fSigmaZ2  = ftime2;

      Int_t patch=fCurrentPatch;
      if(patch==-1) patch=0; //never store negative patch number
      fSpacePointData[counter].fID = counter
	+((fCurrentSlice&0x7f)<<25)+((patch&0x7)<<22);//Uli
#ifdef do_mc
      Int_t trackID[3];
      GetTrackID((Int_t)rint(fpad),(Int_t)rint(ftime),trackID);

      fSpacePointData[counter].fTrackID[0] = trackID[0];
      fSpacePointData[counter].fTrackID[1] = trackID[1];
      fSpacePointData[counter].fTrackID[2] = trackID[2];

      //cout<<"padrow "<<fCurrentRow<<" pad "<<(Int_t)rint(fpad)<<" time "<<(Int_t)rint(ftime)<<" Trackid "<<trackID[0]<<endl;
#endif
      
      fNClusters++;
      counter++;
    }
}

#ifdef do_mc
void AliHLTTPCClustFinderNew::GetTrackID(Int_t pad,Int_t time,Int_t *trackID)
{
  AliHLTTPCDigitRowData *rowPt = (AliHLTTPCDigitRowData*)fDigitRowData;
  
  trackID[0]=trackID[1]=trackID[2]=-2;
  //cout<<"Looking for pad "<<pad<<" time "<<time<<endl;
  for(Int_t i=fFirstRow; i<=fLastRow; i++){
    if(rowPt->fRow < (UInt_t)fCurrentRow){
      AliHLTTPCMemHandler::UpdateRowPointer(rowPt);
      continue;
    }
    AliHLTTPCDigitData *digPt = (AliHLTTPCDigitData*)rowPt->fDigitData;
    for(UInt_t j=0; j<rowPt->fNDigit; j++){
      Int_t cpad = digPt[j].fPad;
      Int_t ctime = digPt[j].fTime;
      if(cpad != pad) continue;
      if(ctime != time) continue;

      trackID[0] = digPt[j].fTrackID[0];
      trackID[1] = digPt[j].fTrackID[1];
      trackID[2] = digPt[j].fTrackID[2];
      
      //cout<<"Reading row "<<fCurrentRow<<" pad "<<cpad<<" time "<<ctime<<" trackID "<<digPt[j].fTrackID[0]<<endl;
      break;
    }
    break;
  }
}
#endif
