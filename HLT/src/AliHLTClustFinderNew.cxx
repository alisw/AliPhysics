// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTClustFinderNew.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTMemHandler.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/** \class AliHLTClustFinderNew
<pre>
//_____________________________________________________________
// AliHLTClustFinderNew
//
// The current cluster finder for HLT
// (Based on STAR L3)
// 
// The cluster finder is initialized with the Init function, 
// providing the slice and patch information to work on. 
// The input is a AliHLTDigitRowData structure using the 
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
// AliHLTFileHandler *file = new AliHLTFileHandler();
// file->SetAliInput(digitfile); //give some input file
// for(int slice=0; slice<=35; slice++){
//   for(int patch=0; pat<6; pat++){
//     file->Init(slice,patch);
//     UInt_t ndigits=0;
//     UInt_t maxclusters=100000;
//     UInt_t pointsize = maxclusters*sizeof(AliHLTSpacePointData);
//     AliHLTSpacePointData *points = (AliHLTSpacePointData*)memory->Allocate(pointsize);
//     AliHLTDigitRowData *digits = (AliHLTDigitRowData*)file->AliAltroDigits2Memory(ndigits,event);
//     AliHLTClustFinderNew *cf = new AliHLTClustFinderNew();
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
//     AliHLTMemHandler *out= new AliHLTMemHandler();
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

ClassImp(AliHLTClustFinderNew)

AliHLTClustFinderNew::AliHLTClustFinderNew()
{
  //constructor
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

AliHLTClustFinderNew::~AliHLTClustFinderNew()
{
  //destructor
  ;
}

void AliHLTClustFinderNew::InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t nmaxpoints)
{
  //init slice
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow = firstrow;
  fLastRow = lastrow;
}

void AliHLTClustFinderNew::InitSlice(Int_t slice,Int_t patch,Int_t nmaxpoints)
{
  //init slice
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow=AliHLTTransform::GetFirstRow(patch);
  fLastRow=AliHLTTransform::GetLastRow(patch);
}

void AliHLTClustFinderNew::SetOutputArray(AliHLTSpacePointData *pt)
{
  //set pointer to output
  fSpacePointData = pt;
}

void AliHLTClustFinderNew::Read(UInt_t ndigits,AliHLTDigitRowData *ptr)
{
  //set input pointer
  fNDigitRowData = ndigits;
  fDigitRowData = ptr;
}

void AliHLTClustFinderNew::ProcessDigits()
{
  //Loop over rows, and call processrow
  AliHLTDigitRowData *tempPt = (AliHLTDigitRowData*)fDigitRowData;
  fNClusters = 0; 
  for(Int_t i=fFirstRow; i<=fLastRow; i++)
    {
      fCurrentRow = i;
      if((Int_t)tempPt->fRow!=fCurrentRow){
	LOG(AliHLTLog::kWarning,"AliHLTClustFinderNew::ProcessDigits","Digits")
	  <<"Row number should match! "<<tempPt->fRow<<" "<<fCurrentRow<<ENDLOG;
	continue;
      }
      ProcessRow(tempPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliHLTDigitRowData) + tempPt->fNDigit*sizeof(AliHLTDigitData);
      tmp += size;
      tempPt = (AliHLTDigitRowData*)tmp;
    }
  LOG(AliHLTLog::kInformational,"AliHLTClustFinderNew::ProcessDigits","Space points")
    <<"Cluster finder found "<<fNClusters<<" clusters in slice "<<fCurrentSlice
    <<" patch "<<fCurrentPatch<<ENDLOG;
}

void AliHLTClustFinderNew::ProcessRow(AliHLTDigitRowData *tempPt)
{
  //process row
  UInt_t lastpad = 123456789;

  AliClusterData *pad1[5000]; //2 lists for internal memory=2pads
  AliClusterData *pad2[5000]; //2 lists for internal memory=2pads
  AliClusterData clusterlist[10000]; //Clusterlist

  AliClusterData **currentPt;  //List of pointers to the current pad
  AliClusterData **previousPt; //List of pointers to the previous pad
  currentPt = pad2;
  previousPt = pad1;
  UInt_t nprevious=0,ncurrent=0,ntotal=0;

  //Loop over sequences of this row:
  for(UInt_t bin=0; bin<tempPt->fNDigit; bin++)
    {
      AliHLTDigitData *data = tempPt->fDigitData;
      if(data[bin].fPad != lastpad)
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
	  nprevious = ncurrent;
	  ncurrent = 0;
	  if(bin[data].fPad != lastpad+1)
	    {
	      //this happens if there is a pad with no signal.
	      nprevious = ncurrent = 0;
	    }
	  lastpad = data[bin].fPad;
	}

      Bool_t newcluster = kTRUE;
      UInt_t seqcharge=0,seqaverage=0,seqerror=0;
      UInt_t lastcharge=0,lastwas_falling=0;
      Int_t newbin=-1;

      if(fDeconvTime)
	{
	redo: //This is a goto.
	  if(newbin > -1)
	    {
	      bin = newbin;
	      newbin = -1;
	    }
	  
	  lastcharge=0;
	  lastwas_falling = 0;
	}
      
      while(1) //Loop over current sequence
	{
	  if(data[bin].fTime >= AliHLTTransform::GetNTimeBins())
	    {
	      LOG(AliHLTLog::kFatal,"AliHLTClustFinderNew::ProcessRow","Digits")
		<<"Timebin out of range "<<(Int_t)data[bin].fTime<<ENDLOG;
	      break;
	    }

	  //Get the current ADC-value
	  UInt_t charge = data[bin].fCharge;
	  
	  if(fDeconvTime)
	    {
	      //Check if the last pixel in the sequence is smaller than this
	      if(charge > lastcharge)
		{
		  if(lastwas_falling)
		    {
		      newbin = bin;
		      break;
		    }
		}
	      else lastwas_falling = 1; //last pixel was larger than this
	      lastcharge = charge;
	    }
	  
	  //Sum the total charge of this sequence
	  seqcharge += charge;
	  seqaverage += data[bin].fTime*charge;
	  seqerror += data[bin].fTime*data[bin].fTime*charge;

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
      Int_t seqmean=0;
      if(seqcharge)
	seqmean = seqaverage/seqcharge;
      else
	{
	  LOG(AliHLTLog::kFatal,"AliHLTClustFinderNew::ProcessRow","Data")
	    <<"Error in data given to the cluster finder"<<ENDLOG;
	  seqmean = 1;
	  seqcharge = 1;
	}
      
      //Calculate mean in pad direction:
      Int_t padmean = seqcharge*data[bin].fPad;
      Int_t paderror = data[bin].fPad*padmean;

      //Compare with results on previous pad:
      for(UInt_t p=0; p<nprevious; p++)
	{
	  //dont merge sequences on the same pad twice
	  if(previousPt[p]->fLastMergedPad==data[bin].fPad) continue;

	  Int_t difference = seqmean - previousPt[p]->fMean;
	  if(difference < -fMatch) break;

	  if(difference <= fMatch) //There is a match here!!
	    {
	      AliClusterData *local = previousPt[p];

	      if(fDeconvPad)
		{
		  if(seqcharge > local->fLastCharge)
		    {
		      if(local->fChargeFalling) //The previous pad was falling
			{			
			  break; //create a new cluster
			}		    
		    }
		  else
		    local->fChargeFalling = 1;
		  local->fLastCharge = seqcharge;
		}
	      
	      //Don't create a new cluster, because we found a match
	      newcluster = kFALSE;
	      
	      //Update cluster on current pad with the matching one:
	      local->fTotalCharge += seqcharge;
	      local->fPad += padmean;
	      local->fPad2 += paderror;
	      local->fTime += seqaverage;
	      local->fTime2 += seqerror;
	      local->fMean = seqmean;
	      local->fFlags++; //means we have more than one pad 
	      local->fLastMergedPad = data[bin].fPad;

	      currentPt[ncurrent] = local;
	      ncurrent++;
	      
	      break;
	    } //Checking for match at previous pad
	} //Loop over results on previous pad.
      
      if(newcluster)
	{
	  //Start a new cluster. Add it to the clusterlist, and update
	  //the list of pointers to clusters in current pad.
	  //current pad will be previous pad on next pad.

	  //Add to the clusterlist:
	  AliClusterData *tmp = &clusterlist[ntotal];
	  tmp->fTotalCharge = seqcharge;
	  tmp->fPad = padmean;
	  tmp->fPad2 = paderror;
	  tmp->fTime = seqaverage;
	  tmp->fTime2 = seqerror;
	  tmp->fMean = seqmean;
	  tmp->fFlags = 0;  //flags for single pad clusters
	  tmp->fLastMergedPad = data[bin].fPad;

	  if(fDeconvPad)
	    {
	      tmp->fChargeFalling = 0;
	      tmp->fLastCharge = seqcharge;
	    }

	  //Update list of pointers to previous pad:
	  currentPt[ncurrent] = &clusterlist[ntotal];
	  ntotal++;
	  ncurrent++;
	}

      if(fDeconvTime)
	if(newbin >= 0) goto redo;
    }//Loop over digits on this padrow
  
  WriteClusters(ntotal,clusterlist);
}

void AliHLTClustFinderNew::WriteClusters(Int_t nclusters,AliClusterData *list)
{
  //write cluster to output pointer
  Int_t thisrow,thissector;
  UInt_t counter = fNClusters;
  
  for(int j=0; j<nclusters; j++)
    {
      if(!list[j].fFlags) continue; //discard single pad clusters
      if(list[j].fTotalCharge < fThreshold) continue; //noise cluster

      Float_t xyz[3];      
      Float_t fpad =(Float_t)list[j].fPad / list[j].fTotalCharge;
      Float_t fpad2=fXYErr*fXYErr; //fixed given error
      Float_t ftime =(Float_t)list[j].fTime / list[j].fTotalCharge;
      Float_t ftime2=fZErr*fZErr;  //fixed given error

      if(fCalcerr) { //calc the errors, otherwice take the fixed error 
	Int_t patch = AliHLTTransform::GetPatch(fCurrentRow);
	UInt_t q2=list[j].fTotalCharge*list[j].fTotalCharge;
	if (!q2) {
	    LOG(AliHLTLog::kError,"AliHLTClustFinderNew::WriteClusters","Total charge")
	      <<"Zero total charge "<<q2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	    continue;

	}
	Float_t sy2=list[j].fPad2 * list[j].fTotalCharge - list[j].fPad * list[j].fPad;
	sy2/=q2;
	if(sy2 < 0) {
	    LOG(AliHLTLog::kError,"AliHLTClustFinderNew::WriteClusters","Cluster width")
	      <<"SigmaY2 negative "<<sy2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	    continue;
	} else {
	  if(!fRawSP){
	    fpad2 = (sy2 + 1./12)*AliHLTTransform::GetPadPitchWidth(patch)*AliHLTTransform::GetPadPitchWidth(patch);
	    if(sy2 != 0){
	      fpad2*=0.108; //constants are from offline studies
	      if(patch<2)
		fpad2*=2.07;
	    }
	  } else fpad2=sy2; //take the width not the error
	}
	Float_t sz2=list[j].fTime2*list[j].fTotalCharge - list[j].fTime*list[j].fTime;
	sz2/=q2;
	if(sz2 < 0){
	  LOG(AliHLTLog::kError,"AliHLTClustFinderNew::WriteClusters","Cluster width")
	    <<"SigmaZ2 negative "<<sz2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	  continue;
	} else {
	  if(!fRawSP){
	    ftime2 = (sz2 + 1./12)*AliHLTTransform::GetZWidth()*AliHLTTransform::GetZWidth();
	    if(sz2 != 0) {
	      ftime2 *= 0.169; //constants are from offline studies
	      if(patch<2)
		ftime2 *= 1.77;
	    }
	  } else ftime2=sz2; //take the width, not the error
	}
      }
      if(fStdout==kTRUE)
	cout<<"WriteCluster: padrow "<<fCurrentRow<<" pad "<<fpad << " +- "<<fpad2<<" time "<<ftime<<" +- "<<ftime2<<" charge "<<list[j].fTotalCharge<<endl;
      
      if(!fRawSP){
	AliHLTTransform::Slice2Sector(fCurrentSlice,fCurrentRow,thissector,thisrow);
	AliHLTTransform::Raw2Local(xyz,thissector,thisrow,fpad,ftime);
	
	if(xyz[0]==0) LOG(AliHLTLog::kError,"AliHLTClustFinder","Cluster Finder")
	  <<AliHLTLog::kDec<<"Zero cluster"<<ENDLOG;
	if(fNClusters >= fMaxNClusters)
	  {
	    LOG(AliHLTLog::kError,"AliHLTClustFinder::WriteClusters","Cluster Finder")
	      <<AliHLTLog::kDec<<"Too many clusters "<<fNClusters<<ENDLOG;
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
void AliHLTClustFinderNew::GetTrackID(Int_t pad,Int_t time,Int_t *trackID)
{
  //get mc id
  AliHLTDigitRowData *rowPt = (AliHLTDigitRowData*)fDigitRowData;
  
  trackID[0]=trackID[1]=trackID[2]=-2;
  //cout<<"Looking for pad "<<pad<<" time "<<time<<endl;
  for(Int_t i=fFirstRow; i<=fLastRow; i++){
    if(rowPt->fRow < (UInt_t)fCurrentRow){
      AliHLTMemHandler::UpdateRowPointer(rowPt);
      continue;
    }
    AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
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
