// @(#) $Id$
// Original: AliHLTClustFinderNew.cxx,v 1.29 2005/06/14 10:55:21 cvetan Exp 

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Anders Vestbo <mailto:vestbo@fi.uib.no>, 		          *
 *          Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>    *
 *          Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>         *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCClusterFinder.cxx
    @author Anders Vestbo <mailto:vestbo@fi.uib.no>, 		     
	    Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
	    Jochen Thaeder <mailto:thaeder@kip.uni-heidelberg.de>     
    @date   
    @brief  Cluster Finder for the TPC
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCClusterFinder.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCPad.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/** \class AliHLTTPCClusterFinder
//
// The current cluster finder for HLT
// (Based on STAR L3)
// 
// The cluster finder is initialized with the Init function, 
// providing the slice and patch information to work on. 
//
// The input is a provided by the AliHLTTPCDigitReader class,
// using the init() funktion, and the next() funktion in order 
// to get the next bin. Either packed or unpacked data can be
// processed, dependent if one uses AliHLTTPCDigitReaderPacked 
// class or AliHLTTPCDigitReaderUnpacked class in the 
// Clusterfinder Component.
// The resulting space points will be in the
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
//     AliHLTTPCClusterFinder *cf = new AliHLTTPCClusterFinder();
//     cf->SetMatchWidth(2);
//     cf->InitSlice( slice, patch, row[0], row[1], maxPoints );
//     cf->SetSTDOutput(kTRUE);    //Some output to standard IO
//     cf->SetRawSP(kFALSE);       //Convert space points to local system
//     cf->SetThreshold(5);        //Threshold of cluster charge
//     cf->SetDeconv(kTRUE);       //Deconv in pad and time direction
//     cf->SetCalcErr(kTRUE);      //Calculate the errors of the spacepoints
//     cf->SetOutputArray(points); //Move the spacepoints to the array
//     cf->Read(iter->fPtr, iter->fSize ); //give the data to the cf
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
*/

ClassImp(AliHLTTPCClusterFinder)

AliHLTTPCClusterFinder::AliHLTTPCClusterFinder()
  :
  fMatch(1),
  fThreshold(10),
  fSignalThreshold(-1),
  fOccupancyLimit(1.0),
  fXYErr(0.2),
  fZErr(0.3),
  fDeconvPad(kTRUE),
  fDeconvTime(kTRUE),
  fStdout(kFALSE),
  fCalcerr(kTRUE),
  fRawSP(kFALSE),
  fFirstRow(0),
  fLastRow(0),
  fDigitReader(NULL)
{
  //constructor
}

AliHLTTPCClusterFinder::AliHLTTPCClusterFinder(const AliHLTTPCClusterFinder& src)
  :
  fMatch(src.fMatch),
  fThreshold(src.fThreshold),
  fSignalThreshold(src.fSignalThreshold),
  fOccupancyLimit(src.fOccupancyLimit),
  fXYErr(src.fXYErr),
  fZErr(src.fZErr),
  fDeconvPad(src.fDeconvPad),
  fDeconvTime(src.fDeconvTime),
  fStdout(src.fStdout),
  fCalcerr(src.fCalcerr),
  fRawSP(src.fRawSP),
  fFirstRow(src.fFirstRow),
  fLastRow(src.fLastRow),
  fDigitReader(src.fDigitReader)
{
}

AliHLTTPCClusterFinder& AliHLTTPCClusterFinder::operator=(const AliHLTTPCClusterFinder& src)
{
  fMatch=src.fMatch;
  fThreshold=src.fThreshold;
  fSignalThreshold=src.fSignalThreshold;
  fOccupancyLimit=src.fOccupancyLimit;
  fXYErr=src.fXYErr;
  fZErr=src.fZErr;
  fDeconvPad=src.fDeconvPad;
  fDeconvTime=src.fDeconvTime;
  fStdout=src.fStdout;
  fCalcerr=src.fCalcerr;
  fRawSP=src.fRawSP;
  fFirstRow=src.fFirstRow;
  fLastRow=src.fLastRow;
  fDigitReader=src.fDigitReader;
  return (*this);
}

AliHLTTPCClusterFinder::~AliHLTTPCClusterFinder()
{
  //destructor
}
 
void AliHLTTPCClusterFinder::InitSlice(Int_t slice,Int_t patch,Int_t firstrow, Int_t lastrow,Int_t nmaxpoints)
{
  //init slice
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow = firstrow;
  fLastRow = lastrow;
}

void AliHLTTPCClusterFinder::InitSlice(Int_t slice,Int_t patch,Int_t nmaxpoints)
{
  //init slice
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow=AliHLTTPCTransform::GetFirstRow(patch);
  fLastRow=AliHLTTPCTransform::GetLastRow(patch);
}

void AliHLTTPCClusterFinder::SetOutputArray(AliHLTTPCSpacePointData *pt)
{
  //set pointer to output
  fSpacePointData = pt;
}

void AliHLTTPCClusterFinder::Read(void* ptr,unsigned long size){
  //set input pointer
  fPtr = (UChar_t*)ptr;
  fSize = size;
}

void AliHLTTPCClusterFinder::ProcessDigits()
{
  bool readValue = true;
  Int_t newRow = 0;    
  Int_t rowOffset = 0;
  UShort_t time=0,newTime=0;
  UInt_t pad=0,newPad=0;
  AliHLTTPCSignal_t charge=0;

  fNClusters = 0;

  // initialize block for reading packed data
  fDigitReader->InitBlock(fPtr,fSize,fFirstRow,fLastRow,fCurrentPatch,fCurrentSlice);
  readValue = fDigitReader->Next();

  // Matthias 08.11.2006 the following return would cause termination without writing the
  // ClusterData and thus would block the component. I just want to have the commented line
  // here for information
  //if (!readValue)return;

  pad = fDigitReader->GetPad();
  time = fDigitReader->GetTime();
  fCurrentRow = fDigitReader->GetRow();

  if ( fCurrentPatch >= 2 ) // Outer sector, patches 2, 3, 4, 5
    rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

  fCurrentRow += rowOffset;

  UInt_t lastpad = 123456789;
  const Int_t kPadArraySize=5000;
  const Int_t kClusterListSize=10000;
  AliClusterData *pad1[kPadArraySize]; //2 lists for internal memory=2pads
  AliClusterData *pad2[kPadArraySize]; //2 lists for internal memory=2pads
  AliClusterData clusterlist[kClusterListSize]; //Clusterlist

  AliClusterData **currentPt;  //List of pointers to the current pad
  AliClusterData **previousPt; //List of pointers to the previous pad
  currentPt = pad2;
  previousPt = pad1;
  UInt_t nprevious=0,ncurrent=0,ntotal=0;

  /* quick implementation of baseline calculation and zero suppression
     open a pad object for each pad and delete it after processing.
     later a list of pad objects with base line history can be used
     The whole thing only works if we really get unprocessed raw data, if
     the data is already zero suppressed, there might be gaps in the time
     bins.
   */
  Int_t gatingGridOffset=50;
  AliHLTTPCPad baseline(gatingGridOffset, AliHLTTPCTransform::GetNTimeBins());
  // just to make later conversion to a list of objects easier
  AliHLTTPCPad* pCurrentPad=NULL;
  if (fSignalThreshold>=0) {
    pCurrentPad=&baseline;
    baseline.SetThreshold(fSignalThreshold);
  }

  while ( readValue ){   // Reads through all digits in block

    if(pad != lastpad){
      //This is a new pad
      
      //Switch the lists:
      if(currentPt == pad2){
	currentPt = pad1;
	previousPt = pad2;
      }
      else {
	currentPt = pad2;
	previousPt = pad1;
      }
      nprevious = ncurrent;
      ncurrent = 0;
      if(pad != lastpad+1){
	//this happens if there is a pad with no signal.
	nprevious = ncurrent = 0;
      }
      lastpad = pad;
    }

    Bool_t newcluster = kTRUE;
    UInt_t seqcharge=0,seqaverage=0,seqerror=0;
    UInt_t lastcharge=0,lastwas_falling=0;
    Int_t newbin=-1;


    if(fDeconvTime){
      redo: //This is a goto.
      
      if(newbin > -1){
	//bin = newbin;
	newbin = -1;
      }
	  
      lastcharge=0;
      lastwas_falling = 0;
    }

    while(1){ //Loop over time bins of current pad
      // read all the values for one pad at once to calculate the base line
      if (pCurrentPad) {
	if (!pCurrentPad->IsStarted()) {
	  //HLTDebug("reading data for pad %d, padrow %d", fDigitReader->GetPad(), fDigitReader->GetRow()+rowOffset);
	  pCurrentPad->SetID(fDigitReader->GetRow()+rowOffset,fDigitReader->GetPad());
	  if ((pCurrentPad->StartEvent())>=0) {
	    do {
	      if ((fDigitReader->GetRow()+rowOffset)!=pCurrentPad->GetRowNumber()) break;
	      if (fDigitReader->GetPad()!=pCurrentPad->GetPadNumber()) break;
	      pCurrentPad->SetRawData(fDigitReader->GetTime(), fDigitReader->GetSignal());
	      //HLTDebug("set raw data to pad: bin %d charge %d", fDigitReader->GetTime(), fDigitReader->GetSignal());
	    } while ((readValue = fDigitReader->Next())!=0);
	  }
	  pCurrentPad->CalculateBaseLine(AliHLTTPCTransform::GetNTimeBins()/2);
	  if (pCurrentPad->Next(kTRUE/*do zero suppression*/)==0) {
	    //HLTDebug("no data available after zero suppression");
	    pCurrentPad->StopEvent();
	    pCurrentPad->ResetHistory();
	    break;
	  }
	  time=pCurrentPad->GetCurrentPosition();
	  if (time>pCurrentPad->GetSize()) {
	    HLTError("invalid time bin for pad");
	    break;
	  }
	}
      }

      if (pCurrentPad) {
	Float_t occupancy=pCurrentPad->GetOccupancy();
	//HLTDebug("pad %d occupancy level: %f", pCurrentPad->GetPadNumber(), occupancy);
	if ( occupancy < fOccupancyLimit ) {
	  charge = pCurrentPad->GetCorrectedData();
	} else {
	  charge = 0;
	  //HLTDebug("ignoring pad %d with occupancy level %f", pCurrentPad->GetPadNumber(), occupancy);
	}
      } else {
	charge = fDigitReader->GetSignal();
      }
      //HLTDebug("get next charge value: position %d charge %d", time, charge);


      // CHARGE DEBUG
      if (fDigitReader->GetRow() == 90){
/////	  LOG(AliHLTTPCLog::kFatal,"AliHLTTPCClusterFinder::Row","row90")  << "PAD=" <<  fDigitReader->GetPad() << "  TIME=" <<  fDigitReader->GetTime() 
	  //					   << "  SIGNAL=" <<  fDigitReader->GetSignal() << ENDLOG;

      }

      if(time >= AliHLTTPCTransform::GetNTimeBins()){
	HLTWarning("Timebin (%d) out of range (%d)", time, AliHLTTPCTransform::GetNTimeBins());
	break;
      }


      //Get the current ADC-value
      if(fDeconvTime){

	//Check if the last pixel in the sequence is smaller than this
	if(charge > lastcharge){
	  if(lastwas_falling){
	    newbin = 1;
	    break;
	  }
	}
	else lastwas_falling = 1; //last pixel was larger than this
	lastcharge = charge;
      }
	  
      //Sum the total charge of this sequence
      seqcharge += charge;
      seqaverage += time*charge;
      seqerror += time*time*charge;
      
      if (pCurrentPad) {
	
	if((pCurrentPad->Next(kTRUE/*do zero suppression*/))==0) {
	  pCurrentPad->StopEvent();
	  pCurrentPad->ResetHistory();
	  if(readValue) {
	    newPad = fDigitReader->GetPad();
	    newTime = fDigitReader->GetTime();
	    newRow = fDigitReader->GetRow() + rowOffset;
	  }
	  break;
	}

	newPad=pCurrentPad->GetPadNumber();
	newTime=pCurrentPad->GetCurrentPosition();
	newRow=pCurrentPad->GetRowNumber();
      } else {
      readValue = fDigitReader->Next();
      //Check where to stop:
      if(!readValue) break; //No more value

      newPad = fDigitReader->GetPad();
      newTime = fDigitReader->GetTime();
      newRow = fDigitReader->GetRow() + rowOffset;
      }

      if(newPad != pad)break; //new pad
      if(newTime != time+1) break; //end of sequence

      // pad = newpad;    is equal
      time = newTime;

    }//end loop over sequence

    //HLTDebug("ended time bin sequence loop: seqcharge=%d readValue=%d", seqcharge, readValue);
    //HLTDebug("pad=%d newpad=%d current row=%d newrow=%d", pad, newPad, fCurrentRow, newRow);
    if (seqcharge<=0) {
      // with active zero suppression zero values are possible
      continue;
    }

    //Calculate mean of sequence:
    Int_t seqmean=0;
    if(seqcharge)
      seqmean = seqaverage/seqcharge;
    else{
      LOG(AliHLTTPCLog::kFatal,"AliHLTTPCClusterFinder::ProcessRow","Data")
	<<"Error in data given to the cluster finder"<<ENDLOG;
      seqmean = 1;
      seqcharge = 1;
    }

    //Calculate mean in pad direction:
    Int_t padmean = seqcharge*pad;
    Int_t paderror = pad*padmean;


    //Compare with results on previous pad:
    for(UInt_t p=0; p<nprevious && p<kPadArraySize && ncurrent<kPadArraySize; p++){
      
      //dont merge sequences on the same pad twice
      if(previousPt[p]->fLastMergedPad==pad) continue;

      Int_t difference = seqmean - previousPt[p]->fMean;
      if(difference < -fMatch) break;

      if(difference <= fMatch){ //There is a match here!!
	AliClusterData *local = previousPt[p];
	
	if(fDeconvPad){
	  if(seqcharge > local->fLastCharge){
	    if(local->fChargeFalling){ //The previous pad was falling
	      break; //create a new cluster
	    }		    
	  }
	  else local->fChargeFalling = 1;
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
	local->fLastMergedPad = pad;

	currentPt[ncurrent] = local;
	ncurrent++;
	      
	break;
      } //Checking for match at previous pad
    } //Loop over results on previous pad.


    if(newcluster && ncurrent<kPadArraySize){
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
      tmp->fLastMergedPad = pad;

      if(fDeconvPad){
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
  
    // to prevent endless loop  
    if(time >= AliHLTTPCTransform::GetNTimeBins()){
      HLTWarning("Timebin (%d) out of range (%d)", time, AliHLTTPCTransform::GetNTimeBins());
      break;
    }


    if(!readValue) break; //No more value
    
    if (ntotal>=kClusterListSize || ncurrent>=kPadArraySize) {
      HLTWarning("pad array size exceeded ntotal=%d ncurrent=%d, skip rest of the data", ntotal, ncurrent);
      break;
    }

    if(fCurrentRow != newRow){
      WriteClusters(ntotal,clusterlist);

      lastpad = 123456789;

      currentPt = pad2;
      previousPt = pad1;
      nprevious=0;
      ncurrent=0;
      ntotal=0;
      
      fCurrentRow = newRow;
    }

    pad = newPad;
    time = newTime;

  } // END while(readValue)

  WriteClusters(ntotal,clusterlist);

  HLTInfo("ClusterFinder found %d clusters in slice %d patch %d", fNClusters, fCurrentSlice, fCurrentPatch);

} // ENDEND

void AliHLTTPCClusterFinder::WriteClusters(Int_t nclusters,AliClusterData *list)
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
	Int_t patch = AliHLTTPCTransform::GetPatch(fCurrentRow);
	UInt_t q2=list[j].fTotalCharge*list[j].fTotalCharge;
	Float_t sy2=list[j].fPad2 * list[j].fTotalCharge - list[j].fPad * list[j].fPad;
	sy2/=q2;
	if(sy2 < 0) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCClusterFinder::WriteClusters","Cluster width")
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
	Float_t sz2=list[j].fTime2*list[j].fTotalCharge - list[j].fTime*list[j].fTime;
	sz2/=q2;
	if(sz2 < 0){
	  LOG(AliHLTTPCLog::kError,"AliHLTTPCClusterFinder::WriteClusters","Cluster width")
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
	cout<<"WriteCluster: padrow "<<fCurrentRow<<" pad "<<fpad << " +- "<<fpad2<<" time "<<ftime<<" +- "<<ftime2<<" charge "<<list[j].fTotalCharge<<endl;
      
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

      fSpacePointData[counter].fUsed = kFALSE;         // only used / set in AliHLTTPCDisplay
      fSpacePointData[counter].fTrackN = -1;           // only used / set in AliHLTTPCDisplay

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

// STILL TO FIX  ----------------------------------------------------------------------------

#ifdef do_mc
void AliHLTTPCClusterFinder::GetTrackID(Int_t pad,Int_t time,Int_t *trackID)
{
  //get mc id
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
