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

/*********************************************************************
 *  FIT reconstruction and filling ESD
 *  - reconstruct mean time (interation time) 
 *  - vertex position
 *  -  multiplicity
 ********************************************************************/

#include <AliESDEvent.h>
#include "AliLog.h"
#include "AliFITRecPoint.h"
#include "AliFITDigit.h"
#include "AliFITReconstructor.h"
#include "AliESDFIT.h"
#include "AliLog.h"
#include "AliESDRun.h"
#include "AliFITRawReader.h"

#include <TArrayI.h>
#include <TGraph.h>
#include <TMath.h>
#include <Riostream.h>

using std::cout;
using std::endl;

ClassImp(AliFITReconstructor)

AliFITReconstructor:: AliFITReconstructor(): AliReconstructor(),
					     fESD(NULL),
					     fESDFIT(NULL),
					     fDigits(NULL)

{
 
  //constructor
 
  printf("@@@ AliFITReconstructor:: AliFITReconstructor()\n");
 
}
//_____________________________________________________________________________
void AliFITReconstructor::Init()
{
// initializer

  fESDFIT  = new AliESDFIT;
   printf("@@@ AliFITReconstructor:: Init\n");

 }

//______________________________________________________________________
void AliFITReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
// converts RAW to digits 
  cout<<" AliFITReconstructor::ConvertDigits "<<rawReader<<" "<<digitsTree<<endl;
  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  if (!fDigits)  fDigits = new TClonesArray ("AliFITDigit", 100);
  digitsTree->Branch("FIT", &fDigits);

  AliFITRawReader myrawreader(rawReader);
  if (myrawreader.Next()) 
    {
      Int_t allData[1000];
      for (Int_t i=0; i<1000; i++)  allData[i]=0;
      for (Int_t i=0; i<1000; i++) 
	if(myrawreader.GetData(i)>0)  { 
	  allData[i]=myrawreader.GetData(i);
	}
   
      Int_t timeCFD, timeLED, timeQT1, timeQT0;
      for (Int_t ipmt=0; ipmt<240; ipmt++) {
	if(allData[ipmt]>0) {
	  timeCFD = allData[ipmt];
	  timeLED = allData[ipmt];
	  timeQT0= allData[ipmt+240];
	  timeQT1 = allData[ipmt+480];
	  //add digit
	  new((*fDigits)[fDigits->GetEntriesFast()] )
	    AliFITDigit( ipmt,timeCFD, timeLED, timeQT0, timeQT1, 0);
	}
      }
      
      
      
      digitsTree->Fill(); 
    }

  
}

 //____________________________________________________________
void AliFITReconstructor::FillESD(TTree *digitsTree, TTree * /*clustersTree*/, AliESDEvent *pESD) const
{
  
  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  AliDebug(1,Form("Start FillESD FIT"));
  if(!pESD) {
    AliError("No ESD Event");
    return;
  }
  Float_t channelWidth = 24.4;  
  Float_t c = 0.0299792458; // cm/ps
  Float_t currentVertex, shift=0;
  Int_t ncont=-1;
  const AliESDVertex* vertex = pESD->GetPrimaryVertex();
  if (!vertex)        vertex = pESD->GetPrimaryVertexSPD();
  if (!vertex)        vertex = pESD->GetPrimaryVertexTPC();
  if (!vertex)        vertex = pESD->GetVertex();
  
  if (vertex) {
    AliDebug(2, Form("Got %s (%s) from ESD: %f", 
		     vertex->GetName(), vertex->GetTitle(), vertex->GetZ()));
    currentVertex = vertex->GetZ();
    
    ncont = vertex->GetNContributors();
    if(ncont>0 ) {
      shift = currentVertex/c;
    }
  } //vertex
  
  
  // FIT digits reconstruction
  
  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }
  
  fDigits = new TClonesArray("AliFITDigit", 100);
  // digitsTree->Print();
  TBranch *brDigits=digitsTree->GetBranch("FIT");
  if (brDigits) 
    brDigits->SetAddress(&fDigits);
  else {
    AliError("no FIT branch");
    return;
  }
  const Float_t max_time = 1e7;
  Int_t pmt=-1;
  Float_t  time[240],amp[240];  
  for(Int_t i=0; i<240; i++)   {  time[i]=max_time; amp[i]=0;}
  
  Int_t nEntries = (Int_t)digitsTree->GetEntries();
  for (Int_t iev=0; iev<nEntries; iev++) {
    digitsTree->GetEvent(iev,1);
    brDigits->GetEvent(iev);
    Int_t nDigits = fDigits->GetEntriesFast(); 
    for (Int_t dig=0; dig<nDigits; dig++) {    
      AliFITDigit* digit = (AliFITDigit*) fDigits->At(dig);      
      pmt = digit->NPMT();
      time[pmt] = Float_t (digit->TimeCFD() );
      amp[pmt] = 0.001 * Float_t (digit->TimeQT1() - digit->TimeQT0() );
    } 
    fESDFIT->SetFITtime(time);         // best TOF on each PMT 
    fESDFIT->SetFITamplitude(amp);     // number of particles(MIPs) on each 
     Float_t firsttime[3] = {max_time,max_time,max_time};
    
    Float_t vertexFIT = 9999999;
    for (Int_t ipmt=0; ipmt<96; ipmt++)//timeC
      if(time[ipmt]<firsttime[2]) firsttime[2]=time[ipmt]; 
    
    for ( Int_t ipmt=96; ipmt<240; ipmt++) 
      if(time[ipmt]<firsttime[1]) firsttime[1]=time[ipmt]; 
    if (firsttime[1]<max_time && firsttime[2]<max_time)  {
      firsttime[0] =  channelWidth *(firsttime[1] + firsttime[2])/2;
      vertexFIT =  c*channelWidth*(firsttime[1] - firsttime[2])/2;
    }
    if (firsttime[1]<max_time) 
      firsttime[1] = firsttime[1] * channelWidth + shift; 
    
    if (firsttime[2]<max_time) 
      firsttime[2] = firsttime[2] * channelWidth - shift; 
             
    for(Int_t i=0; i<3; i++) {
      fESDFIT->SetFITT0(i,firsttime[i]);   // interaction time (ns) 
      }
    fESDFIT->SetFITzVertex(vertexFIT); //vertex Z position 
        
    AliDebug(1,Form("FIT: SPDshift %f Vertex %f  FITsignal %f ps FITA %f ps FITC %f ps \n",shift, vertexFIT, firsttime[0], firsttime[1],firsttime[2]));
  }
     if (pESD)    pESD->SetFITData(fESDFIT); 
}
