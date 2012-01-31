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

/*****************************************************************
  AliStarEvent: Event container for star data

  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#include <string.h>
#include <TObjArray.h>
#include "AliStarTrack.h"
#include "AliStarEvent.h"

ClassImp(AliStarEvent)

//______________________________________________________________________________
AliStarEvent::AliStarEvent():
  TObject(),
  fParams(),
  fTracks(NULL)
{
  //ctor
}

//______________________________________________________________________________
AliStarEvent::AliStarEvent( Int_t n ):
  TObject(),
  fParams(),
  fTracks(new TObjArray(n) )
{
  //ctor
}

//______________________________________________________________________________
AliStarEvent::AliStarEvent( const AliStarEvent& event ):
  TObject(event),
  fParams(),
  fTracks((TObjArray*)(event.fTracks)->Clone())
{
  //copy ctor
  memcpy(fParams,event.fParams,fgkNparams*sizeof(Float_t));
}

//______________________________________________________________________________
AliStarEvent& AliStarEvent::operator=( const AliStarEvent& event )
{
  //assignment
  if (this == &event) return *this;
  TObject::operator=(event);
  memcpy(fParams,event.fParams,fgkNparams*sizeof(Float_t));
  if (fTracks) fTracks->Delete();
  delete fTracks;
  if (event.fTracks)
    fTracks = (TObjArray*)(event.fTracks)->Clone(); //deep copy
  return *this;
}

//______________________________________________________________________________
void AliStarEvent::SetParams( const Float_t* params )
{
  //set params
  memcpy(fParams,params,fgkNparams*sizeof(Float_t));
}

//______________________________________________________________________________
AliStarEvent::~AliStarEvent()
{
  //dtor
  if (fTracks) fTracks->Delete();
  delete fTracks;
}

//______________________________________________________________________________
const AliStarTrack* AliStarEvent::GetTrack( const Int_t i ) const
{
  //get track i
  if (!fTracks) return NULL;
  return static_cast<AliStarTrack*>(fTracks->At(i));
}

//______________________________________________________________________________
void AliStarEvent::AddTrack( AliStarTrack* track )
{
  //add a new track to collection
  if (!fTracks) return;
  fTracks->Add(track);
}

//______________________________________________________________________________
void AliStarEvent::Reset()
{
  //remove all tracks if any and zero the event information
  if (!fTracks) 
  {
    fTracks=new TObjArray(1024);
    return;
  }
  fTracks->Delete();
  for (Int_t i=0;i<fgkNparams; i++){fParams[i]=0;}
}

//______________________________________________________________________________
void AliStarEvent::Print( Option_t* option ) const
{
  //print info
  // TNtuple* event: names are documented in the next 2 lines
  // event  = new TNtuple("event","event",
  //   "runId:eventNumber:VtxX:VtxY:VtxZ:BField:refMult:centralityId:numberOfPrimaryTracks:numberOfParticles:h0:h1:h2:h3:h4" ) ;
  //

  Int_t   runId                  = (Int_t)   fParams[0]  ;
  Int_t   eventNumber            = (Int_t)   fParams[1]  ;     
  Float_t   primaryVertexPosition[3]={ (Float_t) fParams[2],
                                       (Float_t) fParams[3],
                                       (Float_t) fParams[4] };
  Float_t magneticField          = (Float_t) fParams[5]  ;  // kilogauss
  Int_t   refMult                = (Int_t)   fParams[6]  ;  // Raw Mult into 0.5 unit: a relative number, not total Mult.
  Int_t   centralityId           = (Int_t)   fParams[7]  ;  // STAR centrality bin # based on refMult
  Int_t   numberOfPrimaryTracks  = (Int_t)   fParams[8]  ;  // # of primaries, including FTPC tracks which are not in ntuple
  Int_t   numberOfParticles      = (Int_t)   fParams[9]  ;  // # of particles passing track cuts, thus in ntuple
  
  printf("\n  runId eventNo    VtxX    VtxY    VtxZ  MagFld  refMult  CentBin  #PrimeTrak  #Tracks \n") ;
  printf("%7d  %6d %7.4f %7.4f %7.3f  %6.3f   %6d   %6d      %6d   %6d \n\n", 
	 runId, eventNumber, primaryVertexPosition[0], primaryVertexPosition[1], primaryVertexPosition[2], 
	 magneticField, refMult, centralityId, numberOfPrimaryTracks, numberOfParticles ) ;
  
  //Int_t newCentralityID ;
  //newCentralityID = Centrality( refMult) ;              // Should be the same as "centralityID" from tape
  //cout << "Test: should be the same " << newCentralityID << " " << centralityId << endl ;  // JT test

  TString optionstr(option);
  if (!optionstr.Contains("all")) return;
  if (!fTracks) return;
  AliStarTrack* track = static_cast<AliStarTrack*>( fTracks->First() );
  if (!track) return;
  track->Print("legend");
  fTracks->Print();
}

//______________________________________________________________________________
Int_t AliStarEvent::CalculateCentrality( Int_t referenceMultiplicity ) const
{
  // Note Carefully:  Centrality is based on refMult.  This is the 'reference' multiplicity that is measured
  // independpently from the TPC.  Selecting the centrality bins according to the refMult is something that
  // is calibrated for each year and each run.  You can get the basic information off the web:
  // For Example .... http://www.star.bnl.gov/protected/common/common2004/trigger2004/200gev/200gevFaq.html
  // An index pointing to FAQs, Trigger and Centrality data, for all years, is available at:
  // http://www.star.bnl.gov/public/all
  //
  // Note: Add 0.5 to the (int) that is returned by centrality() when using it as an argument for a histogram
  // that expects (float) or (double) as input parameters.  This will place the data point in the center of
  // the bin, avoids ambiguities, and is best for plotting scatter plots and contour plots.
  // For example histogram2D[1] -> Fill ( (float)CentralityID + 0.5 , SumData )   ;
  //
  // The refMult quoted in the Centrality bins array is the lower limit on refMult


  Int_t   centralityBins  [] = { 14 , 31 , 57 , 96 , 150 , 222 , 319 , 441 , 520 , 1000 } ;  // Run4 200 GeV
  Int_t   middleBinID     [] = {  0 ,  1 ,  2 ,  3 ,   4 ,   5 ,   6 ,   7 ,   8 ,    9 } ;  // ID Number
  //Info  middleBinPercent[] = { 85., 75., 65., 55.,  45.,  35.,  25.,  15., 7.5 ,  2.5 } ;  // Percent
  Int_t   myCentrality  ;

  if      ( referenceMultiplicity < centralityBins[0] )
  {
    myCentrality = middleBinID[0] ;
  }
  else if ( referenceMultiplicity < centralityBins[1] )
  {
    myCentrality = middleBinID[1] ;
  }
  else if ( referenceMultiplicity < centralityBins[2] )
  {
    myCentrality = middleBinID[2] ;
  }
  else if ( referenceMultiplicity < centralityBins[3] )
  {
    myCentrality = middleBinID[3] ;
  }
  else if ( referenceMultiplicity < centralityBins[4] )
  {
    myCentrality = middleBinID[4] ;
  }
  else if ( referenceMultiplicity < centralityBins[5] )
  {
    myCentrality = middleBinID[5] ;
  }
  else if ( referenceMultiplicity < centralityBins[6] )
  {
    myCentrality = middleBinID[6] ;
  }
  else if ( referenceMultiplicity < centralityBins[7] )
  {
    myCentrality = middleBinID[7] ;
  }
  else if ( referenceMultiplicity < centralityBins[8] )
  {
    myCentrality = middleBinID[8] ;
  }
  else
  {
    myCentrality = middleBinID[9] ;
  }

  return myCentrality ;
}


