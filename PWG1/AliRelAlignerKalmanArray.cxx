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

///////////////////////////////////////////////////////////////////////////////
//
//     Data container for relative ITS-TPC alignment analysis
//     Holds an array of AliRelAlignerKalman objects
//     and takes care of merging when processing data in parallel
//
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TCollection.h>
#include "AliESDEvent.h"
#include "AliRelAlignerKalman.h"
#include "AliRelAlignerKalmanArray.h"

ClassImp(AliRelAlignerKalmanArray)

//______________________________________________________________________________
AliRelAlignerKalmanArray::AliRelAlignerKalmanArray():
    TNamed(),
    fArray(new TObjArray()),
    fSaveInterval(600), //default every 10minutes
    fTimeMatchingTolerance(20),
    fCurrentTimeBin(0),
    fAligner(new AliRelAlignerKalman())
{
  //ctor
}

//______________________________________________________________________________
AliRelAlignerKalmanArray::AliRelAlignerKalmanArray(const char* name):
    TNamed(name, name),
    fArray(new TObjArray()),
    fSaveInterval(600), //default every 10 minutes
    fTimeMatchingTolerance(60),
    fCurrentTimeBin(0),
    fAligner(new AliRelAlignerKalman())
{
  //ctor
}

//______________________________________________________________________________
AliRelAlignerKalmanArray::~AliRelAlignerKalmanArray()
{
  //dtor
  fArray->SetOwner();
  delete fArray;
  delete fAligner;
}

//______________________________________________________________________________
Long64_t AliRelAlignerKalmanArray::Merge( TCollection* list )
{
  //Merge all the arrays
  //the merge is vertical, meaning matching entries in tree are merged

  AliRelAlignerKalmanArray *arrayFromList;
  if (!list) return 0;
  TIter next(list);
  while ( (arrayFromList = dynamic_cast<AliRelAlignerKalmanArray*>(next())) )
  {
    if (arrayFromList==this) continue;

    fArray->AddAll(arrayFromList->fArray); //put all objects in one array

    //do the merge
    fArray->Sort();
    TObjArray* tmpArray = SortedMerge(fArray);
    tmpArray->SetOwner(kTRUE);

    TObjArray* newArray = dynamic_cast<TObjArray*>(tmpArray->Clone());
    delete fArray; //takes care of all loaded objects
    fArray = newArray;

    fArray->AddLast(arrayFromList->fAligner); //add the endofrun aligner
  }

  //TODO: this can be done better!
  //Add own endofrun aligner and clean up
  fArray->AddLast(fAligner);
  fArray->Sort();
  //TObjArray* tmpArray = SortedMerge(fArray);
  //tmpArray->SetOwner(kTRUE);
  //TObjArray* newArray = dynamic_cast<TObjArray*>(tmpArray->Clone());
  //delete fArray; //takes care of all loaded objects
  //fArray = newArray;

  return fArray->GetEntriesFast();
}

//______________________________________________________________________________
TObjArray* AliRelAlignerKalmanArray::SortedMerge( TObjArray* input )
{
  //Merges the adjacent aligners if close enough
  //input needs to be already sorted

  UInt_t timeStampIn;
  AliRelAlignerKalman* alignerIn;
  AliRelAlignerKalman* alignerOut = dynamic_cast<AliRelAlignerKalman*>(input->At(0));
  TObjArray* output = new TObjArray();  //empty array
  output->AddLast(alignerOut);        //first object in: copy of first input element

  timeStampIn = alignerOut->GetTimeStamp();
  SetCurrentTimeBin( timeStampIn );
  TIter next(input);
  while ( (alignerIn = dynamic_cast<AliRelAlignerKalman*>(next())) )
  {
    timeStampIn = alignerIn->GetTimeStamp();
    if ( IsInCurrentTimeBin(timeStampIn) )
    {
      alignerOut->Merge(alignerIn);
    }
    else
    {
      alignerOut = alignerIn;
      output->AddLast(alignerOut);
    }//if
    SetCurrentTimeBin( timeStampIn );
  }
  return output;
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::SetCurrentTimeBin( UInt_t timestamp )
{
  //set the current timebin
  fCurrentTimeBin = TimeBin(timestamp);
}

//______________________________________________________________________________
UInt_t AliRelAlignerKalmanArray::TimeBin( UInt_t timestamp ) const
{
  return (timestamp+(fSaveInterval/2))/fSaveInterval*fSaveInterval; //it's all integers!
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalmanArray::IsInCurrentTimeBin( UInt_t timestamp ) const
{
  //check if timestamp is within the current timebin
  UInt_t timeDiff = (timestamp>=fCurrentTimeBin)?timestamp-fCurrentTimeBin:
                    fCurrentTimeBin-timestamp;
  return (timeDiff < fTimeMatchingTolerance);
}

////______________________________________________________________________________
//void AliRelAlignerKalmanArray::AddESDEvent( AliESDEvent* event )
//{
//  //add an AliESDEvent, take care of bookkeeping
//  if (!fAligner) return;
//  if (event->GetRunNumber() != fAligner->GetRunNumber())
//  {
//    //what to do when a new run starts
//  }
//
//}
//
////______________________________________________________________________________
Bool_t AliRelAlignerKalmanArray::AddCosmicEvent( AliESDEvent* event )
{
  if (!fAligner->AddCosmicEvent(event)) return kFALSE;

  UInt_t currentTimeStamp = event->GetTimeStamp();
  UInt_t timeFromLastBinCentre = currentTimeStamp - fCurrentTimeBin;
  UInt_t binCentre = TimeBin(currentTimeStamp);
  UInt_t timeFromBinCentre = currentTimeStamp-binCentre;
  UInt_t nIntervals = timeFromLastBinCentre/fSaveInterval;

  //////////////////////////////////////////////////////////////////////////////
  //only ONE SAVE PER TIMEBIN!!!!, as close as possible to the bin center
  //////////////////////////////////////////////////////////////////////////////
  if ( (nIntervals == 1) &&                      //is in next time bin passed centre
       (timeFromBinCentre < fTimeMatchingTolerance) )     //and close to it
  {
    AddLast(new AliRelAlignerKalman(*fAligner));
  }
  else if ( (nIntervals > 2) ) //TODO: don't hardwire stuff!
  {
    //if missed a few windows save anyway at current bin centre
    fAligner->SetTimeStamp(binCentre);
    AddLast(new AliRelAlignerKalman(*fAligner));
  }
  //////////////////////////////////////////////////////////////////////////////

  return kTRUE;
}

//______________________________________________________________________________
AliRelAlignerKalmanArray::AliRelAlignerKalmanArray( const AliRelAlignerKalmanArray& in):
    TNamed(in.GetName(), in.GetTitle()),
    fArray(NULL),
    fSaveInterval(in.fSaveInterval),
    fTimeMatchingTolerance(in.fTimeMatchingTolerance),
    fCurrentTimeBin(in.fCurrentTimeBin),
    fAligner(new AliRelAlignerKalman(*in.fAligner))
{
  //copy ctor
  fArray = static_cast<TObjArray*>(in.Clone());
}

//______________________________________________________________________________
AliRelAlignerKalmanArray& AliRelAlignerKalmanArray::operator=(const AliRelAlignerKalmanArray& in)
{
  //assignment operator
  fArray = static_cast<TObjArray*>(in.Clone());
  fSaveInterval = in.fSaveInterval;
  fTimeMatchingTolerance = in.fTimeMatchingTolerance;
  return *this;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalmanArray::SetSaveInterval( const UInt_t s )
{
  //only set if array empty
  if (fArray->GetEntriesFast()) return kFALSE;
  fSaveInterval = s;
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRelAlignerKalmanArray::SetTimeMatchingTolerance( const UInt_t m )
{
  //only set if array empty
  if (fArray->GetEntriesFast()) return kFALSE;
  fTimeMatchingTolerance = m;
  return kTRUE;
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::At( Int_t i ) const
{
  //mimic TObjArray::At( Int_t i )
  return dynamic_cast<AliRelAlignerKalman*>(fArray->At(i));
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::AddLast( AliRelAlignerKalman* al )
{
  //mimic TObjArray::AddLast( TObject* obj )
  fArray->AddLast( al );
  SetCurrentTimeBin(al->GetTimeStamp());
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::Last() const
{
  //mimic TObjArray::Last()
  return dynamic_cast<AliRelAlignerKalman*>(fArray->Last());
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::operator[](Int_t i) const
{
  //mimic TObjArray::operator[](Int_t)
  return dynamic_cast<AliRelAlignerKalman*>(fArray->At(i));
}

