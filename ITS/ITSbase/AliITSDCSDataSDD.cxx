/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////
// Implementation of the class containing SDD DCS data           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//         V. Pospisil, CTU Praguem gdermog@seznam.cz            //
///////////////////////////////////////////////////////////////////

#include "AliITSDCSDataSDD.h"
#include "AliLog.h"

#define AUTORESIZE 50

ClassImp(AliITSDCSDataSDD)

const Float_t AliITSDCSDataSDD::fgkTPrec = 100.0;
const Float_t AliITSDCSDataSDD::fgkMVPrec = 1000.0;

//---------------------------------------------------------------------------
AliITSDCSDataSDD::AliITSDCSDataSDD(): TObject(),
fTempLeft(0),
fTempLeftTimeStamp(0),
fTempLeftMaxPoints(0),
fTempLeftSetPoints(0),
fTempRight(0),
fTempRightTimeStamp(0),
fTempRightMaxPoints(0),
fTempRightSetPoints(0),
fHV(0),
fHVTimeStamp(0),
fHVMaxPoints(0),
fHVSetPoints(0),
fMV(0),
fMVTimeStamp(0),
fMVMaxPoints(0),
fMVSetPoints(0),
fStatus(0),
fStatusTimeStamp(0),
fStatusMaxPoints(0),
fStatusSetPoints(0)
{
// default constructor
} /*AliITSDCSDataSDD::AliITSDCSDataSDD*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetNPointsTempLeft( Int_t npts )
{
  // dimension arrays with left side temperatures
  //

  if( npts < fTempLeftSetPoints)
  {                     // Cannot resize arrays - some elements would be lost
    AliWarning("Attemp to reduce size of full array (SDD DCS _TEMP_L)"); 
    return;
  } /*if*/

  fTempLeft.Set( npts );
  fTempLeftTimeStamp.Set( npts );
                        // Both temperature and tme stamp arrays are resized

  fTempLeftMaxPoints = npts;
                        // New size is stored
} /*AliITSDCSDataSDD::SetNPointsTempLeft*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetNPointsTempRight( Int_t npts )
{
  // dimension arrays with right side temperatures
  //

  if( npts < fTempRightSetPoints)
  {                     // Cannot resize arrays - some elements would be lost
    AliWarning("Attemp to reduce size of full array (SDD DCS _TEMP_R)");
    return;
  } /*if*/

  fTempRight.Set( npts );
  fTempRightTimeStamp.Set( npts );
                        // Both temperature and tme stamp arrays are resized

  fTempRightMaxPoints = npts;
                        // New size is stored
} /*AliITSDCSDataSDD::SetNPointsTempRight*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetNPointsHV( Int_t npts )
{
  // dimension arrays with HV values
  //

  if( npts < fHVSetPoints)
  {                     // Cannot resize arrays - some elements would be lost
    AliWarning("Attemp to reduce size of full array (SDD DCS _HV)");
    return;
  } /*if*/

  fHV.Set( npts );
  fHVTimeStamp.Set( npts );
                        // Both temperature and tme stamp arrays are resized

  fHVMaxPoints = npts;  // New size is stored

}/*AliITSDCSDataSDD::SetNPointsHV*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetNPointsMV( Int_t npts )
{
  // dimension arrays with   MV values
  //

  if( npts < fMVSetPoints)
  {                     // Cannot resize arrays - some elements would be lost
    AliWarning("Attemp to reduce size of full array (SDD DCS _MV)");
    return;
  } /*if*/

  fMV.Set( npts );
  fMVTimeStamp.Set( npts );
                        // Both temperature and tme stamp arrays are resized

  fMVMaxPoints = npts;  // New size is stored 

} /*AliITSDCSDataSDD::SetNPointsMV*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetNPointsStatus( Int_t npts )
{
  // dimension arrays withn DCS channel status
  //

  if( npts < fStatusSetPoints)
  {                     // Cannot resize arrays - some elements would be lost
    AliWarning("Attemp to reduce size of full array (SDD DCS Status)");
    return;
  } /*if*/

  fStatus.Set( npts );
  fStatusTimeStamp.Set( npts );
                        // Both temperature and tme stamp arrays are resized

  fStatusMaxPoints = npts;
                        // New size is stored 

} /*AliITSDCSDataSDD::SetNPointsMV*/

//---------------------------------------------------------------------------
void AliITSDCSDataSDD::SetValueTempLeft(Int_t time, Float_t temperature )
{
  // insert a value for left temperature data point
  //

   if( fTempLeftMaxPoints == fTempLeftSetPoints )
    SetNPointsTempLeft( fTempLeftMaxPoints + AUTORESIZE );
                        // Enlarges arrays if necessary

   Int_t i =  FindIndex( time, fTempLeftTimeStamp, fTempLeftSetPoints );
                        // Finds place where the new value have to be inserted

   if( i < 0 )
    i = 0;              // New value have to be inserted before the first one in the array
   else
    i++;                // New value will be put somewhere in the middle of the array
                        //  or at the end

   if( i < fTempLeftSetPoints )
   {
      Short_t *fromPtrF = fTempLeft.GetArray() + i;
                       // Sets pointer to cell which have to be filled by new value
      Short_t *toPtrF = fromPtrF + 1;
                       // Sets pointer to cell where the array content have to be shifted 

      memmove( toPtrF, fromPtrF, (fTempLeftSetPoints - i)*sizeof(Short_t) );
                       // Shifts array content. Now there is vacant place for new value to be inserted

      Int_t *fromPtrI = fTempLeftTimeStamp.GetArray() + i;
      Int_t *toPtrI = fromPtrI + 1;
      memmove( toPtrI, fromPtrI, (fTempLeftSetPoints - i)*sizeof(Int_t) );
                       // Do the same for time stamp array
   } /*if*/

   UShort_t val = (UShort_t)( temperature * fgkTPrec );
                       // Float value of temperature is stored as UShort_t with given precision

   fTempLeft.AddAt( (Short_t)val, i );
   fTempLeftTimeStamp.AddAt( time, i );
   fTempLeftSetPoints++;
                       // New values are inserted
} /*AliITSDCSDataSDD::SetValueTempLeft*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetValueTempRight(Int_t time, Float_t temperature )
{
  // insert a value for right temperature data point
  //

   if( fTempRightMaxPoints == fTempRightSetPoints )
    SetNPointsTempRight( fTempRightMaxPoints + AUTORESIZE );
                        // Enlarges arrays if necessary

   Int_t i =  FindIndex( time, fTempRightTimeStamp, fTempRightSetPoints );
                        // Finds place where the new value have to be inserted

   if( i < 0 )
    i = 0;              // New value have to be inserted before the first one in the array
   else
    i++;                // New value will be put somewhere in the middle of the array
                        //  or at the end

   if( i < fTempRightSetPoints )
   {                    // Some values have to be moved
      Short_t *fromPtrF = fTempRight.GetArray() + i;
                       // Sets pointer to cell which have to be filled by new value
      Short_t *toPtrF = fromPtrF + 1;
                       // Sets pointer to cell where the array content have to be shifted 

      memmove( toPtrF, fromPtrF, (fTempRightSetPoints - i)*sizeof(Short_t) );
                       // Shifts array content. Now there is vacant place for new value to be inserted

      Int_t *fromPtrI = fTempRightTimeStamp.GetArray() + i;
      Int_t *toPtrI = fromPtrI + 1;
      memmove( toPtrI, fromPtrI, (fTempRightSetPoints - i)*sizeof(Int_t) );
                       // Do the same for time stamp array
   } /*if*/

   UShort_t val = (UShort_t)( temperature * fgkTPrec );
                       // Float value of temperature is stored as UShort_t with given precision

   fTempRight.AddAt( (Short_t)val, i );
   fTempRightTimeStamp.AddAt( time, i );
   fTempRightSetPoints++;
                       // New values are inserted
} /*AliITSDCSDataSDD::SetValueTempRight*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetValueHV(Int_t time, Float_t voltage )
{
  // insert a value for HV data point
  //

   if( fHVMaxPoints == fHVSetPoints )
    SetNPointsHV( fHVMaxPoints + AUTORESIZE );
                        // Enlarges arrays if necessary

   Int_t i =  FindIndex( time, fHVTimeStamp, fHVSetPoints );
                        // Finds place where the new value have to be inserted

   if( i < 0 )
    i = 0;              // New value have to be inserted before the first one in the array
   else
    i++;                // New value will be put somewhere in the middle of the array
                        //  or at the end

   if( i < fHVSetPoints )
   {
      Float_t *fromPtrF = fHV.GetArray() + i;
                       // Sets pointer to cell which have to be filled by new value
      Float_t *toPtrF = fromPtrF + 1;
                       // Sets pointer to cell where the array content have to be shifted 

      memmove( toPtrF, fromPtrF, (fHVSetPoints - i)*sizeof(Float_t) );
                       // Shifts array content. Now there is vacant place for new value to be inserted

      Int_t *fromPtrI = fHVTimeStamp.GetArray() + i;
      Int_t *toPtrI = fromPtrI + 1;
      memmove( toPtrI, fromPtrI, (fHVSetPoints - i)*sizeof(Int_t) );
                       // Do the same for time stamp array
   } /*if*/

   fHV.AddAt( voltage, i );
   fHVTimeStamp.AddAt( time, i );
   fHVSetPoints++;
                       // New values are inserted
} /*AliITSDCSDataSDD::SetValueHV*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetValueMV(Int_t time, Float_t voltage )
{
  // insert a value for MV data point
  //

   if( fMVMaxPoints == fMVSetPoints )
    SetNPointsMV( fMVMaxPoints + AUTORESIZE );
                        // Enlarges arrays if necessary

   Int_t i =  FindIndex( time, fMVTimeStamp, fMVSetPoints );
                        // Finds place where the new value have to be inserted

   if( i < 0 )
    i = 0;              // New value have to be inserted before the first one in the array
   else
    i++;                // New value will be put somewhere in the middle of the array
                        //  or at the end

   if( i < fMVSetPoints )
   {
      Short_t *fromPtrF = fMV.GetArray() + i;
                       // Sets pointer to cell which have to be filled by new value
      Short_t *toPtrF = fromPtrF + 1;
                       // Sets pointer to cell where the array content have to be shifted 

      memmove( toPtrF, fromPtrF, (fMVSetPoints - i)*sizeof(Short_t) );
                       // Shifts array content. Now there is vacant place for new value to be inserted

      Int_t *fromPtrI = fMVTimeStamp.GetArray() + i;
      Int_t *toPtrI = fromPtrI + 1;
      memmove( toPtrI, fromPtrI, (fMVSetPoints - i)*sizeof(Int_t) );
                       // Do the same for time stamp array
   } /*if*/

   UShort_t val = (UShort_t)( voltage * fgkMVPrec );
                       // Float value of temperature is stored as UShort_t with given precision

   fMV.AddAt( (Short_t)val, i );
   fMVTimeStamp.AddAt( time, i );
   fMVSetPoints++;
                       // New values are inserted
} /*AliITSDCSDataSDD::SetValueMV*/

//---------------------------------------------------------------------------

void AliITSDCSDataSDD::SetValueStatus(Int_t time, Char_t status )
{
  // insert a value for channel status
  //

   if( fStatusMaxPoints == fStatusSetPoints )
    SetNPointsStatus( fStatusMaxPoints + AUTORESIZE );
                        // Enlarges arrays if necessary

   Int_t i =  FindIndex( time, fStatusTimeStamp, fStatusSetPoints );
                        // Finds place where the new value have to be inserted

   if( i < 0 )
    i = 0;              // New value have to be inserted before the first one in the array
   else
    i++;                // New value will be put somewhere in the middle of the array
                        //  or at the end

   if( i < fStatusSetPoints )
   {
      Char_t *fromPtrF = fStatus.GetArray() + i;
                       // Sets pointer to cell which have to be filled by new value
      Char_t *toPtrF = fromPtrF + 1;
                       // Sets pointer to cell where the array content have to be shifted 

      memmove( toPtrF, fromPtrF, (fStatusSetPoints - i)*sizeof(Char_t) );
                       // Shifts array content. Now there is vacant place for new value to be inserted

      Int_t *fromPtrI = fStatusTimeStamp.GetArray() + i;
      Int_t *toPtrI = fromPtrI + 1;
      memmove( toPtrI, fromPtrI, (fStatusSetPoints - i)*sizeof(Int_t) );
                       // Do the same for time stamp array
   } /*if*/

   fStatus.AddAt( status, i );
   fStatusTimeStamp.AddAt( time, i );
   fStatusSetPoints++;
                       // New values are inserted

} /*AliITSDCSDataSDD::SetValueStatus*/


//---------------------------------------------------------------------------

void AliITSDCSDataSDD::Compress()
{
// Minimize array sizes

   SetNPointsTempLeft( fTempLeftSetPoints );
   SetNPointsTempRight( fTempRightSetPoints );
   SetNPointsHV( fHVSetPoints );
   SetNPointsMV( fMVSetPoints );
   SetNPointsStatus( fStatusSetPoints );

} /*AliITSDCSDataSDD::Compress*/

//---------------------------------------------------------------------------

Float_t AliITSDCSDataSDD::GetDriftField( Int_t timeStamp ) const
{
// Returns drift field counted for specific time

   Int_t   cathodesNumber = 291;
   Float_t cathodesPitch = 0.0120;

   Float_t hv = GetHV( timeStamp );
   Float_t mv = GetMV( timeStamp );

   if( hv < 0.0 || mv < 0.0 ) return -1.0;
                        // HV or MV is unknown at this time

   return ( hv - mv ) / ( cathodesNumber * cathodesPitch );

} /*AliITSDCSDataSDD::GetDriftField*/

//---------------------------------------------------------------------------


Float_t AliITSDCSDataSDD::GetDriftSpeed( Int_t /*timeStamp*/ ) const
{
// Returns drift speed counted for specific time. Calculation is based on temerature
//  taken  from DCS. This metod is not dedicated for normal usage, it should be used
//  only in cases that the injectors for given module fails.
//
// Presently only a prototype, returns -1.0.

   /* PROTOTYPE */

   return -1.0;

} /*AliITSDCSDataSDD::*/

//---------------------------------------------------------------------------


void AliITSDCSDataSDD:: PrintValues( FILE *output ) const
{
// Prints array contents

    Int_t nTLEntries = GetTempLeftRecords();
    Int_t nTREntries = GetTempRightRecords();
    Int_t nHVEntries = GetHVRecords() ;
    Int_t nMVEntries = GetMVRecords();
    Int_t nStatEntries = GetStatusRecords();

    fprintf( output, "+------------------------------------------------------------------------------------------------------------+\n");
    fprintf( output, "|                                                DCS content                                                 |\n" ); 
    fprintf( output, "+----------------------+-----------------------+---------------------+---------------------+-----------------+\n");
    fprintf( output, "|    %05i  records    |    %05i   records    |    %05i  records   |    %05i  records   |  %05i records  |\n",
                          nHVEntries, nMVEntries, nTLEntries, nTREntries, nStatEntries );
    fprintf( output, "|  time (s)     HV     |  time (s)      MV     |  time (s)     TL    |  time (s)     TR    | time (s)   Stat |\n" );
    fprintf( output, "+----------------------+-----------------------+---------------------+---------------------+-----------------+\n");

    Int_t a = (nHVEntries > nMVEntries ) ? nHVEntries : nMVEntries;
    Int_t b = (nTLEntries > nTREntries ) ? nTLEntries : nTREntries;
    if( a < b ) a = b;
    Int_t loopMax = ( a > nStatEntries ) ? a : nStatEntries ;
                 // Finds maximal entry number

    for( Int_t entryLoop = 0; entryLoop < loopMax; entryLoop++ )
    {

        if( entryLoop < nHVEntries )
         fprintf( output, "| %12i %4.2f | ", GetHVTimeIdx(entryLoop), GetHVIdx(entryLoop) );
        else
         fprintf( output, "|                      | ");

        if( entryLoop < nMVEntries )
         fprintf( output, " %12i  %2.3f | ", GetMVTimeIdx(entryLoop), GetMVIdx(entryLoop) );
        else
         fprintf( output, "                      | ");

        if( entryLoop < nTLEntries )
         fprintf( output, "%12i  %2.2f | ", GetTempLeftTimeIdx(entryLoop), GetTempLeftIdx(entryLoop) );
        else
         fprintf( output, "                    | ");

        if( entryLoop < nTREntries )
         fprintf( output, "%12i  %2.2f | ", GetTempRightTimeIdx(entryLoop), GetTempRightIdx(entryLoop) );
        else
         fprintf( output, "                    | ");

        if( entryLoop < nStatEntries )
         fprintf( output, "%12i  %i |\n", GetStatusTimeIdx(entryLoop), GetStatusIdx(entryLoop) );
        else
         fprintf( output, "                |\n");

    } /*for( entryLoop )*/


} /*AliITSDCSDataSDD::PrintValues()*/

//---------------------------------------------------------------------------

Int_t AliITSDCSDataSDD::FindIndex( Int_t timeStamp, const TArrayI &timeStampArray, Int_t n ) const
{
// Provides binary search in the time array. Returns index in the array of time 
//  stamps by selected value. Returns -1 if the time is less than time stamp in 
//  the timeArray[0]

  if( n < 1 ) return -1;// Empty array or wrong value of array size

  if( timeStamp >= timeStampArray.At(n-1) ) return n-1;
                        // Time is larger than last timestamp - last value in the array have
                        //  to be used. This is the most frequent case, so it have sense
                        //  to check it and avoid searching.

  if( timeStamp < timeStampArray.At(0) ) return -1;
                        // Time is less than all time stamp stored in the array

  Int_t left = 0;
  Int_t right = n-1;
  Int_t middle = (left + right)/2;

  while( !( middle == left || middle == right) )
  {                     // Binary search in the time stamp array

     if( timeStampArray.At(middle) < timeStamp )
      left = middle;
     else
      right = middle;
     middle = (left + right)/2;
  } /*while*/

  if( timeStamp >= timeStampArray.At(right) )
   return right;
  else
   return left;

} /*AliITSDCSDataSDD::FindIndexByTimeStamp*/

//---------------------------------------------------------------------------
