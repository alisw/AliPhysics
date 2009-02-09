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
// Implementation of the class for SDD DCS data analysis         //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//         V.Pospisil, CTU Prague, gdermog@seznam.cz             //
///////////////////////////////////////////////////////////////////


#include "AliITSDCSAnalyzerSDD.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliITSgeomTGeo.h"

ClassImp(AliITSDCSAnalyzerSDD)

//---------------------------------------------------------------
  AliITSDCSAnalyzerSDD::AliITSDCSAnalyzerSDD(): TObject(),
fHVDelay(0),fMVDelay(0),fTLDelay(0),fTRDelay(0),fStTLDelay(0),fStTRDelay(0),fOKDelay(0),
fHVThresholdFrac(0), fMVThresholdFrac(0), fTLThresholdFrac(0), fTRThresholdFrac(0)
{
// Default constructor
  Init();
  SetHVThreshold();
  SetMVThreshold();
  SetTLThreshold();
  SetTRThreshold();
  for( Int_t moduleLoop = 0; moduleLoop < kNmodules; moduleLoop++ ) fDCSData[moduleLoop] = NULL;
} /*AliITSDCSAnalyzerSDD::AliITSDCSAnalyzerSDD*/

//---------------------------------------------------------------

AliITSDCSAnalyzerSDD::AliITSDCSAnalyzerSDD(const AliITSDCSAnalyzerSDD& /* dcsa */): TObject(),
fHVDelay(0),fMVDelay(0),fTLDelay(0),fTRDelay(0),fStTLDelay(0),fStTRDelay(0),fOKDelay(0),
fHVThresholdFrac(0), fMVThresholdFrac(0), fTLThresholdFrac(0), fTRThresholdFrac(0)
{
// Copy constructor
// Copies are not allowed. The method is protected to avoid misuse.
  AliError("Copy constructor not allowed");
} /*AliITSDCSAnalyzerSDD::AliITSDCSAnalyzerSDD*/

//---------------------------------------------------------------

AliITSDCSAnalyzerSDD& AliITSDCSAnalyzerSDD::operator=(const AliITSDCSAnalyzerSDD& /* dcsa */)
{
// Assigment operator
// Assignment is not allowed. The method is protected to avoid misuse.
  AliError("Assignment operator not allowed");
  return *this;
}/*AliITSDCSAnalyzerSDD::operator=*/

//---------------------------------------------------------------

AliITSDCSAnalyzerSDD::~AliITSDCSAnalyzerSDD()
{ 
// Destructor
  for(int j=0; j<kNmodules; j++)
  {
    if( fDCSData[j] ) delete fDCSData[j];
  } /*for( j )*/
} /*AliITSDCSAnalyzerSDD::~AliITSDCSAnalyzerSDD*/

//---------------------------------------------------------------

void AliITSDCSAnalyzerSDD::AnalyzeData(TMap* dcsMap)
{
// Data processing. Takes DCS points from alias map and sorts them into AliITSDCSDataSDD objects.

   Int_t   counter = 0; // Counter of stored DCS records

   Float_t lastTLValUpper;
   Float_t lastTLValLower;
   Float_t lastTRValUpper;
   Float_t lastTRValLower;
   Float_t lastHVValUpper;
   Float_t lastHVValLower;
   Float_t lastMVValUpper;
   Float_t lastMVValLower;
                        // Thresholds for float DCS variables

   Int_t nEntries;      // Number of entries in each TObjArray, that contains DCS variable values
   AliDCSValue *valToProcess;
                        // Pointer to currently processed DCS variable value
   Float_t valToProcessFloat;
                        // Value of currently processed DCS variable


   for( Int_t iLay = 3; iLay < 5; iLay++ )
   {

      Int_t maxLad = ( iLay == 3) ? kNladders3 : kNladders4;
      Int_t maxMod = ( iLay == 3) ? kNmodLad3 : kNmodLad4;

      for(Int_t iLad = 0; iLad < maxLad; iLad++)
      {
         for(Int_t iMod = 0; iMod < maxMod; iMod++)
         {
                        // Loads arrays of DCS variables from map. Variables are 
                        //  searched by names (for ex. SDD_LAYER3_LADDER5_MODULE4_HV)

            Int_t moduleLoop = AliITSgeomTGeo::GetModuleIndex( iLay, iLad + 1, iMod + 1 ) - 240;

            fDCSData[moduleLoop] = new AliITSDCSDataSDD();
                        // DCS data for specific SDD module will be stored in this class

            TObjArray* arrHV = (TObjArray*) dcsMap->GetValue( fHVDPNames[moduleLoop].Data() );
            if(!arrHV) AliWarning( Form("DCS HV alias %s not found!\n", fHVDPNames[moduleLoop].Data()) );

            TObjArray* arrMV = (TObjArray*) dcsMap->GetValue( fMVDPNames[moduleLoop].Data() );
            if(!arrMV) AliWarning( Form("DCS MV alias %s not found!\n", fMVDPNames[moduleLoop].Data()));

            TObjArray* arrOK = (TObjArray*) dcsMap->GetValue( fOKDPNames[moduleLoop].Data() );
            if(!arrOK)  AliWarning( Form("DCS MOD_OK alias %s not found!\n", fOKDPNames[moduleLoop].Data()));

            TObjArray* arrTL = (TObjArray*) dcsMap->GetValue( fTLDPNames[moduleLoop].Data() );
            if(!arrTL) AliWarning( Form("DCS TEMP_L alias %s not found!\n", fTLDPNames[moduleLoop].Data()));

            TObjArray* arrTR = (TObjArray*) dcsMap->GetValue( fTRDPNames[moduleLoop].Data() );
            if(!arrTR) AliWarning( Form("DCS TEMP_R alias %s not found!\n", fTRDPNames[moduleLoop].Data()));

            TObjArray* arrStTL = (TObjArray*) dcsMap->GetValue( fTLStDPNames[moduleLoop].Data() );
            if(!arrStTL) AliWarning( Form("DCS TEMP_L_STATE alias %s not found!\n", fTLStDPNames[moduleLoop].Data()));

            TObjArray* arrStTR = (TObjArray*) dcsMap->GetValue( fTRStDPNames[moduleLoop].Data() );
            if(!arrStTR) AliWarning( Form("DCS TEMP_R_STATE alias %s not found!\n", fTRStDPNames[moduleLoop].Data()));

            lastTLValUpper = -1e-10;
            lastTLValLower = +1e+10;
            lastTRValUpper = -1e-10;
            lastTRValLower = +1e+10;
            lastHVValUpper = -1e-10;
            lastHVValLower = +1e+10;
            lastMVValUpper = -1e-10;
            lastMVValLower = +1e+10;
                        // First value of any DCS variable must be written


            if( arrTL )
            {
               nEntries = arrTL->GetEntries();
               fDCSData[moduleLoop]->SetNPointsTempLeft( nEntries );
                        // Left temperature array size is set

               for( Int_t tlLoop = 0; tlLoop < nEntries; tlLoop++ )
               {        // Left temerature values are copied into the AliITSDCSDataSDD TempLeft array
                  valToProcess = (AliDCSValue *)(arrTL->At(tlLoop));
                  valToProcessFloat = valToProcess->GetFloat();
                        // Value is readed from the input array

                  if( lastTLValLower <= valToProcessFloat && valToProcessFloat <= lastTLValUpper ) continue;
                        // Value did not cross the treshold (upper neither lower),
                        //  it is not necessary to store it.
                  fDCSData[moduleLoop]->SetValueTempLeft( valToProcess->GetTimeStamp() - fTLDelay, valToProcessFloat );
                        // Value is stored
                  lastTLValLower = valToProcessFloat * ( 1.0 - fTLThresholdFrac );
                  lastTLValUpper = valToProcessFloat * ( 1.0 + fTLThresholdFrac );
                        // New tresholds are set
                  counter ++;
               } /*for( tlLoop )*/
            } /*if*/


            if( arrTR )
            {
               nEntries = arrTR->GetEntries();
               fDCSData[moduleLoop]->SetNPointsTempRight( nEntries );
                        // Right temperature array size is set 

               for( Int_t trLoop = 0; trLoop < nEntries; trLoop++ )
               {           // Right temerature values are copied into the AliITSDCSDataSDD TempRight array
                  valToProcess = (AliDCSValue *)(arrTR->At(trLoop));
                  valToProcessFloat = valToProcess->GetFloat();
                         // Value is readed from the input array

                  if( lastTRValLower <= valToProcessFloat && valToProcessFloat <= lastTRValUpper ) continue;
                         // Value did not cross the treshold (upper neither lower),
                         //  it is not necessary to store it.
                  fDCSData[moduleLoop]->SetValueTempRight( valToProcess->GetTimeStamp() - fTRDelay, valToProcessFloat );
                         // Value is stored
                  lastTRValLower = valToProcessFloat * ( 1.0 - fTRThresholdFrac );
                  lastTRValUpper = valToProcessFloat * ( 1.0 + fTRThresholdFrac );
                         // New tresholds are set
                  counter ++;
               } /*for( trLoop )*/
            } /*if*/


            if( arrHV )
            {
               nEntries = arrHV->GetEntries();
               fDCSData[moduleLoop]->SetNPointsHV( nEntries );
                        // HV array size is set 

               for( Int_t hvLoop = 0; hvLoop < nEntries; hvLoop++ )
               {        // HV values are copied into the AliITSDCSDataSDD HV array
                  valToProcess = (AliDCSValue *)(arrHV->At(hvLoop));
                  valToProcessFloat = valToProcess->GetFloat();
                        // Value is readed from the input array
                  if( lastHVValLower <= valToProcessFloat && valToProcessFloat <= lastHVValUpper ) continue;
                        // Value did not cross the treshold (upper neither lower),
                        //  it is not necessary to store it.
                  fDCSData[moduleLoop]->SetValueHV( valToProcess->GetTimeStamp() - fHVDelay, valToProcessFloat );
                        // Value is stored
                  lastHVValLower = valToProcessFloat * ( 1.0 - fHVThresholdFrac );
                  lastHVValUpper = valToProcessFloat * ( 1.0 + fHVThresholdFrac );
                        // New tresholds are set
                  counter ++;
               } /*for( hvLoop )*/

            } /*if*/



            if( arrMV )
            {
               nEntries = arrMV->GetEntries();
               fDCSData[moduleLoop]->SetNPointsMV( nEntries );
                        // MV array size is set 

               for( Int_t mvLoop = 0; mvLoop < nEntries; mvLoop++ )
               {        // MV values are copied into the AliITSDCSDataSDD MV array
                  valToProcess = (AliDCSValue *)(arrMV->At(mvLoop));
                  valToProcessFloat = valToProcess->GetFloat();
                        // Value is readed from the input array
                  if( lastMVValLower <= valToProcessFloat && valToProcessFloat <= lastMVValUpper ) continue;
                        // Value did not cross the treshold (upper neither lower),
                        //  it is not necessary to store it.
                  fDCSData[moduleLoop]->SetValueMV( valToProcess->GetTimeStamp() - fMVDelay, valToProcessFloat );
                        // Value is stored
                  lastMVValLower = valToProcessFloat * ( 1.0 - fMVThresholdFrac );
                  lastMVValUpper = valToProcessFloat * ( 1.0 + fMVThresholdFrac );
                        // New treshold is ser
                  counter ++;
               } /*for( mvLoop )*/

            } /*if*/


/* Following part of the code is responsibile for the condensing of all status information given by DCS
   into one array of Char_t. Each record of this array is in principle a bit map : 

      0. bit ... _OK
      1. bit ... _TEMP_L_STATE
      2. bit ... _TEMP_R_STATE 

   Each record have its own time stamp. Because there are three inputs with independent time stamp,
   some algorithm which assigns new time stamp to bitmap according to three input time stamps is
   necessary.

   Let's vizualize time stamps of the three input arrays. There is time on x-axis :

            +------------+---------------------+------
            |            |                     |          _OK
      +-----++------+----+--------+------------+------
      |      |      |             |            |          _TEMP_L_STATE
   +--+------+---+--+-------+-----+--------+---+------
   |             |          |     |        |   |          _TEMP_R_STATE
   +-------------+----------+-----+--------+---+------

   |   |    | |  |   |   |  |     |        |   |
   V   V    V V  V   V   V  V     V        V   V

   +---+----+-+--+---+---+--+-----+--------+---+------
   |   |    | |  |   |   |  |     |        |   |           Status bitmap
   +---+----+-+--+---+---+--+-----+--------+---+------


   Principle of combining three status records into one is visible from the picture.
   If there are two sequent records with the same status bitmap, they are joined into
   one (with the time stamp of the earliest one).

*/
            Int_t nStTLEntries = 0;
            Int_t nStTREntries = 0;
            Int_t nOKEntries = 0;

            bool arrStTLcreated = false;
            bool arrStTRcreated = false;
            bool arrOKcreated = false;

            if( arrStTL ) 
             nStTLEntries = arrStTL->GetEntries();
            else
             { arrStTL = new TObjArray; arrStTLcreated = true; }

            if( arrStTR ) 
             nStTREntries = arrStTR->GetEntries();
            else
             { arrStTR = new TObjArray; arrStTRcreated = true; }

            if( arrOK ) 
             nOKEntries = arrOK->GetEntries();
            else
             { arrOK = new TObjArray; arrOKcreated = true; }
                        // Gets number of _STAT_L, _STAT_R and _OK values stored in dcsMap. If any array does
                        //  not exist, it must be created (and it will be filled by 0 status later)

            if( nStTLEntries < 1 )
            {           // TObjArray arrStTL is empty. This would cause segmentation violation during
                        //  the condensing, so this case must be handled before algorithm starts
               AliWarning( Form( "%s contains no data!\n", fTLStDPNames[moduleLoop].Data() ) );
               nStTLEntries = 1;
               arrStTL->Add( new AliDCSValue(  (Int_t)0, 0x7FFFFFFF ) );
                        // 0x7FFFFFFF = 2147483647, maximal signed Int_t number. Left temperature
                        //  sensor will be regarded as switched-off during whole run.
            } /*if*/
 
            if( nStTREntries < 1 )
            {           // TObjArray arrStTR is empty. This would cause segmentation violation during
                        //  the condensing, so this case must be handled before algorithm starts
               AliWarning( Form( "%s contains no data!\n", fTRStDPNames[moduleLoop].Data() ) );
               nStTREntries = 1;
               arrStTR->Add( new AliDCSValue(  (Int_t)0, 0x7FFFFFFF ) );
                        // 0x7FFFFFFF = 2147483647, maximal signed Int_t number. Right temperature
                        //  sensor will be regarded as switched-off during whole run.
            } /*if*/

            if( nOKEntries < 1 )
            {           // TObjArray arrOK is empty. This would cause segmentation violation during
                        //  the condensing, so this case must be handled before algorithm starts
               AliWarning( Form( "%s contains no data!\n", fOKDPNames[moduleLoop].Data() ) );
               nOKEntries = 1;
               arrOK->Add( new AliDCSValue(  (Bool_t)0, 0x7FFFFFFF ) );
                        // 0x7FFFFFFF = 2147483647, maximal signed Int_t number.
                        //  Module will be regarded as switched-off during whole run.
            } /*if*/

            arrStTL->Sort();
            arrStTR->Sort();
            arrOK->Sort();
                        // Condensing would work properly only in the case that
                        //  the input arrays are sorted by time

            Int_t nEntriesMax = nStTLEntries + nStTREntries + nOKEntries;
            fDCSData[moduleLoop]->SetNPointsStatus( nEntriesMax );
                        // Determines necessary length of new array and sets its size 
                        //  Lot of space in such defined array will be probably 
                        //  vacant after the analysis, but this will be corrected
                        //  by Compress() method

            Int_t idxStTL = 0;
            Int_t idxStTR = 0;
            Int_t idxOK = 0;
                        // Input arrays indexes

            Int_t tsStTL, tsStTR, tsOK;
                        // Time stamps ofinput arrays
            Int_t tsNew;// New time stamp (output array)

            Char_t bitStatus;   
                        // New status record :
                        // 0. bit ... _OK
                        // 1. bit ... _TEMP_L_STATE
                        // 2. bit ... _TEMP_R_STATE 
            Char_t lastBitStatus = 100;

            AliDCSValue *valStTL, *valStTR, *valOK;
                        // Pointers to input arrays records (input arrays are TObjArrays
                        //  containing objects of type AliDCSValue

            tsStTR = ( (AliDCSValue *)arrStTR->At(0) )->GetTimeStamp() - fStTLDelay;
            tsStTL = ( (AliDCSValue *)arrStTL->At(0) )->GetTimeStamp() - fStTRDelay;
            tsOK = ( (AliDCSValue *)arrOK->At(0) )->GetTimeStamp() - fOKDelay;
                        // Time stamps of first records in input filea are readed (and delays are substracted)

            tsNew = (tsStTR < tsStTL) ? tsStTR : tsStTL;
            if( tsNew > tsOK ) tsNew = tsOK;
                        // Time intervals are "prolonged" to the very eaarliest of time stamps.
                        //  It means that first output time stamp will be the same as the one
                        //  which is first in input arrays. Values of other DCS variables are
                        //  not defined in this time yet, but they will be treated as equal to
                        //  values in first records of input arrays.

            nStTLEntries--; nStTREntries--; nOKEntries--;
                        // Indexes in the input array must not exceed last records.

            while( (idxStTL < nStTLEntries) || (idxStTR < nStTREntries) || (idxOK < nOKEntries) )
            {           // Loop goes throug all three input files

               valStTL = (AliDCSValue *)( arrStTL->At(idxStTL) );
               valStTR = (AliDCSValue *)( arrStTR->At(idxStTR) );
               valOK = (AliDCSValue *)( arrOK->At(idxOK) );
                        // Values are readed from input arrays

               bitStatus = 0;
               if( valOK->GetBool() ) bitStatus += 1;        // 0. bit - _OK
               if( valStTL->GetInt() == 4 ) bitStatus += 2;  // 1. bit - _TEMP_L_STATE
               if( valStTR->GetInt() == 4 ) bitStatus += 4;  // 2. bit - _TEMP_R_STATE
                        // Bit map is created. *TEMP_*_STATE == 4 means "Thermometer is OK"

               if( lastBitStatus != bitStatus )
               {        // If the status bitmap is the same as last one, it would not be stored.
                        //  It will save much space.
                  fDCSData[moduleLoop]->SetValueStatus( tsNew, bitStatus );
                        // Bit map is written into the output array (if different from last value )
                  lastBitStatus = bitStatus;
                  counter += nEntries;
               } /*if*/

               if( idxStTL == nStTLEntries )
                tsStTL = 0x7FFFFFFF;  // = 2147483647, maximal signed Int_t number
               else
                tsStTL = ( (AliDCSValue *)arrStTL->At(idxStTL + 1) )->GetTimeStamp() - fStTLDelay;

               if( idxStTR == nStTREntries )
                tsStTR = 0x7FFFFFFF;  // = 2147483647, maximal signed Int_t number
               else
                tsStTR = ( (AliDCSValue *)arrStTR->At(idxStTR + 1) )->GetTimeStamp() - fStTRDelay;

               if( idxOK == nOKEntries )
                tsOK = 0x7FFFFFFF;    // = 2147483647, maximal signed Int_t number
               else
                tsOK = ( (AliDCSValue *)arrOK->At(idxOK + 1) )->GetTimeStamp() - fOKDelay;
                        // Reads time stamps of folowing records in the input arrays (and substracts delays).
                        // Validity of the last records in the input arrays are prolonged
                        //  to "infinity"

               if( tsStTL == tsOK && tsStTR == tsOK )  { tsNew = tsStTL; idxStTL++; idxStTR++; idxOK++; continue; }
               if( tsStTL == tsStTR && tsStTR < tsOK ) { tsNew = tsStTL; idxStTL++; idxStTR++; continue; }
               if( tsStTL == tsOK && tsOK < tsStTR )   { tsNew = tsStTL; idxStTL++; idxOK++; continue; }
               if( tsStTR == tsOK && tsOK < tsStTL )   { tsNew = tsStTR; idxStTR++; idxOK++; continue; } 
               if( tsOK < tsStTL && tsOK < tsStTR )    { tsNew = tsOK;   idxOK++;   continue; }
               if( tsStTL < tsOK && tsStTL < tsStTR )  { tsNew = tsStTL; idxStTL++; continue; }
               /*Last possibile case*/                 { tsNew = tsStTR; idxStTR++; }

                        // Index of array, whose following record have time stamp closest to just written one,
                        //  is increased. If there are more records with identical time stamps meeting this condition,
                        //  all correspondent indexes are increased. 

            } /*while*/

            fDCSData[moduleLoop]->Compress();
                        // Size taken by data in AliITSDCSDataSDD object is minimalized

            if( arrStTRcreated ) delete arrStTR;
            if( arrStTLcreated ) delete arrStTL;
            if( arrOKcreated ) delete arrOK;

          } /*for( iMod )*/
       } /*for( iLad )*/
    } /*for( iLay )*/


} /*AliITSDCSAnalyzerSDD::AnalyzeData*/


//---------------------------------------------------------------


void AliITSDCSAnalyzerSDD::Init()
{
// Initialization of DCS DP names
  Char_t dpName[50];
  Char_t modName[50];

  for( Int_t iLay = 3; iLay < 5; iLay++ )
  {    
     Int_t maxLad = ( iLay == 3) ? kNladders3 : kNladders4;
     Int_t maxMod = ( iLay == 3) ? kNmodLad3 : kNmodLad4;

     for(Int_t iLad=0; iLad<maxLad; iLad++)
     {
        for(Int_t iMod=0; iMod<maxMod;iMod++)
        {
           sprintf(modName,"SDD_LAYER%i_LADDER%02d_MODULE%d", iLay, iLad, iMod);
           Int_t id = AliITSgeomTGeo::GetModuleIndex( iLay, iLad + 1, iMod + 1 ) - 240;

           sprintf(dpName,"%s_HV",modName);
           fHVDPNames[id]=dpName;
           sprintf(dpName,"%s_MV",modName);
           fMVDPNames[id]=dpName;
           sprintf(dpName,"%s_OK",modName);
           fOKDPNames[id]=dpName;
           sprintf(dpName,"%s_TEMP_L",modName);
           fTLDPNames[id]=dpName;
           sprintf(dpName,"%s_TEMP_R",modName);
           fTRDPNames[id]=dpName;
           sprintf(dpName,"%s_TEMP_L_STATE",modName);
           fTLStDPNames[id]=dpName;
           sprintf(dpName,"%s_TEMP_R_STATE",modName);
           fTRStDPNames[id]=dpName;

        } /*for( iMod )*/
     } /*for( iLad )*/

  } /*for( iLay )*/

  
} /*AliITSDCSAnalyzerSDD::Init*/

//---------------------------------------------------------------
void AliITSDCSAnalyzerSDD::PrintDCSDPNames( FILE *output )
{
// Prints constructed names of DCS variables into specified file (may be even stdout or stderr)
  for( Int_t j = 0; j < kNmodules; j++ )
  {
    fprintf( output, "Module %d      %s   %s   %s   %s\n",j,fHVDPNames[j].Data(),
                          fMVDPNames[j].Data(),fTLDPNames[j].Data(),fTRDPNames[j].Data());
  } /*for( j )*/
} /*AliITSDCSAnalyzerSDD::PrintDCSDPNames*/

//---------------------------------------------------------------

void AliITSDCSAnalyzerSDD::Export( char *outputDCSFileName )
{
// Exports all stored AliITSDCSDataSDD type object into specified root file. Objects are named as
//
// DCSDataSDD_module<number>
//
// where <number> is in range 0..256 and it is obtained by calling
//
// AliITSgeomTGeo::GetModuleIndex( layer, ladder, moduleInLadder ) - 240

   TFile * newFile = new TFile( outputDCSFileName, "RECREATE" );
   if( newFile == NULL )
   {                    // Creates .root file with specified name. if it is not possible,
                        //  warning is displayed and exporting aborted.
     AliWarning( Form( "Cannot create %s - export aborted ", outputDCSFileName ) );
     return;
   } /*if*/

   newFile->cd();

   char buffer[100];

   for( Int_t moduleLoop = 0; moduleLoop < kNmodules; moduleLoop++ )
   {                    // loops through all modules and writes appropriate object into the file
      sprintf( buffer, "DCSDataSDD_module%i", moduleLoop );
      if( fDCSData[moduleLoop] ) fDCSData[moduleLoop]->Write( buffer, TObject::kSingleKey );
   } /*for( moduleLoop )*/

   newFile->Close();
   delete newFile;

} /*AliITSDCSAnalyzerSDD::Export*/
