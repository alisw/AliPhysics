#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTimeStamp.h>
#include <TMap.h>
#include <TMath.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRandom.h>
#include <TSystem.h>
#include <AliDCSValue.h>
#endif

//////////////////////////////////////////////////////////////////////////////////////
// Generator of testing data for AliITSDCSPreprocessorSDD and  AliITSAnalyzerrSDD   //
// Origin: F.Prino, Torino, prino@to.infn.it                                        //
//         V.Pospisil, CTU Prague, gdermog@seznam.cz                                //
//////////////////////////////////////////////////////////////////////////////////////

TMap* CreateRealisticMap( TTimeStamp startTime, TTimeStamp stopTime );
TMap* CreateRandomMap( Int_t maxRecords = 1000, Int_t randomTimeStamps = 0 );
void GenerateOutputFiles( TMap *map, char *dir );

Int_t sizeOfMapContent = 0;

void CreateSDDDCSMap( void )
{

 TTimeStamp startT( 2008, 2, 7, 15, 0, 0 );
 TTimeStamp stopT( 2008, 2, 7, 18, 0, 0 );
                        // The beginning and the end of the simulated run

//  delete gRandom;
//  gRandom = new TRandom2;
 gRandom->SetSeed(0);

/* Macro can generate random purely random data or simulate DCS. Choose
   one of the methods - CreateRandomMap() or CreateRealisticMap() */

// TMap *outMap = CreateRandomMap( 1000 );
 TMap *outMap = CreateRealisticMap( startT, stopT );

 GenerateOutputFiles( outMap, "./maps_3h_SIM" );

 TFile newfile( "DCSAliasMap_3h_SIM.root", "RECREATE" );

 printf( "Writting map into the file ...\n");

 outMap->Write( "DCSAliasMap", TObject::kSingleKey );

 printf( "DCSAliasMap created, it size is %i byte ...\n", sizeOfMapContent );

} /*CreateMap*/


/***************************************************************/
/*     DCS simulation setup - used in CreateRealisticMap()     */
/*   Choose simulation setup by changing following variables   */
/***************************************************************/

// --- Temperature (left side) - behavior of SDD chip -----------

 Float_t mvalTL = 20.0; // Mean value of temperature
 Float_t fluctTL = 0.5; // Fluctuation (gaussian sigma)

// --- Temperature readout (left side) - behavior of DCS --------

 Int_t freqROTL = 300;  // Frequency of fixed readout
 Float_t fluctROTL = 5.0;
                        // Allowed fluctuation (%)

// Temperature readout status (left side) - behavior of SDD chip

  Int_t freqROStTL = 300;// Frequency of fixed readout
  Float_t failureTL = 0.01;
                        // Probability of thermometer
                        //  failure (in one hour)

// --- Temperature (right side) - behavior of SDD chip ----------

 Float_t mvalTR = 20.0; // Mean value of temperature
 Float_t fluctTR = 0.5;// Fluctuation (gaussian sigma)

// --- Temperature readout (right side) - behavior of DCS -------

 Int_t freqROTR = 300;   // Frequency of fixed readout
 Float_t fluctROTR = 5.0;
                        // Allowed fluctuation (%)

// Temperature readout status (right side) - behavior of SDD

  Int_t freqROStTR = 300;// Frequency of fixed readout
  Float_t failureTR = 0.01;
                        // Probability of thermometer
                        //  failure (in one hour) 

// --- High voltage  - behavior of SDD voltage source -----------

 Float_t mvalHV = 1791.0;
                        // Mean value of HV
 Float_t fluctHV = 0.03;// Fluctuation (gaussian sigma)

// --- High voltage readout - behavior of DCS -------------------

 Int_t freqROHV = 300;   // Frequency of fixed readout
 Float_t fluctROHV = 0.01;
                        // Allowed fluctuation (%)

// --- High voltage readout status - behavior of SDD voltage ----

  Int_t freqROOK = 300;  // Frequency of fixed readout
  Float_t failureHV  = 0.005;
                        // Probability of HV source
                        //  failure (in one hour) 

// --- Medium voltage  - behavior of SDD voltage source --------

 Float_t mvalMV = 45.0; // Mean value of MV
 Float_t fluctMV = 0.005;
                        // Fluctuation (gaussian sigma)

// --- Medium voltage readout - behavior of DCS -----------------

 Int_t freqROMV = 300;   // Frequency of fixed readout
 Float_t fluctROMV = 0.1;
                        // Allowed fluctuation (%)

// --- Medium voltage readout status - behavior of SDD voltage --

 Float_t failureMV  = 0.005;
                        // Probability of MV source
                        //  failure (in one hour)

/***************************************************************/
/*   ITS geometry setup                                        */
/***************************************************************/

 Int_t kNladders3=14;
 Int_t kNladders4=22;
 Int_t kNmodLad3=6;
 Int_t kNmodLad4=8;

/***************************************************************/

//----------------------------------------------------------------------------------------

TMap* CreateRealisticMap( TTimeStamp startTime, TTimeStamp stopTime )
{
  // Creates a DCS structure
  // The structure is the following:
  //   TMap (key --> value)
  //     <DCSAlias> --> <valueList>
  //     <DCSAlias> is a string
  //     <valueList> is a TObjArray of AliDCSValue
  //     An AliDCSValue consists of timestamp and a value in form of a AliSimpleValue

  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);
                        // Empty map is created

  sizeOfMapContent += sizeof( TMap );

  TString aliasName;
  Char_t dpName[50];
  Int_t created = 0;
  Int_t created1 = 0;   // DCS records counter

  Double_t failureTLpersec;
  Double_t failureTRpersec;
  Double_t failureHVpersec;
  Double_t failureMVpersec;
                        // Probabilities of failure per one second
  Double_t lambda;
  lambda = ( TMath::Log( 1.0 - failureTL ) ) / -3600.0;
  failureTLpersec = 1.0 - TMath::Exp(-lambda );
  lambda = ( TMath::Log( 1.0 - failureTR ) ) / -3600.0;
  failureTRpersec = 1.0 - TMath::Exp(-lambda );
  lambda = ( TMath::Log( 1.0 - failureHV ) ) / -3600.0;
  failureHVpersec = 1.0 - TMath::Exp(-lambda );
  lambda = ( TMath::Log( 1.0 - failureMV ) ) / -3600.0;
  failureMVpersec = 1.0 - TMath::Exp(-lambda );
                        // Probabilities of failure per one second are counted

  fluctROTL /= 100.0;
  fluctROTR /= 100.0;
  fluctROHV /= 100.0;
  fluctROMV /= 100.0;   // Percents to fractions

  Double_t failureOK = ( failureHV + failureMV ) / 2.0;
                        // This value is used only for noise in timing

  for( Int_t iLay = 3; iLay < 5; iLay++ )
  {

      Int_t maxLad = ( iLay == 3) ? kNladders3 : kNladders4;
      Int_t maxMod = ( iLay == 3) ? kNmodLad3 : kNmodLad4;

      for(Int_t iLad=0; iLad<maxLad; iLad++)
      {
         for(Int_t iMod=0; iMod<maxMod; iMod++)
         {

            fprintf( stderr, "Simulating data taking for layer %i, ladder %i, module %i ... \n", iLay, iLad, iMod );

            TObjArray* valueSetOK = new TObjArray;
            TObjArray* valueSetHV = new TObjArray;
            TObjArray* valueSetMV = new TObjArray;
            TObjArray* valueSetTL = new TObjArray;
            TObjArray* valueSetTR = new TObjArray;
            TObjArray* valueSetStTL = new TObjArray;
            TObjArray* valueSetStTR = new TObjArray;

            sizeOfMapContent += 7 * sizeof( TObjArray );

            valueSetOK->SetOwner(1);
            valueSetHV->SetOwner(1);
            valueSetMV->SetOwner(1);
            valueSetTL->SetOwner(1);
            valueSetTR->SetOwner(1);
            valueSetStTL->SetOwner(1);
            valueSetStTR->SetOwner(1);

            Float_t actualTL, actualTR;
            Float_t actualHV, actualMV;
            Int_t   actualStTL = 1;
            Int_t   actualStTR = 1;
            Int_t   actualStHV = 1;
            Int_t   actualStMV = 1;
                             // Readout devices are alive/dead

            Float_t lastTL = mvalTL;
            Float_t lastTR = mvalTR;
            Float_t lastHV = mvalHV;
            Float_t lastMV = mvalMV;
                             // Last written values - udes for threshold readout

            Int_t counterTL = 0;
            Int_t counterTR = 0;
            Int_t counterHV = 0;
            Int_t counterMV = 0;
            Int_t counterStTL = 0;
            Int_t counterStTR = 0;
            Int_t counterOK = 0;
                              // Periodic readout counters

            Int_t endingTimeStamp = /*1197477000;*/ stopTime.GetSec();

            for( Int_t timeLoop = /*1197470000 */ startTime.GetSec(); timeLoop < endingTimeStamp; timeLoop ++ )
            {                 // Loop goes through period of run second per second and determines
                              //  all values according to rules of DCS

               actualTL = gRandom->Gaus( mvalTL, fluctTL );
               actualTR = gRandom->Gaus( mvalTR, fluctTR );
               actualHV = gRandom->Gaus( mvalHV, fluctHV );
               actualMV = gRandom->Gaus( mvalMV, fluctMV );
                              // Generates random values of temperatures and voltages
                              //  in given ranges

               if( gRandom->Rndm() < failureTLpersec ) 
                if( actualStTL ) actualStTL = 0; else if( gRandom->Rndm() < 0.1 ) actualStTL = 1;
               if( gRandom->Rndm() < failureTRpersec ) 
                if( actualStTR ) actualStTR = 0; else if( gRandom->Rndm() < 0.1 ) actualStTR = 1;
               if( gRandom->Rndm() < failureHVpersec ) 
                if( actualStHV ) actualStHV = 0; else if( gRandom->Rndm() < 0.1 ) actualStHV = 1;
               if( gRandom->Rndm() < failureMVpersec ) 
                if( actualStMV ) actualStMV = 0; else if( gRandom->Rndm() < 0.1 ) actualStMV = 1;
                              // Decides if any thermometer or voltage source becomes 
                              //  dead (or alive again, but this have much 10x les probability )

               if( gRandom->Rndm() < failureTL ) counterTL++;
               if( gRandom->Rndm() < failureTL ) counterTL--;
               if( gRandom->Rndm() < failureTR ) counterTR++;
               if( gRandom->Rndm() < failureTR ) counterTR--;
               if( gRandom->Rndm() < failureHV ) counterHV++;
               if( gRandom->Rndm() < failureHV ) counterHV--;
               if( gRandom->Rndm() < failureMV ) counterMV++;
               if( gRandom->Rndm() < failureMV ) counterMV--;
               if( gRandom->Rndm() < failureTL ) counterStTL++;
               if( gRandom->Rndm() < failureTL ) counterStTL--;
               if( gRandom->Rndm() < failureTR ) counterStTR++;
               if( gRandom->Rndm() < failureTR ) counterStTR--;
               if( gRandom->Rndm() < failureOK ) counterOK++;
               if( gRandom->Rndm() < failureOK ) counterOK--;
                        // Simulating noise in the clock frequency

               if( counterTL >= freqROTL )
                { valueSetTL->Add( new AliDCSValue( actualTL, timeLoop ) ); lastTL = actualTL; counterTL = 0; created++; }
               if( counterTR >= freqROTR )
                { valueSetTR->Add( new AliDCSValue( actualTR, timeLoop ) ); lastTR = actualTR; counterTR = 0; created++; }
               if( counterHV >= freqROHV )
                { valueSetHV->Add( new AliDCSValue( actualHV, timeLoop ) ); lastHV = actualHV; counterHV = 0; created++; }
               if( counterMV >= freqROMV )
                { valueSetMV->Add( new AliDCSValue( actualMV, timeLoop ) ); lastMV = actualMV; counterMV = 0; created++; }
               if( counterStTL >= freqROStTL )
                { valueSetStTL->Add( new AliDCSValue( actualStTL, timeLoop ) ); counterStTL = 0; created++; }
               if( counterStTR >= freqROStTR )
                { valueSetStTR->Add( new AliDCSValue( actualStTR, timeLoop ) ); counterStTR = 0; created++; }
               if( counterOK >= freqROOK )
                { valueSetOK->Add( new AliDCSValue( (Bool_t)(actualStHV & actualStMV), timeLoop ) ); counterOK = 0; created++; }
                        // Periodic readout

               if( TMath::Abs( (lastTL - actualTL) / lastTL ) > fluctROTL )
                { valueSetTL->Add( new AliDCSValue( actualTL, timeLoop ) ); lastTL = actualTL; counterTL = 0; created1++; }
               if( TMath::Abs( (lastTR - actualTR) / lastTR ) > fluctROTR )
                { valueSetTR->Add( new AliDCSValue( actualTR, timeLoop ) ); lastTR = actualTR; counterTR = 0; created1++; }
               if( TMath::Abs( (lastHV - actualHV) / lastHV ) > fluctROHV )
                { valueSetHV->Add( new AliDCSValue( actualHV, timeLoop ) ); lastHV = actualHV; counterHV = 0; created1++; }
               if( TMath::Abs( (lastMV - actualMV) / lastMV ) > fluctROMV )
                { valueSetMV->Add( new AliDCSValue( actualMV, timeLoop ) ); lastMV = actualMV; counterMV = 0; created1++; }
                        // Treshold readout

               counterTL++;
               counterTR++;
               counterHV++;
               counterMV++;
               counterStTL++;
               counterStTR++;
               counterOK++;

            } /*for( timeLoop )*/

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_OK", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetOK );

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_HV", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetHV );

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_MV", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetMV );

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_L", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetTL );

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_R", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetTR );

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_L_STATE", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetStTL );

            sprintf( dpName, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_R_STATE", iLay, iLad, iMod );
            aliasName = dpName;
            aliasMap->Add( new TObjString(aliasName), valueSetStTR );

         } /*for( iMod ) */
      } /*for( iLad )*/
  } /*for( iLay )*/

  fprintf(  stderr, "\nCreated %i objects of type AliDCSValue (%i periodic + %i treshold)... \n", 
            created + created1, created, created1 );

  sizeOfMapContent += (created + created1 ) * ( sizeof( AliDCSValue ) + sizeof( AliDCSValue * ) );

  return aliasMap;

} /*CreateRealisticMap*/

//----------------------------------------------------------------------------------------


TMap* CreateRandomMap( Int_t maxRecords , Int_t randomTimeStamps  )
{
  // Creates a DCS structure
  // The structure is the following:
  //   TMap (key --> value)
  //     <DCSAlias> --> <valueList>
  //     <DCSAlias> is a string
  //     <valueList> is a TObjArray of AliDCSValue
  //     An AliDCSValue consists of timestamp and a value in form of a AliSimpleValue


  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);

  sizeOfMapContent += sizeof( TMap );

  TString aliasName;
  Char_t dpName[50];
  Int_t created = 0;

  for( Int_t iLay = 3; iLay < 5; iLay++ )
  {

      Int_t maxLad = ( iLay == 3) ? kNladders3 : kNladders4;
      Int_t maxMod = ( iLay == 3) ? kNmodLad3 : kNmodLad4;

      for(Int_t iLad = 0; iLad < maxLad; iLad++)
      {
         for(Int_t iMod =0 ; iMod < maxMod;iMod++)
         {

            fprintf( stderr, "Generating data for layer %i, ladder %i, module %i ... \n", iLay, iLad, iMod );

            TObjArray* valueSetOK = new TObjArray;
            TObjArray* valueSetH = new TObjArray;
            TObjArray* valueSetM = new TObjArray;
            TObjArray* valueSetTL = new TObjArray;
            TObjArray* valueSetTR = new TObjArray;
            TObjArray* valueSetStTL = new TObjArray;
            TObjArray* valueSetStTR = new TObjArray;

            sizeOfMapContent += 7 * sizeof( TObjArray );

            valueSetOK->SetOwner(1);
            valueSetH->SetOwner(1);
            valueSetM->SetOwner(1);
            valueSetTL->SetOwner(1);
            valueSetTR->SetOwner(1);
            valueSetStTL->SetOwner(1);
            valueSetStTR->SetOwner(1);

            Int_t nrOfRecords;
            AliDCSValue* dcsVal;
            Int_t timeStamp;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_OK", iLay,iLad,iMod);
            aliasName=dpName;
	    nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ) );
            timeStamp = 1000000000;
	    for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
                timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
               else
                timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               if( 1000*(gRandom->Rndm()) > 50 ) 
                dcsVal = new AliDCSValue((Bool_t)1, timeStamp);
               else
                dcsVal = new AliDCSValue((Bool_t)0, timeStamp);
               valueSetOK->Add(dcsVal);
            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetOK);
            created += nrOfRecords;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_HV", iLay, iLad,iMod);
            aliasName=dpName;
            nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ));
            timeStamp = 1000000000;
            for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
		timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
               else
                timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               dcsVal = new AliDCSValue( (Float_t)( 1600 + 200*(gRandom->Rndm()) ), timeStamp );
               valueSetH->Add(dcsVal);
            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetH);
            created += nrOfRecords;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_MV", iLay, iLad,iMod);
            aliasName=dpName;
	    nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ));
            timeStamp = 1000000000;
            for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
		timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
               else
                timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               dcsVal = new AliDCSValue( (Float_t)( 30 + 20*(gRandom->Rndm()) ), timeStamp );
               valueSetM->Add(dcsVal);
            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetM);
            created += nrOfRecords;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_L", iLay, iLad,iMod);
            aliasName=dpName;
	    nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ));
            timeStamp = 1000000000;
            for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
		timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
                else
		  timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               dcsVal = new AliDCSValue( (Float_t)( 50 + 50*(gRandom->Rndm()) ), timeStamp );
               valueSetTL->Add(dcsVal);
            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetTL);
            created += nrOfRecords;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_R", iLay, iLad,iMod);
            aliasName=dpName;
	    nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ));
	    timeStamp = 1000000000;
            for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
		timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
               else
		 timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               dcsVal = new AliDCSValue( (Float_t)( 50 + 50*(gRandom->Rndm()) ), timeStamp );
               valueSetTR->Add(dcsVal);
            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetTR);
            created += nrOfRecords;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_L_STATE", iLay, iLad,iMod);
            aliasName=dpName;
	    nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ));
            timeStamp = 1000000000;
            for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
		 timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
               else
		 timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               if( 1000*(gRandom->Rndm()) > 50 ) 
                  dcsVal = new AliDCSValue((Int_t)1, timeStamp);
               else
                  dcsVal = new AliDCSValue((Int_t)0, timeStamp);
               valueSetStTL->Add(dcsVal);
            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetStTL);
            created += nrOfRecords;

            sprintf(dpName,"SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_R_STATE", iLay, iLad,iMod);
            aliasName=dpName;
	    nrOfRecords = (Int_t)(maxRecords * ( gRandom->Rndm() ));
            timeStamp = 1000000000;
            for( Int_t recLoop = 0; recLoop < nrOfRecords; recLoop ++ )
            {
               if( randomTimeStamps )
		timeStamp = (Int_t)(1200000000*(gRandom->Rndm()));
               else
		 timeStamp += (Int_t)(50*(gRandom->Rndm()) );
               if( 1000*(gRandom->Rndm()) > 50 ) 
                  dcsVal = new AliDCSValue((Int_t)1, timeStamp);
               else
                  dcsVal = new AliDCSValue((Int_t)0, timeStamp);
               valueSetStTR->Add(dcsVal);

            } /*for( recLoop )*/
            aliasMap->Add(new TObjString(aliasName), valueSetStTR);
            created += nrOfRecords;

         } /*for( iMod ) */
      } /*for( iLad )*/

  } /*for( iLay )*/

  fprintf(  stderr, "\nCreated %i objects of type AliDCSValue ... \n", created);

  sizeOfMapContent += created * ( sizeof( AliDCSValue ) + sizeof( AliDCSValue * ) );

  return aliasMap;
}

// -----------------------------------------------------------------------------

void GenerateOutputFiles( TMap *map, char *dir )
{

  FILE  *outputFile;
  char   buffer[100],cmd[100];
  Int_t  nHVEntries, nMVEntries, nTLEntries, nTREntries;
  Int_t  nOKEntries, nStTLEntries, nStTREntries;

  sprintf(cmd,"ls -l %s >/dev/null 2>&1",dir);
  if(gSystem->Exec(cmd)!=0){
    printf("%s --- NOT EXISTS -- create it\n",dir);
    sprintf(cmd,"mkdir %s",dir);
    gSystem->Exec(cmd);
  }

  for( Int_t iLay = 3; iLay < 5; iLay++ )
  {

     Int_t maxLad = ( iLay == 3) ? kNladders3 : kNladders4;
     Int_t maxMod = ( iLay == 3) ? kNmodLad3 : kNmodLad4;

     for(Int_t iLad = 0; iLad < maxLad; iLad++)
     {


        for(Int_t iMod = 0; iMod < maxMod; iMod++)
        {

           sprintf( buffer, "%s/DCSMapContent_SDD_LAYER%i_LADDER%02d_MODULE%d.txt", dir, iLay, iLad, iMod );

           fprintf( stderr, "Creating file %s ... ", buffer );
           outputFile = fopen( buffer, "w" );
            if( outputFile == NULL )
            {
               fprintf( stderr, "failed!\n" );
               return;
            } /*if*/
            else
             fprintf( stderr, "\n" );


           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_HV", iLay, iLad, iMod );
           TObjArray* arrHV = (TObjArray*) map->GetValue( buffer );
           if( arrHV == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nHVEntries = 0;
           }
           else
            nHVEntries = arrHV->GetEntries();

           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_MV", iLay, iLad, iMod );
           TObjArray* arrMV = (TObjArray*) map->GetValue( buffer );
           if( arrMV == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nMVEntries = 0;
           }
           else
            nMVEntries = arrMV->GetEntries();

           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_L", iLay, iLad, iMod );
           TObjArray* arrTL = (TObjArray*) map->GetValue( buffer );
           if( arrTL == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nTLEntries = 0;
           }
           else
            nTLEntries = arrTL->GetEntries();

           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_R", iLay, iLad, iMod );
           TObjArray* arrTR = (TObjArray*) map->GetValue( buffer );
           if( arrTR == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nTREntries = 0;
           }
           else
            nTREntries = arrTR->GetEntries();

           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_L_STATE", iLay, iLad, iMod );
           TObjArray* arrStTL = (TObjArray*) map->GetValue( buffer );
           if( arrStTL == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nStTLEntries = 0;
           }
           else
            nStTLEntries = arrStTL->GetEntries();

           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_TEMP_R_STATE", iLay, iLad, iMod );
           TObjArray* arrStTR = (TObjArray*) map->GetValue( buffer );
           if( arrStTR == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nStTREntries = 0;
           }
           else
            nStTREntries = arrStTR->GetEntries();

           sprintf( buffer, "SDD_LAYER%i_LADDER%02d_MODULE%d_OK", iLay, iLad, iMod );
           TObjArray* arrOK = (TObjArray*) map->GetValue( buffer );
           if( arrOK == NULL ) 
           {
               fprintf( stderr, "Map record %s does not exist!\n", buffer );
               nOKEntries = 0;
           }
           else
            nOKEntries = arrOK->GetEntries();

           fprintf( outputFile, "+-----------------------------------------------------------------------");
           fprintf( outputFile, "------------------------------------------------------------------------+\n" );
           fprintf( outputFile, "|                                                DCS Map content for SDD_LAYER%i_LADDER%02d_MODULE%d" , iLay, iLad, iMod); 
           fprintf( outputFile, "                                                |\n" );
           fprintf( outputFile, "+----------------------+----------------------+---------------------+----------------------+");
           fprintf( outputFile, "-----------------+-----------------+----------------+\n" );
           fprintf( outputFile, "|    %05i  records    |    %05i  records    |    %05i  records   |    %05i  records   |  %05i records  |",
                                 nHVEntries, nMVEntries, nTLEntries, nTREntries, nStTLEntries );
           fprintf( outputFile, "  %05i records  | %05i  records  |\n", nStTREntries, nOKEntries );

           fprintf( outputFile, "|  time (s)     HV     |  time (s)      MV     |  time (s)     TL    |  time (s)     TR    | time (s)   StTL |" );
           fprintf( outputFile, " time (s)   StTR | time (s)   OK   |\n" );
           fprintf( outputFile, "+----------------------+----------------------+---------------------+---------------------+");
           fprintf( outputFile, "-----------------+-----------------+-----------------+\n" );



           Int_t a = (nHVEntries > nMVEntries ) ? nHVEntries : nMVEntries;
           Int_t b = (nTLEntries > nTREntries ) ? nTLEntries : nTREntries;
           Int_t c = (nStTLEntries > nStTREntries ) ? nStTLEntries : nStTREntries;
           if( a < b ) a = b; if( c < nOKEntries ) c = nOKEntries;
           Int_t loopMax = ( a > c ) ? a : c ;
                        // Finds maximal entry number



           for( Int_t entryLoop = 0; entryLoop < loopMax; entryLoop++ )
           {

               if( entryLoop < nHVEntries )
                fprintf( outputFile, "| %12i %4.2f | ", ((AliDCSValue*)arrHV->At(entryLoop))->GetTimeStamp(),
                                                                   ((AliDCSValue*)arrHV->At(entryLoop))->GetFloat() );
               else
                fprintf( outputFile, "|                      | ");

               if( entryLoop < nMVEntries )
               fprintf( outputFile, " %12i  %2.3f | ", ((AliDCSValue*)arrMV->At(entryLoop))->GetTimeStamp(), 
                                                                   ((AliDCSValue*)arrMV->At(entryLoop))->GetFloat() );
               else
               fprintf( outputFile, "                      | ");

               if( entryLoop < nTLEntries )
               fprintf( outputFile, "%12i  %2.2f | ", ((AliDCSValue*)arrTL->At(entryLoop))->GetTimeStamp(), 
                                                                   ((AliDCSValue*)arrTL->At(entryLoop))->GetFloat() );
               else
               fprintf( outputFile, "                    | ");

               if( entryLoop < nTREntries )
               fprintf( outputFile, "%12i  %2.2f | ", ((AliDCSValue*)arrTR->At(entryLoop))->GetTimeStamp(), 
                                                                   ((AliDCSValue*)arrTR->At(entryLoop))->GetFloat() );
               else
               fprintf( outputFile, "                    | ");

               if( entryLoop < nStTLEntries )
               fprintf( outputFile, "%12i  %i | ", ((AliDCSValue*)arrStTL->At(entryLoop))->GetTimeStamp(),
                                                                   ((AliDCSValue*)arrStTL->At(entryLoop))->GetInt() );
               else
               fprintf( outputFile, "                | ");

               if( entryLoop < nStTREntries )
               fprintf( outputFile, "%12i  %i | ",  ((AliDCSValue*)arrStTR->At(entryLoop))->GetTimeStamp(),
                                                                   ((AliDCSValue*)arrStTR->At(entryLoop))->GetInt() );
               else
               fprintf( outputFile, "                | ");

               if( entryLoop < nOKEntries )
               fprintf( outputFile, "%12i  %i |\n",  ((AliDCSValue*)arrOK->At(entryLoop))->GetTimeStamp(), 
                                                                    ((AliDCSValue*)arrOK->At(entryLoop))->GetBool() );
               else
               fprintf( outputFile, "                |\n");

           } /*for( entryLoop )*/

           fclose( outputFile );


        } /*for(iMod)*/


     } /*for(iLad)*/
   } /*for(iLay)*/


} /*GenerateOutputFiles*/
