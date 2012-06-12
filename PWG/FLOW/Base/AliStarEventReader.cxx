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

// Title:   Class for ccessing data from STAR NTuples
//          produces data encapsulated in AliStarEvent and AliStarTrack classes
//
// Origin:  Jim Thomas,        jhthomas@lbl.gov
//          Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch

#include <Riostream.h>

#include <TSystem.h>
#include <TSystemFile.h>
#include <TFile.h>
#include <TList.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TString.h>

#include "AliStarEventReader.h"
#include "AliStarEvent.h"
#include "AliStarTrack.h"

using std::cout;
using std::endl;
using std::ifstream;
ClassImp(AliStarEventReader)

//______________________________________________________________________________
AliStarEventReader::AliStarEventReader():
  TObject(),
  fFileList(NULL),
  fEvent(NULL)
{
  //ctor
}
//______________________________________________________________________________
AliStarEventReader::AliStarEventReader( const char* inputFileDirectory ):
  TObject(),
  fFileList(NULL),
  fEvent(new AliStarEvent(1024))
{
  //ctor
  MakeFileList ( inputFileDirectory ) ;
}

//______________________________________________________________________________
AliStarEventReader::~AliStarEventReader()
{
  //dtor
  delete fFileList;
  delete fEvent;
}

//______________________________________________________________________________
Bool_t AliStarEventReader::GetNextEvent( )
{
  //gets next event
  static TFile*    nextFile    = NULL ;
  static TNtuple*  ntData      = NULL ;
  static Int_t     doOnce      = 0 ;
  static Long64_t  nextEntry   = 0 ;
  static Long64_t  entries     = 0 ;
  static Long64_t  fileCounter = 0 ;

  if ( doOnce == 0 )
  {
    doOnce = 1     ;
    nextFile = (TFile*)  fFileList->First() ;
    if ( nextFile == 0 ) return false      ;
    ntData        = (TNtuple*) ( nextFile->Get("NTtracks") ) ;
    entries       = ntData->GetEntriesFast() ;
    Int_t columns = ntData->GetNvar() ;
    if ( columns != 15 )
    {
      cout << "Error in reading Ntuple file: no. columns != 15" << endl ;
      return false ;
    }
    fileCounter++ ;
    cout << "Start New File " << fileCounter << endl ;
  }

  while ( nextFile )
  {
    while ( nextEntry < entries )
    {
      Float_t* header = NULL;
      Int_t numberOfParticles =  0 ;                   // Number of particle tracks in the next event
      Long64_t headerEntry    =  0 ;                   // Store position of Header and Set Flag in case of EOF or error
      Long64_t skipEvent      =  0 ;                   // Flag in case of wrong number of tracks in this event

      fEvent->Reset();           //reset the event

      // Search for the first "Event" record

      for ( Long64_t j = nextEntry ; j < entries ; j++ )
      {
        Long64_t BytesRead = ntData->GetEntry(j) ;
        if ( BytesRead < 60 )
        {
          cout << "Warning: error in file or EOF " <<  endl ;
          headerEntry = -1 ;
          break ;
        }
        header = ntData->GetArgs() ;
        if ( (int) header[10] == -1 && (int) header[11] == -1 && (int) header[12] == -1 &&
             (int) header[13] == -1 && (int) header[14] == -1 )
        {
          fEvent->SetParams(header);  //set the event params

          numberOfParticles = (int) header[9]  ;   // # of particles passing track cuts, thus in ntuple
          headerEntry = j ;
          break ;
        }
        cout << "Warning: no header entries found in this file" << endl ;
        headerEntry = -1 ;
      }

      if ( headerEntry == -1 ) break ;                 // Break out of main loop if I/O error

      // Get subsequent "track" data
      for ( Long64_t j = headerEntry + 1 ; j < headerEntry + 1 + numberOfParticles  ; j++ )
      {
        Long64_t BytesRead = ntData->GetEntry(j) ;
        if ( BytesRead < 60 )
        {
          cout << "Warning: error in file sequence or EOF" << endl ;
          nextEntry = -1 ;
          break ;
        }
        header = ntData->GetArgs() ;

        if ( TMath::IsNaN(header[10]) == 1 )
        {
          cout << "IsNan ... dEdx will be zeroed out" << endl ;
          header[10] = 0 ;
          header[11] = 999 ;
          header[12] = 999 ;
          header[13] = 999 ;
          header[14] = 999 ;
          cout << header[0]  << " " << header[1]  << " " << header[2]  << " " << header[3]  << " "
               << header[4]  << " " << header[5]  << " " << header[6]  << " " << header[7]  << " "
               << header[8]  << " " << header[9]  << " " << header[10] << " " << header[11] << " "
               << header[12] << " " << header[13] << " " << header[14] << endl ;  // JT test
        }

        if ( (int) header[10] == -1 && (int) header[11] == -1 && (int) header[12] == -1 &&
             (int) header[13] == -1 && (int) header[14] == -1 )
        {
          cout << "Warning: Header in the wrong place, skipping event" << endl ;
          skipEvent = 1 ;
          nextEntry = j ;          // Skip event and freeze nextEntry counter
          break ;
        }

        fEvent->AddTrack( new AliStarTrack(header) );    //add the new track

        nextEntry = j+1 ;
      }
      if ( nextEntry == -1 ) break ;      // Bad record in file, go to next file in fFileList
      if ( skipEvent ==  1 ) continue ;   // Bad event, go to next event in this file
      return true ;                       // Success: Event read OK, note unusual location for a successful return
    }

    nextEntry = 0 ; // this entry goes before nextFile
    nextFile = (TFile*) fFileList->After(nextFile) ;
    if ( nextFile == 0 ) break ;
    if (ntData) delete ntData;
    ntData        = (TNtuple*) ( nextFile->Get("NTtracks") ) ;
    entries       = ntData->GetEntriesFast() ;
    Int_t columns = ntData->GetNvar() ;
    if ( columns != 15 )
    {
      cout << "Error in reading Ntuple file: no. columns != 15" << endl ;
      break ;
    }
    fileCounter++ ;
    cout << "Start New File " << fileCounter << endl ;
  }

  return false ;  // Failure: Error or EOF
}

//______________________________________________________________________________
Bool_t AliStarEventReader::MakeFileList ( const char* input )
{
  //get the files to process
  TString inputstring(input);
  inputstring = inputstring.Strip(TString::kBoth);
  TSystemFile inputfile(inputstring.Data(),"");
  if (inputfile.IsDirectory())
    return MakeFileListFromDir(inputstring.Data());
  else
    return MakeFileListFromFile(inputstring.Data());
}

//______________________________________________________________________________
Bool_t AliStarEventReader::MakeFileListFromDir ( const char* inputFileDirectory )
{
  //get the files to process
  Int_t  Count        = 0 ;
  static Int_t doOnce = 0 ;
  fFileList =  new TList() ;
  void*   directory = gSystem->OpenDirectory(inputFileDirectory) ;
  const char* entry = gSystem->GetDirEntry(directory) ;

  if ( entry == 0 )
  {
    cout << endl <<  "Error: \"" << inputFileDirectory << "\" does not exist" << endl << endl ;
    return false ;
  }
  else cout << endl ;

  while(entry != 0)
  {
    int len = strlen(entry);
    if( len >= 5 && strcmp( &entry[len - 5], ".root" ) == 0 )
    {
      TString fileName ;
      fileName = inputFileDirectory ;
      if( !fileName.EndsWith("/") ) fileName += "/" ;
      fileName += entry;
      fFileList->Add ( TFile::Open(fileName) ) ;
      if ( doOnce == 0 )
      {
        cout << "Add: " << fileName << endl ;
        doOnce = 1 ;
      }
      Count ++ ;
    }
    entry = gSystem->GetDirEntry(directory) ;
  }

  cout << "Add: " << Count-1 << " more file(s) from this directory for a total of " << Count << " files." << endl ;
  cout << "Finished creating file list ... preparing to open first file." << endl << endl ;
  return true ;
}

//______________________________________________________________________________
Bool_t AliStarEventReader::MakeFileListFromFile ( const char* inputFile )
{
  //get the files to process, from a text file, one file per line
  if (!fFileList) fFileList=new TList();
  ifstream filein;
  filein.open(inputFile);
  if (!filein.good()) 
  {
    printf("problem reading the file list \"%s\"\n",inputFile);
    return kFALSE;
  }
  TString line;
  while (filein.good())
  {
    printf("opening file: ");
    line.ReadLine(filein);
    if (line.Length() == 0) continue;
    TFile* file = TFile::Open(line.Data());
    if (!file) 
    {
      printf("problem opening file \"%s\"\n",line.Data());
      continue;
    }
    fFileList->Add(file);
    printf("%s\n",line.Data());
  }
  if (fFileList->GetEntries()>0) return kTRUE;
  return kFALSE;
}

