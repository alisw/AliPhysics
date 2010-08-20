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

ClassImp(AliStarEventReader)

//______________________________________________________________________________
AliStarEventReader::AliStarEventReader():
  TObject(),
  fFileList(NULL),
  fEventHeader(NULL),
  fTracks(NULL),
  fEvent(NULL)
{
  //ctor
}
//______________________________________________________________________________
AliStarEventReader::AliStarEventReader( const char* inputFileDirectory ):
  TObject(),
  fFileList(NULL),
  fEventHeader( new TNtuple("event","event",
                            "runId:eventNumber:VtxX:VtxY:VtxZ:BField:refMult:centralityId:numberOfPrimaryTracks:numberOfParticles:h0:h1:h2:h3:h4" )),
  fTracks( new TNtuple("tracks","tracks",
                       "ID:Charge:Eta:Phi:Pt:Dca:nHits:nHitsFit:nHitsPoss:nHitsDedx:dEdx:nSigElect:nSigPi:nSigK:nSigProton" )),
  fEvent(new AliStarEvent(1024))
{
  //ctor
  MakeFileList ( inputFileDirectory ) ;
}

//______________________________________________________________________________
AliStarEventReader::~AliStarEventReader()
{
  //dtor
  delete fEventHeader;
  fEventHeader = NULL ;
  delete fTracks;
  fTracks = NULL ;
  delete fFileList;
  fFileList = NULL ;
  delete fEvent;
}

//______________________________________________________________________________
Bool_t AliStarEventReader::GetNextEvent( )
{
  //gets next event
  static TFile*    nextFile    = NULL ;
  static TNtuple*  ntData      = NULL ;
  static Int_t     DoOnce      = 0 ;
  static Long64_t  NextEntry   = 0 ;
  static Long64_t  entries     = 0 ;
  static Long64_t  FileCounter = 0 ;

  if ( DoOnce == 0 )
  {
    DoOnce = 1     ;
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
    FileCounter++ ;
    cout << "Start New File " << FileCounter << endl ;
  }

  while ( nextFile )
  {
    while ( NextEntry < entries )
    {
      Float_t* header = NULL;
      Int_t numberOfParticles =  0 ;                   // Number of particle tracks in the next event
      Long64_t HeaderEntry    =  0 ;                   // Store position of Header and Set Flag in case of EOF or error
      Long64_t SkipEvent      =  0 ;                   // Flag in case of wrong number of tracks in this event

      ////////////////////////////////////////////
      fEvent->Reset();           //reset the event
      ////////////////////////////////////////////

      // Search for the first "Event" record

      delete fEventHeader  ;
      fEventHeader  = NULL ;	           // Note that a first 'event'  ntuple was created in constructor
      fEventHeader  = new TNtuple("event","event",
                                  "runId:eventNumber:VtxX:VtxY:VtxZ:BField:refMult:centralityId:numberOfPrimaryTracks:numberOfParticles:h0:h1:h2:h3:h4" ) ;

      for ( Long64_t j = NextEntry ; j < entries ; j++ )
      {
        Long64_t BytesRead = ntData->GetEntry(j) ;
        if ( BytesRead < 60 )
        {
          cout << "Warning: error in file or EOF " <<  endl ;
          HeaderEntry = -1 ;
          break ;
        }
        header = ntData->GetArgs() ;
        if ( (int) header[10] == -1 && (int) header[11] == -1 && (int) header[12] == -1 &&
             (int) header[13] == -1 && (int) header[14] == -1 )
        {
          fEventHeader->Fill(header) ;

          //////////////////////////////////////////////////
          fEvent->SetParams(header);  //set the event params
          //////////////////////////////////////////////////

          numberOfParticles = (int) header[9]  ;   // # of particles passing track cuts, thus in ntuple
          HeaderEntry = j ;
          break ;
        }
        cout << "Warning: no header entries found in this file" << endl ;
        HeaderEntry = -1 ;
      }

      if ( HeaderEntry == -1 ) break ;                 // Break out of main loop if I/O error

      // Get subsequent "track" data

      delete fTracks ;
      fTracks = NULL ;                 // Note that a first 'tracks' ntuple was created in constructor
      fTracks = new TNtuple("tracks","tracks",
                            "ID:Charge:Eta:Phi:Pt:Dca:nHits:nHitsFit:nHitsPoss:nHitsDedx:dEdx:nSigElect:nSigPi:nSigK:nSigProton" ) ;

      for ( Long64_t j = HeaderEntry + 1 ; j < HeaderEntry + 1 + numberOfParticles  ; j++ )
      {
        Long64_t BytesRead = ntData->GetEntry(j) ;
        if ( BytesRead < 60 )
        {
          cout << "Warning: error in file sequence or EOF" << endl ;
          NextEntry = -1 ;
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
          SkipEvent = 1 ;
          NextEntry = j ;          // Skip event and freeze NextEntry counter
          break ;
        }

        fTracks->Fill(header) ;

        ////////////////////////////////////////////////////////////////////
        fEvent->AddTrack( new AliStarTrack(header) );    //add the new track
        ////////////////////////////////////////////////////////////////////

        NextEntry = j+1 ;
      }
      if ( NextEntry == -1 ) break ;      // Bad record in file, go to next file in fFileList
      if ( SkipEvent ==  1 ) continue ;   // Bad event, go to next event in this file
      return true ;                       // Success: Event read OK, note unusual location for a successful return
    }

    NextEntry = 0 ; // this entry goes before nextFile
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
    FileCounter++ ;
    cout << "Start New File " << FileCounter << endl ;
  }

  return false ;  // Failure: Error or EOF
}

//______________________________________________________________________________
Bool_t AliStarEventReader::AcceptEvent( AliStarEvent* event )
{
  // Cut parameters for each event

  const Float_t VertexXMin  =  -1.0 ;  // cm
  const Float_t VertexXMax  =   1.0 ;  // cm
  const Float_t VertexYMin  =  -1.0 ;  // cm
  const Float_t VertexYMax  =   1.0 ;  // cm
  const Float_t VertexZMin  = -30.0 ;  // cm
  const Float_t VertexZMax  =  30.0 ;  // cm
  const Int_t   MultMin     =    10 ;  // Note: this is a cut on refMult which is not your ordinary multiplicity.
  const Int_t   MultMax     =  1000 ;  // Refmult corresponds to Zhangbu's ascii table of cuts for centrality bins.
  const Int_t   BlackEvent  =  3000 ;  // Maximum number of primary tracks allowed in one event (more is an error)

  Int_t refMult               = event->GetRefMult();  // Reference Multiplicity for Centrality determination
  Int_t NumberOfPrimaryTracks = event->GetNumberOfPrimaryTracks();  // Number of primary tracks in full event (not just those in main TPC)

  if ( refMult < MultMin          || refMult > MultMax )                   return false ;
  if ( NumberOfPrimaryTracks <= 0 || NumberOfPrimaryTracks > BlackEvent )  return false ;

  // Cut on Vertex location
  Float_t vertex[3] = { event->GetVtxX(),
                        event->GetVtxY(),
                        event->GetVtxZ() };

  if ( vertex[0] < VertexXMin || vertex[0] > VertexXMax )    return false ;  // Skip events that fall outside Vtx cuts
  if ( vertex[1] < VertexYMin || vertex[1] > VertexYMax )    return false ;
  if ( vertex[2] < VertexZMin || vertex[2] > VertexZMax )    return false ;

  return true ;
}

//______________________________________________________________________________
Bool_t AliStarEventReader::AcceptTrack( AliStarTrack* track )
{
  // Cut Parameters for individual tracks

  const Float_t dcaCut      =   3.0  ;      // cm
  const Float_t PtMin       =   0.15 ;      // GeV
  const Float_t PtMax       =   2.0  ;      // GeV
  const Float_t EtaMin      =  -1.1  ;
  const Float_t EtaMax      =   1.1  ;
  const Float_t FitRatio    =  0.52  ;      // Number of hits over number of hits possible
  const Int_t   nHitMin     =    15  ;      // 15 is typical but sometimes goes as high as 25
  const Int_t   nHitMax     =   100  ;      // 45 pad rows in the TPC and so anything bigger than 45+Silicon is infinite
  const Int_t   nHitPossMin =     5  ;      // Don't bother to fit tracks if # possible hits is too low, also protect / 0

  if ( track->GetDCA() >  dcaCut ) return false ;                 // magnitude of 3D DCA for global tracks
  if ( track->GetNHitsPoss() <  nHitPossMin ||
       track->GetNHitsPoss() >  nHitMax ) return false ;   // Minimum number of Possible hits, see above.
  if ( track->GetEta() <  EtaMin  || track->GetEta() >  EtaMax  ) return false ;
  if ( track->GetPt()  <  PtMin   || track->GetEta() >  PtMax   ) return false ;
  if ( track->GetNHitsFit() < nHitMin || track->GetNHitsFit() > nHitMax ) return false ;
  if ( ((1.0*track->GetNHitsFit()) /
        (1.0*track->GetNHitsPoss())) < FitRatio )  return false ;

  return true ;
}

//______________________________________________________________________________
Int_t AliStarEventReader::ParticleID( AliStarTrack* track )
{
  // Note: This is a very simple PID selection scheme.  More elaborate methods (with multiple cuts) may be required.
  // When you *are* using dEdx information, you must chose a finite number of good Dedx hits ... but the limit should
  // be about 2/3 of nHitsMin.  This is because some clusters do not form good dEdx hits due to track
  // merging, etc., and so nHitsDedx is always less than nHitsFit.  A rule of thumb says ~2/3 ratio.

  Int_t ID = 0 ;

  const Int_t   nHitDedxMin =    15  ;       // 10 to 20 is typical.  nHitDedxMin is often chosen to be about 2/3 of nHitMin.
  const Float_t nSigmaPID   =    2.0 ;       // Number of Sigma cut to apply to PID bands

  // Test on Number of dE/dx hits required, return 0 if not enough hits
  if ( track->GetNHitsDedx() <  nHitDedxMin ) return ID;

  // Begin PID

  if ( TMath::Abs( track->GetNSigElect() ) >= nSigmaPID )
  {
    if ( TMath::Abs( track->GetNSigK()  ) <= nSigmaPID )
    {
      ID = 321  ;
    }
    if ( TMath::Abs( track->GetNSigProton()  ) <= nSigmaPID )
    {
      ID = 2212 ;
    }
    if ( TMath::Abs( track->GetNSigPi()  ) <= nSigmaPID )
    {
      ID = 211  ;
    }
  }

  // Pion is the default in case of ambiguity because it is most abundent. Don't re-arrange order, above.

  return ID ;
}

//______________________________________________________________________________
Int_t AliStarEventReader::Centrality( Int_t referenceMultiplicity )

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


  Int_t   CentralityBins  [] = { 14 , 31 , 57 , 96 , 150 , 222 , 319 , 441 , 520 , 1000 } ;  // Run4 200 GeV
  Int_t   MiddleBinID     [] = {  0 ,  1 ,  2 ,  3 ,   4 ,   5 ,   6 ,   7 ,   8 ,    9 } ;  // ID Number
  //Info  MiddleBinPercent[] = { 85., 75., 65., 55.,  45.,  35.,  25.,  15., 7.5 ,  2.5 } ;  // Percent
  Int_t   myCentrality  ;

  if      ( referenceMultiplicity < CentralityBins[0] )
  {
    myCentrality = MiddleBinID[0] ;
  }
  else if ( referenceMultiplicity < CentralityBins[1] )
  {
    myCentrality = MiddleBinID[1] ;
  }
  else if ( referenceMultiplicity < CentralityBins[2] )
  {
    myCentrality = MiddleBinID[2] ;
  }
  else if ( referenceMultiplicity < CentralityBins[3] )
  {
    myCentrality = MiddleBinID[3] ;
  }
  else if ( referenceMultiplicity < CentralityBins[4] )
  {
    myCentrality = MiddleBinID[4] ;
  }
  else if ( referenceMultiplicity < CentralityBins[5] )
  {
    myCentrality = MiddleBinID[5] ;
  }
  else if ( referenceMultiplicity < CentralityBins[6] )
  {
    myCentrality = MiddleBinID[6] ;
  }
  else if ( referenceMultiplicity < CentralityBins[7] )
  {
    myCentrality = MiddleBinID[7] ;
  }
  else if ( referenceMultiplicity < CentralityBins[8] )
  {
    myCentrality = MiddleBinID[8] ;
  }
  else
  {
    myCentrality = MiddleBinID[9] ;
  }

  return myCentrality ;
}

//______________________________________________________________________________
void AliStarEventReader::PrintEventHeader ( )
{
  // TNtuple* event: names are documented in the next 2 lines
  // event  = new TNtuple("event","event",
  //   "runId:eventNumber:VtxX:VtxY:VtxZ:BField:refMult:centralityId:numberOfPrimaryTracks:numberOfParticles:h0:h1:h2:h3:h4" ) ;
  //
  Float_t  *eventhd ;
  eventhd = fEventHeader->GetArgs() ;

  Int_t   runId                  = (int)   eventhd[0]  ;
  Int_t   eventNumber            = (int)   eventhd[1]  ;
  Float_t   primaryVertexPosition[3] = { (float) eventhd[2],  // (cm)
                                         (float) eventhd[3],  // (cm)
                                         (float) eventhd[4] };  // (cm)
  Float_t magneticField          = (float) eventhd[5]  ;  // kilogauss
  Int_t   refMult                = (int)   eventhd[6]  ;  // Raw Mult into 0.5 unit: a relative number, not total Mult.
  Int_t   centralityId           = (int)   eventhd[7]  ;  // STAR centrality bin # based on refMult
  Int_t   numberOfPrimaryTracks  = (int)   eventhd[8]  ;  // # of primaries, including FTPC tracks which are not in ntuple
  Int_t   numberOfParticles      = (int)   eventhd[9]  ;  // # of particles passing track cuts, thus in ntuple

  printf("\n  runId eventNo    VtxX    VtxY    VtxZ  MagFld  refMult  CentBin  #PrimeTrak  #Tracks \n") ;
  printf("%7d  %6d %7.4f %7.4f %7.3f  %6.3f   %6d   %6d      %6d   %6d \n\n",
         runId, eventNumber, primaryVertexPosition[0], primaryVertexPosition[1], primaryVertexPosition[2],
         magneticField, refMult, centralityId, numberOfPrimaryTracks, numberOfParticles ) ;

  //Int_t newCentralityID ;
  //newCentralityID = Centrality( refMult) ;              // Should be the same as "centralityID" from tape
  //cout << "Test: should be the same " << newCentralityID << " " << centralityId << endl ;  // JT test
}

//______________________________________________________________________________
void AliStarEventReader::PrintTrack ( Int_t counter )
{
  // TNtuple* track: names are documented in the next 2 lines
  // tracks = new TNtuple("tracks","tracks",
  //   "ID:Charge:Eta:Phi:Pt:Dca:nHits:nHitsFit:nHitsPoss:nHitsDedx:dEdx:nSigElect:nSigPi:nSigK:nSigProton" ) ;
  //
  if ( counter == 0 )
  {
    printf(
      "    id charge     eta     phi      pt     dca  nHits  nFit nPoss ndEdx   dEdx nSigElec nSigPi  nSigK nSigPr\n") ;
  }
  Float_t* data = fTracks -> GetArgs()       ;  // Extract data from the track
  Int_t   id             = (int)   data[0]   ;  // id - a unique integer for each track in this event
  Int_t   charge         = (int)   data[1]   ;  // +1 or -1
  Float_t eta            = (float) data[2]   ;  // Pseudo-rapidity at the vertex
  Float_t phi            = (float) data[3]   ;  // Phi emission angle at the vertexcd
  Float_t pt             = (float) data[4]   ;  // Pt of the track at the vertex
  Float_t dca            = (float) data[5]   ;  // Magnitude of 3D DCA to vertex
  Int_t   nHits          = (int)   data[6]   ;  // Number of clusters available to the Kalman Filter
  Int_t   nHitsFit       = (int)   data[7]   ;  // Number of clusters used in the fit (after cuts)
  Int_t   nHitsPoss      = (int)   data[8]   ;  // Number of possible cluster on track (theoretical max)
  Int_t   nHitsDedx      = (int)   data[9]   ;  // Number of clusters used in the fit (after dEdx cuts)
  Float_t dEdx           = 1.e6*(float)data[10]  ;  // Measured dEdx (Note: GeV/cm so convert to keV/cm!!)
  Float_t nSigmaElectron = (float) data[11]  ;  // Number of sigma from electron Bethe-Bloch curve
  Float_t nSigmaPion     = (float) data[12]  ;  // Number of sigma from pion Bethe-Bloch curve
  Float_t nSigmaKaon     = (float) data[13]  ;  // Number of sigma from kaon Bethe-Bloch curve
  Float_t nSigmaProton   = (float) data[14]  ;  // Number of sigma from proton Bethe-Bloch curve

  // Alternative way to access the data
  nHitsPoss      = (int) ( fTracks->GetLeaf("nHitsPoss")->GetValue() ) ;  // Note alternative method to retrieve data
  // Using the definition of the original NTuple
  // TrackTuple      = new TNtuple("NTtracks","NTtracks",
  // "ID:Charge:Eta:Phi:Pt:Dca:nHits:nHitsFit:nHitsPoss:nHitsDedx:dEdx:nSigElect:nSigPi:nSigK:nSigProton" )

  printf("%6d %4d   %7.3f %7.3f %7.3f %7.4f %6d %5d %5d %5d %6.2f   %6.2f %6.2f %6.2f %6.2f \n",
         id, charge, eta, phi, pt, dca, nHits, nHitsFit, nHitsPoss, nHitsDedx, dEdx,
         nSigmaElectron, nSigmaPion, nSigmaKaon, nSigmaProton ) ;
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
  static Int_t DoOnce = 0 ;
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
      if ( DoOnce == 0 )
      {
        cout << "Add: " << fileName << endl ;
        DoOnce = 1 ;
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

