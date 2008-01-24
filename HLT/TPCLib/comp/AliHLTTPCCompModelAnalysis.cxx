// $Id: AliHLTTPCCompModelAnalysis.cxx,v 1.2 2006/08/10 09:46:51 richterm Exp $

/**************************************************************************
 * TPCCompModelAnalysisright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: J. Wagner <jwagner@cern.ch>                                   *
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

/** @file   AliHLTTPCCompModelAnalysis.cxx
    @author J. Wagner jwagner@cern.ch
    @date   17-11-2007
    @brief  A processing analysis component for the HLT */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelAnalysis.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCCompDataCompressorHelper.h"
#include "TFile.h"
#include "TH1.h"
#include <cerrno>

/** constructor **/
AliHLTTPCCompModelAnalysis::AliHLTTPCCompModelAnalysis(Bool_t modelanalysis, Bool_t trackanalysis, TString dumpfilename, TString graphfilename):
  fModelAnalysis(modelanalysis),
  fTrackAnalysis(trackanalysis),
  fDumpFileName(dumpfilename),
  fGraphFileName(graphfilename),
  fFirstTrackArray(),
  fSecondTrackArray(),
  fFirstTrackList(NULL),
  fSecondTrackList(NULL),
  fFirstTrashTracks(0),
  fSecondTrashTracks(0),
  fTotalComparedTracks(0),
  fMatchedFirstTrashTracks(0),
  fMatchedSecondTrashTracks(0),
  fFirstUnmatchedTracks(0),
  fSecondUnmatchedTracks(0),
  fToleranceDeviation(0.0),
  fTrackListPointer(NULL),
  fTotalDiscardedClusters(0),
  fValuableDiscardedClusters(0),
  fTrashTracks(0)
{
  // see header file for class documentation
}

/** destructor **/
AliHLTTPCCompModelAnalysis::~AliHLTTPCCompModelAnalysis()
    {
    for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	    {
	    if ( fDiscardedClusters[slice][patch]->fSpacePoints )
		{
		delete [] fDiscardedClusters[slice][patch];
		}
	    }
    }

/** initialise arrays for tracks/ discarded clusters depending on model-flags **/
Int_t AliHLTTPCCompModelAnalysis::Init()
{
  // see header file for class documentation

  if(fTrackAnalysis) // track quantities are to be initialised
    {
      fFirstTrackArray.Reset();
      fSecondTrackArray.Reset();

      fFirstTrackList = NULL;
      fSecondTrackList = NULL;

      fFirstTrashTracks = 0;
      fSecondTrashTracks = 0;

      fTotalComparedTracks = 0;
      fMatchedFirstTrashTracks = 0;
      fMatchedSecondTrashTracks = 0;

      fFirstUnmatchedTracks = 0;
      fSecondUnmatchedTracks = 0;

      fToleranceDeviation = 0.001;

    }
  
  if(fModelAnalysis) // cluster array to be initialised
    {
      for ( UInt_t slice=0; slice<36; slice++ )
	{
	  for ( UInt_t patch=0; patch<6; patch++ )
	    {
	      fDiscardedClusters[slice][patch] = NULL;
	    }
	}

      // initialise trash track list to store discarded tracks
      fTrackListPointer = NULL;

      // set all counters to zero:
      fTrashTracks = 0;
      fTotalDiscardedClusters = 0;
      fValuableDiscardedClusters = 0;
    }

  return 0;
}

Int_t AliHLTTPCCompModelAnalysis::DisplayResults()
{
  // see header file for class documentation
  HLTInfo("--------------------DISPLAYING RESULTS---------------------");
  // if model loss analysis, then display these results, else display track results
  if (fModelAnalysis)
    {
      DisplayModelResults();
    };

  if (fTrackAnalysis)
    {
      DisplayTrackResults();
    };

  // error message: if no analysis flag is switched on
  if( (!fModelAnalysis) && (!fTrackAnalysis) )
    {
      HLTError("Error! Display Results called without any analysis flag switched on.");
      return 1;
    }

  return 0;
}

Int_t AliHLTTPCCompModelAnalysis::SetTracks(AliHLTTPCTrackletData* tracklets, Bool_t fillingfirsttracks)
{
  // see header file for class documentation
  // if fFillingFirstTrackArray is true (i.e. first input file in component is processed)
  // first input track array is filled
  // else the second track array is filled

  if(fillingfirsttracks)
    {
      HLTDebug( "Reading %u tracks in first array", (unsigned) tracklets->fTrackletCnt );
      
      if(tracklets->fTrackletCnt == 0)
	{
	  HLTError("Error! No tracklets to fill into first track array!");
	  return EINVAL;
	}
      else
	{
	  fFirstTrackArray.FillTracks(tracklets->fTrackletCnt, tracklets->fTracklets );
	}

    }
  else
    {
      if(tracklets->fTrackletCnt == 0)
	{
	  HLTError("Error! No tracklets to fill into second track array!");
	  return EINVAL;
	}
      else
	{
	  fSecondTrackArray.FillTracks(tracklets->fTrackletCnt, tracklets->fTracklets );
	}
      // read in tracks in second array (tracks from clusters after Vestbo in second tracking process)
      HLTDebug( "Reading %u tracks in second array", (unsigned) tracklets->fTrackletCnt );
    
    }

  return 0;
}

Int_t AliHLTTPCCompModelAnalysis::CompareTracks()
{
  // see header file for class documentation
  // define variables for number of tracks in track arrays:
  Int_t firsttracks = fFirstTrackArray.GetNTracks();
  Int_t secondtracks = fSecondTrackArray.GetNTracks();

  // error checking: if second array has been filled or not:
  if(firsttracks == 0)
    {
      HLTError("No tracks in first track array!");
      return EINVAL;
    };
 

  if(secondtracks == 0)
    {
      HLTError("No tracks in second track array!");
      return EINVAL;
    };

  // take track from first tracking,
  for(Int_t ii=0; ii < firsttracks; ii++)
    {
      // build track list for all tracks in first array
      AliHLTTPCTrackList* currenttrackentry = new AliHLTTPCTrackList;
      currenttrackentry->track = *(fFirstTrackArray.GetCheckedTrack(ii));

      // get its pythia information, 
      currenttrackentry->pythiatrack = GetComparableTrackPythiaInfo(currenttrackentry->track);

      currenttrackentry->matchingindicator = 0;

      // put this element as first in list
      currenttrackentry->next = fFirstTrackList;
      fFirstTrackList = currenttrackentry;

      // count tracks below 0.1GeV
      if(currenttrackentry->track.GetPt()<0.1)
	{
	  ++fFirstTrashTracks;
	}

    }

 // take track from second tracking,
  for(Int_t ii=0; ii < secondtracks; ii++)
    {
      // build track list for all tracks in second array
      AliHLTTPCTrackList* currenttrackentry = new AliHLTTPCTrackList;
      currenttrackentry->track = *(fSecondTrackArray.GetCheckedTrack(ii));

      // get its pythia information, 
      currenttrackentry->pythiatrack = GetComparableTrackPythiaInfo(currenttrackentry->track);

      // put this element as first in list
      currenttrackentry->next = fSecondTrackList;
      fSecondTrackList = currenttrackentry;

      // count tracks below 0.1GeV
      if(currenttrackentry->track.GetPt()<0.1)
	{
	  ++fSecondTrashTracks;
	}

    }

  // search for matching track from secondary tracking
  AliHLTTPCTrackList* firstmatchpointer = fFirstTrackList;
 
  while(firstmatchpointer != NULL)
   {
     AliHLTTPCTrackList* secondmatchpointer = fSecondTrackList;

     while(secondmatchpointer != NULL)
       {
	
	 // compare paramters of the two tracks, 
	 // match only when coincidence >= 50% in fToleranceDeviation range!
	 if ((CompareTrackInfo(firstmatchpointer, secondmatchpointer) > 4))
	   { 
    
	     if((CompareTrackInfo(firstmatchpointer, secondmatchpointer) > firstmatchpointer->matchingindicator))
	       {
		 // look if current better matching track has already been matched before
		 if((CompareTrackInfo(firstmatchpointer, secondmatchpointer) > secondmatchpointer->matchingindicator))
		   {
		     
		     // set previously assigned matchingindicator (if there was one) of secondary track back to zero
		     if(firstmatchpointer->matchingindicator > 0)
		       {
			 firstmatchpointer->matchingtrack->matchingindicator = 0;
			 firstmatchpointer->matchingtrack->matchingtrack = NULL;
		       }

		     if(secondmatchpointer->matchingindicator > 0)
		       {
			 secondmatchpointer->matchingtrack->matchingindicator = 0;
			 secondmatchpointer->matchingtrack->matchingtrack = NULL;
		       }
		     
		     // compare according to tracks themselves (other possibility: compare pythiatracks - better!)
		     secondmatchpointer->matchingindicator = CompareTrackInfo(firstmatchpointer, secondmatchpointer) ;
		     firstmatchpointer->matchingindicator = CompareTrackInfo(firstmatchpointer, secondmatchpointer);
		     
		     // remember which track matches which
		     secondmatchpointer->matchingtrack = firstmatchpointer;
		     firstmatchpointer->matchingtrack = secondmatchpointer;

		   } // end if compare > second matching indicator

	       } // end if compare > first matching indicator
	     
	   }// end if compare > 4
	 	 
	 secondmatchpointer = secondmatchpointer->next;
       }
     
     // go on with next original track
     firstmatchpointer = firstmatchpointer->next;
   }
  
  // count not matched tracks in first and second track list
  AliHLTTPCTrackList* nomatchcounter = fFirstTrackList;

  while(nomatchcounter != NULL)
    {
      if(nomatchcounter->matchingindicator == 0)
	{
	  ++fFirstUnmatchedTracks;
	}
      else
	{
	 ++fTotalComparedTracks;

	 // count matched trash tracks
	 if(nomatchcounter->track.GetPt() < 0.1)
	   {
	     ++fMatchedFirstTrashTracks;
	   };
	}
      nomatchcounter = nomatchcounter->next;
    }

  nomatchcounter = fSecondTrackList;
  while(nomatchcounter != NULL)
    {
      if(nomatchcounter->matchingindicator == 0)
	{
	  ++fSecondUnmatchedTracks;
	}
      else
	{
	  // count matched trash tracks
	 if(nomatchcounter->track.GetPt() < 0.1)
	   {
	     ++fMatchedSecondTrashTracks;
	   };
	}
    
      nomatchcounter = nomatchcounter->next;
    }
 
  // consistency check: fFirstUnmatchedTracks + fTotalComparedTracks = # of tracks in first array
  // ...and analogously for second array:
  if(fFirstUnmatchedTracks + fTotalComparedTracks !=  firsttracks)
    {
      HLTWarning("Warning! Possible inconsistency in original track array: Number of compared and unmatched tracks not equal to total number of tracks!");
    };
  
  if(fSecondUnmatchedTracks + fTotalComparedTracks !=  secondtracks)
    {
      HLTWarning("Warning! Possible inconsistency in second track array: Number of compared and unmatched tracks not equal to total number of tracks!");
    };
  
  return 0;
}

AliHLTTPCTrack AliHLTTPCCompModelAnalysis::GetComparableTrackPythiaInfo(AliHLTTPCTrack comparabletrack)
{
  // see headerfile for class documentation

  AliHLTTPCTrack pythiatrack = comparabletrack;

  return pythiatrack;
}

Int_t AliHLTTPCCompModelAnalysis::MarkTrashTrack(AliHLTTPCTrack *lowpttrack)
{
  // see header file for class documentation

  // save track first in lowpttrack list (all lowpttracks are displayed in display function altogether)
  AliHLTTPCTrackList* tracklistentry =  new AliHLTTPCTrackList;
  tracklistentry->track = *lowpttrack;
  tracklistentry->wronglydiscarded = GetTrashTrackPythiaInfo(lowpttrack);
  tracklistentry->matchingindicator = 0; // not needed here, therefore initialised to zero
  tracklistentry->next = fTrackListPointer;
  tracklistentry->matchingtrack = NULL; // not needed here, therefore initialised to NULL
  fTrackListPointer = tracklistentry;

  ++fTrashTracks;

  return 0;
}

Int_t AliHLTTPCCompModelAnalysis::MarkTrashCluster(AliHLTTPCClusterData *discardedcluster, UInt_t slice, UInt_t patch)
{
  // see header file for class documentation

  // get Pythia information of discarded cluster
  Bool_t wronglydiscarded = GetClusterPythiaInfo(discardedcluster);

  // if cluster has been discarded wrongly, save information
  if(wronglydiscarded)
    {
      // store cluster
      fDiscardedClusters[slice][patch] = discardedcluster;
      // increase number of valuable discarded clusters
      ++fValuableDiscardedClusters;
    }
      
      ++fTotalDiscardedClusters;   
      
  return 0;
}

Bool_t AliHLTTPCCompModelAnalysis::GetTrashTrackPythiaInfo(AliHLTTPCTrack* discardedtrack)
{
  // see header file for class documentation
  // store information from pythia in current track list entry
  // fTrackListPointer.pythiatrack = FillFromPythia...

  return 0;
}

Bool_t AliHLTTPCCompModelAnalysis::GetClusterPythiaInfo(AliHLTTPCClusterData* discardedcluster)
{
  // see header file for class documentation
  // Pythia information can be
  // either: cluster belongs to discarded track with pt < 0.1 Gev (--> cluster correctly discarded)
  //     or: cluster is discarded and does not belong to any pythia track (--> correctly discarded)
  //     or: cluster is discarded but belongs to pythia track (--> cluster WRONGLY discarded!!!)

  return 0;
}
 
Int_t AliHLTTPCCompModelAnalysis::CompareTrackInfo(AliHLTTPCTrackList* firsttracklistelement, AliHLTTPCTrackList* secondtracklistelement)
{
  // see header file for class documentation
  // calculate matching indicator accoring to the track information
  // ++matchingindicator for every paramter that matches

  Int_t currentmatchingindicator = 0;

  // tolerance range of 1 percent deviation for each quantity

  // compare start point (x,y,z)
  if(abs((firsttracklistelement->track.GetFirstPointX() - secondtracklistelement->track.GetFirstPointX()))/firsttracklistelement->track.GetFirstPointX() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->track.GetFirstPointY() - secondtracklistelement->track.GetFirstPointY()))/firsttracklistelement->track.GetFirstPointY() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->track.GetFirstPointZ() - secondtracklistelement->track.GetFirstPointZ()))/firsttracklistelement->track.GetFirstPointZ() <= fToleranceDeviation)
    ++currentmatchingindicator;
  
  // compare end point
  if(abs((firsttracklistelement->track.GetLastPointX() - secondtracklistelement->track.GetLastPointX()))/firsttracklistelement->track.GetLastPointX() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->track.GetLastPointY() - secondtracklistelement->track.GetLastPointY()))/firsttracklistelement->track.GetLastPointY() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->track.GetLastPointZ() - secondtracklistelement->track.GetLastPointZ()))/firsttracklistelement->track.GetLastPointZ() <= fToleranceDeviation)
    ++currentmatchingindicator;

  // compare pt, psi, tgl
  if(abs((firsttracklistelement->track.GetPt() - secondtracklistelement->track.GetPt()))/firsttracklistelement->track.GetPt() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->track.GetPsi() - secondtracklistelement->track.GetPsi()))/firsttracklistelement->track.GetPsi() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->track.GetTgl() - secondtracklistelement->track.GetTgl()))/firsttracklistelement->track.GetTgl() <= fToleranceDeviation)
    ++currentmatchingindicator;

  // compare number of assigned cluster hits
  if(abs((firsttracklistelement->track.GetNHits() - secondtracklistelement->track.GetNHits()))/firsttracklistelement->track.GetNHits() <= fToleranceDeviation)
    ++currentmatchingindicator;

  return currentmatchingindicator;
}

Int_t AliHLTTPCCompModelAnalysis::ComparePythiaTrackInfo(AliHLTTPCTrackList* firsttracklistelement, AliHLTTPCTrackList* secondtracklistelement)
{
  // see header file for class documentation
  // calculate matching indicator accoring to the track information
  // ++matchingindicator for every paramter that matches

  Int_t currentmatchingindicator = 0;

  // tolerance range of 1 percent deviation for each quantity

  // compare start point (x,y,z)
  if(firsttracklistelement->pythiatrack.GetFirstPointX() == secondtracklistelement->pythiatrack.GetFirstPointX())
    ++currentmatchingindicator;

  if(firsttracklistelement->pythiatrack.GetFirstPointY() == secondtracklistelement->pythiatrack.GetFirstPointY())
    ++currentmatchingindicator;

  if(firsttracklistelement->pythiatrack.GetFirstPointZ() == secondtracklistelement->pythiatrack.GetFirstPointZ())
    ++currentmatchingindicator;
  
  // compare end point
  if(firsttracklistelement->pythiatrack.GetLastPointX() == secondtracklistelement->pythiatrack.GetLastPointX())
    ++currentmatchingindicator;

  if(firsttracklistelement->pythiatrack.GetLastPointY() == secondtracklistelement->pythiatrack.GetLastPointY())
    ++currentmatchingindicator;

  if(firsttracklistelement->pythiatrack.GetLastPointZ() == secondtracklistelement->pythiatrack.GetLastPointZ())
    ++currentmatchingindicator;

  // compare pt, psi, tgl
  if(firsttracklistelement->pythiatrack.GetPt() == secondtracklistelement->pythiatrack.GetPt())
    ++currentmatchingindicator;

  if(firsttracklistelement->pythiatrack.GetPsi() == secondtracklistelement->pythiatrack.GetPsi())
    ++currentmatchingindicator;

  if(firsttracklistelement->pythiatrack.GetTgl() == secondtracklistelement->pythiatrack.GetTgl())
    ++currentmatchingindicator;

  // compare number of assigned cluster hits
  if(firsttracklistelement->pythiatrack.GetNHits() == secondtracklistelement->pythiatrack.GetNHits())
    ++currentmatchingindicator;

  return currentmatchingindicator;
}

Int_t AliHLTTPCCompModelAnalysis::CreateGraphs(Bool_t relativedifferences)
{
  // see header file for class documentation
  AliHLTTPCTrackList* tracklistpointer = fFirstTrackList;

  AliHLTTPCTrackList* trackmatchingpointer;

  // set up histogram ranges
  Double_t difffirstxmin, difffirstxmax;
  Double_t difffirstymin, difffirstymax;
  Double_t difffirstzmin, difffirstzmax;
  Int_t difffirstxbins, difffirstybins, difffirstzbins;

  Double_t difflastxmin, difflastxmax;
  Double_t difflastymin, difflastymax;
  Double_t difflastzmin, difflastzmax;
  Int_t difflastxbins, difflastybins, difflastzbins; 

  Double_t diffptmin, diffptmax;
  Double_t diffpsimin, diffpsimax;
  Double_t difftglmin, difftglmax;
  Int_t diffptbins, diffpsibins, difftglbins;

  Double_t diffclustermin, diffclustermax;
  Int_t diffclusterbins;
 
  // resolution histograms (currently not working since pterr = 0 for every track!)
  Double_t diffpterrmin, diffpterrmax;
  Double_t diffpsierrmin, diffpsierrmax;
  Double_t difftglerrmin, difftglerrmax;
  Int_t diffpterrbins, diffpsierrbins, difftglerrbins;

  if(!relativedifferences)
    {
      difffirstxmin = -2;
      difffirstxmax = +2;
      difffirstxbins = (Int_t) ((difffirstxmax - difffirstxmin)/0.0001);

      difffirstymin = -2;
      difffirstymax = +2;
      difffirstybins = (Int_t) ((difffirstymax - difffirstymin)/0.0001);

      difffirstzmin = -2;
      difffirstzmax = +2;
      difffirstzbins = (Int_t) ((difffirstzmax - difffirstzmin)/0.0001);

      difflastxmin = -2;
      difflastxmax = +2;
      difflastxbins = (Int_t) ((difflastxmax - difflastxmin)/0.0001);

      difflastymin = -2;
      difflastymax = +2;
      difflastybins = (Int_t) ((difflastymax - difflastymin)/0.0001);

      difflastzmin = -2;
      difflastzmax = +2;
      difflastzbins = (Int_t) ((difflastzmax - difflastzmin)/0.0001);

      diffptmin = -1;
      diffptmax = +1;
      diffptbins = (Int_t) ((diffptmax - diffptmin)/0.0001);
      
      diffpsimin = -1;
      diffpsimax = +1;
      diffpsibins = (Int_t) ((diffpsimax - diffpsimin)/0.0001);

      difftglmin = -1;
      difftglmax = +1;
      difftglbins = (Int_t) ((difftglmax - difftglmin)/0.0001);

      diffclustermin = -50;
      diffclustermax = +50;
      diffclusterbins = (Int_t) ((diffclustermax - diffclustermin)/1);

#if 0
      diffpterrmin = -1;
      diffpterrmax = 1;
      diffpterrbins = (Int_t) ((diffpterrmax - diffpterrmin)/1);

      diffpsierrmin = -1;
      diffpsierrmax = 1;
      diffpsierrbins = (Int_t) ((diffpsierrmax - diffpsierrmin)/1);

      difftglerrmin = -1;
      difftglerrmax = 1;
      difftglerrbins = (Int_t) ((difftglerrmax - difftglerrmin)/1);
#endif

    }
  else
    {
      difffirstxmin = -1;
      difffirstxmax = +1;
      difffirstxbins = (Int_t) ((difffirstxmax - difffirstxmin)/0.0001);

      difffirstymin = -1;
      difffirstymax = +1;
      difffirstybins = (Int_t) ((difffirstymax - difffirstymin)/0.0001);

      difffirstzmin = -1;
      difffirstzmax = +1;
      difffirstzbins = (Int_t) ((difffirstzmax - difffirstzmin)/0.0001);

      difflastxmin = -1;
      difflastxmax = +1;
      difflastxbins = (Int_t) ((difflastxmax - difflastxmin)/0.0001);

      difflastymin = -1;
      difflastymax = +1;
      difflastybins = (Int_t) ((difflastymax - difflastymin)/0.0001);

      difflastzmin = -1;
      difflastzmax = +1;
      difflastzbins = (Int_t) ((difflastzmax - difflastzmin)/0.0001);

      diffptmin = -1;
      diffptmax = +1;
      diffptbins = (Int_t) ((diffptmax - diffptmin)/0.0001);
      
      diffpsimin = -1;
      diffpsimax = +1;
      diffpsibins = (Int_t) ((diffpsimax - diffpsimin)/0.0001);

      difftglmin = -1;
      difftglmax = +1;
      difftglbins = (Int_t) ((difftglmax - difftglmin)/0.0001);

      diffclustermin = -1;
      diffclustermax = +1;
      diffclusterbins = (Int_t) ((diffclustermax - diffclustermin)/0.0001);

#if 0
      diffpterrmin = -1;
      diffpterrmax = 1;
      diffpterrbins = (Int_t) ((diffpterrmax - diffpterrmin)/0.0001);

      diffpsierrmin = -1;
      diffpsierrmax = 1;
      diffpsierrbins = (Int_t) ((diffpsierrmax - diffpsierrmin)/0.0001);

      difftglerrmin = -1;
      difftglerrmax = 1;
      difftglerrbins = (Int_t) ((difftglerrmax - difftglerrmin)/0.0001);
#endif
    }

  // intialise histogramms
  TH1F* firstxhisto = new TH1F("Differences of first x (original - secondary) track", "Differences of first x (original - secondary) track", difffirstxbins, difffirstxmin, difffirstxmax);
  TH1F* firstyhisto = new TH1F("Differences of first y (original - secondary) track", "Differences of first y (original - secondary) track", difffirstybins, difffirstymin, difffirstymax);
 TH1F* firstzhisto = new TH1F("Differences of first z (original - secondary) track", "Differences of first z (original - secondary) track", difffirstzbins, difffirstzmin, difffirstzmax);
 TH1F* lastxhisto = new TH1F("Differences of last x (original - secondary) track", "Differences of last x (original - secondary) track", difflastxbins, difflastxmin, difflastxmax);
 TH1F* lastyhisto = new TH1F("Differences of last y (original - secondary) track", "Differences of last y (original - secondary) track", difflastybins, difflastymin, difflastymax);
 TH1F* lastzhisto = new TH1F("Differences of last z (original - secondary) track", "Differences of last z (original - secondary) track", difflastzbins, difflastzmin, difflastzmax);
 TH1F* pthisto = new TH1F("Differences of pt (original - secondary) track", "Differences of pt (original - secondary) track", diffptbins, diffptmin, diffptmax);
 TH1F* psihisto = new TH1F("Differences of psi (original - secondary) track", "Differences of psi (original - secondary) track", diffpsibins, diffpsimin, diffpsimax);
 TH1F* tglhisto = new TH1F("Differences of tgl (original - secondary) track", "Differences of tgl (original - secondary) track", difftglbins, difftglmin, difftglmax);
 TH1F* clusterhisto = new TH1F("Differences of asserted clusters (original - secondary) track", "Differences of asserted clusters (original - secondary) track", diffclusterbins, diffclustermin, diffclustermax);

#if 0
 // commented out since pterr is zero for every track!
 TH1F* pterrhisto = new TH1F("Differences of pt error (original - secondary) track", "Differences of pt error (original - secondary) track", diffpterrbins, diffpterrmin, diffpterrmax);
 TH1F* psierrhisto = new TH1F("Differences of psi error (original - secondary) track", "Differences of psi error (original - secondary) track", diffpsierrbins, diffpsierrmin, diffpsierrmax);
 TH1F* tglerrhisto = new TH1F("Differences of tgl error (original - secondary) track", "Differences of tgl error (original - secondary) track", difftglerrbins, difftglerrmin, difftglerrmax);

#endif

 // evaluate quality of fit:
 TH1I* matchinghisto = new TH1I("Matching indicator (5 - 10)", "Matching indicator (5 - 10)", 11, 0, 11);

  while(tracklistpointer != NULL)
    {
      // if currently processed track is matched, store differences of their parameters in histogram
      if(tracklistpointer->matchingindicator > 0)
	{
	  // matching track
	  trackmatchingpointer = tracklistpointer->matchingtrack;

	  if(relativedifferences == 1) // fill histogram with relative differences
	    {
	      firstxhisto->Fill((tracklistpointer->track.GetFirstPointX()-trackmatchingpointer->track.GetFirstPointX())/tracklistpointer->track.GetFirstPointX(),1);
	      firstyhisto->Fill((tracklistpointer->track.GetFirstPointY()-trackmatchingpointer->track.GetFirstPointY())/tracklistpointer->track.GetFirstPointY(),1);
	      firstzhisto->Fill((tracklistpointer->track.GetFirstPointZ()-trackmatchingpointer->track.GetFirstPointZ())/tracklistpointer->track.GetFirstPointZ(),1);
	      lastxhisto->Fill((tracklistpointer->track.GetLastPointX()-trackmatchingpointer->track.GetLastPointX())/tracklistpointer->track.GetLastPointX(),1);
	      lastyhisto->Fill((tracklistpointer->track.GetLastPointY()-trackmatchingpointer->track.GetLastPointY())/tracklistpointer->track.GetLastPointY(),1);
	      lastzhisto->Fill((tracklistpointer->track.GetLastPointZ()-trackmatchingpointer->track.GetLastPointZ())/tracklistpointer->track.GetLastPointZ(),1);
	      pthisto->Fill((tracklistpointer->track.GetPt()-trackmatchingpointer->track.GetPt())/tracklistpointer->track.GetPt(),1);
	      psihisto->Fill((tracklistpointer->track.GetPsi()-trackmatchingpointer->track.GetPsi())/tracklistpointer->track.GetPsi(),1);
	      tglhisto->Fill((tracklistpointer->track.GetTgl()-trackmatchingpointer->track.GetTgl())/tracklistpointer->track.GetTgl(),1);
	      clusterhisto->Fill((tracklistpointer->track.GetNHits()-trackmatchingpointer->track.GetNHits())/tracklistpointer->track.GetNHits(),1);

#if 0
	      pterrhisto->Fill((tracklistpointer->track.GetPterr()-trackmatchingpointer->track.GetPterr())/tracklistpointer->track.GetPterr(),1);
	      psierrhisto->Fill((tracklistpointer->track.GetPsierr()-trackmatchingpointer->track.GetPsierr())/tracklistpointer->track.GetPsierr(),1);
	      tglerrhisto->Fill((tracklistpointer->track.GetTglerr()-trackmatchingpointer->track.GetTglerr())/tracklistpointer->track.GetTglerr(),1);

	      HLTInfo("Pterr: 1st:%f  2nd:%f   value:%f",tracklistpointer->track.GetPterr(),trackmatchingpointer->track.GetPterr(),(tracklistpointer->track.GetPterr()-trackmatchingpointer->track.GetPterr())/tracklistpointer->track.GetPterr());
	      HLTInfo("Psierr: 1st:%f  2nd:%f   value:%f",tracklistpointer->track.GetPsierr(),trackmatchingpointer->track.GetPsierr(),(tracklistpointer->track.GetPsierr()-trackmatchingpointer->track.GetPsierr())/tracklistpointer->track.GetPsierr());
	      HLTInfo("Tglerr: 1st:%f  2nd:%f   value:%f",tracklistpointer->track.GetTglerr(),trackmatchingpointer->track.GetTglerr(),(tracklistpointer->track.GetTglerr()-trackmatchingpointer->track.GetTglerr())/tracklistpointer->track.GetTglerr());
#endif
	    }
	  else // otherwise fill histogram with absolute differences
	    {
	      firstxhisto->Fill(tracklistpointer->track.GetFirstPointX()-trackmatchingpointer->track.GetFirstPointX(),1);
	      firstyhisto->Fill(tracklistpointer->track.GetFirstPointY()-trackmatchingpointer->track.GetFirstPointY(),1);
	      firstzhisto->Fill(tracklistpointer->track.GetFirstPointZ()-trackmatchingpointer->track.GetFirstPointZ(),1);
	      lastxhisto->Fill(tracklistpointer->track.GetLastPointX()-trackmatchingpointer->track.GetLastPointX(),1);
	      lastyhisto->Fill(tracklistpointer->track.GetLastPointY()-trackmatchingpointer->track.GetLastPointY(),1);
	      lastzhisto->Fill(tracklistpointer->track.GetLastPointZ()-trackmatchingpointer->track.GetLastPointZ(),1);
	      pthisto->Fill(tracklistpointer->track.GetPt()-trackmatchingpointer->track.GetPt(),1);
	      psihisto->Fill(tracklistpointer->track.GetPsi()-trackmatchingpointer->track.GetPsi(),1);
	      tglhisto->Fill(tracklistpointer->track.GetTgl()-trackmatchingpointer->track.GetTgl(),1);
	      clusterhisto->Fill(tracklistpointer->track.GetNHits()-trackmatchingpointer->track.GetNHits(),1);
#if 0
	      // commented out since pterr is always zero for every track!
	      pterrhisto->Fill(tracklistpointer->track.GetPterr()-trackmatchingpointer->track.GetPterr(),1);
	      psierrhisto->Fill(tracklistpointer->track.GetPsierr()-trackmatchingpointer->track.GetPsierr(),1);
	      tglerrhisto->Fill(tracklistpointer->track.GetTglerr()-trackmatchingpointer->track.GetTglerr(),1);
#endif
	    }

	  // fill histogram that determines the quality of the fit
	  matchinghisto->Fill(tracklistpointer->matchingindicator,1);
	}

      tracklistpointer = tracklistpointer->next;
    }

  // write histograms to root file specified in command line argument -graphs <filename>.ROOT
  if(!fGraphFileName.IsNull())
    {
      TFile* graphrootfile = new TFile(fGraphFileName, "recreate");

       graphrootfile->WriteObject(firstxhisto,"firstxhistogram");
      //firstxhisto->Write();
      graphrootfile->WriteObject(firstyhisto,"firstyhistogram");
      //firstyhisto->Write();
      graphrootfile->WriteObject(firstzhisto,"firstzhistogram");
      //firstzhisto->Write();
      graphrootfile->WriteObject(lastxhisto,"lastxhistogram");
      //lastxhisto->Write();
      graphrootfile->WriteObject(lastyhisto,"lastyhistogram");
      //lastyhisto->Write();
      graphrootfile->WriteObject(lastzhisto,"lastzhistogram");
      //lastzhisto->Write();
      graphrootfile->WriteObject(pthisto,"pthistogram");
      //pthisto->Write();
      graphrootfile->WriteObject(psihisto,"psihistogram");
      //psihisto->Write();
      graphrootfile->WriteObject(tglhisto,"tglhistogram");
      //tglhisto->Write();
      graphrootfile->WriteObject(clusterhisto,"clusterhistogram");
      //clusterhisto->Write();
      graphrootfile->WriteObject(matchinghisto,"matchinghistogram");
      //matchinghisto->Write();
#if 0
      // errors in tracking (commented out since pterr is always 0 for every track!)
      graphrootfile->WriteObject(pterrhisto, "pterrhistogram");
      //pterrhisto->Write();
      graphrootfile->WriteObject(psierrhisto, "psierrhistogram");
      //psierrhisto->Write();
      graphrootfile->WriteObject(tglerrhisto, "tglerrhistogram");
      //tglerrhisto->Write();
#endif
      graphrootfile->Close();
    }
  else
    {
      HLTError("Error! No file for graphical output specified.");
      return EINVAL;
    }

  return 0;
}

Int_t AliHLTTPCCompModelAnalysis::DisplayModelResults()
{
  // see header file for class documentation

  //print out parameters of discarded track:
  AliHLTTPCTrackList* trackprintpointer = fTrackListPointer;
  AliHLTTPCTrackList* tracklistdeleter= trackprintpointer;

 FILE* dumpfile = NULL;

 if(!fDumpFileName.IsNull())
   {
     // open new file specified by command line argument
     dumpfile = fopen(fDumpFileName.Data(),"a");
     fprintf(dumpfile,"---------------MODEL ANALYSIS--------------- \n");
     fprintf(dumpfile,"---------------DISCARDED TRACKS: %d --------------- \n", fTrashTracks);
   };
 
 if(fTrackListPointer == NULL)
   {
     HLTInfo("No tracks discarded");
     HLTInfo("--------------");
   }
 else
   {
     
     HLTInfo("---------------DISCARDED TRACKS: %d ---------------", fTrashTracks);
     
     Int_t trashtrackcounter = 1;
     
     while(trackprintpointer != NULL)
       {
	 
	 // infos about found track
	 HLTInfo("%d : Discarding track with %d clusters.", trashtrackcounter, trackprintpointer->track.GetNHits());
	 //PYTHIA INFORMATION ABOUT PARTICLE ID
	 HLTInfo("Track parameters of discarded track:");	
	 HLTInfo("First x: %9.6f \t first y: %9.6f \t first z: %9.6f",trackprintpointer->track.GetFirstPointX(),trackprintpointer->track.GetFirstPointY(),trackprintpointer->track.GetFirstPointZ()); 
	 HLTInfo(" Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f", trackprintpointer->track.GetLastPointX(),trackprintpointer->track.GetLastPointY(),trackprintpointer->track.GetLastPointZ());
	 HLTInfo("     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f", trackprintpointer->track.GetPt(), trackprintpointer->track.GetPsi(), trackprintpointer->track.GetTgl());
	 
	 // write results to file if specified by command line argument
	 if(!fDumpFileName.IsNull())
	   {
	     fprintf(dumpfile,"%d : Discarding track with %d clusters. \n", trashtrackcounter, trackprintpointer->track.GetNHits());
	     fprintf(dumpfile,"Track parameters of discarded track: \n");	
	     fprintf(dumpfile,"First x: %9.6f \t first y: %9.6f \t first z: %9.6f \n",trackprintpointer->track.GetFirstPointX(),trackprintpointer->track.GetFirstPointY(),trackprintpointer->track.GetFirstPointZ()); 
	     fprintf(dumpfile," Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f \n", trackprintpointer->track.GetLastPointX(),trackprintpointer->track.GetLastPointY(),trackprintpointer->track.GetLastPointZ());
	     fprintf(dumpfile,"     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f \n", trackprintpointer->track.GetPt(), trackprintpointer->track.GetPsi(), trackprintpointer->track.GetTgl());
	   };
	 
	 // comparison with pythia information
	 if(trackprintpointer->wronglydiscarded)
	   {
	     HLTInfo("Found track has been wrongly discarded according to Pythia information.");
	     
	     // write results to file if specified by command line argument
	     if(!fDumpFileName.IsNull())
	       {
		  fprintf(dumpfile,"Found track has been wrongly discarded accoring to Pythia information. \n");
	       };
	     
	    };
	 
	 // found pt must be in range pt_pythia \pm 10% to be accepted
	 Double_t ptmin = trackprintpointer->pythiatrack.GetPt() - 0.1*trackprintpointer->pythiatrack.GetPt();
	 Double_t ptmax = trackprintpointer->pythiatrack.GetPt() + 0.1*trackprintpointer->pythiatrack.GetPt();
	 
	 if( (trackprintpointer->track.GetPt() < ptmin) ||(trackprintpointer->track.GetPt() > ptmax) )
	   {
	     HLTInfo("Pt of found track %f differs more than 10 %% from pt of pythia track %f.",trackprintpointer->track.GetPt(), trackprintpointer->pythiatrack.GetPt());
	     
	     if(!fDumpFileName.IsNull())
	       {
		 // write result to file if specified by command line argument
		 fprintf(dumpfile,"Pt of found track %f differs more than 10 %% from pt of pythia track %f. \n",trackprintpointer->track.GetPt(), trackprintpointer->pythiatrack.GetPt());
	       };
	     
	   };
	 
	 HLTInfo("--------------");
	 
	 if(!fDumpFileName.IsNull())
	   {
	     fprintf(dumpfile,"-------------- \n");
	   }; 
	 
	 
	  // go to next element in trash track list
	 tracklistdeleter = trackprintpointer;
	 trackprintpointer = trackprintpointer->next;
	 
	 ++trashtrackcounter;
	 
	 // free space 
	 delete tracklistdeleter;
	 tracklistdeleter = NULL;
	 
       } // end of while(trackpointer != NULL)

   } // end of else
 
 // print out number of noise clusters (not assigned to any track and not valid)
 HLTInfo("Number of discarded clusters not assigned to any track: %d", fTotalDiscardedClusters);
 HLTInfo("--------------");

 // write results to file if specified by command line argument
 if(!fDumpFileName.IsNull())
   {
     fprintf(dumpfile,"Number of discarded clusters not assigned to any track: %d \n", fTotalDiscardedClusters);
     fprintf(dumpfile,"-------------- \n");
   };
 
 // print out paramters of discarded valuable clusters
 HLTInfo("Number of discarded VALUABLE clusters: %d", fValuableDiscardedClusters);

 HLTInfo("--------------");
 
 if(!fDumpFileName.IsNull())
   {

     fprintf(dumpfile,"Number of discarded VALUABLE clusters: %d \n", fValuableDiscardedClusters);
     fclose(dumpfile);
   };
 
  return 0;
}

Int_t AliHLTTPCCompModelAnalysis::DisplayTrackResults()
{
  // see header file for class documentation
  // start with comparison
  if(CompareTracks() != 0)
    {
      return EINVAL;
    };

  // if dumptofile is activated, append results to output analysis file
  FILE* dumpfile = NULL;
  
  if(!fDumpFileName.IsNull())
    {
      // open new file specified by command line argument
      dumpfile = fopen(fDumpFileName.Data(),"a");
      
      fprintf(dumpfile,"---------------TRACK ANALYSIS--------------- \n");
      
    }
  
  // print out number of compared tracks
  HLTInfo("---------------ORIGINAL TRACKS: %d ---------------", fFirstTrackArray.GetNTracks());
  HLTInfo("Number of tracks with pt < 0.1 GeV: %d", fFirstTrashTracks);
  HLTInfo("Number of matched tracks with pt < 0.1 GeV: %d", fMatchedFirstTrashTracks);
  //PYTHIA INFORMATION ABOUT PARTICLE IDs
  HLTInfo("---------------2NDARY TRACKS: %d ---------------", fSecondTrackArray.GetNTracks());
  HLTInfo("Number of tracks with pt < 0.1 GeV: %d", fSecondTrashTracks);
  HLTInfo("Number of matched tracks with pt < 0.1 GeV: %d", fMatchedSecondTrashTracks);
  //PYTHIA INFORMATION ABOUT PARTICLE IDs
  HLTInfo("--------------");
  HLTInfo("Comparison of tracks within parameter precision: %f", fToleranceDeviation);
  HLTInfo("Number of compared tracks: %d", fTotalComparedTracks);
  HLTInfo("Number of unmatched original tracks: %d", fFirstUnmatchedTracks);
  HLTInfo("Number of unmatched secondary tracks: %d", fSecondUnmatchedTracks);
  HLTInfo("--------------");
  
  // print results to file
  if(!fDumpFileName.IsNull())
    {
      fprintf(dumpfile, "---------------%d ORIGINAL TRACKS---------------\n", fFirstTrackArray.GetNTracks());
      fprintf(dumpfile,"Number of tracks with pt < 0.1 GeV: %d \n", fFirstTrashTracks);
      fprintf(dumpfile,"Number of matched tracks with pt < 0.1 GeV: %d \n", fMatchedFirstTrashTracks);
      fprintf(dumpfile,"---------------%d 2NDARY TRACKS---------------\n", fSecondTrackArray.GetNTracks());
      fprintf(dumpfile,"Number of tracks with pt < 0.1 GeV: %d \n", fSecondTrashTracks);
      fprintf(dumpfile,"Number of matched tracks with pt < 0.1 GeV: %d \n", fMatchedSecondTrashTracks);
      fprintf(dumpfile,"--------------\n");
      fprintf(dumpfile,"Comparison of tracks within parameter precision: %f \n", fToleranceDeviation);
      fprintf(dumpfile,"Number of compared tracks: %d \n", fTotalComparedTracks); 
      fprintf(dumpfile,"Number of unmatched original tracks: %d \n", fFirstUnmatchedTracks); 
      fprintf(dumpfile,"Number of unmatched secondary tracks: %d \n", fSecondUnmatchedTracks); 
      fprintf(dumpfile,"--------------\n");
    }

  //////////////////////////////////////////////////////////////////////
  // additional files temporarily necessary for output information
  FILE* infofile = fopen("/afsuser/jwagner/TrackerTest_25092007/cosmics/fullanalysis08012008/trackingefficiency.out","a");
  FILE* info2file = fopen("/afsuser/jwagner/TrackerTest_25092007/cosmics/fullanalysis08012008/parameters.out", "a");

  // consistent = 0 if tracks and second tracks match perfectly, i.e. no doubly assigned tracks,etc.
  Int_t consistentfirst = fFirstTrackArray.GetNTracks() - fTotalComparedTracks - fFirstUnmatchedTracks;
  Int_t consistentsecond = fSecondTrackArray.GetNTracks() - fTotalComparedTracks - fSecondUnmatchedTracks;

  //fprintf(infofile, "1st tracks, 2nd tracks, compared, 1st uncompared, 2nd uncompared, cons.1st, cons.2nd \n");
  fprintf(info2file, " %d \t %d \t %d \t %d \t %d \t %d \t %d \n", fFirstTrackArray.GetNTracks(), fSecondTrackArray.GetNTracks(), fTotalComparedTracks, fFirstUnmatchedTracks, fSecondUnmatchedTracks, consistentfirst, consistentsecond);

  //fprintf(infofile, "1st trash tracks, 2nd trash tracks, 1st matched trash tracks, 2nd matched trash tracks \n");
  fprintf(infofile, "%d \t %d \t %d \t %d \n", fFirstTrashTracks, fSecondTrashTracks, fMatchedFirstTrashTracks, fMatchedSecondTrashTracks);

  fclose(infofile);
  fclose(info2file);
  ////////////////////////////////////////////////////////////////////

  // print out deviations
  Int_t tracknumber = 1;
  AliHLTTPCTrackList* listprintpointer = fFirstTrackList;
  while(listprintpointer != NULL)
    {
      // print out parameters of original track in comparison to secondary track:
      
      if(listprintpointer->matchingindicator != 0)
	{
	  
#if 0
	  HLTInfo("Track %d:", tracknumber);
	  HLTInfo("Original track matched to secondary track with matchingindicator %d", listprintpointer->matchingindicator);
	  
	  HLTInfo("Parameter comparison: Original vs. Secondary");

	  HLTInfo("Clusters: %d \t %d ",listprintpointer->track.GetNHits(), listprintpointer->matchingtrack->track.GetNHits());
	  HLTInfo("First x: %9.6f \t %9.6f", listprintpointer->track.GetFirstPointX(), listprintpointer->matchingtrack->track.GetFirstPointX());
	  HLTInfo("First y: %9.6f \t %9.6f", listprintpointer->track.GetFirstPointY(), listprintpointer->matchingtrack->track.GetFirstPointY());
	  HLTInfo("First z: %9.6f \t %9.6f", listprintpointer->track.GetFirstPointZ(), listprintpointer->matchingtrack->track.GetFirstPointZ()); 
	  HLTInfo("Last x: %9.6f \t %9.6f", listprintpointer->track.GetLastPointX(), listprintpointer->matchingtrack->track.GetLastPointX());
	  HLTInfo("Last y: %9.6f \t %9.6f", listprintpointer->track.GetLastPointY(), listprintpointer->matchingtrack->track.GetLastPointY());
	  HLTInfo("Last z: %9.6f \t %9.6f", listprintpointer->track.GetLastPointZ(), listprintpointer->matchingtrack->track.GetLastPointZ());
	  HLTInfo("    Pt: %9.6f \t %9.6f", listprintpointer->track.GetPt(), listprintpointer->matchingtrack->track.GetPt());  
	  HLTInfo("   Psi: %9.6f \t %9.6f", listprintpointer->track.GetPsi(), listprintpointer->matchingtrack->track.GetPsi());
	  HLTInfo("   Tgl: %9.6f \t %9.6f", listprintpointer->track.GetTgl(), listprintpointer->matchingtrack->track.GetTgl());

	  HLTInfo(" Pterr: %9.6f \t %9.6f", listprintpointer->track.GetPterr(), listprintpointer->matchingtrack->track.GetPterr());  
	  HLTInfo("Psierr: %9.6f \t %9.6f", listprintpointer->track.GetPsierr(), listprintpointer->matchingtrack->track.GetPsierr());
	  HLTInfo("Tglerr: %9.6f \t %9.6f", listprintpointer->track.GetTglerr(), listprintpointer->matchingtrack->track.GetTglerr());
	 
	  HLTInfo("--------------");
#endif
	  // print these results to file
	  if(!fDumpFileName.IsNull())
	    {
	      fprintf(dumpfile, "Track %d: \n", tracknumber);
	      fprintf(dumpfile, "Original track matched to secondary track with matchingindicator %d \n", listprintpointer->matchingindicator); 
	      fprintf(dumpfile, "Parameter comparison: Original vs. Secondary \n"); 
	      fprintf(dumpfile, "Clusters: %d \t %d \n ",listprintpointer->track.GetNHits(), listprintpointer->matchingtrack->track.GetNHits());
	      fprintf(dumpfile, "First x: %9.6f \t %9.6f \n", listprintpointer->track.GetFirstPointX(), listprintpointer->matchingtrack->track.GetFirstPointX());
	      fprintf(dumpfile, "First y: %9.6f \t %9.6f \n", listprintpointer->track.GetFirstPointY(), listprintpointer->matchingtrack->track.GetFirstPointY()); 
	      fprintf(dumpfile, "First z: %9.6f \t %9.6f \n", listprintpointer->track.GetFirstPointZ(), listprintpointer->matchingtrack->track.GetFirstPointZ()); 
	      fprintf(dumpfile, "Last x: %9.6f \t %9.6f \n", listprintpointer->track.GetLastPointX(), listprintpointer->matchingtrack->track.GetLastPointX());  
	      fprintf(dumpfile, "Last y: %9.6f \t %9.6f \n", listprintpointer->track.GetLastPointY(), listprintpointer->matchingtrack->track.GetLastPointY());
	      fprintf(dumpfile, "Last z: %9.6f \t %9.6f \n", listprintpointer->track.GetLastPointZ(), listprintpointer->matchingtrack->track.GetLastPointZ());
	      fprintf(dumpfile, "    Pt: %9.6f \t %9.6f \n", listprintpointer->track.GetPt(), listprintpointer->matchingtrack->track.GetPt());  
	      fprintf(dumpfile, "   Psi: %9.6f \t %9.6f \n", listprintpointer->track.GetPsi(), listprintpointer->matchingtrack->track.GetPsi()); 
	      fprintf(dumpfile, "   Tgl: %9.6f \t %9.6f \n", listprintpointer->track.GetTgl(), listprintpointer->matchingtrack->track.GetTgl());
	      fprintf(dumpfile, "--------------\n");
	    }
	  ++tracknumber;
	}
     
      listprintpointer = listprintpointer->next;
    }


  // print out not matched tracks from first track array:
  listprintpointer = fFirstTrackList;
  Int_t notmatchedtracknumber = 1;

  while(listprintpointer != NULL)
    {
      if(listprintpointer->matchingindicator == 0)
	{
#if 0	  
	  HLTInfo("Original Track, not matched with secondary track %d:", notmatchedtracknumber);
	  HLTInfo("Clusters: %d",listprintpointer->track.GetNHits());
	  //PYTHIA INFORMATION ABOUT PARTICLE ID
	  HLTInfo("First x: %9.6f \t first y: %9.6f \t first z: %9.6f",listprintpointer->track.GetFirstPointX(),listprintpointer->track.GetFirstPointY(),listprintpointer->track.GetFirstPointZ()); 
	  HLTInfo(" Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f", listprintpointer->track.GetLastPointX(),listprintpointer->track.GetLastPointY(),listprintpointer->track.GetLastPointZ());
	  HLTInfo("     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f", listprintpointer->track.GetPt(), listprintpointer->track.GetPsi(), listprintpointer->track.GetTgl());
	  HLTInfo("--------------");
#endif
 
	  // print these results to file
	  if(!fDumpFileName.IsNull())
	    {
	      fprintf(dumpfile, "Original Track, not matched with secondary track %d: \n", notmatchedtracknumber); 
	      fprintf(dumpfile, "Clusters: %d \n",listprintpointer->track.GetNHits());
	      fprintf(dumpfile, "First x: %9.6f \t first y: %9.6f \t first z: %9.6f \n",listprintpointer->track.GetFirstPointX(),listprintpointer->track.GetFirstPointY(),listprintpointer->track.GetFirstPointZ());
	      fprintf(dumpfile, " Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f \n", listprintpointer->track.GetLastPointX(),listprintpointer->track.GetLastPointY(),listprintpointer->track.GetLastPointZ());
	      fprintf(dumpfile, "     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f \n", listprintpointer->track.GetPt(), listprintpointer->track.GetPsi(), listprintpointer->track.GetTgl()); 
	      fprintf(dumpfile, "--------------\n");
	    }

	  ++notmatchedtracknumber;
	}

      listprintpointer = listprintpointer->next;
    }

  // print out not matched tracks from second track array:
  listprintpointer = fSecondTrackList;
  notmatchedtracknumber = 1;

  while(listprintpointer != NULL)
    {
      if(listprintpointer->matchingindicator == 0)
	{
#if 0	  
	  HLTInfo("Secondary Track, not matched with original track %d:", notmatchedtracknumber);
	  HLTInfo("Clusters: %d",listprintpointer->track.GetNHits());
	  //PYTHIA INFORMATION ABOUT PARTICLE ID	
	  HLTInfo("First x: %9.6f \t first y: %9.6f \t first z: %9.6f",listprintpointer->track.GetFirstPointX(),listprintpointer->track.GetFirstPointY(),listprintpointer->track.GetFirstPointZ()); 
	  HLTInfo(" Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f", listprintpointer->track.GetLastPointX(),listprintpointer->track.GetLastPointY(),listprintpointer->track.GetLastPointZ());
	  HLTInfo("     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f", listprintpointer->track.GetPt(), listprintpointer->track.GetPsi(), listprintpointer->track.GetTgl());
	  HLTInfo("--------------");
#endif 
	  // print these results to file
	  if(!fDumpFileName.IsNull())
	    {
	      fprintf(dumpfile, "Secondary Track, not matched with original track %d: \n", notmatchedtracknumber);
	      fprintf(dumpfile, "Clusters: %d \n",listprintpointer->track.GetNHits());
	      fprintf(dumpfile, "First x: %9.6f \t first y: %9.6f \t first z: %9.6f \n",listprintpointer->track.GetFirstPointX(),listprintpointer->track.GetFirstPointY(),listprintpointer->track.GetFirstPointZ()); 
	      fprintf(dumpfile, " Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f \n", listprintpointer->track.GetLastPointX(),listprintpointer->track.GetLastPointY(),listprintpointer->track.GetLastPointZ());
	      fprintf(dumpfile, "     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f \n", listprintpointer->track.GetPt(), listprintpointer->track.GetPsi(), listprintpointer->track.GetTgl());
	      fprintf(dumpfile, "--------------\n");
	    }

	  ++notmatchedtracknumber;
	}

      listprintpointer = listprintpointer->next;
    }

  // close output analysis file
 if(!fDumpFileName.IsNull())
   {
     fclose(dumpfile);
   };

 // if results should be written to graphical output:
 if(!fGraphFileName.IsNull())
   {
     CreateGraphs(1); // specifiy if absolute or rel. differences should be saved (CreateGraphs(0) or CreateGraphs()/ CreateGraphs(1)
   };
 
  // free reserved space
  
  return 0;
}
