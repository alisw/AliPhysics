// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: J. Wagner <jwagner@cern.ch>                           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

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
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

/** destructor **/
AliHLTTPCCompModelAnalysis::~AliHLTTPCCompModelAnalysis()
    {
  // see header file for class documentation
    for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	    {
	    if ( fDiscardedClusters[slice][patch]!=NULL )
		{
		delete [] fDiscardedClusters[slice][patch];
		fDiscardedClusters[slice][patch]=NULL;
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

Int_t AliHLTTPCCompModelAnalysis::SetClusters(AliHLTTPCClusterData* clusters, UInt_t slice, UInt_t patch, Bool_t fillingfirstclusters)
{
  // see header file for class documentation
  if(fillingfirstclusters == 1)
    {
      if ( slice>=36 || patch>=6 )
	return EINVAL;
      if ( fOriginalClusters[slice][patch] )
	return EBUSY;
      fOriginalClusters[slice][patch] = clusters;

      //HLTDebug( "Filling %u clusters in first array", (unsigned)clusters->fSpacePointCnt);
    }
  else
    {
      if ( slice>=36 || patch>=6 )
	return EINVAL;
      if ( fSecondaryClusters[slice][patch] )
	return EBUSY;
      fSecondaryClusters[slice][patch] = clusters;

      //HLTDebug( "Filling %u clusters in second array", (unsigned)clusters->fSpacePointCnt);
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
      return -EINVAL;
    };
 

  if(secondtracks == 0)
    {
      HLTError("No tracks in second track array!");
      return -EINVAL;
    };

  // take track from first tracking,
  for(Int_t ii=0; ii < firsttracks; ii++)
    {
      if (fFirstTrackArray.GetCheckedTrack(ii)==NULL) continue;

      // build track list for all tracks in first array
      // FIXME: I can't find the cleanup of the linked list fFirstTrackList
      AliHLTTPCTrackList* currenttrackentry = new AliHLTTPCTrackList;
      if (!currenttrackentry) return -ENOMEM;
      currenttrackentry->fTrack = *(fFirstTrackArray.GetCheckedTrack(ii));

      // get its pythia information, 
      currenttrackentry->fPythiatrack = GetComparableTrackPythiaInfo(currenttrackentry->fTrack);

      currenttrackentry->fMatchingindicator = 0;

      // put this element as first in list
      currenttrackentry->fNext = fFirstTrackList;
      fFirstTrackList = currenttrackentry;

      // count tracks below 0.1GeV
      if(currenttrackentry->fTrack.GetPt()<0.1)
	{
	  ++fFirstTrashTracks;
	}

    }

 // take track from second tracking,
  for(Int_t ii=0; ii < secondtracks; ii++)
    {
      // build track list for all tracks in second array
      // FIXME: I can't find the cleanup of the linked list fSecondTrackArray
      AliHLTTPCTrackList* currenttrackentry = new AliHLTTPCTrackList;
      if (!currenttrackentry) return -ENOMEM;
      currenttrackentry->fTrack = *(fSecondTrackArray.GetCheckedTrack(ii));

      // get its pythia information, 
      currenttrackentry->fPythiatrack = GetComparableTrackPythiaInfo(currenttrackentry->fTrack);

      // put this element as first in list
      currenttrackentry->fNext = fSecondTrackList;
      fSecondTrackList = currenttrackentry;

      // count tracks below 0.1GeV
      if(currenttrackentry->fTrack.GetPt()<0.1)
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
    
	     if((CompareTrackInfo(firstmatchpointer, secondmatchpointer) > firstmatchpointer->fMatchingindicator))
	       {
		 // look if current better matching track has already been matched before
		 if((CompareTrackInfo(firstmatchpointer, secondmatchpointer) > secondmatchpointer->fMatchingindicator))
		   {
		     
		     // set previously assigned matchingindicator (if there was one) of secondary track back to zero
		     if(firstmatchpointer->fMatchingindicator > 0)
		       {
			 firstmatchpointer->fMatchingtrack->fMatchingindicator = 0;
			 firstmatchpointer->fMatchingtrack->fMatchingtrack = NULL;
		       }

		     if(secondmatchpointer->fMatchingindicator > 0)
		       {
			 secondmatchpointer->fMatchingtrack->fMatchingindicator = 0;
			 secondmatchpointer->fMatchingtrack->fMatchingtrack = NULL;
		       }
		     
		     // compare according to tracks themselves (other possibility: compare pythiatracks - better!)
		     secondmatchpointer->fMatchingindicator = CompareTrackInfo(firstmatchpointer, secondmatchpointer) ;
		     firstmatchpointer->fMatchingindicator = CompareTrackInfo(firstmatchpointer, secondmatchpointer);
		     
		     // remember which track matches which
		     secondmatchpointer->fMatchingtrack = firstmatchpointer;
		     firstmatchpointer->fMatchingtrack = secondmatchpointer;

		   } // end if compare > second matching indicator

	       } // end if compare > first matching indicator
	     
	   }// end if compare > 4
	 	 
	 secondmatchpointer = secondmatchpointer->fNext;
       }
     
     // go on with next original track
     firstmatchpointer = firstmatchpointer->fNext;
   }
  
  // count not matched tracks in first and second track list
  AliHLTTPCTrackList* nomatchcounter = fFirstTrackList;

  while(nomatchcounter != NULL)
    {
      if(nomatchcounter->fMatchingindicator == 0)
	{
	  ++fFirstUnmatchedTracks;
	}
      else
	{
	 ++fTotalComparedTracks;

	 // count matched trash tracks
	 if(nomatchcounter->fTrack.GetPt() < 0.1)
	   {
	     ++fMatchedFirstTrashTracks;
	   };
	}
      nomatchcounter = nomatchcounter->fNext;
    }

  nomatchcounter = fSecondTrackList;
  while(nomatchcounter != NULL)
    {
      if(nomatchcounter->fMatchingindicator == 0)
	{
	  ++fSecondUnmatchedTracks;
	}
      else
	{
	  // count matched trash tracks
	 if(nomatchcounter->fTrack.GetPt() < 0.1)
	   {
	     ++fMatchedSecondTrashTracks;
	   };
	}
    
      nomatchcounter = nomatchcounter->fNext;
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

Int_t AliHLTTPCCompModelAnalysis::CompareClusters(Bool_t relativedifferences)
{

  Int_t totaloriginal = 0;
  Int_t totalsecondary = 0;
  Int_t usedclusters = 0;
  Int_t comparedclusters = 0;
  Int_t notcomparedclusters = 0;

  // create graphs out of differences and leave loop
  TFile* clustergraphrootfile = NULL;
  if(!fGraphFileName.IsNull())
    {
      clustergraphrootfile = new TFile(fGraphFileName, "recreate");
    }

  // specifications of histograms
  Double_t clusterdifffxmin, clusterdifffxmax;
  Double_t clusterdifffymin, clusterdifffymax;
  Double_t clusterdifffzmin, clusterdifffzmax;
  Int_t clusterdifffxbins, clusterdifffybins, clusterdifffzbins;

  Double_t clusterdifffsigmay2min, clusterdifffsigmay2max;
  Double_t clusterdifffsigmaz2min, clusterdifffsigmaz2max;
  Int_t clusterdifffsigmay2bins, clusterdifffsigmaz2bins;

  if (!relativedifferences) // not tested yet!
    {
      clusterdifffxmin = -1;
      clusterdifffxmax = +1;
      clusterdifffxbins = (Int_t) ((clusterdifffxmax - clusterdifffxmin)/0.0001);
  
      clusterdifffymin = -1;
      clusterdifffymax = +1;
      clusterdifffybins = (Int_t) ((clusterdifffymax - clusterdifffymin)/0.0001);
  
      clusterdifffzmin = -1;
      clusterdifffzmax = +1;
      clusterdifffzbins = (Int_t) ((clusterdifffzmax - clusterdifffzmin)/0.0001);
  
      clusterdifffsigmay2min = -1;
      clusterdifffsigmay2max = +1;
      clusterdifffsigmay2bins = (Int_t) ((clusterdifffsigmay2max - clusterdifffsigmay2min)/0.0001);
  
      clusterdifffsigmaz2min = -1;
      clusterdifffsigmaz2max = +1;
      clusterdifffsigmaz2bins = (Int_t) ((clusterdifffsigmaz2max - clusterdifffsigmaz2min)/0.0001);
    }
  else
    {
      clusterdifffxmin = -1;
      clusterdifffxmax = +1;
      clusterdifffxbins = (Int_t) ((clusterdifffxmax - clusterdifffxmin)/0.0001);
  
      clusterdifffymin = -1;
      clusterdifffymax = +1;
      clusterdifffybins = (Int_t) ((clusterdifffymax - clusterdifffymin)/0.0001);
  
      clusterdifffzmin = -1;
      clusterdifffzmax = +1;
      clusterdifffzbins = (Int_t) ((clusterdifffzmax - clusterdifffzmin)/0.0001);
  
      clusterdifffsigmay2min = -1;
      clusterdifffsigmay2max = +1;
      clusterdifffsigmay2bins = (Int_t) ((clusterdifffsigmay2max - clusterdifffsigmay2min)/0.0001);
  
      clusterdifffsigmaz2min = -1;
      clusterdifffsigmaz2max = +1;
      clusterdifffsigmaz2bins = (Int_t) ((clusterdifffsigmaz2max - clusterdifffsigmaz2min)/0.0001);
    }
  
  // intialise histogramms
  TH1F* clusterfxhisto = new TH1F("Differences of x (original - secondary) clusters", "Differences of x (original - secondary) clusters", clusterdifffxbins, clusterdifffxmin, clusterdifffxmax);
  TH1F* clusterfyhisto = new TH1F("Differences of y (original - secondary) clusters", "Differences of y (original - secondary) clusters", clusterdifffybins, clusterdifffymin, clusterdifffymax);
  TH1F* clusterfzhisto = new TH1F("Differences of z (original - secondary) clusters", "Differences of z (original - secondary) clusters", clusterdifffzbins, clusterdifffzmin, clusterdifffzmax);
  TH1F* clusterfsigmay2histo = new TH1F("Differences of sigmay2 (original - secondary) clusters", "Differences of sigmay2 (original - secondary) clusters", clusterdifffsigmay2bins, clusterdifffsigmay2min, clusterdifffsigmay2max);
  TH1F* clusterfsigmaz2histo = new TH1F("Differences of sigmaz2 (original - secondary) clusters", "Differences of sigmaz2 (original - secondary) clusters", clusterdifffsigmaz2bins, clusterdifffsigmaz2min, clusterdifffsigmaz2max);
  

  // see headerfile for class documentation
  // compare for each slice and patch the clusters in the original cluster array
  // to the ones of the secondary cluster array
  for(Int_t slicecntr = 0; slicecntr < 36; slicecntr++)
    {
      for (Int_t patchcntr = 0; patchcntr < 6; patchcntr++)
	{
	  if(!fOriginalClusters[slicecntr][patchcntr])
	    {
	      // HLTDebug("No original clusters for slice %d patch %d", slicecntr, patchcntr);
	      continue;
	    }

	  if(!fSecondaryClusters[slicecntr][patchcntr])
	    {
	      //HLTDebug("No secondary clusters for slice %d patch %d", slicecntr, patchcntr);
	      continue;
	    }

	  for ( unsigned long ii=0; ii<fOriginalClusters[slicecntr][patchcntr]->fSpacePointCnt; ii++ )
	    {
	      ++totaloriginal;

	      // search matching secondary cluster by charge and padrow, 
	      // fill histograms if cluster has been used in tracking process
	      Int_t found = 0;
	     
	      for(unsigned long jj=0; jj<fSecondaryClusters[slicecntr][patchcntr]->fSpacePointCnt; jj++)
		{
		  // if fTrackN != -1 -> draw histograms out of used clusters
		  // if fTrackN == -1 -> draw histograms out of unused clusters
		  if (fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fTrackN == -1)
		    {
		      if((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fCharge == fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fCharge) && (fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fPadRow == fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fPadRow) )
			{
			  
			  if (relativedifferences == 1)
			    {
			      // check whether first entries in cluster array are zero
			      if(fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fX == 0)
				{
				  HLTWarning("Warning! x value of original cluster is zero, relative differences cannot be calculated!");
				};

			      if(fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fY == 0)
				{
				  HLTWarning("Warning! y value of original cluster is zero, relative differences cannot be calculated!");
				};

			      if(fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fZ == 0)
				{
				  HLTWarning("Warning! z value of original cluster is zero, relative differences cannot be calculated!");
				};

			      if(fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaY2 == 0)
				{
				  HLTWarning("Warning! sigmay2 value of original cluster is zero, relative differences cannot be calculated!");
				};

			      if(fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaZ2 == 0)
				{
				  HLTWarning("Warning! sigmaz2 value of original cluster is zero, relative differences cannot be calculated!");
				};

			      // fill relative differences in histograms
			      clusterfxhisto->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fX - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fX)/fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fX,1);
			      clusterfyhisto->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fY - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fY)/fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fY,1);
			      clusterfzhisto->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fZ - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fZ)/fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fZ,1);
			      clusterfsigmay2histo->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaY2 - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fSigmaY2)/fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaY2,1);
			      clusterfsigmaz2histo->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaZ2 - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fSigmaZ2)/fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaZ2,1);
			    }
			  else
			    {
			      // fill absolute differences histograms
			      clusterfxhisto->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fX - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fX),1);
			      clusterfyhisto->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fY - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fY),1);
			      clusterfzhisto->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fZ - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fZ),1);
			      clusterfsigmay2histo->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaY2 - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fSigmaY2),1);
			      clusterfsigmaz2histo->Fill((fOriginalClusters[slicecntr][patchcntr]->fSpacePoints[ii].fSigmaZ2 - fSecondaryClusters[slicecntr][patchcntr]->fSpacePoints[jj].fSigmaZ2),1);
			    }
			      
			  found = 1;
			  ++comparedclusters;
			  break;
			}
		    }
		}
	      
	      if(found == 0)
		{
		  ++notcomparedclusters;
		}

	    }
	}
    }

  // write graphs to rootfile
 if(!fGraphFileName.IsNull())
    {
      clustergraphrootfile->WriteObject(clusterfxhisto, "clusterfxhistogram");
      clustergraphrootfile->WriteObject(clusterfyhisto, "clusterfyhistogram");
      clustergraphrootfile->WriteObject(clusterfzhisto, "clusterfzhistogram");
      clustergraphrootfile->WriteObject(clusterfsigmay2histo, "clusterfsigmay2histogram");
      clustergraphrootfile->WriteObject(clusterfsigmaz2histo, "clusterfsigmaz2histogram");
      clustergraphrootfile->Close();
    }
  
  // count clusters used for tracking
  for (Int_t slicecount=0; slicecount<36; slicecount++)
    {
      for(Int_t patchcount=0; patchcount<6; patchcount++)
	{

	  if(!fSecondaryClusters[slicecount][patchcount])
	    {
	      //HLTDebug("No secondary clusters for slice %d patch %d", slicecntr, patchcntr);
	      continue;
	    }

	  for(Int_t count=0; count < (Int_t) fSecondaryClusters[slicecount][patchcount]->fSpacePointCnt; count++)
	    {

	      ++totalsecondary;

	      if(fSecondaryClusters[slicecount][patchcount]->fSpacePoints[count].fTrackN != -1)
		{
		  ++usedclusters;
		};
	    }
	}
    }
  
  // Display results of cluster analysis
  HLTInfo("Number of original clusters: %d", totaloriginal);
  HLTInfo("Number of 2ndary clusters: %d", totalsecondary);
  HLTInfo("Number of 2ndary clusters used for tracking: %d", usedclusters);
  HLTInfo("Number of compared (tracked) original clusters: %d", comparedclusters);
  HLTInfo("Number of uncompared (tracked) original clusters: %d", notcomparedclusters);

  //////////////////////////////////////////////////////
  FILE* clusterfile = fopen("/afsuser/jwagner/TrackerTest_25092007/cosmics/fullanalysis/clusteranalysis.out", "a");

  fprintf(clusterfile, "%d \t %d \t %d \t %d \t %d \t %d \n", totaloriginal, totalsecondary, usedclusters, comparedclusters, notcomparedclusters, totaloriginal-comparedclusters-notcomparedclusters);

  fclose(clusterfile);
  ////////////////////////////////////////////////////////

  // consistency check
  if(comparedclusters + notcomparedclusters != totaloriginal)
    {
      HLTWarning("Possible inconsistency: number of compared clusters %d + number of not compared clusters %d not equal to total number of original clusters %d", comparedclusters, notcomparedclusters, totaloriginal);
    }

  return 0;
}

AliHLTTPCTrack AliHLTTPCCompModelAnalysis::GetComparableTrackPythiaInfo(const AliHLTTPCTrack& comparabletrack) const
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
  tracklistentry->fTrack = *lowpttrack;
  tracklistentry->fWronglydiscarded = GetTrashTrackPythiaInfo(lowpttrack);
  tracklistentry->fMatchingindicator = 0; // not needed here, therefore initialised to zero
  tracklistentry->fNext = fTrackListPointer;
  tracklistentry->fMatchingtrack = NULL; // not needed here, therefore initialised to NULL
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

Bool_t AliHLTTPCCompModelAnalysis::GetTrashTrackPythiaInfo(const AliHLTTPCTrack* /*discardedtrack*/ ) const
{
  // see header file for class documentation
  // store information from pythia in current track list entry
  // fTrackListPointer.pythiatrack = FillFromPythia...

  return 0;
}

Bool_t AliHLTTPCCompModelAnalysis::GetClusterPythiaInfo(const AliHLTTPCClusterData* /*discardedcluster*/) const
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
  if(abs((firsttracklistelement->fTrack.GetFirstPointX() - secondtracklistelement->fTrack.GetFirstPointX()))/firsttracklistelement->fTrack.GetFirstPointX() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->fTrack.GetFirstPointY() - secondtracklistelement->fTrack.GetFirstPointY()))/firsttracklistelement->fTrack.GetFirstPointY() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->fTrack.GetFirstPointZ() - secondtracklistelement->fTrack.GetFirstPointZ()))/firsttracklistelement->fTrack.GetFirstPointZ() <= fToleranceDeviation)
    ++currentmatchingindicator;
  
  // compare end point
  if(abs((firsttracklistelement->fTrack.GetLastPointX() - secondtracklistelement->fTrack.GetLastPointX()))/firsttracklistelement->fTrack.GetLastPointX() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->fTrack.GetLastPointY() - secondtracklistelement->fTrack.GetLastPointY()))/firsttracklistelement->fTrack.GetLastPointY() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->fTrack.GetLastPointZ() - secondtracklistelement->fTrack.GetLastPointZ()))/firsttracklistelement->fTrack.GetLastPointZ() <= fToleranceDeviation)
    ++currentmatchingindicator;

  // compare pt, psi, tgl
  if(abs((firsttracklistelement->fTrack.GetPt() - secondtracklistelement->fTrack.GetPt()))/firsttracklistelement->fTrack.GetPt() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->fTrack.GetPsi() - secondtracklistelement->fTrack.GetPsi()))/firsttracklistelement->fTrack.GetPsi() <= fToleranceDeviation)
    ++currentmatchingindicator;

  if(abs((firsttracklistelement->fTrack.GetTgl() - secondtracklistelement->fTrack.GetTgl()))/firsttracklistelement->fTrack.GetTgl() <= fToleranceDeviation)
    ++currentmatchingindicator;

  // compare number of assigned cluster hits
  if(abs((firsttracklistelement->fTrack.GetNHits() - secondtracklistelement->fTrack.GetNHits()))/firsttracklistelement->fTrack.GetNHits() <= fToleranceDeviation)
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
  if(firsttracklistelement->fPythiatrack.GetFirstPointX() == secondtracklistelement->fPythiatrack.GetFirstPointX())
    ++currentmatchingindicator;

  if(firsttracklistelement->fPythiatrack.GetFirstPointY() == secondtracklistelement->fPythiatrack.GetFirstPointY())
    ++currentmatchingindicator;

  if(firsttracklistelement->fPythiatrack.GetFirstPointZ() == secondtracklistelement->fPythiatrack.GetFirstPointZ())
    ++currentmatchingindicator;
  
  // compare end point
  if(firsttracklistelement->fPythiatrack.GetLastPointX() == secondtracklistelement->fPythiatrack.GetLastPointX())
    ++currentmatchingindicator;

  if(firsttracklistelement->fPythiatrack.GetLastPointY() == secondtracklistelement->fPythiatrack.GetLastPointY())
    ++currentmatchingindicator;

  if(firsttracklistelement->fPythiatrack.GetLastPointZ() == secondtracklistelement->fPythiatrack.GetLastPointZ())
    ++currentmatchingindicator;

  // compare pt, psi, tgl
  if(firsttracklistelement->fPythiatrack.GetPt() == secondtracklistelement->fPythiatrack.GetPt())
    ++currentmatchingindicator;

  if(firsttracklistelement->fPythiatrack.GetPsi() == secondtracklistelement->fPythiatrack.GetPsi())
    ++currentmatchingindicator;

  if(firsttracklistelement->fPythiatrack.GetTgl() == secondtracklistelement->fPythiatrack.GetTgl())
    ++currentmatchingindicator;

  // compare number of assigned cluster hits
  if(firsttracklistelement->fPythiatrack.GetNHits() == secondtracklistelement->fPythiatrack.GetNHits())
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

      //#if 0
      diffpterrmin = -1;
      diffpterrmax = 1;
      diffpterrbins = (Int_t) ((diffpterrmax - diffpterrmin)/1);

      diffpsierrmin = -1;
      diffpsierrmax = 1;
      diffpsierrbins = (Int_t) ((diffpsierrmax - diffpsierrmin)/1);

      difftglerrmin = -1;
      difftglerrmax = 1;
      difftglerrbins = (Int_t) ((difftglerrmax - difftglerrmin)/1);
      //#endif

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

      //#if 0
      diffpterrmin = -1;
      diffpterrmax = 1;
      diffpterrbins = (Int_t) ((diffpterrmax - diffpterrmin)/0.0001);

      diffpsierrmin = -1;
      diffpsierrmax = 1;
      diffpsierrbins = (Int_t) ((diffpsierrmax - diffpsierrmin)/0.0001);

      difftglerrmin = -1;
      difftglerrmax = 1;
      difftglerrbins = (Int_t) ((difftglerrmax - difftglerrmin)/0.0001);
      //#endif
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

 //#if 0
 // commented out since pterr is zero for every track!
 TH1F* pterrhisto = new TH1F("Differences of pt error (original - secondary) track", "Differences of pt error (original - secondary) track", diffpterrbins, diffpterrmin, diffpterrmax);
 TH1F* psierrhisto = new TH1F("Differences of psi error (original - secondary) track", "Differences of psi error (original - secondary) track", diffpsierrbins, diffpsierrmin, diffpsierrmax);
 TH1F* tglerrhisto = new TH1F("Differences of tgl error (original - secondary) track", "Differences of tgl error (original - secondary) track", difftglerrbins, difftglerrmin, difftglerrmax);

 //#endif

 // initialise histograms for tracking efficiency against pt
 // relative p_t resolution: 1.2% -> take 0.1GeV * 1.2 % --> binsize 0.001 sufficient to grant correct resolution for low pt
 TH1F* firsttracks = new TH1F("pt occurrence original", "Occurrence of pt in original tracks", 10000, 0, 10);
 TH1F* matchedtrackeff = new TH1F("matchedtreffeff", "Occurrence of 2ndary tracks with good pt", 10000, 0, 10);
 TH1F* trackeff = new TH1F("tracking efficiency vs. pt", "Tracking efficiency vs. pt", 10000, 0, 10);

 // evaluate quality of fit:
 TH1I* matchinghisto = new TH1I("Matching indicator (5 - 10)", "Matching indicator (5 - 10)", 11, 0, 11);

  while(tracklistpointer != NULL)
    {
      // if currently processed track is matched, store differences of their parameters in histogram
      if(tracklistpointer->fMatchingindicator > 0)
	{
	  // matching track
	  trackmatchingpointer = tracklistpointer->fMatchingtrack;

	  // fill histograms for trackingefficiency vs. pt
	  firsttracks->Fill(tracklistpointer->fTrack.GetPt(),1);
	  
	  if(abs(tracklistpointer->fTrack.GetPt()-trackmatchingpointer->fTrack.GetPt()) < 0.012*tracklistpointer->fTrack.GetPt())
	    {
	      matchedtrackeff->Fill(tracklistpointer->fTrack.GetPt(),1);
	    }

	  //tracklistpointer = tracklistpointer->fNext; // only efficiency is considered...
	  //continue; // only efficiency is considered...

	  if(relativedifferences == 1) // fill histogram with relative differences
	    {

	      // check if first track parameters are not zero!
	      if (tracklistpointer->fTrack.GetFirstPointX()==0)
		{
		  HLTWarning("Warning! First x of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetFirstPointY()==0)
		{
		  HLTWarning("Warning! First y of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetFirstPointZ()==0)
		{
		  HLTWarning("Warning! First z of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetLastPointX()==0)
		{
		  HLTWarning("Warning! Last x of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetLastPointY()==0)
		{
		  HLTWarning("Warning! Last y of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetLastPointZ()==0)
		{
		  HLTWarning("Warning! Last z of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetPt()==0)
		{
		  HLTWarning("Warning! Pt of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetPsi()==0)
		{
		  HLTWarning("Warning! Psi of original track is zero, relative differences cannot be calculated!");
		};
	      if (tracklistpointer->fTrack.GetTgl()==0)
		{
		  HLTWarning("Warning! Tgl of original track is zero, relative differences cannot be calculated!");
		};
	      firstxhisto->Fill((tracklistpointer->fTrack.GetFirstPointX()-trackmatchingpointer->fTrack.GetFirstPointX())/tracklistpointer->fTrack.GetFirstPointX(),1);
	      firstyhisto->Fill((tracklistpointer->fTrack.GetFirstPointY()-trackmatchingpointer->fTrack.GetFirstPointY())/tracklistpointer->fTrack.GetFirstPointY(),1);
	      firstzhisto->Fill((tracklistpointer->fTrack.GetFirstPointZ()-trackmatchingpointer->fTrack.GetFirstPointZ())/tracklistpointer->fTrack.GetFirstPointZ(),1);
	      lastxhisto->Fill((tracklistpointer->fTrack.GetLastPointX()-trackmatchingpointer->fTrack.GetLastPointX())/tracklistpointer->fTrack.GetLastPointX(),1);
	      lastyhisto->Fill((tracklistpointer->fTrack.GetLastPointY()-trackmatchingpointer->fTrack.GetLastPointY())/tracklistpointer->fTrack.GetLastPointY(),1);
	      lastzhisto->Fill((tracklistpointer->fTrack.GetLastPointZ()-trackmatchingpointer->fTrack.GetLastPointZ())/tracklistpointer->fTrack.GetLastPointZ(),1);
	      pthisto->Fill((tracklistpointer->fTrack.GetPt()-trackmatchingpointer->fTrack.GetPt())/tracklistpointer->fTrack.GetPt(),1);
	      psihisto->Fill((tracklistpointer->fTrack.GetPsi()-trackmatchingpointer->fTrack.GetPsi())/tracklistpointer->fTrack.GetPsi(),1);
	      tglhisto->Fill((tracklistpointer->fTrack.GetTgl()-trackmatchingpointer->fTrack.GetTgl())/tracklistpointer->fTrack.GetTgl(),1);
	      clusterhisto->Fill((tracklistpointer->fTrack.GetNHits()-trackmatchingpointer->fTrack.GetNHits())/tracklistpointer->fTrack.GetNHits(),1);

	      //#if 0
	      pterrhisto->Fill((tracklistpointer->fTrack.GetPterr()-trackmatchingpointer->fTrack.GetPterr())/tracklistpointer->fTrack.GetPterr(),1);
	      psierrhisto->Fill((tracklistpointer->fTrack.GetPsierr()-trackmatchingpointer->fTrack.GetPsierr())/tracklistpointer->fTrack.GetPsierr(),1);
	      tglerrhisto->Fill((tracklistpointer->fTrack.GetTglerr()-trackmatchingpointer->fTrack.GetTglerr())/tracklistpointer->fTrack.GetTglerr(),1);

	      //HLTInfo("Pterr: 1st:%f  2nd:%f   value:%f",tracklistpointer->fTrack.GetPterr(),trackmatchingpointer->fTrack.GetPterr(),(tracklistpointer->fTrack.GetPterr()-trackmatchingpointer->fTrack.GetPterr())/tracklistpointer->fTrack.GetPterr());
	      //HLTInfo("Psierr: 1st:%f  2nd:%f   value:%f",tracklistpointer->fTrack.GetPsierr(),trackmatchingpointer->fTrack.GetPsierr(),(tracklistpointer->fTrack.GetPsierr()-trackmatchingpointer->fTrack.GetPsierr())/tracklistpointer->fTrack.GetPsierr());
	      //HLTInfo("Tglerr: 1st:%f  2nd:%f   value:%f",tracklistpointer->fTrack.GetTglerr(),trackmatchingpointer->fTrack.GetTglerr(),(tracklistpointer->fTrack.GetTglerr()-trackmatchingpointer->fTrack.GetTglerr())/tracklistpointer->fTrack.GetTglerr());
	      //#endif
	    }
	  else // otherwise fill histogram with absolute differences
	    {
	      firstxhisto->Fill(tracklistpointer->fTrack.GetFirstPointX()-trackmatchingpointer->fTrack.GetFirstPointX(),1);
	      firstyhisto->Fill(tracklistpointer->fTrack.GetFirstPointY()-trackmatchingpointer->fTrack.GetFirstPointY(),1);
	      firstzhisto->Fill(tracklistpointer->fTrack.GetFirstPointZ()-trackmatchingpointer->fTrack.GetFirstPointZ(),1);
	      lastxhisto->Fill(tracklistpointer->fTrack.GetLastPointX()-trackmatchingpointer->fTrack.GetLastPointX(),1);
	      lastyhisto->Fill(tracklistpointer->fTrack.GetLastPointY()-trackmatchingpointer->fTrack.GetLastPointY(),1);
	      lastzhisto->Fill(tracklistpointer->fTrack.GetLastPointZ()-trackmatchingpointer->fTrack.GetLastPointZ(),1);
	      pthisto->Fill(tracklistpointer->fTrack.GetPt()-trackmatchingpointer->fTrack.GetPt(),1);
	      psihisto->Fill(tracklistpointer->fTrack.GetPsi()-trackmatchingpointer->fTrack.GetPsi(),1);
	      tglhisto->Fill(tracklistpointer->fTrack.GetTgl()-trackmatchingpointer->fTrack.GetTgl(),1);
	      clusterhisto->Fill(tracklistpointer->fTrack.GetNHits()-trackmatchingpointer->fTrack.GetNHits(),1);
	      //#if 0
	      // commented out since pterr is always zero for every track!
	      pterrhisto->Fill(tracklistpointer->fTrack.GetPterr()-trackmatchingpointer->fTrack.GetPterr(),1);
	      psierrhisto->Fill(tracklistpointer->fTrack.GetPsierr()-trackmatchingpointer->fTrack.GetPsierr(),1);
	      tglerrhisto->Fill(tracklistpointer->fTrack.GetTglerr()-trackmatchingpointer->fTrack.GetTglerr(),1);
	      //#endif
	    }

	  // fill histogram that determines the quality of the fit
	  matchinghisto->Fill(tracklistpointer->fMatchingindicator,1);
	}

      tracklistpointer = tracklistpointer->fNext;
    }

  trackeff->Divide(matchedtrackeff, firsttracks,1,1,"");

  // write histograms to root file specified in command line argument -graphs <filename>.ROOT
  if(!fGraphFileName.IsNull())
    {
      TFile* graphrootfile = new TFile(fGraphFileName, "update");
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
      graphrootfile->WriteObject(firsttracks, "firsttrackeff");
      graphrootfile->WriteObject(matchedtrackeff, "secondtrackeff");
      graphrootfile->WriteObject(trackeff, "trackingefficiency");

      // errors in tracking (commented out since pterr is always 0 for every track!)
      graphrootfile->WriteObject(pterrhisto, "pterrhistogram");
      //pterrhisto->Write();
      graphrootfile->WriteObject(psierrhisto, "psierrhistogram");
      //psierrhisto->Write();
      graphrootfile->WriteObject(tglerrhisto, "tglerrhistogram");
      //tglerrhisto->Write();
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
	 HLTInfo("%d : Discarding track with %d clusters.", trashtrackcounter, trackprintpointer->fTrack.GetNHits());
	 //PYTHIA INFORMATION ABOUT PARTICLE ID
	 HLTInfo("Track parameters of discarded track:");	
	 HLTInfo("First x: %9.6f \t first y: %9.6f \t first z: %9.6f",trackprintpointer->fTrack.GetFirstPointX(),trackprintpointer->fTrack.GetFirstPointY(),trackprintpointer->fTrack.GetFirstPointZ()); 
	 HLTInfo(" Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f", trackprintpointer->fTrack.GetLastPointX(),trackprintpointer->fTrack.GetLastPointY(),trackprintpointer->fTrack.GetLastPointZ());
	 HLTInfo("     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f", trackprintpointer->fTrack.GetPt(), trackprintpointer->fTrack.GetPsi(), trackprintpointer->fTrack.GetTgl());
	 
	 // write results to file if specified by command line argument
	 if(!fDumpFileName.IsNull())
	   {
	     fprintf(dumpfile,"%d : Discarding track with %d clusters. \n", trashtrackcounter, trackprintpointer->fTrack.GetNHits());
	     fprintf(dumpfile,"Track parameters of discarded track: \n");	
	     fprintf(dumpfile,"First x: %9.6f \t first y: %9.6f \t first z: %9.6f \n",trackprintpointer->fTrack.GetFirstPointX(),trackprintpointer->fTrack.GetFirstPointY(),trackprintpointer->fTrack.GetFirstPointZ()); 
	     fprintf(dumpfile," Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f \n", trackprintpointer->fTrack.GetLastPointX(),trackprintpointer->fTrack.GetLastPointY(),trackprintpointer->fTrack.GetLastPointZ());
	     fprintf(dumpfile,"     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f \n", trackprintpointer->fTrack.GetPt(), trackprintpointer->fTrack.GetPsi(), trackprintpointer->fTrack.GetTgl());
	   };
	 
	 // comparison with pythia information
	 if(trackprintpointer->fWronglydiscarded)
	   {
	     HLTInfo("Found track has been wrongly discarded according to Pythia information.");
	     
	     // write results to file if specified by command line argument
	     if(!fDumpFileName.IsNull())
	       {
		  fprintf(dumpfile,"Found track has been wrongly discarded accoring to Pythia information. \n");
	       };
	     
	    };
	 
	 // found pt must be in range pt_pythia \pm 10% to be accepted
	 Double_t ptmin = trackprintpointer->fPythiatrack.GetPt() - 0.1*trackprintpointer->fPythiatrack.GetPt();
	 Double_t ptmax = trackprintpointer->fPythiatrack.GetPt() + 0.1*trackprintpointer->fPythiatrack.GetPt();
	 
	 if( (trackprintpointer->fTrack.GetPt() < ptmin) ||(trackprintpointer->fTrack.GetPt() > ptmax) )
	   {
	     HLTInfo("Pt of found track %f differs more than 10 %% from pt of pythia track %f.",trackprintpointer->fTrack.GetPt(), trackprintpointer->fPythiatrack.GetPt());
	     
	     if(!fDumpFileName.IsNull())
	       {
		 // write result to file if specified by command line argument
		 fprintf(dumpfile,"Pt of found track %f differs more than 10 %% from pt of pythia track %f. \n",trackprintpointer->fTrack.GetPt(), trackprintpointer->fPythiatrack.GetPt());
	       };
	     
	   };
	 
	 HLTInfo("--------------");
	 
	 if(!fDumpFileName.IsNull())
	   {
	     fprintf(dumpfile,"-------------- \n");
	   }; 
	 
	 
	  // go to next element in trash track list
	 tracklistdeleter = trackprintpointer;
	 trackprintpointer = trackprintpointer->fNext;
	 
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
  HLTInfo("---------------CLUSTER ANALYSIS---------------");
  if(CompareClusters(1) != 0)
    { 
      return EINVAL;
    }

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
  HLTInfo("---------------TRACK ANALYSIS---------------");
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
  FILE* infofile = fopen("/afsuser/jwagner/TrackerTest_25092007/pp-ca/fullanalysis08012008/trackingefficiency.out","a");
  FILE* info2file = fopen("/afsuser/jwagner/TrackerTest_25092007/pp-ca/fullanalysis08012008/parameters.out", "a");

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
      
      if(listprintpointer->fMatchingindicator != 0)
	{
	  
#if 0
	  HLTInfo("Track %d:", tracknumber);
	  HLTInfo("Original track matched to secondary track with matchingindicator %d", listprintpointer->fMatchingindicator);
	  
	  HLTInfo("Parameter comparison: Original vs. Secondary");

	  HLTInfo("Clusters: %d \t %d ",listprintpointer->fTrack.GetNHits(), listprintpointer->fMatchingtrack->fTrack.GetNHits());
	  HLTInfo("First x: %9.6f \t %9.6f", listprintpointer->fTrack.GetFirstPointX(), listprintpointer->fMatchingtrack->fTrack.GetFirstPointX());
	  HLTInfo("First y: %9.6f \t %9.6f", listprintpointer->fTrack.GetFirstPointY(), listprintpointer->fMatchingtrack->fTrack.GetFirstPointY());
	  HLTInfo("First z: %9.6f \t %9.6f", listprintpointer->fTrack.GetFirstPointZ(), listprintpointer->fMatchingtrack->fTrack.GetFirstPointZ()); 
	  HLTInfo("Last x: %9.6f \t %9.6f", listprintpointer->fTrack.GetLastPointX(), listprintpointer->fMatchingtrack->fTrack.GetLastPointX());
	  HLTInfo("Last y: %9.6f \t %9.6f", listprintpointer->fTrack.GetLastPointY(), listprintpointer->fMatchingtrack->fTrack.GetLastPointY());
	  HLTInfo("Last z: %9.6f \t %9.6f", listprintpointer->fTrack.GetLastPointZ(), listprintpointer->fMatchingtrack->fTrack.GetLastPointZ());
	  HLTInfo("    Pt: %9.6f \t %9.6f", listprintpointer->fTrack.GetPt(), listprintpointer->fMatchingtrack->fTrack.GetPt());  
	  HLTInfo("   Psi: %9.6f \t %9.6f", listprintpointer->fTrack.GetPsi(), listprintpointer->fMatchingtrack->fTrack.GetPsi());
	  HLTInfo("   Tgl: %9.6f \t %9.6f", listprintpointer->fTrack.GetTgl(), listprintpointer->fMatchingtrack->fTrack.GetTgl());

	  HLTInfo(" Pterr: %9.6f \t %9.6f", listprintpointer->fTrack.GetPterr(), listprintpointer->fMatchingtrack->fTrack.GetPterr());  
	  HLTInfo("Psierr: %9.6f \t %9.6f", listprintpointer->fTrack.GetPsierr(), listprintpointer->fMatchingtrack->fTrack.GetPsierr());
	  HLTInfo("Tglerr: %9.6f \t %9.6f", listprintpointer->fTrack.GetTglerr(), listprintpointer->fMatchingtrack->fTrack.GetTglerr());
	 
	  HLTInfo("--------------");
#endif
	  // print these results to file
	  if(!fDumpFileName.IsNull())
	    {
	      fprintf(dumpfile, "Track %d: \n", tracknumber);
	      fprintf(dumpfile, "Original track matched to secondary track with matchingindicator %d \n", listprintpointer->fMatchingindicator); 
	      fprintf(dumpfile, "Parameter comparison: Original vs. Secondary \n"); 
	      fprintf(dumpfile, "Clusters: %d \t %d \n ",listprintpointer->fTrack.GetNHits(), listprintpointer->fMatchingtrack->fTrack.GetNHits());
	      fprintf(dumpfile, "First x: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetFirstPointX(), listprintpointer->fMatchingtrack->fTrack.GetFirstPointX());
	      fprintf(dumpfile, "First y: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetFirstPointY(), listprintpointer->fMatchingtrack->fTrack.GetFirstPointY()); 
	      fprintf(dumpfile, "First z: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetFirstPointZ(), listprintpointer->fMatchingtrack->fTrack.GetFirstPointZ()); 
	      fprintf(dumpfile, "Last x: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetLastPointX(), listprintpointer->fMatchingtrack->fTrack.GetLastPointX());  
	      fprintf(dumpfile, "Last y: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetLastPointY(), listprintpointer->fMatchingtrack->fTrack.GetLastPointY());
	      fprintf(dumpfile, "Last z: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetLastPointZ(), listprintpointer->fMatchingtrack->fTrack.GetLastPointZ());
	      fprintf(dumpfile, "    Pt: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetPt(), listprintpointer->fMatchingtrack->fTrack.GetPt());  
	      fprintf(dumpfile, "   Psi: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetPsi(), listprintpointer->fMatchingtrack->fTrack.GetPsi()); 
	      fprintf(dumpfile, "   Tgl: %9.6f \t %9.6f \n", listprintpointer->fTrack.GetTgl(), listprintpointer->fMatchingtrack->fTrack.GetTgl());
	      fprintf(dumpfile, "--------------\n");
	    }
	  ++tracknumber;
	}
     
      listprintpointer = listprintpointer->fNext;
    }


  // print out not matched tracks from first track array:
  listprintpointer = fFirstTrackList;
  Int_t notmatchedtracknumber = 1;

  while(listprintpointer != NULL)
    {
      if(listprintpointer->fMatchingindicator == 0)
	{
#if 0	  
	  HLTInfo("Original Track, not matched with secondary track %d:", notmatchedtracknumber);
	  HLTInfo("Clusters: %d",listprintpointer->fTrack.GetNHits());
	  //PYTHIA INFORMATION ABOUT PARTICLE ID
	  HLTInfo("First x: %9.6f \t first y: %9.6f \t first z: %9.6f",listprintpointer->fTrack.GetFirstPointX(),listprintpointer->fTrack.GetFirstPointY(),listprintpointer->fTrack.GetFirstPointZ()); 
	  HLTInfo(" Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f", listprintpointer->fTrack.GetLastPointX(),listprintpointer->fTrack.GetLastPointY(),listprintpointer->fTrack.GetLastPointZ());
	  HLTInfo("     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f", listprintpointer->fTrack.GetPt(), listprintpointer->fTrack.GetPsi(), listprintpointer->fTrack.GetTgl());
	  HLTInfo("--------------");
#endif
 
	  // print these results to file
	  if(!fDumpFileName.IsNull())
	    {
	      fprintf(dumpfile, "Original Track, not matched with secondary track %d: \n", notmatchedtracknumber); 
	      fprintf(dumpfile, "Clusters: %d \n",listprintpointer->fTrack.GetNHits());
	      fprintf(dumpfile, "First x: %9.6f \t first y: %9.6f \t first z: %9.6f \n",listprintpointer->fTrack.GetFirstPointX(),listprintpointer->fTrack.GetFirstPointY(),listprintpointer->fTrack.GetFirstPointZ());
	      fprintf(dumpfile, " Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f \n", listprintpointer->fTrack.GetLastPointX(),listprintpointer->fTrack.GetLastPointY(),listprintpointer->fTrack.GetLastPointZ());
	      fprintf(dumpfile, "     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f \n", listprintpointer->fTrack.GetPt(), listprintpointer->fTrack.GetPsi(), listprintpointer->fTrack.GetTgl()); 
	      fprintf(dumpfile, "--------------\n");
	    }

	  ++notmatchedtracknumber;
	}

      listprintpointer = listprintpointer->fNext;
    }

  // print out not matched tracks from second track array:
  listprintpointer = fSecondTrackList;
  notmatchedtracknumber = 1;

  while(listprintpointer != NULL)
    {
      if(listprintpointer->fMatchingindicator == 0)
	{
#if 0	  
	  HLTInfo("Secondary Track, not matched with original track %d:", notmatchedtracknumber);
	  HLTInfo("Clusters: %d",listprintpointer->fTrack.GetNHits());
	  //PYTHIA INFORMATION ABOUT PARTICLE ID	
	  HLTInfo("First x: %9.6f \t first y: %9.6f \t first z: %9.6f",listprintpointer->fTrack.GetFirstPointX(),listprintpointer->fTrack.GetFirstPointY(),listprintpointer->fTrack.GetFirstPointZ()); 
	  HLTInfo(" Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f", listprintpointer->fTrack.GetLastPointX(),listprintpointer->fTrack.GetLastPointY(),listprintpointer->fTrack.GetLastPointZ());
	  HLTInfo("     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f", listprintpointer->fTrack.GetPt(), listprintpointer->fTrack.GetPsi(), listprintpointer->fTrack.GetTgl());
	  HLTInfo("--------------");
#endif
	  // print these results to file
	  if(!fDumpFileName.IsNull())
	    {
	      fprintf(dumpfile, "Secondary Track, not matched with original track %d: \n", notmatchedtracknumber);
	      fprintf(dumpfile, "Clusters: %d \n",listprintpointer->fTrack.GetNHits());
	      fprintf(dumpfile, "First x: %9.6f \t first y: %9.6f \t first z: %9.6f \n",listprintpointer->fTrack.GetFirstPointX(),listprintpointer->fTrack.GetFirstPointY(),listprintpointer->fTrack.GetFirstPointZ()); 
	      fprintf(dumpfile, " Last x: %9.6f \t  last y: %9.6f \t  last z: %9.6f \n", listprintpointer->fTrack.GetLastPointX(),listprintpointer->fTrack.GetLastPointY(),listprintpointer->fTrack.GetLastPointZ());
	      fprintf(dumpfile, "     Pt: %9.6f \t     Psi: %9.6f \t     Tgl: %9.6f \n", listprintpointer->fTrack.GetPt(), listprintpointer->fTrack.GetPsi(), listprintpointer->fTrack.GetTgl());
	      fprintf(dumpfile, "--------------\n");
	    }

	  ++notmatchedtracknumber;
	}

      listprintpointer = listprintpointer->fNext;
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
