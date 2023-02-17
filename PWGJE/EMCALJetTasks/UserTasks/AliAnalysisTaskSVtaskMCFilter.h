#ifndef AliAnalysisTaskSVtaskMCFilter_H
#define AliAnalysisTaskSVtaskMCFilter_H

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Author : Filip Krizek
//-----------------------------------------------------------------------


#include "AliAnalysisTaskEmcal.h"

class TString;
class TClonesArray;
class AliAODMCHeader;
class AliStack;
class AliNormalizationCounter;

class AliAnalysisTaskSVtaskMCFilter : public AliAnalysisTaskEmcal
{

 public:

   AliAnalysisTaskSVtaskMCFilter();
   AliAnalysisTaskSVtaskMCFilter(const Char_t* name);
   virtual ~AliAnalysisTaskSVtaskMCFilter();

   static AliAnalysisTaskSVtaskMCFilter*  AddTaskSVtaskMCFilter(
     const char* trkcontname = "tracks",
	   const char* outtrk      = "mytracks",
	   Bool_t fFilterTracks    = kTRUE,
     const char* clscontname = "caloClusters",
     const char* outcls      = "myclusters",
     Bool_t fFilterClusters  = kFALSE
   );

   void     UserCreateOutputObjects();
   Bool_t   Run();

   void     SetInputTracksName(const Char_t* name){      fInputTracksName      = name;  }
   void     SetInputClustersName(const Char_t* name){    fInputClustersName    = name;  }
   void     SetFilteredTracksName(const Char_t* name){   fFilteredTracksName   = name;  }
   void     SetFilteredClustersName(const Char_t* name){ fFilteredClustersName = name;  }
   void     SetFilterType(Int_t param){                  fFilterType           = param; }
   TString  GetGenerator(Int_t label, AliAODMCHeader* header);

 protected:
   void ExecOnce();


   TClonesArray   *fFilteredTracksArray;   //! output PYTHIA only particles or tracks
   TClonesArray   *fFilteredClustersArray; //! output PYTHIA only clusters
   TString	       fFilteredTracksName;    //  name of output particle or track container
   TString         fFilteredClustersName;  //  name of output cluster container
   TString	       fInputTracksName;       //  name of input particle or track container
   TString         fInputClustersName;     //  name of input cluster container
   Int_t           fFilterType;            //  filtering for tracks=1 or mcParticles=0  //AID

   AliAODEvent    *fAodEvent;              //!
   AliAODMCHeader *fMCHeader;              //!
   TClonesArray   *fMCPartArray;           //! mcparticles


 private:

    AliAnalysisTaskSVtaskMCFilter(const AliAnalysisTaskSVtaskMCFilter &source);
    AliAnalysisTaskSVtaskMCFilter& operator=(const AliAnalysisTaskSVtaskMCFilter& source);

    ClassDef(AliAnalysisTaskSVtaskMCFilter, 4);
};

#endif
