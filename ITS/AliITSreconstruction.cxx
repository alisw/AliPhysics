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
 
/* $Id$ */
/////////////////////////////////////////////////////////////////////////
//                                                                     //
//                                                                     //
//                                                                     //
////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TSeqCollection.h>
#include <TString.h>
#include <TClonesArray.h>

#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliITS.h"
#include "AliITSDetType.h"
#include "AliITSLoader.h"
#include "AliITSreconstruction.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSClusterFinderSPD.h"
#include "AliITSClusterFinderSDD.h"
#include "AliITSClusterFinderSSD.h"
#include "AliITSresponseSDD.h"
#include "AliITSgeom.h"

ClassImp(AliITSreconstruction)

//______________________________________________________________________
AliITSreconstruction::AliITSreconstruction():
 fInit(kFALSE),
 fEnt(0),
 fEnt0(0),
 fITS(0x0),
 fDfArp(kFALSE),
 fLoader(0x0),
 fRunLoader(0x0)
{
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A zero-ed constructed AliITSreconstruction class.
    fDet[0] = fDet[1] = fDet[2] = kTRUE;
}
//______________________________________________________________________

AliITSreconstruction::AliITSreconstruction(AliRunLoader *rl):
 fInit(kFALSE),
 fEnt(0),
 fEnt0(0),
 fITS(0x0),
 fDfArp(kFALSE),
 fLoader(0x0),
 fRunLoader(rl)
{
  fDet[0] = fDet[1] = fDet[2] = kTRUE;
}
//______________________________________________________________________
AliITSreconstruction::AliITSreconstruction(const char* filename):
 fInit(kFALSE),
 fEnt(0),
 fEnt0(0),
 fITS(0x0),
 fDfArp(kFALSE),
 fLoader(0x0),
 fRunLoader(0x0)
{
    // Standard constructor.
    // Inputs:
    //  const char* filename    filename containing the digits to be
    //                          reconstructed. If filename = 0 (nil)
    //                          then no file is opened but a file is
    //                          assumed to already be opened. This 
    //                          already opened file will be used.
    // Outputs:
    //   none.
    // Return:
    //    A standardly constructed AliITSreconstruction class.

    fDet[0] = fDet[1] = fDet[2] = kTRUE;

    fRunLoader = AliRunLoader::Open(filename);
    if (fRunLoader == 0x0)
     {
       Error("AliITSreconstruction","Can not load the session",filename);
       return;
     }
    fRunLoader->LoadgAlice();
    gAlice = fRunLoader->GetAliRun();

    if(!gAlice) {
          Error("AliITSreconstruction","gAlice not found on file. Aborting.");
          fInit = kFALSE;
          return;
      } // end if !gAlice

}
//______________________________________________________________________
AliITSreconstruction::~AliITSreconstruction(){
    //    A destroyed AliITSreconstruction class.
    delete fRunLoader;
    fITS      = 0;
}
//______________________________________________________________________
Bool_t AliITSreconstruction::Init(){
    // Class Initilizer.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    kTRUE if no errors initilizing this class occurse else kFALSE
    Info("Init","");
    if (fRunLoader == 0x0)
     {
       Error("Init","Run Loader is NULL");
       return kFALSE;
     }
    fRunLoader->LoadgAlice();
    fRunLoader->LoadHeader();  

    fLoader = (AliITSLoader*) fRunLoader->GetLoader("ITSLoader");
    if(!fLoader) {
      Error("Init","ITS loader not found");
      fInit = kFALSE;
    }

    //Int_t retcode;
    fITS = (AliITS*) gAlice->GetDetector("ITS");
    if(!fITS){
      cout << "ITS not found aborting. fITS=" << fITS << endl;
      fInit = kFALSE;
      return fInit;
    } // end if !fITS
    if(!(fITS->GetITSgeom())){
      cout << "ITSgeom not found aborting."<< endl;
      fInit = kFALSE;
      return fInit;
    } // end if !GetITSgeom()
    // Now ready to init.

    fDet[0] = fDet[1] = fDet[2] = kTRUE;
    fEnt0 = 0;

    //fEnt  = gAlice->GetEventsPerRun();
    fEnt = Int_t(fRunLoader->TreeE()->GetEntries());

    fLoader->LoadDigits("read");
    fLoader->LoadRecPoints("recreate");
    fLoader->LoadRawClusters("recreate");
    if (fLoader->TreeR() == 0x0) fLoader->MakeTree("R");
    fITS->MakeBranch("R");
    fITS->MakeBranchC();
    fITS->SetTreeAddress();
    fInit = InitRec();

    Info("Init","  Done\n\n\n");

    return fInit;
}
//______________________________________________________________________
Bool_t AliITSreconstruction::InitRec(){
    // Sets up Reconstruction part of AliITSDetType..
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.
    AliITSDetType *idt;

    // SPD
    if(fDet[kSPD]){
      Info("InitRec","SPD");
      idt = fITS->DetType(kSPD);
      AliITSsegmentationSPD *segSPD = (AliITSsegmentationSPD*)idt->GetSegmentationModel();
      TClonesArray *digSPD = fITS->DigitsAddress(kSPD);
      TClonesArray *recpSPD = fITS->ClustersAddress(kSPD);
      Info("InitRec","idt = %#x; digSPD = %#x; recpSPD = %#x",idt,digSPD,recpSPD);
      AliITSClusterFinderSPD *recSPD = new AliITSClusterFinderSPD(segSPD,digSPD,recpSPD);
      fITS->SetReconstructionModel(kSPD,recSPD);
    } // end if fDet[kSPD].
    // SDD
    if(fDet[kSDD]){
      Info("InitRec","SDD");
      idt = fITS->DetType(kSDD);
      AliITSsegmentationSDD *segSDD = (AliITSsegmentationSDD*)
                                         idt->GetSegmentationModel();
      AliITSresponseSDD *resSDD = (AliITSresponseSDD*)
                                         idt->GetResponseModel();
      TClonesArray *digSDD = fITS->DigitsAddress(kSDD);
      TClonesArray *recpSDD = fITS->ClustersAddress(kSDD);
      AliITSClusterFinderSDD *recSDD =new AliITSClusterFinderSDD(segSDD,
                                                   resSDD,
                                                 digSDD,recpSDD);
      fITS->SetReconstructionModel(kSDD,recSDD);
    } // end if fDet[kSDD]
    // SSD
    if(fDet[kSSD]){
      Info("InitRec","SSD");
      idt = fITS->DetType(kSSD);
      AliITSsegmentationSSD *segSSD = (AliITSsegmentationSSD*)
                                       idt->GetSegmentationModel();
      TClonesArray *digSSD = fITS->DigitsAddress(kSSD);
      AliITSClusterFinderSSD *recSSD =new AliITSClusterFinderSSD(segSSD,
                                                   digSSD);
      fITS->SetReconstructionModel(kSSD,recSSD);
    } // end if fDet[kSSD]
    Info("InitRec","    Done\n");
    return kTRUE;
}
//______________________________________________________________________ 
void AliITSreconstruction::Exec(const Option_t *opt){
    // Main reconstruction function.
    // Inputs:
    //      Option_t * opt   list of subdetector to digitize. =0 all.
    // Outputs:
    //      none.
    // Return:
    //      none.
    Option_t *lopt;
    Int_t evnt;

    if(strstr(opt,"All")||strstr(opt,"ALL")||strstr(opt,"ITS")||opt==0){
      fDet[0] = fDet[1] = fDet[2] = kTRUE;
      lopt = "All";
    }else{
      fDet[0] = fDet[1] = fDet[2] = kFALSE;
      if(strstr(opt,"SPD")) fDet[kSPD] = kTRUE;
      if(strstr(opt,"SDD")) fDet[kSDD] = kTRUE;
      if(strstr(opt,"SSD")) fDet[kSSD] = kTRUE;
      if(fDet[kSPD] && fDet[kSDD] && fDet[kSSD]) lopt = "All";
      else lopt = opt;
    } // end if strstr(opt,...)

    if(!fInit){
      cout << "Initilization Failed, Can't run Exec." << endl;
      return;
    } // end if !fInit

    for(evnt=0;evnt<fEnt;evnt++)
     {
      Info("Exec","");
      Info("Exec","Processing Event %d",evnt);
      Info("Exec","");

      fRunLoader->GetEvent(evnt);
      if (fLoader->TreeR() == 0x0) fLoader->MakeTree("R");
      fITS->MakeBranch("R");
      if (fLoader->TreeC() == 0x0) fITS->MakeTreeC();
      fITS->SetTreeAddress();
      fITS->DigitsToRecPoints(evnt,0,lopt);
    } // end for evnt
}
//______________________________________________________________________ 
void AliITSreconstruction::SetOutputFile(TString filename){
  // Set a new file name for recpoints. 
  // It must be called before Init()
  if(!fLoader)fLoader = (AliITSLoader*) fRunLoader->GetLoader("ITSLoader");
  if(fLoader){
    Info("SetOutputFile","name for rec points is %s",filename.Data());
    fLoader->SetRecPointsFileName(filename);
  }
  else {
    Error("SetOutputFile",
    "ITS loader not available. Not possible to set name: %s",filename.Data());
  }
}
