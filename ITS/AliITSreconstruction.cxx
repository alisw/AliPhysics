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
 
/*
$Log$
*/
#include <TROOT.h>
#include <TFile.h>
#include <TSeqCollection.h>
#include <TString.h>
#include <TClonesArray.h>

#include "AliRun.h"

#include "AliITS.h"
#include "AliITSDetType.h"
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
AliITSreconstruction::AliITSreconstruction(){
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A zero-ed constructed AliITSreconstruction class.

    fFilename = "";
    fFile     = 0;
    fITS      = 0;
    fDet[0] = fDet[1] = fDet[2] = kTRUE;
    fInit     = kFALSE;
}
//______________________________________________________________________
AliITSreconstruction::AliITSreconstruction(const char* filename){
    // Standard constructor.
    // Inputs:
    //  const char* filename    filename containing the digits to be
    //                          reconstructed
    // Outputs:
    //   none.
    // Return:
    //    A standardly constructed AliITSreconstruction class.

    fFilename = filename;

    fFile = (TFile*)gROOT->GetListOfFiles()->FindObject(fFilename.Data());
    if(fFile) fFile->Close();
    fFile = new TFile(fFilename.Data(),"UPDATE");
    //
    if(gAlice) delete gAlice;
    gAlice = (AliRun*)fFile->Get("gAlice");
    if(!gAlice) {
	cout << "gAlice not found on file. Aborting." << endl;
	fInit = kFALSE;
	return;
    } // end if !gAlice

    Init();
}
//______________________________________________________________________
AliITSreconstruction::~AliITSreconstruction(){
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A destroyed AliITSreconstruction class.

    if(fFile) fFile->Close();
    fFile     = 0;
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
    Int_t nparticles;

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
    fEnt  = gAlice->GetEventsPerRun();
    fITS->MakeTreeC();
    nparticles = gAlice->GetEvent(fEnt0);
    
    // finished init.
    fInit = InitRec();
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
	idt = fITS->DetType(kSPD);
	AliITSsegmentationSPD *segSPD = (AliITSsegmentationSPD*)
	                                       idt->GetSegmentationModel();
	TClonesArray *digSPD = fITS->DigitsAddress(kSPD);
	TClonesArray *recpSPD = fITS->ClustersAddress(kSPD);
	AliITSClusterFinderSPD *recSPD = new AliITSClusterFinderSPD(segSPD,
								    digSPD,
								    recpSPD);
	fITS->SetReconstructionModel(kSPD,recSPD);
    } // end if fDet[kSPD].
    // SDD
    if(fDet[kSDD]){
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
	idt = fITS->DetType(kSSD);
	AliITSsegmentationSSD *segSSD = (AliITSsegmentationSSD*)
	                                 idt->GetSegmentationModel();
	TClonesArray *digSSD = fITS->DigitsAddress(kSSD);
	AliITSClusterFinderSSD *recSSD =new AliITSClusterFinderSSD(segSSD,
								   digSSD);
	fITS->SetReconstructionModel(kSSD,recSSD);
    } // end if fDet[kSSD]

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
    Int_t nparticles,evnt;

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
    for(evnt=0;evnt<fEnt;evnt++){
	nparticles = gAlice->GetEvent(evnt);
	gAlice->SetEvent(evnt);
	if(!gAlice->TreeR()) gAlice->MakeTree("R");
	fITS->MakeBranch("R");
	fITS->DigitsToRecPoints(evnt,0,lopt);
    } // end for evnt
}
