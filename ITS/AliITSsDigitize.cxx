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

#include <Riostream.h>
 
#include <TROOT.h>
#include <TFile.h>
#include <TSeqCollection.h>
#include <TString.h>
#include <TClonesArray.h>
 
#include "AliHeader.h"
#include "AliRun.h"
 
#include "AliITS.h"
#include "AliITSsDigitize.h"
#include "AliITSgeom.h"

/////////////////////////////////////////////////////////
//                                                     //
//                                                     //
//                                                     //
///////////////////////////////////////////////////////// 
ClassImp(AliITSsDigitize)
//______________________________________________________________________
AliITSsDigitize::AliITSsDigitize():
fInit(kFALSE),
fEnt(0),
fEnt0(0),
fITS(0),
fRunLoader(0x0){
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A zero-ed constructed AliITSsDigitize class.
 
    fDet[0] = fDet[1] = fDet[2] = kTRUE;
}
//______________________________________________________________________
AliITSsDigitize::AliITSsDigitize(const char* filename):
fInit(),
fEnt(0),
fEnt0(0),
fITS(0),
fRunLoader(0x0){
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
    //    A standardly constructed AliITSsDigitize class.
 
    if(gAlice) 
     {
      delete gAlice;
      gAlice = 0;
     }
    fRunLoader = AliRunLoader::Open(filename);
    fRunLoader->LoadgAlice();
    fRunLoader->LoadHeader();
    
 
    Init();
}
//______________________________________________________________________
AliITSsDigitize::~AliITSsDigitize(){
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A destroyed AliITSsDigitize class.
 
    if(fRunLoader) delete fRunLoader;
    fITS      = 0;
 
}
//______________________________________________________________________
Bool_t AliITSsDigitize::Init(){
    // Class Initilizer.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    kTRUE if no errors initilizing this class occurse else kFALSE
    //Int_t nparticles;
   
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
 
    AliLoader* loader = fRunLoader->GetLoader("ITSLoader");
    
    if(!loader->TreeS()){
        cout << "Having to create the SDigits Tree." << endl;
        loader->MakeTree("S");
    } // end if !gAlice->TreeS()
    //make branch
    fITS->MakeBranch("S");
    fITS->SetTreeAddress();
    
    fRunLoader->GetEvent(fEnt0);
    // finished init.
    fInit = InitSDig();
    return fInit;
}
//______________________________________________________________________
Bool_t AliITSsDigitize::InitSDig() const {
    // Sets up SDigitization part of AliITSDetType..
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.
 
    return kTRUE;
}
 
//______________________________________________________________________
void AliITSsDigitize::Exec(const Option_t *opt){
    // Main SDigitization function.
    // Inputs:
    //      Option_t * opt   list of subdetector to digitize. =0 all.
    // Outputs:
    //      none.
    // Return:
    //      none.
    Option_t *lopt;
//    Int_t nparticles,evnt;
 
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

    fITS->HitsToSDigits(fRunLoader->GetHeader()->GetEvent(),0,-1," ",lopt," ");
}
