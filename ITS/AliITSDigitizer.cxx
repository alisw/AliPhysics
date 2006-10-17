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
///////////////////////////////////////////////////////////////////////////
//Piotr.Skowronski@cern.ch :                                             //
//Corrections applied in order to compile (only)                         // 
//   with new I/O and folder structure                                   //
//To be implemented correctly by responsible                             //
//                                                                       //
//  Class used to steer                                                  //
//  the digitization for ITS                                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranch.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunDigitizer.h"
#include "AliITSDigitizer.h"
#include "AliITSgeom.h"
#include "AliITSsimulation.h"

ClassImp(AliITSDigitizer)

//______________________________________________________________________
AliITSDigitizer::AliITSDigitizer() : AliDigitizer(){
    // Default constructor. Assign fITS since it is never written out from
    // here. 
    // Inputs:
    //      Option_t * opt   Not used
    // Outputs:
    //      none.
    // Return:
    //      A blank AliITSDigitizer class.

    fITS      = 0;
    fModActive   = 0;
    fRoif     = -1;
    fRoiifile = 0;
    fInit     = kFALSE;
    fFlagFirstEv =kTRUE;
}
//______________________________________________________________________
AliITSDigitizer::AliITSDigitizer(AliRunDigitizer *mngr) : AliDigitizer(mngr){
    // Standard constructor. Assign fITS since it is never written out from
    // here. 
    // Inputs:
    //      Option_t * opt   Not used
    // Outputs:
    //      none.
    // Return:
    //      An AliItSDigitizer class.

    fITS      = 0;
    fModActive   = 0;
    fRoif     = -1;
    fRoiifile = 0;
    fInit     = kFALSE;
    fFlagFirstEv =kTRUE;
}

//______________________________________________________________________
AliITSDigitizer::AliITSDigitizer(const AliITSDigitizer &/*rec*/):AliDigitizer(/*rec*/){
    // Copy constructor. 

  Error("Copy constructor","Copy constructor not allowed");
  
}
//______________________________________________________________________
AliITSDigitizer& AliITSDigitizer::operator=(const AliITSDigitizer& /*source*/){
    // Assignment operator. This is a function which is not allowed to be
    // done.
    Error("operator=","Assignment operator not allowed\n");
    return *this; 
}

//______________________________________________________________________
AliITSDigitizer::~AliITSDigitizer(){
    // Default destructor. 
    // Inputs:
    //      Option_t * opt   Not used
    // Outputs:
    //      none.
    // Return:
    //      none.
    fITS = 0; // don't delete fITS. Done else where.
    if(fModActive) delete[] fModActive;
}
//______________________________________________________________________
Bool_t AliITSDigitizer::Init(){
    // Initialization. Set up region of interest, if switched on, and
    // loads ITS and ITSgeom.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.

    fInit = kTRUE; // Assume for now init will work.
    if(!gAlice) {
	fITS      = 0;
	fRoiifile = 0;
	fInit     = kFALSE;
	Warning("Init","gAlice not found");
	return fInit;
    } // end if
    fITS = (AliITS *)(gAlice->GetDetector("ITS"));
    if(!fITS){
	fRoiifile = 0;
	fInit     = kFALSE;
	Warning("Init","ITS not found");
	return fInit;
    } else if(fITS->GetITSgeom()){
	//cout << "fRoif,fRoiifile="<<fRoif<<" "<<fRoiifile<<endl;
	fModActive = new Bool_t[fITS->GetITSgeom()->GetIndexMax()];
    } else{
	fRoiifile = 0;
	fInit     = kFALSE;
	Warning("Init","ITS geometry not found");
	return fInit;
    } // end if
    // fModActive needs to be set to a default all kTRUE value
    for(Int_t i=0;i<fITS->GetITSgeom()->GetIndexMax();i++) fModActive[i] = kTRUE;
    return fInit;
}
//______________________________________________________________________
void AliITSDigitizer::Exec(Option_t* opt){
    // Main digitization function. 
    // Inputs:
    //      Option_t * opt   list of sub detector to digitize. =0 all.
    // Outputs:
    //      none.
    // Return:
    //      none.

    char name[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    const char *all;
    const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
                          strstr(opt,"SSD")};
    if( !det[0] && !det[1] && !det[2] ) all = "All";
    else all = 0;
    Int_t nfiles = GetManager()->GetNinputs();
    Int_t event  = GetManager()->GetOutputEventNr();
    AliITSsimulation *sim      = 0;
    if(fFlagFirstEv){
      fITS->SetDefaults();    
      fITS->SetDefaultSimulation();
      fFlagFirstEv=kFALSE;
    }
    if(!fInit){
	Error("Exec","Init not successful, aborting.");
	return;
    } // end if

    sprintf(name,"%s",fITS->GetName());

    Int_t size   = fITS->GetITSgeom()->GetIndexMax();
    Int_t module,id,ifiles,mask;
    Bool_t lmod;
    Int_t *fl = new Int_t[nfiles];
    fl[0] = fRoiifile;
    mask = 1;
    for(id=0;id<nfiles;id++) 
     if(id!=fRoiifile)
      {
       // just in case fRoiifile!=0.
        fl[mask] = id;
        mask++;
      } // end for,if
    TClonesArray * sdig = new TClonesArray( "AliITSpListItem",1000 );
    
    TString loadname(name);
    loadname+="Loader";
    
    AliRunLoader *inRL = 0x0, *outRL = 0x0;
    AliLoader *ingime = 0x0, *outgime = 0x0;    
    
    outRL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());    
    if ( outRL == 0x0)
     {
       Error("Exec","Can not get Output Run Loader");
       return;
     }
    outRL->GetEvent(event);
    outgime = outRL->GetLoader(loadname);
    if ( outgime == 0x0)
     {
       Error("Exec","Can not get Output ITS Loader");
       return;
     }

    outgime->LoadDigits("update");
    if (outgime->TreeD() == 0x0) outgime->MakeTree("D");
    
    // Digitize
    fITS->MakeBranchInTreeD(outgime->TreeD());
    if(fRoif!=0) {
      AliDebug(1,"Region of Interest digitization selected");
    }
    else {
      AliDebug(1,"No Region of Interest selected. Digitizing everything");
    }
    if(fModActive==0) fRoif = 0; // fModActive array must be define for RIO cuts.

    for(ifiles=0; ifiles<nfiles; ifiles++ )
     {
       inRL =  AliRunLoader::GetRunLoader(fManager->GetInputFolderName(fl[ifiles]));
       ingime = inRL->GetLoader(loadname);
       if (ingime->TreeS() == 0x0) ingime->LoadSDigits();
     }

    for(module=0; module<size; module++ )
     {
       if(fModActive && fRoif!=0) if(!fModActive[module]) continue;
       id = fITS->GetITSgeom()->GetModuleType(module);
       if(!all && !det[id]) continue;
       sim      = (AliITSsimulation*)fITS->GetSimulationModel(id);
       if(!sim) {
            Error( "Exec", "The simulation class was not instanciated!" );
            exit(1);
        } // end if !sim
           // Fill the module with the sum of SDigits
        sim->InitSimulationModule(module, event);
	//cout << "Module=" << module;
        for(ifiles=0; ifiles<nfiles; ifiles++ )
         {
           if(fRoif!=0) if(!fModActive[module]) continue;
            
           inRL =  AliRunLoader::GetRunLoader(fManager->GetInputFolderName(fl[ifiles]));
           ingime = inRL->GetLoader(loadname);
           
           TTree *treeS = ingime->TreeS();
           fITS->SetTreeAddress();
           
           if( !(treeS && fITS->GetSDigits()) ) continue; 
           TBranch *brchSDigits = treeS->GetBranch( name );
           if( brchSDigits ) 
            {
                brchSDigits->SetAddress( &sdig ); 
            } else {
                Error( "Exec", "branch ITS not found in TreeS, input file %d ",
                       ifiles );
		return;
            } // end if brchSDigits
            sdig->Clear();
            mask = GetManager()->GetMask(ifiles);
            // add summable digits to module
            brchSDigits->GetEvent( module );
            lmod = sim->AddSDigitsToModule(sdig,mask);
            if(fRegionOfInterest && (ifiles==0))
             {
               fModActive[module] = lmod;
             } // end if
        } // end for ifiles
	//cout << " end ifiles loop" << endl;
        // Digitize current module sum(SDigits)->Digits
	sim->FinishSDigitiseModule();

        // fills all branches - wasted disk space
        outgime->TreeD()->Fill();
	fITS->ResetDigits();
    } // end for module
    
    outgime->TreeD()->AutoSave();
    outgime->WriteDigits("OVERWRITE");
    outgime->UnloadDigits();
    for(ifiles=0; ifiles<nfiles; ifiles++ )
     {
       inRL =  AliRunLoader::GetRunLoader(fManager->GetInputFolderName(fl[ifiles]));
       ingime = inRL->GetLoader(loadname);
       ingime->UnloadSDigits();
     }

    delete[] fl;
    sdig->Clear();
    delete sdig;
    for(Int_t i=0;i<fITS->GetITSgeom()->GetIndexMax();i++) fModActive[i] = kTRUE;
    

    return;
}
//______________________________________________________________________
void AliITSDigitizer::SetByRegionOfInterest(TTree *ts){
    // Scans through the ITS branch of the SDigits tree, ts, for modules
    // which have SDigits in them. For these modules, a flag is set to
    // digitize only these modules. The value of fRoif determines how many
    // neighboring modules will also be turned on. fRoif=0 will turn on only
    // those modules with SDigits in them. fRoif=1 will turn on, in addition,
    // those modules that are +-1 module from the one with the SDigits. And
    // So on. This last feature is not supported yet.
    // Inputs:
    //      TTree *ts  The tree in which the existing SDigits will define the
    //                 region of interest.
    // Outputs:
    //      none.
    // Return:
    //      none.
    Int_t m,nm,i;

    if(fRoif==0) return;
    if(ts==0) return;
    TBranch *brchSDigits = ts->GetBranch(fITS->GetName());
    TClonesArray * sdig = new TClonesArray( "AliITSpListItem",1000 );
    //cout << "Region of Interest ts="<<ts<<" brchSDigits="<<brchSDigits<<" sdig="<<sdig<<endl;

    if( brchSDigits ) {
      brchSDigits->SetAddress( &sdig );
    } else {
      Error( "SetByRegionOfInterest","branch ITS not found in TreeS");
      return;
    } // end if brchSDigits

    nm = fITS->GetITSgeom()->GetIndexMax();
    for(m=0;m<nm;m++){
      //cout << " fModActive["<<m<<"]=";
      fModActive[m] = kFALSE; // Not active by default
      sdig->Clear();
      brchSDigits->GetEvent(m);
      if(sdig->GetLast()>=0) for(i=0;i<sdig->GetLast();i++){
          // activate the necessary modules
          if(((AliITSpList*)sdig->At(m))->GetpListItem(i)->GetSignal()>0.0){ // Must have non zero signal.
            fModActive[m] = kTRUE;
            break;
          } // end if
      } // end if. end for i.
      //cout << fModActive[m];
      //cout << endl;
    } // end for m
    AliDebug(1,"Digitization by Region of Interest selected");
    sdig->Clear();
    delete sdig;
    return;
}
