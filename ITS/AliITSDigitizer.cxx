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
Revision 1.8  2003/03/03 17:00:30  masera
Corrections to comply with coding conventions

Revision 1.7  2002/10/25 18:54:22  barbera
Various improvements and updates from B.S.Nilsen and T. Virgili

Revision 1.6  2002/10/22 14:45:34  alibrary
Introducing Riostream.h

Revision 1.5  2002/10/14 14:57:00  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.3.4.1  2002/06/10 17:51:14  hristov
Merged with v3-08-02

Revision 1.4  2002/04/24 22:08:12  nilsen
New ITS Digitizer/merger with two macros. One to make SDigits (old way) and
one to run the  merger (modified for Jiri).

Revision 1.3  2002/03/25 10:48:55  nilsen
New ITS SDigit merging with region of interest cut. Update for changes in
AliDigitizer. Additional optimization should be done.

Revision 1.2  2002/03/15 17:26:40  nilsen
New SDigit version of ITS Digitizer.

Revision 1.1  2001/11/27 16:27:28  nilsen
Adding AliITSDigitizer class to do merging and digitization . Based on the
TTask method. AliITSDigitizer class added to the Makefile and ITSLinkDef.h
file. The following files required minor changes. AliITS, added functions
SetHitsAddressBranch, MakeBranchInTreeD and modified MakeBranchD.
AliITSsimulationSDD.cxx needed a Tree independent way of returning back to
the original Root Directory in function Compress1D. Now it uses gDirectory.

Revision 1.2  2002/03/01  E. Lopez
Digitization changed to start from SDigits instead of Hits.
The SDigits are reading as TClonesArray of AliITSpListItem
*/
//
//  Class used to steer
//  the digitization for ITS
//
//
#include <stdlib.h>
#include <Riostream.h>
#include <TClonesArray.h> 
#include <TTree.h>
#include <TBranch.h>

#include <AliRun.h>
#include <AliRunDigitizer.h>

#include "AliITSDigitizer.h"
#include "AliITSgeom.h"
#include "AliITSpList.h"
#include "AliITSsimulation.h"
#include "AliITSDetType.h"


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
/*  This will not work from Init. ts is aways returned as zero.
    TTree *ts;
    if(fRoif>=0 && fRoiifile>=0 && fRoiifile<GetManager()->GetNinputs()){
	ts = GetManager()->GetInputTreeS(fRoiifile);
	if(!ts){
	    if(!gAlice) ts = gAlice->TreeS();
	    if(!ts){
		cout <<"The TTree TreeS needed to set by region not found."
		    " No region of interest cut will be applied."<< endl;
		return fInit;
	    } // end if
	} // end if
	cout << "calling SetByReionOfInterest ts="<< ts <<endl;
	SetByRegionOfInterest(ts);
    } // end if
*/
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
    char *all;
    const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
                          strstr(opt,"SSD")};
    if( !det[0] && !det[1] && !det[2] ) all = "All";
    else all = 0;
    AliITSsimulation *sim      = 0;
    AliITSDetType    *iDetType = 0;
    static Bool_t    setDef    = kTRUE;

    if(!fInit){
	Error("Exec","Init not succesfull, aborting.");
	return;
    } // end if

    if( setDef ) fITS->SetDefaultSimulation();
    setDef = kFALSE;
    sprintf(name,"%s",fITS->GetName());

    Int_t nfiles = GetManager()->GetNinputs();
    Int_t event  = GetManager()->GetOutputEventNr();
    Int_t size   = fITS->GetITSgeom()->GetIndexMax();
    Int_t module,id,ifiles,mask;
    Bool_t lmod;
    Int_t *fl = new Int_t[nfiles];
    fl[0] = fRoiifile;
    mask = 1;
    for(id=0;id<nfiles;id++) if(id!=fRoiifile){
	// just in case fRoiifile!=0.
	fl[mask] = id;
	mask++;
    } // end for,if
    TClonesArray * sdig = new TClonesArray( "AliITSpListItem",1000 );
    
    // Digitize
    fITS->MakeBranchInTreeD(GetManager()->GetTreeD());
    if(fRoif!=0) Info("AliITSDigitizer","Region of Interest digitization selected");
    else Info("AliITSDigitizer","No Region of Interest selected. Digitizing everything");
    //cout <<"fModActive="<<fModActive<<" fRoif="<<fRoif;
    if(fModActive==0) fRoif = 0; // fModActive array must be define for RIO cuts.
    //cout <<" fRoif="<<fRoif<<endl;

    for(module=0; module<size; module++ ){
	if(fRoif!=0) if(!fModActive[module]) continue;
        id = fITS->GetITSgeom()->GetModuleType(module);
        if(!all && !det[id]) continue;
        iDetType = fITS->DetType( id );
        sim      = (AliITSsimulation*)iDetType->GetSimulationModel();
        if(!sim) {
            Error( "Exec", "The simulation class was not instanciated!" );
            exit(1);
        } // end if !sim

        // Fill the module with the sum of SDigits
        sim->InitSimulationModule(module, event);
	//cout << "Module=" << module;
        for(ifiles=0; ifiles<nfiles; ifiles++ ){
	    if(fRoif!=0) if(!fModActive[module]) continue;
	    //cout <<" fl[ifiles=" << ifiles << "]=" << fl[ifiles];
            TTree *treeS = GetManager()->GetInputTreeS(fl[ifiles]);
            if( !(treeS && fITS->GetSDigits()) ) continue;   
            TBranch *brchSDigits = treeS->GetBranch( name );
            if( brchSDigits ) {
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
	    if(ifiles==0){
		fModActive[module] = lmod;
	    } // end if
	    //cout << " fModActive["<<module<<"]=";
	    //if(fModActive[module]) cout << "kTRUE";
	    //else cout << "kFALSE";
        } // end for ifiles
	//cout << " end ifiles loop" << endl;
        // Digitize current module sum(SDigits)->Digits
        sim->FinishSDigitiseModule();

        // fills all branches - wasted disk space
        GetManager()->GetTreeD()->Fill();
        fITS->ResetDigits();
    } // end for module
    //cout << "end modules loop"<<endl;

    GetManager()->GetTreeD()->AutoSave();

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
    Info("AliITSDigitizer","Digitization by Region of Interest selected");
    sdig->Clear();
    delete sdig;
    return;
}
