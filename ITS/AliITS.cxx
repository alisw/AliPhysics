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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      An overview of the basic philosophy of the ITS code development      //
// and analysis is show in the figure below.                                 //
//Begin_Html                                                                 //
/*                                               
<img src="picts/ITS/ITS_Analysis_schema.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>Roberto Barbera is in charge of the ITS Offline code (1999).
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//
//  AliITS. Inner Traking System base class.
//  This class contains the base procedures for the Inner Tracking System
//
//Begin_Html
/*
<img src="picts/ITS/AliITS_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the class diagram of the different elements that are part of
the AliITS class.
</font>
<pre>
*/
//End_Html
//
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// Version: 2
// Modified and documented by A. Bologna
// October 18 1999
//
// AliITS is the general base class for the ITS. Also see AliDetector for
// futher information.
//
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include "AliDetector.h"
#include "AliITS.h"
#include "AliITSDetTypeSim.h"
#include "AliITSDDLRawData.h"
#include "AliITSLoader.h"
#include "AliITShit.h"
#include "AliITSmodule.h"
#include "AliITSpListItem.h"
#include "AliITSsimulation.h"
#include "AliITSsimulationFastPoints.h"
#include "AliMC.h"
#include "AliITSDigitizer.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSresponseSDD.h" 
#include "AliRawReader.h"
#include "AliRun.h"
#include "AliLog.h"
#include "AliITSInitGeometry.h"

ClassImp(AliITS)

//______________________________________________________________________
AliITS::AliITS() : AliDetector(),
fDetTypeSim(0),
fEuclidOut(0),
fOpt("All"),
fIdN(0),
fIdSens(0),
fIdName(0),
fITSmodules(0)
{
  // Default initializer for ITS
  //      The default constructor of the AliITS class. In addition to
  // creating the AliITS class it zeros the variables fIshunt (a member
  // of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
  // fIdSens, and fIdName. The AliDetector default constructor
  // is also called.
  
//    SetDetectors(); // default to fOpt="All". This variable not written out.
    SetMarkerColor(kRed);
}
//______________________________________________________________________
AliITS::AliITS(const char *name, const char *title):AliDetector(name,title),
fDetTypeSim(0),
fEuclidOut(0),
fOpt("All"),
fIdN(0),
fIdSens(0),
fIdName(0),
fITSmodules(0)
{
  //     The standard Constructor for the ITS class. 
  // It also zeros the variables
  // fIshunt (a member of AliDetector class), fEuclidOut, and zeros
  // the pointers fIdSens and fIdName. To help in displaying hits via the
  // ROOT macro display.C AliITS also sets the marker color to red. The
  // variables passes with this constructor, const char *name and *title,
  // are used by the constructor of AliDetector class. See AliDetector
  // class for a description of these parameters and its constructor
  // functions.
  
  fHits = new TClonesArray("AliITShit",1560);
  if(gAlice->GetMCApp()) gAlice->GetMCApp()->AddHitList(fHits);
  //fNhits=0;  //done in AliDetector(name,title)

  SetDetectors(); // default to fOpt="All". This variable not written out.
    
  fDetTypeSim   = new AliITSDetTypeSim();
  
  SetMarkerColor(kRed);
  if(!fLoader) MakeLoader(AliConfig::GetDefaultEventFolderName());
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);

}
//______________________________________________________________________
AliITS::~AliITS(){
    // Default destructor for ITS.
    //     The default destructor of the AliITS class. In addition to deleting
    // the AliITS class it deletes the memory pointed to by 
    // fIdSens, fIdName, fDetTypeSim and it's contents.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.

    if (fHits) {
      fHits->Delete();
      delete fHits;
      fHits=0;
    }
    if(fITSmodules) {
        this->ClearModules();
        delete fITSmodules;
	fITSmodules = 0;
    }// end if fITSmodules!=0

    delete[] fIdName;  // Array of TStrings
    delete[] fIdSens;

    if (fDetTypeSim){
      delete fDetTypeSim;
      fDetTypeSim = 0;
    }
}
//______________________________________________________________________
AliITS::AliITS(const AliITS &source) : AliDetector(source),
fDetTypeSim(0),
fEuclidOut(0),
fOpt("All"),
fIdN(0),
fIdSens(0),
fIdName(0),
fITSmodules(0)
{
    // Copy constructor. This is a function which is not allowed to be
    // done to the ITS. It exits with an error.
    // Inputs:
    //      AliITS &source  An AliITS class.
    // Outputs:
    //      none.
    // Return:
    //      none.

    if(this==&source) return;
    Error("Copy constructor",
          "You are not allowed to make a copy of the AliITS");
    exit(1);
}
//______________________________________________________________________
AliITS& AliITS::operator=(const AliITS &source){
    // Assignment operator. This is a function which is not allowed to be
    // done to the ITS. It exits with an error.
    // Inputs:
    //      AliITS &source  An AliITS class.
    // Outputs:
    //      none.
    // Return:
    //      none.

    if(this==&source) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITS");
    exit(1);
    return *this; //fake return
}
//______________________________________________________________________
AliDigitizer* AliITS::CreateDigitizer(AliRunDigitizer* manager)const{
    // Creates the AliITSDigitizer in a standard way for use via AliModule.
    // This function can not be included in the .h file because of problems
    // with the order of inclusion (recursive).
    // Inputs:
    //    AliRunDigitizer *manager  The Manger class for Digitization
    // Output:
    //    none.
    // Return:
    //    A new AliITSRunDigitizer (cast as a AliDigitizer).

     return new AliITSDigitizer(manager);
}
//______________________________________________________________________
void AliITS::Init(){
    // Initializer ITS after it has been built
    //     This routine initializes the AliITS class. It is intended to be
    // called from the Init function in AliITSv?. Besides displaying a banner
    // indicating that it has been called it initializes the array fIdSens
    // and sets the default segmentation, response, digit and raw cluster
    // classes therefore it should be called after a call to CreateGeometry.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.
    Int_t i;

    SetDefaults();
    // Array of TStrings
    if(gMC) for(i=0;i<fIdN;i++) fIdSens[i] = gMC->VolId(fIdName[i]);
 
    WriteGeometry();
}

//______________________________________________________________________
void AliITS::WriteGeometry(){
  
  //Writes ITS geometry on gAlice

  if(!fLoader) MakeLoader(AliConfig::GetDefaultEventFolderName());
  AliRunLoader* rl  = fLoader->GetRunLoader();
  rl->CdGAFile();
  AliITSgeom* geom = GetITSgeom();
  geom->Write();

}

//______________________________________________________________________
void AliITS::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.
    AliInfoClass("AliITS::Setting Defaults");

    if(!fDetTypeSim) { 
     Error("SetDefaults()","fDetTypeSim is 0!"); 
     return;
    }

    fDetTypeSim->SetDefaults();


}
//______________________________________________________________________
void AliITS::SetDefaultSimulation(){
    // sets the default simulation.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.
    if(!fDetTypeSim) { 
     Error("SetDefaultSimulation()","fDetTypeSim is 0!"); 
     return;
    }

    fDetTypeSim->SetDefaultSimulation();

}


//______________________________________________________________________
void AliITS::MakeBranch(Option_t* option){
    // Creates Tree branches for the ITS.
    // Inputs:
    //      Option_t *option    String of Tree types S,D, and/or R.
    //      const char *file    String of the file name where these branches
    //                          are to be stored. If blank then these branches
    //                          are written to the same tree as the Hits were
    //                          read from.
    // Outputs:
    //      none.
    // Return:
    //      none.
  if(!fDetTypeSim) {
    Error("MakeBranch","fDetTypeSim is 0!");
    return;
  }

  Bool_t cH = (strstr(option,"H")!=0);
  Bool_t cS = (strstr(option,"S")!=0);
  Bool_t cD = (strstr(option,"D")!=0);
  
  if(cH && (fHits == 0x0)) fHits = new TClonesArray("AliITShit", 1560);
  AliDetector::MakeBranch(option);
  
  if(cS) MakeBranchS(0);
  if(cD) MakeBranchD(0);


}
//___________________________________________________________________
void AliITS::MakeBranchS(const char* fl){

  // Creates Tree Branch for the ITS summable digits.
  // Inputs:
  //      cont char *fl  File name where SDigits branch is to be written
  //                     to. If blank it write the SDigits to the same
  //                     file in which the Hits were found.

  
  if(!fDetTypeSim){
    Error("MakeBranchS","fDetTypeSim is 0!");
  }
  Int_t buffersize = 4000;
  char branchname[30];

  // only one branch for SDigits.
  sprintf(branchname,"%s",GetName());

  if(fLoader->TreeS()){
    if(fDetTypeSim->GetSDigits()==0x0) fDetTypeSim->SetSDigits(new TClonesArray("AliITSpListItem",1000));
    TClonesArray* sdig = (TClonesArray*)fDetTypeSim->GetSDigits();
    MakeBranchInTree(fLoader->TreeS(),branchname,&sdig,buffersize,fl);
  } 
}
//______________________________________________________________________
void AliITS::MakeBranchD(const char* file){

  //Make branch for digits
  if(!fDetTypeSim) {
    Warning("MakeBranchD","fDetTypeSim is 0!");
    return;
  }
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  MakeBranchInTreeD(fLoader->TreeD(),file);
}

//___________________________________________________________________
void AliITS:: MakeBranchInTreeD(TTree* treeD, const char* file){
  // Creates Tree branches for the ITS.

  if(!fDetTypeSim){
    Error("MakeBranchS","fDetTypeSim is 0!");
  }
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);

  const Char_t *det[3] = {"SPD","SDD","SSD"};
  Char_t* digclass;
  Int_t buffersize = 4000;
  Char_t branchname[30];
  
  if(!fDetTypeSim->GetDigits()){
    fDetTypeSim->SetDigits(new TObjArray(fgkNTYPES));
  }
  for(Int_t i=0;i<fgkNTYPES;i++){
    digclass = fDetTypeSim->GetDigitClassName(i);
    TString classn = digclass;
    if(!((fDetTypeSim->GetDigits())->At(i))){
      (fDetTypeSim->GetDigits())->AddAt(new TClonesArray(classn.Data(),1000),i);
    }
    else ResetDigits(i);  
    if(fgkNTYPES==3) sprintf(branchname,"%sDigits%s",GetName(),det[i]);
    else sprintf(branchname,"%sDigits%d",GetName(),i+1);
    TObjArray* dig = DigitsAddress(i);
    if(GetDigits() && treeD) AliDetector::MakeBranchInTree(treeD,branchname, &dig,buffersize,file);
  }

}
//______________________________________________________________________
void AliITS::SetTreeAddress(){
    // Set branch address for the Trees.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.
    
  if(!fDetTypeSim) {
    Error("SetTreeAddress","fDetTypeSim is 0!");
    return;
  }

  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);

  TTree *treeS = fLoader->TreeS();
  TTree *treeD = fLoader->TreeD();
  if (fLoader->TreeH() && (fHits == 0x0)) {
      fHits = new TClonesArray("AliITShit", 1560);
  }
  AliDetector::SetTreeAddress();

  fDetTypeSim->SetTreeAddressS(treeS, (Char_t*)GetName());
  fDetTypeSim->SetTreeAddressD(treeD, (Char_t*)GetName());
}
//______________________________________________________________________
void AliITS::AddHit(Int_t track, Int_t *vol, Float_t *hits){
    // Add an ITS hit
    //     The function to add information to the AliITShit class. See the
    // AliITShit class for a full description. This function allocates the
    // necessary new space for the hit information and passes the variable
    // track, and the pointers *vol and *hits to the AliITShit constructor
    // function.
    // Inputs:
    //      Int_t   track   Track number which produced this hit.
    //      Int_t   *vol    Array of Integer Hit information. See AliITShit.h
    //      Float_t *hits   Array of Floating Hit information.  see AliITShit.h
    // Outputs:
    //      none.
    // Return:
    //      none.
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliITShit(fIshunt,track,vol,hits);
}

//______________________________________________________________________
void AliITS::FillModules(Int_t evnt,Int_t bgrev,Int_t nmodules,
                         Option_t *option, const char *filename){
  // fill the modules with the sorted by module hits; add hits from
  // background if option=Add.

  static TTree *trH1;                 //Tree with background hits
  static Bool_t first=kTRUE;
  static TFile *file;
  const char *addBgr = strstr(option,"Add");
  
  evnt = nmodules; // Dummy use of variables to remove warnings
  if (addBgr ) {
    if(first) {
      file=new TFile(filename);
    } // end if first
    first=kFALSE;
    file->cd();
    file->ls();
    // Get Hits Tree header from file
    if(trH1) delete trH1;
    trH1=0;
    
    char treeName[20];
    sprintf(treeName,"TreeH%d",bgrev);
    trH1 = (TTree*)gDirectory->Get(treeName);
    if (!trH1) {
      Error("FillModules","cannot find Hits Tree for event:%d",bgrev);
    } // end if !trH1
    // Set branch addresses
  } // end if addBgr
  
  FillModules(fLoader->TreeH(),0); // fill from this file's tree.
    
  if (addBgr ) {
    FillModules(trH1,10000000); // Default mask 10M.
    TTree *fAli=fLoader->GetRunLoader()->TreeK();
    TFile *fileAli=0;
    if (fAli) fileAli =fAli->GetCurrentFile();
    fileAli->cd();
  } // end if add
  
  
}


//______________________________________________________________________
void AliITS::FillModules(TTree *treeH, Int_t mask) {
    // fill the modules with the sorted by module hits; 
    // can be called many times to do a merging
    // Inputs:
    //      TTree *treeH  The tree containing the hits to be copied into
    //                    the modules.
    //      Int_t mask    The track number mask to indecate which file
    //                    this hits came from.
    // Outputs:
    //      none.
    // Return:
    //      none.

    if (treeH == 0x0)
     {
       Error("FillModules","Tree is NULL");
     }
    Int_t lay,lad,det,index;
    AliITShit *itsHit=0;
    AliITSmodule *mod=0;
    char branchname[20];
    sprintf(branchname,"%s",GetName());
    TBranch *branch = treeH->GetBranch(branchname);
    if (!branch) {
        Error("FillModules","%s branch in TreeH not found",branchname);
        return;
    } // end if !branch
    branch->SetAddress(&fHits);
    Int_t nTracks =(Int_t) treeH->GetEntries();
    Int_t iPrimTrack,h;
    for(iPrimTrack=0; iPrimTrack<nTracks; iPrimTrack++){
        ResetHits();
        Int_t nBytes = treeH->GetEvent(iPrimTrack);
        if (nBytes <= 0) continue;
        Int_t nHits = fHits->GetEntriesFast();
        for(h=0; h<nHits; h++){
            itsHit = (AliITShit *)fHits->UncheckedAt(h);
            itsHit->GetDetectorID(lay,lad,det);
            if (GetITSgeom()) {
                index = GetITSgeom()->GetModuleIndex(lay,lad,det);
            } else {
                index=det-1; // This should not be used.
            } // end if [You must have fITSgeom for this to work!]
            mod = GetModule(index);
            itsHit->SetTrack(itsHit->GetTrack()+mask); // Set track mask.
            mod->AddHit(itsHit,iPrimTrack,h);
        } // end loop over hits 
    } // end loop over tracks
}

//______________________________________________________________________
void AliITS::InitModules(Int_t size,Int_t &nmodules){
    // Initialize the modules array.
    // Inputs:
    //      Int_t size  Size of array of the number of modules to be
    //                  created. If size <=0 then the number of modules
    //                  is gotten from AliITSgeom class kept in fITSgeom.
    // Outputs:
    //      Int_t &nmodules The number of modules existing.
    // Return:
    //      none.

    if(fITSmodules){ 
        fITSmodules->Delete();
        delete fITSmodules;
    } // end fir fITSmoudles

    if(!fDetTypeSim) {
      Error("InitModules","fDetTypeSim is null!");
      return;
    }

    Int_t nl,indexMAX,index;

    if(size<=0){ // default to using data stored in AliITSgeom
        if(fDetTypeSim->GetITSgeom()==0) {
            Error("InitModules","fITSgeom not defined");
            return;
        } // end if fITSgeom==0
        nl = fDetTypeSim->GetITSgeom()->GetNlayers();
        indexMAX = fDetTypeSim->GetITSgeom()->GetIndexMax();
        nmodules = indexMAX;
        fITSmodules = new TObjArray(indexMAX);
        for(index=0;index<indexMAX;index++){
            fITSmodules->AddAt( new AliITSmodule(index),index);
        } // end for index
    }else{
        fITSmodules = new TObjArray(size);
        for(index=0;index<size;index++) {
            fITSmodules->AddAt( new AliITSmodule(index),index);
        } // end for index

        nmodules = size;
    } // end i size<=0
}
//______________________________________________________________________
void AliITS::Hits2SDigits(){
    // Standard Hits to summable Digits function.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
  

   if(!fDetTypeSim) {
     Error("Hits2SDigits","fDetTypeSim is null!");
     return;
    
  } 
     
  SetDefaults();
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* rl = fLoader->GetRunLoader(); 
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  for (Int_t iEvent = 0; iEvent < rl->GetNumberOfEvents(); iEvent++) {
        // Do the Hits to Digits operation. Use Standard input values.
        // Event number from file, no background hit merging , use size from
        // AliITSgeom class, option="All", input from this file only.
    rl->GetEvent(iEvent);
    if (!fLoader->TreeS()) fLoader->MakeTree("S");
    MakeBranch("S");
    SetTreeAddress();
    HitsToPreDigits(iEvent,0,-1," ",fOpt," ");
  } // end for iEvent
    
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  
}
//______________________________________________________________________
void AliITS::Hits2Digits(){

  //Conversion from hits to digits
  if(!fDetTypeSim) {
    Error("Hits2SDigits","fDetTypeSim is 0!");
    return;
  }
   
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  SetDefaults();

  fLoader->LoadHits("read");
  fLoader->LoadDigits("recreate");
  AliRunLoader* rl = fLoader->GetRunLoader(); 
  for (Int_t iEvent = 0; iEvent < rl->GetNumberOfEvents(); iEvent++) {
    rl->GetEvent(iEvent);
    if (!fLoader->TreeD()) fLoader->MakeTree("D");
    MakeBranch("D");
    SetTreeAddress();   
    HitsToDigits(iEvent,0,-1," ",fOpt," ");
  } 
  
  fLoader->UnloadHits();
  fLoader->UnloadDigits();
  
}

//______________________________________________________________________
void AliITS::HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                          Option_t *option,Option_t *opt,
                          const char *filename){
    //   Keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
    // the condition below will disappear when the geom class will be
    // initialized for all versions - for the moment it is only for v5 !
    // 7 is the SDD beam test version.
    // Inputs:
    //      Int_t evnt       Event to be processed.
    //      Int_t bgrev      Background Hit tree number.
    //      Int_t nmodules   Not used.
    //      Option_t *option String indicating if merging hits or not. To
    //                       merge hits set equal to "Add". Otherwise no
    //                       background hits are considered.
    //      Test_t *filename File name containing the background hits..
    // Outputs:
    //      none.
    // Return:
    //      none.

  if(!fDetTypeSim) {
    Error("HitsToDigits","fDetTypeSim is null!");
    return;
  }
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  if(!GetITSgeom()) return; // need transformations to do digitization.
  AliITSgeom *geom = GetITSgeom();

  const char *all = strstr(opt,"All");
  const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
			strstr(opt,"SSD")};
  static Bool_t setDef=kTRUE;
  if (setDef) SetDefaultSimulation();
  setDef=kFALSE;
  
  Int_t nmodules;
  InitModules(size,nmodules);
  FillModules(evNumber,bgrev,nmodules,option,filename);

  AliITSsimulation *sim      = 0;
  AliITSmodule     *mod      = 0;
  Int_t id;
  for(Int_t module=0;module<geom->GetIndexMax();module++){
    id       = geom->GetModuleType(module);
    if (!all && !det[id]) continue;
    sim      = (AliITSsimulation*)fDetTypeSim->GetSimulationModel(id);
    if (!sim) {
      Error("HitsToDigits","The simulation class was not "
	    "instanciated for module %d type %s!",module,
	    geom->GetModuleTypeName(module));
      exit(1);
    } // end if !sim
    mod      = (AliITSmodule *)fITSmodules->At(module);
    sim->DigitiseModule(mod,module,evNumber);
    // fills all branches - wasted disk space
    fLoader->TreeD()->Fill(); 
    ResetDigits();
  } // end for module
  
  ClearModules();
  
  fLoader->TreeD()->GetEntries();
  fLoader->TreeD()->AutoSave();
  // reset tree
  fLoader->TreeD()->Reset();
}
//_____________________________________________________________________
void AliITS::Hits2PreDigits(){ 
  // Turn hits into SDigits

  if(!fDetTypeSim) {
    Error("Hits2SDigits","fDetTypeSim is 0!");
    return;
  }
   
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  SetDefaults();
  
  HitsToPreDigits(fLoader->GetRunLoader()->GetEventNumber(),
		  0,-1," ",fOpt," ");
}

//______________________________________________________________________
void AliITS::HitsToPreDigits(Int_t evNumber,Int_t bgrev,Int_t size,
                             Option_t *option,Option_t *opt,
                             const char *filename){
    //   Keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
    // the condition below will disappear when the geom class will be
    // initialized for all versions - for the moment it is only for v5 !
    // 7 is the SDD beam test version.
    // Inputs:
    //      Int_t evnt       Event to be processed.
    //      Int_t bgrev      Background Hit tree number.
    //      Int_t nmodules   Not used.
    //      Option_t *option String indicating if merging hits or not. To
    //                       merge hits set equal to "Add". Otherwise no
    //                       background hits are considered.
    //      Test_t *filename File name containing the background hits..
    // Outputs:
    //      none.
    // Return:
    //      none.

 
  if(!fDetTypeSim) {
    Error("HitsToPreDigits","fDetTypeSim is null!");
    return;
  }
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);

  if(!GetITSgeom()){
    Error("HitsToPreDigits","fGeom is null!");
    return; // need transformations to do digitization.
  }
  AliITSgeom *geom = GetITSgeom();

  const char *all = strstr(opt,"All");
  const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),
			strstr(opt,"SSD")};
  static Bool_t setDef=kTRUE;
  if (setDef) SetDefaultSimulation();
  setDef=kFALSE;
  
  Int_t nmodules;
  InitModules(size,nmodules);
  FillModules(evNumber,bgrev,nmodules,option,filename);
  

  AliITSsimulation *sim      = 0;
  AliITSmodule     *mod      = 0;
  Int_t id,module;
  for(module=0;module<geom->GetIndexMax();module++){
    id       = geom->GetModuleType(module);
    if (!all && !det[id]) continue;
    sim      = (AliITSsimulation*)GetSimulationModel(id);
    if (!sim) {
      Error("HitsToPreDigits","The simulation class was not "
	    "instanciated for module %d type %s!",module,
	    geom->GetModuleTypeName(module));
      exit(1);
    } // end if !sim
    mod      = (AliITSmodule *)fITSmodules->At(module);
    sim->SDigitiseModule(mod,module,evNumber);
    // fills all branches - wasted disk space
    fLoader->TreeS()->Fill(); 
    fDetTypeSim->ResetSDigits();
  } // end for module

  ClearModules();

  
  fLoader->TreeS()->GetEntries();
  fLoader->TreeS()->AutoSave();
  fLoader->WriteSDigits("OVERWRITE");
  // reset tree
  fLoader->TreeS()->Reset();
}

//_____________________________________________________________________
void AliITS::HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
                                  Option_t *opt0,Option_t *opt1,
                                 const char *flnm){
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
    // the condition below will disappear when the geom class will be
    // initialized for all versions - for the moment it is only for v5 !
    // Inputs:
    //      Int_t evnt       Event to be processed.
    //      Int_t bgrev      Background Hit tree number.
    //      Int_t size       Size used by InitModules. See InitModules.
    //      Option_t *opt0   Option passed to FillModules. See FillModules.
    //      Option_t *opt1   String indicating if merging hits or not. To
    //                       merge hits set equal to "Add". Otherwise no
    //                       background hits are considered.
    //      Test_t *flnm     File name containing the background hits..
    // Outputs:
    //      none.
    // Return:
    //      none.



  if(!GetITSgeom()){
    Error("HitsToPreDigits","fGeom is null!");
    return; // need transformations to do digitization.
  }
  AliITSgeom *geom = GetITSgeom();

  AliITSLoader *pITSloader = (AliITSLoader*)fLoader;

  const char *all = strstr(opt1,"All");
  const char *det[3] ={strstr(opt1,"SPD"),strstr(opt1,"SDD"),
		       strstr(opt1,"SSD")};
  Int_t nmodules;
  InitModules(size,nmodules);
  FillModules(evNumber,bgrev,nmodules,opt0,flnm);

  AliITSsimulation *sim      = 0;
  AliITSmodule     *mod      = 0;
  Int_t id,module;

  TTree *lTR = pITSloader->TreeR();
  if(!lTR) {
    pITSloader->MakeTree("R");
    lTR = pITSloader->TreeR();
  }
  
  TClonesArray* ptarray = new TClonesArray("AliITSRecPoint",1000);
  TBranch* branch = (TBranch*)lTR->Branch("ITSRecPointsF",&ptarray);
  branch->SetAddress(&ptarray);
  //m.b. : this change is nothing but a nice way to make sure
  //the CPU goes up !    
  for(module=0;module<geom->GetIndexMax();module++){
    id       = geom->GetModuleType(module);
    if (!all && !det[id]) continue;
    sim      = (AliITSsimulation*)GetSimulationModel(id);
    if (!sim) {
      Error("HitsToFastPoints","The simulation class was not "
	    "instanciated for module %d type %x!",module,
	    geom->GetModuleTypeName(module));
      exit(1);
    } // end if !sim
    mod      = (AliITSmodule *)fITSmodules->At(module);
    sim->CreateFastRecPoints(mod,module,gRandom,ptarray);
    lTR->Fill();
  } // end for module

  ClearModules();
  fLoader->WriteRecPoints("OVERWRITE");
  delete ptarray;
}
//_____________________________________________________________________
Int_t AliITS::Hits2Clusters(TTree *hTree, TTree *cTree) {
  //------------------------------------------------------------
  // This function creates ITS clusters
  //------------------------------------------------------------
  if(!GetITSgeom()){
    Error("HitsToPreDigits","fGeom is null!");
    return 1; // need transformations to do digitization.
  }
  AliITSgeom *geom=GetITSgeom();
  Int_t mmax=geom->GetIndexMax();

  InitModules(-1,mmax);
  FillModules(hTree,0);

  TClonesArray *points = new TClonesArray("AliITSRecPoint",1000);
  TBranch *branch=cTree->GetBranch("ITSRecPoints");
  if (!branch) cTree->Branch("ITSRecPoints",&points);
  else branch->SetAddress(&points);

  AliITSsimulationFastPoints sim;
  Int_t ncl=0;
  for (Int_t m=0; m<mmax; m++) {
    AliITSmodule *mod=GetModule(m);      
    sim.CreateFastRecPoints(mod,m,gRandom,points);      
    ncl+=points->GetEntriesFast();
    cTree->Fill();
    points->Clear();
  }

  Info("Hits2Clusters","Number of found fast clusters : %d",ncl);

  //cTree->Write();

  delete points;
  return 0;
}

//_____________________________________________________________________
void AliITS::CheckLabels(Int_t lab[3]) const {
  //------------------------------------------------------------
  // Tries to find mother's labels
  //------------------------------------------------------------

  if(lab[0]<0 && lab[1]<0 && lab[2]<0) return; // In case of no labels just exit

  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
  for (Int_t i=0;i<3;i++){
    Int_t label = lab[i];
    if (label>=0 && label<ntracks) {
      TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
      if (part->P() < 0.005) {
	Int_t m=part->GetFirstMother();
	if (m<0) {	
	  continue;
	}
	if (part->GetStatusCode()>0) {
	  continue;
	}
	lab[i]=m;       
      }
    }    
  }
  
}

//______________________________________________________________________
void AliITS::SDigitsToDigits(Option_t *opt){
    // Standard Summable digits to Digits function.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    if(!fDetTypeSim) {
      Error("SDigitsToSDigits","fDetTypeSim is 0!");
      return;
    }
   
    fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
    SetDefaults();
    fDetTypeSim->SDigitsToDigits(opt,(Char_t*)GetName());

}

//______________________________________________________________________
void AliITS::ResetDigits(){
    // Reset number of digits and the digits array for the ITS detector.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    if(!fDetTypeSim) {
      Error("ResetDigits","fDetTypeSim is 0!");
      return;
    }
   
    fDetTypeSim->ResetDigits();


}
//______________________________________________________________________
void AliITS::ResetDigits(Int_t branch){
    // Reset number of digits and the digits array for this branch.
    // Inputs:
    //      none.
    // Outputs:
    //      none.

    if(!fDetTypeSim) {
      Error("ResetDigits","fDetTypeSim is 0!");
      return;
    }
   
    fDetTypeSim->ResetDigits(branch);

}
//______________________________________________________________________
void AliITS::AddSumDigit(AliITSpListItem &sdig){
    // Adds the a module full of summable digits to the summable digits tree.
    // Inputs:
    //      AliITSpListItem &sdig   SDigit to be added to SDigits tree.
    // Outputs:
    //      none.
    // Return:
    //      none.

    if(!fDetTypeSim) {
      Error("AddSumDigit","fDetTypeSim is 0!");
      return;
    }
    fDetTypeSim->AddSumDigit(sdig);
    
}
//______________________________________________________________________
void AliITS::AddRealDigit(Int_t branch, Int_t *digits){
    //   Add a real digit - as coming from data.
    // Inputs:
    //      Int_t id        Detector type number.
    //      Int_t *digits   Integer array containing the digits info. See 
    //                      AliITSdigit.h
    // Outputs:
    //      none.
    // Return:
    //      none.

    if(!fDetTypeSim) {
      Error("AddRealDigit","fDetTypeSim is 0!");
      return;
    }
    fDetTypeSim->AddRealDigit(branch,digits);

}
//______________________________________________________________________
void AliITS::AddSimDigit(Int_t branch, AliITSdigit *d){
    //    Add a simulated digit.
    // Inputs:
    //      Int_t id        Detector type number.
    //      AliITSdigit *d  Digit to be added to the Digits Tree. See 
    //                      AliITSdigit.h
    // Outputs:
    //      none.
    // Return:
    //      none.

    if(!fDetTypeSim) {
      Error("AddSimDigit","fDetTypeSim is 0!");
      return;
    }
    fDetTypeSim->AddSimDigit(branch,d);

}
//______________________________________________________________________
void AliITS::AddSimDigit(Int_t branch,Float_t phys,Int_t *digits,Int_t *tracks,
                         Int_t *hits,Float_t *charges){
  //   Add a simulated digit to the list.
  // Inputs:
  //      Int_t id        Detector type number.
  //      Float_t phys    Physics indicator. See AliITSdigits.h
  //      Int_t *digits   Integer array containing the digits info. See 
  //                      AliITSdigit.h
  //      Int_t *tracks   Integer array [AliITSdigitS?D::GetNTracks()] 
  //                      containing the track numbers that contributed to
  //                      this digit.
  //      Int_t *hits     Integer array [AliITSdigitS?D::GetNTracks()]
  //                      containing the hit numbers, from AliITSmodule, that
  //                      contributed to this digit.
  //      Float_t *charge Floating point array of the signals contributed
  //                      to this digit by each track.
  // Outputs:
  //      none.
  // Return:
  //      none.

    if(!fDetTypeSim) {
      Error("AddSimDigit","fDetTypeSim is 0!");
      return;
    }
    fDetTypeSim->AddSimDigit(branch,phys,digits,tracks,hits,charges);

}
//______________________________________________________________________
void AliITS::Digits2Raw(){
    // convert digits of the current event to raw data

  if(!fDetTypeSim) {
    Error("Digits2Raw","fDetTypeSim is 0!");
    return;
  }
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  SetDefaults();
  fDetTypeSim->GetLoader()->LoadDigits();
  TTree* digits = fDetTypeSim->GetLoader()->TreeD();
  if (!digits) {
      Error("Digits2Raw", "no digits tree");
      return;
  }
  fDetTypeSim->SetTreeAddressD(digits,(Char_t*)GetName());

  AliITSDDLRawData rawWriter;
  //Verbose level
  // 0: Silent
  // 1: cout messages
  // 2: txt files with digits 
  //BE CAREFUL, verbose level 2 MUST be used only for debugging and
  //it is highly suggested to use this mode only for debugging digits files
  //reasonably small, because otherwise the size of the txt files can reach
  //quickly several MB wasting time and disk space.
  rawWriter.SetVerbose(0);
    
  //SILICON PIXEL DETECTOR
  Info("Digits2Raw", "Formatting raw data for SPD");
  rawWriter.RawDataSPD(digits->GetBranch("ITSDigitsSPD"));
    
  //SILICON DRIFT DETECTOR
  Info("Digits2Raw", "Formatting raw data for SDD");
  rawWriter.RawDataSDD(digits->GetBranch("ITSDigitsSDD"));
    
  //SILICON STRIP DETECTOR
  Info("Digits2Raw", "Formatting raw data for SSD");
  rawWriter.RawDataSSD(digits->GetBranch("ITSDigitsSSD"));

  fLoader->UnloadDigits();
}
//______________________________________________________________________
AliLoader* AliITS::MakeLoader(const char* topfoldername){ 
    //builds ITSgetter (AliLoader type)
    //if detector wants to use castomized getter, it must overload this method

    AliDebug(1,Form("Creating AliITSLoader. Top folder is %s.",
         topfoldername));
    fLoader = new AliITSLoader(GetName(),topfoldername);
    return fLoader;
}

Bool_t AliITS::Raw2SDigits(AliRawReader* rawReader)
{
  //
  // Converts RAW data to SDigits
  //
  // Get TreeS
  //
    Int_t last   = -1;
    Int_t size   = GetITSgeom()->GetIndexMax();
    TClonesArray** modA = new TClonesArray*[size];
    for (Int_t mod = 0; mod < size; mod++) modA[mod] = new TClonesArray("AliITSpListItem", 10000);
    
    AliLoader* loader =  (gAlice->GetRunLoader())->GetLoader("ITSLoader");
    if (!loader)
    {
	Error("Open","Can not get ITS loader from Run Loader");
	return kFALSE;
    }

    TTree* tree = 0;
    tree = loader->TreeS();
    if (!tree)
    {
	loader->MakeTree("S");
	tree = loader->TreeS();
    }
    //
    // Array for SDigits
    // 
    TClonesArray aSDigits("AliITSpListItem",10000), *itsSDigits=&aSDigits;
    Int_t bufsize = 32000;
    tree->Branch("ITS", &itsSDigits, bufsize);
    Int_t npx = 0;
    //
    // SPD
    //
    AliITSsegmentationSPD* segSPD = (AliITSsegmentationSPD*) fDetTypeSim->GetSegmentationModel(0);
    npx = segSPD->Npx();
    Double_t thr, sigma; 
    
    AliITSRawStreamSPD inputSPD(rawReader);
    while(1){
	Bool_t next  = inputSPD.Next();
	if (!next) break;

	Int_t module = inputSPD.GetModuleID();
	Int_t column = inputSPD.GetColumn();
	Int_t row    = inputSPD.GetRow();
	Int_t index  = npx * column + row;

	if (module >= size) continue;
 
	last = (modA[module])->GetEntries();
	TClonesArray& dum = *modA[module];
	fDetTypeSim->GetCalibrationModel(module)->Thresholds(thr,sigma);
	thr += 1.;
	new (dum[last]) AliITSpListItem(-1, -1, module, index, thr);
    }
    rawReader->Reset();

    //
    // SDD
    // 
    AliITSsegmentationSDD* segSDD = (AliITSsegmentationSDD*) fDetTypeSim->GetSegmentationModel(1);
    npx = segSDD->Npx();
    AliITSRawStreamSDD inputSDD(rawReader);
    while(1){
	Bool_t next  = inputSDD.Next();
	if (!next) break;

	Int_t module = inputSDD.GetModuleID();
	Int_t anode  = inputSDD.GetAnode();
	Int_t time   = inputSDD.GetTime();
	Int_t signal = inputSDD.GetSignal();
	Int_t index  = npx * anode + time;

	if (module >= size) continue;
	// 8bit -> 10 bit
	AliITSresponseSDD *resSDD = (AliITSresponseSDD*) fDetTypeSim->GetResponse(1);
	Int_t signal10 = resSDD->Convert8to10(signal);  // signal is a 8 bit value (if the compression is active)
	
	last = modA[module]->GetEntries();
	TClonesArray& dum = *modA[module];
	new (dum[last]) AliITSpListItem(-1, -1, module, index, Double_t(signal10));
	((AliITSpListItem*) dum.At(last))->AddSignalAfterElect(module, index, Double_t(signal10));
	
    }
    rawReader->Reset();

    //
    // SSD
    // 
    AliITSsegmentationSSD* segSSD = (AliITSsegmentationSSD*) fDetTypeSim->GetSegmentationModel(2);
    npx = segSSD->Npx();
    AliITSRawStreamSSD inputSSD(rawReader);
    while(1){
	Bool_t next  = inputSSD.Next();
	if (!next) break;

	Int_t module  = inputSSD.GetModuleID();
	Int_t side    = inputSSD.GetSideFlag();
	Int_t strip   = inputSSD.GetStrip();
	Int_t signal  = inputSSD.GetSignal();
	Int_t index  = npx * side + strip;

	if (module >= size) continue;
	
	last = modA[module]->GetEntries();
	TClonesArray& dum = *modA[module];
	new (dum[last]) AliITSpListItem(-1, -1, module, index, Double_t(signal));
    }
    rawReader->Reset();
     AliITSpListItem* sdig = 0;
    
    for (Int_t mod = 0; mod < size; mod++)
    {
	Int_t nsdig =  modA[mod]->GetEntries();
	for (Int_t ie = 0; ie < nsdig; ie++) {
	    sdig = (AliITSpListItem*) (modA[mod]->At(ie));
	    new (aSDigits[ie]) AliITSpListItem(-1, -1, mod, sdig->GetIndex(), sdig->GetSignal());
	    Float_t sig = sdig->GetSignalAfterElect();
	    if (sig > 0.) {
		sdig = (AliITSpListItem*)aSDigits[ie];
		sdig->AddSignalAfterElect(mod, sdig->GetIndex(), Double_t(sig));
	    }
	}
	
	tree->Fill();
	aSDigits.Clear();
	modA[mod]->Clear();
    }
    loader->WriteSDigits("OVERWRITE");    
    delete modA;
    return kTRUE;
}


//______________________________________________________________________
void AliITS::UpdateInternalGeometry(){

  //reads new geometry from TGeo 
  Info("UpdateInternalGeometry", "Delete ITSgeom and create a new one reading TGeo");
  AliITSInitGeometry *initgeom = new AliITSInitGeometry("AliITSvPPRasymmFMD",2);
  AliITSgeom* geom = initgeom->CreateAliITSgeom();
  SetITSgeom(geom);

  if(!fLoader) MakeLoader(AliConfig::GetDefaultEventFolderName());
  AliRunLoader* rl  = fLoader->GetRunLoader();
  rl->CdGAFile();
  geom->Write(0,kOverwrite);

}

