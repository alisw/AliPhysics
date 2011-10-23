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
#include <TArrayI.h>
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
#include "AliITSsimulationSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSRawStreamSDD.h"
#include "AliRawReader.h"
#include "AliRun.h"
#include "AliLog.h"
#include "AliITSInitGeometry.h"
#include "AliITSFOSignalsSPD.h"

ClassImp(AliITS)

//______________________________________________________________________
AliITS::AliITS() : AliDetector(),
fDetTypeSim(0),
fEuclidOut(0),
fOpt("All"),
fIdN(0),
fIdSens(0),
fIdName(0),
fITSmodules(0),
fTiming(kFALSE),
fSimuParam(0),
fModA(0),
fpSDigits(0)
{
  // Default initializer for ITS
  //      The default constructor of the AliITS class. In addition to
  // creating the AliITS class it zeros the variables fIshunt (a member
  // of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
  // fIdSens, and fIdName. The AliDetector default constructor
  // is also called.
  
//    SetDetectors(); // default to fOpt="All". This variable not written out.
//PH    SetMarkerColor(kRed);
  for (int i=fgkNTYPES;i--;) fRawID2ClusID[i] = 0;
}
//______________________________________________________________________
AliITS::AliITS(const Char_t *title):
  AliDetector("ITS",title),
  fDetTypeSim(0),
  fEuclidOut(0),
  fOpt("All"),
  fIdN(0),
  fIdSens(0),
  fIdName(0),
  fITSmodules(0),
  fTiming(kFALSE),
  fSimuParam(0),
  fModA(0),
  fpSDigits(0)
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
    // Inputs:
    //   Char_t *title  Simulation title for the ITS
    // Outputs:
    //   none.
    // Return:
    //   none.
  
    fHits = new TClonesArray("AliITShit",1560); // from AliDetector
    if(gAlice->GetMCApp()) gAlice->GetMCApp()->AddHitList(fHits);
    //fNhits=0;  //done in AliDetector(name,title)
    SetDetectors(); // default to fOpt="All". This variable not written out.
    fDetTypeSim   = new AliITSDetTypeSim();
    //PH  SetMarkerColor(kRed);
    if(!fLoader) MakeLoader(AliConfig::GetDefaultEventFolderName());
    fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
    for (int i=fgkNTYPES;i--;) fRawID2ClusID[i] = 0;
}
//______________________________________________________________________
AliITS::AliITS(const char *name, const char *title):
  AliDetector(name,title),
  fDetTypeSim(0),
  fEuclidOut(0),
  fOpt("All"),
  fIdN(0),
  fIdSens(0),
  fIdName(0),
  fITSmodules(0),
  fTiming(kFALSE),
  fSimuParam(0),
  fModA(0),
  fpSDigits(0)
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
  //PH  SetMarkerColor(kRed);
  if(!fLoader) MakeLoader(AliConfig::GetDefaultEventFolderName());
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  for (int i=fgkNTYPES;i--;) fRawID2ClusID[i] = 0;
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
    Int_t size   = AliITSgeomTGeo::GetNModules();
    if (fDetTypeSim){
      delete fDetTypeSim;
      fDetTypeSim = 0;
    }
    if(fSimuParam){
      delete fSimuParam;
      fSimuParam=0;
    }
    if(fModA){
      if(size>0){
	for(Int_t j=0; j<size; j++){
	  fModA[j]->Delete();
	  delete fModA[j];
	}
      }
      delete []fModA;
    }
    if(fpSDigits){
      fpSDigits->Delete();
      delete fpSDigits;
    }
}

//______________________________________________________________________
AliDigitizer* AliITS::CreateDigitizer(AliDigitizationInput* digInput)const{
    // Creates the AliITSDigitizer in a standard way for use via AliModule.
    // This function can not be included in the .h file because of problems
    // with the order of inclusion (recursive).
    // Inputs:
    //    AliDigitizationInput* digInput  The Manger class for Digitization
    // Output:
    //    none.
    // Return:
    //    A new AliITSRunDigitizer (cast as a AliDigitizer).

     return new AliITSDigitizer(digInput);
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
    // Array of TStrings
    if(gMC) for(i=0;i<fIdN;i++) fIdSens[i] = gMC->VolId(fIdName[i]);
 
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
    if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

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
    if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

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
  char branchname[31];

  // only one branch for SDigits.
  snprintf(branchname,30,"%s",GetName());

  if(fLoader->TreeS()){
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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);
  MakeBranchInTreeD(fLoader->TreeD(),file);
}

//___________________________________________________________________
void AliITS:: MakeBranchInTreeD(TTree* treeD, const char* file){
  // Creates Tree branches for the ITS.

  if(!fDetTypeSim){
    Error("MakeBranchS","fDetTypeSim is 0!");
  }
  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

  const Char_t *det[3] = {"SPD","SDD","SSD"};
  const Char_t* digclass;
  Int_t buffersize = 4000;
  Char_t branchname[31];
  
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
    if(fgkNTYPES==3) snprintf(branchname,30,"%sDigits%s",GetName(),det[i]);
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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

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
void AliITS::FillModules(Int_t /* evnt */,Int_t bgrev,Int_t /* nmodules */,
                         Option_t *option, const char *filename){
  // fill the modules with the sorted by module hits; add hits from
  // background if option=Add.

  static TTree *trH1;                 //Tree with background hits
  static Bool_t first=kTRUE;
  static TFile *file;
  const char *addBgr = strstr(option,"Add");
  
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
    
    char treeName[21];
    snprintf(treeName,20,"TreeH%d",bgrev);
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
    if (fAli) {
      fileAli =fAli->GetCurrentFile();
      fileAli->cd();
    }
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
       AliError("Tree H  is NULL");
       return;
     }
    Int_t lay,lad,det,index;
    AliITShit *itsHit=0;
    AliITSmodule *mod=0;
    char branchname[21];
    snprintf(branchname,20,"%s",GetName());
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
Bool_t AliITS::InitModules(Int_t size,Int_t &nmodules){
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
      return kFALSE;
    }
    if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

    Int_t nl,indexMAX,index;

    if(size<=0){ // default to using data stored in AliITSgeom
        if(fDetTypeSim->GetITSgeom()==0) {
            Error("InitModules","fITSgeom not defined");
            return kFALSE;
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
    return kTRUE;
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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);
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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);
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
 
  // Reset Fast-OR signals for this event
  fDetTypeSim->ResetFOSignals();

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
 
  // Add Fast-OR signals to event (only one object per event)
  if (all || det[0]) { // SPD present
    WriteFOSignals();
  }
  
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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);
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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);

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
	    "instantiated for module %d type %s!",module,
	    geom->GetModuleTypeName(module));
      exit(1);
    } // end if !sim
    mod      = (AliITSmodule *)fITSmodules->At(module);
    sim->CreateFastRecPoints(mod,module,gRandom,ptarray);
    lTR->Fill();
    ptarray->Clear();
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

  AliDebug(1,Form("Number of found fast clusters : %d",ncl));

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
  if (!fDetTypeSim) {
    AliError("fDetTypeSim is 0!");
    return;
  }
  const char *all = strstr(opt,"All");
  const char *det[3] ={strstr(opt,"SPD"),strstr(opt,"SDD"),
		       strstr(opt,"SSD")};

  // Reset Fast-OR signals for this event
  fDetTypeSim->ResetFOSignals();

  fDetTypeSim->SetLoader((AliITSLoader*)fLoader);
  SetDefaults();
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);
  fDetTypeSim->SDigitsToDigits(opt,(Char_t*)GetName());

  // Add Fast-OR signals to event (only one object per event)
  if (all || det[0]) { // SPD present
   WriteFOSignals();
  }
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
                         Int_t *hits,Float_t *charges, Int_t sigexpanded){
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
    fDetTypeSim->AddSimDigit(branch,phys,digits,tracks,hits,charges,sigexpanded);

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
  if(fSimuParam) fDetTypeSim->SetSimuParam(fSimuParam);
  fDetTypeSim->GetLoader()->LoadDigits();
  TTree* digits = fDetTypeSim->GetLoader()->TreeD();
  if (!digits) {
      Error("Digits2Raw", "no digits tree");
      return;
  }
  fDetTypeSim->SetTreeAddressD(digits,(Char_t*)GetName());
  
   // Get the FO signals for this event
  AliITSFOSignalsSPD* foSignals = NULL;
  AliRunLoader* runLoader = AliRunLoader::Instance();
  AliITSLoader* itsLoader = (AliITSLoader*) runLoader->GetLoader("ITSLoader");
  if (!itsLoader) {
    AliError("ITS loader is NULL.");
  }
   else {
      if(!itsLoader->TreeD()) AliError("   !!! No TreeD available !!!");
      foSignals = (AliITSFOSignalsSPD*)itsLoader->TreeD()->GetUserInfo()->FindObject("AliITSFOSignalsSPD");
      if(!foSignals) AliError("FO signals not retrieved");
     }
 
  Bool_t deleteFOsignalsLater = kFALSE;
  if (!foSignals) {
    AliError("FO signals not available. No FO bits will be written.");
    foSignals = new AliITSFOSignalsSPD(); // make a temporary dummy signals object
    deleteFOsignalsLater = kTRUE;
  }
  
  
  AliITSDDLModuleMapSDD* ddlsdd=fDetTypeSim->GetDDLModuleMapSDD();
  Char_t rawSDD=fDetTypeSim->GetSimuParam()->GetSDDRawDataFormat();
  AliITSDDLRawData rawWriter;
  
  rawWriter.SetSDDRawFormat(rawSDD);
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
  AliDebug(1,"Formatting raw data for SPD");
  rawWriter.RawDataSPD(digits->GetBranch("ITSDigitsSPD"),foSignals);
  if(deleteFOsignalsLater) delete foSignals;
    
  //SILICON DRIFT DETECTOR
  AliDebug(1,Form("Formatting raw data for SDD - Format code =%d",rawSDD));
  rawWriter.RawDataSDD(digits->GetBranch("ITSDigitsSDD"),ddlsdd);
    
  //SILICON STRIP DETECTOR
  AliDebug(1,"Formatting raw data for SSD");
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
//______________________________________________________________________
Bool_t AliITS::Raw2SDigits(AliRawReader* rawReader)
{
  //
  // Converts RAW data to SDigits
  //
  // Get TreeS
  //
  Int_t last   = -1;
  Int_t size   = AliITSgeomTGeo::GetNModules();
  if(!fModA) {
    fModA = new TClonesArray*[size];
    for (Int_t mod = 0; mod < size; mod++) fModA[mod] = new TClonesArray("AliITSpListItem", 10000);
  }
  AliLoader* loader =  (AliRunLoader::Instance())->GetLoader("ITSLoader");
  if (!loader){
    Error("Open","Can not get ITS loader from Run Loader");
    return kFALSE;
  }

  TTree* tree = 0;
  tree = loader->TreeS();
  if (!tree){
    loader->MakeTree("S");
    tree = loader->TreeS();
  }
  //
  // Array for SDigits
  // 
  if(!fpSDigits){
    fpSDigits = new TClonesArray("AliITSpListItem",10000);
  }
  TClonesArray& aSDigits = *fpSDigits;
  Int_t bufsize = 32000;
  tree->Branch("ITS", &fpSDigits, bufsize);
  Int_t npx = 0;
  //
  // SPD
  //
  AliITSsegmentationSPD* segSPD = (AliITSsegmentationSPD*) fDetTypeSim->GetSegmentationModel(0);
  if(!segSPD){
    AliWarning("Set AliITS defaults");
    SetDefaults();
    segSPD = (AliITSsegmentationSPD*) fDetTypeSim->GetSegmentationModel(0);
  }
  npx = segSPD->Npx();
  Double_t thr, sigma; 

  Int_t countRW = -1; // RS counter for raw -> cluster ID's (used in embedding)
  const TArrayI* rawID2clusID = fRawID2ClusID[kSPD];
  AliITSRawStreamSPD inputSPD(rawReader);
  while(1){
    Bool_t next  = inputSPD.Next();
    if (!next) break;

    countRW++; // RS
    Int_t module = inputSPD.GetModuleID();
    Int_t column = inputSPD.GetColumn();
    Int_t row    = inputSPD.GetRow();
    Int_t index  = npx * column + row;

    if (module >= size) continue;
 
    last = (fModA[module])->GetEntries();
    TClonesArray& dum = *fModA[module];
    fDetTypeSim->GetSimuParam()->SPDThresholds(module,thr,sigma);
    thr += 1.;
    int label = -1;
    if (rawID2clusID) { // RS If the raw->cluster ID is set (filled by cluster finder) store cluster ID's in SDigits
      if (rawID2clusID->GetSize()<=countRW) {AliError(Form("The buffer of rawSPD to clusSPD ID's is shorter than current rawSPD ID=%d",countRW));}
      else label = (*rawID2clusID)[countRW];
    }
    new (dum[last]) AliITSpListItem(label, -1, module, index, thr);
  }
  rawReader->Reset();

  //
  // SDD
  // 
  AliITSsegmentationSDD* segSDD = (AliITSsegmentationSDD*) fDetTypeSim->GetSegmentationModel(1);
  npx = segSDD->Npx();
  Int_t scalef=AliITSsimulationSDD::ScaleFourier(segSDD);
  Int_t firstSDD=AliITSgeomTGeo::GetModuleIndex(3,1,1);
  Int_t firstSSD=AliITSgeomTGeo::GetModuleIndex(5,1,1);
  //
  countRW = -1; // RS
  rawID2clusID = fRawID2ClusID[kSDD];
  AliITSRawStream* inputSDD=AliITSRawStreamSDD::CreateRawStreamSDD(rawReader);
  for(Int_t iMod=firstSDD; iMod<firstSSD; iMod++){
    AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)fDetTypeSim->GetCalibrationModel(iMod);
    Bool_t isZeroSupp=cal->GetZeroSupp();
    if(isZeroSupp){ 
      for(Int_t iSid=0; iSid<2; iSid++) inputSDD->SetZeroSuppLowThreshold(iMod-firstSDD,iSid,cal->GetZSLowThreshold(iSid));
    }else{
      for(Int_t iSid=0; iSid<2; iSid++) inputSDD->SetZeroSuppLowThreshold(iMod-firstSDD,iSid,0);
    }
  }

  AliITSDDLModuleMapSDD* ddlmap=fDetTypeSim->GetDDLModuleMapSDD();
  inputSDD->SetDDLModuleMap(ddlmap);
  while(inputSDD->Next()){
    countRW++; // RS
    if(inputSDD->IsCompletedModule()==kFALSE && 
       inputSDD->IsCompletedDDL()==kFALSE){

      Int_t module = inputSDD->GetModuleID();
      Int_t anode  = inputSDD->GetCoord1()+segSDD->NpzHalf()*inputSDD->GetChannel();
      Int_t time   = inputSDD->GetCoord2();
      Int_t signal10 = inputSDD->GetSignal();
      Int_t index = AliITSpList::GetIndex(anode,time,scalef*npx);

      if (module >= size) continue;
      last = fModA[module]->GetEntries();
      TClonesArray& dum = *fModA[module];
      int label = -1;
      if (rawID2clusID) { // RS If the raw->cluster ID is set (filled by cluster finder) store cluster ID's in SDigits
	if (rawID2clusID->GetSize()<=countRW) {AliError(Form("The buffer of rawSDD to clusSDD ID's is shorter than current rawSDD ID=%d",countRW));}
	else label = (*rawID2clusID)[countRW];
      }
      new (dum[last]) AliITSpListItem(label, -1, module, index, Double_t(signal10));
      ((AliITSpListItem*) dum.At(last))->AddSignalAfterElect(module, index, Double_t(signal10));
    }
  }
  delete inputSDD;
  rawReader->Reset();
    
  //
  // SSD
  // 
  AliITSsegmentationSSD* segSSD = (AliITSsegmentationSSD*) fDetTypeSim->GetSegmentationModel(2);
  npx = segSSD->Npx();
  AliITSRawStreamSSD inputSSD(rawReader);
  countRW = -1;
  rawID2clusID = fRawID2ClusID[kSSD];
  while(1){
    Bool_t next  = inputSSD.Next();
    if (!next) break;
    countRW++; // RS
    Int_t module  = inputSSD.GetModuleID();
    if(module<0)AliError(Form("Invalid SSD  module %d \n",module));
    if(module<0)continue;
    Int_t side    = inputSSD.GetSideFlag();
    Int_t strip   = inputSSD.GetStrip();
    Int_t signal  = inputSSD.GetSignal();
    Int_t index  = npx * side + strip;

    if (module >= size) continue;
	
    last = fModA[module]->GetEntries();
    TClonesArray& dum = *fModA[module];
    int label = -1;
    if (rawID2clusID) { // RS If the raw->cluster ID is set (filled by cluster finder) store cluster ID's in SDigits
      if (rawID2clusID->GetSize()<=countRW) {AliError(Form("The buffer of rawSSD to clusSSD ID's is shorter than current rawSSD ID=%d",countRW));}
      else label = (*rawID2clusID)[countRW];
    }    
    new (dum[last]) AliITSpListItem(label, -1, module, index, Double_t(signal));
  }
  rawReader->Reset();
  AliITSpListItem* sdig = 0;
    
  Int_t firstssd = GetITSgeom()->GetStartDet(kSSD);
  Double_t adcToEv = 1.;
  for (Int_t mod = 0; mod < size; mod++)
    {
      if(mod>=firstssd) {
	AliITSCalibrationSSD* calssd = (AliITSCalibrationSSD*)fDetTypeSim->GetCalibrationModel(mod);
	adcToEv = 1./calssd->GetSSDDEvToADC(1.);
      }
      Int_t nsdig =  fModA[mod]->GetEntries();
      for (Int_t ie = 0; ie < nsdig; ie++) {
	sdig = (AliITSpListItem*) (fModA[mod]->At(ie));
      	Double_t digsig = sdig->GetSignal();
	if(mod>=firstssd) digsig*=adcToEv; // for SSD: convert back charge from ADC to electron
	new (aSDigits[ie]) AliITSpListItem(sdig->GetTrack(0), -1, mod, sdig->GetIndex(), digsig);
	Float_t sig = sdig->GetSignalAfterElect();
	if(mod>=firstssd) sig*=adcToEv;
	if (sig > 0.) {
	  sdig = (AliITSpListItem*)aSDigits[ie];
	  sdig->AddSignalAfterElect(mod, sdig->GetIndex(), Double_t(sig));
	}
      }
	
      tree->Fill();
      aSDigits.Clear();
      fModA[mod]->Clear();
    }
  loader->WriteSDigits("OVERWRITE");    
  return kTRUE;
}

//______________________________________________________________________
void AliITS::UpdateInternalGeometry(){
  //reads new geometry from TGeo 
//   AliDebug(1,"Delete ITSgeom and create a new one reading TGeo");

  AliITSVersion_t version = (AliITSVersion_t)IsVersion();
  Int_t minor = 0;
  AliITSInitGeometry initgeom;
  AliITSgeom* geom = initgeom.CreateAliITSgeom(version,minor);
  SetITSgeom(geom);
}
//______________________________________________________________________
AliTriggerDetector* AliITS::CreateTriggerDetector() const {
  // create an AliITSTrigger object (and set trigger conditions as input)
  return new AliITSTrigger(fDetTypeSim->GetTriggerConditions());
}
//______________________________________________________________________
void AliITS::WriteFOSignals(){
// This method write FO signals in Digits tree both in Hits2Digits
// or SDigits2Digits

  fDetTypeSim->ProcessNoiseForFastOr();
  fDetTypeSim->WriteFOSignals();
}
