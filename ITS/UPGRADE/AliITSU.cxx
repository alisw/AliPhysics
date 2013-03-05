/*
*/

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

/* $Id: AliITSU.cxx $ */


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
//  AliITSU. Inner Traking System base class.
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
// AliITSU is the general base class for the ITS. Also see AliDetector for
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
#include "AliITSU.h"
#include "AliITSLoader.h"
#include "AliITSULoader.h"
#include "AliITSUHit.h"
#include "AliITSUSDigit.h"
#include "AliITSUSimulation.h"
#include "AliITSUSimulationPix.h"
#include "AliITSsimulationFastPoints.h"
#include "AliMC.h"
#include "AliITSUDigitizer.h"
#include "AliITSRecPoint.h"
#include "AliRawReader.h"
#include "AliRun.h"
#include "AliLog.h"
#include "AliITSdigit.h"
#include "AliITSUModule.h"
#include "AliITSUDigitPix.h"
#include "AliITSsegmentation.h"
#include "AliITSUSegmentationPix.h"
#include "AliITSUSimuParam.h"
#include "AliITSFOSignalsSPD.h"
#include "AliITSUParamList.h"
#include "AliCDBManager.h" // tmp! Later the simuparam should be loaded centrally
#include "AliCDBEntry.h"

ClassImp(AliITSU)

//______________________________________________________________________
AliITSU::AliITSU() : 
AliDetector()
  ,fEuclidOut(0)
  ,fNLayers(0)
  ,fIdSens(0)
  ,fLayerName(0)
  ,fTiming(kFALSE)
  ,fGeomTGeo(0)
  ,fSimuParam(0)
  ,fModA(0)
  ,fpSDigits(0)
  ,fSDigits(0)
  ,fDetHits(0)
  ,fModuleHits(0)
  ,fDetDigits(0)
  ,fSensMap(0)
  ,fSimModelLr(0)
  ,fSegModelLr(0)
  ,fResponseLr(0)
  ,fCalibration(0)
  ,fRunNumber(0)
  ,fSimInitDone(kFALSE)
{
  // Default initializer for ITS  
}

//______________________________________________________________________
AliITSU::AliITSU(const Char_t *title, Int_t nlay) :
  AliDetector("ITS",title)
  ,fEuclidOut(0)
  ,fNLayers(nlay)
  ,fIdSens(0)
  ,fLayerName(0)
  ,fTiming(kFALSE)
  ,fGeomTGeo(0)
  ,fSimuParam(0)
  ,fModA(0)
  ,fpSDigits(0)
  ,fSDigits(0)
  ,fDetHits(0)
  ,fModuleHits(0)
  ,fDetDigits(0)
  ,fSensMap(0)
  ,fSimModelLr(0)
  ,fSegModelLr(0)
  ,fResponseLr(0)
  ,fCalibration(0)
  ,fRunNumber(0)
  ,fSimInitDone(kFALSE)
{
  //     The standard Constructor for the ITS class. 
  AliMC* mc = gAlice->GetMCApp();
  if( mc && mc->GetHitLists() ) {
    fHits = new TClonesArray("AliITSUHit",100); // from AliDetector
    mc->AddHitList(fHits);  
  }
}


//______________________________________________________________________
AliITSU::~AliITSU()
{
  // Default destructor for ITS.
  //  
  delete fHits;
  //  delete fSimuParam; // provided by the CDBManager
  delete fSensMap;
  if (fSimModelLr) {
    for (int i=fNLayers;i--;) { // different layers may use the same simulation model
      for (int j=i;j--;) if (fSimModelLr[j]==fSimModelLr[i]) fSimModelLr[j] = 0;
      delete fSimModelLr[i];
    }
    delete[] fSimModelLr;
  }
  if (fSegModelLr) {
    for (int i=fNLayers;i--;) { // different layers may use the same simulation model     
      for (int j=i;j--;) if (fSegModelLr[j]==fSegModelLr[i]) fSegModelLr[j] = 0;
      delete fSegModelLr[i];
    }
    delete[] fSegModelLr;
  }
  //
  delete fResponseLr; // note: the response data is owned by the CDBManager, we don't delete them
  //
  delete[] fLayerName;  // Array of TStrings
  delete[] fIdSens;
  //
  int nmod = fGeomTGeo ? fGeomTGeo->GetNModules() : 0;
  if (fModuleHits) fModuleHits->Delete();
  delete fModuleHits;
  if(fModA){
    for(Int_t j=0; j<nmod; j++){
      fModA[j]->Delete();
      delete fModA[j];
    }
    delete[] fModA;
  }
  if (fpSDigits) { fpSDigits->Delete(); delete fpSDigits; }
  if (fSDigits)  { fSDigits->Delete(); delete fSDigits; }
  delete fGeomTGeo;
  //
}

//______________________________________________________________________
AliDigitizer* AliITSU::CreateDigitizer(AliDigitizationInput* manager) const
{
  // Creates the AliITSDigitizer in a standard way for use via AliModule.
  // This function can not be included in the .h file because of problems
  // with the order of inclusion (recursive).
  // Inputs:
  //    AliDigitizationInput *manager  The Manger class for Digitization
  // Output:
  //    none.
  // Return:
  //    A new AliITSRunDigitizer (cast as a AliDigitizer).
  //
  return new AliITSUDigitizer(manager);
}

//______________________________________________________________________
void AliITSU::Init() 
{
  // Initializer ITS after it has been built
  //     This routine initializes the AliITS class. It is intended to be
  // called from the Init function in AliITSv?. Besides displaying a banner
  // indicating that it has been called it initializes the array fIdSens
  // and sets the default segmentation, response, digit and raw cluster
  // classes therefore it should be called after a call to CreateGeometry.
  //
  if (!fIdSens) fIdSens = new Int_t[fNLayers];
  for(int i=0;i<fNLayers;i++) fIdSens[i] = gMC ? gMC->VolId(fLayerName[i]) : 0;
  fGeomTGeo     = new AliITSUGeomTGeo(kTRUE);
  InitSimulation();
  //
}

//______________________________________________________________________
void AliITSU::MakeBranch(Option_t* option)
{
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
  
  Bool_t cH = (strstr(option,"H")!=0);
  Bool_t cS = (strstr(option,"S")!=0);
  Bool_t cD = (strstr(option,"D")!=0);
  
  if(cH && (fHits == 0x0)) fHits = new TClonesArray("AliITSUHit", 1560);
  AliDetector::MakeBranch(option);
  
  if(cS) MakeBranchS(0);
  if(cD) MakeBranchD(0);
  // 
}

//___________________________________________________________________
void AliITSU::MakeBranchS(const char* fl) 
{
  // Creates Tree Branch for the ITS summable digits.
  // Inputs:
  //      cont char *fl  File name where SDigits branch is to be written
  //                     to. If blank it write the SDigits to the same
  //                     file in which the Hits were found.
  //  
  Int_t buffersize = 4000;
  char branchname[31];
  //
  // only one branch for SDigits.
  snprintf(branchname,30,"%s",GetName());
  if(fLoader->TreeS()) MakeBranchInTree(fLoader->TreeS(),branchname,&fSDigits,buffersize,fl);
  //
}

//______________________________________________________________________
void AliITSU::MakeBranchD(const char* file)
{
  //Make branch for digits
  MakeBranchInTreeD(fLoader->TreeD(),file);
}

//___________________________________________________________________
void AliITSU:: MakeBranchInTreeD(TTree* treeD, const char* file)
{
  // Creates Tree branches for the ITS.
  //
  if (!treeD) {AliFatal("No tree provided");}
  Int_t buffersize = 4000;
  if (!fDetDigits) InitArrays();
  //
  for (Int_t i=0;i<kNDetTypes;i++) {
    ResetDigits(i);
    TClonesArray* darr = (TClonesArray*)fDetDigits->At(i);
    AliDetector::MakeBranchInTree(treeD,Form("%sDigits%s",GetName(),fGeomTGeo->GetDetTypeName(i)),
				  &darr,buffersize,file);
  }
  //
}

//______________________________________________________________________
void AliITSU::InitArrays()
{
  // initialize arrays
  //
  if(!fLoader) MakeLoader(AliConfig::GetDefaultEventFolderName());
  //  
  fDetDigits = new TObjArray(kNDetTypes);
  for (Int_t i=0;i<kNDetTypes;i++) fDetDigits->AddAt(new TClonesArray(GetDigitClassName(i),100),i);
  //
  fSDigits = new TClonesArray("AliITSUSDigit",100);
  //
  fDetHits = new TClonesArray("AliITSUHit",100);
  //
  fModuleHits = new TObjArray(fGeomTGeo->GetNModules());
  for (int i=0;i<fGeomTGeo->GetNModules();i++) fModuleHits->AddLast( new AliITSUModule(i,fGeomTGeo) );
  //
}

//______________________________________________________________________
void AliITSU::SetTreeAddress()
{
  // Set branch address for the Trees.
  TTree *treeS = fLoader->TreeS();
  if (treeS) {
    TBranch* br = treeS->GetBranch(GetName());
    if (br) br->SetAddress(&fSDigits);
  }
  //
  TTree *treeD = fLoader->TreeD();
  if (treeD) {
    if (!fDetDigits) InitArrays();
    for (int i=0;i<kNDetTypes;i++) {
      TString brname = Form("%sDigits%s",GetName(),GetDetTypeName(i));
      TBranch* br = treeD->GetBranch(brname.Data());
      if (!br) continue;
      TClonesArray* darr = (TClonesArray*)fDetDigits->At(i);
      br->SetAddress(&darr);
    }
  }
  if (fLoader->TreeH() && (fHits == 0x0)) fHits = new TClonesArray("AliITSUHit", 1560);
  AliDetector::SetTreeAddress();
  //
}

//______________________________________________________________________
void AliITSU::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  // Add an ITS hit
  //     The function to add information to the AliITSUHit class. See the
  // AliITSUHit class for a full description. This function allocates the
  // necessary new space for the hit information and passes the variable
  // track, and the pointers *vol and *hits to the AliITSUHit constructor
  // function.
  // Inputs:
  //      Int_t   track   Track number which produced this hit.
  //      Int_t   *vol    Array of Integer Hit information. See AliITSUHit.h
  //      Float_t *hits   Array of Floating Hit information.  see AliITSUHit.h
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliITSUHit(fIshunt,track,vol,hits);
  //
}

//______________________________________________________________________
void AliITSU::FillModules(Int_t bgrev, Option_t *option, const char *filename) 
{
  // fill the modules with the sorted by module hits; add hits from
  // background if option=Add.
  //
  static TTree *trH1=0;                 //Tree with background hits
  static Bool_t first=kTRUE;
  static TFile *file = 0;
  const char *addBgr = strstr(option,"Add");
  //
  if (addBgr ) {
    if(first) {
      file = new TFile(filename);    
      first=kFALSE;
    }
    file->cd();
    file->ls();
    // Get Hits Tree header from file
    if (trH1) {delete trH1; trH1=0;}
    //
    char treeName[21];
    snprintf(treeName,20,"TreeH%d",bgrev);
    trH1 = (TTree*)gDirectory->Get(treeName);
    if (!trH1) Error("FillModules","cannot find Hits Tree for event:%d",bgrev);
    // Set branch addresses
  } // end if addBgr
  
  FillModules(fLoader->TreeH(),0); // fill from this file's tree.
  //
  if (addBgr ) {
    FillModules(trH1,10000000); // Default mask 10M.
    TTree *fAli=fLoader->GetRunLoader()->TreeK();
    TFile *fileAli=0;
    if (fAli) {
      fileAli = fAli->GetCurrentFile();
      fileAli->cd();
    }
  } // end if add
  //  
}

//______________________________________________________________________
void AliITSU::FillModules(TTree *treeH, Int_t /*mask*/)
{
  // fill the modules with the sorted by module hits; 
  // can be called many times to do a merging
  // Inputs:
  //      TTree *treeH  The tree containing the hits to be copied into
  //                    the modules.
  //      Int_t mask    The track number mask to indecate which file
  //                    this hits came from.
  //  
  if (treeH == 0x0) { AliError("Tree H  is NULL"); return; }
  //
  Int_t lay,lad,det,index;
  AliITSUHit *itsHit=0;
  char branchname[21];
  snprintf(branchname,20,"%s",GetName());
  TBranch *branch = treeH->GetBranch(branchname);
  if (!branch) {Error("FillModules","%s branch in TreeH not found",branchname); return;} // end if !branch
  //
  branch->SetAddress(&fHits);
  Int_t nTracks =(Int_t) treeH->GetEntries();
  Int_t iPrimTrack,h;
  for (iPrimTrack=0; iPrimTrack<nTracks; iPrimTrack++) {
    ResetHits();
    Int_t nBytes = treeH->GetEvent(iPrimTrack);
    if (nBytes <= 0) continue;
    Int_t nHits = fHits->GetEntriesFast();
    for (h=0; h<nHits; h++){
      itsHit = (AliITSUHit *)fHits->UncheckedAt(h);
      itsHit->GetDetectorID(lay,lad,det);
      index = fGeomTGeo->GetModuleIndex(lay,lad,det); // !!! AliITSHit counts indices from 1!
      itsHit = new( (*fDetHits)[fDetHits->GetEntriesFast()] ) AliITSUHit(*itsHit);
      itsHit->SetUniqueID(h);
      GetModule(index)->AddHit(itsHit);
      // do we need to add a mask?
      // itsHit->SetTrack(itsHit->GetTrack()+mask);
    } // end loop over hits 
  } // end loop over tracks
}

//______________________________________________________________________
void AliITSU::ClearModules()
{
  // clear accumulated hits
  if (!fModuleHits || !fDetHits) AliFatal("Hits accumulation arrays are not defined");
  for (int i=fGeomTGeo->GetNModules();i--;) GetModule(i)->Clear();
  fDetHits->Clear();
}

//______________________________________________________________________
void AliITSU::Hits2SDigits()
{
  // Standard Hits to summable Digits function.
  if (!IsSimInitDone()) InitSimulation();
  fLoader->LoadHits("read");
  fLoader->LoadSDigits("recreate");
  AliRunLoader* rl = fLoader->GetRunLoader(); 
  //
  for (Int_t iEvent = 0; iEvent < rl->GetNumberOfEvents(); iEvent++) {
    rl->GetEvent(iEvent);
    if (!fLoader->TreeS()) fLoader->MakeTree("S");
    MakeBranch("S");
    SetTreeAddress();
    Hits2SDigits(iEvent,0," "," ");
  } // end for iEvent
    //
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  // 
}

//______________________________________________________________________
void AliITSU::Hits2SDigits(Int_t evNumber,Int_t bgrev,Option_t *option,const char *filename)
{
  // Keep galice.root for signal and name differently the file for 
  // background when add! otherwise the track info for signal will be lost !
  // Inputs:
  //      Int_t evnt       Event to be processed.
  //      Int_t bgrev      Background Hit tree number.
  //      Int_t nmodules   Not used.
  //      Option_t *option String indicating if merging hits or not. To
  //                       merge hits set equal to "Add". Otherwise no
  //                       background hits are considered.
  //      Test_t *filename File name containing the background hits..
  //
  if (!IsSimInitDone()) InitSimulation();
  FillModules(bgrev,option,filename);
  //
  Int_t nmodules = fGeomTGeo->GetNModules();
  
  for(int module=0;module<nmodules;module++) {
    int lr = fGeomTGeo->GetLayer(module);
    AliITSUSimulation* sim = GetSimulationModel(lr);
    sim->InitSimulationModule(GetModule(module),evNumber/*,gAlice->GetEvNumber()*/,GetSegmentation(lr),GetResponseParam(lr));
    sim->SDigitiseModule();
    fLoader->TreeS()->Fill();      // fills all branches - wasted disk space
    ResetSDigits();
  } 
  //
  ClearModules();
  //
  fLoader->TreeS()->GetEntries();
  fLoader->TreeS()->AutoSave();
  fLoader->WriteSDigits("OVERWRITE");
  fLoader->TreeS()->Reset();
}

//______________________________________________________________________
void AliITSU::Hits2Digits()
{
  // Standard Hits to Digits function.
  if (!IsSimInitDone()) InitSimulation();
  fLoader->LoadHits("read");
  fLoader->LoadDigits("recreate");
  AliRunLoader* rl = fLoader->GetRunLoader(); 
  //
  for (Int_t iEvent = 0; iEvent < rl->GetNumberOfEvents(); iEvent++) {
    rl->GetEvent(iEvent);
    if (!fLoader->TreeS()) fLoader->MakeTree("S");
    MakeBranch("D");
    SetTreeAddress();
    Hits2Digits(iEvent,0," "," ");
  } // end for iEvent
    //
  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  // 
}

//______________________________________________________________________
void AliITSU::Hits2Digits(Int_t evNumber,Int_t bgrev,Option_t *option,const char *filename)
{
  //   Keep galice.root for signal and name differently the file for 
  // background when add! otherwise the track info for signal will be lost !
  // Inputs:
  //      Int_t evnt       Event to be processed.
  //      Int_t bgrev      Background Hit tree number.
  //      Option_t *option String indicating if merging hits or not. To
  //                       merge hits set equal to "Add". Otherwise no
  //                       background hits are considered.
  //      Test_t *filename File name containing the background hits..
  // Outputs:
  //  
  if (!IsSimInitDone()) InitSimulation();
  FillModules(bgrev,option,filename);
  // 
  Int_t nmodules = fGeomTGeo->GetNModules();
  for (Int_t module=0;module<nmodules;module++) {
    int lr = fGeomTGeo->GetLayer(module);
    AliITSUSimulation* sim = GetSimulationModel(lr);
    sim->InitSimulationModule(GetModule(module),evNumber/*gAlice->GetEvNumber()*/,GetSegmentation(lr),GetResponseParam(lr));
    sim->DigitiseModule();
    // fills all branches - wasted disk space
    fLoader->TreeD()->Fill(); 
    ResetDigits();
  } // end for module
  //
  ClearModules();
  //
  //    WriteFOSignals(); // Add Fast-OR signals to event (only one object per event)
  fLoader->TreeD()->GetEntries();
  fLoader->TreeD()->AutoSave();
  fLoader->TreeD()->Reset(); 
  //
}

//_____________________________________________________________________
void AliITSU::Hits2FastRecPoints(Int_t bgrev,Option_t *opt,const char *flnm)
{
  // keep galice.root for signal and name differently the file for 
  // background when add! otherwise the track info for signal will be lost !
  // Inputs:
  //      Int_t evnt       Event to be processed.
  //      Int_t bgrev      Background Hit tree number.
  //      Option_t *opt    Option passed to FillModules. See FillModules.
  //      Test_t *flnm     File name containing the background hits..
  // Outputs:
  //      none.
  // Return:
  //      none.
  if (!IsSimInitDone()) InitSimulation();
  AliITSULoader *pITSloader = (AliITSULoader*)fLoader;
  Int_t nmodules = fGeomTGeo->GetNModules();
  FillModules(bgrev,opt,flnm);
  //
  TTree *lTR = pITSloader->TreeR();
  if(!lTR) {
    pITSloader->MakeTree("R");
    lTR = pITSloader->TreeR();
  }
  //
  TClonesArray* ptarray = new TClonesArray("AliITSRecPoint",1000);
  TBranch* branch = (TBranch*)lTR->Branch("ITSRecPointsF",&ptarray);
  branch->SetAddress(&ptarray);
  for (int module=0;module<nmodules;module++){
    int id = fGeomTGeo->GetModuleDetTypeID(module);
    AliITSUSimulation* sim = GetSimulationModel(id);
    if (!sim) AliFatal(Form("The sim.class for module %d of DetTypeID %d is missing",module,id));
    sim->CreateFastRecPoints( GetModule(module) ,module,gRandom,ptarray);
    lTR->Fill();
    ptarray->Clear();
  } // end for module
  //
  ClearModules();
  fLoader->WriteRecPoints("OVERWRITE");
  delete ptarray;
}

//_____________________________________________________________________
Int_t AliITSU::Hits2Clusters(TTree */*hTree*/, TTree */*cTree*/)
{
  /* RS: TODO
  // This function creates ITS clusters
  if (!IsSimInitDone()) InitSimulation();
  Int_t mmax = 0;
  FillModules(hTree,0);
  //
  TClonesArray *points = new TClonesArray("AliITSRecPoint",1000);
  TBranch *branch=cTree->GetBranch("ITSRecPoints");
  if (!branch) cTree->Branch("ITSRecPoints",&points);
  else branch->SetAddress(&points);
  //
  AliITSsimulationFastPoints sim;
  Int_t ncl=0;
  for (Int_t m=0; m<mmax; m++) {
    sim.CreateFastRecPoints(GetModule(m),m,gRandom,points);      
    ncl+=points->GetEntriesFast();
    cTree->Fill();
    points->Clear();
  }
  //
  ClearModules();
  //
  AliDebug(1,Form("Number of found fast clusters : %d",ncl));
  //cTree->Write();
  delete points;
  */
  return 0;
}

//_____________________________________________________________________
void AliITSU::CheckLabels(Int_t lab[3]) const //RSDONE
{
  // Tries to find mother's labels
  //
  if(lab[0]<0 && lab[1]<0 && lab[2]<0) return; // In case of no labels just exit
  //
  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
  for (Int_t i=0;i<3;i++){
    Int_t label = lab[i];
    if (label>=0 && label<ntracks) {
      TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
      if (part->P() < 0.005) {
	Int_t m=part->GetFirstMother();
	if (m<0) continue;
	if (part->GetStatusCode()>0) continue;
	lab[i]=m;       
      }
    }    
  }
  //
}

//______________________________________________________________________
void AliITSU::ResetDigits() //RSDONE?
{
  // Reset number of digits and the digits array for the ITS detector.
  if (fDetDigits) for (int i=kNDetTypes;i--;) ResetDigits(i);
  //
}

//______________________________________________________________________
void AliITSU::ResetDigits(Int_t branch)
{
  // Reset number of digits and the digits array for this branch.
  if (fDetDigits) ((TClonesArray*)fDetDigits->At(branch))->Clear();
  //
}

//______________________________________________________________________
void AliITSU::AddSumDigit(AliITSUSDigit &sdig)
{
  // Adds the module summable digits to the summable digits tree.
  new( (*fSDigits)[fSDigits->GetEntriesFast()]) AliITSUSDigit(sdig);
  //  
}

//______________________________________________________________________
void AliITSU::AddSimDigit(Int_t branch, AliITSdigit *d)
{
  //    Add a simulated digit.
  // Inputs:
  //      Int_t id        Detector type number.
  //      AliITSdigit *d  Digit to be added to the Digits Tree. See 
  //                      AliITSdigit.h
  TClonesArray &ldigits = *((TClonesArray*)fDetDigits->At(branch));
  int nd = ldigits.GetEntriesFast();
  switch(branch){
  case AliITSUGeomTGeo::kDetTypePix:
    new(ldigits[nd]) AliITSUDigitPix(*((AliITSUDigitPix*)d));
    break;
  default:
    AliFatal(Form("Unknown digits branch %d",branch));
  }
}

//______________________________________________________________________
void AliITSU::AddSimDigit(Int_t branch,Float_t /*phys*/,Int_t *digits,Int_t *tracks,
			  Int_t *hits,Float_t */*charges*/, Int_t /*sigexpanded*/)
{
  // Add a simulated digit to the list.
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
  TClonesArray &ldigits = *((TClonesArray*)fDetDigits->At(branch));
  int nd = ldigits.GetEntriesFast();
  switch(branch){
  case AliITSUGeomTGeo::kDetTypePix:
    new(ldigits[nd]) AliITSUDigitPix(digits,tracks,hits);
    break;
  default:
    AliFatal(Form("Unknown digits branch %d",branch));
  }  
  //
}

//______________________________________________________________________
void AliITSU::Digits2Raw()
{
  AliError("Not ready");
}

//______________________________________________________________________
AliLoader* AliITSU::MakeLoader(const char* topfoldername)
{ 
  //builds ITSgetter (AliLoader type)
  //if detector wants to use castomized getter, it must overload this method
  
  AliDebug(1,Form("Creating AliITSULoader. Top folder is %s.",topfoldername));
  fLoader = new AliITSULoader(GetName(),topfoldername);
  return fLoader;
}

//______________________________________________________________________
Bool_t AliITSU::Raw2SDigits(AliRawReader* /*rawReader*/)
{
  AliError("Not ready");
  return kFALSE;
}

//______________________________________________________________________
/*
AliTriggerDetector* AliITSU::CreateTriggerDetector() const 
{
  // create an AliITSTrigger object (and set trigger conditions as input)
  return new AliITSTrigger(fDetTypeSim->GetTriggerConditions());
}
*/

//______________________________________________________________________
void AliITSU::WriteFOSignals()
{
  // This method write FO signals in Digits tree both in Hits2Digits
  // or SDigits2Digits
  AliError("Not ready");
  //  fDetTypeSim->ProcessNoiseForFastOr();
}

//_______________________________________________________________________
void AliITSU::SDigits2Digits()
{
  // Standard Summable digits to Digits function.
  //
  if (!IsSimInitDone()) InitSimulation();
  TTree* trees = fLoader->TreeS();
  if( !(trees && fSDigits) ) AliFatal("Error: No trees or SDigits.");
  TBranch* brchSDigits = trees->GetBranch(GetName());
  //
  int nmodules = fGeomTGeo->GetNModules();
  for (int module=0;module<nmodules;module++) {
    int lr = fGeomTGeo->GetLayer(module);
    AliITSUSimulation* sim = GetSimulationModel(lr);
    sim->InitSimulationModule(GetModule(module),gAlice->GetEvNumber(),GetSegmentation(lr),GetResponseParam(lr));
    fSDigits->Clear();
    brchSDigits->GetEvent(module);
    sim->AddSDigitsToModule(fSDigits,0);
    sim->FinishSDigitiseModule();
    fLoader->TreeD()->Fill();
    ResetDigits();
  }
  //  WriteFOSignals(); 
  fLoader->TreeD()->GetEntries();
  fLoader->TreeD()->AutoSave();
  fLoader->TreeD()->Reset();
}

//_______________________________________________________________________
void AliITSU::InitSimulation()
{
  // Initialize arrays, segmentations ets, needed for simulation
  // Equivalent of old AliITSDetTypeSim construction
  //
  if (fSimInitDone) {AliInfo("Already done"); return;}
  //
  AliCDBEntry* cdbEnt = AliCDBManager::Instance()->Get("ITS/Calib/SimuParam"); // tmp: load it centrally
  if (!cdbEnt) {AliFatal("Failed to find ITS/Calib/SimuParam on CDB"); exit(1);}
  fSimuParam    = (AliITSUSimuParam*)cdbEnt->GetObject();
  //
  fSensMap      = new AliITSUSensMap("AliITSUSDigit",0,0);
  fSimModelLr   = new AliITSUSimulation*[fNLayers];
  fSegModelLr   = new AliITSsegmentation*[fNLayers];
  fResponseLr   = new AliITSUParamList*[fNLayers];
  //
  TObjArray arrSeg;
  AliITSUSegmentationPix::LoadSegmentations(&arrSeg, AliITSUGeomTGeo::GetITSsegmentationFileName());
  //
  // add known simulation types used in the setup
  for (int i=fNLayers;i--;) {
    fSimModelLr[i] = 0;
    fSegModelLr[i] = 0;
    fResponseLr[i] = 0;
    int dType = fGeomTGeo->GetLayerDetTypeID(i);           // fine detector type: class + segmentation
    int sType = dType/AliITSUGeomTGeo::kMaxSegmPerDetType; // detector simulation class
    //
    // check if the simulation of this sType was already created for preceeding layers
    AliITSUSimulation* simUpg = 0;
    for (int j=fNLayers-1;j>i;j--) {
      simUpg = GetSimulationModel(j);
      if (simUpg && int(simUpg->GetUniqueID())==sType) break;
      else simUpg = 0;
    }
    //
    if (!simUpg) { // need to create simulation for detector class sType
      switch (sType) 
	{
	case AliITSUGeomTGeo::kDetTypePix : 
	  simUpg = new AliITSUSimulationPix(fSimuParam,fSensMap); 
	  break;
	default: AliFatal(Form("No %d detector type is defined",sType));
	}
    }
    fSimModelLr[i] = simUpg;
    //
    // add segmentations used in the setup
    if (!(fSegModelLr[i]=(AliITSsegmentation*)arrSeg[dType])) {AliFatal(Form("Segmentation for DetType#%d is not found",dType)); exit(1);}
    //
    // add response function for the detectors of this layer
    if ( !(fResponseLr[i]=(AliITSUParamList*)fSimuParam->FindRespFunParams(dType)) ) {AliFatal(Form("Response for DetType#%d is not found in SimuParams",dType)); exit(1);}
  }
  // delete non needed segmentations
  for (int i=fNLayers;i--;) arrSeg.Remove(fSegModelLr[i]);
  arrSeg.Delete();
  //
  InitArrays();
  //
  fSimInitDone = kTRUE;
  //
}
