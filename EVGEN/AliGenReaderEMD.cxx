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


// Class to read events from external (TNtupla) file
// Events -> neutron removal by EM dissociation of Pb nuclei
// Data from RELDIS code (by I. Pshenichov)

#include <TFile.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include "AliGenReaderEMD.h"
#include "AliStack.h"


ClassImp(AliGenReaderEMD)

AliGenReaderEMD::AliGenReaderEMD():
    fStartEvent(0),
    fNcurrent(0),  
    fNparticle(0), 
    fTreeNtuple(0),
    fPcToTrack(0),
    fOffset(0),
    fNnAside(0),
    fEnAside(0),
    fnPDGCode(0),
    fNnCside(0),
    fEnCside(0),
    fNpAside(0),
    fEtapAside(0),
    fpPDGCode(0),
    fNpCside(0),
    fEtapCside(0),
    fNppAside(0),
    fEtappAside(0),
    fppPDGCode(0),
    fNppCside(0),
    fEtappCside(0),
    fNpmAside(0),
    fEtapmAside(0),
    fpmPDGCode(0),
    fNpmCside(0),
    fEtapmCside(0),
    fNp0Aside(0),
    fEtap0Aside(0),
    fp0PDGCode(0),
    fNp0Cside(0),
    fEtap0Cside(0),
    fNetaAside(0),
    fEtaetaAside(0),
    fetaPDGCode(0),
    fNetaCside(0),
    fEtaetaCside(0),
    fNomegaAside(0),
    fEtaomegaAside(0),
    fomegaPDGCode(0),
    fNomegaCside(0),
    fEtaomegaCside(0)
{
// Std constructor
    for(int i=0; i<70; i++){
       fPxnAside[i] = fPynAside[i] = fPznAside[i] = 0.;
       fPxnCside[i] = fPynCside[i] = fPznCside[i] = 0.;
       if(i<50){
         fPxpAside[i] = fPypAside[i] = fPzpAside[i] = 0.;
      	 fPxpCside[i] = fPypCside[i] = fPzpCside[i] = 0.;
         if(i<30){
           fPxppAside[i] = fPyppAside[i] = fPzppAside[i] = 0.;
      	   fPxppCside[i] = fPyppCside[i] = fPzppCside[i] = 0.;
           fPxpmAside[i] = fPypmAside[i] = fPzpmAside[i] = 0.;
      	   fPxpmCside[i] = fPypmCside[i] = fPzpmCside[i] = 0.;
           fPxp0Aside[i] = fPyp0Aside[i] = fPzp0Aside[i] = 0.;
      	   fPxp0Cside[i] = fPyp0Cside[i] = fPzp0Cside[i] = 0.;
	   if(i<15){
             fPxetaAside[i] = fPyetaAside[i] = fPzetaAside[i] = 0.;
      	     fPxetaCside[i] = fPyetaCside[i] = fPzetaCside[i] = 0.;
             fPxomegaAside[i] = fPyomegaAside[i] = fPzomegaAside[i] = 0.;
      	     fPxomegaCside[i] = fPyomegaCside[i] = fPzomegaCside[i] = 0.;
	   }
	 }
       }	
    }
    if(fPcToTrack==kAll) printf("\n\t   *** AliGenReaderEMD will track all produced particles \n\n");
    else if(fPcToTrack==kNotNucleons) printf("\n\t   *** AliGenReaderEMD will track all produced particles except nucleons\n\n");
    else if(fPcToTrack==kOnlyNucleons) printf("\n\t   *** AliGenReaderEMD will track only nucleons\n\n");
}


AliGenReaderEMD::AliGenReaderEMD(const AliGenReaderEMD &reader):
    AliGenReader(reader),
    fStartEvent(0),
    fNcurrent(0),  
    fNparticle(0), 
    fTreeNtuple(0),
    fPcToTrack(0),
    fOffset(0),
    fNnAside(0),
    fEnAside(0),
    fnPDGCode(0),
    fNnCside(0),
    fEnCside(0),
    fNpAside(0),
    fEtapAside(0),
    fpPDGCode(0),
    fNpCside(0),
    fEtapCside(0),
    fNppAside(0),
    fEtappAside(0),
    fppPDGCode(0),
    fNppCside(0),
    fEtappCside(0),
    fNpmAside(0),
    fEtapmAside(0),
    fpmPDGCode(0),
    fNpmCside(0),
    fEtapmCside(0),
    fNp0Aside(0),
    fEtap0Aside(0),
    fp0PDGCode(0),
    fNp0Cside(0),
    fEtap0Cside(0),
    fNetaAside(0),
    fEtaetaAside(0),
    fetaPDGCode(0),
    fNetaCside(0),
    fEtaetaCside(0),
    fNomegaAside(0),
    fEtaomegaAside(0),
    fomegaPDGCode(0),
    fNomegaCside(0),
    fEtaomegaCside(0)
{
    // Copy Constructor
    for(int i=0; i<70; i++){
       fPxnAside[i] = fPynAside[i] = fPznAside[i] = 0.;
       fPxnCside[i] = fPynCside[i] = fPznCside[i] = 0.;
       if(i<50){
         fPxpAside[i] = fPypAside[i] = fPzpAside[i] = 0.;
      	 fPxpCside[i] = fPypCside[i] = fPzpCside[i] = 0.;
         if(i<30){
           fPxppAside[i] = fPyppAside[i] = fPzppAside[i] = 0.;
      	   fPxppCside[i] = fPyppCside[i] = fPzppCside[i] = 0.;
           fPxpmAside[i] = fPypmAside[i] = fPzpmAside[i] = 0.;
      	   fPxpmCside[i] = fPypmCside[i] = fPzpmCside[i] = 0.;
           fPxp0Aside[i] = fPyp0Aside[i] = fPzp0Aside[i] = 0.;
      	   fPxp0Cside[i] = fPyp0Cside[i] = fPzp0Cside[i] = 0.;
	   if(i<15){
             fPxetaAside[i] = fPyetaAside[i] = fPzetaAside[i] = 0.;
      	     fPxetaCside[i] = fPyetaCside[i] = fPzetaCside[i] = 0.;
             fPxomegaAside[i] = fPyomegaAside[i] = fPzomegaAside[i] = 0.;
      	     fPxomegaCside[i] = fPyomegaCside[i] = fPzomegaCside[i] = 0.;
	   }
	 }
       }	
    }
    reader.Copy(*this);
}
  // -----------------------------------------------------------------------------------
AliGenReaderEMD::~AliGenReaderEMD()
{
    delete fTreeNtuple;
}

// -----------------------------------------------------------------------------------
AliGenReaderEMD& AliGenReaderEMD::operator=(const  AliGenReaderEMD& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

// -----------------------------------------------------------------------------------
void AliGenReaderEMD::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}

// -----------------------------------------------------------------------------------
void AliGenReaderEMD::Init() 
{
//
// Reset the existing file environment and open a new root file
    
    TFile *pFile=0;
    if (!pFile) {
	pFile = new TFile(fFileName);
	pFile->cd();
	printf("\n %s file opened to read RELDIS EMD events\n\n", fFileName);
    }
    fTreeNtuple = (TTree*)gDirectory->Get("h2032");
    fNcurrent = fStartEvent;

    TTree *Ntu=fTreeNtuple;
    //
    // Set branch addresses
    // **** neutrons
    Ntu->SetBranchAddress("Nleft",&fNnAside);
    Ntu->SetBranchAddress("Eleft",&fEnAside);
    Ntu->SetBranchAddress("Ipdg_l_n",&fnPDGCode);
    Ntu->SetBranchAddress("Pxl",  fPxnAside);
    Ntu->SetBranchAddress("Pyl",  fPynAside);
    Ntu->SetBranchAddress("Pzl",  fPznAside);
    Ntu->SetBranchAddress("Nright",&fNnCside);
    Ntu->SetBranchAddress("Eright",&fEnCside);
    Ntu->SetBranchAddress("Pxr",   fPxnCside);
    Ntu->SetBranchAddress("Pyr",   fPynCside);
    Ntu->SetBranchAddress("Pzr",   fPznCside);
    // **** protons
    Ntu->SetBranchAddress("Nleft_p",&fNpAside);
    Ntu->SetBranchAddress("Etaleft_p",&fEtapAside);
    Ntu->SetBranchAddress("Ipdg_l_p",&fpPDGCode);
    Ntu->SetBranchAddress("Pxl_p",  fPxpAside);
    Ntu->SetBranchAddress("Pyl_p",  fPypAside);
    Ntu->SetBranchAddress("Pzl_p",  fPzpAside);
    Ntu->SetBranchAddress("Nright_p",&fNpCside);
    Ntu->SetBranchAddress("Etaright_p",&fEtapCside);
    Ntu->SetBranchAddress("Pxr_p",   fPxpCside);
    Ntu->SetBranchAddress("Pyr_p",   fPypCside);
    Ntu->SetBranchAddress("Pzr_p",   fPzpCside);
    // **** pi+
    Ntu->SetBranchAddress("Nleft_pp",&fNppAside);
    Ntu->SetBranchAddress("Etaleft_pp",&fEtappAside);
    Ntu->SetBranchAddress("Ipdg_l_pp",&fppPDGCode);
    Ntu->SetBranchAddress("Pxl_pp",  fPxppAside);
    Ntu->SetBranchAddress("Pyl_pp",  fPyppAside);
    Ntu->SetBranchAddress("Pzl_pp",  fPzppAside);
    Ntu->SetBranchAddress("Nright_pp",&fNppCside);
    Ntu->SetBranchAddress("Etaright_pp",&fEtappCside);
    Ntu->SetBranchAddress("Pxr_pp",   fPxppCside);
    Ntu->SetBranchAddress("Pyr_pp",   fPyppCside);
    Ntu->SetBranchAddress("Pzr_pp",   fPzppCside);
    // **** pi-
    Ntu->SetBranchAddress("Nleft_pm",&fNpmAside);
    Ntu->SetBranchAddress("Etaleft_pm",&fEtapmAside);
    Ntu->SetBranchAddress("Ipdg_l_pm",&fpmPDGCode);
    Ntu->SetBranchAddress("Pxl_pm",  fPxpmAside);
    Ntu->SetBranchAddress("Pyl_pm",  fPypmAside);
    Ntu->SetBranchAddress("Pzl_pm",  fPzpmAside);
    Ntu->SetBranchAddress("Nright_pm",&fNpmCside);
    Ntu->SetBranchAddress("Etaright_pm",&fEtapmCside);
    Ntu->SetBranchAddress("Pxr_pm",   fPxpmCside);
    Ntu->SetBranchAddress("Pyr_pm",   fPypmCside);
    Ntu->SetBranchAddress("Pzr_pm",   fPzpmCside);
    // **** pi0
    Ntu->SetBranchAddress("Nleft_p0",&fNp0Aside);
    Ntu->SetBranchAddress("Etaleft_p0",&fEtap0Aside);
    Ntu->SetBranchAddress("Ipdg_l_p0",&fp0PDGCode);
    Ntu->SetBranchAddress("Pxl_p0",  fPxp0Aside);
    Ntu->SetBranchAddress("Pyl_p0",  fPyp0Aside);
    Ntu->SetBranchAddress("Pzl_p0",  fPzp0Aside);
    Ntu->SetBranchAddress("Nright_p0",&fNp0Cside);
    Ntu->SetBranchAddress("Etaright_p0",&fEtap0Cside);
    Ntu->SetBranchAddress("Pxr_p0",   fPxp0Cside);
    Ntu->SetBranchAddress("Pyr_p0",   fPyp0Cside);
    Ntu->SetBranchAddress("Pzr_p0",   fPzp0Cside);
    // **** eta
    Ntu->SetBranchAddress("Nleft_et",&fNetaAside);
    Ntu->SetBranchAddress("Etaleft_et",&fEtaetaAside);
    Ntu->SetBranchAddress("Ipdg_l_et",&fetaPDGCode);
    Ntu->SetBranchAddress("Pxl_et",  fPxetaAside);
    Ntu->SetBranchAddress("Pyl_et",  fPyetaAside);
    Ntu->SetBranchAddress("Pzl_et",  fPzetaAside);
    Ntu->SetBranchAddress("Nright_et",&fNetaCside);
    Ntu->SetBranchAddress("Etaright_et",&fEtaetaCside);
    Ntu->SetBranchAddress("Pxr_et",   fPxetaCside);
    Ntu->SetBranchAddress("Pyr_et",   fPyetaCside);
    Ntu->SetBranchAddress("Pzr_et",   fPzetaCside);
    // **** omega
    Ntu->SetBranchAddress("Nleft_om",&fNomegaAside);
    Ntu->SetBranchAddress("Etaleft_om",&fEtaomegaAside);
    Ntu->SetBranchAddress("Ipdg_l_om",&fomegaPDGCode);
    Ntu->SetBranchAddress("Pxl_om",  fPxomegaAside);
    Ntu->SetBranchAddress("Pyl_om",  fPyomegaAside);
    Ntu->SetBranchAddress("Pzl_om",  fPzomegaAside);
    Ntu->SetBranchAddress("Nright_om",&fNomegaCside);
    Ntu->SetBranchAddress("Etaright_om",&fEtaomegaCside);
    Ntu->SetBranchAddress("Pxr_om",   fPxomegaCside);
    Ntu->SetBranchAddress("Pyr_om",   fPyomegaCside);
    Ntu->SetBranchAddress("Pzr_om",   fPzomegaCside);
}

// -----------------------------------------------------------------------------------
Int_t AliGenReaderEMD::NextEvent() 
{
    // Read the next event  
    Int_t nTracks=0;
    fNparticle = 0; fOffset=0;
    
    TFile* pFile = fTreeNtuple->GetCurrentFile();
    pFile->cd();
    

    Int_t nentries = (Int_t) fTreeNtuple->GetEntries();
    if(fNcurrent < nentries) {
	fTreeNtuple->GetEvent(fNcurrent);
	if(fNcurrent%100 == 0) printf("\n *** Reading event %d ***\n",fNcurrent);
	//
	if(fPcToTrack==kAll || fPcToTrack==kOnlyNucleons){ // nucleons
	   nTracks = fNnCside+fNnAside+fNpCside+fNpAside;
	}
	if(fPcToTrack==kAll || fPcToTrack==kNotNucleons){ //pions,eta,omega
	    nTracks += fNppCside+fNpmCside+fNppAside+fNpmAside+fNp0Aside+fNp0Cside+
	    	fNetaAside+fNetaCside+fNomegaAside+fNomegaCside;
	}
	fNcurrent++;
	printf("\t #### Putting %d particles in the stack\n", nTracks);
	/*if(fPcToTrack==kAll || fPcToTrack==kOnlyNucleons) printf("\t\t  %d+%d neutrons, %d+%d protons\n", 
		fNnAside,fNnCside, fNpAside,fNpCside);
	if(fPcToTrack==kAll || fPcToTrack==kNotNucleons) printf("\t %d+%d pi+, %d+%d pi-, %d+%d pi0, %d+%d eta, %d+%d omega\n",
	        fNppAside,fNppCside,fNpmAside,fNpmCside, 
		fNp0Aside,fNp0Cside,fNetaAside,fNetaCside, fNomegaAside,fNomegaCside);*/
	return nTracks;
    }

    return 0;
}

// -----------------------------------------------------------------------------------
TParticle* AliGenReaderEMD::NextParticle() 
{
    // Read the next particle
    Float_t p[4]={0.,0.,0.,0.};
    int pdgCode=0;
    
    if(fPcToTrack==kAll || fPcToTrack==kOnlyNucleons){//***********************************************
      if(fNparticle<fNnAside){
        p[0] = fPxnAside[fNparticle];
        p[1] = fPynAside[fNparticle];
        p[2] = fPznAside[fNparticle];  
	pdgCode = fnPDGCode;
//    printf(" pc%d n sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      else if(fNparticle>=fNnAside && fNparticle<(fNnAside+fNnCside)){
        p[0] = fPxnCside[fNparticle];
        p[1] = fPynCside[fNparticle];
        p[2] = fPznCside[fNparticle];
	pdgCode = fnPDGCode;
//    printf(" pc%d n sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      else if(fNparticle>=fNnAside+fNnCside && fNparticle<(fNnAside+fNnCside+fNpAside)){
        p[0] = fPxpAside[fNparticle];
        p[1] = fPypAside[fNparticle];
        p[2] = fPzpAside[fNparticle];
	pdgCode = fpPDGCode;
//    printf(" pc%d p sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      else if(fNparticle>=fNnAside+fNnCside+fNpAside && fNparticle<(fNnAside+fNnCside+fNpCside+fNpAside)){
        p[0] = fPxpCside[fNparticle];
        p[1] = fPypCside[fNparticle];
        p[2] = fPzpCside[fNparticle];
	pdgCode = fpPDGCode;
//    printf(" pc%d p sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      fOffset = fNnAside+fNnCside+fNpCside+fNpAside;
    } //**********************************************************************************************
    if(fPcToTrack==kAll || fPcToTrack==kNotNucleons){
      if(fNparticle>=fOffset && fNparticle<fOffset+fNppAside){ // *** pi +
        p[0] = fPxppAside[fNparticle];
        p[1] = fPyppAside[fNparticle];
        p[2] = fPzppAside[fNparticle];  
	pdgCode = fppPDGCode;
//    printf(" pc%d pi+ sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside && fNparticle<fOffset+fNppAside+fNppCside){
        p[0] = fPxppCside[fNparticle];
        p[1] = fPyppCside[fNparticle];
        p[2] = fPzppCside[fNparticle];
	pdgCode = fppPDGCode;
//    printf(" pc%d pi+ sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside && fNparticle<fOffset+fNppAside+fNppCside+fNpmAside){ // *** pi -
        p[0] = fPxpmAside[fNparticle];
        p[1] = fPypmAside[fNparticle];
        p[2] = fPzpmAside[fNparticle];
	pdgCode = fpmPDGCode;
//    printf(" pc%d pi- sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside && fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside){
        p[0] = fPxpmCside[fNparticle];
        p[1] = fPypmCside[fNparticle];
        p[2] = fPzpmCside[fNparticle];
	pdgCode = fpmPDGCode;
//    printf(" pc%d pi- sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside && 
         fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside){ // *** pi 0
        p[0] = fPxp0Aside[fNparticle];
        p[1] = fPyp0Aside[fNparticle];
        p[2] = fPzp0Aside[fNparticle];
	pdgCode = fp0PDGCode;
//    printf(" pc%d pi0 sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside && 
        fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside){
        p[0] = fPxp0Cside[fNparticle];
        p[1] = fPyp0Cside[fNparticle];
        p[2] = fPzp0Cside[fNparticle];
	pdgCode = fp0PDGCode;
//    printf(" pc%d pi0 sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside && 
         fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside){ // *** eta
        p[0] = fPxetaAside[fNparticle];
        p[1] = fPyetaAside[fNparticle];
        p[2] = fPzetaAside[fNparticle];
	pdgCode = fetaPDGCode;
//    printf(" pc%d eta sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside && 
         fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside+fNetaCside){
        p[0] = fPxetaCside[fNparticle];
        p[1] = fPyetaCside[fNparticle];
        p[2] = fPzetaCside[fNparticle];
	pdgCode = fetaPDGCode;
//    printf(" pc%d eta sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside+fNetaCside && 
         fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside+fNetaCside+fNomegaAside){ // *** omega
        p[0] = fPxomegaAside[fNparticle];
        p[1] = fPyomegaAside[fNparticle];
        p[2] = fPzomegaAside[fNparticle];
	pdgCode = fomegaPDGCode;
//    printf(" pc%d omega sideA: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      if(fNparticle>=fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside+fNetaCside+fNomegaAside
      && fNparticle<fOffset+fNppAside+fNppCside+fNpmAside+fNpmCside+fNp0Aside+fNp0Cside+fNetaAside+fNetaCside+fNomegaAside+fNomegaCside){
        p[0] = fPxomegaCside[fNparticle];
        p[1] = fPyomegaCside[fNparticle];
        p[2] = fPzomegaCside[fNparticle];
	pdgCode = fomegaPDGCode;
//    printf(" pc%d omega sideC: PDG code %d,  momentum (%f, %f, %f) \n", fNparticle, pdgCode, p[0],p[1],p[2]);
      }
      
    } 
   
    Float_t ptot = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    Double_t amass = TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
    p[3] = TMath::Sqrt(ptot*ptot+amass*amass);
    
    if(p[3]<=amass){ 
       Warning("Generate","Particle %d  E = %f GeV mass = %f GeV ",pdgCode,p[3],amass);
    }
    
    //printf("  Pc %d:  PDGcode %d  p(%1.2f, %1.2f, %1.2f, %1.3f)\n",
    //	fNparticle,pdgCode,p[0], p[1], p[2], p[3]);
    
    TParticle* particle = new TParticle(pdgCode, 0, -1, -1, -1, -1, 
    	p[0], p[1], p[2], p[3], 0., 0., 0., 0.);
    if((p[0]*p[0]+p[1]*p[1]+p[2]*p[2])>1e-5) particle->SetBit(kTransportBit);
    fNparticle++;
    return particle;
}
