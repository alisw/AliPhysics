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

//_________________________________________________________________________
// TTask class for TOF Probabilities (ALICE-Offline week June 2002).
// Use case: start root and execute the following macro
// ntuple.root is assumed to be the ntuple with the reconstructed
// data
/*
{
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  // create an instance of the AliTOFProb class
  // You have to pass the ntuple
 AliTOFProb* tofprob=new AliTOFProb("ntuple.root");
 tofprob->Exec("");
}
*/
// 
//-- Author: F. Pierella | pierella@bo.infn.it
//////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH3.h"
//#include "TH3.h"
#include "TCanvas.h"
#include "TFile.h"
#include <TF1.h>
#include <TF2.h>
#include "TROOT.h"
#include "TFolder.h"
#include "TNtuple.h"
#include "TLeaf.h"

#include "AliConst.h"
#include "AliTOFConstants.h"
#include "AliTOFProb.h"

#include <stdlib.h>
#include <Riostream.h>
#include <Riostream.h>

ClassImp(AliTOFProb)

//____________________________________________________________________________ 
  AliTOFProb::AliTOFProb():TTask("AliTOFProb","") 
{
  // default ctor - set the pointer member vars to zero

  fhfile    = 0;
  fNtuple   = 0;
  fgen      = 0;
  foutfileName  = 0;
}
           
//____________________________________________________________________________ 
  AliTOFProb::AliTOFProb(char* headerFile):TTask("AliTOFProb","") 
{

  fhfile = TFile::Open(headerFile); // connect file with ntuple
  foutfileName=headerFile;
  fNtuple   = 0;
  fgen      = 0;

  Init();
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}
//____________________________________________________________________________ 
void AliTOFProb::Init()
{

  fNtuple= (TNtuple*)fhfile->Get("Ntuple"); // get ntuple from file
  Int_t nvar = fNtuple->GetNvar(); cout <<"N of var.="<< nvar << endl;
  fNtuple->GetEvent(0);

}

//____________________________________________________________________________ 
  AliTOFProb::~AliTOFProb()
{
  //
  // dtor (free used memory)
  //

  if (fhfile)
    {
      delete fhfile;
      fhfile = 0;
    }


  if (fNtuple)
    {
      delete fNtuple;
      fNtuple = 0;
    }

  if (fgen)
    {
      delete fgen;
      fgen = 0;
    }

  if (foutfileName)
    {
      delete foutfileName;
      foutfileName = 0;
    }
}


//____________________________________________________________________________
void AliTOFProb::Exec(const Option_t * /*dummyOpt*/) 
{ 
  //
  // Performs Prob for TOF detector
  // 
  fTask=1;
  fhfile->cd();
  // create the histo
  TH3* h3D=new TH3F("h3D", "", 120, -0.1, 1.1, 120, -0.1, 1.1, 120, -0.1, 1.1);
  TH3* h3Ddt=new TH3F("h3Ddt", "", 360, -120., 120., 360, -120., 120., 360, -120., 120.);

  TH1* hpion  =new TH1F("hpion", "", 120, -0.1, 1.1);
  TH1* hkaon  =new TH1F("hkaon", "", 120, -0.1, 1.1);
  TH1* hproton=new TH1F("hproton", "", 120, -0.1, 1.1);

  TH1* hdtpion  =new TH1F("hdtpion", "", 600, -120., 120.);
  TH1* hdtkaon  =new TH1F("hdtkaon", "", 600, -120., 120.);
  TH1* hdtproton=new TH1F("hdtproton", "", 600, -120., 120.);

  const Float_t timeResolution=0.12;

  Int_t nparticles = (Int_t)fNtuple->GetEntries();
  cout << " Number of nparticles =" << nparticles << endl;
  if (nparticles <= 0) return;
  
  for (Int_t i=0; i < nparticles; i++) {
    fNtuple->GetEvent(i);
    //Int_t event=(Int_t)(fNtuple->GetLeaf("event")->GetValue());
    Int_t pdgcode=(Int_t)(fNtuple->GetLeaf("ipart")->GetValue());
    // unused for the time being in the context of probabilities
    // probably it will be included in AliESD approach
    //Float_t mass=fNtuple->GetLeaf("mext")->GetValue(0);
    Int_t matc=(Int_t)(fNtuple->GetLeaf("matc")->GetValue(0));
    Int_t imam=(Int_t)(fNtuple->GetLeaf("imam")->GetValue(0));
    Float_t time=fNtuple->GetLeaf("text")->GetValue(0); // [ns]
    // generate a gaussian around time with sigma timeresolution [ns]
    Float_t lowTime=time-1000.*timeResolution;
    Float_t upTime =time+1000.*timeResolution;
    TF1* gauss=new TF1("tg","gaus",lowTime,upTime);

    Double_t p1=time;
    Double_t p2=timeResolution;
    //Float_t twopi=2*TMath::Pi();
    //Float_t sqrtwopi=TMath::Sqrt(twopi);
    //Double_t p0=1./(p2*sqrtwopi);
    // choising a different normalization
    Double_t p0=1.;
    gauss->SetParameters(p0,p1,p2);

    Float_t px=fNtuple->GetLeaf("pxvtx")->GetValue(0);
    Float_t py=fNtuple->GetLeaf("pyvtx")->GetValue(0);
    Float_t pz=fNtuple->GetLeaf("pzvtx")->GetValue(0);
    // track length [cm]
    Float_t toflen=fNtuple->GetLeaf("leng")->GetValue(0);

    Float_t pvtx=TMath::Sqrt(px*px+py*py+pz*pz);
    Float_t ptvtx=TMath::Sqrt(px*px+py*py);

    Bool_t isSelected=(imam == 0 && pz !=0 && TMath::ATan(TMath::Abs(ptvtx/pz))>TMath::Pi()*45./180.);
    Bool_t trackWithTime=(matc==2 || matc==3 || matc==4);

    Int_t abspdgcode=TMath::Abs(pdgcode);
    switch(abspdgcode){
    case 321:
      break;
    case 2212:
      break;
    case 11:
      break;
    default:
      break;
    }

    if (isSelected && trackWithTime && toflen>370.)
      {//only primary +/-45  
	Float_t beta[3]={0.,0.,0.};
	Float_t timeth[3]={0.,0.,0.};
	Float_t massarray[3]={0.13957,0.493677,0.9382723};
	for (Int_t j=0; j < 3; j++) {
	  Float_t dummy=(pvtx/massarray[j])*(pvtx/massarray[j]);
	  beta[j]=TMath::Sqrt(dummy/(1.+dummy));
	  timeth[j]=toflen/(29.9792458*beta[j]);
	}
	Float_t pPion   =(Float_t)gauss->Eval(timeth[0],0.,0.);
	Float_t pKaon   =(Float_t)gauss->Eval(timeth[1],0.,0.);
	Float_t pProton =(Float_t)gauss->Eval(timeth[2],0.,0.);

	Float_t dtpion=time-timeth[0];
	Float_t dtkaon=time-timeth[1];
	Float_t dtproton=time-timeth[2];
	hdtpion->Fill(dtpion);
	hdtkaon->Fill(dtkaon);
	hdtproton->Fill(dtproton);

	// normalizing probabilities
	Float_t sumProb=pPion+pKaon+pProton;
	pPion=pPion/sumProb;
	pKaon=pKaon/sumProb;
	pProton=pProton/sumProb;
	hpion->Fill(pPion);
	hkaon->Fill(pKaon);
	hproton->Fill(pProton);
	//cout << "pPion " << pPion << " pKaon " << pKaon << " pProton " << pProton << endl;
	h3D->Fill(pPion,pKaon,pProton);
	h3Ddt->Fill(dtpion,dtkaon,dtproton);
	//cout << "pion " << timeth[0] << " kaon " << timeth[1] << " proton " << timeth[2] << endl;
      }// End of cuts appling

    delete gauss;

  }// End of loop over particles

  char outFileName[100];
  strcpy(outFileName,"3DProbHisto.root");
  TFile *houtfile = new TFile(outFileName,"recreate");
  houtfile->cd();
  h3D->Write(0,TObject::kOverwrite);
  h3Ddt->Write(0,TObject::kOverwrite);
  hpion->Write(0,TObject::kOverwrite);
  hkaon->Write(0,TObject::kOverwrite);
  hproton->Write(0,TObject::kOverwrite);
  hdtpion->Write(0,TObject::kOverwrite);
  hdtkaon->Write(0,TObject::kOverwrite);
  hdtproton->Write(0,TObject::kOverwrite);

  houtfile->Close();
  houtfile->Write(0,TObject::kOverwrite);
}


//__________________________________________________________________
Bool_t AliTOFProb::operator==( AliTOFProb const & /*tofrec*/)const
{
  // dummy version of Equal operator.
  // requested by coding conventions
  return kTRUE;

}
