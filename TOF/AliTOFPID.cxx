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
// TTask class for TOF PID.
// Use case: start root and execute the following macro
/*
{
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  // create an instance of the AliTOFPID class
  // You have to pass the ntuple and cuts filenames
 AliTOFPID* tofpid=new AliTOFPID("ntuple.root","cuts.root");

 // option "pp" for pp events (it includes also electron in the analysis)
 // option "visual" to shows interactively histos
 // option "asC" or "asEPS" to save canvas in the current dir in .C or .eps format

 // make a choice: uncomment one of these lines
 // tofpid->Exec("pp","visual","asC");
 // tofpid->Exec("pp","novisual","asC");
 // tofpid->Exec("Pb-Pb","visual","asC");
 // tofpid->Exec("pp","visual","asC");
 // tofpid->Exec("pp","novisual","asEPS");
}
*/
//
//
//
//
//
//
// 
//
//-- Authors: B. Zagreev , F. Pierella
//////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TText.h"
#include "TLine.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFile.h"
#include <TF1.h>
#include <TF2.h>
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TNtuple.h"
#include "TLeaf.h"

#include "AliConst.h"
#include "AliTOFConstants.h"
#include "AliTOFPID.h"

#include <stdlib.h>
#include <Riostream.h>
#include <Riostream.h>

ClassImp(AliTOFPID)

//____________________________________________________________________________ 
  AliTOFPID::AliTOFPID():TTask("AliTOFPID","") 
{
  // default ctor - set the pointer member vars to zero
  felectron = 0;
  fpion     = 0;
  fkaon     = 0;
  fproton   = 0;
  fcut      = 0;
  fhfile    = 0;
  fNtuple   = 0;
  fgen      = 0;
  foutfileName  = 0;
}
           
//____________________________________________________________________________ 
  AliTOFPID::AliTOFPID(char* headerFile, char *cutsFile, const Option_t* opt):TTask("AliTOFPID","") 
{
  felectron = 0;
  fpion     = 0;
  fkaon     = 0;
  fproton   = 0;
  fhfile = TFile::Open(headerFile); // connect file with ntuple
  fcut = TFile::Open(cutsFile); // connect file for cuts
  foutfileName=headerFile;
  fNtuple   = 0;
  fgen      = 0;

  Init(opt);
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}
//____________________________________________________________________________ 
void AliTOFPID::Init(const Option_t* opt)
{
  if(strstr(opt,"pp")){ 
    if(fcut->GetKey("electron")) felectron = (TCutG*)fcut->Get("electron");
    fcut->Print();
    if(fcut->GetKey("pion")) fpion = (TCutG*)fcut->Get("pion");
    fcut->Print();
  }
  if(fcut->GetKey("kaon")) fkaon = (TCutG*)fcut->Get("kaon");
  fcut->Print();
  if(fcut->GetKey("proton")) fproton = (TCutG*)fcut->Get("proton");

  gFile->ls();
  fNtuple= (TNtuple*)fhfile->Get("Ntuple"); // get ntuple from file
  Int_t nvar = fNtuple->GetNvar(); cout <<"N of var.="<< nvar << endl;
  fNtuple->GetEvent(0);

}

//____________________________________________________________________________ 
  AliTOFPID::~AliTOFPID()
{
  //
  // dtor (free used memory)
  //

  if (felectron)
    {
      delete felectron;
      felectron = 0;
    }

  if (fpion)
    {
      delete fpion;
      fpion = 0;
    }

  if (fkaon)
    {
      delete fkaon;
      fkaon = 0;
    }

  if (fproton)
    {
      delete fproton;
      fproton = 0;
    }

  if (fcut)
    {
      delete fcut;
      fcut = 0;
    }


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
void AliTOFPID::Exec(const Option_t *eventType, const Option_t *outputmode, const Option_t *outputsavemode) 
{ 
  //
  // Performs PID for TOF detector
  // 

  fTask=1;
  TAxis *xaxis;
  ////////// Create histograms /////////////////
  // for electron only in pp case
  TH1F* eleff=0;
  TH1F* eleffls=0;
  TH1F *elcon=0;
  TH1F *elid=0;
  TH1F *elmatch=0;
  TH1F *elall=0;

  if(strstr(eventType,"pp")){
    eleff = new TH1F("eleff","",10,0,0.6);
    xaxis=eleff->GetYaxis();
    xaxis->SetLabelSize(.08);
    eleffls = new TH1F("eleffls","",10,0,0.6);
    xaxis=eleffls->GetYaxis();
    xaxis->SetLabelSize(.08);
    elcon = new TH1F("elcon","",10,0,0.6);
    xaxis=elcon->GetXaxis();
    xaxis->SetLabelSize(.09);
    xaxis=elcon->GetYaxis();
    xaxis->SetLabelSize(.08);
    elid = new TH1F("elid","Identified electrons",10,0,0.6);
    elmatch = new TH1F("elmatch","N(e)",10,0,0.6);
    elall = new TH1F("elall","Electrons",10,0,0.6);
  }

  // pions
  TH1F *pit  = new TH1F("pit","",15,0,2.5); //part. with tracks
  TH1F *pig  = new TH1F("pig","",15,0,2.5); //part. in geometry acceptance
  TH1F *pieff = new TH1F("pieff","",15,0,2.5); //efficiency
  TH1F *pieffls = new TH1F("pieffls","",15,0,2.5); //efficiency (last step)
  xaxis=pieff->GetYaxis();
  xaxis->SetLabelSize(.08);
  xaxis=pieffls->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *picon = new TH1F("picon","",15,0,2.5); //contamination
  xaxis=picon->GetXaxis();
  xaxis->SetLabelSize(.09);
  xaxis=picon->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *piid  = new TH1F("piid","Identified pions",15,0,2.5);
  TH1F *piall = new TH1F("piall","Pions",15,0,2.5);
  TH1F *pimatch = new TH1F("pimatch","N(Pions)",15,0,2.5);
  TH1F *pigen = new TH1F("pigen","Pions",15,0,2.5);
  xaxis=pigen->GetXaxis();
  xaxis->SetLabelSize(.09);
  pigen->SetXTitle("P?t! (GeV/c)");
  xaxis->SetTitleSize(.09);
  xaxis=pigen->GetYaxis();
  xaxis->SetLabelSize(.08);
  //pigen->SetYTitle("1/P?t!dN/dP?t! (GeV/c)^-2!");
  xaxis->SetTitleSize(.09);

  // kaons
  TH1F *kat  = new TH1F("kat","",15,0,2.5);
  TH1F *kag  = new TH1F("kag","",15,0,2.5);
  TH1F *kaeff = new TH1F("kaeff","",15,0,2.5);
  xaxis=kaeff->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *kaeffls = new TH1F("kaeffls","",15,0,2.5);
  xaxis=kaeffls->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *kacon = new TH1F("kacon","",15,0,2.5);
  xaxis=kacon->GetXaxis();
  xaxis->SetLabelSize(.09);
  xaxis=kacon->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *kaid  = new TH1F("kaid","Identified kaons",15,0,2.5);
  TH1F *kamatch  = new TH1F("kamatch","N(K)",15,0,2.5);
  TH1F *kaall = new TH1F("kaall","Kaons",15,0,2.5);
  TH1F *kagen = new TH1F("kagen","Kaons",15,0,2.5);
  xaxis=kagen->GetXaxis();
  xaxis->SetLabelSize(.09);
  kagen->SetXTitle("P?t! (GeV/c)");
  xaxis->SetTitleSize(.09);
  xaxis=kagen->GetYaxis();
  xaxis->SetLabelSize(.08);
  //kagen->SetYTitle("1/P?t!dN/dP?t! (GeV/c)^-2!");
  xaxis->SetTitleSize(.09);

  // protons
  TH1F *prt  = new TH1F("prt","",15,0,4.4);
  TH1F *prg  = new TH1F("prg","",15,0,4.4);
  TH1F *preff = new TH1F("preff","",15,0,4.4);
  xaxis=preff->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *preffls = new TH1F("preffls","",15,0,4.4);
  xaxis=preffls->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *prcon = new TH1F("prcon","",15,0,4.4);
  xaxis=prcon->GetXaxis();
  xaxis->SetLabelSize(.09);
  xaxis=prcon->GetYaxis();
  xaxis->SetLabelSize(.08);
  TH1F *prid  = new TH1F("prid","Identified protons",15,0,4.4);
  TH1F *prmatch  = new TH1F("prmatch","N(p)",15,0,4.4);
  TH1F *prall = new TH1F("prall","Protons",15,0,4.4);
  TH1F *prgen = new TH1F("prgen","Protons",15,0,4.4);
  xaxis=prgen->GetXaxis();
  xaxis->SetLabelSize(.09);
  prgen->SetXTitle("P?t! (GeV/c)");
  xaxis->SetTitleSize(.09);
  xaxis=prgen->GetYaxis();
  xaxis->SetLabelSize(.08);
  //prgen->SetYTitle("1/P?t!dN/dP?t! (GeV/c)^-2!");
  xaxis->SetTitleSize(.09);

  // 2-D histos (extrapolated mass vs momentum)
  TH2F* hel=0;
  if(strstr(eventType,"pp")){
    hel = new TH2F("hel","",1000,-.2,1.2,1000,-4.2,0.);
    hel->SetXTitle("Mass (GeV/c^{2})");
    hel->SetYTitle("Momentum (GeV/c)");
  }
  TH2F *hpi = new TH2F("hpi","",1000,-.2,1.2,1000,-4.2,0.);
  hpi->SetXTitle("Mass (GeV/c^{2})");
  hpi->SetYTitle("Momentum (GeV/c)");
  TH2F *hka = new TH2F("hka","",1000,-.2,1.2,1000,-4.2,0.);
  hka->SetXTitle("Mass (GeV/c^{2})");
  hka->SetYTitle("Momentum (GeV/c)");
  TH2F *hpr = new TH2F("hpr","",1000,-.2,1.2,1000,-4.2,0.);
  hpr->SetXTitle("Mass (GeV/c^{2})");
  hpr->SetYTitle("Momentum (GeV/c)");
  

  fhfile->cd();
  Int_t nparticles = (Int_t)fNtuple->GetEntries();
  cout << " Number of nparticles =" << nparticles << endl;
  if (nparticles <= 0) return;

  Float_t ka=0, pi=0, pr=0, kaal=0, pial=0, pral=0;
  Float_t pitrack=0, pimag=0, pigeom=0;
  Float_t katrack=0, kamag=0, kageom=0;
  Float_t prtrack=0, prmag=0, prgeom=0;
  Float_t pif=0, kaf=0, prf=0, pin=0, kan=0, prn=0;
  Float_t px, py, pz, x, y, z, mass;
  Int_t event, matc, imam, pdgcode;
  Int_t indexOfFile=0, numfile=0;
  //////// Loop over tracks (particles)///////////////////////
  
  for (Int_t i=0; i < nparticles; i++) {
    fNtuple->GetEvent(i);
    event=(Int_t)(fNtuple->GetLeaf("event")->GetValue());
    pdgcode=(Int_t)(fNtuple->GetLeaf("ipart")->GetValue());
    mass=fNtuple->GetLeaf("mext")->GetValue(0);
    matc=(Int_t)(fNtuple->GetLeaf("matc")->GetValue(0));
    imam=(Int_t)(fNtuple->GetLeaf("imam")->GetValue(0));
    px=fNtuple->GetLeaf("pxvtx")->GetValue(0);
    py=fNtuple->GetLeaf("pyvtx")->GetValue(0);
    pz=fNtuple->GetLeaf("pzvtx")->GetValue(0);
    x=fNtuple->GetLeaf("xvtx")->GetValue(0);
    y=fNtuple->GetLeaf("yvtx")->GetValue(0);
    z=fNtuple->GetLeaf("zvtx")->GetValue(0);
    Float_t pvtx=TMath::Sqrt(px*px+py*py+pz*pz);
    Float_t ptvtx=TMath::Sqrt(px*px+py*py);
    Float_t mt=0.;
    Bool_t isSelected=(imam == 0 && pz !=0 && TMath::ATan(TMath::Abs(ptvtx/pz))>TMath::Pi()*45./180.);
    Int_t abspdgcode=TMath::Abs(pdgcode);
    switch(abspdgcode){
    case 321:
      if(isSelected && (matc==3 || matc==4)) kamatch->Fill(pvtx);
      mt=TMath::Sqrt(AliTOFConstants::fgkKaonMass*AliTOFConstants::fgkKaonMass+px*px+py*py);
      break;
    case 2212:
      if(isSelected && (matc==2 || matc==3 || matc==4)) prmatch->Fill(pvtx);
      mt=TMath::Sqrt(AliTOFConstants::fgkProtonMass*AliTOFConstants::fgkProtonMass+px*px+py*py);
      break;
    case 11:
      if(strstr(eventType,"pp") && (matc==3 || matc==4)) elmatch->Fill(pvtx); //  as in kaon case
      mt=TMath::Sqrt(AliTOFConstants::fgkElectronMass*AliTOFConstants::fgkElectronMass+px*px+py*py);
      break;
    default:
      if(isSelected && matc>0) pimatch->Fill(pvtx);
      mt=TMath::Sqrt(AliTOFConstants::fgkPionMass*AliTOFConstants::fgkPionMass+px*px+py*py);
      break;
    }

    if (isSelected)
      {//only primary +/-45  
	if (fkaon->IsInside(mass,-pvtx) && matc>2) {
	  ka++;
	  if (fTask!=2) kaid->Fill(pvtx); else {kaid->Fill(ptvtx);}
	  if (TMath::Abs(pdgcode)==321) {kaf++; kaeff->Fill(pvtx); kaeffls->Fill(pvtx);} else {kan++; kacon->Fill(pvtx);}
	} else if (fproton->IsInside(mass,-pvtx) && matc>1) {
	  pr++;
	  if (fTask!=2) prid->Fill(pvtx); else 
	    {prid->Fill(ptvtx);}
	  if (TMath::Abs(pdgcode)==2212) {prf++; preff->Fill(pvtx); preffls->Fill(pvtx);} else {prn++; prcon->Fill(pvtx);}
	} else  if (strstr(eventType,"pp") && felectron->IsInside(mass,-pvtx) && matc>2) {elid->Fill(pvtx);
	if (strstr(eventType,"pp") && TMath::Abs(pdgcode)==11) {eleff->Fill(pvtx); eleffls->Fill(pvtx);} else {elcon->Fill(pvtx);}
	} else if (matc>0) {
	  //||matc==-4&&fpion->IsInside(mass,-pvtx)
	  pi++;
	  if (fTask!=2) piid->Fill(pvtx); else {piid->Fill(ptvtx);}
	  if (TMath::Abs(pdgcode)==211) {pif++; pieff->Fill(pvtx); pieffls->Fill(pvtx);} else {pin++; picon->Fill(pvtx);}
	}

	//////////////// Normalization histograms ////////////////////
	if (strstr(eventType,"pp") && TMath::Abs(pdgcode)==11) {
	  if (fTask!=2) elall->Fill(pvtx); else elall->Fill(ptvtx);
	  if (fTask==1) hel->Fill(mass,-pvtx);
	  
	} else if (TMath::Abs(pdgcode)==211) {
	  pial++;
	  if (matc!=0) {
	    pitrack++; 
	    pit->Fill(pvtx);
	    if (matc!=-1) {
	      pimag++;
	      if (matc>-2 || matc==-4) {
		pigeom++;
		pig->Fill(pvtx);
	      }
	    }
	  }
	  if (fTask!=2) piall->Fill(pvtx); 
	  else {
	    piall->Fill(ptvtx,1/ptvtx);
	  }
	  if (fTask==1) hpi->Fill(mass,-pvtx); 
	  
	} else if (TMath::Abs(pdgcode)==321) {
	  kaal++;
	  if (matc!=0) {
	    katrack++;
	    kat->Fill(pvtx);
	    if (matc!=-1) {
	      kamag++;
	      if (matc>-2 || matc==-4) {
		kageom++;
		kag->Fill(pvtx);
	      }
	    }
	  }
	  if (fTask!=2) kaall->Fill(pvtx); 
	  else {
	    kaall->Fill(ptvtx,1/ptvtx);
	  }
	  if (fTask==1) hka->Fill(mass,-pvtx);
	  
	} else if (TMath::Abs(pdgcode)==2212) {
	  pral++;
	  if (matc!=0) {
	    prtrack++;
	    prt->Fill(pvtx);
	    if (matc!=-1) {
	      prmag++;
	      if (matc>-2 || matc==-4) {
		prgeom++;
		prg->Fill(pvtx);
	      }
	    }
	  }
	  if (fTask!=2) prall->Fill(pvtx); 
	  else {
	    prall->Fill(ptvtx,1/ptvtx);
	  }
	  if (fTask==1) hpr->Fill(mass,-pvtx);}
	
      }// End of cuts appling
  }// End of loop over particles

  // display results
  cout<< "Pions in 45-135 deg.     "<< pial <<" (100%)"<< endl;
  cout<< "Pions that have track    "<< pitrack/pial*100 <<" %"<<endl;
  cout<< "Magnetic field           "<< pimag/pial*100 <<" %"<<endl;
  cout<< "Geometry efficiency      "<< pigeom/pial*100 <<" %"<<endl;
  cout<< "PID procedure            "<< pif/pial*100 <<" %"<<endl;
  cout<< "Contamination            "<< pin/pi*100 <<" %"<<endl;
  cout<<endl;
  cout<< "Kaons in 45-135 deg.     "<< kaal <<" (100%)"<< endl;
  cout<< "Kaons that have track    "<< katrack/kaal*100 <<" %"<<endl;
  cout<< "Magnetic field           "<< kamag/kaal*100 <<" %"<<endl;
  cout<< "Geometry efficiency      "<< kageom/kaal*100 <<" %"<<endl;
  cout<< "PID procedure(+decays)   "<< kaf/kaal*100 <<" %"<<endl;
  cout<< "Contamination            "<< kan/ka*100 <<" %"<<endl;
  cout<<endl;
  cout<< "Protons in 45-135 deg.   "<< pral <<" (100%)"<< endl;
  cout<< "Protons that have track  "<< prtrack/pral*100 <<" %"<<endl;
  cout<< "Magnetic field           "<< prmag/pral*100 <<" %"<<endl;
  cout<< "Geometry efficiency      "<< prgeom/pral*100 <<" %"<<endl;
  cout<< "PID procedure            "<< prf/pral*100 <<" %"<<endl;
  cout<< "Contamination            "<< prn/pr*100 <<" %"<<endl;
  cout<<endl;
  cout<< "All part. in 45-135 deg.  "<< pial+kaal+pral <<" (100%)"<< endl;
  cout<< "All part. that have track "<< (pitrack+katrack+prtrack)/(pial+kaal+pral)*100 <<" %"<<endl;
  cout<< "Magnetic field            "<< (pimag+kamag+prmag)/(pial+kaal+pral)*100 <<" %"<<endl;
  cout<< "Geometry efficiency       "<< (pigeom+kageom+prgeom)/(pial+kaal+pral)*100 <<" %"<<endl;
  cout<< "PID procedure             "<< (pif+kaf+prf)/(pial+kaal+pral)*100 <<" %"<<endl;
  cout<< "Contamination             "<< (pin+kan+prn)/(pi+ka+pr)*100 <<" %"<<endl;
  cout<<endl;
  
  TCanvas *pidCanvas=0;   // overall
  TCanvas *pidCanvasls=0; // last step of PID
  TCanvas *momvsmassCanvas=0;
  // overall Efficiency
  TPad  *tp=0;
  TPad  *pad1=0;
  TPad  *pad2=0;
  TPad  *pad3=0;
  TPad  *pad4=0;
  TPad  *pad5=0;
  TPad  *pad6=0;
  TPad  *pad7=0;
  TPad  *pad8=0;

  // last step Efficiency
  TPad  *tpls=0;
  TPad  *pad1ls=0;
  TPad  *pad2ls=0;
  TPad  *pad3ls=0;
  TPad  *pad4ls=0;
  TPad  *pad5ls=0;
  TPad  *pad6ls=0;
  TPad  *pad7ls=0;
  TPad  *pad8ls=0;

  //////////////////////// For fTask 1 ///////////////////////////
  if (fTask==1) {
    if (strstr(eventType,"pp")){
      eleff->Divide(elall);
      eleffls->Divide(elmatch);
    }
    // overall efficiency
    pieff->Divide(piall);
    kaeff->Divide(kaall);
    preff->Divide(prall);

    // last step efficiency
    pieffls->Divide(pimatch);
    kaeffls->Divide(kamatch);
    preffls->Divide(prmatch);

    pit->Divide(piall);
    kat->Divide(kaall);
    prt->Divide(prall);
    pig->Divide(piall);
    kag->Divide(kaall);
    prg->Divide(prall);
    if (strstr(eventType,"pp")){
      elcon->Divide(elid);
    }
    // contamination
    picon->Divide(piid);
    kacon->Divide(kaid);
    prcon->Divide(prid);

    //Create a canvas, set the view range, show histograms
    // for overall PID
    if (indexOfFile==0) {
      // overall Efficiency canvas
      pidCanvas = new TCanvas("pidCanvas","PID (Overall)",10,100,800,500);
      pidCanvas->SetBorderMode(0);
      pidCanvas->SetBorderSize(0);
      pidCanvas->SetFillColor(0);
      pidCanvas->SetFillStyle(0);

      // last step Efficiency canvas
      pidCanvasls = new TCanvas("pidCanvasls","PID (Last Step)",10,100,800,500);
      pidCanvasls->SetBorderMode(0);
      pidCanvasls->SetBorderSize(0);
      pidCanvasls->SetFillColor(0);
      pidCanvasls->SetFillStyle(0);

      if (strstr(outputmode,"visual")) {
	pidCanvas->Draw();
	pidCanvasls->Draw();
      }

      Float_t pxs=0.25+0.125; //X size of pad
      Float_t pys=0.5+0.055;  //y size of pad
      // overall
      tp   = new TPad("histo","Histograms",.1,.1,.9,.9);
      // last step
      tpls = new TPad("histo","Histograms",.1,.1,.9,.9);

      if (strstr(eventType,"Pb-Pb")){
	// overall efficiency
	//pad1 = new TPad("pad1","electron efficiency",0.,.5-.055,0.+pxs,.5-.055+pys-.00001,0,0,0);
	pad2 = new TPad("pad2","pion efficiency",0.,0.5-.055,0.+pxs,0.5-.055+pys-.00001,0,0,0);
	pad3 = new TPad("pad3","kaon efficiency",0.3,0.5-.055,0.3+pxs,0.5-.055+pys-.00001,0,0,0);
	pad4 = new TPad("pad4","proton efficiency",0.6,0.5-.055,0.6+pxs,0.5-.055+pys-.00001,0,0,0);

	// contamination
	//pad5 = new TPad("pad5","electron contamination",0.,0.,0.+pxs,0.+pys,0,0,0);
	pad6 = new TPad("pad6","pion contamination",0.,0.,0.+pxs,0.+pys,0,0,0);
	pad7 = new TPad("pad7","kaon contamination",.3,0.,0.3+pxs,0.+pys,0,0,0);
	pad8 = new TPad("pad8","proton contamination",.6,0.,0.6+pxs,0.+pys,0,0,0);

	// we repeat the same for the last step of PID
	//pad1ls = new TPad("pad1ls","electron efficiency",0.,.5-.055,0.+pxs,.5-.055+pys-.00001,0,0,0);
	pad2ls = new TPad("pad2ls","pion efficiency",0.,0.5-.055,0.+pxs,0.5-.055+pys-.00001,0,0,0);
	pad3ls = new TPad("pad3ls","kaon efficiency",0.3,0.5-.055,0.3+pxs,0.5-.055+pys-.00001,0,0,0);
	pad4ls = new TPad("pad4ls","proton efficiency",0.6,0.5-.055,0.6+pxs,0.5-.055+pys-.00001,0,0,0);

	// contamination
	//pad5 = new TPad("pad5","electron contamination",0.,0.,0.+pxs,0.+pys,0,0,0);
	pad6ls = new TPad("pad6ls","pion contamination",0.,0.,0.+pxs,0.+pys,0,0,0);
	pad7ls = new TPad("pad7ls","kaon contamination",.3,0.,0.3+pxs,0.+pys,0,0,0);
	pad8ls = new TPad("pad8ls","proton contamination",.6,0.,0.6+pxs,0.+pys,0,0,0);

      }

      if (strstr(eventType,"pp")){
	// overall Efficiency
	pad1 = new TPad("pad1","electron efficiency",0.,.5-.055,0.25+0.045,1.,0,0,0);
	pad2 = new TPad("pad2","pion efficiency",0.25-0.015,0.5-.055,0.5+0.03,1.,0,0,0);
	pad3 = new TPad("pad3","kaon efficiency",0.5-0.03,0.5-.055,0.75+0.015,1.,0,0,0);
	pad4 = new TPad("pad4","proton efficiency",0.75-0.045,0.5-.055,1.,1.,0,0,0);

	// contamination
	pad5 = new TPad("pad5","electron contamination",0.,0.,.25+.045,.5+.055,0,0,0);
	pad6 = new TPad("pad6","pion contamination",.25-.015,0.,.5+.03,.5+.055,0,0,0);
	pad7 = new TPad("pad7","kaon contamination",.5-.03,0.,.75+.015,.5+.055,0,0,0);
	pad8 = new TPad("pad8","proton contamination",.75-.045,0.,1.,.5+.055,0,0,0);


	// we repeat the same for the last step of PID
	pad1ls = new TPad("pad1ls","electron efficiency",0.,.5-.055,0.25+0.045,1.,0,0,0);
	pad2ls = new TPad("pad2ls","pion efficiency",0.25-0.015,0.5-.055,0.5+0.03,1.,0,0,0);
	pad3ls = new TPad("pad3ls","kaon efficiency",0.5-0.03,0.5-.055,0.75+0.015,1.,0,0,0);
	pad4ls = new TPad("pad4ls","proton efficiency",0.75-0.045,0.5-.055,1.,1.,0,0,0);

	// contamination
	pad5ls = new TPad("pad5ls","electron contamination",0.,0.,.25+.045,.5+.055,0,0,0);
	pad6ls = new TPad("pad6ls","pion contamination",.25-.015,0.,.5+.03,.5+.055,0,0,0);
	pad7ls = new TPad("pad7ls","kaon contamination",.5-.03,0.,.75+.015,.5+.055,0,0,0);
	pad8ls = new TPad("pad8ls","proton contamination",.75-.045,0.,1.,.5+.055,0,0,0);


      }

      // last step of PID
      gStyle->SetOptStat(0);
      tpls->SetFillStyle(0);
      tpls->SetFillColor(0);
      tpls->SetBorderSize(0);
      pidCanvasls->cd();
      TText *text1ls= new TText(.1,.2,"Contamination              Efficiency"); 
      text1ls->SetTextAngle(90);
      if (strstr(outputmode,"visual")) text1ls->Draw();
      //tp->DrawText(.3,.0,"p (GeV/c");
      pidCanvasls->cd();
      TText *text2ls= new TText(.8,.0,"p (GeV/c)"); 
      if (strstr(outputmode,"visual")) {
	text2ls->Draw();
	tpls->Draw();
      }

      // overall
      gStyle->SetOptStat(0);
      tp->SetFillStyle(0);
      tp->SetFillColor(0);
      tp->SetBorderSize(0);
      pidCanvas->cd();
      TText *text1= new TText(.1,.2,"Contamination              Efficiency"); 
      text1->SetTextAngle(90);
      if (strstr(outputmode,"visual")) text1->Draw();
      //tp->DrawText(.3,.0,"p (GeV/c");
      pidCanvas->cd();
      TText *text2= new TText(.8,.0,"p (GeV/c)"); 
      if (strstr(outputmode,"visual")) {
	text2->Draw();
	tp->Draw();
      }

    }


    // drawing histos for overall case
    if (strstr(eventType,"pp")){
      pad1->SetFillStyle(0);
      pad1->SetFillColor(10);
      tp->cd();
      if (strstr(outputmode,"visual")) pad1->Draw();
      pad1->cd();
      //eleff->SetLineWidth(15);
      eleff->SetLineWidth(3);
      eleff->SetMaximum(1.);
      if (indexOfFile==0) {
	//eleff->SetFillColor(33);
	//eleff->SetFillColor(30);
	//eleff->SetFillColor(0);
	//eleff->SetLineColor(4);
	eleff->SetLineColor(1);
	if (strstr(outputmode,"visual")) eleff->Draw();
      } else {
	eleff->SetFillStyle(0);
	eleff->SetFillColor(30);
	if (indexOfFile==2) {
	  eleff->SetLineColor(3);
	  eleff->SetLineWidth(3);
	} else {
	  //eleff->SetLineColor(2);
	  eleff->SetLineColor(1);
	  eleff->SetLineStyle(2);}
	if (strstr(outputmode,"visual")) eleff->Draw("same");
      }
      //   eleff->Fit("pol1");
      TPaveLabel *ellab = new TPaveLabel(.42,.85,.52,1.05,"e");
      if (strstr(outputmode,"visual")) ellab->Draw();
    }
    
    pad2->SetFillStyle(0);
    pad2->SetFillColor(10);
    tp->cd();
    if (strstr(outputmode,"visual")) pad2->Draw();
    pad2->cd();
    pieff->SetLineWidth(3);
    pieff->SetMaximum(1.);
    if (indexOfFile==0) {
      pieff->SetLineColor(1);
      if (strstr(outputmode,"visual")) pieff->Draw();
    } else {
      pieff->SetFillStyle(0);
      pieff->SetFillColor(30);
      if (indexOfFile==1) {
	pieff->SetLineStyle(2);
      } else if (indexOfFile==2) {
	pieff->SetLineStyle(3);
      } else {
	pieff->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) pieff->Draw("same");
    }
    TPaveLabel *pilab = new TPaveLabel(1.7,.85,2.2,1.05,"#pi");
    if (strstr(outputmode,"visual")) pilab->Draw();
    
    pad3->SetFillStyle(0);
    pad3->SetFillColor(10);
    tp->cd();
    if (strstr(outputmode,"visual")) pad3->Draw();
    pad3->cd();
    kaeff->SetLineWidth(3);
    kaeff->SetMaximum(1.);
    if (indexOfFile==0) {
      kaeff->SetLineColor(1);
      if (strstr(outputmode,"visual")) kaeff->Draw();
    } else {
      kaeff->SetFillStyle(0);
      kaeff->SetFillColor(30);
      if (indexOfFile==1) { 
	kaeff->SetLineStyle(2);
      } else if (indexOfFile==2) {
	kaeff->SetLineStyle(3);
      } else {
	kaeff->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) kaeff->Draw("same");
    }
    TPaveLabel *kalab = new TPaveLabel(1.7,.85,2.2,1.05,"K");
    if (strstr(outputmode,"visual")) kalab->Draw();
    
    pad4->SetFillStyle(0);
    pad4->SetFillColor(10);
    tp->cd();
    if (strstr(outputmode,"visual")) pad4->Draw();
    pad4->cd();
    preff->SetLineWidth(3);
    preff->SetMaximum(1.);
    if (indexOfFile==0) {
      preff->SetLineColor(1);
      if (strstr(outputmode,"visual")) preff->Draw();
    } else {
      preff->SetFillStyle(0);
      preff->SetFillColor(30);
      if (indexOfFile==1) {
	preff->SetLineStyle(2);
      } else if (indexOfFile==2) {
	preff->SetLineStyle(3);
      } else {
	preff->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) preff->Draw("same");
    }
    TPaveLabel *prlab = new TPaveLabel(3.2,.85,4.1,1.05,"p");
    if (strstr(outputmode,"visual")) prlab->Draw();

    if (strstr(eventType,"pp")){
      pad5->SetFillStyle(0);
      pad5->SetFillColor(10);
      tp->cd();
      if (strstr(outputmode,"visual")) pad5->Draw();
      pad5->cd();
      //elcon->SetLineWidth(5);
      elcon->SetLineWidth(3);
      elcon->SetMaximum(1.);
      if (indexOfFile==0) {
	//elcon->SetFillColor(33);
	//elcon->SetFillColor(30);
	//elcon->SetLineColor(4);
	elcon->SetLineColor(1);
	if (strstr(outputmode,"visual")) elcon->Draw();
      } else {
	elcon->SetFillStyle(4000);
	elcon->SetFillColor(30);
	if (indexOfFile==2) {
	  elcon->SetLineColor(3);
	  elcon->SetLineWidth(3);
	} else {
	  elcon->SetLineColor(2);
	  elcon->SetLineStyle(2);}
	if (strstr(outputmode,"visual")) elcon->Draw("same");
      }
    }


    pad6->SetFillStyle(0);
    pad6->SetFillColor(10);
    tp->cd();
    if (strstr(outputmode,"visual")) pad6->Draw();
    pad6->cd();
    picon->SetLineWidth(3);
    picon->SetMaximum(1.);
    if (indexOfFile==0) {
      picon->SetLineColor(1);
      if (strstr(outputmode,"visual")) picon->Draw();
    } else {
      picon->SetFillStyle(0);
      picon->SetFillColor(30);
      if (indexOfFile==1) { 
	picon->SetLineStyle(2);
      } else if (indexOfFile==2) {
	picon->SetLineStyle(3);
      } else {
	picon->SetLineStyle(4);
	TLine* line;
	line = new TLine(0.2,0.85,0.9,0.85);
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	if (strstr(outputmode,"visual")) line->Draw();
	line = new TLine(0.2,0.65,0.9,0.65);
	line->SetLineWidth(2);
	if (strstr(outputmode,"visual")) line->Draw();
	line = new TLine(0.2,0.45,0.9,0.45);
	line->SetLineStyle(3);
	line->SetLineWidth(1);
	if (strstr(outputmode,"visual")) line->Draw();
	line = new TLine(0.2,0.25,0.9,0.25);
	line->SetLineStyle(4);
	line->SetLineWidth(1);
	if (strstr(outputmode,"visual")) line->Draw();
	TPaveLabel *pl = new TPaveLabel(1.1,0.8,1.9,0.9,"100 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
	pl = new TPaveLabel(1.1,0.6,1.9,0.7,"150 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
	pl = new TPaveLabel(1.1,0.4,1.9,0.5,"200 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
	pl = new TPaveLabel(1.1,0.2,1.9,0.3,"300 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
      }
      if (strstr(outputmode,"visual")) picon->Draw("same");
    }
    
    pad7->SetFillStyle(0);
    pad7->SetFillColor(10);
    tp->cd();
    if (strstr(outputmode,"visual")) pad7->Draw();
    pad7->cd();
    kacon->SetLineWidth(3);
    kacon->SetMaximum(1.);
    if (indexOfFile==0) {
      kacon->SetLineColor(1);
      if (strstr(outputmode,"visual")) kacon->Draw();
    } else {
      kacon->SetFillStyle(0);
      kacon->SetFillColor(30);
      if (indexOfFile==1) {
	kacon->SetLineStyle(2);
      } else if (indexOfFile==2) {
	kacon->SetLineStyle(3);
      } else {
	kacon->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) kacon->Draw("same");
    }
    
    pad8->SetFillStyle(0);
    pad8->SetFillColor(10);
    tp->cd();
    if (strstr(outputmode,"visual")) pad8->Draw();
    pad8->cd();
    prcon->SetLineWidth(3);
    prcon->SetMaximum(1.);
    if (indexOfFile==0) {
      prcon->SetLineColor(1);
      if (strstr(outputmode,"visual")) prcon->Draw();
    } else {
      prcon->SetFillStyle(0);
      prcon->SetFillColor(30);
      if (indexOfFile==1) {
	prcon->SetLineStyle(2);
      } else if (indexOfFile==2) {
	prcon->SetLineStyle(3);
      } else {
	prcon->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) prcon->Draw("same");
    }



    // last step case (it is just a copy of the previous lines)
    // moving to pidCanvasls canvas
    pidCanvasls->cd();
    // drawing histos for overall case
    if (strstr(eventType,"pp")){
      pad1ls->SetFillStyle(0);
      pad1ls->SetFillColor(10);
      tpls->cd();
      if (strstr(outputmode,"visual")) pad1ls->Draw();
      pad1ls->cd();
      //eleff->SetLineWidth(15);
      eleffls->SetLineWidth(3);
      eleffls->SetMaximum(1.);
      if (indexOfFile==0) {
	//eleff->SetFillColor(33);
	//eleff->SetFillColor(30);
	//eleff->SetFillColor(0);
	//eleff->SetLineColor(4);
	eleffls->SetLineColor(1);
	if (strstr(outputmode,"visual")) eleffls->Draw();
      } else {
	eleffls->SetFillStyle(0);
	eleffls->SetFillColor(30);
	if (indexOfFile==2) {
	  eleffls->SetLineColor(3);
	  eleffls->SetLineWidth(3);
	} else {
	  //eleff->SetLineColor(2);
	  eleffls->SetLineColor(1);
	  eleffls->SetLineStyle(2);}
	if (strstr(outputmode,"visual")) eleffls->Draw("same");
      }
      //   eleff->Fit("pol1");
      TPaveLabel *ellabls = new TPaveLabel(.42,.85,.52,1.05,"e");
      if (strstr(outputmode,"visual")) ellabls->Draw();
    }
    
    pad2ls->SetFillStyle(0);
    pad2ls->SetFillColor(10);
    tpls->cd();
    if (strstr(outputmode,"visual")) pad2ls->Draw();
    pad2ls->cd();
    pieffls->SetLineWidth(3);
    pieffls->SetMaximum(1.);
    if (indexOfFile==0) {
      pieffls->SetLineColor(1);
      if (strstr(outputmode,"visual")) pieffls->Draw();
    } else {
      pieffls->SetFillStyle(0);
      pieffls->SetFillColor(30);
      if (indexOfFile==1) {
	pieffls->SetLineStyle(2);
      } else if (indexOfFile==2) {
	pieffls->SetLineStyle(3);
      } else {
	pieffls->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) pieffls->Draw("same");
    }
    TPaveLabel *pilabls = new TPaveLabel(1.7,.85,2.2,1.05,"#pi");
    if (strstr(outputmode,"visual")) pilabls->Draw();
    
    pad3ls->SetFillStyle(0);
    pad3ls->SetFillColor(10);
    tpls->cd();
    if (strstr(outputmode,"visual")) pad3ls->Draw();
    pad3ls->cd();
    kaeffls->SetLineWidth(3);
    kaeffls->SetMaximum(1.);
    if (indexOfFile==0) {
      kaeffls->SetLineColor(1);
      if (strstr(outputmode,"visual")) kaeffls->Draw();
    } else {
      kaeffls->SetFillStyle(0);
      kaeffls->SetFillColor(30);
      if (indexOfFile==1) { 
	kaeffls->SetLineStyle(2);
      } else if (indexOfFile==2) {
	kaeffls->SetLineStyle(3);
      } else {
	kaeffls->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) kaeffls->Draw("same");
    }
    TPaveLabel *kalabls = new TPaveLabel(1.7,.85,2.2,1.05,"K");
    if (strstr(outputmode,"visual")) kalabls->Draw();
    
    pad4ls->SetFillStyle(0);
    pad4ls->SetFillColor(10);
    tpls->cd();
    if (strstr(outputmode,"visual")) pad4ls->Draw();
    pad4ls->cd();
    preffls->SetLineWidth(3);
    preffls->SetMaximum(1.);
    if (indexOfFile==0) {
      preffls->SetLineColor(1);
      if (strstr(outputmode,"visual")) preffls->Draw();
    } else {
      preffls->SetFillStyle(0);
      preffls->SetFillColor(30);
      if (indexOfFile==1) {
	preffls->SetLineStyle(2);
      } else if (indexOfFile==2) {
	preffls->SetLineStyle(3);
      } else {
	preffls->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) preffls->Draw("same");
    }
    TPaveLabel *prlabls = new TPaveLabel(3.2,.85,4.1,1.05,"p");
    if (strstr(outputmode,"visual")) prlabls->Draw();

    if (strstr(eventType,"pp")){
      pad5ls->SetFillStyle(0);
      pad5ls->SetFillColor(10);
      tpls->cd();
      if (strstr(outputmode,"visual")) pad5ls->Draw();
      pad5ls->cd();
      //elcon->SetLineWidth(5);
      elcon->SetLineWidth(3);
      elcon->SetMaximum(1.);
      if (indexOfFile==0) {
	//elcon->SetFillColor(33);
	//elcon->SetFillColor(30);
	//elcon->SetLineColor(4);
	elcon->SetLineColor(1);
	if (strstr(outputmode,"visual")) elcon->Draw();
      } else {
	elcon->SetFillStyle(4000);
	elcon->SetFillColor(30);
	if (indexOfFile==2) {
	  elcon->SetLineColor(3);
	  elcon->SetLineWidth(3);
	} else {
	  elcon->SetLineColor(2);
	  elcon->SetLineStyle(2);}
	if (strstr(outputmode,"visual")) elcon->Draw("same");
      }
    }


    pad6ls->SetFillStyle(0);
    pad6ls->SetFillColor(10);
    tpls->cd();
    if (strstr(outputmode,"visual")) pad6ls->Draw();
    pad6ls->cd();
    picon->SetLineWidth(3);
    picon->SetMaximum(1.);
    if (indexOfFile==0) {
      picon->SetLineColor(1);
      if (strstr(outputmode,"visual")) picon->Draw();
    } else {
      picon->SetFillStyle(0);
      picon->SetFillColor(30);
      if (indexOfFile==1) { 
	picon->SetLineStyle(2);
      } else if (indexOfFile==2) {
	picon->SetLineStyle(3);
      } else {
	picon->SetLineStyle(4);
	TLine* line;
	line = new TLine(0.2,0.85,0.9,0.85);
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	if (strstr(outputmode,"visual")) line->Draw();
	line = new TLine(0.2,0.65,0.9,0.65);
	line->SetLineWidth(2);
	if (strstr(outputmode,"visual")) line->Draw();
	line = new TLine(0.2,0.45,0.9,0.45);
	line->SetLineStyle(3);
	line->SetLineWidth(1);
	if (strstr(outputmode,"visual")) line->Draw();
	line = new TLine(0.2,0.25,0.9,0.25);
	line->SetLineStyle(4);
	line->SetLineWidth(1);
	if (strstr(outputmode,"visual")) line->Draw();
	TPaveLabel *pl = new TPaveLabel(1.1,0.8,1.9,0.9,"100 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
	pl = new TPaveLabel(1.1,0.6,1.9,0.7,"150 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
	pl = new TPaveLabel(1.1,0.4,1.9,0.5,"200 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
	pl = new TPaveLabel(1.1,0.2,1.9,0.3,"300 ps","br");
	pl->SetFillColor(18);
	pl->SetTextSize(0.99);
	if (strstr(outputmode,"visual")) pl->Draw();
      }
      if (strstr(outputmode,"visual")) picon->Draw("same");
    }
    
    pad7ls->SetFillStyle(0);
    pad7ls->SetFillColor(10);
    tpls->cd();
    if (strstr(outputmode,"visual")) pad7ls->Draw();
    pad7ls->cd();
    kacon->SetLineWidth(3);
    kacon->SetMaximum(1.);
    if (indexOfFile==0) {
      kacon->SetLineColor(1);
      if (strstr(outputmode,"visual")) kacon->Draw();
    } else {
      kacon->SetFillStyle(0);
      kacon->SetFillColor(30);
      if (indexOfFile==1) {
	kacon->SetLineStyle(2);
      } else if (indexOfFile==2) {
	kacon->SetLineStyle(3);
      } else {
	kacon->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) kacon->Draw("same");
    }
    
    pad8ls->SetFillStyle(0);
    pad8ls->SetFillColor(10);
    tpls->cd();
    if (strstr(outputmode,"visual")) pad8ls->Draw();
    pad8ls->cd();
    prcon->SetLineWidth(3);
    prcon->SetMaximum(1.);
    if (indexOfFile==0) {
      prcon->SetLineColor(1);
      if (strstr(outputmode,"visual")) prcon->Draw();
    } else {
      prcon->SetFillStyle(0);
      prcon->SetFillColor(30);
      if (indexOfFile==1) {
	prcon->SetLineStyle(2);
      } else if (indexOfFile==2) {
	prcon->SetLineStyle(3);
      } else {
	prcon->SetLineStyle(4);}
      if (strstr(outputmode,"visual")) prcon->Draw("same");
    }







    // momentum vs mass 2-D histos

    if (indexOfFile==0) {
      momvsmassCanvas = new TCanvas("momvsmassCanvas","Momentum vs mass disribution",500,10,700,700);
      momvsmassCanvas->SetFillColor(0);
      momvsmassCanvas->SetBorderMode(0);
      gPad->SetFillStyle(0);
      gPad->SetBorderMode(0);
      gPad->SetFillColor(0);
      //   gStyle->SetOptStat(11);
      if (numfile==4) momvsmassCanvas->Divide(1,2,0,0); 
    } else if (indexOfFile==1 && numfile == 4) {
      momvsmassCanvas->cd(1);
      TPaveLabel *pl = new TPaveLabel(-0.0376218,-3.03586,0.0979277,-2.70158,"100 ps","br");
      pl->SetFillColor(18);
      pl->SetTextSize(0.99);
      if (strstr(outputmode,"visual")) pl->Draw();
    } else if (indexOfFile==3 && numfile == 4) {
      momvsmassCanvas->cd(2);
      TPaveLabel *pl = new TPaveLabel(-0.0591866,-3.17077,0.076363,-2.86857,"300 ps","br");
      pl->SetFillColor(18);
      pl->SetTextSize(0.99);
      if (strstr(outputmode,"visual")) pl->Draw();
    }
    if (numfile !=4) momvsmassCanvas->cd();
    if (numfile !=4 || indexOfFile==1 || indexOfFile==3) {
      //   hpi->PaintStat2(01);
      hpi->SetMarkerColor(5);
      if (strstr(outputmode,"visual")) hpi->Draw();
      if(strstr(eventType,"pp")){
	hel->SetMarkerColor(2);
	if (strstr(outputmode,"visual")) hel->Draw("same");
      }
      hka->SetMarkerColor(4);
      if (strstr(outputmode,"visual")) hka->Draw("same");
      hpr->SetMarkerColor(3);
      if (strstr(outputmode,"visual")) hpr->Draw("same");
      if (strstr(outputmode,"visual")) {
	fkaon->Draw();
	fproton->Draw();
	if(strstr(eventType,"pp")){
	  felectron->Draw();
	  fpion->Draw();
	}
      }
      if(strstr(eventType,"pp")){
	//TPaveText *ep = new TPaveText(-0.05,-0.5,0.05,-0.3);
	//ep->AddText("e");
	TPaveLabel *ep = new TPaveLabel(.42,.85,.52,1.05,"e");
	if (strstr(outputmode,"visual")) ep->Draw();
      }

      TPaveText *pip = new TPaveText(0.15,-1.0,0.25,-0.8);
      pip->AddText("#pi");
      if (strstr(outputmode,"visual")) pip->Draw();
      TPaveText *kp = new TPaveText(0.5,-2.0,0.6,-1.8);
      kp->AddText("K");
      if (strstr(outputmode,"visual")) kp->Draw();
      TPaveText *prp = new TPaveText(0.9,-2.7,1.0,-2.5);
      prp->AddText("p");
      if (strstr(outputmode,"visual")) prp->Draw();
      //   TText *text2= new TText(.59,.06,"Momentum"); 
      //   text1->SetTextAngle(90);
      //   text1->Draw();
      //pidCanvas->DrawText(.1,.2,"Contamination               Efficiency");
      momvsmassCanvas->Update();
      if(strstr(outputsavemode,"asC")) momvsmassCanvas->Print("momvsmassCanvas.C");
      if(strstr(outputsavemode,"asEPS")) momvsmassCanvas->Print("momvsmassCanvas.eps");
      pidCanvas->cd();
      pidCanvas->Update();
      if(strstr(outputsavemode,"asC")) pidCanvas->Print("pidCanvas.C");
      if(strstr(outputsavemode,"asEPS")) pidCanvas->Print("pidCanvas.eps");
      char outFileName[100];
      strcpy(outFileName,"histos");
      strcat(outFileName,foutfileName);
      TFile *houtfile = new TFile(outFileName,"recreate");
      houtfile->cd();
      // saving canvas
      pidCanvas->Write(0,TObject::kOverwrite);
      momvsmassCanvas->Write(0,TObject::kOverwrite);
      // saving 1-D histos
      pit->Write(0,TObject::kOverwrite);
      pig->Write(0,TObject::kOverwrite);
      pieff->Write(0,TObject::kOverwrite);
      pieffls->Write(0,TObject::kOverwrite);
      picon->Write(0,TObject::kOverwrite);
      piid->Write(0,TObject::kOverwrite);
      pimatch->Write(0,TObject::kOverwrite);
      piall->Write(0,TObject::kOverwrite);
      pigen->Write(0,TObject::kOverwrite);
      kat->Write(0,TObject::kOverwrite);
      kag->Write(0,TObject::kOverwrite);
      kaeff->Write(0,TObject::kOverwrite);
      kaeffls->Write(0,TObject::kOverwrite);
      kaid->Write(0,TObject::kOverwrite);
      kamatch->Write(0,TObject::kOverwrite);
      kaall->Write(0,TObject::kOverwrite);
      kagen->Write(0,TObject::kOverwrite);
      kacon->Write(0,TObject::kOverwrite);
      prt->Write(0,TObject::kOverwrite);
      prg->Write(0,TObject::kOverwrite);
      preff->Write(0,TObject::kOverwrite);
      preffls->Write(0,TObject::kOverwrite);
      prcon->Write(0,TObject::kOverwrite);
      prid->Write(0,TObject::kOverwrite);
      prmatch->Write(0,TObject::kOverwrite);
      prall->Write(0,TObject::kOverwrite);
      prgen->Write(0,TObject::kOverwrite);
      // saving 2-D histos
      hpi->Write(0,TObject::kOverwrite);
      hka->Write(0,TObject::kOverwrite);
      hpr->Write(0,TObject::kOverwrite);
      // electron histos
      if (hel && eleff && elcon && elid && elall){
	hel->Write(0,TObject::kOverwrite);
	eleff->Write(0,TObject::kOverwrite);
	elcon->Write(0,TObject::kOverwrite);
	elid->Write(0,TObject::kOverwrite);
	elmatch->Write(0,TObject::kOverwrite);
	elall->Write(0,TObject::kOverwrite);
      }
      cout << "file " << houtfile << " has been created" << endl;
      cout << "it contains PID histos and canvas" << endl;
      houtfile->Close();
      houtfile->Write(0,TObject::kOverwrite);
    }
  }

  if (strstr(outputmode,"novisual")){
    // free used memory
    delete pit  ; pit=0;
    delete pig  ; pig=0;
    delete pieff; pieff=0;
    delete pieffls; pieffls=0;
    delete picon; picon=0;
    delete piid ; piid=0;
    delete pimatch; pimatch=0;
    delete piall; piall=0;
    delete pigen; pigen=0;
    delete kat  ; kat=0;
    delete kag  ; kag=0;
    delete kaeff; kaeff=0;
    delete kaeffls; kaeffls=0;
    delete kaid;  kaid=0;
    delete kamatch; kamatch=0;
    delete kaall; kaall=0;
    delete kagen; kagen=0;
    delete kacon; kacon=0;
    delete prt;   prt=0;
    delete prg;   prg=0;
    delete preff; preff=0;
    delete preffls; preffls=0;
    delete prcon; prcon=0;
    delete prid;  prid=0;
    delete prmatch; prmatch=0;
    delete prall; prall=0;
    delete prgen; prgen=0;
    // 2-D
    delete hpi; hpi=0;
    delete hka; hka=0;
    delete hpr; hpr=0;
    if (hel){ 
      delete hel;
      hel=0;
    }
    if (eleff){ 
      delete eleff;
      eleff=0;
    }
    if (eleffls){ 
      delete eleffls;
      eleffls=0;
    }
    if (elcon){ 
      delete elcon;
      elcon=0;
    }
    if (elid){ 
      delete elid;
      elid=0;
    }
    if (elmatch){ 
      delete elmatch;
      elmatch=0;
    }
    if (elall){ 
      delete elall;
      elall=0;
    }
  }
}


//__________________________________________________________________
Bool_t AliTOFPID::operator==( AliTOFPID const & /*tofrec*/)const
{
  // dummy version of Equal operator.
  // requested by coding conventions
  return kTRUE;

}
