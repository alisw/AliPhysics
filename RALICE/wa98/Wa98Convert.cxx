// $Id$

///////////////////////////////////////////////////////////////////////////
// Class Wa98Convert
// Conversion of Wa98 Hbook ntuple data into Wa98Event physics event structures.
//
// Usage example :
// ---------------
//
// gSystem->Load("libHbook");
// gSystem->Load("ralice");
// gSystem->Load("rwa98");
//
// // Input file with ANALYZ produced Hbook ntuple data
// THbookFile* f=new THbookFile("pb611258.cwn");
// TTree* h999=(TTree*)f->Get(999);
//
// // Output file for the event structures
// TFile* ofile=new TFile("run11258.root","RECREATE","WA98 data in RALICE event structure");
// TTree* otree=new TTree("T","Data of the 1996 Pb+Pb run");
//
// Int_t nentries=h999->GetEntries();
// cout << " Number of entries available : " << nentries << endl;
// cout << endl;
//
// // Limit the number of entries for testing
// nentries=300;
//
// // Print frequency to produce a short summary print every printfreq events
// Int_t printfreq=10;
//
// Convert q(h999);
// q.Loop(otree,nentries,printfreq);
// 
// otree->Print();
//
// // Close output file
// ofile->Write();
// ofile->Close();
//
//--- Author: Nick van Eijndhoven 06-jul-2004 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "Wa98Convert.h"
#include "Riostream.h"

ClassImp(Wa98Convert) // Class implementation to enable ROOT I/O

Wa98Convert::Wa98Convert(TTree* tree)
{
// Default constructor.
// Initialise the input tree (or chain) to be converted.
// By default tree=0;

 fChain=tree;
 if (!fChain) return;

 // Link the variables to the branches of the input tree/chain 
 fChain->SetBranchAddress("Jrun",&Jrun);
 fChain->SetBranchAddress("Jevt",&Jevt);
 fChain->SetBranchAddress("Jdate",&Jdate);
 fChain->SetBranchAddress("Jtime",&Jtime);
 fChain->SetBranchAddress("Jevid",&Jevid);
 fChain->SetBranchAddress("Jwscal",&Jwscal);
 fChain->SetBranchAddress("Itword",&Itword);
 fChain->SetBranchAddress("Zdc",&Zdc);
 fChain->SetBranchAddress("Emir",&Emir);
 fChain->SetBranchAddress("Emire",&Emire);
 fChain->SetBranchAddress("Emirh",&Emirh);
 fChain->SetBranchAddress("Etm",&Etm);
 fChain->SetBranchAddress("Etme",&Etme);
 fChain->SetBranchAddress("Etmh",&Etmh);
 fChain->SetBranchAddress("Nmod",&Nmod);
 fChain->SetBranchAddress("Irowl",Irowl);
 fChain->SetBranchAddress("Icoll",Icoll);
 fChain->SetBranchAddress("Adcl",Adcl);
 fChain->SetBranchAddress("Nclu",&Nclu);
 fChain->SetBranchAddress("Irowc",Irowc);
 fChain->SetBranchAddress("Icolc",Icolc);
 fChain->SetBranchAddress("Adcc",Adcc);
 fChain->SetBranchAddress("Ncluv",&Ncluv);
 fChain->SetBranchAddress("Iadccv",Iadccv);
 fChain->SetBranchAddress("Thetacv",Thetacv);
 fChain->SetBranchAddress("Phicv",Phicv);
}
///////////////////////////////////////////////////////////////////////////
Wa98Convert::~Wa98Convert()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void Wa98Convert::Loop(TTree* otree,Int_t nentries,Int_t printfreq)
{
// Loop over the specified number of entries and convert the 
// ntuple data into the Wa98Event structure.
// The output will be written on the output tree specified as "otree".
// If otree=0, a default standard output tree will be created.
// If nentries<0 (default) all the entries of the input chain
// will be processed.
// Every "printfreq" events a short event summary will be printed.
// The default value is printfreq=1.

 if (fChain==0) return;

 if (nentries<0) nentries=fChain->GetEntriesFast();

 if (!otree) otree=new TTree("T","Data of the 1996 Pb+Pb run");

 Double_t pi=acos(-1.);

 Double_t me=0.51099890221e-3;
 Double_t mpi=0.13956995;
 Double_t mkc=0.493677;
 Double_t mk0=0.497672;
 Double_t mp=0.93827231;
 Double_t mlam=1.115683;

 Wa98Event* evt=new Wa98Event();

 // Branch in the tree for the event structure
 Int_t split=1;
 Int_t bsize=32000;
 otree->Branch("Wa98Event","Wa98Event",&evt,bsize,split); 

 // The LEDA specific output data
 AliCalorimeter* ledaup=new AliCalorimeter(44,144);
 AliCalorimeter* ledalw=new AliCalorimeter(40,144);

 ledaup->SetName("LedaUp");
 ledalw->SetName("LedaDown");

 evt->InitLeda(ledaup);
 evt->InitLeda(ledalw);

 TDatime datim;
 Float_t pos[3],costh;
 AliSignal s;
 s.SetName("CPV signal ADC");

 for (Int_t jentry=0; jentry<nentries; jentry++)
 {
  fChain->GetEntry(jentry);

  // Reset the complete Event structure
  evt->Reset();

  evt->SetRunNumber(Jrun);
  evt->SetEventNumber(Jevt);
  datim.Set(Jdate,Jtime);
  evt->SetDayTime(datim);
  evt->SetProjectile(207,82,158);
  evt->SetTarget(207,82,0);
  evt->SetWeight(Jwscal);
  evt->SetTrig(Itword);
  evt->SetZdc(Zdc*1000.);
  evt->SetMiracE(1000.*Emir,Emire,Emirh);
  evt->SetMiracEt(Etm,Etme,Etmh);
 
  ledaup->Reset();
  ledalw->Reset();
  // Fill calorimeter with module data
  for (Int_t i=0; i<Nmod; i++)
  {
   if (Adcl[i] > 3) // Adc cut of 3 to remove noise
   {
    if (Irowl[i] > 0) ledaup->SetSignal(Irowl[i],Icoll[i],Adcl[i]);
    if (Irowl[i] < 0) ledalw->SetSignal(-Irowl[i],Icoll[i],Adcl[i]);
   }
  }

  // Store associated CPV signals
  for (Int_t j=0; j<Ncluv; j++)
  {
   s.Reset();
   s.SetSignal(Iadccv[j]);
   pos[1]=Thetacv[j]*pi/180.;
   pos[2]=Phicv[j]*pi/180.;
   costh=cos(pos[1]);
   pos[0]=0;
   if (costh) pos[0]=2103./costh;
   s.SetPosition(pos,"sph");
   pos[0]=0.4;
   pos[1]=2.2;
   pos[2]=0;
   s.SetPositionErrors(pos,"car");
   if (Phicv[j]>=0. && Phicv[j]<=180.)
   {
    ledaup->AddVetoSignal(s);
   }
   else
   {
    ledalw->AddVetoSignal(s);
   }
  }

  evt->AddDevice(ledaup);
  evt->AddDevice(ledalw);

  if (!(jentry%printfreq))
  {
//   cout << " Entry : " << jentry << " Run : " << Jrun << " Event : " << Jevt
//        << " Itword : " << Itword << " Etm : " << Etm << endl;
//   cout << " Jdate : " << Jdate << " Jtime : " << Jtime << endl; 
   cout << " Itword : " << Itword << " Nmod : " << Nmod << " Ncluv : " << Ncluv << endl;
   evt->HeaderData();
  }

  // Write the complete structure to the output Tree
  otree->Fill();
 }

 if (evt) delete evt;
 if (ledaup) delete ledaup;
 if (ledalw) delete ledalw;
}
///////////////////////////////////////////////////////////////////////////
