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

//_________________________________________________________________________
//
// This is a TTask that made the calculation of the Time zero using TOF.
// Description: The algorithm used to calculate the time zero of
// interaction using TOF detector is the following.
// We select in the MonteCarlo some primary particles - or tracks in
// the following - that strike the TOF detector (the larger part are
// pions, kaons or protons).
// We choose a set of 10 selected tracks, for each track You have the
// length of the track when the TOF is reached (a standard TOF hit
// does not contain this additional information, this is the reason
// why we implemented a new time zero dedicated TOF hit class
// AliTOFhitT0; in order to store this type of hit You have to use the
// AliTOFv4T0 as TOF class in Your Config.C. In AliTOFv4T0 the
// StepManager was modified in order to fill the TOF hit branch with
// this type of hits; in fact the AliTOF::AddT0Hit is called rather
// that the usual AliTOF::AddHit), the momentum at generation (from
// TreeK) and the time of flight given by the TOF detector.
// (Observe that the ctor of the AliTOF class, when the AliTOFv4T0
// class is used, is called with the "tzero" option: it is in order
// create the fHits TClonesArray filled with AliTOFhitT0 objects,
// rather than with normal AliTOFhit)
// Then Momentum and time of flight for each track are smeared
// according to known experimental resolution (all sources of error
// have been token into account).
// Let consider now only one set of 10 tracks (the algorithm is the
// same for all sets).
// Assuming the (mass) hypothesis that each track can be AUT a pion,
// AUT a kaon, AUT a proton, we consider all the 3 at 10 possible
// cases.
// For each track in each (mass) configuration
// (a configuration can be
// e.g. pion/pion/kaon/proton/pion/proton/kaon/kaon/pion/pion)
// we calculate the time zero (we know in fact the velocity of the
// track after the assumption about its mass, the time of flight given
// by the TOF, and the corresponding path travelled till the TOF
// detector). Then for each mass configuration we have 10 time zero
// and we can calculate the ChiSquare for the current configuration
// using the weighted mean over all 10 time zero.
// We call the best assignment the mass configuration that gives the
// minimum value of the ChiSquare.
// We plot the weighted mean over all 10 time zero for the best
// assignment, the ChiSquare for the best assignment and the
// corresponding confidence level.
// The strong assumption is the MC selection of primary particles. It
// will be introduced in the future also some more realistic
// simulation about this point.
//
// Use case:
// root [0] AliTOFT0 * tzero = new AliTOFT0("galice.root")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root [1] tzero->ExecuteTask()
// root [2] tzero->ExecuteTask("tim")
//             // available parameters:
//             tim - print benchmarking information
//             all - print usefull informations about the number of
//                   misidentified tracks and a comparison about the
//                   true configuration (known from MC) and the best
//                   assignment
//
//-- Author: F. Pierella
//
//_________________________________________________________________________

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TFolder.h>
#include <TFrame.h>
#include <TH1.h>
#include <TParticle.h>
#include <TBenchmark.h>
#include <TTask.h>
#include <TTree.h>
#include <TRandom.h>
#include <TROOT.h>

#include "AliMC.h"
#include "AliRun.h"

#include "AliTOFhitT0.h"
#include "AliTOFT0.h"
#include "AliTOF.h"

ClassImp(AliTOFT0)

//____________________________________________________________________________ 
AliTOFT0::AliTOFT0():
  TTask("AliTOFT0",""),
  fNevents(0),
  fTimeResolution(0),
  fLowerMomBound(0),
  fUpperMomBound(0),
  fT0File(""),
  fHeadersFile("")
{
  // ctor
}
           
//____________________________________________________________________________ 
AliTOFT0::AliTOFT0(char* headerFile, Int_t nEvents):
  TTask("AliTOFT0",""), 
  fNevents(nEvents),
  fTimeResolution(1.2e-10),
  fLowerMomBound(1.5),
  fUpperMomBound(2.),
  fT0File(""),
  fHeadersFile(headerFile)
{
  //
  //
  //

  //fNevents=nEvents ; // Number of events for which calculate the T0, 
                     // default 0: it means all evens in current file
  //fLowerMomBound=1.5; // [GeV/c] default value
  //fUpperMomBound=2. ; // [GeV/c] default value
  //fTimeResolution   = 1.2e-10; // 120 ps by default	
  //fHeadersFile = headerFile ;

  TFile * file = (TFile*) gROOT->GetFile(fHeadersFile.Data() ) ;

  //File was not opened yet
  if(file == 0){
    if(fHeadersFile.Contains("rfio"))
      file =	TFile::Open(fHeadersFile,"update") ;
    else
      file = new TFile(fHeadersFile.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }

  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}

//____________________________________________________________________________ 
AliTOFT0::AliTOFT0(const AliTOFT0 & tzero):
  TTask("AliTOFT0",""),
  fNevents(0),
  fTimeResolution(0),
  fLowerMomBound(0),
  fUpperMomBound(0),
  fT0File(""),
  fHeadersFile("")
{
  // copy ctr

( (AliTOFT0 &)tzero ).Copy(*this);
}

//____________________________________________________________________________ 
  AliTOFT0::~AliTOFT0()
{
  // dtor
}

//____________________________________________________________________________
void AliTOFT0::Exec(Option_t *option) 
{ 
  //
  // calculate T0 distribution for all events using chisquare 
  //
  Int_t ngood=0;
  Int_t nmisidentified=0;
  Int_t nmisidentified0=0;
  Int_t nmisidentified1=0;
  Int_t nmisidentified2=0;
  Int_t nmisidentified3=0;
  Int_t nmisidentified4=0;
  Int_t nmisidentified5=0;
  Int_t nmisidentified6=0;
  Int_t nmisidentified7=0;
  Int_t nmisidentified8=0;
  Int_t nmisidentified9=0;
  Int_t ipartold = -1;
  Int_t ipart;
  Int_t selected=0;
  Int_t istop=0;
  Float_t timeresolutioninns=fTimeResolution*(1.e+9); // convert in [ns]
  const Int_t kUPDATE = 5; // for visual option
  Int_t itimes=0;
  TCanvas* c1=0;
  TCanvas* c2=0;
  TCanvas* c3=0;

  if(strstr(option,"visual")){
    // Create a new canvas.
    //c1 = new TCanvas("c1","Dynamic Visual Filling of time zero histo",10,10,500,500);
    c1 = new TCanvas("c1","Dynamic Visual Filling of time zero histo",10,10,370,370);
    c1->SetFillColor(35);
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(6);
    c1->GetFrame()->SetBorderMode(-1);

    //c2 = new TCanvas("c2","Dynamic Visual Filling of chisquare histo",550,10,500,500);
    c2 = new TCanvas("c2","Dynamic Visual Filling of chisquare histo",380,10,370,370);
    c2->SetFillColor(35);
    c2->GetFrame()->SetFillColor(21);
    c2->GetFrame()->SetBorderSize(6);
    c2->GetFrame()->SetBorderMode(-1);

    //c3 = new TCanvas("c3","Dynamic Visual Filling of confidence level histo",280,550,500,500);
    c3 = new TCanvas("c3","Dynamic Visual Filling of confidence level histo",760,10,370,370);
    c3->SetFillColor(35);
    c3->GetFrame()->SetFillColor(21);
    c3->GetFrame()->SetBorderSize(6);
    c3->GetFrame()->SetBorderMode(-1);
  }

  if(strstr(option,"tim") || strstr(option,"all"))
    gBenchmark->Start("TOFT0");

  TH1F *htzerobest= new TH1F("htzerobest","T0 for best assignment",200,-1.,1.);
  TH1F* hchibest  = new TH1F("hchibest","ChiSquare Min Distribution",80,0.,40.);
  TH1F* hchibestconflevel  = new TH1F("hchibestconflevel","ChiSquare Min Confidence Level",10,0.,1.);

  // setting histo colors
  if(strstr(option,"visual")){
    htzerobest->SetFillColor(48);
    hchibest->SetFillColor(50);
    hchibestconflevel->SetFillColor(52);
  }

  Int_t   assparticle[10]={3,3,3,3,3,3,3,3,3,3};
  Int_t   truparticle[10]={3,3,3,3,3,3,3,3,3,3};
  Float_t t0best=999.;
  Float_t timeofflight[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t momentum[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t timezero[10];
  Float_t weightedtimezero[10];
  Float_t beta[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t sqMomError[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t sqTrackError[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t massarray[3]={0.13957,0.493677,0.9382723};
  Float_t dummychisquare=0.;
  Float_t chisquare=999.;
  Float_t tracktoflen[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  AliTOF *detTOF = (AliTOF *) gAlice->GetDetector ("TOF");

  if (!detTOF) {
    Error("AliTOFT0","TOF not found");
    return;
  }

  if(strstr(option,"all")){
    AliInfo(Form("Selecting primary tracks with momentum between %d GeV/c and %d GeV/c",  fLowerMomBound, fUpperMomBound));
    AliInfo("Memorandum: 0 means PION | 1 means KAON | 2 means PROTON")
  }

  if (fNevents == 0) fNevents = (Int_t) gAlice->TreeE()->GetEntries();

  for (Int_t ievent = 0; ievent < fNevents; ievent++) {
    gAlice->GetEvent(ievent);
    TTree *hitTree = detTOF->TreeH ();
    if (!hitTree)
      return;
    TParticle*    particle;
    AliTOFhitT0*  tofHit;
    TClonesArray* tofHits = detTOF->Hits();

    Int_t lasttrack=-1;
    Int_t nset=0;

    hitTree->SetBranchStatus("*",0); // switch off all branches
    hitTree->SetBranchStatus("TOF*",1); // switch on only TOF

    // Start loop on primary tracks in the hits containers

    Int_t ntracks = static_cast<Int_t>(hitTree->GetEntries());
    for (Int_t track = 0; track < ntracks; track++)
    {
      if(nset>=5) break; // check on the number of set analyzed
      
      gAlice->ResetHits();
      hitTree->GetEvent(track);

      AliMC *mcApplication = (AliMC*)gAlice->GetMCApp();

      particle = mcApplication->Particle(track);
      Int_t nhits = tofHits->GetEntriesFast();

      for (Int_t hit = 0; hit < nhits; hit++)
      {
	tofHit = (AliTOFhitT0 *) tofHits->UncheckedAt(hit);
	ipart    = tofHit->GetTrack();
	// check to discard the case when the same particle is selected more than one
	// time 

	if (ipart != ipartold){
	  
	  particle = (TParticle*)gAlice->GetMCApp()->Particle(ipart);
	  
	  Float_t idealtime=tofHit->GetTof();
	  //	   Float_t time=idealtime;
	  Float_t time   = gRandom->Gaus(idealtime, fTimeResolution);
	  Float_t toflen=tofHit->GetLen();
	  toflen=toflen/100.; // toflen given in m
	  Int_t pdg   = particle->GetPdgCode();
	  Int_t abspdg   =TMath::Abs(pdg);
	  Float_t idealmom  = particle->P();
	  Float_t momres=idealmom*0.025; // 2.5% res token into account for all momenta
	  Float_t mom =gRandom->Gaus(idealmom,momres);

	  Bool_t isgoodpart=(abspdg==211 || abspdg==2212 || abspdg==321);

	  time*=1.E+9; // tof given in nanoseconds	   
	  if (particle->GetFirstMother() < 0 && isgoodpart && mom<=fUpperMomBound && mom>=fLowerMomBound){
	    selected+=1;
	    istop=selected;
	    if(istop>10) break;
	    Int_t index=selected-1;
	    timeofflight[index]=time;
	    tracktoflen[index]=toflen;
	    momentum[index]=mom;
	    //AliInfo(Form(" %d  %d  %d ", timeofflight[index], tracktoflen[index], momentum[index]));
	    switch (abspdg) {
	    case 211:
	      truparticle[index]=0;
	      break ;
	    case 321:
	      truparticle[index]=1;
	      break ;
	    case 2212: 
	      truparticle[index]=2;
	      break ;
	    }
	    
	  }
	  ipartold = ipart;
	  
	  if(istop==10){ // start analysis on current set
	    nset+=1;
	    lasttrack=track;
	    istop=0;
	    selected=0;
	    //AliInfo("starting t0 calculation for current set");
	    for (Int_t i1=0; i1<3;i1++) {
	      beta[0]=momentum[0]/sqrt(massarray[i1]*massarray[i1]+momentum[0]*momentum[0]);
	      for (Int_t i2=0; i2<3;i2++) { 
		beta[1]=momentum[1]/sqrt(massarray[i2]*massarray[i2]+momentum[1]*momentum[1]);
		for (Int_t i3=0; i3<3;i3++) {
		  beta[2]=momentum[2]/sqrt(massarray[i3]*massarray[i3]+momentum[2]*momentum[2]);
		  for (Int_t i4=0; i4<3;i4++) {
		    beta[3]=momentum[3]/sqrt(massarray[i4]*massarray[i4]+momentum[3]*momentum[3]);
		    for (Int_t i5=0; i5<3;i5++) {
		      beta[4]=momentum[4]/sqrt(massarray[i5]*massarray[i5]+momentum[4]*momentum[4]);
		      for (Int_t i6=0; i6<3;i6++) {
			beta[5]=momentum[5]/sqrt(massarray[i6]*massarray[i6]+momentum[5]*momentum[5]);
			for (Int_t i7=0; i7<3;i7++) { 
			  beta[6]=momentum[6]/sqrt(massarray[i7]*massarray[i7]+momentum[6]*momentum[6]);
			  for (Int_t i8=0; i8<3;i8++) {
			    beta[7]=momentum[7]/sqrt(massarray[i8]*massarray[i8]+momentum[7]*momentum[7]);
			    for (Int_t i9=0; i9<3;i9++) {
			      beta[8]=momentum[8]/sqrt(massarray[i9]*massarray[i9]+momentum[8]*momentum[8]);
			      for (Int_t i10=0; i10<3;i10++) { 	
				beta[9]=momentum[9]/sqrt(massarray[i10]*massarray[i10]+momentum[9]*momentum[9]);
				
				Float_t meantzero=0.;
				Float_t sumAllweights=0.;
				for (Int_t itz=0; itz<10;itz++) {
				  sqMomError[itz]=((1.-beta[itz]*beta[itz])*0.025)*((1.-beta[itz]*beta[itz])*0.025)*(tracktoflen[itz]/(0.299792*beta[itz]))*(tracktoflen[itz]/(0.299792*beta[itz])); // this gives the square of the momentum error in nanoseconds
				  sqTrackError[itz]=(timeresolutioninns*timeresolutioninns+sqMomError[itz]); // total error for the current track
				  sumAllweights+=1./sqTrackError[itz];

				  timezero[itz]=(tracktoflen[itz]/(beta[itz]*0.299792))-timeofflight[itz];
				  weightedtimezero[itz]=((tracktoflen[itz]/(beta[itz]*0.299792))-timeofflight[itz])/sqTrackError[itz];// weighted time zero for current track
				  meantzero+=weightedtimezero[itz];
				} // end loop for (Int_t itz=0; itz<10;itz++)
				meantzero=meantzero/sumAllweights; // it is given in [ns]
				
				dummychisquare=0.;
				// calculate the chisquare for the current assignment
				for (Int_t icsq=0; icsq<10;icsq++) {
				  dummychisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/sqTrackError[icsq];
				} // end loop for (Int_t icsq=0; icsq<10;icsq++) 

				if(dummychisquare<=chisquare){
				  assparticle[0]=i1;
				  assparticle[1]=i2;
				  assparticle[2]=i3;
				  assparticle[3]=i4;
				  assparticle[4]=i5;
				  assparticle[5]=i6;
				  assparticle[6]=i7;
				  assparticle[7]=i8;
				  assparticle[8]=i9;
				  assparticle[9]=i10;
				  chisquare=dummychisquare;
				  t0best=meantzero;
				} // close if(dummychisquare<=chisquare)
				
			      } // end loop on i10
			    } // end loop on i9
			  } // end loop on i8
			} // end loop on i7
		      } // end loop on i6
		    } // end loop on i5
		  } // end loop on i4
		} // end loop on i3
	      } // end loop on i2
	    } // end loop on i1

	    if(truparticle[0]==assparticle[0] && truparticle[1]==assparticle[1] && truparticle[2]==assparticle[2]  && truparticle[3]==assparticle[3] && truparticle[4]==assparticle[4]&& truparticle[5]==assparticle[5] && truparticle[6]==assparticle[6] && truparticle[7]==assparticle[7]  && truparticle[8]==assparticle[8] && truparticle[9]==assparticle[9]) ngood+=1;
	    if(truparticle[0]!=assparticle[0]) nmisidentified0+=1;
	    if(truparticle[1]!=assparticle[1]) nmisidentified1+=1;
	    if(truparticle[2]!=assparticle[2]) nmisidentified2+=1;
	    if(truparticle[3]!=assparticle[3]) nmisidentified3+=1;
	    if(truparticle[4]!=assparticle[4]) nmisidentified4+=1;
	    if(truparticle[5]!=assparticle[5]) nmisidentified5+=1;
	    if(truparticle[6]!=assparticle[6]) nmisidentified6+=1;
	    if(truparticle[7]!=assparticle[7]) nmisidentified7+=1;
	    if(truparticle[8]!=assparticle[8]) nmisidentified8+=1;
	    if(truparticle[9]!=assparticle[9]) nmisidentified9+=1;
	    // filling histos
	    htzerobest->Fill(t0best);
	    hchibest->Fill(chisquare);
	    Double_t dblechisquare=(Double_t)chisquare;
	    Float_t confLevel=(Float_t)TMath::Prob(dblechisquare,9); // ndf 10-1=9
	    hchibestconflevel->Fill(confLevel);
	    itimes++;
	    if(strstr(option,"all")){
	      AliInfo(Form("True Assignment %d  %d  %d  %d  %d  %d  %d  %d  %d  %d", truparticle[0], truparticle[1], truparticle[2], truparticle[3], truparticle[4], truparticle[5], truparticle[6], truparticle[7], truparticle[8], truparticle[9]));
	      AliInfo(Form("Best Assignment %d  %d  %d  %d  %d  %d  %d  %d  %d  %d", assparticle[0], assparticle[1], assparticle[2], assparticle[3], assparticle[4], assparticle[5], assparticle[6], assparticle[7], assparticle[8], assparticle[9]));
	      AliInfo(Form("Minimum ChiSquare for current set   %d ", chisquare));
	      AliInfo(Form("Confidence Level (Minimum ChiSquare) %d", confLevel));
	    }
	    if (strstr(option,"visual") && itimes && (itimes%kUPDATE) == 0) {
	      if (itimes == kUPDATE){
		c1->cd();
		htzerobest->Draw();
		c2->cd();
		hchibest->Draw();
		c3->cd();
		hchibestconflevel->Draw();
	      }
	      c1->Modified();
	      c1->Update();
	      c2->Modified();
	      c2->Update();
	      c3->Modified();
	      c3->Update();
	      if (gSystem->ProcessEvents())
		break;
	    }
	    chisquare=999.;
	    t0best=999.;
	    
	  } // end for the current set. close if(istop==5)
	} // end condition on ipartold
      } // end loop on hits for the current track
      if(istop>=10) break;
    } // end loop on ntracks  
  } //event loop
  
  if(strstr(option,"all")){
    nmisidentified=(nmisidentified0+nmisidentified1+nmisidentified2+nmisidentified3+nmisidentified4+nmisidentified5+nmisidentified6+nmisidentified7+nmisidentified8+nmisidentified9);
    AliInfo(Form("total number of tracks token into account  %i", 10*5*fNevents));
    Float_t badPercentage=100.*(Float_t)nmisidentified/(10*5*fNevents);
    AliInfo(Form("total misidentified                       %i (%d %) ", nmisidentified, badPercentage));
    AliInfo(Form("Total Number of set token into account     %i", 5*fNevents));
    Float_t goodSetPercentage=100.*(Float_t)ngood/(5*fNevents);
    AliInfo(Form("Number of set with no misidentified tracks %i (%d %)", ngood, goodSetPercentage));
  }

  // free used memory for canvas
  delete c1; c1=0;
  delete c2; c2=0;
  delete c3; c3=0;

  // generating output filename only if not previously specified using SetTZeroFile
  char outFileName[70];
  strcpy(outFileName,"ht010tr120ps"); // global time resolution has to be converted from Int_t to char
                                      // in order to have in the output filename this parameter
  strcat(outFileName,fHeadersFile);

  if(fT0File.IsNull()) fT0File=outFileName;

  TFile* houtfile = new TFile(fT0File,"recreate");
  houtfile->cd();
  htzerobest->Write(0,TObject::kOverwrite);
  hchibest->Write(0,TObject::kOverwrite);
  hchibestconflevel->Write(0,TObject::kOverwrite);
  houtfile->Close();  
  
  
  if(strstr(option,"tim") || strstr(option,"all")){
    gBenchmark->Stop("TOFT0");
    AliInfo("AliTOFT0:");
    /*
    cout << "   took " << gBenchmark->GetCpuTime("TOFT0") << " seconds in order to calculate T0 " 
	 << gBenchmark->GetCpuTime("TOFT0")/fNevents << " seconds per event " << endl ;
    */
    gBenchmark->Print("TOFT0");
  }
}
 
//__________________________________________________________________
void AliTOFT0::SetTZeroFile(char * file )
{
  //
  //
  //
  printf("Destination file : %s \n", file) ;
  fT0File=file;

}

//__________________________________________________________________
void AliTOFT0::Print(Option_t* /*option*/)const
{
  //
  //
  //
  printf("------------------- %s -------------\n", GetName()) ;
  if(!fT0File.IsNull())
    printf("  Writing T0 Distribution to file  %s \n",(char*) fT0File.Data());

}

//__________________________________________________________________
Bool_t AliTOFT0::operator==( AliTOFT0 const &tzero )const
{
  //
  // Equal operator
  // 

  if( (fTimeResolution==tzero.fTimeResolution)&&(fLowerMomBound==tzero.fLowerMomBound)&&(fUpperMomBound==tzero.fUpperMomBound))
    return kTRUE ;
  else
    return kFALSE ;
}
