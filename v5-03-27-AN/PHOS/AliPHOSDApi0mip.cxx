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

// --
// --
// Implementation for TTree output in PHOS DA
// for calibrating energy by pi0 and MIP.
// --
// -- Author: Hisayuki Torii (Hiroshima Univ.)
// --

#include <stdio.h>
#include <iostream>
#include <math.h>
#include "AliLog.h"
#include "AliPHOSDApi0mip.h"
ClassImp(AliPHOSDApi0mip)
//----------------------------------------------------------------
AliPHOSDApi0mip::AliPHOSDApi0mip(int module,int iterid,const char* fopt) :
  TNamed(), fCreateTree(false), fCreateHist(false), fMod(0), fIterId(0),
  fTFile(0), fTTree(0), fEvent(0), fEventClustered(0), fTime(0),
  fH1Time(0), fH1DigitNum(0), fH1ClusterNum(0),
  fH2EneDigitId(0), fH2MipDigitId(0), fH2Pi0DigitId(0), fH3Pi0AsymPt(0) {
  // Constructor

  char hname[1024], htitle[1024];

  TString shname="AliPHOSDApi0mip_mod%d_iter%d";
  snprintf(hname,shname.Length(),shname.Data(),module,iterid);
  SetName(hname);

  TString shtitle="PHOS MIP/pi0 Calibration DA for Module:%d Iteration:%d";
  snprintf(htitle,shtitle.Length(),shtitle.Data(),module,iterid);
  SetTitle(htitle);

  fMod = module;
  fIterId = iterid;
  fEvent = 0;
  fEventClustered = false;
  fTime = 0;
  fCreateTree = false;
  fCreateHist = false;

  char fname[1024];

  TString sfname="AliPHOSDApi0mip_mod%d.root";
  snprintf(fname,sfname.Length(),sfname.Data(),module);
  fTFile = TFile::Open(fname,fopt);
}
//----------------------------------------------------------------
AliPHOSDApi0mip::AliPHOSDApi0mip(const AliPHOSDApi0mip& da):
  TNamed(da), fCreateTree(false), fCreateHist(false), fMod(da.fMod), fIterId(da.fIterId),
  fTFile(0), fTTree(0), fEvent(0), fEventClustered(0), fTime(0),
  fH1Time(0), fH1DigitNum(0), fH1ClusterNum(0),
  fH2EneDigitId(0), fH2MipDigitId(0), fH2Pi0DigitId(0), fH3Pi0AsymPt(0) {
  // Copy Constructor

  char fname[1024], hname[1024], htitle[1024];;

  TString sfname="%s.root";
  snprintf(fname,sfname.Length(),sfname.Data(),GetName());
  fTFile = TFile::Open(fname,"RECREATE");
  fEvent = 0;
  fEventClustered = false;
  fTime = 0;

  if( da.fCreateHist ){
    fCreateHist = true;
    fH1Time = (TH1I*) da.fH1Time->Clone();
    fH1DigitNum = (TH1F*) da.fH1DigitNum->Clone();
    fH1ClusterNum = (TH1F*) da.fH1ClusterNum->Clone();
    fH2EneDigitId = (TH2F*) da.fH2EneDigitId->Clone();
    fH2MipDigitId = (TH2F*) da.fH2MipDigitId->Clone();
    fH2Pi0DigitId = (TH2F*) da.fH2Pi0DigitId->Clone();
    fH3Pi0AsymPt = (TH3F*) da.fH3Pi0AsymPt->Clone();
  } else {
    fCreateHist = false;
    fH1Time = 0;
    fH1DigitNum = 0;
    fH1ClusterNum = 0;
    fH2EneDigitId = 0;
    fH2MipDigitId = 0;
    fH2Pi0DigitId = 0;
    fH3Pi0AsymPt = 0;
  }
  if( da.fCreateTree ){
    // Create new ttree instead of copy.
    fCreateTree = true;

    TString shname="tevt_mod%d_iter%d";
    snprintf(hname,shname.Length(),shname.Data(),fMod,fIterId);

    TString shtitle="Calibration for Module:%d Iteration:%d";
    snprintf(htitle,shtitle.Length(),shtitle.Data(),fMod,fIterId);

    fTTree = new TTree(hname,htitle);
    fTTree->Branch("AliPHOSDATreeEvent","AliPHOSDATreeEvent",&fEvent);
  } else {
    fCreateTree = false;
    fTTree = 0;
  }
}
//-------------------------------------------------------------------
AliPHOSDApi0mip& AliPHOSDApi0mip::operator= (const AliPHOSDApi0mip&)
{
  // Operator= is not implemented yet
  AliFatal("Operator= is not implemented");  
  return *this;
}
//-------------------------------------------------------------------
AliPHOSDApi0mip::~AliPHOSDApi0mip(){
  // Destructor

  //std::cout<<" DEBUGDEBUG :: AliPHOSDApi0mip:: deconstructor "<<std::endl;
  if( fCreateHist ){
    fH1Time->Write(0,TObject::kOverwrite);
    fH1DigitNum->Write(0,TObject::kOverwrite);
    fH1ClusterNum->Write(0,TObject::kOverwrite);
    fH2EneDigitId->Write(0,TObject::kOverwrite);
    fH2MipDigitId->Write(0,TObject::kOverwrite);
    fH2Pi0DigitId->Write(0,TObject::kOverwrite);
    fH3Pi0AsymPt->Write(0,TObject::kOverwrite);
    delete fH1Time;
    delete fH1DigitNum;
    delete fH1ClusterNum;
    delete fH2EneDigitId;
    delete fH2MipDigitId;
    delete fH2Pi0DigitId;
    delete fH3Pi0AsymPt;
  }
  if( fCreateTree ){
    fTTree->Write();
  }
  if( fEvent ) delete fEvent;
  fTFile->Close();
  delete fTFile;
}
//-------------------------------------------------------------------
void AliPHOSDApi0mip::FillDigit(float adc,int row,int col){
  // Put new digit information in event
  //

  if( fEvent==0 ) fEvent = new AliPHOSDATreeEvent;
  float energy = adc * 0.005; // 5MeV/ch is fixed.
  if( energy > 0.020 ){  // Threshold at 20MeV.
    //std::cout<<" AliPHOSDApi0mip::AppendDigit() DEBUG:: (row,col)=("<<row<<","<<col<<") adc="<<adc<<" energy="<<energy<<std::endl;
    fEvent->Fill(energy,row,col);
  }
}
//-------------------------------------------------------------------
void AliPHOSDApi0mip::NewEvent(){
  // Initiator must be called in the beginning of event
  //

  if(fEvent) fEvent->Reset();
  fEventClustered = false;
}
//----------------------------------------------------------------
bool AliPHOSDApi0mip::CreateTree(){
  if(! fTFile ) {
    return false;
  }
  fTFile->cd();
  char hname[100], htitle[100];
  snprintf(hname,100,"trevt_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"Calibration for Module:%d Iteration:%d",fMod,fIterId);
  fTTree = new TTree(hname,htitle);
  fTTree->Branch("AliPHOSDATreeEvent","AliPHOSDATreeEvent",&fEvent);
  fCreateTree = true;
  return true;
}
//----------------------------------------------------------------
bool AliPHOSDApi0mip::CreateHist(){
  // Create histogram routine called by the constructor only.

  char hname[100], htitle[100];

  snprintf(hname,100,"h1_time_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"Time : Mod:%d Iter:%d",fMod,fIterId);
  fH1Time = (TH1I*) gDirectory->Get(hname);
  if( fH1Time>0 ){
    std::cout<<" AliPHOSDApi0mip:Warning!! Output object already exist : "<<fH1Time->GetName()<<std::endl;
    return false;
  }
  fH1Time = new TH1I(hname,htitle,2,0,2);

  snprintf(hname,100,"h1_digitnum_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"Number of Digits : Mod:%d Iter:%d",fMod,fIterId);
  fH1DigitNum = new TH1F(hname,htitle,100,0,100);

  snprintf(hname,100,"h1_clusternum_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"Number of Clusters : Mod:%d Iter:%d",fMod,fIterId);
  fH1ClusterNum = new TH1F(hname,htitle,100,0,100);

  snprintf(hname,100,"h2_pi0digitid_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"PHOS pi0 mass vs Digit Id : Mod:%d Iter:%d",fMod,fIterId);
  fH2Pi0DigitId = new TH2F(hname,htitle,17920,0,17920,120,0,0.3);
      
  snprintf(hname,100,"h2_mipdigitid_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"PHOS MIP vs Digit Id : Mod:%d Iter:%d",fMod,fIterId);
  fH2MipDigitId = new TH2F(hname,htitle,17920,0,17920,100,0-0.0025,0.5-0.0025);
  //fH2MipDigitId = new TH2F(hname,htitle,17920,0,17920,50,0.5,1.5);

  snprintf(hname,100,"h2_enedigitid_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"PHOS MIP vs Digit Id : Mod:%d Iter:%d",fMod,fIterId);
  //fH2EneDigitId = new TH2F(hname,htitle,17920,0,17920,100,0-0.0025,0.5-0.0025);
  fH2EneDigitId = new TH2F(hname,htitle,17920,0,17920,50,0,5.0);

  snprintf(hname,100,"h3_pi0asympt_mod%d_iter%d",fMod,fIterId);
  snprintf(htitle,100,"PHOS pi0 mass vs pT vs Asym : Mod:%d Iter:%d",fMod,fIterId);
  fH3Pi0AsymPt = new TH3F(hname,htitle,20,0,1,20,0,10,200,0,1);

  fCreateHist = true;
  return true;
}
//-------------------------------------------------------------------
void AliPHOSDApi0mip::FillTree(AliPHOSDATreeEvent* event){
  // Fill all information (AliPHOSDATreeEvent) into TTree
  // if event==NULL, the obtained information (=fEvent) is used.

  if(! fCreateTree ) CreateTree();
  if( event == 0 ){
    if( fEvent == 0 ){
      fEvent = new AliPHOSDATreeEvent();
    }
    // No need to make clustering for tree saving.
    //if(! fEventClustered ){
    //fEvent->ExecuteClustering();
    //fEventClustered = true;
    //}
    fEvent->SetTime(fTime);
    fTTree->Fill();
  } else {
    if( fEvent ) delete fEvent;
    fEvent = event;
    fTTree->Fill();
    fEvent = 0;
  }
}
//-------------------------------------------------------------------
void AliPHOSDApi0mip::FillHist(AliPHOSDATreeEvent* event){
  // Fill all information (AliPHOSDATreeEvent) into histograms
  // if event==NULL, the obtained information (=fEvent) is used.

  if(! fCreateHist ) CreateHist();
  if( event == 0 ){
    if( fEvent==0 ) fEvent = new AliPHOSDATreeEvent;
    if(! fEventClustered ){
      fEvent->ExecuteClustering();
      fEventClustered = true;
    }
    fEvent->SetTime(fTime);
    event = fEvent;
  }

  // -- Constant
  int nclt1, nclt2;
  int ndigits;
  int ndigitsall = 0;
  AliPHOSDATreeDigit digit1, digit2;
  const float kDistanceFromVertex = 460.;
  float dist1, dist2, px, py;
  float mass, cosine, pt, asym;

  // -- Filling into Time Histogram
  double start = fH1Time->GetBinContent(1);
  double stop = fH1Time->GetBinContent(2);
  double currenttime = (double) (event->GetTime());
  if( start == 0 && stop == 0 ) start = stop = currenttime;
  if( currenttime < start ) start = currenttime;
  if( currenttime > stop ) stop = currenttime;
  fH1Time->SetBinContent(1,start);
  fH1Time->SetBinContent(2,stop);

  // -- Start looping clusters
  nclt1 = event->GetNClusters();
  fH1ClusterNum->Fill(nclt1);
  //std::cout<<" AliPHOSDApi0mip::FillHisto() DEBUG:: Number of clusters: "<<nclt1<<std::endl;
  while( nclt1-- ){
    AliPHOSDATreeCluster& clt1 = event->GetCluster(nclt1);
    ndigits = clt1.GetNDigits();
    while( ndigits-- ){
      AliPHOSDATreeDigit& digit = clt1.GetDigit(ndigits);
      fH2EneDigitId->Fill(digit.GetDigitId(),digit.GetEnergy());
      ndigitsall++;
    }
    digit1 = clt1.GetMaxDigit();
    if( clt1.GetEnergy() > 0 && clt1.GetNDigits() >= 2 ){
      fH2MipDigitId->Fill(digit1.GetDigitId(),clt1.GetEnergy());
    }
    if( clt1.GetEnergy() > 0.2 ){
      nclt2 = nclt1;
      while( nclt2-- ){
	AliPHOSDATreeCluster& clt2 = event->GetCluster(nclt2);
	digit2 = clt2.GetMaxDigit();
	if( clt2.GetEnergy() > 0.2 ){
	  // 
	  // Following approximation is applied in order to avoid many calling "sqrt" function.
	  // sqrt( 1 + x^2 ) =~ 1 + 0.5x^2 (for x<<1)
	  //
	  //cosine = ((clt1.GetRow()-31.5)*(clt2.GetRow()-31.5)*4.84 + (clt1.GetCol()-27.5)*(clt2.GetCol()-27.5)*4.84 + kDistanceFromVertex*kDistanceFromVertex)
	  //	    / ( kDistanceFromVertex + (clt1.GetRow()-31.5)*(clt1.GetRow()-31.5)*2.42/kDistanceFromVertex + (clt1.GetCol()-27.5)*(clt1.GetCol()-27.5)*2.42/kDistanceFromVertex )
	  //	    / ( kDistanceFromVertex + (clt2.GetRow()-31.5)*(clt2.GetRow()-31.5)*2.42/kDistanceFromVertex + (clt2.GetCol()-27.5)*(clt2.GetCol()-27.5)*2.42/kDistanceFromVertex );
	  //
	  

	  dist1 = sqrt( kDistanceFromVertex*kDistanceFromVertex + (clt1.GetRow()-31.5)*(clt1.GetRow()-31.5)*4.84 + (clt1.GetCol()-27.5)*(clt1.GetCol()-27.5)*4.84 );
	  dist2 = sqrt( kDistanceFromVertex*kDistanceFromVertex + (clt2.GetRow()-31.5)*(clt2.GetRow()-31.5)*4.84 + (clt2.GetCol()-27.5)*(clt2.GetCol()-27.5)*4.84 );
	  cosine = ((clt1.GetRow()-31.5)*(clt2.GetRow()-31.5)*4.84 + (clt1.GetCol()-27.5)*(clt2.GetCol()-27.5)*4.84 + kDistanceFromVertex*kDistanceFromVertex)
	    / dist1 / dist2;
	  mass = sqrt( 2.0 * clt1.GetEnergy() * clt2.GetEnergy() * (1 - cosine) );
	  //std::cout<<" DEBUG: dist1:"<<dist1<<" dist2:"<<dist2<<" cosine:"<<cosine<<" c1:"<<clt1.GetEnergy()<<" c2:"<<clt2.GetEnergy()<<" mass:"<<mass<<std::endl;
	  //std::cout<<"      : clt1(row,col)=("<<clt1.GetRow()<<","<<clt1.GetCol()<<") clt2(row,col)=("<<clt2.GetRow()<<","<<clt2.GetCol()<<") "<<std::endl;
	  asym = fabs(clt1.GetEnergy() - clt2.GetEnergy())/(clt1.GetEnergy() + clt2.GetEnergy());
	  py = clt1.GetEnergy() * kDistanceFromVertex / dist1 + clt2.GetEnergy() * kDistanceFromVertex / dist2;
	  px = clt1.GetEnergy() * (clt1.GetCol()-27.5)*2.42/dist1 + clt2.GetEnergy() * (clt2.GetCol()-27.5)*2.42/ dist2;
	  pt = sqrt(px*px+py*py);
	  fH3Pi0AsymPt->Fill(asym,pt,mass);
	  if( digit1.IsValid() && digit1.IsValid() &&
	      digit1.GetDigitId() != digit2.GetDigitId() &&
	      //pt > 1.5 && clt1.GetEnergy() > 0.3 && clt2.GetEnergy() > 0.3 ){
	      //pt > 1.3 && clt1.GetEnergy() > 0.25 && clt2.GetEnergy() > 0.25 ){
	      pt > 1.0 && clt1.GetEnergy() > 0.2 && clt2.GetEnergy() > 0.2 ){
	    fH2Pi0DigitId->Fill(digit1.GetDigitId(),mass);
	    fH2Pi0DigitId->Fill(digit2.GetDigitId(),mass);
	  }
	} // End of energy cut
      } // End of second cluster loop (clt2)
    } // End of energy cut
  } // End of first cluster loop (clt1)
  fH1DigitNum->Fill(ndigitsall);
  //
}
//-------------------------------------------------------------------
void AliPHOSDApi0mip::Print(Option_t *option) const
{
  // Print Out

  //fTFile->ls();
  //fTTree->Print();
  if( fEvent ){
    //fEvent->ExecuteClustering();
    fEvent->Print(option);
  }
}
//-------------------------------------------------------------------
