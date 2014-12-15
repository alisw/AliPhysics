/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the base class for dEdx analysis            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TFile.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliITSdEdxAnalyzer.h"
#include "AliExternalTrackParam.h"
//#include "AliITSpidESD.h"
#include "AliITSPIDResponse.h"
#include "AliPID.h"

ClassImp(AliITSdEdxAnalyzer)

const Int_t AliITSdEdxAnalyzer::fgkPdgCode[kNParticles]={211,321,2212};
const Int_t AliITSdEdxAnalyzer::fgkLayerCode[kNValuesdEdx]={3,4,5,6,0};

//______________________________________________________________________
AliITSdEdxAnalyzer::AliITSdEdxAnalyzer() :  
  TObject(),
  fNPBins(10),
  fPMin(0.1),
  fPMax(1.1),
  fHistodEdx(0),
  fHistoDeltadEdx(0),
  fHistodEdxVsP(0),
  fThickness(0.03),
  fDensity(2.33),
  fBBmodel(0),
  fMIP(82.),
  fTPCpidCut(0.5)
{
  // default constructor
  BookHistos();
}
//______________________________________________________________________
AliITSdEdxAnalyzer::AliITSdEdxAnalyzer(const Int_t npbins, const Float_t pmin, const Float_t pmax) :  
  TObject(),
  fNPBins(npbins),
  fPMin(pmin),
  fPMax(pmax),
  fHistodEdx(0),
  fHistoDeltadEdx(0),
  fHistodEdxVsP(0),
  fThickness(0.03),
  fDensity(2.33),
  fBBmodel(0),
  fMIP(82.),
  fTPCpidCut(0.5)
{
  // standard constructor
  BookHistos();
}
//______________________________________________________________________
AliITSdEdxAnalyzer::~AliITSdEdxAnalyzer(){
  // destructor
  DeleteHistos();
}
//______________________________________________________________________
void AliITSdEdxAnalyzer::SetMomentumBins(const Int_t npbins, const Float_t pmin, const Float_t pmax){
  // Kill exisiting histos, set new binning, book new histos
  DeleteHistos();
  fNPBins=npbins;
  fPMin=pmin;
  fPMax=pmax;
  BookHistos();
}
//______________________________________________________________________
void AliITSdEdxAnalyzer::DeleteHistos(){
  // deletes the hitograms
  if(fHistodEdx){
    for(Int_t i=0; i<GetNHistos();i++) delete fHistodEdx[i];
    delete [] fHistodEdx;
    fHistodEdx=0;
  }
  if(fHistoDeltadEdx){
    for(Int_t i=0; i<GetNHistos();i++) delete fHistoDeltadEdx[i];
    delete [] fHistoDeltadEdx;
    fHistoDeltadEdx=0;
  }
  if(fHistodEdxVsP){
    for(Int_t i=0; i<GetNHistos2();i++) delete fHistodEdxVsP[i];
    delete [] fHistodEdxVsP;
    fHistodEdxVsP=0;
  }
}
//______________________________________________________________________
void AliITSdEdxAnalyzer::BookHistos(){
  // defines the output histograms 
  fHistodEdx=new TH1F*[GetNHistos()];
  fHistoDeltadEdx=new TH1F*[GetNHistos()];
  for(Int_t iSpecie=0; iSpecie<kNParticles; iSpecie++){
    for(Int_t iPBin=0; iPBin<fNPBins; iPBin++){
      for(Int_t iVal=0; iVal<kNValuesdEdx; iVal++){
	TString hisnam1=Form("hdEdx%dpbin%dpart%d",fgkLayerCode[iVal],iPBin,fgkPdgCode[iSpecie]);
	TString histit1=Form("dEdx layer %d (keV) pbin%d part%d",fgkLayerCode[iVal],iPBin,fgkPdgCode[iSpecie]);
	if(iVal==kNValuesdEdx-1) histit1=Form("dEdx trunc. mean (keV) pbin%d part%d",iPBin,fgkPdgCode[iSpecie]);
	fHistodEdx[GetIndex(iSpecie,iPBin,iVal)]=new TH1F(hisnam1.Data(),histit1.Data(),100,0.,500.);
	fHistodEdx[GetIndex(iSpecie,iPBin,iVal)]->GetXaxis()->SetTitle(histit1.Data());
	TString hisnam2=Form("hdeltadEdx%dpbin%dpart%d",fgkLayerCode[iVal],iPBin,fgkPdgCode[iSpecie]);
	TString histit2=Form("(dEdx_l%d-BetheBloch)/BetheBloch pbin%d part%d",fgkLayerCode[iVal],iPBin,fgkPdgCode[iSpecie]);
	if(iVal==kNValuesdEdx-1) histit1=Form("(dEdx_truncmean-BetheBloch)/BetheBloch pbin%d part%d",iPBin,fgkPdgCode[iSpecie]);
	fHistoDeltadEdx[GetIndex(iSpecie,iPBin,iVal)]=new TH1F(hisnam2.Data(),histit2.Data(),100,-0.5,0.5);
	fHistoDeltadEdx[GetIndex(iSpecie,iPBin,iVal)]->GetXaxis()->SetTitle(histit2.Data());
      }
    }
  }

  fHistodEdxVsP=new TH2F*[GetNHistos2()];
  for(Int_t iSpecie=0; iSpecie<kNParticles; iSpecie++){
    for(Int_t iVal=0; iVal<kNValuesdEdx; iVal++){
      TString hisnam=Form("hdEdx%dVsPpart%d",fgkLayerCode[iVal],fgkPdgCode[iSpecie]);
      TString histit=Form("dEdx layer %d (keV) vs P part%d",fgkLayerCode[iVal],fgkPdgCode[iSpecie]);
      if(iVal==kNValuesdEdx-1) histit=Form("dEdx trunc. mean (keV) vs P part%d",fgkPdgCode[iSpecie]);      
      fHistodEdxVsP[GetIndex2(iSpecie,iVal)]=new TH2F(hisnam.Data(),histit.Data(),50,fPMin,fPMax,50,0.,500.);
      histit.ReplaceAll(" vs P "," ");
      fHistodEdxVsP[GetIndex2(iSpecie,iVal)]->GetYaxis()->SetTitle(histit.Data());
      fHistodEdxVsP[GetIndex2(iSpecie,iVal)]->GetXaxis()->SetTitle("Momentum (GeV/c)");
    }
  }
}
//______________________________________________________________________
void AliITSdEdxAnalyzer:: WriteHistos(TString filename) const {
  // stores output histogrma to file
  TFile *outfile=new TFile(filename.Data(),"recreate");
  outfile->cd();
  for(Int_t i=0; i<GetNHistos();i++){ 
    fHistodEdx[i]->Write();
    fHistoDeltadEdx[i]->Write();
  }
  for(Int_t i=0; i<GetNHistos2();i++) fHistodEdxVsP[i]->Write();
  outfile->Close();
  delete outfile;
}
//______________________________________________________________________
void AliITSdEdxAnalyzer::ReadEvent(const AliESDEvent* esd, AliStack* stack){
  // Fill histos 
  // if stack is !=0 MC truth is used to define particle specie
  // else TPC pid is used to define particle specie

  if(!esd) AliFatal("Bad ESD event");
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = esd->GetTrack(iTrack);
    Double_t trmean=track->GetITSsignal();
    Double_t sig[4];
    track->GetITSdEdxSamples(sig);
    Double_t preco=track->P();
    Int_t iPBin=GetMomentumBin(preco);
    if(iPBin==-1) continue;
    Int_t iSpecie=-1;
    Double_t dedxbb=0;
    if(stack){
      Int_t lab=track->GetLabel();
      if(lab<0) continue;
      TParticle* part=(TParticle*)stack->Particle(lab);
      Int_t absPdgCode=TMath::Abs(part->GetPdgCode());
      iSpecie=GetSpecieBin(absPdgCode);
      if(iSpecie==-1) continue;
      dedxbb=BetheBloch(part);
    }else{   
      iSpecie=GetPaticleIdFromTPC(track);
      if(iSpecie==-1) continue;
      TParticlePDG* p=TDatabasePDG::Instance()->GetParticle(fgkPdgCode[iSpecie]);    
      dedxbb=BetheBloch(preco,p->Mass());
    }
    for(Int_t ilay=0;ilay<4;ilay++){    
      if(sig[ilay]>0.){
	fHistodEdx[GetIndex(iSpecie,iPBin,ilay)]->Fill(sig[ilay]);
	fHistoDeltadEdx[GetIndex(iSpecie,iPBin,ilay)]->Fill((sig[ilay]-dedxbb)/dedxbb);
	fHistodEdxVsP[GetIndex2(iSpecie,ilay)]->Fill(preco,sig[ilay]);
      }
    }
    if(trmean>0.){
      fHistodEdx[GetIndex(iSpecie,iPBin,4)]->Fill(trmean);
      fHistoDeltadEdx[GetIndex(iSpecie,iPBin,4)]->Fill((trmean-dedxbb)/dedxbb);
      fHistodEdxVsP[GetIndex2(iSpecie,4)]->Fill(preco,trmean);
    }
  }
}
//______________________________________________________________________
Int_t AliITSdEdxAnalyzer::GetPaticleIdFromTPC(const AliESDtrack* track) const {
  // Returns the particle specie as identified by TPC
  Double_t tpcpid[AliPID::kSPECIES];
  track->GetTPCpid(tpcpid);
  Int_t maxPos=-1;
  Float_t maxProb=0.;
  for(Int_t is=0; is<AliPID::kSPECIES; is++){
    if(tpcpid[is]>maxProb){
      maxProb=tpcpid[is];
      maxPos=is;
    }
  }
  Int_t iSpecie=-1;
  if(maxProb>fTPCpidCut){
    if(maxPos==AliPID::kPion) iSpecie=0;
    else if(maxPos==AliPID::kKaon) iSpecie=1;
    else if(maxPos==AliPID::kProton) iSpecie=2;    
  }
  return iSpecie;
}
//______________________________________________________________________
Double_t AliITSdEdxAnalyzer::BetheBloch(const Float_t p, const Float_t m) const {
  // Computes expected dE/dx value for given particle specie and momentum
  static AliITSPIDResponse pidResponse;
  Double_t dedxbb=0.;
  if(fBBmodel==0){
    Double_t betagamma=p/m;
    Double_t conv=fDensity*1E6*fThickness/116.24*fMIP;
    dedxbb=conv*AliExternalTrackParam::BetheBlochSolid(betagamma);
  }else if(fBBmodel==1){
    dedxbb=pidResponse.Bethe(p,m);
  }
  return dedxbb;
}
//______________________________________________________________________
TGraph* AliITSdEdxAnalyzer::GetBetheBlochGraph(const Int_t pdgcode) const {
 // Fills a TGraph with Bethe-Bloch expe
 TGraph* g=new TGraph(0);
 TParticlePDG* part=TDatabasePDG::Instance()->GetParticle(pdgcode);
 Float_t mass=part->Mass();
 for(Int_t ip=0; ip<100;ip++){
   Float_t p=fPMin+(ip+0.5)*(fPMax-fPMin)/100.;
   g->SetPoint(ip,p,BetheBloch(p,mass));
 }
 g->GetXaxis()->SetTitle("Momentum (GeV/c)");
 g->GetYaxis()->SetTitle("dE/dx (keV/300 #mum)");
 return g;
}
