/*************************************************************************
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
//
// Class for impact parameter (DCA) of charged particles
// Study resolution and pull: prepare for beauty study
//
// Authors:
//   Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//   Carlo Bombonati <carlo.bombonati@cern.ch>
//

#include "TMath.h"
#include "TH1F.h"
#include "TList.h"
#include <TParticle.h>
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"

#include "AliHFEdca.h"

ClassImp(AliHFEdca)

//________________________________________________________________________
const Char_t* AliHFEdca::fgkParticles[12] = {
  // particles name
  "electron", "muonMinus","pionMinus", "kaonMinus", "protonMinus", 
  "positron", "muonPlus", "pionPlus", "kaonPlus", "protonPlus",
  "allNegative", "allPositive"
};
//________________________________________________________________________
const Int_t AliHFEdca::fgkColorPart[12] = { 
  // colors assigned to particles
  kRed, kBlue, kGreen+2, kYellow+2, kMagenta, 
  kRed+2, kBlue+2, kGreen+4, kYellow+4, kMagenta+2,
  kBlack, kGray+1
};



//________________________________________________________________________
const Float_t AliHFEdca::fgkPtIntv[44] = {
  // define pT bins
  0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
  1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 
  3.0, 3.3, 3.6, 3.9, 4.2, 4.5, 5.0, 5.5, 6.0, 6.5, 
  7.0, 8.0, 9.0, 10., 11., 12., 13., 14., 15., 16., 
  17., 18., 19., 20.};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkDcaVar[2] = {
  "deltaDcaXY",  "deltaDcaZ"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkDcaVarTitle[2] ={
  ";residual #Delta(d_{xy}) [#mum];couts", ";residual #Delta(d_{z}) [#mum];counts"};


//________________________________________________________________________
const Char_t *AliHFEdca::fgkPullDcaVar[2] = {
  "pullDcaXY", "pullDcaZ"
};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkPullDcaVarTitle[2] = {
  ";(dca_{xy}^{ESD}-dca_{xy}^{MC})/(#sigma_{track}#oplus#sigma_{Vxy});counts",
  ";(dca_{z}^{ESD}-dca_{z}^{MC})/(#sigma_{track}#oplus#sigma_{Vz}); counts"
};

//________________________________________________________________________
AliHFEdca::AliHFEdca():
   fResidualList(0x0)
  , fPullList(0x0)

{
  // default constructor
 
}

//________________________________________________________________________
AliHFEdca::AliHFEdca(const AliHFEdca &dca):
  TObject(dca)
  , fResidualList(0x0)
  , fPullList(0x0)

{
  // copy constructor

}
//_______________________________________________________________________________________________
AliHFEdca&
AliHFEdca::operator=(const AliHFEdca &)
{
  //
  // Assignment operator
  //
  
  Printf("Not yet implemented.");
  return *this;
}

//________________________________________________________________________
AliHFEdca::~AliHFEdca()
{
  // default destructor

  for(Int_t j=0; j<kNParticles; j++){
    for(Int_t i=0; i<kNPtBins; i++){
      if(fHistDcaXYRes[j][i]) delete fHistDcaXYRes[j][i];
      if(fHistDcaZRes[j][i]) delete fHistDcaZRes[j][i];
      if(fHistDcaXYPull[j][i]) delete fHistDcaXYPull[j][i];
      if(fHistDcaXYPull[j][i]) delete fHistDcaZPull[j][i];
    }
  }

  if(fResidualList) delete fResidualList;
  if(fPullList) delete fPullList;

  Printf("analysis done\n");
}

//________________________________________________________________________
void AliHFEdca::InitAnalysis(){
  
  Printf("initialize analysis\n");

}


//________________________________________________________________________
void AliHFEdca::PostAnalysis() const
{
  // do fit
  // moved to dcaPostAnalysis.C

}
//________________________________________________________________________
void AliHFEdca::CreateHistogramsResidual(TList *residualList){
  // define histogram
  // 1. residual

  // for residuals
  fHistDcaXYRes[kNParticles][kNPtBins]=0x0;
  fHistDcaZRes[kNParticles][kNPtBins]=0x0;
  
  const Int_t nBins = 1000;
  const Float_t maxXYBin = 1000.;
  const Float_t maxZBin = 1000.;

  
  for(Int_t k=0; k<kNDcaVar; k++){
    TString histTitle((const char*)fgkDcaVarTitle[k]);
    
    for(Int_t j=0; j<kNParticles; j++){
      for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);
	
	histName += Form("_%s_pT-%.1f-%.1f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistDcaXYRes[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	  fHistDcaXYRes[j][i]->SetLineColor((const int)fgkColorPart[j]);
	}	    
	if(k==1){
	  fHistDcaZRes[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	  fHistDcaZRes[j][i]->SetLineColor((const int)fgkColorPart[j]);
	}   
      } // 43 pt bins
    } //12 nparticles
  } // 2 dca var
  
  //  TList *fResidualList = 0;
  residualList->SetOwner();
  residualList->SetName("residual");
  for(Int_t iPart=0; iPart<kNParticles; iPart++){
    for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
      residualList->Add(fHistDcaXYRes[iPart][iPtBin]);  
      residualList->Add(fHistDcaZRes[iPart][iPtBin]);  
    } // loop over pt bins
  }  // loop over particles (pos, neg)
  

  

}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsPull(TList *pullList){
  // define histogram
  // 2. pull

  const Int_t nBins = 1000;
  const Float_t maxXYBin = 20.;
  const Float_t maxZBin = 20.;

  
  // for pull -----------------------------------------------------------------------
  fHistDcaXYPull[kNParticles][kNPtBins]=0x0;
  fHistDcaZPull[kNParticles][kNPtBins]=0x0;

  
  for(Int_t k=0; k<kNDcaVar; k++){
    TString histTitle((const char*)fgkPullDcaVarTitle[k]);
    
    for(Int_t j=0; j<kNParticles; j++){
      for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);
	
	histName += Form("_%s_pT-%.1f-%.1f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistDcaXYPull[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, 1-maxXYBin, 1+maxXYBin);
	  fHistDcaXYPull[j][i]->SetLineColor((const int)fgkColorPart[j]);
	}	    
	if(k==1){
	  fHistDcaZPull[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, 1-maxZBin, 1+maxZBin);
	  fHistDcaZPull[j][i]->SetLineColor((const int)fgkColorPart[j]);
	}   
      } // 43 pt bins
    } //6 nparticles
  } // 2 dca var
  
  //  TList *fPullList = 0;
  pullList->SetOwner();
  pullList->SetName("pull");
  for(Int_t iPart=0; iPart<kNParticles; iPart++){
    for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
      pullList->Add(fHistDcaXYPull[iPart][iPtBin]);  
      pullList->Add(fHistDcaZPull[iPart][iPtBin]);  
    } // loop over pt bins
  }  // loop over particles (pos, neg)
  
  

}


//_______________________________________________________________________________________________
void AliHFEdca::FillHistograms(AliESDEvent * const esdEvent, AliESDtrack * const track, AliMCEvent * const mcEvent)
{
// filling historgams track by track

// obtaining reconstructed dca ------------------------------------------------------------------
  Float_t esdpx = track->Px();
  Float_t esdpy = track->Py();
  Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  
  Float_t b[2];  // dca in cm
  Float_t bCov[3];  // covariance matrix
  track->GetImpactParameters(b,bCov);
  
// obtaining errors of dca ------------------------------------------------------------------
  const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
  Float_t magneticField = 5;  // initialized as 5kG
  magneticField = esdEvent->GetMagneticField();  // in kG
 
  Double_t dz[2];   // error of dca in cm
  Double_t covardz[3];
  track->PropagateToDCA(primVtx,magneticField, 1000., dz, covardz);

  // calculate mcDca ------------------------------------------------------------------
  
  AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));  
  TParticle *part = mctrack->Particle();
  Int_t pdg = part->GetPdgCode();
  Int_t charge = 1;
  if(pdg==kPDGelectron || pdg==kPDGmuon 
     || pdg==-kPDGpion || pdg==-kPDGkaon || pdg==-kPDGproton) charge = -1;

  Float_t vx = part->Vx();  // in cm
  Float_t vy = part->Vy();  // in cm
  Float_t vz = part->Vz();   // in cm
  
  Float_t vxy = TMath::Sqrt(vx*vx+vy*vy);
  
  Float_t mcpx = part->Px();
  Float_t mcpy = part->Py();
  Float_t mcpt = TMath::Sqrt(mcpx*mcpx+mcpy*mcpy);
  
  const Float_t conv[2] = {1.783/1.6, 2.99792458};
  Float_t radiusMc = mcpt/(TMath::Abs(magneticField)/10.)*conv[0]*conv[1]; // pt in GeV/c, magnetic field in Tesla
  
  Float_t nx = esdpx/mcpt;
  Float_t ny = esdpy/mcpt;
    
  Float_t radius;
  radius = TMath::Abs(radiusMc);
  Float_t mcDcaXY = (radius - TMath::Sqrt(vxy*vxy/100./100. + radius*radius + 2*radius*charge*(vx*ny-vy*nx)/100.)) ;  // in meters

//   printf("magnetic Field = %.3f \t", magneticField);
//   printf("pt=esd  %.3f/mc  %.3f GeV/c, radius=esd  %.3f/mc  %.3f meter \n", esdpt, mcpt, radiusEsd, radiusMc);
//   printf("mcDcaXY=%.5f micron, esdDcaXY=%.5f micron \t ",  mcDcaXY*1.e6, b[0]*1e4);
//   printf("mcDcaZ=%.5f micron, esdDcaZ=%.5f micron\n\n", vz*1e6, b[1]*1e4);
  
  Double_t mcDca[2] = {mcDcaXY*100, vz};  // in cm
  
  Double_t residual[2] = {0, 0};
  Double_t pull[2] = {0, 0};
  Double_t error[2] ={TMath::Sqrt(covardz[0]), TMath::Sqrt(covardz[2])};
  for(Int_t i=0; i<2; i++){
    residual[i] = dz[i] - mcDca[i]; // in centimeters       
    if(error[i]!=0)pull[i] = residual[i]/error[i];   // unitless
    //    printf("error[%d]=%.6f  residual[%d]=%.6f   pull[%d]=%.6f\n", i, error[i], i, residual[i], i, pull[i]);
  }

  Int_t fPdgParticle[10] = { 
    kPDGelectron, kPDGmuon, -kPDGpion, -kPDGkaon, -kPDGproton, 
    -kPDGelectron, -kPDGmuon, kPDGpion, kPDGkaon, kPDGproton};
    
  for(Int_t iPart=0; iPart<kNParticles-2; iPart++){
    
    // identified ones
    for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
      if(pdg==fPdgParticle[iPart] && (esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1])) {
	fHistDcaXYRes[iPart][iPtBin]->Fill(residual[0]*1.0e4);  // in microns
	fHistDcaZRes[iPart][iPtBin]->Fill(residual[1]*1.0e4);   // in microns
	fHistDcaXYPull[iPart][iPtBin]->Fill(pull[0]);
	fHistDcaZPull[iPart][iPtBin]->Fill(pull[1]);
      }
      else
 	continue;
    }
  }
    
  // for charged particles
  for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
    if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
      Int_t iPart = 10;
      if(charge>0) iPart = 11;
      fHistDcaXYRes[iPart][iPtBin]->Fill(residual[0]*1e4);
      fHistDcaZRes[iPart][iPtBin]->Fill(residual[1]*1e4);
      fHistDcaXYPull[iPart][iPtBin]->Fill(pull[0]);
      fHistDcaZPull[iPart][iPtBin]->Fill(pull[1]);
    }
    else
      continue;
  } 
    
}

