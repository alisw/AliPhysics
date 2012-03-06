#include "AliHFEVZEROEventPlane.h"

// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TString.h"
#include "TVector2.h"
#include "AliEventplane.h"


// STL includes
#include <iostream>
using namespace std;


//_____________________________________________________________________________
AliHFEVZEROEventPlane::AliHFEVZEROEventPlane():
  TNamed(),
  fEventPlaneV0A(-100.0),
  fEventPlaneV0C(-100.0),
  fEventPlaneV0(-100.0),
  fESD(0),
  fRun(-1),
  fMultV0(0x0),
  fV0Cpol(100),
  fV0Apol(100),
  fnamefile(TString("")),
  fOutputList(0x0),
  fMultV0Before(0x0),
  fMultV0After(0x0)
{
  //
  // Default constructor (should not be used)
  //

  for(Int_t k = 0; k < fgknCentrBin; k++) {
    for(Int_t iside = 0; iside < 2; iside++) {
      for(Int_t icoord = 0; icoord < 2; icoord++) {
	fQBefore[k][iside][icoord] = 0x0;
	fQAfter[k][iside][icoord] = 0x0;
      }
    }
  }
  
}

//______________________________________________________________________________
AliHFEVZEROEventPlane::AliHFEVZEROEventPlane(const char *name, const Char_t *title):
  TNamed(name,title),
  fEventPlaneV0A(-100.0),
  fEventPlaneV0C(-100.0),
  fEventPlaneV0(-100.0),
  fESD(0),
  fRun(-1),
  fMultV0(0x0),
  fV0Cpol(100),
  fV0Apol(100),
  fnamefile(TString("")),
  fOutputList(0x0),
  fMultV0Before(0x0),
  fMultV0After(0x0)
{
  //
  // Constructor
  //

  char namelist[100];
  sprintf(namelist,"QA_hist_%s",name);
  fOutputList = new TList();
  fOutputList->SetName((const char*)namelist);
  fOutputList->SetOwner(kTRUE);

  // Multiplicity
  fMultV0Before = new TProfile("MultV0Before","",64,0,64);
  fMultV0Before->Sumw2();
  fMultV0After = new TProfile("MultV0After","",64,0,64);
  fMultV0After->Sumw2();
  fOutputList->Add(fMultV0Before);
  fOutputList->Add(fMultV0After);

  // Recentering
  char namecontbefore[100];
  char namecontafter[100];    
  for(Int_t k = 0; k < fgknCentrBin; k++) {
    for(Int_t iside = 0; iside < 2; iside++) {
      for(Int_t icoord = 0; icoord < 2; icoord++) {
	
	if(iside==0 && icoord==0)
	  sprintf(namecontbefore,"hQxc2_%i",k);
	else if(iside==1 && icoord==0)
	  sprintf(namecontbefore,"hQxa2_%i",k);
	else if(iside==0 && icoord==1)
	  sprintf(namecontbefore,"hQyc2_%i",k);
	else if(iside==1 && icoord==1)
	  sprintf(namecontbefore,"hQya2_%i",k);
	//
	if(iside==0 && icoord==0)
	  sprintf(namecontafter,"hQxc2_%i_after",k);
	else if(iside==1 && icoord==0)
	  sprintf(namecontafter,"hQxa2_%i_after",k);
	else if(iside==0 && icoord==1)
	  sprintf(namecontafter,"hQyc2_%i_after",k);
	else if(iside==1 && icoord==1)
	  sprintf(namecontafter,"hQya2_%i_after",k);
	
	//
	fQBefore[k][iside][icoord] = new TH1F(((const char*)namecontbefore),"",800,-400.0,400.0);
	fQBefore[k][iside][icoord]->Sumw2();
	fOutputList->Add(fQBefore[k][iside][icoord]);
	fQAfter[k][iside][icoord] = new TH1F(((const char*)namecontafter),"",800,-400.0,400.0);
	fQAfter[k][iside][icoord]->Sumw2();
	fOutputList->Add(fQAfter[k][iside][icoord]);
	//
      }
    }
  }
  
}
//____________________________________________________________
AliHFEVZEROEventPlane::AliHFEVZEROEventPlane(const AliHFEVZEROEventPlane &ref):
  TNamed(ref),
  fEventPlaneV0A(ref.fEventPlaneV0A),
  fEventPlaneV0C(ref.fEventPlaneV0C),
  fEventPlaneV0(ref.fEventPlaneV0),
  fESD(0),
  fRun(-1),
  fMultV0(0x0),
  fV0Cpol(100),
  fV0Apol(100),
  fnamefile(ref.fnamefile),
  fOutputList(0x0),
  fMultV0Before(0x0),
  fMultV0After(0x0)
{
  //
  // Copy Constructor
  //
  ref.Copy(*this);
}

//____________________________________________________________
AliHFEVZEROEventPlane &AliHFEVZEROEventPlane::operator=(const AliHFEVZEROEventPlane &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}

//____________________________________________________________
void AliHFEVZEROEventPlane::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliHFEVZEROEventPlane &target = dynamic_cast<AliHFEVZEROEventPlane &>(o);

  target.fEventPlaneV0A = fEventPlaneV0A;
  target.fEventPlaneV0C = fEventPlaneV0C;
  target.fEventPlaneV0 = fEventPlaneV0;
  target.fnamefile = fnamefile;
  
}

//__________________________________________________________________
AliHFEVZEROEventPlane::~AliHFEVZEROEventPlane(){
  //
  // Destruktor
  //
  if(fOutputList){
    fOutputList->Delete();
    delete fOutputList;
  }
  fOutputList = 0x0;
  if(fMultV0Before) delete fMultV0Before;
  if(fMultV0After) delete fMultV0After;

  for(Int_t k = 0; k < fgknCentrBin; k++) {
    for(Int_t iside = 0; iside < 2; iside++) {
      for(Int_t icoord = 0; icoord < 2; icoord++) {
	if(fQBefore[k][iside][icoord]) delete fQBefore[k][iside][icoord];
	if(fQAfter[k][iside][icoord]) delete fQAfter[k][iside][icoord];
      }
    }
  }

}
//______________________________________________________________________________
void AliHFEVZEROEventPlane::ProcessEvent(AliESDEvent *event) 
{
  //
  // Process the event
  // 
  
  // Reset
  fEventPlaneV0A = -100.0;
  fEventPlaneV0C = -100.0;
  fEventPlaneV0  = -100.0;
  
  //
  Int_t run = event->GetRunNumber();
  Bool_t ok = kTRUE;
  if(run != fRun){
    if(fMultV0) delete fMultV0;
    fMultV0 = 0x0;
    // Load the calibrations run dependent
    if(!OpenInfoCalbration(run)) ok = kFALSE;
    fRun=run;
  }

  if(ok) Analyze(event);

}

//________________________________________________________________________
void AliHFEVZEROEventPlane::Analyze(AliESDEvent* esdEvent)
{  
  //
  // Do VZERO calibration + centering
  //

  if(!fMultV0) {
    //printf("Did not find calibration VZERO\n");
    return;
  }
  //printf("HERE!!!\n");

  //Centrality
  Float_t v0Centr  = -10.;
  AliCentrality* centrality = esdEvent->GetCentrality();
  if (centrality){
    v0Centr  = centrality->GetCentralityPercentile("V0M");
  }
  
  // Analyse only for 0-80% PbPb collisions
  Int_t iC = -1;    
  if (v0Centr >0 && v0Centr < 80){
    if(v0Centr < 5) iC = 0;
    else if(v0Centr < 10) iC = 1;
    else if(v0Centr < 20) iC = 2;
    else if(v0Centr < 30) iC = 3;
    else if(v0Centr < 40) iC = 4;
    else if(v0Centr < 50) iC = 5;
    else if(v0Centr < 60) iC = 6;
    else if(v0Centr < 70) iC = 7;
    else iC = 8;
    
    //general	
    Double_t qxa2 = 0, qya2 = 0;
    Double_t qxc2 = 0, qyc2 = 0;
    
    //V0 info    
    AliESDVZERO* esdV0 = esdEvent->GetVZEROData();
    
    for (Int_t iv0 = 0; iv0 < 64; iv0++) {
      Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
      Float_t multv0 = esdV0->GetMultiplicity(iv0);
      if (iv0 < 32){
	//printf("Value %f\n",fMultV0->GetBinContent(iv0+1));
	if(fMultV0->GetBinContent(iv0+1)>0.0) {
	  qxc2 += TMath::Cos(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	  qyc2 += TMath::Sin(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	  fMultV0After->Fill(iv0,multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1));
	}
      } else {
	if(fMultV0->GetBinContent(iv0+1)>0.0) {
	  qxa2 += TMath::Cos(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	  qya2 += TMath::Sin(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	  fMultV0After->Fill(iv0,multv0*fV0Apol/fMultV0->GetBinContent(iv0+1));
	}
      }
      fMultV0Before->Fill(iv0,multv0);
    }

    // Fill histos
    fQBefore[iC][1][0] -> Fill(qxa2);
    fQBefore[iC][1][1] -> Fill(qya2);
    fQBefore[iC][0][0] -> Fill(qxc2);
    fQBefore[iC][0][1] -> Fill(qyc2);
    
    //grab for each centrality the proper histo with the Qx and Qy to do the recentering
    Double_t qxamean2 = fMeanQ[iC][1][0];
    Double_t qxarms2  = fWidthQ[iC][1][0];
    Double_t qyamean2 = fMeanQ[iC][1][1];
    Double_t qyarms2  = fWidthQ[iC][1][1];
    
    Double_t qxcmean2 = fMeanQ[iC][0][0];
    Double_t qxcrms2  = fWidthQ[iC][0][0];
    Double_t qycmean2 = fMeanQ[iC][0][1];
    Double_t qycrms2  = fWidthQ[iC][0][1];	

    if((TMath::Abs(qxarms2) < 0.00000001) || (TMath::Abs(qyarms2) < 0.00000001) || (TMath::Abs(qxcrms2) < 0.00000001) || (TMath::Abs(qycrms2) < 0.00000001)) return;    

    Double_t qxaCor2 = (qxa2 - qxamean2)/qxarms2;
    Double_t qyaCor2 = (qya2 - qyamean2)/qyarms2;
    Double_t qxcCor2 = (qxc2 - qxcmean2)/qxcrms2;
    Double_t qycCor2 = (qyc2 - qycmean2)/qycrms2;

    Double_t qxCor2 = qxaCor2 + qxcCor2;
    Double_t qyCor2 = qyaCor2 + qycCor2;

    // Fill histos
    fQAfter[iC][1][0] -> Fill(qxaCor2);
    fQAfter[iC][1][1] -> Fill(qyaCor2);
    fQAfter[iC][0][0] -> Fill(qxcCor2);
    fQAfter[iC][0][1] -> Fill(qycCor2);
    
    Double_t evPlAngV0ACor2 = TVector2::Phi_0_2pi(TMath::ATan2(qyaCor2, qxaCor2))/2.;
    Double_t evPlAngV0CCor2 = TVector2::Phi_0_2pi(TMath::ATan2(qycCor2, qxcCor2))/2.;
    Double_t evPlAngV0Cor2  = TVector2::Phi_0_2pi(TMath::ATan2(qyCor2, qxCor2))/2.;

    fEventPlaneV0A = evPlAngV0ACor2;
    fEventPlaneV0C = evPlAngV0CCor2;
    fEventPlaneV0  = evPlAngV0Cor2;
    
    //printf("Eventplane V0A %f, V0C %f\n",fEventPlaneV0A,fEventPlaneV0C);

  }
}
//_____________________________________________________________________________
Bool_t AliHFEVZEROEventPlane::OpenInfoCalbration(Int_t run)
{
  //
  // Take the calibration coefficients
  //  

  //printf("Name of the file %s\n",(const char*)fnamefile);
 
  TFile *foadb = TFile::Open(fnamefile.Data());
  if(!foadb){
    printf("OADB file %s cannot be opened\n",fnamefile.Data());
    return kFALSE;
  }

  //printf("test\n");

  AliOADBContainer *cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorr");
  if(!cont){
    printf("OADB object hMultV0BefCorr is not available in the file\n");
    return kFALSE;	
  }

  if(!(cont->GetObject(run))){
    printf("OADB object hMultV0BefCorr is not available for run %i\n",run);
    return kFALSE;	
  }
  //TProfile *multV0 = ((TH2F *) cont->GetObject(run))->ProfileX();
  TProfile *multV0 = ((TProfile *) cont->GetObject(run));
  if(!multV0) return kFALSE;
  fMultV0 = (TProfile *) multV0->Clone();
  fMultV0->SetDirectory(0);
  
  TF1 *fpol0 = new TF1("fpol0","pol0"); 
  fMultV0->Fit(fpol0,"Q0","",0,31);
  fV0Cpol = fpol0->GetParameter(0);
  fMultV0->Fit(fpol0,"Q0","",32,64);
  fV0Apol = fpol0->GetParameter(0);
  
  for(Int_t iside=0;iside<2;iside++){
    for(Int_t icoord=0;icoord<2;icoord++){
      for(Int_t i=0;i  < fgknCentrBin;i++){
	char namecont[100];
	if(iside==0 && icoord==0)
	  sprintf(namecont,"hQxc2_%i",i);
	else if(iside==1 && icoord==0)
	  sprintf(namecont,"hQxa2_%i",i);
	else if(iside==0 && icoord==1)
	  sprintf(namecont,"hQyc2_%i",i);
	else if(iside==1 && icoord==1)
	  sprintf(namecont,"hQya2_%i",i);
	
	cont = (AliOADBContainer*) foadb->Get(namecont);
	if(!cont){
	  printf("OADB object %s is not available in the file\n",namecont);
	  delete fpol0;
	  return kFALSE;	
	}
	
	if(!(cont->GetObject(run))){
	  printf("OADB object %s is not available for run %i\n",namecont,run);
	  delete fpol0;
	  return kFALSE;	
	}
	fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
	fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();
      }
    }
  }
  
  delete fpol0;
  foadb->Close();
  
  return kTRUE;

}
