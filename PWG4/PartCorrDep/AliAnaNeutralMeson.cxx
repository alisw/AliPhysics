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
/* $Id: $ */

//_________________________________________________________________________
// a general class to fill two-photon and pi0+gamma invariant mass hisograms 
// to be used to extract pi0, eta and omega raw yield.
// also for PHOS acceptance
//--
//-- by Renzhuo Wan (Iopp-wuhan) May 13,2009
//_________________________________________________________________________

// --- ROOT system ---
#include "TH3F.h"
#include "TH2F.h"
//#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TParticle.h"
//---- AliRoot system ----
#include "AliAnaNeutralMeson.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliStack.h"
#include "AliFidutialCut.h"
#include "AliVEvent.h"
#ifdef __PHOSGEO__
   #include "AliPHOSGeoUtils.h"
#endif

ClassImp(AliAnaNeutralMeson)

//______________________________________________________________________________
AliAnaNeutralMeson::AliAnaNeutralMeson() : AliAnaPartCorrBaseClass(),
fAnaPi0Eta(0), fAnaOmega(0), fNPID(0),fNmaxMixEv(0), fNAsy(0),fNMod(0), 
fNHistos(0), fPerPtBin(0),
fAsyCut(0), fModCut(0),
fDist(0), fAsy(0), fPt(0), fMass(0),
fInvMassCut(0),
fNbinsAsy(0), fMinAsy(0), fMaxAsy(0),
fNbinsPt(0), fMinPt(0), fMaxPt(0),
fNbinsM(0), fMinM(0), fMaxM(0),
fEventsListPi0Eta(0), fEventsListOmega(0),
fPi0Mass(0), fPi0MassPeakWidthCut(0),  fhNClusters(0),
fhRecPhoton(0), fhRecPhotonEtaPhi(0),
fReal2Gamma(0), fMix2Gamma(0),fRealOmega(0),
fMixOmegaA(0),fMixOmegaB(0), fMixOmegaC(0),
fRealTwoGammaAsyPtM(0), fMixTwoGammaAsyPtM(0),
fRealPi0GammaAsyPtM(0), fMixAPi0GammaAsyPtM(0), fMixBPi0GammaAsyPtM(0), fMixCPi0GammaAsyPtM(0),
fhPrimPhotonPt(0), fhPrimPhotonAccPt(0), fhPrimPhotonY(0), fhPrimPhotonAccY(0), 
fhPrimPhotonPhi(0), fhPrimPhotonAccPhi(0),
fhPrimPi0Pt(0), fhPrimPi0AccPt(0), fhPrimPi0Y(0), fhPrimPi0AccY(0), fhPrimPi0Phi(0), fhPrimPi0AccPhi(0), 
fhPrimEtaPt(0), fhPrimEtaAccPt(0), fhPrimEtaY(0), fhPrimEtaAccY(0), fhPrimEtaPhi(0), fhPrimEtaAccPhi(0),
fhPrimOmegaPt(0), fhPrimOmegaAccPt(0), fhPrimOmegaY(0), fhPrimOmegaAccY(0), 
fhPrimOmegaPhi(0), fhPrimOmegaAccPhi(0),
fCalorimeter("")
{
 //Default Ctor
 InitParameters();
}

//______________________________________________________________________________
AliAnaNeutralMeson::AliAnaNeutralMeson(const AliAnaNeutralMeson & ex) : AliAnaPartCorrBaseClass(ex),  
fAnaPi0Eta(ex.fAnaPi0Eta), fAnaOmega(ex.fAnaOmega), fNPID(ex.fNPID),fNmaxMixEv(ex.fNmaxMixEv), fNAsy(ex.fNAsy),fNMod(ex.fNMod),
fNHistos(ex.fNHistos), fPerPtBin(ex.fPerPtBin),
fAsyCut(ex.fAsyCut), fModCut(ex.fModCut),
fDist(ex.fDist), fAsy(ex.fAsy), fPt(ex.fPt), fMass(ex.fMass),
fInvMassCut(ex.fInvMassCut), 
fNbinsAsy(ex.fNbinsAsy), fMinAsy(ex.fMinAsy), fMaxAsy(ex.fMaxAsy),
fNbinsPt(ex.fNbinsPt), fMinPt(ex.fMinPt), fMaxPt(ex.fMaxPt),
fNbinsM(ex.fNbinsM), fMinM(ex.fMinM), fMaxM(ex.fMaxM),
fEventsListPi0Eta(ex.fEventsListPi0Eta), fEventsListOmega(ex.fEventsListOmega),
fPi0Mass(ex.fPi0Mass), fPi0MassPeakWidthCut(ex.fPi0MassPeakWidthCut), fhNClusters(ex.fhNClusters),
fhRecPhoton(ex.fhRecPhoton), fhRecPhotonEtaPhi(ex.fhRecPhotonEtaPhi),
fReal2Gamma(ex.fReal2Gamma), fMix2Gamma(ex.fMix2Gamma),fRealOmega(ex.fRealOmega),
fMixOmegaA(ex.fMixOmegaA),fMixOmegaB(ex.fMixOmegaB), fMixOmegaC(ex.fMixOmegaC),
fRealTwoGammaAsyPtM(ex.fRealTwoGammaAsyPtM), fMixTwoGammaAsyPtM(ex.fMixTwoGammaAsyPtM),
fRealPi0GammaAsyPtM(ex.fRealPi0GammaAsyPtM), fMixAPi0GammaAsyPtM(ex.fMixAPi0GammaAsyPtM),
fMixBPi0GammaAsyPtM(ex.fMixBPi0GammaAsyPtM), fMixCPi0GammaAsyPtM(ex.fMixCPi0GammaAsyPtM),
fhPrimPhotonPt(ex.fhPrimPhotonPt), fhPrimPhotonAccPt(ex.fhPrimPhotonAccPt), fhPrimPhotonY(ex.fhPrimPhotonY), 
fhPrimPhotonAccY(ex.fhPrimPhotonAccY), fhPrimPhotonPhi(ex.fhPrimPhotonPhi), 
fhPrimPhotonAccPhi(ex.fhPrimPhotonAccPhi),
fhPrimPi0Pt(ex.fhPrimPi0Pt), fhPrimPi0AccPt(ex.fhPrimPi0AccPt), fhPrimPi0Y(ex.fhPrimPi0Y),
fhPrimPi0AccY(ex.fhPrimPi0AccY), fhPrimPi0Phi(ex.fhPrimPi0Phi), fhPrimPi0AccPhi(ex.fhPrimPi0AccPhi), 
fhPrimEtaPt(ex.fhPrimEtaPt), fhPrimEtaAccPt(ex.fhPrimEtaAccPt), fhPrimEtaY(ex.fhPrimEtaY), 
fhPrimEtaAccY(ex.fhPrimEtaAccY), fhPrimEtaPhi(ex.fhPrimEtaPhi), fhPrimEtaAccPhi(ex.fhPrimEtaAccPhi),
fhPrimOmegaPt(ex.fhPrimOmegaPt), fhPrimOmegaAccPt(ex.fhPrimOmegaAccPt), fhPrimOmegaY(ex.fhPrimOmegaY),
fhPrimOmegaAccY(ex.fhPrimOmegaAccY), fhPrimOmegaPhi(ex.fhPrimOmegaPhi), fhPrimOmegaAccPhi(ex.fhPrimOmegaAccPhi),
fCalorimeter("")
{
 // cpy ctor
 //Do not need it
}

//______________________________________________________________________________
AliAnaNeutralMeson & AliAnaNeutralMeson::operator = (const AliAnaNeutralMeson & ex)
{
 // assignment operator

 if(this == &ex)return *this;
   ((AliAnaPartCorrBaseClass *)this)->operator=(ex);
	
  fAnaPi0Eta = ex.fAnaPi0Eta;  fAnaOmega = ex.fAnaOmega;
  fNPID=ex.fNPID; fNmaxMixEv=ex.fNmaxMixEv; fNAsy=ex.fNAsy;fNMod=ex.fNMod;
  fNHistos=ex.fNHistos; fPerPtBin=ex.fPerPtBin;
  fAsyCut=ex.fAsyCut; fModCut=ex.fModCut;
  fDist=ex.fDist; fAsy=ex.fAsy; fPt=ex.fPt; fMass=ex.fMass;
  fInvMassCut=ex.fInvMassCut;
  fNbinsAsy=ex.fNbinsAsy; fMinAsy=ex.fMinAsy; fMaxAsy=ex.fMaxAsy;
  fNbinsPt=ex.fNbinsPt; fMinPt=ex.fMinPt; fMaxPt=ex.fMaxPt;
  fNbinsM=ex.fNbinsM; fMinM=ex.fMinM; fMaxM=ex.fMaxM;
  fEventsListPi0Eta=ex.fEventsListPi0Eta; fEventsListOmega=ex.fEventsListOmega;
  fPi0Mass=ex.fPi0Mass; fPi0MassPeakWidthCut=ex.fPi0MassPeakWidthCut; fhNClusters =ex.fhNClusters;
  fhRecPhoton=ex.fhRecPhoton; fhRecPhotonEtaPhi =ex.fhRecPhotonEtaPhi;
  fReal2Gamma=ex.fReal2Gamma; fMix2Gamma=ex.fMix2Gamma;fRealOmega=ex.fRealOmega;
  fMixOmegaA=ex.fMixOmegaA; fMixOmegaB=ex.fMixOmegaB; fMixOmegaC=ex.fMixOmegaC;
  fRealTwoGammaAsyPtM=ex.fRealTwoGammaAsyPtM; fMixTwoGammaAsyPtM =ex.fMixTwoGammaAsyPtM;
  fRealPi0GammaAsyPtM=ex.fRealPi0GammaAsyPtM; fMixAPi0GammaAsyPtM=ex.fMixAPi0GammaAsyPtM;
  fMixBPi0GammaAsyPtM=ex.fMixBPi0GammaAsyPtM; fMixCPi0GammaAsyPtM=ex.fMixCPi0GammaAsyPtM;
  fhPrimPhotonPt=ex.fhPrimPhotonPt; fhPrimPhotonAccPt=ex.fhPrimPhotonAccPt; fhPrimPhotonY=ex.fhPrimPhotonY;
  fhPrimPhotonAccY=ex.fhPrimPhotonAccY; fhPrimPhotonPhi=ex.fhPrimPhotonPhi;
  fhPrimPhotonAccPhi=ex.fhPrimPhotonAccPhi;
  fhPrimPi0Pt=ex.fhPrimPi0Pt; fhPrimPi0AccPt=ex.fhPrimPi0AccPt; fhPrimPi0Y=ex.fhPrimPi0Y;
  fhPrimPi0AccY=ex.fhPrimPi0AccY; fhPrimPi0Phi=ex.fhPrimPi0Phi; fhPrimPi0AccPhi=ex.fhPrimPi0AccPhi;
  fhPrimEtaPt=ex.fhPrimEtaPt; fhPrimEtaAccPt=ex.fhPrimEtaAccPt; fhPrimEtaY=ex.fhPrimEtaY;
  fhPrimEtaAccY=ex.fhPrimEtaAccY; fhPrimEtaPhi=ex.fhPrimEtaPhi; fhPrimEtaAccPhi=ex.fhPrimEtaAccPhi;
  fhPrimOmegaPt=ex.fhPrimOmegaPt; fhPrimOmegaAccPt=ex.fhPrimOmegaAccPt; fhPrimOmegaY=ex.fhPrimOmegaY;
  fhPrimOmegaAccY=ex.fhPrimOmegaAccY; fhPrimOmegaPhi=ex.fhPrimOmegaPhi; fhPrimOmegaAccPhi=ex.fhPrimOmegaAccPhi;
   	
  return *this;
	
}

//______________________________________________________________________________
AliAnaNeutralMeson::~AliAnaNeutralMeson() {

  //dtor
	
 if(fEventsListPi0Eta){
     delete[] fEventsListPi0Eta;
     fEventsListPi0Eta=0 ;
  }
 if(fEventsListOmega){
     delete[] fEventsListOmega;
     fEventsListOmega=0 ;
 }	

#ifdef __PHOSGEO__
    if(fPHOSGeo) delete fPHOSGeo ;
#endif	
}

//______________________________________________________________________________
void AliAnaNeutralMeson::SetCalorimeter(TString det)
{
 //set the detector
 //for PHOS,  we consider three cases with different PHOS number 
 //for EMCAL, we use the default one, I don't know, needs to be confirmed
 Int_t n=0;
 Int_t * nmod=0;
 if(det == "PHOS"){  
     n      = 3;
     nmod= new Int_t [n];
     nmod[0]=1; nmod[1]=3; nmod[2]=5;
  }
  
  
  if(det == "EMCAL") { //the default number of EMCAL modules
     n = 1;
     nmod= new Int_t [n]; 
     nmod[0]=0; 
  }  
  fCalorimeter = det;
  fNMod = n;
  fModCut = nmod;
}

//______________________________________________________________________________
void AliAnaNeutralMeson::InitParameters()
{
//Init parameters when first called the analysis
//Set default parameters
  SetInputAODName("photons");
  AddToHistogramsName("AnaNeutralMesons_");
	
  fInvMassCut = 1.;
 
  fAnaPi0Eta=kTRUE;
  fAnaOmega=kTRUE;
 
  fNPID      = 9;  
  fNAsy      = 4;
  fAsyCut = new Double_t[fNAsy] ;
  fAsyCut[0]=0.7; fAsyCut[1]=0.8; fAsyCut[2]=0.9; fAsyCut[3]=1;
 
  fNHistos=20; // set the max pt range from 0 to 20 GeV
  fPerPtBin=1.;//set the invariant mass spectrum in pt bin

  fNmaxMixEv = 6; //event buffer size
  fPi0Mass=0.135; //nominal pi0 mass in GeV/c^2
  fPi0MassPeakWidthCut=0.015; //the Pi0MassCut should dependent on pt, here we assume it 0.015 
 
  fNbinsAsy =200; 
  fMinAsy =0;
  fMaxAsy =1;
  fNbinsPt =200;
  fMinPt = 0;
  fMaxPt=20;
  fNbinsM = 200;
  fMinM=0;
  fMaxM=1.;
}

//______________________________________________________________________________
void AliAnaNeutralMeson::Init()
{  
//Init some data members needed in analysis
 
  fEventsListPi0Eta = new TList ;
  fEventsListOmega = new TList ;
	
#ifdef __PHOSGEO__
  printf("PHOS geometry initialized!\n");
  fPHOSGeo = new AliPHOSGeoUtils("PHOSgeo") ;
#endif	
	
}

//______________________________________________________________________________
TList * AliAnaNeutralMeson::GetCreateOutputObjects()
{  
 // Create histograms to be saved in output file and 
 // store them in fOutputContainer
	
 TList * outputContainer = new TList() ; 
 outputContainer->SetName(GetName()); 
 
 char key[255] ;
 char title[255] ;
 const char * detector= fCalorimeter;

 fhNClusters = new TH1I("fhNClusters","N clusters per event",100,0,100); //
 outputContainer->Add(fhNClusters);

 fhRecPhoton  = new TH1D*[fNPID*fNMod];
 fhRecPhotonEtaPhi= new TH2F*[fNPID*fNMod];
 // plot the photon cluster pt, eta and phi distribution by using the PID and module cut 
 for(Int_t ipid=0;ipid<fNPID;ipid++){
    for(Int_t jmod=0;jmod<fNMod;jmod++){
       sprintf(key,"fhRecPhoton_%s_pid_%d_nmod_%d",detector,ipid,fModCut[jmod]);
       sprintf(title,"Rec. photon with pid#%d in %d modules with %s",ipid, fModCut[jmod],detector );
       fhRecPhoton[ipid*fNMod+jmod] = new TH1D(key,title,fNbinsPt,fMinPt, fMaxPt);
       fhRecPhoton[ipid*fNMod+jmod]->GetXaxis()->SetTitle("pt (GeV/c)");
       fhRecPhoton[ipid*fNMod+jmod]->GetYaxis()->SetTitle("Entries");
       outputContainer->Add(fhRecPhoton[ipid*fNMod+jmod]);
 
       sprintf(key,"fhRecPhotonEtaPhi_%s_pid_%d_nmod_%d",detector, ipid,fModCut[jmod]);
       sprintf(title,"Rec. photon eta vs phi with pid#%d in %d modules with %s",ipid, fModCut[jmod], detector);
       fhRecPhotonEtaPhi[ipid*fNMod+jmod] = new TH2F(key,title,200,-1,1,200,0,7);
       fhRecPhotonEtaPhi[ipid*fNMod+jmod]->GetXaxis()->SetTitle("#eta");
       fhRecPhotonEtaPhi[ipid*fNMod+jmod]->GetYaxis()->SetTitle("#phi (rad.)");
       outputContainer->Add(fhRecPhotonEtaPhi[ipid*fNMod+jmod]);
    }
 } 


 if(fAnaPi0Eta){
    fReal2Gamma = new TH1F*[fNPID*fNAsy*fNMod*fNHistos];
    fMix2Gamma  = new TH1F*[fNPID*fNAsy*fNMod*fNHistos];
    
    fRealTwoGammaAsyPtM = new TH3F*[fNPID];
    fMixTwoGammaAsyPtM  = new TH3F*[fNPID];
    for(Int_t ipid=0;ipid<fNPID;ipid++){
       sprintf(key,"RealTwoGammaAsyPtM_%s_pid_%d",detector,ipid);
       sprintf(title, "%s Real 2#gamma asy:pt:ivm with pid#%d",detector,ipid);
       fRealTwoGammaAsyPtM[ipid] = new TH3F(key, title,
                               fNbinsAsy,fMinAsy, fMaxAsy, fNbinsPt,fMinPt, fMaxPt,fNbinsM,fMinM, fMaxM);
       fRealTwoGammaAsyPtM[ipid]->GetYaxis()->SetTitle("2#gamma pt");
       fRealTwoGammaAsyPtM[ipid]->GetXaxis()->SetTitle("A=|E_{#gamma1}-E_{#gamma2}|/(E_{#gamma1}+E_{#gamma2})");
       fRealTwoGammaAsyPtM[ipid]->GetZaxis()->SetTitle("M_{2#gamma}");
       outputContainer->Add(fRealTwoGammaAsyPtM[ipid]);
 
       sprintf(key,"MixTwoGammaAsyPtM_%s_pid_%d",detector,ipid);
       sprintf(title, "%s Mix 2#gamma asy:pt:ivm with pid#%d",detector,ipid);
       fMixTwoGammaAsyPtM[ipid] = new TH3F(key, title,
                               fNbinsAsy,fMinAsy, fMaxAsy, fNbinsPt,fMinPt, fMaxPt,fNbinsM,fMinM, fMaxM);
       fMixTwoGammaAsyPtM[ipid]->GetYaxis()->SetTitle("2#gamma pt");
       fMixTwoGammaAsyPtM[ipid]->GetXaxis()->SetTitle("A=|E_{#gamma1}-E_{#gamma2}|/(E_{#gamma1}+E_{#gamma2})");
       fMixTwoGammaAsyPtM[ipid]->GetZaxis()->SetTitle("M_{2#gamma}");
       outputContainer->Add(fMixTwoGammaAsyPtM[ipid]);   
 
      for(Int_t jmod=0; jmod<fNMod; jmod++){
         for(Int_t kasy=0;kasy<fNAsy;kasy++){
            for(Int_t mhist=0;mhist<fNHistos;mhist++){
               Double_t pt1=mhist*fPerPtBin;
               Double_t pt2=(mhist+1)*fPerPtBin;
               Int_t index = ipid*fNMod*fNAsy*fNHistos+jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
 
               sprintf(key,"fReal2Gamma_%s_pid_%d_mod_%d_asy_%d_pt_%d",detector, ipid,fModCut[jmod],kasy,mhist);
               sprintf(title,"real 2#gamma IVM with pid#%d asy<%2.1f in %d modules with %s at %2.2f<pt<%2.2f",
                                                   ipid,fAsyCut[kasy], fModCut[jmod],detector, pt1, pt2);
               fReal2Gamma[index] = new TH1F(key, title,fNbinsM,fMinM,fMaxM);
               fReal2Gamma[index]->GetYaxis()->SetTitle("dN/dM_{#gamma#gamma}");
               fReal2Gamma[index]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/c^{2})");
               outputContainer->Add(fReal2Gamma[index]);
 
               sprintf(key,"fMix2Gamma_%s_pid_%d_mod_%d_asy_%d_pt_%d",detector,ipid,fModCut[jmod],kasy,mhist);
               sprintf(title,"Mix 2#gamma IVM with pid#%d asy<%2.1f in %d modules with %s at %2.2f<pt<%2.2f",
                                                  ipid,fAsyCut[kasy],fModCut[jmod],detector,pt1, pt2);
               fMix2Gamma[index] = new TH1F(key, title,fNbinsM,fMinM,fMaxM);
               fMix2Gamma[index]->GetYaxis()->SetTitle("dN/dM_{#gamma#gamma}");
               fMix2Gamma[index]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/c^{2})");
               outputContainer->Add(fMix2Gamma[index]);
            } 
         }
      }
   }
 }
 
 if(fAnaOmega){
    fRealOmega  = new TH1F*[fNPID*fNAsy*fNMod*fNHistos];
    fMixOmegaA  = new TH1F*[fNPID*fNAsy*fNMod*fNHistos];
    fMixOmegaB  = new TH1F*[fNPID*fNAsy*fNMod*fNHistos];
    fMixOmegaC  = new TH1F*[fNPID*fNAsy*fNMod*fNHistos];
    fRealPi0GammaAsyPtM=new TH3F*[fNPID];
    fMixAPi0GammaAsyPtM=new TH3F*[fNPID];
    fMixBPi0GammaAsyPtM=new TH3F*[fNPID];
    fMixCPi0GammaAsyPtM=new TH3F*[fNPID];
   
    for(Int_t ipid=0;ipid<fNPID;ipid++){
       sprintf(key,"RealPi0GammaAsyPtM_%s_pid_%d",detector,ipid);
       sprintf(title,"%s Real #omega->#pi^{0}#gamma asy:pt:ivm with pid#%d",detector,ipid);
       fRealPi0GammaAsyPtM[ipid] = new TH3F(key, title,
                               fNbinsAsy,fMinAsy, fMaxAsy, fNbinsPt,fMinPt, fMaxPt,fNbinsM,fMinM, fMaxM);
       fRealPi0GammaAsyPtM[ipid]->GetYaxis()->SetTitle("#pi^{0}#gamma pt");
       fRealPi0GammaAsyPtM[ipid]->GetXaxis()->SetTitle("A=|E_{#pi^{0}}-E_{#gamma}|/(E_{#pi^{0}}+E_{#gamma})");
       fRealPi0GammaAsyPtM[ipid]->GetZaxis()->SetTitle("M_{#pi^{0}#gamma}");
       outputContainer->Add(fRealPi0GammaAsyPtM[ipid]);
 
       sprintf(key,"MixAPi0GammaAsyPtM_%s_pid_%d",detector, ipid);
       sprintf(title,"%s MixA #omega->#pi^{0}#gamma asy:pt:ivm with pid#%d",detector,ipid);
       fMixAPi0GammaAsyPtM[ipid] = new TH3F(key,title,
                              fNbinsAsy,fMinAsy, fMaxAsy, fNbinsPt,fMinPt, fMaxPt,fNbinsM,fMinM, fMaxM);
       fMixAPi0GammaAsyPtM[ipid]->GetYaxis()->SetTitle("#pi^{0}#gamma pt");
       fMixAPi0GammaAsyPtM[ipid]->GetXaxis()->SetTitle("A=|E_{#pi^{0}}-E_{#gamma}|/(E_{#pi^{0}}+E_{#gamma})");
       fMixAPi0GammaAsyPtM[ipid]->GetZaxis()->SetTitle("M_{#pi^{0}#gamma}");
       outputContainer->Add(fMixAPi0GammaAsyPtM[ipid]);

       sprintf(key,"MixBPi0GammaAsyPtM_%s_pid_%d",detector,ipid);
       sprintf(title,"%s MixB #omega->#pi^{0}#gamma asy:pt:ivm with pid#%d",detector,ipid); 
       fMixBPi0GammaAsyPtM[ipid] = new TH3F(key, title,
                               fNbinsAsy,fMinAsy, fMaxAsy, fNbinsPt,fMinPt, fMaxPt,fNbinsM,fMinM, fMaxM);
       fMixBPi0GammaAsyPtM[ipid]->GetYaxis()->SetTitle("#pi^{0}#gamma pt");
       fMixBPi0GammaAsyPtM[ipid]->GetXaxis()->SetTitle("A=|E_{#pi^{0}}-E_{#gamma}|/(E_{#pi^{0}}+E_{#gamma})");
       fMixBPi0GammaAsyPtM[ipid]->GetZaxis()->SetTitle("M_{#pi^{0}#gamma}");
       outputContainer->Add(fMixBPi0GammaAsyPtM[ipid]);
 
       sprintf(key,"MixCPi0GammaAsyPtM_%s_pid_%d",detector,ipid);
       sprintf(title,"%s MixC #omega->#pi^{0}#gamma asy:pt:ivm with pid#%d",detector,ipid);
       fMixCPi0GammaAsyPtM[ipid] = new TH3F(key, title,
                               fNbinsAsy,fMinAsy, fMaxAsy, fNbinsPt,fMinPt, fMaxPt,fNbinsM,fMinM, fMaxM);
       fMixCPi0GammaAsyPtM[ipid]->GetYaxis()->SetTitle("#pi^{0}#gamma pt");
       fMixCPi0GammaAsyPtM[ipid]->GetXaxis()->SetTitle("A=|E_{#pi^{0}}-E_{#gamma}|/(E_{#pi^{0}}+E_{#gamma})");
       fMixCPi0GammaAsyPtM[ipid]->GetZaxis()->SetTitle("M_{#pi^{0}#gamma}");
       outputContainer->Add(fMixCPi0GammaAsyPtM[ipid]);

       for(Int_t jmod=0; jmod<fNMod; jmod++){
         for(Int_t kasy=0;kasy<fNAsy;kasy++){
            for(Int_t mhist=0;mhist<fNHistos;mhist++){
               Double_t pt1=mhist*fPerPtBin;
               Double_t pt2=(mhist+1)*fPerPtBin;
               Int_t index = ipid*fNMod*fNAsy*fNHistos+jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
               sprintf(key,"fRealOmega_%s_pid_%d_mod_%d_asy_%d_pt_%d",detector,ipid,fModCut[jmod],kasy,mhist);
               sprintf(title,"Real #omega with pid#%d asy<%2.1f in %d modules with %s at %2.2f<pt<%2.2f",
                                              ipid,fAsyCut[kasy],fModCut[jmod], detector, pt1, pt2);
               fRealOmega[index] = new TH1F(key, title,fNbinsM,fMinM,fMaxM);
               fRealOmega[index]->GetYaxis()->SetTitle("dN/dM_{#pi^{0}#gamma} ");       
               fRealOmega[index]->GetXaxis()->SetTitle("M_{#pi^{0}#gamma} (GeV/c^{2})");
               outputContainer->Add(fRealOmega[index]);
              
               sprintf(key,"fMixOmegaA_%s_pid_%d_mod_%d_asy_%d_pt_%d",detector, ipid,fModCut[jmod],kasy,mhist);
               sprintf(title,"Mix #omega A with pid#%d asy<%2.1f in %d modules with %s at %2.2f<pt<%2.2f",
                                               ipid,fAsyCut[kasy],fModCut[jmod],detector,pt1, pt2);
               fMixOmegaA[index] = new TH1F(key, title,fNbinsM,fMinM,fMaxM);
               fMixOmegaA[index]->GetYaxis()->SetTitle("dN/dM_{#pi^{0}#gamma}");
               fMixOmegaA[index]->GetXaxis()->SetTitle("M_{#pi^{0}#gamma} (GeV/c^{2})");
               outputContainer->Add(fMixOmegaA[index]);
              
               sprintf(key,"fMixOmegaB_%s_pid_%d_mod_%d_asy_%d_pt_%d",detector,ipid,fModCut[jmod],kasy,mhist);
               sprintf(title,"Mix #omega B with pid#%d asy<%2.1f in %d modules with %s at %2.2f<pt<%2.2f",
                                               ipid,fAsyCut[kasy],fModCut[jmod],detector,pt1, pt2);
               fMixOmegaB[index] = new TH1F(key, title,fNbinsM,fMinM,fMaxM);
               fMixOmegaB[index]->GetYaxis()->SetTitle("dN/dM_{#pi^{0}#gamma}");
               fMixOmegaB[index]->GetXaxis()->SetTitle("M_{#pi^{0}#gamma} (GeV/c^{2})");
               outputContainer->Add(fMixOmegaB[index]);
              
               sprintf(key,"fMixOmegaC_%s_pid_%d_mod_%d_asy_%d_pt_%d",detector,ipid,fModCut[jmod],kasy,mhist);
               sprintf(title,"Mix #omega C with pid#%d asy<%2.1f in %d modules with %s at %2.2f<pt<%2.2f",
                                               ipid,fAsyCut[kasy],fModCut[jmod],detector, pt1, pt2);
               fMixOmegaC[index] = new TH1F(key, title,fNbinsM,fMinM,fMaxM);
               fMixOmegaC[index]->GetYaxis()->SetTitle("dN/dM_{#pi^{0}#gamma}");
               fMixOmegaC[index]->GetXaxis()->SetTitle("M_{#pi^{0}#gamma} (GeV/c^{2})");
               outputContainer->Add(fMixOmegaC[index]);
            }
         }
      }
   }
 }

 //Histograms filled only if MC data is requested 	
 //pt distribution
 //currently only PHOS included                                                                 
 if(fCalorimeter=="PHOS" && (IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC)) ){
      fhPrimPhotonPt = new TH1D("hPrimPhotonPt","Primary phton pt at |y|=1",fNbinsPt,fMinPt, fMaxPt) ;
      outputContainer->Add(fhPrimPhotonPt);

      fhPrimPhotonAccPt = new TH1D("hPrimPhotonAccPt","Primary photon pt in acceptance",fNbinsPt,fMinPt, fMaxPt);
      outputContainer->Add(fhPrimPhotonAccPt);
 
      fhPrimPi0Pt = new TH1D("hPrimPi0Pt","Primary pi0 pt at |y|=1",fNbinsPt,fMinPt, fMaxPt) ;
      outputContainer->Add(fhPrimPi0Pt);

      fhPrimPi0AccPt = new TH1D("hPrimPi0AccPt","Primary pi0 pt with both photons in acceptance",fNbinsPt,fMinPt, fMaxPt);
      outputContainer->Add(fhPrimPi0AccPt);
 
      fhPrimEtaPt = new TH1D("hPrimEtaPt","Primary eta pt at |y|=1",fNbinsPt,fMinPt, fMaxPt) ;
      outputContainer->Add(fhPrimEtaPt);

      fhPrimEtaAccPt = new TH1D("hPrimEtaAccPt","Primary eta pt with both photons in acceptance",fNbinsPt,fMinPt, fMaxPt) ;
      outputContainer->Add(fhPrimEtaAccPt);
 
      fhPrimOmegaPt = new TH1D("hPrimOmegaPt","Primary Omega pt at |y|=1",fNbinsPt,fMinPt, fMaxPt) ;
      outputContainer->Add(fhPrimOmegaPt);

      fhPrimOmegaAccPt = new TH1D("hPrimOmegaAccPt","Primary Omega pt with omega->pi0+gamma in acceptance",fNbinsPt,fMinPt, fMaxPt);
      outputContainer->Add(fhPrimOmegaAccPt);
 
      //for photon
      fhPrimPhotonY = new TH1D("hPrimPhotonY","Rapidity of primary pi0",100,-5.,5.) ;
      outputContainer->Add(fhPrimPhotonY) ;
 
      fhPrimPhotonAccY = new TH1D("hPrimPhotonAccY","Rapidity of primary pi0",100,-5.,5.) ;
      outputContainer->Add(fhPrimPhotonAccY) ;
 
      fhPrimPhotonPhi = new TH1D("hPrimPhotonPhi","Azimithal of primary photon",180,0.,360.) ;
      outputContainer->Add(fhPrimPhotonPhi) ;
 
      fhPrimPhotonAccPhi = new TH1D("hPrimPhotonAccPhi","Azimithal of primary photon in Acceptance",180,-0.,360.) ;
      outputContainer->Add(fhPrimPhotonAccPhi) ;
 
      //for pi0
      fhPrimPi0Y = new TH1D("hPrimPi0Y","Rapidity of primary pi0",100,-5.,5.) ;
      outputContainer->Add(fhPrimPi0Y) ;
 
      fhPrimPi0AccY = new TH1D("hPrimPi0AccY","Rapidity of primary pi0",100,-5.,5.) ;
      outputContainer->Add(fhPrimPi0AccY) ;
 
      fhPrimPi0Phi = new TH1D("hPrimPi0Phi","Azimithal of primary pi0",180,0.,360.) ;
      outputContainer->Add(fhPrimPi0Phi) ;
 
      fhPrimPi0AccPhi = new TH1D("hPrimPi0AccPhi","Azimithal of primary pi0 with accepted daughters",180,-0.,360.) ;
      outputContainer->Add(fhPrimPi0AccPhi) ;
 
      //for eta
      fhPrimEtaY = new TH1D("hPrimEtaY","Rapidity of primary Eta",100,-5.,5.) ;
      outputContainer->Add(fhPrimEtaY) ;
 
      fhPrimEtaAccY = new TH1D("hPrimEtaAccY","Rapidity of primary eta",100,-5.,5.) ;
      outputContainer->Add(fhPrimEtaAccY) ;
 
      fhPrimEtaPhi = new TH1D("hPrimEtaPhi","Azimithal of primary eta",180,0.,360.) ;
      outputContainer->Add(fhPrimEtaPhi) ;
 
      fhPrimEtaAccPhi = new TH1D("hPrimEtaAccPhi","Azimithal of primary eta with accepted daughters",180,-0.,360.) ;
      outputContainer->Add(fhPrimEtaAccPhi) ;
 
      //for omega
      fhPrimOmegaY = new TH1D("hPrimOmegaY","Rapidity of primary Omega",100,-5.,5.) ;
      outputContainer->Add(fhPrimOmegaY) ;
 
      fhPrimOmegaAccY = new TH1D("hPrimOmegaAccY","Rapidity of primary omega->pi0+gamma",100,-5.,5.) ;
      outputContainer->Add(fhPrimOmegaAccY) ;
 
      fhPrimOmegaPhi = new TH1D("hPrimOmegaPhi","Azimithal of primary omega",180,0.,360.) ;
      outputContainer->Add(fhPrimOmegaPhi) ;
 
      fhPrimOmegaAccPhi = new TH1D("hPrimOmegaAccPhi","Azimithal of primary omega->pi0+gamma with accepted daughters",180,-0.,360.) ;
      outputContainer->Add(fhPrimOmegaAccPhi);
 }
 return outputContainer;
}

//______________________________________________________________________________
void AliAnaNeutralMeson::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set for the analysis
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  printf("performing the analysis of 2gamma IVM  %d\n", fAnaPi0Eta); //IVM: Invariant Mass
  printf("performing the analysis of pi0+gamma->3gamma IVM  %d\n", fAnaOmega);
  printf("Number of Events to be mixed: %d \n",fNmaxMixEv) ;
  printf("Number of different PID used:  %d \n",fNPID) ;
  printf("Number of different Asy cut:  %d\n",fNAsy);
  printf("Asy bins:  %d  Min:  %2.1f  Max:  %2.1f\n", fNbinsAsy, fMinAsy, fMaxAsy);
  printf("Pt bins:  %d  Min:  %2.1f  Max:  %2.1f  GeV/c\n", fNbinsPt, fMinPt, fMaxPt);
  printf("IVM bins:  %d  Min:  %2.1f  Max:  %2.1f  GeV/c^2\n", fNbinsM, fMinM, fMaxM);
  printf("produce the %d histograms per %2.1fGeV/c pt bin \n", fNHistos,fPerPtBin);
  printf("Vertec, centrality Cuts, RP are not used at present \n") ;
  printf("------------------------------------------------------\n") ;

} 

//______________________________________________________________________________
void AliAnaNeutralMeson::MakeAnalysisFillHistograms() 
{
 //process event from AOD brach
 //extract pi0, eta and omega analysis
 Int_t iRun=(GetReader()->GetInputEvent())->GetRunNumber() ;
 if(IsBadRun(iRun)) return ;	

 //vertex cut not used yet
 //centrality not used yet
 //reaction plane not used yet

 TClonesArray *  array1 = (TClonesArray*)GetInputAODBranch();
 Int_t nPhot = array1 ->GetEntries();
 fhNClusters->Fill(nPhot);

 //fill photon clusters
 for(Int_t i=0;i<nPhot;i++){
     AliAODPWG4Particle * p1 = (AliAODPWG4Particle*) (GetInputAODBranch()->At(i)) ;
     TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
     for(Int_t ipid=0;ipid<fNPID;ipid++){
        for(Int_t jmod=0;jmod<fNMod;jmod++){
           if(IsPhotonSelected(p1, ipid, fModCut[jmod])){
                fhRecPhoton[ipid*fNMod+jmod]->Fill(photon1.Pt());
                Double_t phi = photon1.Phi() ;
                if(phi<0) phi = photon1.Phi()+2*TMath::Pi();
                fhRecPhotonEtaPhi[ipid*fNMod+jmod]->Fill(photon1.Eta(),phi);
           } 
        }
     }
 }

/////////////////////////////////////////////
//real events
 if(nPhot>=2 && fAnaPi0Eta){
//   printf("real nPhot:  %d \n",nPhot);
 
    RealPi0Eta(array1); //real events for pi0 and eta
 
    TList * evMixListPi0Eta=fEventsListPi0Eta ;        //with cluster larger than 1
    Int_t nMixed = evMixListPi0Eta->GetSize() ;
    ////////////////////////////////////////////
    //Mixing events for pi0 and eta to reconstruct the background
    for(Int_t i1=0; i1<nMixed; i1++){
        TClonesArray* array2= (TClonesArray*) (evMixListPi0Eta->At(i1));
        MixPi0Eta(array1, array2);
     }
    ////////////////////////////////////////////
 
   //event buffer for pi0 and eta
    TClonesArray *currentEvent = new TClonesArray(*array1);
 
    if(currentEvent->GetEntriesFast()>=2){
        evMixListPi0Eta->AddFirst(currentEvent) ;
        currentEvent=0 ; //Now list of particles belongs to buffer and it will be deleted with buffer
        if(evMixListPi0Eta->GetSize()>=fNmaxMixEv) {
            TClonesArray * tmp = (TClonesArray*) (evMixListPi0Eta->Last()) ;
            evMixListPi0Eta->RemoveLast() ;
            delete tmp ;
        }
     }
     else{ //empty event
          delete currentEvent ;
          currentEvent=0 ;
     }
 
 }
 
 if(nPhot>=3 && fAnaOmega) {
    RealOmega(array1); //real events for omega
    TList * evMixListOmega = fEventsListOmega ; //with cluster larger than 2
    Int_t nMixedOmega = evMixListOmega->GetSize() ;
 
    ////////////////////////////////////////////
    //mix for omega
    for(Int_t i1=0; i1<nMixedOmega; i1++){
        TClonesArray* array2= (TClonesArray*) (evMixListOmega->At(i1));
        MixOmegaAB(array1, array2);
        //third event
        for(Int_t i2=i1+1;i2<nMixedOmega;i2++){
            TClonesArray * array3 = (TClonesArray*)(evMixListOmega->At(i2));
//            printf("omega  nPhot:  %d   nPhot2:  %d   nPhot3:  %d \n", nPhot,nPhot2,nPhot3);
            MixOmegaC(array1,array2,array3);
        }
     }
 
    //fill event buffer with 2 more clusters events for omega extraction
    TClonesArray *currentEventOmega = new TClonesArray(*array1);
    if(currentEventOmega->GetEntriesFast()>=3){
         evMixListOmega->AddFirst(currentEventOmega) ;
         currentEventOmega=0 ; //Now list of particles belongs to buffer and it will be deleted with buffer
         if(evMixListOmega->GetSize()>=fNmaxMixEv) {
             TClonesArray * tmp = (TClonesArray*) (evMixListOmega->Last()) ;
             evMixListOmega->RemoveLast() ;
             delete tmp ;
          }
     }  
     else{ //empty event
           delete currentEventOmega ;
           currentEventOmega=0 ;
     }      
  
 }          

 if(fCalorimeter=="PHOS" && (IsDataMC() || (GetReader()->GetDataType() == AliCaloTrackReader::kMC)) ){
     AliStack * stack = GetMCStack();
     //photon acceptance
     PhotonAcceptance(stack);
     //pi0 and eta acceptance
     Pi0EtaAcceptance(stack);
     //omega acceptance      
     OmegaAcceptance(stack);
 }

}

//______________________________________________________________________________
Bool_t AliAnaNeutralMeson::IsPhotonSelected(AliAODPWG4Particle *p1, Int_t ipid, Int_t nmod)
{
 //select the photon to be analyzed based on pid and which module 

 TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());

 if(fCalorimeter == "PHOS"){
     Int_t mod = 20 ;
#ifdef __PHOSGEO__
	Double_t vtx[3]={0,0,0};   
//        if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vtx);
	Double_t x = 0 ,z = 0 ;
	fPHOSGeo->ImpactOnEmc(vtx,p1->Theta(),p1->Phi(), mod,z,x) ;
#endif
	  if(mod<=nmod && p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) {return kTRUE;}
	}
	 
  else if(fCalorimeter == "EMCAL" && p1->IsPIDOK(ipid,AliCaloPID::kPhoton) ) {return kTRUE;} //default EMCAL module
	return kFALSE;
}

//______________________________________________________________________________
void AliAnaNeutralMeson::Get2GammaAsyDistPtM(AliAODPWG4Particle *p1, AliAODPWG4Particle *p2,
                                            Double_t &asy, Double_t &dist, Double_t &pt, Double_t &mass)
{
  //get the asy, dist, pair pt and pair inv mass from the 2gammas
  TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
  TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
  if(photon1.Pt()>0 && photon2.Pt()>0){
     asy=TMath::Abs(photon1.E()-photon2.E())/(photon1.E()+photon2.E());
     Double_t phi1 = photon1.Phi();
     Double_t phi2 = photon2.Phi();
     Double_t eta1 = photon1.Eta();
     Double_t eta2 = photon2.Eta();
     dist=TMath::Sqrt((phi1-phi2)*(phi1-phi2)+(eta1-eta2)*(eta1-eta2));
     pt=(photon1+photon2).Pt();
     mass=(photon1+photon2).M();
//   printf("asy:  %2.3f   dist:   %2.3f   pt:  %2.3f   m:  %2.3f\n",asy,dist,pt,mass);
   }
   else { asy=0; dist =0; pt =0; mass=0;}
}

//______________________________________________________________________________
void AliAnaNeutralMeson::GetAsyDistPtM(TLorentzVector photon1, TLorentzVector photon2,
                                            Double_t &asy, Double_t &dist, Double_t &pt, Double_t &mass)
{
  //get the asy, dist, pair pt and pair inv mass from the 2gammas 
  if(photon1.Pt()>0 && photon2.Pt()>0){
     asy=TMath::Abs(photon1.E()-photon2.E())/(photon1.E()+photon2.E());
     Double_t phi1 = photon1.Phi();
     Double_t phi2 = photon2.Phi();
     Double_t eta1 = photon1.Eta();
     Double_t eta2 = photon2.Eta();
     dist=TMath::Sqrt((phi1-phi2)*(phi1-phi2)+(eta1-eta2)*(eta1-eta2));
     pt=(photon1+photon2).Pt();
     mass=(photon1+photon2).M();
//   printf("asy:  %2.3f   dist:   %2.3f   pt:  %2.3f   m:  %2.3f\n",asy,dist,pt,mass);
   }
   else { asy=0; dist =0; pt =0; mass=0;}
}

//______________________________________________________________________________ 
void AliAnaNeutralMeson::RealPi0Eta(TClonesArray * array)
{
 //get the 2gamma invariant mass distribution from one event
 Int_t nclusters =  array->GetEntries();
 for(Int_t ie=0;ie<nclusters;ie++){
     for(Int_t je=ie+1;je<nclusters;je++){
         AliAODPWG4Particle * p1 = (AliAODPWG4Particle *)array->At(ie);
         AliAODPWG4Particle * p2 = (AliAODPWG4Particle *)array->At(je);
         //Asy vs. pt
         Get2GammaAsyDistPtM(p1,p2,fAsy,fDist,fPt,fMass);
         for(Int_t ipid=0;ipid<fNPID;ipid++){

             if((p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton))){ 
                fRealTwoGammaAsyPtM[ipid]->Fill(fAsy, fPt, fMass); //default module with pid
             }

             for(Int_t jmod=0;jmod<fNMod;jmod++){
                 for(Int_t kasy=0;kasy<fNAsy;kasy++){
                     if(IsPhotonSelected(p1,ipid,fModCut[jmod]) && IsPhotonSelected(p2,ipid,fModCut[jmod])){      
                           Get2GammaAsyDistPtM(p1,p2,fAsy,fDist,fPt,fMass);
                           //filling the histo
                           for(Int_t mhist=0;mhist<fNHistos;mhist++){
                               Double_t pt1 =mhist *fPerPtBin;
                               Double_t pt2 =(mhist+1)*fPerPtBin;
                               Int_t index = ipid*fNMod*fNAsy*fNHistos+jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
                               if(fMass<fInvMassCut && fAsy<fAsyCut[kasy] 
                                  && fPt<pt2 && fPt>pt1 &&fAsy && fDist && fPt && fMass)
                               fReal2Gamma[index]->Fill(fMass);
                           }//hist
                      } //photon selection    
                   }//asy
               }//mod
          }//pid
     }//je
 }//ie
 
}

//______________________________________________________________________________ 
void AliAnaNeutralMeson::MixPi0Eta(TClonesArray * array1, TClonesArray * array2)
{
 //mixing events for pi0 and eta invariant mass analysis
 // printf("mixed events for pi0 and eta invariant mass analysis!\n");
 Int_t nclusters1 =  array1->GetEntries();
 Int_t nclusters2 =  array2->GetEntries();
 for(Int_t ie=0;ie<nclusters1;ie++){
     AliAODPWG4Particle * p1 = (AliAODPWG4Particle *)array1->At(ie);
     for(Int_t je=0;je<nclusters2;je++){
         AliAODPWG4Particle * p2 = (AliAODPWG4Particle *)array2->At(je);
         TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
         TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
         Get2GammaAsyDistPtM(p1,p2, fAsy,fDist,fPt,fMass);

         for(Int_t ipid=0;ipid<fNPID;ipid++){
             if((p1->IsPIDOK(ipid,AliCaloPID::kPhoton)) && (p2->IsPIDOK(ipid,AliCaloPID::kPhoton))){
                fMixTwoGammaAsyPtM[ipid]->Fill(fAsy,fPt,fMass); //default module with pid
             }
             for(Int_t jmod=0;jmod<fNMod;jmod++){ 
                for(Int_t kasy=0;kasy<fNAsy;kasy++){
                    if(IsPhotonSelected(p1,ipid,fModCut[jmod]) && IsPhotonSelected(p2,ipid,fModCut[jmod])){
                        for(Int_t mhist=0;mhist<fNHistos;mhist++){
                            Double_t pt1 =mhist *fPerPtBin;
                            Double_t pt2 =(mhist+1)*fPerPtBin;
                            Int_t index = ipid*fNMod*fNAsy*fNHistos+jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
                            if(fMass<fInvMassCut && fAsy<fAsyCut[kasy] 
                               && fPt<pt2 && fPt>pt1 &&fAsy&& fDist && fPt && fMass)
                            fMix2Gamma[index]->Fill(fMass);
                        } //hist
                    } //photon selected
                } //asy
              } //mod
         }//pid
    } //cluster2
 } //cluster1                 
}                   

//______________________________________________________________________________ 
void AliAnaNeutralMeson::RealOmega(TClonesArray * array)
{
 //omega invariant mass extraction
 //
 // 1. select the pi0 candidate from 2gamma 
 // 2. combine the third photon with the pi0 candidate to get the omega invariant mass
 //
 Int_t nclusters =  array->GetEntries();
 for(Int_t ie=0;ie<nclusters;ie++){//ptc1
     for(Int_t je=ie+1;je<nclusters;je++){ //ptc2
         AliAODPWG4Particle * p1 = (AliAODPWG4Particle *)array->At(ie);
         AliAODPWG4Particle * p2 = (AliAODPWG4Particle *)array->At(je);
         //pid, mod and asy cut
         Double_t f2gammaAsy;
         Get2GammaAsyDistPtM(p1,p2,f2gammaAsy,fDist,fPt,fMass);

         if(TMath::Abs(fMass-fPi0Mass)<fPi0MassPeakWidthCut){//select the pi0 candidate
             for(Int_t ke=0;ke<nclusters;ke++){ //ptc3
                 AliAODPWG4Particle * p3 = (AliAODPWG4Particle *)array->At(ke);              
                 if(p3 != p2 && p3 != p1 ){ ////p3!=p1 p3!=p2
                     TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
                     TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
                     Double_t px12 =(photon1+photon2).Px();
                     Double_t py12 =(photon1+photon2).Py();
                     Double_t pz12 =(photon1+photon2).Pz();
                     Double_t e12 =(photon1+photon2).E();
                     TLorentzVector photon12(px12, py12,pz12,e12);
 
                     TLorentzVector photon3(p3->Px(),p3->Py(),p3->Pz(),p3->E());
                     GetAsyDistPtM(photon12, photon3,fAsy,fDist,fPt,fMass);
                     for(Int_t ipid=0;ipid<fNPID;ipid++){ //pid
                         if(p1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&   
                               p2->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                               p3->IsPIDOK(ipid,AliCaloPID::kPhoton)){
                             fRealPi0GammaAsyPtM[ipid]->Fill(fAsy, fPt, fMass); //default module, without pid
                         }

                         for(Int_t jmod=0;jmod<fNMod;jmod++){ //mod
                             for(Int_t kasy=0;kasy<fNAsy;kasy++){ //asy
                                if(IsPhotonSelected(p1,ipid,fModCut[jmod]) && 
                                   IsPhotonSelected(p2,ipid,fModCut[jmod]) &&
                                   IsPhotonSelected(p3,ipid,fModCut[jmod])){ //photon selection
                                 
                                   //filling the histo
                                  for(Int_t mhist=0;mhist<fNHistos;mhist++){ //hist
                                      Double_t pt1 =mhist *fPerPtBin;
                                      Double_t pt2 =(mhist+1)*fPerPtBin;
                                      Int_t index = ipid*fNMod*fNAsy*fNHistos+
                                                    jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
                                      if(fMass<fInvMassCut && f2gammaAsy<fAsyCut[kasy] && 
                                         fPt<pt2 && fPt>pt1&& fAsy<fAsyCut[kasy] &&fAsy&&fDist&&fPt&&fMass)
                                      fRealOmega[index]->Fill(fMass);
                                   }//hist
                                 } //photon selection
                               } //asy
                          }//mod
                      }//pid
                 } //p3!=p1 p3!=p2
             } //ptc3
         } //pi0 candidate
     }//je  ptc2   
 } //ie ptc1
}


//______________________________________________________________________________
void AliAnaNeutralMeson::MixOmegaAB(TClonesArray * array1, TClonesArray * array2)
{
 //omega background
 //three omega background, we classify it A, B and C
 // --A
 // (r1_event1+r2_event1)+r3_event2
 //
 // --B
 // (r1_event1+r2_event2)+r3_event2
 // 
 // --C
 // (r1_event1+r2_event2)+r3_event3
 //
 Int_t nclusters1 =  array1->GetEntries();
 Int_t nclusters2 =  array2->GetEntries();
 for(Int_t ie=0;ie<nclusters1;ie++){
     //mix Omega A
     for(Int_t je=ie+1;je<nclusters1;je++){ //.............
         AliAODPWG4Particle * p1 = (AliAODPWG4Particle *)array1->At(ie);
         AliAODPWG4Particle * p2 = (AliAODPWG4Particle *)array1->At(je);
         Double_t f2gammaAsy ;
         Get2GammaAsyDistPtM(p1,p2,f2gammaAsy,fDist,fPt,fMass);
         if(TMath::Abs(fMass-fPi0Mass)<fPi0MassPeakWidthCut ){ //pi0 candidate 
            for(Int_t ke=0;ke<nclusters2;ke++){ //third photon
                AliAODPWG4Particle * p3 = (AliAODPWG4Particle *)array2->At(ke);
                TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
                TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
                Double_t px12 =(photon1+photon2).Px();
                Double_t py12 =(photon1+photon2).Py();
                Double_t pz12 =(photon1+photon2).Pz();
                Double_t e12 =(photon1+photon2).E();
                TLorentzVector photon12(px12, py12,pz12,e12);
                TLorentzVector photon3(p3->Px(),p3->Py(),p3->Pz(),p3->E());
                GetAsyDistPtM(photon12, photon3,fAsy,fDist,fPt,fMass);
                //pid, mod and asy cut
                for(Int_t ipid=0;ipid<fNPID;ipid++){ //pid
                    if(p1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                            p2->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                            p3->IsPIDOK(ipid,AliCaloPID::kPhoton)){
                         fMixAPi0GammaAsyPtM[ipid]->Fill(fAsy, fPt, fMass); //default module, with pid
                    }
                    for(Int_t jmod=0;jmod<fNMod;jmod++){ //mod
                        for(Int_t kasy=0;kasy<fNAsy;kasy++){  //asy
                           if(IsPhotonSelected(p1,ipid,fModCut[jmod]) && 
                              IsPhotonSelected(p2,ipid,fModCut[jmod]) &&
                              IsPhotonSelected(p3,ipid,fModCut[jmod])){ //photon selection
                               //fill the hist
                               for(Int_t mhist=0;mhist<fNHistos;mhist++){ //hist
                                   Double_t pt1 =mhist *fPerPtBin;
                                   Double_t pt2 =(mhist+1)*fPerPtBin;
                                   Int_t index = ipid*fNMod*fNAsy*fNHistos+                                                                                        jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
                                   if(fMass<fInvMassCut && f2gammaAsy<fAsyCut[kasy] &&
                                      fPt<pt2 && fPt>pt1&& fAsy<fAsyCut[kasy] &&fAsy&&fDist&&fPt&&fMass)
                                   fMixOmegaA[index]->Fill(fMass);
                              }//hist
                            }//photon selection
                        } //asy
                    }//mod
                } //pid
            } //ptc3
         } //pi0 candidate
      } //mixA je ptc2

     //mix omega B
     for(Int_t je=0;je<nclusters2;je++){ //mixb...............
         AliAODPWG4Particle * p1 = (AliAODPWG4Particle *)array1->At(ie);
         AliAODPWG4Particle * p2 = (AliAODPWG4Particle *)array2->At(je);
         Double_t f2gammaAsy;
         Get2GammaAsyDistPtM(p1,p2,f2gammaAsy,fDist,fPt,fMass);
         if(TMath::Abs(fMass-fPi0Mass)<fPi0MassPeakWidthCut){
            for(Int_t ke=0;ke<nclusters2;ke++){ //third photon
                AliAODPWG4Particle * p3 = (AliAODPWG4Particle *)array2->At(ke);
                if(p3 != p2){ //p3 != p2
                    TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
                    TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
                    Double_t px12 =(photon1+photon2).Px();
                    Double_t py12 =(photon1+photon2).Py();
                    Double_t pz12 =(photon1+photon2).Pz();
                    Double_t e12 =(photon1+photon2).E();
                    TLorentzVector photon12(px12, py12,pz12,e12);
                    TLorentzVector photon3(p3->Px(),p3->Py(),p3->Pz(),p3->E());
                    GetAsyDistPtM(photon12, photon3,fAsy,fDist,fPt,fMass);

                    for(Int_t ipid=0;ipid<fNPID;ipid++){ //pid
                        if(p1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                              p2->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                              p3->IsPIDOK(ipid,AliCaloPID::kPhoton)){
                           fMixBPi0GammaAsyPtM[ipid]->Fill(fAsy, fPt, fMass); //default module, with pid
                        }
                        for(Int_t jmod=0;jmod<fNMod;jmod++){ //mod
                            for(Int_t kasy=0;kasy<fNAsy;kasy++){ //asy
                                if(IsPhotonSelected(p1,ipid,fModCut[jmod]) &&
                                   IsPhotonSelected(p2,ipid,fModCut[jmod]) &&
                                   IsPhotonSelected(p3,ipid,fModCut[jmod])){ //photon selectin
                                   //fill the hist
                                   for(Int_t mhist=0;mhist<fNHistos;mhist++){ //hist
                                       Double_t pt1 =mhist *fPerPtBin;
                                       Double_t pt2 =(mhist+1)*fPerPtBin;
                                       Int_t index = ipid*fNMod*fNAsy*fNHistos+                                                                                        jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
                                       if(fMass<fInvMassCut && f2gammaAsy<fAsyCut[kasy] && 
                                          fPt<pt2 && fPt>pt1&& fAsy<fAsyCut[kasy] &&fAsy&&fDist&&fPt&&fMass)
                                       fMixOmegaB[index]->Fill(fMass);
                                  }//hist 
                                } //photon selection
                            }//asy
                        }//mod
                    }//pid
                }    
            }//ptc3          
         } //pi0 candidate
     } //mix B 

 }//ie

}

//______________________________________________________________________________ 
void AliAnaNeutralMeson::MixOmegaC(TClonesArray * array1, TClonesArray * array2, TClonesArray * array3)
{
  //omega background
 //three omega background, we classify it A, B and C
 // --A
 // (r1_event1+r2_event1)+r3_event2
 //
 // --B
 // (r1_event1+r2_event2)+r3_event2
 // 
 // --C
 // (r1_event1+r2_event2)+r3_event3
 Int_t nclusters1 =  array1->GetEntries();
 Int_t nclusters2 =  array2->GetEntries();
 Int_t nclusters3 =  array3->GetEntries();
 for(Int_t ie=0;ie<nclusters1;ie++){
     for(Int_t je=0;je<nclusters2;je++){
         AliAODPWG4Particle * p1 = (AliAODPWG4Particle *)array1->At(ie);
         AliAODPWG4Particle * p2 = (AliAODPWG4Particle *)array2->At(je);
         Double_t f2gammaAsy;
         Get2GammaAsyDistPtM(p1,p2,f2gammaAsy,fDist,fPt,fMass);
         if(TMath::Abs(fMass-fPi0Mass)<fPi0MassPeakWidthCut){
            for(Int_t ke=0;ke<nclusters3;ke++){ 
                AliAODPWG4Particle * p3 = (AliAODPWG4Particle *)array3->At(ke); //third photon
                TLorentzVector photon1(p1->Px(),p1->Py(),p1->Pz(),p1->E());
                TLorentzVector photon2(p2->Px(),p2->Py(),p2->Pz(),p2->E());
                Double_t px12 =(photon1+photon2).Px();
                Double_t py12 =(photon1+photon2).Py();
                Double_t pz12 =(photon1+photon2).Pz();
                Double_t e12 =(photon1+photon2).E();
                TLorentzVector photon12(px12, py12,pz12,e12);
                TLorentzVector photon3(p3->Px(),p3->Py(),p3->Pz(),p3->E());
                GetAsyDistPtM(photon12, photon3,fAsy,fDist,fPt,fMass);

                for(Int_t ipid=0;ipid<fNPID;ipid++){
                    if(p1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                            p2->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                            p3->IsPIDOK(ipid,AliCaloPID::kPhoton)){
                       fMixCPi0GammaAsyPtM[ipid]->Fill(fAsy, fPt, fMass); //default module, with pid
                    }
                    for(Int_t jmod=0;jmod<fNMod;jmod++){
                        for(Int_t kasy=0;kasy<fNAsy;kasy++){
                           if(IsPhotonSelected(p1,ipid,fModCut[jmod]) && 
                              IsPhotonSelected(p2,ipid,fModCut[jmod])&&
                              IsPhotonSelected(p3,ipid,fModCut[jmod])){ 
                               
                              //filling the histo
                               for(Int_t mhist=0;mhist<fNHistos;mhist++){
                                   Double_t pt1 =mhist *fPerPtBin;
                                   Double_t pt2 =(mhist+1)*fPerPtBin;
                                   Int_t index = ipid*fNMod*fNAsy*fNHistos+                                                                                        jmod*fNAsy*fNHistos+kasy*fNHistos+mhist;
                                   if(fMass<fInvMassCut && f2gammaAsy<fAsyCut[kasy] && 
                                      fPt<pt2 && fPt>pt1&& fAsy<fAsyCut[kasy] &&fAsy&&fDist&&fPt&&fMass)
                                   fMixOmegaC[index]->Fill(fMass);
                              }//hist
                            }//photon selection
                        }//asy
                    }//mod
                }//pid
            }//ke
        }//pi0 candidate
     } //je
 }//ie

} 
 
//______________________________________________________________________________
void AliAnaNeutralMeson::PhotonAcceptance(AliStack * stack)
{
 //photon Acceptance
 for(Int_t i=0 ; i<stack->GetNtrack(); i++){
     TParticle * prim = stack->Particle(i) ;
     Int_t pdg = prim->GetPdgCode() ;
//   Int_t ndau = prim->GetNDaughters();
//   Int_t iphot1=prim->GetFirstDaughter() ;
//   Int_t iphot2=prim->GetLastDaughter() ;
     if(pdg ==22 && stack->IsPhysicalPrimary(i)) {
          Double_t photonPt = prim->Pt() ;
          //printf("photon, pt %2.2f\n",photonPt);                           
          Double_t photonY  = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;   
          Double_t phi   = TMath::RadToDeg()*prim->Phi() ;
          if(TMath::Abs(photonY) < 0.5){  fhPrimPhotonPt->Fill(photonPt);  }
              fhPrimPhotonY  ->Fill(photonY) ;
              fhPrimPhotonPhi->Fill(phi) ;
              //Check if both photons hit Calorimeter
              Bool_t inacceptance = kFALSE;
#ifdef __PHOSGEO__
              Int_t mod1 ;
              Double_t x1,z1;
              if(fCalorimeter == "PHOS" && fPHOSGeo->ImpactOnEmc(prim,mod1,z1,x1) )   inacceptance = kTRUE;
              //printf("In REAL PHOS acceptance? %d\n",inacceptance);
#else
              TLorentzVector lv1;
              prim->Momentum(lv1);
              if(GetFidutialCut()->IsInFidutialCut(lv1,fCalorimeter) )    inacceptance = kTRUE ;
              //printf("In %s fFidutial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
#endif
              if(inacceptance){
                  fhPrimPhotonAccPt->Fill(photonPt) ;
                  fhPrimPhotonAccPhi->Fill(phi) ;
                  fhPrimPhotonAccY->Fill(photonY) ;
               }//Accepted
 
  }// Primary photon
 }//loop on primaries  
}
 
//______________________________________________________________________________
void AliAnaNeutralMeson::Pi0EtaAcceptance(AliStack * stack)
{
 //pi0 and eta Acceptance
 for(Int_t i=0 ; i<stack->GetNprimary(); i++){
     TParticle * prim = stack->Particle(i) ;
     Int_t pdg = prim->GetPdgCode() ;
     if(pdg==111 || pdg==221){
        Double_t pi0Pt = prim->Pt() ;
        //printf("  pi0, pt %2.2f\n",pi0Pt);
        Double_t pi0Y  = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
        Double_t phi   = TMath::RadToDeg()*prim->Phi();
        if(pdg==111){
           if(TMath::Abs(pi0Y) < 0.5)   fhPrimPi0Pt->Fill(pi0Pt) ;
              fhPrimPi0Y  ->Fill(pi0Y) ;
              fhPrimPi0Phi->Fill(phi) ;
         } 
        else{
           if(TMath::Abs(pi0Y) < 0.5)  fhPrimEtaPt->Fill(pi0Pt) ;
              fhPrimEtaY  ->Fill(pi0Y) ;
              fhPrimEtaPhi->Fill(phi) ;
        }
        //Check if both photons hit Calorimeter
        Int_t ndau = prim->GetNDaughters();
        Int_t iphot1=prim->GetFirstDaughter() ;
        Int_t iphot2=prim->GetLastDaughter() ;
        if(ndau==2 && ndau>1 && iphot1>-1 && iphot1<stack->GetNtrack() 
                   && iphot2>-1 && iphot2<stack->GetNtrack()){   
            TParticle * phot1 = stack->Particle(iphot1) ;
            TParticle * phot2 = stack->Particle(iphot2) ;
            if(phot1 && phot2 && phot1->GetPdgCode()==22 && phot2->GetPdgCode()==22){
                Bool_t inacceptance = kFALSE;                
#ifdef __PHOSGEO__
                 Int_t mod, mod1,mod2 ;        
                 Double_t x, z, x1,z1,x2,z2;
                 if(fCalorimeter == "PHOS" && fPHOSGeo->ImpactOnEmc(prim, mod, z,x)                                                         && fPHOSGeo->ImpactOnEmc(phot1,mod1,z1,x1) 
                          && fPHOSGeo->ImpactOnEmc(phot2,mod2,z2,x2))  inacceptance = kTRUE;
                                    
                 //printf("In REAL PHOS acceptance? %d\n",inacceptance);
#else 
                 TLorentzVector lv0, lv1, lv3;
                 prim->Momentum(lv0);
                 phot1->Momentum(lv1);
                 phot2->Momentum(lv3);
                 if(GetFidutialCut()->IsInFidutialCut(lv0,fCalorimeter) 
                    &&GetFidutialCut()->IsInFidutialCut(lv1,fCalorimeter) 
                    && GetFidutialCut()->IsInFidutialCut(lv3,fCalorimeter))  inacceptance = kTRUE ;
#endif
//                printf("In %s fFidutial cut acceptance? %d\n",fCalorimeter.Data(),inacceptance);
                 if(inacceptance){
                    if(pdg==111){
                        fhPrimPi0AccPt->Fill(pi0Pt) ;
                        fhPrimPi0AccPhi->Fill(phi) ;                                                                                      fhPrimPi0AccY->Fill(pi0Y) ;
                     }
                     else{
                        fhPrimEtaAccPt->Fill(pi0Pt) ;                                                                                     fhPrimEtaAccPhi->Fill(phi) ;
                        fhPrimEtaAccY->Fill(pi0Y) ;
                    }
                  }//Accepted
           }// 2 photons      
       }//Check daughters exist
    }// Primary pi0 
  }//loop on primaries    
}
 
//______________________________________________________________________________
void AliAnaNeutralMeson::OmegaAcceptance(AliStack * stack)
{
 //acceptance for omega->pi0+gamma       
 for(Int_t i=0 ; i<stack->GetNprimary(); i++){ //aaa
     TParticle * prim = stack->Particle(i) ;
     Int_t pdg = prim->GetPdgCode();
     if(pdg == 223){ //bbb
        Double_t pt = prim->Pt() ;
        Double_t y  = 0.5*TMath::Log((prim->Energy()-prim->Pz())/(prim->Energy()+prim->Pz())) ;
        Double_t phi   = TMath::RadToDeg()*prim->Phi() ;
        if(TMath::Abs(y) < 0.5){  fhPrimOmegaPt->Fill(pt) ;  }
        fhPrimOmegaY  ->Fill(y) ;  
        fhPrimOmegaPhi->Fill(phi) ;
        //check if omega->pi0+gamma->(gamma+gamma)+gamma hit Calorimeter
        Int_t ndau = prim->GetNDaughters();
        Int_t iphot1=prim->GetFirstDaughter() ;
        Int_t iphot2=prim->GetLastDaughter() ;
        if(ndau==2 && iphot1>-1 && iphot1<stack->GetNtrack() && iphot2>-1 && iphot2<stack->GetNtrack()){ //ccc
           TParticle * phot1 = stack->Particle(iphot1) ;
           TParticle * phot2 = stack->Particle(iphot2) ;
           Int_t pdg1 = phot1->GetPdgCode();
           Int_t pdg2 = phot2->GetPdgCode();        
  	   Int_t ndau2 = phot2->GetNDaughters();
           if(ndau2==2 && phot1 && pdg1==22 &&phot2 && pdg2==111){ //ddd
              Int_t jphot1 = phot2->GetFirstDaughter() ;
              Int_t jphot2 = phot2->GetLastDaughter() ;
              if(jphot1>-1 && jphot1<stack->GetNtrack() &&
                              jphot2>-1 && jphot2<stack->GetNtrack()){  //eee                
                  TParticle * phot3 = stack->Particle(jphot1) ;
                  TParticle * phot4 = stack->Particle(jphot2) ;
                  Int_t pdg3 = phot3->GetPdgCode();
                  Int_t pdg4 = phot4->GetPdgCode();
                  if(phot3 && pdg3==22 && phot4 && pdg4==22){ //fff
                     Bool_t inacceptance = kFALSE;
#ifdef __PHOSGEO__
                     Int_t mod, mod1,mod2,mod3 ;
                     Double_t x,z, x1,z1,x2,z2, x3,z3 ;
                     //3 gammas hit Calorimeter
                     if(fCalorimeter == "PHOS" && fPHOSGeo->ImpactOnEmc(prim,mod,z,x) 
                                               && fPHOSGeo->ImpactOnEmc(phot1,mod1,z1,x1) 
                                               && fPHOSGeo->ImpactOnEmc(phot3,mod2,z2,x2) 
                                               && fPHOSGeo->ImpactOnEmc(phot4,mod3,z3,x3))
                      inacceptance = kTRUE;
 
#else
                      TLorentzVector lv0,lv1, lv2, lv3;                                                         
                      prim->Momentum(lv0);
                      phot1->Momentum(lv1);
                      phot3->Momentum(lv2);
                      phot4->Momentum(lv3);
                      if(GetFidutialCut()->IsInFidutialCut(lv0,fCalorimeter) &&
                         GetFidutialCut()->IsInFidutialCut(lv1,fCalorimeter) && 
                         GetFidutialCut()->IsInFidutialCut(lv2,fCalorimeter) && 
                         GetFidutialCut()->IsInFidutialCut(lv3,fCalorimeter))  
                      inacceptance = kTRUE ;
#endif
                      if(inacceptance){ //ggg
                         fhPrimOmegaAccPt->Fill(pt) ;
                         fhPrimOmegaAccPhi->Fill(phi) ;
                         fhPrimOmegaAccY->Fill(y) ;
                       } //ggg
                   }//fff
                }//eee
             }//ddd
         }//ccc
       } //bbb
  }//aaa                                  
}                                                               
 
void AliAnaNeutralMeson::GetMCDistAsy(TParticle *p1,TParticle *p2,Double_t &dist,Double_t &asy)
{   
  //get dist and asy for MC particle                                                            
  Double_t eta1 = p1->Eta();        
  Double_t eta2 = p2->Eta(); 
  Double_t phi1 = p1->Phi();
  Double_t phi2 = p2->Phi();                                            
  Double_t e1 = p1->Energy(); 
  Double_t e2 = p2->Energy();                                                               
  dist = TMath::Sqrt( (eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
  if(e1+e2) asy = TMath::Abs(e1-e2)/(e1+e2);
}                                   

//______________________________________________________________________________
//void AliAnaNeutralMeson::Terminate(TList * outList) 
//{
// //Do some calculations and plots from the final histograms.
//  printf("AliAnaNeutralMeson::Terminate() - Not implemented yet\n");
//
//}


