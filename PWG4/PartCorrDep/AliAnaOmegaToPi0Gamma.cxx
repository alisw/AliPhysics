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
// class to extract omega(782)->pi0+gamma->3gamma
//
//-- Author: Renzhuo Wan (IOPP-Wuhan, China)
//_________________________________________________________________________

// --- ROOT system
class TROOT;

// --- AliRoot system
//class AliVEvent;
// --- ROOT system ---
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TFile.h"
//---- AliRoot system ----
#include "AliAnaOmegaToPi0Gamma.h"
#include "AliCaloTrackReader.h"
#include "AliCaloPID.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
ClassImp(AliAnaOmegaToPi0Gamma)

//______________________________________________________________________________
AliAnaOmegaToPi0Gamma::AliAnaOmegaToPi0Gamma() : AliAnaPartCorrBaseClass(),
fInputAODPi0(0), fInputAODGammaName(""),
fEventsList(0x0),fNVtxZBin(0), fNCentBin(0), fNRpBin(0), fNBadChDistBin(0), fNpid(0),
fNmaxMixEv(0), fVtxZCut(0), fCent(0), fRp(0), 
fPi0Mass(0),fPi0MassWindow(0),fPi0OverOmegaPtCut(0),
fGammaOverOmegaPtCut(0),
fhEtalon(0),
fRealOmega0(0), fMixAOmega0(0),
fMixBOmega0(0), fMixCOmega0(0),
fRealOmega1(0), fMixAOmega1(0),
fMixBOmega1(0), fMixCOmega1(0),
fRealOmega2(0), fMixAOmega2(0),
fMixBOmega2(0), fMixCOmega2(0),
fhOmegaPriPt(0)
{
 //Default Ctor
 InitParameters();
}
/*
//______________________________________________________________________________
AliAnaOmegaToPi0Gamma::AliAnaOmegaToPi0Gamma(const AliAnaOmegaToPi0Gamma & ex) : AliAnaPartCorrBaseClass(ex),
fInputAODPi0(new TClonesArray (*ex.fInputAODPi0)),
fInputAODGammaName(ex.fInputAODGammaName),
fEventsList(0x0), 
fNVtxZBin(ex.fNVtxZBin), fNCentBin(ex.fNCentBin), fNRpBin(ex.fNRpBin),
fNBadChDistBin(ex.fNBadChDistBin),fNpid(ex.fNpid),
fNmaxMixEv(ex.fNmaxMixEv),
fVtxZCut(ex.fVtxZCut), fCent(ex.fCent), fRp(ex.fRp), 
fPi0Mass(ex.fPi0Mass),
fPi0MassWindow(ex.fPi0MassWindow),
fPi0OverOmegaPtCut(ex.fPi0OverOmegaPtCut),
fGammaOverOmegaPtCut(ex.fGammaOverOmegaPtCut),
fhEtalon(ex.fhEtalon),
fRealOmega0(ex.fRealOmega0), fMixAOmega0(ex.fMixAOmega0),
fMixBOmega0(ex.fMixBOmega0), fMixCOmega0(ex.fMixCOmega0),
fRealOmega1(ex.fRealOmega1), fMixAOmega1(ex.fMixAOmega1),
fMixBOmega1(ex.fMixBOmega1), fMixCOmega1(ex.fMixCOmega1),
fRealOmega2(ex.fRealOmega2), fMixAOmega2(ex.fMixAOmega2),
fMixBOmega2(ex.fMixBOmega2), fMixCOmega2(ex.fMixCOmega2),
fhOmegaPriPt(ex.fhOmegaPriPt)
{
 // cpy ctor
 //Do not need it
}
*/
/*
//______________________________________________________________________________
AliAnaOmegaToPi0Gamma & AliAnaOmegaToPi0Gamma::operator = (const AliAnaOmegaToPi0Gamma & ex)
{
 // assignment operator
  
 if(this == &ex)return *this;
   ((AliAnaPartCorrBaseClass *)this)->operator=(ex);
   fInputAODPi0 = new TClonesArray(*ex.fInputAODPi0);
   fInputAODGammaName = ex.fInputAODGammaName;
   fEventsList = ex.fEventsList;

   fNVtxZBin=ex.fNVtxZBin;
   fNCentBin=ex.fNCentBin;
   fNRpBin=ex.fNRpBin;
   fNBadChDistBin=ex.fNBadChDistBin;
   fNpid=ex.fNpid;
   fNmaxMixEv =ex.fNmaxMixEv;

   fVtxZCut=ex.fVtxZCut;
   fCent=ex.fCent;
   fRp=ex.fRp;

   fPi0Mass=ex.fPi0Mass;
   fPi0MassWindow=ex.fPi0MassWindow;
   fPi0OverOmegaPtCut=ex.fPi0OverOmegaPtCut;
   fGammaOverOmegaPtCut=ex.fGammaOverOmegaPtCut;

   fhEtalon=ex.fhEtalon;
   fRealOmega0=ex.fRealOmega0;
   fMixAOmega0=ex.fMixAOmega0;
   fMixBOmega0=ex.fMixBOmega0;
   fMixCOmega0=ex.fMixCOmega0;
   fRealOmega1=ex.fRealOmega1;
   fMixAOmega1=ex.fMixAOmega1;
   fMixBOmega1=ex.fMixBOmega1;
   fMixCOmega1=ex.fMixCOmega1;
   fRealOmega2=ex.fRealOmega2;
   fMixAOmega2=ex.fMixAOmega2;
   fMixBOmega2=ex.fMixBOmega2;
   fMixCOmega2=ex.fMixCOmega2;
   fhOmegaPriPt=ex.fhOmegaPriPt;
  return *this;
	
}
*/

//______________________________________________________________________________
AliAnaOmegaToPi0Gamma::~AliAnaOmegaToPi0Gamma() {

  //dtor
//  Done by the maker  
//  if(fInputAODPi0){
//    fInputAODPi0->Clear();
//    delete fInputAODPi0;
//  }  

  if(fEventsList){
     for(Int_t i=0;i<fNVtxZBin;i++){
        for(Int_t j=0;j<fNCentBin;j++){
           for(Int_t k=0;k<fNRpBin;k++){
               fEventsList[i*fNCentBin*fNRpBin+j*fNRpBin+k]->Clear();
               delete fEventsList[i*fNCentBin*fNRpBin+j*fNRpBin+k];
           }
        }
     }
  }
  delete [] fEventsList;
  fEventsList=0;

 delete [] fVtxZCut;
 delete [] fCent;
 delete [] fRp;

}

//______________________________________________________________________________
void AliAnaOmegaToPi0Gamma::InitParameters()
{
//Init parameters when first called the analysis
//Set default parameters
 fInputAODGammaName="PhotonsDetector";  
 fNVtxZBin=1;              
 fNCentBin=1;               
 fNRpBin=1;                 
 fNBadChDistBin=3;          
 fNpid=9;                   
 fNmaxMixEv=8;              
 
 fPi0Mass=0.1348;             
 fPi0MassWindow=0.015;       
 fPi0OverOmegaPtCut=0.8;   
 fGammaOverOmegaPtCut=0.2; 
}


//______________________________________________________________________________
TList * AliAnaOmegaToPi0Gamma::GetCreateOutputObjects()
{
  //
  fVtxZCut = new Double_t [fNVtxZBin];
  for(Int_t i=0;i<fNVtxZBin;i++) fVtxZCut[i]=10*(i+1);

  fCent=new Double_t[fNCentBin];
  for(int i = 0;i<fNCentBin;i++)fCent[i]=0;

  fRp=new Double_t[fNRpBin];
  for(int i = 0;i<fNRpBin;i++)fRp[i]=0;
  //
  Int_t nptbins   = GetHistoPtBins();
  Float_t ptmax   = GetHistoPtMax();
  Float_t ptmin  = GetHistoPtMin();

  Int_t nmassbins = GetHistoMassBins();
  Float_t massmin = GetHistoMassMin();
  Float_t massmax = GetHistoMassMax();

  fhEtalon = new TH2F("hEtalon","Histo with binning parameters", nptbins,ptmin,ptmax,nmassbins,massmin,massmax) ;
  fhEtalon->SetXTitle("P_{T} (GeV)") ;
  fhEtalon->SetYTitle("m_{inv} (GeV)") ;
 
 // store them in fOutputContainer
 fEventsList = new TList*[fNVtxZBin*fNCentBin*fNRpBin];
 for(Int_t i=0;i<fNVtxZBin;i++){
     for(Int_t j=0;j<fNCentBin;j++){
         for(Int_t k=0;k<fNRpBin;k++){
             fEventsList[i*fNCentBin*fNRpBin+j*fNRpBin+k]=new TList();
         }
     }
 }
	
 TList * outputContainer = new TList() ; 
 outputContainer->SetName(GetName());
 char key[255] ;
 char title[255] ;
 const char * detector= fInputAODGammaName.Data();
 Int_t ndim=fNVtxZBin*fNCentBin*fNRpBin*fNBadChDistBin*fNpid;

 fRealOmega0 =new TH2F*[ndim];
 fMixAOmega0 =new TH2F*[ndim];
 fMixBOmega0 =new TH2F*[ndim];
 fMixCOmega0 =new TH2F*[ndim];
 
 fRealOmega1 =new TH2F*[ndim];
 fMixAOmega1 =new TH2F*[ndim];
 fMixBOmega1 =new TH2F*[ndim];
 fMixCOmega1 =new TH2F*[ndim];

 fRealOmega2 =new TH2F*[ndim];
 fMixAOmega2 =new TH2F*[ndim];
 fMixBOmega2 =new TH2F*[ndim];
 fMixCOmega2 =new TH2F*[ndim];

 for(Int_t i=0;i<fNVtxZBin;i++){
     for(Int_t j=0;j<fNCentBin;j++){
         for(Int_t k=0;k<fNRpBin;k++){ //at event level
             Int_t idim=i*fNCentBin*fNRpBin+j*fNRpBin+k;
             for(Int_t ipid=0;ipid<fNpid;ipid++){
		 for(Int_t idist=0;idist<fNBadChDistBin;idist++){ //at particle level

                     Int_t index=idim*fNpid*fNBadChDistBin+ipid*fNBadChDistBin+idist;

		     sprintf(key,"RealToPi0Gamma_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s Real Pi0GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fRealOmega0[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fRealOmega0[index]->SetName(key) ;
                     fRealOmega0[index]->SetTitle(title);
                     outputContainer->Add(fRealOmega0[index]);

		     sprintf(key,"MixAToPi0Gamma_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixA Pi0GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixAOmega0[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixAOmega0[index]->SetName(key) ;
                     fMixAOmega0[index]->SetTitle(title);
                     outputContainer->Add(fMixAOmega0[index]);

		     sprintf(key,"MixBToPi0Gamma_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixB Pi0GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixBOmega0[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixBOmega0[index]->SetName(key) ;
                     fMixBOmega0[index]->SetTitle(title);
                     outputContainer->Add(fMixBOmega0[index]);

		     sprintf(key,"MixCToPi0Gamma_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixC Pi0GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixCOmega0[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixCOmega0[index]->SetName(key) ;
                     fMixCOmega0[index]->SetTitle(title);
                     outputContainer->Add(fMixCOmega0[index]);

		     sprintf(key,"RealToPi0Gamma1_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s Real Pi0(A<0.7)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fRealOmega1[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fRealOmega1[index]->SetName(key) ;
                     fRealOmega1[index]->SetTitle(title);
                     outputContainer->Add(fRealOmega1[index]);

                     sprintf(key,"MixAToPi0Gamma1_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixA Pi0(A<0.7)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixAOmega1[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixAOmega1[index]->SetName(key) ;
                     fMixAOmega1[index]->SetTitle(title);
                     outputContainer->Add(fMixAOmega1[index]);

                     sprintf(key,"MixBToPi0Gamma1_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixB Pi0(A<0.7)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixBOmega1[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixBOmega1[index]->SetName(key) ;
                     fMixBOmega1[index]->SetTitle(title);
                     outputContainer->Add(fMixBOmega1[index]);

                     sprintf(key,"MixCToPi0Gamma1_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixC Pi0(A<0.7)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixCOmega1[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixCOmega1[index]->SetName(key) ;
                     fMixCOmega1[index]->SetTitle(title);
                     outputContainer->Add(fMixCOmega1[index]);

		     sprintf(key,"RealToPi0Gamma2_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s Real Pi0(A<0.8)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fRealOmega2[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fRealOmega2[index]->SetName(key) ;
                     fRealOmega2[index]->SetTitle(title);
                     outputContainer->Add(fRealOmega2[index]);

                     sprintf(key,"MixAToPi0Gamma2_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixA Pi0(A<0.8)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixAOmega2[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixAOmega2[index]->SetName(key) ;
                     fMixAOmega2[index]->SetTitle(title);
                     outputContainer->Add(fMixAOmega2[index]);

                     sprintf(key,"MixBToPi0Gamma2_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixB Pi0(A<0.8)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixBOmega2[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixBOmega2[index]->SetName(key) ;
                     fMixBOmega2[index]->SetTitle(title);
                     outputContainer->Add(fMixBOmega2[index]);

                     sprintf(key,"MixCToPi0Gamma2_Vz%dC%dRp%dPid%dDist%d",i,j,k,ipid,idist);
                     sprintf(title, "%s MixC Pi0(A<0.8)GammaIVM vz_%2.1f_ct_%2.1f_Rp_%2.1f_pid_%d_dist_%d",detector,fVtxZCut[i],fCent[j],fRp[k],ipid,idist);
                     fhEtalon->Clone(key);
                     fMixCOmega2[index]=(TH2F*)fhEtalon->Clone(key) ;
                     fMixCOmega2[index]->SetName(key) ;
                     fMixCOmega2[index]->SetTitle(title);
                     outputContainer->Add(fMixCOmega2[index]);
                 }
             }
         }
     }  
 }

 if(IsDataMC()){
    sprintf(key, "%sOmegaPri",detector);
    sprintf(title,"primary #omega in %s",detector);
    fhOmegaPriPt=new TH1F(key, title,nptbins,ptmin,ptmax);
    fhOmegaPriPt->GetXaxis()->SetTitle("P_{T}");
    fhOmegaPriPt->GetYaxis()->SetTitle("dN/P_{T}");
    outputContainer->Add(fhOmegaPriPt);
 }

 delete fhEtalon;
 return outputContainer;
}

//______________________________________________________________________________
void AliAnaOmegaToPi0Gamma::Print(const Option_t * /*opt*/) const
{
  //Print some relevant parameters set in the analysis
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");
  printf("Omega->pi0+gamma->3gamma\n");
  printf("Cuts at event level:            \n");
  printf("Bins of vertex Z:                     %d \n", fNVtxZBin);
  printf("Bins of centrality:                   %d \n",fNCentBin);
  printf("Bins of Reaction plane:               %d\n",fNRpBin);
  printf("Cuts at AOD particle level:\n");
  printf("Number of PID:                        %d \n", fNpid);
  printf("Number of DistToBadChannel cuts:      %d\n", fNBadChDistBin);
  printf("number of events buffer to be mixed:  %d\n",fNmaxMixEv);
} 

//______________________________________________________________________________
void AliAnaOmegaToPi0Gamma::MakeAnalysisFillHistograms() 
{
 //fill the MC AOD if needed first
 //-----------
 //need to be further implemented
 AliStack * stack = 0x0;
// TParticle * primary = 0x0;
 TClonesArray * mcparticles0 = 0x0;
 TClonesArray * mcparticles1 = 0x0;
 AliAODMCParticle * aodprimary = 0x0;
 Int_t pdg=0;
 Double_t pt=0;
 Double_t eta=0;

 if(IsDataMC()){
     if(GetReader()->ReadStack()){
        stack =  GetMCStack() ;
        if(!stack)
        printf("AliAnaAcceptance::MakeAnalysisFillHistograms() - There is no stack!\n");
        for(Int_t i=0 ; i<stack->GetNtrack(); i++){
            TParticle * prim = stack->Particle(i) ;
            pdg = prim->GetPdgCode() ;
            eta=prim->Eta();
            pt=prim->Pt();
            if(TMath::Abs(eta)<0.5) {
               if(pdg==223) fhOmegaPriPt->Fill(pt);
            }
        }
      }
      else if(GetReader()->ReadAODMCParticles()){
          //Get the list of MC particles
          mcparticles0 = GetReader()->GetAODMCParticles(0);
          if(!mcparticles0 && GetDebug() > 0)     {
             printf("AliAnaAcceptance::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
           }
           if(GetReader()->GetSecondInputAODTree()){
               mcparticles1 = GetReader()->GetAODMCParticles(1);
               if(!mcparticles1 && GetDebug() > 0)     {
                   printf("AliAnaAcceptance::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
                }
           }
           for(Int_t i=0;i<mcparticles0->GetEntries();i++){
               aodprimary =(AliAODMCParticle*)mcparticles0->At(i);
               pdg = aodprimary->GetPdgCode() ;
               eta=aodprimary->Eta();
               pt=aodprimary->Pt();
               if(TMath::Abs(eta)<0.5) {
                  if(pdg==223) fhOmegaPriPt->Fill(pt);
                }

           }
     }
 }// is data and MC
  

 //process event from AOD brach 
 //extract pi0, eta and omega analysis
 Int_t iRun=(GetReader()->GetInputEvent())->GetRunNumber() ;
 if(IsBadRun(iRun)) return ;	

 //vertex z
 Double_t vert[]={0,0,0} ;
 GetReader()->GetVertex(vert);
 Int_t curEventBin =0;

 Int_t ivtxzbin=(Int_t)TMath::Abs(vert[2])/10;
 if(ivtxzbin>=fNVtxZBin)return;

 //centrality
 Int_t icentbin=0;

 //reaction plane
 Int_t irpbin=0;

 if(ivtxzbin==-1) return; 
 curEventBin = ivtxzbin*fNCentBin*fNRpBin + icentbin*fNRpBin + irpbin;
  TClonesArray *aodGamma = (TClonesArray*) GetAODBranch(fInputAODGammaName); //photon array
//  TClonesArray *aodGamma = (TClonesArray *) GetReader()->GetOutputEvent()->FindListObject(fInputAODGammaName); //photon array
 Int_t nphotons =0;
 if(aodGamma) nphotons= aodGamma->GetEntries(); 
 else return;

 fInputAODPi0 = (TClonesArray*)GetInputAODBranch();  //pi0 array
 Int_t npi0s = 0;
 if(fInputAODPi0) npi0s= fInputAODPi0 ->GetEntries();
 else return;

 if(nphotons<3 || npi0s<1)return; //for pi0, eta and omega->pi0+gamma->3gamma reconstruction

 //reconstruction of omega(782)->pi0+gamma->3gamma
 //loop for pi0 and photon
 if(GetDebug() > 0) printf("omega->pi0+gamma->3gamma invariant mass analysis ! This event have %d photons and %d pi0 \n", nphotons, npi0s);
 for(Int_t i=0;i<npi0s;i++){
     AliAODPWG4Particle * pi0 = (AliAODPWG4Particle*) (fInputAODPi0->At(i)) ; //pi0
     TLorentzVector vpi0(pi0->Px(),pi0->Py(),pi0->Pz(),pi0->E());
     Int_t lab1=pi0->GetCaloLabel(0);  // photon1 from pi0 decay
     Int_t lab2=pi0->GetCaloLabel(1);  // photon2 from pi0 decay
     //for omega->pi0+gamma, it needs at least three photons per event
     //Get the two decay photons from pi0
     AliAODPWG4Particle * photon1 =0;
     AliAODPWG4Particle * photon2 =0;
     for(Int_t d1=0;d1<nphotons;d1++){
         for(Int_t d2=0;d2<nphotons;d2++){
             AliAODPWG4Particle * dp1 = (AliAODPWG4Particle*) (aodGamma->At(d1));
             AliAODPWG4Particle * dp2 = (AliAODPWG4Particle*) (aodGamma->At(d2));
             Int_t dlab1=dp1->GetCaloLabel(0);
             Int_t dlab2=dp2->GetCaloLabel(0);
             if(dlab1==lab1 && dlab2==lab2){
                photon1=dp1;
                photon2=dp2;
             }
             else continue;
         }
    }
    //caculate the asy and dist of the two photon from pi0 decay
    TLorentzVector dph1(photon1->Px(),photon1->Py(),photon1->Pz(),photon1->E());
    TLorentzVector dph2(photon2->Px(),photon2->Py(),photon2->Pz(),photon2->E());

    Double_t pi0asy= TMath::Abs(dph1.E()-dph2.E())/(dph1.E()+dph2.E());
//    Double_t phi1=dph1.Phi();
//    Double_t phi2=dph2.Phi();
//    Double_t eta1=dph1.Eta();
//    Double_t eta2=dph2.Eta();
//    Double_t pi0dist=TMath::Sqrt((phi1-phi2)*(phi1-phi2)+(eta1-eta2)*(eta1-eta2));

    if(pi0->GetPdg()==111  && nphotons>2 && npi0s
           && TMath::Abs(vpi0.M()-fPi0Mass)<fPi0MassWindow) { //pi0 candidates

         //avoid the double counting
         Int_t * dc1= new Int_t[nphotons];
         Int_t * dc2= new Int_t[nphotons];
         Int_t index1=0;
         Int_t index2=0;
         for(Int_t k=0;k<i;k++){
             AliAODPWG4Particle * p3=(AliAODPWG4Particle*)(fInputAODPi0->At(k));
             Int_t lab4=p3->GetCaloLabel(0);
             Int_t lab5=p3->GetCaloLabel(1);
             if(lab1==lab4){ dc1[index1]=lab5;  index1++;  }
             if(lab2==lab5){ dc2[index2]=lab4;  index2++;  }
         }


	 //loop the pi0 with third gamma
  	 for(Int_t j=0;j<nphotons;j++){
             AliAODPWG4Particle *photon3 = (AliAODPWG4Particle*) (aodGamma->At(j));
	     TLorentzVector dph3(photon3->Px(),photon3->Py(),photon3->Pz(),photon3->E());
	     Int_t lab3=photon3->GetCaloLabel(0);
             Double_t pi0gammapt=(vpi0+dph3).Pt();
             Double_t pi0gammamass=(vpi0+dph3).M();
             Double_t pi0OverOmegaPtRatio =vpi0.Pt()/pi0gammapt; 
             Double_t gammaOverOmegaPtRatio= dph3.Pt()/pi0gammapt;

             //pi0, gamma pt cut             
             if(pi0OverOmegaPtRatio>fPi0OverOmegaPtCut || 
                gammaOverOmegaPtRatio<fGammaOverOmegaPtCut) continue;

    	     for(Int_t l=0;l<index1;l++) if(lab3==dc1[l]) lab3=-1;
             for(Int_t l=0;l<index2;l++) if(lab3==dc2[l]) lab3=-1;
    
	     if(lab3>0 && lab3!=lab1 && lab3!=lab2){
	        for(Int_t ipid=0;ipid<fNpid;ipid++){
		    for(Int_t idist=0;idist<fNBadChDistBin;idist++){
			Int_t index=curEventBin*fNpid*fNBadChDistBin+ipid*fNBadChDistBin+idist;
	                if(photon1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
		    	   photon2->IsPIDOK(ipid,AliCaloPID::kPhoton) && 
			   photon3->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
			   photon1->DistToBad()>=idist &&
			   photon2->DistToBad()>=idist &&
			   photon3->DistToBad()>=idist ){
			   //fill the histograms
                           if(GetDebug() > 2) printf("Real: index  %d  pt  %2.3f  mass   %2.3f \n", index, pi0gammapt, pi0gammamass);
			   fRealOmega0[index]->Fill(pi0gammapt,pi0gammamass);
			   if(pi0asy<0.7) fRealOmega1[index]->Fill(pi0gammapt,pi0gammamass);
			   if(pi0asy<0.8) fRealOmega2[index]->Fill(pi0gammapt,pi0gammamass);
 			}
	             }
	        }
             }
         }	
         delete []dc1;
	 delete []dc2;
 if(GetDebug() > 0) printf("MixA: (r1_event1+r2_event1)+r3_event2 \n");
	 //-------------------------
	 //background analysis
	 //three background
	 // --A   (r1_event1+r2_event1)+r3_event2
         Int_t nMixed = fEventsList[curEventBin]->GetSize();
         for(Int_t im=0;im<nMixed;im++){
             TClonesArray* ev2= (TClonesArray*) (fEventsList[curEventBin]->At(im));
	     for(Int_t mix1=0;mix1<ev2->GetEntries();mix1++){
                AliAODPWG4Particle *mix1ph = (AliAODPWG4Particle*) (ev2->At(mix1));     
		TLorentzVector vmixph(mix1ph->Px(),mix1ph->Py(),mix1ph->Pz(),mix1ph->E());
                Double_t pi0gammapt=(vpi0+vmixph).Pt();
                Double_t pi0gammamass=(vpi0+vmixph).M();
                Double_t pi0OverOmegaPtRatio =vpi0.Pt()/pi0gammapt;
                Double_t gammaOverOmegaPtRatio= vmixph.Pt()/pi0gammapt;

              //pi0, gamma pt cut             
               if(pi0OverOmegaPtRatio>fPi0OverOmegaPtCut || 
                  gammaOverOmegaPtRatio<fGammaOverOmegaPtCut) continue;

	        for(Int_t ipid=0;ipid<fNpid;ipid++){
		   for(Int_t idist=0;idist<fNBadChDistBin;idist++){
		       Int_t index=curEventBin*fNpid*fNBadChDistBin+ipid*fNBadChDistBin+idist;
                       if(photon1->IsPIDOK(ipid,AliCaloPID::kPhoton)&&
                          photon2->IsPIDOK(ipid,AliCaloPID::kPhoton)&&
                          mix1ph->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                          photon1->DistToBad()>=idist &&
                          photon2->DistToBad()>=idist &&
                          mix1ph->DistToBad()>=idist ){
                          if(GetDebug() > 2) printf("MixA: index  %d   pt  %2.3f  mass   %2.3f \n",index, pi0gammapt, pi0gammamass);
		 	  //fill the histograms
			  fMixAOmega0[index]->Fill(pi0gammapt,pi0gammamass);
                          if(pi0asy<0.7)fMixAOmega1[index]->Fill(pi0gammapt,pi0gammamass);
                          if(pi0asy<0.8)fMixAOmega2[index]->Fill(pi0gammapt,pi0gammamass);
//printf("mix A  %d  %2.2f \n", index, pi0gammamass);
			  
                       }
		   }
                }
	     }
         }
     }
 }

 //
 // --B   (r1_event1+r2_event2)+r3_event2
 //
 if(GetDebug() >0)printf("MixB:  (r1_event1+r2_event2)+r3_event2 \n");
 for(Int_t i=0;i<nphotons;i++){
     AliAODPWG4Particle *ph1 = (AliAODPWG4Particle*) (aodGamma->At(i)); 
     TLorentzVector vph1(ph1->Px(),ph1->Py(),ph1->Pz(),ph1->E());

     Int_t nMixed = fEventsList[curEventBin]->GetSize();
     for(Int_t ie=0;ie<nMixed;ie++){
         TClonesArray* ev2= (TClonesArray*) (fEventsList[curEventBin]->At(ie));
         for(Int_t mix1=0;mix1<ev2->GetEntries();mix1++){
             AliAODPWG4Particle *ph2 = (AliAODPWG4Particle*) (ev2->At(mix1));
	     TLorentzVector vph2(ph2->Px(),ph2->Py(),ph2->Pz(),ph2->E());
	     Double_t pi0asy = TMath::Abs(vph1.E()-vph2.E())/(vph1.E()+vph2.E()); 	     
	     Double_t pi0mass=(vph1+vph2).M();
		
	     if(TMath::Abs(pi0mass-fPi0Mass)<fPi0MassWindow){//for pi0 selection
		for(Int_t mix2=(mix1+1);mix2<ev2->GetEntries();mix2++){
                    AliAODPWG4Particle *ph3 = (AliAODPWG4Particle*) (ev2->At(mix2));
                    TLorentzVector vph3(ph3->Px(),ph3->Py(),ph3->Pz(),ph3->E());
		    
                    Double_t pi0gammapt=(vph1+vph2+vph3).Pt();
                    Double_t pi0gammamass=(vph1+vph2+vph3).M(); 
                    Double_t pi0OverOmegaPtRatio =(vph1+vph2).Pt()/pi0gammapt;
                    Double_t gammaOverOmegaPtRatio= vph3.Pt()/pi0gammapt;
                   //pi0, gamma pt cut             
                   if(pi0OverOmegaPtRatio>fPi0OverOmegaPtCut ||
                       gammaOverOmegaPtRatio<fGammaOverOmegaPtCut) continue;

		    for(Int_t ipid=0;ipid<fNpid;ipid++){
			 for(Int_t idist=0;idist<fNBadChDistBin;idist++){
                            Int_t index=curEventBin*fNpid*fNBadChDistBin+ipid*fNBadChDistBin+idist;
  	                    if(ph1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
		               ph2->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                               ph3->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
			       ph1->DistToBad()>=idist &&
                               ph2->DistToBad()>=idist &&
                               ph3->DistToBad()>=idist ){
                               if(GetDebug() > 2) printf("MixB: index  %d   pt  %2.3f  mass   %2.3f \n", index, pi0gammapt, pi0gammamass);
			       //fill histograms
			       fMixBOmega0[index]->Fill(pi0gammapt,pi0gammamass);
                               if(pi0asy<0.7) fMixBOmega1[index]->Fill(pi0gammapt,pi0gammamass);
                               if(pi0asy<0.8) fMixBOmega2[index]->Fill(pi0gammapt,pi0gammamass);
//printf("mix B  %d  %2.2f \n", index, pi0gammamass);
			    }
		        }		    
		    }
		}
	
		//
	        // --C   (r1_event1+r2_event2)+r3_event3
                //
                if(GetDebug() >0)printf("MixC: (r1_event1+r2_event2)+r3_event3\n");
		for(Int_t je=(ie+1);je<nMixed;je++){
		    TClonesArray* ev3= (TClonesArray*) (fEventsList[curEventBin]->At(je));
		    for(Int_t mix3=0;mix3<ev3->GetEntries();mix3++){
	                AliAODPWG4Particle *ph3 = (AliAODPWG4Particle*) (ev3->At(mix3));
            		TLorentzVector vph3(ph3->Px(),ph3->Py(),ph3->Pz(),ph3->E());
			
                        Double_t pi0gammapt=(vph1+vph2+vph3).Pt();
	                Double_t pi0gammamass=(vph1+vph2+vph3).M();
		        Double_t pi0OverOmegaPtRatio =(vph1+vph2).Pt()/pi0gammapt;
                        Double_t gammaOverOmegaPtRatio= vph3.Pt()/pi0gammapt;
                        //pi0, gamma pt cut             
                        if(pi0OverOmegaPtRatio>fPi0OverOmegaPtCut ||
                          gammaOverOmegaPtRatio<fGammaOverOmegaPtCut) continue;
	
			for(Int_t ipid=0;ipid<fNpid;ipid++){
			    for(Int_t idist=0;idist<fNBadChDistBin;idist++){
                                Int_t index=curEventBin*fNpid*fNBadChDistBin+ipid*fNBadChDistBin+idist;
                                if(ph1->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                                   ph2->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
                                   ph3->IsPIDOK(ipid,AliCaloPID::kPhoton) &&
				   ph1->DistToBad()>=idist &&
                                   ph2->DistToBad()>=idist &&
                                   ph3->DistToBad()>=idist ){
                                   if(GetDebug() > 2) printf("MixC: index  %d  pt  %2.3f  mass   %2.3f \n", index, pi0gammapt, pi0gammamass);
                  	           //fill histograms
		  	          fMixCOmega0[index]->Fill(pi0gammapt,pi0gammamass);
                                  if(pi0asy<0.7) fMixCOmega1[index]->Fill(pi0gammapt,pi0gammamass);
				  if(pi0asy<0.8) fMixCOmega2[index]->Fill(pi0gammapt,pi0gammamass);
//printf("mix C  %d  %2.2f \n", index, pi0gammamass);
                               }
			   }
                        }
		    }
		}
	     } //for pi0 selecton		
	 }
     }
 }


 //event buffer 
 TClonesArray *currentEvent = new TClonesArray(*aodGamma);
 if(currentEvent->GetEntriesFast()>0){
        fEventsList[curEventBin]->AddFirst(currentEvent) ;
        currentEvent=0 ; 
        if(fEventsList[curEventBin]->GetSize()>=fNmaxMixEv) {
            TClonesArray * tmp = (TClonesArray*) (fEventsList[curEventBin]->Last()) ;
            fEventsList[curEventBin]->RemoveLast() ;
            delete tmp ;
        }
 }
 else{ 
       delete currentEvent ;
       currentEvent=0 ;
 }
 
}

//______________________________________________________________________________
void AliAnaOmegaToPi0Gamma::ReadHistograms(TList * outputList)
{
 //read the histograms 
 //for the finalization of the terminate analysis

 Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"RealToPi0Gamma_Vz0C0Rp0Pid0Dist0"));

  Int_t ndim=fNVtxZBin*fNCentBin*fNRpBin*fNBadChDistBin*fNpid;

 if(!fRealOmega0) fRealOmega0 =new TH2F*[ndim];
 if(!fMixAOmega0) fMixAOmega0 =new TH2F*[ndim];
 if(!fMixBOmega0) fMixBOmega0 =new TH2F*[ndim];
 if(!fMixCOmega0) fMixCOmega0 =new TH2F*[ndim];

 if(!fRealOmega1) fRealOmega1 =new TH2F*[ndim];
 if(!fMixAOmega1) fMixAOmega1 =new TH2F*[ndim];
 if(!fMixBOmega1) fMixBOmega1 =new TH2F*[ndim];
 if(!fMixCOmega1) fMixCOmega1 =new TH2F*[ndim];

 if(!fRealOmega2) fRealOmega2 =new TH2F*[ndim];
 if(!fMixAOmega2) fMixAOmega2 =new TH2F*[ndim];
 if(!fMixBOmega2) fMixBOmega2 =new TH2F*[ndim];
 if(!fMixCOmega2) fMixCOmega2 =new TH2F*[ndim];

  for(Int_t i=0;i<fNVtxZBin;i++){
     for(Int_t j=0;j<fNCentBin;j++){
         for(Int_t k=0;k<fNRpBin;k++){ //at event level
             Int_t idim=i*fNCentBin*fNRpBin+j*fNRpBin+k;
             for(Int_t ipid=0;ipid<fNpid;ipid++){ 
                for(Int_t idist=0;idist<fNBadChDistBin;idist++){ //at particle
                    Int_t ind=idim*fNpid*fNBadChDistBin+ipid*fNBadChDistBin+idist;
                    fRealOmega0[ind]= (TH2F*) outputList->At(index++);
                    fMixAOmega0[ind]= (TH2F*) outputList->At(index++);
                    fMixBOmega0[ind]= (TH2F*) outputList->At(index++);
                    fMixCOmega0[ind]= (TH2F*) outputList->At(index++);

                    fRealOmega1[ind]= (TH2F*) outputList->At(index++);
                    fMixAOmega1[ind]= (TH2F*) outputList->At(index++);
                    fMixBOmega1[ind]= (TH2F*) outputList->At(index++);
                    fMixCOmega1[ind]= (TH2F*) outputList->At(index++);

                    fRealOmega2[ind]= (TH2F*) outputList->At(index++);
                    fMixAOmega2[ind]= (TH2F*) outputList->At(index++);
                    fMixBOmega2[ind]= (TH2F*) outputList->At(index++);
                    fMixCOmega2[ind]= (TH2F*) outputList->At(index++);
                    
                 
                }
              }
          }
      }
  }
  
  if(IsDataMC()){
     fhOmegaPriPt  = (TH1F*)  outputList->At(index++);
  }

}

//______________________________________________________________________________
void AliAnaOmegaToPi0Gamma::Terminate(TList * outputList) 
{
// //Do some calculations and plots from the final histograms.
  if(GetDebug() >= 0) printf("AliAnaOmegaToPi0Gamma::Terminate() \n");
  ReadHistograms(outputList);
  char cvs1[128];  
  sprintf(cvs1, "Neutral_%s_IVM",fInputAODGammaName.Data());

  TCanvas * cvsIVM = new TCanvas(cvs1, cvs1, 400, 10, 600, 700) ;
  cvsIVM->Divide(2, 2);

  cvsIVM->cd(1);
  char dec[128];
  sprintf(dec,"h2Real_%s",fInputAODGammaName.Data());
  TH2F * h2Real= (TH2F*)fRealOmega0[0]->Clone(dec);
  h2Real->GetXaxis()->SetRangeUser(4,6);
  TH1F * hRealOmega = (TH1F*) h2Real->ProjectionY();
  hRealOmega->SetTitle("RealPi0Gamma 4<pt<6");
  hRealOmega->SetLineColor(2);
  hRealOmega->Draw();

  cvsIVM->cd(2);
  sprintf(dec,"hMixA_%s",fInputAODGammaName.Data());
  TH2F *h2MixA= (TH2F*)fMixAOmega0[0]->Clone(dec);
  h2MixA->GetXaxis()->SetRangeUser(4,6);
  TH1F * hMixAOmega = (TH1F*) h2MixA->ProjectionY();
  hMixAOmega->SetTitle("MixA 4<pt<6");
  hMixAOmega->SetLineColor(2);
  hMixAOmega->Draw();

  cvsIVM->cd(3);
  sprintf(dec,"hMixB_%s",fInputAODGammaName.Data());
  TH2F * h2MixB= (TH2F*)fMixBOmega0[0]->Clone(dec);
  h2MixB->GetXaxis()->SetRangeUser(4,6);
  TH1F * hMixBOmega = (TH1F*) h2MixB->ProjectionY();
  hMixBOmega->SetTitle("MixB 4<pt<6");
  hMixBOmega->SetLineColor(2);
  hMixBOmega->Draw();

  cvsIVM->cd(4);
  sprintf(dec,"hMixC_%s",fInputAODGammaName.Data());
  TH2F *h2MixC= (TH2F*)fMixCOmega0[0]->Clone(dec);
  h2MixC->GetXaxis()->SetRangeUser(4,6);
  TH1F * hMixCOmega = (TH1F*) h2MixC->ProjectionY();
  hMixCOmega->SetTitle("MixC 4<pt<6");
  hMixCOmega->SetLineColor(2);
  hMixCOmega->Draw();

  char eps[128 ];
  sprintf(eps,"CVS_%s_IVM.eps",fInputAODGammaName.Data());
  cvsIVM->Print(eps);
  cvsIVM->Modified();
 
}
