#define AliFlowAnalysisWithCumulants_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"
#include "TFile.h"
#include "TParticle.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"

class TH1;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TRandom3;
class TObjArray;
class TList;
class TCanvas;
class TSystem;
class TROOT;

class AliFlowVector;
class TVector;

//******************************* 
// flow analysis with cumulants *   
// author: Ante Bilandzic       * 
// email: anteb@nikhef.nl       *
//******************************* 

ClassImp(AliFlowAnalysisWithCumulants)

//________________________________________________________________________

AliFlowAnalysisWithCumulants::AliFlowAnalysisWithCumulants():  
  fTrack(NULL),
  fAvM(0),
  fR0(0),
  fPtMax(0),
  fPtMin(0),
  fBinWidth(0),
  fAvQx(0),
  fAvQy(0),
  fAvQ2x(0),
  fAvQ2y(0),
  fCommonHists(NULL),
  fCommonHistsRes2(NULL),
  fCommonHistsRes4(NULL),
  fCommonHistsRes6(NULL),
  fCommonHistsRes8(NULL)
  {
   //constructor 
   fR0=AliFlowCumuConstants::fgR0;
   fPtMax=AliFlowCommonConstants::GetPtMax(); 
   fPtMin=AliFlowCommonConstants::GetPtMin();
   fBinWidth=(fPtMax-fPtMin)/fgknBins;
  
   for(Int_t n=0;n<fgknBins;n++){
    fBinEventEntries[n]=0;
    fBinNoOfParticles[n]=0;
    fBinMeanPt[n]=0;
    for(Int_t p=0;p<fgkPmax;p++){
     for(Int_t q=0;q<fgkQmax;q++){
      fAvG[p][q]=0;
      fBinEventDRe[n][p][q]=0; 
      fBinEventDIm[n][p][q]=0;
     }
    }
   }
  }


//___________________________________________________________________________
void AliFlowAnalysisWithCumulants::CreateOutputObjects(){
 //output histograms
 TString fHistFileName = "cumulants.root";
 TFile* fHistFile;
 fHistFile = new TFile(fHistFileName.Data(),"RECREATE");
 fCommonHists = new AliFlowCommonHist("Cumulants");//control histograms
 fCommonHistsRes2 = new AliFlowCommonHistResults("Cumulants2");
 fCommonHistsRes4 = new AliFlowCommonHistResults("Cumulants4");
 fCommonHistsRes6 = new AliFlowCommonHistResults("Cumulants6");
 fCommonHistsRes8 = new AliFlowCommonHistResults("Cumulants8");
}

//________________________________________________________________________
void AliFlowAnalysisWithCumulants::Exec(AliFlowEventSimple* anEvent) {
  //running over data
 
  fCommonHists->FillControlHistograms(anEvent);   
  
  Double_t fG[fgkPmax][fgkQmax];//generating function for integrated flow
  for(Int_t p=0;p<fgkPmax;p++){
   for(Int_t q=0;q<fgkQmax;q++){
    fG[p][q]=1.; 
   }
  }
  
  //---------------------------------------------------------
  //Q vector stuff 
  AliFlowVector fQVector;
  fQVector.Set(0.,0.);
  fQVector.SetMult(0);
  
  fQVector=anEvent->GetQ();//get the Q vector for this event
  
  fAvQx+=fQVector.X();
  fAvQy+=fQVector.Y();
  fAvQ2x+=pow(fQVector.X(),2.);
  fAvQ2y+=pow(fQVector.Y(),2.);
  //----------------------------------------------------------
    
  Int_t nPrim = anEvent->NumberOfTracks();
  Int_t fEventNSelTracksIntFlow = anEvent->GetEventNSelTracksIntFlow();
  Int_t fSelTracksIntFlow = 0;
    
  cout<<"Number of input tracks for cumulant analysis: "<<nPrim<<endl;
  cout<<"Number of selected tracks for cumulant analysis: "<<fEventNSelTracksIntFlow<<endl;
  
    //------------------------------------------------------------------------------------
    //STARTING THE FIRST LOOP (CALCULATING THE GENERATING FUNCTION FOR INTEGRATED FLOW)
    for(Int_t i=0;i<nPrim;i++){
    fTrack=anEvent->GetTrack(i);
     if(fTrack&&fTrack->UseForIntegratedFlow()){
      fSelTracksIntFlow++;
      for(Int_t p=0;p<fgkPmax;p++){
       for(Int_t q=0;q<fgkQmax;q++){
	fG[p][q]*=(1.+(2.*fR0*sqrt(p+1.)/fEventNSelTracksIntFlow)*cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)); 
       }
      }
     }
    }
    // ENDING THE FIRST LOOP OVER TRACKS
    //------------------------------------------------------------------------------------
          
  
  //------------------------------------------------------------------------------------
  //avarage multiplicity
  fAvM+=fSelTracksIntFlow;
  //avarage of the generating function for integrated flow
  for(Int_t p=0;p<fgkPmax;p++){
   for(Int_t q=0;q<fgkQmax;q++){
    fAvG[p][q]+=1.*fG[p][q];
   } 
  }  
  //------------------------------------------------------------------------------------
  
  
  //STARTING WITH DIFFERENTIAL FLOW...  
  Int_t fBinEntries[fgknBins]={0};//stores the number of particles per bin for the current event
  Double_t fNumDRe[fgknBins][fgkPmax][fgkQmax]={0.};//real part of the numerator of D (see relation (10) in PG) 
  Double_t fNumDIm[fgknBins][fgkPmax][fgkQmax]={0.};//imaginary part of the numerator D   
  
  //------------------------------------------------------------------------------------------------
  //STARTING THE SECOND LOOP OVER TRACKS (CALCULATING THE GENERATING FUNCTION FOR DIFFERENTIAL FLOW)
  for(Int_t i=0;i<nPrim;i++){
    fTrack=anEvent->GetTrack(i);
    if (fTrack && fTrack->UseForDifferentialFlow()) {
      Int_t fBin=TMath::Nint(floor(fTrack->Pt()/fBinWidth));
      if(fBin>=fgknBins)continue;//ignoring the particles with pt>ptMax
      fBinNoOfParticles[fBin]++;
      fBinEntries[fBin]++;
      fBinMeanPt[fBin]+=fTrack->Pt();
      for(Int_t p=0;p<fgkPmax;p++){
	for(Int_t q=0;q<fgkQmax;q++){
	  fNumDRe[fBin][p][q]+=fG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/
	    (1.+(2.*fR0*sqrt(p+1.)/fSelTracksIntFlow) *
	     cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)); 
	  fNumDIm[fBin][p][q]+=fG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/
	    (1.+(2.*fR0*sqrt(p+1.)/fSelTracksIntFlow) * 
	     cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)); 
	}
      }
    }
  } 
  //ENDING THE SECOND LOOP OVER TRACKS 
  //----------------------------------------------------------------------------------------------- 
   
   
   
  //----------------------------------------------------------
  //AVARAGING OVER ALL pt BINS WITHIN ONE EVENT 
  for(Int_t b=0;b<fgknBins;b++){
    if(fBinEntries[b]==0)continue;
    fBinEventEntries[b]++;
    for(Int_t p=0;p<fgkPmax;p++){
      for(Int_t q=0;q<fgkQmax;q++){
	fBinEventDRe[b][p][q]+=fNumDRe[b][p][q]/fBinEntries[b];
	fBinEventDIm[b][p][q]+=fNumDIm[b][p][q]/fBinEntries[b];
      }
    }
  }
  //----------------------------------------------------------
}

//________________________________________________________________________
void AliFlowAnalysisWithCumulants::Terminate(Int_t nEvents){
  //final results
  cout<<""<<endl;
  cout<<"***************************************"<<endl;
  cout<<"**** results of cumulant analysis: ****"<<endl;
  cout<<"***************************************"<<endl;
  cout<<""<<endl;
  cout<<"number of events = "<<nEvents<<endl;
  
  //final avarage multiplicity
  fAvM/=nEvents;
  
  //final avarage of generating function for the integrated flow
  for(Int_t p=0;p<fgkPmax;p++){
   for(Int_t q=0;q<fgkQmax;q++){
    fAvG[p][q]/=nEvents;
   }
  }    
  
  //final avarage of the Q vector stuff
  fAvQx/=nEvents;
  fAvQy/=nEvents;
  fAvQ2x/=nEvents;
  fAvQ2y/=nEvents;
  
  /////////////////////////////////////////////////////////////////////////////      
  //////////////////gen. function for the cumulants////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  Double_t fC[fgkPmax][fgkQmax]={0.};
  for (Int_t p=0;p<fgkPmax;p++){
    for (Int_t q=0;q<fgkQmax;q++){
      fC[p][q]=1.*fAvM*(pow(fAvG[p][q],(1./fAvM))-1.);
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////
  ///////avaraging the gen. function for the cumulants over azimuth////////////
  /////////////////////////////////////////////////////////////////////////////
  
  Double_t fCAv[fgkPmax]={0.};
  for (Int_t p=0;p<fgkPmax;p++){ 
    Double_t fTempHere=0.; 
    for (Int_t q=0;q<fgkQmax;q++){
      fTempHere+=1.*fC[p][q];
    } 
    fCAv[p]=1.*fTempHere/fgkQmax;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////final results//////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  Double_t fCumulant[fgkPmax];//c_{iFlow}\{2(p+1)\}
  
  fCumulant[0]=(1./(fR0*fR0)) * (8.*fCAv[0] - 14.*fCAv[1] + (56./3.)*fCAv[2] - (35./2.)*fCAv[3] + 
			      (56./5.)*fCAv[4] - (14./3.)*fCAv[5] + (8./7.)*fCAv[6] - (1./8.)*fCAv[7]);

  fCumulant[1]=(1./pow(fR0,4.)) * ((-1924./35.)*fCAv[0] + (621./5.)*fCAv[1] - (8012./45.)*fCAv[2] + 
				 (691./4.)*fCAv[3] - (564./5.)*fCAv[4] + (2143./45.)*fCAv[5] - 
				 (412./35.)*fCAv[6] + (363./280.)*fCAv[7]);

  fCumulant[2]=(1./pow(fR0,6.)) * (349.*fCAv[0] - (18353./20.)*fCAv[1] + (7173./5.)*fCAv[2] - 
				 1457.*fCAv[3] + (4891./5.)*fCAv[4] - (1683./4.)*fCAv[5] + 
				 (527./5.)*fCAv[6] - (469./40.)*fCAv[7]);

  fCumulant[3]=(1./pow(fR0,8.)) * ((-10528./5.)*fCAv[0] + (30578./5.)*fCAv[1] - (51456./5.)*fCAv[2] + 
				 10993.*fCAv[3] - (38176./5.)*fCAv[4] + (16818./5.)*fCAv[5] - 
				 (4288./5.)*fCAv[6] + (967./10.)*fCAv[7]);

  fCumulant[4]=(1./pow(fR0,10.)) * (11500.*fCAv[0] - 35800.*fCAv[1] + 63900.*fCAv[2] - 71600.*fCAv[3] + 
				  51620.*fCAv[4] - 23400.*fCAv[5] + 6100.*fCAv[6] - 700.*fCAv[7]);

  fCumulant[5]=(1./pow(fR0,12.)) * (-52560.*fCAv[0] + 172080.*fCAv[1] - 321840.*fCAv[2] + 376200.*fCAv[3] - 
				  281520.*fCAv[4] + 131760.*fCAv[5] - 35280.*fCAv[6] + 4140.*fCAv[7]);

  fCumulant[6]=(1./pow(fR0,14.)) * (176400.*fCAv[0] - 599760.*fCAv[1] + 1164240.*fCAv[2] - 1411200.*fCAv[3] + 
				  1093680.*fCAv[4] - 529200.*fCAv[5] + 146160.*fCAv[6] - 17640.*fCAv[7]);

  fCumulant[7]=(1./pow(fR0,16.)) * (-322560*fCAv[0] + 1128960.*fCAv[1] - 2257920.*fCAv[2] + 2822400.*fCAv[3] - 
				  2257920.*fCAv[4] + 1128960.*fCAv[5] - 322560.*fCAv[6] + 40320.*fCAv[7]);
  
  cout<<""<<endl;
  cout<<"***************************"<<endl;
  cout<<"cumulants:"<<endl;
  
  cout<<" c_"<<fgkFlow<<"{2} = "<<fCumulant[0]<<endl; 
  cout<<" c_"<<fgkFlow<<"{4} = "<<fCumulant[1]<<endl;
  cout<<" c_"<<fgkFlow<<"{6} = "<<fCumulant[2]<<endl;
  cout<<" c_"<<fgkFlow<<"{8} = "<<fCumulant[3]<<endl; 
  cout<<"c_"<<fgkFlow<<"{10} = "<<fCumulant[4]<<endl; 
  cout<<"c_"<<fgkFlow<<"{12} = "<<fCumulant[5]<<endl;
  cout<<"c_"<<fgkFlow<<"{14} = "<<fCumulant[6]<<endl; 
  cout<<"c_"<<fgkFlow<<"{16} = "<<fCumulant[7]<<endl;  
  cout<<""<<endl;
  cout<<"integrated flow: "<<endl;
  
  Double_t fV2=0.,fV4=0.,fV6=0.,fV8=0.;
  Double_t fSdQ[4]={0.};
  Double_t fChiQ[4]={0.};
   if (fCumulant[0]>=0.){ 
    fV2=sqrt(fCumulant[0]);    
    fChiQ[0]=fAvM*fV2/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV2*fAvM,2.),0.5);
    fSdQ[0]=pow(((1./(2.*fAvM*nEvents))*((1.+1.*pow(fChiQ[0],2))/(1.*pow(fChiQ[0],2)))),0.5);
    cout<<" v_"<<fgkFlow<<"{2} = "<<100.*fV2<<"%, chi{2} = "<<fChiQ[0]<<", sd{2} = "<<100.*fSdQ[0]<<"%"<<endl;
    fCommonHistsRes2->FillIntegratedFlow(100.*fV2,100.*fSdQ[0]);
    fCommonHistsRes2->FillChi(fChiQ[0]);
   } else {
    cout<<" v_"<<fgkFlow<<"{2} = Im"<<endl;  
  }
  if (fCumulant[1]<=0.){
    fV4=pow(-fCumulant[1],(1./4.));
    fChiQ[1]=fAvM*fV4/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV4*fAvM,2.),0.5);
    fSdQ[1]=(1./(pow(2.*fAvM*nEvents,.5)))*pow((1.+2.*pow(fChiQ[1],2)+(1./4.)*pow(fChiQ[1],4.)+(1./4.)*pow(fChiQ[1],6.))/((1./4.)*pow(fChiQ[1],6.)),.5);
    cout<<" v_"<<fgkFlow<<"{4} = "<<100.*fV4<<"%, chi{4} = "<<fChiQ[1]<<", sd{4} = "<<100.*fSdQ[1]<<"%"<<endl;
    fCommonHistsRes4->FillIntegratedFlow(100.*fV4,100.*fSdQ[1]);
    fCommonHistsRes4->FillChi(fChiQ[1]);
  } else {
    cout<<" v_"<<fgkFlow<<"{4} = Im"<<endl;  
  } 
  if (fCumulant[2]>=0.){
    fV6=pow((1./4.)*fCumulant[2],(1./6.));
    fChiQ[2]=fAvM*fV6/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV6*fAvM,2.),0.5);
    fSdQ[2]=(1./(pow(2.*fAvM*nEvents,.5)))*pow((3.+18.*pow(fChiQ[2],2)+9.*pow(fChiQ[2],4.)+28.*pow(fChiQ[2],6.)+12.*pow(fChiQ[2],8.)+24.*pow(fChiQ[2],10.))/(24.*pow(fChiQ[2],10.)),.5);
    cout<<" v_"<<fgkFlow<<"{6} = "<<100.*fV6<<"%, chi{6} = "<<fChiQ[2]<<", sd{6} = "<<100.*fSdQ[2]<<"%"<<endl;
    fCommonHistsRes6->FillIntegratedFlow(100.*fV6,100.*fSdQ[2]);
    fCommonHistsRes6->FillChi(fChiQ[2]);
  } else {
    cout<<" v_"<<fgkFlow<<"{6} = Im"<<endl;  
  }
  if (fCumulant[3]<=0.){
    fV8=pow(-(1./33.)*fCumulant[3],(1./8.));
    fChiQ[3]=fAvM*fV8/pow(fAvQ2x+fAvQ2y-pow(fAvQx,2.)-pow(fAvQy,2.)-pow(fV8*fAvM,2.),0.5);
    fSdQ[3]=(1./(pow(2.*fAvM*nEvents,.5)))*pow((12.+96.*pow(fChiQ[3],2)+72.*pow(fChiQ[3],4.)+304.*pow(fChiQ[3],6.)+257.*pow(fChiQ[3],8.)+804.*pow(fChiQ[3],10.)+363.*pow(fChiQ[3],12.)+726.*pow(fChiQ[3],14.))/(726.*pow(fChiQ[3],14.)),.5);
    cout<<" v_"<<fgkFlow<<"{8} = "<<100.*fV8<<"%, chi{8} = "<<fChiQ[3]<<", sd{8} = "<<100.*fSdQ[3]<<"%"<<endl;
     fCommonHistsRes8->FillIntegratedFlow(100.*fV8,100.*fSdQ[3]);
     fCommonHistsRes8->FillChi(fChiQ[3]);
  } else {
    cout<<" v_"<<fgkFlow<<"{8} = Im"<<endl;
     
  }
  if (fCumulant[4]>=0.){
    cout<<"v_"<<fgkFlow<<"{10} = "<<100.*pow((1./456.)*fCumulant[4],(1./10.))<<"%"<<endl;
  } else {
    cout<<"v_"<<fgkFlow<<"{10} = Im"<<endl;  
  }
  if (fCumulant[5]<=0.){
    cout<<"v_"<<fgkFlow<<"{12} = "<<100.*pow(-(1./9460.)*fCumulant[5],(1./12.))<<"%"<<endl;
  } else {
    cout<<"v_"<<fgkFlow<<"{12} = Im"<<endl;  
  }
  if (fCumulant[6]>=0.){
    cout<<"v_"<<fgkFlow<<"{14} = "<<100.*pow((1./274800.)*fCumulant[6],(1./14.))<<"%"<<endl;
  } else {
    cout<<"v_"<<fgkFlow<<"{14} = Im"<<endl;  
  }
  if (fCumulant[7]<=0.){
    cout<<"v_"<<fgkFlow<<"{16} = "<<100.*pow(-(1./10643745.)*fCumulant[7],(1./16.))<<"%"<<endl;
  } else {
    cout<<"v_"<<fgkFlow<<"{16} = Im"<<endl;  
  }
  cout<<"***************************"<<endl;
  
  cout<<""<<endl;
  cout<<"continuing with calculations for differential flow..."<<endl;
  cout<<""<<endl;
 
  Double_t fBinEventDReAv[fgknBins][fgkPmax][fgkQmax]={0.};
  Double_t fBinEventDImAv[fgknBins][fgkPmax][fgkQmax]={0.};
  
  for(Int_t b=0;b<fgknBins;b++){
    if(fBinEventEntries[b]==0) continue;
    for(Int_t p=0;p<fgkPmax;p++){
      for(Int_t q=0;q<fgkQmax;q++){
	fBinEventDReAv[b][p][q]=fBinEventDRe[b][p][q]/fBinEventEntries[b];//avarage of the real part of numerator in relation (10) in PG 
	fBinEventDImAv[b][p][q]=fBinEventDIm[b][p][q]/fBinEventEntries[b];//avarage of the imaginary part of numerator in relation (10) in PG
      }
    }                                                                                                                                                                     
  }

  Double_t fX[fgknBins][fgkPmax][fgkQmax]={0.};//see the text bellow relation (11) in PG
  Double_t fY[fgknBins][fgkPmax][fgkQmax]={0.};
  
  for(Int_t b=0;b<fgknBins;b++){
    for(Int_t p=0;p<fgkPmax;p++){
      for(Int_t q=0;q<fgkQmax;q++){
	fX[b][p][q]=fBinEventDReAv[b][p][q]/fAvG[p][q];
	fY[b][p][q]=fBinEventDImAv[b][p][q]/fAvG[p][q];
      }
    }   
  } 
  cout<<""<<endl;
  cout<<"I have calculated X and Y."<<endl;
  cout<<""<<endl;
  
  Double_t fD[fgknBins][fgkPmax]={0.};//implementing relation (11) from PG
  
  for (Int_t b=0;b<fgknBins;b++){
    for (Int_t p=0;p<fgkPmax;p++){ 
      Double_t fTempHere3=0.; 
      for (Int_t q=0;q<fgkQmax;q++){
	fTempHere3+=cos(fgkMltpl*2.*q*TMath::Pi()/fgkQmax)*fX[b][p][q] + sin(fgkMltpl*2.*q*TMath::Pi()/fgkQmax)*fY[b][p][q];
      } 
      fD[b][p]=1.*(pow(fR0*pow(p+1.,.5),fgkMltpl)/fgkQmax)*fTempHere3;
    }
  } 
  
  cout<<""<<endl;
  cout<<"calculating differential cumulants now..."<<endl;
  cout<<""<<endl;
  
  Double_t fDiffCumulant2[fgknBins]={0.};//implementing relation (12) from PG
  Double_t fDiffCumulant4[fgknBins]={0.};
  Double_t fDiffCumulant6[fgknBins]={0.};
  Double_t fDiffCumulant8[fgknBins]={0.};
  
  for (Int_t b=0;b<fgknBins;b++){
    fDiffCumulant2[b]=(1./(fR0*fR0))*(4.*fD[b][0]-3.*fD[b][1]+(4./3.)*fD[b][2]-(1./4.)*fD[b][3]);
    fDiffCumulant4[b]=(1./pow(fR0,4.))*((-26./3.)*fD[b][0]+(19./2.)*fD[b][1]-(14./3.)*fD[b][2]+(11./12.)*fD[b][3]);
    fDiffCumulant6[b]=(1./pow(fR0,6.))*(18.*fD[b][0]-24.*fD[b][1]+14.*fD[b][2]-3.*fD[b][3]);
    fDiffCumulant8[b]=(1./pow(fR0,8.))*(-24.*fD[b][0]+36.*fD[b][1]-24.*fD[b][2]+6.*fD[b][3]);
  }
  
  Double_t fv2[fgknBins],fv4[fgknBins],fv6[fgknBins],fv8[fgknBins];
  Double_t fAvPt[fgknBins];
  Double_t fSddiff2[fgknBins],fSddiff4[fgknBins];

  cout<<"number of pt bins: "<<fgknBins<<endl;
  cout<<"****************************************"<<endl;
  for (Int_t b=0;b<fgknBins;b++){ 
    if(fBinNoOfParticles[b]==0)continue;
    fAvPt[b]=fBinMeanPt[b]/fBinNoOfParticles[b];
    cout<<"pt bin: "<<b*fBinWidth<<"-"<<(b+1)*fBinWidth<<" GeV"<<endl;
    cout<<"number of particles in this pt bin: "<<fBinNoOfParticles[b]<<endl;
    cout<<"mean pt in this bin: "<<fAvPt[b]<<" GeV"<<endl;
    if(fCumulant[0]>=0){
      fv2[b]=100.*fDiffCumulant2[b]/pow(fCumulant[0],.5);
      fSddiff2[b]=pow((1./(2.*fBinNoOfParticles[b]))*((1.+pow(fChiQ[0],2.))/pow(fChiQ[0],2.)),0.5);
      cout<<"v'_2/2{2} = "<<fv2[b]<<"%, "<<" "<<"sd{2} = "<<100.*fSddiff2[b]<<"%"<<endl;
      fCommonHistsRes2->FillDifferentialFlow(b+1,fv2[b],100.*fSddiff2[b]);
    }else{
      cout<<"v'_2/2{2} = Im"<<endl;
    } 
    if(fCumulant[1]<=0){
      fv4[b]=-100.*fDiffCumulant4[b]/pow(-fCumulant[1],.75);
      fSddiff4[b]=pow((1./(2.*fBinNoOfParticles[b]))*((2.+6.*pow(fChiQ[1],2.)+pow(fChiQ[1],4.)+pow(fChiQ[1],6.))/pow(fChiQ[1],6.)),0.5);
      cout<<"v'_2/2{4} = "<<fv4[b]<<"%, "<<" "<<"sd{4} = "<<100.*fSddiff4[b]<<"%"<<endl;
      fCommonHistsRes4->FillDifferentialFlow(b+1,fv4[b],100.*fSddiff4[b]);
    }else{
      cout<<"v'_2/2{4} = Im"<<endl;
    }  
    if(fCumulant[2]>=0){
      cout<<"v'_2/2{6} = "<<100.*fDiffCumulant6[b]/(4.*pow((1./4.)*fCumulant[2],(5./6.)))<<"%"<<endl;
      fv6[b]=100.*fDiffCumulant6[b]/(4.*pow((1./4.)*fCumulant[2],(5./6.)));
      fCommonHistsRes6->FillDifferentialFlow(b+1,fv6[b],0.);
    }else{
      cout<<"v'_2/2{6} = Im"<<endl;
    }     
    if(fCumulant[3]<=0){
      cout<<"v'_2/2{8} = "<<-100.*fDiffCumulant8[b]/(33.*pow(-(1./33.)*fCumulant[3],(7./8.)))<<"%"<<endl;
      fv8[b]=-100.*fDiffCumulant8[b]/(33.*pow(-(1./33.)*fCumulant[3],(7./8.))); 
      fCommonHistsRes8->FillDifferentialFlow(b+1,fv8[b],0.);
    }else{
      cout<<"v'_2/2{8} = Im"<<endl;
    }       
    cout<<"****************************************"<<endl;
  }  
}



