#define AliFlowAnalysisWithCumulants_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h" //NEW
#include "TParticle.h"







#include "TProfile.h"
#include "TProfile2D.h" 
#include "TProfile3D.h"






#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"

class TH1;
class TH3;
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
  //fHistFileName("cumulants.root"),//NEW: In the old version theis was not commented out
  //fHistFilefHistFile(NULL),//NEW: In the old version theis was not commented out
  fHistList(NULL), //NEW
  fAvM(0),
  fR0(0),
  fPtMax(0),
  fPtMin(0),
  fBinWidth(0),
  fAvQx(0),
  fAvQy(0),
  fAvQ2x(0),
  fAvQ2y(0),
  fNumberOfEvents(0),
  
  
  
  fHistProAvM(NULL),
  fIntFlowResults(NULL),
  fDiffFlowResults2(NULL),//to be improved
  fDiffFlowResults4(NULL),
  fDiffFlowResults6(NULL),
  fDiffFlowResults8(NULL),
  fIntFlowGenFun(NULL),
  fDiffFlowGenFunRe(NULL),
  fDiffFlowGenFunIm(NULL),
  
  
  
  fCommonHists(NULL)
  //fCommonHistsRes2(NULL),//to be improved
  //fCommonHistsRes4(NULL),//to be improved
  //fCommonHistsRes6(NULL),//to be improved
  //fCommonHistsRes8(NULL)//to be improved
 {
   //constructor 
   
   fHistList = new TList(); //NEW
   
   fR0=AliFlowCumuConstants::fgR0;
   fPtMax=AliFlowCommonConstants::GetPtMax(); 
   fPtMin=AliFlowCommonConstants::GetPtMin();
   fBinWidth=(fPtMax-fPtMin)/fgknBins;
   
 }

AliFlowAnalysisWithCumulants::~AliFlowAnalysisWithCumulants(){
  //desctructor
   delete fHistList; //NEW
}

  
//___________________________________________________________________________
void AliFlowAnalysisWithCumulants::CreateOutputObjects(){
 //output histograms
 
  fHistProAvM = new TProfile("fHistProAvM","Avarage Multiplicity",1,0,1,0,100000);
  fHistProAvM->SetXTitle("");
  fHistProAvM->SetYTitle("Avarage Multiplicity");
  fHistList->Add(fHistProAvM);
  
  fIntFlowResults = new TH1D("fIntFlowResults","Integrated Flow From Cumulants",8,0,8);
  fIntFlowResults->SetXTitle("");
  fIntFlowResults->SetYTitle("Integrated Flow [%]");
  fHistList->Add(fIntFlowResults);
 
  //to be improved, I should store the results as CommonHistResults
  fDiffFlowResults2 = new TH1D("fDiffFlowResults2","v'_2/2{2}",fgknBins,fPtMin,fPtMax);
  fDiffFlowResults2->SetXTitle("pt [GeV]");
  fDiffFlowResults2->SetYTitle("Differential Flow [%]");
  fHistList->Add(fDiffFlowResults2);
 
  //to be improved, I should store the results as CommonHistResults
  fDiffFlowResults4 = new TH1D("fDiffFlowResults4","v'_2/2{4}",fgknBins,fPtMin,fPtMax);
  fDiffFlowResults4->SetXTitle("pt [GeV]");
  fDiffFlowResults4->SetYTitle("Differential Flow [%]");
  fHistList->Add(fDiffFlowResults4);
  
  //to be improved, I should store the results as CommonHistResults
  fDiffFlowResults6 = new TH1D("fDiffFlowResults6","v'_2/2{6}",fgknBins,fPtMin,fPtMax);
  fDiffFlowResults6->SetXTitle("pt [GeV]");
  fDiffFlowResults6->SetYTitle("Differential Flow [%]");
  fHistList->Add(fDiffFlowResults6);
  
  //to be improved, I should store the results as CommonHistResults
  fDiffFlowResults8 = new TH1D("fDiffFlowResults8","v'_2/2{8}",fgknBins,fPtMin,fPtMax);
  fDiffFlowResults8->SetXTitle("pt [GeV]");
  fDiffFlowResults8->SetYTitle("Differential Flow [%]");
  fHistList->Add(fDiffFlowResults8);
 
  fIntFlowGenFun = new TProfile2D("fIntFlowGenFun","G[p][q]",8,0.,8.,17,0.,17.);//to be improved (z_down and z_up)
  fIntFlowGenFun->SetXTitle("p");
  fIntFlowGenFun->SetYTitle("q");
  fHistList->Add(fIntFlowGenFun);
     
  fDiffFlowGenFunRe = new TProfile3D("fDiffFlowGenFunRe","Re(D[p][q])",fgknBins,fPtMin/fBinWidth,fPtMax/fBinWidth,8,0.,8.,17,0.,17.);
  fHistList->Add(fDiffFlowGenFunRe);
 
  fDiffFlowGenFunIm = new TProfile3D("fDiffFlowGenFunIm","Im(D[p][q])",fgknBins,fPtMin/fBinWidth,fPtMax/fBinWidth,8,0.,8.,17,0.,17.);
  fHistList->Add(fDiffFlowGenFunIm);
  

 /* 
 fCommonHistsRes2 = new AliFlowCommonHistResults("Cumulants2");//to be improved
 fCommonHistsRes4 = new AliFlowCommonHistResults("Cumulants4");//to be improved
 fCommonHistsRes6 = new AliFlowCommonHistResults("Cumulants6");//to be improved
 fCommonHistsRes8 = new AliFlowCommonHistResults("Cumulants8");//to be improved
 
 
 fHistList->Add(fCommonHistsRes2->GetHistList()); //NEW //to be improved
 fHistList->Add(fCommonHistsRes4->GetHistList()); //NEW //to be improved
 fHistList->Add(fCommonHistsRes6->GetHistList()); //NEW //to be improved 
 fHistList->Add(fCommonHistsRes8->GetHistList()); //NEW //to be improved
 */

 //control histograms
 fCommonHists = new AliFlowCommonHist("Cumulants");
 fHistList->Add(fCommonHists->GetHistList()); 
 
}

//________________________________________________________________________
void AliFlowAnalysisWithCumulants::Make(AliFlowEventSimple* anEvent) {
  //running over data
  
  //---------------------------------------------------------  
  //fill the common control histograms
  fCommonHists->FillControlHistograms(anEvent);   
  //---------------------------------------------------------
  
  
  //---------------------------------------------------------
  //initializing the generating function for integrated flow
  Double_t fG[fgkPmax][fgkQmax];
  for(Int_t p=0;p<fgkPmax;p++){
   for(Int_t q=0;q<fgkQmax;q++){
    fG[p][q]=1.; 
   }
  }
  //---------------------------------------------------------
  
  
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
  }// ENDING THE FIRST LOOP OVER TRACKS
  //------------------------------------------------------------------------------------ 
  
  
  fHistProAvM->Fill(0.,fSelTracksIntFlow,1);
  fAvM=fHistProAvM->GetBinContent(1);

  
  //------------------------------------------------------------------------------------
  //STORING THE VALUE OF GENERATING FUNCTION FOR INTEGRATED FLOW INTO 2D PROFILE
  for(Int_t p=0;p<fgkPmax;p++){
   for(Int_t q=0;q<fgkQmax;q++){
    fIntFlowGenFun->Fill((Double_t)p,(Double_t)q,fG[p][q],1);
   } 
  }
  //------------------------------------------------------------------------------------
   
                
  //------------------------------------------------------------------------------------
  //STARTING THE SECOND LOOP (CALCULATING THE GENERATING FUNCTION FOR DIFFERENTIAL FLOW)
  // Remark 0: generating function for diff. flow is complex number, I need to calcuate separately real and imaginary part
  // Remark 1: note that I need here fG[p][q], the value of generating function for integrated flow for the CURRENT event
  // Remark 2: results are immediately stored in two 3D profiles, one for real and one for imaginary part
  for(Int_t i=0;i<nPrim;i++){
   fTrack=anEvent->GetTrack(i);
   if (fTrack && fTrack->UseForDifferentialFlow()){
    for(Int_t p=0;p<fgkPmax;p++){
     for(Int_t q=0;q<fgkQmax;q++){
      fDiffFlowGenFunRe->Fill(fTrack->Pt()/fBinWidth,(Double_t)p,(Double_t)q,fG[p][q]*cos(fgkMltpl*fgkFlow*fTrack->Phi())/
	    (1.+(2.*fR0*sqrt(p+1.)/fSelTracksIntFlow) *
	     cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);
      fDiffFlowGenFunIm->Fill(fTrack->Pt()/fBinWidth,(Double_t)p,(Double_t)q,fG[p][q]*sin(fgkMltpl*fgkFlow*fTrack->Phi())/
	    (1.+(2.*fR0*sqrt(p+1.)/fSelTracksIntFlow) * 
	     cos(fgkFlow*fTrack->Phi()-2.*q*TMath::Pi()/fgkQmax)),1.);        
     }
    }
   }  
  }// ENDING THE SECOND LOOP OVER TRACKS
  //------------------------------------------------------------------------------------
  
}//end of Make()

//________________________________________________________________________
void AliFlowAnalysisWithCumulants::Finish(){
  //final results
}


















