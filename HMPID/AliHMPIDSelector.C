#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>  //Terminate()
#include <TChain.h>
#include <TBenchmark.h>
#include <fstream>    //caf()      
#include <TProof.h>   //caf()
#include <AliSelector.h>      //base class
#include <AliESD.h>           


class AliHMPIDSelector : public AliSelector {
 public :
           AliHMPIDSelector():AliSelector(),fChain(0),fEsd(0),fCkovP(0),fMipXY(0),fDifXY(0),fSigP(0) {for(Int_t i=0;i<5;i++) fProb[i]=0;}
  virtual ~AliHMPIDSelector()                                                                      {delete fEsd;}


  virtual Int_t   Version        () const {return 2;}
  virtual void    Begin          (TTree *) {}
  virtual void    SlaveBegin     (TTree *tree);
  virtual void    Init           (TTree *tree);
  virtual Bool_t  Notify         ()    {return kTRUE;}
  virtual Bool_t  Process        (Long64_t entry);
  virtual void    SetOption      (const char *option) { fOption = option; }
  virtual void    SetObject      (TObject *obj) { fObject = obj; }
  virtual void    SetInputList   (TList *input) {fInput = input;}
  virtual TList  *GetOutputList  () const { return fOutput; }
  virtual void    SlaveTerminate ();
  virtual void    Terminate      ();

 private: 
  TTree          *fChain ;   //!pointer to the analyzed TTree or TChain
  AliESD         *fEsd ;     //!

  TH2F           *fCkovP,*fMipXY,*fDifXY,*fSigP; //!
  TH1F           *fProb[5];                      //!

  ClassDef(AliHMPIDSelector,0);  
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDSelector::SlaveBegin(TTree *tree)
{
// The SlaveBegin() function is called after the Begin() function.  When running with PROOF SlaveBegin() is called on each slave server.
// The tree argument is deprecated (on PROOF 0 is passed).

   Init(tree);

   TString option = GetOption();

   // create histograms on each slave server
   fCkovP    = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]", 150,   0,  7  ,100, -3, 1); 
   fSigP     = new TH2F("SigP"  ,"#sigma_{#theta_c}"          , 150,   0,  7  ,100, 0, 1e20);
   fMipXY    = new TH2F("MipXY" ,"mip position"               , 260,   0,130  ,252,0,126); 
   fDifXY    = new TH2F("DifXY" ,"diff"                       , 260, -10, 10  ,252,-10,10); 

   fProb[0] = new TH1F("PidE" ,"PID: e yellow #mu magenta"  ,100,0,1); fProb[0]->SetLineColor(kYellow);
   fProb[1] = new TH1F("PidMu","pid of #mu"                 ,100,0,1); fProb[1]->SetLineColor(kMagenta);
   fProb[2] = new TH1F("PidPi","PID: #pi red K green p blue",100,0,1); fProb[2]->SetLineColor(kRed);
   fProb[3] = new TH1F("PidK" ,"pid of K"                   ,100,0,1); fProb[3]->SetLineColor(kGreen);
   fProb[4] = new TH1F("PidP" ,"pid of p"                   ,100,0,1); fProb[4]->SetLineColor(kBlue);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDSelector::Init(TTree *pTr)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if ( !pTr )  return ;
   fChain = pTr ;
   fChain->SetBranchAddress("ESD", &fEsd) ;
   fChain->SetBranchStatus("*", 0);
   fChain->SetBranchStatus("fTracks.*", 1);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDSelector::Process(Long64_t entry)
{

  fChain->GetTree()->GetEntry(entry);

  for(Int_t iTrk=0;iTrk<fEsd->GetNumberOfTracks();iTrk++){
     AliESDtrack *pTrk=fEsd->GetTrack(iTrk);
     
//     if(pTrk->GetHMPIDsignal()<0) continue;
     
     fCkovP->Fill(pTrk->GetP(),pTrk->GetHMPIDsignal()) ; 
     fSigP ->Fill(pTrk->GetP(),TMath::Sqrt(pTrk->GetHMPIDchi2()));
     
//     Float_t xm,ym; Int_t q,np;  pTrk->GetHMPIDmip(xm,ym,q,np);  fMipXY->Fill(xm,ym); //mip info
//     Float_t xd,yd,th,ph;        pTrk->GetHMPIDtrk(xd,yd,th,ph); fDifXY->Fill(xd,yd); //track info 
     
     Double_t pid[5];  pTrk->GetHMPIDpid(pid); for(Int_t i =0;i<5;i++) fProb[i]->Fill(pid[i]);
  }//tracks loop 
     
  return kTRUE;
}//Process()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  // Add the histograms to the output on each slave server
  
  fOutput->Add(fCkovP);
  fOutput->Add(fSigP); 
  fOutput->Add(fMipXY);
  fOutput->Add(fDifXY);
  
  for(Int_t i=0;i<5;i++) fOutput->Add(fProb[i]);
}//SlaveTerminate()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  fCkovP   = dynamic_cast<TH2F*>(fOutput->FindObject("CkovP")) ;
  fSigP    = dynamic_cast<TH2F*>(fOutput->FindObject("SigP")) ; 
  fMipXY   = dynamic_cast<TH2F*>(fOutput->FindObject("MipXY")) ;
  fDifXY   = dynamic_cast<TH2F*>(fOutput->FindObject("DifXY")) ;
  
  fProb[0] = dynamic_cast<TH1F*>(fOutput->FindObject("PidE")) ;
  fProb[1] = dynamic_cast<TH1F*>(fOutput->FindObject("PidMu")) ;
  fProb[2] = dynamic_cast<TH1F*>(fOutput->FindObject("PidPi")) ;
  fProb[3] = dynamic_cast<TH1F*>(fOutput->FindObject("PidK")) ;
  fProb[4] = dynamic_cast<TH1F*>(fOutput->FindObject("PidP")) ;

  Float_t n=1.292; //mean freon ref idx 
  TF1 *pPi=new TF1("RiPiTheo","acos(sqrt(x*x+[0]*[0])/(x*[1]))",1.2,7); pPi->SetLineWidth(1); pPi->SetParameter(1,n); 
  AliPID ppp;                 pPi->SetLineColor(kRed);   pPi->SetParameter(0,AliPID::ParticleMass(AliPID::kPion));    //mass
  TF1 *pK=(TF1*)pPi->Clone(); pK ->SetLineColor(kGreen); pK ->SetParameter(0,AliPID::ParticleMass(AliPID::kKaon)); 
  TF1 *pP=(TF1*)pPi->Clone(); pP ->SetLineColor(kBlue);  pP ->SetParameter(0,AliPID::ParticleMass(AliPID::kProton)); 

  TCanvas *pC=new TCanvas("c1","ESD QA");pC->SetFillColor(10); pC->SetHighLightColor(10); pC->Divide(3,2);
  pC->cd(1); fCkovP->Draw(); pPi->Draw("same"); pK->Draw("same"); pP->Draw("same");   pC->cd(2); fMipXY->Draw();   pC->cd(3); fProb[0]->Draw(); fProb[1]->Draw("same"); 
  pC->cd(4); fSigP ->Draw();                                                          pC->cd(5); fDifXY->Draw();   pC->cd(6); fProb[2]->Draw(); fProb[3]->Draw("same"); fProb[4]->Draw("same"); 
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void loc()
{
  TChain* pChain =new TChain("esdTree");
  pChain->Add("AliESDs.root");

  pChain->Process("AliHMPIDSelector.C+");	
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void caf()
{
  gBenchmark->Start("PRooF exec");
  TChain* pChain =new TChain("esdTree");
  
  ifstream list; list.open("list.txt");

  TString file;
  while(list.good()) {
    list>>file;
    if (!file.Contains("root")) continue; //it's wrong file name
    pChain->Add(file.Data());
  }
  list.close();
  
  pChain->GetListOfFiles()->Print();
  
  TVirtualProof *pProof=TProof::Open("kir@lxb6046.cern.ch");	
  pProof->UploadPackage("ESD.par");
  pProof->EnablePackage("ESD");
  
  pChain->SetProof(pProof);
  pChain->Process("AliHMPIDSelector.C+");
  
  gBenchmark->Show("PRooF exec");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

