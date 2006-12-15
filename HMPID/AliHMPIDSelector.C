#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>  //Terminate()
#include <TChain.h>
#include <TBenchmark.h>
#include <TFile.h>    //docosmic()    
#include <fstream>    //caf()      
#include <TProof.h>   //caf()
#include <AliSelector.h>      //base class
#include <AliESD.h>           
#include <AliBitPacking.h> //HmpidPayload()
#include "AliHMPIDDigit.h" 
#include "AliHMPIDCluster.h" 
#include "AliHMPIDReconstructor.h" //docosmic()

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
Int_t DateHeader(ifstream *pFile,Bool_t isPrint=0)
{
  Int_t iSize=-1;
  pFile->read((char*)&iSize,4);
  if(!isPrint)     
    pFile->seekg(16*4,ios::cur);
  else{
    Int_t w32=-1; 
                                Printf("");
                                Printf("Event size        %i bytes",iSize);                           //1  common DATE header 17 words
    pFile->read((char*)&w32,4); Printf("Event magic       0x%x"    ,w32);                             //2
    pFile->read((char*)&w32,4); Printf("Event head size   %i bytes",w32);                             //3  
    pFile->read((char*)&w32,4); Printf("Event version     0x%x"    ,w32);                             //4
    pFile->read((char*)&w32,4); Printf("Event type        %i (%s)" ,w32,(w32==7)? "physics":"SOR");   //5 
    pFile->read((char*)&w32,4); Printf("Run number        %i"      ,w32);                             //6
  
    pFile->read((char*)&w32,4); Printf("Event ID 1        %i"      ,w32);                             //7
    pFile->read((char*)&w32,4); Printf("Event ID 2        %i"      ,w32);                             //8
  
    pFile->read((char*)&w32,4); Printf("Trigger pattern 1 %i"      ,w32);                             //9
    pFile->read((char*)&w32,4); Printf("Trigger pattern 2 %i"      ,w32);                             //10
  
    pFile->read((char*)&w32,4); Printf("Detector pattern  %i"      ,w32);                             //11
  
    pFile->read((char*)&w32,4); Printf("Type attribute 1  %i"      ,w32);                             //12
    pFile->read((char*)&w32,4); Printf("Type attribute 2  %i"      ,w32);                             //13
    pFile->read((char*)&w32,4); Printf("Type attribute 3  %i"      ,w32);                             //14
  
    pFile->read((char*)&w32,4); Printf("LDC ID            %i"      ,w32);                             //15
    pFile->read((char*)&w32,4); Printf("GDC ID            %i"      ,w32);                             //16
    pFile->read((char*)&w32,4); TDatime time(w32); time.Print();                                      //17
  
    Printf("");
  }
  return iSize;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpidHeader(ifstream *pFile,Bool_t isPrint=kFALSE)
{
// Prints hmpid trailer and returns number of words for this trailer  
  if(!isPrint) {pFile->seekg(15*4,ios::cur);return;}
  Int_t w32=-1;
  Printf("\nHMPID Header:");//private HEADER is 15 words
  for(Int_t i=1;i<=11;i++) { pFile->read((char*)&w32,4); Printf("Word #%2i=%12i meaningless",i,w32);}
                             pFile->read((char*)&w32,4); Printf("Word #12=%12i event counter",w32);   
  for(Int_t i=13;i<=15;i++){ pFile->read((char*)&w32,4); Printf("Word #%2i=%12i empty",i,w32);}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HmpidPayload(ifstream *pFile,Int_t iDdl,TObjArray *pDigAll)
{
// payload is 8 sequences with structure: WC36A8 then WC number of w32
  
  TClonesArray *pDig1=(TClonesArray*)pDigAll->At(iDdl/2); //get list of digits for requested chamber
  
  UInt_t w32=0;
  Int_t iDigCnt=pDig1->GetEntriesFast();
  for(Int_t row=1;row<=8;row++){  
    pFile->read((char*)&w32,4);  Int_t wc=AliBitPacking::UnpackWord(w32,16,31); Int_t ma=AliBitPacking::UnpackWord(w32, 0,15);    
    if(ma!=0x36a8) Printf("ERROR ERROR ERROR WRONG Marker=0x%x ERROR ERROR ERROR",ma);    
    for(Int_t i=1;i<=wc;i++){//words loop
      pFile->read((char*)&w32,4);
      if(w32&0x08000000) continue; //it's DILOGIC CW
      AliHMPIDDigit *pDig=new AliHMPIDDigit;
      pDig->Raw(iDdl,w32);   
      new ((*pDig1)[iDigCnt++]) AliHMPIDDigit(*pDig);
    }//words loop 
  }//rows loop
}//HmpidPayload()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void docosmic(const char* name,Int_t iMaxEvt=9999999)
{
  gBenchmark->Start("Cosmic converter");
  ifstream in(name);

  
  Bool_t isPrint=kFALSE;
  
  TString rooName=name; rooName.Replace(1+rooName.First('.'),4,"root");
  Int_t ddl=0;    
  
  TFile *pOut=new TFile(rooName,"recreate");
  TTree *pTr=new TTree("cosmic","some for time being");
  
  
  TObjArray *pDigAll=new TObjArray(7); pDigAll->SetOwner(); for(Int_t i=0;i<7;i++) pDigAll->AddAt(new TClonesArray("AliHMPIDDigit")  ,i); pTr->Branch("Digs",&pDigAll,64000,0); 
  TObjArray *pCluAll=new TObjArray(7); pCluAll->SetOwner(); for(Int_t i=0;i<7;i++) pCluAll->AddAt(new TClonesArray("AliHMPIDCluster"),i); pTr->Branch("Clus",&pCluAll,64000,0);

  Int_t iEvt=0;
  while(1){      
    Int_t iSize=DateHeader(&in,      isPrint);  if(iSize==68) continue;  //DATE header 
    if(in.eof()) break;
    HmpidHeader (&in,      isPrint);    //HMPID header 
    HmpidPayload(&in,ddl+1,pDigAll);    //HMPID payload
    HmpidHeader (&in,      isPrint);    //HMPID header 
    HmpidPayload(&in,ddl  ,pDigAll);    //HMPID payload
    
    AliHMPIDReconstructor::Dig2Clu(pDigAll,pCluAll);
    pTr->Fill();
    for(Int_t i=0;i<7;i++){
      ((TClonesArray*)pDigAll->At(i))->Clear(); 
      ((TClonesArray*)pCluAll->At(i))->Clear(); 
    }
    
    if(!(iEvt%200)) Printf("Event %i processed",iEvt);
    iEvt++;
    if(iEvt==iMaxEvt) break;
  }//events loop
  
  pTr->Write();  pOut->Close();
  pDigAll->Delete(); pCluAll->Delete();
  in.close();
  
  Printf("Total %i events processed",iEvt);
  gBenchmark->Show("Cosmic converter");
}//docosmic()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cosmic()
{
  TChain *pCh=new TChain("cosmic");  pCh->Add("test.root");

  
  TH1F *pDigQ=new TH1F("digQ","Digits QDC",500,0,4100);
  TH1F *pDigO=new TH1F("digO","Digits # per event",500,0,8000);
  TH2F *pDigM=new TH2F("digM","Digits map",500,0,131,500,0,127);
    
  TH1F *pCluQ   =new TH1F("cluQ","Clusters QDC",500,0,12100);
  TH1F *pCluQmax=new TH1F("cluQmax","Clusters Max QDC",500,0,12100); pCluQmax->SetLineColor(kRed);
  TH1F *pCluO=new TH1F("cluO","Clusters # per event",500,0,5000);
  TH2F *pCluM=new TH2F("cluM","Clusters map",500,0,131,500,0,127);
  
  TObjArray *pDigAll=0; pCh->SetBranchAddress("Digs",&pDigAll);
  TObjArray *pCluAll=0; pCh->SetBranchAddress("Clus",&pCluAll);

    
  for(Int_t iEvt=0;iEvt<pCh->GetEntries();iEvt++){
    pCh->GetEntry(iEvt);
    
    TClonesArray *pDigCh=(TClonesArray*)pDigAll->At(0);
    TClonesArray *pCluCh=(TClonesArray*)pCluAll->At(0);
    pDigO->Fill(pDigCh->GetEntriesFast());
    pCluO->Fill(pCluCh->GetEntriesFast());
    
    for(Int_t iDig=0;iDig<pDigCh->GetEntriesFast();iDig++){//digits loop
      AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigCh->UncheckedAt(iDig);
      pDigQ->Fill(pDig->Q());
      pDigM->Fill(pDig->LorsX(),pDig->LorsY());
    }//digits loop
    Int_t qmax=0;    
    for(Int_t iClu=0;iClu<pCluCh->GetEntriesFast();iClu++){//clusters loop
      AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluCh->UncheckedAt(iClu);
      pCluQ->Fill(pClu->Q());
      if(pClu->Q()>qmax) qmax=pClu->Q();
      pCluM->Fill(pClu->X(),pClu->Y());
    }//digits loop
    pCluQmax->Fill(qmax);
  }//entries loop
  
  TCanvas *pC=new TCanvas("comic","cosmic"); pC->Divide(2,3);
  
  pC->cd(1); pDigM->Draw(); pC->cd(2); pCluM->Draw();
  pC->cd(3); gPad->SetLogy(); pDigQ->Draw(); pC->cd(4); gPad->SetLogy(); pCluQ->Draw(); pCluQmax->Draw("same");
  pC->cd(5); gPad->SetLogy(); pDigO->Draw(); pC->cd(6); gPad->SetLogy(); pCluO->Draw();
}//cosmic()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

