/*
#if !defined(__CINT__) || defined(__MAKECINT__)
*/
#include <Riostream.h>     

#include <TSystem.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1I.h>

#include "AliHMPIDDigit.h"
#include "AliHMPIDRawStream.h"
/*
#endif
*/


TH2F hgPedMapMean[14][6]; 
TH2F hgPedMapSigma[14][6];
TH1F hgPedMapMean1D[14][6]; 
TH1F hgPedMapSigma1D[14][6]; 
TH2F hgPedMapMeanSigma[14][6];

TH1I hgDdlErr[14];

Int_t fgRunNum;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void Convert(Int_t ddl,Int_t r,Int_t d,Int_t a,Int_t &ch, Int_t &pc, Int_t &px, Int_t &py)
{
 
  
  Int_t a2y[6]={3,2,4,1,5,0};//pady for a given address (for single DILOGIC chip)
  
  Int_t ch=ddl/2;
  Int_t tmp=(24-r)/8;              Int_t pc=(ddl%2)? 5-2*tmp:2*tmp; 
                                  Int_t px=(d-1)*8+a/6;
        tmp=(ddl%2)?r-1:(24-r);   Int_t py=6*(tmp%8)+a2y[a%6];
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ProcPed(Int_t nDDL)
{

    if(gSystem->IsFileInIncludePath(Form("./HmpidPedDdl%02i.txt",nDDL)))
          ifstream infile(Form("./HmpidPedDdl%02i.txt",nDDL));
    else return;
    
    for(Int_t ni=0;ni<6;ni++)
      {
        hgPedMapMean[nDDL][ni]  = new TH2F(Form("hPedMapMean_DDL%d_PC_%d" ,nDDL,ni),         Form("Pedestal: Mean, DDL(0-13): %d PC(0-5): %d;padx;pady" ,nDDL,ni),80,0,80,48,0,48);
        hgPedMapSigma[nDDL][ni] = new TH2F(Form("hPedMapSigma_DDL%d_PC_%d",nDDL,ni),         Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;padx;pady",nDDL,ni),80,0,80,48,0,48);
        hgPedMapMean1D[nDDL][ni]  = new TH1F(Form("hPedMapMean1D_DDL%d_PC_%d" ,nDDL,ni),     Form("Pedestal: Mean, DDL(0-13): %d PC(0-5): %d;pad;Mean" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapSigma1D[nDDL][ni]  = new TH1F(Form("hPedMapSigma1D_DDL%d_PC_%d" ,nDDL,ni),   Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;pad;Sigma" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapMeanSigma[nDDL][ni] = new TH2F(Form("hPedMapMeanSigma_DDL%d_PC_%d" ,nDDL,ni),Form("Pedestal, DDL(0-13): %d PC(0-5): %d;Mean;Sigma"),300,0,300,50,0,5);
       }
      
    Int_t nSigCut,r,d,a,hard;  Float_t mean,sigma;
    Int_t ch=0,pc=0,px=0,py=0;
    Int_t cnt=0;
    Int_t runNumber=-99999;
    Char_t tName[10];
    Int_t  ldcId;
    Int_t  timeStamp;
    Int_t  nEv,nDdlEv;
    Int_t nBadEv;Float_t nBadEvPer;
    Printf("Start reading DDL: %d ...",nDDL);  
    infile>>tName>>runNumber;
    infile>>tName>>ldcId;
    infile>>tName>>timeStamp;
    infile>>tName>>nEv;  
    infile>>tName>>nDdlEv;  
    infile>>tName>>nBadEv;
    infile>>tName>>nBadEvPer;
    infile>>tName>>nSigCut;
    
    Printf("RunNumber: %d, LdcId: %d TimeStamp: %d nEv: %d",runNumber,ldcId,timeStamp,nEv);
    
    while(!infile.eof()){
      infile>>dec>>r>>d>>a>>mean>>sigma>>hex>>hard;
      Convert(nDDL,r,d,a,ch,pc,px,py); 
      hgPedMapMean[nDDL][pc].Fill(px,py,mean);  
      hgPedMapSigma[nDDL][pc].Fill(px,py,sigma);
      hgPedMapMean1D[nDDL][pc].Fill(cnt%3840,mean);
      hgPedMapSigma1D[nDDL][pc].Fill(cnt%3840,sigma);
      hgPedMapMeanSigma[nDDL][pc].Fill(mean,sigma);
      cnt++;      
    }
  infile.close();
  Printf("Stop reading DDL: %d ...",nDDL);
  
  /* fill the overall histos */
  
  fgRunNum=runNumber;    
  gStyle->SetPalette(1);    
  
  TCanvas *cped = new TCanvas("cped","cped",800,800);     cped->Divide(2,2);
  TCanvas *cped2 = new TCanvas("cped2","cped2",800,800);  cped2->Divide(2,2);
  
  for(Int_t npc=0;npc<6;npc++) {
    if(hgPedMapMean[nDDL][npc].GetEntries()==0) continue;
    cped->cd(1); hgPedMapMean[nDDL][npc].Draw("surf1");
    cped->cd(2); hgPedMapSigma[nDDL][npc].Draw("surf1");
    
    cped->cd(3); hgPedMapMean[nDDL][npc].SetStats(0);hgPedMapMean[nDDL][npc].Draw("colz");
    cped->cd(4); hgPedMapSigma[nDDL][npc].SetStats(0);hgPedMapSigma[nDDL][npc].Draw("colz");
    cped->SaveAs(Form("PedMap1_DDL%d_PC_%d.eps",nDDL,npc));
    cped->SaveAs(Form("PedMap1_DDL%d_PC_%d.gif",nDDL,npc));
    //cped->SaveAs(Form("PedMap1_DDL%d_PC_%d.pdf",nDDL,npc));
        
    cped2->cd(1); hgPedMapMean1D[nDDL][npc].Draw();
    cped2->cd(2); hgPedMapSigma1D[nDDL][npc].Draw();
    cped2->cd(3); hgPedMapMeanSigma[nDDL][npc].Draw("colz");
    cped2->SaveAs(Form("PedMap2_DDL%d_PC_%d.eps",nDDL,npc));           
    cped2->SaveAs(Form("PedMap2_DDL%d_PC_%d.gif",nDDL,npc));           
    //cped2->SaveAs(Form("PedMap2_DDL%d_PC_%d.pdf",nDDL,npc));           
    
  }
    
       
}//ProcPed()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void ProcErr(Int_t nDDL)
{
  if(gSystem->IsFileInIncludePath(Form("./HmpidErrorsDdl%02i.txt",nDDL)))
          ifstream infile(Form("./HmpidErrorsDdl%02i.txt",nDDL));
    else return;
    
    hgDdlErr[nDDL] = new TH1I(Form("hPedErr_DDL%d",nDDL),Form("DDL(%d) Decoding Errors",nDDL),AliHMPIDRawStream::kSumErr+1,-0.5,AliHMPIDRawStream::kSumErr+0.5);
     for(Int_t ilabel=0; ilabel< AliHMPIDRawStream::kSumErr; ilabel++) {
      hgDdlErr[nDDL].SetStats(0);
      hgDdlErr[nDDL].GetXaxis()->CenterLabels(kTRUE);
      hgDdlErr[nDDL].GetXaxis()->SetBinLabel((ilabel+1),Form("%i  %s",ilabel+1,AliHMPIDRawStream::GetErrName(ilabel)));
      hgDdlErr[nDDL].SetYTitle("Error #");
     }
    hgDdlErr[nDDL].SetFillColor(5);

    
    Int_t runNumber=-99999;
    Char_t tName[10];
    Int_t  ldcId;
    Int_t  timeStamp;
    Int_t  nEv,nDdlEv;
    Int_t rerr;Int_t nBadEv;Float_t nBadEvPer;
    Printf("Start reading Error File DDL: %d ...",nDDL);  
    infile>>tName>>runNumber;
    infile>>tName>>ldcId;
    infile>>tName>>timeStamp;
    infile>>tName>>nEv;
    infile>>tName>>nDdlEv;
    infile>>tName>>nBadEv;
    infile>>tName>>nBadEvPer;
    Printf("RunNumber: %d, LdcId: %d TimeStamp: %d nEv: %d",runNumber,ldcId,timeStamp,nEv);
    for(Int_t ierr=0;ierr<AliHMPIDRawStream::kSumErr;ierr++)
    {
       infile>>rerr;hgDdlErr[nDDL].SetBinContent(ierr+1,rerr);
     }
    
    
    TCanvas *cerr = new TCanvas("cped","cped");
    cerr.cd();
    hgDdlErr[nDDL]->Draw();
    cerr->SaveAs(Form("PedError_DDL%d.eps",nDDL)); 
    cerr->SaveAs(Form("PedError_DDL%d.gif",nDDL)); 
      
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void PlotErrAllCh()
{
  
  TCanvas *chmpid1= new TCanvas("chmpid1","chmpid1",1024,768);
  TPaveLabel *pl = new TPaveLabel(1,16.3,24,17.5,Form("HMPID Pedestal values Run #: %d",fgRunNum),"br");
   pl->SetFillColor(18);
   pl->SetTextFont(32);
   pl->SetTextColor(49);
   pl->Draw();
  chmpid1->Divide(6,9);
  chmpid1->cd(1);
  hgPedMapMean[0][0].Draw("colz");
  
  chmpid1->SaveAs(Form("hmpid_pedestal_run_%d.eps",fgRunNum));
  
  
  TCanvas *chmpid3= new TCanvas("chmpid3","chmpid3",1280,960);
           chmpid3->Divide(6,3);

  /* chamber layout */        
           
          chmpid3->cd(1); hgDdlErr[12]->Draw(); chmpid3->cd(2); hgDdlErr[13]->Draw();  /* */ chmpid3->cd(3);  hgDdlErr[10]->Draw(); chmpid3->cd(4);  hgDdlErr[11]->Draw();  /* */ /* empty empty */
          chmpid3->cd(7); hgDdlErr[8]->Draw();  chmpid3->cd(8); hgDdlErr[9]->Draw();  /* */  chmpid3->cd(9);  hgDdlErr[6]->Draw();  chmpid3->cd(10); hgDdlErr[7]->Draw();  /* */ chmpid3->cd(11); hgDdlErr[4]->Draw();  chmpid3->cd(12); hgDdlErr[5]->Draw(); 
          /* empty empty */                                                           /* */  chmpid3->cd(15); hgDdlErr[2]->Draw();  chmpid3->cd(16); hgDdlErr[3]->Draw();  /* */ chmpid3->cd(17); hgDdlErr[0]->Draw();  chmpid3->cd(18); hgDdlErr[1]->Draw(); 
        
      chmpid3->SaveAs(Form("hmpid_pedestal_run_%d_errors.eps",fgRunNum));          
      chmpid3->SaveAs(Form("hmpid_pedestal_run_%d_errors.gif",fgRunNum));          
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void HplotDA()
{
  for(Int_t i=0;i<14;i++) 
  {
   ProcPed(i);  
   ProcErr(i);
  }
  
  PlotErrAllCh();
  
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
