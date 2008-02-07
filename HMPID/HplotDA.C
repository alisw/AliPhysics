#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>     

#include <TSystem.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1I.h>

#include "AliHMPIDRawStream.h"
#endif



TH2F hgPedMapMean[6]; 
TH2F hgPedMapSigma[6];
TH1F hgPedMapMean1D[6]; 
TH1F hgPedMapSigma1D[6]; 
TH2F hgPedMapMeanSigma[6];
TH1I hgDdlErr;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Convert(Int_t ddl,Int_t r,Int_t d,Int_t a,Int_t &ch, Int_t &pc, Int_t &px, Int_t &py)
{
  
  Int_t a2y[6]={3,2,4,1,5,0};//pady for a given address (for single DILOGIC chip)
  
  Int_t ch=ddl/2;
  Int_t tmp=(r-1)/8;              Int_t pc=(ddl%2)? 5-2*tmp:2*tmp; 
                                  Int_t px=(d-1)*8+a/6;
        tmp=(ddl%2)?(24-r):r-1;   Int_t py=6*(tmp%8)+a2y[a%6];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void ProcPed(Int_t nDDL)
{

    if(gSystem->IsFileInIncludePath(Form("./HmpidPedDdl%02i.txt",nDDL)))
          ifstream infile(Form("./HmpidPedDdl%02i.txt",nDDL));
    else return;
    
    for(Int_t ni=0;ni<6;ni++)
      {
        hgPedMapMean[ni]  = new TH2F(Form("hPedMapMean_DDL%d_PC_%d" ,nDDL,ni),         Form("Pedestal: Mean, DDL(0-13): %d PC(0-5): %d;padx;pady" ,nDDL,ni),80,0,80,48,0,48);
        hgPedMapSigma[ni] = new TH2F(Form("hPedMapSigma_DDL%d_PC_%d",nDDL,ni),         Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;padx;pady",nDDL,ni),80,0,80,48,0,48);
        hgPedMapMean1D[ni]  = new TH1F(Form("hPedMapMean1D_DDL%d_PC_%d" ,nDDL,ni),     Form("Pedestal: Mean, DDL(0-13): %d PC(0-5): %d;pad;Mean" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapSigma1D[ni]  = new TH1F(Form("hPedMapSigma1D_DDL%d_PC_%d" ,nDDL,ni),   Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;pad;Sigma" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapMeanSigma[ni] = new TH2F(Form("hPedMapMeanSigma_DDL%d_PC_%d" ,nDDL,ni),Form("Pedestal, DDL(0-13): %d PC(0-5): %d;Mean;Sigma"),300,0,300,50,0,5);
       }
      
    Int_t nSigCut,r,d,a,hard;  Float_t mean,sigma;
    Int_t ch=0,pc=0,px=0,py=0;
    Int_t cnt=0;
    Int_t runNumber=-99999;
    Char_t tName[10];
    Int_t  ldcId;
    Int_t  timeStamp;
    Int_t  nEv;
    Printf("Start reading DDL: %d ...",nDDL);  
    infile>>tName>>runNumber;
    infile>>tName>>ldcId;
    infile>>tName>>timeStamp;
    infile>>tName>>nEv;  
    infile>>tName>>nSigCut;
    Printf("RunNumber: %d, LdcId: %d TimeStamp: %d nEv: %d",runNumber,ldcId,timeStamp,nEv);
    while(!infile.eof()){
      infile>>dec>>r>>d>>a>>mean>>sigma>>hex>>hard;
      Convert(nDDL,r,d,a,ch,pc,px,py); 
      hgPedMapMean[pc].Fill(px,py,mean);  
      hgPedMapSigma[pc].Fill(px,py,sigma);
      hgPedMapMean1D[pc].Fill(cnt%3840,mean);
      hgPedMapSigma1D[pc].Fill(cnt%3840,sigma);
      hgPedMapMeanSigma[pc].Fill(mean,sigma);
      cnt++;
      
    }
  infile.close();
  Printf("Stop reading DDL: %d ...",nDDL);
      
  gStyle->SetPalette(1);    
  
  TCanvas *cped = new TCanvas("cped","cped",800,800);     cped->Divide(2,2);
  TCanvas *cped2 = new TCanvas("cped2","cped2",800,800);  cped2->Divide(2,2);
  
  for(Int_t npc=0;npc<6;npc++) {
    if(hgPedMapMean[npc].GetEntries()==0) continue;
    cped->cd(1); hgPedMapMean[npc].Draw("surf1");
    cped->cd(2); hgPedMapSigma[npc].Draw("surf1");
    
    cped->cd(3); hgPedMapMean[npc].SetStats(0);hgPedMapMean[npc].Draw("colz");
    cped->cd(4); hgPedMapSigma[npc].SetStats(0);hgPedMapSigma[npc].Draw("colz");
    cped->SaveAs(Form("PedMap1_DDL%d_PC_%d.eps",nDDL,npc));
    cped->SaveAs(Form("PedMap1_DDL%d_PC_%d.gif",nDDL,npc));
    //cped->SaveAs(Form("PedMap1_DDL%d_PC_%d.pdf",nDDL,npc));
        
    cped2->cd(1); hgPedMapMean1D[npc].Draw();
    cped2->cd(2); hgPedMapSigma1D[npc].Draw();
    cped2->cd(3); hgPedMapMeanSigma[npc].Draw("colz");
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
    
    hgDdlErr = new TH1I(Form("hPedErr_DDL%d",nDDL),Form("DDL(%d) Decoding Errors",nDDL),AliHMPIDRawStream::kSumErr+1,0,AliHMPIDRawStream::kSumErr);
     for(Int_t ilabel=0; ilabel< AliHMPIDRawStream::kSumErr; ilabel++) {
      hgDdlErr->SetStats(0);
      hgDdlErr->GetXaxis()->CenterLabels(kTRUE);
      hgDdlErr->GetXaxis()->SetBinLabel((ilabel+1),Form("%i  %s",ilabel+1,AliHMPIDRawStream::GetErrName(ilabel)));
      }
    hgDdlErr->SetFillColor(5);

    
    Int_t runNumber=-99999;
    Char_t tName[10];
    Int_t  ldcId;
    Int_t  timeStamp;
    Int_t  nEv;Int_t rerr;
     Printf("Start reading Error File DDL: %d ...",nDDL);  
    infile>>tName>>runNumber;
    infile>>tName>>ldcId;
    infile>>tName>>timeStamp;
    infile>>tName>>nEv;  
    Printf("RunNumber: %d, LdcId: %d TimeStamp: %d nEv: %d",runNumber,ldcId,timeStamp,nEv);
    for(Int_t ierr=0;ierr<AliHMPIDRawStream::kSumErr;ierr++)
    {
       infile>>rerr;hgDdlErr->SetBinContent(ierr+1,rerr);
     }
    
    
    TCanvas *cerr = new TCanvas("cped","cped");
    cerr->cd();
    hgDdlErr->Draw();
    cerr->SaveAs(Form("PedError_DDL%d.eps",nDDL)); 
    cerr->SaveAs(Form("PedError_DDL%d.gif",nDDL)); 
      
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void HplotDA()
{
  for(Int_t i=0;i<14;i++) 
  {
   ProcPed(i);  
   ProcErr(i);
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
