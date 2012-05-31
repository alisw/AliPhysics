#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>     

#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1I.h>

#include "AliHMPIDDigit.h"
#include "AliHMPIDRawStream.h"

#endif



TH2F hgPedMapMean[14][6]; 
TH2F hgPedMapSigma[14][6];
TH1F hgPedMapMean1D[14][6]; 
TH1F hgPedMapSigma1D[14][6]; 
TH1F hgPedMapSigma1Db[14][6]; 


TH1I hgDdlErr[14];
TH1I *hgtmp,*hgtmp2;


Int_t fgRunNum;
TFile *fgin[14];
TFile *fgout=0x0;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void Convert(Int_t ddl,Int_t r,Int_t d,Int_t a,Int_t &ch, Int_t &pc, Int_t &px, Int_t &py)
{
  Int_t a2y[6]={3,2,4,1,5,0};//pady for a given address (for single DILOGIC chip)
  
  Int_t ch=ddl/2;
  Int_t tmp=(24-r)/8;              Int_t pc=(ddl%2)? 5-2*tmp:2*tmp; 
  Int_t px=dil*8-pad/6-1;  //flip according to Paolo (26-3-2008)
        tmp=(ddl%2)?r-1:(24-r);   Int_t py=6*(tmp%8)+a2y[a%6];
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PlotPedHisto(Char_t *tag,Int_t runNum, Int_t ddl=0,Int_t row=1,Int_t dil=1,Int_t pad=0)
{
 
  if( 0<=ddl && ddl<=13 && 1<=dil && dil<=10 && 1<=row && row<=24 && 0<=pad && pad<=47) 
  { 
 
  fgin[ddl]->cd();
  
  hgtmp=(TH1I*)(fgin[ddl]->Get(Form("hDDL_%d_Row_%d_Dil_%d_Pad_%d",ddl,row,dil,pad)))->Clone();
  hgtmp2=(TH1I*)(fgin[ddl]->Get(Form("hDDL_%d_Row_%d_Dil_%d_Pad_%d",ddl,row,dil,pad)))->Clone();
  
  TCanvas *c1=new TCanvas(Form("hDDL%s_%d_Row_%d_Dil_%d_Pad_%d",tag,ddl,row,dil,pad),Form("hDDL%s_%d_Row_%d_Dil_%d_Pad_%d",tag,ddl,row,dil,pad));
  
  hgtmp->SetXTitle("ADC");
  hgtmp->SetYTitle("Entries");
  hgtmp->Draw();
  hgtmp->SetAxisRange(0,350);
  hgtmp->Draw("hist same");
  hgtmp2->SetFillColor(5);
  hgtmp2->Draw("hist same");
  hgtmp->Draw("same");
  
  c1->SaveAs(Form("hDDL%s_%d_Row_%d_Dil_%d_Pad_%d.eps",tag,ddl,row,dil,pad));
  c1->SaveAs(Form("hDDL%s_%d_Row_%d_Dil_%d_Pad_%d.gif",tag,ddl,row,dil,pad));
  
  if(fgout!=0x0) 
  {
    fgout->cd();
    hgtmp2->Write();  
  }
  
  hgtmp->Reset();
  hgtmp2->Reset();
  
  }
  else 
  {
   if( ddl < 0 || ddl > 13) {Printf("Not a valid DDL, exiting %d ...",ddl); } 
   if( dil < 1 || dil > 10) {Printf("Not a valid DIL, exiting %d ...",dil); } 
   if( row < 1 || row > 24) {Printf("Not a valid ROW, exiting %d ...",row); } 
   if( pad < 0 || pad > 48) {Printf("Not a valid PAD, exiting %d ...",pad); } 
  }   
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PlotHisto1D(TH1& hin,Char_t* name)
{
  TCanvas *c1=new TCanvas("c1","c1");
  hin.Draw();
  c1->SaveAs(Form("%s.eps",name));
  c1->SaveAs(Form("%s.gif",name));  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void PlotHisto2D(TH2 &hin,Char_t* name,Char_t* opt)
{
 TCanvas *c1=new TCanvas("c1","c1");
 hin.Draw(opt);
 c1->SaveAs(Form("%s.eps",name));
 c1->SaveAs(Form("%s.gif",name));  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ProcPed(Int_t nDDL,Int_t mode)
{

    if(gSystem->IsFileInIncludePath(Form("./HmpidPedDdl%02i.txt",nDDL)))
          ifstream infile(Form("./HmpidPedDdl%02i.txt",nDDL));
    else return;
    if(mode==1) {
    for(Int_t ni=0;ni<6;ni++)
      {
        hgPedMapMean[nDDL][ni]  = new TH2F(Form("hPedMapMean_DDL%d_PC_%d" ,nDDL,ni),         Form("Pedestal: Mean, DDL(0-13): %d PC(0-5): %d;padx;pady" ,nDDL,ni),80,0,80,48,0,48);
        hgPedMapSigma[nDDL][ni] = new TH2F(Form("hPedMapSigma_DDL%d_PC_%d",nDDL,ni),         Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;padx;pady",nDDL,ni),80,0,80,48,0,48);
        hgPedMapMean1D[nDDL][ni]  = new TH1F(Form("hPedMapMean1D_DDL%d_PC_%d" ,nDDL,ni),     Form("Pedestal: Mean, DDL(0-13): %d PC(0-5): %d;pad;Mean" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapSigma1D[nDDL][ni]  = new TH1F(Form("hPedMapSigma1D_DDL%d_PC_%d" ,nDDL,ni),   Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;pad;Sigma" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapSigma1Db[nDDL][ni]  = new TH1F(Form("hPedMapSigma1Db_DDL%d_PC_%d" ,nDDL,ni), Form("Pedestal: Sigma, DDL(0-13): %d PC(0-5): %d;pad;Sigma" ,nDDL,ni),3841,-0.5,3840.5);
        hgPedMapSigma1Db[nDDL][ni]->SetMaximum(2.5);
        }
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
    infile>>tName>>runNumber;fgRunNum=runNumber;
    infile>>tName>>ldcId;
    infile>>tName>>timeStamp;
    infile>>tName>>nEv;  
    infile>>tName>>nDdlEv;  
    infile>>tName>>nBadEv;
    infile>>tName>>nBadEvPer;
    infile>>tName>>nSigCut;
    
    Printf("RunNumber: %d, LdcId: %d TimeStamp: %d nEv: %d",fgRunNum,ldcId,timeStamp,nEv);
    
    while(!infile.eof()){
      infile>>dec>>r>>d>>a>>mean>>sigma>>hex>>hard;
      switch(mode)
      {
      case 1:
        Convert(nDDL,r,d,a,ch,pc,px,py); 
        hgPedMapMean[nDDL][pc].Fill(px,py,mean);  
        hgPedMapSigma[nDDL][pc].Fill(px,py,sigma);
        hgPedMapMean1D[nDDL][pc].Fill(cnt%3840,mean);
        hgPedMapSigma1D[nDDL][pc].Fill(cnt%3840,sigma);
        hgPedMapSigma1Db[nDDL][pc].Fill(cnt%3840,sigma);
      break;
      case 2:
        if(sigma>3 && sigma!=1000)  PlotPedHisto("Sigma",runNumber,nDDL,r,d,a);    
      break;
      case 3:
        if(sigma==1000)  PlotPedHisto("NoEntry",runNumber,nDDL,r,d,a);    
      break;
      
      }
        
      cnt++;      
    }
  infile.close();
  Printf("Stop reading DDL: %d ...",nDDL);
  
  /* fill the overall histos */
  
  if(mode==1) 
  {
    for(Int_t npc=0;npc<6;npc++) 
      {
          /* plot pedestal mean values */
         if(hgPedMapMean[nDDL][npc].GetEntries()!=0) {
           hgPedMapMean[nDDL][npc].SetStats(0);
           PlotHisto2D(hgPedMapMean[nDDL][npc],Form("Run%d_PedMeanA_DDL%d_PC%d",fgRunNum,nDDL,npc),"colz");
           PlotHisto2D(hgPedMapMean[nDDL][npc],Form("Run%d_PedMeanB_DDL%d_PC%d",fgRunNum,nDDL,npc),"surf1");
         }
         /* plot pedestal sigma values */
         if(hgPedMapSigma[nDDL][npc].GetEntries()!=0) {
           hgPedMapSigma[nDDL][npc].SetStats(0);
           PlotHisto2D(hgPedMapSigma[nDDL][npc],Form("Run%d_PedSigmaA_DDL%d_PC%d",fgRunNum,nDDL,npc),"colz");
           PlotHisto2D(hgPedMapSigma[nDDL][npc],Form("Run%d_PedSigmaB_DDL%d_PC%d",fgRunNum,nDDL,npc),"surf1");
         }
         /* plot pedestal sigma values */
         if(hgPedMapSigma[nDDL][npc].GetEntries()!=0) {
           hgPedMapSigma[nDDL][npc].SetStats(0);
           hgPedMapSigma[nDDL][npc].SetMaximum(2.5);
           PlotHisto2D(hgPedMapSigma[nDDL][npc],Form("Run%d_PedSigmaC_DDL%d_PC%d",fgRunNum,nDDL,npc),"colz");
           PlotHisto2D(hgPedMapSigma[nDDL][npc],Form("Run%d_PedSigmaD_DDL%d_PC%d",fgRunNum,nDDL,npc),"surf1");
         }
         /* plot pedestal mean values */
         if(hgPedMapMean1D[nDDL][npc].GetEntries()!=0) {
            hgPedMapMean1D[nDDL][npc].SetStats(0);
           PlotHisto1D(hgPedMapMean1D[nDDL][npc],Form("Run%d_PedMean1DA_DDL%d_PC%d",fgRunNum,nDDL,npc));
         }
         /* plot pedestal sigma values */
         if(hgPedMapSigma1D[nDDL][npc].GetEntries()!=0) {
           hgPedMapSigma1D[nDDL][npc].SetStats(0);
           PlotHisto1D(hgPedMapSigma1D[nDDL][npc],Form("Run%d_PedSigma1DA_DDL%d_PC%d",fgRunNum,nDDL,npc));
         }
         /* plot pedestal sigma values */
         if(hgPedMapSigma1Db[nDDL][npc].GetEntries()!=0) {
           hgPedMapSigma1Db[nDDL][npc].SetStats(0);
           PlotHisto1D(hgPedMapSigma1D[nDDL][npc],Form("Run%d_PedSigma1DA_DDL%d_PC%d",fgRunNum,nDDL,npc));
         }       
      }  
  }
  
       
}//ProcPed()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void ProcErr(Int_t nDDL,Int_t mode)
{
  if(gSystem->IsFileInIncludePath(Form("./HmpidErrorsDdl%02i.txt",nDDL)))
          ifstream infile(Form("./HmpidErrorsDdl%02i.txt",nDDL));
    else return;
    if(mode==1)
    {
    hgDdlErr[nDDL] = new TH1I(Form("hPedErr_DDL%d",nDDL),Form("DDL(%d) Decoding Errors",nDDL),AliHMPIDRawStream::kSumErr+1,-0.5,AliHMPIDRawStream::kSumErr+0.5);
     for(Int_t ilabel=0; ilabel< AliHMPIDRawStream::kSumErr; ilabel++) {
      hgDdlErr[nDDL].SetStats(0);
      hgDdlErr[nDDL].GetXaxis()->CenterLabels(kTRUE);
      hgDdlErr[nDDL].GetXaxis()->SetBinLabel((ilabel+1),Form("%i  %s",ilabel+1,AliHMPIDRawStream::GetErrName(ilabel)));
      hgDdlErr[nDDL].SetYTitle("Error #");
     }
    hgDdlErr[nDDL].SetFillColor(5);
  }  
    
    Int_t runNumber=-99999;
    Char_t tName[10];
    Int_t  ldcId;
    Int_t  timeStamp;
    Int_t  nEv,nDdlEv;
    Int_t rerr;Int_t nBadEv;Float_t nBadEvPer;
    Int_t row,dil,pad,nzero;
    Int_t ch,pc,px,py;
    Printf("Start reading Error File DDL: %d ...",nDDL);  
    infile>>tName>>runNumber;fgRunNum=runNumber;
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
    
    switch(mode)
    {
      case 1:          
          TCanvas *cerr = new TCanvas("cped","cped");
          cerr.cd();
          hgDdlErr[nDDL]->Draw();
          cerr->SaveAs(Form("PedError_DDL%d.eps",nDDL)); 
          cerr->SaveAs(Form("PedError_DDL%d.gif",nDDL)); 
    
        break;
        
      case 2:
        while(!infile.eof()){
        infile>>dec>>row>>dil>>pad>>nzero; 
        if( 0<=nDDL && nDDL<=13 && 1<=dil && dil<=10 && 1<=row && row<=24 && 0<=pad && pad<=47) 
          { 
           PlotPedHisto("ZeroQ",runNumber,nDDL,row,dil,pad);
         }
       } 
        
    break;
      
    default:
        break;  
    }
   infile.close();     
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
void HplotDA(Int_t runNumber)
{
  gStyle->SetPalette(1);
  fgRunNum=runNumber;
  
  Printf("********************************************************************");
  Printf("********************** HplotDA *************************************");
  Printf("***                                                              ***");
  Printf("*** Choose from the following options:                           ***");
  Printf("***                                                              ***");
  Printf("*** 1.) Draw errors for all DDL                                  ***");
  Printf("*** 2.) Draw pedestals for all DDL                               ***");
  Printf("*** 3.) Draw ADC histo for pads with sigma > 3                   ***");
  Printf("*** 4.) Draw ADC histo for pads with zero charge                 ***");
  Printf("*** 5.) Draw ADC histo for pads with no valid charge             ***");
  Printf("*** 6.) Run all                                                  ***");
  Printf("***                                                              ***");
  Printf("********************************************************************");
  Printf("*** Please select: ");
  
  Int_t set=0;
  cin>>set;
  
  
 for(Int_t i=0;i<14;i++) 
    {
    switch(set){
  
    case 1:  
     ProcErr(i,1);
     break;
     case 2:
      ProcPed(i,1);  
     break;
     case 3:
       fgin[i]=new TFile(Form("Run%d_DDL%d.root",runNumber,i),"read");
       ProcPed(i,2);  
       fgin[i]->Close();
       break;
     case 4:
       fgin[i]=new TFile(Form("Run%d_DDL%d.root",runNumber,i),"read");
       ProcErr(i,2);
       fgin[i]->Close();
     break;  
     case 5:
       fgin[i]=new TFile(Form("Run%d_DDL%d.root",runNumber,i),"read");
       ProcPed(i,3);  
       fgin[i]->Close();
     break; 
     case 6:
      if(i==0) fgout=new TFile(Form("SummaryOfRun%d.root",fgRunNum),"recreate");
  
      ProcErr(i,1);
      ProcPed(i,1);
      fgin[i]=new TFile(Form("Run%d_DDL%d.root",runNumber,i),"read");
      ProcPed(i,2);  
      ProcErr(i,2);
      ProcPed(i,3);  
      fgin[i]->Close(); 
       
      if(i==13) {
        
        for(Int_t ipc=0;ipc<6;ipc++) 
        {
          hgPedMapMean[i][ipc].Write(); 
          hgPedMapSigma[i][ipc].Write(); 
          hgPedMapMean1D[i][ipc].Write(); 
          hgPedMapSigma1D[i][ipc].Write(); 
          hgPedMapSigma1Db[i][ipc].Write(); 
        }
        hgDdlErr[i].Write();
      }
     break;
     default:
         Printf("Not a valid selection bye-bye....");
     break;
    }
  }
 
  if(fgout!=0x0)fgout->Close();    
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
