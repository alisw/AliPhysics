#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TFileMerger.h>
#include <TAlienFile.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TClass.h>
#include <TKey.h>
#include <TFile.h>
#include <TFileMerger.h>
#include <TAlienFile.h>
//#include <TExec.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TClass.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include "AliRecoParam.h"
#endif

enum rawexpert_t{pattern,patternnorm,layer3norm,layer4norm,layer3occ,layer4occ,rawtot};
enum nonrawexpert_t{layer3,layer4,nonrawtot};

enum rpexpert_t{rlocaldistro,rlayer3norm,rlayer4norm,rmodpattern,rmodpatternnorm,rphilayer3,rphilayer4,rrlayer3,rrlayer4,occlayer3,occlayer4,recraw2dlayer3,recraw2dlayer4,recrawlayer3,recrawlayer4,dedxlayer3,dedxlayer4,rtot};
enum rnonexpert_t{rphizlayer3,rphizlayer4,rglobalrz,rglobalxy,rchargelayer3,rchargelayer4,rtimelayer3,rtimelayer4,rlayer3,rlayer4,nrtot};

enum canvsname_t{canvname0,canvname1,canvname2,canvname3,canvname4,canvname5,canvname6,canvname7,canvname8,canvname9,canvname10,canvname11,canvname12,canvname13,canvname14,canvname15,canvname16,canvname17,canvname18,canvname19,canvname20,canvname21,canvname22,canvname23};//



void PlotQASDD(Char_t fileName[100]="File.QA.111333.2010.LHC09b.pass2.root",Char_t eventspecie[25]="LowMultiplicity",Bool_t kDoEps=kFALSE)
{
  const TString Rawsexpertname[]={"SDDModPattern","SDDModPatternNORM","SDDphizL3","SDDphizL4","SDDL3_RelativeOccupancy","SDDL4_RelativeOccupancy"};//6
  const TString Rawsnonexpertname[]={"SDDphizL3NORM","SDDphizL4NORM"};//2

  const TString RecPointsexpertname[]={"SDDLocalCoordDistrib","SDDModPatternL3RP","SDDModPatternL4RP","SDDModPatternRP","SDDModPatternRPNORM","SDDphidistrib_Layer3","SDDphidistrib_Layer4","SDDrdistrib_Layer3","SDDrdistrib_Layer4","SDDL3_RelativeOccupancy","SDDL4_RelativeOccupancy","SDDL3_Rec2Raw_2D","SDDL4_Rec2Raw_2D","SDDL3_Rec2Raw","SDDL4_Rec2Raw","SDDL3_dedx","SDDL4_dedx"};//17
  const TString RecPointsnonexpertname[]={"SDDGlobalCoordDistribL3PHIZ","SDDGlobalCoordDistribL4PHIZ","SDDGlobalCoordDistribRZ","SDDGlobalCoordDistribYX","SDDLay3TotCh","SDDLay4TotCh","SDDdrifttime_Layer3","SDDdrifttime_Layer4","SDDModPatternL3RPNORM","SDDModPatternL4RPNORM"};//10

  const TString canvassavedname[]={"RawLayers","Rawpatterns","RawLayersNORM","RawLayersRelOcc","RawLayersRelOccBoth","RecPointsphiz","RecPointsGlobalDistributions","RecPointsCharge","RecPointsDriftTime","RecPointsLocalDistribution","RecPointsLayers","RecPointsLayersNORM","RecPointsLayersPatterns","RecPointsPhiDistribution","RecPointsRDistribution","RecPointsChargeBothLayers","RecPointsDriftTimeBothLayers","RecPointsRelLayersOcc","RecPointsLayersRelOccBoth","RecPointsLayersRecRaw2D","RecPointsLayersRecRaw","RecPointsLayersRecRaw","RecPointsLayersRecRawBoth","RecPointsLayersdedx","RecPointsLayersdedxBoth"};//15

  Char_t trendfile[25]="FileQAtrend.root";


  gStyle->SetPalette(1);
  //  gStyle->SetOptStat("ne");
  if(gSystem->Exec(Form("ls %s >/dev/null 2>&1", fileName))!=0)
    {
      printf(" No file created --- Exit");
      return;
    }
  TFile mergedfile(fileName);
  if(mergedfile.GetKey("ITS")==0x0){
    printf("In this run ITS QA has not been executed.\n\nExit macro \n\n");
    return;
  }
  Char_t filepath[100];
  Char_t namecanvas[50];
  TDirectory *directory=NULL;
  TDirectory *directory2=NULL;
  Char_t histoname[200];
  TH2D *historaw=NULL;
  TH1F *historaw2=NULL;
  TH2F *histodraw=NULL;
  
  TH1F *histo2save[9];
  for(Int_t i=0;i<9;i++){histo2save[i]=NULL;}
  
  TF1 *f1 = new TF1("f1","1",-100000,100000); // this is the function defined to multiply histograms
  Char_t layer[10];
  float fmax=0.;
  float fmaxold=0.;
  float fmaxmargin;
  
  Int_t cannum=0;
  Int_t histo2savenumber[9]={4,5,6,7,7,8,3,5,1};
  TCanvas *canvas1 = new TCanvas("canvas1","SDD QA Plot",1000,600);
  //   TLegend *legend=new TLegend(0.83,0.8,0.97,0.7);
  TLegend *legend=new TLegend(0.81,0.895,0.95,0.995);
  Bool_t updatecanvas[2];
  for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
  
  // -------------This part is to read the number of chunks that were merged
  Float_t fChnknmbr=0.;
  FILE * pChunkNumber = fopen ("ChunkNumber.txt","r");
  if(pChunkNumber){
  Int_t rv=fscanf(pChunkNumber, "%f", &fChnknmbr);
  fclose (pChunkNumber);
  }
  else fChnknmbr=1.;

  gStyle->SetPalette(1);
  float fCNinv=1./fChnknmbr;
  //   printf("\n====================>%f\n\n", fCNinv);
  
  if(!mergedfile.IsOpen()){return;}else{printf("file is open\n");}
  for(Int_t ispecie=0;ispecie<AliRecoParam::kNSpecies;ispecie++){
    //__________________________________________________________________
    //raw data
    sprintf(filepath,"ITS/Raws/%s",AliRecoParam::GetEventSpecieName(ispecie));
    printf("%s",filepath);
    TString especie(filepath);
    if(!especie.Contains(eventspecie)){printf("...Found and Skipped\n");continue;}
    canvas1->Print("SDDQAPlot.ps[");
    directory=(TDirectory*)mergedfile.Get(filepath);
    if(!directory){printf("...Not Found\n ");}//faccio l'istogramma
    else{
      printf("...Found: The histograms of this EventSpecie will be displayed\n");
      printf("1\n");
      canvas1->Clear();
      canvas1->Divide(2,1);
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      for(Int_t iraws=2;iraws<4;iraws++){//non expert raws
	canvas1->cd(iraws-1);
	
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[iraws].Data());
	printf("histo name %s ",Rawsexpertname[iraws].Data());
	historaw=(TH2D*)mergedfile.Get(histoname);
	if(historaw){printf("...Found\n");historaw->DrawCopy("colz");}
	else{updatecanvas[iraws-2]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw=NULL;
      }//end for	
      
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps");
	}
      //else{delete canvas1; canvas1=NULL;}
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("2 \n");
      cannum++;
      for(Int_t inraws=0;inraws<2;inraws++){//non expert raws
	canvas1->cd(inraws+1);
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[inraws].Data());
	printf("histo name %s ",Rawsexpertname[inraws].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	
	
	// 	  -----------------This is the part were I divide between the number of chunks to normalize the histogram----------------
	
	
	if(historaw2)
	  {
	    printf("...Found\n");
	    if (inraws==1)
	      {
		historaw2->Multiply(f1,fCNinv);
		
	      }
	    historaw2->DrawCopy();
	  }
	else{updatecanvas[inraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas2; canvas2=NULL;}
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      printf("3\n");
      cannum++;
      canvas1->Clear();
      canvas1->Divide(2,1);
      for(Int_t iraws=0;iraws<2;iraws++){//non expert raws
	canvas1->cd(iraws+1);
	
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsnonexpertname[iraws].Data());
	printf("histo name %s",Rawsnonexpertname[iraws].Data());
	historaw=(TH2D*)mergedfile.Get(histoname);
	if(historaw){
	  printf("...Found\n");
	  historaw->Multiply(f1,fCNinv);
	  historaw->DrawCopy("colz");
	}
	else
	  {updatecanvas[iraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");} 
	historaw=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas3; canvas3=NULL;}
      
      //--------------------- new plots
      printf("================= 4\n");
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      cannum++;
      canvas1->Clear();
      canvas1->Divide(2,1);
      for(Int_t iraws=4;iraws<rawtot;iraws++){//non expert raws
	canvas1->cd(iraws-3);
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[iraws].Data());
	printf("histo name %s",Rawsexpertname[iraws].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->DrawCopy();
	  }else
	  {
	    updatecanvas[iraws-4]=kFALSE;
	    printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");
	  } 
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas3; canvas3=NULL;}
      
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("14\n");
      cannum++;
      for (Int_t i=4;i<rawtot;i++){
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==4){
	    fmaxold=fmax;}
	  if (i!=4){
	    if(fmaxold<fmax){
	      fmaxold=fmax;
	    }
	  }
	}
	historaw2=NULL;
      }
      fmaxmargin=1.1*fmaxold;
      for(Int_t irrpp=4;irrpp<rawtot;irrpp++){//non expert raws
	
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[irrpp].Data());
	printf("histo name %s",Rawsexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-1);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("Raw Relative Occupancy");
	    if (irrpp==4) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=4) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-4]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      
    }//end if directory
    
    
    
    //_______________________________________________________________________________
    //rec point
    sprintf(filepath,"ITS/RecPoints/%s",AliRecoParam::GetEventSpecieName(ispecie));
    printf("%s\n",filepath);
    directory2=(TDirectory*)mergedfile.Get(filepath);
    if(directory2){
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("4\n");
      cannum++;
      for(Int_t irp=0;irp<2;irp++){//non expert rec point
	canvas1->cd(irp+1);
	
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[irp].Data());
	printf("histo name %s ",RecPointsnonexpertname[irp].Data());
	histodraw=(TH2F*)mergedfile.Get(histoname);
	if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[irp]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	histodraw=NULL;
      }
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //	else{delete canvas4; canvas4=NULL;}
      
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("5\n");
      cannum++;
      for(Int_t irrp=2;irrp<4;irrp++){//non expert raws
	canvas1->cd(irrp-1);
	
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[irrp].Data());
	printf("histo name %s",RecPointsnonexpertname[irrp].Data() );
	histodraw=(TH2F*)mergedfile.Get(histoname);
	if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[irrp-2]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	histodraw=NULL;
      }
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas5; canvas5=NULL;}
      
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("6\n");
      cannum++;
      for(Int_t irrpp=4;irrpp<6;irrpp++){//non expert raws
	canvas1->cd(irrpp-3);
	
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[irrpp].Data());
	printf("histo name %s",RecPointsnonexpertname[irrpp].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	
	if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[irrpp-4]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas6; canvas6=NULL;}
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("7\n");
      cannum++;
      for(Int_t irrpp=6;irrpp<8;irrpp++){//non expert raws
	canvas1->cd(irrpp-5);
	//	  printf("histo name\n");
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[irrpp].Data());
	printf("histo name %s",RecPointsnonexpertname[irrpp].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[irrpp-6]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas7; canvas7=NULL;}

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      printf("8\n");
      cannum++;
      //canvas1->Divide(2,1);
      sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[rlocaldistro].Data());
      printf("histo name %s",RecPointsexpertname[rlocaldistro].Data());
      histodraw=(TH2F*)mergedfile.Get(histoname);

      if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[0]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
      histodraw=NULL;
      if(updatecanvas[0]==kTRUE){
	canvas1->Update();
	if(kDoEps){
	  sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	  canvas1->SaveAs(namecanvas);
	}
	canvas1->Print("SDDQAPlot.ps"); 
      }
      //else{delete canvas8; canvas8=NULL;}
       
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("9\n");
      cannum++;
      for(Int_t i=1;i<3;i++)
	{
	  canvas1->cd(i);

	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  histodraw=(TH2F*)mergedfile.Get(histoname);

	  if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[i-1]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  histodraw=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas9; canvas9=NULL;}
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("10\n");
      cannum++;
      for(Int_t i=8;i<10;i++)
	{
	  canvas1->cd(i-7);

	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[i].Data());
	  histodraw=(TH2F*)mergedfile.Get(histoname);
	  printf("histo name %s",RecPointsnonexpertname[i].Data());

	  if(histodraw){
	    printf("...Found\n");histodraw->Multiply(f1,fCNinv);
	    histodraw->DrawCopy("colz");
	  }else{updatecanvas[i-8]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  histodraw=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas10; canvas10=NULL;}

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("11\n");
      cannum++;
      for(Int_t i=3;i<5;i++)
	{
	  canvas1->cd(i-2);

	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2)
	    {
	      printf("...Found\n");
	      if (i==4){historaw2->Multiply(f1,fCNinv);}
	      historaw2->DrawCopy();
	    }
	  else{updatecanvas[i-3]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas11; canvas11=NULL;}

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("12\n");
      cannum++;
      for(Int_t i=5;i<7;i++)
	{
	  canvas1->cd(i-4);

	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-5]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //else{delete canvas12; canvas12=NULL;}

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("13\n");
      cannum++;
      for(Int_t i=7;i<9;i++)
	{
	  canvas1->cd(i-6);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-7]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
	  

      //superimposed
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("14\n");
      cannum++;
      for (Int_t i=4;i<6;i++){
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==4){
	    fmaxold=fmax;}
	  if (i!=4){
	    if(fmaxold<fmax){
	      fmaxold=fmax;
	    }
	  }
	}
	historaw2=NULL;
      }
      fmaxmargin=1.1*fmaxold;
      for(Int_t irrpp=4;irrpp<6;irrpp++){//non expert raws
    
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[irrpp].Data());
	printf("histo name %s",RecPointsnonexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-1);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("Charge");
	    historaw2->SetFillColor(0);
	    if (irrpp==4) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=4) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-4]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
	  
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("15\n");
      cannum++;
      for (Int_t i=6;i<8;i++){
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==6){
	    fmaxold=fmax;}
	  if (i!=6){
	    if(fmaxold<fmax){
	      fmaxold=fmax;}}}
	fmaxmargin=1.1*fmaxold;
      }
      for(Int_t irrpp=6;irrpp<8;irrpp++){//non expert raws
    
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[irrpp].Data());
	printf("histo name %s",RecPointsnonexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-3);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("Drift Time");
	    historaw2->SetFillColor(0);
	    if (irrpp==6) {historaw2->SetStats(0); historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=6) {historaw2->SetStats(0); historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-6]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	printf("%s\n%s\n",historaw2->GetName(),historaw2->GetTitle());
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}

      //------------------------------------------- new plot

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("13\n");
      cannum++;
      for(Int_t i=9;i<11;i++)
	{
	  canvas1->cd(i-8);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-9]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
	  
      //------------------------------------
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("14\n");
      cannum++;
      for (Int_t i=9;i<11;i++){
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==4){
	    fmaxold=fmax;}
	  if (i!=4){
	    if(fmaxold<fmax){
	      fmaxold=fmax;
	    }
	  }
	}
	historaw2=NULL;
      }
      fmaxmargin=1.1*fmaxold;
      for(Int_t irrpp=9;irrpp<11;irrpp++){//non expert raws
	  
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[irrpp].Data());
	printf("histo name %s",RecPointsexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-6);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("RecPoint Relative Occupancy");
	    if (irrpp==9) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=9) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-9]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}


      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("13\n");
      cannum++;
      for(Int_t i=11;i<13;i++)
	{
	  canvas1->cd(i-10);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw=(TH2D*)mergedfile.Get(histoname);

	  if(historaw){printf("...Found\n");historaw->DrawCopy("colz");}else{updatecanvas[i-11]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
	  


      //-----------------
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("13\n");
      cannum++;
      for(Int_t i=13;i<15;i++)
	{
	  canvas1->cd(i-12);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-13]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      //--------------------------------------------

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("14\n");
      cannum++;
      for (Int_t i=13;i<15;i++){
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==13){
	    fmaxold=fmax;}
	  if (i!=13){
	    if(fmaxold<fmax){
	      fmaxold=fmax;
	    }
	  }
	}
	historaw2=NULL;
      }
      fmaxmargin=1.1*fmaxold;
      for(Int_t irrpp=13;irrpp<15;irrpp++){//non expert raws
    
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[irrpp].Data());
	printf("histo name %s",RecPointsexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-10);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("Rec2Raw Ratio");
	    if (irrpp==13) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=13) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-13]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}

	  
      //--------------------------------------
	  
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("13\n");
      cannum++;
      for(Int_t i=15;i<17;i++)
	{
	  canvas1->cd(i-14);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);
		
	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-15]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
	  

      //--------------------------------------------

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("14\n");
      cannum++;
      for (Int_t i=15;i<17;i++){
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==15){
	    fmaxold=fmax;}
	  if (i!=15){
	    if(fmaxold<fmax){fmaxold=fmax;}
	  }
	}
	historaw2=NULL;
      }
      fmaxmargin=1.1*fmaxold;
      for(Int_t irrpp=15;irrpp<17;irrpp++){//non expert raws
    
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[irrpp].Data());
	printf("histo name %s",RecPointsexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-12);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("RecPoint dEdx");
	    if (irrpp==15) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=15) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-15]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print("SDDQAPlot.ps"); 
	}
      canvas1->Print("SDDQAPlot.ps]"); 


      for(Int_t isave=0;isave<9;isave++){
	if(isave<4)sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[histo2savenumber[isave]].Data());
	if(isave>4&&isave<7)sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[histo2savenumber[isave]].Data());
	if(isave>7)sprintf(histoname,"ITS/Raws/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[histo2savenumber[isave]].Data());
	printf("file2save name:\t %s\n",histoname);
	histo2save[isave]=(TH1F*)mergedfile.Get(histoname);
      }

      TFile file2savefortrend(trendfile,"recreate");
      file2savefortrend.cd();
      for(Int_t iss=0;iss<9;iss++){printf("Saved %d\n",iss); histo2save[iss]->Write();}
      file2savefortrend.Close();
	
      //else{delete canvas13; canvas13=NULL;}	
      //directory2=NULL;
    }//end directory



  }//end for   




  delete  directory;
  directory=NULL;
  delete directory2;
  directory2=NULL;

}//end macro


