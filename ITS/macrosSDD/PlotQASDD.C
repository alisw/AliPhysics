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
//#include "AliRecoParam.h"
#endif

enum rawexpert_t{pattern,patternnorm,layer3norm,layer4norm,rawtot};
enum nonrawexpert_t{layer3,layer4,nonrawtot};

enum rpexpert_t{rlocaldistro,rlayer3norm,rlayer4norm,rmodpattern,rmodpatternnorm,rphilayer3,rphilayer4,rrlayer3,rrlayer4,rtot};
enum rnonexpert_t{rphizlayer3,rphizlayer4,rglobalrz,rglobalxy,rchargelayer3,rchargelayer4,rtimelayer3,rtimelayer4,rlayer3,rlayer4,nrtot};

enum canvsname_t{canvname0,canvname1,canvname2,canvname3,canvname4,canvname5,canvname6,canvname7,canvname8,canvname9,canvname10,canvname11,canvname12};

void PlotQASDD(Char_t fileName[100]="File.QA.111333.2010.LHC09b.pass2.root",Char_t eventspecie[25]="LowMultiplicity")
{
  const TString Rawsexpertname[]={"SDDModPattern","SDDModPatternNORM","SDDphizL3NORM","SDDphizL4NORM"};//4
  const TString Rawsnonexpertname[]={"SDDphizL3","SDDphizL4"};//2

  const TString RecPointsexpertname[]={"SDDLocalCoordDistrib","SDDModPatternL3RPNORM","SDDModPatternL4RPNORM","SDDModPatternRP","SDDModPatternRPNORM","SDDphidistrib_Layer3","SDDphidistrib_Layer4","SDDrdistrib_Layer3","SDDrdistrib_Layer4"};//9
  const TString RecPointsnonexpertname[]={"SDDGlobalCoordDistribL3PHIZ","SDDGlobalCoordDistribL4PHIZ","SDDGlobalCoordDistribRZ","SDDGlobalCoordDistribYX","SDDLay3TotCh","SDDLay4TotCh","SDDdrifttime_Layer3","SDDdrifttime_Layer4","SDDModPatternL3RP","SDDModPatternL4RP"};//1

  const TString canvassavedname[]={"RawLayers","Rawpatterns","RawLayersNORM","RecPointsphiz","RecPointsGlobalDistributions","RecPointsCharge","RecPointsDriftTime","RecPointsLocalDistribution","RecPointsLayers","RecPointsLayersNORM","RecPointsLayersPatterns","RecPointsPhiDistribution","RecPointsRDistribution"};//13

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
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

  TCanvas *canvas1 = new TCanvas("canvas1","SDD QA Plot",1000,600);
//   canvas1->Divide(2,1);
//   TCanvas *canvas3 = new TCanvas("canvas3","RawLayersNORM",1000,600);
//   canvas3->Divide(2,1);
//   TCanvas *canvas2 = new TCanvas("canvas2","Rawpatterns",1000,600);
//   canvas2->Divide(1,2);

//   TCanvas *canvas4 = new TCanvas("canvas4","RecPointsphiz",1000,600);
//   canvas4->Divide(2,1);
//   TCanvas *canvas5 = new TCanvas("canvas5","RecPointsGlobalDistributions",1000,600);
//   canvas5->Divide(2,1);
//   TCanvas *canvas6 = new TCanvas("canvas6","RecPointsCharge",1000,600);
//   canvas6->Divide(1,2);
//   TCanvas *canvas7 = new TCanvas("canvas7","RecPointsDriftTime",1000,600);
//   canvas7->Divide(1,2);

//   TCanvas *canvas8 = new TCanvas("canvas8","RecPointsLocalDistribution",1000,600);
//   TCanvas *canvas9 = new TCanvas("canvas9","RecPointsLayers",1000,600);
//   canvas9->Divide(2,1);
  
//   TCanvas *canvas10 = new TCanvas("canvas10","RecPointsLayersNORM",1000,600);
//   canvas10->Divide(2,1);
//   TCanvas *canvas11 = new TCanvas("canvas11","RecPointsLayersPatterns",1000,600);
//   canvas11->Divide(1,2);
//   TCanvas *canvas12 = new TCanvas("canvas12","RecPointsPhiDistribution",1000,600);
//   canvas12->Divide(2,1);
//   TCanvas *canvas13 = new TCanvas("canvas13","RecPointsRDistribution",1000,600);
//   canvas13->Divide(2,1);
  

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
      Bool_t updatecanvas[2];
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      for(Int_t iraws=0;iraws<nonrawtot;iraws++){//non expert raws
	canvas1->cd(iraws+1);

	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsnonexpertname[iraws].Data());
	printf("histo name %s ",Rawsnonexpertname[iraws].Data());
	historaw=(TH2D*)mergedfile.Get(histoname);
	if(historaw){printf("...Found\n");historaw->DrawCopy("colz");}
	else{updatecanvas[iraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw=NULL;
      }//end for	

      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  sprintf(namecanvas,"%s.eps",canvassavedname[canvname0].Data());
	  canvas1->SaveAs(namecanvas);
	  canvas1->Print("SDDQAPlot.ps");
	}
      //else{delete canvas1; canvas1=NULL;}
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	canvas1->Divide(1,2);
	printf("2 \n");
	for(Int_t inraws=0;inraws<2;inraws++){//non expert raws
	  canvas1->cd(inraws+1);
	  sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[inraws].Data());
	  printf("histo name %s ",Rawsexpertname[inraws].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);
	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[inraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}//end for
	if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	  {
	    canvas1->Update();
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname1].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas2; canvas2=NULL;}
     
	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	printf("3\n");
	canvas1->Clear();
	canvas1->Divide(2,1);
	for(Int_t iraws=2;iraws<rawtot;iraws++){//non expert raws
	  canvas1->cd(iraws-1);
	  
	  sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[iraws].Data());
	  printf("histo name %s",Rawsexpertname[iraws].Data());
	  historaw=(TH2D*)mergedfile.Get(histoname);
	  if(historaw){printf("...Found\n");historaw->DrawCopy("colz");}else{updatecanvas[iraws-2]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");} 
	  historaw=NULL;
	}//end for
	if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	  {
	    canvas1->Update();
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname2].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas3; canvas3=NULL;}
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
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname3].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//	else{delete canvas4; canvas4=NULL;}
       

	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	canvas1->Divide(2,1);
	printf("5\n");
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
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname4].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas5; canvas5=NULL;}
     

	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	canvas1->Divide(1,2);
	printf("6\n");
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
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname5].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas6; canvas6=NULL;}
	
	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	canvas1->Divide(1,2);
	printf("7\n");
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
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname6].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas7; canvas7=NULL;}

	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	printf("8\n");
	//canvas1->Divide(2,1);
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[rlocaldistro].Data());
	printf("histo name %s",RecPointsexpertname[rlocaldistro].Data());
	histodraw=(TH2F*)mergedfile.Get(histoname);

	if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[0]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	histodraw=NULL;
	if(updatecanvas[0]==kTRUE){
	  canvas1->Update();
	  sprintf(namecanvas,"%s.eps",canvassavedname[canvname7].Data());
	  canvas1->SaveAs(namecanvas);
	  canvas1->Print("SDDQAPlot.ps"); 
	}
	//else{delete canvas8; canvas8=NULL;}
       
	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	canvas1->Divide(2,1);
	printf("9\n");
	for(Int_t i=8;i<10;i++)
	  {
	    canvas1->cd(i-7);

	    sprintf(histoname,"ITS/RecPoints/%s/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[i].Data());
	    printf("histo name %s",RecPointsnonexpertname[i].Data());
	    histodraw=(TH2F*)mergedfile.Get(histoname);

	    if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[i-8]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	    histodraw=NULL;
	  }
	if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	  {
	    canvas1->Update();
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname8].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas9; canvas9=NULL;}
	    for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	    canvas1->Clear();
	    canvas1->Divide(2,1);
	    printf("10\n");
	for(Int_t i=1;i<3;i++)
	  {
	    canvas1->cd(i);

	    sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	    histodraw=(TH2F*)mergedfile.Get(histoname);
	    printf("histo name %s",RecPointsexpertname[i].Data());

	    if(histodraw){printf("...Found\n");histodraw->DrawCopy("colz");}else{updatecanvas[i-1]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	    histodraw=NULL;
	  }
	if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	  {
	    canvas1->Update();
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname9].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas10; canvas10=NULL;}

	for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	canvas1->Clear();
	canvas1->Divide(1,2);
	printf("11\n");
	for(Int_t i=3;i<5;i++)
	  {
	    canvas1->cd(i-2);

	    sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	    printf("histo name %s",RecPointsexpertname[i].Data());
	    historaw2=(TH1F*)mergedfile.Get(histoname);

	    if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-3]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	    historaw2=NULL;
	  }
	if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	  {
	    canvas1->Update();
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname10].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas11; canvas11=NULL;}

	    for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	    canvas1->Clear();
	    canvas1->Divide(2,1);
	    printf("12\n");
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
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname11].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	//else{delete canvas12; canvas12=NULL;}

	    for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
	    canvas1->Clear();
	    canvas1->Divide(2,1);
	    printf("13\n");
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
	    sprintf(namecanvas,"%s.eps",canvassavedname[canvname12].Data());
	    canvas1->SaveAs(namecanvas);
	    canvas1->Print("SDDQAPlot.ps"); 
	  }
	canvas1->Print("SDDQAPlot.ps]"); 
	//else{delete canvas13; canvas13=NULL;}	
	    directory2=NULL;
      }//end directory
  }//end for   

  delete  directory;
  directory=NULL;
  delete directory2;
  directory2=NULL;

}//end macro


