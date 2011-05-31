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
#include <TString.h>
#include "AliRecoParam.h"
#endif

enum rawexpert_t{pattern,patternnorm,layer3norm,layer4norm,layer3occ,layer4occ,rawcheck,rawtot};
enum nonrawexpert_t{layer3,layer4,calibl3,calibl4,nonrawtot};

enum rpexpert_t{rlocaldistro,rlayer3norm,rlayer4norm,rmodpattern,rmodpatternnorm,rphilayer3,rphilayer4,rrlayer3,rrlayer4,occlayer3,occlayer4,recraw2dlayer3,recraw2dlayer4,recrawlayer3,recrawlayer4,dedxlayer3,dedxlayer4,rlayer3,rlayer4,recpointcheck,rtot};
enum rnonexpert_t{rphizlayer3,rphizlayer4,rglobalrz,rglobalxy,rchargelayer3,rchargelayer4,rtimelayer3,rtimelayer4,nrtot};

enum canvsname_t{canvname0,canvname1,canvname2,canvname3,canvname4,canvname5,canvname6,canvname7,canvname8,canvname9,canvname10,canvname11,canvname12,canvname13,canvname14,canvname15,canvname16,canvname17,canvname18,canvname19,canvname20,canvname21,canvname22,canvname23,canvasname24,canvasname25,canvasname26};//27



void PlotQASDD(Char_t fileName[100]="File.QA.111333.2010.LHC09b.pass2.root",Bool_t kDoEps=kFALSE,Char_t eventspecie[25]="LowMultiplicity")
{
  const TString Rawsexpertname[]={"SDDModPattern",//0     1
				  "SDDModPatternNORM",//1  //1
				  "SDDphizL3NORM",    //2  //0
				  "SDDphizL4NORM",    //3  //0
				  "SDDL3_RelativeOccupancy", //4     4  5
				  "SDDL4_RelativeOccupancy",//5      4  5
				  "SDDRawDataCheck"//6              6
                                  };//tot 7
  const TString Rawsnonexpertname[]={"SDDphizL3",//0   2         
				     "SDDphizL4",//1    2
				     "SDDphizCalibL3",//2  //   3
				     "SDDphizCalibL4"//3     3
                                     };//tot 4

  const TString RecPointsexpertname[]={"SDDLocalCoordDistrib",//0    //11
				       "SDDModPatternL3RP",//1       //12
				       "SDDModPatternL4RP",//2       //12
				       "SDDModPatternL3RPNORM",//3   //13
				       "SDDModPatternL4RPNORM",//4   //13
				       "SDDModPatternRP",//5         //14
				       "SDDModPatternRPNORM",//6     //14
				       "SDDphidistrib_Layer3",//7    //15
				       "SDDphidistrib_Layer4",//8    //15
				       "SDDrdistrib_Layer3",//9      //16
				       "SDDrdistrib_Layer4",//10     //16
				       "SDDL3_RelativeOccupancy",//11 //19  20
				       "SDDL4_RelativeOccupancy",//12 //19  20
				       "SDDL3_Rec2Raw_2D",       //13 //21
				       "SDDL4_Rec2Raw_2D",       //14 //21
				       "SDDL3_Rec2Raw",          //15 //22  23 
				       "SDDL4_Rec2Raw",          //16 //22  23 
				       "SDDL3_dedx",             //17 //24  25
				       "SDDL4_dedx",             //18   24  25
				       "SDDRecPointCheck"//19           26
                                       };//tot 20

  const TString RecPointsnonexpertname[]={"SDDGlobalCoordDistribL3PHIZ",//0      7
					  "SDDGlobalCoordDistribL4PHIZ",//1      7
					  "SDDGlobalCoordDistribRZ",//2          8
					  "SDDGlobalCoordDistribYX",//3          8
					  "SDDLay3TotCh",           //4          9    17
					  "SDDLay4TotCh",           //5          9    17
					  "SDDdrifttime_Layer3",    //6          10   18
					  "SDDdrifttime_Layer4"//7               10   18
                                         };//tot 8

  const TString canvassavedname[]={"RawLayersNORM",                 //0
				   "Rawpatterns",                   //1
				   "RawLayers",                     //2
				   "CalibConfiguration",            //3
				   "RawLayersRelOcc",               //4
				   "RawLayersRelOccBoth",           //5
				   "CheckRawdata",                  //6
				   "RecPointsphiz",                 //7
				   "RecPointsGlobalDistributions",  //8
				   "RecPointsCharge",               //9
				   "RecPointsDriftTime",            //10
				   "RecPointsLocalDistribution",    //11
				   "RecPointsLayers",               //12
				   "RecPointsLayersNORM",           //13
				   "RecPointsLayersPatterns",       //14
				   "RecPointsPhiDistribution",      //15
				   "RecPointsRDistribution",        //16
				   "RecPointsChargeBothLayers",     //17
				   "RecPointsDriftTimeBothLayers",  //18
				   "RecPointsRelLayersOcc",         //19
				   "RecPointsLayersRelOccBoth",     //20
				   "RecPointsLayersRecRaw2D",       //21
				   "RecPointsLayersRecRaw",         //22
				   "RecPointsLayersRecRawBoth",     //23
				   "RecPointsLayersdedx",           //24
				   "RecPointsLayersdedxBoth",       //25
				   "CheckRecPoints"                 //26
				  };//tot 27

  Char_t trendfile[25]="FileQAtrend.root";

  TString psfilestart=Form("SDDQAPlot.ps[");
  TString psfile     =Form("SDDQAPlot.ps");
  TString psfileend  =Form("SDDQAPlot.ps]");


  gStyle->SetTitleFont(62,"");
  gStyle->SetTitleFontSize(0.025);
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
  
  //  TH1F *histo2save[9];
  //for(Int_t i=0;i<9;i++){histo2save[i]=NULL;}
  
  TF1 *f1 = new TF1("f1","1",-100000,100000); // this is the function defined to multiply histograms
  Char_t layer[10];
  float fmax=0.;
  float fmaxold=0.;
  float fmaxmargin;
  
  Int_t cannum=0;
  //Int_t histo2savenumber[9]={4,5,6,7,7,8,3,5,1};
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
    canvas1->Print(psfilestart.Data());
    directory=(TDirectory*)mergedfile.Get(filepath);
    if(!directory){printf("...Not Found\n ");}//faccio l'istogramma
    else{
      
      printf("...Found: The histograms of this EventSpecie will be displayed\n");
      printf("0\n");
      canvas1->Clear();
      canvas1->Divide(2,1);
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      for(Int_t iraws=2;iraws<4;iraws++){//non expert raws
	canvas1->cd(iraws-1);
	
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[iraws].Data());
	printf("histo name %s ",Rawsexpertname[iraws].Data());
	historaw=(TH2D*)mergedfile.Get(histoname);
	if(historaw){printf("...Found\n"); 	  historaw->Multiply(f1,fCNinv);  historaw->DrawCopy("colz");}
	else
	  {
	    sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[iraws].Data());
	    printf("...Not found\nSecond try for histo name %s ",Rawsexpertname[iraws].Data());
	    historaw=(TH2D*)mergedfile.Get(histoname);
	    if(historaw){printf("...Found\n");	  historaw->Multiply(f1,fCNinv);historaw->DrawCopy("colz");} 
	    else
	      {
		updatecanvas[iraws-2]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");
	      }
	  }
	historaw=NULL;
      }//end for	
      
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data());
	}
      //else{delete canvas1; canvas1=NULL;}
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("1\n");
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
	else {updatecanvas[inraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}
      //else{delete canvas2; canvas2=NULL;}
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      printf("2\n");
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
	  //	  historaw->SetTitleSize(0.02);

	  historaw->DrawCopy("colz");
	}
	else
	  {
	    sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsnonexpertname[iraws].Data());
	    printf("...Not Found.\n Second Try for histo name %s",Rawsnonexpertname[iraws].Data());
	    historaw=(TH2D*)mergedfile.Get(histoname);
	    if(historaw){
	      printf("...Found\n");
	      //	      historaw->SetTitleSize(0.02);
	      historaw->Multiply(f1,fCNinv);
	      historaw->DrawCopy("colz");
	    }
	    else
	      {
		updatecanvas[iraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");
	      }
	  } 
	historaw=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}
      //else{delete canvas3; canvas3=NULL;}
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      printf("3\n");
      cannum++;
      canvas1->Clear();
      canvas1->Divide(2,1);
      for(Int_t iraws=2;iraws<4;iraws++){//non expert raws
	canvas1->cd(iraws-1);
	
	sprintf(histoname,"%s/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsnonexpertname[iraws].Data());
	printf("histo name %s",Rawsnonexpertname[iraws].Data());
	historaw=(TH2D*)mergedfile.Get(histoname);
	if(historaw){
	  printf("...Found\n");
	  //	  historaw->SetTitleSize(0.02);
	  //	  historaw->Multiply(f1,fCNinv);
	  historaw->DrawCopy("colz");
	}
	else
	  {
	    updatecanvas[iraws]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");
	  } 
	historaw=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}

      //--------------------- new plots
      printf("================= 4\n");
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      cannum++;
      canvas1->Clear();
      canvas1->Divide(2,1);
      for(Int_t iraws=4;iraws<rawtot-1;iraws++){//non expert raws
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
	  canvas1->Print(psfile.Data()); 
	}
      //else{delete canvas3; canvas3=NULL;}
      
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("5\n");
      cannum++;
      for (Int_t i=4;i<rawtot-1;i++){
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
      for(Int_t irrpp=4;irrpp<rawtot-1;irrpp++){//non expert raws
	
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
	  canvas1->Print(psfile.Data()); 
	}
      //---------------------------- summary plot ------------------------//
      printf("6\n");
      canvas1->Clear();
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[rawtot-1].Data());
      printf("histo name %s",Rawsexpertname[rawtot-1].Data());
      cannum++;
      historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->SetStats(0);
	    historaw2->DrawCopy();
	    canvas1->Update();
	  }
	else{updatecanvas[0]=kFALSE;updatecanvas[1]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
	if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	  {
	    canvas1->Update();
	    if(kDoEps){
	      sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	      canvas1->SaveAs(namecanvas);
	    }
	    canvas1->Print(psfile.Data()); 
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
      printf("7\n");
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
	  canvas1->Print(psfile.Data()); 
	}
      //	else{delete canvas4; canvas4=NULL;}
      
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("8\n");
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
	  canvas1->Print(psfile.Data()); 
	}
      //else{delete canvas5; canvas5=NULL;}
      
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("9\n");
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
	  canvas1->Print(psfile.Data()); 
	}
      //else{delete canvas6; canvas6=NULL;}
      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);
      printf("10\n");
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
	  canvas1->Print(psfile.Data()); 
	}
      //else{delete canvas7; canvas7=NULL;}

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      printf("11\n");
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
	canvas1->Print(psfile.Data()); 
      }
      //else{delete canvas8; canvas8=NULL;}
       
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("12\n");
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
	  canvas1->Print(psfile.Data()); 
	}
      
      //else{delete canvas9; canvas9=NULL;}
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("13\n");
      cannum++;
      for(Int_t i=3;i<5;i++)
	{
	  canvas1->cd(i-2);

	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  histodraw=(TH2F*)mergedfile.Get(histoname);
	  printf("histo name %s",RecPointsexpertname[i].Data());

	  if(histodraw){
	    printf("...Found\n");histodraw->Multiply(f1,fCNinv);
	    histodraw->DrawCopy("colz");
	  }else{updatecanvas[i-3]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  histodraw=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}

      
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(1,2);

      printf("14\n");
      cannum++;
      for(Int_t i=5;i<7;i++)
	{
	  canvas1->cd(i-4);

	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2)
	    {
	      printf("...Found\n");
	      if (i==4){historaw2->Multiply(f1,fCNinv);}
	      historaw2->DrawCopy();
	    }
	  else{updatecanvas[i-5]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}


      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("15\n");
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
	  canvas1->Print(psfile.Data()); 
	}
 
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("16\n");
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
	  canvas1->Print(psfile.Data()); 
	}
	  

      //superimposed
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("17\n");
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
	  canvas1->Print(psfile.Data()); 
	}
	  
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("18\n");
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
	  canvas1->Print(psfile.Data()); 
	}

      //------------------------------------------- new plot

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("19\n");
      cannum++;
      for(Int_t i=11;i<13;i++)
	{
	  canvas1->cd(i-10);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-11]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}
	  
      //------------------------------------
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("20\n");
      cannum++;
      for (Int_t i=11;i<13;i++){
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
      for(Int_t irrpp=11;irrpp<13;irrpp++){//non expert raws
	  
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[irrpp].Data());
	printf("histo name %s",RecPointsexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-8);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("RecPoint Relative Occupancy");
	    if (irrpp==11) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=11) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
	    legend->AddEntry(historaw2,layer,"l");
	    legend->Draw();
	    canvas1->Update();
	  }
	else{updatecanvas[irrpp-11]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	historaw2=NULL;
      }//end for
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}


      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("21\n");
      cannum++;
      for(Int_t i=13;i<15;i++)
	{
	  canvas1->cd(i-12);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw=(TH2D*)mergedfile.Get(histoname);

	  if(historaw){printf("...Found\n");  historaw->DrawCopy("colz");}else{updatecanvas[i-13]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}
	  


      //-----------------
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("22\n");
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
	  canvas1->Print(psfile.Data()); 
	}
      //--------------------------------------------

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("23\n");
      cannum++;
      for (Int_t i=15;i<17;i++){
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==15){
	    fmaxold=fmax;}
	  if (i!=15){
	    if(fmaxold<fmax){
	      fmaxold=fmax;
	    }
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
	    historaw2->SetTitle("Rec2Raw Ratio");
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
	  canvas1->Print(psfile.Data()); 
	}

	  
      //--------------------------------------
	  
      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      canvas1->Divide(2,1);
      printf("24\n");
      cannum++;
      for(Int_t i=17;i<19;i++)
	{
	  canvas1->cd(i-16);
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	  printf("histo name %s",RecPointsexpertname[i].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);
		
	  if(historaw2){printf("...Found\n");historaw2->DrawCopy();}else{updatecanvas[i-17]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;
	}
      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}
	  

      //--------------------------------------------

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();
      legend->Clear();
      printf("25\n");
      cannum++;
      for (Int_t i=17;i<19;i++){
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[i].Data());
	historaw2=(TH1F*)mergedfile.Get(histoname);
	if(historaw2){
	  fmax=historaw2->GetMaximum();
	  if (i==17){
	    fmaxold=fmax;}
	  if (i!=17){
	    if(fmaxold<fmax){fmaxold=fmax;}
	  }
	}
	historaw2=NULL;
      }
      fmaxmargin=1.1*fmaxold;
      for(Int_t irrpp=17;irrpp<19;irrpp++){//non expert raws
    
	sprintf(histoname,"%s/Expert/%s_%s",filepath,AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[irrpp].Data());
	printf("histo name %s",RecPointsexpertname[irrpp].Data());
	sprintf(layer, "layer %d",irrpp-14);
	historaw2=(TH1F*)mergedfile.Get(histoname);
	gStyle->SetOptStat(0);
	  
	if(historaw2)
	  {
	    printf("...Found\n");
	    historaw2->GetYaxis()->SetRangeUser(0,fmaxmargin);
	    historaw2->SetTitle("RecPoint dEdx");
	    if (irrpp==17) {historaw2->SetStats(0);historaw2->SetLineColor(2);historaw2->DrawCopy();}
	    if (irrpp!=17) {historaw2->SetStats(0);historaw2->SetLineColor(4);historaw2->DrawCopy("same");}
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
	  canvas1->Print(psfile.Data()); 

      //------------------------summary plot ---------------//

      for(Int_t i=0;i<2;i++){updatecanvas[i]=kTRUE;}
      canvas1->Clear();

      printf("26\n");
      cannum++;


	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[rtot-1].Data());
	  printf("histo name %s",RecPointsexpertname[rtot-1].Data());
	  historaw2=(TH1F*)mergedfile.Get(histoname);

	  if(historaw2)
	    {
	      printf("...Found\n");
	
	      historaw2->DrawCopy();
	    }
	  else{updatecanvas[0]=kFALSE;updatecanvas[1]=kFALSE;printf("...Not Found: the histogram or this QA has been done before the histograms was added to QA\n");}
	  historaw2=NULL;

      if(updatecanvas[0]==kTRUE||updatecanvas[1]==kTRUE)
	{
	  canvas1->Update();
	  if(kDoEps){
	    sprintf(namecanvas,"%s.eps",canvassavedname[cannum].Data());
	    canvas1->SaveAs(namecanvas);
	  }
	  canvas1->Print(psfile.Data()); 
	}
    



	}
      canvas1->Print(psfileend.Data()); 

      TFile file2savefortrend(trendfile,"recreate");
      file2savefortrend.cd();
      
      for(Int_t isave=0;isave<nonrawtot;isave++)
	{
	  sprintf(histoname,"ITS/Raws/%s/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),Rawsnonexpertname[isave].Data());
	  historaw=(TH2D*)mergedfile.Get(histoname);
	  historaw->Write();
	  printf("Saved  %s\n",histoname);
	  historaw=NULL;
	}
      
      for(Int_t isave1=0;isave1<rawtot;isave1++)
	{
	  sprintf(histoname,"ITS/Raws/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),Rawsexpertname[isave1].Data());
	  if(isave1==2||isave1==3)
	    {
	      historaw=(TH2D*)mergedfile.Get(histoname);
	      historaw->Write();
	      historaw=NULL;
	    }
	  else
	    {
	      historaw2=(TH1F*)mergedfile.Get(histoname);
	      historaw2->Write();
	      historaw2=NULL;
	    }
	  printf("Saved %s\n",histoname);
	}
      
      for(Int_t isave2=0;isave2<nrtot;isave2++)
	{
	  sprintf(histoname,"ITS/RecPoints/%s/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsnonexpertname[isave2].Data());
	  if(isave2<4||isave2==17||isave2==18)
	    {
	      historaw=(TH2D*)mergedfile.Get(histoname);
	      historaw->Write();
	      historaw=NULL;
	    }
	  else
	    {
	      historaw2=(TH1F*)mergedfile.Get(histoname);
	      historaw2->Write();
	      historaw2=NULL;
	    }
	  printf("Saved %s\n",histoname);
	}
      for(Int_t isave3=0;isave3<rtot;isave3++)
	{
	  sprintf(histoname,"ITS/RecPoints/%s/Expert/%s_%s",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[isave3].Data());
	  //printf(histoname,"ITS/RecPoints/%s/Expert/%s_%s \n",AliRecoParam::GetEventSpecieName(ispecie),AliRecoParam::GetEventSpecieName(ispecie),RecPointsexpertname[isave3].Data());
	  if(isave3<3||isave3==11||isave3==12||isave3==17||isave3==18)
	    {
	      histodraw=(TH2F*)mergedfile.Get(histoname);
	      histodraw->Write();
	      histodraw=NULL;
	    }
	  else
	    {
	      historaw2=(TH1F*)mergedfile.Get(histoname);
	      historaw2->Write();
	      historaw2=NULL;
	    }
	  printf("Saved %s\n",histoname);
	}
      
      //for(Int_t iss=0;iss<9;iss++){printf("Saved %d\n",iss); histo2save[iss]->Write();}
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


