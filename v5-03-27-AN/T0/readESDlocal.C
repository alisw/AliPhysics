#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TMap.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TGLabel.h"
#include "TGrid.h"
#include "TFitResult.h"

#include "TControlBar.h"
#include <Riostream.h>
#include "AliCDBManager.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliT0LookUpValue.h"
#include "AliT0LookUpKey.h"
#include "AliT0Parameters.h"
#include "AliT0RawReader.h"
#include "AliT0CalibLaserData.h"
#include "AliESDTZERO.h"
#include "AliESDVertex.h"
#include "AliESDVZERO.h"
#include "AliESDEvent.h"

void readESDlocal(Int_t run) 
{
  //read T0 ESD from reconstructed T0 data 
  TH1F *hMean= new TH1F("hMean"," T0 ", 1000,2000,3000);
  // TH1F *hMeanTOF= new TH1F("hMeanTOF"," T0 ", 1000,-10,10);
  // TH1F *hOrA= new TH1F("hOrA"," OrA T0 ",1000,-10,10);
  //  TH1F *hOrC= new TH1F("hOrC"," OrC T0 ",1000,-10,10);
  TH1F *hVertex = new TH1F("hVertex","Z position of vertex", 180,-100,100);
  TH1F *hVertexcalc = new TH1F("hVertexcalc","Z position of vertex", 180,-100,100);
  TH1F *hVertexT0 = new TH1F("hVertexT0","Z position of vertex", 1000,-10,10);
 TH1F *hVertexT0calc = new TH1F("hVertexT0calc","Z position of vertex", 1000,-10,10);
  TH1F *hVertexComp = new TH1F("hVertexComparison ","Z position of vertex", 180,-30,30);
  TH1F *hVertexT0only = new TH1F("hVertexT0only","Z position of vertex without SPD vertex", 180,-30,30);
  TH1F *hVertexSPD = new TH1F("hVertexSPD","Z position of vertex SPD", 180,-30,30);

  TH2F *hVertexSPDT0 = new TH2F("hVertexSPDT0","Z position of vertex SPD-T0",   180,-30,30, 180,-30,30);
  //
   TH1F *hMeanTOF= new TH1F("hMeanTOF"," T0 ", 8000, -100,100);
    TH1F *hOrA= new TH1F("hOrA"," OrA T0 ",8000, -100,100);
    TH1F *hOrC= new TH1F("hOrC"," OrC T0 ",8000, -100,100);
    // TH1F *hMeanTOF= new TH1F("hMeanTOF"," T0 ", 2000, -8950,-8900);
    // TH1F *hOrA= new TH1F("hOrA"," OrA T0 ",2000, -8950,-8900);
    //   TH1F *hOrC= new TH1F("hOrC"," OrC T0 ",2000, -8950,-8900);


    TH1F *hMeanTOFcalc= new TH1F("hMeanTOFcalc"," T0 ", 2000, -8950,-8900);
    TH1F *hOrAcalc= new TH1F("hOrAcalc"," OrA T0 ",2000, -8950,-8900);
    TH1F *hOrCcalc= new TH1F("hOrCcalc"," OrC T0 ",2000, -8950,-8900);
    //  TH1F *hMeanTOFcalc= new TH1F("hMeanTOFcalc"," T0 ", 8000, -100,100);
    //   TH1F *hOrAcalc= new TH1F("hOrAcalc"," OrA T0 ",8000, -100,100);
    //   TH1F *hOrCcalc= new TH1F("hOrCcalc"," OrC T0 ",8000, -100,100);
   
  TH1F *hAmp[24];  TH1F * hTime[24]; TH1F * hTimeDiff[24]; TH1F * hTimecorr[24];
  //  TProfile * hAmpTime[24];
  TH2F * hAmpTime[24];
  
  for(Int_t ic=0; ic<24; ic++) 
    {
      hAmp[ic] = new TH1F(Form("hAmp%i",ic+1),"Amp",100, -10, 20 );
      //      hAmp[ic] = new TH1F(Form("hAmp%i",ic),"Amp",400, 300, 700 );
      hTime[ic] = new TH1F(Form("hTime%i",ic+1)," time",  5000,0,10000);
      hTimecorr[ic] = new TH1F(Form("hTimecorr%i",ic+1)," time",  5000,0,10000);
      hTimeDiff[ic] = new TH1F(Form("hTimeDiff%i",ic+1)," time", 300,-300,300);
      //   hAmpTime[ic] = new TProfile(Form("hAmpTime%i",ic)," time",200, 300,700, 49000, 53000);
      hAmpTime[ic] = new TH2F(Form("hAmpTime%i",ic+1)," time",100, 0,10, 100,3000,4000);
	// hAmpTime[ic] = new TH2F(Form("hAmpTime%i",ic)," time",250, 350,600, 400, 49400, 49800);
    }
  
  Float_t channelWidth=24.4;
  Float_t c = 0.0299792458; // cm/ps
  Float_t shift=0;
  Float_t fLatencyHPTDC = 9000;


 //run 124702
  // Float_t fLatencyL1 = 8.91358e+03;
  // Float_t fLatencyL1A = 8.91352e+03;
  // Float_t fLatencyL1C =8.91361e+03;

  //run 125097
  //   Float_t fLatencyL1 =   8.914520e+03 ;
  // Float_t fLatencyL1A =  8.914860e+03 ;
  // Float_t fLatencyL1C =  8.914180e+03;

//run 125295
  Float_t fLatencyL1 =  8.91406e+03  ;
  Float_t fLatencyL1A = 8.91401e+03;
  Float_t fLatencyL1C = 8.91412e+03 ;


  //  Float_t fLatencyL1 = 0;
  //  Float_t fLatencyL1A = 0;
  // Float_t fLatencyL1C =0;
     //125842
  //  Float_t fLatencyL1 = 8.91306e+03;
  //  Float_t fLatencyL1A =8.91338e+03;
  // Float_t fLatencyL1C =8.91274e+03;
 //
     //126407
     //  Float_t fLatencyL1 = 8.91345e+03;
     //  Float_t fLatencyL1A =8.91378e+03;
     // Float_t fLatencyL1C =8.91311e+03;
 //
  const Double32_t *amp, *time, *mean;
  Double32_t time1[24];
  TString filenam[100];
 
 

  //   Float_t equal[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  //124358
  //Float_t equal[24]={22, 0 ,0 ,1 ,3 ,1 ,0 ,0 ,-4 ,1 ,1 ,-1 ,0 ,-3 ,0 ,-6 ,0 ,-3 ,2 ,0 ,-1 ,-1 ,-3 ,0 };

     //124702
  // Float_t equal[24]={20, 0 ,15, 50, 37, 36, 23, -20, 77, 54, 44, -12, 
  //		     0, 19, -25, -13, 35, 1, 23, 46, 20, 77, 50, 122};


  //125097
   Float_t equal[24]={8,0 ,15 ,50 ,36 ,34 ,24 ,-21 ,81 ,53 ,42 ,-12 ,0 ,19 ,-20 ,-7 ,37 ,5 ,25 ,49 ,19 ,78 ,52 ,124};

  //  125295

  //Float_t equal[24]={16, 0,  2, 2, 30, 3, -5, -8, 5, -4, 16, -5, 
   //		 0, 10, 10, 12, 13, -24, 15, 26, -2, 17, 10, -30};

 

   //125842
  // Float_t equal[24]={16, 0, 1, 3,  32,   4, -2, -8, 5, -3, 15, -5, 
  //		      0, 10, 9, 12, 14, -23, 16, 26, 0, 19, 9, -31 };
 
  //126407
  //  Float_t equal[24]={16, 0,   1, 4, 33, 4, 0, -7, 6, -3, 15, -4,
    //                    0 , 10, 10, 13, 13, -23, 16, 28, 0, 19, 11, -29};

   //for ( Int_t indexfile=filestart; indexfile < filestop+1;indexfile++ ) 
    for ( Int_t indexfile=0; indexfile < 1; indexfile++ ) 
    {
      TString  fFileName=Form("%i/noeq/AliESDs.root", run);
      //  cout<<" filenam "<<fFileName<<endl;
      //  fFileName=filenam[indexfile];
       
      TFile *file = TFile::Open(fFileName);
      //  cout<<file<<endl;  
      
      //   TFile *file = new TFile("AliESDs.root");
      
      TTree *esdTree =  (TTree*)file->Get("esdTree");
      
      AliESDEvent *esdevent = new AliESDEvent;
      esdevent->ReadFromTree(esdTree);
      cout<<"esdevent "<<esdevent<<endl;
      
      TBranch *brSPD=esdTree->GetBranch("SPDVertex.");
      AliESDVertex *fverSPD = new AliESDVertex();
      if (brSPD) {
	brSPD->SetAddress(&fverSPD);
      }else{
	cerr<<"EXEC Branch SPDvertex  not found"<<endl;
	return;
      }
      
      TBranch *brRec=esdTree->GetBranch("AliESDTZERO.");
      AliESDTZERO *fesdT0 = new AliESDTZERO();
      if (brRec) {
	brRec->SetAddress(&fesdT0);
      }else{
	cerr<<"EXEC Branch T0  not found"<<endl;
	return;
      }
      cout<<" branches "<<brRec<<endl;
      // Event ------------------------- LOOP  
      for (Int_t ievent=0; ievent<brRec->GetEntries(); ievent++){
	
	Double32_t fOrA, besttimeA=9999999;
	Double32_t fOrC, besttimeC=9999999;
	for (Int_t ip=0; ip<24; ip++) time1[ip]=0;
	Float_t timecorr=0;
	Int_t ncont=0;
	Float_t shift=0;
	brRec->GetEntry(ievent);
	brSPD->GetEntry(ievent);
	esdTree->GetEvent(ievent);
	
	TString triggers = esdevent->GetFiredTriggerClasses();
	//	printf("Event %d: trigger classes: \"%s\"\n",
	//     esdevent->GetEventNumberInFile(),
	//     esdevent->GetFiredTriggerClasses().Data());
	if ( !triggers.Contains("CINT1B-ABCE-NOPF-ALL") ) {
	  //  Printf("Skip event with trigger class \"%s\"",triggers.Data());
	  continue;
	}
	Double_t spdver = fverSPD->GetZ();
	ncont = fverSPD->GetNContributors();
	// cout<<" spdver "<<spdver<<" ncont "<<ncont<<endl;
	if(ncont>2) {
	  hVertexSPD->Fill(spdver);
	  shift = spdver/(c*channelWidth );
	  //	  cout<<ievent<<" vertex shif "<<shift<<" vertex "<<spdver<<" IsFromVertexer3D  "<<fverSPD->IsFromVertexer3D()<<endl;
	}
	
	Int_t trig=fesdT0->GetT0Trig();
	Float_t   meanTOF = fesdT0->GetT0();
	mean = fesdT0->GetT0TOF();
	Double32_t vertex= fesdT0->GetT0zVertex();
	Double32_t orA=0.001* mean[1]  ;
	Double32_t orC= 0.001*mean[2] ;
	if (orA<99){
	  hOrA->Fill(orA);
	  //	  cout<<ievent<<" ORA "<< orA<<endl;
	}
	if (orC<99) {
	  hOrC->Fill(orC);
	  //	  cout<<ievent<<" ORC "<< orC<<endl;
	}

	if(orA<99 && orC<99) hVertexT0->Fill((orA-orC)/2.);
	
	if (vertex<99990) {
	  hMeanTOF->Fill(0.001*mean[0] );
	  hMean->Fill(meanTOF);
	  hVertex->Fill(vertex);
	  if(ncont>1) {
	    hVertexComp->Fill(vertex-spdver);
	    hVertexSPDT0->Fill(vertex,spdver);
	  }
	  if(ncont<0) hVertexT0only->Fill(vertex); 
	  
	}
	amp=fesdT0->GetT0amplitude();
	time=fesdT0->GetT0time();
	for (Int_t i=0; i<24; i++){ 
	  if( time[i]>0){
	    hAmp[i]->Fill(amp[i]);
	    hTime[i]->Fill(time[i]);
	    time1[i]=time[i] -equal[i];
	    
	    if(i<12){
	      if(time[1] >0) hTimeDiff[i]->Fill((time1[i] - time1[1]));
	      if(ncont>2) {	timecorr = time[i] - shift;
		hTimecorr[i]->Fill(timecorr);
		hAmpTime[i]->Fill(amp[i],timecorr);
	      }
	    }
	    
	    
	    if(i>11) {
	      if(time[12] >0)  hTimeDiff[i]->Fill((time1[i]- time1[12]));
	      if(ncont>2) {	timecorr = time[i] + shift;
		hTimecorr[i]->Fill(timecorr); 
		hAmpTime[i]->Fill(amp[i],timecorr);
	      }
	    } 
	    //	      cout<<ievent<<" "<<i<<" "<<<endl;
	   }
	} 
      	
	for (Int_t ipmt=0; ipmt<12; ipmt++){
	  
	    if( time1[ipmt]>1 && time1[ipmt] < besttimeC ){
	      besttimeC=time1[ipmt]; //timeC
	      
	  }
	}
	for ( Int_t ipmt=12; ipmt<23; ipmt++){
	    if( time1[ipmt]>1 && time1[ipmt]<besttimeA ) {
	      besttimeA=time1[ipmt]; //timeA
	  }
	}
	if(besttimeA<9999999 ) {
	  fOrA=0.001*besttimeA * channelWidth - fLatencyHPTDC +  fLatencyL1A +0.001*shift*channelWidth;
	  hOrAcalc->Fill(fOrA);
	  //	  cout<<ievent<<" calc OrA "<<fOrA<<endl;
	}
	if( besttimeC < 9999999 ) {
	  fOrC=0.001*besttimeC * channelWidth - fLatencyHPTDC +fLatencyL1C -0.001*shift*channelWidth;
	  hOrCcalc->Fill(fOrC);
	  //	  cout<<ievent<<" calc OrC "<<fOrC<<endl;
	}
	if(besttimeA <9999999 && besttimeC < 9999999 ) 
	  {
	    Double_t meanCalc= 0.001*channelWidth * Float_t( besttimeA+besttimeC)/2.- fLatencyHPTDC + fLatencyL1;
	    hMeanTOFcalc->Fill(meanCalc ); 
	    hVertexT0calc->Fill((fOrA-fOrC)/2.);
	    Float_t  timeDiff = ( besttimeA - besttimeC)* 0.001* channelWidth + fLatencyL1A - fLatencyL1C;
	 Float_t vertexcal =   - 1000*c*(timeDiff)/2. ; //+ (fdZonA - fdZonC)/2; 

	    hVertexcalc->Fill(vertexcal);


	    //	    cout<<ievent<<" calc meanTOF "<<meanCalc<<endl;
	  }
	
      }
      
    } 
    TFile* Hfile = new TFile(Form("FigESD_%i.noeq.root",run),"RECREATE","Histograms for T0 digits");
  printf("Writting histograms to root file \n");
  Hfile->cd();
  //Create a canvas, set the view range, show histograms
  //  
  Float_t differ;
  hMean->Write();
  hMeanTOF->Write();
  hMeanTOFcalc->Write();
  hVertex->Write();
  hVertexT0->Write();
  hVertexT0calc->Write();
  hVertexSPD->Write();
  hVertexComp->Write();
  hVertexT0only->Write();
  hOrC->Write();
  hOrA->Write();
  hOrCcalc->Write();
  hOrAcalc->Write();
  hVertexSPDT0->Write();
  hVertexcalc->Write();
  for (Int_t i=0; i<24; i++){ 
    hAmp[i]->Write();
    hTime[i]->Write();
     TFitResultPtr r = hTimeDiff[i]->Fit("gaus","SQ");
      Double_t par0;
       if ((Int_t)r == 0) par0 = r->Parameters()[1];
       else par0 = 99999;
    hTimeDiff[i]->Write();
    cout<<Int_t(par0)<<", ";
      //      cout<<Int_t(hTimeDiff[i]->GetMean())<<" fit  "<<par0<<endl;
      //  cout<<Int_t(hTimeDiff[i]->GetMean())<<" , ";
    hAmpTime[i]->Write();
   hTimecorr[i]->Write();
  }
  cout<<endl;
 for (Int_t i=1; i<24; i++)  cout<<i<<" ,";
  cout<<endl;
  for (Int_t i=1; i<24; i++)  cout<<Int_t(hTimeDiff[i]->GetMean())<<" ,";
  cout<<endl;

}
