#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH2F.h>
#include <fstream.h>
#include <TMath.h>
#include <TProfile.h>

#endif

/* Macro to produce LHC file for the luminous region.
   Author: caffarri@pd.infn.it
*/ 


void VtxInfoLHC(TString vtxtype="TRKnc",
		TString fname="Vertex.Performance.root",
		TString ntname="fNtupleBeamSpot"){
  
  // const Int_t nbinsmult=11;
  // Float_t limmult[nbinsmult+1]={0.,1.5, 2.5, 5.5,10.5,14.5,24.5,30.5,44.5,60.5,74.5, 99999.}; 
  
  Float_t totev=0,totevtriggered=0;

  Float_t resolVtx=510.;
  Float_t errResolVtx = 10; 
  Float_t p2= 1.3;
  Float_t errp2 = 0.1;
  Float_t meanMult=38;
  Float_t errTotResol =(errResolVtx*errResolVtx)+(0.3*(errp2*errp2)/TMath::Power(meanMult,p2+1)); 

  Float_t timeBinning = 120.; 
  

  Char_t hisnam[20];
  Double_t xRange = 1.;
  Int_t xBins = 125.;
  Double_t zRange = 40.;
  Int_t zBins = 40.;

  /*
    Int_t binsHRes=250*0.5;
    Float_t rangeHRes=10000.*0.5;
    Float_t rangeGrRms = 5000; // micron
    Float_t rangeGrPull = 15; 
    Float_t rangeHPull=10.; 
    Float_t rangeHResS = 1.;
    Float_t rangeHResS = 10000*0.5.;
  */
  
  TFile *f=new TFile(fname.Data());
  TList *cOutput = (TList*)f->Get("cOutputVtxESD");
  TNtuple *nt=(TNtuple*)cOutput->FindObject(ntname.Data());
  // TNtuple *nt=(TNtuple*)f->FindObject(ntname.Data());
  Int_t nnnev=nt->GetEntries();
  printf("Events = %d\n",nnnev);
  
  Float_t xVtx,xVtxSc,xerrVtx;
  Float_t yVtx,yVtxSc,yerrVtx;
  Float_t zVtx,zVtxSc,zerrVtx;
  Float_t ntrklets,ncontrVtx,triggered;
  Float_t run, tstamp;
  
  TString sxx="x"; sxx.Append(vtxtype.Data());
  nt->SetBranchAddress(sxx.Data(),&xVtx);
  TString syy="y"; syy.Append(vtxtype.Data());
  nt->SetBranchAddress(syy.Data(),&yVtx);
  TString szz="z"; szz.Append(vtxtype.Data());
  nt->SetBranchAddress(szz.Data(),&zVtx);
  
  TString xerr="xerr"; xerr.Append(vtxtype.Data());
  nt->SetBranchAddress(xerr.Data(),&xerrVtx);
  TString yerr="yerr"; yerr.Append(vtxtype.Data());
  nt->SetBranchAddress(yerr.Data(),&yerrVtx);
  TString zerr="zerr"; zerr.Append(vtxtype.Data());
  nt->SetBranchAddress(zerr.Data(),&zerrVtx);
  
  
  TCanvas *stampCanvas = new TCanvas ("stampCanvas", "stampCanvas");
  stampCanvas->Divide(3);
  TCanvas *stampCanvasSig = new TCanvas ("stampCanvasSig", "stampCanvasSig");
  stampCanvasSig->Divide(3);
  TCanvas *diamCanvas = new TCanvas ("diamCanvas", "diamCanvas");
  diamCanvas->Divide(2);

  TString trkstitle;
  nt->SetBranchAddress("ntrklets",&ntrklets);
  trkstitle="Tracklets Multiplicity";
  
  TString ntrks="ntrks"; ntrks.Append(vtxtype.Data());
  nt->SetBranchAddress(ntrks.Data(),&ncontrVtx);
  // nt->SetBranchAddress(ntrks.Data(),&ntrklets);
  nt->SetBranchAddress("triggered",&triggered);
  nt->SetBranchAddress("run", &run);
  nt->SetBranchAddress("tstamp",&tstamp);
 
  Float_t minTstamp=10E+13;
  Float_t maxTstamp=-1;
  Int_t subevents;
  
  Float_t minRun = 10E+7;
  Float_t maxRun = 1;
  
  for (Int_t ientries=0; ientries<nt->GetEntriesFast(); ientries++){
    nt->GetEntry(ientries);
    
    if (tstamp<minTstamp) minTstamp=tstamp;
    if (tstamp>maxTstamp) maxTstamp=tstamp;
    
    if (run<minRun) minRun=run-10;
    if (run>maxRun) maxRun=run+10;    
  }

  Int_t nTimeBins = (Int_t)(maxTstamp-minTstamp) / timeBinning;
  
  printf ("number of time bins %d\n", nTimeBins);
  subevents=0;
  
  TH1F **hxVtx = new TH1F*[nTimeBins+1]; //"hxVtx", "hxVtx", 125,-1,1);
  TH1F **hyVtx = new TH1F*[nTimeBins+1]; //("hyVtx", "hyVtx", 125,-1,1);
  TH1F **hzVtx = new TH1F*[nTimeBins+1]; //("hzVtx", "hzVtx", 40,-40,40);
  
  TH1F **hxVtxMult = new TH1F*[nTimeBins+1]; //
  TH1F **hyVtxMult = new TH1F*[nTimeBins+1]; //

  for(Int_t kTimeBin=0;kTimeBin<nTimeBins+1;kTimeBin++){
    printf ("here we are, kTimeBin = %d\n", kTimeBin);
    
    sprintf(hisnam,"hxVtx%d",kTimeBin);
    hxVtx[kTimeBin]=new TH1F(hisnam,"",xBins,-xRange,xRange);
    sprintf(hisnam,"hyVtx%d",kTimeBin);
    hyVtx[kTimeBin]=new TH1F(hisnam,"",xBins,-xRange,xRange);
    sprintf(hisnam,"hzVtx%d",kTimeBin);
    hzVtx[kTimeBin]=new TH1F(hisnam,"",zBins,-zRange,zRange);
    sprintf(hisnam,"hxVtxMult%d",kTimeBin);
    hxVtxMult[kTimeBin]=new TH1F(hisnam,"",xBins,-xRange,xRange);
    sprintf(hisnam,"hyVtxMult%d",kTimeBin);
    hyVtxMult[kTimeBin]=new TH1F(hisnam,"",xBins,-xRange,xRange);

  }

  TH1F *htstampX = new TH1F("tstamp Vx","tstamp Vx", nTimeBins+1, minTstamp, maxTstamp);
  TH1F *htstampY = new TH1F("tstamp Vy","tstamp Vy", nTimeBins+1, minTstamp, maxTstamp);
  TH1F *htstampZ = new TH1F("tstamp Vz","tstamp Vz", nTimeBins+1, minTstamp, maxTstamp); 
  
  TH1F *htstampSigX = new TH1F("tstamp Sig x","tstamp Sig x", nTimeBins+1, minTstamp, maxTstamp);
  TH1F *htstampSigY = new TH1F("tstamp Sig y","tstamp Sig y", nTimeBins+1, minTstamp, maxTstamp);
  TH1F *htstampSigZ = new TH1F("tstamp Sig z","tstamp Sig z", nTimeBins+1, minTstamp, maxTstamp); 
  
  TH1F *htstampXdiam = new TH1F("tstamp Diam X","tstamp Diam X", nTimeBins+1, minTstamp, maxTstamp);
  TH1F *htstampYdiam = new TH1F("tstamp Diam Y","tstamp Diam Y", nTimeBins+1, minTstamp, maxTstamp);
  
  TCanvas *cVtx = new TCanvas("cVtx", "Vtx Canvas", 1000,700);
  cVtx->Divide(3);
  
  Int_t *fill = new Int_t [nTimeBins+1];
  
  for(Int_t iev=0;iev<nnnev;iev++){
    nt->GetEvent(iev);
    
    Int_t rightTimeBin =(Int_t)(tstamp-minTstamp)/ 120;
    printf ("right Time Bin %d \n", rightTimeBin);
    
    if ((run>=114783) && (run<=114798)) fill[rightTimeBin]=1005;
    if ((run>=114916) && (run<=114931)) fill[rightTimeBin]=1013;
    if ((run>=115186) && (run<=115193)) fill[rightTimeBin]=1019;
    if ((run>=115310) && (run<=115345)) fill[rightTimeBin]=1022;
    if ((run>=115393) && (run<=115414)) fill[rightTimeBin]=1023;
    if ((run>=115514) && (run<=115521)) fill[rightTimeBin]=1026;
    if ((run>=115880) && (run<=115892)) fill[rightTimeBin]=1031;
    if ((run>=116079) && (run<=116081)) fill[rightTimeBin]=1033;
    if ((run>=116102) && (run<=116134)) fill[rightTimeBin]=1034;
    if ((run>=116197) && (run<=116204)) fill[rightTimeBin]=1035;
    if ((run>=116287) && (run<=116288)) fill[rightTimeBin]=1038;
    if ((run>=116351) && (run<=116354)) fill[rightTimeBin]=1042;
    if ((run>=116401) && (run<=116403)) fill[rightTimeBin]=1044;
    if ((run>=116559) && (run<=116574)) fill[rightTimeBin]=1045;
    if ((run>=116609) && (run<=116611)) fill[rightTimeBin]=1046; 
    if ((run>=116640) && (run<=116645)) fill[rightTimeBin]=1047;
    if ((run>=116681) && (run<=116684)) fill[rightTimeBin]=1049;


    //      if((tstamp<stopTime)&&(tstamp>=iTstamp)){
    subevents++;
    // ncontrVtx = ntrklets;
    
    //printf ("ntrksTRKnc = %d ncontr = %d \n ", ntrklets, ncontrVtx);
    
    xVtxSc = 10000.*xVtx;
    yVtxSc = 10000.*yVtx;
    zVtxSc = 10000.*zVtx;
    
    if (ncontrVtx>0){
      hxVtx[rightTimeBin]->Fill(xVtx);
      hyVtx[rightTimeBin]->Fill(yVtx);
      hzVtx[rightTimeBin]->Fill(zVtx);
   
      if((ntrklets>30)&&(ntrklets<45)){
	hxVtxMult[rightTimeBin]->Fill(xVtx);
	hyVtxMult[rightTimeBin]->Fill(yVtx);
      }
    }

  }

  Double_t xMeanVtx, xMeanVtxErr, xSigmaVtx, xSigmaVtxErr;
  Double_t yMeanVtx, yMeanVtxErr, ySigmaVtx, ySigmaVtxErr;
  Double_t zMeanVtx, zMeanVtxErr, zSigmaVtx, zSigmaVtxErr;
  Double_t xVtxRms, yVtxRms, zVtxRms;
  Double_t xVtxMult, yVtxMult;  
  Double_t xDiam,  xDiamErr, yDiam,  yDiamErr, zDiam, zDiamErr;
  
  /*
    TString outFileName ="aliceDataVertex"; 
    TString outFileExt = ".txt";
    TString outFileFill ;
  */
  
  ofstream outFile("1045_lumireg_ALICE.txt");
  
  for (Int_t iTimeBin=0; iTimeBin<nTimeBins+1; iTimeBin++){

    if (fill[iTimeBin] != 1045) continue;
    
    Int_t corrTstamp = iTimeBin*120+minTstamp;

    if (iTimeBin == nTimeBins+1) {
      cVtx->cd(1);
      hxVtxMult[iTimeBin]->Draw();
      cVtx->cd(2);
      hyVtxMult[iTimeBin]->Draw();
    
    }

    TF1 *fitVtx, *fitVtxMult;
    if((hxVtx[iTimeBin]->GetEffectiveEntries()<500)){// && (hxVtx[iTimeBin]->GetEffectiveEntries()<500)){
      //outFile << "0." <<" " << "0."<<" ";
      xMeanVtx =0.;
      xMeanVtxErr = 0.;
      xSigmaVtx = 0.;
      xSigmaVtxErr = 0.;
      xDiam = 0.;
      xDiamErr =0.;
      xVtxRms =0.;
    }
    else{
      hxVtx[iTimeBin]->Fit("gaus", "M");
      fitVtx = hxVtx[iTimeBin]->GetFunction("gaus");
      xMeanVtx = fitVtx->GetParameter(1);
      xMeanVtxErr = fitVtx->GetParError(1);
      xSigmaVtx = fitVtx->GetParameter(2);
      xSigmaVtxErr = fitVtx->GetParError(2);
      xVtxRms = hxVtx[iTimeBin]->GetRMS();

      hxVtxMult[iTimeBin]->Fit("gaus", "M");
      fitVtxMult= hxVtxMult[iTimeBin]->GetFunction("gaus");
      xVtxMult = fitVtxMult->GetParameter(2)*10000;
      xDiam=TMath::Sqrt(xVtxMult*xVtxMult - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
      if(xDiam >= 500 ) xDiam=0.;
      xDiamErr = TMath::Sqrt(fitVtxMult->GetParError(2)*10000+ errTotResol);
      printf ("xDiam = %f errDiam = %f \n", xDiam, xDiamErr);
      
      }

    
    if((hyVtx[iTimeBin] ->GetEffectiveEntries()<500)){// && (hyVtxMult[iTimeBin]->GetEffectiveEntries()<500)){
      //outFile << "0." <<" " << "0."<<" ";
      yMeanVtx=0.;
      yMeanVtxErr = 0.;
      ySigmaVtx = 0.;
      ySigmaVtxErr = 0.;
      yDiam = 0;
      yDiamErr =0;
      yVtxRms=0.;
    }
    else{
      hyVtx[iTimeBin]->Fit("gaus", "M");
      fitVtx = hyVtx[iTimeBin]->GetFunction("gaus");
      yMeanVtx = fitVtx->GetParameter(1);
      yMeanVtxErr = fitVtx->GetParError(1);
      ySigmaVtx = fitVtx->GetParameter(2);
      ySigmaVtxErr = fitVtx->GetParError(2);
      yVtxRms = hyVtx[iTimeBin]->GetRMS();
      
      hyVtxMult[iTimeBin]->Fit("gaus", "M");
      fitVtxMult= hyVtxMult[iTimeBin]->GetFunction("gaus");
      yVtxMult = fitVtxMult->GetParameter(2)*10000;
      yDiam=TMath::Sqrt(yVtxMult*yVtxMult - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
      if(yDiam >= 500 ) yDiam=0.;
      yDiamErr = TMath::Sqrt(fitVtxMult->GetParError(2)*10000 + errTotResol);
      printf("yDiam = %f  errDiam= %f \n", yDiam, yDiamErr);
    }
    
    if((hzVtx[iTimeBin]->GetEffectiveEntries()<500)){
      //  outFile << "0." <<" " << "0."<<" ";
      zMeanVtx=0.;
      zMeanVtxErr = 0.;
      zSigmaVtx = 0.;
      zSigmaVtxErr = 0.;
      zDiam = 0.;
      zDiamErr =0.;  
      zVtxRms =0.;
    }
    else{
      hzVtx[iTimeBin]->Fit("gaus", "M");
      fitVtx = hzVtx[iTimeBin]->GetFunction("gaus");
      zMeanVtx = fitVtx->GetParameter(1);
      zMeanVtxErr = fitVtx->GetParError(1);
      zSigmaVtx = fitVtx->GetParameter(2);
      zSigmaVtxErr = fitVtx->GetParError(2);
      zVtxRms = hzVtx[iTimeBin]->GetRMS();

      zDiam = zSigmaVtx;
      zDiamErr = zSigmaVtxErr;
      
     
    }
      
    //outfile<<

    Int_t iBin = htstampX->FindBin(corrTstamp);
    htstampX->SetBinContent(iBin, xMeanVtx);
    htstampX->SetBinError(iBin, xMeanVtxErr); 

    iBin = htstampY->FindBin(corrTstamp);
    htstampY->SetBinContent(iBin, yMeanVtx);
    htstampY->SetBinError(iBin, yMeanVtxErr); 

    iBin = htstampZ->FindBin(corrTstamp);
    htstampZ->SetBinContent(iBin, zMeanVtx);
    htstampZ->SetBinError(iBin, zMeanVtxErr);

    iBin = htstampSigX->FindBin(corrTstamp);
    htstampSigX->SetBinContent(iBin, xSigmaVtx);
    htstampSigX->SetBinError(iBin, xSigmaVtxErr);

    iBin = htstampSigY->FindBin(corrTstamp);
    htstampSigY->SetBinContent(iBin, ySigmaVtx);
    htstampSigY->SetBinError(iBin, ySigmaVtxErr);

    iBin = htstampSigZ->FindBin(corrTstamp);
    htstampSigZ->SetBinContent(iBin, zSigmaVtx);
    htstampSigZ->SetBinError(iBin, zSigmaVtxErr);

    iBin = htstampXdiam->FindBin(corrTstamp);
    htstampXdiam->SetBinContent(iBin, xDiam);
    htstampXdiam->SetBinError(iBin, xDiamErr);

    iBin = htstampYdiam->FindBin(corrTstamp);
    htstampYdiam->SetBinContent(iBin,yDiam);
    htstampYdiam->SetBinError(iBin, yDiamErr);
  
    outFile <<"2"<<" "<<fill[iTimeBin]<<" "<<corrTstamp<<" ";
    outFile <<xMeanVtx*10<<" "<<xMeanVtxErr*10<<" ";
    outFile <<yMeanVtx*10<<" "<<yMeanVtxErr*10<<" ";
    outFile <<zMeanVtx*10<<" "<<zMeanVtxErr*10<<" ";
    outFile <<xSigmaVtx*10<<" "<<xSigmaVtxErr*10<<" ";
    outFile <<ySigmaVtx*10<<" "<<ySigmaVtxErr*10<<" ";
    outFile <<zSigmaVtx*10<<" "<<zSigmaVtxErr*10<<" ";
    outFile <<xDiam*0.001<<" "<<xDiamErr*0.001<<" ";
    outFile <<yDiam*0.001<<" "<<yDiamErr*0.001<<" ";
    outFile <<zDiam*10<<" "<<zDiamErr*10<<" ";
    outFile <<xVtxRms*10 <<" "<<yVtxRms*10<<" "<<zVtxRms*10<<endl;
    



  }
  printf ("Err Resol = %f \n", errTotResol);
  stampCanvas->cd(1);
  htstampX->Draw();
  stampCanvas->cd(2);
  htstampY->Draw();
  stampCanvas->cd(3);
  htstampZ->Draw();
  
  stampCanvasSig->cd(1);
  htstampSigX->Draw();
  stampCanvasSig->cd(2);
  htstampSigY->Draw();
  stampCanvasSig->cd(3);
  htstampSigZ->Draw(); 
  
  diamCanvas->cd(1);
  htstampXdiam->Draw();
  diamCanvas->cd(2);
  htstampYdiam->Draw();
}
