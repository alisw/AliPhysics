#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <iostream>
  #include <fstream>
  
  #include <TStopwatch.h>
  #include <TStyle.h>
  #include <TFile.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TString.h>
  #include <TCanvas.h>
  #include <TChain.h>
  #include <TSystem.h>
  #include <TGrid.h>
  #include <TMath.h>

  #include "AliPHOSFEEMapRun2.h"
  
#endif

using namespace std;

Bool_t CreatePedestalTable(Int_t runnumber)
{
  TFile *rootfile = TFile::Open(Form("ped%d.root",runnumber),"READ");
  if (rootfile)
    printf("File ped%d.root is opened successfully\n",runnumber);
  else {
    printf("Cannot open ped%d.root\n",runnumber);
    return kFALSE;
  }
  TFile *outfile = new TFile(Form("PedestalTable_%d.root",runnumber),"RECREATE");

  const Int_t Nmod=5;
  const Int_t Nsru=4;
  const Int_t Ndtc=40;
  const Int_t Nx=64;
  const Int_t Nz=56;

  TH2F *hALTROCHMapHigh[Nmod];
  TH2F *hALTROCHMapLow [Nmod];

  TH2F *hPedHiMean[Nmod];
  TH2F *hPedLoMean[Nmod];

  AliPHOSFEEMapRun2 *calib = new AliPHOSFEEMapRun2();
  for(Int_t imod=0;imod<Nmod;imod++){
    hPedHiMean[imod] = (TH2F*)rootfile->Get(Form("hPedHiMeanm%d",imod));
    hPedLoMean[imod] = (TH2F*)rootfile->Get(Form("hPedLoMeanm%d",imod));

    hALTROCHMapHigh[imod] = new TH2F(Form("hALTROCHMapHighM%d",imod),"hALTROCHMapHigh",64,0,64,56,0,56);
    hALTROCHMapLow[imod]  = new TH2F(Form("hALTROCHMapLowM%d",imod),"hALTROCHMapLow"  ,64,0,64,56,0,56);

    if(!hPedHiMean[imod] || !hPedLoMean[imod]) continue;

    outfile->WriteTObject(hPedHiMean[imod]);
    outfile->WriteTObject(hPedLoMean[imod]);

    for(Int_t isru=0;isru<Nsru;isru++){
      if(imod==1 && (isru==0 || isru==1)) continue;//SRU M1-0 and M1-1 do not exist in Run2.

      for(Int_t idtc=0;idtc<Ndtc;idtc++){

        if(idtc==0 || idtc==20) continue;//for TRU.
        if((14 < idtc && idtc<20) || (34 < idtc && idtc<40)) continue;//empty dtc port

        fstream ftxt;
        Char_t name[255];
	Int_t iDTC=idtc;
	if (idtc==1) iDTC=15; // FEE 1 is connected to SRU port 15
        sprintf(name,"%d/Pedestal_%d_SRUM%d-%d_DTC%d.txt",runnumber,runnumber,imod,isru,iDTC);
	gSystem->Exec(Form("mkdir -p %d",runnumber));
        ftxt.open(name,ios::out);

        for(Int_t ix=0;ix<Nx;ix++){// high gain
          for(Int_t iz=0;iz<Nz;iz++){
            Int_t SRU = calib->CellToSRUID(ix);
            Int_t DTC = calib->CellToFEEID(iz);

            if(isru != SRU) continue;
            if(idtc != DTC) continue;

            Int_t ALTRO = calib->CellToALTRO(ix,iz);
            Int_t CSP = calib->CellToCSPID(imod,ix,iz);
            Int_t ALTROCH_Hi = calib->CSPToALTROChannel(CSP,"High");
            hALTROCHMapHigh[imod]->SetBinContent(ix+1,iz+1,ALTROCH_Hi);

            Int_t PED = round(hPedHiMean[imod]->GetBinContent(ix+1,iz+1));


						ftxt <<"ALTRO: " << ALTRO <<" CHANNEL: "<< ALTROCH_Hi << " PED: "<< PED << endl;

          }// end of cell z loop
        }// end of cell x loop

        for(Int_t ix=0;ix<Nx;ix++){// low gain
          for(Int_t iz=0;iz<Nz;iz++){
            Int_t SRU = calib->CellToSRUID(ix);
            Int_t DTC = calib->CellToFEEID(iz);

            if(isru != SRU) continue;
            if(idtc != DTC) continue;

            Int_t ALTRO = calib->CellToALTRO(ix,iz);
            Int_t CSP = calib->CellToCSPID(imod,ix,iz);
            Int_t ALTROCH_Lo = calib->CSPToALTROChannel(CSP,"Low");
            hALTROCHMapLow[imod] ->SetBinContent(ix+1,iz+1,ALTROCH_Lo);

            Int_t PED = round(hPedLoMean[imod]->GetBinContent(ix+1,iz+1));

						ftxt <<"ALTRO: " << ALTRO <<" CHANNEL: "<< ALTROCH_Lo << " PED: "<< PED << endl;

          }// end of cell z loop
        }// end of cell x loop



        ftxt.close();

      }//end of dtc loop

    }//end of sru loop

    hALTROCHMapHigh[imod]->SetMinimum(-0.5);
    hALTROCHMapHigh[imod]->SetMaximum(15.5);
    hALTROCHMapHigh[imod]->SetContour(16);

    hALTROCHMapLow[imod] ->SetMinimum(-0.5);
    hALTROCHMapLow[imod] ->SetMaximum(15.5);
    hALTROCHMapLow[imod]->SetContour(16);


    outfile->WriteTObject(hALTROCHMapHigh[imod]);
    outfile->WriteTObject(hALTROCHMapLow[imod]);

  }//end of module loop


  outfile->Close();
  return kTRUE;

}

