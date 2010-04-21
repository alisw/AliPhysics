/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONGain.h"
#include "AliMUONPedestal.h"
#include "AliMUONErrorCounter.h"
#include "AliMUON2DMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMpIntPair.h"

#include <TString.h>
#include <THashList.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <Riostream.h>

#include <sstream>

#define  NFITPARAMS 4

//-----------------------------------------------------------------------------
/// \class AliMUONGain
///
/// Implementation of calibration computing
///
/// add
/// 
///
/// \author Alberto Baldisseri, JL Charvet 
//-----------------------------------------------------------------------------




// functions

namespace {
  
  //______________________________________________________________________________
  Double_t funcLin (const Double_t *x, const Double_t *par)
  {
    /// Linear function
    return par[0] + par[1]*x[0];
  }
  
  //______________________________________________________________________________
  Double_t funcParabolic (const Double_t *x, const Double_t *par)
  {
    /// Parabolic function
    return par[0]*x[0]*x[0];
  }
  
  //______________________________________________________________________________
  Double_t funcCalib (const Double_t *x, const Double_t *par)  
  {
    /// Calibration function
    Double_t xLim= par[3];
    
    if(x[0] <= xLim) return par[0] + par[1]*x[0];
    
    Double_t yLim = par[0]+ par[1]*xLim;
    return yLim + par[1]*(x[0] - xLim) + par[2]*(x[0] - xLim)*(x[0] - xLim);
  }
  
}

/// \cond CLASSIMP
ClassImp(AliMUONGain)
/// \endcond

//______________________________________________________________________________
AliMUONGain::AliMUONGain()
: AliMUONPedestal(),
fInjCharge(0), 
fRootDataFileName(),
fnInit(1),
fnEntries(11),
fnbpf1(6),
fPrintLevel(0), 
fPlotLevel(0)
{
/// Default constructor

}
  
//______________________________________________________________________________
AliMUONGain::~AliMUONGain()
{
/// Destructor
}

//______________________________________________________________________________
TString AliMUONGain::WriteDummyHeader(void) 
{
///

  ostringstream stream;
  stream<<"//DUMMY FILE (to prevent Shuttle failure)"<< endl;
  stream<<"//================================================" << endl;
  stream<<"//       MUONTRKGAINda: Calibration run  " << endl;
  stream<<"//================================================" << endl;
  stream<<"//   * Run           : " << fRunNumber << endl; 
  stream<<"//   * Date          : " << fDate->AsString("l") <<endl;
  stream<<"//   * DAC           : " << fInjCharge << endl;
  stream<<"//-------------------------------------------------" << endl;

  return TString(stream.str().c_str());
}

//______________________________________________________________________________
void AliMUONGain::MakePedStoreForGain(TString shuttleFile)
{
/// Store Pedmean and sigma to pedestal-like ascii file

  ofstream fileout;
  TString tempstring;
  TString flatFile;
  TString outputFile;
  
  // Store pedestal map in root file
  TTree* tree = 0x0;

  // write dummy ascii file -> Shuttle
  if(fIndex<fnEntries)
    {  
      FILE *pfilew=0;
      pfilew = fopen (shuttleFile,"w");

      fprintf(pfilew,"//DUMMY FILE (to prevent Shuttle failure)\n");
      fprintf(pfilew,"//================================================\n");
      fprintf(pfilew,"//       MUONTRKGAINda: Calibration run  \n");
      fprintf(pfilew,"//=================================================\n");
      fprintf(pfilew,"//   * Run           : %d \n",fRunNumber); 
      fprintf(pfilew,"//   * Date          : %s \n",fDate->AsString("l"));
      fprintf(pfilew,"//   * DAC           : %d \n",fInjCharge);
      fprintf(pfilew,"//-------------------------------------------------\n");
      fclose(pfilew);
    }



  Finalize();
  MakeControlHistos();
  if(fPrintLevel>0)
    {
      // compute and store mean DAC values (like pedestals)
      flatFile = Form("%s.ped",fPrefixDA.Data());
      outputFile=flatFile;
      cout << "\n" << fPrefixDA.Data() << " : Flat file  generated  : " << flatFile.Data() << "\n";
      if (!outputFile.IsNull())  
      {
        ofstream out(outputFile.Data());
	MakeASCIIoutput(out);
	out.close();
      }      
    }

  TString mode("UPDATE");

  if (fIndex==1) {
    mode = "RECREATE";
  }
  TFile* histoFile = new TFile(fRootDataFileName.Data(), mode.Data(), "MUON Tracking Gains");

  // second argument should be the injected charge, taken from config crocus file
  // put also info about run number could be usefull
  AliMpIntPair* pair   = new AliMpIntPair(fRunNumber,fInjCharge );

  if (mode.CompareTo("UPDATE") == 0) {
    tree = (TTree*)histoFile->Get("t");
    tree->SetBranchAddress("run",&pair);
    tree->SetBranchAddress("ped",&fPedestalStore);

  } else {
    tree = new TTree("t","Pedestal tree");
    tree->Branch("run", "AliMpIntPair",&pair);
    tree->Branch("ped", "AliMUON2DMap",&fPedestalStore);
    tree->SetBranchAddress("run",&pair);
    tree->SetBranchAddress("ped",&fPedestalStore);

  }

  tree->Fill();
  tree->Write("t", TObject::kOverwrite); // overwrite the tree
  histoFile->Close();

  delete pair;
}

//______________________________________________________________________________
TString AliMUONGain::WriteGainHeader(Int_t nInit, Int_t nEntries, Int_t nbpf2, Int_t *numrun, Double_t *injCharge) 
{
/// Header of the calibration output file

  ostringstream stream;


  stream<<"//=======================================================" << endl;
  stream<<"//      Calibration file calculated by " << fPrefixDA.Data() <<endl;
  stream<<"//=======================================================" << endl;
  stream<<"//   * Run           : " << fRunNumber << endl; 
  stream<<"//   * Date          : " << fDate->AsString("l") <<endl;
  stream<<"//   * Statictics    : " << fNEvents << endl;
  if(fConfig)
  stream<<"//   * # of MANUS    : " << fNManuConfig << " read in the Det. config. " << endl;
  stream<<"//   * # of MANUS    : " << fNManu << " read in raw data " << endl;
  stream<<"//   * # of MANUS    : " << fNChannel/64 << " written in calibration file " << endl;
  stream<<"//   * # of channels : " << fNChannel << endl;
  stream<<"//-------------------------------------------------------" << endl;

  if(nInit==0)
    stream<<"//  "<< nEntries <<" DAC values  fit: "<< fnbpf1 << " pts (1st order) " << nbpf2 << " pts (2nd order)" << endl;
  if(nInit==1)
    stream<<"//  "<< nEntries <<" DAC values  fit: "<< fnbpf1 << " pts (1st order) " << nbpf2 << " pts (2nd order) DAC=0 excluded" << endl;
  stream<<"//   *  nInit = " << nInit << "  *  f1nbp = " << fnbpf1 << "  *  f2nbp = " <<  nbpf2 << endl; 

  stream<<"//   RUN     DAC   " << endl;
  stream<<"//-----------------" << endl;
  for (Int_t i = 0; i < nEntries; ++i) {
 	stream<<Form("//   %d    %5.0f \n",numrun[i],injCharge[i]);
  }
  stream<<"//=======================================" << endl;
  stream<<"// BP MANU CH.   a1      a2     thres. q " << endl;
  stream<<"//=======================================" << endl;

  return TString(stream.str().c_str());
}

//______________________________________________________________________________
TString AliMUONGain::WriteGainData(Int_t BP, Int_t Manu, Int_t ch, Double_t p1, Double_t p2, Int_t threshold, Int_t q) 
{
/// Write calibration parameters per channel

  ostringstream stream("");
  stream << Form("%4i %5i %2i %7.4f %10.3e %4i %2x\n",BP,Manu,ch,p1,p2,threshold,q);
  return TString(stream.str().c_str());

}

//_______________________________________________________________________________
void AliMUONGain::MakeGainStore(TString shuttleFile)
{
  /// Store gains in ASCII files
  ofstream fileout;
  ofstream filcouc;
  TString tempstring;	
  TString filename; 

  Double_t goodA1Min =  0.5;
  Double_t goodA1Max =  2.;
//   Double_t goodA2Min = -0.5E-03;
//   Double_t goodA2Max =  1.E-03;
  Double_t goodA2Min = -0.5E-01; // changed 28/10/2009 (JLC) <=> enlarged condition on a2
  Double_t goodA2Max =  1.E-01;
  // Table for uncalibrated  buspatches and manus
  THashList* uncalBuspatchManuTable = new THashList(1000,2);

  Int_t numrun[11];

  // open file MUONTRKda_gain_data.root
  // read again the pedestal for the calibration runs (11 runs)
  // need the injection charge from config file (to be done)
  // For each channel make a TGraphErrors (mean, sigma) vs injected charge
  // Fit with a polynomial fct
  // store the result in a flat file.

  if(fIndex==0)cout << " Root data file = " << fRootDataFileName.Data() << endl;  
  TFile*  histoFile = new TFile(fRootDataFileName.Data());

  AliMUON2DMap* map[11];
  AliMUONVCalibParam* ped[11];
  AliMpIntPair* run[11];

  //read back from root file
  TTree* tree = (TTree*)histoFile->Get("t");
  Int_t nEntries = tree->GetEntries();

  Int_t nbpf2 = nEntries - (fnInit + fnbpf1) + 1; // nb pts used for 2nd order fit

  // read back info
  for (Int_t i = 0; i < nEntries; ++i) {
    map[i] = 0x0;
    run[i] = 0x0;
    tree->SetBranchAddress("ped",&map[i]);
    tree->SetBranchAddress("run",&run[i]);
    tree->GetEvent(i);
    //        std::cout << map[i] << " " << run[i] << std::endl;
  }
  // RunNumber extracted from Root data fil
  if(fIndex==0)fRunNumber=(UInt_t)run[nEntries-1]->GetFirst();
  //     sscanf(getenv("DATE_RUN_NUMBER"),"%d",&gAliRunNumber);

  Double_t pedMean[11];
  Double_t pedSigma[11];
  for ( Int_t i=0 ; i<11 ; i++) {pedMean[i]=0.;pedSigma[i]=1.;};
  Double_t injCharge[11];
  Double_t injChargeErr[11];
  for ( Int_t i=0 ; i<11 ; i++) {injCharge[i]=0.;injChargeErr[i]=1.;};

  // some print
  cout<<"\n ********  MUONTRKGAINda for Gain computing (Last Run = " << fRunNumber << ") ********\n" << endl;
  cout<<" * Date          : " << fDate->AsString("l") << "\n" << endl;
  cout << " Entries = " << nEntries << " DAC values \n" << endl; 
  for (Int_t i = 0; i < nEntries; ++i) {
    cout<< " Run = " << run[i]->GetFirst() << "    DAC = " << run[i]->GetSecond() << endl;
    numrun[i] = run[i]->GetFirst();
    injCharge[i] = run[i]->GetSecond();
    injChargeErr[i] = 0.01*injCharge[i];
    if(injChargeErr[i] <= 1.) injChargeErr[i]=1.;
  }
  cout << "" << endl;

  //  print out in .log file

  (*fFilcout)<<"\n\n//=================================================" << endl;
  (*fFilcout)<<"//    MUONTRKGAINda: Gain Computing  Run = " << fRunNumber << endl;
  (*fFilcout)<<"//    RootDataFile  = "<< fRootDataFileName.Data() << endl;
  (*fFilcout)<<"//=================================================" << endl;
  (*fFilcout)<<"//* Date          : " << fDate->AsString("l") << "\n" << endl;



  // why 2 files ? (Ch. F.)  => second file contains detailed results
    FILE *pfilen = 0;
    if(fPrintLevel>1)
      {
        filename=Form("%s.param",fPrefixDA.Data());
        cout << " Second fit parameter file        = " << filename.Data() << "\n";
        pfilen = fopen (filename.Data(),"w");

        fprintf(pfilen,"//===================================================================\n");
        fprintf(pfilen,"//  BP MANU CH. par[0]     [1]     [2]     [3]      xlim          P(chi2) p1        P(chi2)2  p2\n");
        fprintf(pfilen,"//===================================================================\n");
        fprintf(pfilen,"//   * Run           : %d \n",fRunNumber); 
        fprintf(pfilen,"//===================================================================\n");
      }


  // file outputs for gain

  ofstream pfilew;
  pfilew.open(shuttleFile.Data());
  // Write Header Data of the .par file
  pfilew << WriteGainHeader(fnInit,nEntries,nbpf2,numrun,injCharge);

  // print mean and sigma values in file
  FILE *pfilep = 0;
  if(fPrintLevel>1)
    {
      filename=Form("%s.peak",fPrefixDA.Data());
      cout << " File containing Peak mean values = " << filename << "\n";
      pfilep = fopen (filename,"w");

      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"//   * Run           : %d \n",fRunNumber); 
      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"// BP  MANU  CH.    Ped.     <0>      <1>      <2>      <3>      <4>      <5>      <6>      <7>      <8>      <9>     <10> \n"); 
      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"//                 DAC= %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f  fC\n",injCharge[0],injCharge[1],injCharge[2],injCharge[3],injCharge[4],injCharge[5],injCharge[6],injCharge[7],injCharge[8],injCharge[9],injCharge[10]);
      fprintf(pfilep,"//                      %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f\n",injChargeErr[0],injChargeErr[1],injChargeErr[2],injChargeErr[3],injChargeErr[4],injChargeErr[5],injChargeErr[6],injChargeErr[7],injChargeErr[8],injChargeErr[9],injChargeErr[10]);
      fprintf(pfilep,"//==============================================================================================================================\n");
    }

  Double_t chi2    = 0.;
  Double_t chi2P2  = 0.;
  Double_t prChi2  = 0; 
  Double_t prChi2P2 =0;
  Double_t a0=0.,a1=1.,a2=0.;
  Int_t busPatchId ;
  Int_t manuId     ;
  Int_t channelId ;
  Int_t threshold = 0;
  Int_t q = 0;
  Int_t p1 =0;
  Int_t p2 =0;
  Double_t gain=10.; // max value (= bad value)
  Double_t capa=0.2; // internal capacitor (pF)

  //  plot out 

  TFile* gainFile = 0x0;
  TTree* tg = 0x0;
  if(fPlotLevel>0)
    {
      fHistoFileName=Form("%s.root",fPrefixDA.Data());
      gainFile = new TFile(fHistoFileName.Data(),"RECREATE","MUON Tracking gains");
      tg = new TTree("tg","TTree avec class Manu_DiMu");

      tg->Branch("bp",&busPatchId, "busPatchId/I");
      tg->Branch("manu",&manuId, "manuId/I");
      tg->Branch("channel",&channelId, "channelId/I");

      tg->Branch("a0",&a0, "a0/D");
      tg->Branch("a1",&a1, "a1/D");
      tg->Branch("a2",&a2, "a2/D");
      tg->Branch("Pchi2",&prChi2, "prChi2/D");
      tg->Branch("Pchi2_2",&prChi2P2, "prChi2P2/D");
      tg->Branch("Threshold",&threshold, "threshold/I");
      tg->Branch("q",&q, "q/I");
      tg->Branch("p1",&p1, "p1/I");
      tg->Branch("p2",&p2, "p2/I");
      tg->Branch("gain",&gain, "gain/D");
    }

  char graphName[256];

  // iterates over the first pedestal run
  TIter next(map[0]->CreateIterator());
  AliMUONVCalibParam* p;

  Int_t    nmanu         = 0;
  Int_t    nGoodChannel   = 0;
  Int_t    nBadChannel   = 0;
  Int_t    noFitChannel   = 0;
  Int_t    nplot=0;
  Double_t sumProbChi2   = 0.;
  Double_t sumA1         = 0.;
  Double_t sumProbChi2P2 = 0.;
  Double_t sumA2         = 0.;

  Double_t x[11], xErr[11], y[11], yErr[11];
  Double_t xp[11], xpErr[11], yp[11], ypErr[11];

  Int_t uncalcountertotal=0 ;
  Int_t unparabolicfit=0;

  while ( ( p = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
    {
      ped[0]  = p;
      unparabolicfit=0;

      busPatchId = p->ID0();
      manuId     = p->ID1();

      // read back pedestal from the other runs for the given (bupatch, manu)
      for (Int_t i = 1; i < nEntries; ++i) {
	ped[i] = static_cast<AliMUONVCalibParam*>(map[i]->FindObject(busPatchId, manuId));
      }

      // compute for each channel the gain parameters
      for ( channelId = 0; channelId < ped[0]->Size() ; ++channelId ) 
	{

	  Int_t n = 0;
	  for (Int_t i = 0; i < nEntries; ++i) {

	    if (!ped[i]) continue; //shouldn't happen.
	    pedMean[i]      = ped[i]->ValueAsDouble(channelId, 0);
	    pedSigma[i]     = ped[i]->ValueAsDouble(channelId, 1);

	    if (pedMean[i] < 0) continue; // not connected

	    if (pedSigma[i] <= 0) pedSigma[i] = 1.; // should not happen.
	    n++;
	  }

	  // print_peak_mean_values
	  if(fPrintLevel>1)
	    {

	      fprintf(pfilep,"%4i%5i%5i%10.3f",busPatchId,manuId,channelId,0.);
	      fprintf(pfilep,"%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f \n",pedMean[0],pedMean[1],pedMean[2],pedMean[3],pedMean[4],pedMean[5],pedMean[6],pedMean[7],pedMean[8],pedMean[9],pedMean[10]);
	      fprintf(pfilep,"                   sig= %9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f \n",pedSigma[0],pedSigma[1],pedSigma[2],pedSigma[3],pedSigma[4],pedSigma[5],pedSigma[6],pedSigma[7],pedSigma[8],pedSigma[9],pedSigma[10]);
	    }

	  // makegain 


	  // Fit Method:  Linear fit over gAlinbpf1 points + parabolic fit  over nbpf2  points) 
	  // nInit=1 : 1st pt DAC=0 excluded

	  // 1. - linear fit over gAlinbpf1 points

	  Double_t par[4] = {0.,0.5,0.,ADCMax()};
	  Int_t nbs   = nEntries - fnInit;
	  if(nbs < fnbpf1)fnbpf1=nbs;

	  Int_t fitproceed=1;
	  Int_t nbpf2Dynamic=nbpf2;
	  Int_t ADClimit=4090; // when RMS < 0.5 (in other cases mean values forced to 4095, see DA_PED)
	  for (Int_t j = 0; j < nbs; ++j)
	    {
	      Int_t k = j + fnInit;
	      x[j]    = pedMean[k];
	      if(x[j]<=0.){fitproceed=0; break;}
	      //	      if(x[j]>= ADCMax())
	      if(x[j]>= ADClimit)
		{
		  if(j < nbs-1){fitproceed=0; break;}
		  else { nbpf2Dynamic=nbpf2-1; break;}
		}
	      xErr[j] = pedSigma[k];
	      y[j]    = injCharge[k];
	      yErr[j] = injChargeErr[k];

	    }

	  TGraphErrors *graphErr;
	  if(!fitproceed) { p1=0; p2=0; noFitChannel++;}

	  if(fitproceed)
	    {
		      
	      TF1 *f1 = new TF1("f1",funcLin,0.,ADCMax(),2);
	      graphErr = new TGraphErrors(fnbpf1, x, y, xErr, yErr);

	      f1->SetParameters(0,0);

	      graphErr->Fit("f1","RQ");

	      chi2 = f1->GetChisquare();
	      f1->GetParameters(par);

	      delete graphErr;
	      graphErr=0;
	      delete f1;

	      prChi2 = TMath::Prob(chi2, fnbpf1 - 2);

	      Double_t xLim = pedMean[fnInit + fnbpf1 - 1];
	      Double_t yLim = par[0]+par[1] * xLim;

	      a0 = par[0];
	      a1 = par[1];

	      // 2. - Translation : new origin (xLim, yLim) + parabolic fit over nbf2 points
	      //checking:         if(busPatchId ==1841 && manuId==4)nbpf2Dynamic=2;
	      if(nbpf2Dynamic > 2)
		{
		  for (Int_t j = 0; j < nbpf2Dynamic; j++)
		    {
		      Int_t k  = j + (fnInit + fnbpf1) - 1;
		      xp[j]    = pedMean[k] - xLim;
		      xpErr[j] = pedSigma[k];

		      yp[j]    = injCharge[k] - yLim - par[1]*xp[j];
		      ypErr[j] = injChargeErr[k];
		    }

		  TF1 *f2 = new TF1("f2",funcParabolic,0.,ADCMax(),1);
		  graphErr = new TGraphErrors(nbpf2Dynamic, xp, yp, xpErr, ypErr);

		  graphErr->Fit(f2,"RQ");
		  chi2P2 = f2->GetChisquare();
		  f2->GetParameters(par);

		  delete graphErr;
		  graphErr=0;
		  delete f2;

		  prChi2P2 = TMath::Prob(chi2P2, nbpf2Dynamic-1);
		  a2 = par[0];
		}
	      else 
		{ 
		  unparabolicfit++;
		  (*fFilcout) << " Warning : BP = " << busPatchId << " Manu = " << manuId <<  " Channel = " << channelId <<": parabolic fit not possible (nbpf2=" <<  nbpf2Dynamic  << ") => a2=0 and linear fit OK" << std::endl;
		  if(unparabolicfit==1) std::cout << " Warning : BP = " << busPatchId << " Manu = " << manuId <<  ": no parabolic fit for some channels (nbpf2=" <<  nbpf2Dynamic  << "), linear fit is OK (see .log for details)" << std::endl;
		  a2=0. ; prChi2P2=0. ;
		}

	      par[0] = a0;
	      par[1] = a1;
	      par[2] = a2;
	      par[3] = xLim;

	      if(prChi2>0.999999)prChi2=0.999999 ; if(prChi2P2>0.999999)prChi2P2=0.9999999; // avoiding Pr(Chi2)=1 value
	      p1 = TMath::Nint(floor(prChi2*15))+1;    // round down value : floor(2.8)=2.
	      p2 = TMath::Nint(floor(prChi2P2*15))+1;
	      q  = p1*16 + p2;  // fit quality 

	      Double_t x0 = -par[0]/par[1]; // value of x corresponding to à 0 fC 
	      threshold = TMath::Nint(ceil(par[3]-x0)); // linear if x < threshold

	      if(fPrintLevel>1)
		{
		  fprintf(pfilen,"%4i %4i %2i",busPatchId,manuId,channelId);
		  fprintf(pfilen," %6.2f %6.4f %10.3e %4.2f %4i          %8.6f %8.6f   %x          %8.6f  %8.6f   %x\n",
			  par[0], par[1], par[2], par[3], threshold, prChi2, floor(prChi2*15), p1,  prChi2P2, floor(prChi2P2*15),p2);
		}


	      // tests
	      if(par[1]< goodA1Min ||  par[1]> goodA1Max) p1=0;
	      if(par[2]< goodA2Min ||  par[2]> goodA2Max) p2=0;

	    } // fitproceed

	  if(fitproceed && p1>0 && p2>0) 
	    {
	      nGoodChannel++;
	      sumProbChi2   += prChi2;
	      sumA1         += par[1];
	      sumProbChi2P2   += prChi2P2;
	      sumA2         += par[2];
	    }
	  else // bad calibration
	    {
	      nBadChannel++;
	      q=0;  
	      par[1]=0.5; a1=0.5; p1=0;
	      par[2]=0.;  a2=0.;  p2=0;
	      threshold=ADCMax();	

	      // bad calibration counter
	      char bpmanuname[256];
	      AliMUONErrorCounter* uncalcounter;

	      sprintf(bpmanuname,"bp%dmanu%d",busPatchId,manuId);
	      if (!(uncalcounter = (AliMUONErrorCounter*)uncalBuspatchManuTable->FindObject(bpmanuname)))
		{
		  // New buspatch_manu name
		  uncalcounter= new AliMUONErrorCounter (busPatchId,manuId);
		  uncalcounter->SetName(bpmanuname);
		  uncalBuspatchManuTable->Add(uncalcounter);
		}
	      else
		{
		  // Existing buspatch_manu name
		  uncalcounter->Increment();
		}
	      //			    uncalcounter->Print_uncal()
	      uncalcountertotal ++;
	    }
	  gain=1./(par[1]*capa); // mv/fC

	  if(fPlotLevel>0)
	    {if(fPlotLevel>1)
		{
		  //		      if(q==0  and  nplot < 100)
		  // 	  if(p1>1 && p2==0  and  nplot < 100)
		  //		  	    if(p1>10 && p2>10  and  nplot < 100)
		  //	if(p1>=1 and p1<=2  and  nplot < 100)
//		  if((p1==1 || p2==1) and  nplot < 100)
		  	    if(nbpf2Dynamic<nbpf2  and  nplot < 100)
		    {
		      nplot++;
		      // 	      cout << " nplot = " << nplot << endl;
		      TF1 *f2Calib = new TF1("f2Calib",funcCalib,0.,ADCMax(),NFITPARAMS);

		      graphErr = new TGraphErrors(nEntries,pedMean,injCharge,pedSigma,injChargeErr);

		      sprintf(graphName,"BusPatch_%d_Manu_%d_Ch_%d",busPatchId, manuId,channelId);

		      graphErr->SetTitle(graphName);
		      graphErr->SetMarkerColor(3);
		      graphErr->SetMarkerStyle(12);
		      graphErr->Write(graphName);

		      sprintf(graphName,"f2_BusPatch_%d_Manu_%d_Ch_%d",busPatchId, manuId,channelId);
		      f2Calib->SetTitle(graphName);
		      f2Calib->SetLineColor(4);
		      f2Calib->SetParameters(par);
		      f2Calib->Write(graphName);

		      delete graphErr;
		      graphErr=0;
		      delete f2Calib;
		    }
		}
	      tg->Fill();
	    }
	  pfilew << WriteGainData(busPatchId,manuId,channelId,par[1],par[2],threshold,q);
	}
      nmanu++;
      Int_t step=500;
      if(nmanu % step == 0)std::cout << " Nb manu = " << nmanu << std::endl;
    }

  //      print in logfile
  if (uncalBuspatchManuTable->GetSize())
    {
      uncalBuspatchManuTable->Sort();  // use compare
      TIterator* iter = uncalBuspatchManuTable->MakeIterator();
      AliMUONErrorCounter* uncalcounter;
      (*fFilcout) << "\n List of problematic BusPatch and Manu " << endl;
      (*fFilcout) << " ========================================" << endl;
      (*fFilcout) << "        BP       Manu        Nb Channel  " << endl ;
      (*fFilcout) << " ========================================" << endl;
      while((uncalcounter = (AliMUONErrorCounter*) iter->Next()))
	{
	  (*fFilcout) << "\t" << uncalcounter->BusPatch() << "\t " << uncalcounter->ManuId() << "\t\t"   << uncalcounter->Events() << endl;
	}
      (*fFilcout) << " ========================================" << endl;

      (*fFilcout) << " Number of bad calibrated Manu    = " << uncalBuspatchManuTable->GetSize() << endl ;
      (*fFilcout) << " Number of bad calibrated channel = " << uncalcountertotal << endl;
	
    }


  (*fFilcout) << "\n Nb of channels in raw data = " << nmanu*64 << " (" << nmanu << " Manu)" <<  endl;
  (*fFilcout) << " Nb of calibrated channel   = " << nGoodChannel << " (" << goodA1Min << "<a1<" << goodA1Max 
	      << " and " << goodA2Min << "<a2<" << goodA2Max << ") " << endl;
  (*fFilcout) << " Nb of uncalibrated channel = " << nBadChannel << " (" << noFitChannel << " unfitted)" << endl;

  cout << "\n Nb of channels in raw data = " << nmanu*64 << " (" << nmanu << " Manu)" <<  endl;
  cout << " Nb of calibrated channel   = " << nGoodChannel << " (" << goodA1Min << "<a1<" << goodA1Max 
       << " and " << goodA2Min << "<a2<" << goodA2Max << ") " << endl;
  cout << " Nb of uncalibrated channel = " << nBadChannel << " (" << noFitChannel << " unfitted)" << endl;

  Double_t meanA1         = sumA1/(nGoodChannel);
  Double_t meanProbChi2   = sumProbChi2/(nGoodChannel);
  Double_t meanA2         = sumA2/(nGoodChannel);
  Double_t meanProbChi2P2 = sumProbChi2P2/(nGoodChannel);

  Double_t capaManu = 0.2; // pF
  (*fFilcout) << "\n linear fit   : <a1> = " << meanA1 << "\t  <gain>  = " <<  1./(meanA1*capaManu) 
	      << " mV/fC (capa= " << capaManu << " pF)" << endl;
  (*fFilcout) <<   "        Prob(chi2)>  = " <<  meanProbChi2 << endl;
  (*fFilcout) << "\n parabolic fit: <a2> = " << meanA2  << endl;
  (*fFilcout) <<   "        Prob(chi2)>  = " <<  meanProbChi2P2 << "\n" << endl;

  cout << "\n  <gain>  = " <<  1./(meanA1*capaManu) 
       << " mV/fC (capa= " << capaManu << " pF)" 
       <<  "  Prob(chi2)>  = " <<  meanProbChi2 << endl;
  
  pfilew.close();

  if(fPlotLevel>0){tg->Write();histoFile->Close();}
  if(fPrintLevel>1){fclose(pfilep); fclose(pfilen);}
}
