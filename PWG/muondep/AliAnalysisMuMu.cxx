/**************************************************************************
 * Copyright(c) 1996-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
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

#include "AliAnalysisMuMu.h"
#include "AliAnalysisMuMuBase.h" //Just to have avaliable the MCInputPrefix() method

#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuConfig.h"
#include "AliAnalysisMuMuFnorm.h"
#include "AliAnalysisMuMuGraphUtil.h"
#include "AliAnalysisMuMuJpsiResult.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliAnalysisMuMuSpectraCapsulePbPb.h"
#include "AliAnalysisMuMuSpectraCapsulePbP.h"
#include "AliAnalysisMuMuSpectraCapsulePP.h"
#include "AliAnalysisTriggerScalers.h"
#include "AliCounterCollection.h"
#include "AliHistogramCollection.h"
#include "AliLog.h"
#include "AliMergeableCollection.h"
#include "Riostream.h"
#include "TArrayL64.h"
#include "TASImage.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TBox.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGrid.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "THashList.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TMap.h"
#include "TMath.h"
#include "TMethodCall.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TPaveText.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cassert>
#include <map>
#include <set>
#include <vector>
#include "TLatex.h"

using std::cout;
using std::endl;
using std::ifstream;

ClassImp(AliAnalysisMuMu)

//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu(const char* filename, AliAnalysisMuMuConfig& config) : TObject(),
fFilename(),
fDirectory(""),
fCounterCollection(0x0),
fBinning(0x0),
fMergeableCollection(0x0),
fRunNumbers(),
fCorrectionPerRun(0x0),
fAssociatedSimulation(0x0),
fAssociatedSimulation2(0x0),
fParticleName(""),
fConfig(new AliAnalysisMuMuConfig(config))
{
  GetFileNameAndDirectory(filename);

  GetCollections(fFilename,fDirectory,fMergeableCollection,fCounterCollection,fBinning,fRunNumbers);

  if ( IsSimulation() ) SetParticleNameFromFileName(fFilename);
}


//_____________________________________________________________________________
AliAnalysisMuMu::AliAnalysisMuMu(const char* filename, const char* associatedSimFileName, const char* associatedSimFileName2, const char* configurationFile) : TObject(),
fFilename(filename),
fDirectory(""),
fCounterCollection(0x0),
fBinning(0x0),
fMergeableCollection(0x0),
fRunNumbers(),
fCorrectionPerRun(0x0),
fAssociatedSimulation(0x0),
fAssociatedSimulation2(0x0),
fParticleName(""),
fConfig(0x0)
{
  /// ctor

  GetFileNameAndDirectory(filename);

  GetCollections(fFilename,fDirectory,fMergeableCollection,fCounterCollection,fBinning,fRunNumbers);

  if ( IsSimulation() ) SetParticleNameFromFileName(fFilename);

  if ( fCounterCollection ) {

    if ( strlen(associatedSimFileName) ) fAssociatedSimulation = new AliAnalysisMuMu(associatedSimFileName);
    if ( strlen(associatedSimFileName2) ) fAssociatedSimulation2 = new AliAnalysisMuMu(associatedSimFileName2);

    fConfig = new AliAnalysisMuMuConfig;

    fConfig->ReadFromFile(configurationFile);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMu::~AliAnalysisMuMu()
{
  /// dtor

  if ( fAssociatedSimulation ) fAssociatedSimulation->Update();
  if ( fAssociatedSimulation2 ) fAssociatedSimulation2->Update();

  Update();

  delete fCounterCollection;
  delete fBinning;
  delete fMergeableCollection;
  delete fCorrectionPerRun;
  delete fAssociatedSimulation;
  delete fAssociatedSimulation2;
  delete fConfig;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::BasicCounts(Bool_t detailTriggers,
                                  ULong64_t* totalNmb,
                                  ULong64_t* totalNmsl,
                                  ULong64_t* totalNmul)
{
/**
 * @brief Report of some basic trigger counts (for MB,MUL and MSL) both before and after physics selection.
 * @details Amount all the triggers available in our counter collection,we only consider the triggers that are defined in the configuration.
 * If detailTriggers is kTRUE, will show the detail, including Physics Selection fraction,
 * for each trigger found (as opposed to just showing info for MB,MSL and MUL triggers)
 *
 * @param detailTriggers
 * @param totalNmb
 * @param totalNmsl
 * @param totalNmul
 */

  if (!fMergeableCollection || !fCounterCollection) return;

  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);

  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);

  TObjArray* events = fCounterCollection->GetKeyWords("event").Tokenize(",");

  Bool_t doPS = (events->FindObject("PSALL") != 0x0);

  TObjString* srun;
  TObjString* strigger;

  ULong64_t localNmb(0);
  ULong64_t localNmsl(0);
  ULong64_t localNmul(0);

  if ( totalNmb)  *totalNmb = 0;
  if ( totalNmsl) *totalNmsl = 0;
  if ( totalNmul )*totalNmul = 0;

  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    std::cout << Form("RUN %09d ",srun->String().Atoi());

    TString details;
    ULong64_t nmb(0);
    ULong64_t nmsl(0);
    ULong64_t nmul(0);

    nextTrigger.Reset();

    Int_t nofPS(0);

    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {

      if ( !Config()->Has(Config()->MinbiasTriggerKey(),strigger->String().Data(),IsSimulation()) &&
           !Config()->Has(Config()->DimuonTriggerKey(),strigger->String().Data(),IsSimulation()) &&
           !Config()->Has(Config()->MuonTriggerKey(),strigger->String().Data(),IsSimulation()) ) continue;

      ULong64_t n = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d/centrality:all",
                                                    strigger->String().Data(),"ALL",srun->String().Atoi())));

      details += TString::Format("\n%50s %10lld",strigger->String().Data(),n);


      if ( doPS )
      {
        ULong64_t nps = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d/centrality:all",
                                                                    strigger->String().Data(),"PSALL",srun->String().Atoi())));

        details += TString::Format(" PS %5.1f %%",nps*100.0/n);

        if (nps) ++nofPS;
      }

      if ( Config()->Has(Config()->MinbiasTriggerKey(),strigger->String(),IsSimulation() ) )
      {
        nmb             += n;
        if ( totalNmb) (*totalNmb) += n;
        localNmb        += n;
      }
      else if ( Config()->Has(Config()->MuonTriggerKey(),strigger->String(),IsSimulation()))
      {
        nmsl             += n;
        if ( totalNmsl) (*totalNmsl) += n;
        localNmsl        += n;
      }
      else if ( Config()->Has(Config()->DimuonTriggerKey(),strigger->String(),IsSimulation()) )
      {
        nmul              += n;
        if ( totalNmul ) (*totalNmul) += n;
        localNmul         += n;
      }
    }

    std::cout << Form("MB %10lld MSL %10lld MUL %10lld %s",
                 nmb,nmsl,nmul,(nofPS == 0 ? "(NO PS AVAIL)": ""));



    if ( detailTriggers ) std::cout << details.Data();

    std::cout << std::endl;
  }

  if ( !totalNmul && !totalNmsl && !totalNmb ) std::cout << std::endl << Form("%13s MB %10lld MSL %10lld MUL %10lld ","TOTAL",localNmb,localNmsl,localNmul) << std::endl;

  delete runs;
  delete triggers;
  delete events;
}



//_____________________________________________________________________________
void AliAnalysisMuMu::CleanAllSpectra()
{
  /**
  * @brief Delete all the spectra we may have
  */
  OC()->RemoveByType("AliAnalysisMuMuSpectra");
  Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::CleanFNorm()
{
  /**
  * @brief Delete the FNorm related results we may have
  */

  OC()->Prune("/FNORM");
  Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::CleanMix()
{
  /**
   * @brief Delete all the spectra we may have
   */
  OC()->Prune("/MIX");
  Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DrawFitResults(const char* what,
                                     const char* spectraName,
                                     const char* subresults,
                                     Bool_t mix,
                                     Bool_t AccEffCorr
                                     )const
{
  /**
  * @brief Draw all fit results
  * @details Draw all results/subresults (i.e fit functions) spectras on a single canvas for every combination of Trigger/eventType/centrality/Pair cut.
  *
  * @param particle    particle name
  * @param binType     [description]
  * @param subresults  If one wants to draw one/several specific results. Tokenized by ','
  * @param AccEffCorr  Spectra type
  */

  if (!OC() || !CC()){
    AliError("No mergeable/counter collection. Consider Upgrade()");
    return ;
  }

  // Get configuration settings

  TObjArray* pairCutArray    = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
  TObjArray* centralityArray = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());

  TObjArray* triggerArray(0x0);
  if (!mix) triggerArray = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
  else                              triggerArray = Config()->GetListElements(Config()->MixTriggerKey(),IsSimulation());

  TObjArray* eventTypeArray(0x0);
  if (!mix) eventTypeArray = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
  else                              eventTypeArray = Config()->GetListElements(Config()->EventSelectionMixKey(),IsSimulation());

  TString refTrigger(Form("%s",Config()->First(Config()->RefMixTriggerKey(),IsSimulation()).Data()));
  TString refEvent(Form("%s",Config()->First(Config()->RefMixEventSelectionKey(),IsSimulation()).Data()));

  // Iterater for loops
  TIter nextTrigger(triggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextCentrality(centralityArray);

  // Strings
  TObjString* strigger;
  TObjString* seventType;
  TObjString* spairCut;
  TObjString* scentrality;

  AliAnalysisMuMuSpectra * spectra=0x0;

  //Loop on particle
  nextEventType.Reset();

  // Loop on each envenType (see MuMuConfig)
  while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
  {
    AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
    nextTrigger.Reset();

    // Loop on each trigger (see MuMuConfig)
    while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
    {
      AliDebug(1,Form("-TRIGGER %s",strigger->String().Data()));
      nextCentrality.Reset();

      // Loop on each centrality (not the ones in MuMuConfig but the ones set)
      while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
      {
        AliDebug(1,Form("--CENTRALITY %s",scentrality->String().Data()));
        nextPairCut.Reset();

        // Loop on each paircut (not the ones in MuMuConfig but the ones set)
        while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
        {
          AliDebug(1,Form("---PAIRCUT %s",spairCut->String().Data()));

          cout << "---------------------" << endl;
          cout << "Looking for spectras ..."<< endl;

          // Get spectra
          TString spectraPath ="";
          if(mix)
            spectraPath = Form("/FitResults/%s_%s/%s/%s/%s/%s/%s",refEvent.Data(),refTrigger.Data(),seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName);
          else
            spectraPath = Form("/FitResults/%s/%s/%s/%s/%s",seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName);
          if (AccEffCorr)spectraPath+="-AccEffCorr";

          AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));
          if(!spectra){
            AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
            return;
          }

          // Create pointer on fitted spectra. Any kind of capsule do the job
          AliAnalysisMuMuSpectraCapsulePbPb * capsule = new AliAnalysisMuMuSpectraCapsulePbPb(spectra,spectraPath,"","");
          if(!capsule){
            AliError("Could not find spetra !");
            return;
          }
          // Draw results
          TString sspectraName(spectraName);
          if(sspectraName.Contains("PSI-"))capsule->DrawResults(what,"PSI",subresults);
          delete capsule;
        }
      }
    }
  }


  delete eventTypeArray ;
  delete triggerArray ;
  delete pairCutArray ;
  delete centralityArray ;

  return ;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::PrintFitParam(TString spectraName, const char* subresult,const char* param)const
{
/**
 * @brief Draw fit parameters
 *
 * @param particle   particle name
 * @param param      Tokenized by ','
 * @param binType    [description]
 * @param subresult  If one wants to draw one/several specific results. Tokenized by ','
 * @param AccEffCorr  Spectra type
 */
  if (!OC() || !CC()){
    AliError("No mergeable/counter collection. Consider Upgrade()");
    return ;
  }

  // Get configuration settings
  TObjArray* eventTypeArray   = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
  TObjArray* triggerArray     = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
  TObjArray* fitfunctionArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());
  TObjArray* pairCutArray     = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
  TObjArray* centralityArray  = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
  TObjArray* paramArray       = TString(param).Tokenize(",");
  TObjArray* subresultArray   = TString(subresult).Tokenize(",");
  TObjArray* bins;

  // Iterater for loops
  TIter nextTrigger(triggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextparam(paramArray);
  TIter nextSubResult(subresultArray);
  TIter nextCentrality(centralityArray);

  // Strings
  TObjString* strigger;
  TObjString* seventType;
  TObjString* spairCut;
  TObjString* sparam;
  TObjString* ssubresult;
  TObjString* scentrality;

  AliAnalysisMuMuSpectra* spectra=0x0;
  TH1*       h          = 0x0;
  TH1*       hcent      = 0x0;
  TH1*       href       = 0x0;



  nextEventType.Reset();
  // Loop on each envenType (see MuMuConfig)
  while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
  {
    AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
    nextTrigger.Reset();

    // Loop on each trigger (see MuMuConfig)
    while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
    {
      AliDebug(1,Form("-TRIGGER %s",strigger->String().Data()));

      nextPairCut.Reset();
      // Loop on each paircut (not the ones in MuMuConfig but the ones set)
      while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
      {
        AliDebug(1,Form("---PAIRCUT %s",spairCut->String().Data()));
        nextCentrality.Reset();

        // Loop on each centrality (not the ones in MuMuConfig but the ones set)
        while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
        {
          AliDebug(1,Form("--CENTRALITY %s",scentrality->String().Data()));


          // Output message
          cout << "---------------------" << endl;
          cout << "Looking for spectras ..."<< endl;

          // Get spectra
          TString spectraPath= Form("/FitResults/%s/%s/%s/%s/%s",seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName.Data());

          AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));
          if(!spectra) {
            AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
            return;
          }

          TCanvas    * c = new TCanvas;
          int nparam = paramArray->GetEntries();
          c->Divide(2,nparam);

          // Loop on param
          int k = 1;
          nextparam.Reset();
          while ( ( sparam = static_cast<TObjString*>(nextparam()) ) )
          {
            AliDebug(1,Form("param %s",sparam->String().Data()));

            if(c)
            {
              // --- First canvas ---
              c->cd(k);

              TLegend*leg = new TLegend(0.1,0.7,0.48,0.9);
              leg->SetHeader(Form("Fit Parameters %s ",sparam->String().Data()));
              leg->SetTextSize(0.03);

              printf("going to subcanvas %d\n",k );
              int i = 1;
              nextSubResult.Reset();
              while ( ( ssubresult = static_cast<TObjString*>(nextSubResult()) ) )
              {
                //Loop over subresults
                AliDebug(1,Form("-----SubResults %s",ssubresult->String().Data()));
                if  ( !spectraName.Contains("VS") )
                  h = spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
                else if ( spectraName.Contains("YVSPT") )
                  h = static_cast<TH2*>(spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();
                else if ( spectraName.Contains("PTVSY") )
                  h = static_cast<TH2*>(spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();

                if(!h) {
                  AliError(Form("Cannot find histo for SubResults %s",ssubresult->String().Data()));
                  return;
                }

                // beautifull histo
                if( i!=3 && i!=5 && i!=10 && i!=11 && i!=12 && i!=13 && i!=14 ) h->SetMarkerColor(i); //nobody likes green and yellow
                else               h->SetMarkerColor(i+5);

                h->SetMarkerSize(1.);
                h->SetMarkerStyle(20+i);
                if(i==11)h->SetMarkerStyle(20+i+3);
                if(i==1)
                {
                  h->GetYaxis()->SetTitleSize(0.05);
                  h->GetYaxis()->SetLabelSize(0.05);
                  h->GetXaxis()->SetLabelSize(0.05);
                  h->GetXaxis()->SetTitleSize(0.05);
                  h->SetTitle(Form(" %s for bin %s",sparam->String().Data(),spectraName.Data()));
                  h->GetYaxis()->SetTitle(sparam->String().Data());
                  if (spectraName.Contains("YVSPT") )
                    h->GetXaxis()->SetTitle("PT");
                  else if (spectraName.Contains("PTVSY") )
                    h->GetXaxis()->SetTitle("Y");
                }

                if(! sparam->String().Contains("FitChi2PerNDF"))
                {
                  if(i==1)h->DrawCopy();
                  else    h->DrawCopy("same");
                }
                else
                {
                  if(i==1)h->DrawCopy("p");
                  else    h->DrawCopy("samep");
                }

                leg->AddEntry(h,Form("%s with %s",sparam->String().Data(),ssubresult->String().Data()),"p");
                i++;
              }

              leg->Draw("same");

              // --- ratio ---
              c->cd(++k);

              TLegend * leg2 = new TLegend(0.1,0.7,0.48,0.9);
              leg2->SetHeader(Form("Ratio"));
              leg2->SetTextSize(0.03);

              nextSubResult.Reset();
              int j= 1;
              TString refName;
              while ( ( ssubresult = static_cast<TObjString*>(nextSubResult()) ) )
              {
                AliDebug(1,Form("-----SubResults %s",ssubresult->String().Data()));

                if(j==1)
                {
                  href    = spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
                  if  ( !spectraName.Contains("VS") ) href = spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
                  else if ( spectraName.Contains("YVSPT") )
                    href = static_cast<TH2*>(spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();
                  else if ( spectraName.Contains("PTVSY") )
                    href = static_cast<TH2*>(spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();

                  refName = href->GetName();
                  j++;
                  continue;
                }

                if  ( !spectraName.Contains("VS") )
                  h = spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
                else if ( spectraName.Contains("YVSPT") )
                  h = static_cast<TH2*>(spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();
                else if ( spectraName.Contains("PTVSY") )
                  h = static_cast<TH2*>(spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE))->ProjectionX();

                if(!h || !href )
                {
                  AliError(Form("Cannot find histos for SubResults  ratio "));
                  return;
                }

                if( j!=3 && j!=5 && j!= 10 && j!=11 && j!=12 && j!=13 && j!=14 ) h->SetMarkerColor(j); //nobody likes green and yellow
                else               h->SetMarkerColor(j+5);
                h->SetMarkerSize(1.);
                h->SetMarkerStyle(20+j);
                if(j==11)h->SetMarkerStyle(20+j+3);
                if(j==2)
                {
                  h->GetYaxis()->SetTitleSize(0.05);
                  h->GetYaxis()->SetLabelSize(0.05);
                  h->GetXaxis()->SetLabelSize(0.05);
                  h->GetXaxis()->SetTitleSize(0.05);
                  h->SetTitle(Form(" %s Ratio over %s for %s",sparam->String().Data(),refName.Data(),spectraName.Data()));
                  if (spectraName.Contains("YVSPT") ) h->GetXaxis()->SetTitle("PT");
                  else if (spectraName.Contains("PTVSY") ) h->GetXaxis()->SetTitle("Y");
                }
                h->Divide(href);
                if(! sparam->String().Contains("FitChi2PerNDF"))
                {
                  if(j==2)h->DrawCopy();
                  else    h->DrawCopy("same");
                }
                else
                {
                  if(j==2)h->DrawCopy("p");
                  else    h->DrawCopy("samep");
                }
                leg2->AddEntry(h,Form("Results %d over Result 1",j),"pe");

                j++;
              }

              leg2->Draw("same");
              ++k;
            } else {
              nextSubResult.Reset();
              int i= 1;
              while ( ( ssubresult = static_cast<TObjString*>(nextSubResult()) ) ){
                //Loop over subresults
                AliDebug(1,Form("-----SubResults %s",ssubresult->String().Data()));
                h = spectra->Plot(sparam->String().Data(),ssubresult->String().Data(),kFALSE);
                if(!h) {
                  AliError(Form("Cannot find histo for SubResults %s",ssubresult->String().Data()));
                  return;
                }
              }
            }
          }
        }
      }
    }
  }
  delete eventTypeArray ;
  delete triggerArray ;
  delete pairCutArray ;
  delete centralityArray ;
  delete paramArray ;

  return ;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DivideRawMixHisto(const char* binType, const char* particle, const char* flavour, Bool_t corrected, Double_t Mmin, Double_t Mmax)
{
/**
 * @brief Divide raw histo with the related mix histo
 * @details To be run after NormMixEvent (run-by-run) on the merged file.
 *
 * Mixed Histograms are selected according to Mix...Key() in config. file. See AliAnalysisMuMuConfig class.
 *
 * @param binType    [description]
 * @param particle   particle name
 * @param flavour    [description]
 * @param corrected  If Minv are already AccxEff corrected or not
 */
    if(!OC())
    {
        AliError("No mergeable. Consider Upgrade()");
        return;
    }
    else
    {
        cout <<      " ================================================================ " << endl;
        cout <<      "                       DivideRawMixHisto                  " << endl;
        cout <<      " ================================================================ " << endl;
    }

    // Get configuration settings
    TObjArray* eventTypeArray    = Config()->GetListElements(Config()->RefMixEventSelectionKey(),IsSimulation());
    TObjArray* triggerArray      = Config()->GetListElements(Config()->RefMixTriggerKey(),IsSimulation());
    TObjArray* eventTypeMixArray = Config()->GetListElements(Config()->EventSelectionMixKey(),IsSimulation());
    TObjArray* triggerMixArray   = Config()->GetListElements(Config()->MixTriggerKey(),IsSimulation());
    TObjArray* pairCutArray      = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
    TObjArray* centralityArray   = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
    TObjArray* binTypeArray      = TString(binType).Tokenize(",");

    TObjArray* bins;

    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextTriggerMix(triggerMixArray);
    TIter nextEventType(eventTypeArray);
    TIter nextEventTypeMix(eventTypeMixArray);
    TIter nextPairCut(pairCutArray);
    TIter nextbinType(binTypeArray);
    TIter nextCentrality(centralityArray);

    // Strings
    TObjString* strigger;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* sbinType;
    TObjString* scentrality;
    TObjString* sTriggerMix;
    TObjString* seventTypeMix;


    // For the loop comming
    TString signFlagMinv[3] ={"","PlusPlus","MinusMinus"};
    TString signFlagDist[3] ={"","PP","MM"};

    THnSparse* n        =0x0;
    TObject* o          =0x0;

    // Get Binning
    AliAnalysisMuMuBinning* binning(0x0);

    // Loop on each envenTypeMix (see MuMuConfig)
    while ( ( seventTypeMix = static_cast<TObjString*>(nextEventTypeMix())) ){
      AliDebug(1,Form("EVENTTYPEMIX %s",seventTypeMix->String().Data()));
      nextTriggerMix.Reset();

      // Loop on each triggerMix (see MuMuConfig)
      while ( ( sTriggerMix = static_cast<TObjString*>(nextTriggerMix())) ){
        AliDebug(1,Form("-TRIGGERMIX %s",sTriggerMix->String().Data()));
        nextEventType.Reset();

        // Loop on each envenType (see MuMuConfig)
        while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
        {
          AliDebug(1,Form("--REFEVENTTYPE %s",seventType->String().Data()));
          nextTrigger.Reset();

          // Loop on each trigger (see MuMuConfig)
          while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
          {
            AliDebug(1,Form("---REFTRIGGER %s",strigger->String().Data()));
            nextPairCut.Reset();

            // Loop on each paircut (see MuMuConfig)
            while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
            {
              AliDebug(1,Form("----PAIRCUT %s",spairCut->String().Data()));
              nextbinType.Reset();

              // Loop on each type (see MuMuConfig)
              while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
              {
                AliDebug(1,Form("-----TYPE %s",sbinType->String().Data()));
                nextCentrality.Reset();

                // Loop on each centrality (see MuMuConfig)
                while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
                {
                  AliDebug(1,Form("------CENTRALITY %s",scentrality->String().Data()));
                  nextbinType.Reset();

                  // Loop on each centrality (see MuMuConfig)
                  while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
                  {
                    AliDebug(1,Form("-------BIN %s",sbinType->String().Data()));

                    // Pointers
                    TH1* hTableMinv[6] = {0x0,0x0,0x0,0x0,0x0,0x0};
    //                      {hMinvPM,hMinvPP,hMinvMM,hMinvMixPM,hMinvMixPP,hMinvMixMM};
                    TH1* hTableDistRaw[6] = {0x0,0x0,0x0,0x0,0x0,0x0};
    //                      {hpT,hpTPP,hpTMM,hY,hYPP,hYMM};
                    TH1* hTableDistMix[6] = {0x0,0x0,0x0,0x0,0x0,0x0};
    //                      {hpTMix,hpTMixPP,hpTMixMM,hYMix,hYMixPP,hYMixMM};

                    // Get binning
                    if ( fBinning && sbinType->String().Length() > 0 ) binning = fBinning->Project(particle,sbinType->String().Data(),flavour);
                    else  {
                      binning = new AliAnalysisMuMuBinning;
                      binning->AddBin(particle,sbinType->String().Data());
                    }
                    if (!binning) {
                      AliError("oups. binning is NULL");
                      continue;
                    }

                    // Create array
                    TObjArray* bins = binning->CreateBinObjArray(particle);
                    if (!bins){
                      AliError(Form("Did not get any bin for particle %s",particle));
                      delete binning;
                      continue;
                    }

                    // Create ID for the fit which will be used to name results
                    TString idMix(Form("/%s/%s/%s/%s",seventTypeMix->String().Data(),sTriggerMix->String().Data(),scentrality->String().Data(),spairCut->String().Data()));
                    AliDebug(1,Form("idMix = %s\n",idMix.Data() ));
                    TString idMinv(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()));
                    AliDebug(1,Form("idMinv = %s\n",idMinv.Data() ));
                    TString idDist(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()));
                    AliDebug(1,Form("idDist = %s\n",idDist.Data() ));

                    // The binning pointer, which point at Pt binning, Y binning etc.
                    AliAnalysisMuMuBinning::Range* bin;
                    TIter next(bins);

                    //MAIN PART : Loop on every binning range
                    while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
                    {
                      // Get Minv Histo
                      for (int j = 0; j < 3; ++j) {
                        TString hnameraw = corrected ? Form("MinvUS+%s%s_AccEffCorr",bin->AsString().Data(),signFlagMinv[j].Data()) : Form("MinvUS+%s%s",bin->AsString().Data(),signFlagMinv[j].Data());
                        TString hnamemix = corrected ? Form("MinvUS+%s%sMix_AccEffCorr",bin->AsString().Data(),signFlagMinv[j].Data()) : Form("MinvUS+%s%sMix",bin->AsString().Data(),signFlagMinv[j].Data());
                        // Pointer to the histo from histo collection (Yes it is discusting )
                        if ( OC()->Histo(idMinv.Data(),hnameraw.Data()) ) hTableMinv[j] = static_cast<TH1*>(OC()->Histo(idMinv.Data(),hnameraw.Data())->Clone());
                        else { AliError(Form("Could not find histo %s/%s",idMinv.Data(),hnameraw.Data())); continue ; }

                        if ( OC()->Histo(idMinv.Data(),hnamemix.Data()) ) hTableMinv[j+3] = static_cast<TH1*>(OC()->Histo(idMinv.Data(),hnamemix.Data())->Clone());
                        else { AliError(Form("Could not find histo %s/%s",idMinv.Data(),hnamemix.Data())); continue ; }
                      }

                      // Get Dist Histo
                      for (int j = 0; j < 3; ++j) {
                        TString hnamePt = corrected ? Form("Pt%s_AccEffCorr_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax) : Form("Pt%s_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax);
                        // TString hnamePt = corrected ? Form("Pt%s_AccEffCorr",signFlagDist[j].Data()) : Form("Pt%s",signFlagDist[j].Data());
                        TString hnameY  = corrected ? Form("Y%s_AccEffCorr_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax) : Form("Y%s_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax);
                        // TString hnameY  = corrected ? Form("Y%s_AccEffCorr",signFlagDist[j].Data()) : Form("Y%s",signFlagDist[j].Data());

                        TString hnamePtMix = corrected ? Form("PtMix%s_AccEffCorr_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax) : Form("PtMix%s_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax);
                        // TString hnamePtMix = corrected ? Form("PtMix%s_AccEffCorr",signFlagDist[j].Data()) : Form("PtMix%s",signFlagDist[j].Data());
                        TString hnameYMix  = corrected ? Form("YMix%s_AccEffCorr_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax) : Form("YMix%s_proj_0_%.2f_%.2f",signFlagDist[j].Data(),Mmin,Mmax);
                        // TString hnameYMix  = corrected ? Form("YMix%s_AccEffCorr",signFlagDist[j].Data()) : Form("YMix%s",signFlagDist[j].Data());

                        // Pointer to the histo from histo collection (Yes it is discusting )
                        if ( OC()->GetObject(idDist.Data(),hnamePt.Data()) ){
                          hTableDistRaw[j] = static_cast<TH1*>(OC()->Histo(idDist.Data(),hnamePt.Data())->Clone());
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idDist.Data(),hnamePt.Data()));
                          continue ;
                        }

                        if ( OC()->GetObject(idDist.Data(),hnameY.Data()) ) {
                          hTableDistRaw[j+3] = static_cast<TH1*>(OC()->Histo(idDist.Data(),hnameY.Data())->Clone());
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idDist.Data(),hnameY.Data()));
                          continue ;
                        }

                        if ( OC()->GetObject(idDist.Data(),hnamePtMix.Data()) ){
                          hTableDistMix[j] = static_cast<TH1*>(OC()->Histo(idDist.Data(),hnamePtMix.Data())->Clone());
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idDist.Data(),hnamePtMix.Data()));
                          continue ;
                        }

                        if ( OC()->GetObject(idDist.Data(),hnameYMix.Data()) ) {
                          hTableDistMix[j+3] = static_cast<TH1*>(OC()->Histo(idDist.Data(),hnameYMix.Data())->Clone());
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idDist.Data(),hnameYMix.Data()));
                          continue ;
                        }
                      }
                      // check if we have all the histo
                      Bool_t okAllHisto = kTRUE;
                      for (int i = 0; i < 3; ++i)
                      {
                        if ( !hTableMinv[i]    || !hTableMinv[i+3])   okAllHisto =kFALSE;
                        if ( !hTableDistRaw[i] || !hTableDistRaw[i+3])okAllHisto =kFALSE;
                        if ( !hTableDistMix[i] || !hTableDistMix[i+3])okAllHisto =kFALSE;
                      }
                      if(!okAllHisto) continue;

                      // --- remove background and make ratio ---

                      for (int i = 0; i < 3; ++i) {
                        if ( i==0 ) {
                          // Make errors by hands in case of empty entries
                          for (int j = 0; j < hTableMinv[i]->GetEntries(); ++j){
                            Double_t binContent = hTableMinv[i]->GetBinContent(j);
                            Double_t binError   = hTableMinv[i]->GetBinError(j);
                            if(binContent  == 0. && binError==0. ){
                              hTableMinv[i]->SetBinContent(j,0.);
                              hTableMinv[i]->SetBinError(j,1.);
                            }
                          }
                          hTableMinv[i]->Add(hTableMinv[i+3],-1); // Norm MinvMix histo
                          hTableMinv[i]->SetName(Form("%s_wbck",hTableMinv[i]->GetName())); // Norm MinvMix histo
                        }
                        else {
                          hTableMinv[i]->Divide(hTableMinv[i+3]); // Norm pt and rapidy mix histo
                          hTableMinv[i]->SetName(Form("%s_ratio",hTableMinv[i]->GetName())); // Norm MinvMix histo
                        }
                      }

                      for (int i = 0; i < 6; ++i) {
                        hTableDistRaw[i]->Divide(hTableDistMix[i]); // Norm MinvMix histo
                        hTableDistRaw[i]->SetName(Form("%s_ratio",hTableDistMix[i]->GetName())); // Norm MinvMix histo
                      }

                      // save results in mergeable collection
                      for (int i = 0; i < 6; ++i){
                        o = 0x0;
                        // Store mix histo
                        if ( i<3 ) {
                          o = fMergeableCollection->GetObject(idMinv.Data(),hTableMinv[i]->GetName());
                          AliDebug(1,Form("----idMinv=%s o=%p",idMinv.Data(),o));

                          if (o) fMergeableCollection->Remove(Form("%s/%s",idMinv.Data(),hTableMinv[i]->GetName()));

                          Bool_t adoptOK = fMergeableCollection->Adopt(idMinv.Data(),hTableMinv[i]);

                          if ( adoptOK ) std::cout << "\n" << "+++histo " << hTableMinv[i]->GetName() << " adopted in " << idMinv.Data() << std::endl;
                          else AliError(Form("Could not adopt histo %s",hTableMinv[i]->GetName()));
                        }

                        o = 0x0;
                        if ( hTableDistRaw[i] ) {
                          o = fMergeableCollection->GetObject(idDist.Data(),hTableDistRaw[i]->GetName());
                          AliDebug(1,Form("----idDist=%s o=%p",idDist.Data(),o));

                          if (o) fMergeableCollection->Remove(Form("%s/%s",idDist.Data(),hTableDistRaw[i]->GetName()));

                          Bool_t adoptOK = fMergeableCollection->Adopt(idDist.Data(),hTableDistRaw[i]);

                          if ( adoptOK ) std::cout << "\n" << "+++histo " << hTableDistRaw[i]->GetName() << " adopted in " << idDist.Data() << std::endl;
                          else AliError(Form("Could not adopt histo %s",hTableDistRaw[i]->GetName()));
                        }
                      }

                       // clean pointers
                       for (int i = 0; i < 6; ++i)
                       {
                          if(hTableMinv[i])hTableMinv[i] = 0x0;
                          if(hTableDistRaw[i])hTableDistRaw[i] = 0x0;
                          if(hTableDistMix[i])hTableDistMix[i] = 0x0;
                       }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    Update();

    delete eventTypeArray;
    delete eventTypeMixArray;
    delete triggerMixArray;
    delete pairCutArray;
    delete centralityArray;
    delete triggerArray;

    return;
}


//_____________________________________________________________________________
void AliAnalysisMuMu::DrawMinv(const char* type,const char* particle,const char* trigger,const char* eventType,const char* pairCut,const char* centrality,const char* subresultname,const char* flavour) const
{
  /**
  * @brief   Draw minv spectra
  *
  * @param type          'minv','mpt'...
  * @param particle      particle name
  * @param trigger       [description]
  * @param eventType     [description]
  * @param pairCut       [description]
  * @param centrality    [description]
  * @param subresultname If one wants to draw one/several specific results.
  * @param flavour       [description]
  */

  if (!OC() || !BIN()) return;

  TObjArray* bins = BIN()->CreateBinObjArray(particle,type,flavour);
  if (!bins)
  {
    AliError(Form("Could not get %s bins",type));
    return;
  }

  Double_t xmin(-1);
  Double_t xmax(-1);

  TString sparticle(particle);
  if ( sparticle=="PSI" )
  {
    xmin = 2;
    xmax = 6;
  }

  Int_t nx(1);
  Int_t ny(1);

  Int_t n = bins->GetEntries();

  if ( n == 2 )
  {
    nx = 2;
  }
  else if ( n > 2 )
  {
    ny = TMath::Nint(TMath::Sqrt(n));
    nx = n/ny;
  }

  TString stype(type);
  stype.ToUpper();

  TString spectraName(Form("/%s/%s/%s/%s/%s-%s",eventType,trigger,centrality,pairCut,particle,stype.Data()));

  if ( strlen(flavour))
  {
    spectraName += "-";
    spectraName += flavour;
  }

  AliAnalysisMuMuSpectra* spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraName.Data()));

  AliDebug(1,Form("spectraName=%s spectra=%p",spectraName.Data(),spectra));

  TObjArray* spectraBins(0x0);
  if ( spectra )
  {
    spectraBins = spectra->BinContentArray();
  }

  TCanvas* c = new TCanvas;
  c->Divide(nx,ny);
  c->Draw();
  gStyle->SetOptFit(1112);

  c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "AliAnalysisMuMu",
              (void*)this, "ExecuteCanvasEvent(Int_t,Int_t,Int_t,TObject*)");


  TIter next(bins);
  AliAnalysisMuMuBinning::Range* r;
  Int_t ci(0);

  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    TString name(Form("/%s/%s/%s/%s/MinvUS+%s",eventType,trigger,centrality,pairCut,r->AsString().Data()));

    AliDebug(1,name.Data());

    AliAnalysisMuMuJpsiResult* spectraBin(0x0);

    if ( spectraBins )
    {
      AliAnalysisMuMuResult* sr = static_cast<AliAnalysisMuMuResult*>(spectraBins->At(ci));

      spectraBin = static_cast<AliAnalysisMuMuJpsiResult*>(sr->SubResult(subresultname));

      AliDebug(1,Form("spectraBin(%s)=%p",subresultname,spectraBin));
    }

    TH1* h = OC()->Histo(name.Data());

    if ( spectraBin )
    {
      h = spectraBin->Histo();
    }

    if (h)
    {
      ++ci;
      c->cd(ci);
      gPad->SetLogy();
      if (xmin>0)
      {
        h->GetXaxis()->SetRangeUser(xmin,xmax);
      }
      h->Draw("histes");

      TObject* f1 = h->GetListOfFunctions()->FindObject("fitTotal");
      if (f1)
      {
        f1->Draw("same");
      }

      gPad->Modified();
      gPad->Update();

      TObject* stats = h->FindObject("stats");
      if (stats)
      {
        stats->Draw("same");
      }
    }
  }

  delete bins;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::DrawMinv(const char* type, const char* particle, const char* flavour, const char* subresultname) const
{
  /**
  * @brief   Draw minv spectra for binning of given type
  * @details Need to be reimplemented !
  *
  * @param type          'minv','mpt'...
  * @param particle      particle name
  * @param flavour       [description]
  * @param subresultname If one wants to draw one/several specific results.
  */

  if (!fConfig)
  {
    AliError("No configuration available yet. Don't know what to draw");
    return;
  }

  const AliAnalysisMuMuConfig& c = *(Config());

  DrawMinv(type,particle,
           c.First(c.DimuonTriggerKey(),IsSimulation()),
           c.First(c.EventSelectionKey(),IsSimulation()),
           c.First(c.PairSelectionKey(),IsSimulation()),
           c.First(c.CentralitySelectionKey(),IsSimulation()),
           subresultname,
           flavour);
}

//___________________________________________________________________
void AliAnalysisMuMu::ExecuteCanvasEvent(Int_t event, Int_t /*px*/, Int_t /*py*/, TObject *sel)
{
    // Actions in reponse to mouse button events.

    TCanvas* c = static_cast<TCanvas*>(gTQSender);
    TPad* pad = static_cast<TPad*>(c->GetSelectedPad());
    if (!pad) return;

    //  if ((event == kButton1Down) ||
    if (event == kButton1Double)
        {

        //    Float_t x = pad->AbsPixeltoX(px);
        //    Float_t y = pad->AbsPixeltoY(py);
        //    x = pad->PadtoX(x);
        //    y = pad->PadtoY(y);

        //    std::cout << "event=" << event << " px=" << px << " py=" << py << " ";

        if ( sel && sel->InheritsFrom("TH1") )
            {
            TCanvas* clocal = new TCanvas;
            clocal->SetLogy();
            clocal->Draw();
            sel->Draw();
            }
        else
            {
            TList* list = pad->GetListOfPrimitives();
            TIter next(list);
            TObject* h;

            while ( ( h = next() ) )
                {
                if ( h->InheritsFrom("TH1") )
                    {
                    TCanvas* clocal = new TCanvas;
                    clocal->SetLogy();
                    clocal->Draw();
                    h->Draw();
                    break;
                    }
                }

            }

        //      std::cout  << std::endl;

        pad->Modified();
        }

}

//_____________________________________________________________________________
void AliAnalysisMuMu::RAAasGraphic(const char* particle,const char* binType,const char* externfile,const char* externfile2,const char* RefCent,Bool_t print,Bool_t AccEffCorr)const
{
  /**
   * @brief Compute, store and print R_AA
   * @details Should be used after a fit process (FitJpsi() for instance). Work is delegated to a AliAnalysisMuMuSpectraCapsule class according to beam year.
   *
   * @param particle    particle name
   * @param binType     [description]
   * @param externfile  Config. file readed by AliAnalysisMuMuSpectraCapsule classes. See AliAnalysisMuMuSpectraCapsule documentation
   * @param externfile2 Config. file readed by AliAnalysisMuMuSpectraCapsule classes. See AliAnalysisMuMuSpectraCapsule documentation
   * @param RefCent     Centrality bin from which we get the number of trigger in the counter collection. V0M_00.00_00.90 by default
   * @param print       Print error details
   * @param AccEffCorr  Spectra type
  */
  if (!OC() || !CC())
      {
      AliError("No mergeable/counter collection. Consider Upgrade()");
      return ;
      }
  else
      {
      cout <<      " ================================================================ " << endl;
      cout <<      "                             Computing RAA                             " << endl;
      cout <<      " ================================================================ " << endl;
      }

  // Get configuration settings
  TObjArray* eventTypeArray   = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
  TObjArray* triggerArray     = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
  TObjArray* fitfunctionArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());// to add here an entry
  TObjArray* pairCutArray     = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
  TObjArray* centralityArray  = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
  TObjArray* particleArray    = TString(particle).Tokenize(",");
  TObjArray* binTypeArray     = TString(binType).Tokenize(",");
  TObjArray* bins;

  // Iterator for loops
  TIter nextTrigger(triggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextParticle(particleArray);
  TIter nextbinType(binTypeArray);
  TIter nextCentrality(centralityArray);

  // Strings
  TObjString* strigger;
  TObjString* seventType;
  TObjString* spairCut;
  TObjString* sparticle;
  TObjString* sbinType;
  TObjString* scentrality;

  // Pointers
  TGraphErrors* graph=0x0;
  TGraphErrors* graphErr=0x0;
  TGraphErrors* graphErrCorr=0x0;
  TGraphErrors* graphCent=0x0;
  TGraphErrors* graphCentErr=0x0;
  TList* list=0x0;
  AliAnalysisMuMuSpectra spectra=0x0;

  //Loop on particle type
  while ( ( sparticle = static_cast<TObjString*>(nextParticle()) ) )
      {
      AliDebug(1,Form("particle %s",sparticle->String().Data()));
      nextEventType.Reset();
      // Loop on each envenType (see MuMuConfig)
      //==============================================================================
      while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
          {
          AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
          nextTrigger.Reset();
          // Loop on each trigger (see MuMuConfig)
          //==============================================================================
          while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
              {
              AliDebug(1,Form("-TRIGGER %s",strigger->String().Data()));
              nextPairCut.Reset();
              // Loop on each paircut (not the ones in MuMuConfig but the ones set)
              //==============================================================================
              while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
                  {
                  AliDebug(1,Form("--PAIRCUT %s",spairCut->String().Data()));
                  nextbinType.Reset();
                  // Loop on each type (pt or y)
                  //==============================================================================
                  while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
                      {
                      AliDebug(1,Form("---TYPE %s",sbinType->String().Data()));

                      //canvas
                      TCanvas *c1 = new TCanvas;
                      c1->Draw();
                      //Divide canvas for pt and y bin
                      if (!sbinType->String().Contains("INTEGRATED"))
                      {
                        Int_t nx(1);
                        Int_t ny(1);
                        Int_t nofGraph = centralityArray->GetEntries(); // # of histo
                        if ( nofGraph == 2 ){
                          nx=2;
                          ny=0;
                        }
                        else if ( nofGraph > 2 ){
                          ny = TMath::Nint(TMath::Sqrt(nofGraph));
                          nx = TMath::Nint((nofGraph/ny) +0.6);
                        }
                        c1->Divide(nx,ny);
                      }
                      else if (sbinType->String().Contains("INTEGRATED")){
                        graphCent = new TGraphErrors(centralityArray->GetEntries());
                        graphCent->SetMinimum(0.);
                        graphCent->SetMaximum(1.2);
                        graphCentErr = new TGraphErrors(centralityArray->GetEntries());
                        graphCentErr->SetMinimum(0.);
                        graphCentErr->SetMaximum(1.2);
                      }
                      gStyle->SetOptStat(0);

                      Int_t n=1; // counter
                      nextCentrality.Reset();
                      // Loop on each centrality (not the ones in MuMuConfig but the ones set)
                      //==============================================================================
                      while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) ) {
                          AliDebug(1,Form("---CENTRALITY %s",scentrality->String().Data()));

                          //________Get spectra
                          TString spectraPath= Form("/%s/%s/%s/%s/%s-%s",seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data(),sparticle->String().Data(),sbinType->String().Data());
                          if (AccEffCorr)spectraPath+="-AccEffCorr";

                          AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));
                          if(!spectra){//protection
                            AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
                            return;
                          }
                          //________

                          //________Get Trigger sum
                          Int_t NofMUL = TMath::Nint(CC()->GetSum(Form("trigger:%s/event:%s/centrality:%s",strigger->String().Data(),seventType->String().Data(),RefCent)));
                          //________

                          AliAnalysisMuMuSpectraCapsulePbPb * capsule = new AliAnalysisMuMuSpectraCapsulePbPb(spectra,spectraPath,externfile,externfile2);
                          if(!capsule) continue;
                          AliDebug(1,Form("Spectra = %p",capsule));
                          if(print)capsule->SetPrintFlag();

                          // Get Graph with RAA results
                          list = capsule->RAAasGraphic(NofMUL);

                          if(!list) continue;
                          AliDebug(1,Form("list = %p",list));

                          // Draw Graph according to bin type
                          if (!sbinType->String().Contains("INTEGRATED")){ // PT/Y/single centrality
                            // Select subcanvas
                            c1->cd(n);
                            //legend
                            TLegend * leg = new TLegend(0.2,0.8,0.90,0.9);
                            leg->SetTextSize(0.04);
                            leg->SetHeader(Form("ALICE, Pb-Pb #sqrt{s_{NN}}=5.02 TeV, L_{int}=222 #mub^{-1}, %s",scentrality->String().Data()));
                            //Draw it
                            graph = static_cast<TGraphErrors*>(list->At(0)->Clone());
                            graphErr = static_cast<TGraphErrors*>(list->At(1)->Clone());
                            graphErrCorr = static_cast<TGraphErrors*>(list->At(2)->Clone());
                            if (!graph || ! graphErr) {
                                printf("Did not find graph in the list");
                                return;
                            }

                            leg->AddEntry(graph,"Inclusive J/#psi","pe");
                            leg->AddEntry(graphErr,"Syst. uncertainty","f");

                            //Graph1
                            graph->GetYaxis()->SetRangeUser(0,1.41);
                            graph->Draw("ap");
                            // Graph Syst.
                            graphErr->SetFillColorAlpha(4,0.);
                            graphErr->Draw("same5");
                            leg->Draw();

                            // Global box
                            Double_t global = graphErrCorr->GetErrorY(0)/100.;
                            TBox *globalBox =0x0;
                            if(sbinType->String().Contains("PT"))globalBox= new TBox(11.5,1.-global,12,1+global);
                            if(sbinType->String().Contains("Y"))globalBox= new TBox(-2.6,1.-global,-2.5,1+global);
                            globalBox->SetFillColor(4);
                            globalBox->Draw("F");
                            printf("Global syst : %f\n",global );

                            TLine *l = 0x0;
                            if(sbinType->String().Contains("PT"))l = new TLine(0.,1.,12.,1.);
                            if(sbinType->String().Contains("PT"))l = new TLine(-4.0,1.,-2.5,1.);
                            l->SetLineStyle(2);
                            l->Draw("same");

                          }

                          else {// Integrated case

                            if(!sbinType->String().Contains("INTEGRATED")) {// Protection
                              cout << "Cannot plot INTEGRATED  ! Check it please :) " << endl;
                              delete c1;
                              return;
                            }

                            Double_t x=0;
                            Double_t y=0;
                            graph = static_cast<TGraphErrors*>(list->At(0)->Clone());
                            graphErr = static_cast<TGraphErrors*>(list->At(1)->Clone());
                            graphErrCorr = static_cast<TGraphErrors*>(list->At(2)->Clone());

                            if (!graph || ! graphErr) {// Protection
                                printf("Did not find graph in the list");
                                return;
                            }

                            // Get point for each centrality
                            Double_t dx = graph->GetErrorX(0);
                            Double_t dy = graph->GetErrorY(0);
                            Double_t dysys = graphErr->GetErrorY(0);
                            graph->GetPoint(0,x,y);

                            // Set them to a new graphic
                            graphCent->SetPoint(n-1,x,y);
                            graphCent->SetPointError(n-1,dx,dy);
                            graphCent->SetTitle(graph->GetTitle());
                            graphCentErr->SetPoint(n-1,x,y);
                            graphCentErr->SetPointError(n-1,2.5,dysys);
                          }
                        n++;
                        delete capsule;
                      }
                      cout << "" << endl;
                      if (sbinType->String().Contains("INTEGRATED")){ //Print

                        graphCent->GetXaxis()->SetTitle("<N_{part}>");
                        graphCent->GetYaxis()->SetTitle("R_{AA}");
                        graphCent->GetYaxis()->SetRangeUser(0,1.5);
                        graphCent->SetMarkerColor(4);
                        graphCent->SetMarkerStyle(21);
                        graphCentErr->SetFillColorAlpha(4,0.);
                        graphCent->SetTitle(Form("%s/%s/%s/%s/%s",seventType->String().Data(),
                                                 strigger->String().Data(),
                                                 spairCut->String().Data(),
                                                 sparticle->String().Data(),
                                                 sbinType->String().Data()));
                        TLegend * leg = new TLegend(0.2,0.8,0.90,0.9);
                        leg->SetHeader("ALICE, Pb-Pb #sqrt{s_{NN}}=5.02 TeV, L_{int}=222 #mub^{-1}, PT/Y integrated");
                        leg->SetTextSize(0.03);
                        leg->AddEntry(graphCent,"Inclusive J/#psi","pe");
                        leg->AddEntry(graphCentErr,"Syst. uncertainty","f");

                        graphCent->Draw("ap");
                        graphCentErr->Draw("same5");
                        leg->Draw();

                        // Global box
                        Double_t global = graphErrCorr->GetErrorY(0)/100.;
                        TBox *globalBox = globalBox= new TBox(420,1.-global,430,1+global);
                        globalBox->SetFillColor(4);
                        globalBox->Draw("F");
                        printf("Global syst : %f\n",global );

                        TLine *l = new TLine(0.,1.,12.,1.);
                        l->SetLineStyle(2);
                        l->Draw("same");
                      }
                      //________ Update results in Mergeable collection
                      TString id(Form("/RAA-%s/%s/%s/%s/%s",strigger->String().Data(),seventType->String().Data(),spairCut->String().Data(),sbinType->String().Data(),sparticle->String().Data()));
                      TObject* o = 0x0;

                      if (graph){// first graph

                        o = fMergeableCollection->GetObject(Form("%s/%s",id.Data(),graph->GetName()));
                        if (o){
                          AliWarning(Form("Replacing %s/%s",id.Data(),graph->GetName()));
                          fMergeableCollection->Remove(Form("%s/%s",id.Data(),graph->GetName()));
                        }

                        Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),graph);

                        if ( adoptOK ) std::cout << "+++RAA graph " << graph->GetName() << " adopted" << std::endl;
                        else AliError(Form("Could not adopt RAA grap %s",graph->GetName()));
                      }
                      if (graphCent){// second graph
                        o = fMergeableCollection->GetObject(Form("%s/%s",id.Data(),graphCent->GetName()));
                        if (o){
                          AliWarning(Form("Replacing %s/%s",id.Data(),graphCent->GetName()));
                          fMergeableCollection->Remove(Form("%s/%s",id.Data(),graphCent->GetName()));
                        }

                        Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),graphCent);

                        if ( adoptOK ) std::cout << "+++RAA graph " << graphCent->GetName() << " adopted" << std::endl;
                        else AliError(Form("Could not adopt RAA graph %s",graphCent->GetName()));
                      }
                      //________
                      }
                  }
              }
          }
      }
  delete list;
  delete eventTypeArray ;
  delete triggerArray ;
  delete fitfunctionArray ;
  delete pairCutArray ;
  delete centralityArray ;
  delete particleArray ;
  delete binTypeArray ;

  return ;
}

//_____________________________________________________________________________
TString AliAnalysisMuMu::ExpandPathName(const char* file)
{
/**
 * @brief An expand method that lives alien URL as they are
 *
 * @param file [description]
 */
  TString sfile;

  if ( !sfile.BeginsWith("alien://") )
  {
    return gSystem->ExpandPathName(file);
  }
  else
  {
    if (!gGrid) TGrid::Connect("alien://");
    if (!gGrid) return "";
  }

  return file;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::TwikiOutputFnorm(const char* series) const
{
/**
 * @brief  Make a twiki-compatible output of the Fnorm factor(s) ( Outdated ?)
 *
 * @param series [description]
 */
  TObjArray* what = TString(series).Tokenize(",");
  TObjString* s;
  TObjArray graphs;
  TIter next(what);

  std::cout << "| *Run* |";
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TGraph* g = static_cast<TGraph*>(OC()->GetObject(Form("/FNORM/GRAPHS/%s",s->String().Data())));
    if (!g)
    {
      AliError(Form("Could not find graph for %s",s->String().Data()));
      continue;
    }
    std::cout << " *" << s->String().Data();
    if ( s->String().BeginsWith("RelDif") ) std::cout << " %";
    std::cout << "*|";
    graphs.Add(g);
  }

  std::cout << std::endl;

  TGraphErrors* g0 = static_cast<TGraphErrors*>(graphs.First());
  if (!g0) return;

  for ( Int_t i = 0; i < g0->GetN(); ++i )
  {
    TString msg;

    msg.Form("|%6d|",TMath::Nint(g0->GetX()[i]));

    for ( Int_t j = 0; j < graphs.GetEntries(); ++j )
    {
      TGraphErrors* g = static_cast<TGraphErrors*>(graphs.At(j));

      msg += TString::Format(" %6.2f +- %6.2f |",g->GetY()[i],g->GetEY()[i]);
    }

    std::cout << msg.Data() << std::endl;
  }

  next.Reset();

  std::cout << "|*Weigthed mean (*)*|";

  AliAnalysisMuMuResult* r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/Fnorm"));

  if (!r)
  {
    AliError("Could not find Fnorm result !");
    return;
  }


  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TString var("Fnorm");
    TString unit;

    if ( s->String().BeginsWith("Fnorm") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/Fnorm"));
    }
    else if ( s->String().BeginsWith("RelDif") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/RelDif"));
      unit = "%";
    }

    r->Exclude("*");
    r->Include(s->String().Data());

    std::cout << Form(" * %5.2f +- %5.2f %s * |",
                      r->GetValue(var.Data()),
                      r->GetErrorStat(var.Data()),
                      unit.Data());
  }

  next.Reset();

  std::cout << std::endl;

  std::cout << "|*RMS*|";

  while ( ( s = static_cast<TObjString*>(next())) )
  {
    TString var("Fnorm");

    if ( s->String().BeginsWith("Fnorm") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/Fnorm"));
    }
    else if ( s->String().BeginsWith("RelDif") )
    {
      r = static_cast<AliAnalysisMuMuResult*>(OC()->GetObject("/FNORM/RESULTS/RelDif"));
    }

    r->Exclude("*");
    r->Include(s->String().Data());

    Double_t d = 100.0*r->GetRMS(var.Data())/r->GetValue(var.Data());

    std::cout << Form(" * %5.2f (%5.2f %%) * |",
                      r->GetRMS(var.Data()),d);
  }

  std::cout << std::endl;
  std::cout << "(*) weight is the number of CMUL7-B-NOPF-MUON triggers (physics-selected and pile-up corrected) in each run" << std::endl;

  delete what;
}

//_____________________________________________________________________________
TFile* AliAnalysisMuMu::FileOpen(const char* file)
{
/**
 * @brief Open a file after expansion of its name
 *
 * @param file [description]
 * @return [description]
 */
    return TFile::Open(ExpandPathName(file).Data());
}

//_____________________________________________________________________________
TString AliAnalysisMuMu::First(const TString& list) const
{
/**
 * @brief [brief description]
 * @details [long description]
 *
 * @param list [description]
 * @return [description]
 */
    TObjArray* a = list.Tokenize(",");
    if ( a->GetLast() < 0 ) return "";

    TString rv = static_cast<TObjString*>(a->First())->String();

    delete a;

    return rv;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::FitParticle(const char* particle,const char* trigger,const char* eventType,const char* pairCut,const char* centrality,const AliAnalysisMuMuBinning& binning,Bool_t corrected,const TString* fitMethod,const char* flavour,const char* histoType)
{
/**
 * @brief Fit minv histo
 * @details Create and store an AliAnalysisMuMuResult classes ( for the moment, only AliAnalysisMuMuJpsiResult is used and implemented)
 * which is the class who actually does process the fit.
 *
 * @param particle     particle name
 * @param trigger     [description]
 * @param eventType   [description]
 * @param pairCut     [description]
 * @param centrality  [description]
 * @param binning     [description]
 * @param corrected   [description]
 * @param fitMethod   '' or 'mix'
 * @return            AliAnalysisMuMuSpectra to be handled by the owner
 */
  // To avoid bins with error=0 due to low statstics
  TProfile::Approximate();

  static int n(0);

  Bool_t mix = kFALSE;
  TString* id(0x0);
  TString* refTrigger(0x0);
  TString* refEvent(0x0);
  TObjString* fitType;

  // Get number of runs and store it in nruns
  TObjArray   * runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  Int_t nruns = runs->GetEntries();
  delete runs;

  // ---- Some cross check  with the counter collection----

  // Check Binning list
  TObjArray* bins = binning.CreateBinObjArray(particle);
  if (!bins){
    AliError(Form("Did not get any bin for particle %s",particle));
    return 0x0;
  }

  // to avoid case sensitive stuff
  TString seventType(eventType);
  if(seventType.Contains("PSINT7inMUON")) seventType = "PSINT7INMUON";

  // Check trigger list
  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  if ( !triggers->FindObject(trigger) ){
    AliError(Form("Did not find trigger %s",trigger));
    delete bins;
    delete triggers;
    return 0x0;
  }
  delete triggers;

  // Check event list
  TObjArray* events = fCounterCollection->GetKeyWords("event").Tokenize(",");
  if ( !events->FindObject(seventType.Data()) ){
    AliError(Form("Did not find eventType %s",seventType.Data()));
    delete bins;
    delete events;
    return 0x0;
  }
  delete events;

  Int_t ntrigger = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s",trigger,seventType.Data())));

  // Check trigger
  if  ( ntrigger<=0 ){
    AliError(Form("No trigger for trigger:%s/event:%s",trigger,seventType.Data()));
    delete bins;
    return 0x0;
  }

  // ---- Here we select triggers and the fit method (mix or not) ----

  // Set mix flag
  if ( fitMethod->Contains("mix") && !IsSimulation() ) mix = kTRUE;

  // Create ID for the fit which will be used to name results
  if ( !mix ) id = new TString(Form("/%s/%s/%s/%s",eventType,trigger,centrality,pairCut));
  else {
    refTrigger = new TString(Form("%s",Config()->First(Config()->RefMixTriggerKey(),IsSimulation()).Data()));
    refEvent   = new TString(Form("%s",Config()->First(Config()->RefMixEventSelectionKey(),IsSimulation()).Data()));
    id         = new TString(Form("/MIX/%s_%s/%s/%s/%s/%s",refEvent->Data(),refTrigger->Data(),eventType,trigger,centrality,pairCut));
  }

  // The result pointer, will be return at the end
  AliAnalysisMuMuSpectra* spectra(0x0);
  TString spectraName(binning.GetName());
  TString sFlavour(flavour);
  if( !sFlavour.IsNull()) spectraName += Form("-%s",flavour);
  if ( corrected )   spectraName += "-AccEffCorr";

  // ---- MAIN PART : Loop on every binning range ----

  AliAnalysisMuMuBinning::Range* bin;
  TIter next(bins);
  next.Reset();
  while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
  {
    Int_t added(0);
    AliAnalysisMuMuJpsiResult* r    = 0x0;
    Bool_t adoptOk           = kFALSE;
    Bool_t adoptMix          = kFALSE;

    TH1* histo(0x0);

    // ---- Here we select the histo name we want/need to proceed the fit ----

    // Get fitType as a string
    TString sHistoType(histoType);

    // Select name histo
    TString hname;
    TString mixflag1 = mix ? "_wbck" : "" ;
    if( sHistoType.Contains("minv"))
      hname = corrected ? Form("MinvUS_AccEffCorr+%s%s",bin->AsString().Data(),mixflag1.Data())  : Form("MinvUS+%s%s",bin->AsString().Data(),mixflag1.Data());
    else if( sHistoType.Contains("mpt") && !sHistoType.Contains("mpt2") )
      hname = corrected ? Form("MeanPtVsMinvUS_AccEffCorr+%s%s",bin->AsString().Data(),mixflag1.Data()) : Form("MeanPtVsMinvUS+%s%s",bin->AsString().Data(),mixflag1.Data());
    else if( sHistoType.Contains("mpt2") )
      hname = corrected ? Form("MeanPtSquareVsMinvUS_AccEffCorr+%s%s",bin->AsString().Data(),mixflag1.Data()) : Form("MeanPtSquareVsMinvUS+%s%s",bin->AsString().Data(),mixflag1.Data());
    else {
      AliError("Wrong spectra type choice: Possibilities are: 'minv' or 'mpt' ");
      continue;
    }

    // Print the fitting process on the terminal
    TString isCorr(corrected ? " AccEffCorr " : " ");
    isCorr += ( mix )  ? " with mixing event method " : " ";
    std::cout << "---------------------------------//---------------------------------" << std::endl;
    std::cout << "Fitting" << isCorr.Data() << sHistoType.Data() << " spectra in " << id->Data() << std::endl;

    // Finally gets it
    if ( OC()->Histo(id->Data(),hname.Data()) ) histo = static_cast<TH1*>(OC()->Histo(id->Data(),hname.Data())->Clone(Form("%s%d",sHistoType.Data(),n++)));
    if ( !histo ) {
      AliError(Form("Could not find histo %s/%s",id->Data(),hname.Data()));
      continue;
    }

    // At some point particleTmp should become particle (but for now particle is always = "psi")
    const char* particleTmp = IsSimulation() ? GetParticleName() : "JPsi";
    cout << "particleTmp =" << particleTmp << endl;
    TString sparticleTmp(particleTmp);

    // Create the AliAnalysisMuMuResults that will fit the JPsi
    r = new AliAnalysisMuMuJpsiResult(particleTmp,*histo,trigger,eventType,pairCut,centrality,*bin);
    if ( !r ){
      AliError("Cannot create a AliAnalysisMuMuJpsiResult");
      continue;
    }
    r->SetNofTriggers(ntrigger);
    r->SetNofRuns(nruns);

    // Create an array (fitTypeArray) pointing on AliAnalysisMuMuConfig and store kFitTypeList.  Also create pointers and strings for several pointers
    TObjArray* fitTypeArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());
    TIter nextFitType(fitTypeArray);  // Iterater for every fit types, i.e fitting functions and their config.
    nextFitType.Reset();

    // Loop on every fittype and create a subresult inside the spectra.
    while ( ( fitType = static_cast<TObjString*>(nextFitType())) )
    {
      AliDebug(1,Form("<<<<<< fitType=%s bin=%s",fitType->String().Data(),bin->Flavour().Data()));

      std::cout << "" << std::endl;
      std::cout << "---------------" << "Fit " << added + 1 << "------------------" << std::endl;
      if(!mix) std::cout << "Fitting " << hname.Data() << " with " << fitType->String().Data() << std::endl;
      else     std::cout << "Fitting " << hname.Data() << " with " << fitType->String().Data() << " and after remmoving backround from mixing " << std::endl;
      std::cout << "" << std::endl;

      // Look for specific fit param.
      TObjArray* fitSingle = Config()->GetListElements(Config()->FitSingleKey(),IsSimulation());
      TIter nextFitSingle(fitSingle);
      TObjString* specifit;
      nextFitSingle.Reset();
      while ( ( specifit = static_cast<TObjString*>(nextFitSingle())) ){
        // must find a better way to do ...
        // Here we tokenize the initial FitType string with ":" and see if function names and ranges match

        if ( !specifit->String().Contains(bin->AsString().Data()) ) continue; // Check binning

        TString func   ="";
        TString range  ="";
        TString Weight ="";
        TString mctails="";

        TObjArray* oldFitParam = fitType->String().Tokenize(":");
        TIter nextoldFitParam(oldFitParam);
        TObjString* param;
        nextoldFitParam.Reset();
        while ( ( param = static_cast<TObjString*>(nextoldFitParam())) ){
          if ( param->String().Contains("func=") )    func    = param->String().Data();
          if ( param->String().Contains("range=") )   range   = param->String().Data();
          if ( param->String().Contains("weight=") )  Weight  = param->String().Data();
          if ( param->String().Contains("mctails") )  mctails = param->String().Data();
        }
        if ( specifit->String().Contains(func.Data())
          && specifit->String().Contains(range.Data())
          && specifit->String().Contains(Weight.Data())
          && specifit->String().Contains(mctails.Data()) ) fitType->String() = specifit->String().Data();

        delete oldFitParam;
      }
      delete fitSingle;

      if(  mix != fitType->String().Contains("mix")  )  {
        printf("skip %s because inconsistant with FitMethod \n",fitType->String().Data() );
        continue;
      }

      // Conf. for MC Tails (see function type)
      if ( fitType->String().Contains("mctails",TString::kIgnoreCase) || fitType->String().Contains("mctails2",TString::kIgnoreCase) ){

        TString sbin                          = bin->AsString();
        TString spectraMCName                 = spectraName;
        AliAnalysisMuMuBinning::Range* binMC  = bin;

        // Javier's Legacy
        if( (sbin.Contains("MULT") || sbin.Contains("NCH") || sbin.Contains("DNCHDETA") || sbin.Contains("V0A") || sbin.Contains("V0ACENT") || sbin.Contains("V0C") || sbin.Contains("V0M") || sbin.Contains("NTRCORR")|| sbin.Contains("RELNTRCORR")) && !sbin.Contains("NTRCORRPT") && !sbin.Contains("NTRCORRY")){

          //-------has to have a better way to do it
          AliAnalysisMuMuBinning* b = new AliAnalysisMuMuBinning;
          b->AddBin("psi","INTEGRATED");

          binMC = static_cast<AliAnalysisMuMuBinning::Range*>(b->CreateBinObjArray()->At(0));

          spectraMCName = b->GetName();
          delete b;

          if ( corrected ){
            spectraMCName += "-";
            spectraMCName += "AccEffCorr";
          }
        }

        Bool_t okMCtails = kFALSE;
        if( !fitType->String().Contains("momo",TString::kIgnoreCase) )
          okMCtails = GetParametersFromMC(fitType->String(),Form("/%s/%s",centrality,pairCut),spectraMCName.Data(),binMC);
        else
          okMCtails = GetParametersFromMC(fitType->String(),Form("/PP/%s",pairCut),spectraMCName.Data(),binMC);

        if(!okMCtails) continue;

        added += ( r->AddFit(fitType->String().Data()) == kTRUE );
      }

      // Config. for mpt (see function type)
      else if ( fitType->String().Contains("histoType=mpt",TString::kIgnoreCase) && !fitType->String().Contains("histoType=minv",TString::kIgnoreCase) && !fitType->String().Contains("MPTPSI_HFUNCTION",TString::kIgnoreCase) ){

        std::cout << "++The Minv parameters will be taken from " << spectraName.Data() << std::endl;
        std::cout << "" << std::endl;

        AliAnalysisMuMuSpectra* minvSpectra = dynamic_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/FitResults%s",id->Data()),spectraName.Data()));

        if ( !minvSpectra ){
          AliError(Form("Cannot fit mean pt: could not get the minv spectra for /FitResults%s",id->Data()));
          continue;
        }

        AliAnalysisMuMuJpsiResult* minvResult = static_cast<AliAnalysisMuMuJpsiResult*>(minvSpectra->GetResultForBin(*bin));

        if ( !minvResult ){
          AliError(Form("Cannot fit mean pt: could not get the minv result for bin %s in /FitResults%s",bin->AsString().Data(),id->Data()));
          continue;
        }

        TObjArray* minvSubResults = minvResult->SubResults();
        TIter nextSubResult(minvSubResults);
        AliAnalysisMuMuJpsiResult* fitMinv;
        TString subResultName;

        Int_t nSubFit(0);
        while ( ( fitMinv = static_cast<AliAnalysisMuMuJpsiResult*>(nextSubResult())) ){

          TString fitMinvName(fitMinv->GetName());
          fitMinvName.Remove(fitMinvName.First("_"),fitMinvName.Sizeof()-fitMinvName.First("_"));

          if ( !fitType->String().Contains(fitMinvName) ) continue;

          std::cout << "" << std::endl;
          std::cout <<  "      /-- SubFit " << nSubFit + 1 << " --/ " << std::endl;
          std::cout << "" << std::endl;

          TString sMinvfitType(fitType->String());

          GetParametersFromResult(sMinvfitType,fitMinv);//FIXME: Think about if this is necessary

          added += ( r->AddFit(sMinvfitType.Data()) == kTRUE );

          nSubFit++;
        }
      }

      // Config. for mpt (see function type)
      else if ( fitType->String().Contains("histoType=mpt2",TString::kIgnoreCase) && !fitType->String().Contains("histoType=minv",TString::kIgnoreCase) && !fitType->String().Contains("MPTPSI_HFUNCTION",TString::kIgnoreCase) ){

        std::cout << "++The Minv parameters will be taken from " << spectraName.Data() << std::endl;
        std::cout << "" << std::endl;

        AliAnalysisMuMuSpectra* minvSpectra = dynamic_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("/FitResults%s",id->Data()),spectraName.Data()));

        if ( !minvSpectra ){
          AliError(Form("Cannot fit mean pt: could not get the minv spectra for /FitResults%s",id->Data()));
          continue;
        }

        AliAnalysisMuMuJpsiResult* minvResult = static_cast<AliAnalysisMuMuJpsiResult*>(minvSpectra->GetResultForBin(*bin));

        if ( !minvResult ){
          AliError(Form("Cannot fit mean pt: could not get the minv result for bin %s in /FitResults%s",bin->AsString().Data(),id->Data()));
          continue;
        }

        TObjArray* minvSubResults = minvResult->SubResults();
        TIter nextSubResult(minvSubResults);
        AliAnalysisMuMuJpsiResult* fitMinv;
        TString subResultName;

        Int_t nSubFit(0);
        while ( ( fitMinv = static_cast<AliAnalysisMuMuJpsiResult*>(nextSubResult())) ){

          TString fitMinvName(fitMinv->GetName());
          fitMinvName.Remove(fitMinvName.First("_"),fitMinvName.Sizeof()-fitMinvName.First("_"));

          if ( !fitType->String().Contains(fitMinvName) ) continue;

          std::cout << "" << std::endl;
          std::cout <<  "      /-- SubFit " << nSubFit + 1 << " --/ " << std::endl;
          std::cout << "" << std::endl;

          TString sMinvfitType(fitType->String());

          GetParametersFromResult(sMinvfitType,fitMinv);//FIXME: Think about if this is necessary

          added += ( r->AddFit(sMinvfitType.Data()) == kTRUE );

          nSubFit++;
        }
      }

      //Config. for the rest (see function type)
      else {

        if ( fitType->String().Contains("PSICB2",TString::kIgnoreCase) || fitType->String().Contains("PSINA60NEW",TString::kIgnoreCase))
          std::cout << "+Free tails fit... " << std::endl;
        else if ( fitType->String().Contains("PSICOUNT",TString::kIgnoreCase) )
          std::cout << Form("+Just counting %s...",GetParticleName()) << std::endl;
        else
          std::cout << "+Using predefined tails... " << std::endl;

        if ( fitType->String().Contains("minvJPsi") && !sparticleTmp.Contains("JPsi") ){
          std::cout << "This fitting funtion is set to fit JPsi: Skipping fit..." << std::endl;
          continue;
        }

        if ( fitType->String().Contains("minvPsiP") && !sparticleTmp.Contains("PsiP") ){
          std::cout << "This fitting funtion is set to fit PsiP: Skipping fit..." << std::endl;
          continue;
        }
        // Here we call  FINALLY the fit functions
        added += ( r->AddFit(fitType->String().Data()) == kTRUE );
      }

      std::cout << "-------------------------------------" << std::endl;
      std::cout << "" << std::endl;
    }
    if ( !added )
    {
      delete fitTypeArray;
      continue;
    }

    // Get <flavour>
    flavour = bin->Flavour();

    // Implement <spectra> and set its name
    if (!spectra){

      TString spectraSaveName = spectraName;

      // Check if we fit meanPt
      nextFitType.Reset();
      Bool_t meanptVSminvFlag = kFALSE;
      Bool_t meanpt2VSminvFlag = kFALSE;
      while ( ( fitType = static_cast<TObjString*>(nextFitType())) ){
        meanpt2VSminvFlag = fitType->String().Contains("histoType=mpt2");
        if(!meanpt2VSminvFlag)meanptVSminvFlag  = fitType->String().Contains("histoType=mpt");
      }
      if ( meanptVSminvFlag){
        spectraSaveName += "-";
        spectraSaveName += "MeanPtVsMinvUS";
      }
      if ( meanpt2VSminvFlag){
        spectraSaveName += "-";
        spectraSaveName += "MeanPtSquareVsMinvUS";
      }
      spectra = new AliAnalysisMuMuSpectra(spectraSaveName.Data());
    }

    // We adopt the Result for current bin into the spectra
    if ( r ) adoptOk = spectra->AdoptResult(*bin,r);

    if ( r && spectra && adoptOk ) printf("Result %s adopted in spectra %s  \n",r->GetName(),spectra->GetName());
    else AliError(Form("Error adopting result "));

    if ( IsSimulation() ) {
      std::cout << "Computing AccEff Value Spectra " << std::endl;
      SetNofInputParticles(*r,eventType,trigger,centrality);
    }
    delete fitTypeArray;
  }

  delete bins;
  if (refTrigger) delete  refTrigger;
  if (refEvent)   delete  refEvent;
  delete  id;

  return spectra;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::GetParametersFromMC(TString& fitType, const char* pathCentrPairCut, const char* spectraName, AliAnalysisMuMuBinning::Range* bin) const
{
/**
 * @brief   Get fit tails parameters from the associated simulations file.
 * @details Intended to be used in minv fits, where we need the tails from JPsi (and Psi').
 *
 * @param fitType          'NA60NEW','CB2','INDEPTAILS' implemented at the moment
 * @param pathCentrPairCut [description]
 * @param spectraName      [description]
 * @param bin              [description]
 */
    if ( !SIM() && !SIM2() )
    {
      AliError("Cannot get MC tails without associated simulation(s) file(s) !");
      fitType = "";
      return kFALSE;
    }

    TString subResultName("");
    if ( fitType.Contains("NA60NEW",TString::kIgnoreCase) )  subResultName = "PSINA60NEW";//_1.5_4.2
    else if ( fitType.Contains("CB2",TString::kIgnoreCase) ) subResultName = "PSICB2";//_2.2_3.9
    else {
      AliError("I don't know from which MC subresult take the tails");
      return kFALSE;
    }

    TObjArray* simArr = new TObjArray;
    if ( SIM() && !fitType.Contains("mctails2",TString::kIgnoreCase)){
      printf("Using first mc file ...\n");
      simArr->Add(SIM());
    }
    else if ( SIM2() && fitType.Contains("mctails2",TString::kIgnoreCase)){
      printf("Using second mc file \n");
      simArr->Add(SIM2());
    }
    else if ( SIM2() && fitType.Contains("INDEPTAILS",TString::kIgnoreCase)   )
      simArr->Add(SIM2()); // If we have set the key to get the JPsi ans PsiP tails
    else {
      AliError("Don't add any SIM file for this test, see inside for details\n");
      return kFALSE;
    }

    TIter nextSim(simArr);
    AliAnalysisMuMu* currentSIM;

    TString spath(pathCentrPairCut);

    spath.Prepend(Form("/%s",Config()->First(Config()->DimuonTriggerKey(),kTRUE).Data()));//FIXME: Care with this when there is more than one selection in the list
    spath.Prepend(Form("/%s",Config()->First(Config()->EventSelectionKey(),kTRUE).Data()));

    Bool_t okMCtails = kFALSE;
    while ( (currentSIM = static_cast<AliAnalysisMuMu*>(nextSim())) )
    {
      TString sspectraName(spectraName);
      if (sspectraName.EndsWith("-AccEffCorr"))
      {
        sspectraName.ReplaceAll("-AccEffCorr","");
        sspectraName.Remove(sspectraName.Length());
      }
      AliAnalysisMuMuSpectra* minvMCSpectra = 0x0;
      minvMCSpectra = currentSIM->SPECTRA(Form("/FitResults%s/%s",spath.Data(),sspectraName.Data()));
      if (!minvMCSpectra){
       printf("######## LOOKING FOR MOHAMAD4S TAIL ########\n");
       minvMCSpectra = currentSIM->SPECTRA(Form("/FitResults/ALL/ANY/PP/pALLPAIRYPAIRPTIN0.0-12.0RABSMATCHLOWETA/%s",sspectraName.Data()));
      }
      if (!minvMCSpectra)
      {
        AliError(Form("Could not find spectra /FitResults%s/%s for associated simulation",spath.Data(),sspectraName.Data()));
        currentSIM->OC()->Print("*:Ali*");
        fitType = "";
        continue;
      }

      AliAnalysisMuMuJpsiResult* minvMCResult = static_cast<AliAnalysisMuMuJpsiResult*>(minvMCSpectra->GetResultForBin(*bin));
      if ( !minvMCResult )
      {
        AliError(Form("Cannot get MC tails cause the minv result for bin %s in %s/%s was not found",bin->AsString().Data(),spath.Data(),sspectraName.Data()));
        fitType = "";
        continue;
      }

      AliAnalysisMuMuJpsiResult* r = dynamic_cast<AliAnalysisMuMuJpsiResult*>(minvMCResult->SubResult(subResultName.Data())); //FIXME: Find an independet way of naming results
      if  ( r && subResultName.Contains("PSICB2") )
      {
        fitType += Form(":al%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("al%s",currentSIM->GetParticleName())));
        fitType += Form(":nl%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("nl%s",currentSIM->GetParticleName())));
        fitType += Form(":au%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("au%s",currentSIM->GetParticleName())));
        fitType += Form(":nu%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("nu%s",currentSIM->GetParticleName())));

        std::cout << " Using MC " << currentSIM->GetParticleName() << " CB2 tails... " << std::endl;
        std::cout << std::endl;
        okMCtails =kTRUE;
      }
      else if ( r && subResultName.Contains("PSINA60NEW") )
      {
        fitType += Form(":p1L%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p1L%s",currentSIM->GetParticleName())));
        fitType += Form(":p2L%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p2L%s",currentSIM->GetParticleName())));
        fitType += Form(":p3L%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p3L%s",currentSIM->GetParticleName())));
        fitType += Form(":p1R%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p1R%s",currentSIM->GetParticleName())));
        fitType += Form(":p2R%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p2R%s",currentSIM->GetParticleName())));
        fitType += Form(":p3R%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("p3R%s",currentSIM->GetParticleName())));

        fitType += Form(":aL%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("aL%s",currentSIM->GetParticleName())));
        fitType += Form(":aR%s=%f",currentSIM->GetParticleName(),r->GetValue(Form("aR%s",currentSIM->GetParticleName())));

        std::cout << " Using MC " << currentSIM->GetParticleName() << " NA60New tails... " << std::endl;
        std::cout << std::endl;
        okMCtails =kTRUE;
      }
      else
      {
        AliError(Form("Cannot get MC tails. MC Subresult %s not found",minvMCResult->GetName()));
        fitType = "";
      }
    }
    delete simArr;
    return okMCtails;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::GetParametersFromResult(TString& fitType, AliAnalysisMuMuJpsiResult* minvResult) const
{
/**
 * @brief   Get parameters from a result
 * @details intended to be used for mean pt fits where we need the signal and backgroud parameters
 *
 * @param fitType    'NA60NEW','CB2','INDEPTAILS' implemented at the moment
 * @param minvResult [description]
 */
    AliWarning("Re-implement me !!!"); //FIXME: The parameters to get will depend on the fit function and also in this way is not suitable for other particles (ie Upsilon)(Find a way to get the particle(s) name)

    TString msg("");
    if ( minvResult )
    {
      // We take the usual parameters (the ones from JPsi and the normalization of the Psi')
      fitType += minvResult->HasValue("Weight") ? Form(":Weight=%f",minvResult->GetValue("Weight")) : Form(":Weight=%f",1.); //FIXME: Names are not correct
      fitType += minvResult->HasValue("kJPsi") ? Form(":kJPsi=%f",minvResult->GetValue("kJPsi")) : Form(":kJPsi=%f",1.);
      fitType += minvResult->HasValue("mJPsi") ? Form(":mJPsi=%f",minvResult->GetValue("mJPsi")) : Form(":mJPsi=%f",1.);
      fitType += minvResult->HasValue("sJPsi") ? Form(":sJPsi=%f",minvResult->GetValue("sJPsi")) : Form(":sJPsi=%f",1.);

      fitType += minvResult->HasValue("NofJPsi") ? Form(":NofJPsi=%f",minvResult->GetValue("NofJPsi")) :  Form(":NofJPsi=%f",1.) ;
      fitType += minvResult->HasValue("NofJPsi") ? Form(":ErrStatNofJPsi=%f",minvResult->GetErrorStat("NofJPsi")) : Form(":ErrStatNofJPsi=%f",1.) ;

      fitType += minvResult->HasValue("kPsiP") ? Form(":kPsiP=%f",minvResult->GetValue("kPsiP")) :  Form(":kPsiP=%f",1.) ;

      TString minvName(minvResult->GetName());

      TString minvRangeParam = minvName;
      minvRangeParam.Remove(0,minvRangeParam.First("_") + 1);
      fitType += Form(":MinvRS=%s",minvRangeParam.Data());

      fitType += minvResult->HasValue("FSigmaPsiP") ? Form(":FSigmaPsiP=%f",minvResult->GetValue("FSigmaPsiP")) :  Form(":FSigmaPsiP=%f",1.) ;

      if ( fitType.Contains("CB2",TString::kIgnoreCase) )
      {
        fitType += minvResult->HasValue("alJPsi") ? Form(":alJPsi=%f",minvResult->GetValue("alJPsi")) :  Form(":alJPsi=%f",1.) ;
        fitType += minvResult->HasValue("nlJPsi") ? Form(":nlJPsi=%f",minvResult->GetValue("nlJPsi")) :  Form(":nlJPsi=%f",1.) ;
        fitType += minvResult->HasValue("auJPsi") ? Form(":auJPsi=%f",minvResult->GetValue("auJPsi")) :  Form(":auJPsi=%f",1.) ;
        fitType += minvResult->HasValue("nuJPsi") ? Form(":nuJPsi=%f",minvResult->GetValue("nuJPsi")) :  Form(":nuJPsi=%f",1.) ;

        msg += "JPsi CB2 signal parameters";
        //    fitType += minvResult->HasValue("NofPsiP") ? Form(":NofPsiP=%f",minvResult->GetValue("NofPsiP")) :  ? Form(":NofPsiP=%f",1.) ;
        //    fitType += minvResult->HasValue() ? Form(":ErrStatNofPsiP=%f",minvResult->GetErrorStat("NofPsiP"));

        if ( fitType.Contains("INDEPTAILS") )
        {
          //        minvName = minvResult->GetName();
          if ( minvName.Contains("INDEPTAILS") )
          {
            // In case we use independent parameters tails for JPsi and Psi' we take also the Psi' ones
            fitType += minvResult->HasValue("alPsiP") ? Form(":alPsiP=%f",minvResult->GetValue("alPsiP")) :  Form(":alPsiP=%f",1.) ;
            fitType += minvResult->HasValue("nlPsiP") ? Form(":nlPsiP=%f",minvResult->GetValue("nlPsiP")) :  Form(":nlPsiP=%f",1.) ;
            fitType += minvResult->HasValue("auPsiP") ? Form(":auPsiP=%f",minvResult->GetValue("auPsiP")) :  Form(":auPsiP=%f",1.) ;
            fitType += minvResult->HasValue("nuPsiP") ? Form(":nuPsiP=%f",minvResult->GetValue("nuPsiP")) :  Form(":nuPsiP=%f",1.) ;
            fitType += minvResult->HasValue("mPsiP") ? Form(":mPsiP=%f",minvResult->GetValue("mPsiP")) :  Form(":mPsiP=%f",1.) ;
            fitType += minvResult->HasValue("sPsiP") ? Form(":sPsiP=%f",minvResult->GetValue("sPsiP")) :  Form(":sPsiP=%f",1.) ;

            msg += " + PsiP CB2 signal parameters";
          }
          else
          {
            AliError(Form("Cannot get PsiP tails from result. Result %s does not contain PsiP tails info => Fit will fail ",minvResult->GetName()));
            fitType = "";
            return;
          }
        }
      }
      else if ( fitType.Contains("NA60NEW",TString::kIgnoreCase) )
      {
        fitType += minvResult->HasValue("p1LJPsi") ? Form(":p1LJPsi=%f",minvResult->GetValue("p1LJPsi")) :  Form(":p1LJPsi=%f",1.) ;
        fitType += minvResult->HasValue("p2LJPsi") ? Form(":p2LJPsi=%f",minvResult->GetValue("p2LJPsi")) :  Form(":p2LJPsi=%f",1.) ;
        fitType += minvResult->HasValue("p3LJPsi") ? Form(":p3LJPsi=%f",minvResult->GetValue("p3LJPsi")) :  Form(":p3LJPsi=%f",1.) ;
        fitType += minvResult->HasValue("p1RJPsi") ? Form(":p1RJPsi=%f",minvResult->GetValue("p1RJPsi")) :  Form(":p1RJPsi=%f",1.) ;
        fitType += minvResult->HasValue("p2RJPsi") ? Form(":p2RJPsi=%f",minvResult->GetValue("p2RJPsi")) :  Form(":p2RJPsi=%f",1.) ;
        fitType += minvResult->HasValue("p3RJPsi") ? Form(":p3RJPsi=%f",minvResult->GetValue("p3RJPsi")) :  Form(":p3RJPsi=%f",1.) ;

        fitType += minvResult->HasValue("aLJPsi") ? Form(":aLJPsi=%f",minvResult->GetValue("aLJPsi")) :  Form(":aLJPsi=%f",1.) ;
        fitType += minvResult->HasValue("aRJPsi") ? Form(":aRJPsi=%f",minvResult->GetValue("aRJPsi")) :  Form(":aRJPsi=%f",1.) ;

        msg += "JPsi NA60New signal parameters";

        if ( fitType.Contains("INDEPTAILS") )
        {
          //        TString minvName(minvResult->GetName());
          if ( minvName.Contains("INDEPTAILS") )
          {
            // In case we use independent parameters tails for JPsi and Psi' we take also the Psi' ones
            fitType += minvResult->HasValue("p1LPsiP") ? Form(":p1LPsiP=%f",minvResult->GetValue("p1LPsiP")) :  Form(":p1LPsiP=%f",1.) ;
            fitType += minvResult->HasValue("p2LPsiP") ? Form(":p2LPsiP=%f",minvResult->GetValue("p2LPsiP")) :  Form(":p2LPsiP=%f",1.) ;
            fitType += minvResult->HasValue("p3LPsiP") ? Form(":p3LPsiP=%f",minvResult->GetValue("p3LPsiP")) :  Form(":p3LPsiP=%f",1.) ;
            fitType += minvResult->HasValue("p1RPsiP") ? Form(":p1RPsiP=%f",minvResult->GetValue("p1RPsiP")) :  Form(":p1RPsiP=%f",1.) ;
            fitType += minvResult->HasValue("p2RPsiP") ? Form(":p2RPsiP=%f",minvResult->GetValue("p2RPsiP")) :  Form(":p2RPsiP=%f",1.) ;
            fitType += minvResult->HasValue("p3RPsiP") ? Form(":p3RPsiP=%f",minvResult->GetValue("p3RPsiP")) :  Form(":p3RPsiP=%f",1.) ;

            fitType += minvResult->HasValue("aLPsiP") ? Form(":aLPsiP=%f",minvResult->GetValue("aLPsiP")) :  Form(":aLPsiP=%f",1.) ;
            fitType += minvResult->HasValue("aRPsiP") ? Form(":aRPsiP=%f",minvResult->GetValue("aRPsiP")) :  Form(":aRPsiP=%f",1.) ;

            msg += " + PsiP NA60New signal parameters";

          }
          else
          {
            AliError(Form("Cannot get PsiP tails from result. Result %s does not contain PsiP tails info => Fit will fail ",minvResult->GetName()));
            fitType = "";
            return;
          }
        }
      }
      else
      {
        AliError(Form("Cannot get the parameters from %s",minvResult->GetName()));
        fitType = "";
        return;
      }
      // Now we take the background parameters
      if ( fitType.Contains("VWG_") || fitType.Contains("VWGINDEPTAILS") ) //FIXME: Check that cannot be misunderstood(like Exp x Pol2..). In fact it can be misunderstood since the meanpt function name has also the name of the function to fit the bkg (free parameters). Also add the rest of the BKG functions
      {
        fitType += minvResult->HasValue("kVWG") ? Form(":kVWG=%f",minvResult->GetValue("kVWG")) :  Form(":kVWG=%f",1.) ;
        fitType += minvResult->HasValue("mVWG") ? Form(":mVWG=%f",minvResult->GetValue("mVWG")) :  Form(":mVWG=%f",1.) ;
        fitType += minvResult->HasValue("sVWG1") ? Form(":sVWG1=%f",minvResult->GetValue("sVWG1")) :  Form(":sVWG1=%f",1.) ;
        fitType += minvResult->HasValue("sVWG2") ? Form(":sVWG2=%f",minvResult->GetValue("sVWG2")) :  Form(":sVWG2=%f",1.) ;

        msg += " + VWG Bkg parameters";
      }
      else if ( fitType.Contains("POL2EXP_") || fitType.Contains("POL2EXPINDEPTAILS") )
      {
        fitType += minvResult->HasValue("kPol2Exp") ? Form(":kPol2Exp=%f",minvResult->GetValue("kPol2Exp")) :  Form(":kPol2Exp=%f",1.) ;
        fitType += minvResult->HasValue("pol0") ? Form(":pol0=%f",minvResult->GetValue("pol0")) :  Form(":pol0=%f",1.) ;
        fitType += minvResult->HasValue("pol1") ? Form(":pol1=%f",minvResult->GetValue("pol1")) :  Form(":pol1=%f",1.) ;
        fitType += minvResult->HasValue("pol2") ? Form(":pol2=%f",minvResult->GetValue("pol2")) :  Form(":pol2=%f",1.) ;
        fitType += minvResult->HasValue("exp") ? Form(":exp=%f",minvResult->GetValue("exp")) :  Form(":exp=%f",1.) ;

        msg += " + Pol2xExp Bkg parameters";
      }
      else if ( fitType.Contains("POL4EXP_") || fitType.Contains("POL4EXPINDEPTAILS") )
      {
        fitType += minvResult->HasValue("kPol4Exp") ? Form(":kPol4Exp=%f",minvResult->GetValue("kPol4Exp")) :  Form(":kPol4Exp=%f",1.) ;
        fitType += minvResult->HasValue("pol0") ? Form(":pol0=%f",minvResult->GetValue("pol0")) :  Form(":pol0=%f",1.) ;
        fitType += minvResult->HasValue("pol1") ? Form(":pol1=%f",minvResult->GetValue("pol1")) :  Form(":pol1=%f",1.) ;
        fitType += minvResult->HasValue("pol2") ? Form(":pol2=%f",minvResult->GetValue("pol2")) :  Form(":pol2=%f",1.) ;
        fitType += minvResult->HasValue("pol3") ? Form(":pol3=%f",minvResult->GetValue("pol3")) :  Form(":pol3=%f",1.) ;
        fitType += minvResult->HasValue("pol4") ? Form(":pol4=%f",minvResult->GetValue("pol4")) :  Form(":pol4=%f",1.) ;
        fitType += minvResult->HasValue("exp") ? Form(":exp=%f",minvResult->GetValue("exp")) :  Form(":exp=%f",1.) ;

        msg += " + Pol4xExp Bkg parameters";
      }
      else if ( fitType.Contains("POL1POL2_") )
      {
        fitType += minvResult->HasValue("a") ? Form(":a=%f",minvResult->GetValue("a")) :  Form(":a=%f",1.) ;
        fitType += minvResult->HasValue("b") ? Form(":b=%f",minvResult->GetValue("b")) :  Form(":b=%f",1.) ;
        fitType += minvResult->HasValue("a'") ? Form(":a'=%f",minvResult->GetValue("a'")) :  Form(":a'=%f",1.) ;
        fitType += minvResult->HasValue("b'") ? Form(":b'=%f",minvResult->GetValue("b'")) :  Form(":b'=%f",1.) ;
        fitType += minvResult->HasValue("c'") ? Form(":c'=%f",minvResult->GetValue("c'")) :  Form(":c'=%f",1.) ;

        msg += " + pol1/pl2 Bkg parameters";
      }

      std::cout << "Using " << msg.Data() << " from " << minvResult->GetName() <<  " inv mass result" << std::endl;
      std::cout << "" << std::endl;
    }
    else
    {
      AliError(Form("Cannot get tails from result. Result %s not found",minvResult->GetName()));
      fitType = "";
      return;
    }
}

//_____________________________________________________________________________
ULong64_t AliAnalysisMuMu::GetTriggerScalerCount(const char* triggerList, Int_t runNumber)
{
/**
 * @brief Get the expected (from OCDB scalers) trigger count
 *
 * @param triggerList Tokenize by ','
 * @param runNumber   [description]
 *
 * @return Counts
 */
    AliAnalysisTriggerScalers ts(runNumber,Config()->OCDBPath());

    TObjArray* triggers = TString(triggerList).Tokenize(",");
    TObjString* trigger;
    TIter next(triggers);
    ULong64_t n(0);

    while ( ( trigger = static_cast<TObjString*>(next()) ) )
        {
        AliAnalysisTriggerScalerItem* item = ts.GetTriggerScaler(runNumber,"L2A",trigger->String().Data());
        if (item)
            {
            n += item->Value();
            }
        delete item;
        }
    delete triggers;

    return n;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::GetSpectraFromConfig(const char* binType, const char* flavour) const
{
/**
 * @brief Pointer to given spectra
 * @details The spectra path is selected according to first keys in the config. file. Spectra name should be 'PSI-<binType>'
 *
 * @param binType    [description]
 * @param flavour    [description]
 *
 * @return Pointer to the spectra
 */
    TString sbinType(binType);
    TString sflavour(flavour);
    sbinType.ToUpper();
    sflavour.ToUpper();

    TString spectraName(Form("/FitResults/%s/%s/%s/%s/PSI-%s",
                             Config()->First(Config()->EventSelectionKey(),IsSimulation()).Data(),
                             Config()->First(Config()->DimuonTriggerKey(),IsSimulation()).Data(),
                             Config()->First(Config()->CentralitySelectionKey(),IsSimulation()).Data(),
                             Config()->First(Config()->PairSelectionKey(),IsSimulation()).Data(),
                             sbinType.Data()));

    cout << "spectraName : " << spectraName.Data() << endl;

    if (sflavour.Length()>0){
      spectraName += "-";
      spectraName += sflavour.Data();
    }

    return SPECTRA(spectraName.Data());
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::GetSpectra(const char* binType , const char* EventSelection,  const char* DimuonTrigger,
  const char* Centrality, const char* PairSelectionKey, const char* flavour) const
{
/**
 * @brief Pointer to a given spectra
 * @details More specific method
 *
 * @param binType [description]
 * @param flavour [description]
 *
 * @return AliAnalysisMuMuSpectra
 */
    TString sbinType(binType);
    TString sflavour(flavour);
    sbinType.ToUpper();
    sflavour.ToUpper();

    TString spectraName(Form("/FitResults/%s/%s/%s/%s/PSI-%s",
                             EventSelection,
                             DimuonTrigger,
                             Centrality,
                             PairSelectionKey,
                             sbinType.Data()));

    cout << "spectraName MC : " << spectraName.Data() << endl;

    if (sflavour.Length()>0)
    {
      spectraName += "-";
      spectraName += sflavour.Data();
    }

    return SPECTRA(spectraName.Data());
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::MergeSpectra(const char* spectraPath1,const char* spectraPath2,const char* spectraname)
{


  if (!OC() || !CC())
  {
    AliError("No mergeable/counter collection. Consider Upgrade()");
    return 0x0;
  }
  else
  {
    cout <<      " ================================================================ " << endl;
    cout <<      "                             MergeSpectra RAA                             " << endl;
    cout <<      " ================================================================ " << endl;
  }

  TString spectraName1(Form("/%s/%s",spectraPath1,spectraname));
  TString spectraName2(Form("/%s/%s",spectraPath2,spectraname));

  cout << "spectraName 1 : " << spectraName1.Data() << endl;
  cout << "spectraName 2 : " << spectraName2.Data() << endl;

  AliAnalysisMuMuSpectra* mspectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraName1.Data())->Clone());
  AliAnalysisMuMuSpectra* spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraName2.Data()));
  if(!mspectra || !spectra){
    AliError("Could not get spectras");
    return 0x0 ;
  }

  mspectra->Merge(spectra);
  mspectra->SetName(Form("%s_merged",spectraname));

  TObject* o = fMergeableCollection->GetObject(spectraPath1,mspectra->GetName());
  if (o) {
    AliWarning(Form("Replacing %s/%s",spectraPath1,mspectra->GetName()));
    fMergeableCollection->Remove(Form("/%s/%s",spectraPath1,mspectra->GetName()));

  }
  Bool_t adoptOK = fMergeableCollection->Adopt(spectraPath1,mspectra);

  if ( adoptOK ) std::cout << "+++Spectra " << mspectra->GetName() << " adopted" << std::endl;
  else AliError(Form("Could not adopt spectra %s",mspectra->GetName()));
  Update();
  return mspectra;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::MergeHistoInCentalityBins(const char* binToMerge, const char* newBin, Bool_t mix, const char* refbin) const
{

  if (!OC() || !CC())
  {
    AliError("No mergeable/counter collection. Consider Upgrade()");
    return;
  }
  else
  {
    cout <<      " ================================================================ " << endl;
    cout <<      "    MergeHistoInCentalityBins "<< binToMerge << " in "<< newBin << endl;
    cout <<      "    Warning : spectra must have been deleted first (CleanAllSpectra()) 'cause not mergeable  " << endl;
    cout <<      " ================================================================ " << endl;
  }
    // Get configuration settings

    TObjArray* pairCutArray    = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());

    TObjArray* triggerArray(0x0);
    if (!mix) triggerArray = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
    else                              triggerArray = Config()->GetListElements(Config()->MixTriggerKey(),IsSimulation());

    TObjArray* eventTypeArray(0x0);
    if (!mix) eventTypeArray = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
    else                              eventTypeArray = Config()->GetListElements(Config()->EventSelectionMixKey(),IsSimulation());

    TString refTrigger(Form("%s",Config()->First(Config()->RefMixTriggerKey(),IsSimulation()).Data()));
    TString refEvent(Form("%s",Config()->First(Config()->RefMixEventSelectionKey(),IsSimulation()).Data()));



    TObjArray* bins;

    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextEventType(eventTypeArray);
    TIter nextPairCut(pairCutArray);

    // Strings
    TObjString* strigger;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* sparticle;
    TObjString* scentrality;

    nextEventType.Reset();
    // Loop on each envenType (see MuMuConfig)
    //==============================================================================
    while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
    {
      AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
      nextTrigger.Reset();
      // Loop on each trigger (see MuMuConfig)
      //==============================================================================
      while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
      {
      AliDebug(1,Form("-TRIGGER %s",strigger->String().Data()));
          TString identifier = Form("/%s/%s",seventType->String().Data(),strigger->String().Data());

        nextPairCut.Reset();
        // Loop on each paircut (not the ones in MuMuConfig but the ones set)
        //==============================================================================
        while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
        {
          AliDebug(1,Form("---PAIRCUT %s",spairCut->String().Data()));

          TList* list(0x0);
          if(!mix)     list = OC()->CreateListOfObjectNames(Form("%s/%s/%s/",identifier.Data(),refbin,spairCut->String().Data()));
          else if (mix)list = OC()->CreateListOfObjectNames(Form("/MIX/%s_%s%s/%s/%s/",refEvent.Data(),refTrigger.Data(),identifier.Data(),refbin,spairCut->String().Data()));

          // printf("%s/%s/%s/",identifier.Data(),refbin,spairCut->String().Data());
          // printf("list : %p\n",list );
          // list->Print("");
          if(list->IsEmpty()) continue;

          TObject* p =0x0;

          for (int j = 0; j < list->GetEntries(); ++j) {

                const char* objName =list->At(j)->GetName();

                TObject * o =0x0;
                if(!mix) o = OC()->GetSum(Form("%s/%s/%s/%s",identifier.Data(),binToMerge,spairCut->String().Data(),objName));
                else     o = OC()->GetSum(Form("/MIX/%s_%s%s/%s/%s/%s",refEvent.Data(),refTrigger.Data(),identifier.Data(),binToMerge,spairCut->String().Data(),objName));

                if ( o ) {
                   if(!mix) p = OC()->GetObject(Form("%s/%s/%s",identifier.Data(),newBin,spairCut->String().Data()),o->GetName());
                   else     p = OC()->GetObject(Form("/MIX/%s_%s%s/%s/%s",refEvent.Data(),refTrigger.Data(),identifier.Data(),newBin,spairCut->String().Data()),o->GetName());

                    if (p && !mix) OC()->Remove(Form("%s/%s/%s/%s",identifier.Data(),newBin,spairCut->String().Data(),o->GetName()));
                    if (p &&  mix) OC()->Remove(Form("/MIX/%s_%s%s/%s/%s/%s",refEvent.Data(),refTrigger.Data(),identifier.Data(),newBin,spairCut->String().Data(),objName));

                    Bool_t adoptOK = kFALSE;
                    if(!mix) adoptOK = OC()->Adopt(Form("%s/%s/%s",identifier.Data(),newBin,spairCut->String().Data()),o->Clone());
                    else     adoptOK = OC()->Adopt(Form("/MIX/%s_%s%s/%s/%s",refEvent.Data(),refTrigger.Data(),identifier.Data(),newBin,spairCut->String().Data()),o->Clone());

                    if ( adoptOK && !mix) printf(" --- Merge done for %s/%s/%s/%s \n",identifier.Data(),newBin,spairCut->String().Data(),o->GetName());
                    else if ( adoptOK &&  mix) printf(" --- Merge done for /MIX/%s_%s%s/%s/%s/%s \n",refEvent.Data(),refTrigger.Data(),identifier.Data(),newBin,spairCut->String().Data(),o->GetName());
                    else printf(" --- A problem occurs \n");
                } else {
                  printf("Cannot get object\n");
                  continue;
                }
            }

            delete list;
        }
    }
  }


  delete eventTypeArray ;
  delete triggerArray ;
  delete pairCutArray ;

  return ;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMu::PlotAccEfficiency(const char* SpectraName, const char* subresultname )
{
/**
 * @brief Plot AccxEff
 * @details Search for an AliAnalysisMuMuSpectra according to first entry of each keys in the config. file. FIXME : make it general.
 *
 * @param SpectraName [description]
 * @return TH1 to be handle by the user
 */
    if ( !IsSimulation() ){
      AliError("Could not get AccxEff histo: Not simulation file");
      return 0x0;
    }

    TString path(Form("/FitResults/%s/%s/%s/%s",
                      Config()->First(Config()->EventSelectionKey(),kTRUE).Data(),
                      Config()->First(Config()->DimuonTriggerKey(),kTRUE).Data(),
                      Config()->First(Config()->CentralitySelectionKey(),kTRUE).Data(),
                      Config()->First(Config()->PairSelectionKey(),kTRUE).Data()));

    AliAnalysisMuMuSpectra* s = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(Form("%s/%s",path.Data(),SpectraName)));
    if ( !s )
        {
        AliError(Form("No AccEff spectra %s found in %s",SpectraName,path.Data()));
        return 0x0;
        }

    if(!TString(SpectraName).Contains("VS"))return s->Plot(Form("AccEff%s",GetParticleName()),subresultname,kFALSE);//_2.2_3.9
    else return static_cast<TH2*>(s->Plot(Form("AccEff%s",GetParticleName()),subresultname,kFALSE))->ProjectionX();//_2.2_3.9

}

//_____________________________________________________________________________
UInt_t AliAnalysisMuMu::GetSum(AliCounterCollection& cc, const char* triggerList,
                               const char* eventSelection, Int_t runNumber)
{
  TObjArray* ktrigger = cc.GetKeyWords("trigger").Tokenize(",");
  TObjArray* kevent = cc.GetKeyWords("event").Tokenize(",");
  TObjArray* a = TString(triggerList).Tokenize(" ");
  TIter next(a);
  TObjString* str;
  UInt_t n(0);

  TString sEventSelection(eventSelection);
  sEventSelection.ToUpper();

  if ( kevent->FindObject(sEventSelection.Data()) )
  {
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      if ( ktrigger->FindObject(str->String().Data()) )
      {
        if ( runNumber < 0 )
        {
          n +=  static_cast<UInt_t>(cc.GetSum(Form("trigger:%s/event:%s",str->String().Data(),eventSelection)));
        }
        else
        {
          n +=  static_cast<UInt_t>(cc.GetSum(Form("trigger:%s/event:%s/run:%d",str->String().Data(),eventSelection,runNumber)));
        }
      }
    }
  }

  delete a;
  delete ktrigger;
  delete kevent;
  return n;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::GetCollectionsFromAnySubdir(TDirectory& dir,
                                                    AliMergeableCollection*& oc,
                                                    AliCounterCollection*& cc,
                                                    AliAnalysisMuMuBinning*& bin)
{
/**
 * @brief Find, within dir and its sub-directories, the objects OC,CC and BIN
 *
 * @param dir [description]
 * @param oc [description]
 * @param cc [description]
 * @param bin [description]
 */
  TList* keys = dir.GetListOfKeys();
  TIter next(keys);

  TKey* k;

  while ( ( k = static_cast<TKey*>(next()) ) )
  {
    TObject* object = k->ReadObj();

    if ( object->InheritsFrom("TDirectory") )
    {
      TDirectory* d = static_cast<TDirectory*>(object);
      GetCollectionsFromAnySubdir(*d,oc,cc,bin);
      continue;
    }

    if ( ( object->InheritsFrom("AliMergeableCollection") ) &&
        ( strcmp(object->GetName(),"OC")==0 ) )
    {
      oc = dynamic_cast<AliMergeableCollection*>(object);
      fDirectory = dir.GetName();
    }

    if ( ( object->InheritsFrom("AliCounterCollection") ) &&
        ( strcmp(object->GetName(),"CC")==0 ) )
    {
      cc = dynamic_cast<AliCounterCollection*>(object);
      fDirectory = dir.GetName();
    }

    if ( ( object->InheritsFrom("AliAnalysisMuMuBinning") ) &&
        ( strncmp(object->GetName(),"BIN",3)==0 ) )
    {
      bin = dynamic_cast<AliAnalysisMuMuBinning*>(object);
      fDirectory = dir.GetName();
    }

  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::GetCollections(const char* rootfile,const char* subdir,AliMergeableCollection*& oc,AliCounterCollection*& cc,AliAnalysisMuMuBinning*& bin,std::set<int>& runnumbers)
{
/**
 * @brief Get access to the mergeable collection, counter collection and binning within file rootfile.
 * @details   rootfile is a filename, with an optional directory (with the syntax filename.root:directory) where the collections are to be found.
 *
 * The rootfile should be the outfile of the AliAnalysisTaskMuMu class.
 *
 */
  oc = 0x0;
  cc = 0x0;
  bin = 0x0;


  TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(rootfile));
  if (!f)
  {
    f = TFile::Open(rootfile);
  }

  if ( !f || f->IsZombie() )
  {
    return kFALSE;
  }

  TString dir = subdir;

  if (dir.Length())
  {
    dir.Remove(TString::kBoth,'/');
    dir += "/";
  }

  f->GetObject(Form("%sOC",dir.Data()),oc);
  if (!oc)
  {
    f->GetObject(Form("%sMC",dir.Data()),oc);
  }
  f->GetObject(Form("%sCC",dir.Data()),cc);

  TIter next(f->GetListOfKeys());
  TKey* key;

  while ( ( key = static_cast<TKey*>(next())) && !bin )
  {
    if ( strcmp(key->GetClassName(),"AliAnalysisMuMuBinning")==0 )
    {
      bin = dynamic_cast<AliAnalysisMuMuBinning*>(key->ReadObj());
    }
  }

  if ( (!oc || !cc) && fDirectory.Length()==0 )
  {
    // one more try, searching in subdirectories as well
    GetCollectionsFromAnySubdir(*f,oc,cc,bin);
  }

  if (!oc || !cc)
  {
    AliError("Could not get OC, CC and BIN. Is that an old file ? Try to upgrade it or check it's the right file...");
    return kFALSE;
  }

  // get run list
  TObjArray* runs = cc->GetKeyWords("run").Tokenize(",");
  runs->Sort();
  TIter nextRun(runs);
  TObjString* srun;

  runnumbers.clear();

  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    runnumbers.insert(srun->String().Atoi());
  }

  delete runs;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::IsSimulation() const
{
/**
 * @brief whether or not we have MC information
 */
    if (!fMergeableCollection) return kFALSE;

    TList* list = fMergeableCollection->CreateListOfKeys(0);
    TIter next(list);
    TObjString* str;
    Bool_t ok(kFALSE);

    while ( ( str = static_cast<TObjString*>(next()) ) ){
      if ( str->String().Contains(AliAnalysisMuMuBase::MCInputPrefix()) ) ok = kTRUE;
    }
    delete list;

    return ok;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::NormMixedMinv(const char* binType, const char* particle, const char* flavour, Bool_t corrected, Double_t Mmin, Double_t Mmax)
{
/**
 * @brief Create normalized mixed unlike-sign minv spectrum.
 * @details The normalization is made via : int{N_mix_+-} = \int{2*R*\sqrt{N_raw_++ * N_raw_--}}. Raw histograms are selected according to RefMix...Key() in the config. file.
 *
 * Mixed Histograms are selected according to Mix...Key() in config. file. See AliAnalysisMuMuConfig class.
 *
 * @param binType    [description]
 * @param particle   particle name
 * @param flavour    [description]
 * @param corrected  If Minv are already AccxEff corrected or not
 */
    if(!OC())
    {
        AliError("No mergeable. Consider Upgrade()");
        return;
    }
    else
    {
        cout <<      " ================================================================ " << endl;
        cout <<      "                       NormMixedMinv                  " << endl;
        cout <<      " ================================================================ " << endl;
    }

    // Get configuration settings
    TObjArray* eventTypeArray    = Config()->GetListElements(Config()->RefMixEventSelectionKey(),IsSimulation());
    TObjArray* triggerArray      = Config()->GetListElements(Config()->RefMixTriggerKey(),IsSimulation());
    TObjArray* eventTypeMixArray = Config()->GetListElements(Config()->EventSelectionMixKey(),IsSimulation());
    TObjArray* triggerMixArray   = Config()->GetListElements(Config()->MixTriggerKey(),IsSimulation());
    TObjArray* pairCutArray      = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
    TObjArray* centralityArray   = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
    TObjArray* binTypeArray      = TString(binType).Tokenize(",");

    TObjArray* bins;

    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextTriggerMix(triggerMixArray);
    TIter nextEventType(eventTypeArray);
    TIter nextEventTypeMix(eventTypeMixArray);
    TIter nextPairCut(pairCutArray);
    TIter nextbinType(binTypeArray);
    TIter nextCentrality(centralityArray);

    // Strings
    TObjString* strigger;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* sbinType;
    TObjString* scentrality;
    TObjString* sTriggerMix;
    TObjString* seventTypeMix;

    // For the loop comming
    TString signFlagMinv[3] ={"","PlusPlus","MinusMinus"};
    TString signFlagDist[3] ={"","PP","MM"};

    THnSparse* n        =0x0;
    TObject* o          =0x0;

    // Loop on each envenTypeMix (see MuMuConfig)
    while ( ( seventTypeMix = static_cast<TObjString*>(nextEventTypeMix())) ){
      AliDebug(1,Form("EVENTTYPEMIX %s",seventTypeMix->String().Data()));
      nextTriggerMix.Reset();

      // Loop on each triggerMix (see MuMuConfig)
      while ( ( sTriggerMix = static_cast<TObjString*>(nextTriggerMix())) ){
        AliDebug(1,Form("-TRIGGERMIX %s",sTriggerMix->String().Data()));
        nextEventType.Reset();

        // Loop on each envenType (see MuMuConfig)
        while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
        {
          AliDebug(1,Form("--REFEVENTTYPE %s",seventType->String().Data()));
          nextTrigger.Reset();

          // Loop on each trigger (see MuMuConfig)
          while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
          {
            AliDebug(1,Form("---REFTRIGGER %s",strigger->String().Data()));
            nextPairCut.Reset();

            // Loop on each paircut (see MuMuConfig)
            while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
            {
              AliDebug(1,Form("----PAIRCUT %s",spairCut->String().Data()));
              nextbinType.Reset();

              // Loop on each type (see MuMuConfig)
              while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
              {
                AliDebug(1,Form("-----TYPE %s",sbinType->String().Data()));
                nextCentrality.Reset();

                // Loop on each centrality (see MuMuConfig)
                while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
                {
                  AliDebug(1,Form("------CENTRALITY %s",scentrality->String().Data()));
                  nextbinType.Reset();

                  // Loop on each centrality (see MuMuConfig)
                  while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
                  {
                    AliDebug(1,Form("-------BIN %s",sbinType->String().Data()));

                    // Get Binning
                    AliAnalysisMuMuBinning* binning(0x0);
                    if ( fBinning && sbinType->String().Length() > 0 ) binning = fBinning->Project(particle,sbinType->String().Data(),flavour);
                    else  {
                      binning = new AliAnalysisMuMuBinning;
                      binning->AddBin(particle,sbinType->String().Data());
                    }
                    if (!binning) {
                      AliError("oups. binning is NULL");
                      continue;
                    }

                    // Create array
                    TObjArray* bins = binning->CreateBinObjArray(particle);
                    if (!bins){
                      AliError(Form("Did not get any bin for particle %s",particle));
                      delete binning;
                      continue;
                    }

                    // Create ID for the fit which will be used to name results
                    TString idraw(Form("/%s/%s/%s/%s",seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data()));
                    AliDebug(1,Form("idraw = %s\n",idraw.Data() ));
                    TString idMix(Form("/%s/%s/%s/%s",seventTypeMix->String().Data(),sTriggerMix->String().Data(),scentrality->String().Data(),spairCut->String().Data()));
                    AliDebug(1,Form("idMix = %s\n",idMix.Data() ));

                    // The binning pointer, which point at Pt binning, Y binning etc.
                    AliAnalysisMuMuBinning::Range* bin;
                    TIter next(bins);

                    // Add some element to ID
                    TString spectraName(binning->GetName());
                    if ( strcmp(flavour,"") != 0 ){
                      spectraName += "-";
                      spectraName += flavour;
                    }
                    if ( corrected ){
                      spectraName += "-";
                      spectraName += "AccEffCorr";
                    }

                    //MAIN PART : Loop on every binning range
                    while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
                    {
                      // ---- Here we get all the histos we need ----
                      TH1* hTableMinv[6]    = {0x0,  0x0,  0x0,  0x0,  0x0,  0x0};
                      //                      {hMinvPM,  hMinvPP,  hMinvMM,  hMinvMixPM,  hMinvMixPP,  hMinvMixMM};
                      TH1* hTableDistRaw[6] = {0x0,  0x0,  0x0,  0x0,  0x0,  0x0};
                      //                      {hpT,  hpTPP,  hpTMM,  hY,  hYPP,  hYMM};
                      TH1* hTableDistMix[6] = {0x0,  0x0,  0x0,  0x0,  0x0,  0x0};
                      //                      {hpTMix,  hpTMixPP,  hpTMixMM,  hYMix,  hYMixPP,  hYMixMM};
                      TH1* hRCoef;

                      // Get Minv Histo
                      for (int j = 0; j < 3; ++j) {
                        TString hnameraw = corrected ? Form("MinvUS+%s%s_AccEffCorr",bin->AsString().Data(),signFlagMinv[j].Data()) : Form("MinvUS+%s%s",bin->AsString().Data(),signFlagMinv[j].Data());
                        TString hnamemix = corrected ? Form("MinvUS+%s%sMix_AccEffCorr",bin->AsString().Data(),signFlagMinv[j].Data()) : Form("MinvUS+%s%sMix",bin->AsString().Data(),signFlagMinv[j].Data());
                        // Pointer to the histo from histo collection (Yes it is discusting )
                        if ( OC()->Histo(idraw.Data(),hnameraw.Data()) ) hTableMinv[j] = static_cast<TH1*>(OC()->Histo(idraw.Data(),hnameraw.Data())->Clone());
                        else {
                          AliError(Form("Could not find histo %s/%s",idraw.Data(),hnameraw.Data()));
                          continue ;
                        }

                        if ( OC()->Histo(idMix.Data(),hnamemix.Data()) ) hTableMinv[j+3] = static_cast<TH1*>(OC()->Histo(idMix.Data(),hnamemix.Data())->Clone());
                        else {
                          AliError(Form("Could not find histo %s/%s",idMix.Data(),hnamemix.Data()));
                          continue;
                        }
                      }

                      // Get Dist Histo
                      for (int j = 0; j < 3; ++j) {
                        TString hnamePt = corrected ? Form("Pt%s_AccEffCorr",signFlagDist[j].Data()) : Form("Pt%s",signFlagDist[j].Data());
                        TString hnameY  = corrected ? Form("Y%s_AccEffCorr",signFlagDist[j].Data()) : Form("Y%s",signFlagDist[j].Data());

                        TString hnamePtMix = corrected ? Form("PtMix%s_AccEffCorr",signFlagDist[j].Data()) : Form("PtMix%s",signFlagDist[j].Data());
                        TString hnameYMix  = corrected ? Form("YMix%s_AccEffCorr",signFlagDist[j].Data()) : Form("YMix%s",signFlagDist[j].Data());

                        // Pointer to the histo from histo collection (Yes it is discusting )
                        if ( OC()->GetObject(idraw.Data(),hnamePt.Data()) ){
                          // hTableDistRaw[j] = OC()->Histo(idraw.Data(),hnamePt.Data());
                          n = static_cast<THnSparse*>(OC()->GetObject(idraw.Data(),hnamePt.Data()));
                          Int_t binmin = n->GetAxis(1)->FindBin(Mmin);
                          Int_t binmax = n->GetAxis(1)->FindBin(Mmax);
                          n->GetAxis(1)->SetRange(binmin,binmax);
                          hTableDistRaw[j] = n->Projection(0,"e");
                          hTableDistRaw[j]->SetName(Form("%s_%.2f_%.2f",hTableDistRaw[j]->GetName(),Mmin,Mmax)) ;
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idraw.Data(),hnamePt.Data()));
                          continue ;
                        }

                        if ( OC()->GetObject(idraw.Data(),hnameY.Data()) ) {
                          // hTableDistRaw[j+3] = OC()->Histo(idraw.Data(),hnameY.Data());
                          n = static_cast<THnSparse*>(OC()->GetObject(idraw.Data(),hnameY.Data()));
                          Int_t binmin = n->GetAxis(1)->FindBin(Mmin);
                          Int_t binmax = n->GetAxis(1)->FindBin(Mmax);
                          n->GetAxis(1)->SetRange(binmin,binmax);
                          hTableDistRaw[j+3] = n->Projection(0,"e");
                          hTableDistRaw[j+3]->SetName(Form("%s_%.2f_%.2f",hTableDistRaw[j+3]->GetName(),Mmin,Mmax)) ;
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idraw.Data(),hnameY.Data()));
                          continue ;
                        }

                        if ( OC()->GetObject(idMix.Data(),hnamePtMix.Data()) ){
                          // hTableDistMix[j] = OC()->Histo(idMix.Data(),hnamePtMix.Data());
                          n = static_cast<THnSparse*>(OC()->GetObject(idMix.Data(),hnamePtMix.Data()));
                          Int_t binmin = n->GetAxis(1)->FindBin(Mmin);
                          Int_t binmax = n->GetAxis(1)->FindBin(Mmax);
                          n->GetAxis(1)->SetRange(binmin,binmax);
                          hTableDistMix[j] = n->Projection(0,"e");
                          hTableDistMix[j]->SetName(Form("%s_%.2f_%.2f",hTableDistMix[j]->GetName(),Mmin,Mmax)) ;
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idMix.Data(),hnamePtMix.Data()));
                          continue ;
                        }

                        if ( OC()->GetObject(idMix.Data(),hnameYMix.Data()) ) {
                          // hTableDistMix[j+3] = OC()->Histo(idMix.Data(),hnameYMix.Data());
                          n = static_cast<THnSparse*>(OC()->GetObject(idMix.Data(),hnameYMix.Data()));
                          Int_t binmin = n->GetAxis(1)->FindBin(Mmin);
                          Int_t binmax = n->GetAxis(1)->FindBin(Mmax);
                          n->GetAxis(1)->SetRange(binmin,binmax);
                          hTableDistMix[j+3] = n->Projection(0,"e");
                          hTableDistMix[j+3]->SetName(Form("%s_%.2f_%.2f",hTableDistMix[j+3]->GetName(),Mmin,Mmax)) ;
                        } else {
                          AliError(Form("Could not find GetObject %s/%s",idMix.Data(),hnameYMix.Data()));
                          continue ;
                        }
                      }

                      // Check if we have all the histo we need
                      Bool_t missAnHisto = kFALSE;
                      for (int i = 0; i < 6; ++i){
                        if(!hTableMinv[i] || !hTableDistRaw[i] || !hTableDistMix[i]) missAnHisto=kTRUE;
                      }
                      if(missAnHisto){
                        AliError(Form("One histo is missing ... go to next bin\n"));
                        for (int i=0; i<6; i++){ if (hTableMinv[i]!=0x0) delete hTableMinv[i];  }
                        for (int i=0; i<6; i++){ if (hTableDistRaw[i]!=0x0) delete hTableDistRaw[i]; }
                        for (int i=0; i<6; i++){ if (hTableDistMix[i]!=0x0) delete hTableDistMix[i]; }
                        continue;
                      }

                      // ---- First part : compute the R Factor ----
                      hRCoef = static_cast<TH1*>(hTableMinv[3]->Clone());
                      hRCoef->SetName(Form("RCoefficient_%s",bin->AsString().Data()));
                      hRCoef->SetTitle("R Coefficient versus M");
                      hRCoef->GetYaxis()->SetTitle("R");

                      // Copy and Multiply like sign histograms
                      TH1* hMinvMixPPCopy = static_cast<TH1*>(hTableMinv[4]->Clone());
                      hMinvMixPPCopy->Multiply(hTableMinv[5]);

                      // Compute and fill square roots of like-signs histo
                      Int_t nEntries = hMinvMixPPCopy->GetEntries();
                      for (int j = 0; j < nEntries ; j++) {
                        Double_t binContent = hMinvMixPPCopy->GetBinContent(j);
                        Double_t binError   = hMinvMixPPCopy->GetBinError(j);
                        if( binContent !=0 && binError !=0 ){
                          hMinvMixPPCopy->SetBinContent(j,TMath::Sqrt(binContent));
                          hMinvMixPPCopy->SetBinError(j,TMath::Sqrt(binContent)*(binError/binContent));
                        }/* else {
                          hMinvMixPPCopy->SetBinContent(j,0.);
                          hMinvMixPPCopy->SetBinError(j,0.);
                        }*/
                      }

                      // Final R Factor
                      hRCoef->Divide(hMinvMixPPCopy);
                      hRCoef->Scale(0.5);

                      // Save results in mergeable collection
                      if ( hRCoef ) {
                        o = fMergeableCollection->GetObject(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hRCoef->GetName());
                        AliDebug(1,Form("----idMix=MIX_%s_%s%s o=%p",seventType->String().Data(),strigger->String().Data(),idMix.Data(),o));

                        if (o) fMergeableCollection->Remove(Form("/MIX/%s_%s%s/%s",seventType->String().Data(),strigger->String().Data(),idMix.Data(),hRCoef->GetName()));

                        Bool_t adoptOK = fMergeableCollection->Adopt(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hRCoef);
                        if(!adoptOK) AliError(Form("Could not adopt histo %s",hRCoef->GetName()));
                      } else {
                        for (int i=0; i<6; i++){ if (hTableMinv[i]!=0x0) delete hTableMinv[i]; }
                        for (int i=0; i<6; i++){ if (hTableDistRaw[i]!=0x0) delete hTableDistRaw[i]; }
                        for (int i=0; i<6; i++){ if (hTableDistMix[i]!=0x0) delete hTableDistMix[i]; }
                        delete hMinvMixPPCopy;
                        continue;
                      }

                      // ---- normalize unlike-sign mixed minv histo ----

                      // Copy and Multiply like sign histograms
                      TH1* hMinvPPCopy = static_cast<TH1*>(hTableMinv[1]->Clone());
                      hMinvPPCopy->Multiply(hTableMinv[2]);
                      nEntries = hMinvPPCopy->GetEntries();

                      // Compute and fill square roots of like-signs histo
                      for (int i = 0; i < nEntries ; i++){
                        Double_t binContent = hMinvPPCopy->GetBinContent(i);
                        Double_t binError   = hMinvPPCopy->GetBinError(i);
                        if( binContent !=0 && binError !=0 ){
                          hMinvPPCopy->SetBinContent(i,TMath::Sqrt(binContent));
                          hMinvPPCopy->SetBinError(i,TMath::Sqrt(binContent)*(binError/binContent));
                        } /*else {
                          hMinvPPCopy->SetBinContent(i,0.);
                          hMinvPPCopy->SetBinError(i,0.);
                        }*/
                      }

                      // Multiply by 2R
                      hMinvPPCopy->Multiply(hRCoef);
                      hMinvPPCopy->Scale(2.);

                      Double_t intMix = 0;
                      Double_t intRaw = 0;

                      Int_t bmin  = hTableMinv[3]->GetXaxis()->FindBin(Mmin);
                      Int_t bmax  = hTableMinv[3]->GetXaxis()->FindBin(Mmax);
                      intMix      = hTableMinv[3]->Integral(bmin,bmax);

                      Int_t bmin2 = hMinvPPCopy->GetXaxis()->FindBin(Mmin);
                      Int_t bmax2 = hMinvPPCopy->GetXaxis()->FindBin(Mmax);
                      intRaw      = hMinvPPCopy->Integral(bmin2,bmax2);

                      if ( intMix !=0. && intRaw !=0. && intRaw/intMix !=0. ) {
                        for (int i = 0; i < 6; ++i) {
                          if ( i<3 ) hTableMinv[i+3]->Scale(intRaw/intMix); // Norm MinvMix histo
                          hTableDistMix[i]->Scale(intRaw/intMix); // Norm pt and rapidy mix histo
                        }
                      } else {
                        AliError(Form("\n Cannot compute integral from one of the histos since intRaw = %f and intMix = %f \n",intRaw,intMix));
                        for (int i=0; i<6; i++){ if (hTableMinv[i]!=0x0) delete hTableMinv[i]; }
                        for (int i=0; i<6; i++){ if (hTableDistRaw[i]!=0x0) delete hTableDistRaw[i]; }
                        for (int i=0; i<6; i++){ if (hTableDistMix[i]!=0x0) delete hTableDistMix[i]; }
                        delete hMinvMixPPCopy;
                        delete hMinvPPCopy;
                        continue;
                      }

                      // ---- save results in mergeable collection ----

                      for (int i = 0; i < 6; ++i){
                        o = 0x0;
                        if ( hTableMinv[i] ) {
                          o = fMergeableCollection->GetObject(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hTableMinv[i]->GetName());
                          AliDebug(1,Form("----idMix=MIX_%s_%s%s o=%p",seventType->String().Data(),strigger->String().Data(),idMix.Data(),o));

                          if (o) fMergeableCollection->Remove(Form("/MIX/%s_%s%s/%s",seventType->String().Data(),strigger->String().Data(),idMix.Data(),hTableMinv[i]->GetName()));

                          Bool_t adoptOK = fMergeableCollection->Adopt(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hTableMinv[i]);
                          if(! adoptOK) AliError(Form("Cannot adopt histo %s ",hTableMinv[i]->GetName()));
                        }
                        o = 0x0;
                        if ( hTableDistRaw[i] ) {
                          o = fMergeableCollection->GetObject(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hTableDistRaw[i]->GetName());
                          AliDebug(1,Form("----idMix=MIX_%s_%s%s o=%p",seventType->String().Data(),strigger->String().Data(),idMix.Data(),o));

                          if (o) fMergeableCollection->Remove(Form("/MIX/%s_%s%s/%s",seventType->String().Data(),strigger->String().Data(),idMix.Data(),hTableDistRaw[i]->GetName()));

                          Bool_t adoptOK = fMergeableCollection->Adopt(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hTableDistRaw[i]);
                          if(! adoptOK) AliError(Form("Cannot adopt histo %s ",hTableDistRaw[i]->GetName()));
                        }
                        o = 0x0;
                        if ( hTableDistMix[i] ) {
                          o = fMergeableCollection->GetObject(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hTableDistMix[i]->GetName());
                          AliDebug(1,Form("----idMix=MIX_%s_%s%s o=%p",seventType->String().Data(),strigger->String().Data(),idMix.Data(),o));

                          if (o) fMergeableCollection->Remove(Form("/MIX/%s_%s%s/%s",seventType->String().Data(),strigger->String().Data(),idMix.Data(),hTableDistMix[i]->GetName()));

                          Bool_t adoptOK = fMergeableCollection->Adopt(Form("/MIX/%s_%s%s",seventType->String().Data(),strigger->String().Data(),idMix.Data()),hTableDistMix[i]);
                          if(! adoptOK) AliError(Form("Cannot adopt histo %s ",hTableDistMix[i]->GetName()));
                        }
                      }
                    }
                    delete binning;
                    delete bins;
                  }
                }
              }
            }
          }
        }
      }
    }

    Update();

    delete eventTypeArray;
    delete eventTypeMixArray;
    delete triggerMixArray;
    delete pairCutArray;
    delete centralityArray;
    delete triggerArray;

    return;
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMu::FitJpsi(const char* binType, const char* flavour, const TString fitMethod, const char* histoType )
{
/**
 * @brief  Fit the J/psi and Psi(2S) pick.
 * @details Fit process delegated to AliAnalysisMuMu::FitParticle().
 *
 * @param binType  [description]
 * @param flavour  [description]
 * @param fitMethod '' or 'mix'
 * @return number of proceeded fits
 */
  TStopwatch timer;

  if (!fMergeableCollection) {
    AliError("No mergeable collection. Consider Upgrade()");
    return 0;
  } else {
      cout <<      " ================================================================ " << endl;
      cout <<      "                       FitJpsi                  " << endl;
      cout <<      " ================================================================ " << endl;
  }

  const char* particle = "psi";

  Int_t nfits(0);

  TObjArray* pairCutArray    = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
  TObjArray* centralityArray = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
  TObjArray* binTypeArray    = TString(binType).Tokenize(",");

  TObjArray* TriggerArray(0x0);
  if ( !fitMethod.Contains("mix") ) TriggerArray = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
  else                              TriggerArray = Config()->GetListElements(Config()->MixTriggerKey(),IsSimulation());

  TObjArray* eventTypeArray(0x0);
  if ( !fitMethod.Contains("mix") ) eventTypeArray = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
  else                              eventTypeArray = Config()->GetListElements(Config()->EventSelectionMixKey(),IsSimulation());

  TString refTrigger(Form("%s",Config()->First(Config()->RefMixTriggerKey(),IsSimulation()).Data()));
  TString refEvent(Form("%s",Config()->First(Config()->RefMixEventSelectionKey(),IsSimulation()).Data()));

  TIter nextTrigger(TriggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextbinType(binTypeArray);
  TIter nextCentrality(centralityArray);

  TObjString* trigger;
  TObjString* eventType;
  TObjString* pairCut;
  TObjString* sbinType;
  TObjString* centrality;

  while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
  {
    AliAnalysisMuMuBinning* binning(0x0);

    if ( fBinning && sbinType->String().Length() > 0 ) {
      binning = fBinning->Project(particle,sbinType->String().Data(),flavour);
    } else {
      binning = new AliAnalysisMuMuBinning;
      binning->AddBin(particle,sbinType->String().Data());
    }

    StdoutToAliDebug(1,std::cout << "++++++++++++ binning=" << sbinType->String().Data() << std::endl;);

    std::cout << "" << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "+++++++++++++++++++ binning  = " << sbinType->String().Data() << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    if (!binning) {
      AliError("oups. binning is NULL");
      continue;
    }

    StdoutToAliDebug(1,binning->Print(););

    // Loop over all trigger
    while ( ( trigger = static_cast<TObjString*>(nextTrigger())) ) {
      AliDebug(1,Form("--TRIGGER %s",trigger->String().Data()));
      nextEventType.Reset();

      // Loop over all evenType
      while ( ( eventType = static_cast<TObjString*>(nextEventType())) ){
        AliDebug(1,Form("----EVENTTYPE %s",eventType->String().Data()));
        nextPairCut.Reset();

        // Loop over all paircut
        while ( ( pairCut = static_cast<TObjString*>(nextPairCut())) ) {
          AliDebug(1,Form("------PAIRCUT %s",pairCut->String().Data()));
          nextCentrality.Reset();

          // Loop over all centrality
          while ( ( centrality = static_cast<TObjString*>(nextCentrality()) ) ) {
            AliDebug(1,"------Fitting...");

            // Select stored path
            TObject* o;
            TString id = "";
            if(fitMethod.Contains("mix"))
              id = Form("/FitResults/%s_%s/%s/%s/%s/%s",refEvent.Data(),refTrigger.Data(),eventType->String().Data(),trigger->String().Data(),centrality->String().Data(),pairCut->String().Data());
            else
              id = Form("/FitResults/%s/%s/%s/%s",eventType->String().Data(),trigger->String().Data(),centrality->String().Data(),pairCut->String().Data());

            printf("\n ----- id:%s -----\n",id.Data() );
            AliAnalysisMuMuSpectra* spectra(0x0);
            AliAnalysisMuMuSpectra* spectraCorr(0x0);

            // ---- The main part. The fit method is called ----

            AliDebug(1,"------Fitting spectra...");
            spectra = FitParticle(particle,trigger->String().Data(),eventType->String().Data(),pairCut->String().Data(),centrality->String().Data(),*binning,kFALSE,&fitMethod,flavour,histoType);
            AliDebug(1,Form("------fitting done spectra = %p",spectra));

            // AliDebug(1,"------Fitting corrected spectra...");
            // spectraCorr = FitParticle(particle,trigger->String().Data(),eventType->String().Data(),pairCut->String().Data(),centrality->String().Data(),*binning,kTRUE,&fitMethod,flavour,histoType);
            // AliDebug(1,Form("------fitting done corrected spectra = %p",spectraCorr));

            // ----

            // --- save results in mergeable collection ---
            if ( spectra ) {
              ++nfits;

              o = fMergeableCollection->GetObject(id.Data(),spectra->GetName());
              AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));

              if (o) {
                AliWarning(Form("Replacing %s/%s",id.Data(),spectra->GetName()));
                fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectra->GetName()));
              }

              Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),spectra);

              if ( adoptOK ) std::cout << "+++Spectra " << spectra->GetName() << " adopted" << std::endl;
              else AliError(Form("Could not adopt spectra %s",spectra->GetName()));

              StdoutToAliDebug(1,spectra->Print(););
            } else AliError("Error creating spectra");

            o = 0x0;
            if ( spectraCorr ) {
              ++nfits;

              o = fMergeableCollection->GetObject(id.Data(),spectraCorr->GetName());
              AliDebug(1,Form("----nfits=%d id=%s o=%p",nfits,id.Data(),o));

              if (o) {
                AliWarning(Form("Replacing %s/%s",id.Data(),spectraCorr->GetName()));
                fMergeableCollection->Remove(Form("%s/%s",id.Data(),spectraCorr->GetName()));
              }

              Bool_t adoptOK = fMergeableCollection->Adopt(id.Data(),spectraCorr);

              if ( adoptOK ) std::cout << "+++Spectra " << spectraCorr->GetName() << " adopted" << std::endl;
              else AliError(Form("Could not adopt spectra %s",spectraCorr->GetName()));

              StdoutToAliDebug(1,spectraCorr->Print(););
            }
            else AliError("Error creating spectra");
          }
        }
      }
    }
  }

  delete eventTypeArray;
  delete TriggerArray;
  delete pairCutArray;
  delete centralityArray;

  StdoutToAliDebug(1,timer.Print(););

  if (nfits) Update();

  return nfits;
}

//_____________________________________________________________________________
TGraph* AliAnalysisMuMu::PlotEventSelectionEvolution(const char* trigger1, const char* event1,
                                                     const char* trigger2, const char* event2,
                                                     Bool_t drawFills,
                                                     Bool_t asRejection) const
{
/**
 * @brief Plot relative Event Selection Evolution
 *
 * @param trigger1    [description]
 * @param event1      [description]
 * @param trigger2    [description]
 * @param event2      [description]
 * @param drawFills   [description]
 * @param asRejection [description]
 * @return TGraph to be handle by the user
 */
    if (!CC()) return 0x0;

    const std::set<int>& runnumbers = RunNumbers();

    TGraphErrors* g = new TGraphErrors(runnumbers.size());

    std::set<int>::const_iterator it;
    Int_t i(0);

    Double_t ymin(TMath::Limits<double>::Max());
    Double_t ymax(TMath::Limits<double>::Min());

    for ( it = runnumbers.begin(); it != runnumbers.end(); ++it )
        {
        Int_t runNumber = *it;
        Double_t n = CC()->GetSum(Form("trigger:%s/event:%s/run:%d",trigger1,event1,runNumber));
        Double_t d = CC()->GetSum(Form("trigger:%s/event:%s/run:%d",trigger2,event2,runNumber));
        if (n>0 && d>0)
            {
            Double_t y = n/d;

            if ( fCorrectionPerRun )
                {
                Double_t xcorr,ycorr;
                fCorrectionPerRun->GetPoint(i,xcorr,ycorr); // note that the fact that xcorr==runNumber has been checked by the SetCorrectionPerRun method
                y *= ycorr;
                // FIXME: should get the correction error here
                }

            if ( asRejection ) y = 100*(1.0 - y);
            ymin = TMath::Min(ymin,y);
            ymax = TMath::Max(ymax,y);
            Double_t yerr = y*AliAnalysisMuMuResult::ErrorAB(n,TMath::Sqrt(n),d,TMath::Sqrt(d));
            g->SetPoint(i,runNumber,y);
            g->SetPointError(i,0.5,yerr);

            ++i;
            }

        }

    TH2* hframe = new TH2F(Form("%s %s-%s / %s-%s",(asRejection ? "1 - ":""),trigger1,event1,trigger2,event2),
                           Form("%s %s-%s / %s-%s",(asRejection ? "1 - ":""),trigger1,event1,trigger2,event2),
                           runnumbers.size()+50,
                           *(runnumbers.begin())-25,
                           *(runnumbers.rbegin())+25,100,0,ymax*1.3);

    gStyle->SetOptStat(0);

    hframe->Draw();

    hframe->GetXaxis()->SetNoExponent();

    hframe->GetYaxis()->SetTitle(asRejection ? "Rejection (%)" : "Ratio");

    g->Set(i);
    g->SetTitle(Form("%s %s-%s / %s-%s",(asRejection ? "1 - ":""),trigger1,event1,trigger2,event2));
    g->GetXaxis()->SetNoExponent();
    g->Draw("lp");

    AliAnalysisTriggerScalers ts(RunNumbers(),Config()->OCDBPath());

    if ( drawFills )
        {
        ts.DrawFills(ymin,ymax);
        g->Draw("lp");
        }


    std::map<std::string, std::pair<int,int> > periods;

    ts.GetLHCPeriodBoundaries(periods);

    TLegend* legend = new TLegend(0.15,0.82,0.90,0.92);
    legend->SetFillColor(0);
    Int_t n(0);


    for ( std::map<std::string, std::pair<int,int> >::const_iterator pit = periods.begin(); pit != periods.end(); ++pit )
        {
        std::string period = pit->first;
        int run1 = (pit->second).first;
        int run2 = (pit->second).second;
        int nruns(0);
        for ( std::set<int>::const_iterator rit = RunNumbers().begin(); rit != RunNumbers().end(); ++ rit )
            {
            if ( (*rit) >= run1 && (*rit) <= run2 )
                {
                ++nruns;
                }
            }
        AliInfo(Form("Period %s runs %6d-%6d ; %d actual runs",period.c_str(),run1,run2,nruns));

        g->Fit("pol0","+Q","",run1,run2);
        TF1* func = static_cast<TF1*>(g->GetListOfFunctions()->Last());
        if (func)
            {
            func->SetLineColor(2+n);
            legend->AddEntry(func,Form("%s %5.2f #pm %5.2f %s (rel. error %5.2f %%)",period.c_str(),func->GetParameter(0),func->GetParError(0),
                                       (asRejection ? "%":""),100*func->GetParError(0)/func->GetParameter(0)));
            ++n;
            }
        }

    legend->SetNColumns(3);

    Double_t mean = TMath::Mean(g->GetN(),g->GetY());
    Double_t rms = TMath::RMS(g->GetN(),g->GetY());

    legend->AddEntry("",Form("Mean %5.2f RMS %5.2f (%5.2f %%)",mean,rms,(mean) ? 100.0*rms/mean : 0.0),"");

    legend->Draw();

    return g;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ShowList(const char* title, const TString& list, const char separator) const
{
  /**
   * @brief convenient method to show a string as a list providing the separator.
   */
  TObjArray* parts = list.Tokenize(separator);

  std::cout << title << " (" << parts->GetEntries() << ") : " << std::endl;

  TIter next(parts);
  TObjString* str;

  while ( ( str = static_cast<TObjString*>(next()) ) ) std::cout << "    " << str->String().Data() << std::endl;

  if ( parts->GetEntries()==0) std::cout << endl;

  delete parts;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Print(Option_t* opt) const
{
  /// printout

  std::cout << "Reading from file : " << fFilename.Data() << std::endl;

  TString copt(opt);

  if (IsSimulation() || SIM() )copt += "SIM";
  if ( !IsSimulation() )copt += "REAL";

  Config()->Print(copt.Data());

  if ( RunNumbers().size() > 1 )std::cout << RunNumbers().size() << " runs";
  else std::cout << RunNumbers().size() << " run";


  if ( fCorrectionPerRun ) std::cout << " with correction factors";

  std::cout << std::endl;
  Int_t i(0);
  for ( std::set<int>::const_iterator it = RunNumbers().begin(); it != RunNumbers().end(); ++it ){
    std::cout << (*it);
    if ( fCorrectionPerRun ) std::cout << Form("(%e)",fCorrectionPerRun->GetY()[i]);
    std::cout << ",";
    ++i;
  }
  std::cout << std::endl;

  TString sopt(opt);
  sopt.ToUpper();

  if ( sopt.Contains("BIN") && BIN() ){
    std::cout << "Binning : " << std::endl;
    TString topt(sopt);
    topt.ReplaceAll("BIN","");
    BIN()->Print(topt.Data());
  }

  if ( sopt.Contains("MC") && OC() ){
    TString topt(sopt);
    topt.ReplaceAll("MC","");
    OC()->Print(topt.Data());
  }

  if ( sopt.Contains("CC") && CC() )CC()->Print("trigger/event");


  if ( sopt.Contains("SIZE") ){
    TFile* f = ReOpen(fFilename,"READ");
    TIter next(f->GetListOfKeys());
    TKey* key;

    while ( ( key = static_cast<TKey*>(next()) ) )std::cout << key->GetName() << " " << key->GetNbytes() << " " << key->GetObjlen() << std::endl;
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMu::PrintNofWhat(const char* what, const char* spectraName, Bool_t mix, Bool_t AccEffCorr) const
{
    /**
    * @brief            Print number of particle. Method to use after FitJPsi() or any other fit process.
    * @details          Delegate to AliAnalysisMuMuCapsulePbPb::PrintNumberOfWhat().
    * @param particle   particle name
    * @param what Quantity stored in the AliAnalysisMuMuSpectra
    */


    if (!OC() || !CC()){
      AliError("No me rgeable/counter collection. Consider Upgrade()");
      return ;
    } else {
      cout <<      " ================================================================ " << endl;
      cout <<      "                       Number of "<< what << endl;
      cout <<      " ================================================================ " << endl;
    }

    // Get configuration settings

    TObjArray* pairCutArray    = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
    TObjArray* centralityArray = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());

    TObjArray* triggerArray(0x0);
    if (!mix) triggerArray = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
    else                              triggerArray = Config()->GetListElements(Config()->MixTriggerKey(),IsSimulation());

    TObjArray* eventTypeArray(0x0);
    if (!mix) eventTypeArray = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
    else                              eventTypeArray = Config()->GetListElements(Config()->EventSelectionMixKey(),IsSimulation());

    TString refTrigger(Form("%s",Config()->First(Config()->RefMixTriggerKey(),IsSimulation()).Data()));
    TString refEvent(Form("%s",Config()->First(Config()->RefMixEventSelectionKey(),IsSimulation()).Data()));

    TObjArray* bins;

    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextEventType(eventTypeArray);
    TIter nextPairCut(pairCutArray);
    TIter nextCentrality(centralityArray);

    // Strings
    TObjString* strigger;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* sparticle;
    TObjString* scentrality;

    nextEventType.Reset();
    // Loop on each envenType (see MuMuConfig)
    //==============================================================================
    while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
    {
      AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
      nextTrigger.Reset();
      // Loop on each trigger (see MuMuConfig)
      //==============================================================================
      while ( ( strigger = static_cast<TObjString*>(nextTrigger())) )
      {
      AliDebug(1,Form("-TRIGGER %s",strigger->String().Data()));
      nextCentrality.Reset();
      // Loop on each centrality (not the ones in MuMuConfig but the ones set)
      //==============================================================================
      while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
      {
        AliDebug(1,Form("--CENTRALITY %s",scentrality->String().Data()));
        nextPairCut.Reset();
        // Loop on each paircut (not the ones in MuMuConfig but the ones set)
        //==============================================================================
        while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
        {
          AliDebug(1,Form("---PAIRCUT %s",spairCut->String().Data()));

          //________Get spectra
          TString spectraPath = "";
          if(mix)
            spectraPath = Form("/FitResults/%s_%s/%s/%s/%s/%s/%s",refEvent.Data(),refTrigger.Data(),seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName);
          else
            spectraPath = Form("/FitResults/%s/%s/%s/%s/%s",seventType->String().Data(),strigger->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName);
          if (AccEffCorr)spectraPath+="-AccEffCorr";

          AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));

          if(!spectra){
            AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
            continue;
          }
          //________

          // Create pointer on fitted spectra
          AliAnalysisMuMuSpectraCapsulePbPb * capsule = new AliAnalysisMuMuSpectraCapsulePbPb(spectra,spectraPath,"","");
          AliDebug(1,Form("capsule = %p",capsule));
          if(!capsule){
            AliError("Could not find spetra !");
            continue;
          }

          capsule->PrintNofWhat(what);
          delete capsule;
        }
      }
    }
  }


  delete eventTypeArray ;
  delete triggerArray ;
  delete pairCutArray ;
  delete centralityArray ;

  return ;

}

//_____________________________________________________________________________
TFile* AliAnalysisMuMu::ReOpen(const char* filename, const char* mode) const
{
  /**
   * @brief Tries to reopen the file with a new mode
   */

  TFile* f = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(filename));
  if (f)delete f;

  f = TFile::Open(filename,mode);

  if ( !f || !f->IsOpen() ){
    AliError(Form("Cannot open file %s in mode %s",filename,mode));
    return 0x0;
  }

  return f;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::SetCorrectionPerRun(const TGraph& corr, const char* formula)
{
/**
 * @brief Sets the graph used to correct values per run
 *
 * @param corr    [description]
 * @param formula [description]
 *
 */
    delete fCorrectionPerRun;
    fCorrectionPerRun=0x0;

    // check that corr has the same runs as we do

    Int_t i(0);

    for ( std::set<int>::const_iterator it = RunNumbers().begin(); it != RunNumbers().end(); ++it ){
      Int_t corrRun = TMath::Nint(corr.GetX()[i]);

      if (corrRun != *it){
        AliError(Form("%d-th run mistmatch %d vs %d",i,corrRun,*it));
        return kFALSE;
      }
      ++i;
    }

    fCorrectionPerRun = new TGraphErrors(corr.GetN());

    TFormula* tformula(0x0);
    if ( strlen(formula) > 0 ) tformula = new TFormula("SetCorrectionPerRunFormula",formula);

    i = 0;

    for ( std::set<int>::const_iterator it = RunNumbers().begin(); it != RunNumbers().end(); ++it ){
      Double_t y = corr.GetY()[i];

      if ( tformula ) y = tformula->Eval(y);
      fCorrectionPerRun->SetPoint(i,corr.GetX()[i],y);
      ++i;
    }

    delete formula;

    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SetNofInputParticles(AliAnalysisMuMuJpsiResult& r,const char* event, const char* trigger, const char* centrality)
{
  /// Set the "NofInput" variable(s) of one result

  TString hname(Form("MinvUS+%s",r.Bin().AsString().Data()));

  TH1* hinput = fMergeableCollection->Histo(Form("/%s/%s/%s/%s/INYRANGE",AliAnalysisMuMuBase::MCInputPrefix(),event,trigger,centrality),hname.Data());

  if (!hinput)
  {
    AliError(Form("Got a simulation file where I did not find histogram /%s/%s/%s/%s/INYRANGE/%s",AliAnalysisMuMuBase::MCInputPrefix(),event,trigger,centrality,hname.Data()));
  }
  else
  {
    r.SetNofInputParticles(*hinput);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::SPECTRA(const char* fullpath) const
{
  /**
   * @brief Shortcut method to get to a spectra
   */
  if (!OC()) return 0x0;

  return static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(fullpath));
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SelectRunByTrigger(const char* triggerList)
{
  /**
  *   [AliAnalysisMuMu::SelectRunByTrigger description]
  *   @brief   Old function with a lot of hard coded things, not sure what it does... Maybe outdated
  *   @details [long   description]
  *
  *   @param   triggerList [description]
  */
  if (!fMergeableCollection || !fCounterCollection) return;

  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);

  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);

  TObjString* srun;
  TObjString* strigger;

  TString striggerList(triggerList);

  TList* runList = new TList();

  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {

    TList* runList = new TList();

    while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
    {

      nextTrigger.Reset();

      while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
      {
        if ( !striggerList.Contains(strigger->String().Data()) ) continue;

        ULong64_t n = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                                  strigger->String().Data(),"PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00",srun->String().Atoi())));
        if ( n > 0 ) runList->Add(srun);

      }
    }
  }
    runList->Sort();
    TIter nextRunOK(runList);
    while ( ( srun = static_cast<TObjString*>(nextRunOK()) ) )
    {
      std::cout << srun->String().Atoi() << std::endl;
    }
    delete runList;

    delete triggers;
    delete runs;
}
//_____________________________________________________________________________
void AliAnalysisMuMu::TriggerCountCoverage(const char* triggerList,
                                           Bool_t compact,
                                           Bool_t orderByTriggerCount)
{
/**
 *   [AliAnalysisMuMu::TriggerCountCoverage description]
 *   @brief   Give the fraction of triggers (in triggerList) relative to what is expected in the scalers
 *
 *   @param   triggerList         [description]
 *   @param   compact             [description]
 *   @param   orderByTriggerCount [description]
 */
  TGrid::Connect("alien://"); // to insure the "Trying to connect to server... message does not pollute our output later on...

  AliLog::EType_t oldLevel = static_cast<AliLog::EType_t>(AliLog::GetGlobalLogLevel());

  AliLog::SetGlobalLogLevel(AliLog::kFatal);

  if (!fMergeableCollection || !fCounterCollection) return;

  TObjArray* runs = fCounterCollection->GetKeyWords("run").Tokenize(",");
  TIter nextRun(runs);

  TObjArray* triggers = fCounterCollection->GetKeyWords("trigger").Tokenize(",");
  TIter nextTrigger(triggers);

  TObjString* srun;
  TObjString* strigger;

  TString striggerList(triggerList);

  ULong64_t total(0);
  ULong64_t totalExpected(0);
  TString msg;
  std::multimap<ULong64_t,std::string> messages;

  while ( ( srun = static_cast<TObjString*>(nextRun()) ) )
  {
    msg.Form("RUN %09d ",srun->String().Atoi());

    if (!compact) msg += "\n";

    ULong64_t nmax(0);

    nextTrigger.Reset();

    while ( ( strigger = static_cast<TObjString*>(nextTrigger()) ) )
    {
      if ( !striggerList.Contains(strigger->String().Data()) ) continue;

      ULong64_t n = TMath::Nint(fCounterCollection->GetSum(Form("trigger:%s/event:%s/run:%d",
                                                            strigger->String().Data(),"ALL",srun->String().Atoi())));

      ULong64_t expected = GetTriggerScalerCount(strigger->String().Data(),srun->String().Atoi());


      nmax = TMath::Max(n,nmax);

      total += n;
      totalExpected += expected;

      msg += TString::Format("%30s %9lld expected %9lld [%s] ",strigger->String().Data(),n,expected,
                             (n>expected ? "!" : " "));

      if ( expected > 0 ) {
        msg += TString::Format("fraction %5.1f %%",n*100.0/expected);
      }

      if (!compact) msg += "\n";
    }

    if (nmax>0){
      if (!orderByTriggerCount) std::cout << msg.Data() << std::endl;
      else messages.insert(std::make_pair(nmax,static_cast<std::string>(msg.Data())));
    }
  }

  std::multimap<ULong64_t,std::string>::const_reverse_iterator it;

  ULong64_t current(0);
  Int_t n(0);

  for ( it = messages.rbegin(); it != messages.rend(); ++it ){
    ++n;
    current += it->first;
    Double_t percent = ( total > 0.0 ? current*100.0/total : 0.0);
    std::cout << Form("%10lld",it->first) << " " << it->second << " percentage of total = " << Form("%7.2f %% %3d",percent,n ) << std::endl;
  }

  std::cout << Form("--- TOTAL %lld expected %lld fraction %5.1f %%",
                    total,totalExpected,totalExpected ? total*100.0/totalExpected : 0.0) << std::endl;


    for ( it = messages.rbegin(); it != messages.rend(); ++it ) {
      ++n;
      current += it->first;
      Double_t percent = ( total > 0.0 ? current*100.0/total : 0.0);
      std::cout << Form("%10lld",it->first) << " " << it->second << " percentage of total = " << Form("%7.2f %% %3d",percent,n ) << std::endl;
    }

    std::cout << Form("--- TOTAL %lld expected %lld fraction %5.1f %%",
                      total,totalExpected,totalExpected ? total*100.0/totalExpected : 0.0) << std::endl;



    AliLog::SetGlobalLogLevel(oldLevel);
    delete triggers;
    delete runs;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::UnsetCorrectionPerRun()
{
    // drop the correction factors
    delete fCorrectionPerRun;
    fCorrectionPerRun=0x0;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::Update()
{
  /// update the current file with memory

  if (!CC() || !OC()) return;

  ReOpen(fFilename,"UPDATE");

  if (OC()){
    if (fDirectory.Length()) gDirectory->cd(fDirectory.Data());
    OC()->Write("OC",TObject::kSingleKey|TObject::kOverwrite);
  }

  ReOpen(fFilename,"READ");

  GetCollections(fFilename,fDirectory,fMergeableCollection,fCounterCollection,fBinning,fRunNumbers);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::Upgrade(const char* filename)
{
    /// Upgrade a file
    AliAnalysisMuMu m(filename);

    return m.Upgrade();
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::Upgrade()
{
    /// Upgrade the current file
    /// - from single list to one key per object, if needed
    /// - from histogramCollection to mergeableCollection, if needed

    AliWarning("Out of date method");

    TFile* f = ReOpen(fFilename,"UPDATE");

    TList* list = static_cast<TList*>(f->Get("chist"));

    if (list) {
      // really old file where everything was in a single list

      AliHistogramCollection* hc = static_cast<AliHistogramCollection*>(list->At(0));
      AliCounterCollection* cc = static_cast<AliCounterCollection*>(list->At(1));

      AliMergeableCollection* mc = hc->Convert();

    AliWarning("Out of date method");

      mc->Write("MC",TObject::kSingleKey);
      cc->Write("CC",TObject::kSingleKey);

      f->Delete("chist;*");

      f->Write();

    } else {
      AliHistogramCollection* hc = static_cast<AliHistogramCollection*>(f->Get("HC"));

      if ( hc ) {
        // old file with histogram collection instead of mergeable collection

        AliMergeableCollection* mc = hc->Convert();

        f->cd();

        mc->Write("MC",TObject::kSingleKey);
        f->Delete("HC;*");
        f->Write();
      }
    }

  delete f;
  return kTRUE;
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMu::PrintDistribution(const char* spectraName, const char* what, const char* subresultname, Bool_t divideByBinWidth, Bool_t AccEffCorr)
{
    /**
     * @brief    Print distribution vs binType.
     * @details  Delegate procedure to AliAnalysisSpectra object. DimuonTriggerKey is used to select the trigger.
     *
     *  @param   binType          [description]
     *  @param   what               Quantity stored in the AliAnalysisMuMuSpectra
     *  @param   subresultname    If one wants to draw one/several specific results.
     *  @param   divideByBinWidth [description]
     *  @param   AccEffCorr       Spectra type
    **/
    if (!OC() || !CC()){
      AliError("No mergeable/counter collection. Consider Upgrade()");
      return 0x0;
    } else {
      cout <<      " ================================================================ " << endl;
      cout <<      "                       PrintDistribution                          " << endl;
      cout <<      " ================================================================ " << endl;
    }

    // Get configuration settings
    TObjArray* eventTypeArray   = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
    TObjArray* triggerArray     = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
    TObjArray* fitfunctionArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());// to add here an entry
    TObjArray* pairCutArray     = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
    TObjArray* centralityArray  = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
    TObjArray* whatArray        = TString(what).Tokenize(",");
    TObjArray* spectraNameArray     = TString(spectraName).Tokenize(",");
    TObjArray* bins;

    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextEventType(eventTypeArray);
    TIter nextPairCut(pairCutArray);
    TIter nextWhat(whatArray);
    TIter nextspectraName(spectraNameArray);
    TIter nextCentrality(centralityArray);

    // Strings
    TObjString* striggerDimuon;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* swhat;
    TObjString* sspectraName;
    TObjString* scentrality;


    // Pointers

    TObjArray* h= new TObjArray();
    h->SetOwner(kTRUE);
    int i_TH1 = 1 ;

    AliAnalysisMuMuSpectra spectra=0x0;

    //Loop on what type
    while ( ( swhat = static_cast<TObjString*>(nextWhat()) ) )
    {
      AliDebug(1,Form("what %s",swhat->String().Data()));
      nextEventType.Reset();
      // Loop on each envenType (see MuMuConfig)
      //==============================================================================
      while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
      {
        AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
        nextTrigger.Reset();
        // Loop on each trigger (see MuMuConfig)
        //==============================================================================
        while ( ( striggerDimuon = static_cast<TObjString*>(nextTrigger())) )
        {
          AliDebug(1,Form("-TRIGGER %s",striggerDimuon->String().Data()));
          nextPairCut.Reset();
          // Loop on each paircut (not the ones in MuMuConfig but the ones set)
          //==============================================================================
          while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
          {
            AliDebug(1,Form("--PAIRCUT %s",spairCut->String().Data()));
            nextspectraName.Reset();
            // Loop on each type (pt or y)
            //==============================================================================
            while ( ( sspectraName = static_cast<TObjString*>(nextspectraName()) ) )
            {
              AliDebug(1,Form("---TYPE %s",sspectraName->String().Data()));

              // //canvas
              // TCanvas *c1 = new TCanvas;
              // c1->Draw();
              // gStyle->SetOptStat(0);

              nextCentrality.Reset();
              // Loop on each centrality (not the ones in MuMuConfig but the ones set)
              //==============================================================================
              while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
              {
                AliDebug(1,Form("---CENTRALITY %s",scentrality->String().Data()));

                //________Get spectra
                TString spectraPath= Form("/FitResults/%s/%s/%s/%s/%s",seventType->String().Data(),striggerDimuon->String().Data(),scentrality->String().Data(),spairCut->String().Data(),sspectraName->String().Data());
                if (AccEffCorr)spectraPath+="-AccEffCorr";
                AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));
                if(!spectra)
                    {
                    AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
                    return 0x0;
                    }
                //________

                h->Add(((TObject*)spectra->Plot(swhat->String().Data(),subresultname,divideByBinWidth)));
                TCanvas *c = new TCanvas;
                // c->SetLogy();
                spectra->Plot(swhat->String().Data(),subresultname,divideByBinWidth)->Draw();
                ++i_TH1;
              }
            }
          }
        }
      }
    }

    return h;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeDimuonRawCount(const Double_t rlow, const Double_t rhigh, const char* binType, const char* binRangeExluded, const char* flavour, Bool_t corrected)
{
/**
 *   [AliAnalysisMuMu::ComputeDimuonRawCount description]
 *   @brief   Compute the raw count of dimuon pair
 *
 *   @param   rlow            Low intervale limit in bin for the raw count
 *   @param   rhigh           High intervale limit in bin for the raw count
 *   @param   binType         [description]
 *   @param   binRangeExluded Excluded range inside interval
 *   @param   flavour         [description]
 *   @param   corrected       histo type
 */
    if(!OC())
    {
        AliError("No mergeable. Consider Upgrade()");
        return;
    }
    else
    {
        cout <<      " ================================================================ " << endl;
        cout <<      "                       ComputeDimuonRawCount                          " << endl;
        cout <<      " ================================================================ " << endl;
    }

     // Get configuration settings
    TObjArray* eventTypeArray   = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
    TObjArray* triggerArray     = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
    TObjArray* fitfunctionArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());// to add here an entry
    TObjArray* pairCutArray     = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
    TObjArray* centralityArray  = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
    TObjArray* binTypeArray     = TString(binType).Tokenize(",");
    TObjArray* binRangeExcludedArray     = TString(binRangeExluded).Tokenize(",");
    TObjArray* bins;

    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextEventType(eventTypeArray);
    TIter nextPairCut(pairCutArray);
    TIter nextbinType(binTypeArray);
    TIter nextCentrality(centralityArray);

    // Strings
    TObjString* striggerDimuon;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* sbinType;
    TObjString* scentrality;

    TString striggerMB  = Config()->First(Config()->MinbiasTriggerKey(),IsSimulation());

    // Pointers
    TH1* h= 0x0;
    AliAnalysisMuMuSpectra spectra=0x0;


    // Loop on each envenType (see MuMuConfig)
    while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
    {
        AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
        nextTrigger.Reset();
        // Loop on each trigger (see MuMuConfig)
        while ( ( striggerDimuon = static_cast<TObjString*>(nextTrigger())) )
        {
            AliDebug(1,Form("-TRIGGER %s",striggerDimuon->String().Data()));
            nextPairCut.Reset();
            // Loop on each paircut (not the ones in MuMuConfig but the ones set)
            while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
            {
                AliDebug(1,Form("--PAIRCUT %s",spairCut->String().Data()));
                nextbinType.Reset();
                // Loop on each type (pt or y)
                while ( ( sbinType = static_cast<TObjString*>(nextbinType()) ) )
                {
                    AliDebug(1,Form("---TYPE %s",sbinType->String().Data()));
                    nextCentrality.Reset();


                     AliAnalysisMuMuBinning* binning(0x0);

                    if ( fBinning && sbinType->String().Length() > 0 )
                    {
                        binning = fBinning->Project("psi",sbinType->String().Data(),flavour);
                    }
                    else
                    {
                        binning = new AliAnalysisMuMuBinning;
                        binning->AddBin("psi",sbinType->String().Data());
                    }

                    // Check Binning list
                    TObjArray* bins = binning->CreateBinObjArray("psi");
                    if (!bins)
                    {
                        AliError(Form("Did not get any bin for particle psi"));
                        return;
                    }
                    // Loop on each centrality (not the ones in MuMuConfig but the ones set)
                    while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
                    {
                        AliDebug(1,Form("----CENTRALITY %s",scentrality->String().Data()));
                        // The binning pointer, which point at Pt binning, Y binning etc.
                        AliAnalysisMuMuBinning::Range* sbin;
                        TIter next(bins);
                        next.Reset();

                        // Create ID for the fit which will be used to name results
                        TString id(Form("/%s/%s/%s/%s",seventType->String().Data(),striggerDimuon->String().Data(),scentrality->String().Data(),spairCut->String().Data()));

                        // Histo where we put the count
                        const Double_t * binArrayX = binning->CreateBinArrayX();
                        Int_t nBinX = binning->GetNBinsX();

                        TH1* hraw = new TH1D(Form("hRawCountVS%s_%f-%f",sbinType->String().Data(),rlow,rhigh),Form("raw count of dimuon pairs for %s",id.Data()),nBinX,binArrayX);
                        hraw->GetYaxis()->SetTitle(Form("raw count of dimuon pairs in [%0.2f;%0.2f] GeV/c",rlow,rhigh));
                        hraw->GetXaxis()->SetTitle(sbinType->String().Data());

                        // Loop on each range in bin
                        while ( ( sbin = static_cast<AliAnalysisMuMuBinning::Range*>(next())) )
                        {
                            AliDebug(1,Form("-----Bin range %s",sbin->AsString().Data()));

                            if(binRangeExcludedArray->FindObject(sbin->AsString().Data())) {
                               AliDebug(1,Form("-----Bin range is excluded"));
                               continue;
                            }
                            // Create ID for the fit which will be used to name results
                            TString id(Form("/%s/%s/%s/%s",seventType->String().Data(),striggerDimuon->String().Data(),scentrality->String().Data(),spairCut->String().Data()));
                            TString hname = corrected ? Form("MinvUS+%s_AccEffCorr",sbin->AsString().Data()) : Form("MinvUS+%s",sbin->AsString().Data());

                            // Pointer to the histo from histo collection
                            h = OC()->Histo(id.Data(),hname.Data());
                            if (!h){
                                AliError(Form("Could not find histo %s",hname.Data()));
                                continue;
                            }

                            // Get X
                            Double_t xmin = sbin->Xmin();
                            Double_t xmax = sbin->Xmax();
                            Double_t x = xmin + (xmax-xmin)/2;
                            AliDebug(1,Form("x = %f \n", x));

                            //find bin
                            Int_t binLow   = h->FindBin(rlow);
                            Int_t binHight = h->FindBin(rhigh);

                            if(binLow==0 || binHight ==0)continue;
                            // Fill
                            Double_t rawCount =0.;
                            for (Int_t i = binLow; i < binHight; ++i){
                                rawCount = rawCount + h->GetBinContent(i);
                                AliDebug(1,Form("rawCount for %s in bin [%d]= %f \n",h->GetTitle(),i,rawCount));
                            }
                            if(rawCount!=0. )hraw->Fill(x,rawCount);

                        }
                        TH1* o = OC()->Histo(id.Data(),hraw->GetName());

                        if (o){
                            AliWarning(Form("Replacing %s/%s",id.Data(),hraw->GetName()));
                            OC()->Remove(Form("%s/%s",id.Data(),hraw->GetName()));
                        }

                        //Adopt
                        Bool_t adoptOK = OC()->Adopt(Form("%s",id.Data()),hraw);

                        if ( adoptOK ) std::cout << "+++raw histo " << hraw->GetName() << " adopted" << std::endl;
                        else AliError(Form("Could not adopt Yield histo %s",hraw->GetName()));
                        new TCanvas;
                        hraw->DrawCopy("e0");
                    }
                    delete binning;
                }
            }
        }
    }
    return;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputePPCrossSection(const char* spectraName,const char* externfile, const char* externfile2, Bool_t print, const char* what)
{
/**
 *   [AliAnalysisMuMu::ComputePPCrossSection description]
 *   @brief   Compute the PP Cross section. At the moment, only implemented when <what> == CorrNofJPsi
 *   @details Delegate the process to AliAnalysisMuMuCapsulePP.
 *
 *   @param   spectraName spectra Name
 *   @param   what        Quantity stored in the AliAnalysisMuMuSpectra used for the cross-section. Should always be already Accxeff corrected ("CorrNofJPsi" for instance)
 *   @param   externfile  Config. file readed by AliAnalysisMuMuSpectraCapsule classes. See AliAnalysisMuMuSpectraCapsule documentation
 *   @param   externfile2 Config. file readed by AliAnalysisMuMuSpectraCapsule classes. See AliAnalysisMuMuSpectraCapsule documentation
 *   @param   print       Print more details
 */

  if(!OC()){
    AliError("No mergeable. Consider Upgrade()");
    return;
  } else {
    cout <<      " ================================================================ " << endl;
    cout <<      "                       ComputePPCrossSection                          " << endl;
    cout <<      " ================================================================ " << endl;
  }

   // Get configuration settings
  TObjArray* eventTypeArray   = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
  TObjArray* triggerArray     = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
  TObjArray* fitfunctionArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());// to add here an entry
  TObjArray* pairCutArray     = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
  TObjArray* centralityArray  = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
  TObjArray* bins;

  // Iterator for loops
  TIter nextTrigger(triggerArray);
  TIter nextEventType(eventTypeArray);
  TIter nextPairCut(pairCutArray);
  TIter nextCentrality(centralityArray);

  // Strings
  TObjString* striggerDimuon;
  TObjString* seventType;
  TObjString* spairCut;
  TObjString* scentrality;

  // Pointers
  TGraphErrors* graph=0x0;
  TGraphErrors* graphErr=0x0;

  TList* list=0x0;
  AliAnalysisMuMuSpectra spectra=0x0;


  // Loop on each envenType (see MuMuConfig)
  while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
  {
    AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
    nextTrigger.Reset();

    // Loop on each trigger (see MuMuConfig)
    while ( ( striggerDimuon = static_cast<TObjString*>(nextTrigger())) )
    {
      AliDebug(1,Form("-TRIGGER %s",striggerDimuon->String().Data()));
      nextPairCut.Reset();

      // Loop on each paircut (not the ones in MuMuConfig but the ones set)
      while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
      {
        AliDebug(1,Form("--PAIRCUT %s",spairCut->String().Data()));
        nextCentrality.Reset();

        // Loop on each centrality (not the ones in MuMuConfig but the ones set)
        while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
        {
          AliDebug(1,Form("----CENTRALITY %s",scentrality->String().Data()));
          // The binning pointer, which point at Pt binning, Y binning etc.


          //________Get spectra
          TString spectraPath= Form("/FitResults/%s/%s/%s/%s/%s",seventType->String().Data(),striggerDimuon->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName);

          AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));
          if(!spectra){
            AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
            return;
          }

          // Create capsule who will compute the cross section
          AliAnalysisMuMuSpectraCapsulePP * capsule = new AliAnalysisMuMuSpectraCapsulePP(spectra,spectraPath,externfile,externfile2);
          if(!capsule) continue;
          AliDebug(1,Form("Spectra = %p",capsule));

          if(print)capsule->SetPrintFlag();

          // Get Graph with RAA results
          list = capsule->ComputeJpsiPPCrossSection(what);

          AliDebug(1,Form("list = %p",list));
          if(!list) continue;

          graph = static_cast<TGraphErrors*>(list->At(0)->Clone());
          graphErr = static_cast<TGraphErrors*>(list->At(1)->Clone());

          // --- playground for a beautyfull plot ! ---

          Config()->LoadAliceStyles();
          //points
          graph->SetMarkerColor(2);
          graph->SetMarkerStyle(8);
          graph->SetLineColor(2);
          graph->SetLineWidth(2);
          graph->SetMarkerSize(1.5);
          //syst.
          graphErr->SetMarkerColor(2);
          graphErr->SetFillColor(0);
          graphErr->SetLineColor(2);
          graphErr->SetLineWidth(2);
          graphErr->SetTitle(Form("J/#psi cross section"));

          TCanvas *c = new TCanvas("c","c",4,132,1024,768);
          c->SetFillColor(0);
          c->SetBorderMode(0);
          c->SetBorderSize(2);
          c->SetLeftMargin(0.15);
          c->SetRightMargin(0.03);
          c->SetTopMargin(0.03);
          c->SetBottomMargin(0.13);
          c->SetFrameBorderMode(0);

          TLegend * leg = new TLegend(0.2,0.2,0.50,0.4);
          leg->SetHeader(Form("ALICE, pp #sqrt{#it{s}} = 5.02 TeV"));

          TString sspectraName(spectraName);

          if(sspectraName.Contains("PT") && !sspectraName.Contains("VS"))
          {
            graphErr->GetXaxis()->SetTitle(Form("#it{p}_{T}(GeV/#it{c})"));
            graphErr->GetYaxis()->SetTitle(Form("d^{2}#sigma/d#it{y}d#it{p}_{T} (#mub(GeV/#it{c})^{-1})"));
            graphErr->GetYaxis()->SetTitleSize(0.06);
            graphErr->GetYaxis()->SetRangeUser(10e-4,2);
            gPad->SetLogy();
            leg->AddEntry(graph,"Inclusive J/#psi , 2.5 < #it{y} < 4 ","pe");

          }
          else if(sspectraName.Contains("Y") && !sspectraName.Contains("VS"))
          {
            graphErr->GetXaxis()->SetTitle(Form(" #it{y}"));
            graphErr->GetXaxis()->SetRangeUser(2.49,4.01);
            graphErr->GetYaxis()->SetTitle(Form("d^{2}#sigma/d#it{y}(#muB)"));
            graphErr->GetYaxis()->SetTitleSize(0.06);
            leg->AddEntry(graph,"Inclusive J/#psi , 0 < #it{p}_{T} < 12 GeV/#it{c}","pe");
          }
          else if(sspectraName.Contains("PTVSY-INT_PT_CUT"))
          {
            graphErr->GetXaxis()->SetTitle(Form(" #it{y}"));
            graphErr->GetXaxis()->SetRangeUser(2.49,4.01);
            graphErr->GetYaxis()->SetTitle(Form("d^{2}#sigma/d#it{y}(#muB)"));
            graphErr->GetYaxis()->SetTitleSize(0.06);
            leg->AddEntry(graph,"Inclusive J/#psi , 0.3 < #it{p}_{T} < 12 GeV/#it{c}","pe");
          }
          else if(sspectraName.Contains("PTVSY-FULL_PT_CUT"))
          {
            graphErr->GetXaxis()->SetTitle(Form(" #it{y}"));
            graphErr->GetXaxis()->SetRangeUser(2.49,4.01);
            graphErr->GetYaxis()->SetTitle(Form("d^{2}#sigma/d#it{y}(#muB)"));
            graphErr->GetYaxis()->SetTitleSize(0.06);
            leg->AddEntry(graph,"Inclusive J/#psi , 0.3 < #it{p}_{T} < 12 GeV/#it{c}","pe");
          }
          else if(sspectraName.Contains("YVSPT-2DBIN1"))
          {
            graphErr->GetXaxis()->SetTitle(Form("#it{p}_{T}(GeV/#it{c})"));
            graphErr->GetYaxis()->SetTitle(Form("d^{2}#sigma/d#it{y}d#it{p}_{T} (#mub(GeV/#it{c})^{-1})"));
            graphErr->GetYaxis()->SetTitleSize(0.06);
            graphErr->GetYaxis()->SetRangeUser(10e-4,2);
            gPad->SetLogy();
            leg->AddEntry(graph,"Inclusive J/#psi ,3.25 < #it{y} < 4","pe");
          }
          else if(sspectraName.Contains("YVSPT-2DBIN2"))
          {
            graphErr->GetXaxis()->SetTitle(Form("#it{p}_{T}(GeV/#it{c})"));
            graphErr->GetYaxis()->SetTitle(Form("d^{2}#sigma/d#it{y}d#it{p}_{T} (#mub(GeV/#it{c})^{-1})"));
            graphErr->GetYaxis()->SetTitleSize(0.06);
            graphErr->GetYaxis()->SetRangeUser(10e-4,2);
            gPad->SetLogy();
            leg->AddEntry(graph,"Inclusive J/#psi ,2.5 < #it{y} < 3.25","pe");
          }
          graphErr->SetMarkerSize(1.7);

          leg->AddEntry((TObject*)0,"L_{int} = 106.3 #pm 2.1% nb^{-1}","");
          leg->SetTextSize(0.04);
          TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
          header->SetTextSize(.055);

          graphErr->DrawClone("A5.[]");

          // --- Playground for a fit ---

          // TF1* fit=0x0;
          // TFitResultPtr Fitpoint=0x0;
          // Pt distribution
          // fit = new TF1("fit","[0]*x/TMath::Power([1] + TMath::Power(x,[2]),[3])",0,12);
          // fit->SetParameters(5.42*4654.3/1.5,12.8133,1.9647,3.66641);
          // fit->DrawCopy("same");
          // Y distribution
          // fit= new TF1("fit","[0] * TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2.))", -10, 2.5);
          // fit->SetParameters(2.8,0.00,2.9271);
          // fit->FixParameter(1,0.00);

          // if(fit) Fitpoint = graph->Fit("fit","WMR");
          // fit->DrawCopy("same");
          // if(sspectraName.Contains("PT") && !sspectraName.Contains("YVSPT") && static_cast<int>(Fitpoint)==0) printf("cross section = %f +/- %f \n",fit->Integral(0.,8.),fit->IntegralError(0.,8.));
          // if(sspectraName.Contains("Y")  && static_cast<int>(Fitpoint)==0)  printf("cross section = %f +/- %f \n",fit->Integral(-4.,-2.5),fit->IntegralError(-4.,-2.5));

          graph->DrawClone("Psame");
          leg->DrawClone("same");
          gPad->RedrawAxis();
        }
      }
    }
  }
  return;
}

//_____________________________________________________________________________
TH2* AliAnalysisMuMu::ComputeSPDCorrection(const char* type, const char* eventSel, const char* triggerSel, Bool_t bkgReject)
{
/**
 *   [AliAnalysisMuMu::ComputeSPDCorrection description]
 *   @brief   Old method with a lot of hard coded lines, not sure of what it does... Needs to be more general
 *   @details make it general
 *
 *   @param   type       [description]
 *   @param   eventSel   [description]
 *   @param   triggerSel [description]
 *   @param   bkgReject  [description]
 *   @return             [description]
 */
    TString stype(type);
    TString evtype(eventSel);
    TString trigtype(triggerSel);
    // >Add centrality and title in the path<
    TH2* hNch = static_cast<TH2*>(OC()->Histo(Form("/MCINPUT/%s/%s/V0A/NchVsZVertexVsEta",evtype.Data(),
                                                   trigtype.Data()))); // Input Nch // //"/INPUT/QASPDZPSALL/NchVSEtaVSZVertMC"
    if ( !hNch )
        {
        AliError("No Nch histo found");
        return 0x0;
        }
    // >Add centrality and title in the path<
    TH2* hNtr = static_cast<TH2*>(OC()->Histo(Form("/%s/%s/V0A/TrackletsVsZVertexVsEta",evtype.Data(),
                                                   trigtype.Data()))); // Reco tracklets //  //"/RECO/QASPDZPSALL/MB1/NtrVSEtaVSZVertMC"
    if ( !hNtr )
        {
        AliError("No tracklets histo found");
        return 0x0;
        }
    // >Add centrality and title in the path<
    TH2* hNtrBkg = static_cast<TH2*>(OC()->Histo(Form("/MCINPUT/%s/%s/V0A/NBkgTrackletsVsZVertexVsEta",evtype.Data(),
                                                      trigtype.Data()))); // Reco tracklets //  //"/RECO/QASPDZPSALL/MB1/NtrVSEtaVSZVertMC"
    if ( !hNtrBkg )
        {
        AliError("No background tracklets histo found");
        return 0x0;
        }


    TH2D* hSPDCorr = static_cast<TH2D*>(hNtr->Clone("SPDCorr"));
    TString title("");\
    if ( stype.Contains("oneOverAccEff")) hSPDCorr->SetTitle("SPD 1/AccxEff correction");
    else if ( stype.Contains("AccEffOnly")) hSPDCorr->SetTitle("SPD AccxEff correction");
    else if ( stype.Contains("statOneOverAccEff")) hSPDCorr->SetTitle("SPD 1/AccxEff correction stat. unc.");

    for (Int_t i = 1 ; i < hNch->GetNbinsX() ; i++)
        {
        for (Int_t j = 1 ; j < hNch->GetNbinsY() ; j++)
            {
            Int_t n = hNch->GetBin(i,j);
            Double_t nch = hNch->GetBinContent(n);
            Double_t ntr = hNtr->GetBinContent(n);
            Double_t nBkgtr(0.);
            if ( bkgReject ) nBkgtr = hNtrBkg->GetBinContent(n);

            Double_t corr(0.),corrErr(0.);
            if ( nch != 0. )
                {
                corr = (ntr - nBkgtr)/nch;
                corrErr = TMath::Max( 1./nch,TMath::Sqrt( corr*(1.-corr)/nch ) );
                }

            if ( stype.Contains("oneOverAccEff"))
                {
                if ( corr > 0. )
                    {
                    hSPDCorr->SetBinContent(n,1./corr);
                    hSPDCorr->SetBinError(n,corrErr/TMath::Power(corr,2.));
                    }
                else
                    {
                    hSPDCorr->SetBinContent(n,0.);
                    hSPDCorr->SetBinError(n,1.);
                    }

                }
            else if ( stype.Contains("AccEffOnly"))
                {
                hSPDCorr->SetBinContent(n,corr);
                hSPDCorr->SetBinError(n,corrErr);
                }
            else if ( stype.Contains("statOneOverAccEff"))
                {
                if ( corr != 0. )
                    {
                    hSPDCorr->SetBinContent(n,(corrErr/TMath::Power(corr,2.)*100)/(1./corr));
                    }
                else
                    {
                    hSPDCorr->SetBinContent(n,-1);
                    }

                }
            }
        }

    return hSPDCorr;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeNumberOfEvent()
{
/**
 *   [AliAnalysisMuMu::ComputeNumberOfEvent description]
 *   @brief   Compute the CMUL to CINT ratio(s)
 *
 */
    if (!OC() || !CC()|| !Config())
    {
        AliError("No mergeable/counter collection. Consider Upgrade()");
        return ;
    }
    else
    {
        cout <<      " ================================================================ " << endl;
        cout <<      "                     ComputeNumberOfEvent                 " << endl;
        cout <<      " ================================================================ " << endl;
    }

    const AliAnalysisMuMuFnorm::ETriggerType triggerTypes[] = { AliAnalysisMuMuFnorm::kMB, AliAnalysisMuMuFnorm::kMUL, AliAnalysisMuMuFnorm::kMSL, AliAnalysisMuMuFnorm::kMSH };
    const Bool_t trueFalse[] = { kTRUE, kFALSE };

    //Delete precedent Fnorm
    OC()->Prune("/FNORM/Nevent");

    AliAnalysisMuMuFnorm computer(*(CC()),*(Config()),AliAnalysisMuMuFnorm::kMUL,Config()->OCDBPath(),Config()->CompactGraphs());// here trigger type doesn't count

    for ( Int_t i = 0; i < 4; ++i )
    {
        for ( Int_t pileup = 0; pileup < 2; ++pileup )
        {
          for ( Int_t ps = 0; ps < 2; ++ps )
          {
            computer.ComputeNofEvents(triggerTypes[i],trueFalse[pileup],ps);
          }
        }
    }

    AliMergeableCollection* fnorm = computer.DetachMC();
    //Set the new ownership
    OC()->Attach(fnorm,"/FNORM/Nevent");
    //Update
    Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeFnormWeightedMeanGraphs(AliAnalysisMuMuFnorm::ETriggerType refTrigger,const char* patternOrList, const char* graphName )
{
/**
 *   [AliAnalysisMuMu::ComputeFnormWeightedMeanGraphs description]
 *   @brief   Compute the FNorm Weighted Mean Graph.
 *   @details Delegate to AliAnalysisMuMuFNorm.
 *
 *   @param   refTrigger    kMUL, kMB ... See AliAnalysisMuMuFNorm
 *   @param   patternOrList Name of the graphs included in the mean process. Two mean graphs are separated by ":" and graphics used to compute a mean graphs are
 *                          separated by ','. Example : "FnormOffline2PUPS,FnormOffline1PUPS:FnormOffline2PUPS,FnormOffline1PUPS,FnormScalers1PUPS" computes
 *                          mean graph for "FnormOffline2PUPS - FnormOffline1PUPS", another one for "FnormOffline2PUPS - FnormOffline1PUPS - FnormScalers1PUPS".
 *   @param   graphName     Name of the output mean graph.
 */

    if (!OC() || !CC()|| !Config())
    {
      AliError("No mergeable/counter collection. Consider Upgrade()");
      return ;
    } else {
      cout <<      " ================================================================ " << endl;
      cout <<      "                     ComputeFnormWeightedMeanGraphs                 " << endl;
      cout <<      " ================================================================ " << endl;
    }

    //Delete precedent Fnorm
    OC()->Prune("/FNORM/WeightedMeanGraphs");

    AliAnalysisMuMuFnorm computer(*(CC()),*(Config()),refTrigger,Config()->OCDBPath(),Config()->CompactGraphs());

    TString spattern(patternOrList);
    TObjArray* slist(0x0);
    slist = spattern.Tokenize(":");

    TIter next(slist);
    TObjString* str(0x0);

    while ( ( str = static_cast<TObjString*>(next()) ) ) computer.WeightedMeanGraphs(str->String().Data(),graphName,OC());

    AliMergeableCollection* fnorm = computer.DetachMC();
    //Set the new ownership
    OC()->Attach(fnorm,"/FNORM/WeightedMeanGraphs");
    //Update
    Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeFnormScalers(AliAnalysisMuMuFnorm::ETriggerType refTrigger, Bool_t PileUpCorr)
{
/**
 *   [AliAnalysisMuMu::ComputeFnormScalers description]
 *   @brief   Compute the MB to REF ratio using the scalers method (from OCDB)
 *   @details Delegate to AliAnalysisMuMuFnorm.
 *
 *   @param   refTrigger [description]
 *   @param   PileUpCorr [description]
 */

    if (!OC() || !CC()|| !Config()){
      AliError("No mergeable/counter collection. Consider Upgrade()");
      return ;
    } else {
      cout <<      " ================================================================ " << endl;
      cout <<      "                     ComputeFnormScalers                 " << endl;
      cout <<      " ================================================================ " << endl;
    }

    //Delete precedent Fnorm
    OC()->Prune("/FNORM/Scaler");

    AliAnalysisMuMuFnorm computer(*(CC()),*(Config()),refTrigger,Config()->OCDBPath(),Config()->CompactGraphs());

    computer.ComputeFnormScalers(kFALSE,0);
    if(PileUpCorr){
        computer.ComputeFnormScalers(kTRUE,0);
        computer.ComputeFnormScalers(kTRUE,1);
    }

    AliMergeableCollection* fnorm = computer.DetachMC();
    //Set the new ownership
    OC()->Attach(fnorm,"/FNORM/Scaler");
    //Update
    Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::ComputeIntFnormFromCounters(AliAnalysisMuMuFnorm::ETriggerType refTrigger, Bool_t PileUpCorr)
{
 /// Compute the CMUL to CINT ratio(s) in 1 or 2 steps (Offline method) from the CC(), in bins.
 ///
 /// Important considerations:
 ///
 ///   - The analysed file must contain the CMUL, CINT, CMSL, CMSL&0MUL and CINT&0MSL triggers to work correctly.
 ///
 ///    FIX ME : quantity not define in CC() for "pt" and "y".


    if (!OC() || !CC()|| !Config()){
      AliError("No mergeable/counter collection. Consider Upgrade()");
      return ;
    } else {
        cout <<      " ================================================================ " << endl;
        cout <<      "                     ComputeIntFnormFromCounters                 " << endl;
        cout <<      " ================================================================ " << endl;
    }

    //Delete precedent Fnorm
    OC()->Prune("/FNORM/Offline");

    AliAnalysisMuMuFnorm computer(*(CC()),*(Config()),refTrigger,Config()->OCDBPath(),Config()->CompactGraphs());

    computer.ComputeFnormOffline(1,kFALSE,0); // CINT/CINT&0MUL
    computer.ComputeFnormOffline(1,kFALSE,1); // CINT/CINT&0MUL + Event selection corrected
    if(PileUpCorr) computer.ComputeFnormOffline(1,kTRUE,1); // CINT/CINT&0MUL + Event selection corrected + pileup correction

    computer.ComputeFnormOffline(2,kFALSE,0); // CMSL/CMSL&0MUL x CINT/CINT&0MSL
    computer.ComputeFnormOffline(2,kFALSE,1); // CMSL/CMSL&0MUL x CINT/CINT&0MSL + Event selection corrected
    if(PileUpCorr) computer.ComputeFnormOffline(2,kTRUE,1); // CMSL/CMSL&0MUL x CINT/CINT&0MSL+ Event selection corrected + pileup correction

    AliMergeableCollection* fnorm = computer.DetachMC();
    //Set the new ownership
    OC()->Attach(fnorm,"/FNORM/Offline/");
    Update();
}

//_____________________________________________________________________________
void AliAnalysisMuMu::PlotJpsiYield(const char* spectraName, const char* subresultname, const char* beamYear, int NofMUL, const char* externfile, const char* externfile2 )
{
/**
 *   [AliAnalysisMuMu::ComputeJpsiYield description]
 *   @brief   Compute the Jpsi yield, i.e NofJPsi/AccxEff/FNorm/MUL.
 *   @details Delegate to AliAnalysisMuMuCapsule. It can be compute for an specific subresult (fit with an specific background shape, signal, fitting range... combination) or from the mean of all the subresults.
 *
 *   @param binType        [description]
 *   @param externfile     Config. file readed by AliAnalysisMuMuSpectraCapsule classes. See AliAnalysisMuMuSpectraCapsule documentation
 *   @param externfile2    Config. file readed by AliAnalysisMuMuSpectraCapsule classes. See AliAnalysisMuMuSpectraCapsule documentation
 *   @param AccEffCorr     Spectra type
 *   @param subresultname  If one wants to draw one/several specific results.
 */

    if (!OC() || !CC()) {
      AliError("No mergeable/counter collection. Consider Upgrade()");
      return ;
    } else {
      cout <<      " ================================================================ " << endl;
      cout <<      "                       ComputeJpsiYield                           " << endl;
      cout <<      " ================================================================ " << endl;
    }

    // Get configuration settings
    TObjArray* eventTypeArray   = Config()->GetListElements(Config()->EventSelectionKey(),IsSimulation());
    TObjArray* triggerArray     = Config()->GetListElements(Config()->DimuonTriggerKey(),IsSimulation());
    TObjArray* fitfunctionArray = Config()->GetListElements(Config()->FitTypeKey(),IsSimulation());// to add here an entry
    TObjArray* pairCutArray     = Config()->GetListElements(Config()->PairSelectionKey(),IsSimulation());
    TObjArray* centralityArray  = Config()->GetListElements(Config()->CentralitySelectionKey(),IsSimulation());
    TObjArray* bins;


    // Iterator for loops
    TIter nextTrigger(triggerArray);
    TIter nextEventType(eventTypeArray);
    TIter nextPairCut(pairCutArray);
    TIter nextCentrality(centralityArray);

    // Strings
    TObjString* striggerDimuon;
    TObjString* seventType;
    TObjString* spairCut;
    TObjString* scentrality;

    TString striggerMB  = Config()->First(Config()->MinbiasTriggerKey(),IsSimulation());
    const TString syear(beamYear);

    // Pointers
    TH1* h                         = 0x0;
    TGraphErrors* graph            = 0x0;
    TGraphErrors* graphCent        = 0x0;
    AliAnalysisMuMuSpectra spectra = 0x0;

    // Loop on each envenType (see MuMuConfig)
    while ( ( seventType = static_cast<TObjString*>(nextEventType())) )
    {
      AliDebug(1,Form("EVENTTYPE %s",seventType->String().Data()));
      nextTrigger.Reset();
      // Loop on each trigger (see MuMuConfig)
      while ( ( striggerDimuon = static_cast<TObjString*>(nextTrigger())) )
      {
        AliDebug(1,Form("-TRIGGER %s",striggerDimuon->String().Data()));
        nextPairCut.Reset();

        // Loop on each paircut (not the ones in MuMuConfig but the ones set)
        while ( ( spairCut = static_cast<TObjString*>(nextPairCut())) )
        {
          AliDebug(1,Form("--PAIRCUT %s",spairCut->String().Data()));
          //canvas
          TCanvas *c1 = new TCanvas;
          c1->Draw();
          gStyle->SetOptStat(0);

          nextCentrality.Reset();

          // Loop on each centrality (not the ones in MuMuConfig but the ones set)
          while ( ( scentrality = static_cast<TObjString*>(nextCentrality()) ) )
          {
            AliDebug(1,Form("---CENTRALITY %s",scentrality->String().Data()));

            //________Get Fnorm Histo
            TString id(Form("/FNORM-%s/%s/%s/%s",striggerDimuon->String().Data(),seventType->String().Data(),scentrality->String().Data(),syear.Data())); // Path to save the Fnorm and EqNofMB histos in the mergeable collection

            TString idHisto="";
            TString sspectraName(spectraName);
            if (!sspectraName.Contains("INTEGRATED")) idHisto= Form("hNofEqMBVS%s",spectraName);
            else                                      idHisto= Form("hFNormInt_%s",striggerMB.Data());

            h = OC()->Histo(Form("%s/%s",id.Data(),idHisto.Data()));
            if (!h) AliError(Form("Could not find histo in %s/%s. Take FNorm from the capsule",id.Data(),idHisto.Data()));
            //________

            //________Get spectra
            TString spectraPath= Form("/FitResults/%s/%s/%s/%s/%s",seventType->String().Data(),striggerDimuon->String().Data(),scentrality->String().Data(),spairCut->String().Data(),spectraName);
            AliAnalysisMuMuSpectra * spectra = static_cast<AliAnalysisMuMuSpectra*>(OC()->GetObject(spectraPath.Data()));
            if(!spectra)
            {
              AliError(Form("Cannot find spectra with name %s",spectraPath.Data()));
              continue;
            }
            //________

            //________Select methods
            if(syear.Contains("PbPb"))
            {
              printf("there !!\n");
              AliAnalysisMuMuSpectraCapsulePbPb * capsule = new AliAnalysisMuMuSpectraCapsulePbPb(spectra,spectraPath,externfile,externfile2);
              AliDebug(1,Form("Spectra = %p",capsule));

              const char* what ="NofJPsi";

              // Get Graph with Yield results
              printf("coucou !! !!\n");
              graph = capsule->ComputeYield(what,h,subresultname,NofMUL);

              TLegend * leg = new TLegend(0.4,0.7,0.90,0.9);
              leg->SetHeader(Form("ALICE, Pb-Pb #sqrt{s_{NN}}=2.76 TeV, L_{int}=70 #mub^{-1}, %s",scentrality->String().Data()));
              leg->AddEntry(graph,"Inclusive J/#psi Yield","pe");
              graph->Draw("ap");
              leg->Draw();

              delete capsule;
            }
            else if(syear.Contains("pPb") || syear.Contains("Pbp"))
            {
              AliError("Not implemented yet ! You are welcome to do so :D ");
              return;
              // AliAnalysisMuMuSpectraCapsulePbP * capsule = new AliAnalysisMuMuSpectraCapsulePbP(spectra,spectraPath,externfile,externfile2);
              // AliDebug(1,Form("Spectra = %p",capsule));

              // Int_t NofMUL= TMath::Nint(CC()->GetSum(Form("trigger:%s/event:%s/centrality:%s",striggerDimuon->String().Data(),seventType->String().Data(),"V0M_00.00_90.00")));
              // AliDebug(1,Form("Reference centrality for NofMUL = V0M_00.00_90.00"));

              // // Get Graph with Yield results
              // graph = capsule->ComputeYield(what,h,subresultname,NofMUL);

              // delete capsule;
            }
            else if(syear.Contains("pp"))
            {
              const char* what ="CorrNofJPsi";

              AliAnalysisMuMuSpectraCapsulePP * capsule = new AliAnalysisMuMuSpectraCapsulePP(spectra,spectraPath,externfile,externfile2);
              AliDebug(1,Form("Spectra = %p",capsule));

              Int_t NofMUL= TMath::Nint(CC()->GetSum(Form("trigger:%s/event:%s/centrality:%s",striggerDimuon->String().Data(),seventType->String().Data(),"PP")));
              AliDebug(1,Form("Reference centrality for NofMUL = PP"));

              // Get Graph with Yield results
              graph = capsule->ComputeYield(what,h,subresultname,NofMUL);

              TLegend * leg = new TLegend(0.4,0.7,0.90,0.9);
              leg->SetHeader(Form("ALICE, pp #sqrt{s}=5.02 TeV, L_{int}=116.3 #mub^{-1}, %s",scentrality->String().Data()));
              leg->AddEntry(graph,"Inclusive J/#psi Yield","pe");
              graph->Draw("ap");
              leg->Draw();

              delete capsule;
            }

            //________ Update resultes in Mergeable collection
            TString id2(Form("/JpsiYield-%s/%s/%s",striggerDimuon->String().Data(),seventType->String().Data(),spairCut->String().Data()));

            TObject* o = 0x0;

            if (graph)// first graph
            {
              o = fMergeableCollection->GetObject(Form("%s/%s",id2.Data(),graph->GetName()));

              if (o)
              {
                AliWarning(Form("Replacing %s/%s",id2.Data(),graph->GetName()));
                fMergeableCollection->Remove(Form("%s/%s",id2.Data(),graph->GetName()));
              }

              Bool_t adoptOK = fMergeableCollection->Adopt(id2.Data(),graph);

              if ( adoptOK ) std::cout << "+++JpsiYield graph " << graph->GetName() << " adopted" << std::endl;
              else AliError(Form("Could not adopt JpsiYield grap %s",graph->GetName()));
            }
            //________
          }
        }
      }
    }

    return;
}

//_____________________________________________________________________________

void AliAnalysisMuMu::ComputeMBXSectionFractionInBins(const char* filePileUpCorr, const char* eventSelection, const char* what,const char* quantity
                                                      ,const char* flavour)
{
  ////////////////
  // OLD METHOD //
  ////////////////
/**
 *   [AliAnalysisMuMu::ComputeMBXSectionFractionInBins description]
 *   @brief   Compute the CMUL to CINT ratio(s) in 2 steps (Offline method) from the CC(), in bins.
 *   @details Important considerations:
 *
 *   - The analysed file must contain the CMUL, CINT, CMSL, CMSL&0MUL and CINT&0MSL triggers to work correctly.
 *
 *
 *   - The analysed file must contain the event selection PSALL and the one used in the yield analysis (i.e. PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00) to work correctly. (the first to compute the Fnorm and the second to get the NofCMUL used in the yield analysis to get the correct NofEqMB = Fnorm*NofCMUL)
 *
 *   OLD METHOD WITH A LOT OF HARD CODED LINES, MAY BE OUTDATED
 *
 *   @param   filePileUpCorr txt file with the pile up correction run by run. Each line in the file must have the format : RUN 195681 PERIOD LHC13d PILE-UP CORRECTION FACTOR (mu/(1-exp(-mu)) =  1.0015
 *   @param   eventSelection desired event selection
 *   @param   what           what the binning range is about (J/psi, event...). By default is "psi".
 *   @param   quantity       binning type.
 *   @param   flavour        binning flavour.
 */
    //________Decoding of the pileup correction file
    Bool_t corrPU(kFALSE);
    TObjArray* pUCorr = new TObjArray();
    if ( strlen(filePileUpCorr) > 0 )
        {
        std::cout << "Extracting Pile-Up correction factors from " << filePileUpCorr << std::endl;
        char line[1024];
        ifstream in(filePileUpCorr);

        while ( in.getline(line,1024,'\n'))
            {
            TString lrun(line);
            TString lvalue(line);

            lrun.Remove(0,4);
            lrun.Remove(6,67);

            lvalue.Remove(0,lvalue.First("=")+1);

            std::cout << "RUN: " << lrun.Data() << " PUFactor = " << lvalue.Data() << std::endl;

            pUCorr->Add(new TParameter<Double_t>(lrun.Data(),lvalue.Atof()));
            }
        corrPU = kTRUE;
        }
    //________


    TString seventSelection(eventSelection);
    TString sQuantity(quantity);
    TString sruns = CC()->GetKeyWords("run");
    TObjArray* runs = sruns.Tokenize(",");

    TIter nextRun(runs);
    TObjString* s;

    AliAnalysisMuMuBinning* binning = BIN()->Project(what,sQuantity.Data(),flavour);
    if ( !binning )
        {
        AliError(Form("%s-%s-%s binning does not exist",what,sQuantity.Data(),flavour));
        return;
        }
    TObjArray* bin = binning->CreateBinObjArray(what,sQuantity.Data(),flavour);
    Int_t nEntries = bin->GetEntries();

    Double_t* nCINTBin = new Double_t[nEntries];
    Double_t nCINTTot = 0. ;

    TIter nextBin(bin);
    AliAnalysisMuMuBinning::Range* r;
    Int_t i(0); // Bin number
    while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ) //Bin loop
        {
        nCINTBin[i] = 0;

        std::cout << "______________________________" << std::endl;
        std::cout << "Bin: " << r->AsString().Data() << std::endl;

        nextRun.Reset();
        while ( ( s = static_cast<TObjString*>(nextRun())) ) //Run loop
            {
            Double_t nCINT = CC()->GetSum(Form("/event:%s/trigger:CINT7-B-NOPF-ALLNOTRD/centrality:V0A/run:%s/bin:%s",
                                               seventSelection.Data(),s->GetName(),r->AsString().Data()));
            Double_t pUfactor = 1.;
            if ( nCINT !=0. )
                {
                if (corrPU)
                    {
                    TParameter<Double_t>* p = static_cast<TParameter<Double_t>*>(pUCorr->FindObject(s->GetName()));
                    if ( p ) pUfactor = p->GetVal();
                    else
                        {
                        AliError(Form("Run %s not found in pile-up correction list",s->GetName()));
                        }
                    }

                nCINT = nCINT*pUfactor;
                }
            else
                {
                std::cout << " Warning: Run " << s->GetName() << " has no MB trigger in this bin" << std::endl;
                continue;
                }
            nCINTBin[i] += nCINT;

            std::cout << "Run " << s->GetName() <<  " ; " << "Nof CINT = " << nCINT << std::endl;

            }

        nCINTTot += nCINTBin[i];

        std::cout << std::endl;

        std::cout << "Nof CINT in bin = " << nCINTBin[i]  << std::endl;

        i++;
        }


    nextBin.Reset();
    i=0;
    Double_t totXSectFract(0.);
    while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
        {
        Double_t xSecFract = (nCINTBin[i] / nCINTTot)*100.;

        totXSectFract += xSecFract;

        std::cout << "Cross section Fraction in bin " << r->AsString().Data() << " = " << xSecFract  << std::endl;

        i++;
        }

    std::cout << std::endl;
    std::cout << "Total xSection in bins = " << totXSectFract  << std::endl;

    delete binning;
    delete runs;
    delete bin;
    delete[] nCINTBin;

    return;
}

//_____________________________________________________________________________
Double_t AliAnalysisMuMu::ErrorPropagationAxBoverCxD(Double_t a,Double_t b,Double_t c,Double_t d)
{
    //Just valid for counts
    Double_t error2 = TMath::Power(b/(c*d),2.)*a + TMath::Power(a/(c*d),2.)*b + TMath::Power(a*b*d,2.)*(c/TMath::Power(c*d,4.)) + TMath::Power(a*b*c,2.)*(d/TMath::Power(c*d,4.));

    return TMath::Sqrt(error2);
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::CorrectSpectra(const char* binType, const char* flavour,const char* accEffSubResultName)
{
/**
 *   [AliAnalysisMuMu::CorrectSpectra description]
 *   @brief   Correct an AliAnalysisMuMuSpectra
 *   @details Compute AccxEff corrected yield of a 'real' spectra with a 'mc' spectra. Delegate to AliAnalysisMuMuSpectra.
 *
 *   Spectra selection made with first entries of config. file keys.
 *
 *   DimuonTriggerKey() used to select trigger.
 *
 *   @param   binType             [description]
 *   @param   flavour             [description]
 *   @param   accEffSubResultName Results from which we take the AccxEff values
 *   @return  The 'real' spectra with corrected yield values stored in it.
 */
    if (!SIM()){
      AliError("Cannot compute corrected yield without associated MC file !");
      return 0x0;
    }
    else{
      cout <<      " ================================================================ " << endl;
      cout <<      "                       CorrectSpectra                           " << endl;
      cout <<      " ================================================================ " << endl;
    }

    TString Event         =Config()->First(Config()->EventSelectionKey(),kTRUE);
    TString Centrality    =Config()->First(Config()->CentralitySelectionKey(),kTRUE);
    TString DimuonTrigger =Config()->First(Config()->DimuonTriggerKey(),kTRUE);
    TString PairSelection =Config()->First(Config()->PairSelectionKey(),kTRUE);

    AliAnalysisMuMuSpectra* realSpectra = GetSpectraFromConfig(binType,flavour);
    AliAnalysisMuMuSpectra* simSpectra = SIM()->GetSpectra(binType,Event.Data(),DimuonTrigger.Data(),Centrality.Data(),PairSelection.Data(),flavour);

    if ( !realSpectra ){
      AliError("could not get real spectra");
      return 0x0;
    }

    if ( !simSpectra){
      AliError("could not get sim spectra");
      return 0x0;
    }

    realSpectra->Correct(*simSpectra,"JPsi",accEffSubResultName);

    Update();

    return realSpectra;
}

//_____________________________________________________________________________
AliAnalysisMuMuSpectra* AliAnalysisMuMu::ComputeYield(const char* binType, const char* flavour,const char* accEffSubResultName)
{
/**
 *   [AliAnalysisMuMu::ComputeYield description]
 *   @brief   Compute yield from simulation file
 *   @details Delegate to AliAnalysisMuMuSpectra. Compute yield using a 'real' spectra and a 'mc' spectra. Spectra are selected using first key entries of the config. file.
 *
 *   FIX ME: Make it general.
 *
 *   @param   binType             [description]
 *   @param   flavour             [description]
 *   @param   accEffSubResultName Results from which we take the AccxEff values
 *   @return                      The 'real' spectra with a corrected yield entry
 */

    TString Event         =Config()->First(Config()->EventSelectionKey(),kTRUE);
    TString Centrality    =Config()->First(Config()->CentralitySelectionKey(),kTRUE);
    TString DimuonTrigger =Config()->First(Config()->DimuonTriggerKey(),kTRUE);
    TString PairSelection =Config()->First(Config()->PairSelectionKey(),kTRUE);

    if (!SIM()){
      AliError("Cannot compute corrected yield without associated MC file !");
      return 0x0;
    }

    AliAnalysisMuMuSpectra* realSpectra = GetSpectraFromConfig(binType,flavour);

    if ( !realSpectra ){
      AliError("could not get real spectra");
      return 0x0;
    }

    if (!realSpectra->HasValue("CoffNofJPsi")){
      if (!CorrectSpectra(binType,flavour,accEffSubResultName)){
        AliError("Could not get corrected spectra");
        return 0x0;
      }
    }

    AliAnalysisMuMuSpectra* simSpectra = SIM()->GetSpectra(binType,Event.Data(),DimuonTrigger.Data(),Centrality.Data(),PairSelection.Data(),flavour);

    if ( !simSpectra){
      AliErrorClass("could not get sim spectra");
      return 0x0;
    }

    AliAnalysisMuMuBinning* binning = realSpectra->Binning();
    TObjArray* bins = binning->CreateBinObjArray();
    TIter nextBin(bins);
    AliAnalysisMuMuBinning::Range* bin;
    Int_t i(0);
    AliAnalysisMuMuJpsiResult* r;

    while ( ( bin = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) ){

      r = static_cast<AliAnalysisMuMuJpsiResult*>(realSpectra->BinContentArray()->At(i));
      StdoutToAliDebug(1,std::cout << "bin=";r->Print(););

      AliAnalysisMuMuJpsiResult* rsim = static_cast<AliAnalysisMuMuJpsiResult*>(simSpectra->BinContentArray()->At(i));

      r->Set("NofInputJPsi",rsim->GetValue("NofInputJPsi",accEffSubResultName),rsim->GetErrorStat("NofInputJPsi",accEffSubResultName));
      r->Set("AccEffJPsi",rsim->GetValue("AccEffJPsi",accEffSubResultName),rsim->GetErrorStat("AccEffJPsi",accEffSubResultName));

      ++i;
    }

    delete bins;
    Update();
    return realSpectra;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::SetConfig(const AliAnalysisMuMuConfig& config)
{
  /// (re)set the config
  delete fConfig;
  fConfig = new AliAnalysisMuMuConfig(config);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMu::SetParticleNameFromFileName(const char* filename)
{
  /// Try to get the simulated particle name from the filename

  TString sFilename(filename);

  if ( sFilename.Contains("JPSI",TString::kIgnoreCase) )
  {
    SetParticleName("JPsi");
  }
  else if ( sFilename.Contains("PSIP",TString::kIgnoreCase) )
  {
    SetParticleName("PsiP");
  }
  else
  {
    AliError("Unknown Particle Name in simulation");
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisMuMu::GetFileNameAndDirectory(const char* filename)
{
  /// Split the filename[:subdir] into a file name and a subdir

  fFilename = gSystem->ExpandPathName(filename);

  fDirectory = "";

  if ( fFilename.CountChar(':') )
  {
    fDirectory = filename;
    Int_t colon = fFilename.Index(':');
    fFilename = fFilename(0,colon);
    fDirectory = fDirectory(colon+1,strlen(filename)-colon);
  }
}
