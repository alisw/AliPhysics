/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// Root
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TH1F.h>
#include <TList.h>
#include <TProfile.h>

// AliRoot
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"

#include "AliAnaWeights.h"

/// \cond CLASSIMP
ClassImp(AliAnaWeights) ;
/// \endcond

//______________________________________________________________
/// Constructor.
//______________________________________________________________
AliAnaWeights::AliAnaWeights() 
: TObject(), fDebug(0),
  fhCentralityWeight(0),
  fCentrality(0),
  fUseCentralityWeight(0),
  fDoMCParticlePtWeights(1),
  fEtaFunction(0), fPi0Function(0), 
  fCheckGeneratorName(0),
  fMCWeight(1.),
  fCurrFileName(0),
  fCheckMCCrossSection(kFALSE),
  fJustFillCrossSecHist(0),
  fhXsec(0),
  fhTrials(0),
  fPyEventHeader(0),
  fCheckPythiaEventHeader(0)
{
}

//______________________________________________________________
/// Destructor.
//______________________________________________________________
AliAnaWeights::~AliAnaWeights() 
{ 
  if ( fhCentralityWeight ) delete fhCentralityWeight ; 
  
  if ( fEtaFunction ) delete fEtaFunction ;
 
  if ( fPi0Function ) delete fPi0Function ;
}

//____________________________________________________
/// Init histogram pointers.
/// \return TList containing histograms.
//____________________________________________________
TList * AliAnaWeights::GetCreateOutputHistograms()
{
  TList * outputContainer = new TList() ;
  outputContainer->SetName("MCWeightHistogram") ;
    
  if(!fCheckMCCrossSection) return outputContainer;
  
  outputContainer->SetOwner(kFALSE);

  fhXsec = new TH1F("hXsec","xsec from pyxsec.root",1,0,1);
  fhXsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  outputContainer->Add(fhXsec);
    
  fhTrials = new TH1F("hTrials","trials root file",1,0,1);
  fhTrials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  outputContainer->Add(fhTrials);
  
  return outputContainer ;
}

//__________________________________
/// \return weight to be applied to the event.
/// The weight can be right now:
///  * pT-hard bin in PYTHIA MC events
///  * centrality dependent weight
//_________________________________
Double_t AliAnaWeights::GetWeight()
{
  Double_t weight = 1.;
      
  if ( fCheckMCCrossSection )
  {
      Double_t temp = GetPythiaCrossSection() ;
      
      AliDebug(1,Form("MC pT-hard weight: %e",temp));
      weight*=temp;
  }
    
  if ( fUseCentralityWeight )
  {
      Double_t temp = fhCentralityWeight->GetBinContent(fhCentralityWeight->FindBin(fCentrality));
      
      AliDebug(1,Form("Centrality %2.1f, weight: %2.2f",fCentrality,temp));
      
      weight*=temp;
  }
    
  AliDebug(1,Form("Event weight %e",weight));
    
  return weight ;
}

//_____________________________________________
//
/// \return the particle pT weight depending on
///
/// \param pt: input particle transverse momentum 
/// \param pdg: input particle type pdg code
/// \param genName: TString with generator name
/// \param igen: index of generator (not used yet).
///
/// currently, only pi0 and eta cases are considered. The parametrizations are passed via
/// SetEtaFunction(TF1* fun) and SetPi0Function(TF1* fun)
/// ex param: TF1* PowerEta0 = new TF1("PowerEta0","[0]*pow([1]/(([1]*exp(-[3]*x)+x)),[2])",4,25);
/// Only active when SwitchOnMCParticlePtWeights() 
//_____________________________________________
Double_t AliAnaWeights::GetParticlePtWeight(Float_t pt, Int_t pdg, TString genName, Int_t igen) const 
{
  Double_t weight = 1.;

  if ( !fDoMCParticlePtWeights ) return weight ;

  //printf("Get particle MC weight %p %p\n",fPi0Function,fEtaFunction);
  
  if ( !fCheckGeneratorName )
  {
    if      (pdg == 111 && fPi0Function)
    {
      weight = fPi0Function->Eval(pt);
      //printf("GetParticlePtWeights:: Pi0 w %2.3f, pt %2.3f\n",weight,pt);
    }
    else if (pdg == 221 && fEtaFunction)
      weight = fEtaFunction->Eval(pt);
  }
  else // Check particular generator names, to be better done
  {
    if      (pdg == 111 && fPi0Function && 
             ( genName.Contains("Pi0") || genName.Contains("pi0") || genName.Contains("PARAM") || genName.Contains("BOX")))
      weight = fPi0Function->Eval(pt);
    else if (pdg == 221 && fEtaFunction && 
             ( genName.Contains("Eta") || genName.Contains("eta") || genName.Contains("PARAM") || genName.Contains("BOX")))
      weight = fEtaFunction->Eval(pt);
  }
  
  AliDebug(1,Form("MC particle pdg %d, pt %2.2f, generator %s with index %d: weight %f",pdg,pt,genName.Data(),igen, weight));
  
  return weight ;
}

//_____________________________________________
//
/// Read the cross sections and number of trials 
/// from pyxsec.root (ESD) or pysec_hists.root (AODs), 
/// values stored in specific histograms-trees.
/// If no file available, get information from Pythia event header
///
/// Fill the control histograms on number of trials and xsection
/// optionally, with fJustFillCrossSecHist, just do that and do not return a weight.
/// For cross section obtained from pythia event header, this is already the case.
//_____________________________________________
Double_t AliAnaWeights::GetPythiaCrossSection()
{
  Float_t xsection  = 0;
  Float_t trials    = 1;
  Float_t avgTrials = 0;
  
  if ( !fhXsec || !fhTrials ) return 1;
  
  // -----------------------------------------------
  // Check if info is already in Pythia event header
  // Do not apply the weight per event, too much number of trial variation
  // -----------------------------------------------
  if ( fPyEventHeader && fCheckPythiaEventHeader )
  {    
    AliDebug(fDebug,Form("Pythia event header: xs %2.2e, trial %d", 
                         fPyEventHeader->GetXsection(),fPyEventHeader->Trials()));
    
    fhXsec  ->Fill("<#sigma>"     ,fPyEventHeader->GetXsection());
    fhTrials->Fill("#sum{ntrials}",fPyEventHeader->Trials());
    
    if ( !fJustFillCrossSecHist )
    {
      fMCWeight =  fPyEventHeader->GetXsection() / fPyEventHeader->Trials() ;
      AliDebug(1,Form("MC Weight: %e",fMCWeight));
    }
    else fMCWeight = 1;
    
    return  fMCWeight ;
  }

  // -----------------------------------------------
  // Get cross section from corresponding files
  // -----------------------------------------------

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if ( !tree    ) return 1;
  
  TFile *curfile = tree->GetCurrentFile();
  
  if ( !curfile ) return 1;
  
  // Check if file not accessed previously, if so
  // return the previously calculated weight
  if(fCurrFileName == curfile->GetName()) return fMCWeight;
  
  fCurrFileName = TString(curfile->GetName());
    
  Bool_t ok = GetPythiaInfoFromFile(fCurrFileName,xsection,trials);
  
  if ( !ok )
  {
      AliWarning("Parameters from file not recovered properly");
      return 1;
  }
    
  fhXsec->Fill("<#sigma>",xsection);

  // average number of trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();

  if(trials >= nEntries && nEntries > 0.) avgTrials = trials/nEntries;
  
  fhTrials->Fill("#sum{ntrials}",avgTrials);
  
  AliInfo(Form("xs %2.2e, trial %e, avg trials %2.2e, events per file %e",
               xsection,trials,avgTrials,nEntries));
  
  AliDebug(1,Form("Reading File %s",curfile->GetName()));
    
  if(avgTrials > 0.)
  {
      if(!fJustFillCrossSecHist) 
      {
        fMCWeight =  xsection / avgTrials ;
        
        AliInfo(Form("MC Weight: %e",fMCWeight));
      }
      else  fMCWeight = 1; // do not add weight to histograms
  }
  else
  {
      AliWarning(Form("Average number of trials is NULL!! Set weight to 1: xs : %e, trials %e, entries %e",
                      xsection,trials,nEntries));
      
      fMCWeight = 1;
  }
    
  return fMCWeight ;
}

//_______________________________________________________________________________________
/// This method gets and returns the 
/// \param file : either pyxsec.root (ESDs) or pysec_hists.root (AODs) files
/// \param xsec : cross section 
/// \param trials : number of event trials 
/// that should be located where the main data file is.
/// This is called in Notify and should provide the path to the AOD/ESD file
//_______________________________________________________________________________________
Bool_t AliAnaWeights::GetPythiaInfoFromFile(TString file,Float_t & xsec,Float_t & trials)
{
  xsec   = 0;
  trials = 1;
  
  if(file.Contains("root_archive.zip#"))
  {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos  = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  }
  else
  {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  
  //Printf("%s",file.Data());
  
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec)
  {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec)
    {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else
    {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if(!key)
      {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list)
      {
        fxsec->Close();
        return kFALSE;
      }
      
      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else
  {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree)
    {
      fxsec->Close();
      return kFALSE;
    }
    
    UInt_t   ntrials  = 0;
    Double_t xsection = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
    fxsec->Close();
  }
  
  return kTRUE;
}

//_______________________________________________________________________________________
/// This method intializes the histogram containing the centrality weights.
/// Use in case of different binning than default.
/// \param nbins : number of bins in histogram, default 100
/// \param minCen : minimum value of centrality, default 0
/// \param maxCen : maximum value of centrality, default 100
//_______________________________________________________________________________________
void AliAnaWeights::InitCentralityWeightsHistogram(Int_t nbins, Int_t minCen, Int_t maxCen)
{
  if ( fhCentralityWeight ) delete fhCentralityWeight ;
  
  fhCentralityWeight = new TH1F("hCentralityWeights","Centrality weights",nbins,minCen,maxCen);
    
  for(Int_t ibin = 0; ibin < fhCentralityWeight->GetNbinsX(); ibin++)
      fhCentralityWeight->SetBinContent(ibin,1.) ;
}

//_________________________________________________
/// \return histogram with centrality weights.
//_________________________________________________
TH1F* AliAnaWeights::GetCentralityWeightsHistogram()
{
  if ( !fhCentralityWeight ) InitCentralityWeightsHistogram();
    
  return fhCentralityWeight ;
}

