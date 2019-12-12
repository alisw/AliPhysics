/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake.                                                      *
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

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGrid.h>
#include <TFile.h>
#include <TSystem.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  #include <TPython.h>
#endif

#include "AliAnalysisManager.h"
#include "AliRhoParameter.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskEmcalJetCorrection.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetCorrection)
/// \endcond


//________________________________________________________________________
AliAnalysisTaskEmcalJetCorrection::AliAnalysisTaskEmcalJetCorrection() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetCorrection", kTRUE),
  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  fPythonCLI(),
  #endif
  fCustomPackages(),
  fPythonModulePath("lib/python3.6/site-packages/"),
  fJetsCont(),
  fTracksCont(),
  fJetOutputArray(),
  fBackgroundModelFileName(""),
  fBackgroundModelInputParameters(""),
  fModelName("Model")
{
  // nothing to do
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetCorrection::AliAnalysisTaskEmcalJetCorrection(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  fPythonCLI(),
  #endif
  fCustomPackages(),
  fPythonModulePath("lib/python3.6/site-packages/"),
  fJetsCont(),
  fTracksCont(),
  fJetOutputArray(),
  fBackgroundModelFileName(""),
  fBackgroundModelInputParameters(""),
  fModelName("Model")
  {
  // nothing to do
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetCorrection::~AliAnalysisTaskEmcalJetCorrection()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCorrection::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont           = GetJetContainer(0);
  if(!fJetsCont)
    AliFatal("Jet input container not found!");
  fTracksCont       = static_cast<AliParticleContainer*>(fJetsCont->GetParticleContainer());
  if(!fTracksCont)
    AliFatal("Particle input container not found attached to jets!");

  AddHistogram2D<TH2D>("hJetPtCorrelation", "Jet p_{T} (model vs. area-based)", "COLZ", 400, -100., 300., 400, -100., 300., "p_{T, ML} (GeV/c)", "p_{T} (GeV/c)", "dN2^{Jets}/d2p_{T}");

  PostData(1, fOutput);

}


//________________________________________________________________________
void AliAnalysisTaskEmcalJetCorrection::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();

  // ### Need to explicitly tell jet container to load rho mass object
  fJetsCont->LoadRhoMass(InputEvent());

  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    fPythonCLI = new TPython();
    // ### Install custom packages locally on worker node
    if(!fCustomPackages.IsNull())
    {
      std::ofstream outFile;
      outFile.open("./myScript.sh");
      outFile << "#!/bin/bash" << std::endl;
      outFile << "export PYTHONUSERBASE=./my-local-python" << std::endl;
      outFile << "export PATH=\"$PYTHONUSERBASE/bin:$PATH\"" << std::endl;
      TObjArray* packages = fCustomPackages.Tokenize(",");
      for(Int_t iPackage = 0; iPackage < packages->GetEntries(); iPackage++)
      {
        TString package = ((TObjString *)(packages->At(iPackage)))->String();
        outFile << "pip install --user " << package.Data() << std::endl;
      }
      outFile.close();
      gSystem->Exec("bash ./myScript.sh");

      packages->SetOwner();
      delete packages;
    }
    // ### Get background model from alien and load it
    TGrid::Connect("alien://");
    if (gSystem->AccessPathName(fBackgroundModelFileName.Data()))
      AliFatal(Form("Background model %s does not exist!", fBackgroundModelFileName.Data()));
    TFile::Cp(fBackgroundModelFileName.Data(), "./Model.joblib");
    fPythonCLI->Exec("import sys");
    fPythonCLI->Exec(Form("sys.path.insert(0, './my-local-python/%s')", fPythonModulePath.Data()));
    fPythonCLI->Exec("import sklearn, numpy, joblib");
    fPythonCLI->Exec("estimator = joblib.load(\"./Model.joblib\")");
  #endif

  // ### Prepare jet output array
  if (!(InputEvent()->FindListObject(Form("%s_RhoMVA_%s", fJetsCont->GetArrayName().Data(), fModelName.Data())))) {
    fJetOutputArray = new TClonesArray("AliEmcalJet");
    fJetOutputArray->SetName(Form("%s_RhoMVA_%s", fJetsCont->GetArrayName().Data(), fModelName.Data()));
    InputEvent()->AddObject(fJetOutputArray);
  }
  else {
    AliFatal(Form("Jet output array with name %s_RhoMVA_%s already in event!", fJetsCont->GetArrayName().Data(), fModelName.Data()));
  }

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetCorrection::Run()
{

  // ################################### MAIN JET LOOP
  fJetsCont->ResetCurrentID();
  Int_t jetCount = 0;
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    Float_t pt_ML = 0;
    GetPtFromModel(jet, pt_ML);
    AliEmcalJet *outJet = new ((*fJetOutputArray)[jetCount]) AliEmcalJet(*jet);
    outJet->SetPtEtaPhi(pt_ML, jet->Eta(), jet->Phi());
    FillHistogram("hJetPtCorrelation", pt_ML, jet->Pt()-jet->Area()*fJetsCont->GetRhoVal());
    jetCount++;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCorrection::GetPtFromModel(AliEmcalJet* jet, Float_t& pt_ML)
{
  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    TString arrayStr = GetBackgroundModelArrayString(jet);
    // ####### Run Python script that does inference on jet
    fPythonCLI->Exec(Form("data_inference = numpy.array([[%s]])", arrayStr.Data()));
    fPythonCLI->Exec("result = estimator.predict(data_inference)");
    pt_ML   = fPythonCLI->Eval("result[0]");
  #endif
}

//________________________________________________________________________
TString AliAnalysisTaskEmcalJetCorrection::GetBackgroundModelArrayString(AliEmcalJet* jet)
{
  // ####### Calculate inference input parameters
  std::vector<PWG::JETFW::AliEmcalParticleJetConstituent> tracks_sorted = jet->GetParticleConstituents();
  std::sort(tracks_sorted.rbegin(), tracks_sorted.rend());
  Double_t trackPts[100] = {0};
  for(Int_t i = 0; i < TMath::Min(Int_t(tracks_sorted.size()), 100); i++)
  {
    const AliVParticle* particle = tracks_sorted[i].GetParticle();
    trackPts[i] = particle->Pt();
  }

  // Calculate jet shapes that could be demanded
  Double_t leSub_noCorr = 0;
  Double_t angularity = 0;
  Double_t momentumDispersion = 0;
  Double_t trackPtMean = 0;
  Double_t trackPtMedian = 0;
  CalculateJetShapes(jet, leSub_noCorr, angularity, momentumDispersion, trackPtMean, trackPtMedian);

  TString resultStr = "";
  TObjArray* data_tokens = fBackgroundModelInputParameters.Tokenize(",");
  for(Int_t iToken = 0; iToken < data_tokens->GetEntries(); iToken++)
  {
    TString token = ((TObjString *)(data_tokens->At(iToken)))->String();
    if(token == "Jet_Pt")
      resultStr += Form("%E", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area());
    else if(token == "Jet_NumTracks")
      resultStr += Form("%E", (Double_t)tracks_sorted.size());
    else if(token == "Jet_Shape_Mass")
      resultStr += Form("%E", jet->M());
    else if(token == "Jet_Shape_Mass_DerivCorr_1")
      resultStr += Form("%E", jet->GetShapeProperties()->GetFirstOrderSubtracted());
    else if(token == "Jet_Shape_Mass_DerivCorr_2")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtracted());
    else if(token == "Jet_Shape_pTD_DerivCorr_1")
      resultStr += Form("%E", jet->GetShapeProperties()->GetFirstOrderSubtractedpTD());
    else if(token == "Jet_Shape_pTD_DerivCorr_2")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtractedpTD());
    else if(token == "Jet_Shape_LeSub")
      resultStr += Form("%E", leSub_noCorr);
    else if(token == "Jet_Shape_LeSub_DerivCorr")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub());
    else if(token == "Jet_Shape_Sigma2_DerivCorr_1")
      resultStr += Form("%E", jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2());
    else if(token == "Jet_Shape_Sigma2_DerivCorr_2")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2());
    else if(token == "Jet_Shape_Angularity")
      resultStr += Form("%E", angularity);
    else if(token == "Jet_Shape_Angularity_DerivCorr_1")
      resultStr += Form("%E", jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity());
    else if(token == "Jet_Shape_Angularity_DerivCorr_2")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity());
    else if(token == "Jet_Shape_Circularity_DerivCorr_1")
      resultStr += Form("%E", jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity());
    else if(token == "Jet_Shape_Circularity_DerivCorr_2")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity());
    else if(token == "Jet_Shape_NumTracks_DerivCorr")
      resultStr += Form("%E", jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent());
    else if(token == "Jet_Shape_MomentumDispersion")
      resultStr += Form("%E", momentumDispersion);
    else if(token == "Jet_Shape_TrackPtMean")
      resultStr += Form("%E", trackPtMean);
    else if(token == "Jet_Shape_TrackPtMedian")
      resultStr += Form("%E", trackPtMedian);
    else if(token == "Event_BackgroundDensity")
      resultStr += Form("%E", fJetsCont->GetRhoVal());
    else if(token == "Event_BackgroundDensityMass")
      resultStr += Form("%E", fJetsCont->GetRhoMassVal());
    else if(token == "Jet_Area")
      resultStr += Form("%E", jet->Area());
    else if(token.BeginsWith("Jet_TrackPt"))
    {
      TString num = token(11,(token.Length()-11));
      resultStr += Form("%E", trackPts[num.Atoi()]);
    }

    // Add comma after numbers
    if((iToken < data_tokens->GetEntries()-1))
      resultStr += ", ";
  }

  data_tokens->SetOwner();
  delete data_tokens;
  return resultStr;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCorrection::CalculateJetShapes(AliEmcalJet* jet, Double_t& leSub_noCorr, Double_t& angularity, Double_t& momentumDispersion, Double_t& trackPtMean, Double_t& trackPtMedian)
{
  // #### Calculate mean, median of constituents, radial moment (angularity), momentum dispersion, leSub (no correction)
  Double_t jetLeadingHadronPt = -999.;
  Double_t jetSubleadingHadronPt = -999.;
  Double_t jetSummedPt = 0;
  Double_t jetSummedPt2 = 0;
  trackPtMean = 0;
  trackPtMedian = 0;
  angularity = 0;
  momentumDispersion = 0;
  std::vector<PWG::JETFW::AliEmcalParticleJetConstituent> tracks_sorted = jet->GetParticleConstituents();
  std::sort(tracks_sorted.rbegin(), tracks_sorted.rend());
  Int_t     numTracks = tracks_sorted.size();
  if(!numTracks) return;
  Double_t* trackPts = new Double_t[numTracks];

  // Loop over all constituents and do jet shape calculations
  for (Int_t i=0;i<numTracks;i++)
  {
    const AliVParticle* particle = tracks_sorted[i].GetParticle();
    trackPtMean += particle->Pt();
    trackPts[i] = particle->Pt();
    if(particle->Pt() > jetLeadingHadronPt)
    {
      jetSubleadingHadronPt = jetLeadingHadronPt;
      jetLeadingHadronPt = particle->Pt();
    }
    else if(particle->Pt() > jetSubleadingHadronPt)
      jetSubleadingHadronPt = particle->Pt();

    Double_t deltaR = GetDistance(particle->Eta(), jet->Eta(), particle->Phi(), jet->Phi());
    jetSummedPt += particle->Pt();
    jetSummedPt2 += particle->Pt()*particle->Pt();
    angularity += particle->Pt() * deltaR;
  }

  if(numTracks)
  {
    trackPtMean   /= numTracks;
    trackPtMedian = TMath::Median(numTracks, trackPts);
  }

  if(numTracks > 1)
    leSub_noCorr = jetLeadingHadronPt - jetSubleadingHadronPt;
  else
    leSub_noCorr = jetLeadingHadronPt;

  if(jetSummedPt)
  {
    momentumDispersion = TMath::Sqrt(jetSummedPt2)/jetSummedPt;
    angularity /= jetSummedPt;
  }
}


//########################################################################
// HELPERS
//########################################################################

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetCorrection::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetCorrection::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetCorrection::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
{
  TH2* tmpHist = static_cast<TH2*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}

//________________________________________________________________________
inline void AliAnalysisTaskEmcalJetCorrection::FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add)
{
  TH3* tmpHist = static_cast<TH3*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  if(add)
    tmpHist->Fill(x,y,z,add);
  else
    tmpHist->Fill(x,y,z);
}


//________________________________________________________________________
template <class T> T* AliAnalysisTaskEmcalJetCorrection::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax);

  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskEmcalJetCorrection::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskEmcalJetCorrection::AddHistogram3D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, Int_t zBins, Double_t zMin, Double_t zMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax, zBins, zMin, zMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetCorrection::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

// ### ADDTASK MACRO
//________________________________________________________________________
AliAnalysisTaskEmcalJetCorrection* AliAnalysisTaskEmcalJetCorrection::AddTaskEmcalJetCorrection(TString modelName, TString trackArray, TString jetArray, TString rhoObject, Double_t jetRadius, const char* taskNameSuffix)
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Double_t    minJetEta          = 0.5;
  Double_t    minJetPt           = 0.15;
  Double_t    minTrackPt         = 0.15;
  Double_t    minJetAreaPerc     = 0.557;
  TString     suffix             = "";
  if(taskNameSuffix != 0)
    suffix = taskNameSuffix;

  // ###### Task name
  TString name("AliAnalysisTaskEmcalJetCorrection");
  if (jetArray != "") {
    name += "_";
    name += jetArray;
  }
  if (rhoObject != "") {
    name += "_";
    name += rhoObject;
  }
  if(modelName != "") {
    name += "_";
    name += modelName;
  }
  if (suffix != "") {
    name += "_";
    name += suffix;
  }

  // ###### Setup task with default settings
  AliAnalysisTaskEmcalJetCorrection* myTask = new AliAnalysisTaskEmcalJetCorrection(name);
  myTask->SetVzRange(-10.,10.);
  myTask->SetModelName(modelName);

  // Particle container and track pt cut
  AliParticleContainer* trackCont = 0;
  if(trackArray == "mcparticles")
    trackCont = myTask->AddMCParticleContainer(trackArray);
  else if(trackArray =="mctracks")
    trackCont = myTask->AddParticleContainer(trackArray);
  else
    trackCont = myTask->AddTrackContainer(trackArray);
  trackCont->SetParticlePtCut(minTrackPt);

  // Jet container
  AliJetContainer *jetCont = myTask->AddJetContainer(jetArray,6,jetRadius);
  if (jetCont) {
    jetCont->SetRhoName(rhoObject);
    jetCont->SetPercAreaCut(minJetAreaPerc);
    jetCont->SetJetPtCut(minJetPt);
    jetCont->SetLeadingHadronType(0);
    jetCont->SetPtBiasJetTrack(minTrackPt);
    jetCont->SetJetEtaLimits(-minJetEta, +minJetEta);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->SetMaxTrackPt(1000);
  }

  mgr->AddTask(myTask);

  // ###### Connect inputs/outputs
  mgr->ConnectInput  (myTask, 0,  mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (myTask, 1, mgr->CreateContainer(Form("%s", name.Data()), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", mgr->GetCommonFileName())) );

  return myTask;
}
