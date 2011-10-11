// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Per-Ivar Lønne <perivarlonne@gmail.com>               *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************/

/// @file   AliHLTTPCdEdxMonitoringComponent.cxx
/// @author Per-Ivar Lønne, Jochen Thaeder, Matthias Richter, Alexander Kalweit
/// @date   21.08.2011
/// @brief  Component for reading ESD from chain and produce a dEdx monitoring plot
///

#if __GNUC__>= 3
using namespace std;
#endif


#include "TSystem.h"
#include "AliESDtrackCuts.h"
#include <AliHLTDAQ.h>
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "TMap.h"
#include "TObjString.h"
#include "AliHLTDataTypes.h"
#include "AliHLTComponentBenchmark.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TH2F.h"

#include "AliHLTTPCdEdxMonitoringComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp( AliHLTTPCdEdxMonitoringComponent )
  
AliHLTTPCdEdxMonitoringComponent::AliHLTTPCdEdxMonitoringComponent() 
  :
  AliHLTProcessor(),
  fESDTrackCuts(NULL),
  fHist(),
  fxbins(),
  fxmin(),
  fxmax(),
  fybins(),
  fymin(),
  fymax()
{
  //Constructor
}

AliHLTTPCdEdxMonitoringComponent::~AliHLTTPCdEdxMonitoringComponent()
{
  //Destructor
}

// ######################################################################### //
const char* AliHLTTPCdEdxMonitoringComponent::GetComponentID()
{
  // get component id
  return "TPCdEdxMonitoring";
}

// ######################################################################### //
void AliHLTTPCdEdxMonitoringComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // get list of input data types
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
}

// ######################################################################### //
AliHLTComponentDataType AliHLTTPCdEdxMonitoringComponent::GetOutputDataType()
{
  // get the output data size
  return kAliHLTDataTypeHistogram|kAliHLTDataOriginHLT;
}

// ######################################################################### //
void AliHLTTPCdEdxMonitoringComponent::GetOutputDataSize( unsigned long& constBase, Double_t& inputMultiplier )
{
  // get output size estimator
  constBase = 4096;
  inputMultiplier = 0;
}

// ######################################################################### //
void AliHLTTPCdEdxMonitoringComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  
    if (!targetMap) return;
    targetMap->Add(new TObjString("HLT/ConfigTPC/TPCdEdxMonitoring"),
    new TObjString("configuration object"));
    targetMap->Add(new TObjString("GRP/GRP/Data"),
    new TObjString("GRP object"));
  
}

// ######################################################################### //
AliHLTComponent* AliHLTTPCdEdxMonitoringComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTTPCdEdxMonitoringComponent;
}

// ######################################################################### //
Int_t AliHLTTPCdEdxMonitoringComponent::DoInit( Int_t argc, const char** argv )
{
  // init the component
  Int_t iResult=0;
  fESDTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","HLT");

  fxbins=Int_t(300);
  fxmin=Double_t(-2);
  fxmax=Double_t(2);
  fybins=Int_t(500);
  fymin=Double_t(0);
  fymax=Double_t(500);
  
if (!fESDTrackCuts)
    {
      iResult=-ENOMEM;
    }

  if (iResult<0)
    {
      if (fESDTrackCuts)
	delete fESDTrackCuts;
      fESDTrackCuts = NULL;
    }      
       
  if (iResult>=0) 
    {
      SetDefaultConfiguration();
      TString cdbPath="HLT/ConfigTPC/"; 
      cdbPath+=GetComponentID();
      iResult=ConfigureFromCDBTObjString(cdbPath);
    
  
      if (iResult>=0) 
	{
	  iResult=ConfigureFromArgumentString(argc, argv);	
	}
    }

  if (iResult>=0) {
    HLTInfo("ESD track cuts : %s",fESDTrackCuts->GetTitle() );
  }
  fHist = new TH2F("hHLT", "HLT", fxbins, fxmin, fxmax, fybins, fymin, fymax);
  fHist->GetXaxis()->SetTitle("momentum/charge #frac{p}{z}  (GeV/c)");
  fHist->GetYaxis()->SetTitle("dE/dx in TPC (a.u.)");
  Plotstyle();
  return iResult;
}
  

// ######################################################################### //
Int_t AliHLTTPCdEdxMonitoringComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here
  if (fESDTrackCuts)
    delete fESDTrackCuts;
  fESDTrackCuts = NULL;
  if (fHist)
    delete fHist;
  fHist=NULL;
  return 0;
}


// ######################################################################### //
Int_t AliHLTTPCdEdxMonitoringComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					      AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function
  float sig;
  Double_t ptot;

  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;

  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");

  // input objects are not supposed to be changed by the component, so they
  // are defined const. However, the implementation of AliESDEvent does not
  // support this and we need the const_cast
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  if (esd != NULL) {
    esd->GetStdContent();
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      AliESDtrack* track = esd->GetTrack(i);
      sig=track->GetTPCsignal();
      if(!fESDTrackCuts->AcceptTrack(track)) continue;
      if (!track->GetInnerParam()) continue;
      ptot = track->GetInnerParam()->GetP()*track->GetSign();
      fHist->Fill(ptot, sig);
    }
  }
  // publish the histogram
  PushBack(fHist, kAliHLTDataTypeHistogram | kAliHLTDataOriginHLT);

  return 0;
}

// ######################################################################### //
Int_t AliHLTTPCdEdxMonitoringComponent::ScanConfigurationArgument(Int_t argc, const char** argv)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.
  
  if (argc<=0) return 0;
  Int_t ii =0;
  TString argument=argv[ii];
  
  if (argument.IsNull()) return 0;
  
  if( !fESDTrackCuts){
    HLTError("No ESD track cuts availible");
    return -ENOMEM;
  }
  

  //**********************************//
  //        Histogram Binning         //
  //**********************************//


  // -xbins
  if (argument.CompareTo("-xbins")==0) 
    {
      if (++ii>=argc) return -EPROTO;
      argument=argv[ii];

      fxbins = argument.Atoi();

      return 2;
    }

  // -xmin
  if (argument.CompareTo("-xmin")==0) 
    {
      if (++ii>=argc) return -EPROTO;
      argument=argv[ii];

      fxmin = argument.Atoi();

      return 2;
    }

  if (argument.CompareTo("-xmax")==0) 
    {
      if (++ii>=argc) return -EPROTO;
      argument=argv[ii];

      fxmax = argument.Atoi();

      return 2;
    }

  // -xbins
  if (argument.CompareTo("-ybins")==0) 
    {
      if (++ii>=argc) return -EPROTO;
      argument=argv[ii];

      fybins = argument.Atoi();

      return 2;
    }

  // -xmin
  if (argument.CompareTo("-ymin")==0) 
    {
      if (++ii>=argc) return -EPROTO;
      argument=argv[ii];

      fymin = argument.Atoi();

      return 2;
    }

  if (argument.CompareTo("-ymax")==0) 
    {
      if (++ii>=argc) return -EPROTO;
      argument=argv[ii];

      fymax = argument.Atoi();

      return 2;
    }


  //**********************************//
  //           Track Cuts             //
  //**********************************//

 // -maxpt
  if (argument.CompareTo("-maxpt")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    Float_t minPt, maxPt;
    fESDTrackCuts->GetPtRange(minPt,maxPt);
    maxPt = argument.Atof(); 
    fESDTrackCuts->SetPtRange(minPt,maxPt);

    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("p_t < %f", maxPt);
    fESDTrackCuts->SetTitle(title);
    return 2;
  }    

  // -minpt
  if (argument.CompareTo("-minpt")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    Float_t minPt, maxPt;
    fESDTrackCuts->GetPtRange(minPt,maxPt);
    minPt = argument.Atof(); 
    fESDTrackCuts->SetPtRange(minPt,maxPt);

    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("p_t > %f", minPt);
    fESDTrackCuts->SetTitle(title);
    return 2;
  }    

  // -min-ldca
  // minimum longitudinal dca to vertex
  if (argument.CompareTo("-min-ldca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMinDCAToVertexZ(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAz > %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }
  
  // -max-ldca
  // maximum longitudinal dca to vertex
  if (argument.CompareTo("-max-ldca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMaxDCAToVertexZ(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAz < %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -min-tdca
  // minimum transverse dca to vertex
  if (argument.CompareTo("-min-tdca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMinDCAToVertexXY(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAr > %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }
  
  // -max-tdca
  // maximum transverse dca to vertex
  if (argument.CompareTo("-max-tdca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMaxDCAToVertexXY(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAr < %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -etarange
  // +/- eta 
  if (argument.CompareTo("-etarange")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t eta = argument.Atof();

    fESDTrackCuts->SetEtaRange(-eta,eta);     
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("Eta[%f,%f]", argument.Atof(),argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -minNClsTPC
  // minimum clusters in TPC
  if (argument.CompareTo("-minNClsTPC")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    Int_t ncls;
    ncls = Int_t(argument.Atof()); 
    fESDTrackCuts->SetMinNClustersTPC(ncls);  

    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("minNClsTPC < %i", ncls);
    fESDTrackCuts->SetTitle(title);
    return 2;
  }    



  // unknown argument
  return -EINVAL;
  
  return 0;
  
}

// ######################################################################### //
Int_t AliHLTTPCdEdxMonitoringComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  HLTInfo("reconfigure '%s' from entry %s", chainId, cdbEntry);

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry)
    {
      cdbPath=cdbEntry;
    } else {
      cdbPath="HLT/ConfigTPC/";
      cdbPath+=GetComponentID();
    }

  iResult=ConfigureFromCDBTObjString(cdbPath); //// Or use return 0, and skip this line?

  return iResult;
}

// ######################################################################### //
Int_t AliHLTTPCdEdxMonitoringComponent::ReadPreprocessorValues(const char* modules)
{
  // read the preprocessor values for the detectors in the modules list
  Int_t iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult; //Done differently in AliHLTMultiplicityCorrelationsComponent...
}



// ######################################################################### //
void AliHLTTPCdEdxMonitoringComponent::SetDefaultConfiguration() 
{
  if (fESDTrackCuts)
    {
      //fESDTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
      fESDTrackCuts->SetEtaRange(-0.8,+0.8);
      fESDTrackCuts->SetPtRange(0.15,1e10);
      fESDTrackCuts->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);  // BEWARE STANDARD VALUES ARE: 2, 2, 0.5, 0.5, 2
      fESDTrackCuts->SetMaxNsigmaToVertex(3);
      fESDTrackCuts->SetRequireSigmaToVertex(kTRUE);
      fESDTrackCuts->SetAcceptKinkDaughters(kFALSE);
      fESDTrackCuts->SetMinNClustersTPC(70);
      fESDTrackCuts->SetMaxChi2PerClusterTPC(4);
      fESDTrackCuts->SetMaxDCAToVertexXY(3);
      fESDTrackCuts->SetMaxDCAToVertexZ(3);
      fESDTrackCuts->SetRequireTPCRefit(kTRUE);
      //fESDTrackCuts->SetRequireITSRefit(kTRUE); //Kills HLT simulated reconstructions?
      fESDTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //TEMPORARY <-> REMOVE
      fESDTrackCuts->SetMinNClustersITS(3);
      
    }
 fxbins=300;
 fxmin=-2;
 fxmax=2;
 fybins=500;  
 fymin=0;
 fymax=500;
  return;
}



// ######################################################################### //
void AliHLTTPCdEdxMonitoringComponent::Plotstyle()
{
  gROOT->SetStyle("Plain");
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasDefH(550);
  gStyle->SetCanvasDefW(575);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetStatColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPalette(1,0);
  
  //
  
  //gStyle->SetStatX(0.7);
  //gStyle->SetStatW(0.2);
  //gStyle->SetLabelOffset(1.2);
  //gStyle->SetLabelFont(72);
  //gStyle->SetLabelSize(0.6);
  //gStyle->SetTitleOffset(1.2);
  gStyle->SetTitleFontSize(0.04);
  
  
  gStyle->SetOptStat(10);
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetTextSize(0.04);
  gStyle->SetTitleSize(0.04,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.02,"xyz");
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleOffset(1.3,"x");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleOffset(1.6,"z");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleColor(1,"xyz");
  //gStyle->SetPadTopMargin(0.1);
  //gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadLeftMargin(0.2);
  
  
  const Int_t NCont=255;
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont); 
}
