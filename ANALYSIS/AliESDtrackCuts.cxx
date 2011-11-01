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

/* $Id: AliESDtrackCuts.cxx 24534 2008-03-16 22:22:11Z fca $ */

#include "AliESDtrackCuts.h"

#include <AliESDtrack.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliMultiplicity.h>
#include <AliLog.h>

#include <TTree.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TH2F.h>
#include <TF1.h>

//____________________________________________________________________
ClassImp(AliESDtrackCuts)

// Cut names
const Char_t* AliESDtrackCuts::fgkCutNames[kNCuts] = {
 "require TPC refit",
 "require TPC standalone",
 "require ITS refit",
 "n clusters TPC",
 "n clusters ITS",
 "#Chi^{2}/cluster TPC",
 "#Chi^{2}/cluster ITS",
 "cov 11",
 "cov 22",
 "cov 33",
 "cov 44",
 "cov 55",
 "trk-to-vtx",
 "trk-to-vtx failed",
 "kink daughters",
 "p",
 "p_{T}",
 "p_{x}",
 "p_{y}",
 "p_{z}",
 "eta",
 "y",
 "trk-to-vtx max dca 2D absolute",
 "trk-to-vtx max dca xy absolute",
 "trk-to-vtx max dca z absolute",
 "trk-to-vtx min dca 2D absolute",
 "trk-to-vtx min dca xy absolute",
 "trk-to-vtx min dca z absolute",
 "SPD cluster requirement",
 "SDD cluster requirement",
 "SSD cluster requirement",
 "require ITS stand-alone",
 "rel 1/pt uncertainty",
 "TPC n shared clusters",
 "TPC rel shared clusters",
 "require ITS Pid",
 "n crossed rows TPC",
 "n crossed rows / n findable clusters",
 "missing ITS points",
 "#Chi^{2} TPC constrained vs. global"
};

AliESDtrackCuts* AliESDtrackCuts::fgMultEstTrackCuts[AliESDtrackCuts::kNMultEstTrackCuts] = { 0, 0, 0, 0 };

//____________________________________________________________________
AliESDtrackCuts::AliESDtrackCuts(const Char_t* name, const Char_t* title) : AliAnalysisCuts(name,title),
  fCutMinNClusterTPC(0),
  fCutMinNClusterITS(0),
  fCutMinNCrossedRowsTPC(0),
  fCutMinRatioCrossedRowsOverFindableClustersTPC(0),
  f1CutMinNClustersTPCPtDep(0x0),
  fCutMaxPtDepNClustersTPC(0),
  fCutMaxChi2PerClusterTPC(0),
  fCutMaxChi2PerClusterITS(0),
  fCutMaxChi2TPCConstrainedVsGlobal(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(kVertexTracks | kVertexSPD),
  fCutMaxMissingITSPoints(0),
  fCutMaxC11(0),
  fCutMaxC22(0),
  fCutMaxC33(0),
  fCutMaxC44(0),
  fCutMaxC55(0),
  fCutMaxRel1PtUncertainty(0),
  fCutAcceptKinkDaughters(0),
  fCutAcceptSharedTPCClusters(0),
  fCutMaxFractionSharedTPCClusters(0),
  fCutRequireTPCRefit(0),
  fCutRequireTPCStandAlone(0),
  fCutRequireITSRefit(0), 
  fCutRequireITSPid(0),
  fCutRequireITSStandAlone(0),
  fCutRequireITSpureSA(0),
  fCutNsigmaToVertex(0),
  fCutSigmaToVertexRequired(0),
  fCutMaxDCAToVertexXY(0),
  fCutMaxDCAToVertexZ(0),
  fCutMinDCAToVertexXY(0),
  fCutMinDCAToVertexZ(0),
  fCutMaxDCAToVertexXYPtDep(""),
  fCutMaxDCAToVertexZPtDep(""),
  fCutMinDCAToVertexXYPtDep(""),
  fCutMinDCAToVertexZPtDep(""),
  f1CutMaxDCAToVertexXYPtDep(0x0),
  f1CutMaxDCAToVertexZPtDep(0x0),
  f1CutMinDCAToVertexXYPtDep(0x0),
  f1CutMinDCAToVertexZPtDep(0x0),
  fCutDCAToVertex2D(0),
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapMin(0),
  fRapMax(0),
  fHistogramsOn(0),
  ffDTheoretical(0),
  fhCutStatistics(0),         
  fhCutCorrelation(0)
{
  //
  // constructor
  //

  Init();

  //##############################################################################
  // setting default cuts
  SetMinNClustersTPC();
  SetMinNClustersITS();
  SetMinNCrossedRowsTPC();
  SetMinRatioCrossedRowsOverFindableClustersTPC();
  SetMaxChi2PerClusterTPC();
  SetMaxChi2PerClusterITS();
  SetMaxChi2TPCConstrainedGlobal();
  SetMaxChi2TPCConstrainedGlobalVertexType();
  SetMaxNOfMissingITSPoints();
  SetMaxCovDiagonalElements();
  SetMaxRel1PtUncertainty();
  SetRequireTPCRefit();
  SetRequireTPCStandAlone();
  SetRequireITSRefit();
  SetRequireITSPid(kFALSE);
  SetRequireITSStandAlone(kFALSE);
  SetRequireITSPureStandAlone(kFALSE);
  SetAcceptKinkDaughters();
  SetAcceptSharedTPCClusters();
  SetMaxFractionSharedTPCClusters();
  SetMaxNsigmaToVertex();
  SetMaxDCAToVertexXY();
  SetMaxDCAToVertexZ();
  SetDCAToVertex2D();
  SetMinDCAToVertexXY();
  SetMinDCAToVertexZ();
  SetPRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();
  SetEtaRange();
  SetRapRange();
  SetClusterRequirementITS(kSPD);
  SetClusterRequirementITS(kSDD);
  SetClusterRequirementITS(kSSD);

  SetHistogramsOn();
}

//_____________________________________________________________________________
AliESDtrackCuts::AliESDtrackCuts(const AliESDtrackCuts &c) : AliAnalysisCuts(c),
  fCutMinNClusterTPC(0),
  fCutMinNClusterITS(0),
  fCutMinNCrossedRowsTPC(0),
  fCutMinRatioCrossedRowsOverFindableClustersTPC(0),
  f1CutMinNClustersTPCPtDep(0x0),
  fCutMaxPtDepNClustersTPC(0),
  fCutMaxChi2PerClusterTPC(0),
  fCutMaxChi2PerClusterITS(0),
  fCutMaxChi2TPCConstrainedVsGlobal(0),
  fCutMaxChi2TPCConstrainedVsGlobalVertexType(kVertexTracks | kVertexSPD),
  fCutMaxMissingITSPoints(0),
  fCutMaxC11(0),
  fCutMaxC22(0),
  fCutMaxC33(0),
  fCutMaxC44(0),
  fCutMaxC55(0),
  fCutMaxRel1PtUncertainty(0),
  fCutAcceptKinkDaughters(0),
  fCutAcceptSharedTPCClusters(0),
  fCutMaxFractionSharedTPCClusters(0),
  fCutRequireTPCRefit(0),
  fCutRequireTPCStandAlone(0),
  fCutRequireITSRefit(0),
  fCutRequireITSPid(0),
  fCutRequireITSStandAlone(0),
  fCutRequireITSpureSA(0),
  fCutNsigmaToVertex(0),
  fCutSigmaToVertexRequired(0),
  fCutMaxDCAToVertexXY(0),
  fCutMaxDCAToVertexZ(0),
  fCutMinDCAToVertexXY(0),
  fCutMinDCAToVertexZ(0),
  fCutMaxDCAToVertexXYPtDep(""),
  fCutMaxDCAToVertexZPtDep(""),
  fCutMinDCAToVertexXYPtDep(""),
  fCutMinDCAToVertexZPtDep(""),
  f1CutMaxDCAToVertexXYPtDep(0x0),
  f1CutMaxDCAToVertexZPtDep(0x0),
  f1CutMinDCAToVertexXYPtDep(0x0),
  f1CutMinDCAToVertexZPtDep(0x0),
  fCutDCAToVertex2D(0),
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fRapMin(0),
  fRapMax(0),
  fHistogramsOn(0),
  ffDTheoretical(0),				     
  fhCutStatistics(0),         
  fhCutCorrelation(0)
{
  //
  // copy constructor
  //

  ((AliESDtrackCuts &) c).Copy(*this);
}

AliESDtrackCuts::~AliESDtrackCuts()
{
  //
  // destructor
  //

  for (Int_t i=0; i<2; i++) {
    
    if (fhNClustersITS[i])
      delete fhNClustersITS[i];            
    if (fhNClustersTPC[i])
      delete fhNClustersTPC[i];                
    if (fhNSharedClustersTPC[i])
      delete fhNSharedClustersTPC[i];                
    if (fhNCrossedRowsTPC[i])
      delete fhNCrossedRowsTPC[i];                
    if (fhRatioCrossedRowsOverFindableClustersTPC[i])
      delete fhRatioCrossedRowsOverFindableClustersTPC[i];                
    if (fhChi2PerClusterITS[i])
      delete fhChi2PerClusterITS[i];       
    if (fhChi2PerClusterTPC[i])
      delete fhChi2PerClusterTPC[i];    
    if (fhChi2TPCConstrainedVsGlobal[i])
      delete fhChi2TPCConstrainedVsGlobal[i];
    if(fhNClustersForITSPID[i])
      delete fhNClustersForITSPID[i];
    if(fhNMissingITSPoints[i])
      delete fhNMissingITSPoints[i];
    if (fhC11[i])
      delete fhC11[i];                     
    if (fhC22[i])
      delete fhC22[i];                     
    if (fhC33[i])
      delete fhC33[i];                     
    if (fhC44[i])
      delete fhC44[i];                     
    if (fhC55[i])
      delete fhC55[i];

    if (fhRel1PtUncertainty[i])
      delete fhRel1PtUncertainty[i];
    
    if (fhDXY[i])
      delete fhDXY[i];                     
    if (fhDZ[i])
      delete fhDZ[i];
    if (fhDXYDZ[i])
      delete fhDXYDZ[i];
    if (fhDXYvsDZ[i])
      delete fhDXYvsDZ[i];

    if (fhDXYNormalized[i])
      delete fhDXYNormalized[i];           
    if (fhDZNormalized[i])
      delete fhDZNormalized[i];
    if (fhDXYvsDZNormalized[i])
      delete fhDXYvsDZNormalized[i];
    if (fhNSigmaToVertex[i])
      delete fhNSigmaToVertex[i];
    if (fhPt[i])
      delete fhPt[i];
    if (fhEta[i])
      delete fhEta[i];
  }

  if(f1CutMaxDCAToVertexXYPtDep)delete f1CutMaxDCAToVertexXYPtDep;
  f1CutMaxDCAToVertexXYPtDep = 0;
  if( f1CutMaxDCAToVertexZPtDep) delete  f1CutMaxDCAToVertexZPtDep;
  f1CutMaxDCAToVertexZPtDep = 0;
  if( f1CutMinDCAToVertexXYPtDep)delete  f1CutMinDCAToVertexXYPtDep;
  f1CutMinDCAToVertexXYPtDep = 0;
  if(f1CutMinDCAToVertexZPtDep)delete  f1CutMinDCAToVertexZPtDep; 
  f1CutMinDCAToVertexZPtDep = 0;


  if (ffDTheoretical)
    delete ffDTheoretical;

  if (fhCutStatistics)
    delete fhCutStatistics;             
  if (fhCutCorrelation)
    delete fhCutCorrelation;    

  if(f1CutMinNClustersTPCPtDep)
    delete f1CutMinNClustersTPCPtDep;
        
}

void AliESDtrackCuts::Init()
{
  //
  // sets everything to zero
  //

  fCutMinNClusterTPC = 0;
  fCutMinNClusterITS = 0;

  fCutMaxChi2PerClusterTPC = 0;
  fCutMaxChi2PerClusterITS = 0;
  fCutMaxChi2TPCConstrainedVsGlobal = 0;
  fCutMaxChi2TPCConstrainedVsGlobalVertexType = kVertexTracks | kVertexSPD;
  fCutMaxMissingITSPoints  = 0;
  
  for (Int_t i = 0; i < 3; i++)
  	fCutClusterRequirementITS[i] = kOff;

  fCutMaxC11 = 0;
  fCutMaxC22 = 0;
  fCutMaxC33 = 0;
  fCutMaxC44 = 0;
  fCutMaxC55 = 0;
  
  fCutMaxRel1PtUncertainty = 0;

  fCutAcceptKinkDaughters = 0;
  fCutAcceptSharedTPCClusters = 0;
  fCutMaxFractionSharedTPCClusters = 0;
  fCutRequireTPCRefit = 0;
  fCutRequireTPCStandAlone = 0;
  fCutRequireITSRefit = 0;
  fCutRequireITSPid = 0;
  fCutRequireITSStandAlone = 0;
  fCutRequireITSpureSA = 0;

  fCutNsigmaToVertex = 0;
  fCutSigmaToVertexRequired = 0;
  fCutMaxDCAToVertexXY = 0;
  fCutMaxDCAToVertexZ = 0;
  fCutDCAToVertex2D = 0;
  fCutMinDCAToVertexXY = 0;
  fCutMinDCAToVertexZ = 0;
  fCutMaxDCAToVertexXYPtDep = "";
  fCutMaxDCAToVertexZPtDep = "";
  fCutMinDCAToVertexXYPtDep = "";
  fCutMinDCAToVertexZPtDep = "";

  if(f1CutMaxDCAToVertexXYPtDep)delete f1CutMaxDCAToVertexXYPtDep;
  f1CutMaxDCAToVertexXYPtDep = 0;
  if( f1CutMaxDCAToVertexXYPtDep) delete  f1CutMaxDCAToVertexXYPtDep;
  f1CutMaxDCAToVertexXYPtDep = 0;
  if( f1CutMaxDCAToVertexZPtDep) delete  f1CutMaxDCAToVertexZPtDep;
  f1CutMaxDCAToVertexZPtDep = 0;
  if( f1CutMinDCAToVertexXYPtDep)delete  f1CutMinDCAToVertexXYPtDep;
  f1CutMinDCAToVertexXYPtDep = 0;
  if(f1CutMinDCAToVertexZPtDep)delete f1CutMinDCAToVertexZPtDep;
  f1CutMinDCAToVertexZPtDep = 0;

  
  fPMin = 0;
  fPMax = 0;
  fPtMin = 0;
  fPtMax = 0;
  fPxMin = 0;
  fPxMax = 0;
  fPyMin = 0;
  fPyMax = 0;
  fPzMin = 0;
  fPzMax = 0;
  fEtaMin = 0;
  fEtaMax = 0;
  fRapMin = 0;
  fRapMax = 0;

  fHistogramsOn = kFALSE;

  for (Int_t i=0; i<2; ++i)
  {
    fhNClustersITS[i] = 0;
    fhNClustersTPC[i] = 0;
    fhNSharedClustersTPC[i] = 0;
    fhNCrossedRowsTPC[i] = 0;
    fhRatioCrossedRowsOverFindableClustersTPC[i] = 0;

    fhChi2PerClusterITS[i] = 0;
    fhChi2PerClusterTPC[i] = 0;
    fhChi2TPCConstrainedVsGlobal[i] = 0;
    fhNClustersForITSPID[i] = 0;
    fhNMissingITSPoints[i] = 0;

    fhC11[i] = 0;
    fhC22[i] = 0;
    fhC33[i] = 0;
    fhC44[i] = 0;
    fhC55[i] = 0;

    fhRel1PtUncertainty[i] = 0;

    fhDXY[i] = 0;
    fhDZ[i] = 0;
    fhDXYDZ[i] = 0;
    fhDXYvsDZ[i] = 0;

    fhDXYNormalized[i] = 0;
    fhDZNormalized[i] = 0;
    fhDXYvsDZNormalized[i] = 0;
    fhNSigmaToVertex[i] = 0;
    
    fhPt[i] = 0;
    fhEta[i] = 0;
  }
  ffDTheoretical = 0;

  fhCutStatistics = 0;
  fhCutCorrelation = 0;
}

//_____________________________________________________________________________
AliESDtrackCuts &AliESDtrackCuts::operator=(const AliESDtrackCuts &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliESDtrackCuts &) c).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliESDtrackCuts::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliESDtrackCuts& target = (AliESDtrackCuts &) c;

  target.Init();

  target.fCutMinNClusterTPC = fCutMinNClusterTPC;
  target.fCutMinNClusterITS = fCutMinNClusterITS;
  target.fCutMinNCrossedRowsTPC = fCutMinNCrossedRowsTPC;
  target.fCutMinRatioCrossedRowsOverFindableClustersTPC = fCutMinRatioCrossedRowsOverFindableClustersTPC;
  if(f1CutMinNClustersTPCPtDep){
    target.f1CutMinNClustersTPCPtDep = (TFormula*) f1CutMinNClustersTPCPtDep->Clone("f1CutMinNClustersTPCPtDep");
  }
  target.fCutMaxPtDepNClustersTPC =   fCutMaxPtDepNClustersTPC;

  target.fCutMaxChi2PerClusterTPC = fCutMaxChi2PerClusterTPC;
  target.fCutMaxChi2PerClusterITS = fCutMaxChi2PerClusterITS;
  target.fCutMaxChi2TPCConstrainedVsGlobal = fCutMaxChi2TPCConstrainedVsGlobal;
  target.fCutMaxChi2TPCConstrainedVsGlobalVertexType = fCutMaxChi2TPCConstrainedVsGlobalVertexType;
  target.fCutMaxMissingITSPoints = fCutMaxMissingITSPoints;

  for (Int_t i = 0; i < 3; i++)
    target.fCutClusterRequirementITS[i] = fCutClusterRequirementITS[i];

  target.fCutMaxC11 = fCutMaxC11;
  target.fCutMaxC22 = fCutMaxC22;
  target.fCutMaxC33 = fCutMaxC33;
  target.fCutMaxC44 = fCutMaxC44;
  target.fCutMaxC55 = fCutMaxC55;

  target.fCutMaxRel1PtUncertainty = fCutMaxRel1PtUncertainty;

  target.fCutAcceptKinkDaughters = fCutAcceptKinkDaughters;
  target.fCutAcceptSharedTPCClusters = fCutAcceptSharedTPCClusters;
  target.fCutMaxFractionSharedTPCClusters = fCutMaxFractionSharedTPCClusters;
  target.fCutRequireTPCRefit = fCutRequireTPCRefit;
  target.fCutRequireTPCStandAlone = fCutRequireTPCStandAlone;
  target.fCutRequireITSRefit = fCutRequireITSRefit;
  target.fCutRequireITSPid = fCutRequireITSPid;
  target.fCutRequireITSStandAlone = fCutRequireITSStandAlone;
  target.fCutRequireITSpureSA = fCutRequireITSpureSA;

  target.fCutNsigmaToVertex = fCutNsigmaToVertex;
  target.fCutSigmaToVertexRequired = fCutSigmaToVertexRequired;
  target.fCutMaxDCAToVertexXY = fCutMaxDCAToVertexXY;
  target.fCutMaxDCAToVertexZ = fCutMaxDCAToVertexZ;
  target.fCutDCAToVertex2D = fCutDCAToVertex2D;
  target.fCutMinDCAToVertexXY = fCutMinDCAToVertexXY;
  target.fCutMinDCAToVertexZ = fCutMinDCAToVertexZ;

  target.fCutMaxDCAToVertexXYPtDep = fCutMaxDCAToVertexXYPtDep;
  target.SetMaxDCAToVertexXYPtDep(fCutMaxDCAToVertexXYPtDep.Data());

  target.fCutMaxDCAToVertexZPtDep = fCutMaxDCAToVertexZPtDep;
  target.SetMaxDCAToVertexZPtDep(fCutMaxDCAToVertexZPtDep.Data());

  target.fCutMinDCAToVertexXYPtDep = fCutMinDCAToVertexXYPtDep;
  target.SetMinDCAToVertexXYPtDep(fCutMinDCAToVertexXYPtDep.Data());

  target.fCutMinDCAToVertexZPtDep = fCutMinDCAToVertexZPtDep;
  target.SetMinDCAToVertexZPtDep(fCutMinDCAToVertexZPtDep.Data());

  target.fPMin = fPMin;
  target.fPMax = fPMax;
  target.fPtMin = fPtMin;
  target.fPtMax = fPtMax;
  target.fPxMin = fPxMin;
  target.fPxMax = fPxMax;
  target.fPyMin = fPyMin;
  target.fPyMax = fPyMax;
  target.fPzMin = fPzMin;
  target.fPzMax = fPzMax;
  target.fEtaMin = fEtaMin;
  target.fEtaMax = fEtaMax;
  target.fRapMin = fRapMin;
  target.fRapMax = fRapMax;

  target.fHistogramsOn = fHistogramsOn;

  for (Int_t i=0; i<2; ++i)
  {
    if (fhNClustersITS[i]) target.fhNClustersITS[i] = (TH1F*) fhNClustersITS[i]->Clone();
    if (fhNClustersTPC[i]) target.fhNClustersTPC[i] = (TH1F*) fhNClustersTPC[i]->Clone();
    if (fhNSharedClustersTPC[i]) target.fhNSharedClustersTPC[i] = (TH1F*) fhNSharedClustersTPC[i]->Clone();
    if (fhNCrossedRowsTPC[i]) target.fhNCrossedRowsTPC[i] = (TH1F*) fhNCrossedRowsTPC[i]->Clone();
    if (fhRatioCrossedRowsOverFindableClustersTPC[i]) target.fhRatioCrossedRowsOverFindableClustersTPC[i] = (TH1F*) fhRatioCrossedRowsOverFindableClustersTPC[i]->Clone();

    if (fhChi2PerClusterITS[i]) target.fhChi2PerClusterITS[i] = (TH1F*) fhChi2PerClusterITS[i]->Clone();
    if (fhChi2PerClusterTPC[i]) target.fhChi2PerClusterTPC[i] = (TH1F*) fhChi2PerClusterTPC[i]->Clone();
    if (fhChi2TPCConstrainedVsGlobal[i]) target.fhChi2TPCConstrainedVsGlobal[i] = (TH1F*) fhChi2TPCConstrainedVsGlobal[i]->Clone();
    if (fhNClustersForITSPID[i]) target.fhNClustersForITSPID[i] = (TH1F*) fhNClustersForITSPID[i]->Clone();
    if (fhNMissingITSPoints[i]) target.fhNMissingITSPoints[i] = (TH1F*) fhNMissingITSPoints[i]->Clone();

    if (fhC11[i]) target.fhC11[i] = (TH1F*) fhC11[i]->Clone();
    if (fhC22[i]) target.fhC22[i] = (TH1F*) fhC22[i]->Clone();
    if (fhC33[i]) target.fhC33[i] = (TH1F*) fhC33[i]->Clone();
    if (fhC44[i]) target.fhC44[i] = (TH1F*) fhC44[i]->Clone();
    if (fhC55[i]) target.fhC55[i] = (TH1F*) fhC55[i]->Clone();

    if (fhRel1PtUncertainty[i]) target.fhRel1PtUncertainty[i] = (TH1F*) fhRel1PtUncertainty[i]->Clone();

    if (fhDXY[i]) target.fhDXY[i] = (TH1F*) fhDXY[i]->Clone();
    if (fhDZ[i]) target.fhDZ[i] = (TH1F*) fhDZ[i]->Clone();
    if (fhDXYDZ[i]) target.fhDXYDZ[i] = (TH1F*) fhDXYDZ[i]->Clone();
    if (fhDXYvsDZ[i]) target.fhDXYvsDZ[i] = (TH2F*) fhDXYvsDZ[i]->Clone();

    if (fhDXYNormalized[i]) target.fhDXYNormalized[i] = (TH1F*) fhDXYNormalized[i]->Clone();
    if (fhDZNormalized[i]) target.fhDZNormalized[i] = (TH1F*) fhDZNormalized[i]->Clone();
    if (fhDXYvsDZNormalized[i]) target.fhDXYvsDZNormalized[i] = (TH2F*) fhDXYvsDZNormalized[i]->Clone();
    if (fhNSigmaToVertex[i]) target.fhNSigmaToVertex[i] = (TH1F*) fhNSigmaToVertex[i]->Clone();
    
    if (fhPt[i]) target.fhPt[i] = (TH1F*) fhPt[i]->Clone();
    if (fhEta[i]) target.fhEta[i] = (TH1F*) fhEta[i]->Clone();
  }
  if (ffDTheoretical) target.ffDTheoretical = (TF1*) ffDTheoretical->Clone();

  if (fhCutStatistics) target.fhCutStatistics = (TH1F*) fhCutStatistics->Clone();
  if (fhCutCorrelation) target.fhCutCorrelation = (TH2F*) fhCutCorrelation->Clone();

  TNamed::Copy(c);
}

//_____________________________________________________________________________
Long64_t AliESDtrackCuts::Merge(TCollection* list) {
  // Merge a list of AliESDtrackCuts objects with this (needed for PROOF)
  // Returns the number of merged objects (including this)
  if (!list)
    return 0;
  if (list->IsEmpty())
    return 1;
  if (!fHistogramsOn)
    return 0;
  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collection of measured and generated histograms
  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliESDtrackCuts* entry = dynamic_cast<AliESDtrackCuts*>(obj);
    if (entry == 0)
      continue;

    if (!entry->fHistogramsOn)
      continue;

    for (Int_t i=0; i<2; i++) {
      
      fhNClustersITS[i]      ->Add(entry->fhNClustersITS[i]     );      
      fhNClustersTPC[i]      ->Add(entry->fhNClustersTPC[i]     ); 
      if (fhNSharedClustersTPC[i])
        fhNSharedClustersTPC[i]      ->Add(entry->fhNSharedClustersTPC[i]     ); 
      if (fhNCrossedRowsTPC[i])
        fhNCrossedRowsTPC[i]   ->Add(entry->fhNCrossedRowsTPC[i]     ); 
      if (fhRatioCrossedRowsOverFindableClustersTPC[i])
        fhRatioCrossedRowsOverFindableClustersTPC[i]      ->Add(entry->fhRatioCrossedRowsOverFindableClustersTPC[i]     ); 
      					  			    
      fhChi2PerClusterITS[i] ->Add(entry->fhChi2PerClusterITS[i]); 
      fhChi2PerClusterTPC[i] ->Add(entry->fhChi2PerClusterTPC[i]); 
      if (fhChi2TPCConstrainedVsGlobal[i])
	fhChi2TPCConstrainedVsGlobal[i]->Add(entry->fhChi2TPCConstrainedVsGlobal[i]);
      if (fhNClustersForITSPID[i])
	fhNClustersForITSPID[i]->Add(entry->fhNClustersForITSPID[i]);
      if (fhNMissingITSPoints[i])
	fhNMissingITSPoints[i] ->Add(entry->fhNMissingITSPoints[i]);

      fhC11[i]               ->Add(entry->fhC11[i]              ); 
      fhC22[i]               ->Add(entry->fhC22[i]              ); 
      fhC33[i]               ->Add(entry->fhC33[i]              ); 
      fhC44[i]               ->Add(entry->fhC44[i]              ); 
      fhC55[i]               ->Add(entry->fhC55[i]              );

      fhRel1PtUncertainty[i] ->Add(entry->fhRel1PtUncertainty[i]);
      					  			    
      fhDXY[i]               ->Add(entry->fhDXY[i]              ); 
      fhDZ[i]                ->Add(entry->fhDZ[i]               ); 
      fhDXYDZ[i]             ->Add(entry->fhDXYDZ[i]          );
      fhDXYvsDZ[i]           ->Add(entry->fhDXYvsDZ[i]          );

      fhDXYNormalized[i]     ->Add(entry->fhDXYNormalized[i]    );
      fhDZNormalized[i]      ->Add(entry->fhDZNormalized[i]     );
      fhDXYvsDZNormalized[i] ->Add(entry->fhDXYvsDZNormalized[i]);
      fhNSigmaToVertex[i]    ->Add(entry->fhNSigmaToVertex[i]); 

      fhPt[i]                ->Add(entry->fhPt[i]); 
      fhEta[i]               ->Add(entry->fhEta[i]); 
    }      

    fhCutStatistics  ->Add(entry->fhCutStatistics);        
    fhCutCorrelation ->Add(entry->fhCutCorrelation);      

    count++;
  }
  return count+1;
}

void AliESDtrackCuts::SetMinNClustersTPCPtDep(TFormula *f1, Float_t ptmax) 
{
  //
  // Sets the pT dependent NClustersTPC cut
  //

  if(f1){ 
    delete f1CutMinNClustersTPCPtDep;
    f1CutMinNClustersTPCPtDep = (TFormula*)f1->Clone("f1CutMinNClustersTPCPtDep"); 
  }
  fCutMaxPtDepNClustersTPC=ptmax; 
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardTPCOnlyTrackCuts()
{
  // creates an AliESDtrackCuts object and fills it with standard (pre data-taking) values for TPC-only cuts
  
  AliInfoClass("Creating track cuts for TPC-only.");
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  
  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  
  esdTrackCuts->SetMaxDCAToVertexZ(3.2);
  esdTrackCuts->SetMaxDCAToVertexXY(2.4);
  esdTrackCuts->SetDCAToVertex2D(kTRUE);
  
  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS-TPC cuts for pp 2009 data
  
  AliInfoClass("Creating track cuts for ITS+TPC (2009 definition).");
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;

  // TPC  
  esdTrackCuts->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  if(selPrimaries) {
    // 7*(0.0050+0.0060/pt^0.9)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
  }
  esdTrackCuts->SetMaxDCAToVertexZ(1.e6);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //esdTrackCuts->SetEtaRange(-0.8,+0.8);
  
  esdTrackCuts->SetMaxChi2PerClusterITS(36);
  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  
  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(Bool_t selPrimaries,Int_t clusterCut)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS-TPC cuts for pp 2010 data  
  // if clusterCut = 1, the cut on the number of clusters is replaced by
  // a cut on the number of crossed rows and on the ration crossed
  // rows/findable clusters

  AliInfoClass("Creating track cuts for ITS+TPC (2010 definition).");
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;

  // TPC  
  if(clusterCut == 0)  esdTrackCuts->SetMinNClustersTPC(70);
  else if (clusterCut == 1) {
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  }
  else {
    AliWarningClass(Form("Wrong value of the clusterCut parameter (%d), using cut on Nclusters",clusterCut));
    esdTrackCuts->SetMinNClustersTPC(70);
  }
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  if(selPrimaries) {
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  }
  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  
  esdTrackCuts->SetMaxChi2PerClusterITS(36);
  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);

  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSPureSATrackCuts2009(Bool_t selPrimaries, Bool_t useForPid)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS pure SA tracks
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetRequireITSPureStandAlone(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE); 
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2PerClusterITS(1.);

  if(selPrimaries) {
    // 7*(0.0085+0.0026/pt^1.55)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55");
  }
  if(useForPid){
    esdTrackCuts->SetRequireITSPid(kTRUE);
  }
  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(Bool_t selPrimaries, Bool_t useForPid)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS pure SA tracks - pp 2010
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetRequireITSPureStandAlone(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE); 
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2PerClusterITS(2.5);

  if(selPrimaries) {
    // 7*(0.0033+0.0045/pt^1.3)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0231+0.0315/pt^1.3");
  }
  if(useForPid){
    esdTrackCuts->SetRequireITSPid(kTRUE);
  }
  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSSATrackCuts2009(Bool_t selPrimaries, Bool_t useForPid)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS pure SA tracks
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetRequireITSStandAlone(kTRUE);
  esdTrackCuts->SetRequireITSPureStandAlone(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE); 
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2PerClusterITS(1.);

  if(selPrimaries) {
    // 7*(0.0085+0.0026/pt^1.55)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55");
  }
  if(useForPid){
    esdTrackCuts->SetRequireITSPid(kTRUE);
  }
  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSSATrackCuts2010(Bool_t selPrimaries, Bool_t useForPid)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS pure SA tracks --pp 2010
  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetRequireITSStandAlone(kTRUE);
  esdTrackCuts->SetRequireITSPureStandAlone(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE); 
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2PerClusterITS(2.5);

  if(selPrimaries) {
    // 7*(0.0033+0.0045/pt^1.3)
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0231+0.0315/pt^1.3");
  }
  if(useForPid){
    esdTrackCuts->SetRequireITSPid(kTRUE);
  }
  return esdTrackCuts;
}

//____________________________________________________________________
AliESDtrackCuts* AliESDtrackCuts::GetStandardITSSATrackCutsPbPb2010(Bool_t selPrimaries, Bool_t useForPid)
{
  // creates an AliESDtrackCuts object and fills it with standard values for ITS pure SA tracks -- PbPb 2010
  
  AliESDtrackCuts* esdTrackCuts = GetStandardITSSATrackCuts2010(selPrimaries, useForPid);
  esdTrackCuts->SetMaxNOfMissingITSPoints(1);

  return esdTrackCuts;
}
//____________________________________________________________________

AliESDtrackCuts* AliESDtrackCuts::GetStandardV0DaughterCuts()
{
  // creates a AliESDtrackCuts object and fills it with standard cuts for V0 daughters
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  return esdTrackCuts;
}

//____________________________________________________________________
Int_t AliESDtrackCuts::GetReferenceMultiplicity(const AliESDEvent* esd, Bool_t tpcOnly)
{
  // Gets reference multiplicity following the standard cuts and a defined fiducial volume
  // tpcOnly = kTRUE -> consider TPC-only tracks
  //         = kFALSE -> consider global tracks
  //
  // DEPRECATED Use GetReferenceMultiplicity with the enum as second argument instead
  
  if (!tpcOnly)
  {
    AliErrorClass("Not implemented for global tracks!");
    return -1;
  }
  
  static AliESDtrackCuts* esdTrackCuts = 0;
  if (!esdTrackCuts)
  {
    esdTrackCuts = GetStandardTPCOnlyTrackCuts();
    esdTrackCuts->SetEtaRange(-0.8, 0.8);
    esdTrackCuts->SetPtRange(0.15);
  }
  
  Int_t nTracks = esdTrackCuts->CountAcceptedTracks(esd);
  
  return nTracks;
}

//____________________________________________________________________
Float_t AliESDtrackCuts::GetSigmaToVertex(const AliESDtrack* const esdTrack)
{
  // Calculates the number of sigma to the vertex.

  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebugClass(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-d**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0)
    return -1;

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // work around precision problem
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  // 1e-15 corresponds to nsigma ~ 7.7
  if (TMath::Exp(-d * d / 2) < 1e-15)
    return 1000;

  Float_t nSigma = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return nSigma;
}

void AliESDtrackCuts::EnableNeededBranches(TTree* tree)
{
  // enables the branches needed by AcceptTrack, for a list see comment of AcceptTrack

  tree->SetBranchStatus("fTracks.fFlags", 1);
  tree->SetBranchStatus("fTracks.fITSncls", 1);
  tree->SetBranchStatus("fTracks.fTPCncls", 1);
  tree->SetBranchStatus("fTracks.fITSchi2", 1);
  tree->SetBranchStatus("fTracks.fTPCchi2", 1);
  tree->SetBranchStatus("fTracks.fC*", 1);
  tree->SetBranchStatus("fTracks.fD", 1);
  tree->SetBranchStatus("fTracks.fZ", 1);
  tree->SetBranchStatus("fTracks.fCdd", 1);
  tree->SetBranchStatus("fTracks.fCdz", 1);
  tree->SetBranchStatus("fTracks.fCzz", 1);
  tree->SetBranchStatus("fTracks.fP*", 1);
  tree->SetBranchStatus("fTracks.fR*", 1);
  tree->SetBranchStatus("fTracks.fKinkIndexes*", 1);
}

//____________________________________________________________________
Bool_t AliESDtrackCuts::AcceptTrack(const AliESDtrack* esdTrack, const AliESDEvent* esdEvent) 
{
  // 
  // figure out if the tracks survives all the track cuts defined
  //
  // the different quality parameter and kinematic values are first
  // retrieved from the track. then it is found out what cuts the
  // track did not survive and finally the cuts are imposed.

  // this function needs the following branches:
  // fTracks.fFlags
  // fTracks.fITSncls
  // fTracks.fTPCncls
  // fTracks.fITSchi2
  // fTracks.fTPCchi2
  // fTracks.fC   //GetExternalCovariance
  // fTracks.fD   //GetImpactParameters
  // fTracks.fZ   //GetImpactParameters
  // fTracks.fCdd //GetImpactParameters
  // fTracks.fCdz //GetImpactParameters
  // fTracks.fCzz //GetImpactParameters
  // fTracks.fP   //GetPxPyPz
  // fTracks.fR   //GetMass
  // fTracks.fP   //GetMass
  // fTracks.fKinkIndexes
  //
  // esdEvent is only required for the MaxChi2TPCConstrainedVsGlobal

  UInt_t status = esdTrack->GetStatus();

  // getting quality parameters from the ESD track
  Int_t nClustersITS = esdTrack->GetITSclusters(0);
  Int_t nClustersTPC = -1;
  if(fCutRequireTPCStandAlone) {
    nClustersTPC = esdTrack->GetTPCNclsIter1();
  }
  else {
    nClustersTPC = esdTrack->GetTPCclusters(0);
  }

  //Pt dependent NClusters Cut
  if(f1CutMinNClustersTPCPtDep) {
    if(esdTrack->Pt()<fCutMaxPtDepNClustersTPC)
      fCutMinNClusterTPC = f1CutMinNClustersTPCPtDep->Eval(esdTrack->Pt());
    else
      fCutMinNClusterTPC = f1CutMinNClustersTPCPtDep->Eval(fCutMaxPtDepNClustersTPC);
  }

  Float_t nCrossedRowsTPC = esdTrack->GetTPCCrossedRows();
  Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  if (esdTrack->GetTPCNclsF()>0) {
    ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / esdTrack->GetTPCNclsF();
  }
  
  Int_t nClustersTPCShared = esdTrack->GetTPCnclsS();
  Float_t fracClustersTPCShared = -1.;

  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = esdTrack->GetITSchi2()/Float_t(nClustersITS);
  if (nClustersTPC!=0) {
    if(fCutRequireTPCStandAlone) {
      chi2PerClusterTPC = esdTrack->GetTPCchi2Iter1()/Float_t(nClustersTPC);
    } else {
      chi2PerClusterTPC = esdTrack->GetTPCchi2()/Float_t(nClustersTPC);
    }
    fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
  }

  Double_t extCov[15];
  esdTrack->GetExternalCovariance(extCov);

  Float_t b[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }


  // set pt-dependent DCA cuts, if requested
  SetPtDepDCACuts(esdTrack->Pt());


  Float_t dcaToVertexXY = b[0];
  Float_t dcaToVertexZ = b[1];

  Float_t dcaToVertex = -1;
  
  if (fCutDCAToVertex2D)
  {
    dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY/fCutMaxDCAToVertexXY/fCutMaxDCAToVertexXY + dcaToVertexZ*dcaToVertexZ/fCutMaxDCAToVertexZ/fCutMaxDCAToVertexZ);
  }
  else
    dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY + dcaToVertexZ*dcaToVertexZ);
    
  // getting the kinematic variables of the track
  // (assuming the mass is known)
  Double_t p[3];
  esdTrack->GetPxPyPz(p);

  Float_t momentum = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2) + TMath::Power(p[2],2));
  Float_t pt       = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2));
  Float_t energy   = TMath::Sqrt(TMath::Power(esdTrack->GetMass(),2) + TMath::Power(momentum,2));

  //y-eta related calculations
  Float_t eta = -100.;
  Float_t y   = -100.;
  if((momentum != TMath::Abs(p[2]))&&(momentum != 0))
    eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  if((energy != TMath::Abs(p[2]))&&(momentum != 0))
    y = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));
    
  if (extCov[14] < 0) 
  {
    AliWarning(Form("GetSigma1Pt2() returns negative value for external covariance matrix element fC[14]: %f. Corrupted track information, track will not be accepted!", extCov[14]));
    return kFALSE;
  }
  Float_t relUncertainty1Pt = TMath::Sqrt(extCov[14])*pt;
  
  //########################################################################
  // cut the track?
  
  Bool_t cuts[kNCuts];
  for (Int_t i=0; i<kNCuts; i++) cuts[i]=kFALSE;
  
  // track quality cuts
  if (fCutRequireTPCRefit && (status&AliESDtrack::kTPCrefit)==0)
    cuts[0]=kTRUE;
  if (fCutRequireTPCStandAlone && (status&AliESDtrack::kTPCin)==0)
    cuts[1]=kTRUE;
  if (fCutRequireITSRefit && (status&AliESDtrack::kITSrefit)==0)
    cuts[2]=kTRUE;
  if (nClustersTPC<fCutMinNClusterTPC)
    cuts[3]=kTRUE;
  if (nClustersITS<fCutMinNClusterITS) 
    cuts[4]=kTRUE;
  if (chi2PerClusterTPC>fCutMaxChi2PerClusterTPC) 
    cuts[5]=kTRUE; 
  if (chi2PerClusterITS>fCutMaxChi2PerClusterITS) 
    cuts[6]=kTRUE;
  if (extCov[0]  > fCutMaxC11) 
    cuts[7]=kTRUE;  
  if (extCov[2]  > fCutMaxC22) 
    cuts[8]=kTRUE;  
  if (extCov[5]  > fCutMaxC33) 
    cuts[9]=kTRUE;  
  if (extCov[9]  > fCutMaxC44) 
    cuts[10]=kTRUE;  
  if (extCov[14]  > fCutMaxC55) 
    cuts[11]=kTRUE;  
  
  // cut 12 and 13 see below
  
  if (!fCutAcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0)
    cuts[14]=kTRUE;
  // track kinematics cut
  if((momentum < fPMin) || (momentum > fPMax)) 
    cuts[15]=kTRUE;
  if((pt < fPtMin) || (pt > fPtMax)) 
    cuts[16] = kTRUE;
  if((p[0] < fPxMin) || (p[0] > fPxMax)) 
    cuts[17] = kTRUE;
  if((p[1] < fPyMin) || (p[1] > fPyMax)) 
    cuts[18] = kTRUE;
  if((p[2] < fPzMin) || (p[2] > fPzMax))
    cuts[19] = kTRUE;
  if((eta < fEtaMin) || (eta > fEtaMax))
    cuts[20] = kTRUE;
  if((y < fRapMin) || (y > fRapMax)) 
    cuts[21] = kTRUE;
  if (fCutDCAToVertex2D && dcaToVertex > 1)
    cuts[22] = kTRUE;
  if (!fCutDCAToVertex2D && TMath::Abs(dcaToVertexXY) > fCutMaxDCAToVertexXY)
    cuts[23] = kTRUE;
  if (!fCutDCAToVertex2D && TMath::Abs(dcaToVertexZ) > fCutMaxDCAToVertexZ)
    cuts[24] = kTRUE;
  if (fCutDCAToVertex2D && fCutMinDCAToVertexXY > 0 && fCutMinDCAToVertexZ > 0 && dcaToVertexXY*dcaToVertexXY/fCutMinDCAToVertexXY/fCutMinDCAToVertexXY + dcaToVertexZ*dcaToVertexZ/fCutMinDCAToVertexZ/fCutMinDCAToVertexZ < 1)
    cuts[25] = kTRUE;
  if (!fCutDCAToVertex2D && TMath::Abs(dcaToVertexXY) < fCutMinDCAToVertexXY)
    cuts[26] = kTRUE;
  if (!fCutDCAToVertex2D && TMath::Abs(dcaToVertexZ) < fCutMinDCAToVertexZ)
    cuts[27] = kTRUE;
  
  for (Int_t i = 0; i < 3; i++)
    cuts[28+i] = !CheckITSClusterRequirement(fCutClusterRequirementITS[i], esdTrack->HasPointOnITSLayer(i*2), esdTrack->HasPointOnITSLayer(i*2+1));
  
  if(fCutRequireITSStandAlone || fCutRequireITSpureSA){
    if ((status & AliESDtrack::kITSin) == 0 || (status & AliESDtrack::kTPCin)){
      // TPC tracks
      cuts[31] = kTRUE; 
    }else{
      // ITS standalone tracks
      if(fCutRequireITSStandAlone && !fCutRequireITSpureSA){
	if(status & AliESDtrack::kITSpureSA) cuts[31] = kTRUE;
      }else if(fCutRequireITSpureSA){
	if(!(status & AliESDtrack::kITSpureSA)) cuts[31] = kTRUE;
      }
    }
  }

  if (relUncertainty1Pt > fCutMaxRel1PtUncertainty)
     cuts[32] = kTRUE;

  if (!fCutAcceptSharedTPCClusters && nClustersTPCShared!=0)
    cuts[33] = kTRUE;

  if (fracClustersTPCShared > fCutMaxFractionSharedTPCClusters)
    cuts[34] = kTRUE;  

  Int_t nITSPointsForPid=0;
  UChar_t clumap=esdTrack->GetITSClusterMap();
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nITSPointsForPid;
  }
  if(fCutRequireITSPid && nITSPointsForPid<3) cuts[35] = kTRUE;
  

  if (nCrossedRowsTPC<fCutMinNCrossedRowsTPC)
    cuts[36]=kTRUE;
  if (ratioCrossedRowsOverFindableClustersTPC<fCutMinRatioCrossedRowsOverFindableClustersTPC) 
    cuts[37]=kTRUE;

  Int_t nMissITSpts=0;
  Int_t idet,statusLay;
  Float_t xloc,zloc;
  for(Int_t iLay=0; iLay<6; iLay++){
    Bool_t retc=esdTrack->GetITSModuleIndexInfo(iLay,idet,statusLay,xloc,zloc);
    if(retc && statusLay==5) ++nMissITSpts;
  }
  if(nMissITSpts>fCutMaxMissingITSPoints) cuts[38] = kTRUE;
  
  Bool_t cut=kFALSE;
  for (Int_t i=0; i<kNCuts; i++) 
    if (cuts[i]) {cut = kTRUE;}

  // for performance evaluate the CPU intensive cuts only when the others have passed, and when they are requested
  Double_t chi2TPCConstrainedVsGlobal = -2;
  Float_t nSigmaToVertex = -2;
  if (!cut)
  {
    // getting the track to vertex parameters
    if (fCutSigmaToVertexRequired)
    {
      nSigmaToVertex = GetSigmaToVertex(esdTrack);
      if (nSigmaToVertex > fCutNsigmaToVertex && fCutSigmaToVertexRequired)
      {
	cuts[12] = kTRUE;
	cut = kTRUE;
      }
      // if n sigma could not be calculated
      if (nSigmaToVertex<0 && fCutSigmaToVertexRequired)
      {
	cuts[13] = kTRUE;
	cut = kTRUE;
      }
    }
      
    // max chi2 TPC constrained vs global track only if track passed the other cut
    if (fCutMaxChi2TPCConstrainedVsGlobal < 1e9)
    {
      if (!esdEvent)
	AliFatal("fCutMaxChi2TPCConstrainedVsGlobal set but ESD event not provided.");
      
      // get vertex
      const AliESDVertex* vertex = 0;
      if (fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexTracks)
	vertex = esdEvent->GetPrimaryVertexTracks();
      
      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexSPD)
	vertex = esdEvent->GetPrimaryVertexSPD();
	
      if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexTPC)
	vertex = esdEvent->GetPrimaryVertexTPC();

      if (vertex->GetStatus())
	chi2TPCConstrainedVsGlobal = esdTrack->GetChi2TPCConstrainedVsGlobal(vertex);
      
      if (chi2TPCConstrainedVsGlobal < 0 || chi2TPCConstrainedVsGlobal > fCutMaxChi2TPCConstrainedVsGlobal)
      {
	cuts[39] = kTRUE;
	cut = kTRUE;
      }
    }
  }
  
  //########################################################################
  // filling histograms
  if (fHistogramsOn) {
    fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n tracks")));
    if (cut)
      fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n cut tracks")));
    
    for (Int_t i=0; i<kNCuts; i++) {
      if (fhCutStatistics->GetXaxis()->FindBin(fgkCutNames[i]) < 1)
        AliFatal(Form("Inconsistency! Cut %d with name %s not found", i, fgkCutNames[i]));
    
      if (cuts[i])
        fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin(fgkCutNames[i])));

      for (Int_t j=i; j<kNCuts; j++) {
        if (cuts[i] && cuts[j]) {
          Float_t xC = fhCutCorrelation->GetXaxis()->GetBinCenter(fhCutCorrelation->GetXaxis()->FindBin(fgkCutNames[i]));
          Float_t yC = fhCutCorrelation->GetYaxis()->GetBinCenter(fhCutCorrelation->GetYaxis()->FindBin(fgkCutNames[j]));
          fhCutCorrelation->Fill(xC, yC);
        }
      }
    }
  }
  
  // now we loop over the filling of the histograms twice: once "before" the cut, once "after"
  // the code is not in a function due to too many local variables that would need to be passed

  for (Int_t id = 0; id < 2; id++)
  {
    // id = 0 --> before cut
    // id = 1 --> after cut

    if (fHistogramsOn)
    {
      fhNClustersITS[id]->Fill(nClustersITS);
      fhNClustersTPC[id]->Fill(nClustersTPC);
      fhNSharedClustersTPC[id]->Fill(nClustersTPCShared);
      fhNCrossedRowsTPC[id]->Fill(nCrossedRowsTPC);
      fhRatioCrossedRowsOverFindableClustersTPC[id]->Fill(ratioCrossedRowsOverFindableClustersTPC);
      fhChi2PerClusterITS[id]->Fill(chi2PerClusterITS);
      fhChi2PerClusterTPC[id]->Fill(chi2PerClusterTPC);
      fhChi2TPCConstrainedVsGlobal[id]->Fill(chi2TPCConstrainedVsGlobal);
      fhNClustersForITSPID[id]->Fill(nITSPointsForPid);
      fhNMissingITSPoints[id]->Fill(nMissITSpts);

      fhC11[id]->Fill(extCov[0]);
      fhC22[id]->Fill(extCov[2]);
      fhC33[id]->Fill(extCov[5]);
      fhC44[id]->Fill(extCov[9]);
      fhC55[id]->Fill(extCov[14]);

      fhRel1PtUncertainty[id]->Fill(relUncertainty1Pt);

      fhPt[id]->Fill(pt);
      fhEta[id]->Fill(eta);

      Float_t bRes[2];
      bRes[0] = TMath::Sqrt(bCov[0]);
      bRes[1] = TMath::Sqrt(bCov[2]);

      fhDZ[id]->Fill(b[1]);
      fhDXY[id]->Fill(b[0]);
      fhDXYDZ[id]->Fill(dcaToVertex);
      fhDXYvsDZ[id]->Fill(b[1],b[0]);

      if (bRes[0]!=0 && bRes[1]!=0) {
        fhDZNormalized[id]->Fill(b[1]/bRes[1]);
        fhDXYNormalized[id]->Fill(b[0]/bRes[0]);
        fhDXYvsDZNormalized[id]->Fill(b[1]/bRes[1], b[0]/bRes[0]);
        fhNSigmaToVertex[id]->Fill(nSigmaToVertex);
      }
    }

    // cut the track
    if (cut)
      return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________
Bool_t AliESDtrackCuts::CheckITSClusterRequirement(ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2)
{
  // checks if the cluster requirement is fullfilled (in this case: return kTRUE)
  
  switch (req)
  {
  	case kOff:        return kTRUE;
  	case kNone:       return !clusterL1 && !clusterL2;
  	case kAny:        return clusterL1 || clusterL2;
  	case kFirst:      return clusterL1;
  	case kOnlyFirst:  return clusterL1 && !clusterL2;
  	case kSecond:     return clusterL2;
  	case kOnlySecond: return clusterL2 && !clusterL1;
  	case kBoth:       return clusterL1 && clusterL2;
  }
  
  return kFALSE;
}

//____________________________________________________________________
AliESDtrack* AliESDtrackCuts::GetTPCOnlyTrack(const AliESDEvent* esd, Int_t iTrack)
{
  // Utility function to create a TPC only track from the given esd track
  // 
  // IMPORTANT: The track has to be deleted by the user
  //
  // NB. most of the functionality to get a TPC only track from an ESD track is in AliESDtrack, where it should be
  // there are only missing propagations here that are needed for old data
  // this function will therefore become obsolete
  //
  // adapted from code provided by CKB

  if (!esd->GetPrimaryVertexTPC())
    return 0; // No TPC vertex no TPC tracks

  if(!esd->GetPrimaryVertexTPC()->GetStatus())
    return 0; // TPC Vertex is created by default in AliESDEvent, do not use in this case

  AliESDtrack* track = esd->GetTrack(iTrack);
  if (!track)
    return 0;

  AliESDtrack *tpcTrack = new AliESDtrack();

  // only true if we have a tpc track
  if (!track->FillTPCOnlyTrack(*tpcTrack))
  {
    delete tpcTrack;
    return 0;
  }

  return tpcTrack;
}

//____________________________________________________________________
TObjArray* AliESDtrackCuts::GetAcceptedTracks(const AliESDEvent* esd, Bool_t bTPC)
{
  //
  // returns an array of all tracks that pass the cuts
  // or an array of TPC only tracks (propagated to the TPC vertex during reco)
  // tracks that pass the cut
  //
  // NOTE: List has to be deleted by the user

  TObjArray* acceptedTracks = new TObjArray();

  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    if(bTPC){
      if(!esd->GetPrimaryVertexTPC())return acceptedTracks; // No TPC vertex no TPC tracks
      if(!esd->GetPrimaryVertexTPC()->GetStatus())return acceptedTracks; // No proper TPC vertex, only the default

      AliESDtrack *tpcTrack = GetTPCOnlyTrack(esd, iTrack);
      if (!tpcTrack)
        continue;

      if (AcceptTrack(tpcTrack, esd)) {
        acceptedTracks->Add(tpcTrack);
      }
      else
        delete tpcTrack;
    }
    else
    {
      AliESDtrack* track = esd->GetTrack(iTrack);
      if(AcceptTrack(track, esd))
        acceptedTracks->Add(track);
    }
  } 
  if(bTPC)acceptedTracks->SetOwner(kTRUE);
  return acceptedTracks;
}

//____________________________________________________________________
Int_t AliESDtrackCuts::CountAcceptedTracks(const AliESDEvent* const esd)
{
  //
  // returns an the number of tracks that pass the cuts
  //

  Int_t count = 0;

  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = esd->GetTrack(iTrack);
    if (AcceptTrack(track, esd))
      count++;
  }

  return count;
}

//____________________________________________________________________
 void AliESDtrackCuts::DefineHistograms(Int_t color) {
   // 
   // diagnostics histograms are defined
   // 

   fHistogramsOn=kTRUE;
   
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   
   //###################################################################################
   // defining histograms

   fhCutStatistics = new TH1F("cut_statistics","cut statistics",kNCuts+4,-0.5,kNCuts+3.5);

   fhCutStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
   fhCutStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");

   fhCutCorrelation = new TH2F("cut_correlation","cut correlation",kNCuts,-0.5,kNCuts-0.5,kNCuts,-0.5,kNCuts-0.5);;
  
   for (Int_t i=0; i<kNCuts; i++) {
     fhCutStatistics->GetXaxis()->SetBinLabel(i+4,fgkCutNames[i]);
     fhCutCorrelation->GetXaxis()->SetBinLabel(i+1,fgkCutNames[i]);
     fhCutCorrelation->GetYaxis()->SetBinLabel(i+1,fgkCutNames[i]);
   }

  fhCutStatistics  ->SetLineColor(color);
  fhCutCorrelation ->SetLineColor(color);
  fhCutStatistics  ->SetLineWidth(2);
  fhCutCorrelation ->SetLineWidth(2);

  for (Int_t i=0; i<2; i++) {
    fhNClustersITS[i]        = new TH1F("nClustersITS"    ,"",8,-0.5,7.5);
    fhNClustersTPC[i]        = new TH1F("nClustersTPC"     ,"",165,-0.5,164.5);
    fhNSharedClustersTPC[i]  = new TH1F("nSharedClustersTPC"     ,"",165,-0.5,164.5);
    fhNCrossedRowsTPC[i]     = new TH1F("nCrossedRowsTPC"  ,"",165,-0.5,164.5);
    fhRatioCrossedRowsOverFindableClustersTPC[i]     = new TH1F("ratioCrossedRowsOverFindableClustersTPC"  ,"",60,0,1.5);
    fhChi2PerClusterITS[i]   = new TH1F("chi2PerClusterITS","",500,0,10);
    fhChi2PerClusterTPC[i]   = new TH1F("chi2PerClusterTPC","",500,0,10);
    fhChi2TPCConstrainedVsGlobal[i]   = new TH1F("chi2TPCConstrainedVsGlobal","",600,-2,50);
    fhNClustersForITSPID[i]  = new TH1F("nPointsForITSpid","",5,-0.5,4.5);
    fhNMissingITSPoints[i]   = new TH1F("nMissingITSClusters","",7,-0.5,6.5);

    fhC11[i]                 = new TH1F("covMatrixDiagonal11","",2000,0,20);
    fhC22[i]                 = new TH1F("covMatrixDiagonal22","",2000,0,20);
    fhC33[i]                 = new TH1F("covMatrixDiagonal33","",1000,0,0.1);
    fhC44[i]                 = new TH1F("covMatrixDiagonal44","",1000,0,0.1);
    fhC55[i]                 = new TH1F("covMatrixDiagonal55","",1000,0,5);

    fhRel1PtUncertainty[i]   = new TH1F("rel1PtUncertainty","",1000,0,5);

    fhDXY[i]                 = new TH1F("dXY"    ,"",500,-10,10);
    fhDZ[i]                  = new TH1F("dZ"     ,"",500,-10,10);
    fhDXYDZ[i]               = new TH1F("dXYDZ"  ,"",500,0,10);
    fhDXYvsDZ[i]             = new TH2F("dXYvsDZ","",200,-10,10,200,-10,10);

    fhDXYNormalized[i]       = new TH1F("dXYNormalized"    ,"",500,-10,10);
    fhDZNormalized[i]        = new TH1F("dZNormalized"     ,"",500,-10,10);
    fhDXYvsDZNormalized[i]   = new TH2F("dXYvsDZNormalized","",200,-10,10,200,-10,10);

    fhNSigmaToVertex[i]      = new TH1F("nSigmaToVertex","",500,0,10);

    fhPt[i]                  = new TH1F("pt"     ,"p_{T} distribution;p_{T} (GeV/c)", 800, 0.0, 10.0);
    fhEta[i]                 = new TH1F("eta"     ,"#eta distribution;#eta",40,-2.0,2.0);
    
    fhNClustersITS[i]->SetTitle("n ITS clusters");
    fhNClustersTPC[i]->SetTitle("n TPC clusters");
    fhNSharedClustersTPC[i]->SetTitle("n TPC shared clusters");
    fhChi2PerClusterITS[i]->SetTitle("#Chi^{2} per ITS cluster");
    fhChi2PerClusterTPC[i]->SetTitle("#Chi^{2} per TPC cluster");
    fhChi2TPCConstrainedVsGlobal[i]->SetTitle("#Chi^{2} TPC constrained track vs global track");
    fhNClustersForITSPID[i]->SetTitle("n ITS points for PID");
    fhNMissingITSPoints[i]->SetTitle("n ITS layers with missing cluster");

    fhC11[i]->SetTitle("cov 11 : #sigma_{y}^{2} [cm^{2}]");
    fhC22[i]->SetTitle("cov 22 : #sigma_{z}^{2} [cm^{2}]");
    fhC33[i]->SetTitle("cov 33 : #sigma_{sin(#phi)}^{2}");
    fhC44[i]->SetTitle("cov 44 : #sigma_{tan(#theta_{dip})}^{2}");
    fhC55[i]->SetTitle("cov 55 : #sigma_{1/p_{T}}^{2} [(c/GeV)^2]");

    fhRel1PtUncertainty[i]->SetTitle("rel. uncertainty of 1/p_{T}");

    fhDXY[i]->SetXTitle("transverse impact parameter (cm)");
    fhDZ[i]->SetXTitle("longitudinal impact parameter (cm)");
    fhDXYDZ[i]->SetTitle("absolute impact parameter;sqrt(dXY**2 + dZ**2)  (cm)");
    fhDXYvsDZ[i]->SetXTitle("longitudinal impact parameter (cm)");
    fhDXYvsDZ[i]->SetYTitle("transverse impact parameter (cm)");

    fhDXYNormalized[i]->SetTitle("normalized trans impact par (n#sigma)");
    fhDZNormalized[i]->SetTitle("normalized long impact par (n#sigma)");
    fhDXYvsDZNormalized[i]->SetTitle("normalized long impact par (n#sigma)");
    fhDXYvsDZNormalized[i]->SetYTitle("normalized trans impact par (n#sigma)");
    fhNSigmaToVertex[i]->SetTitle("n #sigma to vertex");

    fhNClustersITS[i]->SetLineColor(color);   fhNClustersITS[i]->SetLineWidth(2);
    fhNClustersTPC[i]->SetLineColor(color);   fhNClustersTPC[i]->SetLineWidth(2);
    fhNSharedClustersTPC[i]->SetLineColor(color);  fhNSharedClustersTPC[i]->SetLineWidth(2);
    fhChi2PerClusterITS[i]->SetLineColor(color);   fhChi2PerClusterITS[i]->SetLineWidth(2);
    fhChi2PerClusterTPC[i]->SetLineColor(color);   fhChi2PerClusterTPC[i]->SetLineWidth(2);
    fhChi2TPCConstrainedVsGlobal[i]->SetLineColor(color);   fhChi2TPCConstrainedVsGlobal[i]->SetLineWidth(2);
    fhNClustersForITSPID[i]->SetLineColor(color);  fhNClustersForITSPID[i]->SetLineWidth(2);
    fhNMissingITSPoints[i]->SetLineColor(color);   fhNMissingITSPoints[i]->SetLineWidth(2);

    fhC11[i]->SetLineColor(color);   fhC11[i]->SetLineWidth(2);
    fhC22[i]->SetLineColor(color);   fhC22[i]->SetLineWidth(2);
    fhC33[i]->SetLineColor(color);   fhC33[i]->SetLineWidth(2);
    fhC44[i]->SetLineColor(color);   fhC44[i]->SetLineWidth(2);
    fhC55[i]->SetLineColor(color);   fhC55[i]->SetLineWidth(2);

    fhRel1PtUncertainty[i]->SetLineColor(color); fhRel1PtUncertainty[i]->SetLineWidth(2);

    fhDXY[i]->SetLineColor(color);   fhDXY[i]->SetLineWidth(2);
    fhDZ[i]->SetLineColor(color);    fhDZ[i]->SetLineWidth(2);
    fhDXYDZ[i]->SetLineColor(color); fhDXYDZ[i]->SetLineWidth(2);

    fhDXYNormalized[i]->SetLineColor(color);   fhDXYNormalized[i]->SetLineWidth(2);
    fhDZNormalized[i]->SetLineColor(color);    fhDZNormalized[i]->SetLineWidth(2);
    fhNSigmaToVertex[i]->SetLineColor(color);  fhNSigmaToVertex[i]->SetLineWidth(2); 
  }

  // The number of sigmas to the vertex is per definition gaussian
  ffDTheoretical = new TF1("nSigmaToVertexTheoretical","([0]/2.506628274)*exp(-(x**2)/2)",0,50);
  ffDTheoretical->SetParameter(0,1);

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
Bool_t AliESDtrackCuts::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms from a file
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  ffDTheoretical = dynamic_cast<TF1*> (gDirectory->Get("nSigmaToVertexTheory"));

  fhCutStatistics = dynamic_cast<TH1F*> (gDirectory->Get("cut_statistics"));
  fhCutCorrelation = dynamic_cast<TH2F*> (gDirectory->Get("cut_correlation"));

  for (Int_t i=0; i<2; i++) {
    if (i==0)
    {
      gDirectory->cd("before_cuts");
    }
    else
      gDirectory->cd("after_cuts");

    fhNClustersITS[i]      = dynamic_cast<TH1F*> (gDirectory->Get("nClustersITS"     ));
    fhNClustersTPC[i]      = dynamic_cast<TH1F*> (gDirectory->Get("nClustersTPC"     ));
    fhNSharedClustersTPC[i] = dynamic_cast<TH1F*> (gDirectory->Get("nSharedClustersTPC"     ));
    fhNCrossedRowsTPC[i]   = dynamic_cast<TH1F*> (gDirectory->Get("nCrossedRowsTPC"  ));
    fhRatioCrossedRowsOverFindableClustersTPC[i]   = dynamic_cast<TH1F*> (gDirectory->Get("ratioCrossedRowsOverFindableClustersTPC"  ));
    fhChi2PerClusterITS[i] = dynamic_cast<TH1F*> (gDirectory->Get("chi2PerClusterITS"));
    fhChi2PerClusterTPC[i] = dynamic_cast<TH1F*> (gDirectory->Get("chi2PerClusterTPC"));
    fhChi2TPCConstrainedVsGlobal[i] = dynamic_cast<TH1F*> (gDirectory->Get("fhChi2TPCConstrainedVsGlobal"));
    fhNClustersForITSPID[i] = dynamic_cast<TH1F*> (gDirectory->Get("nPointsForITSpid"));
    fhNMissingITSPoints[i] = dynamic_cast<TH1F*> (gDirectory->Get("nMissingITSClusters"));

    fhC11[i] = dynamic_cast<TH1F*> (gDirectory->Get("covMatrixDiagonal11"));
    fhC22[i] = dynamic_cast<TH1F*> (gDirectory->Get("covMatrixDiagonal22"));
    fhC33[i] = dynamic_cast<TH1F*> (gDirectory->Get("covMatrixDiagonal33"));
    fhC44[i] = dynamic_cast<TH1F*> (gDirectory->Get("covMatrixDiagonal44"));
    fhC55[i] = dynamic_cast<TH1F*> (gDirectory->Get("covMatrixDiagonal55"));

    fhRel1PtUncertainty[i] = dynamic_cast<TH1F*> (gDirectory->Get("rel1PtUncertainty"));

    fhDXY[i] =     dynamic_cast<TH1F*> (gDirectory->Get("dXY"    ));
    fhDZ[i] =      dynamic_cast<TH1F*> (gDirectory->Get("dZ"     ));
    fhDXYDZ[i] =   dynamic_cast<TH1F*> (gDirectory->Get("dXYDZ"));
    fhDXYvsDZ[i] = dynamic_cast<TH2F*> (gDirectory->Get("dXYvsDZ"));

    fhDXYNormalized[i] =     dynamic_cast<TH1F*> (gDirectory->Get("dXYNormalized"    ));
    fhDZNormalized[i] =      dynamic_cast<TH1F*> (gDirectory->Get("dZNormalized"     ));
    fhDXYvsDZNormalized[i] = dynamic_cast<TH2F*> (gDirectory->Get("dXYvsDZNormalized"));
    fhNSigmaToVertex[i] = dynamic_cast<TH1F*> (gDirectory->Get("nSigmaToVertex"));

    fhPt[i] = dynamic_cast<TH1F*> (gDirectory->Get("pt"));
    fhEta[i] = dynamic_cast<TH1F*> (gDirectory->Get("eta"));

    gDirectory->cd("../");
  }

  gDirectory->cd("..");

  return kTRUE;
}

//____________________________________________________________________
void AliESDtrackCuts::SaveHistograms(const Char_t* dir) {
  //
  // saves the histograms in a directory (dir)
  //

  if (!fHistogramsOn) {
    AliDebug(0, "Histograms not on - cannot save histograms!!!");
    return;
  }

  if (!dir)
    dir = GetName();

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");
 
  // a factor of 2 is needed since n sigma is positive
  ffDTheoretical->SetParameter(0,2*fhNSigmaToVertex[0]->Integral("width"));
  ffDTheoretical->Write("nSigmaToVertexTheory");

  fhCutStatistics->Write();
  fhCutCorrelation->Write();

  for (Int_t i=0; i<2; i++) {
    if (i==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");

    fhNClustersITS[i]        ->Write();
    fhNClustersTPC[i]        ->Write();
    fhNSharedClustersTPC[i]  ->Write();
    fhNCrossedRowsTPC[i]     ->Write();
    fhRatioCrossedRowsOverFindableClustersTPC[i]     ->Write();
    fhChi2PerClusterITS[i]   ->Write();
    fhChi2PerClusterTPC[i]   ->Write();
    fhChi2TPCConstrainedVsGlobal[i]   ->Write();
    fhNClustersForITSPID[i]  ->Write();
    fhNMissingITSPoints[i]   ->Write();

    fhC11[i]                 ->Write();
    fhC22[i]                 ->Write();
    fhC33[i]                 ->Write();
    fhC44[i]                 ->Write();
    fhC55[i]                 ->Write();

    fhRel1PtUncertainty[i]   ->Write();

    fhDXY[i]                 ->Write();
    fhDZ[i]                  ->Write();
    fhDXYDZ[i]               ->Write();
    fhDXYvsDZ[i]             ->Write();

    fhDXYNormalized[i]       ->Write();
    fhDZNormalized[i]        ->Write();
    fhDXYvsDZNormalized[i]   ->Write();
    fhNSigmaToVertex[i]      ->Write();

    fhPt[i]                  ->Write();
    fhEta[i]                 ->Write();
    
    gDirectory->cd("../");
  }

  gDirectory->cd("../");
}

//____________________________________________________________________
void AliESDtrackCuts::DrawHistograms()
{
  // draws some histograms

  TCanvas* canvas1 = new TCanvas(Form("%s_1", GetName()), "Track Quality Results1", 800, 800);
  canvas1->Divide(2, 2);

  canvas1->cd(1);
  fhNClustersTPC[0]->SetStats(kFALSE);
  fhNClustersTPC[0]->Draw();

  canvas1->cd(2);
  fhChi2PerClusterTPC[0]->SetStats(kFALSE);
  fhChi2PerClusterTPC[0]->Draw();

  canvas1->cd(3);
  fhNSigmaToVertex[0]->SetStats(kFALSE);
  fhNSigmaToVertex[0]->GetXaxis()->SetRangeUser(0, 10);
  fhNSigmaToVertex[0]->Draw();

  canvas1->SaveAs(Form("%s_%s.gif", GetName(), canvas1->GetName()));

  TCanvas* canvas2 = new TCanvas(Form("%s_2", GetName()), "Track Quality Results2", 1200, 800);
  canvas2->Divide(3, 2);

  canvas2->cd(1);
  fhC11[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhC11[0]->Draw();

  canvas2->cd(2);
  fhC22[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhC22[0]->Draw();

  canvas2->cd(3);
  fhC33[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhC33[0]->Draw();

  canvas2->cd(4);
  fhC44[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhC44[0]->Draw();

  canvas2->cd(5);
  fhC55[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhC55[0]->Draw();

  canvas2->cd(6);
  fhRel1PtUncertainty[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhRel1PtUncertainty[0]->Draw();

  canvas2->SaveAs(Form("%s_%s.gif", GetName(), canvas2->GetName()));

  TCanvas* canvas3 = new TCanvas(Form("%s_3", GetName()), "Track Quality Results3", 1200, 800);
  canvas3->Divide(3, 2);

  canvas3->cd(1);
  fhDXY[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhDXY[0]->Draw();

  canvas3->cd(2);
  fhDZ[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhDZ[0]->Draw();

  canvas3->cd(3);
  fhDXYvsDZ[0]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetRightMargin(0.15);
  fhDXYvsDZ[0]->Draw("COLZ");

  canvas3->cd(4);
  fhDXYNormalized[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhDXYNormalized[0]->Draw();

  canvas3->cd(5);
  fhDZNormalized[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhDZNormalized[0]->Draw();

  canvas3->cd(6);
  fhDXYvsDZNormalized[0]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetRightMargin(0.15);
  fhDXYvsDZNormalized[0]->Draw("COLZ");

  canvas3->SaveAs(Form("%s_%s.gif", GetName(), canvas3->GetName()));

  TCanvas* canvas4 = new TCanvas(Form("%s_4", GetName()), "Track Quality Results4", 800, 500);
  canvas4->Divide(2, 1);

  canvas4->cd(1);
  fhCutStatistics->SetStats(kFALSE);
  fhCutStatistics->LabelsOption("v");
  gPad->SetBottomMargin(0.3);
  fhCutStatistics->Draw();

  canvas4->cd(2);
  fhCutCorrelation->SetStats(kFALSE);
  fhCutCorrelation->LabelsOption("v");
  gPad->SetBottomMargin(0.3);
  gPad->SetLeftMargin(0.3);
  fhCutCorrelation->Draw("COLZ");

  canvas4->SaveAs(Form("%s_%s.gif", GetName(), canvas4->GetName()));

  /*canvas->cd(1);
  fhDXYvsDZNormalized[0]->SetStats(kFALSE);
  fhDXYvsDZNormalized[0]->DrawCopy("COLZ");

  canvas->cd(2);
  fhNClustersTPC[0]->SetStats(kFALSE);
  fhNClustersTPC[0]->DrawCopy();

  canvas->cd(3);
  fhChi2PerClusterITS[0]->SetStats(kFALSE);
  fhChi2PerClusterITS[0]->DrawCopy();
  fhChi2PerClusterITS[1]->SetLineColor(2);
  fhChi2PerClusterITS[1]->DrawCopy("SAME");

  canvas->cd(4);
  fhChi2PerClusterTPC[0]->SetStats(kFALSE);
  fhChi2PerClusterTPC[0]->DrawCopy();
  fhChi2PerClusterTPC[1]->SetLineColor(2);
  fhChi2PerClusterTPC[1]->DrawCopy("SAME");*/
}
//--------------------------------------------------------------------------
void AliESDtrackCuts::SetPtDepDCACuts(Double_t pt) {
  //
  // set the pt-dependent DCA cuts
  //

  if(f1CutMaxDCAToVertexXYPtDep) {
     fCutMaxDCAToVertexXY=f1CutMaxDCAToVertexXYPtDep->Eval(pt);
  }

  if(f1CutMaxDCAToVertexZPtDep) {
    fCutMaxDCAToVertexZ=f1CutMaxDCAToVertexZPtDep->Eval(pt);
  }

  if(f1CutMinDCAToVertexXYPtDep) {
    fCutMinDCAToVertexXY=f1CutMinDCAToVertexXYPtDep->Eval(pt);
  }

  if(f1CutMinDCAToVertexZPtDep) {
    fCutMinDCAToVertexZ=f1CutMinDCAToVertexZPtDep->Eval(pt);
  }


  return;
}



//--------------------------------------------------------------------------
Bool_t AliESDtrackCuts::CheckPtDepDCA(TString dist,Bool_t print) const {
  //
  // Check the correctness of the string syntax
  //
  Bool_t retval=kTRUE;

  if(!dist.Contains("pt")) {
    if(print) AliError("string must contain \"pt\"");
    retval= kFALSE;
  } 
  return retval;
}

 void AliESDtrackCuts::SetMaxDCAToVertexXYPtDep(const char *dist){

   if(f1CutMaxDCAToVertexXYPtDep){
     delete f1CutMaxDCAToVertexXYPtDep;
     // resetiing both
     f1CutMaxDCAToVertexXYPtDep = 0;
     fCutMaxDCAToVertexXYPtDep = "";
   }
   if(!CheckPtDepDCA(dist,kTRUE)){
     return;
   }  
   fCutMaxDCAToVertexXYPtDep = dist;
   TString tmp(dist);
   tmp.ReplaceAll("pt","x");
   f1CutMaxDCAToVertexXYPtDep = new TFormula("f1CutMaxDCAToVertexXYPtDep",tmp.Data());
 
}

 void AliESDtrackCuts::SetMaxDCAToVertexZPtDep(const char *dist){


   if(f1CutMaxDCAToVertexZPtDep){
     delete f1CutMaxDCAToVertexZPtDep;
     // resetiing both
     f1CutMaxDCAToVertexZPtDep = 0;
     fCutMaxDCAToVertexZPtDep = "";
   }
   if(!CheckPtDepDCA(dist,kTRUE))return;
     
   fCutMaxDCAToVertexZPtDep = dist;
   TString tmp(dist);
   tmp.ReplaceAll("pt","x");
   f1CutMaxDCAToVertexZPtDep = new TFormula("f1CutMaxDCAToVertexZPtDep",tmp.Data());

   
}


 void AliESDtrackCuts::SetMinDCAToVertexXYPtDep(const char *dist){


   if(f1CutMinDCAToVertexXYPtDep){
     delete f1CutMinDCAToVertexXYPtDep;
     // resetiing both
     f1CutMinDCAToVertexXYPtDep = 0;
     fCutMinDCAToVertexXYPtDep = "";
   }
   if(!CheckPtDepDCA(dist,kTRUE))return;

   fCutMinDCAToVertexXYPtDep = dist;
   TString tmp(dist);
   tmp.ReplaceAll("pt","x");
   f1CutMinDCAToVertexXYPtDep = new TFormula("f1CutMinDCAToVertexXYPtDep",tmp.Data());

}


 void AliESDtrackCuts::SetMinDCAToVertexZPtDep(const char *dist){

   

   if(f1CutMinDCAToVertexZPtDep){
     delete f1CutMinDCAToVertexZPtDep;
     // resetiing both
     f1CutMinDCAToVertexZPtDep = 0;
     fCutMinDCAToVertexZPtDep = "";
   }
   if(!CheckPtDepDCA(dist,kTRUE))return;
   fCutMinDCAToVertexZPtDep = dist;
   TString tmp(dist);
   tmp.ReplaceAll("pt","x");
   f1CutMinDCAToVertexZPtDep = new TFormula("f1CutMinDCAToVertexZPtDep",tmp.Data());
}

AliESDtrackCuts* AliESDtrackCuts::GetMultEstTrackCuts(MultEstTrackCuts cut)
{
  // returns the multiplicity estimator track cuts objects to allow for user configuration
  // upon first call the objects are created
  //
  // the cut defined here correspond to GetStandardITSTPCTrackCuts2010 (apart from the one for without SPD)
  
  if (!fgMultEstTrackCuts[kMultEstTrackCutGlobal])
  {
    // quality cut on ITS+TPC tracks
    fgMultEstTrackCuts[kMultEstTrackCutGlobal] = new AliESDtrackCuts();
    fgMultEstTrackCuts[kMultEstTrackCutGlobal]->SetMinNClustersTPC(70);
    fgMultEstTrackCuts[kMultEstTrackCutGlobal]->SetMaxChi2PerClusterTPC(4);
    fgMultEstTrackCuts[kMultEstTrackCutGlobal]->SetAcceptKinkDaughters(kFALSE);
    fgMultEstTrackCuts[kMultEstTrackCutGlobal]->SetRequireTPCRefit(kTRUE);
    fgMultEstTrackCuts[kMultEstTrackCutGlobal]->SetRequireITSRefit(kTRUE);
    //multiplicity underestimate if we use global tracks with |eta| > 0.9
    fgMultEstTrackCuts[kMultEstTrackCutGlobal]->SetEtaRange(-0.9, 0.9);
    
    // quality cut on ITS_SA tracks (complementary to ITS+TPC)
    fgMultEstTrackCuts[kMultEstTrackCutITSSA] = new AliESDtrackCuts();
    fgMultEstTrackCuts[kMultEstTrackCutITSSA]->SetRequireITSRefit(kTRUE);
    
    // primary selection for tracks with SPD hits
    fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD] = new AliESDtrackCuts();
    fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD]->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD]->SetMaxDCAToVertexZ(2);
    
    // primary selection for tracks w/o SPD hits
    fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD] = new AliESDtrackCuts();
    fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD]->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
    fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD]->SetMaxDCAToVertexXYPtDep("1.5*(0.0182+0.0350/pt^1.01)");
    fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD]->SetMaxDCAToVertexZ(2);
  }
  
  return fgMultEstTrackCuts[cut];
}

Int_t AliESDtrackCuts::GetReferenceMultiplicity(const AliESDEvent* esd, MultEstTrackType trackType, Float_t etaRange) 
{
  // Get multiplicity estimate based on TPC/ITS tracks and tracklets
  // Adapted for AliESDtrackCuts from a version developed by: Ruben Shahoyan, Anton Alkin, Arvinder Palaha
  //
  // Returns a negative value if no reliable estimate can be provided:
  //   -1 SPD vertex missing
  //   -2 SPD VertexerZ dispersion too large
  //   -3 Track vertex missing (not checked for kTracklets)
  //   -4 Distance between SPD and track vertices too large (not checked for kTracklets)
  //
  // WARNING This functions does not cut on the z vtx. Depending on the eta range requested, you need to restrict your z vertex range!
  //
  // Strategy for combined estimators
  //   1. Count global tracks and flag them
  //   2. Count ITSSA as complementaries for ITSTPC+ or as main tracks
  //   3. Count the complementary SPD tracklets
  
  const AliESDVertex* vertices[2];
  vertices[0] = esd->GetPrimaryVertexSPD();
  vertices[1] = esd->GetPrimaryVertexTracks();
  
  if (!vertices[0]->GetStatus())
  {
    AliDebugClass(AliLog::kDebug, "No SPD vertex. Not able to make a reliable multiplicity estimate.");
    return -1;
  }
  
  if (vertices[0]->IsFromVertexerZ() && vertices[0]->GetDispersion() > 0.02)
  {
    AliDebugClass(AliLog::kDebug, "Vertexer z dispersion > 0.02. Not able to make a reliable multiplicity estimate.");
    return -2;
  }
  
  Int_t multiplicityEstimate = 0;
  
  // SPD tracklet-only estimate
  if (trackType == kTracklets)
  {
    const AliMultiplicity* spdmult = esd->GetMultiplicity();    // spd multiplicity object
    for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i) 
    {
      if (TMath::Abs(spdmult->GetEta(i)) > etaRange) 
	continue; // eta selection for tracklets
      multiplicityEstimate++;
    }
    return multiplicityEstimate;
  }

  if (!vertices[1]->GetStatus())
  {
    AliDebugClass(AliLog::kDebug, "No track vertex. Not able to make a reliable multiplicity estimate.");
    return -3;
  }

  // TODO value of displacement to be studied
  const Float_t maxDisplacement = 0.5;
  //check for displaced vertices
  Double_t displacement = TMath::Abs(vertices[0]->GetZ() - vertices[1]->GetZ());
  if (displacement > maxDisplacement) 
  {
    AliDebugClass(AliLog::kDebug, Form("Displaced vertices %f > %f",displacement,maxDisplacement));
    return -4;
  }
  
  // update eta range in track cuts
  GetMultEstTrackCuts(kMultEstTrackCutITSSA)->SetEtaRange(-etaRange, etaRange);
  GetMultEstTrackCuts(kMultEstTrackCutDCAwSPD)->SetEtaRange(-etaRange, etaRange);
  GetMultEstTrackCuts(kMultEstTrackCutDCAwoSPD)->SetEtaRange(-etaRange, etaRange);
  
  //*******************************************************************************************************
  //set counters to initial zeros
  Int_t tracksITSTPC = 0;     //number of global tracks for a given event
  Int_t tracksITSSA = 0;      //number of ITS standalone tracks for a given event
  Int_t tracksITSTPCSA_complementary = 0; //number of ITS standalone tracks complementary to TPC for a given event
  Int_t trackletsITSTPC_complementary = 0;//number of SPD tracklets complementary to global/ITSSA tracks for a given events
  Int_t trackletsITSSA_complementary = 0; //number of SPD tracklets complementary to ITSSA tracks for a given event

  const Int_t nESDTracks = esd->GetNumberOfTracks();
  Int_t highestID = 0;

  // flags for secondary and rejected tracks
  const Int_t kRejBit = BIT(15); // set this bit in global tracks if it is rejected by a cut
  const Int_t kSecBit = BIT(16); // set this bit in global tracks if it is secondary according to a cut
  
  for(Int_t itracks=0; itracks < nESDTracks; itracks++) {
    if(esd->GetTrack(itracks)->GetLabel() > highestID) highestID = esd->GetTrack(itracks)->GetLabel();
    esd->GetTrack(itracks)->ResetBit(kSecBit|kRejBit); //reset bits used for flagging secondaries and rejected tracks in case they were changed before this analysis
  }
  const Int_t maxid = highestID+1; // used to define bool array for check multiple associations of tracklets to one track. array starts at 0.
  
  // bit mask for esd tracks, to check if multiple tracklets are associated to it
  Bool_t globalBits[maxid], pureITSBits[maxid];
  for(Int_t i=0; i<maxid; i++){ // set all bools to false
      globalBits[i]=kFALSE;
      pureITSBits[i]=kFALSE;
  }

  //*******************************************************************************************************
  // get multiplicity from global tracks
  for(Int_t itracks = 0; itracks < nESDTracks; itracks++) { // flag the tracks
    AliESDtrack* track = esd->GetTrack(itracks);

    // if track is a secondary from a V0, flag as a secondary
    if (track->IsOn(AliESDtrack::kMultInV0)) {
      track->SetBit(kSecBit);
      continue;
    }
    
    //secondary?
    if (track->IsOn(AliESDtrack::kMultSec)) {
      track->SetBit(kSecBit);
      continue;
    }

    // check tracks with ITS part
    //*******************************************************************************************************
    if (track->IsOn(AliESDtrack::kITSin) && !track->IsOn(AliESDtrack::kITSpureSA) && trackType == kTrackletsITSTPC) { // track has ITS part but is not an ITS_SA
      //*******************************************************************************************************
      // TPC+ITS
      if (track->IsOn(AliESDtrack::kTPCin)) {  // Global track, has ITS and TPC contributions
          if (fgMultEstTrackCuts[kMultEstTrackCutGlobal]->AcceptTrack(track)) { // good ITSTPC track
            if (fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD]->AcceptTrack(track) || fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD]->AcceptTrack(track)) {
              tracksITSTPC++; //global track counted
              globalBits[itracks] = kTRUE;
            }
            else track->SetBit(kSecBit); // large DCA -> secondary, don't count either track not associated tracklet
          }
          else track->SetBit(kRejBit); // bad quality, don't count the track, but may count tracklet if associated
      }
      //*******************************************************************************************************
      // ITS complementary
      else if (fgMultEstTrackCuts[kMultEstTrackCutITSSA]->AcceptTrack(track)) { // good ITS complementary track
        if (fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD]->AcceptTrack(track) || fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD]->AcceptTrack(track)) {
          tracksITSTPCSA_complementary++;
          globalBits[itracks] = kTRUE;
        }
        else track->SetBit(kSecBit); // large DCA -> secondary, don't count either track not associated tracklet
      }
      else track->SetBit(kRejBit); // bad quality, don't count the track, but may count tracklet if associated
    }
    //*******************************************************************************************************
    // check tracks from ITS_SA_PURE
    if (track->IsOn(AliESDtrack::kITSin) && track->IsOn(AliESDtrack::kITSpureSA) && trackType == kTrackletsITSSA){
      if (fgMultEstTrackCuts[kMultEstTrackCutITSSA]->AcceptTrack(track)) { // good ITSSA track
        if (fgMultEstTrackCuts[kMultEstTrackCutDCAwSPD]->AcceptTrack(track) || fgMultEstTrackCuts[kMultEstTrackCutDCAwoSPD]->AcceptTrack(track)) {
          tracksITSSA++;
          pureITSBits[itracks] = kTRUE;
        }
        else track->SetBit(kRejBit);
      }
      else track->SetBit(kRejBit);
    }
  }//ESD tracks counted

  //*******************************************************************************************************
  // get multiplicity from ITS tracklets to complement TPC+ITS, and ITSpureSA
  const AliMultiplicity* spdmult = esd->GetMultiplicity();    // spd multiplicity object
  for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i) {
    if (TMath::Abs(spdmult->GetEta(i)) > etaRange) continue; // eta selection for tracklets
    
    // if counting tracks+tracklets, check if clusters were already used in tracks
    Int_t id1,id2,id3,id4;
    spdmult->GetTrackletTrackIDs(i,0,id1,id2); // references for eventual Global/ITS_SA tracks
    AliESDtrack* tr1 = id1>=0 ? esd->GetTrack(id1) : 0;
    spdmult->GetTrackletTrackIDs(i,1,id3,id4); // references for eventual ITS_SA_pure tracks
    AliESDtrack* tr3 = id3>=0 ? esd->GetTrack(id3) : 0;

    // are both clusters from the same tracks? If not, skip the tracklet (shouldn't change things much)
    if ((id1!=id2 && id1>=0 && id2>=0) || (id3!=id4 && id3>=0 && id4>=0)) continue;

    Bool_t bUsedInGlobal = (id1 != -1) ? globalBits[id1] : 0;// has associated global track been associated to a previous tracklet?
    Bool_t bUsedInPureITS = (id3 != -1) ? pureITSBits[id3] : 0;// has associated pure ITS track been associated to a previous tracklet?
    //*******************************************************************************************************
    if (trackType == kTrackletsITSTPC) {
      // count tracklets towards global+complementary tracks
      if (  (tr1 && !tr1->TestBit(kSecBit)) &&   // reject as secondary
        (tr1 &&  tr1->TestBit(kRejBit)) ) {  // count tracklet as bad quality track
          if(!bUsedInGlobal){
            ++trackletsITSTPC_complementary;
            if(id1>0) globalBits[id1] = kTRUE; // mark global track linked to this tracklet as "associated"
          }
      }
      else if(id1<0) {
        ++trackletsITSTPC_complementary; // if no associated track, count the tracklet
      }
    } else {
      // count tracklets towards ITS_SA_pure tracks
      if (  (tr3 && !tr3->TestBit(kSecBit)) &&   // reject as secondary
        (tr3 &&  tr3->TestBit(kRejBit)) ) { // count tracklet as bad quality track
        if(!bUsedInPureITS) {
          ++trackletsITSSA_complementary;
          if(id3>0) pureITSBits[id3] = kTRUE; // mark global track linked to this tracklet as "associated"
        }
      }
      else if(id3<0) {
        ++trackletsITSSA_complementary; // if no associated track, count the tracklet
      }
    }
  }
  
  //*******************************************************************************************************
  if (trackType == kTrackletsITSTPC)
    multiplicityEstimate = tracksITSTPC + tracksITSTPCSA_complementary + trackletsITSTPC_complementary;
  else
    multiplicityEstimate = tracksITSSA + trackletsITSSA_complementary;

  return multiplicityEstimate;
}
