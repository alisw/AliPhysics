#ifndef __CINT__

#include "TString.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskEmcalJetCDF.h"

#endif

/// \file AddTaskEmcalJetCDF.C
/// \brief Adds a AliAnalysisTaskEmcalJetCDF analysis task and coresponding containers
///
/// Analysis of jets (all and leading) distribution of pt and multiplicity, R distribution, N{70,75,80,85,90} and Pt{70,75,80,85,90}
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 19, 2015

/// Add a AliAnalysisTaskEmcalJetCDF task - detailed signature
/// \param const char* ntracks : name of tracks collection
/// \param const char* nclusters : name of clusters collection
/// \param const char* ncells : name of EMCAL cell collection
/// \param const char* tag
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF* AddTaskEmcalJetCDF (
  const char* ntracks                      = "usedefault",
  const char* nclusters                    = "usedefault",
  const char* ncells                       = "usedefault",
  const char* tag                          = "CDF"
)
  {
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if ( !mgr ) { ::Error ( "AddTaskEmcalJetCDF", "No analysis manager to connect to." );  return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) { ::Error ( "AddTaskEmcalJetCDF", "This task requires an input event handler" ); return NULL; }

  enum EDataType_t { kUnknown, kESD, kAOD }; EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) { dataType = kESD; }
  else
  if (handler->InheritsFrom("AliAODInputHandler")) { dataType = kAOD; }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString suffix   ( tag );
  TString tracks   ( ntracks );
  TString clusters ( nclusters );
  TString cells    ( ncells );

  if ( tracks.EqualTo("usedefault") )
    {
    if ( dataType == kESD ) { tracks = "Tracks"; }
    else
    if ( dataType == kAOD ) { tracks = "tracks"; }
    else
      { tracks = ""; }
    }

  if ( clusters.EqualTo("usedefault") )
    {
    if ( dataType == kESD ) { clusters = "CaloClusters"; }
    else
    if ( dataType == kAOD ) { clusters = "caloClusters"; }
    else
      { clusters = ""; }
    }

  if ( cells.EqualTo("usedefault") )
    {
    if (dataType == kESD) { cells = "EMCALCells"; }
    else
    if (dataType == kAOD) { cells = "emcalCells"; }
    else
      { cells = ""; }
    }

  TString name("JetCDF");
  if (!tracks.IsNull())   { name += "_" + tracks; }
  if (!clusters.IsNull()) { name += "_" + clusters; }
  if (!cells.IsNull())    { name += "_" + cells; }
  if (!suffix.IsNull())   { name += "_" + suffix; }

  AliAnalysisTaskEmcalJetCDF* cdfTask = new AliAnalysisTaskEmcalJetCDF ( name.Data() );
  cdfTask->SetVzRange(-10,10);
  cdfTask->SetCaloCellsName(cells.Data());

  if ( tracks.EqualTo("mcparticles") )
    { AliMCParticleContainer* mcpartCont = cdfTask->AddMCParticleContainer ( tracks.Data() ); }
  else
  if ( tracks.EqualTo("tracks") || tracks.EqualTo("Tracks") )
    { AliTrackContainer* trackCont = cdfTask->AddTrackContainer( tracks.Data() ); }
  else
  if ( !tracks.IsNull())
    { cdfTask->AddParticleContainer(tracks.Data()); }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask ( cdfTask );

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();

  TString contname = name + "_histos";
  TString outfile (Form("%s", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer ( contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data() );

  mgr->ConnectInput  ( cdfTask, 0,  cinput1 );
  mgr->ConnectOutput ( cdfTask, 1, coutput1 );

  return cdfTask;
  }


// kate: indent-mode none; indent-width 2; replace-tabs on;
