// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/**
 * Display HLT TPC clusters in AliEVE.
 * The cluster data must be stored in an HLTOUT collection originating
 * either from an HLT simulation or real run with a chain adding
 * HLT TPC cluster structures of type AliHLTTPCDefinitions::fgkClustersDataType
 * to the HLTOUT.
 *
 * The HLTOUT data is read either from the RawReader which AliEve has been
 * initialized or the HLT.Digits.root if the RawReader is not available.
 * Please note: As of Nov 2008, AliEve can handle raw data noy through
 * AliRawReaderFile. It can only open a RawReader on a single file, by default
 * raw.root. If the file is not existing, the RawReader is not initialized.
 *
 * Usage:
 * <pre>
 *   alieve $ALICE_ROOT/EVE/macros/event_next.C \
 *          $ALICE_ROOT/EVE/macros/alieve_init.C \
 *          $ALICE_ROOT/EVE/macros/geom_simple.C \
 *          $ALICE_ROOT/HLT/TPCLib/EVE/macros/hlt_tpc_clusters.C
 * </pre>
 * Display is changed to next event by executing event_next()
 * from the root prompt.
 * <pre>
 *   event_next(); hlt_tpc_clusters();
 * </pre>
 *
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 * @date   2008-11-22
 */
TEvePointSet* hlt_tpc_clusters(const char* digitfile=NULL, TEveElement* cont=0, Float_t maxR=270)
{
  if (!TClass::GetClass("AliEveEventManager")) {
    Error("hlt_tpc_clusters.C", "EVE library not loaded, please start alieve correctly");
    return NULL;
  }

  AliEveEventManager* eveManager=AliEveEventManager::GetMaster();
  if (!eveManager) {
    Error("hlt_tpc_clusters.C", "EVE manager not initialized");
    return NULL;
  }

  eveManager->AssertGeometry();

  TClass* pCl=NULL;
  int iLibResult=0;
  gSystem->Load("libAliHLTUtil");
  gSystem->Load("libAliHLTRCU");
  do {
    pCl=TClass::GetClass("AliHLTTPDefinitions");
  } while (!pCl && (iLibResult=gSystem->Load("libAliHLTTPC"))==0);
  do {
    pCl=TClass::GetClass("AliHLTTPCEVE");
  } while (!pCl && (iLibResult=gSystem->Load("libAliHLTTPCEVE"))==0);

  AliHLTTPCEVE hlttpceve;
  TEvePointSet* clusters = NULL;

  AliESDEvent* pESD=eveManager->AssertESD();
  // extract from RawReader if available and no digit file has been specified
  if (digitfile==NULL && eveManager->HasRawReader()) {
    AliRawReader* pRawReader=eveManager->AssertRawReader();
    if (pRawReader) {
      Info("hlt_tpc_clusters.C", "extracting HLT TPC clusters from RawReader");
      clusters=hlttpceve.MakePointSetFromHLTOUT(pRawReader, cont, maxR);
    }
  } else {
    Info("hlt_tpc_clusters.C", Form("extracting HLT TPC clusters from digit file %s", digitfile!=NULL?digitfile:""));
    clusters=hlttpceve.MakePointSetFromHLTDigits(digitfile, eveManager->GetEventId(), cont, maxR);
  }
  if (!clusters) return NULL;
  
  if (clusters->Size() == 0) {
    Info("hlt_tpc_clusters.C", "No TPC clusters");
  }

  //clusters->ApplyVizTag(clusters->GetName());

  gEve->AddElement(clusters, cont);

  gEve->Redraw3D();

  return clusters;
}
