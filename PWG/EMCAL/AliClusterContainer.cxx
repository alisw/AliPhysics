//
// Cluster Container
//
// Author: M. Verweij

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliLog.h"

#include "AliClusterContainer.h"

ClassImp(AliClusterContainer)

//________________________________________________________________________
AliClusterContainer::AliClusterContainer():
  AliEmcalContainer("AliClusterContainer"),
  fClusPtCut(0.15),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fClusterBitMap(0),
  fMCClusterBitMap(0),
  fMinMCLabel(0)
{
  // Default constructor.

}

//________________________________________________________________________
AliClusterContainer::AliClusterContainer(const char *name):
  AliEmcalContainer(name),
  fClusPtCut(0.15),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fClusterBitMap(0),
  fMCClusterBitMap(0),
  fMinMCLabel(0)
{
  // Standard constructor.

}

//________________________________________________________________________
void AliClusterContainer::SetClusterArray(AliVEvent *event) 
{
  // Set jet array

  SetArray(event, "AliVCluster");
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetLeadingCluster(const char* opt) const
{
  // Get the leading cluster; use e if "e" is contained in opt (otherwise et)

  TString option(opt);
  option.ToLower();

  AliVCluster *clusterMax = GetNextAcceptCluster(0);
  AliVCluster *cluster = 0;

  if (option.Contains("e")) {
    while ((cluster = GetNextAcceptCluster())) {
      if (cluster->E() > clusterMax->E()) clusterMax = cluster;
    }
  }
  else {
    Double_t et = 0;
    Double_t etmax = 0;
    while ((cluster = GetNextAcceptCluster())) {
      TLorentzVector mom;
      cluster->GetMomentum(mom,const_cast<Double_t*>(fVertex));
      et = mom.Et();
      if (et > etmax) { 
	clusterMax = cluster;
	etmax = et;
      }
    }
  }

  return clusterMax;
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetCluster(Int_t i) const {

  //Get i^th cluster in array

  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliVCluster *vp = static_cast<AliVCluster*>(fClArray->At(i));
  return vp;

}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetAcceptCluster(Int_t i) const {
  //return pointer to cluster if cluster is accepted

  AliVCluster *vc = GetCluster(i);
  if(!vc) return 0;

  if(AcceptCluster(vc))
    return vc;
  else {
    AliDebug(2,"Cluster not accepted.");
    return 0;
  }
}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetNextAcceptCluster(Int_t i) const {

  //Get next accepted cluster; if i >= 0 (re)start counter from i; return 0 if no accepted particle could be found

  static Int_t counter = -1;
  if (i>=0) counter = i;

  const Int_t n = GetNEntries();
  AliVCluster *c = 0;
  while (counter < n && !c) { 
    c = GetAcceptCluster(counter);
    counter++;
  }

  return c;
}

//________________________________________________________________________
void AliClusterContainer::GetMomentum(TLorentzVector &mom, Int_t i) const
{
  //Get momentum of the i^th cluster in array

  AliVCluster *vc = GetCluster(i);
  if(vc) vc->GetMomentum(mom,const_cast<Double_t*>(fVertex));
}

//________________________________________________________________________
Bool_t AliClusterContainer::AcceptCluster(AliVCluster *clus) const
{
  // Return true if cluster is accepted.

  if (!clus)
    return kFALSE;
      
  if (!clus->IsEMCAL())
    return kFALSE;

  if (clus->GetLabel() > fMinMCLabel) {
    if (clus->TestBits(fMCClusterBitMap) != (Int_t)fMCClusterBitMap) {
      AliDebug(2,"MC Cluster not accepted because of MC bit map.");
      return kFALSE;
    }
  }
  else {
    if (clus->TestBits(fClusterBitMap) != (Int_t)fClusterBitMap) {
      AliDebug(2,"Cluster not accepted because of bit map.");
      return kFALSE;
    }
  }

  if (clus->GetTOF() > fClusTimeCutUp || clus->GetTOF() < fClusTimeCutLow)
    return kFALSE;

  TLorentzVector nPart;
  clus->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

  if (nPart.Et() < fClusPtCut)
    return kFALSE;
  
  return kTRUE;

}
