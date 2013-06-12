//
// Cluster Container
//
// Author: M. Verweij

#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include <TChain.h>
#include <TClonesArray.h>
#include <TObject.h>
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
  fMinMCLabel(0),
  fVVertex(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

}

//________________________________________________________________________
AliClusterContainer::AliClusterContainer(const char *name):
  AliEmcalContainer(name),
  fClusPtCut(0.15),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fClusterBitMap(0),
  fMCClusterBitMap(0),
  fMinMCLabel(0),
  fVVertex(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

}

//________________________________________________________________________
AliClusterContainer::~AliClusterContainer()
{
  // Destructor.
}

//________________________________________________________________________
void AliClusterContainer::SetClusterArray(AliVEvent *event) 
{
  // Set jet array

  SetArray(event, "AliVCluster");

  fVVertex = event->GetPrimaryVertex();
  if (fVVertex) {
    fVVertex->GetXYZ(fVertex);
  }

}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetCluster(Int_t i) const {

  //Get i^th jet in array

  if(i<0 || i>fClArray->GetEntriesFast()) return 0;
  AliVCluster *vp = static_cast<AliVCluster*>(fClArray->At(i));
  return vp;

}

//________________________________________________________________________
AliVCluster* AliClusterContainer::GetAcceptCluster(Int_t i) const {
  //return pointer to particle if particle is accepted

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
