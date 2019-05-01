// #include "iterator"

// --- Custom header files ---
#include "AliPP13EfficiencySelectionSPMC.h"

// --- ROOT system ---
#include <TParticle.h>
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliLog.h>
#include <AliVCluster.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13EfficiencySelectionSPMC);


//________________________________________________________________
void AliPP13EfficiencySelectionSPMC::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	TLorentzVector p1 = ClusterMomentum(c1, eflags);
	TLorentzVector p2 = ClusterMomentum(c2, eflags);
	TLorentzVector psum = p1 + p2;

	// Pair cuts can be applied here
	if (psum.M2() < 0)  return;

	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok


	Double_t ma12 = psum.M();
	Double_t pt12 = psum.Pt();

	Double_t w = fWeights->Weights(pt12, eflags);
	TH2 * hist = dynamic_cast<TH2 *> (fInvMass[eflags.isMixing]);
	hist->Fill(ma12, pt12, w);
}
