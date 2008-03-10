////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticleMomRes - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticleMomRes.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include "AliFemtoModelHiddenInfo.h"

AliFemtoCutMonitorParticleMomRes::AliFemtoCutMonitorParticleMomRes():
  fMomRes3D(0),
  fMomResXvsP(0),
  fMomResYvsP(0),
  fMomResZvsP(0),
  fImpactXY(0),
  fImpactZ(0),
  fSigma(0),
  fMass(0.13957)
{
  // Default constructor
  fMomRes3D = new TH3D("MomRes3D", "Momentum resolution", 100, -0.05, 0.05, 100, -0.05, 0.05, 100, -0.05, 0.05);
  fMomResXvsP = new TH2D("MomResXvsP", "X momentum resolution vs P", 100, 0.1, 2.0, 100, -0.05, 0.05);
  fMomResYvsP = new TH2D("MomResYvsP", "Y momentum resolution vs P", 100, 0.1, 2.0, 100, -0.05, 0.05);
  fMomResZvsP = new TH2D("MomResZvsP", "Z momentum resolution vs P", 100, 0.1, 2.0, 100, -0.05, 0.05);
  fImpactXY   = new TH2D("ImpactXY", "XY impact parameter vs P", 100, 0.1, 2.0, 200, -1.0, 1.0);
  fImpactZ    = new TH2D("ImpactZ",  "Z impact parameter vs P" , 100, 0.1, 2.0, 200, -1.0, 1.0);
  fSigma      = new TH2D("Sigma",     "Sigma to vertex vs P" , 100, 0.1, 2.0, 200, -5.0, 5.0);
}

AliFemtoCutMonitorParticleMomRes::AliFemtoCutMonitorParticleMomRes(const char *aName, float aMass):
  fMomRes3D(0),
  fMomResXvsP(0),
  fMomResYvsP(0),
  fMomResZvsP(0),
  fImpactXY(0),
  fImpactZ(0),
  fSigma(0),
  fMass(aMass)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "MomRes3D%s", aName);
  fMomRes3D = new TH3D(name, "Momentum resolution", 100, -0.05, 0.05, 100, -0.05, 0.05, 100, -0.05, 0.05);
  snprintf(name, 200, "MomResXvsP%s", aName);
  fMomResXvsP = new TH2D(name, "X momentum resolution vs P", 100, 0.1, 2.0, 100, -0.05, 0.05);
  snprintf(name, 200, "MomResYvsP%s", aName);
  fMomResYvsP = new TH2D(name, "Y momentum resolution vs P", 100, 0.1, 2.0, 100, -0.05, 0.05);
  snprintf(name, 200, "MomResZvsP%s", aName);
  fMomResZvsP = new TH2D(name, "Z momentum resolution vs P", 100, 0.1, 2.0, 100, -0.05, 0.05);
  snprintf(name, 200, "ImpactXY%s", aName);
  fImpactXY   = new TH2D(name, "XY impact parameter vs P", 100, 0.1, 2.0, 200, -1.0, 1.0);
  snprintf(name, 200, "ImpactZ%s", aName);
  fImpactZ    = new TH2D(name,  "Z impact parameter vs P" , 100, 0.1, 2.0, 200, -1.0, 1.0);
  snprintf(name, 200, "Sigma%s", aName);
  fSigma    = new TH2D(name,  "Z impact parameter vs P" , 100, 0.1, 2.0, 200, -5.0, 5.0);
}

AliFemtoCutMonitorParticleMomRes::AliFemtoCutMonitorParticleMomRes(const AliFemtoCutMonitorParticleMomRes &aCut):
  AliFemtoCutMonitor(),
  fMomRes3D(0),
  fMomResXvsP(0),
  fMomResYvsP(0),
  fMomResZvsP(0),
  fImpactXY(0),
  fImpactZ(0),
  fSigma(0),
  fMass(0.13957)
{
  // copy constructor
  if (fMomRes3D) delete fMomRes3D;
  fMomRes3D = new TH3D(*aCut.fMomRes3D);
  if (fMomResXvsP) delete fMomResXvsP;
  fMomResXvsP = new TH2D(*aCut.fMomResXvsP);
  if (fMomResYvsP) delete fMomResYvsP;
  fMomResYvsP = new TH2D(*aCut.fMomResYvsP);
  if (fMomResZvsP) delete fMomResZvsP;
  fMomResZvsP = new TH2D(*aCut.fMomResZvsP);
  if (fImpactXY) delete fImpactXY;
  fImpactXY = new TH2D(*aCut.fImpactXY);
  if (fImpactZ) delete fImpactZ;
  fImpactZ = new TH2D(*aCut.fImpactZ);
  if (fSigma) delete fSigma;
  fSigma = new TH2D(*aCut.fSigma);
  fMass = aCut.fMass; 
}

AliFemtoCutMonitorParticleMomRes::~AliFemtoCutMonitorParticleMomRes()
{
  // Destructor
  delete fMomRes3D;
  delete fMomResXvsP;
  delete fMomResYvsP;
  delete fMomResZvsP;
  delete fImpactXY;
  delete fImpactZ;
  delete fSigma;
}

AliFemtoCutMonitorParticleMomRes& AliFemtoCutMonitorParticleMomRes::operator=(const AliFemtoCutMonitorParticleMomRes& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fMomRes3D) delete fMomRes3D;
  fMomRes3D = new TH3D(*aCut.fMomRes3D);
  if (fMomResXvsP) delete fMomResXvsP;
  fMomResXvsP = new TH2D(*aCut.fMomResXvsP);
  if (fMomResYvsP) delete fMomResYvsP;
  fMomResYvsP = new TH2D(*aCut.fMomResYvsP);
  if (fMomResZvsP) delete fMomResZvsP;
  fMomResZvsP = new TH2D(*aCut.fMomResZvsP);
  if (fImpactXY) delete fImpactXY;
  fImpactXY = new TH2D(*aCut.fImpactXY);
  if (fImpactZ) delete fImpactZ;
  fImpactZ = new TH2D(*aCut.fImpactZ);
  if (fSigma) delete fSigma;
  fSigma = new TH2D(*aCut.fSigma);
  
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticleMomRes::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorParticleMomRes report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorParticleMomRes::Fill(const AliFemtoTrack* aTrack)
{
  // Fill momentum resolution histograms for the particle
  AliFemtoModelHiddenInfo *tInf = ( AliFemtoModelHiddenInfo *) aTrack->GetHiddenInfo();
  fMomRes3D->Fill(tInf->GetTrueMomentum()->x() - aTrack->P().x(),
		  tInf->GetTrueMomentum()->y() - aTrack->P().y(),
		  tInf->GetTrueMomentum()->z() - aTrack->P().z());

  fMomResXvsP->Fill(aTrack->P().mag(),
		    tInf->GetTrueMomentum()->x() - aTrack->P().x());
  fMomResYvsP->Fill(aTrack->P().mag(),
		    tInf->GetTrueMomentum()->y() - aTrack->P().y());
  fMomResZvsP->Fill(aTrack->P().mag(),
		    tInf->GetTrueMomentum()->z() - aTrack->P().z());
  fImpactXY->Fill(aTrack->P().mag(),
		  aTrack->ImpactD());
  fImpactZ->Fill(aTrack->P().mag(),
		 aTrack->ImpactZ());
  fSigma->Fill(aTrack->P().mag(),
	       aTrack->SigmaToVertex());
}

void AliFemtoCutMonitorParticleMomRes::Write()
{
  // Write out the relevant histograms
  fMomRes3D->Write();
  fMomResXvsP->Write();
  fMomResYvsP->Write();
  fMomResZvsP->Write();
  fImpactXY->Write();
  fImpactZ->Write();
  fSigma->Write();
}

TList *AliFemtoCutMonitorParticleMomRes::GetOutputList()
{
  // Get the list of histograms to write
  TList *tOutputList = new TList();
  tOutputList->Add(fMomRes3D);
  tOutputList->Add(fMomResXvsP);
  tOutputList->Add(fMomResYvsP);
  tOutputList->Add(fMomResZvsP);
  tOutputList->Add(fImpactXY);
  tOutputList->Add(fImpactZ);
  tOutputList->Add(fSigma);

  return tOutputList;
}
