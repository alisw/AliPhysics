////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorV0                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorV0.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include "AliFemtoModelHiddenInfo.h"

AliFemtoCutMonitorV0::AliFemtoCutMonitorV0():
  fLambdaMass(0),
  fAntiLambdaMass(0),
  fK0ShortMass(0),
  fDcaDaughters(0), 
  fDcaV0ToPrimVertex(0),
  fCosPointingAngle(0),
  fEtaV0(0),
  fPtV0(0),
  fPtPosDaughter(0),
  fPtNegDaughter(0)
{
  // Default constructor
  fLambdaMass = new TH1F("LambdaMass", "Mass Assuming Lambda Hypothesis", 10000, 0, 5);
  fAntiLambdaMass = new TH1F("AntiLambdaMass", "Mass Assuming AntiLambda Hypothesis", 10000, 0, 5);
  fK0ShortMass= new TH1F("K0ShortMass", "Mass Assuming K0 short Hypothesis", 500, 0, 5);
  fDcaDaughters = new TH1F("DcaDaughters", "DCA Daughters", 500, 0, 2);
  fDcaV0ToPrimVertex = new TH1F("DcaV0ToPrimVertex", "DCA V0 to primary vertex", 500, 0, 3);
  fCosPointingAngle = new TH1F("CosPointingAngle","Cosinus Pointing Angle",500,0,1);
  fEtaV0 = new TH1F("EtaV0", "|Eta| distribution of V0s", 500, 0.0, 8.);
  fPtV0 = new TH1F("PtV0", "Pt distribution of V0s", 500, 0.0, 8.);
  fPtPosDaughter = new TH1F("PtPosDaughter", "Pt distribution of positive daughters", 500, 0.0, 5.);
  fPtNegDaughter = new TH1F("PtNegDaughter", "Pt distribution of negative daughters", 500, 0.0, 5.);


  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();

}

AliFemtoCutMonitorV0::AliFemtoCutMonitorV0(const char *aName):
  AliFemtoCutMonitor(),
  fLambdaMass(0),
  fAntiLambdaMass(0),
  fK0ShortMass(0),
  fDcaDaughters(0), 
  fDcaV0ToPrimVertex(0),
  fCosPointingAngle(0),
  fEtaV0(0),
  fPtV0(0),
  fPtPosDaughter(0),
  fPtNegDaughter(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "LambdaMass%s", aName);
  fLambdaMass = new TH1F(name, "Mass Assuming Lambda Hypothesis", 10000, 0, 5);
  snprintf(name, 200, "AntiLambdaMass%s", aName);
  fAntiLambdaMass = new TH1F(name, "Mass Assuming AntiLambda Hypothesis", 10000, 0, 5);
  snprintf(name, 200, "K0ShortMass%s", aName);
  fK0ShortMass = new TH1F(name, "Mass Assuming K0 short Hypothesis", 500, 0, 5);
  snprintf(name, 200, "DcaDaughters%s", aName);
  fDcaDaughters = new TH1F(name, "DCA Daughters", 500, 0, 2);
  snprintf(name, 200, "DcaV0ToPrimVertex%s", aName);
  fDcaV0ToPrimVertex = new TH1F(name, "DCA V0 to primary vertex", 500, 0, 3);
  snprintf(name, 200, "CosPointingAngle%s", aName);
  fCosPointingAngle = new TH1F(name,"Cosinus Pointing Angle",500,0,1);
  snprintf(name, 200, "EtaV0%s", aName);
  fEtaV0 = new TH1F(name, "|Eta| distribution of V0s", 500, 0.0, 1.); 
  snprintf(name, 200, "PtV0%s", aName);
  fPtV0 = new TH1F(name, "Pt distribution of V0s", 500, 0.0, 8.); 
  snprintf(name, 200, "fPtPosDaughter%s", aName);
  fPtPosDaughter = new TH1F(name, "Pt distribution of positive daughters", 500, 0.0, 5.);
  snprintf(name, 200, "fPtNegDaughter%s", aName);
  fPtNegDaughter = new TH1F(name, "Pt distribution of negative daughters", 500, 0.0, 5.);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();
}

AliFemtoCutMonitorV0::AliFemtoCutMonitorV0(const AliFemtoCutMonitorV0 &aCut):
  AliFemtoCutMonitor(),
  fLambdaMass(0),
  fAntiLambdaMass(0),
  fK0ShortMass(0),
  fDcaDaughters(0), 
  fDcaV0ToPrimVertex(0),
  fCosPointingAngle(0),
  fEtaV0(0),
  fPtV0(0),
  fPtPosDaughter(0),
  fPtNegDaughter(0)
{
  // copy constructor
  if (fLambdaMass) delete fLambdaMass;
  fLambdaMass = new TH1F(*aCut.fLambdaMass);
  if (fAntiLambdaMass) delete fAntiLambdaMass;
  fAntiLambdaMass = new TH1F(*aCut.fAntiLambdaMass);
  if (fK0ShortMass) delete fK0ShortMass;
  fK0ShortMass = new TH1F(*aCut.fK0ShortMass);
  if (fDcaDaughters) delete fDcaDaughters;
  fDcaDaughters = new TH1F(*aCut.fDcaDaughters);
  if (fDcaV0ToPrimVertex) delete fDcaV0ToPrimVertex;
  fDcaV0ToPrimVertex = new TH1F(*aCut.fDcaV0ToPrimVertex);
  if(fCosPointingAngle) delete fCosPointingAngle;
  fCosPointingAngle = new TH1F(*aCut.fCosPointingAngle);
  if(fEtaV0) delete fEtaV0;
  fEtaV0 = new TH1F(*aCut.fEtaV0);
  if(fPtV0) delete fPtV0;
  fPtV0 = new TH1F(*aCut.fPtV0);
  if(fPtPosDaughter) delete fPtPosDaughter;
  fPtPosDaughter = new TH1F(*aCut.fPtPosDaughter);
  if(fPtNegDaughter) delete fPtNegDaughter;
  fPtNegDaughter = new TH1F(*aCut.fPtNegDaughter);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();
}

AliFemtoCutMonitorV0::~AliFemtoCutMonitorV0()
{
  // Destructor
  delete fLambdaMass;
  delete fAntiLambdaMass;
  delete fK0ShortMass;
  delete fDcaDaughters; 
  delete fDcaV0ToPrimVertex;
  delete fCosPointingAngle;
  delete fEtaV0;
  delete fPtV0;
  delete fPtPosDaughter;
  delete fPtNegDaughter;
}

AliFemtoCutMonitorV0& AliFemtoCutMonitorV0::operator=(const AliFemtoCutMonitorV0& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fLambdaMass) delete fLambdaMass;
  fLambdaMass = new TH1F(*aCut.fLambdaMass);
  if (fAntiLambdaMass) delete fAntiLambdaMass;
  fAntiLambdaMass = new TH1F(*aCut.fAntiLambdaMass);
  if (fK0ShortMass) delete fK0ShortMass;
  fK0ShortMass = new TH1F(*aCut.fK0ShortMass);
  if (fDcaDaughters) delete fDcaDaughters;
  fDcaDaughters = new TH1F(*aCut.fDcaDaughters);
  if (fDcaV0ToPrimVertex) delete fDcaV0ToPrimVertex;
  fDcaV0ToPrimVertex = new TH1F(*aCut.fDcaV0ToPrimVertex);
  if(fCosPointingAngle) delete fCosPointingAngle;
  fCosPointingAngle = new TH1F(*aCut.fCosPointingAngle);
  if(fEtaV0) delete fEtaV0;
  fEtaV0 = new TH1F(*aCut.fEtaV0);
  if(fPtV0) delete fPtV0;
  fPtV0 = new TH1F(*aCut.fPtV0);
  if(fPtPosDaughter) delete fPtPosDaughter;
  fPtPosDaughter = new TH1F(*aCut.fPtPosDaughter);
  if(fPtNegDaughter) delete fPtNegDaughter;
  fPtNegDaughter = new TH1F(*aCut.fPtNegDaughter);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();

  return *this;
}

AliFemtoString AliFemtoCutMonitorV0::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorV0 report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorV0::Fill(const AliFemtoV0* aV0)
{
  // Fill momentum resolution histograms for the particle
  fLambdaMass->Fill(aV0->MassLambda());
  fAntiLambdaMass->Fill(aV0->MassAntiLambda());
  fK0ShortMass->Fill(aV0->MassK0Short());
  fDcaDaughters->Fill(aV0->DcaV0Daughters()); 
  fDcaV0ToPrimVertex->Fill(aV0->DcaV0ToPrimVertex());
  fCosPointingAngle->Fill(aV0->CosPointingAngle());
  fEtaV0->Fill(aV0->EtaV0());
  fPtV0->Fill(aV0->PtV0());
  fPtPosDaughter->Fill(aV0->PtPos());
  fPtNegDaughter->Fill(aV0->PtNeg());
}

void AliFemtoCutMonitorV0::Write()
{
  // Write out the relevant histograms
  fLambdaMass->Write();
  fAntiLambdaMass->Write();
  fK0ShortMass->Write();
  fDcaDaughters->Write();
  fDcaV0ToPrimVertex->Write();
  fCosPointingAngle->Write();
  fEtaV0->Write();
  fPtV0->Write();
  fPtPosDaughter->Write();
  fPtNegDaughter->Write();
}

TList *AliFemtoCutMonitorV0::GetOutputList()
{
  // Get the list of histograms to write
  TList *tOutputList = new TList();
  tOutputList->Add(fLambdaMass);
  tOutputList->Add(fAntiLambdaMass);
  tOutputList->Add(fK0ShortMass);
  tOutputList->Add(fDcaDaughters);
  tOutputList->Add(fDcaV0ToPrimVertex);
  tOutputList->Add(fCosPointingAngle);
  tOutputList->Add(fEtaV0);
  tOutputList->Add(fPtV0);
  tOutputList->Add(fPtPosDaughter);
  tOutputList->Add(fPtNegDaughter);

  return tOutputList;
}
