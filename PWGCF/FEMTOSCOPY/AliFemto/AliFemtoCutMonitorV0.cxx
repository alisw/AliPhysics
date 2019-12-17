/// \class AliFemtoCutMonitorV0
/// \brief AliFemtoCutMonitorV0

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
  fDcaPosToPrimVertex(0),
  fDcaNegToPrimVertex(0),
  fCosPointingAngle(0),
  fDecayLength(0),
  fEtaV0(0),
  fPtV0(0),
  fPtPosDaughter(0),
  fPtNegDaughter(0),
  fdEdxPosDaughter(0),
  fdEdxNegDaughter(0),
  fTOFtimePosDaughter(0),
  fTOFtimeNegDaughter(0),
  fMINVvsPt(0),
  fnsigmaPosL(0),
  fnsigmaNegL(0),
  fnsigmaPosAL(0),
  fnsigmaNegAL(0),
  fParticleOrigin(0),
  fParticleId(0)
{
  /// Default constructor

  fLambdaMass = new TH1F("LambdaMass", "Mass Assuming Lambda Hypothesis", 10000, 0, 5);
  fAntiLambdaMass = new TH1F("AntiLambdaMass", "Mass Assuming AntiLambda Hypothesis", 10000, 0, 5);
  fK0ShortMass= new TH1F("K0ShortMass", "Mass Assuming K0 short Hypothesis", 500, 0, 5);
  fDcaDaughters = new TH1F("DcaDaughters", "DCA Daughters", 500, 0, 2);
  fDcaV0ToPrimVertex = new TH1F("DcaV0ToPrimVertex", "DCA V0 to primary vertex", 500, 0, 3);
  fDcaPosToPrimVertex = new TH1F("DcaPosToPrimVertex", "DCA V0 to primary vertex", 500, 0, 3);
  fDcaNegToPrimVertex = new TH1F("DcaNegToPrimVertex", "DCA V0 to primary vertex", 500, 0, 3);
  fCosPointingAngle = new TH1F("CosPointingAngle","Cosinus Pointing Angle",500,0,1);
  fDecayLength = new TH1F("DecayLength","V0 Decay Length",100,0,100);
  fEtaV0 = new TH1F("EtaV0", "|Eta| distribution of V0s", 500, 0.0, 8.);
  fPtV0 = new TH1F("PtV0", "Pt distribution of V0s", 500, 0.0, 8.);
  fPtPosDaughter = new TH1F("PtPosDaughter", "Pt distribution of positive daughters", 500, 0.0, 5.);
  fPtNegDaughter = new TH1F("PtNegDaughter", "Pt distribution of negative daughters", 500, 0.0, 5.);

  fdEdxPosDaughter = new TH2D("dEdxPosDaughter","dEdx of positive daughters",200, 0.1, 4.0, 250, 0.0, 500.0);
  fdEdxNegDaughter = new TH2D("dEdxNegDaughter","dEdx of negative daughters",200, 0.1, 4.0, 250, 0.0, 500.0);
  fTOFtimePosDaughter = new TH2D("TOFtimePosDaughter","TOF time of positive daughters",100,0.,1.1,100,0.,3.0);
  fTOFtimeNegDaughter = new TH2D("TOFtimeNegDaughter","TOF time of negative daughters",100,0.,1.1,100,0.,3.0);
  fMINVvsPt = new TH2D("fMINVvsPt", "Invarint mass V0 vs Pt",50, 0., 5., 1000, 0., 2.);

  fnsigmaPosL = new TH1D("fnsigmaPosL","Number of sigmas of positive Lambda daughters",200,-8,8);
  fnsigmaNegL = new TH1D("fnsigmaNegL","Number of sigmas of negative Lambda daughters",200,-8,8);
  fnsigmaPosAL = new TH1D("fnsigmaPosAL","Number of sigmas of positive AntiLambda daughters",200,-8,8);
  fnsigmaNegAL = new TH1D("fnsigmaNegAL","Number of sigmas of negative AntiLambda daughters",200,-8,8);

  fParticleOrigin =  new TH1D("POrigin", "Mothers PDG Codes", 6000, 0.0, 6000.0);
  fParticleId =  new TH1D("PId", "Particle PDG Codes", 6000, 0.0, 6000.0);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fDcaPosToPrimVertex->Sumw2();
  fDcaNegToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fDecayLength->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();
  fdEdxPosDaughter->Sumw2();
  fdEdxNegDaughter->Sumw2();
  fTOFtimePosDaughter->Sumw2();
  fTOFtimeNegDaughter->Sumw2();
  fMINVvsPt->Sumw2();
  fnsigmaPosL->Sumw2();
  fnsigmaNegL->Sumw2();
  fnsigmaPosAL->Sumw2();
  fnsigmaNegAL->Sumw2();
}

AliFemtoCutMonitorV0::AliFemtoCutMonitorV0(const char *aName):
  AliFemtoCutMonitor(),
  fLambdaMass(0),
  fAntiLambdaMass(0),
  fK0ShortMass(0),
  fDcaDaughters(0),
  fDcaV0ToPrimVertex(0),
  fDcaPosToPrimVertex(0),
  fDcaNegToPrimVertex(0),
  fCosPointingAngle(0),
  fDecayLength(0),
  fEtaV0(0),
  fPtV0(0),
  fPtPosDaughter(0),
  fPtNegDaughter(0),
  fdEdxPosDaughter(0),
  fdEdxNegDaughter(0),
  fTOFtimePosDaughter(0),
  fTOFtimeNegDaughter(0),
  fMINVvsPt(0),
  fnsigmaPosL(0),
  fnsigmaNegL(0),
  fnsigmaPosAL(0),
  fnsigmaNegAL(0),
  fParticleOrigin(0),
  fParticleId(0)
{
  /// Normal constructor

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
  snprintf(name, 200, "DcaPosToPrimVertex%s", aName);
  fDcaPosToPrimVertex = new TH1F(name, "DCA pos. daughter V0 to primary vertex", 500, 0, 3);
  snprintf(name, 200, "DcaNegToPrimVertex%s", aName);
  fDcaNegToPrimVertex = new TH1F(name, "DCA neg. daughter V0 to primary vertex", 500, 0, 3);
  snprintf(name, 200, "CosPointingAngle%s", aName);
  fCosPointingAngle = new TH1F(name,"Cosinus Pointing Angle",500,0,1);
  snprintf(name, 200, "DecayLength%s", aName);
  fDecayLength = new TH1F(name,"Decay Length",100,0,100);
  snprintf(name, 200, "EtaV0%s", aName);
  fEtaV0 = new TH1F(name, "|Eta| distribution of V0s", 500, 0.0, 1.);
  snprintf(name, 200, "PtV0%s", aName);
  fPtV0 = new TH1F(name, "Pt distribution of V0s", 500, 0.0, 8.);
  snprintf(name, 200, "fPtPosDaughter%s", aName);
  fPtPosDaughter = new TH1F(name, "Pt distribution of positive daughters", 500, 0.0, 5.);
  snprintf(name, 200, "fPtNegDaughter%s", aName);
  fPtNegDaughter = new TH1F(name, "Pt distribution of negative daughters", 500, 0.0, 5.);
  snprintf(name, 200, "fdEdxPosDaughter%s", aName);
  fdEdxPosDaughter = new TH2D(name,"dEdx of positive daughters",200, 0.1, 4.0, 250, 0.0, 500.0);
  snprintf(name, 200, "fdEdxNegDaughter%s", aName);
  fdEdxNegDaughter = new TH2D(name,"dEdx of negative daughters",200, 0.1, 4.0, 250, 0.0, 500.0);
  snprintf(name, 200, "fTOFtimePosDaughter%s", aName);
  fTOFtimePosDaughter = new TH2D(name,"TOF time of positive daughters",190, 0.1, 2.0, 400, -4000.0, 4000.0);
  snprintf(name, 200, " fTOFtimeNegDaughter%s", aName);
  fTOFtimeNegDaughter = new TH2D(name,"TOF time of negative daughters",190, 0.1, 2.0, 400, -4000.0, 4000.0);
  snprintf(name, 200, "fMINVvsPt%s", aName);
  fMINVvsPt = new  TH2D(name, "Invarint mass V0 vs Pt", 50, 0., 5., 1000, 0., 2.);
  snprintf(name, 200, " fnsigmaPosL%s", aName);
  fnsigmaPosL = new TH1D(name,"Number of sigmas of positive Lambda daughters",200,-10,10);
  snprintf(name, 200, " fnsigmaNegL%s", aName);
  fnsigmaNegL = new TH1D(name,"Number of sigmas of negative Lambda daughters",200,-10,10);
  snprintf(name, 200, " fnsigmaPosAL%s", aName);
  fnsigmaPosAL = new TH1D(name,"Number of sigmas of positive AntiLambda daughters",200,-10,10);
  snprintf(name, 200, " fnsigmaNegAL%s", aName);
  fnsigmaNegAL = new TH1D(name,"Number of sigmas of negative AntiLambda daughters",200,-10,10);

  snprintf(name, 200, "POrigin%s", aName);
  fParticleOrigin =  new TH1D(name, "Mothers PDG Codes", 6000, 0.0, 6000.0);

  snprintf(name, 200, "PId%s", aName);
  fParticleId =  new TH1D(name, "Particle PDG Codes", 6000, 0.0, 6000.0);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fDcaPosToPrimVertex->Sumw2();
  fDcaNegToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fDecayLength->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();
  fdEdxPosDaughter->Sumw2();
  fdEdxNegDaughter->Sumw2();
  fTOFtimePosDaughter->Sumw2();
  fTOFtimeNegDaughter->Sumw2();
  fMINVvsPt->Sumw2();
  fnsigmaPosL->Sumw2();
  fnsigmaNegL->Sumw2();
  fnsigmaPosAL->Sumw2();
  fnsigmaNegAL->Sumw2();
}

AliFemtoCutMonitorV0::AliFemtoCutMonitorV0(const AliFemtoCutMonitorV0 &aCut):
  AliFemtoCutMonitor(),
  fLambdaMass(0),
  fAntiLambdaMass(0),
  fK0ShortMass(0),
  fDcaDaughters(0),
  fDcaV0ToPrimVertex(0),
  fDcaPosToPrimVertex(0),
  fDcaNegToPrimVertex(0),
  fCosPointingAngle(0),
  fDecayLength(0),
  fEtaV0(0),
  fPtV0(0),
  fPtPosDaughter(0),
  fPtNegDaughter(0),
  fdEdxPosDaughter(0),
  fdEdxNegDaughter(0),
  fTOFtimePosDaughter(0),
  fTOFtimeNegDaughter(0),
  fMINVvsPt(0),
  fnsigmaPosL(0),
  fnsigmaNegL(0),
  fnsigmaPosAL(0),
  fnsigmaNegAL(0),
  fParticleOrigin(0),
  fParticleId(0)
{
  /// copy constructor

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
  if (fDcaPosToPrimVertex) delete fDcaPosToPrimVertex;
  fDcaPosToPrimVertex = new TH1F(*aCut.fDcaPosToPrimVertex);
  if (fDcaNegToPrimVertex) delete fDcaNegToPrimVertex;
  fDcaNegToPrimVertex = new TH1F(*aCut.fDcaNegToPrimVertex);
  if(fCosPointingAngle) delete fCosPointingAngle;
  fCosPointingAngle = new TH1F(*aCut.fCosPointingAngle);
  if(fDecayLength) delete fDecayLength;
  fDecayLength = new TH1F(*aCut.fDecayLength);
  if(fEtaV0) delete fEtaV0;
  fEtaV0 = new TH1F(*aCut.fEtaV0);
  if(fPtV0) delete fPtV0;
  fPtV0 = new TH1F(*aCut.fPtV0);
  if(fPtPosDaughter) delete fPtPosDaughter;
  fPtPosDaughter = new TH1F(*aCut.fPtPosDaughter);
  if(fPtNegDaughter) delete fPtNegDaughter;
  fPtNegDaughter = new TH1F(*aCut.fPtNegDaughter);
  if(fdEdxPosDaughter) delete fdEdxPosDaughter;
  fdEdxNegDaughter = new TH2D(*aCut.fdEdxNegDaughter);
  if(fdEdxNegDaughter) delete fdEdxNegDaughter;
  fdEdxNegDaughter = new TH2D(*aCut.fdEdxNegDaughter);
  if(fTOFtimePosDaughter) delete fTOFtimePosDaughter;
  fTOFtimePosDaughter = new TH2D(*aCut.fTOFtimePosDaughter);
  if(fTOFtimeNegDaughter) delete fTOFtimeNegDaughter;
  fTOFtimeNegDaughter = new TH2D(*aCut.fTOFtimeNegDaughter);
  if(fMINVvsPt) delete fMINVvsPt;
  fMINVvsPt = new TH2D(*aCut.fMINVvsPt);
  if(fnsigmaPosL) delete fnsigmaPosL;
  fnsigmaPosL = new TH1D(*aCut.fnsigmaPosL);
  if(fnsigmaNegL) delete fnsigmaNegL;
  fnsigmaNegL = new TH1D(*aCut.fnsigmaNegL);
  if(fnsigmaPosAL) delete fnsigmaPosAL;
  fnsigmaPosAL = new TH1D(*aCut.fnsigmaPosAL);
  if(fnsigmaNegAL) delete fnsigmaNegAL;
  fnsigmaNegAL = new TH1D(*aCut.fnsigmaNegAL);

  if (fParticleOrigin) delete fParticleOrigin;
  fParticleOrigin= new TH1D(*aCut.fParticleOrigin);
  if (fParticleId) delete fParticleId;
  fParticleId= new TH1D(*aCut.fParticleId);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fDcaPosToPrimVertex->Sumw2();
  fDcaNegToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fDecayLength->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();
  fdEdxPosDaughter->Sumw2();
  fdEdxNegDaughter->Sumw2();
  fTOFtimePosDaughter->Sumw2();
  fTOFtimeNegDaughter->Sumw2();
  fMINVvsPt->Sumw2();
  fnsigmaPosL->Sumw2();
  fnsigmaNegL->Sumw2();
  fnsigmaPosAL->Sumw2();
  fnsigmaNegAL->Sumw2();
}

AliFemtoCutMonitorV0::~AliFemtoCutMonitorV0()
{
  /// Destructor

  delete fLambdaMass;
  delete fAntiLambdaMass;
  delete fK0ShortMass;
  delete fDcaDaughters;
  delete fDcaV0ToPrimVertex;
  delete fDcaPosToPrimVertex;
  delete fDcaNegToPrimVertex;
  delete fCosPointingAngle;
  delete fDecayLength;
  delete fEtaV0;
  delete fPtV0;
  delete fPtPosDaughter;
  delete fPtNegDaughter;
  delete fdEdxPosDaughter;
  delete fdEdxNegDaughter;
  delete fTOFtimePosDaughter;
  delete fTOFtimeNegDaughter;
  delete fMINVvsPt;
  delete fnsigmaPosL;
  delete fnsigmaNegL;
  delete fnsigmaPosAL;
  delete fnsigmaNegAL;

  delete fParticleOrigin;
  delete fParticleId;
}

AliFemtoCutMonitorV0& AliFemtoCutMonitorV0::operator=(const AliFemtoCutMonitorV0& aCut)
{
  /// assignment operator

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
  if (fDcaPosToPrimVertex) delete fDcaPosToPrimVertex;
  fDcaPosToPrimVertex = new TH1F(*aCut.fDcaPosToPrimVertex);
  if (fDcaNegToPrimVertex) delete fDcaNegToPrimVertex;
  fDcaNegToPrimVertex = new TH1F(*aCut.fDcaNegToPrimVertex);
  if(fCosPointingAngle) delete fCosPointingAngle;
  fCosPointingAngle = new TH1F(*aCut.fCosPointingAngle);
  if(fDecayLength) delete fDecayLength;
  fDecayLength = new TH1F(*aCut.fDecayLength);
  if(fEtaV0) delete fEtaV0;
  fEtaV0 = new TH1F(*aCut.fEtaV0);
  if(fPtV0) delete fPtV0;
  fPtV0 = new TH1F(*aCut.fPtV0);
  if(fPtPosDaughter) delete fPtPosDaughter;
  fPtPosDaughter = new TH1F(*aCut.fPtPosDaughter);
  if(fPtNegDaughter) delete fPtNegDaughter;
  fPtNegDaughter = new TH1F(*aCut.fPtNegDaughter);
  if(fdEdxPosDaughter) delete fdEdxPosDaughter;
  fdEdxNegDaughter = new TH2D(*aCut.fdEdxNegDaughter);
  if(fdEdxNegDaughter) delete fdEdxNegDaughter;
  fdEdxNegDaughter = new TH2D(*aCut.fdEdxNegDaughter);
  if(fTOFtimePosDaughter) delete fTOFtimePosDaughter;
  fTOFtimePosDaughter = new TH2D(*aCut.fTOFtimePosDaughter);
  if(fTOFtimeNegDaughter) delete fTOFtimeNegDaughter;
  fTOFtimeNegDaughter = new TH2D(*aCut.fTOFtimeNegDaughter);
  if(fMINVvsPt) delete fMINVvsPt;
  fMINVvsPt = new TH2D(*aCut.fMINVvsPt);
  if(fnsigmaPosL) delete fnsigmaPosL;
  fnsigmaPosL = new TH1D(*aCut.fnsigmaPosL);
  if(fnsigmaNegL) delete fnsigmaNegL;
  fnsigmaNegL = new TH1D(*aCut.fnsigmaNegL);
  if(fnsigmaPosAL) delete fnsigmaPosAL;
  fnsigmaPosAL = new TH1D(*aCut.fnsigmaPosAL);
  if(fnsigmaNegAL) delete fnsigmaNegAL;
  fnsigmaNegAL = new TH1D(*aCut.fnsigmaNegAL);

  if (fParticleOrigin) delete fParticleOrigin;
  fParticleOrigin= new TH1D(*aCut.fParticleOrigin);
  if (fParticleId) delete fParticleId;
  fParticleId= new TH1D(*aCut.fParticleId);

  fLambdaMass->Sumw2();
  fAntiLambdaMass->Sumw2();
  fK0ShortMass->Sumw2();
  fDcaDaughters->Sumw2();
  fDcaV0ToPrimVertex->Sumw2();
  fDcaPosToPrimVertex->Sumw2();
  fDcaNegToPrimVertex->Sumw2();
  fCosPointingAngle->Sumw2();
  fDecayLength->Sumw2();
  fEtaV0->Sumw2();
  fPtPosDaughter->Sumw2();
  fPtNegDaughter->Sumw2();
  fdEdxPosDaughter->Sumw2();
  fdEdxNegDaughter->Sumw2();
  fTOFtimePosDaughter->Sumw2();
  fTOFtimeNegDaughter->Sumw2();
  fMINVvsPt->Sumw2();
  fnsigmaPosL->Sumw2();
  fnsigmaNegL->Sumw2();
  fnsigmaPosAL->Sumw2();
  fnsigmaNegAL->Sumw2();

  return *this;
}

AliFemtoString AliFemtoCutMonitorV0::Report(){
  /// Prepare report from the execution

  string stemp = "*** AliFemtoCutMonitorV0 report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoCutMonitorV0::Fill(const AliFemtoV0* aV0)
{
  /// Fill momentum resolution histograms for the particle

  fLambdaMass->Fill(aV0->MassLambda());
  fAntiLambdaMass->Fill(aV0->MassAntiLambda());
  fK0ShortMass->Fill(aV0->MassK0Short());
  fDcaDaughters->Fill(aV0->DcaV0Daughters());
  fDcaV0ToPrimVertex->Fill(aV0->DcaV0ToPrimVertex());
  fDcaPosToPrimVertex->Fill(aV0->DcaPosToPrimVertex());
  fDcaNegToPrimVertex->Fill(aV0->DcaNegToPrimVertex());
  fCosPointingAngle->Fill(aV0->CosPointingAngle());
  fDecayLength->Fill(aV0->DecayLengthV0());
  fEtaV0->Fill(aV0->EtaV0());
  fPtV0->Fill(aV0->PtV0());
  fPtPosDaughter->Fill(aV0->PtPos());
  fPtNegDaughter->Fill(aV0->PtNeg());
  fdEdxPosDaughter->Fill(aV0->GetTPCMomentumPos(),aV0->DedxPos());
  fdEdxNegDaughter->Fill(aV0->GetTPCMomentumNeg(),aV0->DedxNeg());

  fTOFtimePosDaughter->Fill(aV0->PtPos(),aV0->TOFProtonTimePos()); //true only for lambdas
  fTOFtimeNegDaughter->Fill(aV0->PtNeg(),aV0->TOFPionTimeNeg());
  fMINVvsPt->Fill(aV0->PtV0(),aV0->MassLambda());
  fnsigmaPosL->Fill(aV0->PosNSigmaTPCP());
  fnsigmaNegL->Fill(aV0->NegNSigmaTPCPi());
  fnsigmaNegAL->Fill(aV0->NegNSigmaTPCP());
  fnsigmaPosAL->Fill(aV0->PosNSigmaTPCPi());

  AliFemtoModelHiddenInfo *tInfo = (AliFemtoModelHiddenInfo*)aV0->GetHiddenInfo();
  if(tInfo!=NULL) {
    Int_t partID = TMath::Abs(tInfo->GetPDGPid());
    Int_t motherID = TMath::Abs(tInfo->GetMotherPdgCode());

    fParticleId->Fill(partID);
    fParticleOrigin->Fill(motherID);
  }
}

void AliFemtoCutMonitorV0::Write()
{
  /// Write out the relevant histograms

  fLambdaMass->Write();
  fAntiLambdaMass->Write();
  fK0ShortMass->Write();
  fDcaDaughters->Write();
  fDcaV0ToPrimVertex->Write();
  fDcaPosToPrimVertex->Write();
  fDcaNegToPrimVertex->Write();
  fCosPointingAngle->Write();
  fDecayLength->Write();
  fEtaV0->Write();
  fPtV0->Write();
  fPtPosDaughter->Write();
  fPtNegDaughter->Write();
  fdEdxPosDaughter->Write();
  fdEdxNegDaughter->Write();
  fTOFtimePosDaughter->Write();
  fTOFtimeNegDaughter->Write();
  fMINVvsPt->Write();
  fnsigmaPosL->Write();
  fnsigmaNegL->Write();
  fnsigmaPosAL->Write();
  fnsigmaNegAL->Write();

  fParticleId->Write();
  fParticleOrigin->Write();
}

TList *AliFemtoCutMonitorV0::GetOutputList()
{
  /// Get the list of histograms to write

  TList *tOutputList = new TList();
  tOutputList->Add(fLambdaMass);
  tOutputList->Add(fAntiLambdaMass);
  tOutputList->Add(fK0ShortMass);
  tOutputList->Add(fDcaDaughters);
  tOutputList->Add(fDcaV0ToPrimVertex);
  tOutputList->Add(fDcaPosToPrimVertex);
  tOutputList->Add(fDcaNegToPrimVertex);
  tOutputList->Add(fCosPointingAngle);
  tOutputList->Add(fDecayLength);
  tOutputList->Add(fEtaV0);
  tOutputList->Add(fPtV0);
  tOutputList->Add(fPtPosDaughter);
  tOutputList->Add(fPtNegDaughter);
  tOutputList->Add(fdEdxPosDaughter);
  tOutputList->Add(fdEdxNegDaughter);
  tOutputList->Add(fTOFtimePosDaughter);
  tOutputList->Add(fTOFtimeNegDaughter);
  tOutputList->Add(fMINVvsPt);
  tOutputList->Add(fnsigmaPosL);
  tOutputList->Add(fnsigmaNegL);
  tOutputList->Add(fnsigmaPosAL);
  tOutputList->Add(fnsigmaNegAL);

  tOutputList->Add(fParticleId);
  tOutputList->Add(fParticleOrigin);

  return tOutputList;
}
