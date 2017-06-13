///
/// \file AliFemtoV0PurityBgdEstimator.cxx
///

#include "AliFemtoV0PurityBgdEstimator.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoV0PurityBgdEstimator);
  /// \endcond
#endif

const float DFLT_MinvMin = 0.0, DFLT_MinvMax = 1.0;
const int DFLT_NbinsMinv = 200;



//____________________________
AliFemtoV0PurityBgdEstimator::AliFemtoV0PurityBgdEstimator():
  AliFemtoV0PurityBgdEstimator("V0PurityBgdEst", DFLT_NbinsMinv, DFLT_MinvMin, DFLT_MinvMax)
{
// no-op
}


//____________________________
AliFemtoV0PurityBgdEstimator::AliFemtoV0PurityBgdEstimator(const char* title,
                                                           const int nbins,
                                                           const float MinvLo,
                                                           const float MinvHi):
  AliFemtoCorrFctn(),
  fTitle(title),
  fNbinsMinv(nbins),
  fMinvLow(MinvLo),
  fMinvHigh(MinvHi),
  fNumeratorWithoutSharedDaughterCut(nullptr),
  fDenominatorWithoutSharedDaughterCut(nullptr),
  fRatioWithoutSharedDaughterCut(nullptr),
  fNumerator(nullptr),
  fDenominator(nullptr),
  fRatio(nullptr),
  fFemtoV0(nullptr),
  fFemtoV0TrackCut(nullptr),
  fCurrentPrimVtx(),
  fRealV0Collection(nullptr),
  fMixedV0Collection(nullptr)
{
  fNumeratorWithoutSharedDaughterCut = new TH1D(TString::Format("NumWithoutSharedDaughterCut%s", fTitle.Data()),
                        "V0Purity - Sig+Bgd; M_{inv}(GeV/c^{2});",
                        nbins, MinvLo, MinvHi);
  fDenominatorWithoutSharedDaughterCut = new TH1D(TString::Format("DenWithoutSharedDaughterCut%s", fTitle.Data()),
                          "V0Purity - Bgd; M_{inv}(GeV/c^{2});",
                          nbins, MinvLo, MinvHi);
  fRatioWithoutSharedDaughterCut = new TH1D(TString::Format("RatWithoutSharedDaughterCut%s", fTitle.Data()),
                    "V0Purity - Ratio; M_{inv}(GeV/c^{2});",
                    nbins, MinvLo, MinvHi);
  fNumerator = new TH1D(TString::Format("Num%s", fTitle.Data()),
                        "V0Purity - Sig+Bgd; M_{inv}(GeV/c^{2});",
                        nbins, MinvLo, MinvHi);
  fDenominator = new TH1D(TString::Format("Den%s", fTitle.Data()),
                          "V0Purity - Bgd; M_{inv}(GeV/c^{2});",
                          nbins, MinvLo, MinvHi);
  fRatio = new TH1D(TString::Format("Rat%s", fTitle.Data()),
                    "V0Purity - Ratio; M_{inv}(GeV/c^{2});",
                    nbins, MinvLo, MinvHi);

  // to enable error bar calculation...
  fNumeratorWithoutSharedDaughterCut->Sumw2();
  fDenominatorWithoutSharedDaughterCut->Sumw2();
  fRatioWithoutSharedDaughterCut->Sumw2();
  fNumerator->Sumw2();
  fDenominator->Sumw2();
  fRatio->Sumw2();


  fFemtoV0 = new AliFemtoV0();
  fFemtoV0TrackCut = new AliFemtoV0TrackCutNSigmaFilter();

  fCurrentPrimVtx = AliFmThreeVector<double>(0.,0.,0.);

  fRealV0Collection = new AliFemtoV0Collection();
  fMixedV0Collection = new AliFemtoV0Collection();
}

//____________________________
AliFemtoV0PurityBgdEstimator::AliFemtoV0PurityBgdEstimator(const AliFemtoV0PurityBgdEstimator& aCorrFctn):
  AliFemtoCorrFctn(aCorrFctn),
  fTitle(aCorrFctn.fTitle),
  fNbinsMinv(aCorrFctn.fNbinsMinv),

  fNumeratorWithoutSharedDaughterCut(aCorrFctn.fNumeratorWithoutSharedDaughterCut ? new TH1D(*aCorrFctn.fNumeratorWithoutSharedDaughterCut) : nullptr),
  fDenominatorWithoutSharedDaughterCut(aCorrFctn.fDenominatorWithoutSharedDaughterCut ? new TH1D(*aCorrFctn.fDenominatorWithoutSharedDaughterCut) : nullptr),
  fNumerator(aCorrFctn.fNumerator ? new TH1D(*aCorrFctn.fNumerator) : nullptr),
  fDenominator(aCorrFctn.fDenominator ? new TH1D(*aCorrFctn.fDenominator) : nullptr),

  fFemtoV0(aCorrFctn.fFemtoV0 ? new AliFemtoV0(*aCorrFctn.fFemtoV0) : nullptr),
  fFemtoV0TrackCut(aCorrFctn.fFemtoV0TrackCut ? new AliFemtoV0TrackCutNSigmaFilter(*aCorrFctn.fFemtoV0TrackCut) : nullptr),

  fRealV0Collection(aCorrFctn.fRealV0Collection ? new AliFemtoV0Collection(*aCorrFctn.fRealV0Collection) : nullptr),
  fMixedV0Collection(aCorrFctn.fMixedV0Collection ? new AliFemtoV0Collection(*aCorrFctn.fMixedV0Collection) : nullptr),

  fMinvLow(aCorrFctn.fMinvLow),
  fMinvHigh(aCorrFctn.fMinvHigh)
{
  // copy constructor
  fRatioWithoutSharedDaughterCut = (aCorrFctn.fRatioWithoutSharedDaughterCut) ? new TH1D(*aCorrFctn.fRatioWithoutSharedDaughterCut) : nullptr;
  fRatio = (aCorrFctn.fRatio) ? new TH1D(*aCorrFctn.fRatio) : nullptr;
}

//_________________________
AliFemtoV0PurityBgdEstimator& AliFemtoV0PurityBgdEstimator::operator=(const AliFemtoV0PurityBgdEstimator& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(aCorrFctn);

  fNbinsMinv = aCorrFctn.fNbinsMinv;
  fMinvLow = aCorrFctn.fMinvLow;
  fMinvHigh = aCorrFctn.fMinvHigh;

  fTitle = aCorrFctn.fTitle;

  if(fNumeratorWithoutSharedDaughterCut) delete fNumeratorWithoutSharedDaughterCut;
    fNumeratorWithoutSharedDaughterCut = new TH1D(*aCorrFctn.fNumeratorWithoutSharedDaughterCut);
  if(fDenominatorWithoutSharedDaughterCut) delete fDenominatorWithoutSharedDaughterCut;
    fDenominatorWithoutSharedDaughterCut = new TH1D(*aCorrFctn.fDenominatorWithoutSharedDaughterCut);
  if(fRatioWithoutSharedDaughterCut) delete fRatioWithoutSharedDaughterCut;
    fRatioWithoutSharedDaughterCut = new TH1D(*aCorrFctn.fRatioWithoutSharedDaughterCut);
  if(fNumerator) delete fNumerator;
    fNumerator = new TH1D(*aCorrFctn.fNumerator);
  if(fDenominator) delete fDenominator;
    fDenominator = new TH1D(*aCorrFctn.fDenominator);
  if(fRatio) delete fRatio;
    fRatio = new TH1D(*aCorrFctn.fRatio);

  if(fFemtoV0) delete fFemtoV0;
    fFemtoV0 = new AliFemtoV0(*aCorrFctn.fFemtoV0);

  if(fFemtoV0TrackCut) delete fFemtoV0TrackCut;
    fFemtoV0TrackCut = new AliFemtoV0TrackCutNSigmaFilter(*aCorrFctn.fFemtoV0TrackCut);

  if(fRealV0Collection) delete fRealV0Collection;
    fRealV0Collection = new AliFemtoV0Collection(*aCorrFctn.fRealV0Collection);

  if(fMixedV0Collection) delete fMixedV0Collection;
    fMixedV0Collection = new AliFemtoV0Collection(*aCorrFctn.fMixedV0Collection);

  return *this;
}

//____________________________
AliFemtoV0PurityBgdEstimator* AliFemtoV0PurityBgdEstimator::Clone()
{
  return(new AliFemtoV0PurityBgdEstimator(*this));
}

//____________________________
AliFemtoV0PurityBgdEstimator::~AliFemtoV0PurityBgdEstimator()
{
  // destructor
  delete fNumeratorWithoutSharedDaughterCut;
  delete fDenominatorWithoutSharedDaughterCut;
  delete fRatioWithoutSharedDaughterCut;
  delete fNumerator;
  delete fDenominator;
  delete fRatio;

  delete fFemtoV0;
  delete fFemtoV0TrackCut;
}

//____________________________
void AliFemtoV0PurityBgdEstimator::ClearV0Collections()
{
  for(AliFemtoV0Iterator iter = fRealV0Collection->begin(); iter != fRealV0Collection->end(); iter++) delete *iter;
  fRealV0Collection->clear();
  for(AliFemtoV0Iterator iter = fMixedV0Collection->begin(); iter != fMixedV0Collection->end(); iter++) delete *iter;
  fMixedV0Collection->clear();
}

//____________________________
void AliFemtoV0PurityBgdEstimator::FillHistograms()
{
  double weight = 1.0;
  double tMinv;
  short tV0Type = fFemtoV0TrackCut->GetParticleType();

  AliFemtoV0SharedDaughterCut shared_daughter_cut;
  AliFemtoV0Collection tCorrectedRealV0Collection = shared_daughter_cut.AliFemtoV0SharedDaughterCutCollection(fRealV0Collection, (AliFemtoV0Cut*)fFemtoV0TrackCut);
  AliFemtoV0Collection tCorrectedMixedV0Collection = shared_daughter_cut.AliFemtoV0SharedDaughterCutCollection(fMixedV0Collection, (AliFemtoV0Cut*)fFemtoV0TrackCut);

  AliFemtoV0* tV0 = NULL;
  for(AliFemtoV0Iterator iter = tCorrectedRealV0Collection.begin(); iter != tCorrectedRealV0Collection.end(); iter++)
  {
    tV0 = *iter;
    if(tV0Type==AliFemtoV0TrackCut::kLambda || tV0Type==AliFemtoV0TrackCut::kLambdaMC) tMinv = tV0->MassLambda();
    else if(tV0Type==AliFemtoV0TrackCut::kAntiLambda || tV0Type==AliFemtoV0TrackCut::kAntiLambdaMC) tMinv = tV0->MassAntiLambda();
    else if(tV0Type==AliFemtoV0TrackCut::kK0s || tV0Type==AliFemtoV0TrackCut::kK0sMC) tMinv = tV0->MassK0Short();
    else tMinv = 0.;

    fNumerator->Fill(tMinv);
  }
  for(AliFemtoV0Iterator iter = tCorrectedMixedV0Collection.begin(); iter != tCorrectedMixedV0Collection.end(); iter++)
  {
    tV0 = *iter;
    if(tV0Type==AliFemtoV0TrackCut::kLambda || tV0Type==AliFemtoV0TrackCut::kLambdaMC) tMinv = tV0->MassLambda();
    else if(tV0Type==AliFemtoV0TrackCut::kAntiLambda || tV0Type==AliFemtoV0TrackCut::kAntiLambdaMC) tMinv = tV0->MassAntiLambda();
    else if(tV0Type==AliFemtoV0TrackCut::kK0s || tV0Type==AliFemtoV0TrackCut::kK0sMC) tMinv = tV0->MassK0Short();
    else tMinv = 0.;

    fDenominator->Fill(tMinv,weight);
  }
}

//____________________________
void AliFemtoV0PurityBgdEstimator::IsSameEvent()
{
  AliFmThreeVector<double> tPrimVtx = fFemtoV0->PrimaryVertex();
  if(tPrimVtx.x() != fCurrentPrimVtx.x() ||
     tPrimVtx.y() != fCurrentPrimVtx.y() ||
     tPrimVtx.z() != fCurrentPrimVtx.z())
  {
    fCurrentPrimVtx.SetX(tPrimVtx.x());
    fCurrentPrimVtx.SetY(tPrimVtx.y());
    fCurrentPrimVtx.SetZ(tPrimVtx.z());

    FillHistograms();
    ClearV0Collections();
  }
}

//____________________________
AliFmThreeVector<double> AliFemtoV0PurityBgdEstimator::ShiftMomentumToDCA(const AliFmHelix& tHelix, const AliFmThreeVector<double>& tMomentum, double tPathLengthToDCA)
{
  //Based on AliFmPhysicalHelix::MomentumAt(double S, double B)
  double xc = tHelix.XCenter();
  double yc = tHelix.YCenter();
  double rx = (tHelix.Y(tPathLengthToDCA)-yc)/(tHelix.Origin().y()-yc);
  double ry = (tHelix.X(tPathLengthToDCA)-xc)/(tHelix.Origin().x()-xc);
  return tMomentum.PseudoProduct(rx,ry,1.0);
}

//____________________________
void AliFemtoV0PurityBgdEstimator::UseCorrectedDaughterHelices(const AliFemtoTrack* tTrackPos, const AliFemtoTrack* tTrackNeg)
{
  //Set the positive and negative daughter momenta to the correct value at their DCA to each other
  //Set the DCA of the daughters to each other
  //Set the decay vertex of the V0
  //The UpdateV0 call below in BuildV0 will handle setting a number of other related attributes

  const AliFmPhysicalHelixD tHelixPosOG = tTrackPos->Helix();
  const AliFmPhysicalHelixD tHelixNegOG = tTrackNeg->Helix();

  //When the helix was built, it assumed that particles come directly from event vertex
  //This is not useful for me, so I build new ones with more correct (hopefully) origins
  //Also, define relative to event vertex, this is necessary in the mixed event case
  const AliFmThreeVector<double> tOriginPos(tTrackPos->XatDCA()-tHelixPosOG.Origin().x(), 
                                            tTrackPos->YatDCA()-tHelixPosOG.Origin().y(), 
                                            tTrackPos->ZatDCA()-tHelixPosOG.Origin().z());
  const AliFmThreeVector<double> tOriginNeg(tTrackNeg->XatDCA()-tHelixNegOG.Origin().x(), 
                                            tTrackNeg->YatDCA()-tHelixNegOG.Origin().y(), 
                                            tTrackNeg->ZatDCA()-tHelixNegOG.Origin().z());

  const AliFmPhysicalHelixD tHelixPos(tHelixPosOG.Curvature(), tHelixPosOG.DipAngle(), tHelixPosOG.Phase(), tOriginPos, tHelixPosOG.H());
  const AliFmPhysicalHelixD tHelixNeg(tHelixNegOG.Curvature(), tHelixNegOG.DipAngle(), tHelixNegOG.Phase(), tOriginNeg, tHelixNegOG.H());

  pair<double, double> tPathLengthsAtDca = tHelixPos.PathLengths(tHelixNeg);

  AliFmThreeVector<double> tPositionAtDcaPos = tHelixPos.At(tPathLengthsAtDca.first);
  AliFmThreeVector<double> tPositionAtDcaNeg = tHelixNeg.At(tPathLengthsAtDca.second);

  AliFmThreeVector<double> tMomPosAtOrigin = tTrackPos->P();
  AliFmThreeVector<double> tMomNegAtOrigin = tTrackNeg->P();

  AliFmThreeVector<double> tMomPosAtDCA = ShiftMomentumToDCA(tHelixPos,tMomPosAtOrigin,tPathLengthsAtDca.first);
  fFemtoV0->SetmomPos(tMomPosAtDCA);
    fFemtoV0->SetmomPosX(tMomPosAtDCA.x());
    fFemtoV0->SetmomPosY(tMomPosAtDCA.y());
    fFemtoV0->SetmomPosZ(tMomPosAtDCA.z());

  AliFmThreeVector<double> tMomNegAtDCA = ShiftMomentumToDCA(tHelixNeg,tMomNegAtOrigin,tPathLengthsAtDca.second);
  fFemtoV0->SetmomNeg(tMomNegAtDCA);
    fFemtoV0->SetmomNegX(tMomNegAtDCA.x());
    fFemtoV0->SetmomNegY(tMomNegAtDCA.y());
    fFemtoV0->SetmomNegZ(tMomNegAtDCA.z());

  AliFmThreeVector<double> tDecayVertex = tPositionAtDcaPos;
  tDecayVertex += tPositionAtDcaNeg;
  tDecayVertex /= 2.0;

  double tDca = tHelixPos.Distance(tPositionAtDcaNeg);

  fFemtoV0->SetdcaV0Daughters(tDca);
  fFemtoV0->SetdecayVertexV0(tDecayVertex);
}


//____________________________
void AliFemtoV0PurityBgdEstimator::BuildV0(AliFemtoPair* aPair)
{
  AliFemtoTrack* tTrackPos = aPair->Track1()->Track();
  AliFemtoTrack* tTrackNeg = aPair->Track2()->Track();

  //Make sure positive and negative daughters have correct charge,
  //If not, switch them
  if(tTrackPos->Charge() < 0 && tTrackNeg->Charge()> 0)
  {
    AliFemtoTrack* tTmpTrack = tTrackPos;
    tTrackPos = tTrackNeg;
    tTrackNeg = new AliFemtoTrack(*tTmpTrack);
    delete tTmpTrack;
  }

  UseCorrectedDaughterHelices(tTrackPos,tTrackNeg);

  //Due to relative shift above in UseCorrectedDaughterHelices, use the primary vertex = {0,0,0} in
  //decay length etc. calculations below.
  //  However, still call AliFemtoV0::SetPrimaryVertex to set primary vertex to actual value, 
  //  as, for now, the position of the event primary vertex is how we distinguish events in IsSameEvent function

  double tEventPrimVtx[3];
  tTrackPos->GetPrimaryVertex(tEventPrimVtx);
  fFemtoV0->SetprimaryVertex(AliFemtoThreeVector(tEventPrimVtx));
  double tPrimVtx[3] = {0., 0., 0.};

  double tDecayLengthV0 = (fFemtoV0->DecayVertexV0X()-tPrimVtx[0])*(fFemtoV0->DecayVertexV0X()-tPrimVtx[0]) +
                          (fFemtoV0->DecayVertexV0Y()-tPrimVtx[1])*(fFemtoV0->DecayVertexV0Y()-tPrimVtx[1]) +
                          (fFemtoV0->DecayVertexV0Z()-tPrimVtx[2])*(fFemtoV0->DecayVertexV0Z()-tPrimVtx[2]);
  tDecayLengthV0 = sqrt(tDecayLengthV0);
  fFemtoV0->SetdecayLengthV0(tDecayLengthV0);

  fFemtoV0->UpdateV0();

  TVector3 tPrimToSecVtxVec;
  tPrimToSecVtxVec.SetX(fFemtoV0->DecayVertexV0X()-tPrimVtx[0]);
  tPrimToSecVtxVec.SetY(fFemtoV0->DecayVertexV0Y()-tPrimVtx[1]);
  tPrimToSecVtxVec.SetZ(fFemtoV0->DecayVertexV0Z()-tPrimVtx[2]);

  TVector3 tMomV0;
  tMomV0.SetX(fFemtoV0->MomV0X());
  tMomV0.SetY(fFemtoV0->MomV0Y());
  tMomV0.SetZ(fFemtoV0->MomV0Z());

  double tDcaV0 = TMath::Abs(tPrimToSecVtxVec.Perp(tMomV0));
  double tCosPointing = TMath::Abs(tPrimToSecVtxVec.Dot(tMomV0)/(tPrimToSecVtxVec.Mag()*tMomV0.Mag()));

  fFemtoV0->SetdcaV0ToPrimVertex(tDcaV0);
  fFemtoV0->SetCosPointingAngle(tCosPointing);

//  fFemtoV0->SetdcaPosToPrimVertex(sqrt(pow(tTrackPos->ImpactD(),2)+pow(tTrackPos->ImpactZ(),2)));  //V0s apparently only use DcaXY here
  fFemtoV0->SetdcaPosToPrimVertex(tTrackPos->ImpactD());
  fFemtoV0->SetPosNSigmaTPCK(tTrackPos->NSigmaTPCK());
  fFemtoV0->SetPosNSigmaTPCPi(tTrackPos->NSigmaTPCPi());
  fFemtoV0->SetPosNSigmaTPCP(tTrackPos->NSigmaTPCP());
  fFemtoV0->SetPosNSigmaTOFK(tTrackPos->NSigmaTOFK());
  fFemtoV0->SetPosNSigmaTOFPi(tTrackPos->NSigmaTOFPi());
  fFemtoV0->SetPosNSigmaTOFP(tTrackPos->NSigmaTOFP());
  fFemtoV0->SetTPCNclsPos(tTrackPos->TPCncls());
  fFemtoV0->SetNdofPos(tTrackPos->TPCchi2());
  fFemtoV0->SetStatusPos(tTrackPos->Flags());
  fFemtoV0->SetidPos(tTrackPos->TrackId());
  fFemtoV0->SetEtaPos(tTrackPos->P().PseudoRapidity());

//  fFemtoV0->SetdcaNegToPrimVertex(sqrt(pow(tTrackNeg->ImpactD(),2)+pow(tTrackNeg->ImpactZ(),2)));  //V0s apparently only use DcaXY here
  fFemtoV0->SetdcaNegToPrimVertex(tTrackNeg->ImpactD());
  fFemtoV0->SetNegNSigmaTPCK(tTrackNeg->NSigmaTPCK());
  fFemtoV0->SetNegNSigmaTPCPi(tTrackNeg->NSigmaTPCPi());
  fFemtoV0->SetNegNSigmaTPCP(tTrackNeg->NSigmaTPCP());
  fFemtoV0->SetNegNSigmaTOFK(tTrackNeg->NSigmaTOFK());
  fFemtoV0->SetNegNSigmaTOFPi(tTrackNeg->NSigmaTOFPi());
  fFemtoV0->SetNegNSigmaTOFP(tTrackNeg->NSigmaTOFP());
  fFemtoV0->SetTPCNclsNeg(tTrackNeg->TPCncls());
  fFemtoV0->SetNdofNeg(tTrackNeg->TPCchi2());
  fFemtoV0->SetStatusNeg(tTrackNeg->Flags());
  fFemtoV0->SetidNeg(tTrackNeg->TrackId());
  fFemtoV0->SetEtaNeg(tTrackNeg->P().PseudoRapidity());
}


//____________________________
AliFemtoString AliFemtoV0PurityBgdEstimator::Report()
{
  // construct report
  TString report = "V0 Purity Background Estimator Report:\n";
  report += TString::Format("Number of entries in numeratorWithoutSharedDaughterCut:\t%E\n", fNumeratorWithoutSharedDaughterCut->GetEntries());
  report += TString::Format("Number of entries in denominatorWithoutSharedDaughterCut:\t%E\n", fDenominatorWithoutSharedDaughterCut->GetEntries());
  report += TString::Format("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += TString::Format("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());
  return AliFemtoString(report);
}

//______________________________
TList* AliFemtoV0PurityBgdEstimator::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumeratorWithoutSharedDaughterCut);
  tOutputList->Add(fDenominatorWithoutSharedDaughterCut);
  tOutputList->Add(fRatioWithoutSharedDaughterCut);
  tOutputList->Add(fNumerator);
  tOutputList->Add(fDenominator);
  tOutputList->Add(fRatio);

  return tOutputList;
}

//_________________________
void AliFemtoV0PurityBgdEstimator::Finish()
{
  fRatioWithoutSharedDaughterCut->Divide(fNumeratorWithoutSharedDaughterCut, fDenominatorWithoutSharedDaughterCut, 1.0, 1.0);
  fRatio->Divide(fNumerator, fDenominator, 1.0, 1.0);
}

//____________________________
void AliFemtoV0PurityBgdEstimator::Write()
{
  // Write out neccessary objects
  fNumeratorWithoutSharedDaughterCut->Write();
  fDenominatorWithoutSharedDaughterCut->Write();
  fRatioWithoutSharedDaughterCut->Write();
  fNumerator->Write();
  fDenominator->Write();
  fRatio->Write();
}


//____________________________
void AliFemtoV0PurityBgdEstimator::AddRealPair(AliFemtoPair* aPair)
{
  // add true pair
  BuildV0(aPair);
  if (!fFemtoV0TrackCut->Pass(fFemtoV0)) return;

  IsSameEvent();
  double tMinv;
  short tV0Type = fFemtoV0TrackCut->GetParticleType();
  if(tV0Type==AliFemtoV0TrackCut::kLambda || tV0Type==AliFemtoV0TrackCut::kLambdaMC) tMinv = fFemtoV0->MassLambda();
  else if(tV0Type==AliFemtoV0TrackCut::kAntiLambda || tV0Type==AliFemtoV0TrackCut::kAntiLambdaMC) tMinv = fFemtoV0->MassAntiLambda();
  else if(tV0Type==AliFemtoV0TrackCut::kK0s || tV0Type==AliFemtoV0TrackCut::kK0sMC) tMinv = fFemtoV0->MassK0Short();
  else tMinv = 0.;

  fRealV0Collection->push_back(new AliFemtoV0(*fFemtoV0));
  fNumeratorWithoutSharedDaughterCut->Fill(tMinv);
}

//____________________________
void AliFemtoV0PurityBgdEstimator::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
  BuildV0(aPair);
  if (!fFemtoV0TrackCut->Pass(fFemtoV0)) return;

  IsSameEvent();
  double weight = 1.0;
  double tMinv;
  short tV0Type = fFemtoV0TrackCut->GetParticleType();
  if(tV0Type==AliFemtoV0TrackCut::kLambda || tV0Type==AliFemtoV0TrackCut::kLambdaMC) tMinv = fFemtoV0->MassLambda();
  else if(tV0Type==AliFemtoV0TrackCut::kAntiLambda || tV0Type==AliFemtoV0TrackCut::kAntiLambdaMC) tMinv = fFemtoV0->MassAntiLambda();
  else if(tV0Type==AliFemtoV0TrackCut::kK0s || tV0Type==AliFemtoV0TrackCut::kK0sMC) tMinv = fFemtoV0->MassK0Short();
  else tMinv = 0.;

  fMixedV0Collection->push_back(new AliFemtoV0(*fFemtoV0));
  fDenominatorWithoutSharedDaughterCut->Fill(tMinv,weight);

}



