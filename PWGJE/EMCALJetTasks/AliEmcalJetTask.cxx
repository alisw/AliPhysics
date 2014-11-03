//
// Emcal jet finder task.
//
// Authors: C.Loizides, S.Aiola, M. Verweij

#include <vector>
#include "AliEmcalJetTask.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"
#include "AliMCEvent.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliRhoParameter.h"
using std::cout;
using std::endl;
using std::cerr;

ClassImp(AliEmcalJetTask)

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask() :
  AliAnalysisTaskSE("AliEmcalJetTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fJetsSubName(""),
  fJetType(kNone),
  fConstSel(0),
  fMCConstSel(0),
  fMarkConst(kFALSE),
  fRadius(0.4),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-0.9),
  fEtaMax(+0.9),
  fMinJetArea(0.005),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fMinMCLabel(0),
  fRecombScheme(fastjet::pt_scheme),
  fTrackEfficiency(1.),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fIsInit(0),
  fLocked(0),
  fIsPSelSet(0),
  fIsMcPart(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fCodeDebug(kFALSE),
  fPionMassClusters(kFALSE),
  fDoGenericSubtractionJetMass(kFALSE),
  fDoGenericSubtractionGR(kFALSE),
  fDoGenericSubtractionExtraJetShapes(kFALSE),
  fDoConstituentSubtraction(kFALSE),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(0),
  fRhom(0),
  fRMax(0.4),
  fDRStep(0.04),
  fPtMinGR(40.),
  fJets(0),
  fJetsSub(0),
  fEvent(0),
  fTracks(0),
  fClus(0),
  fRhoParam(0),
  fRhomParam(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliEmcalJetTask::AliEmcalJetTask(const char *name) :
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fJetsSubName(""),
  fJetType(kAKT|kFullJet|kRX1Jet),
  fConstSel(0),
  fMCConstSel(0),
  fMarkConst(kFALSE),
  fRadius(0.4),
  fMinJetTrackPt(0.15),
  fMinJetClusPt(0.15),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-0.9),
  fEtaMax(+0.9),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10),
  fJetPhiMax(+10),
  fJetEtaMin(-1),
  fJetEtaMax(+1),
  fGhostArea(0.005),
  fMinMCLabel(0),
  fRecombScheme(fastjet::pt_scheme),
  fTrackEfficiency(1.),
  fMCFlag(0),
  fGeneratorIndex(-1),
  fIsInit(0),
  fLocked(0),
  fIsPSelSet(0),
  fIsMcPart(0),
  fIsEmcPart(0),
  fLegacyMode(kFALSE),
  fCodeDebug(kFALSE),
  fPionMassClusters(kFALSE),
  fDoGenericSubtractionJetMass(kFALSE),
  fDoGenericSubtractionGR(kFALSE),
  fDoGenericSubtractionExtraJetShapes(kFALSE),
  fDoConstituentSubtraction(kFALSE),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(0),
  fRhom(0),
  fRMax(0.4),
  fDRStep(0.04),
  fPtMinGR(40.),
  fJets(0),
  fJetsSub(0),
  fEvent(0),
  fTracks(0),
  fClus(0),
  fRhoParam(0),
  fRhomParam(0)
{
  // Standard constructor.

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliEmcalJetTask::~AliEmcalJetTask()
{
  // Destructor
}

//________________________________________________________________________
void AliEmcalJetTask::UserCreateOutputObjects()
{
  // Create user objects.

  fJets = new TClonesArray("AliEmcalJet");
  fJets->SetName(fJetsName);

  if(!fJetsSubName.IsNull()) {
    fJetsSub = new TClonesArray("AliEmcalJet");
    fJetsSub->SetName(fJetsSubName);
  }
}

//________________________________________________________________________
void AliEmcalJetTask::UserExec(Option_t *)
{
  // Main loop, called for each event.
  if (!fIsInit) {
    if (!DoInit())
      return;
    fIsInit = kTRUE;
  }

  // clear the jet array (normally a null operation)
  fJets->Delete();
  if(fJetsSub) fJetsSub->Delete();

  FindJets();
}

//________________________________________________________________________
void AliEmcalJetTask::Terminate(Option_t *)
{
  // Called once at the end of the analysis.
}

//________________________________________________________________________
void AliEmcalJetTask::FindJets()
{
  // Find jets.
  if (!fTracks && !fClus){
    cout << "WARNING NO TRACKS OR CLUSTERS:"  <<endl;
    return;
  }

 if (fRhoParam)
    fRho = fRhoParam->GetVal();
  if (fRhomParam)
    fRhom = fRhomParam->GetVal();

  TString name("kt");
  fastjet::JetAlgorithm jalgo(fastjet::kt_algorithm);
  if ((fJetType & kAKT) != 0) {
    name  = "antikt";
    jalgo = fastjet::antikt_algorithm;
    AliDebug(1,"Using AKT algorithm");
  }
  else {
    AliDebug(1,"Using KT algorithm");
  }

  if ((fJetType & kR020Jet) != 0)
    fRadius = 0.2;
  else if ((fJetType & kR030Jet) != 0)
    fRadius = 0.3;
  else if ((fJetType & kR040Jet) != 0)
    fRadius = 0.4;

  // setup fj wrapper
  AliFJWrapper fjw(name, name);
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetGhostArea(fGhostArea);
  fjw.SetR(fRadius);
  fjw.SetAlgorithm(jalgo);
  fjw.SetRecombScheme(static_cast<fastjet::RecombinationScheme>(fRecombScheme));
  fjw.SetMaxRap(fEtaMax);
  fjw.Clear();

  // get primary vertex
  Double_t vertex[3] = {0, 0, 0};
  if(fEvent->GetPrimaryVertex()) fEvent->GetPrimaryVertex()->GetXYZ(vertex);

  AliDebug(2,Form("Jet type = %d", fJetType));

  if ((fIsMcPart || ((fJetType & kFullJet) != 0) || ((fJetType & kChargedJet) != 0)) && fTracks) {
    const Int_t Ntracks = fTracks->GetEntries();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVParticle *t = static_cast<AliVParticle*>(fTracks->At(iTracks));
      if (!t)
        continue;
      if (fIsMcPart) {
	if (((fJetType & kChargedJet) != 0) && (t->Charge() == 0)) {
	  AliDebug(2,Form("Skipping track %d because it is neutral.", iTracks));
	  continue;
	}
	if (((fJetType & kNeutralJet) != 0) && (t->Charge() != 0)) {
	  AliDebug(2,Form("Skipping track %d because it is charged.", iTracks));
	  continue;
	}
      }
      if (fIsMcPart || TMath::Abs(t->GetLabel()) > fMinMCLabel) {
	if (t->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	  AliDebug(2,Form("Skipping track %d because it does not match the bit mask (%d, %d)", iTracks, fMCConstSel, t->TestBits(TObject::kBitMask)));
	  continue;
	}
	else {
	  AliDebug(2,Form("Track %d matches the bit mask (%d, %d)", iTracks, fMCConstSel, t->TestBits(TObject::kBitMask)));
	}
      }
      else {
	if (t->TestBits(fConstSel) != (Int_t)fConstSel) {
	  AliDebug(2,Form("Skipping track %d because it does not match the bit mask (%d, %d)", iTracks, fConstSel, t->TestBits(TObject::kBitMask)));
	  continue;
	}
	else {
	  AliDebug(2,Form("Track %d matches the bit mask (%d, %d)", iTracks, fConstSel, t->TestBits(TObject::kBitMask)));
	}
      }
      if ((t->GetFlag() & fMCFlag) != fMCFlag) {
	AliDebug(2,Form("Skipping track %d because it does not match the MC flags", iTracks));
	continue;
      }
      if (fGeneratorIndex >= 0 && t->GetGeneratorIndex() != fGeneratorIndex) {
	AliDebug(2,Form("Skipping track %d because it does not match the MC generator index", iTracks));
	continue;
      }
      if (t->Pt() < fMinJetTrackPt)
        continue;
      Double_t eta = t->Eta();
      Double_t phi = t->Phi();
      if ((eta<fEtaMin) || (eta>fEtaMax) ||
          (phi<fPhiMin) || (phi>fPhiMax))
        continue;

      // artificial inefficiency
      if (fTrackEfficiency < 1.) {
	Double_t rnd = gRandom->Rndm();
	if (fTrackEfficiency < rnd) {
	  AliDebug(2,Form("Track %d rejected due to artificial tracking inefficiency", iTracks));
	  continue;
	}
      }

      // offset of 100 for consistency with cluster ids
      AliDebug(2,Form("Track %d accepted (label = %d, pt = %f)", iTracks, t->GetLabel(), t->Pt()));
      fjw.AddInputVector(t->Px(), t->Py(), t->Pz(), t->E(), iTracks + 100);
    }
  }

  if ((((fJetType & kFullJet) != 0) || ((fJetType & kNeutralJet) != 0)) && fClus) {
    const Int_t Nclus = fClus->GetEntries();
    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = 0;
      Double_t cEta=0,cPhi=0,cPt=0;
      Double_t cPx=0,cPy=0,cPz=0;
      if (fIsEmcPart) {
	AliEmcalParticle *ep = static_cast<AliEmcalParticle*>(fClus->At(iClus));
	if (!ep)
	  continue;

	c = ep->GetCluster();
	if (!c)
	  continue;

	if (c->GetLabel() > fMinMCLabel) {
	  if (ep->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	    AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fMCConstSel, ep->TestBits(TObject::kBitMask)));
	    continue;
	  }
	  else {
	    AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fMCConstSel, ep->TestBits(TObject::kBitMask)));
	  }
	}
	else {
	  if (ep->TestBits(fConstSel) != (Int_t)fConstSel) {
	    AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fConstSel, ep->TestBits(TObject::kBitMask)));
	    continue;
	  }
	  else {
	    AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fConstSel, ep->TestBits(TObject::kBitMask)));
	  }
	}

	cEta = ep->Eta();
	cPhi = ep->Phi();
	cPt  = ep->Pt();
	cPx  = ep->Px();
	cPy  = ep->Py();
	cPz  = ep->Pz();
      } else {
	c = static_cast<AliVCluster*>(fClus->At(iClus));
	if (!c)
	  continue;

	if (c->GetLabel() > fMinMCLabel) {
	  if (c->TestBits(fMCConstSel) != (Int_t)fMCConstSel) {
	    AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fMCConstSel, c->TestBits(TObject::kBitMask)));
	    continue;
	  }
	  else {
	    AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fMCConstSel, c->TestBits(TObject::kBitMask)));
	  }
	}
	else {
	  if (c->TestBits(fConstSel) != (Int_t)fConstSel) {
	    AliDebug(2,Form("Skipping cluster %d because it does not match the bit mask (%d, %d)", iClus, fConstSel, c->TestBits(TObject::kBitMask)));
	    continue;
	  }
	  else {
	    AliDebug(2,Form("Cluster %d matches the bit mask (%d, %d)", iClus, fConstSel, c->TestBits(TObject::kBitMask)));
	  }
	}

	TLorentzVector nP;
	c->GetMomentum(nP, vertex);
	cEta = nP.Eta();
	cPhi = nP.Phi();
	cPt  = nP.Pt();
	cPx  = nP.Px();
	cPy  = nP.Py();
	cPz  = nP.Pz();
      }
      if (!c->IsEMCAL())
	continue;
      if (cPt < fMinJetClusPt)
	continue;
      if ((cEta<fEtaMin) || (cEta>fEtaMax) ||
	  (cPhi<fPhiMin) || (cPhi>fPhiMax))
	continue;
      // offset of 100 to skip ghost particles uid = -1
      AliDebug(2,Form("Cluster %d accepted (label = %d)", iClus, c->GetLabel()));
      Double_t e = TMath::Sqrt(cPx*cPx+cPy*cPy+cPz*cPz);
      if(fPionMassClusters) e = TMath::Sqrt(cPx*cPx+cPy*cPy+cPz*cPz + 0.13957*0.13957); //MV: dirty, need better solution
      fjw.AddInputVector(cPx, cPy, cPz, e, -iClus - 100);
      //      fjw.AddInputVector(cPx, cPy, cPz, TMath::Sqrt(cPx*cPx+cPy*cPy+cPz*cPz), -iClus - 100);
    }
  }

  // setting legacy mode
  if (fLegacyMode) {
    fjw.SetLegacyMode(kTRUE);
  }

  // run jet finder
  fjw.Run();

  //run generic subtractor
  if(fDoGenericSubtractionJetMass) {
    fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
    fjw.DoGenericSubtractionJetMass();
  }

  if(fDoGenericSubtractionExtraJetShapes) {
    fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
    fjw.DoGenericSubtractionJetAngularity();
    fjw.DoGenericSubtractionJetpTD();
    fjw.DoGenericSubtractionJetCircularity();
    fjw.DoGenericSubtractionJetSigma2();
    fjw.DoGenericSubtractionJetConstituent();
    fjw.DoGenericSubtractionJetLeSub();
  }

  //run constituent subtractor
  if(fDoConstituentSubtraction) {
    fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
    fjw.DoConstituentSubtraction();
  }

  // get geometry
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom) {
    AliFatal(Form("%s: Can not create geometry", GetName()));
    return;
  }

  // loop over fastjet jets
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  AliDebug(1,Form("%d jets found", (Int_t)jets_incl.size()));
  for (UInt_t ijet=0, jetCount=0; ijet<jets_incl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    AliDebug(3,Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

    if (jets_incl[ij].perp()<fMinJetPt)
      continue;
    if (fjw.GetJetArea(ij)<fMinJetArea)
      continue;
    if ((jets_incl[ij].eta()<fJetEtaMin) || (jets_incl[ij].eta()>fJetEtaMax) ||
	(jets_incl[ij].phi()<fJetPhiMin) || (jets_incl[ij].phi()>fJetPhiMax))
      continue;

    AliEmcalJet *jet = new ((*fJets)[jetCount])
      AliEmcalJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());
    jet->SetLabel(ij);

    //do generic subtraction if requested
#ifdef FASTJET_VERSION
    if(fDoGenericSubtractionJetMass) {
      std::vector<fastjet::contrib::GenericSubtractorInfo> jetMassInfo = fjw.GetGenSubtractorInfoJetMass();
      Int_t n = (Int_t)jetMassInfo.size();
      if(n>ij && n>0) {
        jet->SetFirstDerivative(jetMassInfo[ij].first_derivative());
        jet->SetSecondDerivative(jetMassInfo[ij].second_derivative());
        jet->SetFirstOrderSubtracted(jetMassInfo[ij].first_order_subtracted());
        jet->SetSecondOrderSubtracted(jetMassInfo[ij].second_order_subtracted());
      }
    }

    //here do generic subtraction for angular structure function
    Double_t ptcorr = jets_incl[ij].perp()-fjw.GetJetArea(ij)*fRho;
    fRMax = fRadius+0.2;
    if(fDoGenericSubtractionGR && ptcorr>fPtMinGR) {
      fjw.SetUseExternalBkg(fUseExternalBkg,fRho,fRhom);
      fjw.SetRMaxAndStep(fRMax,fDRStep);
      fjw.DoGenericSubtractionGR(ij);
      std::vector<double> num = fjw.GetGRNumerator();
      std::vector<double> den = fjw.GetGRDenominator();
      std::vector<double> nums = fjw.GetGRNumeratorSub();
      std::vector<double> dens = fjw.GetGRDenominatorSub();
      //pass this to AliEmcalJet
      jet->SetGRNumSize(num.size());
      jet->SetGRDenSize(den.size());
      jet->SetGRNumSubSize(nums.size());
      jet->SetGRDenSubSize(dens.size());
      Int_t nsize = (Int_t)num.size();
      for(Int_t g = 0; g<nsize; ++g) {
        jet->AddGRNumAt(num[g],g);
        jet->AddGRNumSubAt(nums[g],g);
      }
      Int_t dsize = (Int_t)den.size();
      for(Int_t g = 0; g<dsize; ++g) {
        jet->AddGRDenAt(den[g],g);
        jet->AddGRDenSubAt(dens[g],g);
      }
    }

   if(fDoGenericSubtractionExtraJetShapes){
      std::vector<fastjet::contrib::GenericSubtractorInfo> jetAngularityInfo = fjw.GetGenSubtractorInfoJetAngularity();
      Int_t na = (Int_t)jetAngularityInfo.size();
      if(na>ij && na>0) {
	jet->SetFirstDerivativeAngularity(jetAngularityInfo[ij].first_derivative());
	jet->SetSecondDerivativeAngularity(jetAngularityInfo[ij].second_derivative());
	jet->SetFirstOrderSubtractedAngularity(jetAngularityInfo[ij].first_order_subtracted());
	jet->SetSecondOrderSubtractedAngularity(jetAngularityInfo[ij].second_order_subtracted());
      }

      std::vector<fastjet::contrib::GenericSubtractorInfo> jetpTDInfo = fjw.GetGenSubtractorInfoJetpTD();
      Int_t np = (Int_t)jetpTDInfo.size();
      if(np>ij && np>0) {
	jet->SetFirstDerivativepTD(jetpTDInfo[ij].first_derivative());
	jet->SetSecondDerivativepTD(jetpTDInfo[ij].second_derivative());
	jet->SetFirstOrderSubtractedpTD(jetpTDInfo[ij].first_order_subtracted());
	jet->SetSecondOrderSubtractedpTD(jetpTDInfo[ij].second_order_subtracted());
      }

      std::vector<fastjet::contrib::GenericSubtractorInfo> jetCircularityInfo = fjw.GetGenSubtractorInfoJetCircularity();
      Int_t nc = (Int_t)jetCircularityInfo.size();
      if(nc>ij && nc>0) {
	jet->SetFirstDerivativeCircularity(jetCircularityInfo[ij].first_derivative());
	jet->SetSecondDerivativeCircularity(jetCircularityInfo[ij].second_derivative());
	jet->SetFirstOrderSubtractedCircularity(jetCircularityInfo[ij].first_order_subtracted());
	jet->SetSecondOrderSubtractedCircularity(jetCircularityInfo[ij].second_order_subtracted());
      }

      std::vector<fastjet::contrib::GenericSubtractorInfo> jetSigma2Info = fjw.GetGenSubtractorInfoJetSigma2();
      Int_t ns = (Int_t)jetSigma2Info.size();
      if(ns>ij && ns>0) {
	jet->SetFirstDerivativeSigma2(jetSigma2Info[ij].first_derivative());
	jet->SetSecondDerivativeSigma2(jetSigma2Info[ij].second_derivative());
	jet->SetFirstOrderSubtractedSigma2(jetSigma2Info[ij].first_order_subtracted());
	jet->SetSecondOrderSubtractedSigma2(jetSigma2Info[ij].second_order_subtracted());
      }


      std::vector<fastjet::contrib::GenericSubtractorInfo> jetConstituentInfo = fjw.GetGenSubtractorInfoJetConstituent();
      Int_t nco= (Int_t)jetConstituentInfo.size();
      if(nco>ij && nco>0) {
	jet->SetFirstDerivativeConstituent(jetConstituentInfo[ij].first_derivative());
	jet->SetSecondDerivativeConstituent(jetConstituentInfo[ij].second_derivative());
	jet->SetFirstOrderSubtractedConstituent(jetConstituentInfo[ij].first_order_subtracted());
	jet->SetSecondOrderSubtractedConstituent(jetConstituentInfo[ij].second_order_subtracted());
      }

      std::vector<fastjet::contrib::GenericSubtractorInfo> jetLeSubInfo = fjw.GetGenSubtractorInfoJetLeSub();
      Int_t nlsub= (Int_t)jetLeSubInfo.size();
      if(nlsub>ij && nlsub>0) {
	jet->SetFirstDerivativeLeSub(jetLeSubInfo[ij].first_derivative());
	jet->SetSecondDerivativeLeSub(jetLeSubInfo[ij].second_derivative());
	jet->SetFirstOrderSubtractedLeSub(jetLeSubInfo[ij].first_order_subtracted());
	jet->SetSecondOrderSubtractedLeSub(jetLeSubInfo[ij].second_order_subtracted());
      }
   }
#endif

    // Fill constituent info
    std::vector<fastjet::PseudoJet> constituents(fjw.GetJetConstituents(ij));
    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfClusters(constituents.size());

    Int_t nt            = 0;
    Int_t nc            = 0;
    Double_t neutralE   = 0;
    Double_t maxCh      = 0;
    Double_t maxNe      = 0;
    Int_t gall          = 0;
    Int_t gemc          = 0;
    Int_t cemc          = 0;
    Int_t ncharged      = 0;
    Int_t nneutral      = 0;
    Double_t mcpt       = 0;
    Double_t emcpt      = 0;

    FillJetConstituents(constituents,jet,vertex,jetCount,nt,nc,maxCh,maxNe,ncharged,nneutral,neutralE,mcpt,cemc,emcpt,gall,gemc,constituents,0);
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SortConstituents();
    jet->SetMaxNeutralPt(maxNe);
    jet->SetMaxChargedPt(maxCh);
    jet->SetNEF(neutralE / jet->E());
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());
    jet->SetNumberOfCharged(ncharged);
    jet->SetNumberOfNeutrals(nneutral);
    jet->SetMCPt(mcpt);
    jet->SetNEmc(cemc);
    jet->SetPtEmc(emcpt);

    if (gall > 0)
      jet->SetAreaEmc(fjw.GetJetArea(ij) * gemc / gall);
    else
      jet->SetAreaEmc(-1);
    if ((jet->Phi() > geom->GetArm1PhiMin() * TMath::DegToRad()) &&
	(jet->Phi() < geom->GetArm1PhiMax() * TMath::DegToRad()) &&
	(jet->Eta() > geom->GetArm1EtaMin()) &&
	(jet->Eta() < geom->GetArm1EtaMax()))
      jet->SetAxisInEmcal(kTRUE);

    AliDebug(2,Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), (Int_t)constituents.size()));
    jetCount++;
  }
  //fJets->Sort();

  //run constituent subtractor if requested
  if(fDoConstituentSubtraction) {
#ifdef FASTJET_VERSION
    if(!fJetsSub) AliWarning(Form("No jet branch to write to for subtracted jets. fJetsSubName: %s",fJetsSubName.Data()));
    else {
      std::vector<fastjet::PseudoJet> jets_sub;
      jets_sub = fjw.GetConstituentSubtrJets();
      AliDebug(1,Form("%d constituent subtracted jets found", (Int_t)jets_sub.size()));
      for (UInt_t ijet=0, jetCount=0; ijet<jets_sub.size(); ++ijet) {
	//Only storing 4-vector and jet area of unsubtracted jet
	AliEmcalJet *jet_sub = new ((*fJetsSub)[jetCount])
	  AliEmcalJet(jets_sub[ijet].perp(), jets_sub[ijet].eta(), jets_sub[ijet].phi(), jets_sub[ijet].m());
	jet_sub->SetLabel(ijet);
	 // Fill constituent info
        std::vector<fastjet::PseudoJet> constituents_unsub(fjw.GetJetConstituents(ijet));
	std::vector<fastjet::PseudoJet> constituents_sub = jets_sub[ijet].constituents();
        jet_sub->SetNumberOfTracks(constituents_sub.size());
	jet_sub->SetNumberOfClusters(constituents_sub.size());
	Int_t nt            = 0;
	Int_t nc            = 0;
	Double_t neutralE   = 0;
	Double_t maxCh      = 0;
	Double_t maxNe      = 0;
	Int_t gall          = 0;
	Int_t gemc          = 0;
	Int_t cemc          = 0;
	Int_t ncharged      = 0;
	Int_t nneutral      = 0;
	Double_t mcpt       = 0;
	Double_t emcpt      = 0;

	FillJetConstituents(constituents_sub,jet_sub,vertex,jetCount,nt,nc,maxCh,maxNe,ncharged,nneutral,neutralE,mcpt,cemc,emcpt,gall,gemc,constituents_unsub,1);
	jet_sub->SetNumberOfTracks(nt);
	jet_sub->SetNumberOfClusters(nc);
	jet_sub->SortConstituents();

	fastjet::PseudoJet area(fjw.GetJetAreaVector(ijet));
	jet_sub->SetArea(area.perp());
	jet_sub->SetAreaEta(area.eta());
	jet_sub->SetAreaPhi(area.phi());
	jet_sub->SetAreaEmc(area.perp());
	jetCount++;
      }
    }
#endif
  } //constituent subtraction
}

//________________________________________________________________________
Bool_t AliEmcalJetTask::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  // Get the leading jets.

  static Float_t pt[9999] = {0};

  const Int_t n = (Int_t)array.size();

  if (n < 1)
    return kFALSE;

  for (Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalJetTask::DoInit()
{
  // Init. Return true if successful.

  if (fTrackEfficiency < 1.) {
    if (gRandom) delete gRandom;
    gRandom = new TRandom3(0);
  }

  // get input collections
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  // get the event
  fEvent = InputEvent();
  if (!fEvent) {
    AliError(Form("%s: Could not retrieve event! Returning", GetName()));
    return 0;
  }

  if (!(fEvent->FindListObject(fJetsSubName)) && fJetsSub)
    fEvent->AddObject(fJetsSub);

  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");
  if (!fTracks && !fTracksName.IsNull()) {
    fTracks = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Pointer to tracks %s == 0", GetName(), fTracksName.Data()));
      return 0;
    }
  }
  if (fTracks) {
    TClass cls(fTracks->GetClass()->GetName());
    if (cls.InheritsFrom("AliMCParticle") || cls.InheritsFrom("AliAODMCParticle"))
      fIsMcPart = 1;
  }

  if (fCaloName == "CaloClusters")
    am->LoadBranch("CaloClusters");
  if (!fClus && !fCaloName.IsNull()) {
    fClus = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fCaloName));
    if (!fClus) {
      AliError(Form("%s: Pointer to clus %s == 0", GetName(), fCaloName.Data()));
      return 0;
    }
  }
  if (fClus) {
    TClass cls(fClus->GetClass()->GetName());
    if (cls.InheritsFrom("AliEmcalParticle"))
      fIsEmcPart = 1;
  }

  if (!fRhoName.IsNull() && !fRhoParam) { // get rho from the event
    fRhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
    if (!fRhoParam) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return 0;
    }
  }
  if (!fRhomName.IsNull() && !fRhomParam) { // get rhom from the event
    fRhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhomName));
    if (!fRhomParam) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhomName.Data()));
      return 0;
    }
  }

  // add jets to event if not yet there
  if (!(fEvent->FindListObject(fJetsName)))
    fEvent->AddObject(fJets);
  else {
    AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fJetsName.Data()));
    return 0;
  }

  return 1;
}

//___________________________________________________________________________________________________________________
void  AliEmcalJetTask::FillJetConstituents(std::vector<fastjet::PseudoJet>& constituents,AliEmcalJet *jet,Double_t vertex[3],UInt_t jetCount,Int_t& nt,Int_t& nc,Double_t& maxCh,Double_t& maxNe,Int_t& ncharged,Int_t& nneutral,Double_t& neutralE,Double_t& mcpt,Int_t& cemc,Double_t& emcpt,Int_t& gall,Int_t& gemc,std::vector<fastjet::PseudoJet>& constituentsunsub,Int_t flagsub){
    nt            = 0;
    nc            = 0;
    neutralE   = 0;
    maxCh      = 0;
    maxNe      = 0;
    gall          = 0;
    gemc       =0;
    cemc          = 0;
    ncharged      = 0;
    nneutral      = 0;
    mcpt       = 0;
    emcpt      = 0;
    Int_t uid   = -1;
   AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
    for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
      if(flagsub==0) uid = constituents[ic].user_index();
      if(flagsub!=0) uid=GetIndexSub(constituents[ic].perp(),constituentsunsub);
      AliDebug(3,Form("Processing constituent %d", uid));
      if ((uid == -1) /*&& (constituents[ic].kt2() < 1e-25)*/) { //ghost particle
        ++gall;
        Double_t gphi = constituents[ic].phi();
        if (gphi<0)
          gphi += TMath::TwoPi();
        gphi *= TMath::RadToDeg();
        Double_t geta = constituents[ic].eta();
        if ((gphi > geom->GetArm1PhiMin()) && (gphi < geom->GetArm1PhiMax()) &&
            (geta > geom->GetArm1EtaMin()) && (geta < geom->GetArm1EtaMax()))
          ++gemc;
      }	else if ((uid > 0) && fTracks) { // track constituent
	Int_t tid = uid - 100;
        AliVParticle *t = static_cast<AliVParticle*>(fTracks->At(tid));
        if (!t) {
	  AliError(Form("Could not find track %d",tid));
          continue;
	}
	if (jetCount < fMarkConst) {
	  AliDebug(2,Form("Marking track %d with bit map %d", tid, fJetType));
	  t->SetBit(fJetType);
	}
        Double_t cEta = t->Eta();
        Double_t cPhi = t->Phi();
        Double_t cPt  = t->Pt();
        Double_t cP   = t->P();
        if (t->Charge() == 0) {
          neutralE += cP;
          ++nneutral;
          if (cPt > maxNe)
            maxNe = cPt;
        } else {
          ++ncharged;
          if (cPt > maxCh)
            maxCh = cPt;
        }
        if (fIsMcPart || TMath::Abs(t->GetLabel()) > fMinMCLabel) // check if MC particle
          mcpt += cPt;

        if (cPhi<0)
          cPhi += TMath::TwoPi();
        cPhi *= TMath::RadToDeg();
        if ((cPhi > geom->GetArm1PhiMin()) && (cPhi < geom->GetArm1PhiMax()) &&
            (cEta > geom->GetArm1EtaMin()) && (cEta < geom->GetArm1EtaMax())) {
          emcpt += cPt;
          ++cemc;
        }
	// cout<<"Adding the track"<<endl;
        jet->AddTrackAt(tid, nt);
        ++nt;
      } else if (fClus) { // cluster constituent
	Int_t cid = -uid - 100;
	AliVCluster *c = 0;
        Double_t cEta=0,cPhi=0,cPt=0,cP=0;
        if (fIsEmcPart) {
          AliEmcalParticle *ep = static_cast<AliEmcalParticle*>(fClus->At(cid));
          if (!ep)
            continue;
          c = ep->GetCluster();
          if (!c)
            continue;
	  if (jetCount < fMarkConst)
	    ep->SetBit(fJetType);
          cEta = ep->Eta();
          cPhi = ep->Phi();
          cPt  = ep->Pt();
          cP   = ep->P();
        } else {
          c = static_cast<AliVCluster*>(fClus->At(cid));
          if (!c)
            continue;
	  if (jetCount < fMarkConst)
	    c->SetBit(fJetType);
          TLorentzVector nP;
          c->GetMomentum(nP, vertex);
          cEta = nP.Eta();
          cPhi = nP.Phi();
          cPt  = nP.Pt();
          cP   = nP.P();
        }

        neutralE += cP;
        if (cPt > maxNe)
          maxNe = cPt;

        if (c->GetLabel() > fMinMCLabel) // MC particle
          mcpt += c->GetMCEnergyFraction() > 1e-6 ? cPt * c->GetMCEnergyFraction() : cPt;

        if (cPhi<0)
          cPhi += TMath::TwoPi();
        cPhi *= TMath::RadToDeg();
        if ((cPhi > geom->GetArm1PhiMin()) && (cPhi < geom->GetArm1PhiMax()) &&
            (cEta > geom->GetArm1EtaMin()) && (cEta < geom->GetArm1EtaMax())) {
          emcpt += cPt;
          ++cemc;
        }

        jet->AddClusterAt(cid, nc);
        ++nc;
        ++nneutral;
      } else {
        AliError(Form("%s: No logical way to end up here.", GetName()));
        continue;
      }
    }
}
//______________________________________________________________________________________
Int_t AliEmcalJetTask::GetIndexSub(Double_t ptsub,std::vector<fastjet::PseudoJet>& constituentsunsub){

    for(UInt_t ii=0;ii<constituentsunsub.size();ii++){
    Double_t ptunsub=constituentsunsub[ii].perp();
    //cout<<ptunsub<<" "<<ptsub<<endl;
    if(ptsub==ptunsub) return constituentsunsub[ii].user_index();

    } return -1;}
//______________________________________________________________________________________
