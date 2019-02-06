/**
 * @Author: Pascal Dillenseger <pascaldillenseger>
 * @Date:   2017-08-09, 17:39:09
 * @Email:  pdillens@cern.ch
 * @Filename: AliDielectronEvtVsTrkHist.cxx
 * @Last modified by:   pascaldillenseger
 * @Last modified time: 2017-10-02, 11:21:56
 */



// Class to create histograms plotting event properties versus track properties e.g. ITS-TPC matching efficiency vs. pT
// Since these are special case histograms the filling of each histogram is kind of treated individually, therefore check if your desired histogram exists in the filling procedure
#include <TAxis.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include <TH1.h>

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliESDEvent.h"
#include "AliDielectronTrackCuts.h"

#include "AliDielectronEvtVsTrkHist.h"

#define ADVM AliDielectronVarManager
#define ADEVTH AliDielectronEvtVsTrkHist

ClassImp(AliDielectronEvtVsTrkHist)


//______________________________________________
AliDielectronEvtVsTrkHist::AliDielectronEvtVsTrkHist(const char *name, const char *title) :
TNamed(name, title),
fEventNumber(-1),
fIsAODEvent(kFALSE),
fSameEvent(kFALSE),
fNObjs(-1),
fNHistos(0),
fNSparses(0),
fNtracksITS(0),
fNtracksTPC(0),
fHistoList(0x0),
fSparseObjs{0x0},
fPIDResponse(0x0),
fMatchEffITS(0x0),
fMatchEffTPC(0x0)
{
  // Default named constructor
}

//______________________________________________
AliDielectronEvtVsTrkHist::~AliDielectronEvtVsTrkHist()
{
  // Default destructor
  if(fMatchEffITS) delete fMatchEffITS;
  if(fMatchEffTPC) delete fMatchEffTPC;
}

//______________________________________________
void AliDielectronEvtVsTrkHist::Init()
{
  if(!fHistoList){
    AliFatal("Histo list from AliDielectronHistos object not setted!");
    return;
  }
  fNObjs = fHistoList->GetEntries();
  TObject *obj;
  for (Int_t iObj = 0; iObj < fNObjs; iObj++) {
    obj = fHistoList->At(iObj);

    if (obj->InheritsFrom(TH1::Class())) {
      if(!ADEVTH::SetHistogramObj(obj))  continue;
      else fNHistos++;
    }

    if (obj->InheritsFrom(THnBase::Class())) {
      if(!ADEVTH::SetSparseObj(obj)) continue;
      else fNSparses++;
    }
  }
}

//______________________________________________
void AliDielectronEvtVsTrkHist::FillHistograms(const AliVEvent *ev)
{
  if (ev->IsA() == AliAODEvent::Class()) fIsAODEvent = kTRUE;
  if(!fIsAODEvent){
    // Fill information for ESDs
    // Currently not implemented therefore an error is printed and the function is returned
    AliError("Event vs. track histograms not implemented for ESDs!");
    return;
  }
    if(fIsAODEvent){
    // Fill information for AODs
      AliAODEvent *event = (AliAODEvent*) ev;
      AliAODHeader *header = (AliAODHeader*) event->GetHeader();

      fSameEvent = fEventNumber == header->GetEventNumberESDFile() ? kTRUE : kFALSE;

      if(fSameEvent)  return;

      fEventNumber = header->GetEventNumberESDFile();
      Int_t nTracks = event->GetNumberOfTracks();
      ADEVTH::SetEventplaneAngles(event);


    //Needed Vars for matching eff
      const Int_t nDimMatchEff = fSparseObjs[ADEVTH::kSparseMatchEffITSTPC]->GetNdimensions();

      Double_t *bin = new Double_t[nDimMatchEff];
      for (Int_t iDim = 0; iDim < nDimMatchEff; iDim++) {
        if(!fSparseIsTrackVar[ADEVTH::kSparseMatchEffITSTPC][iDim]){
          bin[iDim] = ADEVTH::GetVarValueEvent(event, fSparseVars[ADEVTH::kSparseMatchEffITSTPC][iDim]);
        }
      }
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliVTrack *track        = static_cast<AliVTrack*>(event->GetTrack(iTrack));
      if (!track) continue;


      // Fill Matching efficiecy Sparse:
      if (fSparseObjs[kSparseMatchEffITSTPC]) {
        for (Int_t iDim = 0; iDim < nDimMatchEff; iDim++) {
          if(fSparseIsTrackVar[ADEVTH::kSparseMatchEffITSTPC][iDim]){
            bin[iDim] = ADEVTH::GetTrackValue(track, fSparseVars[ADEVTH::kSparseMatchEffITSTPC][iDim]);
          }
        }
        if(ADEVTH::IsSelectedKinematics(track)){
          if(ADEVTH::IsSelectedTPC(track)){
            fMatchEffTPC->Fill(bin);
            fNtracksTPC++;
            if(ADEVTH::IsSelectedITS(track)){
              fMatchEffITS->Fill(bin);
              fNtracksITS++;
            }
            else fMatchEffITS->GetBin(bin); // make sure bin is allocated
          }
        }
      }
    }
    delete [] bin;
  }
}

//______________________________________________
Float_t AliDielectronEvtVsTrkHist::GetVarValueEvent(const AliVEvent *ev, Int_t var)
{
  AliMultSelection *multSelection = 0x0; // no var initialisation in switch

  switch (var) {
    case ADVM::kCentrality:
    case ADVM::kCentralityNew:
      multSelection = (AliMultSelection*) ev->FindListObject("MultSelection");
      if(!multSelection) return -999.;
      else return multSelection->GetMultiplicityPercentile("V0M",kFALSE);
    break;
    case ADVM::kMatchEffITSTPC:
      return 1;
    break;

    default:
      AliError(Form("Variable %s not defined for AliDielectronEvtVsTrkHist", ADVM::GetValueName(var)));
      return -1.;
    break;
  }
}

//______________________________________________
Float_t AliDielectronEvtVsTrkHist::GetTrackValue(AliVTrack *trk, Int_t var)
{
  // Call only after the eventplaneangles have been set
  Float_t rValue = -1.;
  switch (var) {
    case ADVM::kPhi:
      rValue = trk->Phi();
    break;
    case ADVM::kPt:
      rValue = trk->Pt();
    break;
    case ADVM::kQnDeltaPhiTrackTPCrpH2:
      rValue = TVector2::Phi_mpi_pi((trk->Phi() - fEventPlaneAngle[0]));
    break;
    case ADVM::kQnDeltaPhiTrackV0CrpH2:
      rValue = TVector2::Phi_mpi_pi((trk->Phi() - fEventPlaneAngle[1]));
    break;
    default:
      AliError(Form("Variable %d %s not defined for AliDielectronEvtVsTrkHist", var, ADVM::GetValueName(var)));
    break;
  }
  return rValue;
}

//______________________________________________
Bool_t AliDielectronEvtVsTrkHist::SetHistogramObj(TObject *obj)
{
  TString name = obj->GetName();

  // for (Int_t i = 0; i < ADEVTH::kHistoNMaxValues; i++) {
  //   if(name.Contains(ADEVTH::Get3DHistoUniqueNameInfo(i))){
  //     f3DHisoObjs[i] = (THnSparseD*) obj;
  //     return kTRUE;
  //   }
  // }
  AliError(Form("No histogram defined for %s in AliDielectronEvtVsTrkHist!", name.Data()));
  return kFALSE;
}

//______________________________________________
Bool_t AliDielectronEvtVsTrkHist::SetSparseObj(TObject *obj)
{
  TString name = obj->GetName();
  printf("ADEVTH - Setting sparse object %s\n", name.Data());
  for (Int_t i = 0; i < ADEVTH::kSparseNMaxValues; i++) {
    if(name.Contains(ADEVTH::GetSparseUniqueNameInfo(i))){
      fSparseObjs[i] = (THnSparseD*) obj;
      for (Int_t j = 0; j < fSparseObjs[i]->GetNdimensions(); j++) {
        ADEVTH::SetSparseVars(i,j,fSparseObjs[i]->GetAxis(j)->GetUniqueID());
      }
      if(i == ADEVTH::kSparseMatchEffITSTPC){
        fMatchEffITS = (THnSparseD*) fSparseObjs[i]->Clone();
        fMatchEffTPC = (THnSparseD*) fSparseObjs[i]->Clone();
      }
      return kTRUE;
    }
  }
  AliError(Form("No sparse defined for %s in AliDielectronEvtVsTrkHist!", name.Data()));
  return kFALSE;
}

//______________________________________________
const char* AliDielectronEvtVsTrkHist::GetSparseUniqueNameInfo( Int_t iName)
{
  switch (iName) {
    case ADEVTH::kSparseMatchEffITSTPC:
      return ADVM::GetValueName(ADVM::kMatchEffITSTPC);
    break;
    case ADEVTH::kSparseNMaxValues:
      return "NotExistingSparse";
    break;
  }
  return "NotExistingSparse";
}

//______________________________________________
const char* AliDielectronEvtVsTrkHist::Get3DHistoUniqueNameInfo( Int_t iName)
{
  switch (iName) {
    case ADEVTH::k3DHistoNMaxValues:
      return "NotExisting3DHistogram";
    break;
  }
  return "NotExisting3DHistogram";
}

//______________________________________________
void AliDielectronEvtVsTrkHist::SetSparseVars(Int_t nSparse, Int_t nAxis, Int_t type)
{
  switch (type) {
    case ADVM::kCentrality:
    case ADVM::kCentralityNew:
    case ADVM::kMatchEffITSTPC:
      fSparseIsTrackVar[nSparse][nAxis] = kFALSE;
    break;
    case ADVM::kPhi:
    case ADVM::kPt:
    case ADVM::kQnDeltaPhiTrackTPCrpH2:
    case ADVM::kQnDeltaPhiTrackV0CrpH2:
        fSparseIsTrackVar[nSparse][nAxis] = kTRUE;
    break;

    default:
      AliFatal(Form("Variable %s not defined yet please check the code!", ADVM::GetValueName(type)));
      return;
    break;
  }
  // printf("Variable added to sparse %d %s\n", nSparse, ADVM::GetValueName(type));
  fSparseVars[nSparse][nAxis] = type;
}

//______________________________________________
void AliDielectronEvtVsTrkHist::SetEventplaneAngles(const AliVEvent *ev)
{
  TString epDetector[2] = {"TPC","VZEROC"};
  fEventPlaneAngle[0] = -999.;
  fEventPlaneAngle[1] = -999.;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man)
  if( AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
      dynamic_cast<AliAnalysisTaskFlowVectorCorrections*> (man->GetTask("FlowQnVectorCorrections")) )
    if(flowQnVectorTask != NULL){
      AliQnCorrectionsManager *flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
      TList *qnlist = flowQnVectorMgr->GetQnVectorList();
      if(qnlist != NULL){
        for (Int_t i = 0; i < 2; i++) {
          AliQnCorrectionsQnVector *qVecQnFramework = AliDielectronQnEPcorrection::GetQnVectorFromList( qnlist, epDetector[i].Data(), "latest", "latest" );
          if(qVecQnFramework != NULL){
            TVector2 qVector( qVecQnFramework->Qx(2), qVecQnFramework->Qy(2) );
            fEventPlaneAngle[i] = TVector2::Phi_mpi_pi(qVector.Phi())/2;
          }
        }
      }
    }
}

//______________________________________________
Bool_t AliDielectronEvtVsTrkHist::IsSelectedKinematics(AliVTrack *trk)
{
  Float_t impactParXY = -99.;
  Float_t impactParZ = -99.;
  Float_t eta = -99.;

  if(trk->IsA() == AliESDtrack::Class()) trk->GetImpactParameters(impactParXY, impactParZ);
  else{
    AliAODTrack *trackAOD = (AliAODTrack*) trk;
    Double_t xyz[2] = {-100.,-100.};
    Double_t dcaRes[3] = {-100.,-100.,-100.};
    AliDielectronVarManager::GetDCA(trackAOD, xyz, dcaRes);
    impactParXY = xyz[0];
    impactParZ = xyz[1];
  }

  if(TMath::Abs(impactParXY) > 1.0) return kFALSE;
  if(TMath::Abs(impactParZ) > 3.0) return kFALSE;
  eta = trk->Eta();
  if(TMath::Abs(eta) > 0.9) return kFALSE;

  return kTRUE;
}

//______________________________________________
Bool_t AliDielectronEvtVsTrkHist::IsSelectedITS(AliVTrack *trk)
{
  // TODO the cuts should be made available to be set by a setter
  if(!trk) return kFALSE;

  AliDielectronTrackCuts trkCutsITS;
  trkCutsITS.SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  trkCutsITS.SetRequireITSRefit(kTRUE);

  if(!trkCutsITS.IsSelected(trk)) return kFALSE;

  Float_t itsChi2Cl = -99.;

  if(trk->IsA() == AliAODTrack::Class()){
    AliAODTrack *aodtrk = (AliAODTrack*) trk;
    itsChi2Cl = ( aodtrk->GetITSNcls()>0 ) ? aodtrk->GetITSchi2() / aodtrk->GetITSNcls() : 0;
  }
  if(trk->IsA() == AliESDtrack::Class()){
    AliESDtrack *esdtrk = (AliESDtrack*) trk;
    itsChi2Cl = ( esdtrk->GetNcls(0)>0 ) ? esdtrk->GetITSchi2() / esdtrk->GetITSNcls() : 0;
  }

  if(0. >= itsChi2Cl || itsChi2Cl >= 5.) return kFALSE;

  return kTRUE;
}

//______________________________________________
Bool_t AliDielectronEvtVsTrkHist::IsSelectedTPC(AliVTrack *trk)
{
  // TODO the cuts should be made available to be set by a setter
  if(!trk) return kFALSE;

  AliDielectronTrackCuts trkCutsTPC;
  trkCutsTPC.SetRequireTPCRefit(kTRUE);

  if(!trkCutsTPC.IsSelected(trk)) return kFALSE;

  Float_t nClsTPC = -99.;
  Float_t tpcChi2Cl = -99.;

  nClsTPC = trk->GetTPCNcls();
  if(nClsTPC < 50. || nClsTPC > 160.) return kFALSE;
  tpcChi2Cl = nClsTPC > 0 ? trk->GetTPCchi2() / nClsTPC : -1.;
  if(tpcChi2Cl < 0. || tpcChi2Cl > 3.) return kFALSE;

  Float_t nSigmaEleTPC = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kElectron); // Use only raw tpc pid not corrected with the dq-framework
  if( TMath::Abs(nSigmaEleTPC) > 3. ) return kFALSE;

  return kTRUE;
}

//______________________________________________
void AliDielectronEvtVsTrkHist::CalculateMatchingEfficiency()
{
  const Int_t nBins = fMatchEffTPC->GetNbins();
  printf("NbinsTPC %lld NbinsITS %lld\n", fMatchEffTPC->GetNbins(), fMatchEffITS->GetNbins());
  Double_t matchEff;
  Double_t matchEffErr;
  Int_t matchEffDim = -1;

  Double_t nCountsITS;
  Double_t nCountsTPC;

  const Int_t nDimMatchEff = fSparseObjs[ADEVTH::kSparseMatchEffITSTPC]->GetNdimensions();
  for (Int_t iDim = 0; iDim < nDimMatchEff; iDim++) {
    if(fSparseVars[ADEVTH::kSparseMatchEffITSTPC][iDim] == ADVM::kMatchEffITSTPC) matchEffDim = iDim;
  }

  Int_t *coord = new Int_t[nDimMatchEff];
  Long64_t errBin = -1;
  for (Int_t iBin = 0; iBin < nBins; iBin++) {
    // Get the number of tracks in the bin from the TPC and ITS histograms
    nCountsTPC = fMatchEffTPC->GetBinContent(iBin, coord);
    nCountsITS = fMatchEffITS->GetBinContent(iBin);

    if (nCountsTPC < 0.9) {
      nCountsITS = 0.;
      nCountsTPC = 1.;
    }

    // Calculate the matching efficiency for the given bin
    matchEff = nCountsITS/nCountsTPC;

    // Set the matching efficiency bin in the coordinates
    coord[matchEffDim] = fMatchEffTPC->GetAxis(matchEffDim)->FindBin(matchEff);

    Double_t err1 = fMatchEffTPC->GetBinError(iBin) * nCountsTPC;
    Double_t err2 = fMatchEffITS->GetBinError(iBin) * nCountsITS;
    Double_t b22  = nCountsITS * nCountsTPC;
    if(b22 > 0) matchEffErr = (err1 * err1 + err2 * err2) / (b22 * b22);
    else matchEffErr = 0.;
    matchEffErr = TMath::Sqrt(nCountsITS);

// Now set the correct match efficiency value to the bin coordinates
    fSparseObjs[ADEVTH::kSparseMatchEffITSTPC]->SetBinContent(coord, nCountsITS);
    errBin = fSparseObjs[ADEVTH::kSparseMatchEffITSTPC]->GetBin(coord);

    fSparseObjs[ADEVTH::kSparseMatchEffITSTPC]->SetBinError2(errBin, matchEffErr);
  }
  delete [] coord;
}
