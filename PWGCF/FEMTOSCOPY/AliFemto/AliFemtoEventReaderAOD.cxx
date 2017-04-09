///
/// \file AliFemtoEventReaderAOD.cxx
///

#include "AliFemtoEventReaderAOD.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliPID.h"

#include "AliAODpidUtil.h"
#include "AliAnalysisUtils.h"
#include "AliGenHijingEventHeader.h"

#include "AliExternalTrackParam.h"

#include <algorithm>
#include <cassert>
#include <map>


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventReaderAOD);
  /// \endcond
#endif

#if !(ST_NO_NAMESPACES)
using namespace units;
#endif

using namespace std;

double fV1[3];

//____________________________
//constructor with 0 parameters , look at default settings
AliFemtoEventReaderAOD::AliFemtoEventReaderAOD():
  fNumberofEvent(0),
  fCurEvent(0),
  fEvent(NULL),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fFilterMask(0),
  //  fPWG2AODTracks(0x0),
  fReadMC(0),
  fReadV0(0),
  fReadCascade(0),
  fUsePreCent(0),
  fEstEventMult(kCentrality),
  fAODpidUtil(NULL),
  fAODheader(NULL),
  fAnaUtils(NULL),
  fInputFile(""),
  fTree(NULL),
  fAodFile(NULL),
  fMagFieldSign(1),
  fisEPVZ(kTRUE),
  fpA2013(kFALSE),
  fisPileUp(kFALSE),
  fMVPlp(kFALSE),
  fMinVtxContr(0),
  fMinPlpContribMV(0),
  fMinPlpContribSPD(0),
  fDCAglobalTrack(0),
  fFlatCent(kFALSE),
  fShiftPosition(0.),
  fPrimaryVertexCorrectionTPCPoints(kFALSE),
  f1DcorrectionsPions(0),
  f1DcorrectionsKaons(0),
  f1DcorrectionsProtons(0),
  f1DcorrectionsPionsMinus(0),
  f1DcorrectionsKaonsMinus(0),
  f1DcorrectionsProtonsMinus(0),
  f1DcorrectionsAll(0),
  f1DcorrectionsLambdas(0),
  f1DcorrectionsLambdasMinus(0),
  fIsKaonAnalysis(kFALSE),
  fIsProtonAnalysis(kFALSE),
  fIsPionAnalysis(kFALSE),
  fIsElectronAnalysis(kFALSE)
{
  // default constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fCentRange[0] = 0;
  fCentRange[1] = 1000;
}

AliFemtoEventReaderAOD::AliFemtoEventReaderAOD(const AliFemtoEventReaderAOD &aReader):
  AliFemtoEventReader(),
  fNumberofEvent(aReader.fNumberofEvent),
  fCurEvent(aReader.fCurEvent),
  fEvent(new AliAODEvent()),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(aReader.fFilterBit),
  fFilterMask(aReader.fFilterMask),
  //  fPWG2AODTracks(0x0),
  fReadMC(aReader.fReadMC),
  fReadV0(aReader.fReadV0),
  fReadCascade(aReader.fReadCascade),
  fUsePreCent(aReader.fUsePreCent),
  fEstEventMult(aReader.fEstEventMult),
  fAODpidUtil(aReader.fAODpidUtil),
  fAODheader(aReader.fAODheader),
  fAnaUtils(aReader.fAnaUtils),
  fInputFile(aReader.fInputFile),
  fTree(NULL),
  fAodFile(new TFile(aReader.fAodFile->GetName())),
  fMagFieldSign(aReader.fMagFieldSign),
  fisEPVZ(aReader.fisEPVZ),
  fpA2013(aReader.fpA2013),
  fisPileUp(aReader.fisPileUp),
  fMVPlp(aReader.fMVPlp),
  fMinVtxContr(aReader.fMinVtxContr),
  fMinPlpContribMV(aReader.fMinPlpContribMV),
  fMinPlpContribSPD(aReader.fMinPlpContribSPD),
  fDCAglobalTrack(aReader.fDCAglobalTrack),
  fFlatCent(aReader.fFlatCent),
  fShiftPosition(aReader.fShiftPosition),
  fPrimaryVertexCorrectionTPCPoints(aReader.fPrimaryVertexCorrectionTPCPoints),
  f1DcorrectionsPions(aReader.f1DcorrectionsPions),
  f1DcorrectionsKaons(aReader.f1DcorrectionsKaons),
  f1DcorrectionsProtons(aReader.f1DcorrectionsProtons),
  f1DcorrectionsPionsMinus(aReader.f1DcorrectionsPionsMinus),
  f1DcorrectionsKaonsMinus(aReader.f1DcorrectionsKaonsMinus),
  f1DcorrectionsProtonsMinus(aReader.f1DcorrectionsProtonsMinus),
  f1DcorrectionsAll(aReader.f1DcorrectionsAll),
  f1DcorrectionsLambdas(aReader.f1DcorrectionsLambdas),
  f1DcorrectionsLambdasMinus(aReader.f1DcorrectionsLambdasMinus),
  fIsKaonAnalysis(aReader.fIsKaonAnalysis),
  fIsProtonAnalysis(aReader.fIsProtonAnalysis),
  fIsPionAnalysis(aReader.fIsPionAnalysis),
  fIsElectronAnalysis(aReader.fIsElectronAnalysis)

{
  // copy constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);

  //  fPWG2AODTracks = aReader.fPWG2AODTracks;

  fCentRange[0] = aReader.fCentRange[0];
  fCentRange[1] = aReader.fCentRange[1];
}
//__________________
AliFemtoEventReaderAOD::~AliFemtoEventReaderAOD()
{ // destructor
  delete fTree;
  delete fEvent;
  delete fAodFile;
//   if (fPWG2AODTracks) {
//     fPWG2AODTracks->Delete();
//     delete fPWG2AODTracks;
//   }
}

//__________________
AliFemtoEventReaderAOD &AliFemtoEventReaderAOD::operator=(const AliFemtoEventReaderAOD &aReader)
{ // assignment operator
  if (this == &aReader) {
    return *this;
  }

  fInputFile = aReader.fInputFile;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  if (fTree) delete fTree;
  if (fEvent) delete fEvent;
  fEvent = new AliAODEvent();
  if (fAodFile) delete fAodFile;
  fAodFile = new TFile(aReader.fAodFile->GetName());
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fFilterBit = aReader.fFilterBit;
  fFilterMask = aReader.fFilterMask;
  //  fPWG2AODTracks = aReader.fPWG2AODTracks;
  fAODpidUtil = aReader.fAODpidUtil;
  fAODheader = aReader.fAODheader;
  fAnaUtils = aReader.fAnaUtils;
  fCentRange[0] = aReader.fCentRange[0];
  fCentRange[1] = aReader.fCentRange[1];
  fUsePreCent = aReader.fUsePreCent;
  fEstEventMult = aReader.fEstEventMult;
  fMagFieldSign = aReader.fMagFieldSign;
  fpA2013 = aReader.fpA2013;
  fisPileUp = aReader.fisPileUp;
  fMVPlp = aReader.fMVPlp;
  fMinVtxContr = aReader.fMinVtxContr;
  fMinPlpContribMV = aReader.fMinPlpContribMV;
  fMinPlpContribSPD = aReader.fMinPlpContribSPD;
  fDCAglobalTrack = aReader.fDCAglobalTrack;
  fFlatCent = aReader.fFlatCent;
  fShiftPosition = aReader.fShiftPosition;
  fPrimaryVertexCorrectionTPCPoints = aReader.fPrimaryVertexCorrectionTPCPoints;
  f1DcorrectionsPions = aReader.f1DcorrectionsPions;
  f1DcorrectionsKaons = aReader.f1DcorrectionsKaons;
  f1DcorrectionsProtons = aReader.f1DcorrectionsProtons;
  f1DcorrectionsPionsMinus = aReader.f1DcorrectionsPionsMinus;
  f1DcorrectionsKaonsMinus = aReader.f1DcorrectionsKaonsMinus;
  f1DcorrectionsProtonsMinus = aReader.f1DcorrectionsProtonsMinus;
  f1DcorrectionsAll = aReader.f1DcorrectionsAll;
  f1DcorrectionsLambdas = aReader.f1DcorrectionsLambdas;
  f1DcorrectionsLambdasMinus = aReader.f1DcorrectionsLambdasMinus;
  fIsKaonAnalysis = aReader.fIsKaonAnalysis;
  fIsProtonAnalysis = aReader.fIsProtonAnalysis;
  fIsPionAnalysis = aReader.fIsPionAnalysis;
  fIsElectronAnalysis = aReader.fIsElectronAnalysis;

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderAOD::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderAOD\n";
  return temp;
}

//__________________
void AliFemtoEventReaderAOD::SetInputFile(const char *inputFile)
{
  /// Reads a list of filenames from the file 'inputFile'. Each filename is
  /// checked for a AOD TTree, if present it is added to the fTree chain.
  fInputFile = inputFile;
  ifstream infile(inputFile);

  delete fTree;
  fTree = new TChain("aodTree");

  for (std::string line; std::getline(infile, line);) {
    const char *filename = line.c_str();
    TFile *aodFile = TFile::Open(filename, "READ");
    if (aodFile) {
      TTree *tree = (TTree *) aodFile->Get("aodTree");
      if (tree) {
        fTree->AddFile(filename);
        delete tree;
      }
      aodFile->Close();
    }
    delete aodFile;
  }
}

AliFemtoEvent *AliFemtoEventReaderAOD::ReturnHbtEvent()
{
  /// Reads in the next event from the chain and converts it to an AliFemtoEvent
  AliFemtoEvent *hbtEvent = NULL;

  // We have hit the end of our range -> open the next file
  if (fCurEvent == fNumberofEvent) {
    // We haven't loaded anything yet - open
    if (fNumberofEvent == 0) {
      // cout << "fEvent: " << fEvent << "\n";
      fEvent = new AliAODEvent();
      fEvent->ReadFromTree(fTree);

      fNumberofEvent = fTree->GetEntries();
      // cout << "Number of entries in file " << fNumberofEvent << endl;
      fCurEvent = 0;
    } else { //no more data to read
      fReaderStatus = 1;
      return NULL;
    }
  }

  // cout << "starting to read event " << fCurEvent << endl;
  fTree->GetEvent(fCurEvent);
  // cout << "Read event " << fEvent << " from file " << fTree << endl;

  hbtEvent = CopyAODtoFemtoEvent();
  fCurEvent++;

  return hbtEvent;
}

AliFemtoEvent *AliFemtoEventReaderAOD::CopyAODtoFemtoEvent()
{

  // A function that reads in the AOD event
  // and transfers the neccessary information into
  // the internal AliFemtoEvent

  AliFemtoEvent *tEvent = new AliFemtoEvent();

  // setting global event characteristics
  tEvent->SetRunNumber(fEvent->GetRunNumber());
  tEvent->SetMagneticField(fEvent->GetMagneticField()*kilogauss);//to check if here is ok
  tEvent->SetZDCN1Energy(fEvent->GetZDCN1Energy());
  tEvent->SetZDCP1Energy(fEvent->GetZDCP1Energy());
  tEvent->SetZDCN2Energy(fEvent->GetZDCN2Energy());
  tEvent->SetZDCP2Energy(fEvent->GetZDCP2Energy());
  tEvent->SetZDCEMEnergy(fEvent->GetZDCEMEnergy(0));
  tEvent->SetZDCParticipants(0);
  tEvent->SetTriggerMask(fEvent->GetTriggerMask());
  tEvent->SetTriggerCluster(fEvent->GetTriggerCluster());

  // Attempt to access MC header
  AliAODMCHeader *mcH = NULL;
  TClonesArray *mcP = NULL;

  if (fReadMC) {
    mcH = (AliAODMCHeader *) fEvent->FindListObject(AliAODMCHeader::StdBranchName());
    if (!mcH) {
      cout << "AOD MC information requested, but no header found!" << endl;
    }

    mcP = (TClonesArray *) fEvent->FindListObject(AliAODMCParticle::StdBranchName());
    if (!mcP) {
      cout << "AOD MC information requested, but no particle array found!" << endl;
    }
  }

  AliAODHeader *header = dynamic_cast<AliAODHeader *>(fEvent->GetHeader());
  assert(header && "Not a standard AOD");

  tEvent->SetReactionPlaneAngle(header->GetQTheta(0) / 2.0);
  // Int_t *motherids=0;
  // if (mcP) {
  //   const int motherTabSize = ((AliAODMCParticle *) mcP->At(mcP->GetEntries()-1))->GetLabel();
  //   motherids = new int[motherTabSize+1];
  //   for (int ip=0; ip<motherTabSize+1; ip++) motherids[ip] = 0;

  //   // Read in mother ids
  //   AliAODMCParticle *motherpart;
  //   for (int ip=0; ip<mcP->GetEntries(); ip++) {
  //     motherpart = (AliAODMCParticle *) mcP->At(ip);
  //     if (motherpart->GetDaughter(0) > 0)
  //       motherids[motherpart->GetDaughter(0)] = ip;
  //     if (motherpart->GetDaughter(1) > 0)
  //       motherids[motherpart->GetDaughter(1)] = ip;
  //   }
  // }

  // AliAnalysisUtils
  if (fisPileUp || fpA2013) {
    fAnaUtils = new AliAnalysisUtils();
    if (fMinVtxContr) {
      fAnaUtils->SetMinVtxContr(fMinVtxContr);
    }
    if (fpA2013 && fAnaUtils->IsVertexSelected2013pA(fEvent) == kFALSE) {
      delete fAnaUtils;
      delete tEvent;
      return NULL;  // Vertex rejection for pA analysis.
    }
    fAnaUtils->SetUseMVPlpSelection(fMVPlp);

    if (fMinPlpContribMV) {
      fAnaUtils->SetMinPlpContribMV(fMinPlpContribMV);
    }
    if (fMinPlpContribSPD) {
      fAnaUtils->SetMinPlpContribSPD(fMinPlpContribSPD);
    }
    if (fisPileUp && fAnaUtils->IsPileUpEvent(fEvent)) {
      delete fAnaUtils;
      delete tEvent;
      return NULL;  // Pile-up rejection.
    }
    delete fAnaUtils;
  }

  // Primary Vertex position
  const AliAODVertex *aodvertex = (AliAODVertex *) fEvent->GetPrimaryVertex();
  if (!aodvertex || aodvertex->GetNContributors() < 1) {
    delete tEvent;  // Bad vertex, skip event.
    return NULL;
  }

  aodvertex->GetPosition(fV1);
  AliFmThreeVectorF vertex(fV1[0], fV1[1], fV1[2]);
  tEvent->SetPrimVertPos(vertex);

  //starting to reading tracks
  int nofTracks = 0; //number of reconstructed tracks in event

  nofTracks = fEvent->GetNumberOfTracks();
  //cout<<"NEW: ======> nofTracks "<<nofTracks<<endl;
  AliEventplane *ep = fEvent->GetEventplane();
  if (ep) {
    const float event_plane_angle = (fisEPVZ)
                                  ? ep->GetEventplane("V0", fEvent, 2)
                                  : ep->GetEventplane("Q");
    tEvent->SetEP(ep);
    tEvent->SetReactionPlaneAngle(event_plane_angle);
  }

  AliCentrality *cent = fEvent->GetCentrality();

  if (!fEstEventMult && cent && fUsePreCent) {
    if ((cent->GetCentralityPercentile("V0M") * 10 < fCentRange[0]) ||
        (cent->GetCentralityPercentile("V0M") * 10 > fCentRange[1])) {
      delete tEvent;
      return NULL;
    }
  }

  const Float_t percent = cent->GetCentralityPercentile("V0M");

  // Flatten centrality distribution
  if (percent < 9 && fFlatCent) {
    bool reject_event = RejectEventCentFlat(fEvent->GetMagneticField(), percent);
    if (reject_event) {
      delete tEvent;
      return NULL;
    }
  }

  int realnofTracks = 0; // number of track which we use in a analysis
  int tracksPrim = 0;

  // 'labels' maps a track's id to the track's index in the Event
  // i.e. labels[Event->GetTrack(x)->GetID()] == x
  std::map<int, int> labels;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i = 0; i < nofTracks; i++) {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(fEvent->GetTrack(i));
    assert(aodtrack && "Not a standard AOD");
    if (!aodtrack->TestFilterBit(fFilterBit)) {
      // Skip TPC-only tracks
      if (aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }

  int tNormMult = 0;
  for (int i = 0; i < nofTracks; i++) {

    //  const AliAODTrack *aodtrack=dynamic_cast<AliAODTrack*>(fEvent->GetTrack(i));
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(fEvent->GetTrack(i));
    assert(aodtrack && "Not a standard AOD"); // Getting the AODtrack directly

    if (aodtrack->IsPrimaryCandidate()) tracksPrim++;

    if ((fFilterBit && !aodtrack->TestFilterBit(fFilterBit)) ||
        (fFilterMask && !aodtrack->TestFilterBit(fFilterMask))) {
      continue;
    }

    // Check the sanity of the tracks - reject zero momentum tracks
    if (aodtrack->P() == 0.0) {
      continue;
    }

    // Counting particles to set multiplicity
    if (fEstEventMult == kGlobalCount) {
      //if (aodtrack->IsPrimaryCandidate()) //? instead of kinks?
      if (aodtrack->Chi2perNDF() < 4.0)
        if (0.15 <= aodtrack->Pt() && aodtrack->Pt() < 20)
          if (aodtrack->GetTPCNcls() > 70)
            if (aodtrack->Eta() < 0.8)
              tNormMult++;
    }

    AliFemtoTrack *trackCopy = CopyAODtoFemtoTrack(aodtrack);

    // copying PID information from the correspondent track
    //  const AliAODTrack *aodtrackpid = fEvent->GetTrack(labels[-1-fEvent->GetTrack(i)->GetID()]);


    // For TPC Only tracks we have to copy PID information from corresponding global tracks
    const Int_t pid_track_id = (fFilterBit == (1 << 7) || fFilterMask == 128)
                             ? labels[-1 - fEvent->GetTrack(i)->GetID()]
                             : i;
    AliAODTrack *aodtrackpid = dynamic_cast<AliAODTrack *>(fEvent->GetTrack(pid_track_id));
    assert(aodtrackpid && "Not a standard AOD");

    CopyPIDtoFemtoTrack(aodtrackpid, trackCopy);

    if (mcP) {

      // Fill the hidden information with the simulated data
      Int_t track_label = aodtrack->GetLabel();
      const AliAODMCParticle *tPart = (track_label > -1)
                                    ? (AliAODMCParticle *)mcP->At(track_label)
                                    : NULL;
      AliFemtoModelGlobalHiddenInfo *tInfo = new AliFemtoModelGlobalHiddenInfo();
      double fpx = 0.0, fpy = 0.0, fpz = 0.0, fpt = 0.0;
      if (!tPart) {
        fpx = fV1[0];
        fpy = fV1[1];
        fpz = fV1[2];
        tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);
        tInfo->SetPDGPid(0);
        tInfo->SetTrueMomentum(0.0, 0.0, 0.0);
        tInfo->SetEmissionPoint(0.0, 0.0, 0.0, 0.0);
        tInfo->SetMass(0);
      } else {
        // Check the mother information

        // Using the new way of storing the freeze-out information
        // Final state particle is stored twice on the stack
        // one copy (mother) is stored with original freeze-out information
        //   and is not tracked
        // the other one (daughter) is stored with primary vertex position
        //   and is tracked

        // Freeze-out coordinates
        fpx = tPart->Xv() - fV1[0];
        fpy = tPart->Yv() - fV1[1];
        fpz = tPart->Zv() - fV1[2];
        //    fpt = tPart->T();

        tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);

        fpx *= 1e13;
        fpy *= 1e13;
        fpz *= 1e13;
        //    fpt *= 1e13;

        //      cout << "Looking for mother ids " << endl;

        //if (motherids[TMath::Abs(aodtrack->GetLabel())]>0) {
        if (tPart->GetMother() > -1) { //MC particle has a mother
          //  cout << "Got mother id" << endl;
          //          AliAODMCParticle *mother = GetParticleWithLabel(mcP, motherids[TMath::Abs(aodtrack->GetLabel())]);
          AliAODMCParticle *mother = (AliAODMCParticle *)mcP->At(tPart->GetMother());
          // Check if this is the same particle stored twice on the stack
          if (mother) {
            if ((mother->GetPdgCode() == tPart->GetPdgCode() || (mother->Px() == tPart->Px()))) {
              // It is the same particle
              // Read in the original freeze-out information
              // and convert it from to [fm]

              // EPOS style
              //    fpx = mother->Xv()*1e13*0.197327;
              //    fpy = mother->Yv()*1e13*0.197327;
              //    fpz = mother->Zv()*1e13*0.197327;
              //    fpt = mother->T() *1e13*0.197327*0.5;


              // Therminator style
              fpx = mother->Xv() * 1e13;
              fpy = mother->Yv() * 1e13;
              fpz = mother->Zv() * 1e13;
              //        fpt = mother->T() *1e13*3e10;

            } else { //particle's mother exists and the information about it can be added to hiddeninfo:
              tInfo->SetMotherPdgCode(mother->GetPdgCode());
            }
          }
        }

        //       if (fRotateToEventPlane) {
        //  double tPhi = TMath::ATan2(fpy, fpx);
        //  double tRad = TMath::Hypot(fpx, fpy);

        //  fpx = tRad*TMath::Cos(tPhi - tReactionPlane);
        //  fpy = tRad*TMath::Sin(tPhi - tReactionPlane);
        //       }

        tInfo->SetPDGPid(tPart->GetPdgCode());

        //    if (fRotateToEventPlane) {
        //      double tPhi = TMath::ATan2(tPart->Py(), tPart->Px());
        //      double tRad = TMath::Hypot(tPart->Px(), tPart->Py());

        //      tInfo->SetTrueMomentum(tRad*TMath::Cos(tPhi - tReactionPlane),
        //           tRad*TMath::Sin(tPhi - tReactionPlane),
        //           tPart->Pz());
        //    }
        //       else
        tInfo->SetTrueMomentum(tPart->Px(), tPart->Py(), tPart->Pz());
        Double_t mass2 = (tPart->E() * tPart->E()
                          - tPart->Px() * tPart->Px()
                          - tPart->Py() * tPart->Py()
                          - tPart->Pz() * tPart->Pz());
        if (mass2 > 0.0)
          tInfo->SetMass(TMath::Sqrt(mass2));
        else
          tInfo->SetMass(0.0);

        tInfo->SetEmissionPoint(fpx, fpy, fpz, fpt);

	if(tPart->IsPhysicalPrimary())
	  tInfo->SetOrigin(0);
	else if(tPart->IsSecondaryFromWeakDecay())
	  tInfo->SetOrigin(1);
	else if(tPart->IsSecondaryFromMaterial())
	  tInfo->SetOrigin(2);
	else
	  tInfo->SetOrigin(-1);
		

        // // fillDCA
        // //if (TMath::Abs(impact[0]) > 0.001) {
        // if (tPart->IsPhysicalPrimary()){
        //   tInfo->SetPartOrigin(0);
        //   // trackCopy->SetImpactDprim(impact[0]);
        //   //cout << "Read prim" << endl;
        // }
        // else if (tPart->IsSecondaryFromWeakDecay()) {
        //   tInfo->SetPartOrigin(1);
        //   // trackCopy->SetImpactDweak(impact[0]);
        //   //cout << "Read wea" << endl;
        // }
        // else if (tPart->IsSecondaryFromMaterial()) {
        //   tInfo->SetPartOrigin(2);
        //   // trackCopy->SetImpactDmat(impact[0]);
        //   //cout << "Read mat" << endl;
        // }
        // //}
        // //  end fillDCA


      }
      trackCopy->SetHiddenInfo(tInfo);
    }

    //AliExternalTrackParam *param = new AliExternalTrackParam(*aodtrack->GetInnerParam());
    trackCopy->SetInnerMomentum(aodtrack->GetTPCmomentum());

    //Special MC analysis for pi,K,p,e slected by PDG code -->
    if(fIsKaonAnalysis || fIsProtonAnalysis || fIsPionAnalysis || fIsElectronAnalysis) {
      Int_t pdg = ((AliFemtoModelHiddenInfo*)trackCopy->GetHiddenInfo())->GetPDGPid();
      Double_t ptrue = ((AliFemtoModelHiddenInfo*)trackCopy->GetHiddenInfo())->GetTrueMomentum()->Mag();
  
      //if(fIsKaonAnalysis&&TMath::Abs(pdg) == 321)cout<<"AOD REader pdg cod "<<pdg<<" ptrue "<<ptrue<< endl;
      
      Bool_t trackAccept = true;

      if (fIsKaonAnalysis == true && TMath::Abs(pdg) != 321) 	trackAccept = false;
      if (fIsProtonAnalysis == true && TMath::Abs(pdg) != 2212) trackAccept = false;
      if (fIsPionAnalysis == true && TMath::Abs(pdg) != 211)    trackAccept = false;
      if (fIsElectronAnalysis == true && TMath::Abs(pdg) != 11) trackAccept = false;

      if (trackAccept == true && ptrue > 0) {    
	tEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
	realnofTracks++;//real number of tracks
      }
      else {
	// cout<<"bad track : AOD REader pdg cod"<<pdg<<" ptrue "<<ptrue<<endl;
	delete  trackCopy;
      }
      //Special MC analysis for pi,K,p,e slected by PDG code <--
    }
    else {
      tEvent->TrackCollection()->push_back(trackCopy);  // Adding track to analysis
      realnofTracks++; // Real number of tracks
    }
  }
  //cout<<"======================> realnofTracks"<<realnofTracks<<endl;
  tEvent->SetNumberOfTracks(realnofTracks); // Setting number of track which we read in event

  if (cent) {
    tEvent->SetCentralityV0(cent->GetCentralityPercentile("V0M"));
    tEvent->SetCentralityV0A(cent->GetCentralityPercentile("V0A"));
    tEvent->SetCentralityV0C(cent->GetCentralityPercentile("V0C"));
    tEvent->SetCentralityZNA(cent->GetCentralityPercentile("ZNA"));
    tEvent->SetCentralityZNC(cent->GetCentralityPercentile("ZNC"));
    tEvent->SetCentralityCL1(cent->GetCentralityPercentile("CL1"));
    tEvent->SetCentralityCL0(cent->GetCentralityPercentile("CL0"));
    tEvent->SetCentralityTKL(cent->GetCentralityPercentile("TKL"));
    tEvent->SetCentralityFMD(cent->GetCentralityPercentile("FMD"));
    tEvent->SetCentralityFMD(cent->GetCentralityPercentile("NPA"));
    //    tEvent->SetCentralitySPD1(cent->GetCentralityPercentile("CL1"));
    tEvent->SetCentralityTrk(cent->GetCentralityPercentile("TRK"));
    tEvent->SetCentralityCND(cent->GetCentralityPercentile("CND"));
  }


  Int_t norm_mult = tracksPrim;

  if (cent) {
    switch (fEstEventMult) {
      case kCentrality:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("V0M"));
        break;
      case kCentralityV0A:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("V0A"));
        break;
      case kCentralityV0C:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("V0C"));
        break;
      case kCentralityZNA:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("ZNA"));
        break;
      case kCentralityZNC:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("ZNC"));
        break;
      case kCentralityCL1:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("CL1"));
        break;
      case kCentralityCL0:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("CL0"));
        break;
      case kCentralityTRK:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("TRK"));
        break;
      case kCentralityTKL:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("TKL"));
        break;
      case kCentralityCND:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("CND"));
        break;
      case kCentralityNPA:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("NPA"));
        break;
      case kCentralityFMD:
        norm_mult = lrint(10 * cent->GetCentralityPercentile("FMD"));
        break;
    }
  } else {
    switch (fEstEventMult) {
      case kGlobalCount:
        norm_mult = tNormMult; // Particles counted in the loop, trying to reproduce GetReferenceMultiplicity. If better (default) method appears it should be changed
        break;
      case kReference:
        norm_mult = fAODheader->GetRefMultiplicity();
        break;
      case kTPCOnlyRef:
        norm_mult = fAODheader->GetTPConlyRefMultiplicity();
        break;
      case kVZERO:
        Float_t multV0 = 0.0;
        for (Int_t i = 0; i < 64; i++)
          multV0 += fEvent->GetVZEROData()->GetMultiplicity(i);
        norm_mult = lrint(multV0);
    }
  }

  tEvent->SetNormalizedMult(norm_mult);

  if (fReadV0) {
    int count_pass = 0;
    for (Int_t i = 0; i < fEvent->GetNumberOfV0s(); i++) {
      AliAODv0 *aodv0 = fEvent->GetV0(i);
      if (!aodv0) continue;
      if (aodv0->GetNDaughters() > 2) continue;
      if (aodv0->GetNProngs() > 2) continue;
      if (aodv0->GetCharge() != 0) continue;
      if (aodv0->ChargeProng(0) == aodv0->ChargeProng(1)) continue;
      if (aodv0->CosPointingAngle(fV1) < 0.98) continue;

      AliAODTrack *daughterTrackPos = (AliAODTrack *)aodv0->GetDaughter(0),   // getting positive daughter track
                  *daughterTrackNeg = (AliAODTrack *)aodv0->GetDaughter(1);   // getting negative daughter track
      if (daughterTrackPos == NULL || daughterTrackNeg == NULL) continue;     // daughter tracks must exist
      if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) continue; // and have different charge

      // ensure that pos and neg are pointing to the correct children
      if (daughterTrackPos->Charge() < 0 && daughterTrackNeg->Charge() > 0) {
        // std::swap(daughterTrackPos, daughterTrackNeg); // can we use this?
        AliAODTrack *tmp = daughterTrackPos;
        daughterTrackPos = daughterTrackNeg;
        daughterTrackNeg = tmp;
      }

      AliFemtoV0 *trackCopyV0 = CopyAODtoFemtoV0(aodv0);
      if (mcP) {
        daughterTrackPos->SetAODEvent(fEvent);
        daughterTrackNeg->SetAODEvent(fEvent);
        if (daughterTrackPos->GetLabel() > 0 && daughterTrackNeg->GetLabel() > 0) {
          // get the MC data for the two daughter particles
          const AliAODMCParticle *mcParticlePos = static_cast<AliAODMCParticle*>(mcP->At(daughterTrackPos->GetLabel())),
                                 *mcParticleNeg = static_cast<AliAODMCParticle*>(mcP->At(daughterTrackNeg->GetLabel()));

          // They daughter info MUST exist for both
          if ((mcParticlePos != NULL) && (mcParticleNeg != NULL)) {
            // Get the mother ID of the two daughters
            const int motherOfPosID = mcParticlePos->GetMother(),
                      motherOfNegID = mcParticleNeg->GetMother();

            // If both daughter tracks refer to the same mother, we can continue
            if ((motherOfPosID > -1) && (motherOfPosID == motherOfNegID)) {
              // Create the MC data store
              AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();

              // Our V0 particle
              const AliAODMCParticle *v0 = static_cast<AliAODMCParticle*>(mcP->At(motherOfPosID));

              if (v0 == NULL) {
                tInfo->SetPDGPid(0);
                tInfo->SetTrueMomentum(0.0, 0.0, 0.0);
                tInfo->SetMass(0);
              } else {
                //-----v0 particle-----
                const int v0MotherId = v0->GetMother();
                if (v0MotherId > -1) { //V0 particle has a mother
                  AliAODMCParticle *motherOfV0 = static_cast<AliAODMCParticle*>(mcP->At(v0MotherId));
                  tInfo->SetMotherPdgCode(motherOfV0->GetPdgCode());
                }
                tInfo->SetPDGPid(v0->GetPdgCode());
                tInfo->SetMass(v0->GetCalcMass());
                tInfo->SetTrueMomentum(v0->Px(), v0->Py(), v0->Pz());
                tInfo->SetEmissionPoint(v0->Xv(), v0->Yv(), v0->Zv(), v0->T());

		if(v0->IsPhysicalPrimary())
		  tInfo->SetOrigin(0);
		else if(v0->IsSecondaryFromWeakDecay())
		  tInfo->SetOrigin(1);
		else if(v0->IsSecondaryFromMaterial())
		  tInfo->SetOrigin(2);
		else
		  tInfo->SetOrigin(-1);
		
                //-----Positive daughter of v0-----
                tInfo->SetPDGPidPos(mcParticlePos->GetPdgCode());
                tInfo->SetMassPos(mcParticlePos->GetCalcMass());
                tInfo->SetTrueMomentumPos(mcParticlePos->Px(), mcParticlePos->Py(), mcParticlePos->Pz());
                tInfo->SetEmissionPointPos(mcParticlePos->Xv(), mcParticlePos->Yv(), mcParticlePos->Zv(), mcParticlePos->T());

                //-----Negative daughter of v0-----
                tInfo->SetPDGPidNeg(mcParticleNeg->GetPdgCode());
                tInfo->SetMassNeg(mcParticleNeg->GetCalcMass());
                tInfo->SetTrueMomentumNeg(mcParticleNeg->Px(), mcParticleNeg->Py(), mcParticleNeg->Pz());
                tInfo->SetEmissionPointNeg(mcParticleNeg->Xv(), mcParticleNeg->Yv(), mcParticleNeg->Zv(), mcParticleNeg->T());
              }
              trackCopyV0->SetHiddenInfo(tInfo);
            }
          }
        }
      }
      tEvent->V0Collection()->push_back(trackCopyV0);
      count_pass++;
    }
  }


  if (fReadCascade) {
    int count_pass = 0;
    for (Int_t i = 0; i < fEvent->GetNumberOfCascades(); i++) {
      AliAODcascade *aodxi = fEvent->GetCascade(i);
      if (!aodxi) continue;
      //if (aodxi->GetNDaughters() > 2) continue;
      //if (aodxi->GetNProngs() > 2) continue;
      //if (aodxi->GetCharge() != 0) continue;
      if (aodxi->ChargeProng(0) == aodxi->ChargeProng(1)) continue;
      if (aodxi->CosPointingAngle(fV1) < 0.9) continue;
      if (aodxi->CosPointingAngleXi(fV1[0],fV1[1],fV1[2]) < 0.98) continue;

      const AliAODTrack *daughterTrackPos = (AliAODTrack *)aodxi->GetDaughter(0), // getting positive daughter track
                        *daughterTrackNeg = (AliAODTrack *)aodxi->GetDaughter(1); // getting negative daughter track
      if (daughterTrackPos == NULL || daughterTrackNeg == NULL) continue;         // daughter tracks must exist
      if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) continue;     // and have different charge

      AliFemtoXi *trackCopyXi = CopyAODtoFemtoXi(aodxi);
      tEvent->XiCollection()->push_back(trackCopyXi);
      count_pass++;
    }
  }

  return tEvent;
}

AliFemtoTrack *AliFemtoEventReaderAOD::CopyAODtoFemtoTrack(AliAODTrack *tAodTrack)
{
  // Copy the track information from the AOD into the internal AliFemtoTrack
  // If it exists, use the additional information from the PWG2 AOD
  AliFemtoTrack *tFemtoTrack = new AliFemtoTrack();

  // Primary Vertex position

  fEvent->GetPrimaryVertex()->GetPosition(fV1);
  //  fEvent->GetPrimaryVertex()->GetXYZ(fV1);
  tFemtoTrack->SetPrimaryVertex(fV1);

  tFemtoTrack->SetCharge(tAodTrack->Charge());

  double pxyz[3];
  tAodTrack->PxPyPz(pxyz);//reading noconstrained momentum

  AliFemtoThreeVector v(pxyz[0], pxyz[1], pxyz[2]);
  tFemtoTrack->SetP(v);//setting momentum
  tFemtoTrack->SetPt(sqrt(pxyz[0]*pxyz[0] + pxyz[1]*pxyz[1]));
  const AliFmThreeVectorD kOrigin(fV1[0], fV1[1], fV1[2]);
  //setting track helix
  const AliFmThreeVectorD ktP(pxyz[0], pxyz[1], pxyz[2]);
  AliFmPhysicalHelixD helix(ktP, kOrigin, (double)(fEvent->GetMagneticField())*kilogauss, (double)(tFemtoTrack->Charge()));
  tFemtoTrack->SetHelix(helix);


  // Flags
  tFemtoTrack->SetTrackId(tAodTrack->GetID());
  tFemtoTrack->SetFlags(tAodTrack->GetFlags());
  tFemtoTrack->SetLabel(tAodTrack->GetLabel());

  // Track quality information
  float covmat[6];
  tAodTrack->GetCovMatrix(covmat);

  if (fDCAglobalTrack==0) {

    tFemtoTrack->SetImpactD(tAodTrack->DCA());
    tFemtoTrack->SetImpactZ(tAodTrack->ZAtDCA());

    tFemtoTrack->SetXatDCA(tAodTrack->XAtDCA());
    tFemtoTrack->SetYatDCA(tAodTrack->YAtDCA());
    tFemtoTrack->SetZatDCA(tAodTrack->ZAtDCA());
  }

  tFemtoTrack->SetCdd(covmat[0]);
  tFemtoTrack->SetCdz(covmat[1]);
  tFemtoTrack->SetCzz(covmat[2]);
  tFemtoTrack->SetITSchi2(tAodTrack->Chi2perNDF());
  tFemtoTrack->SetITSncls(tAodTrack->GetITSNcls());
  tFemtoTrack->SetTPCchi2(tAodTrack->Chi2perNDF());
  tFemtoTrack->SetTPCncls(tAodTrack->GetTPCNcls());
  tFemtoTrack->SetTPCnclsF(tAodTrack->GetTPCNcls());
  tFemtoTrack->SetTPCsignalN(1);
  tFemtoTrack->SetTPCsignalS(1);
  tFemtoTrack->SetTPCsignal(tAodTrack->GetTPCsignal());
  tFemtoTrack->SetTPCClusterMap(tAodTrack->GetTPCClusterMap());
  tFemtoTrack->SetTPCSharedMap(tAodTrack->GetTPCSharedMap());

  float globalPositionsAtRadii[9][3];
  float bfield = 5 * fMagFieldSign;

  GetGlobalPositionAtGlobalRadiiThroughTPC(tAodTrack, bfield, globalPositionsAtRadii);
  double tpcEntrance[3] = {globalPositionsAtRadii[0][0], globalPositionsAtRadii[0][1], globalPositionsAtRadii[0][2]};
  double **tpcPositions;
  tpcPositions = new double*[9];

  for (int i = 0; i < 9; i++) {
    tpcPositions[i] = new double[3];
  }

  double tpcExit[3] = {globalPositionsAtRadii[8][0], globalPositionsAtRadii[8][1], globalPositionsAtRadii[8][2]};
  for (int i = 0; i < 9; i++) {
    tpcPositions[i][0] = globalPositionsAtRadii[i][0];
    tpcPositions[i][1] = globalPositionsAtRadii[i][1];
    tpcPositions[i][2] = globalPositionsAtRadii[i][2];
  }

  if (fPrimaryVertexCorrectionTPCPoints) {
    tpcEntrance[0] -= fV1[0];
    tpcEntrance[1] -= fV1[1];
    tpcEntrance[2] -= fV1[2];

    tpcExit[0] -= fV1[0];
    tpcExit[1] -= fV1[1];
    tpcExit[2] -= fV1[2];

    for (int i = 0; i < 9; i++) {
      tpcPositions[i][0] -= fV1[0];
      tpcPositions[i][1] -= fV1[1];
      tpcPositions[i][2] -= fV1[2];
    }
  }

  tFemtoTrack->SetNominalTPCEntrancePoint(tpcEntrance);
  tFemtoTrack->SetNominalTPCPoints(tpcPositions);
  tFemtoTrack->SetNominalTPCExitPoint(tpcExit);

  if (fShiftPosition > 0.) {
    Float_t posShifted[3];
    SetShiftedPositions(tAodTrack, bfield, posShifted, fShiftPosition);
    tFemtoTrack->SetNominalTPCPointShifted(posShifted);
  }

  for (int i = 0; i < 9; i++) {
    delete [] tpcPositions[i];
  }
  delete [] tpcPositions;


  int indexes[3];
  for (int ik = 0; ik < 3; ik++) {
    indexes[ik] = 0;
  }
  tFemtoTrack->SetKinkIndexes(indexes);


  for (int ii = 0; ii < 6; ii++) {
    tFemtoTrack->SetITSHitOnLayer(ii, tAodTrack->HasPointOnITSLayer(ii));
  }

  //Corrections
  if(f1DcorrectionsPions){
    tFemtoTrack->SetCorrectionPion(f1DcorrectionsPions->GetBinContent(f1DcorrectionsPions->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionPion(1.0);

  if(f1DcorrectionsKaons){
    tFemtoTrack->SetCorrectionKaon(f1DcorrectionsKaons->GetBinContent(f1DcorrectionsKaons->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionKaon(1.0);

  if(f1DcorrectionsProtons){
    tFemtoTrack->SetCorrectionProton(f1DcorrectionsProtons->GetBinContent(f1DcorrectionsProtons->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionProton(1.0);

  if(f1DcorrectionsPionsMinus){
    tFemtoTrack->SetCorrectionPionMinus(f1DcorrectionsPionsMinus->GetBinContent(f1DcorrectionsPionsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionPionMinus(1.0);

  if(f1DcorrectionsKaonsMinus){
    tFemtoTrack->SetCorrectionKaonMinus(f1DcorrectionsKaonsMinus->GetBinContent(f1DcorrectionsKaonsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionKaonMinus(1.0);

  if(f1DcorrectionsProtonsMinus){
    tFemtoTrack->SetCorrectionProtonMinus(f1DcorrectionsProtonsMinus->GetBinContent(f1DcorrectionsProtonsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionProtonMinus(1.0);
  

  if(f1DcorrectionsAll){
    tFemtoTrack->SetCorrectionAll(f1DcorrectionsAll->GetBinContent(f1DcorrectionsAll->FindFixBin(tAodTrack->Pt())));
  }
  else tFemtoTrack->SetCorrectionAll(1.0);


  return tFemtoTrack;
}

AliFemtoV0 *AliFemtoEventReaderAOD::CopyAODtoFemtoV0(AliAODv0 *tAODv0)
{
  AliFemtoV0 *tFemtoV0 = new AliFemtoV0();

  tFemtoV0->SetdecayLengthV0(tAODv0->DecayLength(fV1));
  tFemtoV0->SetdecayVertexV0X(tAODv0->DecayVertexV0X());
  tFemtoV0->SetdecayVertexV0Y(tAODv0->DecayVertexV0Y());
  tFemtoV0->SetdecayVertexV0Z(tAODv0->DecayVertexV0Z());
  AliFemtoThreeVector decayvertex(tAODv0->DecayVertexV0X(), tAODv0->DecayVertexV0Y(), tAODv0->DecayVertexV0Z());
  tFemtoV0->SetdecayVertexV0(decayvertex);
  tFemtoV0->SetdcaV0Daughters(tAODv0->DcaV0Daughters());
  tFemtoV0->SetdcaV0ToPrimVertex(tAODv0->DcaV0ToPrimVertex());
  tFemtoV0->SetdcaPosToPrimVertex(tAODv0->DcaPosToPrimVertex());
  tFemtoV0->SetdcaNegToPrimVertex(tAODv0->DcaNegToPrimVertex());
  tFemtoV0->SetmomPosX(tAODv0->MomPosX());
  tFemtoV0->SetmomPosY(tAODv0->MomPosY());
  tFemtoV0->SetmomPosZ(tAODv0->MomPosZ());
  AliFemtoThreeVector mompos(tAODv0->MomPosX(), tAODv0->MomPosY(), tAODv0->MomPosZ());
  tFemtoV0->SetmomPos(mompos);
  tFemtoV0->SetmomNegX(tAODv0->MomNegX());
  tFemtoV0->SetmomNegY(tAODv0->MomNegY());
  tFemtoV0->SetmomNegZ(tAODv0->MomNegZ());
  AliFemtoThreeVector momneg(tAODv0->MomNegX(), tAODv0->MomNegY(), tAODv0->MomNegZ());
  tFemtoV0->SetmomNeg(momneg);
  tFemtoV0->SetradiusV0(tAODv0->RadiusV0());
  tFemtoV0->SetprimaryVertex(fV1);

  //jest cos takiego w AliFemtoV0.h czego nie ma w AliAODv0.h
  //void SettpcHitsPos(const int& i);
  //void SettpcHitsNeg(const int& i);

  //void SetTrackTopologyMapPos(unsigned int word, const unsigned long& m);
  //void SetTrackTopologyMapNeg(unsigned int word, const unsigned long& m);


  tFemtoV0->SetmomV0X(tAODv0->MomV0X());
  tFemtoV0->SetmomV0Y(tAODv0->MomV0Y());
  tFemtoV0->SetmomV0Z(tAODv0->MomV0Z());
  AliFemtoThreeVector momv0(tAODv0->MomV0X(), tAODv0->MomV0Y(), tAODv0->MomV0Z());
  tFemtoV0->SetmomV0(momv0);
  tFemtoV0->SetalphaV0(tAODv0->AlphaV0());
  tFemtoV0->SetptArmV0(tAODv0->PtArmV0());
  tFemtoV0->SeteLambda(tAODv0->ELambda());
  tFemtoV0->SeteK0Short(tAODv0->EK0Short());
  tFemtoV0->SetePosProton(tAODv0->EPosProton());
  tFemtoV0->SeteNegProton(tAODv0->ENegProton());
  tFemtoV0->SetmassLambda(tAODv0->MassLambda());
  tFemtoV0->SetmassAntiLambda(tAODv0->MassAntiLambda());
  tFemtoV0->SetmassK0Short(tAODv0->MassK0Short());
  tFemtoV0->SetrapLambda(tAODv0->RapLambda());
  tFemtoV0->SetrapK0Short(tAODv0->RapK0Short());

  //void SetcTauLambda( float x);
  //void SetcTauK0Short( float x);

  //tFemtoV0->SetptV0(::sqrt(tAODv0->Pt2V0())); //!
  tFemtoV0->SetptV0(tAODv0->Pt());

  if(f1DcorrectionsLambdas){
    tFemtoV0->SetCorrectionLambdas(f1DcorrectionsLambdas->GetBinContent(f1DcorrectionsLambdas->FindFixBin(tAODv0->Pt())));
  }
  else {
    tFemtoV0->SetCorrectionLambdas(1.0);
  }

    if(f1DcorrectionsLambdasMinus){
    tFemtoV0->SetCorrectionLambdasMinus(f1DcorrectionsLambdasMinus->GetBinContent(f1DcorrectionsLambdasMinus->FindFixBin(tAODv0->Pt())));
  }
  else {
    tFemtoV0->SetCorrectionLambdasMinus(1.0);
  }


  tFemtoV0->SetptotV0(::sqrt(tAODv0->Ptot2V0()));
  //tFemtoV0->SetptPos(::sqrt(tAODv0->MomPosX()*tAODv0->MomPosX()+tAODv0->MomPosY()*tAODv0->MomPosY()));
  //tFemtoV0->SetptotPos(::sqrt(tAODv0->Ptot2Pos()));
  //tFemtoV0->SetptNeg(::sqrt(tAODv0->MomNegX()*tAODv0->MomNegX()+tAODv0->MomNegY()*tAODv0->MomNegY()));
  //tFemtoV0->SetptotNeg(::sqrt(tAODv0->Ptot2Neg()));

  tFemtoV0->SetidNeg(tAODv0->GetNegID());
  //cout<<"tAODv0->GetNegID(): "<<tAODv0->GetNegID()<<endl;
  //cout<<"tFemtoV0->IdNeg(): "<<tFemtoV0->IdNeg()<<endl;
  tFemtoV0->SetidPos(tAODv0->GetPosID());

  tFemtoV0->SetEtaV0(tAODv0->Eta());
  tFemtoV0->SetPhiV0(tAODv0->Phi());
  tFemtoV0->SetCosPointingAngle(tAODv0->CosPointingAngle(fV1));
  //tFemtoV0->SetYV0(tAODv0->Y());


  //void SetdedxNeg(float x);
  //void SeterrdedxNeg(float x);//Gael 04Fev2002
  //void SetlendedxNeg(float x);//Gael 04Fev2002
  //void SetdedxPos(float x);
  //void SeterrdedxPos(float x);//Gael 04Fev2002
  //void SetlendedxPos(float x);//Gael 04Fev2002

  //tFemtoV0->SetEtaPos(tAODv0->PseudoRapPos());
  //tFemtoV0->SetEtaNeg(tAODv0->PseudoRapNeg());

  AliAODTrack *trackpos = (AliAODTrack *)tAODv0->GetDaughter(0);
  AliAODTrack *trackneg = (AliAODTrack *)tAODv0->GetDaughter(1);

  // ensure that trackpos and trackneg are pointing to the correct children
  // This confusion seems to arise when fOnFlyStatusV0 = true
  if(trackpos->Charge() < 0 && trackneg->Charge() > 0)
  {
    AliAODTrack *tmp = trackpos;
    trackpos = trackneg;
    trackneg = tmp;
  }

  if (trackpos && trackneg) {
    tFemtoV0->SetEtaPos(trackpos->Eta());
    tFemtoV0->SetEtaNeg(trackneg->Eta());
    tFemtoV0->SetptotPos(tAODv0->PProng(0));
    tFemtoV0->SetptotNeg(tAODv0->PProng(1));
    tFemtoV0->SetptPos(trackpos->Pt());
    tFemtoV0->SetptNeg(trackneg->Pt());

    //tFemtoV0->SetEtaPos(trackpos->Eta()); //tAODv0->PseudoRapPos()
    //tFemtoV0->SetEtaNeg(trackneg->Eta()); //tAODv0->PseudoRapNeg()
    tFemtoV0->SetTPCNclsPos(trackpos->GetTPCNcls());
    tFemtoV0->SetTPCNclsNeg(trackneg->GetTPCNcls());
    tFemtoV0->SetTPCclustersPos(trackpos->GetTPCClusterMap());
    tFemtoV0->SetTPCclustersNeg(trackneg->GetTPCClusterMap());
    tFemtoV0->SetTPCsharingPos(trackpos->GetTPCSharedMap());
    tFemtoV0->SetTPCsharingNeg(trackneg->GetTPCSharedMap());
    tFemtoV0->SetNdofPos(trackpos->Chi2perNDF());
    tFemtoV0->SetNdofNeg(trackneg->Chi2perNDF());
    tFemtoV0->SetStatusPos(trackpos->GetStatus());
    tFemtoV0->SetStatusNeg(trackneg->GetStatus());

    tFemtoV0->SetPosNSigmaTPCK(fAODpidUtil->NumberOfSigmasTPC(trackpos, AliPID::kKaon));
    tFemtoV0->SetNegNSigmaTPCK(fAODpidUtil->NumberOfSigmasTPC(trackneg, AliPID::kKaon));
    tFemtoV0->SetPosNSigmaTPCP(fAODpidUtil->NumberOfSigmasTPC(trackpos, AliPID::kProton));
    tFemtoV0->SetNegNSigmaTPCP(fAODpidUtil->NumberOfSigmasTPC(trackneg, AliPID::kProton));
    tFemtoV0->SetPosNSigmaTPCPi(fAODpidUtil->NumberOfSigmasTPC(trackpos, AliPID::kPion));
    tFemtoV0->SetNegNSigmaTPCPi(fAODpidUtil->NumberOfSigmasTPC(trackneg, AliPID::kPion));


    float bfield = 5 * fMagFieldSign;
    float globalPositionsAtRadiiPos[9][3];
    GetGlobalPositionAtGlobalRadiiThroughTPC(trackpos, bfield, globalPositionsAtRadiiPos);
    double tpcEntrancePos[3] = {globalPositionsAtRadiiPos[0][0], globalPositionsAtRadiiPos[0][1], globalPositionsAtRadiiPos[0][2]};
    double tpcExitPos[3] = {globalPositionsAtRadiiPos[8][0], globalPositionsAtRadiiPos[8][1], globalPositionsAtRadiiPos[8][2]};

    float globalPositionsAtRadiiNeg[9][3];
    GetGlobalPositionAtGlobalRadiiThroughTPC(trackneg, bfield, globalPositionsAtRadiiNeg);
    double tpcEntranceNeg[3] = {globalPositionsAtRadiiNeg[0][0], globalPositionsAtRadiiNeg[0][1], globalPositionsAtRadiiNeg[0][2]};
    double tpcExitNeg[3] = {globalPositionsAtRadiiNeg[8][0], globalPositionsAtRadiiNeg[8][1], globalPositionsAtRadiiNeg[8][2]};

    if (fPrimaryVertexCorrectionTPCPoints) {
      tpcEntrancePos[0] -= fV1[0];
      tpcEntrancePos[1] -= fV1[1];
      tpcEntrancePos[2] -= fV1[2];

      tpcExitPos[0] -= fV1[0];
      tpcExitPos[1] -= fV1[1];
      tpcExitPos[2] -= fV1[2];

      tpcEntranceNeg[0] -= fV1[0];
      tpcEntranceNeg[1] -= fV1[1];
      tpcEntranceNeg[2] -= fV1[2];

      tpcExitNeg[0] -= fV1[0];
      tpcExitNeg[1] -= fV1[1];
      tpcExitNeg[2] -= fV1[2];
    }

    AliFemtoThreeVector tmpVec;
    tmpVec.SetX(tpcEntrancePos[0]);
    tmpVec.SetY(tpcEntrancePos[1]);
    tmpVec.SetZ(tpcEntrancePos[2]);
    tFemtoV0->SetNominalTpcEntrancePointPos(tmpVec);

    tmpVec.SetX(tpcExitPos[0]);
    tmpVec.SetY(tpcExitPos[1]);
    tmpVec.SetZ(tpcExitPos[2]);
    tFemtoV0->SetNominalTpcExitPointPos(tmpVec);

    tmpVec.SetX(tpcEntranceNeg[0]);
    tmpVec.SetY(tpcEntranceNeg[1]);
    tmpVec.SetZ(tpcEntranceNeg[2]);
    tFemtoV0->SetNominalTpcEntrancePointNeg(tmpVec);

    tmpVec.SetX(tpcExitNeg[0]);
    tmpVec.SetY(tpcExitNeg[1]);
    tmpVec.SetZ(tpcExitNeg[2]);
    tFemtoV0->SetNominalTpcExitPointNeg(tmpVec);


    AliFemtoThreeVector vecTpcPos[9];
    AliFemtoThreeVector vecTpcNeg[9];
    for (int i = 0; i < 9; i++) {
      vecTpcPos[i].SetX(globalPositionsAtRadiiPos[i][0]);
      vecTpcPos[i].SetY(globalPositionsAtRadiiPos[i][1]);
      vecTpcPos[i].SetZ(globalPositionsAtRadiiPos[i][2]);
      vecTpcNeg[i].SetX(globalPositionsAtRadiiNeg[i][0]);
      vecTpcNeg[i].SetY(globalPositionsAtRadiiNeg[i][1]);
      vecTpcNeg[i].SetZ(globalPositionsAtRadiiNeg[i][2]);
    }

    if (fPrimaryVertexCorrectionTPCPoints) {
      AliFemtoThreeVector tmpVertexVec;
      tmpVertexVec.SetX(fV1[0]);
      tmpVertexVec.SetY(fV1[1]);
      tmpVertexVec.SetZ(fV1[2]);

      for (int i = 0; i < 9; i ++) {
        vecTpcPos[i] -= tmpVertexVec;
        vecTpcNeg[i] -= tmpVertexVec;
      }
    }

    tFemtoV0->SetNominalTpcPointPos(vecTpcPos);
    tFemtoV0->SetNominalTpcPointNeg(vecTpcNeg);

    if (fShiftPosition > 0.) {
      Float_t posShiftedPos[3];
      Float_t posShiftedNeg[3];
      SetShiftedPositions(trackpos, bfield, posShiftedPos, fShiftPosition);
      SetShiftedPositions(trackneg, bfield, posShiftedNeg, fShiftPosition);
      AliFemtoThreeVector tmpVecPos;
      AliFemtoThreeVector tmpVecNeg;
      tmpVecPos.SetX(posShiftedPos[0]);
      tmpVecPos.SetY(posShiftedPos[1]);
      tmpVecPos.SetZ(posShiftedPos[2]);
      tmpVecNeg.SetX(posShiftedNeg[0]);
      tmpVecNeg.SetY(posShiftedNeg[1]);
      tmpVecNeg.SetZ(posShiftedNeg[2]);
      tFemtoV0->SetNominalTpcPointPosShifted(tmpVecPos);
      tFemtoV0->SetNominalTpcPointNegShifted(tmpVecNeg);
    }

    tFemtoV0->SetTPCMomentumPos(trackpos->GetTPCmomentum());
    tFemtoV0->SetTPCMomentumNeg(trackneg->GetTPCmomentum());

    tFemtoV0->SetdedxPos(trackpos->GetTPCsignal());
    tFemtoV0->SetdedxNeg(trackneg->GetTPCsignal());


    Float_t probMisPos = 1.0;
    Float_t probMisNeg = 1.0;

    if (((tFemtoV0->StatusPos() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoV0->StatusPos() & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      // if (tFemtoV0->StatusPos() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000
      probMisPos = fAODpidUtil->GetTOFMismatchProbability(trackpos);
    }
    if (((tFemtoV0->StatusNeg() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoV0->StatusNeg() & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      // if (tFemtoV0->StatusNeg() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000
      probMisNeg = fAODpidUtil->GetTOFMismatchProbability(trackneg);
    }

    // if(// (tFemtoV0->StatusPos()& AliESDtrack::kTOFpid)==0 ||
    //    (tFemtoV0->StatusPos()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTOFout)==0 || probMisPos > 0.01)

    if (!(((tFemtoV0->StatusPos() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoV0->StatusPos() & AliVTrack::kTIME) == AliVTrack::kTIME)) || probMisPos > 0.01) {
      // if(// (tFemtoV0->StatusNeg()&AliESDtrack::kTOFpid)==0 ||
      //    (tFemtoV0->StatusNeg()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTOFout)==0 || probMisNeg > 0.01)
      if (!(((tFemtoV0->StatusNeg() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoV0->StatusNeg() & AliVTrack::kTIME) == AliVTrack::kTIME)) || probMisNeg > 0.01) {
        tFemtoV0->SetPosNSigmaTOFK(-1000);
        tFemtoV0->SetNegNSigmaTOFK(-1000);
        tFemtoV0->SetPosNSigmaTOFP(-1000);
        tFemtoV0->SetNegNSigmaTOFP(-1000);
        tFemtoV0->SetPosNSigmaTOFPi(-1000);
        tFemtoV0->SetNegNSigmaTOFPi(-1000);

        tFemtoV0->SetTOFProtonTimePos(-1000);
        tFemtoV0->SetTOFPionTimePos(-1000);
        tFemtoV0->SetTOFKaonTimePos(-1000);
        tFemtoV0->SetTOFProtonTimeNeg(-1000);
        tFemtoV0->SetTOFPionTimeNeg(-1000);
        tFemtoV0->SetTOFKaonTimeNeg(-1000);
      }
    } else {
      if (((tFemtoV0->StatusPos() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoV0->StatusPos() & AliVTrack::kTIME) == AliVTrack::kTIME) && probMisPos < 0.01) {

        // if(trackpos->IsOn(AliESDtrack::kTOFout & AliESDtrack::kTIME)) {
        tFemtoV0->SetPosNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kKaon));
        tFemtoV0->SetPosNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kProton));
        tFemtoV0->SetPosNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kPion));
      }
      if (((tFemtoV0->StatusNeg() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoV0->StatusNeg() & AliVTrack::kTIME) == AliVTrack::kTIME) && probMisNeg < 0.01) {

        // if(trackneg->IsOn(AliESDtrack::kTOFout & AliESDtrack::kTIME)) {
        tFemtoV0->SetNegNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kKaon));
        tFemtoV0->SetNegNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kProton));
        tFemtoV0->SetNegNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kPion));
      }
      double TOFSignalPos = trackpos->GetTOFsignal();
      double TOFSignalNeg = trackneg->GetTOFsignal();
      TOFSignalPos -= fAODpidUtil->GetTOFResponse().GetStartTime(trackpos->P());
      TOFSignalNeg -= fAODpidUtil->GetTOFResponse().GetStartTime(trackneg->P());
      double pidPos[5];
      double pidNeg[5];
      trackpos->GetIntegratedTimes(pidPos);
      trackneg->GetIntegratedTimes(pidNeg);

      tFemtoV0->SetTOFPionTimePos(TOFSignalPos - pidPos[2]);
      tFemtoV0->SetTOFKaonTimePos(TOFSignalPos - pidPos[3]);
      tFemtoV0->SetTOFProtonTimePos(TOFSignalPos - pidPos[4]);
      tFemtoV0->SetTOFPionTimeNeg(TOFSignalNeg - pidNeg[2]);
      tFemtoV0->SetTOFKaonTimeNeg(TOFSignalNeg - pidNeg[3]);
      tFemtoV0->SetTOFProtonTimeNeg(TOFSignalNeg - pidNeg[4]);
    }

  } else {
    tFemtoV0->SetStatusPos(999);
    tFemtoV0->SetStatusNeg(999);
  }

  tFemtoV0->SetOnFlyStatusV0(tAODv0->GetOnFlyStatus());
  return tFemtoV0;
}



AliFemtoXi *AliFemtoEventReaderAOD::CopyAODtoFemtoXi(AliAODcascade *tAODxi)
{

  AliFemtoXi *tFemtoXi = NULL;

  { // this is to keep tmpV0 in its own scope
    AliFemtoV0 *tmpV0 = CopyAODtoFemtoV0(tAODxi);
    tFemtoXi = new AliFemtoXi(tmpV0);
    delete tmpV0;
  }

  //The above lines set the V0 attributes.  However, some of these are now set incorrectly
  //For instance:
  //  tFemtoV0->SetdecayLengthV0(tAODv0->DecayLength(fV1)); in CopyAODtoFemtoXi now sets the decay length
  //    of the V0 roughly equal to the decay length of the cascade + decay length of V0
  //  This is fixed by:
  //    tFemtoXi->SetdecayLengthV0(tAODv0->DecayLengthV0());
  //  Or:
  //    double temp[3];
  //    temp[0] = tAODxi->DecayVertexXiX();
  //    temp[1] = tAODxi->DecayVertexXiY();
  //    temp[2] = tAODxi->DecayVertexXiZ();
  //    tFemtoXi->SetdecayLengthV0(tAODxi->DecayLength(temp));

  //-- Include any fixes to V0 attributes set in CopyAODtoFemtoXi
  tFemtoXi->SetdecayLengthV0(tAODxi->DecayLengthV0());


  //xi
  tFemtoXi->SetmassXi(tAODxi->MassXi());
  tFemtoXi->SetdecayLengthXi(tAODxi->DecayLengthXi(fV1[0],fV1[1],fV1[2]));
  tFemtoXi->SetdecayVertexXiX(tAODxi->DecayVertexXiX());
  tFemtoXi->SetdecayVertexXiY(tAODxi->DecayVertexXiY());
  tFemtoXi->SetdecayVertexXiZ(tAODxi->DecayVertexXiZ());
  tFemtoXi->SetdcaXiDaughters(tAODxi->DcaXiDaughters());
  //tFemtoXi->SetdcaXiToPrimVertex(tAODxi->DcaXiToPrimVertex());  //This doesn't work because fDcaXiToPrimVertex was not set
  tFemtoXi->SetdcaXiToPrimVertex(tAODxi->DcaXiToPrimVertex(fV1[0],fV1[1],fV1[2]));
  tFemtoXi->SetdcaBacToPrimVertex(tAODxi->DcaBachToPrimVertex());

  tFemtoXi->SetmomBacX(tAODxi->MomBachX());
  tFemtoXi->SetmomBacY(tAODxi->MomBachY());
  tFemtoXi->SetmomBacZ(tAODxi->MomBachZ());
  tFemtoXi->SetmomXiX(tAODxi->MomXiX());
  tFemtoXi->SetmomXiY(tAODxi->MomXiY());
  tFemtoXi->SetmomXiZ(tAODxi->MomXiZ());
  AliFemtoThreeVector momxi(tAODxi->MomXiX(), tAODxi->MomXiY(), tAODxi->MomXiZ());
  tFemtoXi->SetmomXi(momxi);

  tFemtoXi->SetidBac(tAODxi->GetBachID());

  double tEtaXi = 0.5*TMath::Log((TMath::Sqrt(tAODxi->Ptot2Xi())+tAODxi->MomXiZ())/(TMath::Sqrt(tAODxi->Ptot2Xi())-tAODxi->MomXiZ()+1.e-13));
  tFemtoXi->SetEtaXi(tEtaXi);
  double tPhiXi = TMath::Pi()+TMath::ATan2(-tAODxi->MomXiY(),-tAODxi->MomXiX());
  tFemtoXi->SetPhiXi(tPhiXi);

  tFemtoXi->SetCosPointingAngleXi(tAODxi->CosPointingAngleXi(fV1[0],fV1[1],fV1[2]));
  tFemtoXi->SetCosPointingAngleV0toXi(tAODxi->CosPointingAngle(tAODxi->GetDecayVertexXi()));

  tFemtoXi->SetChargeXi(tAODxi->ChargeXi());
  tFemtoXi->SetptXi(std::sqrt(tAODxi->Pt2Xi()));

  AliAODTrack *trackbac = (AliAODTrack *)tAODxi->GetDecayVertexXi()->GetDaughter(0);

  if (trackbac) {
    tFemtoXi->SetptBac(trackbac->Pt()); //setting pt? px and py was set!

    tFemtoXi->SetEtaBac(trackbac->Eta());//bac!
    tFemtoXi->SetTPCNclsBac(trackbac->GetTPCNcls()); //bac!
    tFemtoXi->SetNdofBac(trackbac->Chi2perNDF());//bac!
    tFemtoXi->SetStatusBac(trackbac->GetStatus()); //bac!

    tFemtoXi->SetBacNSigmaTPCK(fAODpidUtil->NumberOfSigmasTPC(trackbac, AliPID::kKaon));
    tFemtoXi->SetBacNSigmaTPCP(fAODpidUtil->NumberOfSigmasTPC(trackbac, AliPID::kProton));
    tFemtoXi->SetBacNSigmaTPCPi(fAODpidUtil->NumberOfSigmasTPC(trackbac, AliPID::kPion));


    //NEED TO ADD:  
    //    tFemtoXi->SetNominalTpcEntrancePointBac(tmpVec);
    //    tFemtoXi->SetNominalTpcExitPointBac(tmpVec);
    //    tFemtoXi->SetNominalTpcPointBac(vecTpcPos);
    //    tFemtoXi->SetTPCMomentumBac(trackbac->GetTPCmomentum());
    //    if (fShiftPosition > 0.)

    tFemtoXi->SetdedxBac(trackbac->GetTPCsignal());

    Float_t probMisBac = 1.0;

    if (((tFemtoXi->StatusBac() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoXi->StatusBac() & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      // if (tFemtoXi->StatusBac() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000
      probMisBac = fAODpidUtil->GetTOFMismatchProbability(trackbac);
    }

    // if(// (tFemtoXi->StatusPos()& AliESDtrack::kTOFpid)==0 ||
    //    (tFemtoXi->StatusPos()&AliESDtrack::kTIME)==0 || (tFemtoXi->StatusPos()&AliESDtrack::kTOFout)==0 || probMisPos > 0.01)


    if (!(((tFemtoXi->StatusBac() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoXi->StatusBac() & AliVTrack::kTIME) == AliVTrack::kTIME)) || probMisBac > 0.01) 
    {
        tFemtoXi->SetBacNSigmaTOFK(-1000);
	tFemtoXi->SetBacNSigmaTOFP(-1000);
	tFemtoXi->SetBacNSigmaTOFPi(-1000);

	tFemtoXi->SetTOFProtonTimeBac(-1000);
        tFemtoXi->SetTOFPionTimeBac(-1000);
        tFemtoXi->SetTOFKaonTimeBac(-1000);
    }
    else   
    {
      if (((tFemtoXi->StatusBac() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((tFemtoXi->StatusBac() & AliVTrack::kTIME) == AliVTrack::kTIME) && probMisBac < 0.01)
      {
        tFemtoXi->SetBacNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackbac, AliPID::kKaon));
        tFemtoXi->SetBacNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackbac, AliPID::kProton));
        tFemtoXi->SetBacNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackbac, AliPID::kPion));
      }

      double TOFSignalBac = trackbac->GetTOFsignal();
      TOFSignalBac -= fAODpidUtil->GetTOFResponse().GetStartTime(trackbac->P());
      double pidBac[5];
      trackbac->GetIntegratedTimes(pidBac);

      tFemtoXi->SetTOFPionTimeBac(TOFSignalBac - pidBac[2]);
      tFemtoXi->SetTOFKaonTimeBac(TOFSignalBac - pidBac[3]);
      tFemtoXi->SetTOFProtonTimeBac(TOFSignalBac - pidBac[4]);
    }
  }
  else 
  {
    tFemtoXi->SetStatusBac(999);
  }

  return tFemtoXi;
}



void AliFemtoEventReaderAOD::SetFilterBit(UInt_t ibit)
{
  fFilterBit = (1 << (ibit));
}


void AliFemtoEventReaderAOD::SetFilterMask(int ibit)
{
  fFilterMask = ibit;
}

void AliFemtoEventReaderAOD::SetReadMC(unsigned char a)
{
  fReadMC = a;
}


void AliFemtoEventReaderAOD::SetReadV0(unsigned char a)
{
  fReadV0 = a;
}

void AliFemtoEventReaderAOD::SetReadCascade(unsigned char a)
{
  fReadCascade = a;
}

void AliFemtoEventReaderAOD::SetUseMultiplicity(EstEventMult aType)
{
  fEstEventMult = aType;
}

AliAODMCParticle *AliFemtoEventReaderAOD::GetParticleWithLabel(TClonesArray *mcP, Int_t aLabel)
{
  if (aLabel < 0) return NULL;

  Int_t posstack = TMath::Min(aLabel, mcP->GetEntries());
  AliAODMCParticle *aodP = (AliAODMCParticle *) mcP->At(posstack);

  if (aodP->GetLabel() > posstack) {
    do {
      aodP = (AliAODMCParticle *) mcP->At(posstack);
      if (aodP->GetLabel() == aLabel) return aodP;
      posstack--;
    } while (posstack > 0);
  } else {
    do {
      aodP = (AliAODMCParticle *) mcP->At(posstack);
      if (aodP->GetLabel() == aLabel) return aodP;
      posstack++;
    } while (posstack < mcP->GetEntries());
  }

  return NULL;
}

void AliFemtoEventReaderAOD::CopyPIDtoFemtoTrack(AliAODTrack *tAodTrack, AliFemtoTrack *tFemtoTrack)
{
  // A cache which maps vertices to the number of tracks used to determine the vertex
  // Added due to slow calculation in AliAODVertex::GetNContributors - if that changes, remove this.
  static std::map<Short_t, Int_t> _vertex_NContributors_cache;

  if (fDCAglobalTrack==1) {

    // code from Michael and Prabhat from AliAnalysisTaskDptDptCorrelations
    const AliAODVertex *vertex = (AliAODVertex *) fEvent->GetPrimaryVertex();

    float vertexX = -999.;
    float vertexY = -999.;
    float vertexZ = -999.;

    if (vertex) {
      Double_t fCov[6] = {0.0};
      vertex->GetCovarianceMatrix(fCov);
      if (fCov[5] != 0) {
        vertexX = vertex->GetX();
        vertexY = vertex->GetY();
        vertexZ = vertex->GetZ();
      }
    }

    Double_t pos[3];
    tAodTrack->GetXYZ(pos);

    Double_t DCAX = pos[0] - vertexX;
    Double_t DCAY = pos[1] - vertexY;
    Double_t DCAZ = pos[2] - vertexZ;

    Double_t DCAXY = TMath::Sqrt((DCAX * DCAX) + (DCAY * DCAY));

    tFemtoTrack->SetImpactD(DCAXY);
    tFemtoTrack->SetImpactZ(DCAZ);

    tFemtoTrack->SetXatDCA(pos[0]);
    tFemtoTrack->SetYatDCA(pos[1]);
    tFemtoTrack->SetZatDCA(pos[2]);


  }
  else if(fDCAglobalTrack==2)
    {
      Double_t DCAXY = -999;
      Double_t DCAZ = -999;
      
      //DCA for TPC only - from PropagateToDCA method
      AliExternalTrackParam aliextparam;
      aliextparam.CopyFromVTrack(tAodTrack);
      if (aliextparam.GetX() > 3.0)
	{
	  DCAXY = -999;
	  DCAZ = -999;
	}
      else
	{
	  Double_t covar[3]={0,0,0};
	  Double_t DCA[2]={0,0};
	  if (!aliextparam.PropagateToDCA(fEvent->GetPrimaryVertex(),fEvent->GetMagneticField(),99999.0,DCA,covar))
	    {
	      DCAXY = -999;
	      DCAZ = -999;
	    }
	  else
	    {
	      DCAXY = DCA[0];
	      DCAZ = DCA[1];
	    }
	}
 
      tFemtoTrack->SetImpactD(DCAXY);
      tFemtoTrack->SetImpactZ(DCAZ);    
    }

  double aodpid[10];
  tAodTrack->GetPID(aodpid);
  tFemtoTrack->SetPidProbElectron(aodpid[0]);
  tFemtoTrack->SetPidProbMuon(aodpid[1]);
  tFemtoTrack->SetPidProbPion(aodpid[2]);
  tFemtoTrack->SetPidProbKaon(aodpid[3]);
  tFemtoTrack->SetPidProbProton(aodpid[4]);

  aodpid[0] = -100000.0;
  aodpid[1] = -100000.0;
  aodpid[2] = -100000.0;
  aodpid[3] = -100000.0;
  aodpid[4] = -100000.0;

  double tTOF = 0.0;
  Float_t probMis = 1.0;

  //what is that code? for what do we need it? nsigma values are not enough?
  // if (tAodTrack->GetStatus() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000

  ULong_t status = tAodTrack->GetStatus();

  if (((status & AliVTrack::kTOFout) == AliVTrack::kTOFout)
      && ((status & AliVTrack::kTIME) == AliVTrack::kTIME)) {
    tTOF = tAodTrack->GetTOFsignal();
    tAodTrack->GetIntegratedTimes(aodpid);

    tTOF -= fAODpidUtil->GetTOFResponse().GetStartTime(tAodTrack->P());

    probMis = fAODpidUtil->GetTOFMismatchProbability(tAodTrack);
  }



  tFemtoTrack->SetTofExpectedTimes(tTOF - aodpid[2], tTOF - aodpid[3], tTOF - aodpid[4]);

  //////  TPC ////////////////////////////////////////////

  const float nsigmaTPCK = fAODpidUtil->NumberOfSigmasTPC(tAodTrack, AliPID::kKaon);
  const float nsigmaTPCPi = fAODpidUtil->NumberOfSigmasTPC(tAodTrack, AliPID::kPion);
  const float nsigmaTPCP = fAODpidUtil->NumberOfSigmasTPC(tAodTrack, AliPID::kProton);
  const float nsigmaTPCE = fAODpidUtil->NumberOfSigmasTPC(tAodTrack, AliPID::kElectron);

  tFemtoTrack->SetNSigmaTPCPi(nsigmaTPCPi);
  tFemtoTrack->SetNSigmaTPCK(nsigmaTPCK);
  tFemtoTrack->SetNSigmaTPCP(nsigmaTPCP);
  tFemtoTrack->SetNSigmaTPCE(nsigmaTPCE);

  tFemtoTrack->SetTPCchi2(tAodTrack->Chi2perNDF());
  tFemtoTrack->SetTPCncls(tAodTrack->GetTPCNcls());
  tFemtoTrack->SetTPCnclsF(tAodTrack->GetTPCNcls());

  tFemtoTrack->SetTPCsignalN(1);
  tFemtoTrack->SetTPCsignalS(1);
  tFemtoTrack->SetTPCsignal(tAodTrack->GetTPCsignal());


  //////  TOF ////////////////////////////////////////////

  float vp = -1000.;
  float nsigmaTOFPi = -1000.;
  float nsigmaTOFK = -1000.;
  float nsigmaTOFP = -1000.;
  float nsigmaTOFE = -1000.;

  if (((status & AliVTrack::kTOFout) == AliVTrack::kTOFout)
      && ((status & AliVTrack::kTIME) == AliVTrack::kTIME)
      && probMis < 0.01) {

    nsigmaTOFPi = fAODpidUtil->NumberOfSigmasTOF(tAodTrack, AliPID::kPion);
    nsigmaTOFK = fAODpidUtil->NumberOfSigmasTOF(tAodTrack, AliPID::kKaon);
    nsigmaTOFP = fAODpidUtil->NumberOfSigmasTOF(tAodTrack, AliPID::kProton);
    nsigmaTOFE = fAODpidUtil->NumberOfSigmasTOF(tAodTrack, AliPID::kElectron);

    Double_t len = 200; // esdtrack->GetIntegratedLength(); !!!!!
    Double_t tof = tAodTrack->GetTOFsignal();
    if (tof > 0.) vp = len / tof / 0.03;
  }

  tFemtoTrack->SetVTOF(vp);
  tFemtoTrack->SetNSigmaTOFPi(nsigmaTOFPi);
  tFemtoTrack->SetNSigmaTOFK(nsigmaTOFK);
  tFemtoTrack->SetNSigmaTOFP(nsigmaTOFP);
  tFemtoTrack->SetNSigmaTOFE(nsigmaTOFE);

  //////////////////////////////////////
}

void AliFemtoEventReaderAOD::SetCentralityPreSelection(double min, double max)
{
  fCentRange[0] = min;
  fCentRange[1] = max;
  fUsePreCent = 1;
  fEstEventMult = kCentrality;
}

void AliFemtoEventReaderAOD::SetNoCentrality(bool anocent)
{
  if (anocent == false) {
    fEstEventMult = kCentrality;
  } else {
    fEstEventMult = kReference;
    fUsePreCent = 0;
  }
}

void AliFemtoEventReaderAOD::SetAODpidUtil(AliAODpidUtil *aAODpidUtil)
{
  fAODpidUtil = aAODpidUtil;
  //  printf("fAODpidUtil: %x\n",fAODpidUtil);
}

void AliFemtoEventReaderAOD::SetAODheader(AliAODHeader *aAODheader)
{
  fAODheader = aAODheader;
}


void AliFemtoEventReaderAOD::SetMagneticFieldSign(int s)
{
  if (s > 0)
    fMagFieldSign = 1;
  else if (s < 0)
    fMagFieldSign = -1;
  else
    fMagFieldSign = 0;
}

void AliFemtoEventReaderAOD::SetEPVZERO(Bool_t iepvz)
{
  fisEPVZ = iepvz;
}

void AliFemtoEventReaderAOD::GetGlobalPositionAtGlobalRadiiThroughTPC(AliAODTrack *track, Float_t bfield, Float_t globalPositionsAtRadii[9][3])
{
  // Gets the global position of the track at nine different radii in the TPC
  // params:
  //   track - the track to propagate
  //   bfield - magnetic field of event
  //   globalPositionsAtRadii - Output array of global positions in the radii and xyz
  const Float_t DEFAULT_VALUE = -9999.0;

  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  const Float_t Rwanted[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  // index of global position we are filling
  //  - first we use AliExternalTrackParam, then just default value
  Int_t radius_index = 0;

  // loop over the array of radii
  for (; radius_index < 9; radius_index++) {

    // extracted radius
    const Float_t radius = Rwanted[radius_index];
    // buffer to store position
    Double_t pos_buffer[3] = {0};

    // get the global position of the track at this radial location
    bool good = etp.GetXYZatR(radius, bfield, pos_buffer, NULL);

    // if value is not good, break loading loop
    if (!good || fabs(AliFemtoThreeVector(pos_buffer).Perp() - radius) > 0.5) {
      radius_index--; // decrement to fill current location with default value
      break;
    }

    // store the global position
    globalPositionsAtRadii[radius_index][0] = pos_buffer[0];
    globalPositionsAtRadii[radius_index][1] = pos_buffer[1];
    globalPositionsAtRadii[radius_index][2] = pos_buffer[2];
  }

  // Fill any remaining positions with the default value
  for (; radius_index < 9; radius_index++) {
    std::fill_n(globalPositionsAtRadii[radius_index], 3, DEFAULT_VALUE);
  }
}

//________________________________________________________________________
void AliFemtoEventReaderAOD::SetShiftedPositions(const AliAODTrack *track
						 ,const Float_t bfield
						 ,Float_t posShifted[3]
						 ,const Double_t radius) {
  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system, code adapted from Hans Beck analysis

  // Initialize the array to something indicating there was no propagation
  posShifted[0]=-9999.; // THIS IS THE DATA MEMBER OF YOUR FEMTOTRACK
  posShifted[1]=-9999.;
  posShifted[2]=-9999.;
   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};

  // The radius in cm we want to propagate to, squared
  const Float_t RSquaredWanted(radius*radius*1e4);

  // Propagation is done in local x of the track
  for (Float_t x = 58.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent.
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared.
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-fV1[0])*(xyz[0]-fV1[0])
                                 + (xyz[1]-fV1[1])*(xyz[1]-fV1[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted){

      // Bigger loop has bad precision, we're nearly one centimeter too far,
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted){
	// Propagate a mm inwards
	x-=.1;
	if(!etp.PropagateTo(x,bfield)){
	  // Propagation failed but we're already with a
	  // cm precision at R=1.25m so we only break the
	  // inner loop
	  break;
	}
	// Get the global position
	etp.GetXYZ(xyz);
	// Calculate shifted radius, squared
	shiftedRadiusSquared = (xyz[0]-fV1[0])*(xyz[0]-fV1[0])
	                     + (xyz[1]-fV1[1])*(xyz[1]-fV1[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      posShifted[0]=xyz[0]-fV1[0];
      posShifted[1]=xyz[1]-fV1[1];
      posShifted[2]=xyz[2]-fV1[2];
      // Done
      return;
    } // End of if roughly reached radius
  } // End of coarse propagation loop
}

void AliFemtoEventReaderAOD::SetpA2013(Bool_t pa2013)
{
  fpA2013 = pa2013;
}

void AliFemtoEventReaderAOD::SetUseMVPlpSelection(Bool_t mvplp)
{
  fMVPlp = mvplp;
}

void AliFemtoEventReaderAOD::SetIsPileUpEvent(Bool_t ispileup)
{
  fisPileUp = ispileup;
}

void AliFemtoEventReaderAOD::SetDCAglobalTrack(Int_t dcagt)
{
  fDCAglobalTrack = dcagt;
}

bool AliFemtoEventReaderAOD::RejectEventCentFlat(float MagField, float CentPercent)
{
  // Flattens the centrality distribution

  // Setting 0 as seed ensures random seed every time.
  TRandom3 RNG(0); // for 3D, random sign switching

  float kCentWeight[2][9] = {
    {.878, .876, .860, .859, .859, .880, .873, .879, .894},
    {.828, .793, .776, .772, .775, .796, .788, .804, .839}
  };

  int weightBinCent = (int) CentPercent,
      weightBinSign = (MagField > 0) ? 0 : 1;

  bool rejectEvent = RNG.Rndm() > kCentWeight[weightBinSign][weightBinCent];
  return rejectEvent;
}

void AliFemtoEventReaderAOD::SetCentralityFlattening(Bool_t dcagt)
{
  fFlatCent = dcagt;
}

void AliFemtoEventReaderAOD::SetShiftPosition(Double_t dcagt)
{
  fShiftPosition = dcagt;
}

void AliFemtoEventReaderAOD::SetPrimaryVertexCorrectionTPCPoints(bool correctTpcPoints)
{
  fPrimaryVertexCorrectionTPCPoints = correctTpcPoints;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsPions(TH1D *h1)
{
  f1DcorrectionsPions = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsKaons(TH1D *h1)
{
  f1DcorrectionsKaons = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsProtons(TH1D *h1)
{
  f1DcorrectionsProtons = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsPionsMinus(TH1D *h1)
{
  f1DcorrectionsPionsMinus = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsKaonsMinus(TH1D *h1)
{
  f1DcorrectionsKaonsMinus = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsProtonsMinus(TH1D *h1)
{
  f1DcorrectionsProtonsMinus = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsAll(TH1D *h1)
{
  f1DcorrectionsAll = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsLambdas(TH1D *h1)
{
  f1DcorrectionsLambdas = h1;
}

void AliFemtoEventReaderAOD::Set1DCorrectionsLambdasMinus(TH1D *h1)
{
  f1DcorrectionsLambdasMinus = h1;
}

//Special MC analysis for pi,K,p,e selected by PDG code -->
void AliFemtoEventReaderAOD::SetPionAnalysis(Bool_t aSetPionAna)
{
  fIsPionAnalysis = aSetPionAna;
}
void AliFemtoEventReaderAOD::SetKaonAnalysis(Bool_t aSetKaonAna)
{
  fIsKaonAnalysis = aSetKaonAna;
}
void AliFemtoEventReaderAOD::SetProtonAnalysis(Bool_t aSetProtonAna)
{
  fIsProtonAnalysis = aSetProtonAna;
}
void AliFemtoEventReaderAOD::SetElectronAnalysis(Bool_t aSetElectronAna)
{
  fIsElectronAnalysis = aSetElectronAna;
}
//Special MC analysis for pi,K,p,e selected by PDG code <--
