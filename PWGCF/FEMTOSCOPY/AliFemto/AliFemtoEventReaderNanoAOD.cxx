///
/// \file AliFemtoEventReaderNanoAOD.cxx
///

#include "AliFemtoEventReaderNanoAOD.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"

#include "AliNanoAODHeader.h"
//#include "/home/wfpw/alice5/AliPhysics/PWG/DevNanoAOD/AliNanoAODTrack.h"
#include "AliNanoAODTrack.h"

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

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include <algorithm>
#include <cassert>
#include <map>


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventReaderNanoAOD);
  /// \endcond
#endif

#if !(ST_NO_NAMESPACES)
using namespace units;
#endif

using namespace std;

double fVe1[3];

//____________________________
//constructor with 0 parameters , look at default settings
AliFemtoEventReaderNanoAOD::AliFemtoEventReaderNanoAOD():
  fNumberofEvent(0),
  fCurEvent(0),
  fEvent(nullptr),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fFilterMask(0),
  //  fPWG2AODTracks(0x0),
  fReadMC(0),
  fReadV0(0),
  fReadCascade(0),
  fUsePreCent(0),
  fEstEventMult("RefMult08"),
  fAODpidUtil(nullptr),
  fAODheader(nullptr),
  fAnaUtils(nullptr),
  fEventCuts(nullptr),
  fUseAliEventCuts(0),
  fReadFullMCData(false),
  fInputFile(""),
  fTree(nullptr),
  fAodFile(nullptr),
  fMagFieldSign(1),
  fisEPVZ(kTRUE),
  fpA2013(kFALSE),
  fisPileUp(kFALSE),
  fCascadePileUpRemoval(kFALSE),
  fV0PileUpRemoval(kFALSE),
  fTrackPileUpRemoval(kFALSE),
  fMVPlp(kFALSE),
  fOutOfBunchPlp(kFALSE),
  fMinVtxContr(0),
  fMinPlpContribMV(0),
  fMinPlpContribSPD(0),
  fDCAglobalTrack(0),
  fFlatCent(kFALSE),
  fPrimaryVertexCorrectionTPCPoints(kFALSE),
  fShiftPosition(0.),
  fCovMatPresent(kTRUE),
  f1DcorrectionsPions(0),
  f1DcorrectionsKaons(0),
  f1DcorrectionsProtons(0),
  f1DcorrectionsDeuterons(0),
  f1DcorrectionsPionsMinus(0),
  f1DcorrectionsKaonsMinus(0),
  f1DcorrectionsProtonsMinus(0),
  f1DcorrectionsDeuteronsMinus(0),
  f1DcorrectionsAll(0),
  f1DcorrectionsLambdas(0),
  f1DcorrectionsLambdasMinus(0),
  f1DcorrectionsXiPlus(0),
  f1DcorrectionsXiMinus(0),
  f4DcorrectionsPions(0),
  f4DcorrectionsKaons(0),
  f4DcorrectionsProtons(0),
  f4DcorrectionsDeuterons(0),
  f4DcorrectionsPionsMinus(0),
  f4DcorrectionsKaonsMinus(0),
  f4DcorrectionsProtonsMinus(0),
  f4DcorrectionsDeuteronsMinus(0),
  f4DcorrectionsAll(0),
  f4DcorrectionsLambdas(0),
  f4DcorrectionsLambdasMinus(0)
{
  // default constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fCentRange[0] = 0;
  fCentRange[1] = 1000;
}

AliFemtoEventReaderNanoAOD::AliFemtoEventReaderNanoAOD(const AliFemtoEventReaderNanoAOD &aReader):
  AliFemtoEventReader(),
  fNumberofEvent(aReader.fNumberofEvent),
  fCurEvent(aReader.fCurEvent),
  //fEvent(new AliAODEvent()),
  fEvent(0),
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
  fEventCuts(aReader.fEventCuts),
  fUseAliEventCuts(aReader.fUseAliEventCuts),
  fReadFullMCData(aReader.fReadFullMCData),
  fInputFile(aReader.fInputFile),
  fTree(nullptr),
  fAodFile(new TFile(aReader.fAodFile->GetName())),
  fMagFieldSign(aReader.fMagFieldSign),
  fisEPVZ(aReader.fisEPVZ),
  fpA2013(aReader.fpA2013),
  fisPileUp(aReader.fisPileUp),
  fCascadePileUpRemoval(aReader.fCascadePileUpRemoval),
  fV0PileUpRemoval(aReader.fV0PileUpRemoval),
  fTrackPileUpRemoval(aReader.fTrackPileUpRemoval),
  fMVPlp(aReader.fMVPlp),
  fOutOfBunchPlp(aReader.fOutOfBunchPlp),
  fMinVtxContr(aReader.fMinVtxContr),
  fMinPlpContribMV(aReader.fMinPlpContribMV),
  fMinPlpContribSPD(aReader.fMinPlpContribSPD),
  fDCAglobalTrack(aReader.fDCAglobalTrack),
  fFlatCent(aReader.fFlatCent),
  fPrimaryVertexCorrectionTPCPoints(aReader.fPrimaryVertexCorrectionTPCPoints),
  fShiftPosition(aReader.fShiftPosition),
  fCovMatPresent(kTRUE),
  f1DcorrectionsPions(aReader.f1DcorrectionsPions),
  f1DcorrectionsKaons(aReader.f1DcorrectionsKaons),
  f1DcorrectionsProtons(aReader.f1DcorrectionsProtons),
  f1DcorrectionsDeuterons(aReader.f1DcorrectionsDeuterons),
  f1DcorrectionsPionsMinus(aReader.f1DcorrectionsPions),
  f1DcorrectionsKaonsMinus(aReader.f1DcorrectionsKaonsMinus),
  f1DcorrectionsProtonsMinus(aReader.f1DcorrectionsProtonsMinus),
  f1DcorrectionsDeuteronsMinus(aReader.f1DcorrectionsDeuteronsMinus),
  f1DcorrectionsAll(aReader.f1DcorrectionsAll),
  f1DcorrectionsLambdas(aReader.f1DcorrectionsLambdas),
  f1DcorrectionsLambdasMinus(aReader.f1DcorrectionsLambdasMinus),
  f1DcorrectionsXiPlus(aReader.f1DcorrectionsXiPlus),
  f1DcorrectionsXiMinus(aReader.f1DcorrectionsXiMinus),
  f4DcorrectionsPions(aReader.f4DcorrectionsPions),
  f4DcorrectionsKaons(aReader.f4DcorrectionsKaons),
  f4DcorrectionsProtons(aReader.f4DcorrectionsProtons),
  f4DcorrectionsDeuterons(aReader.f4DcorrectionsDeuterons),
  f4DcorrectionsPionsMinus(aReader.f4DcorrectionsPionsMinus),
  f4DcorrectionsKaonsMinus(aReader.f4DcorrectionsKaonsMinus),
  f4DcorrectionsProtonsMinus(aReader.f4DcorrectionsProtonsMinus),
  f4DcorrectionsDeuteronsMinus(aReader.f4DcorrectionsDeuteronsMinus),
  f4DcorrectionsAll(aReader.f4DcorrectionsAll),
  f4DcorrectionsLambdas(aReader.f4DcorrectionsLambdas),
  f4DcorrectionsLambdasMinus(aReader.f4DcorrectionsLambdasMinus)
{
  // copy constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);

  //  fPWG2AODTracks = aReader.fPWG2AODTracks;

  fCentRange[0] = aReader.fCentRange[0];
  fCentRange[1] = aReader.fCentRange[1];
}
//__________________
AliFemtoEventReaderNanoAOD::~AliFemtoEventReaderNanoAOD()
{ // destructor
  delete fEventCuts;
  delete fTree;
  delete fEvent;
  delete fAodFile;
  //   if (fPWG2AODTracks) {
  //     fPWG2AODTracks->Delete();
  //     delete fPWG2AODTracks;
  //   }
}

//__________________
AliFemtoEventReaderNanoAOD &AliFemtoEventReaderNanoAOD::operator=(const AliFemtoEventReaderNanoAOD &aReader)
{ // assignment operator
  if (this == &aReader) {
    return *this;
  }

  fInputFile = aReader.fInputFile;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  delete fTree;
  delete fEvent;
  //fEvent = new AliAODEvent();
  //fEvent = new AliVEvent();
  delete fAodFile;
  fAodFile = new TFile(aReader.fAodFile->GetName());
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fFilterBit = aReader.fFilterBit;
  fFilterMask = aReader.fFilterMask;
  fReadMC = aReader.fReadMC;
  fReadV0 = aReader.fReadV0;
  fReadCascade = aReader.fReadCascade;
  fReadFullMCData = aReader.fReadFullMCData;
  //  fPWG2AODTracks = aReader.fPWG2AODTracks;
  fAODpidUtil = aReader.fAODpidUtil;
  fAODheader = aReader.fAODheader;
  fAnaUtils = aReader.fAnaUtils;
  fEventCuts = aReader.fEventCuts;
  fUseAliEventCuts = aReader.fUseAliEventCuts;
  fCentRange[0] = aReader.fCentRange[0];
  fCentRange[1] = aReader.fCentRange[1];
  fUsePreCent = aReader.fUsePreCent;
  fEstEventMult = aReader.fEstEventMult;
  fMagFieldSign = aReader.fMagFieldSign;
  fpA2013 = aReader.fpA2013;
  fisPileUp = aReader.fisPileUp;
  fCascadePileUpRemoval = aReader.fCascadePileUpRemoval;
  fV0PileUpRemoval = aReader.fV0PileUpRemoval;
  fTrackPileUpRemoval = aReader.fTrackPileUpRemoval;
  fMVPlp = aReader.fMVPlp;
  fOutOfBunchPlp = aReader.fOutOfBunchPlp;
  fMinVtxContr = aReader.fMinVtxContr;
  fMinPlpContribMV = aReader.fMinPlpContribMV;
  fMinPlpContribSPD = aReader.fMinPlpContribSPD;
  fDCAglobalTrack = aReader.fDCAglobalTrack;
  fFlatCent = aReader.fFlatCent;
  fPrimaryVertexCorrectionTPCPoints = aReader.fPrimaryVertexCorrectionTPCPoints;
  fShiftPosition = aReader.fShiftPosition;
  fCovMatPresent = aReader.fCovMatPresent;

  f1DcorrectionsPions = aReader.f1DcorrectionsPions;
  f1DcorrectionsKaons = aReader.f1DcorrectionsKaons;
  f1DcorrectionsProtons = aReader.f1DcorrectionsProtons;
  f1DcorrectionsDeuterons = aReader.f1DcorrectionsDeuterons;
  f1DcorrectionsPionsMinus = aReader.f1DcorrectionsPionsMinus;
  f1DcorrectionsKaonsMinus = aReader.f1DcorrectionsKaonsMinus;
  f1DcorrectionsProtonsMinus = aReader.f1DcorrectionsProtonsMinus;
  f1DcorrectionsDeuteronsMinus = aReader.f1DcorrectionsDeuteronsMinus;
  f1DcorrectionsAll = aReader.f1DcorrectionsAll;
  f1DcorrectionsLambdas = aReader.f1DcorrectionsLambdas;
  f1DcorrectionsLambdasMinus = aReader.f1DcorrectionsLambdasMinus;
  f1DcorrectionsXiPlus = aReader.f1DcorrectionsXiPlus;
  f1DcorrectionsXiMinus = aReader.f1DcorrectionsXiMinus;
  
  f4DcorrectionsPions = aReader.f4DcorrectionsPions;
  f4DcorrectionsKaons = aReader.f4DcorrectionsKaons;
  f4DcorrectionsProtons = aReader.f4DcorrectionsProtons;
  f4DcorrectionsDeuterons = aReader.f4DcorrectionsDeuterons;
  f4DcorrectionsPionsMinus = aReader.f4DcorrectionsPionsMinus;
  f4DcorrectionsKaonsMinus = aReader.f4DcorrectionsKaonsMinus;
  f4DcorrectionsProtonsMinus = aReader.f4DcorrectionsProtonsMinus;
  f4DcorrectionsDeuteronsMinus = aReader.f4DcorrectionsDeuteronsMinus;
  f4DcorrectionsAll = aReader.f4DcorrectionsAll;
  f4DcorrectionsLambdas = aReader.f4DcorrectionsLambdas;
  f4DcorrectionsLambdasMinus = aReader.f4DcorrectionsLambdasMinus;

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderNanoAOD::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderNanoAOD\n";
  return temp;
}

//__________________
void AliFemtoEventReaderNanoAOD::SetInputFile(const char *inputFile)
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
      TTree *tree = static_cast<TTree*>(aodFile->Get("aodTree"));
      if (tree) {
        fTree->AddFile(filename);
        delete tree;
      }
      aodFile->Close();
    }
    delete aodFile;
  }
}

AliFemtoEvent *AliFemtoEventReaderNanoAOD::ReturnHbtEvent()
{
  /// Reads in the next event from the chain and converts it to an AliFemtoEvent
  AliFemtoEvent *hbtEvent = nullptr;

  // We have hit the end of our range -> open the next file
  if (fCurEvent == fNumberofEvent) {
    // We haven't loaded anything yet - open
    if (fNumberofEvent == 0) {
      // cout << "fEvent: " << fEvent << "\n";
      //fEvent = new AliVEvent();
      //fEvent->ReadFromTree(fTree);

      fNumberofEvent = fTree->GetEntries();
      // cout << "Number of entries in file " << fNumberofEvent << endl;
      fCurEvent = 0;
    } else {
      //no more data to read
      fReaderStatus = 1;
      return nullptr;
    }
  }

  // cout << "starting to read event " << fCurEvent << endl;
  fTree->GetEvent(fCurEvent);
  // cout << "Read event " << fEvent << " from file " << fTree << endl;

  hbtEvent = CopyAODtoFemtoEvent();
  fCurEvent++;

  return hbtEvent;
}

AliFemtoEvent *AliFemtoEventReaderNanoAOD::CopyAODtoFemtoEvent()
{

  // A function that reads in the AOD event
  // and transfers the neccessary information into
  // the internal AliFemtoEvent

  AliFemtoEvent *tEvent = new AliFemtoEvent();



  //AliNanoAODHeader *header = dynamic_cast<AliNanoAODHeader *>(fEvent->GetHeader());
  //assert(header && "Not a standard AOD");

  //tEvent->SetReactionPlaneAngle(header->GetQTheta(0) / 2.0);


  // Primary Vertex position
  const auto *aodvertex = static_cast<const AliAODVertex *>(fEvent->GetPrimaryVertex());

  if (!aodvertex || aodvertex->GetNContributors() < 1) {
    delete tEvent;  // Bad vertex, skip event.
    return nullptr;
  }

  aodvertex->GetPosition(fVe1);
  AliFmThreeVectorF vertex(fVe1[0], fVe1[1], fVe1[2]);
  tEvent->SetPrimVertPos(vertex);

  //starting to reading tracks
  int nofTracks = 0; //number of reconstructed tracks in event

  nofTracks = fEvent->GetNumberOfTracks();


  AliCentrality *cent = fEvent->GetCentrality();


  //if (!fEstEventMult && cent && fUsePreCent) {
  //if ((cent->GetCentralityPercentile("V0M") * 10 < fCentRange[0]) ||
  //    (cent->GetCentralityPercentile("V0M") * 10 > fCentRange[1])) {
  //  delete tEvent;
  //  return nullptr;
  //}
  //}

  //const Float_t percent = cent->GetCentralityPercentile("V0M");

  // Flatten centrality distribution
  //if (percent < 9 && fFlatCent) {
  //bool reject_event = RejectEventCentFlat(fEvent->GetMagneticField(), percent);
  //if (reject_event) {
  //  delete tEvent;
  //  return nullptr;
  //}
    //}

  int realnofTracks = 0; // number of track which we use in a analysis
  int tracksPrim = 0;

  // constant indicating label has been unset
  const int UNDEFINED_LABEL = -1;

  // 'labels' maps a track's id to the track's index in the Event
  // i.e. labels[Event->GetTrack(x)->GetID()] == x
  std::vector<int> labels(nofTracks, UNDEFINED_LABEL);

  // looking for global tracks and saving their numbers to copy from
  // them PID information to TPC-only tracks in the main loop over tracks
  /*
  for (int i = 0; i < nofTracks; i++) {
    const auto *aodtrack = static_cast<AliAODTrack *>(fEvent->GetTrack(i));

    if (!aodtrack->TestFilterBit(fFilterBit)) {
      // Skip TPC-only tracks
      const int id = aodtrack->GetID();
      if (id < 0) {
        continue;
      }

      // Resize labels vector if "id" is larger than mapping allows
      if (static_cast<size_t>(id) >= labels.size()) {
        labels.resize(id + 1024, UNDEFINED_LABEL);
      }
      labels[id] = i;
    }
  }
  */

  int tNormMult = 0;
  double norm_mult = 0;
  for (int i = 0; i < nofTracks; i++) {
    //  const AliAODTrack *aodtrack=dynamic_cast<AliAODTrack*>(fEvent->GetTrack(i));
    AliNanoAODTrack *aodtrack = static_cast<AliNanoAODTrack *>(fEvent->GetTrack(i));


    assert(aodtrack && "Not a standard AOD"); // Getting the AODtrack directly

    if (aodtrack->IsPrimary()) {
      tracksPrim++;
    }

    if ((fFilterBit && !aodtrack->TestFilterBit(fFilterBit)) ||
        (fFilterMask && !aodtrack->TestFilterBit(fFilterMask))) {
      continue;
    }

    // Check the sanity of the tracks - reject zero momentum tracks
    if (aodtrack->P() == 0.0) {
      continue;
    }

    // Counting particles to set multiplicity
    if ((fEstEventMult == "kGlobalCount")
      //&& (aodtrack->IsPrimaryCandidate()) //? instead of kinks?
        && (aodtrack->Chi2perNDF() < 4.0)
        && (0.15 <= aodtrack->Pt() && aodtrack->Pt() < 20)
        && (aodtrack->GetTPCNcls() > 70)
        && (aodtrack->Eta() < 0.8)) {
      tNormMult++;
    }

    norm_mult = tracksPrim;


    char fEstEventMult_cstr[fEstEventMult.size() + 1];
    strcpy(fEstEventMult_cstr, fEstEventMult.c_str());
    if (fEstEventMult.find("MultSelection.") != std::string::npos)
      {
	static const Int_t kMult = fAODheader->GetVarIndex(fEstEventMult_cstr);
	norm_mult  = fAODheader->GetVar(kMult);
      }
    else if(fEstEventMult == "TRK")
      {
	norm_mult = 10*fAODheader->GetCentr("TRK");
      }
    else if(fEstEventMult == "CL0")
      {
	norm_mult = 10*fAODheader->GetCentr("CL0");
      }
    else if(fEstEventMult == "CL1")
      {
	norm_mult = 10*fAODheader->GetCentr("CL0");
      }
    else if(fEstEventMult == "V0M")
      {
	norm_mult = 10*fAODheader->GetCentr("V0M");
      }
    else
      {
	norm_mult = 10*fAODheader->GetCentr(fEstEventMult_cstr);
      }

    tEvent->SetNormalizedMult(norm_mult);

    AliFemtoTrack *trackCopy = CopyAODtoFemtoTrack(aodtrack);

    trackCopy->SetMultiplicity(norm_mult);
    trackCopy->SetZvtx(fVe1[2]);

    // copying PID information from the correspondent track
    // not possible in NanoAOD?
    /*const AliAODTrack *aodtrackpid = fEvent->GetTrack(labels[-1-fEvent->GetTrack(i)->GetID()]);

    // For TPC Only tracks we have to copy PID information from corresponding global tracks
    const Int_t pid_track_id = (fFilterBit == (1 << 7) || fFilterMask == 128)
                           ? labels[-1 - fEvent->GetTrack(i)->GetID()]
                           : i;
    AliNanoAODTrack *aodtrackpid = static_cast<AliNanoAODTrack *>(fEvent->GetTrack(pid_track_id));
    assert(aodtrackpid && "Not a standard AOD");*/

    //CopyPIDtoFemtoTrack(aodtrackpid, trackCopy); //not needed anymore
    //now way simpler PID
    static Bool_t bPIDAvailable = AliNanoAODTrack::InitPIDIndex();


    if (aodtrack && bPIDAvailable)
      {
      	//////  TPC ////////////////////////////////////////////
	static const Int_t kcstNSigmaTPCPi  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kPion);
	static const Int_t kcstNSigmaTPCK  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kKaon);
	static const Int_t kcstNSigmaTPCPr  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kProton);
	static const Int_t kcstNSigmaTPCE  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kElectron);
	static const Int_t kcstNSigmaTPCD  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kDeuteron);

	const float nsigmaTPCK = aodtrack->GetVar(kcstNSigmaTPCK);
	const float nsigmaTPCPi = aodtrack->GetVar(kcstNSigmaTPCPi);
	const float nsigmaTPCP = aodtrack->GetVar(kcstNSigmaTPCPr);
	const float nsigmaTPCE = aodtrack->GetVar(kcstNSigmaTPCE);
	const float nsigmaTPCD = aodtrack->GetVar(kcstNSigmaTPCD);

	trackCopy->SetNSigmaTPCPi(nsigmaTPCPi);
	trackCopy->SetNSigmaTPCK(nsigmaTPCK);
	trackCopy->SetNSigmaTPCP(nsigmaTPCP);
	trackCopy->SetNSigmaTPCE(nsigmaTPCE);
	trackCopy->SetNSigmaTPCD(nsigmaTPCD);

 	trackCopy->SetTPCsignal(aodtrack->GetTPCsignal());
//	trackCopy->SetTPCsignalS(1);
//	trackCopy->SetTPCsignalN(aodtrack->GetTPCsignalN());


	trackCopy->SetITSchi2(aodtrack->GetITSchi2());
	//trackCopy->SetITSncls(aodtrack->GetITSNcls()); //not implemented

	for (int ii = 0; ii < 6; ii++) {
	  trackCopy->SetITSHitOnLayer(ii, aodtrack->HasPointOnITSLayer(ii));
	}

  	//////  TOF ////////////////////////////////////////////
	static const Int_t kcstNSigmaTOFPi  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kPion);
	static const Int_t kcstNSigmaTOFK  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kKaon);
	static const Int_t kcstNSigmaTOFPr  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kProton);
        static const Int_t kcstNSigmaTOFE  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kElectron);
	static const Int_t kcstNSigmaTOFD  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kDeuteron);

	const float nsigmaTOFK = aodtrack->GetVar(kcstNSigmaTOFK);
	const float nsigmaTOFPi = aodtrack->GetVar(kcstNSigmaTOFPi);
	const float nsigmaTOFP = aodtrack->GetVar(kcstNSigmaTOFPr);
	const float nsigmaTOFE = aodtrack->GetVar(kcstNSigmaTOFE);
	const float nsigmaTOFD = aodtrack->GetVar(kcstNSigmaTOFD);

	trackCopy->SetNSigmaTOFPi(nsigmaTOFPi);
	trackCopy->SetNSigmaTOFK(nsigmaTOFK);
	trackCopy->SetNSigmaTOFP(nsigmaTOFP);
	trackCopy->SetNSigmaTOFE(nsigmaTOFE);
	trackCopy->SetNSigmaTOFD(nsigmaTOFD);

//	trackCopy->SetTOFsignal(aodtrack->GetTOFsignal());

      }



  /*************************************************************************************/


    //AliExternalTrackParam *param = new AliExternalTrackParam(*aodtrack->GetInnerParam());
    trackCopy->SetInnerMomentum(aodtrack->GetTPCmomentum());

    tEvent->TrackCollection()->push_back(trackCopy); // Adding track to analysis
    realnofTracks++;                                 // Real number of tracks
  }

  //cout<<"======================> realnofTracks"<<realnofTracks<<endl;
  tEvent->SetNumberOfTracks(realnofTracks); // Setting number of track which we read in event

  /*
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
  */


  if (fReadV0) {
    AliAODEvent* aod = dynamic_cast<AliAODEvent*> (fEvent);
    int count_pass = 0;
    if (aod->GetV0s()) {
    for (Int_t i = 0; i < aod->GetNumberOfV0s(); i++) {
      AliAODv0 *aodv0 = aod->GetV0(i);
      // ensure a "good" v0 particle passes these conditions
      if (!aodv0
          || aodv0->GetNDaughters() > 2
          || aodv0->GetNProngs() > 2
          || aodv0->GetCharge() != 0
          || aodv0->ChargeProng(0) == aodv0->ChargeProng(1)
          || aodv0->CosPointingAngle(fVe1) < 0.98) {
        continue;
      }

      AliNanoAODTrack *daughterTrackPos = (AliNanoAODTrack *)aodv0->GetDaughter(0), // getting positive daughter track
                      *daughterTrackNeg = (AliNanoAODTrack *)aodv0->GetDaughter(1); // getting negative daughter track

     if (daughterTrackPos == nullptr || daughterTrackNeg == nullptr) {
        continue; // daughter tracks must exist
      }

      if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) {
        continue; // and have different charge
      }

      // ensure that pos and neg are pointing to the correct children
      if (daughterTrackPos->Charge() < 0 && daughterTrackNeg->Charge() > 0) {
        std::swap(daughterTrackPos, daughterTrackNeg); // can we use this?
      }


      AliFemtoV0 *trackCopyV0 = CopyAODtoFemtoV0(aodv0);
      trackCopyV0->SetMultiplicity(norm_mult);
      trackCopyV0->SetZvtx(fVe1[2]);

      tEvent->V0Collection()->push_back(trackCopyV0);
      count_pass++;

    }
  }
}


  if (fReadCascade) {
    int count_pass = 0;
    AliAODEvent* aod = dynamic_cast<AliAODEvent*> (fEvent);
    for (Int_t i = 0; i < aod->GetNumberOfCascades(); i++) {
      AliAODcascade *aodxi = aod->GetCascade(i);
      if (!aodxi) {
        continue;
      }
      //if (aodxi->GetNDaughters() > 2) continue;
      //if (aodxi->GetNProngs() > 2) continue;
      //if (aodxi->GetCharge() != 0) continue;
      if ((aodxi->ChargeProng(0) == aodxi->ChargeProng(1))
          || (aodxi->CosPointingAngle(fVe1) < 0.9)
          || (aodxi->CosPointingAngleXi(fVe1[0], fVe1[1], fVe1[2]) < 0.98)) {
        continue;
      }

      AliAODTrack *daughterTrackPos = (AliAODTrack *)aodxi->GetDaughter(0), // getting positive daughter track
                  *daughterTrackNeg = (AliAODTrack *)aodxi->GetDaughter(1), // getting negative daughter track
                  *bachTrack = (AliAODTrack *)aodxi->GetDecayVertexXi()->GetDaughter(0);

      if (daughterTrackPos == nullptr || daughterTrackNeg == nullptr || bachTrack == nullptr) {
        continue; // daughter tracks must exist
      }
      if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) {
        continue; // and have different charge
      }


      AliFemtoXi *trackCopyXi = CopyAODtoFemtoXi(aodxi);
      tEvent->XiCollection()->push_back(trackCopyXi);
      count_pass++;
    }
  }

  return tEvent;
}

AliFemtoTrack *AliFemtoEventReaderNanoAOD::CopyAODtoFemtoTrack(AliNanoAODTrack *tAodTrack)
{
  // Copy the track information from the AOD into the internal AliFemtoTrack
  // If it exists, use the additional information from the PWG2 AOD
  AliFemtoTrack *tFemtoTrack = new AliFemtoTrack();

  // Primary Vertex position
  const AliVVertex* vertex = fEvent->GetPrimaryVertex();
  // if (vertex->GetNContributors() < 1)
  //return;
  vertex->GetXYZ(fVe1);

  //fEvent->GetPrimaryVertex()->GetPosition(fVe1);
  //  fEvent->GetPrimaryVertex()->GetXYZ(fVe1);
  tFemtoTrack->SetPrimaryVertex(fVe1);

  tFemtoTrack->SetCharge(tAodTrack->Charge());

  double pxyz[3];
  tAodTrack->PxPyPz(pxyz); //reading noconstrained momentum


  AliFemtoThreeVector v(pxyz[0], pxyz[1], pxyz[2]);
  tFemtoTrack->SetP(v); //setting momentum
  tFemtoTrack->SetPt(sqrt(pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1]));
  const AliFmThreeVectorD kOrigin(fVe1[0], fVe1[1], fVe1[2]);
  //setting track helix
  const AliFmThreeVectorD ktP(pxyz[0], pxyz[1], pxyz[2]);
  AliFmPhysicalHelixD helix(ktP, kOrigin, (double)(fEvent->GetMagneticField()) * kilogauss, (double)(tFemtoTrack->Charge()));
  tFemtoTrack->SetHelix(helix);

  // Flags
  tFemtoTrack->SetTrackId(tAodTrack->GetID());
  tFemtoTrack->SetFlags(tAodTrack->GetFlag());
  tFemtoTrack->SetLabel(tAodTrack->GetLabel());

  // Track quality information
  //Double_t covmat[21];
  //tAodTrack->GetCovarianceXYZPxPyPz(covmat);

  if (fDCAglobalTrack == 0) {
    tFemtoTrack->SetImpactD(tAodTrack->DCA());
    tFemtoTrack->SetImpactZ(tAodTrack->ZAtDCA());

    tFemtoTrack->SetXatDCA(tAodTrack->XAtDCA());
    tFemtoTrack->SetYatDCA(tAodTrack->YAtDCA());
    tFemtoTrack->SetZatDCA(tAodTrack->ZAtDCA());
  }
  else if (fDCAglobalTrack == 1) {

    auto *vertex = static_cast<const AliAODVertex *>(fEvent->GetPrimaryVertex());
//    const Int_t pid_track_id = (fFilterBit == (1 << 7) || fFilterMask == 128)
 //                          ? labels[-1 - fEvent->GetTrack(i)->GetID()]
 //                          : i;
//    const auto *aodtrackpid = static_cast<AliAODTrack *>(fEvent->GetTrack(pid_track_id));
    //assert(aodtrackpid && "Not a standard AOD");

    float vertexX = -999.;
    float vertexY = -999.;
    float vertexZ = -999.;

    if (vertex) {
      Double_t fCov[6] = {0.0};
      vertex->GetCovarianceMatrix(fCov);
      if (fCov[5] != 0.0) {
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
  //  tFemtoTrack->SetCdd(covmat[0]);
  //  tFemtoTrack->SetCdz(covmat[1]);
  //  tFemtoTrack->SetCzz(covmat[2]);

  tFemtoTrack->SetTPCchi2(tAodTrack->GetTPCchi2());
  tFemtoTrack->SetTPCncls(tAodTrack->GetTPCNcls());
  tFemtoTrack->SetTPCnclsF(tAodTrack->GetTPCNclsF());

  tFemtoTrack->SetTPCsignal(tAodTrack->GetTPCsignal());
  //tFemtoTrack->SetTPCClusterMap(tAodTrack->GetTPCClusterMap());
  //tFemtoTrack->SetTPCSharedMap(tAodTrack->GetTPCSharedMap());

  float globalPositionsAtRadii[9][3];
  float bfield = 5 * fMagFieldSign;

  if(fCovMatPresent)
     GetGlobalPositionAtGlobalRadiiThroughTPC(tAodTrack, bfield, globalPositionsAtRadii);

  AliFemtoThreeVector tpcPositions[9];
  if(fCovMatPresent)
    std::copy_n(globalPositionsAtRadii, 9, tpcPositions);

  if (fPrimaryVertexCorrectionTPCPoints) {
    for (int i = 0; i < 9; i++) {
      tpcPositions[i] -= kOrigin;
    }
  }

  tFemtoTrack->SetNominalTPCEntrancePoint(tpcPositions[0]);
  tFemtoTrack->SetNominalTPCPoints(tpcPositions);
  tFemtoTrack->SetNominalTPCExitPoint(tpcPositions[8]);

  //if (fShiftPosition > 0.) {
  //Float_t posShifted[3];
  //SetShiftedPositions(tAodTrack, bfield, posShifted, fShiftPosition);
  //tFemtoTrack->SetNominalTPCPointShifted(posShifted);
  //}

  int kink_indexes[3] = { 0, 0, 0 };
  tFemtoTrack->SetKinkIndexes(kink_indexes);

//Corrections
  if (f1DcorrectionsPions) {
    tFemtoTrack->SetCorrectionPion(f1DcorrectionsPions->GetBinContent(f1DcorrectionsPions->FindFixBin(tAodTrack->Pt())));
    //cout<<"Setting pion correction to: "<<f1DcorrectionsPions->GetBinContent(f1DcorrectionsPions->FindFixBin(tAodTrack->Pt()))<<endl;
  }
  else if (f4DcorrectionsPions) {
   

    Int_t idx[4] = {f4DcorrectionsPions->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsPions->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsPions->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsPions->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    // cout<<"Track with pT "<<tAodTrack->Pt()<<" eta: "<<tAodTrack->Eta()<<" zv: "<<tAodTrack->Zv()<<" phi: "<<tAodTrack->Phi()<<endl;
    //cout<<"Pion bin: "<<idx[0]<<" "<<idx[1]<<" "<<idx[2]<<" "<<idx[3]<<" val: "<<f4DcorrectionsPions->GetBinContent(idx)<<endl;
    double correction = f4DcorrectionsPions->GetBinContent(idx);
    tFemtoTrack->SetCorrectionPion(correction == 0.0 ? 1.0 : correction);
    
  }
  else {
    tFemtoTrack->SetCorrectionPion(1.0);
  }

  if (f1DcorrectionsKaons) {
    tFemtoTrack->SetCorrectionKaon(f1DcorrectionsKaons->GetBinContent(f1DcorrectionsKaons->FindFixBin(tAodTrack->Pt())));
    //cout<<"Setting kaon correction to: "<<f1DcorrectionsKaons->GetBinContent(f1DcorrectionsKaons->FindFixBin(tAodTrack->Pt()))<<endl;
  }
  else if (f4DcorrectionsKaons) {
    Int_t idx[4] = {f4DcorrectionsKaons->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsKaons->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsKaons->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsKaons->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsKaons->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionKaon(f4DcorrectionsKaons->GetBinContent(idx));
    }else {
      tFemtoTrack->SetCorrectionKaon(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionKaon(1.0);
  }

  if (f1DcorrectionsProtons) {
    tFemtoTrack->SetCorrectionProton(f1DcorrectionsProtons->GetBinContent(f1DcorrectionsProtons->FindFixBin(tAodTrack->Pt())));
  }
  else if (f4DcorrectionsProtons) {
    Int_t idx[4] = {f4DcorrectionsProtons->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsProtons->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsProtons->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsProtons->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsProtons->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionProton(f4DcorrectionsProtons->GetBinContent(idx));
    } else {
      tFemtoTrack->SetCorrectionProton(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionProton(1.0);
  }

  if (f1DcorrectionsDeuterons) {
    tFemtoTrack->SetCorrectionDeuteron(f1DcorrectionsDeuterons->GetBinContent(f1DcorrectionsDeuterons->FindFixBin(tAodTrack->Pt())));
  }
  else if (f4DcorrectionsDeuterons) {
    Int_t idx[4] = {f4DcorrectionsDeuterons->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsDeuterons->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsDeuterons->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsDeuterons->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsDeuterons->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionDeuteron(f4DcorrectionsDeuterons->GetBinContent(idx));
    } else {
      tFemtoTrack->SetCorrectionDeuteron(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionDeuteron(1.0);
  }

  if (f1DcorrectionsPionsMinus) {
    tFemtoTrack->SetCorrectionPionMinus(f1DcorrectionsPionsMinus->GetBinContent(f1DcorrectionsPionsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else if (f4DcorrectionsPionsMinus) {
    Int_t idx[4] = {f4DcorrectionsPionsMinus->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsPionsMinus->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsPionsMinus->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsPionsMinus->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsPionsMinus->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionPionMinus(f4DcorrectionsPionsMinus->GetBinContent(idx));
    } else {
      tFemtoTrack->SetCorrectionPionMinus(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionPionMinus(1.0);
  }

  if (f1DcorrectionsKaonsMinus) {
    tFemtoTrack->SetCorrectionKaonMinus(f1DcorrectionsKaonsMinus->GetBinContent(f1DcorrectionsKaonsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else if (f4DcorrectionsKaonsMinus) {
    Int_t idx[4] = {f4DcorrectionsKaonsMinus->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsKaonsMinus->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsKaonsMinus->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsKaonsMinus->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsKaonsMinus->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionKaonMinus(f4DcorrectionsKaonsMinus->GetBinContent(idx));
    } else {
      tFemtoTrack->SetCorrectionKaonMinus(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionKaonMinus(1.0);
  }

  if (f1DcorrectionsProtonsMinus) {
    tFemtoTrack->SetCorrectionProtonMinus(f1DcorrectionsProtonsMinus->GetBinContent(f1DcorrectionsProtonsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else if (f4DcorrectionsProtonsMinus) {
    Int_t idx[4] = {f4DcorrectionsProtonsMinus->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsProtonsMinus->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsProtonsMinus->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsProtonsMinus->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsProtonsMinus->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionProtonMinus(f4DcorrectionsProtonsMinus->GetBinContent(idx));
    } else {
      tFemtoTrack->SetCorrectionProtonMinus(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionProtonMinus(1.0);
  }

  if (f1DcorrectionsDeuteronsMinus) {
    tFemtoTrack->SetCorrectionDeuteronMinus(f1DcorrectionsDeuteronsMinus->GetBinContent(f1DcorrectionsDeuteronsMinus->FindFixBin(tAodTrack->Pt())));
  }
  else if (f4DcorrectionsDeuteronsMinus) {
    Int_t idx[4] = {f4DcorrectionsDeuteronsMinus->GetAxis(0)->FindFixBin(tAodTrack->Eta()),
                    f4DcorrectionsDeuteronsMinus->GetAxis(1)->FindFixBin(tAodTrack->Pt()),
                    f4DcorrectionsDeuteronsMinus->GetAxis(2)->FindFixBin(0.0),
                    f4DcorrectionsDeuteronsMinus->GetAxis(3)->FindFixBin(tAodTrack->Phi())};
    if (f4DcorrectionsDeuteronsMinus->GetBinContent(idx) != 0) {
      tFemtoTrack->SetCorrectionDeuteronMinus(f4DcorrectionsDeuteronsMinus->GetBinContent(idx));
    } else {
      tFemtoTrack->SetCorrectionDeuteronMinus(1.0);
    }
  }
  else {
    tFemtoTrack->SetCorrectionDeuteronMinus(1.0);
  }

  if (f1DcorrectionsAll) {
    tFemtoTrack->SetCorrectionAll(f1DcorrectionsAll->GetBinContent(f1DcorrectionsAll->FindFixBin(tAodTrack->Pt())));
  }
  else {
    tFemtoTrack->SetCorrectionAll(1.0);
  }
  

  //
  /*******************************************************************/
  return tFemtoTrack;
}

AliFemtoV0 *AliFemtoEventReaderNanoAOD::CopyAODtoFemtoV0(AliAODv0 *tAODv0)
{
  AliFemtoV0 *tFemtoV0 = new AliFemtoV0();
  tFemtoV0->SetdecayLengthV0(tAODv0->DecayLength(fVe1));
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
  tFemtoV0->SetprimaryVertex(fVe1);
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
  tFemtoV0->SetptV0(tAODv0->Pt());
  tFemtoV0->SetptotV0(::sqrt(tAODv0->Ptot2V0()));
  tFemtoV0->SetidNeg(tAODv0->GetNegID());
  tFemtoV0->SetidPos(tAODv0->GetPosID());
  tFemtoV0->SetEtaV0(tAODv0->Eta());
  tFemtoV0->SetPhiV0(tAODv0->Phi());
  tFemtoV0->SetCosPointingAngle(tAODv0->CosPointingAngle(fVe1));

  AliNanoAODTrack *trackpos = (AliNanoAODTrack *)tAODv0->GetDaughter(0);
  AliNanoAODTrack *trackneg = (AliNanoAODTrack *)tAODv0->GetDaughter(1);

  // ensure that trackpos and trackneg are pointing to the correct children
  // This confusion seems to arise when fOnFlyStatusV0 = true
  if (trackpos->Charge() < 0 && trackneg->Charge() > 0) {
    AliNanoAODTrack *tmp = trackpos;
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
     tFemtoV0->SetTPCNclsPos(trackpos->GetTPCNcls());
     tFemtoV0->SetTPCNclsNeg(trackneg->GetTPCNcls());
     tFemtoV0->SetNdofPos(trackpos->Chi2perNDF());
     tFemtoV0->SetNdofNeg(trackneg->Chi2perNDF());

     static const Int_t kcstNSigmaTPCPi  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kPion);
     static const Int_t kcstNSigmaTPCK  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kKaon);
     static const Int_t kcstNSigmaTPCPr = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kProton);

     const float nsigmaTPCPipos = trackpos->GetVar(kcstNSigmaTPCPi);
     const float nsigmaTPCPineg = trackneg->GetVar(kcstNSigmaTPCPi);
     const float nsigmaTPCKpos = trackpos->GetVar(kcstNSigmaTPCK);
     const float nsigmaTPCKneg = trackneg->GetVar(kcstNSigmaTPCK);
     const float nsigmaTPCPpos = trackpos->GetVar(kcstNSigmaTPCPr);
     const float nsigmaTPCPneg = trackneg->GetVar(kcstNSigmaTPCPr);

     tFemtoV0->SetPosNSigmaTPCPi(nsigmaTPCPipos);
     tFemtoV0->SetNegNSigmaTPCPi(nsigmaTPCPineg);
     tFemtoV0->SetPosNSigmaTPCK(nsigmaTPCKpos);
     tFemtoV0->SetNegNSigmaTPCK(nsigmaTPCKneg);
     tFemtoV0->SetPosNSigmaTPCP(nsigmaTPCPpos);
     tFemtoV0->SetNegNSigmaTPCP(nsigmaTPCPneg);


     float bfield = 5 * fMagFieldSign;
     float globalPositionsAtRadiiPos[9][3];
     for(int i=0; i<9; i++)
        {
         for(int j=0; j<3; j++)
          {
            globalPositionsAtRadiiPos[i][j] = 0;
          }
        }
     if(fCovMatPresent)
       GetGlobalPositionAtGlobalRadiiThroughTPC(trackpos, bfield, globalPositionsAtRadiiPos);
     double tpcEntrancePos[3] = {globalPositionsAtRadiiPos[0][0], globalPositionsAtRadiiPos[0][1], globalPositionsAtRadiiPos[0][2]};
     double tpcExitPos[3] = {globalPositionsAtRadiiPos[8][0], globalPositionsAtRadiiPos[8][1], globalPositionsAtRadiiPos[8][2]};

     float globalPositionsAtRadiiNeg[9][3];
     for(int i=0; i<9; i++)
        {
        for(int j=0; j<3; j++)
          {
            globalPositionsAtRadiiNeg[i][j] = 0;
          }
        }
     if(fCovMatPresent)
       GetGlobalPositionAtGlobalRadiiThroughTPC(trackneg, bfield, globalPositionsAtRadiiNeg);
     double tpcEntranceNeg[3] = {globalPositionsAtRadiiNeg[0][0], globalPositionsAtRadiiNeg[0][1], globalPositionsAtRadiiNeg[0][2]};
     double tpcExitNeg[3] = {globalPositionsAtRadiiNeg[8][0], globalPositionsAtRadiiNeg[8][1], globalPositionsAtRadiiNeg[8][2]};

     if (fPrimaryVertexCorrectionTPCPoints) {
        tpcEntrancePos[0] -= fVe1[0];
        tpcEntrancePos[1] -= fVe1[1];
        tpcEntrancePos[2] -= fVe1[2];

        tpcExitPos[0] -= fVe1[0];
        tpcExitPos[1] -= fVe1[1];
        tpcExitPos[2] -= fVe1[2];

        tpcEntranceNeg[0] -= fVe1[0];
        tpcEntranceNeg[1] -= fVe1[1];
        tpcEntranceNeg[2] -= fVe1[2];

        tpcExitNeg[0] -= fVe1[0];
        tpcExitNeg[1] -= fVe1[1];
        tpcExitNeg[2] -= fVe1[2];
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

  tFemtoV0->SetNominalTpcPointPos(vecTpcPos);
  tFemtoV0->SetNominalTpcPointNeg(vecTpcNeg);

   tFemtoV0->SetTPCMomentumPos(trackpos->GetTPCmomentum());
   tFemtoV0->SetTPCMomentumNeg(trackneg->GetTPCmomentum());
  // tFemtoV0->SetdedxPos(trackpos->GetTPCsignal());
  // tFemtoV0->SetdedxNeg(trackneg->GetTPCsignal());

  static const Int_t kcstNSigmaTOFPi  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kPion);
 	static const Int_t kcstNSigmaTOFK  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kKaon);
 	static const Int_t kcstNSigmaTOFPr  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kProton);

  const float nsigmaTOFPipos = trackpos->GetVar(kcstNSigmaTOFPi);
  const float nsigmaTOFKpos = trackpos->GetVar(kcstNSigmaTOFK);
  const float nsigmaTOFPpos = trackpos->GetVar(kcstNSigmaTOFPr);
  const float nsigmaTOFPineg = trackneg->GetVar(kcstNSigmaTOFPi);
  const float nsigmaTOFKneg = trackneg->GetVar(kcstNSigmaTOFK);
  const float nsigmaTOFPneg = trackneg->GetVar(kcstNSigmaTOFPr);

  tFemtoV0->SetPosNSigmaTOFPi(nsigmaTOFPipos);
  tFemtoV0->SetPosNSigmaTOFK(nsigmaTOFKpos);
  tFemtoV0->SetPosNSigmaTOFP(nsigmaTOFPpos);
  tFemtoV0->SetNegNSigmaTOFPi(nsigmaTOFPineg);
  tFemtoV0->SetNegNSigmaTOFK(nsigmaTOFKneg);
  tFemtoV0->SetNegNSigmaTOFP(nsigmaTOFPneg);


  /*double TOFSignalPos = trackpos->GetTOFsignal();
  double TOFSignalNeg = trackneg->GetTOFsignal();
  tFemtoV0->SetTOFPionTimePos(TOFSignalPos);
  tFemtoV0->SetTOFKaonTimePos(TOFSignalPos);
  tFemtoV0->SetTOFProtonTimePos(TOFSignalPos);
  tFemtoV0->SetTOFPionTimeNeg(TOFSignalNeg);
  tFemtoV0->SetTOFKaonTimeNeg(TOFSignalNeg);
  tFemtoV0->SetTOFProtonTimeNeg(TOFSignalNeg);
*/
   }

  tFemtoV0->SetOnFlyStatusV0(tAODv0->GetOnFlyStatus());


  //corrections
  if (f1DcorrectionsLambdas) {
    tFemtoV0->SetCorrectionLambdas(f1DcorrectionsLambdas->GetBinContent(f1DcorrectionsLambdas->FindFixBin(tAODv0->Pt())));
  }
  else if (f4DcorrectionsLambdas) {
    Int_t idx[4] = {f4DcorrectionsLambdas->GetAxis(0)->FindFixBin(tAODv0->Eta()),
                    f4DcorrectionsLambdas->GetAxis(1)->FindFixBin(tAODv0->Pt()),
                    f4DcorrectionsLambdas->GetAxis(2)->FindFixBin(tAODv0->Zv()),
                    f4DcorrectionsLambdas->GetAxis(3)->FindFixBin(tAODv0->Phi())};
    if (f4DcorrectionsLambdas->GetBinContent(idx) != 0) {
      tFemtoV0->SetCorrectionLambdas(1. / f4DcorrectionsLambdas->GetBinContent(idx));
    } else {
      tFemtoV0->SetCorrectionLambdas(1.0);
    }
  }
  else {
    tFemtoV0->SetCorrectionLambdas(1.0);
  }

  if (f1DcorrectionsLambdasMinus) {
    tFemtoV0->SetCorrectionLambdasMinus(f1DcorrectionsLambdasMinus->GetBinContent(f1DcorrectionsLambdasMinus->FindFixBin(tAODv0->Pt())));
  }
   else if (f4DcorrectionsLambdasMinus) {
    Int_t idx[4] = {f4DcorrectionsLambdasMinus->GetAxis(0)->FindFixBin(tAODv0->Eta()),
                    f4DcorrectionsLambdasMinus->GetAxis(1)->FindFixBin(tAODv0->Pt()),
                    f4DcorrectionsLambdasMinus->GetAxis(2)->FindFixBin(tAODv0->Zv()),
                    f4DcorrectionsLambdasMinus->GetAxis(3)->FindFixBin(tAODv0->Phi())};
    if (f4DcorrectionsLambdasMinus->GetBinContent(idx) != 0) {
      tFemtoV0->SetCorrectionLambdasMinus(1. / f4DcorrectionsLambdasMinus->GetBinContent(idx));
    } else {
      tFemtoV0->SetCorrectionLambdasMinus(1.0);
    }
  }
  else {
    tFemtoV0->SetCorrectionLambdasMinus(1.0);
  }
  //*****************

  return tFemtoV0;
}


AliFemtoXi *AliFemtoEventReaderNanoAOD::CopyAODtoFemtoXi(AliAODcascade *tAODxi)
{
  AliFemtoXi *tFemtoXi = nullptr;

  { // this is to keep tmpV0 in its own scope
    AliFemtoV0 *tmpV0 = CopyAODtoFemtoV0(tAODxi);
    tFemtoXi = new AliFemtoXi(tmpV0);
    delete tmpV0;
  }

  //The above lines set the V0 attributes.  However, some of these are now set incorrectly
  //For instance:
  //  tFemtoV0->SetdecayLengthV0(tAODv0->DecayLength(fVe1)); in CopyAODtoFemtoXi now sets the decay length
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
  tFemtoXi->SetmassOmega(tAODxi->MassOmega());
  tFemtoXi->SetdecayLengthXi(tAODxi->DecayLengthXi(fVe1[0], fVe1[1], fVe1[2]));
  tFemtoXi->SetdecayVertexXiX(tAODxi->DecayVertexXiX());
  tFemtoXi->SetdecayVertexXiY(tAODxi->DecayVertexXiY());
  tFemtoXi->SetdecayVertexXiZ(tAODxi->DecayVertexXiZ());
  tFemtoXi->SetdcaXiDaughters(tAODxi->DcaXiDaughters());
  //tFemtoXi->SetdcaXiToPrimVertex(tAODxi->DcaXiToPrimVertex());  //This doesn't work because fDcaXiToPrimVertex was not set
  tFemtoXi->SetdcaXiToPrimVertex(tAODxi->DcaXiToPrimVertex(fVe1[0], fVe1[1], fVe1[2]));
  tFemtoXi->SetdcaBacToPrimVertex(tAODxi->DcaBachToPrimVertex());

  tFemtoXi->SetmomBacX(tAODxi->MomBachX());
  tFemtoXi->SetmomBacY(tAODxi->MomBachY());
  tFemtoXi->SetmomBacZ(tAODxi->MomBachZ());
  tFemtoXi->SetmomXiX(tAODxi->MomXiX());
  tFemtoXi->SetmomXiY(tAODxi->MomXiY());
  tFemtoXi->SetmomXiZ(tAODxi->MomXiZ());
  AliFemtoThreeVector momxi(tAODxi->MomXiX(), tAODxi->MomXiY(), tAODxi->MomXiZ());
  tFemtoXi->SetmomXi(momxi);
  tFemtoXi->SetRadiusXi(TMath::Sqrt(tAODxi->DecayVertexXiX() * tAODxi->DecayVertexXiX()
                                    + tAODxi->DecayVertexXiY() * tAODxi->DecayVertexXiY()));

  tFemtoXi->SetidBac(tAODxi->GetBachID());

  double tEtaXi = 0.5 * TMath::Log((TMath::Sqrt(tAODxi->Ptot2Xi()) + tAODxi->MomXiZ()) / (TMath::Sqrt(tAODxi->Ptot2Xi()) - tAODxi->MomXiZ() + 1.e-13));
  tFemtoXi->SetEtaXi(tEtaXi);
  double tPhiXi = TMath::Pi() + TMath::ATan2(-tAODxi->MomXiY(), -tAODxi->MomXiX());
  tFemtoXi->SetPhiXi(tPhiXi);

  tFemtoXi->SetCosPointingAngleXi(tAODxi->CosPointingAngleXi(fVe1[0], fVe1[1], fVe1[2]));
  tFemtoXi->SetCosPointingAngleV0toXi(tAODxi->CosPointingAngle(tAODxi->GetDecayVertexXi()));

  tFemtoXi->SetChargeXi(tAODxi->ChargeXi());
  tFemtoXi->SetptXi(std::sqrt(tAODxi->Pt2Xi()));

  if (f1DcorrectionsXiPlus) {
    tFemtoXi->SetCorrectionXiPlus(f1DcorrectionsXiPlus->GetBinContent(f1DcorrectionsXiPlus->FindFixBin(tAODxi->Pt())));
  }
  else {
    tFemtoXi->SetCorrectionXiPlus(1.0);
  }

  if (f1DcorrectionsXiMinus) {
    tFemtoXi->SetCorrectionXiMinus(f1DcorrectionsXiMinus->GetBinContent(f1DcorrectionsXiMinus->FindFixBin(tAODxi->Pt())));
  }
  else {
    tFemtoXi->SetCorrectionXiMinus(1.0);
  }

  AliNanoAODTrack *trackbac = (AliNanoAODTrack *)tAODxi->GetDecayVertexXi()->GetDaughter(0);

  if (trackbac) {
    tFemtoXi->SetptBac(trackbac->Pt()); //setting pt? px and py was set!
    tFemtoXi->SetEtaBac(trackbac->Eta());            //bac!
    tFemtoXi->SetTPCNclsBac(trackbac->GetTPCNcls()); //bac!
    tFemtoXi->SetNdofBac(trackbac->Chi2perNDF());    //bac!
  //  tFemtoXi->SetStatusBac(trackbac->GetStatus());   //bac!

    static const Int_t kcstNSigmaTPCPi  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kPion);
    static const Int_t kcstNSigmaTPCK  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kKaon);
    static const Int_t kcstNSigmaTPCPr = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, AliPID::kProton);

    const float nsigmaTPCPibac = trackbac->GetVar(kcstNSigmaTPCPi);
    const float nsigmaTPCKbac = trackbac->GetVar(kcstNSigmaTPCK);
    const float nsigmaTPCPbac = trackbac->GetVar(kcstNSigmaTPCPr);

    tFemtoXi->SetBacNSigmaTPCK(nsigmaTPCKbac);
    tFemtoXi->SetBacNSigmaTPCP(nsigmaTPCPbac);
    tFemtoXi->SetBacNSigmaTPCPi(nsigmaTPCPibac);

    //NEED TO ADD:
    //    tFemtoXi->SetNominalTpcEntrancePointBac(tmpVec);
    //    tFemtoXi->SetNominalTpcExitPointBac(tmpVec);
    //    tFemtoXi->SetNominalTpcPointBac(vecTpcPos);
    //    tFemtoXi->SetTPCMomentumBac(trackbac->GetTPCmomentum());
    //    if (fShiftPosition > 0.)

    float bfield = 5 * fMagFieldSign;
    float globalPositionsAtRadiiBac[9][3];
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<3; j++)
      {
        globalPositionsAtRadiiBac[i][j] = 0;
      }
    }
    if(fCovMatPresent)
      GetGlobalPositionAtGlobalRadiiThroughTPC(trackbac, bfield, globalPositionsAtRadiiBac);
    double tpcEntranceBac[3] = {globalPositionsAtRadiiBac[0][0], globalPositionsAtRadiiBac[0][1], globalPositionsAtRadiiBac[0][2]};
    double tpcExitBac[3] = {globalPositionsAtRadiiBac[8][0], globalPositionsAtRadiiBac[8][1], globalPositionsAtRadiiBac[8][2]};

   if (fPrimaryVertexCorrectionTPCPoints) {
      tpcEntranceBac[0] -= fVe1[0];
      tpcEntranceBac[1] -= fVe1[1];
      tpcEntranceBac[2] -= fVe1[2];

      tpcExitBac[0] -= fVe1[0];
      tpcExitBac[1] -= fVe1[1];
      tpcExitBac[2] -= fVe1[2];
    }

    AliFemtoThreeVector tmpVec;
    tmpVec.SetX(tpcEntranceBac[0]);
    tmpVec.SetY(tpcEntranceBac[1]);
    tmpVec.SetZ(tpcEntranceBac[2]);
    tFemtoXi->SetNominalTpcEntrancePointBac(tmpVec);

    tmpVec.SetX(tpcExitBac[0]);
    tmpVec.SetY(tpcExitBac[1]);
    tmpVec.SetZ(tpcExitBac[2]);
    tFemtoXi->SetNominalTpcExitPointBac(tmpVec);

    AliFemtoThreeVector vecTpcBac[9];
    //cout<<"La"<<endl;
  //  for (int i = 0; i < 9; i++) {
      //cout<<globalPositionsAtRadiiBac[i][0]<<" "<<globalPositionsAtRadiiBac[i][1]<<" "<<globalPositionsAtRadiiBac[i][2]<<endl;
      //vecTpcBac[i].SetX(globalPositionsAtRadiiBac[i][0]);
      //vecTpcBac[i].SetY(globalPositionsAtRadiiBac[i][1]);
      //vecTpcBac[i].SetZ(globalPositionsAtRadiiBac[i][2]);
  //  }

    if (fPrimaryVertexCorrectionTPCPoints) {
      AliFemtoThreeVector tmpVertexVec;
      tmpVertexVec.SetX(fVe1[0]);
      tmpVertexVec.SetY(fVe1[1]);
      tmpVertexVec.SetZ(fVe1[2]);

      for (int i = 0; i < 9; i++) {
        vecTpcBac[i] -= tmpVertexVec;
      }
    }

    //tFemtoXi->SetNominalTpcPointBac(vecTpcBac);
    tFemtoXi->SetTPCMomentumBac(trackbac->GetTPCmomentum());
    //tFemtoXi->SetdedxBac(trackbac->GetTPCsignal());

    static const Int_t kcstNSigmaTOFPi  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kPion);
    static const Int_t kcstNSigmaTOFK  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kKaon);
    static const Int_t kcstNSigmaTOFPr  = AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, AliPID::kProton);

    const float nsigmaTOFPibac = trackbac->GetVar(kcstNSigmaTOFPi);
    const float nsigmaTOFKbac = trackbac->GetVar(kcstNSigmaTOFK);
    const float nsigmaTOFPbac = trackbac->GetVar(kcstNSigmaTOFPr);

    tFemtoXi->SetBacNSigmaTOFK(nsigmaTOFKbac);
    tFemtoXi->SetBacNSigmaTOFP(nsigmaTOFPbac);
    tFemtoXi->SetBacNSigmaTOFPi(nsigmaTOFPibac);

  /*  double TOFSignalBac = trackbac->GetTOFsignal();  //tego nie było, moze nie dzialac
    tFemtoXi->SetTOFPionTimeBac(TOFSignalBac);
    tFemtoXi->SetTOFKaonTimeBac(TOFSignalBac);
    tFemtoXi->SetTOFProtonTimeBac(TOFSignalBac);*/
  }
   return tFemtoXi;
}

void AliFemtoEventReaderNanoAOD::SetFilterBit(UInt_t ibit)
{
  fFilterBit = (1 << (ibit));
}

void AliFemtoEventReaderNanoAOD::SetFilterMask(int ibit)
{
  fFilterMask = ibit;
}

void AliFemtoEventReaderNanoAOD::SetReadMC(unsigned char a)
{
  fReadMC = a;
}

void AliFemtoEventReaderNanoAOD::SetReadV0(unsigned char a)
{
  fReadV0 = a;
}

void AliFemtoEventReaderNanoAOD::SetReadCascade(unsigned char a)
{
  fReadCascade = a;
}

void AliFemtoEventReaderNanoAOD::SetUseMultiplicity(string aType)
{
  fEstEventMult = aType;
}



void AliFemtoEventReaderNanoAOD::SetAODpidUtil(AliAODpidUtil *aAODpidUtil)
{
  fAODpidUtil = aAODpidUtil;
  //  printf("fAODpidUtil: %x\n",fAODpidUtil);
}

void AliFemtoEventReaderNanoAOD::SetAODheader(AliNanoAODHeader *aAODheader)
{
  fAODheader = aAODheader;
}

void AliFemtoEventReaderNanoAOD::SetMagneticFieldSign(int s)
{
  if (s > 0)
    fMagFieldSign = 1;
  else if (s < 0)
    fMagFieldSign = -1;
  else
    fMagFieldSign = 0;
}

void AliFemtoEventReaderNanoAOD::SetEPVZERO(Bool_t iepvz)
{
  fisEPVZ = iepvz;
}

void AliFemtoEventReaderNanoAOD
  ::GetGlobalPositionAtGlobalRadiiThroughTPC(AliNanoAODTrack *track,
                                             Float_t bfield,
                                             Float_t globalPositionsAtRadii[9][3])
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
    Double_t pos_buffer[3] = {0.0, 0.0, 0.0};

    // get the global position of the track at this radial location
    bool good = etp.GetXYZatR(radius, bfield, pos_buffer, nullptr);

    // if value is not good, break loading loop
    if (!good || fabs(AliFemtoThreeVector(pos_buffer).Perp() - radius) > 0.5) {
      radius_index--; // decrement to fill current location with default value
      break;
    }

    // store the global position
    std::copy_n(pos_buffer, 3, globalPositionsAtRadii[radius_index]);
  }

  // Fill any remaining positions with the default value
  for (; radius_index < 9; radius_index++) {
    std::fill_n(globalPositionsAtRadii[radius_index], 3, DEFAULT_VALUE);
  }
}

//________________________________________________________________________
/*void AliFemtoEventReaderNanoAOD::SetShiftedPositions(const AliNanoAODTrack *track,
                                                 const Float_t bfield,
                                                 Float_t posShifted[3],
                                                 const Double_t radius)
{
  // Sets the spatial position of the track at the radius R=1.25m in
  // the shifted coordinate system, code adapted from Hans Beck analysis

  // Initialize the array to something indicating there was no propagation
  posShifted[0] = -9999.; // THIS IS THE DATA MEMBER OF YOUR FEMTOTRACK
  posShifted[1] = -9999.;
  posShifted[2] = -9999.;

  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);

  // The global position of the the track
  Double_t xyz[3] = {-9999., -9999., -9999.};

  // The radius in cm we want to propagate to, squared
  const Float_t RSquaredWanted(radius * radius * 1e4);

  // Propagation is done in local x of the track
  for (Float_t x = 58.; x < 247.; x += 1.) {
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent.
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if (!etp.PropagateTo(x, bfield)) {
      break;
    }
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared.
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0] - fVe1[0]) * (xyz[0] - fVe1[0])
                                 + (xyz[1] - fVe1[1]) * (xyz[1] - fVe1[1]);

    // Roughly reached the radius we want
    if (shiftedRadiusSquared > RSquaredWanted) {

      // Bigger loop has bad precision, we're nearly one centimeter too far,
      // go back in small steps.
      while (shiftedRadiusSquared > RSquaredWanted) {
        // Propagate a mm inwards
        x -= 0.1;
        if (!etp.PropagateTo(x, bfield)) {
          // Propagation failed but we're already with a
          // cm precision at R=1.25m so we only break the
          // inner loop
          break;
        }
        // Get the global position
        etp.GetXYZ(xyz);
        // Calculate shifted radius, squared
        shiftedRadiusSquared = (xyz[0] - fVe1[0]) * (xyz[0] - fVe1[0])
                             + (xyz[1] - fVe1[1]) * (xyz[1] - fVe1[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      posShifted[0] = xyz[0] - fVe1[0];
      posShifted[1] = xyz[1] - fVe1[1];
      posShifted[2] = xyz[2] - fVe1[2];
      // Done
      return;
    } // End of if roughly reached radius
  }   // End of coarse propagation loop
}
*/
void AliFemtoEventReaderNanoAOD::SetUseAliEventCuts(Bool_t useAliEventCuts)
{
  fUseAliEventCuts = useAliEventCuts;
  fEventCuts = new AliEventCuts();
}

void AliFemtoEventReaderNanoAOD::SetpA2013(Bool_t pa2013)
{
  fpA2013 = pa2013;
}

void AliFemtoEventReaderNanoAOD::SetUseMVPlpSelection(Bool_t mvplp)
{
  fMVPlp = mvplp;
}

void AliFemtoEventReaderNanoAOD::SetUseOutOfBunchPlpSelection(Bool_t outOfBunchPlp)
{
  fOutOfBunchPlp = outOfBunchPlp;
}

void AliFemtoEventReaderNanoAOD::SetIsPileUpEvent(Bool_t ispileup)
{
  fisPileUp = ispileup;
}

void AliFemtoEventReaderNanoAOD::SetCascadePileUpRemoval(Bool_t cascadePileUpRemoval)
{
  fCascadePileUpRemoval = cascadePileUpRemoval;
}

void AliFemtoEventReaderNanoAOD::SetV0PileUpRemoval(Bool_t v0PileUpRemoval)
{
  fV0PileUpRemoval = v0PileUpRemoval;
}

void AliFemtoEventReaderNanoAOD::SetTrackPileUpRemoval(Bool_t trackPileUpRemoval)
{
  fTrackPileUpRemoval = trackPileUpRemoval;
}

void AliFemtoEventReaderNanoAOD::SetDCAglobalTrack(Int_t dcagt)
{
  fDCAglobalTrack = dcagt;
}

bool AliFemtoEventReaderNanoAOD::RejectEventCentFlat(float MagField, float CentPercent)
{
  // Flattens the centrality distribution

  // Setting 0 as seed ensures random seed every time.
  TRandom3 RNG(0); // for 3D, random sign switching

  float kCentWeight[2][9] = {
      {0.878, .876, .860, .859, .859, .880, .873, .879, .894},
      {0.828, .793, .776, .772, .775, .796, .788, .804, .839}};

  int weightBinCent = (int)CentPercent,
      weightBinSign = (MagField > 0) ? 0 : 1;

  bool rejectEvent = RNG.Rndm() > kCentWeight[weightBinSign][weightBinCent];
  return rejectEvent;
}

void AliFemtoEventReaderNanoAOD::SetCentralityFlattening(Bool_t dcagt)
{
  fFlatCent = dcagt;
}

void AliFemtoEventReaderNanoAOD::SetShiftPosition(Double_t dcagt)
{
  fShiftPosition = dcagt;
}

void AliFemtoEventReaderNanoAOD::SetPrimaryVertexCorrectionTPCPoints(bool correctTpcPoints)
{
  fPrimaryVertexCorrectionTPCPoints = correctTpcPoints;
}


void AliFemtoEventReaderNanoAOD::Set1DCorrectionsPions(TH1D *h1)
{
  f1DcorrectionsPions = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsKaons(TH1D *h1)
{
  f1DcorrectionsKaons = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsProtons(TH1D *h1)
{
  f1DcorrectionsProtons = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsDeuterons(TH1D *h1)
{
  f1DcorrectionsDeuterons = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsPionsMinus(TH1D *h1)
{
  f1DcorrectionsPionsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsKaonsMinus(TH1D *h1)
{
  f1DcorrectionsKaonsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsProtonsMinus(TH1D *h1)
{
  f1DcorrectionsProtonsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsDeuteronsMinus(TH1D *h1)
{
  f1DcorrectionsDeuteronsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsAll(TH1D *h1)
{
  f1DcorrectionsAll = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsLambdas(TH1D *h1)
{
  f1DcorrectionsLambdas = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsLambdasMinus(TH1D *h1)
{
  f1DcorrectionsLambdasMinus = h1;
}
void AliFemtoEventReaderNanoAOD::Set1DCorrectionsXiPlus(TH1D *h1)
{
  f1DcorrectionsXiPlus = h1;
}

void AliFemtoEventReaderNanoAOD::Set1DCorrectionsXiMinus(TH1D *h1)
{
  f1DcorrectionsXiMinus = h1;
}


void AliFemtoEventReaderNanoAOD::Set4DCorrectionsPions(THnSparse *h1)
{
  f4DcorrectionsPions = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsKaons(THnSparse *h1)
{
  f4DcorrectionsKaons = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsProtons(THnSparse *h1)
{
  f4DcorrectionsProtons = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsDeuterons(THnSparse *h1)
{
  f4DcorrectionsDeuterons = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsPionsMinus(THnSparse *h1)
{
  f4DcorrectionsPionsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsKaonsMinus(THnSparse *h1)
{
  f4DcorrectionsKaonsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsProtonsMinus(THnSparse *h1)
{
  f4DcorrectionsProtonsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsDeuteronsMinus(THnSparse *h1)
{
  f4DcorrectionsDeuteronsMinus = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsAll(THnSparse *h1)
{
  f4DcorrectionsAll = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsLambdas(THnSparse *h1)
{
  f4DcorrectionsLambdas = h1;
}

void AliFemtoEventReaderNanoAOD::Set4DCorrectionsLambdasMinus(THnSparse *h1)
{
  f4DcorrectionsLambdasMinus = h1;
}
