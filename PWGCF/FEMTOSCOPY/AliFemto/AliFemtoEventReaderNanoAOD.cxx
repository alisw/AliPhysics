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
  fShiftPosition(0.)
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
  fShiftPosition(aReader.fShiftPosition)
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



  AliNanoAODHeader *header = dynamic_cast<AliNanoAODHeader *>(fEvent->GetHeader());
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

    string multString = "MultSelection."+fEstEventMult;
    int nMultString = multString.length(); 
    char multChar[nMultString + 1]; 
    strcpy(multChar, multString.c_str()); 
    
    norm_mult  = header->GetVarIndex(multChar);
    //norm_mult  = header->GetVarIndex("MultSelection.RefMult08");
    cout<<"header: "<<header<<endl;
    cout<<"multChar: "<<multChar<<endl;
    cout<<"norm_mult: "<<norm_mult<<endl;
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

	
	const float nsigmaTPCK = aodtrack->GetVar(kcstNSigmaTPCK); 
	const float nsigmaTPCPi = aodtrack->GetVar(kcstNSigmaTPCPi); 
	const float nsigmaTPCP = aodtrack->GetVar(kcstNSigmaTPCPr);
	const float nsigmaTPCE = aodtrack->GetVar(kcstNSigmaTPCE);


	trackCopy->SetNSigmaTPCPi(nsigmaTPCPi);
	trackCopy->SetNSigmaTPCK(nsigmaTPCK);
	trackCopy->SetNSigmaTPCP(nsigmaTPCP);
	trackCopy->SetNSigmaTPCE(nsigmaTPCE);

	trackCopy->SetTPCsignal(aodtrack->GetTPCsignal());
	trackCopy->SetTPCsignalS(1);
	trackCopy->SetTPCsignalN(aodtrack->GetTPCsignalN());


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

	
	const float nsigmaTOFK = aodtrack->GetVar(kcstNSigmaTOFK); 
	const float nsigmaTOFPi = aodtrack->GetVar(kcstNSigmaTOFPi); 
	const float nsigmaTOFP = aodtrack->GetVar(kcstNSigmaTOFPr);
	const float nsigmaTOFE = aodtrack->GetVar(kcstNSigmaTOFE);

	trackCopy->SetNSigmaTOFPi(nsigmaTOFPi);
	trackCopy->SetNSigmaTOFK(nsigmaTOFK);
	trackCopy->SetNSigmaTOFP(nsigmaTOFP);
	trackCopy->SetNSigmaTOFE(nsigmaTOFE);

	trackCopy->SetTOFsignal(aodtrack->GetTOFsignal());
	
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

/*
  if (fReadV0) {
    int count_pass = 0;
    for (Int_t i = 0; i < fEvent->GetNumberOfV0s(); i++) {
      AliAODv0 *aodv0 = fEvent->GetV0(i);
      // ensure a "good" v0 particle passes these conditions
      if (!aodv0
          || aodv0->GetNDaughters() > 2
          || aodv0->GetNProngs() > 2
          || aodv0->GetCharge() != 0
          || aodv0->ChargeProng(0) == aodv0->ChargeProng(1)
          || aodv0->CosPointingAngle(fVe1) < 0.98) {
        continue;
      }

      AliAODTrack *daughterTrackPos = (AliAODTrack *)aodv0->GetDaughter(0), // getting positive daughter track
                  *daughterTrackNeg = (AliAODTrack *)aodv0->GetDaughter(1); // getting negative daughter track

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
*/
/*
  if (fReadCascade) {
    int count_pass = 0;
    for (Int_t i = 0; i < fEvent->GetNumberOfCascades(); i++) {
      AliAODcascade *aodxi = fEvent->GetCascade(i);
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
*/
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
  Double_t covmat[21];
  tAodTrack->GetCovarianceXYZPxPyPz(covmat);

  if (fDCAglobalTrack == 0) {
    tFemtoTrack->SetImpactD(tAodTrack->DCA());
    tFemtoTrack->SetImpactZ(tAodTrack->ZAtDCA());

    tFemtoTrack->SetXatDCA(tAodTrack->XAtDCA());
    tFemtoTrack->SetYatDCA(tAodTrack->YAtDCA());
    tFemtoTrack->SetZatDCA(tAodTrack->ZAtDCA());
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

  GetGlobalPositionAtGlobalRadiiThroughTPC(tAodTrack, bfield, globalPositionsAtRadii);

  AliFemtoThreeVector tpcPositions[9];
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
  tFemtoV0->SetCosPointingAngle(tAODv0->CosPointingAngle(fVe1));
  //tFemtoV0->SetYV0(tAODv0->Y());

  //void SetdedxNeg(float x);
  //void SeterrdedxNeg(float x);//Gael 04Fev2002
  //void SetlendedxNeg(float x);//Gael 04Fev2002
  //void SetdedxPos(float x);
  //void SeterrdedxPos(float x);//Gael 04Fev2002
  //void SetlendedxPos(float x);//Gael 04Fev2002

  //tFemtoV0->SetEtaPos(tAODv0->PseudoRapPos());
  //tFemtoV0->SetEtaNeg(tAODv0->PseudoRapNeg());

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

    //tFemtoV0->SetEtaPos(trackpos->Eta()); //tAODv0->PseudoRapPos()
    //tFemtoV0->SetEtaNeg(trackneg->Eta()); //tAODv0->PseudoRapNeg()
    tFemtoV0->SetTPCNclsPos(trackpos->GetTPCNcls());
    tFemtoV0->SetTPCNclsNeg(trackneg->GetTPCNcls());
    //tFemtoV0->SetTPCclustersPos(trackpos->GetTPCClusterMap());
    //tFemtoV0->SetTPCclustersNeg(trackneg->GetTPCClusterMap());
    //tFemtoV0->SetTPCsharingPos(trackpos->GetTPCSharedMap());
    //tFemtoV0->SetTPCsharingNeg(trackneg->GetTPCSharedMap());
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

    if (fPrimaryVertexCorrectionTPCPoints) {
      AliFemtoThreeVector tmpVertexVec;
      tmpVertexVec.SetX(fVe1[0]);
      tmpVertexVec.SetY(fVe1[1]);
      tmpVertexVec.SetZ(fVe1[2]);

      for (int i = 0; i < 9; i++) {
        vecTpcPos[i] -= tmpVertexVec;
        vecTpcNeg[i] -= tmpVertexVec;
      }
    }

    tFemtoV0->SetNominalTpcPointPos(vecTpcPos);
    tFemtoV0->SetNominalTpcPointNeg(vecTpcNeg);

    /*
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
    */
    tFemtoV0->SetTPCMomentumPos(trackpos->GetTPCmomentum());
    tFemtoV0->SetTPCMomentumNeg(trackneg->GetTPCmomentum());

    tFemtoV0->SetdedxPos(trackpos->GetTPCsignal());
    tFemtoV0->SetdedxNeg(trackneg->GetTPCsignal());

    Float_t probMisPos = 1.0;
    Float_t probMisNeg = 1.0;

    if (((tFemtoV0->StatusPos() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
        && ((tFemtoV0->StatusPos() & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      // if (tFemtoV0->StatusPos() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000
      probMisPos = fAODpidUtil->GetTOFMismatchProbability(trackpos);
    }
    if (((tFemtoV0->StatusNeg() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
        && ((tFemtoV0->StatusNeg() & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      // if (tFemtoV0->StatusNeg() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000
      probMisNeg = fAODpidUtil->GetTOFMismatchProbability(trackneg);
    }

    // if(// (tFemtoV0->StatusPos()& AliESDtrack::kTOFpid)==0 ||
    //    (tFemtoV0->StatusPos()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTOFout)==0 || probMisPos > 0.01)

    if (!(((tFemtoV0->StatusPos() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
          && ((tFemtoV0->StatusPos() & AliVTrack::kTIME) == AliVTrack::kTIME))
        || probMisPos > 0.01) {
      // if(// (tFemtoV0->StatusNeg()&AliESDtrack::kTOFpid)==0 ||
      //    (tFemtoV0->StatusNeg()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTOFout)==0 || probMisNeg > 0.01)
      if (!(((tFemtoV0->StatusNeg() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
            && ((tFemtoV0->StatusNeg() & AliVTrack::kTIME) == AliVTrack::kTIME))
          || probMisNeg > 0.01) {
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
      if (((tFemtoV0->StatusPos() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
          && ((tFemtoV0->StatusPos() & AliVTrack::kTIME) == AliVTrack::kTIME)
          && probMisPos < 0.01) {
        // if(trackpos->IsOn(AliESDtrack::kTOFout & AliESDtrack::kTIME)) {
        tFemtoV0->SetPosNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kKaon));
        tFemtoV0->SetPosNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kProton));
        tFemtoV0->SetPosNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kPion));
      }
      if (((tFemtoV0->StatusNeg() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
          && ((tFemtoV0->StatusNeg() & AliVTrack::kTIME) == AliVTrack::kTIME)
          && probMisNeg < 0.01) {
        // if(trackneg->IsOn(AliESDtrack::kTOFout & AliESDtrack::kTIME)) {
        tFemtoV0->SetNegNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kKaon));
        tFemtoV0->SetNegNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kProton));
        tFemtoV0->SetNegNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kPion));
      }

      double TOFSignalPos = trackpos->GetTOFsignal();
      double TOFSignalNeg = trackneg->GetTOFsignal();
      TOFSignalPos -= fAODpidUtil->GetTOFResponse().GetStartTime(trackpos->P());
      TOFSignalNeg -= fAODpidUtil->GetTOFResponse().GetStartTime(trackneg->P());
      //double pidPos[5];
      //double pidNeg[5];
      //trackpos->GetIntegratedTimes(pidPos);
      //trackneg->GetIntegratedTimes(pidNeg);

      //tFemtoV0->SetTOFPionTimePos(TOFSignalPos - pidPos[2]);
      //tFemtoV0->SetTOFKaonTimePos(TOFSignalPos - pidPos[3]);
      //tFemtoV0->SetTOFProtonTimePos(TOFSignalPos - pidPos[4]);
      //tFemtoV0->SetTOFPionTimeNeg(TOFSignalNeg - pidNeg[2]);
      //tFemtoV0->SetTOFKaonTimeNeg(TOFSignalNeg - pidNeg[3]);
      //tFemtoV0->SetTOFProtonTimeNeg(TOFSignalNeg - pidNeg[4]);
    }
  } else {
    tFemtoV0->SetStatusPos(999);
    tFemtoV0->SetStatusNeg(999);
  }

  tFemtoV0->SetOnFlyStatusV0(tAODv0->GetOnFlyStatus());
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

  double tEtaXi = 0.5 * TMath::Log((TMath::Sqrt(tAODxi->Ptot2Xi()) + tAODxi->MomXiZ())
                                    / (TMath::Sqrt(tAODxi->Ptot2Xi()) - tAODxi->MomXiZ() + 1.e-13));
  tFemtoXi->SetEtaXi(tEtaXi);
  double tPhiXi = TMath::Pi() + TMath::ATan2(-tAODxi->MomXiY(), -tAODxi->MomXiX());
  tFemtoXi->SetPhiXi(tPhiXi);

  tFemtoXi->SetCosPointingAngleXi(tAODxi->CosPointingAngleXi(fVe1[0], fVe1[1], fVe1[2]));
  tFemtoXi->SetCosPointingAngleV0toXi(tAODxi->CosPointingAngle(tAODxi->GetDecayVertexXi()));

  tFemtoXi->SetChargeXi(tAODxi->ChargeXi());
  tFemtoXi->SetptXi(std::sqrt(tAODxi->Pt2Xi()));

  AliNanoAODTrack *trackbac = (AliNanoAODTrack *)tAODxi->GetDecayVertexXi()->GetDaughter(0);

  if (trackbac) {
    tFemtoXi->SetptBac(trackbac->Pt()); //setting pt? px and py was set!

    tFemtoXi->SetEtaBac(trackbac->Eta());            //bac!
    tFemtoXi->SetTPCNclsBac(trackbac->GetTPCNcls()); //bac!
    tFemtoXi->SetNdofBac(trackbac->Chi2perNDF());    //bac!
    tFemtoXi->SetStatusBac(trackbac->GetStatus());   //bac!

    tFemtoXi->SetBacNSigmaTPCK(fAODpidUtil->NumberOfSigmasTPC(trackbac, AliPID::kKaon));
    tFemtoXi->SetBacNSigmaTPCP(fAODpidUtil->NumberOfSigmasTPC(trackbac, AliPID::kProton));
    tFemtoXi->SetBacNSigmaTPCPi(fAODpidUtil->NumberOfSigmasTPC(trackbac, AliPID::kPion));

    //NEED TO ADD:
    //    tFemtoXi->SetNominalTpcEntrancePointBac(tmpVec);
    //    tFemtoXi->SetNominalTpcExitPointBac(tmpVec);
    //    tFemtoXi->SetNominalTpcPointBac(vecTpcPos);
    //    tFemtoXi->SetTPCMomentumBac(trackbac->GetTPCmomentum());
    //    if (fShiftPosition > 0.)

    float bfield = 5 * fMagFieldSign;
    float globalPositionsAtRadiiBac[9][3];
    GetGlobalPositionAtGlobalRadiiThroughTPC(trackbac, bfield, globalPositionsAtRadiiBac);
    double tpcEntranceBac[3] = {globalPositionsAtRadiiBac[0][0],
                                globalPositionsAtRadiiBac[0][1],
                                globalPositionsAtRadiiBac[0][2]};
    double tpcExitBac[3] = {globalPositionsAtRadiiBac[8][0],
                            globalPositionsAtRadiiBac[8][1],
                            globalPositionsAtRadiiBac[8][2]};

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
    for (int i = 0; i < 9; i++) {
      vecTpcBac[i].SetX(globalPositionsAtRadiiBac[i][0]);
      vecTpcBac[i].SetY(globalPositionsAtRadiiBac[i][1]);
      vecTpcBac[i].SetZ(globalPositionsAtRadiiBac[i][2]);
    }

    if (fPrimaryVertexCorrectionTPCPoints) {
      AliFemtoThreeVector tmpVertexVec;
      tmpVertexVec.SetX(fVe1[0]);
      tmpVertexVec.SetY(fVe1[1]);
      tmpVertexVec.SetZ(fVe1[2]);

      for (int i = 0; i < 9; i++) {
        vecTpcBac[i] -= tmpVertexVec;
      }
    }

    tFemtoXi->SetNominalTpcPointBac(vecTpcBac);

    /*    if (fShiftPosition > 0.) {
      Float_t posShiftedBac[3];
      SetShiftedPositions(trackbac, bfield, posShiftedBac, fShiftPosition);
      AliFemtoThreeVector tmpVecBac;
      tmpVecBac.SetX(posShiftedBac[0]);
      tmpVecBac.SetY(posShiftedBac[1]);
      tmpVecBac.SetZ(posShiftedBac[2]);
      tFemtoXi->SetNominalTpcPointBacShifted(tmpVecBac);
    }
    */
    tFemtoXi->SetTPCMomentumBac(trackbac->GetTPCmomentum());
    tFemtoXi->SetdedxBac(trackbac->GetTPCsignal());

    Float_t probMisBac = 1.0;

    if (((tFemtoXi->StatusBac() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
         && ((tFemtoXi->StatusBac() & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      // if (tFemtoXi->StatusBac() & AliESDtrack::kTOFout & AliESDtrack::kTIME) {  //AliESDtrack::kTOFpid=0x8000
      probMisBac = fAODpidUtil->GetTOFMismatchProbability(trackbac);
    }

    // if(// (tFemtoXi->StatusPos()& AliESDtrack::kTOFpid)==0 ||
    //    (tFemtoXi->StatusPos()&AliESDtrack::kTIME)==0 || (tFemtoXi->StatusPos()&AliESDtrack::kTOFout)==0 || probMisPos > 0.01)

    if (!(((tFemtoXi->StatusBac() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
          && ((tFemtoXi->StatusBac() & AliVTrack::kTIME) == AliVTrack::kTIME))
        || probMisBac > 0.01) {
      tFemtoXi->SetBacNSigmaTOFK(-1000);
      tFemtoXi->SetBacNSigmaTOFP(-1000);
      tFemtoXi->SetBacNSigmaTOFPi(-1000);

      tFemtoXi->SetTOFProtonTimeBac(-1000);
      tFemtoXi->SetTOFPionTimeBac(-1000);
      tFemtoXi->SetTOFKaonTimeBac(-1000);
    } else {
      if (((tFemtoXi->StatusBac() & AliVTrack::kTOFout) == AliVTrack::kTOFout)
          && ((tFemtoXi->StatusBac() & AliVTrack::kTIME) == AliVTrack::kTIME)
          && probMisBac < 0.01) {
        tFemtoXi->SetBacNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackbac, AliPID::kKaon));
        tFemtoXi->SetBacNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackbac, AliPID::kProton));
        tFemtoXi->SetBacNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackbac, AliPID::kPion));
      }

      double TOFSignalBac = trackbac->GetTOFsignal();
      TOFSignalBac -= fAODpidUtil->GetTOFResponse().GetStartTime(trackbac->P());
      //double pidBac[5];
      //trackbac->GetIntegratedTimes(pidBac);

      //tFemtoXi->SetTOFPionTimeBac(TOFSignalBac - pidBac[2]);
      //tFemtoXi->SetTOFKaonTimeBac(TOFSignalBac - pidBac[3]);
      //tFemtoXi->SetTOFProtonTimeBac(TOFSignalBac - pidBac[4]);
    }
  } else {
    tFemtoXi->SetStatusBac(999);
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

