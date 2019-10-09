///
/// \file PWGCF/FEMTOSCOPY/AliFemto/AliFemtoEventReaderAOD.h
///

#ifndef ALIFEMTOEVENTREADERNANOAOD_H
#define ALIFEMTOEVENTREADERNANOAOD_H
#include "AliFemtoEventReader.h"
#include "AliFemtoEnumeration.h"

#include <string>

#include "TTree.h"
#include "TChain.h"
#include "TBits.h"
#include "THnSparse.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliFemtoV0.h"
#include "AliFemtoXi.h"
#include "AliAODpidUtil.h"
#include "AliNanoAODHeader.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"


//#include "/home/wfpw/alice5/AliPhysics/PWG/DevNanoAOD/AliNanoAODTrack.h"
#include "AliNanoAODTrack.h"

class AliFemtoEvent;
class AliFemtoTrack;

class AliFemtoEventReaderNanoAOD : public AliFemtoEventReader {
public:
  AliFemtoEventReaderNanoAOD();
  AliFemtoEventReaderNanoAOD(const AliFemtoEventReaderNanoAOD &aReader);
  virtual ~AliFemtoEventReaderNanoAOD();

  AliFemtoEventReaderNanoAOD &operator=(const AliFemtoEventReaderNanoAOD &aReader);

  virtual AliFemtoEvent *ReturnHbtEvent();
  AliFemtoString Report();
  void SetInputFile(const char *inputfile);
  void SetFilterBit(UInt_t ibit);
  void SetFilterMask(int ibit);
  void SetReadMC(unsigned char a);
  void SetReadV0(unsigned char a);
  void SetReadCascade(unsigned char a);
  void SetAODpidUtil(AliAODpidUtil *aAODpidUtil);
  void SetAODheader(AliNanoAODHeader *aAODheader);
  void SetMagneticFieldSign(int s);
  void SetEPVZERO(Bool_t);
  void GetGlobalPositionAtGlobalRadiiThroughTPC(AliNanoAODTrack *track, Float_t bfield, Float_t globalPositionsAtRadii[9][3]);
  void SetUseMultiplicity(std::string aType);
  void SetpA2013(Bool_t pa2013); ///< set vertex configuration for pA (2013): IsVertexSelected2013pA
  void SetUseMVPlpSelection(Bool_t mvplp);
  void SetUseOutOfBunchPlpSelection(Bool_t outOfBunchPlp);
  void SetIsPileUpEvent(Bool_t ispileup);
  void SetCascadePileUpRemoval(Bool_t cascadePileUpRemoval);
  void SetV0PileUpRemoval(Bool_t v0PileUpRemoval);
  void SetTrackPileUpRemoval(Bool_t trackPileUpRemoval);

  void SetMinVtxContr(Int_t contr = 1) {
    fMinVtxContr = contr;
  }
  void SetMinPlpContribMV(Int_t minPlpContribMV) {
    fMinPlpContribMV = minPlpContribMV;
  }
  void SetMinPlpContribSPD(Int_t minPlpContribSPD) {
    fMinPlpContribSPD = minPlpContribSPD;
  }
  void SetDCAglobalTrack(Int_t dcagt);

  bool RejectEventCentFlat(float MagField, float CentPercent);
  void SetCentralityFlattening(Bool_t flat);
  void SetShiftPosition(Double_t rad);

  void SetPrimaryVertexCorrectionTPCPoints(bool correctTpcPoints);
  //void SetShiftedPositions(const AliAODTrack *track ,const Float_t bfield, Float_t posShifted[3], const Double_t radius=1.25);

  void SetUseAliEventCuts(Bool_t useAliEventCuts);
  void SetReadFullMCData(Bool_t should_read=true);
  bool GetReadFullMCData() const;

  void SetInputEvent(AliVEvent *inputEvent){fEvent = inputEvent;};

  void SetCovMatPresent(Bool_t pres){fCovMatPresent = pres;}


  void Set1DCorrectionsPions(TH1D *h1);
  void Set1DCorrectionsKaons(TH1D *h1);
  void Set1DCorrectionsProtons(TH1D *h1);
  void Set1DCorrectionsPionsMinus(TH1D *h1);
  void Set1DCorrectionsKaonsMinus(TH1D *h1);
  void Set1DCorrectionsProtonsMinus(TH1D *h1);
  void Set1DCorrectionsAll(TH1D *h1);
  void Set1DCorrectionsLambdas(TH1D *h1);
  void Set1DCorrectionsLambdasMinus(TH1D *h1);
  void Set1DCorrectionsXiPlus(TH1D *h1);
  void Set1DCorrectionsXiMinus(TH1D *h1);

  void Set4DCorrectionsPions(THnSparse *h1);
  void Set4DCorrectionsKaons(THnSparse *h1);
  void Set4DCorrectionsProtons(THnSparse *h1);
  void Set4DCorrectionsPionsMinus(THnSparse *h1);
  void Set4DCorrectionsKaonsMinus(THnSparse *h1);
  void Set4DCorrectionsProtonsMinus(THnSparse *h1);
  void Set4DCorrectionsAll(THnSparse *h1);
  void Set4DCorrectionsLambdas(THnSparse *h1);
  void Set4DCorrectionsLambdasMinus(THnSparse *h1);

protected:
  virtual AliFemtoEvent *CopyAODtoFemtoEvent();
  virtual AliFemtoTrack *CopyAODtoFemtoTrack(AliNanoAODTrack *tAodTrack
      //            AliPWG2AODTrack *tPWG2AODTrack
                                            );
  virtual AliFemtoV0 *CopyAODtoFemtoV0(AliAODv0 *tAODv0);
  virtual AliFemtoXi *CopyAODtoFemtoXi(AliAODcascade *tAODxi);


  int            fNumberofEvent;    ///< number of Events in AOD file
  int            fCurEvent;         ///< number of current event
  //AliAODEvent   *fEvent;            ///< AOD event
  AliVEvent   *fEvent;            ///< V event
  TBits          fAllTrue;          ///< Bit set with all true bits
  TBits          fAllFalse;         ///< Bit set with all false bits
  UInt_t         fFilterBit;        ///< Bitmap bit for AOD filters
  UInt_t         fFilterMask;
  //  TClonesArray*  fPWG2AODTracks;    // Link to PWG2 specific AOD information (if it exists)

  unsigned char  fReadMC;           ///< Attempt to read the MC information from the AOD
  unsigned char  fReadV0;           ///< Read V0 information from the AOD and put it into V0Collection
  unsigned char  fReadCascade;      ///< Read Cascade information from the AOD and put it into V0Collection
  unsigned char  fUsePreCent;       ///< Use centrality pre-selection to speed up analysis
  string   fEstEventMult;     ///< Type of the event multiplicity estimator
  double         fCentRange[2];     ///< Centrality pre-selection range
  AliAODpidUtil *fAODpidUtil;
  AliNanoAODHeader *fAODheader;
  AliAnalysisUtils *fAnaUtils;
  AliEventCuts     *fEventCuts;
  Bool_t           fUseAliEventCuts;

  /// Read generated particle info of "low quality" MC tracks
  /// (i.e. tracks with negative labels)
  Bool_t           fReadFullMCData;

private:


  string fInputFile;       ///< name of input file with AOD filenames
  TChain *fTree;           ///< AOD tree
  TFile *fAodFile;         ///< AOD file
  int fMagFieldSign;       ///< Magnetic field sign
  Bool_t fisEPVZ;          ///< to get event plane angle from VZERO
  Bool_t fpA2013;          ///< analysis on pA 2013 data
  Bool_t fisPileUp;        ///< pile up rejection on?
  Bool_t fCascadePileUpRemoval;//pile-up removal for cascades (its+tof hits for pos, neg and bac tracks)
  Bool_t fV0PileUpRemoval;//pile-up removal for V0s
  Bool_t fTrackPileUpRemoval;//pile-up removal for tracks (its+tof hits of tracks)
  Bool_t fMVPlp;           ///< multi-vertex pileup rejection?
  Bool_t fOutOfBunchPlp;   ///out-of-bunch pileup rejection
  Int_t fMinVtxContr;      ///< no of contributors for pA 2013 data
  Int_t fMinPlpContribMV;  ///< no of contributors for multivertex pile-up rejection
  Int_t fMinPlpContribSPD; ///< no of contributors for SPD pile-up rejection
  Int_t fDCAglobalTrack;  ///< to get DCA from global tracks instead of TPC-only
  Bool_t fFlatCent;        ///< Boolean determining if the user should flatten the centrality
  Bool_t fPrimaryVertexCorrectionTPCPoints; ///< Boolean determining if the reader should shift all TPC points to be relative to event vertex
  Double_t fShiftPosition; ///< radius at which the spatial position of the track in the shifted coordinate system is calculated
  Bool_t fCovMatPresent; /// flag if covariance matrix is not present in NanoAOD

  TH1D *f1DcorrectionsPions;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsKaons;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsProtons;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsPionsMinus;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsKaonsMinus;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsProtonsMinus;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsAll;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsLambdas;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsLambdasMinus;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsXiPlus;    ///<file with corrections, pT dependant
  TH1D *f1DcorrectionsXiMinus;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsPions;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsKaons;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsProtons;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsPionsMinus;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsKaonsMinus;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsProtonsMinus;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsAll;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsLambdas;    ///<file with corrections, pT dependant
  THnSparse *f4DcorrectionsLambdasMinus;    ///<file with corrections, pT dependant

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventReaderNanoAOD, 13);
  /// \endcond
#endif

};


inline void AliFemtoEventReaderNanoAOD::SetReadFullMCData(bool read)
  { fReadFullMCData = read; }

inline bool AliFemtoEventReaderNanoAOD::GetReadFullMCData() const
  { return fReadFullMCData; }

#endif
