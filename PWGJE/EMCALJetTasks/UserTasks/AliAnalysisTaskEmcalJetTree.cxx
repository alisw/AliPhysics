/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliLog.h"
#include "AliAnalysisTaskEmcalJetTree.h"

// Definitions of class AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryBase

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryBase)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryBase::AliEmcalJetInfoSummaryBase(const AliEmcalJetInfo& source) :
  fPt(0),
  fEta(0),
  fPhi(0),
  fNEF(0),
  fZLeading(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryBase::Reset()
{
  fPt = 0;
  fEta = 0;
  fPhi = 0;
  fNEF = 0;
  fZLeading = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryBase::Set(const AliEmcalJetInfo& source)
{
  fPt = source.Pt();
  fEta = source.Eta();
  fPhi = source.Phi_0_2pi();
  fNEF = source.fNEF;
  fZLeading = source.fZ;
}

// Definitions of class AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPP

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPP)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPP::AliEmcalJetInfoSummaryPP(const AliEmcalJetInfo& source) :
  AliEmcalJetInfoSummaryBase(),
  fNConstituents(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPP::Reset()
{
  AliEmcalJetInfoSummaryBase::Reset();
  fNConstituents = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPP::Set(const AliEmcalJetInfo& source)
{
  AliEmcalJetInfoSummaryBase::Set(source);
  fNConstituents = Char_t(source.fNConstituents);
}

// Definitions of class AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb::AliEmcalJetInfoSummaryPbPb(const AliEmcalJetInfo& source) :
  AliEmcalJetInfoSummaryBase(),
  fCent(0),
  fEP(0),
  fArea(0),
  fNConstituents(0),
  fCorrPt(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb::Reset()
{
  AliEmcalJetInfoSummaryBase::Reset();
  fCent = 0;
  fEP = 0;
  fArea = 0;
  fNConstituents = 0;
  fCorrPt = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb::Set(const AliEmcalJetInfo& source)
{
  AliEmcalJetInfoSummaryBase::Set(source);
  fCent = Char_t(source.fCent);
  fEP = source.fEP;
  fArea = source.fArea;
  fNConstituents = Short_t(source.fNConstituents);
  fCorrPt = source.fCorrPt;
}

// Definitions of class AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryEmbedding)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryEmbedding::AliEmcalJetInfoSummaryEmbedding(const AliEmcalJetInfo& source) :
  AliEmcalJetInfoSummaryPbPb(),
  fMCPt(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryEmbedding::Reset()
{
  AliEmcalJetInfoSummaryPbPb::Reset();
  fMCPt = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryEmbedding::Set(const AliEmcalJetInfo& source)
{
  AliEmcalJetInfoSummaryPbPb::Set(source);
  fMCPt = source.fMCPt;
}

// Definitions of class AliAnalysisTaskEmcalJetTreeBase

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetTreeBase)
/// \endcond

AliAnalysisTaskEmcalJetTreeBase::AliAnalysisTaskEmcalJetTreeBase(const char *name) :
  AliAnalysisTaskEmcalJetSpectraQA(name),
  fTree(0)
{
  DefineOutput(2, TTree::Class());
  SetHistoType(kTTree);
}

/// Static method to create a specialized instance of AliAnalysisTaskEmcalJetTree
///
/// \param name Name of the task
/// \param type Analysis type (see enum definition)
/// \return pointer to a new instance of AliAnalysisTaskEmcalJetTree<T>
AliAnalysisTaskEmcalJetTreeBase* AliAnalysisTaskEmcalJetTreeBase::CreateInstance(const char* name, EAnalisysType_t type)
{
  switch (type) {
  case kJetPP:
    ::Info("AliAnalysisTaskEmcalJetTreeBase::CreateInstance", "Created an instance of AliAnalysisTaskEmcalJetTree<AliEmcalJetInfoSummaryPP>");
    return new AliAnalysisTaskEmcalJetTree<AliEmcalJetInfoSummaryPP>(name);
    break;
  case kJetPbPb:
    ::Info("AliAnalysisTaskEmcalJetTreeBase::CreateInstance", "Created an instance of AliAnalysisTaskEmcalJetTree<AliEmcalJetInfoSummaryPbPb>");
    return new AliAnalysisTaskEmcalJetTree<AliEmcalJetInfoSummaryPbPb>(name);
    break;
  case kJetEmbedding:
    ::Info("AliAnalysisTaskEmcalJetTreeBase::CreateInstance", "Created an instance of AliAnalysisTaskEmcalJetTree<AliEmcalJetInfoSummaryEmbedding>");
    return new AliAnalysisTaskEmcalJetTree<AliEmcalJetInfoSummaryEmbedding>(name);
    break;
  default:
    ::Error("AliAnalysisTaskEmcalJetTreeBase::CreateInstance", "Type %d not implemented!", type);
    return 0;
  }
}

// Definitions of class AliAnalysisTaskEmcalJetTree

/// \cond CLASSIMP
templateClassImp(AliAnalysisTaskEmcalJetTree)
/// \endcond

/// Default constructor for ROOT I/O purposes
template <class T>
AliAnalysisTaskEmcalJetTree<T>::AliAnalysisTaskEmcalJetTree() :
  AliAnalysisTaskEmcalJetTreeBase("AliAnalysisTaskEmcalJetTree"),
  fCurrentOutput(0)
{
}

/// Specialized default constructor (AliEmcalJetInfoSummaryPP) for ROOT I/O purposes
/// This is needed to address a "feature" (aka bug) of ROOT CINT, to be checked (maybe fixed in newer versions of ROOT)
template <>
AliAnalysisTaskEmcalJetTree<AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPP>::AliAnalysisTaskEmcalJetTree() :
  AliAnalysisTaskEmcalJetTreeBase("AliAnalysisTaskEmcalJetTree"),
  fCurrentOutput(0)
{
}

/// Specialized default constructor (AliEmcalJetInfoSummaryPbPb) for ROOT I/O purposes
/// This is needed to address a "feature" (aka bug) of ROOT CINT, to be checked (maybe fixed in newer versions of ROOT)
template <>
AliAnalysisTaskEmcalJetTree<AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryPbPb>::AliAnalysisTaskEmcalJetTree() :
  AliAnalysisTaskEmcalJetTreeBase("AliAnalysisTaskEmcalJetTree"),
  fCurrentOutput(0)
{
}

/// Specialized default constructor (AliEmcalJetInfoSummaryEmbedding) for ROOT I/O purposes
/// This is needed to address a "feature" (aka bug) of ROOT CINT, to be checked (maybe fixed in newer versions of ROOT)
template <>
AliAnalysisTaskEmcalJetTree<AliAnalysisTaskEmcalJetTreeBase::AliEmcalJetInfoSummaryEmbedding>::AliAnalysisTaskEmcalJetTree() :
  AliAnalysisTaskEmcalJetTreeBase("AliAnalysisTaskEmcalJetTree"),
  fCurrentOutput(0)
{
}

/// Standard named constructor
///
/// \param name Name of the task
template <class T>
AliAnalysisTaskEmcalJetTree<T>::AliAnalysisTaskEmcalJetTree(const char *name) :
  AliAnalysisTaskEmcalJetTreeBase(name),
  fCurrentOutput(0)
{
}

/// Allocate output TTree for a jet container
///
/// \param jets Valid pointer to an AliJetContainer object
template <class T>
void AliAnalysisTaskEmcalJetTree<T>::AllocateTTree(const AliJetContainer* jets)
{
  typename std::map<std::string,std::vector<T> >::iterator it = (fCurrentOutput->insert(std::pair<std::string,std::vector<T> >(jets->GetName(), std::vector<T>()))).first;
  fTree->Branch(jets->GetName(), &((*it).second));
}

/// Overloads base class method. Creates output objects
template <class T>
void AliAnalysisTaskEmcalJetTree<T>::UserCreateOutputObjects()
{
  fCurrentOutput = new std::map<std::string, std::vector<T> >();
  TString treeName = TString::Format("%s_jets", GetName());
  fTree = new TTree(treeName, treeName);

  AliAnalysisTaskEmcalJetSpectraQA::UserCreateOutputObjects();

  PostData(2, fTree);
}

/// Fill tree with jet info
///
/// \param jet  Jet containing the information to be sent to the tree/histograms
template <class T>
void AliAnalysisTaskEmcalJetTree<T>::FillTTree(const AliEmcalJetInfo& jet, const AliJetContainer* jets)
{
  static typename std::map<std::string, std::vector<T> >::iterator it = fCurrentOutput->end();

  if (it == fCurrentOutput->end() || TString(jets->GetName()) != it->first) {
    it = fCurrentOutput->find(jets->GetName());
    if (it == fCurrentOutput->end()) return;
  }

  it->second.push_back(T(jet));
}

/// Overloads base class method. Fills the output histograms
///
/// \return kTRUE if successful
template <class T>
Bool_t AliAnalysisTaskEmcalJetTree<T>::FillHistograms()
{
  typedef typename std::map<std::string, std::vector<T> >::iterator iterator_type;

  for (iterator_type it = fCurrentOutput->begin(); it != fCurrentOutput->end(); it++) {
    it->second.clear();
  }

  Bool_t r = AliAnalysisTaskEmcalJetSpectraQA::FillHistograms();
  if (!r) return kFALSE;
  fTree->Fill();
  PostData(2, fTree);
  return kTRUE;
}
