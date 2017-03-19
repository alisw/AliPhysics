/// \file AliQnCorrectionsProfileChannelizedIngress.cxx
/// \brief Implementation of the multidimensional ingress channelized profile 

#include "TList.h"

#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfileChannelizedIngress.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsProfileChannelizedIngress);
/// \endcond

/// Default constructor
AliQnCorrectionsProfileChannelizedIngress::AliQnCorrectionsProfileChannelizedIngress() :
    AliQnCorrectionsHistogramBase() {

  fValues = NULL;
  fGroupValues = NULL;
  fValidated = NULL;
  fUsedChannel = NULL;
  fChannelGroup = NULL;
  fNoOfChannels = 0;
  fActualNoOfChannels = 0;
  fChannelMap = NULL;
  fUseGroups = kFALSE;
  fUsedGroup = NULL;
  fActualNoOfGroups = 0;
  fNoOfGroups = 0;
  fGroupMap = NULL;
}

/// Normal constructor
///
/// Stores the set of variables that identify the
/// different event classes passing them to its parent
/// and prepares the object for actual histogram
/// attachment
///
/// \param name base for the name of the histograms
/// \param title base for the title of the histograms
/// \param ecvs the event classes variables set
/// \param nNoOfChannels the number of channels associated
/// \param option option for errors computation
///     ' '  (Default) the bin errors are the standard error on the mean of the
///          bin values
///
///     's'            the bin are the standard deviation of of the bin values
AliQnCorrectionsProfileChannelizedIngress::AliQnCorrectionsProfileChannelizedIngress(const char *name,
    const char *title,
    AliQnCorrectionsEventClassVariablesSet &ecvs,
    Int_t nNoOfChannels,
    Option_t *option) : AliQnCorrectionsHistogramBase(name, title, ecvs, option) {

  fValues = NULL;
  fGroupValues = NULL;
  fValidated = NULL;
  fUsedChannel = NULL;
  fChannelGroup = NULL;
  fNoOfChannels = nNoOfChannels;
  fActualNoOfChannels = 0;
  fChannelMap = NULL;
  fUseGroups = kFALSE;
  fUsedGroup = NULL;
  fActualNoOfGroups = 0;
  fNoOfGroups = 0;
  fGroupMap = NULL;
}

/// Default destructor
/// Releases the memory taken
AliQnCorrectionsProfileChannelizedIngress::~AliQnCorrectionsProfileChannelizedIngress() {

  if (fUsedChannel != NULL) delete [] fUsedChannel;
  if (fChannelGroup != NULL) delete [] fChannelGroup;
  if (fChannelMap != NULL) delete [] fChannelMap;
  if (fValues != NULL) delete fValues;
  if (fGroupValues != NULL) delete fGroupValues;
  if (fValidated != NULL) delete fValidated;
  if (fUsedGroup != NULL) delete [] fUsedGroup;
  if (fGroupMap != NULL) delete [] fGroupMap;
}



/// Attaches existing histograms as the support histograms for the profile function
///
/// The histograms are located in the passed list and if found and with the
/// proper dimensions their references are stored in member variables.
///
/// Channel information is used to build internal structures such as
/// the channel map and the actual number of channels and the channels groups
/// and the actual number of groups. The information
/// is matched with the found histogram to validate it. If
/// bUsedChannel is NULL all channels
/// within fNoOfChannels are assigned to this profile.
/// If nChannelGroup is NULL all channels assigned to this
/// profile are allocated to the same group.
///
/// Once the histograms are found and validated, a unique value / error channel histogram
/// is created for efficient access and a potential unique value / error channels group
/// histogram is created.
/// \param histogramList list where the histograms have to be located
/// \param bUsedChannel array of booleans one per each channel
/// \param nChannelGroup array of group number for each channel
/// \return true if properly attached else false
Bool_t AliQnCorrectionsProfileChannelizedIngress::AttachHistograms(TList *histogramList, const Bool_t *bUsedChannel, const Int_t *nChannelGroup) {
  /* let's build the histograms names */
  TString histoName = GetName();
  TString entriesHistoName = GetName(); entriesHistoName += szEntriesHistoSuffix;

  /* initialize. Remember we own the histograms */
  if (fValues != NULL) delete fValues;
  if (fGroupValues != NULL) delete fGroupValues;
  if (fValidated != NULL) delete fValidated;
  fValues = NULL;
  fGroupValues = NULL;
  fValidated = NULL;
  if (fUsedChannel != NULL) delete [] fUsedChannel;
  if (fChannelGroup != NULL) delete [] fChannelGroup;
  if (fChannelMap != NULL) delete [] fChannelMap;
  if (fUsedGroup != NULL) delete [] fUsedGroup;
  if (fGroupMap != NULL) delete [] fGroupMap;


  /* lets consider now the channel information */
  fUsedChannel = new Bool_t[fNoOfChannels];
  fChannelGroup = new Int_t[fNoOfChannels];
  fChannelMap = new Int_t[fNoOfChannels];
  for (Int_t i = 0; i < fNoOfChannels; i++) {
    fUsedChannel[i] = kFALSE;
    fChannelGroup[i] = 0;
    fChannelMap[i] = -1;
  }

  Int_t nMinGroup = 0xFFFF;
  Int_t nMaxGroup = 0x0000;
  fActualNoOfChannels = 0;
  for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
    if (bUsedChannel != NULL) {
      fUsedChannel[ixChannel] = bUsedChannel[ixChannel];
    }
    else {
      fUsedChannel[ixChannel] = kTRUE;
    }
    if (fUsedChannel[ixChannel]) {
      if (nChannelGroup != NULL) {
        fChannelGroup[ixChannel] = nChannelGroup[ixChannel];
        /* update min max group number */
        if (nChannelGroup[ixChannel] < nMinGroup)
          nMinGroup = nChannelGroup[ixChannel];
        if (nMaxGroup < nChannelGroup[ixChannel])
          nMaxGroup = nChannelGroup[ixChannel];
      }
      else {
        fChannelGroup[ixChannel] = 0;
        nMinGroup = 0;
        nMaxGroup = 0;
      }
      fChannelMap[ixChannel] = fActualNoOfChannels;
      fActualNoOfChannels++;
    }
  }
  fUseGroups = (nChannelGroup != NULL) && (nMinGroup != nMaxGroup);

  if (fUseGroups) {
    /* let's build the groups support structures */
    fNoOfGroups = nMaxGroup + 1; /* just in case group number starts from zero */
    fUsedGroup = new Bool_t[fNoOfGroups];
    fGroupMap = new Int_t[fNoOfGroups];
    for (Int_t i = 0; i < fNoOfGroups; i++) {
      fUsedGroup[i] = kFALSE;
      fGroupMap[i] = -1;
    }
    fActualNoOfGroups = 0;
    for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
      if (fUsedChannel[ixChannel]) {
        if (fUsedGroup[fChannelGroup[ixChannel]]) {
          /* group already considered */
          continue;
        }
        else {
          /* new group number */
          fUsedGroup[fChannelGroup[ixChannel]] = kTRUE;
          fGroupMap[fChannelGroup[ixChannel]] = fActualNoOfGroups;
          fActualNoOfGroups++;
        }
      }
    }
  }

  /* let's first try the Values / Entries structure */
  THnI *origEntries = (THnI *) histogramList->FindObject((const char*) entriesHistoName);
  if (origEntries != NULL && origEntries->GetEntries() != 0) {
    /* so we get it! */
    /* let's check the channel axis */
    if (fActualNoOfChannels != origEntries->GetAxis(fEventClassVariables.GetEntriesFast())->GetNbins())
      return kFALSE;
    THnF *origValues = (THnF *) histogramList->FindObject((const char *)histoName);
    if (origValues == NULL)
      return kFALSE;
    /* let's check the channel axis */
    if (fActualNoOfChannels != origValues->GetAxis(fEventClassVariables.GetEntriesFast())->GetNbins())
      return kFALSE;

    /* so we got the original histograms */
    /* now we should build the definitive histogram value / error */
    /* and the group value / error histogram if applicable */

    /* let's first prepare the bin validation information */
    /* we open space for channel variable as well */
    Int_t nVariables = fEventClassVariables.GetEntriesFast();
    Double_t *minvals = new Double_t[nVariables+1];
    Double_t *maxvals = new Double_t[nVariables+1];
    Int_t *nbins = new Int_t[nVariables+1];
    /* get the multidimensional structure */
    fEventClassVariables.GetMultidimensionalConfiguration(nbins,minvals,maxvals);
    minvals[nVariables] = -0.5;
    maxvals[nVariables] = -0.5 + fActualNoOfChannels;
    nbins[nVariables] = fActualNoOfChannels;
    /* create the values and entries multidimensional histograms */
    fValidated = new THnC(Form("%s_Validated",(const char *) histoName), Form("%s_Validated",(const char *) histoName),nVariables+1,nbins,minvals,maxvals);
    /* and now the definitive histogram value /error getting validation information */
    fValues = DivideTHnF(origValues, origEntries,fValidated);

    if (fUseGroups) {
      /* let's then build the groups histogram */
      TString histoGroupName = szGroupHistoPrefix;
      histoGroupName += GetName();
      TString histoGroupTitle = szGroupHistoPrefix;
      histoGroupTitle += GetTitle();

      /* There will be a wrong external view of the channel number especially */
      /* manifested when there are holes in the channel assignment */
      /* so, lets complete the dimension information */
      /* WARNING: be aware that ROOT does not keep label information when projecting THn */
      minvals[nVariables] = -0.5;
      maxvals[nVariables] = -0.5 + fActualNoOfGroups;
      nbins[nVariables] = fActualNoOfGroups;

      /* create the values and entries multidimensional histograms */
      fGroupValues = new THnF((const char *) histoGroupName, (const char *) histoGroupTitle,nVariables+1,nbins,minvals,maxvals);

      /* now let's set the proper binning and label on each axis */
      for (Int_t var = 0; var < nVariables; var++) {
        fGroupValues->GetAxis(var)->Set(fEventClassVariables.At(var)->GetNBins(),fEventClassVariables.At(var)->GetBins());
        fGroupValues->GetAxis(var)->SetTitle(fEventClassVariables.At(var)->GetVariableLabel());
      }

      /* and now the channel axis */
      fGroupValues->GetAxis(nVariables)->SetTitle(szGroupAxisTitle);
      /* and the proper group number labels if needed */
      if (fActualNoOfGroups != fNoOfGroups) {
        for (Int_t ixGroup = 0; ixGroup < fNoOfGroups; ixGroup++) {
          if (fUsedGroup[ixGroup]) {
            fGroupValues->GetAxis(nVariables)->SetBinLabel(fGroupMap[ixGroup]+1, Form("%d", ixGroup));
          }
        }
      }

      fGroupValues->Sumw2();

      /* now let's build its content */
      /* the procedure is as follows: we will project and add together the values histogram */
      /* of the channels corresponding to a group number and the number of entries within these */
      /* channels then we divide both sums and then store the result in the corresponding group values */
      Int_t *dimToProject = new Int_t[nVariables];
      for (Int_t var = 0; var < nVariables; var++)
        dimToProject[var] = var;
      Int_t *binsArray = new Int_t[nVariables+1];

      for (Int_t ixGroup = 0; ixGroup < fNoOfGroups; ixGroup++) {
        if (fUsedGroup[ixGroup]) {
          /* new group, let's add its contributor channels */
          THnF *hCumProjected = NULL;
          THnI *hCumProjectedEntries = NULL;
          for (Int_t ixChannel = 0; ixChannel < fNoOfChannels; ixChannel++) {
            if (fUsedChannel[ixChannel]) {
              if (fChannelGroup[ixChannel] == ixGroup) {
                /* channel within the group found let's add its content */
                /* first filter the projection */
                origValues->GetAxis(nVariables)->SetRange(fChannelMap[ixChannel]+1, fChannelMap[ixChannel]+1);
                origEntries->GetAxis(nVariables)->SetRange(fChannelMap[ixChannel]+1, fChannelMap[ixChannel]+1);
                if (hCumProjected != NULL) {
                  /* let's accumulate the new channel */
                  THnF *hProjected = (THnF *) origValues->Projection(nVariables,dimToProject, "E");
                  THnI *hProjectedEntries = (THnI *) origEntries->Projection(nVariables,dimToProject);
                  hCumProjected->Add(hProjected);
                  hCumProjectedEntries->Add(hProjectedEntries);
                  delete hProjected;
                  delete hProjectedEntries;
                }
                else {
                  /* first channel in the group */
                  hCumProjected = (THnF *) origValues->Projection(nVariables,dimToProject, "E");
                  hCumProjectedEntries = (THnI *) origEntries->Projection(nVariables,dimToProject);
                }
              }
            }
          }
          /* let's build the final group weight */
          THnF *hChannelsGroupWeights = DivideTHnF(hCumProjected, hCumProjectedEntries);

          /* let's store the channels contribution to the group values for the group */
          /* the group for which we are storing values */
          for (Int_t var = 0; var < nVariables; var++)
            binsArray[var] = 0;
          binsArray[nVariables] = fGroupMap[ixGroup] + 1;
          /* remember, we copy the values because we are visiting each group once */
          CopyTHnF(fGroupValues, hChannelsGroupWeights, binsArray);
          delete hCumProjected;
          delete hCumProjectedEntries;
          delete hChannelsGroupWeights;
        }
      }
      /* reset the ranges */
      origValues->GetAxis(nVariables)->SetRange(0, 0);
      origEntries->GetAxis(nVariables)->SetRange(0, 0);
      delete [] dimToProject;
      delete [] binsArray;
    }
    /* we finished here with this stuff */
    delete [] minvals;
    delete [] maxvals;
    delete [] nbins;
  }
  else
    return kFALSE;

  return kTRUE;
}

/// Get the bin number for the current variable content and passed channel
///
/// The bin number identifies the event class the current
/// variable content points to under the passed channel.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsProfileChannelizedIngress::GetBin(const Float_t *variableContainer, Int_t nChannel) {

  /* store also the channel number */
  FillBinAxesValues(variableContainer, fChannelMap[nChannel]);
  return fValues->GetBin(fBinAxesValues);
}

/// Check the validity of the content of the passed bin
/// For the time being this kind of histograms cannot check
/// bin content validity so, kTRUE is returned.
/// \param bin the bin to check its content validity
/// \return kTRUE if the content is valid kFALSE otherwise
Bool_t AliQnCorrectionsProfileChannelizedIngress::BinContentValidated(Long64_t bin) {

  if (fValidated->GetBinContent(bin) < 0.5) {
    return kFALSE;
  }
  else {
    return kTRUE;
  }
}

/// Get the bin content for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// is requested. The bin average stored as bin content is returned
///
/// \param bin the interested bin number
/// \return the bin number content
Float_t AliQnCorrectionsProfileChannelizedIngress::GetBinContent(Long64_t bin) {

  return fValues->GetBinContent(bin);
}

/// Get the bin content error for the passed bin number
///
/// The bin number identifies a desired event class whose content
/// error is requested. The bin error stored as bin error is returned.
///
/// \param bin the interested bin number
/// \return the bin number content error
Float_t AliQnCorrectionsProfileChannelizedIngress::GetBinError(Long64_t bin) {

  return fValues->GetBinError(bin);
}

/// Get the bin number for the current variable content and passed channel group number
///
/// The bin number identifies the event class the current
/// variable content points to under the passed channel group number.
///
/// \param variableContainer the current variables content addressed by var Id
/// \param nChannel the interested external channel number which group number is asked
/// \return the associated bin to the current variables content
Long64_t AliQnCorrectionsProfileChannelizedIngress::GetGrpBin(const Float_t *variableContainer, Int_t nChannel) {

  /* check the groups structures are in place */
  if (fUseGroups) {
    /* store also the group number */
    FillBinAxesValues(variableContainer, fGroupMap[fChannelGroup[nChannel]]);
    return fGroupValues->GetBin(fBinAxesValues);
  }
  return -1;
}

/// Get the group bin content for the passed bin number
///
/// The bin number identifies a desired event class whose group content
/// is requested.
///
/// \param bin the interested group bin number
/// \return the group bin number content
Float_t AliQnCorrectionsProfileChannelizedIngress::GetGrpBinContent(Long64_t bin) {

  /* check the groups structures are in place */
  if (fUseGroups) {
    return fGroupValues->GetBinContent(bin);
  }
  return 1.0;
}

/// Get the group bin content error for the passed bin number
///
/// The group bin number identifies a desired event class whose content
/// error is requested.
///
/// \param bin the interested group bin number
/// \return the bin number content error
Float_t AliQnCorrectionsProfileChannelizedIngress::GetGrpBinError(Long64_t bin) {

  /* check the groups structures are in place */
  if (fUseGroups) {
    return fGroupValues->GetBinError(bin);
  }
  return 1.0;
}


