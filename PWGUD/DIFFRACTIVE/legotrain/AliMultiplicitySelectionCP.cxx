#include "TMath.h"
#include "TObject.h"
#include "TIterator.h"
#include "TList.h"

#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliSPDUtils.h"
#include "AliITSsegmentationSPD.h"
#include "AliESDv0.h"

#include "AliMultiplicitySelectionCP.h"


ClassImp(AliMultiplicitySelectionCP)

AliMultiplicitySelectionCP::AliMultiplicitySelectionCP():TObject(),
  fkCheckReferenceMultiplicity(0)
{

  fTrackCutListPrim = new TList();
  fTrackCutListPrim->SetOwner();
  fTrackCutListPrim->SetName("PrimaryTrackCut");
 
  SetTPCnclsS();
  SetTrackDCAz();
  SetTrackEtaRange();

  IgnoreV0s();
  for(Int_t i = 0; i< fkNtrackMax; i++)
    fkIsTrackSec[i]= kFALSE;



}



AliMultiplicitySelectionCP::~AliMultiplicitySelectionCP()
{
  
  if(fTrackCutListPrim)
    {
      fTrackCutListPrim->Delete();
      delete fTrackCutListPrim;
    }
  fTrackCutListPrim = 0;
}


void AliMultiplicitySelectionCP::InitDefaultTrackCuts(Int_t clusterCut)
{
  /*
    Alexander Kalweit 
    Email to PWG conveners on 22 Apr 2014
    LHC10b&c (pass2):
    ==================

    Default cut which is currently recommended: AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE/kFALSE, 1)
    Important is the second argument (=1) which replaces the cut on 70 clusters with a crossed rows cuts.  
    !!! Please note, that a cut on 70 clusters is strongly discouraged in LHC10b&c pass2 data analysis!!! 
    Changing to number of clusters (=0) and variations of the cut to 60 or 80 should be included in the systematic studies

    LHC10deh (pass2):
    ==================
    Default cut which is currently recommended: AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE/kFALSE, 0)
    In this period, a cut on 70 clusters should be okay, however, changing to a crossed rows cut and lowering the cut to 60 clusters should be included in the systematic error study.
  */

  AliESDtrackCuts *fcutITSTPC_P = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, clusterCut);
  fcutITSTPC_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  fcutITSTPC_P->SetName("ITSTPC");
  AddPrimaryTrackCut(fcutITSTPC_P);
  AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
  fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  fcutITSSA_P->SetName("ITSSA");
  AddPrimaryTrackCut(fcutITSSA_P);

  return;
}

void AliMultiplicitySelectionCP::AddPrimaryTrackCut(AliESDtrackCuts *cut)
{
  fTrackCutListPrim->Add(cut);
}

Int_t AliMultiplicitySelectionCP::GetNumberOfITSTPCtracks(AliESDEvent *esd)
{
  TArrayI indices;
  return GetNumberOfITSTPCtracks(esd, indices);
}


Bool_t AliMultiplicitySelectionCP::InitV0Daughters(AliESDEvent *esd)
{


  if(fkNtrackMax < esd->GetNumberOfTracks() )
    {
      AliFatal(" fkNtrackMax < esd->GetNumberOfTracks() !!!\n");
    }


  for(Int_t i=0; i< esd->GetNumberOfTracks(); i++)
    {
      fkIsTrackSec[i] = kFALSE;
    }

  //  if(!fkIgnoreV0s) return kTRUE;

  Int_t Nv0  = esd->GetNumberOfV0s();

  for(Int_t iv0 = 0; iv0<Nv0; iv0++)
    {
      AliESDv0 *v0 = esd->GetV0(iv0);
      if(!v0) continue;

      fkIsTrackSec[v0->GetPindex()] = kTRUE;
      fkIsTrackSec[v0->GetNindex()] = kTRUE;
    }

  return kTRUE;
}


Int_t AliMultiplicitySelectionCP::GetNumberOfITSTPCtracks(AliESDEvent *esd, TArrayI &indices)
{
  //  return -1000;
  indices.Set(esd->GetNumberOfTracks());
  indices.Reset(-1);


  fIndicesN.Set(esd->GetNumberOfTracks());
  fIndicesN.Reset(-1);
  fIndicesP.Set(esd->GetNumberOfTracks());
  fIndicesP.Reset(-1);

  const AliESDVertex *vtxESD = esd->GetPrimaryVertex();

  Int_t NtracksSel = 0;
  Int_t NpureITStracks = 0;


  Int_t NtracksSelN = 0;
  Int_t NtracksSelP = 0;

  Double_t bfield = esd->GetMagneticField();
  Double_t dca[2], cov[3];


  esd->ConnectTracks();

  if(fkIgnoreV0s) InitV0Daughters(esd); 




  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack* track = esd->GetTrack(iTrack);
      track->SetESDEvent(esd);

      if(track->GetTPCnclsS()>fTPCnclsS) return -1;

      if(!track->PropagateToDCA(vtxESD, bfield, 500., dca, cov))
	continue;

      if(fkIgnoreV0s && fkIsTrackSec[iTrack])
        continue;

      Bool_t isITSpureSA = ((track->GetStatus() & AliESDtrack::kITSpureSA) != 0);
      if(isITSpureSA) 
	{
	  NpureITStracks++;
	  continue;
	}

      if(TMath::Abs(track->Zv() - vtxESD->GetZ())>fTrackDCAz) 
	continue;

      indices.AddAt(iTrack, NtracksSel);

      NtracksSel++;
      
      if(track->GetSign()<0)
	{
	  fIndicesN.AddAt(iTrack, NtracksSelN);
	  NtracksSelN++;
	}
      else if(track->GetSign()>0)
	{
	  fIndicesP.AddAt(iTrack, NtracksSelP);
	  NtracksSelP++;
	}
      

    }

  indices.Set(NtracksSel);

//  printf("NtracksSelN = %d   NtracksSelP = %d   ***************\n",NtracksSelN,NtracksSelP);


  fIndicesN.Set(NtracksSelN);
  fIndicesP.Set(NtracksSelP);


  for(Int_t i = 0; i< NtracksSel; i++)
    {
      AliESDtrack* tr = esd->GetTrack(indices.At(i));
      if(tr->Eta() < fTrackEtaMin || tr->Eta() > fTrackEtaMax)
	return -2;
      if(!AcceptTrack(tr, kTRUE))
	return -3;
    }

 

  const AliMultiplicity *mult = esd->GetMultiplicity();

  if(NpureITStracks>NtracksSel || mult->GetNumberOfTracklets() > NtracksSel)
    return -4;


  if(!TestFiredChips(esd, indices))
    return -5;


  if(fkCheckReferenceMultiplicity)
    {
      Int_t NRefMult = AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTrackletsITSTPC, 3);

      if(NRefMult > NtracksSel)
	return -6;
    }

  return NtracksSel; 

}



Bool_t AliMultiplicitySelectionCP::AcceptTrack(AliESDtrack *track, Bool_t asPrimary)
{
  if(asPrimary)
    {
      TIter next(fTrackCutListPrim);
      AliESDtrackCuts *cut;
      while ((cut=(AliESDtrackCuts*)next()))
	{
	  if(cut->AcceptTrack(track))
	    return kTRUE;
	}
      return kFALSE;
    }

  else{
    Bool_t isITSrefit = ((track->GetStatus() & AliESDtrack::kITSrefit) != 0);
    Bool_t isTPCrefit = ((track->GetStatus() & AliESDtrack::kTPCrefit) != 0);

    if(isITSrefit || isTPCrefit) return kTRUE;
    else return kFALSE;
  }

  return kFALSE;
}


Bool_t AliMultiplicitySelectionCP::IsTrackSelected(Int_t index)
{

  for(Int_t i = 0; i< fIndicesN.GetSize(); i++)
    {
      if(fIndicesN.At(i)==index) 
	return kTRUE;
    }

  for(Int_t i = 0; i< fIndicesP.GetSize(); i++)
    {
      if(fIndicesP.At(i)==index) 
	return kTRUE;
    }

  return  kFALSE;
}


Bool_t AliMultiplicitySelectionCP::TestFiredChips(AliESDEvent *esd, TArrayI indices)
{
  const AliMultiplicity *mult = esd->GetMultiplicity();
  Int_t Ntracks = indices.GetSize();
  UInt_t *Modules = new UInt_t[2*Ntracks];

  for(Int_t iT = 0; iT< Ntracks; iT++)
    {

//      printf("AliMultiplicitySelectionCP::TestFiredChips:  indices.At(%d) = %d \n", iT, indices.At(iT));

      Int_t statusLay;
      Int_t idet = -1;
      Float_t xloc,zloc;
      AliESDtrack* track = esd->GetTrack(indices.At(iT));
      Bool_t retc=track->GetITSModuleIndexInfo(0,idet,statusLay,xloc,zloc);
      if(retc && statusLay!=5) Modules[2*iT] = idet;
      retc=track->GetITSModuleIndexInfo(1,idet,statusLay,xloc,zloc);
      if(retc && statusLay!=5) Modules[2*iT+1] = idet;
    }

  UInt_t eq, hs, chip;
  for (Int_t i=0; i<1200; i++)
    {
      if (!mult->TestFiredChipMap(i)) continue;
      AliSPDUtils::GetOnlineFromOfflineChipKey(i, eq, hs,  chip);
      UInt_t module = AliSPDUtils::GetOfflineModuleFromOnline(eq, hs, chip);

      Bool_t ktmp = kFALSE;
      for(Int_t iM = 0; iM<2*Ntracks; iM++)
	{
	  if(Modules[iM]==module)
	    ktmp=kTRUE;
	}
      if(!ktmp) 
	{
	  delete[] Modules;
	  return kFALSE;
	}
    }

  delete[] Modules;
  return kTRUE;
}
