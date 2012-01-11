/**************************************************************************
  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////////
//
//     Data container for relative ITS-TPC alignment analysis
//     Holds an array of AliRelAlignerKalman objects
//     and takes care of merging when processing data in parallel
//
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TTree.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TCollection.h>
#include <AliESDEvent.h>
#include <AliRelAlignerKalman.h>
#include "AliRelAlignerKalmanArray.h"
#include "TList.h"
#include "TBrowser.h"

ClassImp(AliRelAlignerKalmanArray)

//______________________________________________________________________________
AliRelAlignerKalmanArray::AliRelAlignerKalmanArray():
    TNamed("alignerArray","array of aligners"),
    fT0(0),
    fTimebinWidth(0),
    fSize(0),
    fOutRejSigmaOnMerge(10.),
    fOutRejSigmaOnSmooth(1.),
    fAlignerTemplate(),
    fPArray(NULL),
    fListOfGraphs(NULL)
{
  //ctor
}

//______________________________________________________________________________
AliRelAlignerKalmanArray::AliRelAlignerKalmanArray(Int_t t0, Int_t tend, Int_t slotwidth):
    TNamed("alignerArray","array of aligners"),
    fT0(t0),
    fTimebinWidth(slotwidth),
    fSize(0),
    fOutRejSigmaOnMerge(10.),
    fOutRejSigmaOnSmooth(1.),
    fAlignerTemplate(),
    fPArray(NULL),
    fListOfGraphs(new TList)
{
  //ctor
  if (slotwidth==0) fSize = 1;
  else fSize = (tend-t0)/slotwidth;
  fPArray = new AliRelAlignerKalman* [fSize];
  if (fPArray) for (Int_t i=0;i<fSize;i++){fPArray[i]=NULL;}//fill with zeros
  else fSize=0;
  fListOfGraphs->SetName("graphs");
  fListOfGraphs->SetOwner(kTRUE);
}

//______________________________________________________________________________
AliRelAlignerKalmanArray::AliRelAlignerKalmanArray( const AliRelAlignerKalmanArray& in):
    TNamed(in.GetName(), in.GetTitle()),
    fT0(in.fT0),
    fTimebinWidth(in.fTimebinWidth),
    fSize(in.fSize),
    fOutRejSigmaOnMerge(in.fOutRejSigmaOnMerge),
    fOutRejSigmaOnSmooth(in.fOutRejSigmaOnSmooth),
    fAlignerTemplate(in.fAlignerTemplate),
    fPArray(NULL),
    fListOfGraphs(new TList)
{
  //copy ctor
  fPArray = new AliRelAlignerKalman* [fSize];
  if (!fPArray) {fSize=0;return;} //if fail
  for (Int_t i=0;i<fSize;i++)
  {
    if (in.fPArray[i]) fPArray[i]=new AliRelAlignerKalman(*(in.fPArray[i]));
    else fPArray[i]=NULL;
  }
  fListOfGraphs->SetName("graphs");
  fListOfGraphs->SetOwner(kTRUE);
}

//______________________________________________________________________________
AliRelAlignerKalmanArray::~AliRelAlignerKalmanArray()
{
  //dtor
  ClearContents();
  delete [] fPArray;
  delete fListOfGraphs;
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::SetupArray(Int_t t0, Int_t tend, Int_t slotwidth)
{
  //setup array - clears old content
  ClearContents();
  fT0=t0;
  fTimebinWidth=slotwidth;
  Int_t tmpsize;
  if (slotwidth==0) tmpsize = 1;
  else tmpsize = (tend-t0)/slotwidth;
  if (tmpsize!=fSize)
  {
    delete [] fPArray;
    fPArray=new AliRelAlignerKalman* [tmpsize];
    if (fPArray) fSize=tmpsize;
    else fSize=0;
  }
  for (Int_t i=0;i<fSize;i++){fPArray[i]=NULL;}//fill with zeros
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::GetAlignerTemplate()
{
  //get aligner template
  return &fAlignerTemplate;
}

//______________________________________________________________________________
Long64_t AliRelAlignerKalmanArray::Merge( TCollection* list )
{
  //Merge all the arrays
  //tlihe merge is vertical, meaning matching entries in tree are merged

  AliRelAlignerKalmanArray *arrayFromList;
  if (!list) return 0;
  TIter next(list);
  while ( (arrayFromList = dynamic_cast<AliRelAlignerKalmanArray*>(next())) )
  {
    if (fT0 != arrayFromList->fT0) continue;
    if (fTimebinWidth != arrayFromList->fTimebinWidth) continue;
    if (fSize != arrayFromList->fSize) continue;

    for (Int_t i=0; i<GetSize(); i++)
    {
      AliRelAlignerKalman* a1 = fPArray[i];
      AliRelAlignerKalman* a2 = arrayFromList->fPArray[i];
      if (a1 && a2)
      {
        a1->SetRejectOutliers(kFALSE);
        a1->SetOutRejSigma(fOutRejSigmaOnMerge); a1->Merge(a2);
      }
      else
        if (!a1 && a2) fPArray[i] = new AliRelAlignerKalman(*a2);
    }
  }
  fListOfGraphs->Delete();
  return 0;
}

//______________________________________________________________________________
Int_t AliRelAlignerKalmanArray::Timebin( UInt_t timestamp ) const
{
  //calculate binnumber for given timestamp
  if (fTimebinWidth==0) return 0;
  Int_t slot = (timestamp-fT0)/fTimebinWidth; //it's all integers!
  return slot;
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::GetAligner(UInt_t ts)
{
  //get the aligner for specified timestamp
  Int_t tb = Timebin(ts);
  if (tb<0) return NULL;
  if (tb>=fSize) return NULL;
  if (!fPArray) return NULL;
  if (!fPArray[tb])
  {
    fPArray[tb] = new AliRelAlignerKalman( fAlignerTemplate );
    fPArray[tb]->SetTimeStamp(fT0+tb*fTimebinWidth);
  }
  return fPArray[tb];
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::GetAligner(AliESDEvent* event)
{
  //get the aligner for this event and set relevant info
  AliRelAlignerKalman* a = GetAligner(event->GetTimeStamp());
  if (a) a->SetRunNumber(event->GetRunNumber());
  if (a) a->SetMagField(event->GetMagneticField());
  return a;
}

//______________________________________________________________________________
AliRelAlignerKalmanArray& AliRelAlignerKalmanArray::operator=(const AliRelAlignerKalmanArray& in)
{
  //assignment operator
  if (fSize!=in.fSize)
  {
    //if sizes different, delete curent and make a new one with new size
    ClearContents();
    delete [] fPArray;
    fPArray = new AliRelAlignerKalman* [in.fSize];
  }

  fOutRejSigmaOnMerge=in.fOutRejSigmaOnMerge;
  fOutRejSigmaOnSmooth=in.fOutRejSigmaOnSmooth;

  if (!fPArray) fSize=0;
  else fSize = in.fSize;
  fTimebinWidth = in.fTimebinWidth;
  fT0 = in.fT0;
  for (Int_t i=0;i<fSize;i++)
  {
    if (in.fPArray[i]) fPArray[i] = new AliRelAlignerKalman(*(in.fPArray[i]));
    else fPArray[i]=NULL;
  }
  return *this;
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::ClearContents()
{
  //clear contents of array
  for (Int_t i=0;i<fSize;i++)
  {
    if (fPArray[i]) delete fPArray[i];
  }
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::At( Int_t i ) const
{
  //mimic TObjArray::At( Int_t i )
  if (i>=fSize) {printf("AliRelAlignerKalmanArray::At index %i out of bounds, size=%i\n",i,fSize); return NULL;}
  if (i<0) return NULL;
  return fPArray[i];
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::operator[](Int_t i) const
{
  //mimic TObjArray::operator[](Int_t)
  if (i>=fSize) {printf("AliRelAlignerKalmanArray::operator[] index %i out of bounds, size=%i\n",i,fSize); return NULL;}
  if (i<0) return NULL;
  return fPArray[i];
}

//______________________________________________________________________________
AliRelAlignerKalman*& AliRelAlignerKalmanArray::operator[](Int_t i)
{
  //mimic TObjArray::operator[](Int_t) can be used as lvalue
  if (i>=fSize||i<0) {Error("operator[]", "index %i out of bounds, size=%i\n",i,fSize); return fPArray[0];}
  return fPArray[i];
}

//______________________________________________________________________________
AliRelAlignerKalman* AliRelAlignerKalmanArray::Last() const
{
  //return aligner in last filled slot
  if (fSize==0) return NULL;
  return fPArray[fSize-1];
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::Print(Option_t* option) const
{
  // print some info
  TString optionStr(option);
  printf( "%s: t0: %i, tend: %i, width: %i, size: %i, entries: %i\n",
          GetName(),
          fT0, (fT0+(fSize)*fTimebinWidth), fTimebinWidth,
          fSize, GetEntries() );
  if (optionStr.Contains("a"))
    for (Int_t i=0; i<fSize; i++)
    {
      AliRelAlignerKalman* al = fPArray[i];
      if (!al) continue;
      al->Print();
    }
}

//______________________________________________________________________________
Int_t AliRelAlignerKalmanArray::GetEntries() const
{
  //get number of filled slots
  if (!fPArray) return 0;
  Int_t entries=0;
  for (Int_t i=0; i<fSize; i++)
  {
    if (fPArray[i]) entries++;
  }
  return entries;
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::FillTree( TTree* tree ) const
{
  AliRelAlignerKalman* al = NULL;
  tree->Branch("aligner","AliRelAlignerKalman",&al);
  //fill a tree with filled slots
  for (Int_t i=0; i<fSize; i++)
  {
    al = fPArray[i];
    if (al) tree->Fill();
  }
}

//______________________________________________________________________________
TGraphErrors* AliRelAlignerKalmanArray::MakeGraph(Int_t iparam) const
{
  //return a graph for the iparam-th parameter
  if (iparam>8 || iparam<0)
  {
    printf("wrong parameter number. must be from 0-8");
    return NULL;
  }

  Int_t n=GetEntries();
  TVectorD vx(n);
  TVectorD vy(n);
  TVectorD vex(n);
  TVectorD vey(n);
  Int_t entry=0;
  for (Int_t i=0; i<fSize; i++)
  {
    AliRelAlignerKalman* al = fPArray[i];
    if (!al) continue;
    vx(entry) = al->GetTimeStamp();
    vy(entry) = al->GetStateArr()[iparam];
    TMatrixDSym* cm = al->GetStateCov();
    vey(entry) = TMath::Sqrt((*cm)(iparam,iparam));
    entry++;
  }

  TGraphErrors* graph = new TGraphErrors(vx,vy,vex,vey);
  TString graphtitle;
  TString graphtitley;
  switch (iparam)
  {
  case 0:
    graphtitle="rotation \\psi";
    graphtitley="\\psi [rad]";
    break;
  case 1:
    graphtitle="rotation \\theta";
    graphtitley="\\theta [rad]";
    break;
  case 2:
    graphtitle="rotation \\phi";
    graphtitley="\\phi [rad]";
    break;
  case 3:
    graphtitle="shift x";
    graphtitley="x [cm]";
    break;
  case 4:
    graphtitle="shift y";
    graphtitley="y [cm]";
    break;
  case 5:
    graphtitle="shift z";
    graphtitley="z [cm]";
    break;
  case 6:
    graphtitle="TPC vd correction";
    graphtitley="correction factor []";
    break;
  case 7:
    graphtitle="TPC t0 correction";
    graphtitley="t0 correction [\\micros]";
    break;
  case 8:
    graphtitle="TPC dv/dy";
    graphtitley="dv/dy [cm/\\micros/m]";
    break;
  }
  graph->SetName(graphtitle);
  graph->SetTitle(graphtitle);
  TAxis* xas = graph->GetXaxis();
  TAxis* yas = graph->GetYaxis();
  xas->SetTitle("time");
  xas->SetTimeDisplay(1);
  yas->SetTitle(graphtitley);
  return graph;
}

//______________________________________________________________________________
AliRelAlignerKalmanArray* AliRelAlignerKalmanArray::MakeSmoothArray() const
{
  //return a smoothed version of the data
  AliRelAlignerKalmanArray* outputarr = new AliRelAlignerKalmanArray(fT0,
                                        fT0+fSize*fTimebinWidth,fTimebinWidth);

  AliRelAlignerKalman* tmpaligner = outputarr->GetAlignerTemplate();
  tmpaligner->SetOutRejSigma(fOutRejSigmaOnSmooth);
  //copy the first filled slot
  Int_t n=0;
  while (!fPArray[n]) {n++;}
  if (n==fSize) return NULL;
  *tmpaligner   = *fPArray[n];
  if (fPArray[n]->GetNUpdates()>10)
  (*outputarr)[n] = new AliRelAlignerKalman(*(fPArray[n]));
  //for the rest of slots use merge
  for (Int_t i=n+1;i<fSize;i++)
  {
    if (!fPArray[i]) continue;
    PropagateToTime(tmpaligner, fPArray[i]->GetTimeStamp());
    if (tmpaligner->Merge(fPArray[i]))
    (*outputarr)[i] = new AliRelAlignerKalman(*tmpaligner);
    else
    (*outputarr)[i] = new AliRelAlignerKalman(*(fPArray[i]));
  }
  
  tmpaligner->Reset(); //clean up

  return outputarr;
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::PropagateToTime(AliRelAlignerKalman* al, UInt_t timestamp ) const
{
  //propagate an aligner in time
  TMatrixDSym* cov = al->GetStateCov();
  TMatrixDSym corr(TMatrixDSym::kZero,*cov);
  UInt_t dt = (timestamp>al->GetTimeStamp())?
    timestamp-al->GetTimeStamp():al->GetTimeStamp()-timestamp;
  //the propagation matrix
  //to be added to the covariance matrix of kalman filter
  corr(6,6) = dt*1e-8;  //vdrift
  (*cov) += corr;
}

//______________________________________________________________________________
void AliRelAlignerKalmanArray::Browse(TBrowser* b )
{
   if (!b) return;
   if (!fListOfGraphs) 
   {
     fListOfGraphs=new TList();
     fListOfGraphs->SetName("graphs");
     fListOfGraphs->SetOwner(kTRUE);
   }
   for (Int_t i=0;i<9;i++)
   {
     TGraphErrors* gr = dynamic_cast<TGraphErrors*>(fListOfGraphs->At(i));
     if (gr) continue;
     fListOfGraphs->AddAt(MakeGraph(i),i);
   }
   b->Add(fListOfGraphs);
}

