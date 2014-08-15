#include <TMath.h>
#include "AliITSSAPLayer.h"
#include "AliITSSAPAux.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliVertex.h"
#include <TRandom.h>
#include <TStopwatch.h>
#include <TString.h>


using namespace AliITSSAPAux;

//_________________________________________________________________
AliITSSAPLayer::AliITSSAPLayer() :
  fClusters(0)
  ,fLrID(-1)
  ,fVIDOffset(0)
  ,fNClusters(0)
  ,fZMin(0)
  ,fZMax(0)
  ,fDZInv(-1)
  ,fDPhiInv(-1)
  ,fNZBins(20)
  ,fNPhiBins(20)
  ,fQueryZBmin(-1)
  ,fQueryZBmax(-1)
  ,fQueryPhiBmin(-1)
  ,fQueryPhiBmax(-1)
  ,fBins(0)
  ,fOccBins(0)
  ,fNOccBins(0)
  ,fNFoundClusters(0)
  ,fFoundClusterIterator(0)
  ,fFoundBinIterator(0)
  ,fFoundBins(0)
  ,fSortedClInfo(0)
{
  // def. c-tor
}

//_________________________________________________________________
AliITSSAPLayer::AliITSSAPLayer(int id, float zspan,int nzbins,int nphibins, int buffer) :
  fClusters(0)
  ,fLrID(id)
  ,fVIDOffset((id+1)*2048)
  ,fNClusters(0)
  ,fZMin(-zspan)
  ,fZMax(zspan)
  ,fDZInv(-1)
  ,fDPhiInv(-1)
  ,fNZBins(nzbins)
  ,fNPhiBins(nphibins)
  ,fQueryZBmin(-1)
  ,fQueryZBmax(-1)
  ,fQueryPhiBmin(-1)
  ,fQueryPhiBmax(-1)
  ,fBins(0)
  ,fOccBins(0)
  ,fNOccBins(0)
  ,fNFoundClusters(0)
  ,fFoundClusterIterator(0)
  ,fFoundBinIterator(0)
  ,fFoundBins(0)
  ,fSortedClInfo(0)
{
  // c-tor
  Init(buffer);
}

//_________________________________________________________________
AliITSSAPLayer::~AliITSSAPLayer()
{
  // d-tor
  delete[] fBins;
  delete[] fOccBins;
  delete fClusters;
}

//_________________________________________________________________
void AliITSSAPLayer::Init(int buffer)
{
  if (fClusters) {
    printf("Already initialized\n");
    return;
  }
  if (fNZBins<1)   fNZBins = 2;
  if (fNPhiBins<1) fNPhiBins = 1;
  fDZInv   = fNZBins/(fZMax-fZMin);
  fDPhiInv = fNPhiBins/TMath::TwoPi();
  //
  fBins = new ClBinInfo_t[fNZBins*fNPhiBins];
  fOccBins = new int[fNZBins*fNPhiBins];
  if (buffer<100) buffer = 100;
  fClusters = new TObjArray(buffer);
  fSortedClInfo.reserve(buffer);
  //
  // prepare detectors info
  int id1 = fLrID+1;
  Int_t nlad=AliITSgeomTGeo::GetNLadders(id1);
  Int_t ndet=AliITSgeomTGeo::GetNDetectors(id1);
  int detID = 0;
  for (Int_t j=1; j<nlad+1; j++) {
    for (Int_t k=1; k<ndet+1; k++) { //Fill this layer with detectors
      ITSDetInfo_t det;
      det.index = detID++;
      //
      TGeoHMatrix m; AliITSgeomTGeo::GetOrigMatrix(id1,j,k,m);
      const TGeoHMatrix *tm=AliITSgeomTGeo::GetTracking2LocalMatrix(id1,j,k);
      m.Multiply(tm);
      Double_t txyz[3] = {0.}, xyz[3] = {0.};
      m.LocalToMaster(txyz,xyz);
      det.xTF = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
      det.phiTF = TMath::ATan2(xyz[1],xyz[0]);
      //BringTo02Pi(det.phiTF);
      det.sinTF = TMath::Sin(det.phiTF);
      det.cosTF = TMath::Cos(det.phiTF);
      //
      // compute the real radius (with misalignment)
      TGeoHMatrix mmisal(*(AliITSgeomTGeo::GetMatrix(id1,j,k)));
      mmisal.Multiply(tm);
      xyz[0]=0.;xyz[1]=0.;xyz[2]=0.;
      mmisal.LocalToMaster(txyz,xyz);
      det.xTFmisal=TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);      
      //
      fDetectors.push_back(det);
    } // end loop on detectors
  } // end loop on ladders 

}

//_________________________________________________________________
void AliITSSAPLayer::SortClusters(const AliVertex* vtx)
{
  // sort clusters and build fast lookup table
  //
  ClearSortedInfo();
  fSortedClInfo.reserve(fNClusters);
  //
  ClsInfo_t cl;
  for (int icl=fNClusters;icl--;) {
    AliITSRecPoint* cluster = (AliITSRecPoint*)fClusters->UncheckedAt(icl);
    cluster->GetGlobalXYZ( (float*)&cl );
    //
    if (vtx) { // phi and r will be computed wrt vertex
      cl.x -= vtx->GetX();
      cl.y -= vtx->GetY();
    }
    //
    cl.r = TMath::Sqrt(cl.x*cl.x + cl.y*cl.y);
    cl.phi = TMath::ATan2(cl.y,cl.x);
    BringTo02Pi(cl.phi);
    cl.index = icl;
    cl.zphibin = GetBinIndex(GetZBin(cl.z),GetPhiBin(cl.phi));
    cl.detid = cluster->GetVolumeId() - fVIDOffset;
    //
    fSortedClInfo.push_back( cl );
    //
  }
  sort(fSortedClInfo.begin(), fSortedClInfo.end()); // sort in phi, z
  //
  // fill cells in phi,z
  int currBin = -1;
  for (int icl=0;icl<fNClusters;icl++) {
    ClsInfo_t &t = fSortedClInfo[icl]; 
    if (t.zphibin>currBin) { // register new occupied bin
      currBin = t.zphibin;
      fBins[currBin].first = icl;
      fBins[currBin].index = fNOccBins;
      fOccBins[fNOccBins++] = currBin;
    }
    fBins[currBin].ncl++;
  }
  //  Print("clb"); //RS
}

//_________________________________________________________________
void AliITSSAPLayer::Clear(Option_t *)
{
  // clear cluster info
  ClearSortedInfo();
  fNClusters = 0;
  if (fClusters) fClusters->Clear();
  //
}

//_________________________________________________________________
void AliITSSAPLayer::ClearSortedInfo()
{
  // clear cluster info
  fSortedClInfo.clear();
  memset(fBins,0,fNZBins*fNPhiBins*sizeof(ClBinInfo_t));
  memset(fOccBins,0,fNZBins*fNPhiBins*sizeof(int));
  fNOccBins = 0;
}

//_________________________________________________________________
void AliITSSAPLayer::Print(Option_t *opt) const
{
  // dump cluster bins info
  TString opts = opt;
  opts.ToLower();
  printf("Stored %d clusters in %d occupied bins\n",fNClusters,fNOccBins);
  //
  if (opts.Contains("c")) {
    printf("\nCluster info\n");
    for (int i=0;i<fNClusters;i++) {
      const ClsInfo_t &t = fSortedClInfo[i];
      printf("#%5d Bin(phi/z):%03d/%03d Z:%+8.3f Phi:%+6.3f R:%7.3f Ind:%d ",
	     i,t.zphibin/fNZBins,t.zphibin%fNZBins,t.z,t.phi,t.r,t.index);
      if (opts.Contains("l")) { // mc labels
	AliITSRecPoint* rp = (AliITSRecPoint*)fClusters->UncheckedAt(t.index);
	for (int l=0;l<3;l++) if (rp->GetLabel(l)>=0) printf("| %d ",rp->GetLabel(l));
      }
      printf("\n");
    }
  }
  //
  if (opts.Contains("b")) {
    printf("\nBins info (occupied only)\n");
    for (int i=0;i<fNOccBins;i++) {
      printf("%4d %5d(phi/z: %03d/%03d) -> %3d cl from %d\n",i,fOccBins[i],fOccBins[i]/fNZBins,fOccBins[i]%fNZBins,
	     fBins[fOccBins[i]].ncl,fBins[fOccBins[i]].first);
    }
  }
  //
}
 
//_____________________________________________________________
int AliITSSAPLayer::SelectClusters(float zmin,float zmax,float phimin,float phimax)
{
  // prepare occupied bins in the requested region
  //printf("Select: Z %f %f | Phi: %f %f\n",zmin,zmax,phimin,phimax);
  if (!fNOccBins) return 0;
  if (zmax<fZMin || zmin>fZMax || zmin>zmax) return 0;
  fFoundBins.clear();

  fQueryZBmin = GetZBin(zmin);
  if (fQueryZBmin<0) fQueryZBmin = 0;
  fQueryZBmax = GetZBin(zmax);
  if (fQueryZBmax>=fNZBins) fQueryZBmax = fNZBins-1;
  BringTo02Pi(phimin);
  BringTo02Pi(phimax);
  fQueryPhiBmin = GetPhiBin(phimin);
  fQueryPhiBmax = GetPhiBin(phimax);
  int dbz=0;
  fNFoundClusters = 0;
  int nbcheck = fQueryPhiBmax - fQueryPhiBmin + 1;
  if (nbcheck>0) { // no wrapping around 0-2pi, fast case
    for (int ip=fQueryPhiBmin;ip<=fQueryPhiBmax;ip++) {
      int binID = GetBinIndex(fQueryZBmin,ip);
      if ( !(dbz=(fQueryZBmax-fQueryZBmin)) ) { // just one Z bin in the query range 
	ClBinInfo_t& binInfo = fBins[binID];
	if (!binInfo.ncl) continue;
	fNFoundClusters += binInfo.ncl;
	fFoundBins.push_back(binID);
	continue;
      }
      int binMax = binID+dbz;
      for (;binID<=binMax;binID++) {
	ClBinInfo_t& binInfo = fBins[binID];
	if (!binInfo.ncl) continue;
	fNFoundClusters += binInfo.ncl;
	fFoundBins.push_back(binID);
      }      
    }
  }
  else {  // wrapping
    nbcheck += fNPhiBins;
    for (int ip0=0;ip0<=nbcheck;ip0++) {
      int ip = fQueryPhiBmin+ip0;
      if (ip>=fNPhiBins) ip -= fNPhiBins;
      int binID = GetBinIndex(fQueryZBmin,ip);
      if ( !(dbz=(fQueryZBmax-fQueryZBmin)) ) { // just one Z bin in the query range 
	ClBinInfo_t& binInfo = fBins[binID];
	if (!binInfo.ncl) continue;
	fNFoundClusters += binInfo.ncl;
	fFoundBins.push_back(binID);
	continue;
      }
      int binMax = binID+dbz;
      for (;binID<=binMax;binID++) {
	ClBinInfo_t& binInfo = fBins[binID];
	if (!binInfo.ncl) continue;
	fNFoundClusters += binInfo.ncl;
	fFoundBins.push_back(binID);
      }
    }
  }
  fFoundClusterIterator = fFoundBinIterator = 0;
  /*
  //printf("Selected -> %d cl in %d bins\n",fNFoundClusters,(int)fFoundBins.size());
  for (int i=0;i<(int)fFoundBins.size();i++) {
  int bn = fFoundBins[i];
  ClBinInfo_t& bin=fBins[bn];
  printf("#%d b:%d 1st: %3d Ncl:%d\n",i,bn,bin.first,bin.ncl);
  }
  printf("\n");
  */
  return fNFoundClusters;
}

//_____________________________________________________________
int AliITSSAPLayer::GetNextClusterInfoID()
{
  if (fFoundBinIterator<0) return 0;
  int currBin = fFoundBins[fFoundBinIterator];
  if (fFoundClusterIterator<fBins[currBin].ncl) { // same bin
    return fBins[currBin].first+fFoundClusterIterator++;
  }
  if (++fFoundBinIterator<int(fFoundBins.size())) {  // need to change bin
    currBin = fFoundBins[fFoundBinIterator];
    fFoundClusterIterator = 1;
    return fBins[currBin].first;
  }
  fFoundBinIterator = -1;
  return -1;
}

//_____________________________________________________________
void AliITSSAPLayer::ResetFoundIterator()
{
  // prepare for a new loop over found clusters
  if (fNFoundClusters)  fFoundClusterIterator = fFoundBinIterator = 0;
}
