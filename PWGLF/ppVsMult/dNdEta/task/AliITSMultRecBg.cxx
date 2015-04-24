#include "AliITSMultRecBg.h"
#include "AliGeomManager.h"
#include "AliMultiplicity.h"
#include "AliITSgeomTGeo.h"
#include <TH2F.h>
#include <TTree.h>
#include <TRandom.h>
#include <TBits.h>
#include <TClonesArray.h>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for generating combinatorial backgroung                             //
// for the SPD tracklets                                                     //
//                                                                           //
// Modes for "insertion", "rotation" and "mixing" are supported              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


ClassImp(AliITSMultRecBg)

//_________________________________________________________________
AliITSMultRecBg::AliITSMultRecBg() 
: fRecType(kData),
  fInjLr(0),
  fInjStave(0),
  fInjModule(0),
  fInjModInStave(0),
  fInjX(0),
  fInjZ(0),
  fInjHitsN(0),
  fInjScale(1),
  fInjCurrTrial(0),
  fInjCluster(),
  fInjBuffer(20*kTrNPar)
  //
{
  // constructor
  fCreateClustersCopy = kTRUE;
  for (int i=0;i<2;i++) {
    fAssociations[i] = 0;
    fInjSPDOcc[i]    = 0;
    fInjSPDPatt[i]   = 0;
    fInjNTrials[i]   = 0;
    fInjNSuccess[i]  = 0;    
    for (int s=4;s--;) fInjModuleClStart[i][s] = fInjModuleClN[i][s] = 0;
  }
  //
}

//_________________________________________________________________
/*
AliITSMultRecBg::AliITSMultRecBg(const AliITSMultRecBg &src) 
  : AliITSMultReconstructor(src),
  fRecType(kData),
  fInjLr(0),
  fInjStave(0),
  fInjModule(0),
  fInjModInStave(0),
  fInjX(0),
  fInjZ(0),
  fInjHitsN(0),
  fInjScale(1),
  fInjCurrTrial(0),
  fInjCluster(),
  fInjBuffer(20*kTrNPar)
  //
{
  // dummy copy c-tor
  fCreateClustersCopy = kTRUE;
  for (int i=0;i<2;i++) {
    fAssociations[i] = 0;
    fInjSPDOcc[i]    = 0;
    fInjSPDPatt[i]   = 0;
    fInjNTrials[i]   = 0;
    fInjNSuccess[i]  = 0;    
    for (int s=4;s--;) fInjModuleClStart[i][s] = fInjModuleClN[i][s] = 0;
  }
  //
}
*/
//_________________________________________________________________
AliITSMultRecBg::~AliITSMultRecBg()
{
  // destructor
  for (int i=0;i<2;i++) {
    delete[] fAssociations[i];
    delete   fInjSPDOcc[i];
    delete   fInjSPDPatt[i];
    for (int s=4;s--;) {
      delete[] fInjModuleClStart[i][s];
      delete[] fInjModuleClN[i][s];
    }
  }
  //
}

//____________________________________________________________________
void AliITSMultRecBg::CreateMultiplicityObject()
{
  // create AliMultiplicity object
  //
  if (fRecType==kData || fRecType==kBgRot || fRecType==kBgMix) {
    AliITSMultReconstructor::CreateMultiplicityObject();
  }
  else if (fRecType==kBgInj) {
    if (fMult) delete fMult;
    int ntr = GetNInjSuccsess();
    TBits fastOrFiredMap(0);
    fMult = new AliMultiplicity(ntr,0,fNFiredChips[0],fNFiredChips[1],fastOrFiredMap);
    fMult->SetMultTrackRefs(kTRUE);
    fMult->SetScaleDThetaBySin2T(fScaleDTBySin2T);
    fMult->SetDPhiWindow2(fDPhiWindow2);
    fMult->SetDThetaWindow2(fDThetaWindow2);
    fMult->SetDPhiShift(fDPhiShift);
    fMult->SetNStdDev(fNStdDev);
    //
    for (int itr=ntr;itr--;)  {
      Float_t *tlInfo = fInjBuffer.GetArray() + itr*kTrNPar;
      fMult->SetTrackletData(itr,tlInfo);
    }
  }
  //

}

//_________________________________________________________________________
void AliITSMultRecBg::Run(TTree* tree, Float_t* vtx, TTree* treeMix)
{
  // reconstruct with current settings
  if (!tree) AliFatal("The RP Tree is missing");
  if (fRecType==kBgMix && !treeMix) AliFatal("Mixed Mode requested but 2nd RP Tree is missing");
  if (!vtx)  return;
  //
  if      (fRecType==kData)  Reconstruct(tree, vtx);
  else if (fRecType==kBgRot) Reconstruct(tree, vtx);
  else if (fRecType==kBgMix) ReconstructMix(tree, treeMix, vtx);
  else if (fRecType==kBgInj) GenInjBgSample(tree,vtx); // if needed, the reco will be done internally
  else {  AliError(Form("Unknown running mode %d",fRecType)); }
  //
  CreateMultiplicityObject();
}

//_________________________________________________________________________
Bool_t AliITSMultRecBg::PrepareInjBgGenerator(Float_t *vtx)
{
  // prepare histograms/patterns for bg generation
  //
  if (!fRecoDone) Reconstruct(fTreeRP,vtx);
  //
  float xshift = fSPDSeg.Dx()/2, zshift = fSPDSeg.Dz()/2;
  int  maxRow1 = fSPDSeg.Npx(), maxCol1 = fSPDSeg.Npz();
  int nChipsPerLadder = fSPDSeg.GetNumberOfChips();  
  //
  fInjCluster.SetLabel(kInjFakeLabel,0);
  fInjCluster.SetLabel(-2,1);
  fInjCluster.SetLabel(-2,2);
  //
  for (int ilr=0;ilr<2;ilr++) {
    int nLadPerStave = AliITSgeomTGeo::GetNDetectors(ilr+1);
    int nStaves = AliITSgeomTGeo::GetNLadders(ilr+1);
    for (int ild=nLadPerStave;ild--;) {
      fInjModuleTotClN[ilr][ild] = 0;
      for (int is=nStaves;is--;) {
	fInjModuleClStart[ilr][ild][is] = -1;
	fInjModuleClN[ilr][ild][is] = 0;
      }
    }
    //
    fInjSPDOcc[ilr]->Reset();
    fInjSPDPatt[ilr]->ResetAllBits();
    TClonesArray* clusters = fClArr[ilr];
    Int_t nclus = fNClustersLay[ilr];
    //
    for (int icl=0;icl<nclus;icl++) {
      AliITSRecPoint* clus = (AliITSRecPoint*)clusters->UncheckedAt(icl); 
      if (!clus) continue;
      int ladder = clus->GetDetectorIndex();   // ladder id within the layer
      int stave = ladder/nLadPerStave;
      ladder = nLadPerStave-1 - ladder%nLadPerStave; // ladder id within the stave !!!! invert to get human readble histos
      // the clusters are packed per modules, register the beginning and n.clusters of each module
      if (fInjModuleClStart[ilr][ladder][stave]<0)  fInjModuleClStart[ilr][ladder][stave] = icl;
      fInjModuleClN[ilr][ladder][stave]++;
      fInjModuleTotClN[ilr][ladder]++; 
      //
      int chip = fSPDSeg.GetChipFromLocal(clus->GetDetLocalX(),clus->GetDetLocalZ()); // chip within the ladder
      if (ilr==1) chip = nChipsPerLadder - 1 - chip; //!!!! invert to get human readble histos
      chip += ladder*nChipsPerLadder; // chip within the stave
      ladder %= nLadPerStave; // ladder id within the stave
      fInjSPDOcc[ilr]->Fill(chip, stave);                      // register cluster for hit density profile
      //
      float xloc = clus->GetDetLocalX()*1e4+xshift;
      float zloc = clus->GetDetLocalZ()*1e4+zshift;
      int row,col;
      fSPDSeg.GetPadIxz(xloc,zloc,row,col); // row,col here stats from 1
      row--;
      col--;
      //
      // generate bit pattern according to hit type
      int npix = clus->GetType();
      int nrows = clus->GetNy();
      int ncols = clus->GetNz();
      float cx=0,cz=0;
      UInt_t *patt = GenClusterPattern(npix,nrows,ncols,cx,cz);
      row = int(row - cx); if (row<0) row=0; else if (row+nrows>maxRow1) row = maxRow1-nrows;
      col = int(col - cz); if (col<0) col=0; else if (col+ncols>maxCol1) col = maxCol1-ncols;
      for (int icol=ncols;icol--;) {
	UInt_t hcol = patt[icol];
	for (int irow=nrows;irow--;) {
	  if (!(hcol&(1<<irow))) continue; // no hit here
	  int pbit = GetPixelBitL(stave,ladder,col+icol,row+irow);
	  fInjSPDPatt[ilr]->SetBitNumber( pbit );
	}
      }
      //
    } // loop over clusters of layer ilr  
  } // loop over layers
  //
  if (fNClustersLay[0]==0||fNClustersLay[1]==0) {
    AliInfo(Form("Trackleting is impossible: Ncl1=%d Ncl2=%d",fNClustersLay[0],fNClustersLay[1]));
    return kFALSE;
  }
  for (int i=0;i<2;i++) {
    if (fAssociations[i]) delete[] fAssociations[i];
    int* tmpi = fAssociations[i] = new Int_t[ fNClustersLay[i] ];
    for (int ic=fNClustersLay[i];ic--;) tmpi[ic] = -1;
  }
  for (int itr=fNTracklets;itr--;) {
    Float_t* tracklet = GetTracklet(itr);
    int cll1 = (int)tracklet[kClID1];
    int cll2 = (int)tracklet[kClID2];
    fAssociations[0][cll1] = cll2;
    fAssociations[1][cll2] = cll1;
  }
  //
  fInjBuffer.Set(fNTracklets*kTrNPar);
  //
  return kTRUE;
  //
}

//_________________________________________________________________________
Bool_t AliITSMultRecBg::GenClusterToInject()
{
  // generate bg cluster on layer lr
  //
  int nLadPerStave = AliITSgeomTGeo::GetNDetectors(fInjLr+1);
  int nChipsPerModule = fSPDSeg.GetNumberOfChips();
  //
  int clID;
  //
  //RRR  printf("On Lr %d | %d %d\n",fInjLr,fNClustersLay[0],fNClustersLay[1]);
  do {
    fInjSPDOcc[fInjLr]->GetRandom2(fInjZ, fInjX); 
    //    printf("raw: %f %f\n",fInjZ,fInjX);
    fInjStave = int(fInjX);
    int chipInStave  = int(fInjZ);    // chip in the stave
    //
    fInjX = (fInjX - fInjStave)*fSPDSeg.Dx();     // local x coordinate
    fInjModInStave = chipInStave/nChipsPerModule;
    fInjModInStave = nLadPerStave - 1 - fInjModInStave;  //!!! go from human readable to formal one
    fInjModule = fInjStave*nLadPerStave + fInjModInStave; // formal module number
    //
    fInjZ = (fInjZ-chipInStave)*fSPDSeg.Dz();     // local z coordinate
    //    printf("Z %e X %e | MinSt%d Mod%d\n",fInjZ, fInjX,fInjModInStave,fInjModule);
    //
    clID = PickClusterPrototype(fInjLr, fInjModInStave, fInjStave);    
  } while(clID<0);
  //
  //  printf("clID: %d %d %d %d\n",clID,fNClustersLay[0],fNClustersLay[1],fInjLr);
  AliITSRecPoint* rClus = (AliITSRecPoint*)fClArr[fInjLr]->UncheckedAt(clID);
  fInjCluster.SetLayer(fInjLr);
  fInjCluster.SetType(TMath::Min(kInjMaxNY*kInjMaxNZ,rClus->GetType()));
  fInjCluster.SetNy(TMath::Min(kInjMaxNY,rClus->GetNy()));
  fInjCluster.SetNz(TMath::Min(kInjMaxNZ,rClus->GetNz()));
  fInjCluster.SetDetectorIndex(fInjModule);
  fInjCluster.SetVolumeId(AliGeomManager::LayerToVolUID(AliGeomManager::kSPD1+fInjLr,fInjModule));
  //
  PlaceInjCluster();
  return kTRUE;
  //
}

//_________________________________________________________________________
void AliITSMultRecBg::PlaceInjCluster()
{
  // place the injected cluster on the selected module, 
  // avoiding overlaps with real clusters
  int npix = fInjCluster.GetType(), nrows = fInjCluster.GetNy(), ncols = fInjCluster.GetNz();
  Float_t cx=0,cz=0;
  UInt_t* pattern = GenClusterPattern(npix, nrows, ncols, cx,cz);
  //
  // try to embedd on top of real clusters
  int maxRow1 = fSPDSeg.Npx(), maxCol1 = fSPDSeg.Npz();
  int maxRow  = maxRow1-1,     maxCol  = maxCol1-1;
  TBits &bits = *fInjSPDPatt[fInjLr];   // hits pattern of selected layer
  int row0=0,col0=0;
  Bool_t failed;
  do {
    failed = kFALSE;
    fSPDSeg.GetPadIxz(fInjX,fInjZ,row0,col0);
    row0--; col0--;                  // row,col here start from 1
    if (row0+nrows > maxRow1) row0 = maxRow1-nrows;
    if (col0+ncols > maxCol1) col0 = maxCol1-ncols;
    //
    //    printf("Cluster at %f %f: col:%d row:%d\n",fInjZ,fInjX,col0,row0);
    // 
    // check if generated pattern is mergable with data clusters
    fInjHitsN = 0;
    for (int ic=ncols;ic--;) {
      int colt = col0+ic; 
      for (int ir=nrows;ir--;) {
	if ( !(pattern[ic]&(1<<ir)) ) continue; // not used
	int rowt = row0+ir;
	int bitID = GetPixelBitL(fInjStave,fInjModInStave,colt,rowt);
	if ( bits.TestBitNumber(bitID) ||	// is pixel free?
	     (colt>0      && bits.TestBitNumber(bitID-maxRow1)) || // pixel  1 column below
	     (colt<maxCol && bits.TestBitNumber(bitID+maxRow1)) || // pixel  1 column above
	     (rowt>0      && bits.TestBitNumber(bitID-1))       || // pixel  1 row below
	     (rowt<maxRow && bits.TestBitNumber(bitID+1)))         // pixel in 1 row above
	  {failed=kTRUE; break;}
	// ok for this pixel
	fInjHits[fInjHitsN++] = bitID;
      }
      if (failed) break;
    }
    if (failed) { // generate new x,z
      //      printf("Conflict found, retry\n");
      fInjX = gRandom->Rndm()*fSPDSeg.Dx();
      fInjZ = gRandom->Rndm()*fSPDSeg.Dz();
      continue;
    }
  } while(failed);
  //
  // caclulate cluster coordinates
  float x=0,z=0;
  fInjX = fInjZ = 0;
  for (int pix=fInjHitsN;pix--;) {
    ExtractPixelL(fInjHits[pix], fInjStave, fInjModInStave, col0, row0);
    fSPDSeg.GetPadCxz(row0+1,col0,x,z); // !!! Note: here row starts from 1 but col from 0!!!
    fInjX += x;
    fInjZ += z;
  }
  fInjX = (fInjX/fInjHitsN-fSPDSeg.Dx()/2)*1e-4;
  fInjZ = (fInjZ/fInjHitsN-fSPDSeg.Dz()/2)*1e-4;
  const TGeoHMatrix *mT2L = AliITSgeomTGeo::GetTracking2LocalMatrix(fInjLr+1,fInjStave+1,fInjModInStave+1);
  Double_t loc[3]={fInjX,0.,fInjZ},trk[3]={0.,0.,0.};
  mT2L->MasterToLocal(loc,trk);
  //
  fInjCluster.SetX(0);
  fInjCluster.SetY(trk[1]);
  fInjCluster.SetZ(trk[2]);
  //
  // printf("ClCoord: %+e %+e %+e\n",fInjCluster.GetX(),fInjCluster.GetY(),fInjCluster.GetZ());
}

//_________________________________________________________________________
void AliITSMultRecBg::InitInjBg()
{
  // initialize data for bg injection
  char buff[100];
  for (int ilr=0;ilr<2;ilr++) {
    sprintf(buff,"occL%d",ilr);
    int nLaddersStave = AliITSgeomTGeo::GetNDetectors(ilr+1);
    int nChipsStave   = fSPDSeg.GetNumberOfChips()*nLaddersStave;
    int nStaves       = AliITSgeomTGeo::GetNLadders(ilr+1);
    int nColsStave    = fSPDSeg.Npz();
    int nRowsStave    = fSPDSeg.Npz();
    fInjSPDOcc[ilr]  = new TH2F(buff,buff,nChipsStave,0,nChipsStave, nStaves,0,nStaves);
    fInjSPDPatt[ilr] = new TBits(nStaves*nLaddersStave*nColsStave*nRowsStave);
    for (int is=0;is<nStaves;is++) {
      for (int il=0;il<nLaddersStave;il++) {
	fInjModuleClStart[ilr][il] = new Int_t[nStaves];
	fInjModuleClN[ilr][il]     = new Int_t[nStaves];
      }
    }
    //
  }
}

//_________________________________________________________________________
UInt_t* AliITSMultRecBg::GenClusterPattern(Int_t &npix, Int_t &ny, Int_t &nz, Float_t cy,Float_t &cz)
{
  // generate random pattern for the cluster
  static UInt_t hitPattern[160];
  if (ny>kInjMaxNY) ny = kInjMaxNY;
  if (nz>kInjMaxNZ) nz = kInjMaxNZ;
  for (int iz=nz;iz--;) hitPattern[iz]=0;
  //
  // handle explicitly easy cases: () is for pixels in the same column ...
  if      (npix==1) hitPattern[0] = 0x1; // type (|)
  else if (npix==2) {
    if (nz==1) hitPattern[0] = 0x3;      // tpye (||)
    else       hitPattern[0] = hitPattern[1] = 0x1; // tpye (|)(|)
  }
  else if (npix==3) {
    if      (nz==1) hitPattern[0] = 0x7;      // tpye (|||)
    else if (ny==1) hitPattern[0] = hitPattern[1] = hitPattern[2] = 0x1; // type (|)(|)(|)
    else { hitPattern[0] = 0x3; hitPattern[1] = 0x1;}    // type (||)(|)
  }
  else if (npix==4) {
    if      (nz==1) hitPattern[0] = 0xf;          // type (||||)
    else if (ny==1) hitPattern[0] = hitPattern[1] = hitPattern[2] = hitPattern[3] = 0x1; // type (|)(|)(|)(|)
    else if (ny==2) {
      if (nz==2) hitPattern[0] = hitPattern[1] = 0x3; // tpye (||)(||)
      else { 
	hitPattern[0] = hitPattern[1] = hitPattern[2] = 0x1; // type (||)(|)(|) or (|)(||)(|)
	hitPattern[gRandom->Rndm()>0.5 ? 0 : 1] |= 0x2; 
      } 
    }
    else { 
      hitPattern[0] = 0x7; 
      hitPattern[1] = (gRandom->Rndm()>0.8) ? 0x1 : 0x2;  // type (|||) + (_|_) or (|||) + (|__) 
    }
  }
  // more complex topologies
  else {
    UInt_t mask = 0xffffffff>>(32-ny);
    for (int i=nz;i--;) hitPattern[i] = mask;
    int toSup = ny*nz - npix;
    int trial = toSup*5;
    int colToTouch = nz-1; // at least 1 column should not be touched
    while(toSup>0 && (trial--)>0) { // suppress random pixel until needed npix, ny and nz is obtained
      int col = gRandom->Integer(nz);
      if (hitPattern[col]==0x1) continue; // don't lose this column
      if (hitPattern[col]==mask) {
	if (!colToTouch) continue; // this is the only remaining column with ny rows hit
	colToTouch--; // there are other columns of ny rows hit, may touch it
      } 
      hitPattern[col] >>= 1;
      toSup--;
    }
    if (toSup) npix += toSup; // failed to suppress exact number of pixels 
  }
  //
  cy = cz = 0; // get centroid
  for (int col=nz;col--;) {
    int npxr = 0;
    for (int row=ny;row--;) if (hitPattern[col] & (0x1<<row)) {npxr++; cy+=(1+row);}
    cz += npxr*(1+col);      
  }
  cz = cz/npix - 0.5;
  cy = cy/npix - 0.5;
  //
  return (UInt_t*)&hitPattern[0];
}

//_________________________________________________________________________
Int_t AliITSMultRecBg::PickClusterPrototype(Int_t lr, Int_t ladInStave, Int_t stave2avoid)
{
  // pick random cluster to use as a prototype for injection. It should come 
  // from ladders at given Z belt and should not be from the stave2avoid
  // pick random cluster to use as a prototype for injection. It should come 
  // from ladders at given Z belt and should not be from the stave2avoid
  static int tried = 0;
  static int ladInStaveOrig = 0;
  //
  int ncl = fInjModuleTotClN[lr][ladInStave];
  if (stave2avoid>=0) ncl -= fInjModuleClN[lr][ladInStave][stave2avoid];
  if (ncl<1) {
    int totLad = AliITSgeomTGeo::GetNDetectors(lr+1);
    if (!tried) ladInStaveOrig = ladInStave; // before starting the resque solution, save the orig.ladder
    if (++tried >= totLad) { // failed to find cluster avoiding a stave2avoid
      tried = 0;
      return PickClusterPrototype(lr,ladInStaveOrig,-1);
    }
    int useLad = ladInStave+1;
    if (useLad>=AliITSgeomTGeo::GetNDetectors(lr+1)) useLad = 0;
    return PickClusterPrototype(lr,useLad,stave2avoid); // look in the neigbouring ladder
  }
  //
  int pick = gRandom->Integer(ncl);
  int nst = AliITSgeomTGeo::GetNLadders(lr+1);
  int stave = 0;
  for (stave=0;stave<nst;stave++) {
    if (stave==stave2avoid) continue;
    if (pick<fInjModuleClN[lr][ladInStave][stave]) break;
    pick -= fInjModuleClN[lr][ladInStave][stave];
  }
  //
  tried = 0;
  return fInjModuleClStart[lr][ladInStave][stave]+pick;
}

//_________________________________________________________________________
Int_t AliITSMultRecBg::SearchInjTracklet(const Float_t *vtx)
{
  // reconstruct tracklets which may involve injected cluster
  // fake arrays to be used for injected cluster in MultReco
  Float_t clustersLayInj[kClNPar];
  //  Int_t   detectorIndexClustersLayInj[1];
  //  Bool_t  overlapFlagClustersLayInj[1];
  //  UInt_t  usedClusLayInj[1];
  //
  Bool_t kUseOrig = kFALSE;//kTRUE;
  if (kUseOrig) {
    // to try with unused clusters
    if ( fAssociations[fInjLr][fInjCurrTrial]>=0 ) return 0; // associated
    float *origCl = GetClusterOfLayer(fInjLr,fInjCurrTrial);
    for (int i=0;i<kClNPar;i++) clustersLayInj[i] = origCl[i];
  }
  else {
    // >> fill cluster data: equivavlent of data fetched in AliITSMultReconstructor::LoadClusterArrays
    //    detectorIndexClustersLayInj[0] = fInjCluster.GetDetectorIndex();
    //    overlapFlagClustersLayInj[0]   = kFALSE;
    //    usedClusLayInj[0]              = 0;
    fInjCluster.GetGlobalXYZ( clustersLayInj );
    for (int i=3;i--;) clustersLayInj[kClMC0+i] = fInjCluster.GetLabel(i);
    ClusterPos2Angles(clustersLayInj, vtx); // convert to angles
  }
  //
  // compare injected cluster with real ones
  int partnerLr = 1 - fInjLr;

  double bestChiUsed = fNStdDev*2, bestChiFree = fNStdDev*2;
  int bestPartnerUsed=-1,bestPartnerFree=-1;
  //
  for (int icl=fNClustersLay[partnerLr];icl--;) {
    Float_t *partnerCl = GetClusterOfLayer(partnerLr,icl);
    Double_t dTheta = clustersLayInj[kClTh] - partnerCl[kClTh]; 
    Double_t dPhi   = clustersLayInj[kClPh] - partnerCl[kClPh];
    if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;     // take into account boundary condition
    Float_t d = CalcDist(dPhi,fInjLr==0 ? -dTheta:dTheta, fInjLr==0 ? clustersLayInj[kClTh]:partnerCl[kClTh]);
    if (d>fNStdDev) continue;
    //
    int competitor = fAssociations[partnerLr][icl]; // is the cluster of partner layer already used by some tracklet?
    if (competitor>=0) { if (d<bestChiUsed && d < fMinDists[ partnerLr==1 ? icl : competitor ]) { bestChiUsed = d; bestPartnerUsed = icl; } }
    else               { if (d<bestChiFree) { bestChiFree = d; bestPartnerFree = icl; } }
  }
  //
  int winner = -1;
  //
  if (bestPartnerUsed>=0) {
    if (bestPartnerFree>=0) { winner = bestChiFree<bestChiUsed ? bestPartnerFree : bestPartnerUsed;}  // shall we subtract old real tracklet if the winner cluster is used one?
    else { winner = bestPartnerUsed; }
  }
  else winner = bestPartnerFree;

  //  printf("%d\n",winner);

  if (winner<0) return 0;
  //
  int nCurrTr = GetNInjSuccsess();
  if (nCurrTr >= fInjBuffer.GetSize()/kTrNPar) fInjBuffer.Set((20+nCurrTr*2)*kTrNPar);
  //
  Float_t *tracklet = fInjBuffer.GetArray() + nCurrTr*kTrNPar;
  Float_t *clPar1=0,*clPar2=0;
  //  AliInfo(Form("Size: %d NCurrTr: %d El: %d Winner: %d trackler: %p",fInjBuffer.GetSize(), nCurrTr, nCurrTr*kTrNPar, winner,tracklet));
  if (fInjLr) {
    clPar1 = GetClusterOfLayer(partnerLr,winner);
    clPar2 = clustersLayInj;
    tracklet[kClID1] = winner;
    tracklet[kClID2] = -1;
  }
  else {
    clPar2 = GetClusterOfLayer(partnerLr,winner);
    clPar1 = clustersLayInj;
    tracklet[kClID1] = -1;
    tracklet[kClID2] =  winner;
  }
  tracklet[kTrTheta] = clPar1[kClTh];    // use the theta from the clusters in the first layer
  tracklet[kTrPhi]   = clPar1[kClPh];    // use the phi from the clusters in the first layer
  tracklet[kTrDPhi]  = clPar1[kClPh] - clPar2[kClPh];  // store the difference between phi1 and phi2 
  //
  // define dphi in the range [0,pi] with proper sign (track charge correlated)
  if (tracklet[kTrDPhi] > TMath::Pi())   tracklet[kTrDPhi] = tracklet[kTrDPhi]-2.*TMath::Pi();
  if (tracklet[kTrDPhi] < -TMath::Pi())  tracklet[kTrDPhi] = tracklet[kTrDPhi]+2.*TMath::Pi();
  tracklet[kTrDTheta] = clPar1[kClTh] - clPar2[kClTh]; // store the theta1-theta2
  tracklet[kTrLab1] = clPar1[kClMC0];
  tracklet[kTrLab2] = clPar2[kClMC0];
  //
  //  printf("Got Tracklet from lr%d: %f %f\n",fInjLr,fInjNSuccess[0],fInjNSuccess[1]);
  //
  return 1;
  // 
}

//_________________________________________________________________________
void AliITSMultRecBg::GenInjBgSample(TTree* treeRP, Float_t *vtx)
{
  // generate a sample of tracklets from injected bg
  //
  SetTreeRP(treeRP);
  InitInjBg();
  if (!PrepareInjBgGenerator(vtx)) return;
  //
  fInjNSuccess[0] = fInjNSuccess[1] = 0;
  for (int i=0;i<2;i++) {
    fInjLr = i;
    fInjNTrials[i] = TMath::Max(1,int(fInjScale*fNClustersLay[i]));
    for (fInjCurrTrial=0;fInjCurrTrial<fInjNTrials[i];fInjCurrTrial++) {
      if (!GenClusterToInject()) break;
      fInjNSuccess[i] += SearchInjTracklet(vtx);
    }
  }
  printf("Successes/Trials: %d/%d %d/%d\n",fInjNSuccess[0],fInjNTrials[0],fInjNSuccess[1],fInjNTrials[1]);
  //  
}

