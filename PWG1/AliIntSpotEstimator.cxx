#include "AliIntSpotEstimator.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliLog.h"
#include <TH1.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TStyle.h>


ClassImp(AliIntSpotEstimator)

//______________________________________________________________________________________
AliIntSpotEstimator::AliIntSpotEstimator(const char* name,Double_t outcut,Int_t ntrIP,
					 Int_t nPhiBins,Int_t nestb,
					 Double_t estmin,Double_t estmax,
					 Int_t ntrBins,Int_t ntMn,Int_t ntMx,
					 Int_t nPBins,Double_t pmn,Double_t pmx,Bool_t ntuple) 
: TNamed(name,""),
  fEvProc(0),fIPCenterStat(0),fMinTracksForIP(ntrIP>2?ntrIP:2),fOutlierCut(outcut),fEstimIP(0),
  fEstimVtx(0),fEstimTrc(0),fHVtxXY(0),fNtuple(0),fVertexer(0),fTracks(0)
{
  InitEstimators(nPhiBins,nestb,estmin,estmax,ntrBins,ntMn,ntMx,nPBins,pmn,pmx,ntuple);
}

//______________________________________________________________________________________                                    
AliIntSpotEstimator::AliIntSpotEstimator(Bool_t initDef) 
  : TNamed("IPEstimator",""),
    fEvProc(0),fIPCenterStat(0),fMinTracksForIP(2),fOutlierCut(1e-4),
    fEstimIP(0),fEstimVtx(0),fEstimTrc(0),fHVtxXY(0),fNtuple(0),fVertexer(0),fTracks(0)
{
  if (initDef) InitEstimators();
  //
}

//______________________________________________________________________________________
AliIntSpotEstimator::~AliIntSpotEstimator()
{
  delete fEstimIP;
  delete fEstimVtx;
  delete fEstimTrc;
  delete fVertexer;
  delete fNtuple;
}

//______________________________________________________________________________________
Bool_t AliIntSpotEstimator::ProcessIPCenter(const AliESDVertex* vtx)
{
  // account the vertex in IP center estimation
  if (!vtx) return kFALSE;
  double xyz[3];
  vtx->GetXYZ(xyz);
  double r = xyz[0]*xyz[0] + xyz[1]*xyz[1];
  if (r>2.) return kFALSE;
  for (int i=3;i--;) {
    fIPCenter[i] += xyz[i];
    fIPCen2[i] += xyz[i]*xyz[i];
  }
  fIPCenterStat++;
  fHVtxXY->Fill(xyz[0],xyz[1]);
  return kTRUE;
  //
}

//______________________________________________________________________________________
Bool_t AliIntSpotEstimator::ProcessEvent(const AliESDEvent* esd, const AliESDVertex* vtx)
{
  // Note: this method will modify the indices of the provided vertex
  // Account the vertex in widths estimation
  // if vtx is provided its tracks will be used otherwise use the primary vertex of the event
  if (!IsValid()) {AliError("Not initialized yet"); return -999;}
  if (!fTracks)   fTracks   = new TObjArray(50);
  if (!fVertexer) fVertexer = new AliVertexerTracks();
  //
  fEvProc++;
  int nTracks;
  if (!vtx) vtx = (AliESDVertex*) esd->GetPrimaryVertex();
  if (!vtx || (nTracks=vtx->GetNIndices())<GetMinTracks()+1) return kFALSE;
  //
  UShort_t *trackID = (UShort_t*)vtx->GetIndices(); 
  for (int itr=0;itr<nTracks;itr++) fTracks->Add(esd->GetTrack(trackID[itr]));
  //
  fVertexer->SetFieldkG( esd->GetMagneticField() );
  return ProcessTracks();
  //
}

//______________________________________________________________________________________
Bool_t AliIntSpotEstimator::ProcessEvent(const TObjArray* tracks)
{
  // This method allows to process the vertex made of arbitrary tracks
  // Account the vertex in widths estimation with externally provided tracks 
  // The indices of the vertex must correspond to array of tracks
  //
  if (!IsValid()) {AliError("Not initialized yet"); return -999;}
  if (!fVertexer) fVertexer = new AliVertexerTracks();
  //
  fEvProc++;
  int nTracks = tracks->GetEntriesFast();
  if ( nTracks<GetMinTracks()+1 ) return kFALSE;
  for (int itr=0;itr<nTracks;itr++) fTracks->Add(tracks->At(itr));
  //
  return ProcessTracks();
}

//______________________________________________________________________________________
Bool_t AliIntSpotEstimator::ProcessTracks()
{
  // account the vertex made of provided tracks
  //
  int nTracks = fTracks->GetEntriesFast();
  UShort_t *trackID = new UShort_t[nTracks];
  for (int i=nTracks;i--;) trackID[i] = i;
  double xyzDCA[3],ddca[2],covdca[3];
  int nTracks1 = nTracks - 1;
  AliExternalTrackParam *selTrack,*movTrack=0;
  UShort_t selTrackID,movTrackID=0;
  //
  AliESDVertex* recNewVtx = fVertexer->VertexForSelectedTracks(fTracks,trackID,kTRUE,kFALSE,kFALSE);
  if (!recNewVtx || ((nTracks=recNewVtx->GetNIndices())<GetMinTracks())) {
    if (recNewVtx) delete recNewVtx;
    fTracks->Clear();
    delete[] trackID;
    return kFALSE;
  }
  if (nTracks>=fMinTracksForIP) ProcessIPCenter(recNewVtx);
  //
  double pmn = GetTrackMinP();
  double pmx = GetTrackMaxP();
  double fieldVal = fVertexer->GetFieldkG();
  for (int itr=0;itr<nTracks;itr++) {
    selTrack   = (AliExternalTrackParam*) (*fTracks)[itr]; // track to probe
    double pTrack = selTrack->GetP();
    if (!IsZero(fieldVal) && (pTrack<pmn || pTrack>pmx)) continue;
    selTrackID = trackID[itr];
    //
    if (itr<nTracks1) {
      movTrack   = (AliExternalTrackParam*) (*fTracks)[nTracks1];   // save the track
      movTrackID = trackID[nTracks1];
      (*fTracks)[itr] = movTrack;  
      trackID[itr] = movTrackID;
    }
    fTracks->RemoveAt(nTracks1);     // move the last track in the place of the probed one
    //
    // refit the vertex w/o the probed track
    recNewVtx = fVertexer->VertexForSelectedTracks(fTracks,trackID,kTRUE,kFALSE,kFALSE);
    if (recNewVtx) {
      double told = selTrack->GetX();  // store the original track position
      selTrack->PropagateToDCA(recNewVtx,fieldVal,1e4,ddca,covdca); // in principle, done in the vertexer
      selTrack->GetXYZ(xyzDCA);
      //
      double phiTrack = selTrack->Phi();
      double cs = TMath::Cos(phiTrack);
      double sn = TMath::Sin(phiTrack);
      double trDCA = (xyzDCA[0]-fIPCenIni[0])         *sn - (xyzDCA[1]-fIPCenIni[1])         *cs;  // track signed DCA to origin
      double vtDCA = (recNewVtx->GetXv()-fIPCenIni[0])*sn - (recNewVtx->GetYv()-fIPCenIni[1])*cs;  // vertex signed DCA to origin
      UpdateEstimators(vtDCA,trDCA, nTracks1, pTrack, phiTrack);
      selTrack->PropagateTo(told,fieldVal);    // restore the track
      if (fNtuple) {
	static float ntf[8];
	ntf[0] = float(nTracks1);
	ntf[1] = recNewVtx->GetXv();
	ntf[2] = recNewVtx->GetYv();
	ntf[3] = recNewVtx->GetZv();
	ntf[4] = xyzDCA[0];
	ntf[5] = xyzDCA[1];
	ntf[6] = phiTrack;
	ntf[7] = pTrack;
	fNtuple->Fill(ntf);
      }
    }
    delete recNewVtx;
    // restore the track indices
    (*fTracks)[itr] = selTrack; 
    trackID[itr] = selTrackID;
    if (itr<nTracks1) {
      (*fTracks)[nTracks1] = movTrack; 
      trackID[nTracks1] = movTrackID;
    }
  }
  //
  fTracks->Clear();
  delete[] trackID;
  return kTRUE;
  //
}

//______________________________________________________________________________________
void AliIntSpotEstimator::UpdateEstimators(double rvD, double rtD, double nTracks, double pTrack, double phiTrack)
{
  // update the estimator values
  double estIP  = rvD*rtD;
  double estVtx = rvD*(rvD - rtD);
  double estTrc = rtD*(rtD - rvD);
  //
  if (nTracks >= fMinTracksForIP) fEstimIP->Fill(phiTrack, estIP);
  fEstimVtx->Fill(nTracks, estVtx);
  if (pTrack<1e-6) pTrack = GetTrackMinP()+1e6;
  fEstimTrc->Fill(1./pTrack,estTrc);
  //
}
 
//______________________________________________________________________________________
void AliIntSpotEstimator::InitEstimators(Int_t nPhiBins,Int_t nestb,Double_t estmin,Double_t estmax,
					 Int_t ntrBins,Int_t ntMn,Int_t ntMx,
					 Int_t nPBins,Double_t pmn,Double_t pmx,Bool_t ntuple)
{
  Clear();
  // regularize binning/limits for DCA resolution
  nPBins = nPBins<1 ? 1: nPBins;
  pmn    = pmn>0.1  ? pmn : 0.1;
  pmx    = pmx>pmn  ? pmx : pmn+0.1;
  //
  // regularize binning/limits for vertex resolution
  ntMn = ntMn>2 ? ntMn : 2;
  ntMx = ntMx>ntMn ? ntMx : ntMn;
  ntrBins = ntrBins<1 ? 1:ntrBins;
  int step = (ntMx-ntMn)/ntrBins;
  if (step<1) step = 1;
  ntrBins = (ntMx-ntMn)/step;
  //
  // regularize binning for IP sigma vs angle
  nPhiBins = nPhiBins>1 ? nPhiBins : 1;
  //
  // regularize binning/limits for estimators
  nestb  = nestb>300 ? nestb:300;
  estmin = estmin<-2.e-2 ? estmin : -2.e-2;
  estmax = estmax> 4.e-2 ? estmax :  4.e-2;
  //
  TString nm;
  nm = GetName();  
  nm += "diamondEst";
  fEstimIP  = new TH2F(nm.Data(),nm.Data(),nPhiBins,0.,2.*TMath::Pi(),nestb,estmin,estmax);
  //  
  nm = GetName();  
  nm += "VResEst";  
  fEstimVtx = new TH2F(nm.Data(),nm.Data(),ntrBins,ntMn,ntMx,         nestb,estmin,estmax);
  // 
  nm = GetName();  
  nm += "dcaEst";   
  fEstimTrc = new TH2F(nm.Data(),nm.Data(),nPBins,1./pmx,1./pmn,      nestb,estmin,estmax);
  // 
  nm = GetName();  
  nm += "VtxXY";   
  fHVtxXY = new TH2F(nm.Data(),nm.Data(),200, -1,1, 200,-1,1);
  //
  if (ntuple) {
    nm = GetName();  
    nm += "ntuple";
    fNtuple = new TNtuple(nm.Data(),nm.Data(),"ntrack:xv:yv:zv:xt:yt:phi:p");
  }
  //
  fVertexer = new AliVertexerTracks(); // the field is set dynamically
  fVertexer->SetConstraintOff();
  fTracks = new TObjArray(50);
  //

}

//______________________________________________________________________________________                                                                                             
Double_t AliIntSpotEstimator::GetIPSigma(Int_t bin, Double_t *err) const
{                                                                                                                                                                                    
  // get estimate for the IP sigma
  if (!IsValid()) {AliError("Not initialized yet"); return -999;}
  double cxe,cye;
  double cx  = GetIPCenter(0,&cxe) - GetIPCenIni(0);
  double cy  = GetIPCenter(1,&cye) - GetIPCenIni(1);
  TH1* proj = fEstimIP->ProjectionY("ipProj",bin<1 ? 1:bin, bin<1 ? GetNPhiBins():bin,"e");
  double merr = 0;
  double est = CalcMean(proj, fOutlierCut, &merr) - (cx*cx + cy*cy)/2.;
  if (est>0) {
    est = TMath::Sqrt(est);
    if (err) {
      *err = 0;
      *err = merr*merr;
      *err += cx*cx*cxe*cxe + cy*cy*cye*cye;
      *err = TMath::Sqrt(*err)/est/2.;
    }
  }
  else {
    est = 0;
    if (err) *err = 0;
  }
  delete proj;
  return est;
}

//______________________________________________________________________________________                                                                                             
Double_t AliIntSpotEstimator::GetVtxSigma(int ntr, double* err) const
{                                                                                                                                                                                    
  // get estimate for the IP sigma^2
  if (!IsValid()) {AliError("Not initialized yet"); return -999;}
  int bin = fEstimVtx->GetXaxis()->FindBin(ntr);
  if (bin<1 || bin>GetNTrackBins()) {
    AliError(Form("Requested vertex multiplicity %d out of defined %d-%d range",
		  ntr,GetMinTracks(),GetMaxTracks()));
    return -1;
  }
  TH1* proj = fEstimVtx->ProjectionY("vrProj",bin,bin,"e");
  double est = CalcMean(proj, fOutlierCut, err);
  delete proj;
  if (est>0) {
    est = TMath::Sqrt(est);
    if (err) *err /= 2*est;
  }
  else {
    est = 0;
    if (err) *err = 0;
  }
  return est;
  //
}

//______________________________________________________________________________________                                                                                             
Double_t AliIntSpotEstimator::GetDCASigma(double pt, double *err) const
{                                                                                                                                                                                    
  // get estimate for the IP sigma^2
  if (!IsValid()) {AliError("Not initialized yet"); return -999;}
  if (pt<1e-6) pt = GetTrackMinP()+1e6;
  pt = 1./pt;
  int bin = fEstimTrc->GetXaxis()->FindBin(pt);
  if (bin<1 || bin>GetNPBins()) {
    AliError(Form("Requested track P %.2f out of defined %.2f-%.2f range",1/pt,GetTrackMinP(),GetTrackMaxP()));
    return -1;
  }
  TH1* proj = fEstimTrc->ProjectionY("trProj",bin,bin,"e");
  double est = CalcMean(proj, fOutlierCut, err);
  delete proj;
   if (est>0) {
    est = TMath::Sqrt(est);
    if (err) *err /= 2*est;
  }
  else {
    est = 0;
    if (err) *err = 0;
  }
  return est;
  //
}

//______________________________________________________________________________________                                                                                             
Double_t AliIntSpotEstimator::CalcMean(TH1* histo, Double_t ctfact, Double_t *err)
{
  // calculate mean applying the cut on the outliers
  //
  double max = histo->GetMaximum();
  double cut = (ctfact>0&&ctfact<1.) ? max*ctfact : 0;//TMath::Max(1.0,max*ctfact) : 0;
  int nb = histo->GetNbinsX();
  double mean = 0.,cumul = 0, rms = 0;
  for (int i=1;i<=nb;i++) {
    double vl = histo->GetBinContent(i) - cut;
    if (vl<1e-10) continue;
    double x = histo->GetBinCenter(i);
    mean  += vl*x;
    rms   += vl*x*x;
    cumul += vl;
  }
  //
  mean = cumul>0 ? mean/cumul : 0;
  rms -= mean*mean*cumul;
  if (err) {
    *err = cumul > 1 ? rms/(cumul-1) : 0;
    if (*err>0) *err = TMath::Sqrt(*err/cumul);
  }
  //
  return mean;
}

//______________________________________________________________________________________
void AliIntSpotEstimator::Print(Option_t *) const                                                                                                                              
{   
  if (!IsValid()) {printf("Not initialized yet\n"); return;}
  //
  double cx,cy,cz,cxe,cye,cze;
  cx = GetIPCenter(0,&cxe);
  cy = GetIPCenter(1,&cye);
  cz = GetIPCenter(2,&cze);
  printf("Processed %d events\n",fEvProc);
  printf("Estimator for IP center: %+.4e+-%.3e | %+.4e+-%.3e | %+.4e+-%.3e\n",
	 cx,cxe,cy,cye,cz,cze);
  printf("Initial IP center was  : %+.4e %+.4e %+.4e\n",
	 GetIPCenIni(0),GetIPCenIni(1),GetIPCenIni(2));
  double sgIP = GetIPSigma(0,&cxe);
  printf("Estimator for IP sigma : %.4e+-%.3e\n",sgIP,cxe);
  //
  printf("Estimators for vertex resolution vs Ntracks:\n");
  
  for (int i=1;i<=GetNTrackBins();i++) {
    double sig = GetVtxSigma( (int)fEstimVtx->GetXaxis()->GetBinCenter(i), &cxe );
    if (IsZero(sig)) continue;
    int tmin = TMath::Nint(fEstimVtx->GetXaxis()->GetBinLowEdge(i));
    int tmax = tmin + int(fEstimVtx->GetXaxis()->GetBinWidth(i));
    printf("%3d-%3d : %.4e+-%.3e\n",tmin,tmax,sig,cxe);
  }
  //
  printf("Estimators for track DCA resolution vs P:\n");
  for (int i=1;i<=GetNPBins();i++) {
    double sig = GetDCASigma( 1./fEstimTrc->GetXaxis()->GetBinCenter(i), &cxe );
    if (IsZero(sig)) continue;
    double pmax = 1./fEstimTrc->GetXaxis()->GetBinLowEdge(i);
    double pmin = 1./(fEstimTrc->GetXaxis()->GetBinLowEdge(i)+fEstimTrc->GetXaxis()->GetBinWidth(i));
    printf("%.2f-%.2f : %.4e+-%.3e\n",pmin,pmax,sig, cxe);
  }
}

//______________________________________________________________________________________
void AliIntSpotEstimator::Clear(Option_t *)                                                                                   
{  
  // reset all
  fEvProc = 0;
  fIPCenterStat = 0;
  for (int i=3;i--;) fIPCenter[i] = fIPCenIni[i] = 0.;
  if (fEstimIP)  fEstimIP->Reset();
  if (fEstimVtx) fEstimVtx->Reset();
  if (fEstimTrc) fEstimTrc->Reset();
}

//_____________________________________________________________________
AliIntSpotEstimator &AliIntSpotEstimator::operator += (const AliIntSpotEstimator &src)
{
  fEvProc        += src.fEvProc;
  fIPCenterStat  += src.fIPCenterStat;
  for (int i=3;i--;) fIPCenter[i] += src.fIPCenter[i];
  if (fEstimIP  && src.fEstimIP ) fEstimIP->Add(src.fEstimIP);
  if (fEstimVtx && src.fEstimVtx) fEstimVtx->Add(src.fEstimVtx);
  if (fEstimTrc && src.fEstimTrc) fEstimTrc->Add(src.fEstimTrc);
  if (fHVtxXY   && src.fHVtxXY)   fHVtxXY->Add(src.fHVtxXY);
  //
  return *this;
}

//_____________________________________________________________________
TCanvas* AliIntSpotEstimator::CreateReport(const char* outname)
{
  // prepare canvas with report
  TCanvas *cnv = new TCanvas(GetName(), GetName(),5,5,700,1000);
  gStyle->SetOptStat(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleW(0.7);
  gStyle->SetTitleY(1);
  gStyle->SetTitleX(0.2);
  //
  const Int_t nc=200;
  char buff[nc];
  //
  cnv->Divide(2,2);
  cnv->cd(1);
  //
  TPaveText *pt = new TPaveText(0.05,0.05,0.95,0.95,"blNDC");
  snprintf(buff,nc,"%s | Outliers Cut : %.2e",GetName(),fOutlierCut);
  pt->AddText(buff);
  //
  snprintf(buff,nc,"Processed Events:\n%d",fEvProc);
  pt->AddText(buff);
  //
  snprintf(buff,nc,"Accepted  Events\n%d",fIPCenterStat);
  pt->AddText(buff);
  //
  double cx,cy,cz,cxe,cye,cze;
  cx = GetIPCenter(0,&cxe);
  cy = GetIPCenter(1,&cye);
  cz = GetIPCenter(2,&cze);
  //  
  snprintf(buff,nc,"Int.Spot (#mum)\n%+d#pm%d\n%+d#pm%d\n%+d#pm%d",
	  int(cx*1e4),int(cxe*1e4),int(cy*1e4),int(cye*1e4),int(cz*1e4),int(cze*1e4));
  pt->AddText(buff);
  //
  cx = GetIPSigma(0,&cxe);
  snprintf(buff,nc,"Int.Spot #sigma (#mum):\n%d#pm%d",int(cx*1e4),int(cxe*1e4));
  pt->AddText(buff);
  pt->Draw();
  gPad->Modified();
  //
  cnv->cd(2);
  gPad->SetLeftMargin(0.2);
  TH1* iph = fEstimIP->ProjectionX();
  iph->Reset();
  iph->SetTitle("Int.Spot size vs #phi");
  for (int i=1;i<=iph->GetNbinsX();i++) {
    cx = GetIPSigma(i,&cxe);
    iph->SetBinContent(i,cx*1e4);
    iph->SetBinError(i,cxe*1e4);
  }
  iph->GetXaxis()->SetTitle("#phi");
  iph->GetYaxis()->SetTitle("#sigma_{IP} [#mum]");
  iph->SetMarkerStyle(20);
  iph->Draw("p");
  gPad->Modified();
  iph->SetTitleOffset(2.5,"Y");
  //
  cnv->cd(3);
  gPad->SetLeftMargin(0.2);
  TH1* vrh = fEstimVtx->ProjectionX();
  vrh->Reset();
  vrh->SetTitle("Vertex resolution vs N tracks");
  for (int i=1;i<=vrh->GetNbinsX();i++) {
    cx = GetVtxSigma( TMath::Nint(vrh->GetBinCenter(i)), &cxe);
    vrh->SetBinContent(i,cx*1e4);
    vrh->SetBinError(i,cxe*1e4);
  }
  vrh->GetXaxis()->SetTitle("n tracks");
  vrh->GetYaxis()->SetTitle("#sigma_{VTX} [#mum]");
  vrh->SetMarkerStyle(20);
  vrh->Draw("p");
  gPad->Modified();
  vrh->SetTitleOffset(2.5,"Y");
  //
  cnv->cd(4);
  gPad->SetLeftMargin(0.2);
  TH1* trh = fEstimTrc->ProjectionX();
  trh->Reset();
  trh->SetTitle("Track DCA resolution vs 1/P");
  for (int i=1;i<=trh->GetNbinsX();i++) {
    cx = GetDCASigma(1./trh->GetBinCenter(i), &cxe);
    trh->SetBinContent(i,cx*1e4);
    trh->SetBinError(i,cxe*1e4);
  }
  trh->GetXaxis()->SetTitle("1/p [GeV]");
  trh->GetYaxis()->SetTitle("#sigma_{DCA} [#mum]");
  trh->SetMarkerStyle(20);
  trh->Draw("p");
  gPad->Modified();
  trh->SetTitleOffset(2.5,"Y");
  //
  cnv->cd();
  gPad->Modified();
  //
  if (outname) cnv->Print(outname);
  //
  return cnv;
}


//_____________________________________________________________________
Long64_t AliIntSpotEstimator::Merge(TCollection *coll)
{
  // Merge estimators: used to put together the results of parallel processing
  if(!coll) return 0;
  if (coll->IsEmpty()) return 1;
  TIterator* iter = coll->MakeIterator();
  TObject* obj;
  int count = 0;
  while ((obj = iter->Next())) {
    AliIntSpotEstimator* entry = dynamic_cast<AliIntSpotEstimator*>(obj);
    if (!entry) continue;
    (*this) += *entry;
    count++;
  }
  return count;
  //
}

//_____________________________________________________________________
Double_t AliIntSpotEstimator::GetIPCenter(Int_t id,Double_t *err) const
{
  // calculate IP center in axis id and error
  double cen = fIPCenterStat>0 ? fIPCenter[id]/fIPCenterStat:0;
  if (err) {
    *err = fIPCenterStat>1 ? (fIPCen2[id] - cen*cen*fIPCenterStat)/(fIPCenterStat-1) : 0;
    *err = *err > 0 ? TMath::Sqrt(*err/fIPCenterStat) : 0;
  }
  return cen;
}
