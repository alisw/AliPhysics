///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber                                                  //
//  This class contains the basic functions for the Time Projection Chamber  //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTPCClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBS.h>
#include <TObjectTable.h>
#include "GParticle.h"
#include "AliTPC.h"
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMC.h"

ClassImp(AliTPC) 

//_____________________________________________________________________________
AliTPC::AliTPC()
{
  //
  // Default constructor
  //
  fIshunt   = 0;
  fClusters = 0;
  fHits     = 0;
  fDigits   = 0;
  fTracks   = 0;
  fNsectors = 0;
  fNtracks  = 0;
  fNclusters= 0;
}
 
//_____________________________________________________________________________
AliTPC::AliTPC(const char *name, const char *title)
      : AliDetector(name,title)
{
  //
  // Standard constructor
  //

  //
  // Initialise arrays of hits and digits 
  fHits     = new TClonesArray("AliTPChit",  176);
  fDigits   = new TClonesArray("AliTPCdigit",10000);
  //
  // Initialise counters
  fClusters = 0;
  fTracks   = 0;
  fNsectors = 72;
  fNtracks  = 0;
  fNclusters= 0;
  fDigitsIndex = new Int_t[fNsectors+1];
  fClustersIndex = new Int_t[fNsectors+1];
  //
  fIshunt     =  0;
  //
  // Initialise color attributes
  SetMarkerColor(kYellow);
}

//_____________________________________________________________________________
AliTPC::~AliTPC()
{
  //
  // TPC destructor
  //
  fIshunt   = 0;
  delete fHits;
  delete fDigits;
  delete fClusters;
  delete fTracks;
  if (fDigitsIndex)   delete [] fDigitsIndex;
  if (fClustersIndex) delete [] fClustersIndex;
}

//_____________________________________________________________________________
void AliTPC::AddCluster(Float_t *hits, Int_t *tracks)
{
  //
  // Add a simulated cluster to the list
  //
  if(!fClusters) fClusters=new TClonesArray("AliTPCcluster",10000);
  TClonesArray &lclusters = *fClusters;
  new(lclusters[fNclusters++]) AliTPCcluster(hits,tracks);
}
 
//_____________________________________________________________________________
void AliTPC::AddCluster(const AliTPCcluster &c)
{
  //
  // Add a simulated cluster copy to the list
  //
  if(!fClusters) fClusters=new TClonesArray("AliTPCcluster",10000);
  TClonesArray &lclusters = *fClusters;
  new(lclusters[fNclusters++]) AliTPCcluster(c);
}
 
//_____________________________________________________________________________
void AliTPC::AddDigit(Int_t *tracks, Int_t *digits)
{
  //
  // Add a TPC digit to the list
  //
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliTPCdigit(tracks,digits);
}
 
//_____________________________________________________________________________
void AliTPC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit to the list
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTPChit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliTPC::AddTrack(Float_t *hits)
{
  //
  // Add a track to the list of tracks
  //
  TClonesArray &ltracks = *fTracks;
  new(ltracks[fNtracks++]) AliTPCtrack(hits);
}

//_____________________________________________________________________________
void AliTPC::AddTrack(const AliTPCtrack& t)
{
  //
  // Add a track copy to the list of tracks
  //
  if(!fTracks) fTracks=new TClonesArray("AliTPCtrack",10000);
  TClonesArray &ltracks = *fTracks;
  new(ltracks[fNtracks++]) AliTPCtrack(t);
}

//_____________________________________________________________________________
void AliTPC::BuildGeometry()
{
  //
  // Build TPC ROOT TNode geometry for the event display
  //
  TNode *Node, *Top;
  TTUBS *tubs;
  Int_t i;
  const int kColorTPC=19;
  char name[5], title[18];
  const Double_t kDegrad=TMath::Pi()/180;
  const Double_t loAng=30;
  const Double_t hiAng=15;
  const Int_t nLo = Int_t (360/loAng+0.5);
  const Int_t nHi = Int_t (360/hiAng+0.5);
  const Double_t loCorr = 1/TMath::Cos(0.5*loAng*kDegrad);
  const Double_t hiCorr = 1/TMath::Cos(0.5*hiAng*kDegrad);
  //
  // Get ALICE top node
  Top=gAlice->GetGeometry()->GetNode("alice");
  //
  // Inner sectors
  for(i=0;i<nLo;i++) {
    sprintf(name,"LS%2.2d",i);
    sprintf(title,"TPC low sector %d",i);
    tubs = new TTUBS(name,title,"void",88*loCorr,136*loCorr,250,loAng*(i-0.5),loAng*(i+0.5));
    tubs->SetNumberOfDivisions(1);
    Top->cd();
    Node = new TNode(name,title,name,0,0,0,"");
    Node->SetLineColor(kColorTPC);
    fNodes->Add(Node);
  }
  // Outer sectors
  for(i=0;i<nHi;i++) {
    sprintf(name,"US%2.2d",i);
    sprintf(title,"TPC upper sector %d",i);
    tubs = new TTUBS(name,title,"void",142*hiCorr,250*hiCorr,250,hiAng*(i-0.5),hiAng*(i+0.5));
    tubs->SetNumberOfDivisions(1);
    Top->cd();
    Node = new TNode(name,title,name,0,0,0,"");
    Node->SetLineColor(kColorTPC);
    fNodes->Add(Node);
  }
}
 

//_____________________________________________________________________________
void AliTPC::CreateList(Int_t *tracks,Float_t signal[][MAXTPCTBK+1],
			Int_t ntr,Int_t time)
{
  //
  // Creates list of tracks contributing to a given digit
  // Only the 3 most significant tracks are taken into account
  //
  
  Int_t i,j;
  
  for(i=0;i<3;i++) tracks[i]=-1;
  
  //
  //  Loop over signals, only 3 times
  //
  
  Float_t qmax;
  Int_t jmax;
  Int_t jout[3] = {-1,-1,-1};

  for(i=0;i<3;i++){
    qmax=0.;
    jmax=0;
    
    for(j=0;j<ntr;j++){
      
      if((i == 1 && j == jout[i-1]) 
	 ||(i == 2 && (j == jout[i-1] || j == jout[i-2]))) continue;
      
      if(signal[j][time] > qmax) {
	qmax = signal[j][time];
	jmax=j;
      }       
    } 
    
    if(qmax > 0.) {
      tracks[i]=jmax; 
      jout[i]=jmax;
    }
    
  } 
  
  for(i=0;i<3;i++){
    if(tracks[i] < 0){
      tracks[i]=0;
    }
    else {
      tracks[i]=(Int_t) signal[tracks[i]][0]; // first element is a track number
    }
  }
}

//_____________________________________________________________________________
void AliTPC::DigSignal(Int_t isec,Int_t irow,TObjArray *pointer)
{
  //
  // Digitalise TPC signal
  //
  Int_t pad_c,n_of_pads;
  Int_t pad_number;
  
  n_of_pads = (isec < 25) ? npads_low[irow] : npads_up[irow];
  pad_c=(n_of_pads+1)/2; // this is the "central" pad for a row

  Int_t track,idx;
  Int_t entries;
  TVector *pp;
  TVector *ppp;
  Float_t y,yy,z;
  Float_t pad_signal = 0;
  Float_t signal[MAXTPCTBK]; // Integrated signal over all tracks
  Float_t signal_tr[100][MAXTPCTBK+1]; // contribution from less than 50 tracks
  Int_t flag; // flag indicating a track contributing to a pad signal
  
  //

  Int_t ntracks = pointer->GetEntriesFast();

  if(ntracks == 0) return; // no signal on this row!
  
  //--------------------------------------------------------------
  // For each electron calculate the pad number and the avalanche
  //             This is only once per pad row
  //--------------------------------------------------------------
  
  TVector **el = new TVector* [ntracks]; // each track is a new vector

  TObjArray *arr = new TObjArray; // array of tracks for this row
  
  for(track=0;track<ntracks;track++) {
    pp = (TVector*) pointer->At(track);
    entries = pp->GetNrows();
    el[track] = new TVector(entries-1);
    TVector &v1 = *el[track];
    TVector &v2 = *pp;
    
    for(idx=0;idx<entries-2;idx+=2)
      {
	y=v2(idx+1);
	yy=TMath::Abs(y);
	
	Float_t range=((n_of_pads-1)/2 + 0.5)*pad_pitch_w;
	//
	// Pad number and pad range
	//
	if(yy < 0.5*pad_pitch_w){
	  pad_number=pad_c;
	}
	else if (yy < range){
	  pad_number=(Int_t) ((yy-0.5*pad_pitch_w)/pad_pitch_w +1.);
	  pad_number=(Int_t) (pad_number*TMath::Sign(1.,(double) y)+pad_c);
	}  
	else{
	  pad_number=0;
	}
	
	v1(idx) = (Float_t) pad_number;
	
	// Avalanche, taking the fluctuations into account
	
	Int_t gain_fluct = (Int_t) (-gas_gain*TMath::Log(gRandom->Rndm()));
	v1(idx+1)= (Float_t) gain_fluct;
	
      } // end of loop over electrons
    
    arr->Add(el[track]); // add the vector to the array
    
  }// end of loop over tracks
  
  delete [] el; //  delete an array of pointers
  
  //-------------------------------------------------------------
  //  Calculation of signal for every pad 
  //-------------------------------------------------------------

  //-------------------------------------------------------------
  // Loop over pads
  //-------------------------------------------------------------


  for(Int_t np=1;np<n_of_pads+1;np++)
    {
      for(Int_t l =0;l<MAXTPCTBK;l++) signal[l]=0.; // set signals for this pad to 0
      
      for(Int_t k1=0;k1<100;k1++){
	for(Int_t k2=0;k2<MAXTPCTBK+1;k2++) signal_tr[k1][k2]=0.;
      }
      Int_t track_counter=0; 
      //
      
      //---------------------------------------------------------
      // Loop over tracks
      // --------------------------------------------------------
      
      for(track=0;track<ntracks;track++)
	{
          flag = 0;
          pp = (TVector*) pointer->At(track);
          ppp = (TVector*) arr->At(track);
	  
          TVector &v1 = *pp;
          TVector &v2 = *ppp;
          
          entries = pp->GetNrows();
	  

	  //----------------------------------------------------
	  // Loop over electrons
	  //----------------------------------------------------
	  
          for(idx=0;idx<entries-2;idx+=2)
	    {
	      
	      pad_number = (Int_t) v2(idx);
	      
	      if(pad_number == 0) continue; // skip electrons outside range
	      
	      Int_t pad_dist = pad_number-np;
	      Int_t abs_dist = TMath::Abs(pad_dist);
	      
	      if(abs_dist > 3) continue; // beyond signal range
	      
	      y=  v1(idx+1);
	      z = v1(idx+2);
	      
	      Float_t dist = y-(pad_number-pad_c)*pad_pitch_w;
	      
	      //----------------------------------------------
	      // Calculate the signal induced on a pad "np"
	      //----------------------------------------------
	      
	      if(pad_dist < 0) dist = -dist;
	      
	      switch((Int_t) abs_dist){
	      case 0 : pad_signal = P4(dist); // electron is on pad "np"
                          break;
	      case 1 : pad_signal = P3(dist); // electron is 1 pad away
		break;
	      case 2 : pad_signal = P2(dist); // electron is 2 pads away
		break;
	      case 3 : pad_signal = P1(dist); // electron is 3 pads away
	      }
	      
	      //---------------------------------
              //  Multiply by a gas gain
	      //---------------------------------
	      
              pad_signal=pad_signal*v2(idx+1);
	      
              flag = 1;
	      
	      
	      //-----------------------------------------------
	      //  Sample the pad signal in time
	      //-----------------------------------------------
	      
	      Float_t t_drift = (z_end-TMath::Abs(z))/v_drift; // drift time
	      
	      Float_t t_offset = t_drift-t_sample*(Int_t)(t_drift/t_sample);   
	      Int_t first_bucket = (Int_t) (t_drift/t_sample+1.); 
	      
	      for(Int_t bucket = 1;bucket<6;bucket++){
		Float_t time = (bucket-1)*t_sample+t_offset;
		Int_t time_idx = first_bucket+bucket-1;
		Float_t ampl = pad_signal*TimeRes(time);
		if (time_idx > MAXTPCTBK) break; //Y.Belikov
		if (track_counter >=100) break;  //Y.Belikov

		signal_tr[track_counter][time_idx] += ampl; // single track only
		signal[time_idx-1] += ampl;  // fill a signal array for this pad
		
	      } // end of time sampling
	      
	    } // end of loop over electrons for a current track
	  
	  //-----------------------------------------------
	  //  add the track number and update the counter
	  //  if it gave a nonzero contribution to the pad
	  //-----------------------------------------------
	  if(flag != 0 && track_counter < 100){
          signal_tr[track_counter][0] = v1(0);
          track_counter++;      // counter is looking at the NEXT track!
	  }
	  
	} // end of loop over tracks
      
      //----------------------------------------------
      //  Fill the Digits for this pad
      //----------------------------------------------
      
      Int_t tracks[3];
      Int_t digits[5];
      
      digits[0]=isec; // sector number
      digits[1]=irow+1; // row number
      digits[2]=np;   // pad number
      
      Float_t q;
      
      for(Int_t time = 0;time<MAXTPCTBK;time++){
	digits[3] = time+1; // time bucket
	
	q = signal[time];
	q = gRandom->Gaus(q,sigma_noise); // apply noise
	
	q *= (q_el*1.e15); // convert to fC
	q *= chip_gain; // convert to mV   
	q *= (adc_sat/dyn_range); // convert to ADC counts
	
	if(q < zero_sup) continue; // do not fill "zeros"
	if(q > adc_sat) q = adc_sat;  // saturation
	digits[4] = (Int_t) q;                   // ADC counts
	
	//--------------------------------------------------
	//    "Real signal" or electronics noise
	//--------------------------------------------------
	
	if(signal[time] > 0.){
	  
	  //-----------------------------------------------------
	  //  Create a list of tracks contributing to this digit
	  //  If the digit results from a noise, track number is 0
	  //-----------------------------------------------------
	  
	  CreateList(tracks,signal_tr,track_counter,time);
	}
	else {
	  for(Int_t ii=0;ii<3;ii++) tracks[ii]=0;
	}
	
	AddDigit(tracks,digits);
	
	
      } // end of digits for this pad
      
    } // end of loop over pads
  
  arr->Delete(); // delete objects in this array
  
  delete arr; // deletes the TObjArray itselves
  
}

//_____________________________________________________________________________
Int_t AliTPC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Calculate distance from TPC to mouse on the display
  // Dummy procedure
  //
  return 9999;
}

//_____________________________________________________________________________
const int MAX_CLUSTER=nrow_low+nrow_up; 
const int S_MAXSEC=24;
const int L_MAXSEC=48;
const int ROWS_TO_SKIP=21;
const Float_t MAX_CHI2=15.;
const Float_t THRESHOLD=8*zero_sup;

//_____________________________________________________________________________
static Double_t SigmaY2(Double_t r, Double_t tgl, Double_t pt)
{
  //
  // Calculate spread in Y
  //
  pt=TMath::Abs(pt)*1000.;
  Double_t x=r/pt;
  tgl=TMath::Abs(tgl);
  Double_t s=a_rphi - b_rphi*r*tgl + c_rphi*x*x + d_rphi*x;
  if (s<0.4e-3) s=0.4e-3;
  return s;
}

//_____________________________________________________________________________
static Double_t SigmaZ2(Double_t r, Double_t tgl) 
{
  //
  // Calculate spread in Z
  //
  tgl=TMath::Abs(tgl);
  Double_t s=a_z - b_z*r*tgl + c_z*tgl*tgl;
  if (s<0.4e-3) s=0.4e-3;
  return s;
}

//_____________________________________________________________________________
inline Double_t f1(Double_t x1,Double_t y1,   //C
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //
  // Function f1
  //
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));
  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -xr*yr/sqrt(xr*xr+yr*yr); 
}


//_____________________________________________________________________________
inline Double_t f2(Double_t x1,Double_t y1,  //eta=C*x0
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //
  // Function f2
  //
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));
  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
inline Double_t f3(Double_t x1,Double_t y1,  //tgl
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //
  // Function f3
  //
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//_____________________________________________________________________________
static int FindProlongation(AliTPCtrack& t, const AliTPCSector *sec, 
			    int s, int ri, int rf=0) 
{
  //
  // Propagate track
  //
  int try_again=ROWS_TO_SKIP;
  Double_t alpha=sec->GetAlpha();
  int ns=int(2*TMath::Pi()/alpha)+1;
  for (int nr=ri; nr>=rf; nr--) {
    Double_t x=sec[s].GetX(nr), ymax=sec[s].GetMaxY(nr);
    if (!t.PropagateTo(x)) break;
    
    const AliTPCcluster *cl=0;
    Double_t max_chi2=MAX_CHI2;
    const AliTPCRow& row=sec[s][nr];
    Double_t sy2=SigmaY2(t.GetX(),t.GetTgl(),t.GetPt());
    Double_t sz2=SigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=3.*sqrt(t.GetSigmaY2() + sy2), y=t.GetY(), z=t.GetZ();
    
    if (road>30) {
      if (t>3) cerr<<t<<" AliTPCtrack warning: Too broad road !\n"; 
      break;
    }
    
    if (row) 
      for (int i=row.Find(y-road); i<row; i++) {
	AliTPCcluster* c=(AliTPCcluster*)(row[i]);
	if (c->fY > y+road) break;
	if (c->IsUsed()) continue;
	if ((c->fZ - z)*(c->fZ - z) > 9.*(t.GetSigmaZ2() + sz2)) continue;
	Double_t chi2=t.GetPredictedChi2(c);
	if (chi2 > max_chi2) continue;
	max_chi2=chi2;
	cl=c;       
      }
    if (cl) {
      t.Update(cl,max_chi2);
      try_again=ROWS_TO_SKIP;
    } else {
      if (try_again--) {
	if (y > ymax) {
	  s = (s+1) % ns;
	  if (!t.Rotate(alpha)) break;
	} else
	  if (y <-ymax) {
	    s = (s-1+ns) % ns;
	    if (!t.Rotate(-alpha)) break;
	  };
	continue;
      }
      break;
    }
  }
  return s;
}

//_____________________________________________________________________________
static void MakeSeeds(TObjArray& seeds,const AliTPCSector* sec,int i1,int i2)
{
  //
  // Find seed for tracking
  //
  TMatrix C(5,5); TVector x(5);
  int max_sec=L_MAXSEC/2;
  for (int ns=0; ns<max_sec; ns++) {
    int nl=sec[(ns-1+max_sec)%max_sec][i2];
    int nm=sec[ns][i2];
    int nu=sec[(ns+1)%max_sec][i2];
    Double_t alpha=sec[ns].GetAlpha();
    const AliTPCRow& r1=sec[ns][i1];
    for (int is=0; is < r1; is++) {
      Double_t x1=sec[ns].GetX(i1), y1=r1[is]->fY, z1=r1[is]->fZ;
      for (int js=0; js < nl+nm+nu; js++) {
	const AliTPCcluster *cl;
	Double_t cs,sn;
	//int ks;
	
	if (js<nl) {
	  //ks=(ns-1+max_sec)%max_sec;
	  const AliTPCRow& r2=sec[(ns-1+max_sec)%max_sec][i2];
	  cl=r2[js];
	  cs=cos(alpha); sn=sin(alpha);
	} else 
	  if (js<nl+nm) {
	    //ks=ns;
	    const AliTPCRow& r2=sec[ns][i2];
	    cl=r2[js-nl];
	    cs=1; sn=0.;
	  } else {
	    //ks=(ns+1)%max_sec;
	    const AliTPCRow& r2=sec[(ns+1)%max_sec][i2];
	    cl=r2[js-nl-nm];
	    cs=cos(alpha); sn=-sin(alpha);
	  }
	Double_t x2=sec[ns].GetX(i2), y2=cl->fY, z2=cl->fZ;
	//if (z1*z2 < 0) continue;
	//if (TMath::Abs(z1) < TMath::Abs(z2)) continue;
	
	Double_t tmp= x2*cs+y2*sn;
	y2 =-x2*sn+y2*cs;
	x2=tmp;     
	
	x(0)=y1;
	x(1)=z1;
	x(2)=f1(x1,y1,x2,y2,0.,0.);
	x(3)=f2(x1,y1,x2,y2,0.,0.);
	x(4)=f3(x1,y1,x2,y2,z1,z2);
	
	if (TMath::Abs(x(2)*x1-x(3)) >= 0.999) continue;
	
	if (TMath::Abs(x(4)) > 1.2) continue;
	
	Double_t a=asin(x(3));
	/*
	  Double_t tgl1=z1*x(2)/(a+asin(x(2)*x1-x(3)));
	  Double_t tgl2=z2*x(2)/(a+asin(x(2)*x2-x(3)));
	  Double_t ratio=2*(tgl1-tgl2)/(tgl1+tgl2);
	  if (TMath::Abs(ratio)>0.0170) continue; //or > 0.005
	*/
	Double_t zv=z1 - x(4)/x(2)*(a+asin(x(2)*x1-x(3)));
	if (TMath::Abs(zv)>33.) continue; 
	
	
	TMatrix X(6,6); X=0.; 
	X(0,0)=r1[is]->fSigmaY2; X(1,1)=r1[is]->fSigmaZ2;
	X(2,2)=cl->fSigmaY2;     X(3,3)=cl->fSigmaZ2;
	X(4,4)=3./12.; X(5,5)=3./12.;
	TMatrix F(5,6); F.UnitMatrix();
	Double_t sy=sqrt(X(0,0)), sz=sqrt(X(1,1));
	F(2,0)=(f1(x1,y1+sy,x2,y2,0.,0.)-x(2))/sy;
	F(2,2)=(f1(x1,y1,x2,y2+sy,0.,0.)-x(2))/sy;
	F(2,4)=(f1(x1,y1,x2,y2,0.,0.+sy)-x(2))/sy;
	F(3,0)=(f2(x1,y1+sy,x2,y2,0.,0.)-x(3))/sy;
	F(3,2)=(f2(x1,y1,x2,y2+sy,0.,0.)-x(3))/sy;
	F(3,4)=(f2(x1,y1,x2,y2,0.,0.+sy)-x(3))/sy;
	F(4,0)=(f3(x1,y1+sy,x2,y2,z1,z2)-x(4))/sy;
	F(4,1)=(f3(x1,y1,x2,y2,z1+sz,z2)-x(4))/sz;
	F(4,2)=(f3(x1,y1,x2,y2+sy,z1,z2)-x(4))/sy;
	F(4,3)=(f3(x1,y1,x2,y2,z1,z2+sz)-x(4))/sz;
	F(4,4)=0;
	F(3,3)=0;
	
	TMatrix t(F,TMatrix::kMult,X);
	C.Mult(t,TMatrix(TMatrix::kTransposed,F));
	
	TrackSeed *track=new TrackSeed(*(r1[is]),x,C);
	FindProlongation(*track,sec,ns,i1-1,i2);
	int ii=(i1-i2)/2;
	if (*track >= ii) {seeds.AddLast(track); continue;}
	else delete track;
      }
    }
  }
}

//_____________________________________________________________________________
void AliTPC::Clusters2Tracks()
{
  //
  // TPC Track finder from clusters.
  //
  if (!fClusters) return;
  AliTPCSSector ssec[S_MAXSEC/2];
  AliTPCLSector lsec[L_MAXSEC/2];
  int ncl=fClusters->GetEntriesFast();
  while (ncl--) {
    AliTPCcluster *c=(AliTPCcluster*)fClusters->UncheckedAt(ncl);
    int sec=int(c->fSector), row=int(c->fPadRow);
    
    if (sec<24) {
      if (row<0 || row>nrow_low) {cerr<<"low !!!"<<row<<endl; continue;}
      ssec[sec%12][row].InsertCluster(c);
    } else {
      if (row<0 || row>nrow_up ) {cerr<<"up  !!!"<<row<<endl; continue;}
      sec -= 24;
      lsec[sec%24][row].InsertCluster(c);
    }
  }
  
  
  TObjArray seeds(20000);
  MakeSeeds(seeds,lsec,nrow_up-1,nrow_up-1-8);
  MakeSeeds(seeds,lsec,nrow_up-1-4,nrow_up-1-4-8);
  
  seeds.Sort();
  
  int found=0;
  int nseed=seeds.GetEntriesFast();
  
  for (int s=0; s<nseed; s++) {
    AliTPCtrack& t=*((AliTPCtrack*)seeds.UncheckedAt(s));
    Double_t alpha=t.GetAlpha();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
    int ns=int(alpha/lsec->GetAlpha() + 0.5);
    
    Double_t x=t.GetX();
    int nr;
    if (x<pad_row_up[nrow_up-1-4-7]) nr=nrow_up-1-4-8;
    else if (x<pad_row_up[nrow_up-1-7]) nr=nrow_up-1-8;
    else {cerr<<x<<" =x !!!\n"; continue;}
    
    int ls=FindProlongation(t,lsec,ns,nr-1);
    
    //   if (t < 25) continue;
    
    x=t.GetX(); alpha=lsec[ls].GetAlpha();          //
    Double_t phi=ls*alpha + atan(t.GetY()/x);       // Find S-sector
    int ss=int(0.5*(phi/alpha+1));                  //
    alpha *= 2*(ss-0.5*ls);                         // and rotation angle
    if (!t.Rotate(alpha)) continue;                 //
    ss %= (S_MAXSEC/2);                             //
    
    ss=FindProlongation(t,ssec,ss,nrow_low-1);
    if (t < 30) continue;
    
    AddTrack(t);
    t.UseClusters();
    cerr<<found++<<'\r';
  }  
}

//_____________________________________________________________________________
void AliTPC::CreateMaterials()
{
  //
  // Create Materials for for TPC
  // Origin M.Kowalski
  //
  
  AliMC* pMC = AliMC::GetMC();

  Int_t ISXFLD=gAlice->Field()->Integ();
  Float_t SXMGMX=gAlice->Field()->Max();
  
  Float_t absl, radl, a, d, z;
  Float_t dg;
  Float_t x0ne;
  Float_t buf[1];
  Int_t nbuf;
  
  // --- Methane (CH4) --- 
  Float_t am[2] = { 12.,1. };
  Float_t zm[2] = { 6.,1. };
  Float_t wm[2] = { 1.,4. };
  Float_t dm    = 7.17e-4;
  // --- The Neon CO2 90/10 mixture --- 
  Float_t ag[2] = { 20.18 };
  Float_t zg[2] = { 10. };
  Float_t wg[2] = { .8,.2 };
  Float_t dne   = 9e-4;        // --- Neon density in g/cm3 ---     
  // --- Mylar (C5H4O2) --- 
  Float_t amy[3] = { 12.,1.,16. };
  Float_t zmy[3] = { 6.,1.,8. };
  Float_t wmy[3] = { 5.,4.,2. };
  Float_t dmy    = 1.39;
  // --- CO2 --- 
  Float_t ac[2] = { 12.,16. };
  Float_t zc[2] = { 6.,8. };
  Float_t wc[2] = { 1.,2. };
  Float_t dc    = .001977;
  // --- Carbon density and radiation length --- 
  Float_t densc = 2.265;
  Float_t radlc = 18.8;
  // --- Silicon --- 
  Float_t asi   = 28.09;
  Float_t zsi   = 14.;
  Float_t desi  = 2.33;
  Float_t radsi = 9.36;
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(0, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  x0ne = 28.94 / dne;
  AliMaterial(1, "Ne $", 20.18, 10., dne, x0ne, 999.);
  
  // --  Methane, defined by the proportions of atoms 
  
  AliMixture(2, "Methane$", am, zm, dm, -2, wm);
  
  // --- CO2, defined by the proportion of atoms 
  
  AliMixture(7, "CO2$", ac, zc, dc, -2, wc);
  
  // --  Get A,Z etc. for CO2 
  
  char namate[21];
  pMC->Gfmate((*fIdmate)[7], namate, a, z, d, radl, absl, buf, nbuf);

  ag[1] = a;
  zg[1] = z;
  dg = dne * .9 + dc * .1;

  //-------------------------------------
  
  // --  Create Ne/CO2 90/10 mixture 
  
  AliMixture(3, "Gas-mixt $", ag, zg, dg, 2, wg);
  AliMixture(4, "Gas-mixt $", ag, zg, dg, 2, wg);
  
  AliMaterial(5, "G10$", 20., 10., 1.7, 19.4, 999.);
  AliMixture(6, "Mylar$", amy, zmy, dmy, -3, wmy);
  
  a = ac[0];
  z = zc[0];
  AliMaterial(8, "Carbon", a, z, densc, radlc, 999.);
  
  AliMaterial(9, "Silicon", asi, zsi, desi, radsi, 999.);
  AliMaterial(99, "Air$", 14.61, 7.3, .001205, 30420., 67500.);
  
  AliMedium(400, "Al wall$",  0, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1);
  AliMedium(402, "Gas mix1$", 3, 0, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(403, "Gas mix2$", 3, 0, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(404, "Gas mix3$", 4, 1, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(405, "G10 pln$",  5, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
  AliMedium(406, "Mylar  $",  6, 0, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(407, "CO2    $",  7, 0, ISXFLD, SXMGMX, 10., .01,.1, .01,  .01);
  AliMedium(408, "Carbon $",  8, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
  AliMedium(409, "Silicon$",  9, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
  AliMedium(499, "Air gap$", 99, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
}

//_____________________________________________________________________________
struct Bin {
   const AliTPCdigit *dig;
   int idx;
   Bin() {dig=0; idx=-1;}
};

struct PreCluster : public AliTPCcluster {
   const AliTPCdigit* summit;
   int idx;
   int cut;
   int npeaks;
   PreCluster() : AliTPCcluster() {cut=npeaks=0;}
};

//_____________________________________________________________________________
static void FindCluster(int i, int j, Bin bins[][MAXTPCTBK+2], PreCluster &c) 
{
  //
  // Find clusters
  //
  Bin& b=bins[i][j];
  int q=b.dig->fSignal;
  
  if (q<0) { // digit is at the edge of the pad row
    q=-q;
    c.cut=1;
  } 
  if (b.idx >= 0 && b.idx != c.idx) {
    c.idx=b.idx;
    c.npeaks++;
  }
  
  if (q > TMath::Abs(c.summit->fSignal)) c.summit=b.dig;
  
  c.fY += i*q;
  c.fZ += j*q;
  c.fSigmaY2 += i*i*q;
  c.fSigmaZ2 += j*j*q;
  c.fQ += q;
  
  b.dig = 0;  b.idx = c.idx;
  
  if (bins[i-1][j].dig) FindCluster(i-1,j,bins,c);
  if (bins[i][j-1].dig) FindCluster(i,j-1,bins,c);
  if (bins[i+1][j].dig) FindCluster(i+1,j,bins,c);
  if (bins[i][j+1].dig) FindCluster(i,j+1,bins,c);
  
}

//_____________________________________________________________________________
void AliTPC::Digits2Clusters()
{
  //
  // simple TPC cluster finder from digits.
  //
  const Int_t MAX_PAD=200+2, MAX_BUCKET=MAXTPCTBK+2;
  const Int_t Q_min=200;//75;
  
  TTree *t=gAlice->TreeD();
  t->GetBranch("TPC")->SetAddress(&fDigits);
  Int_t sectors_by_rows=(Int_t)t->GetEntries();
  
  int ncls=0;
  
  for (Int_t n=0; n<sectors_by_rows; n++) {
    if (!t->GetEvent(n)) continue;
    Bin bins[MAX_PAD][MAX_BUCKET];
    AliTPCdigit *dig=(AliTPCdigit*)fDigits->UncheckedAt(0);
    Int_t nsec=dig->fSector, nrow=dig->fPadRow;
    Int_t ndigits=fDigits->GetEntriesFast();
    
    int npads;  int sign_z;
    if (nsec<25) {
      sign_z=(nsec<13) ? 1 : -1;
      npads=npads_low[nrow-1];
    } else {
      sign_z=(nsec<49) ? 1 : -1;
      npads=npads_up[nrow-1];
    }
    
    int ndig;
    for (ndig=0; ndig<ndigits; ndig++) {
      dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
      int i=dig->fPad, j=dig->fTime;
      if (dig->fSignal >= THRESHOLD) bins[i][j].dig=dig;
      if (i==1 || i==npads || j==1 || j==MAXTPCTBK) dig->fSignal*=-1;
    }
    
    int ncl=0;
    int i,j;
    for (i=1; i<MAX_PAD-1; i++) {
      for (j=1; j<MAX_BUCKET-1; j++) {
	if (bins[i][j].dig == 0) continue;
	PreCluster c; c.summit=bins[i][j].dig; c.idx=ncls;
	FindCluster(i, j, bins, c);
	//if (c.fQ <= Q_min) continue; //noise cluster
	c.fY /= c.fQ;
	c.fZ /= c.fQ;
	c.fSigmaY2 = c.fSigmaY2/c.fQ - c.fY*c.fY + 1./12.;
	c.fSigmaZ2 = c.fSigmaZ2/c.fQ - c.fZ*c.fZ + 1./12.;
	c.fSigmaY2 *= pad_pitch_w*pad_pitch_w;
	c.fSigmaZ2 *= z_end/MAXTPCTBK*z_end/MAXTPCTBK;
	c.fSigmaY2 *= 0.022*8;
	c.fSigmaZ2 *= 0.068*4;
	c.fY = (c.fY - 0.5 - 0.5*npads)*pad_pitch_w;
	c.fZ = z_end/MAXTPCTBK*c.fZ; 
	c.fZ -= 3.*fwhm/2.35482*v_drift; // PASA delay 
	c.fZ = sign_z*(z_end - c.fZ);
	c.fSector=nsec-1;
	c.fPadRow=nrow-1;
	c.fTracks[0]=c.summit->fTracks[0];
	c.fTracks[1]=c.summit->fTracks[1];
	c.fTracks[2]=c.summit->fTracks[2];
	
	if (c.cut) {
	  c.fSigmaY2 *= 25.;
	  c.fSigmaZ2 *= 4.;
	}
	
	AddCluster(c); ncls++; ncl++;
      }
    }
    
    for (ndig=0; ndig<ndigits; ndig++) {
      dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
      if (TMath::Abs(dig->fSignal) >= THRESHOLD/3) 
	bins[dig->fPad][dig->fTime].dig=dig;
    }
    
    for (i=1; i<MAX_PAD-1; i++) {
      for (j=1; j<MAX_BUCKET-1; j++) {
	if (bins[i][j].dig == 0) continue;
	PreCluster c; c.summit=bins[i][j].dig; c.idx=ncls;
	FindCluster(i, j, bins, c);
	if (c.fQ <= Q_min) continue; //noise cluster
	if (c.npeaks>1) continue;    //overlapped cluster
	c.fY /= c.fQ;
	c.fZ /= c.fQ;
	c.fSigmaY2 = c.fSigmaY2/c.fQ - c.fY*c.fY + 1./12.;
	c.fSigmaZ2 = c.fSigmaZ2/c.fQ - c.fZ*c.fZ + 1./12.;
	c.fSigmaY2 *= pad_pitch_w*pad_pitch_w;
	c.fSigmaZ2 *= z_end/MAXTPCTBK*z_end/MAXTPCTBK;
	c.fSigmaY2 *= 0.022*4;
	c.fSigmaZ2 *= 0.068*4;
	c.fY = (c.fY - 0.5 - 0.5*npads)*pad_pitch_w;
	c.fZ = z_end/MAXTPCTBK*c.fZ; 
	c.fZ -= 3.*fwhm/2.35482*v_drift; // PASA delay 
	c.fZ = sign_z*(z_end - c.fZ);
	c.fSector=nsec-1;
	c.fPadRow=nrow-1;
	c.fTracks[0]=c.summit->fTracks[0];
	c.fTracks[1]=c.summit->fTracks[1];
	c.fTracks[2]=c.summit->fTracks[2];
	
	if (c.cut) {
	  c.fSigmaY2 *= 25.;
	  c.fSigmaZ2 *= 4.;
	}
	
	if (c.npeaks==0) {AddCluster(c); ncls++; ncl++;}
	else {
	  new ((*fClusters)[c.idx]) AliTPCcluster(c);
	}
      }
    }
    
    cerr<<"sector, row, digits, clusters: "
	<<nsec<<' '<<nrow<<' '<<ndigits<<' '<<ncl<<"                  \r";
    
    fDigits->Clear();
    
  }
}

//_____________________________________________________________________________
void AliTPC::ElDiff(Float_t *xyz)
{
  //
  // calculates the diffusion of a single electron
  //
  
  Float_t driftl;
  //
  Float_t z0=xyz[2];
  driftl=z_end-TMath::Abs(xyz[2]);
  if(driftl<0.01) driftl=0.01;
  driftl=TMath::Sqrt(driftl);
  Float_t sig_t = driftl*diff_t;
  Float_t sig_l = driftl*diff_l;
  //
  xyz[0]=gRandom->Gaus(xyz[0],sig_t);
  xyz[1]=gRandom->Gaus(xyz[1],sig_t);
  xyz[2]=gRandom->Gaus(xyz[2],sig_l);
  //
  if (TMath::Abs(xyz[2])>z_end){
    xyz[2]=z_end*TMath::Sign(1.,(double) z0);
  }
  if(xyz[2]*z0 < 0.){
    xyz[2]=0.0001*TMath::Sign(1.,(double) z0);
  } 
}

//_____________________________________________________________________________
void AliTPC::Hits2Clusters()
{
  //
  // TPC simple cluster generator from hits
  // obtained from the TPC Fast Simulator
  //
  Float_t sigma_rphi,sigma_z,cl_rphi,cl_z;
  //
  GParticle *particle; // pointer to a given particle
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  TClonesArray *Particles; //pointer to the particle list
  Int_t sector,nhits;
  Int_t ipart;
  Float_t xyz[5];
  Float_t pl,pt,tanth,rpad,ratio;
  Float_t rot_angle;
  Float_t cph,sph;
  
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *TH = gAlice->TreeH();
  Stat_t ntracks = TH->GetEntries();
  
  //------------------------------------------------------------
  // Loop over all sectors (72 sectors)
  // Sectors 1-24 are lower sectors, 1-12 z>0, 13-24 z<0
  // Sectors 25-72 are upper sectors, 25-48 z>0, 49-72 z<0
  //
  // First cluster for sector 1 starts at "0"
  //------------------------------------------------------------
  
  
  fClustersIndex[0] = 0;
  
  //
  for(Int_t isec=1;isec<fNsectors+1;isec++){
    //
    if(isec < 25){
      rot_angle = (isec < 13) ? (isec-1)*alpha_low : (isec-13)*alpha_low;
    }
    else {
      rot_angle = (isec < 49) ? (isec-25)*alpha_up : (isec-49)*alpha_up;
    }
    
    cph=TMath::Cos(rot_angle);
    sph=TMath::Sin(rot_angle);
    
    //------------------------------------------------------------
    // Loop over tracks
    //------------------------------------------------------------
    
    for(Int_t track=0;track<ntracks;track++){
      ResetHits();
      TH->GetEvent(track);
      //
      //  Get number of the TPC hits and a pointer
      //  to the particles
      //
      nhits=fHits->GetEntriesFast();
      Particles=gAlice->Particles();
      //
      // Loop over hits
      //
      for(Int_t hit=0;hit<nhits;hit++){
	tpcHit=(AliTPChit*)fHits->UncheckedAt(hit);
	sector=tpcHit->fSector; // sector number
	if(sector != isec) continue; //terminate iteration
	ipart=tpcHit->fTrack;
	particle=(GParticle*)Particles->UncheckedAt(ipart);
	pl=particle->GetPz();
	pt=particle->GetPT();
	if(pt < 1.e-9) pt=1.e-9;
	tanth=pl/pt;
	tanth = TMath::Abs(tanth);
	rpad=TMath::Sqrt(tpcHit->fX*tpcHit->fX + tpcHit->fY*tpcHit->fY);
	ratio=0.001*rpad/pt; // pt must be in MeV/c - historical reason
	
	//   space-point resolutions
	
	sigma_rphi=SigmaY2(rpad,tanth,pt);
	sigma_z   =SigmaZ2(rpad,tanth   );
	
	//   cluster widths
	
	cl_rphi=ac_rphi-bc_rphi*rpad*tanth+cc_rphi*ratio*ratio;
	cl_z=ac_z-bc_z*rpad*tanth+cc_z*tanth*tanth;
	
	// temporary protection
	
	if(sigma_rphi < 0.) sigma_rphi=0.4e-3;
	if(sigma_z < 0.) sigma_z=0.4e-3;
	if(cl_rphi < 0.) cl_rphi=2.5e-3;
	if(cl_z < 0.) cl_z=2.5e-5;
	
	//
	
	//
	// smearing --> rotate to the 1 (13) or to the 25 (49) sector,
	// then the inaccuracy in a X-Y plane is only along Y (pad row)!
	//
	Float_t xprim= tpcHit->fX*cph + tpcHit->fY*sph;
	Float_t yprim=-tpcHit->fX*sph + tpcHit->fY*cph;
	xyz[0]=gRandom->Gaus(yprim,TMath::Sqrt(sigma_rphi));   // y
	xyz[1]=gRandom->Gaus(tpcHit->fZ,TMath::Sqrt(sigma_z)); // z 
	xyz[2]=tpcHit->fQ;                                     // q
	xyz[3]=sigma_rphi;                                     // fSigmaY2
	xyz[4]=sigma_z;                                        // fSigmaZ2
	
	//find row number
	int row;
	if (xprim > 0.5*(pad_row_up[0]+pad_row_low[nrow_low-1])) {
          for (row=0; row<nrow_up; row++) if (xprim < pad_row_up[row]) break;
	} else {
          for (row=0; row<nrow_low; row++) if (xprim < pad_row_low[row]) break;
	}
	
	// and finally add the cluster
	Int_t tracks[5]={tpcHit->fTrack+1, 0, 0, sector-1, row-1};
	AddCluster(xyz,tracks);
	
      } // end of loop over hits
    }   // end of loop over tracks 
    
    fClustersIndex[isec] = fNclusters; // update clusters index
    
  } // end of loop over sectors
  
  fClustersIndex[fNsectors]--; // set end of the clusters buffer
  
}

//_____________________________________________________________________________
void AliTPC::Hits2Digits()
{
  //
  // TPC conversion from hits to digits.
  //
  
  Int_t i;
  //
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  //
  Float_t xyz[3];
  Float_t rot_angle;
  //-------------------------------------------------------
  //  Get the access to the tracks 
  //-------------------------------------------------------
  TTree *TH = gAlice->TreeH();
  Stat_t ntracks = TH->GetEntries();
  
  //----------------------------------------------------
  // Loop over all sectors (72 sectors)
  // Sectors 1-24 are lower sectors, 1-12 z>0, 13-24 z<0
  // Sectors 25-72 are upper sectors, 25-48 z>0, 49-72 z<0
  //----------------------------------------------------
  
  for(Int_t isec=1;isec<fNsectors+1;isec++){ 
    
    //
    printf("*** Processing sector number %d ***\n",isec);
    
    if(isec < 25){
      rot_angle = (isec < 13) ? (isec-1)*alpha_low : (isec-13)*alpha_low;
    }
    else {
      rot_angle = (isec < 49) ? (isec-25)*alpha_up : (isec-49)*alpha_up;
    }
    
    Int_t nrows = (isec<25) ? 23 : 52;
    
    
    Float_t cph=TMath::Cos(rot_angle);
    Float_t sph=TMath::Sin(rot_angle);
    
    
    
    //----------------------------------------------
    // Create TObjArray-s, one for each row
    //---------------------------------------------- 
    
    TObjArray **row = new TObjArray* [nrows];
    for(i=0; i<nrows; i++){
      row[i] = new TObjArray;
    }
    
    //----------------------------------------------
    // Loop over tracks
    //----------------------------------------------
    for(Int_t track=0;track<ntracks;track++){  
      ResetHits();
      TH->GetEvent(track); 
      
      //------------------------------------------------
      //  Get number of the TPC hits and a pointer
      //  to the particles
      //------------------------------------------------
      Int_t nhits=fHits->GetEntriesFast();
      if(nhits == 0) continue; 
      //-----------------------------------------------
      //  Create vectors for storing the track information,
      //  one vector per track per row,
      //  first element is a track number and then
      //  there are (y,z) pairs * number of electrons
      //----------------------------------------------
      
      TVector **tr = new TVector* [nrows];
      Int_t *n_of_electrons= new int [nrows]; // electron counter 
      for(i=0;i<nrows;i++){
	tr[i] = new TVector(241); // 120 electrons initialy
	n_of_electrons[i]=0;
      } 
      //-----------------------------------------------------
      // Loop over hits
      //------------------------------------------------------
      for(Int_t hit=0;hit<nhits;hit++){
	tpcHit=(AliTPChit*)fHits->UncheckedAt(hit);
	Int_t sector=tpcHit->fSector; // sector number
	if(sector != isec) continue; //terminate iteration
	
	xyz[0]=tpcHit->fX;
	xyz[1]=tpcHit->fY;
	xyz[2]=tpcHit->fZ;
	Int_t QI = (Int_t) (tpcHit->fQ); // energy loss (number of electrons)
	
	//-----------------------------------------------
	//  Rotate the electron cluster to sector 1,13,25,49
	//-----------------------------------------------
	Float_t xprim=xyz[0]*cph+xyz[1]*sph;
	Float_t yprim=-xyz[0]*sph+xyz[1]*cph;
	Float_t zprim=xyz[2]; 
	
	//-------------------------------------
	// Loop over electrons
	//-------------------------------------
	for(Int_t nel=0;nel<QI;nel++){
	  xyz[0]=xprim;  //
	  xyz[1]=yprim;  // Keep the initial cluster position!
	  xyz[2]=zprim;  //
	  
	  ElDiff(xyz); // Appply the diffusion
	  
	  Float_t row_first; 
	  Int_t row_number;
	  row_first = (isec<25) ? pad_row_low[0] : pad_row_up[0];
	  
	  row_number=(Int_t) ((xyz[0]-row_first+0.5*pad_pitch_l)/pad_pitch_l);
	  
	  // Check if out of range
	  
          if((xyz[0]-row_first+0.5*pad_pitch_l) < 0 
	     || row_number > (nrows-1)) continue;
          
	  n_of_electrons[row_number]++;
	  
	  //
	  // Expand vector if necessary
	  //
	  if(n_of_electrons[row_number]>120){
	    Int_t range = tr[row_number]->GetNrows();
	    if(n_of_electrons[row_number] > (range-1)/2){
	      tr[row_number]->ResizeTo(range+30); // Add 15 electrons
	    }
	  }
	  
	  //---------------------------------
	  //  E x B effect at the wires
	  //---------------------------------
	  Int_t nw;
	  nw=(nwires+1)/2;
	  Float_t xx,dx;
	  for (Int_t nwire=1;nwire<=nwires;nwire++){
	    xx=(nwire-nw)*ww_pitch+
	      ((isec<13) ? pad_row_low[row_number]:pad_row_up[row_number]);
	    dx=xx-xyz[0];
	    if(TMath::Abs(dx) < 0.5*ww_pitch) {
	      xyz[1]=dx*omega_tau+xyz[1];
	      break;
	    }
	  } // end of loop over the wires
	  
	  TVector &v = *tr[row_number];
	  Int_t idx = 2*n_of_electrons[row_number]-1;
	  v(idx)=xyz[1];
	  v(idx+1)=xyz[2];
	} // end of loop over electrons         
      } // end of loop over hits  
      
      //
      //  The track is finished 
      //     
      int trk=((AliTPChit*)fHits->UncheckedAt(0))->fTrack; //Y.Belikov
      for(i=0;i<nrows;i++){
	TVector &v = *tr[i];
	if(n_of_electrons[i] >0) {
	  //        v(0)=(Float_t)(track+1); // track number
          v(0)=(Float_t)(trk+1); // Y.Belikov 
	  tr[i]->ResizeTo(2*n_of_electrons[i]+1); // shrink if necessary
	  row[i]->Add(tr[i]); // add to the row-array
	}
	else{
	  delete tr[i]; // delete TVector if empty
	}
      }   
      delete [] tr; // delete track pointers  
      delete [] n_of_electrons; // delete n_of_electrons array     
    } //end of loop over tracks 
    //---------------------------------------------------
    //  Digitize the sector data, row by row
    //---------------------------------------------------  
    
    printf("*** Digitizing sector number %d ***\n",isec);
    
    for(i=0;i<nrows;i++){
      
      DigSignal(isec,i,row[i]); 
      gAlice->TreeD()->Fill();      
      int ndig=fDigits->GetEntriesFast();
      printf("*** Sector, row, digits %d %d %d ***\n",isec,i+1,ndig);
      
      ResetDigits();
      
      row[i]->Delete(); // delete objects in this array
      delete row[i]; // delete the TObjArray itselves
      
    } // end of digitization
    
    delete [] row; // delete vectors of pointers to sectors
    
  } //end of loop over sectors 
  
}

//_____________________________________________________________________________
void AliTPC::Init()
{
  //
  // Initialise TPC detector after definition of geometry
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TPC_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliTPC::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the TPC.
  //
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *D = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  }	

  char *R = strstr(option,"R");

  if (fClusters && gAlice->TreeR() && R) {
    gAlice->TreeR()->Branch(branchname,&fClusters, buffersize);
    printf("Making Branch %s for Clusters\n",branchname);
  }	
}
 
//_____________________________________________________________________________
void AliTPC::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  // reset clusters
  //
  fNdigits   = 0;
  if (fDigits)   fDigits->Clear();
  fNclusters = 0;
  if (fClusters) fClusters->Clear();
}

//_____________________________________________________________________________
void AliTPC::SetSecAL(Int_t sec)
{
  //
  // Activate/deactivate selection for lower sectors
  //
  fSecAL = sec;
}

//_____________________________________________________________________________
void AliTPC::SetSecAU(Int_t sec)
{
  //
  // Activate/deactivate selection for upper sectors
  //
  fSecAU = sec;
}

//_____________________________________________________________________________
void AliTPC::SetSecLows(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6)
{
  //
  // Select active lower sectors
  //
  fSecLows[0] = s1;
  fSecLows[1] = s2;
  fSecLows[2] = s3;
  fSecLows[3] = s4;
  fSecLows[4] = s5;
  fSecLows[5] = s6;
}

//_____________________________________________________________________________
void AliTPC::SetSecUps(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6,
                       Int_t s7, Int_t s8 ,Int_t s9 ,Int_t s10, 
                       Int_t s11 , Int_t s12)
{
  // 
  // Select active upper sectors
  //
  fSecUps[0] = s1;
  fSecUps[1] = s2;
  fSecUps[2] = s3;
  fSecUps[3] = s4;
  fSecUps[4] = s5;
  fSecUps[5] = s6;
  fSecUps[6] = s7;
  fSecUps[7] = s8;
  fSecUps[8] = s9;
  fSecUps[9] = s10;
  fSecUps[10] = s11;
  fSecUps[11] = s12;
}

//_____________________________________________________________________________
void AliTPC::SetSens(Int_t sens)
{
  fSens = sens;
}

//_____________________________________________________________________________
void AliTPC::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliTPC.
  //
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      if (R__v < 2) return;
      R__b >> fNsectors;
      R__b >> fNclusters;
      R__b >> fNtracks;
      fClustersIndex = new Int_t[fNsectors+1];
      fDigitsIndex   = new Int_t[fNsectors+1];
   } else {
      R__b.WriteVersion(AliTPC::IsA());
      AliDetector::Streamer(R__b);
      R__b << fNsectors;
      R__b << fNclusters;
      R__b << fNtracks;
   }
}
  
ClassImp(AliTPCcluster)
 
//_____________________________________________________________________________
AliTPCcluster::AliTPCcluster(Float_t *hits, Int_t *lab)
{
  //
  // Creates a simulated cluster for the TPC
  //
  fTracks[0]  = lab[0];
  fTracks[1]  = lab[1];
  fTracks[2]  = lab[2];
  fSector     = lab[3];
  fPadRow     = lab[4];
  fY          = hits[0];
  fZ          = hits[1];
  fQ          = hits[2];
  fSigmaY2    = hits[3];
  fSigmaZ2    = hits[4];
}
 
//_____________________________________________________________________________
void AliTPCcluster::GetXYZ(Double_t& x, Double_t& y, Double_t& z) const 
{
  //
  // Returns centroid for of a simulated TPC cluster
  //
  Double_t alpha,xrow;
  if (fSector<24) {
    alpha=(fSector<12) ? fSector*alpha_low : (fSector-12)*alpha_low;
    xrow=pad_row_low[fPadRow];
  } else {
    alpha=(fSector<48) ? (fSector-24)*alpha_up : (fSector-48)*alpha_up;
    xrow=pad_row_up[fPadRow];
  }
  x=xrow*cos(alpha) - fY*sin(alpha);
  y=xrow*sin(alpha) + fY*cos(alpha);
  z=fZ;
}
 
ClassImp(AliTPCdigit)
 
//_____________________________________________________________________________
AliTPCdigit::AliTPCdigit(Int_t *tracks, Int_t *digits):
  AliDigit(tracks)
{
  //
  // Creates a TPC digit object
  //
  fSector     = digits[0];
  fPadRow     = digits[1];
  fPad        = digits[2];
  fTime       = digits[3];
  fSignal     = digits[4];
}

 
ClassImp(AliTPChit)
 
//_____________________________________________________________________________
AliTPChit::AliTPChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
AliHit(shunt,track)
{
  //
  // Creates a TPC hit object
  //
  fSector     = vol[0];
  fPadRow     = vol[1];
  fX          = hits[0];
  fY          = hits[1];
  fZ          = hits[2];
  fQ          = hits[3];
}
 
 
ClassImp(AliTPCtrack)
 
//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(Float_t *hits)
{
  //
  // Default creator for a TPC reconstructed track object
  //
  ref=hits[0]; // This is dummy code !
}

AliTPCtrack::AliTPCtrack(const AliTPCcluster& c,const TVector& xx,
			 const TMatrix& CC):
  x(xx),C(CC),clusters(MAX_CLUSTER)
{
  //
  // Standard creator for a TPC reconstructed track object
  //
  chi2=0.;
  int sec=c.fSector, row=c.fPadRow;
  if (sec<24) { 
    fAlpha=(sec%12)*alpha_low; ref=pad_row_low[0] + row*pad_pitch_l;
  } else { 
    fAlpha=(sec%24)*alpha_up;  ref=pad_row_up[0]  + row*pad_pitch_l;
  }
  clusters.AddLast((AliTPCcluster*)(&c));
}

//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(const AliTPCtrack& t) : x(t.x), C(t.C),
  clusters(t.clusters.GetEntriesFast()) 
{
  //
  // Copy creator for a TPC reconstructed track
  //
  ref=t.ref;
  chi2=t.chi2;
  fAlpha=t.fAlpha;
  int n=t.clusters.GetEntriesFast();
  for (int i=0; i<n; i++) clusters.AddLast(t.clusters.UncheckedAt(i));
}

//_____________________________________________________________________________
Double_t AliTPCtrack::GetY(Double_t xk) const 
{
  //
  //
  //
  Double_t c2=x(2)*xk - x(3);
  if (TMath::Abs(c2) >= 0.999) {
    if (*this>10) cerr<<*this<<" AliTPCtrack warning: No y for this x !\n";
    return 0.;
  }
  Double_t c1=x(2)*ref - x(3);
  Double_t r1=sqrt(1.-c1*c1), r2=sqrt(1.-c2*c2);
  Double_t dx=xk-ref;
  return x(0) + dx*(c1+c2)/(r1+r2);
}

//_____________________________________________________________________________
int AliTPCtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho,Double_t pm)
{
  //
  // Propagate a TPC reconstructed track
  //
  if (TMath::Abs(x(2)*xk - x(3)) >= 0.999) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Propagation failed !\n";
    return 0;
  }
  
  Double_t x1=ref, x2=x1+0.5*(xk-x1), dx=x2-x1, y1=x(0), z1=x(1);
  Double_t c1=x(2)*x1 - x(3), r1=sqrt(1.- c1*c1);
  Double_t c2=x(2)*x2 - x(3), r2=sqrt(1.- c2*c2);
  
  x(0) += dx*(c1+c2)/(r1+r2);
  x(1) += dx*(c1+c2)/(c1*r2 + c2*r1)*x(4);
  
  TMatrix F(5,5); F.UnitMatrix();
  Double_t rr=r1+r2, cc=c1+c2, xx=x1+x2;
  F(0,2)= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  F(0,3)=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  Double_t cr=c1*r2+c2*r1;
  F(1,2)= dx*x(4)*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);
  F(1,3)=-dx*x(4)*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  F(1,4)= dx*cc/cr; 
  TMatrix tmp(F,TMatrix::kMult,C);
  C.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  
  ref=x2;
  
  //Multiple scattering******************
  Double_t ey=x(2)*ref - x(3);
  Double_t ex=sqrt(1-ey*ey);
  Double_t ez=x(4);
  TMatrix Q(5,5); Q=0.;
  Q(2,2)=ez*ez+ey*ey;   Q(2,3)=-ex*ey;       Q(2,4)=-ex*ez;
  Q(3,2)=Q(2,3);        Q(3,3)= ez*ez+ex*ex; Q(3,4)=-ey*ez;
  Q(4,2)=Q(2,4);        Q(4,3)= Q(3,4);      Q(4,4)=1.;
  
  F=0;
  F(2,2)=-x(2)*ex;          F(2,3)=-x(2)*ey;
  F(3,2)=-ex*(x(2)*ref-ey); F(3,3)=-(1.+ x(2)*ref*ey - ey*ey);
  F(4,2)=-ez*ex;            F(4,3)=-ez*ey;           F(4,4)=1.;
  
  tmp.Mult(F,Q);
  Q.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  
  Double_t p2=GetPt()*GetPt()*(1.+x(4)*x(4));
  Double_t beta2=p2/(p2 + pm*pm);
  Double_t d=sqrt((x1-ref)*(x1-ref)+(y1-x(0))*(y1-x(0))+(z1-x(1))*(z1-x(1)));
  d*=2.;
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho;
  Q*=theta2;
  C+=Q;
  
  //Energy losses************************
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d*rho;
  if (x1 < x2) dE=-dE;
  x(2)*=(1.- sqrt(p2+pm*pm)/p2*dE);
  
  x1=ref; x2=xk; y1=x(0); z1=x(1);
  c1=x(2)*x1 - x(3); r1=sqrt(1.- c1*c1);
  c2=x(2)*x2 - x(3); r2=sqrt(1.- c2*c2);
  
  x(0) += dx*(c1+c2)/(r1+r2);
  x(1) += dx*(c1+c2)/(c1*r2 + c2*r1)*x(4);
  
  F.UnitMatrix();
  rr=r1+r2; cc=c1+c2; xx=x1+x2;
  F(0,2)= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  F(0,3)=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  cr=c1*r2+c2*r1;
  F(1,2)= dx*x(4)*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);
  F(1,3)=-dx*x(4)*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  F(1,4)= dx*cc/cr; 
  tmp.Mult(F,C);
  C.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  
  ref=x2;
  
  return 1;
}

//_____________________________________________________________________________
void AliTPCtrack::PropagateToVertex(Double_t x0,Double_t rho,Double_t pm) 
{
  //
  // Propagate a reconstructed track from the vertex
  //
  Double_t c=x(2)*ref - x(3);
  Double_t tgf=-x(3)/(x(2)*x(0) + sqrt(1-c*c));
  Double_t snf=tgf/sqrt(1.+ tgf*tgf);
  Double_t xv=(x(3)+snf)/x(2);
  PropagateTo(xv,x0,rho,pm);
}

//_____________________________________________________________________________
void AliTPCtrack::Update(const AliTPCcluster *c, Double_t chisq)
{
  //
  // Update statistics for a reconstructed TPC track
  //
  TMatrix H(2,5); H.UnitMatrix();
  TMatrix Ht(TMatrix::kTransposed,H);
  TVector m(2);   m(0)=c->fY; m(1)=c->fZ; 
  TMatrix V(2,2); V(0,0)=c->fSigmaY2; V(0,1)=0.; V(1,0)=0.; V(1,1)=c->fSigmaZ2;

  TMatrix tmp(H,TMatrix::kMult,C);
  TMatrix R(tmp,TMatrix::kMult,Ht); R+=V;
  
  Double_t det=(Double_t)R(0,0)*R(1,1) - (Double_t)R(0,1)*R(1,0);
  R(0,1)=R(0,0); R(0,0)=R(1,1); R(1,1)=R(0,1); 
  R(1,0)*=-1; R(0,1)=R(1,0);
  R*=1./det;
  
  //R.Invert();
  
  TMatrix K(C,TMatrix::kMult,Ht); K*=R;
  
  TVector savex=x;
  x*=H; x-=m; x*=-1; x*=K; x+=savex;
  if (TMath::Abs(x(2)*ref-x(3)) >= 0.999) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Filtering failed !\n";
    x=savex;
    return;
  }
  
  TMatrix saveC=C;
  C.Mult(K,tmp); C-=saveC; C*=-1;
  
  clusters.AddLast((AliTPCcluster*)c);
  chi2 += chisq;
}

//_____________________________________________________________________________
int AliTPCtrack::Rotate(Double_t alpha)
{
  //
  // Rotate a reconstructed TPC track
  //
  fAlpha += alpha;
  
  Double_t x1=ref, y1=x(0);
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=x(2)*ref - x(3);
  
  ref = x1*ca + y1*sa;
  x(0)=-x1*sa + y1*ca;
  x(3)=x(3)*ca + (x(2)*y1 + sqrt(1.- r1*r1))*sa;
  
  Double_t r2=x(2)*ref - x(3);
  if (TMath::Abs(r2) >= 0.999) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Rotation failed !\n";
    return 0;
  }
  
  Double_t y0=x(0) + sqrt(1.- r2*r2)/x(2);
  if ((x(0)-y0)*x(2) >= 0.) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Rotation failed !!!\n";
    return 0;
  }
  
  TMatrix F(5,5); F.UnitMatrix();
  F(0,0)=ca;
  F(3,0)=x(2)*sa; 
  F(3,2)=(y1 - r1*x1/sqrt(1.- r1*r1))*sa; 
  F(3,3)= ca + sa*r1/sqrt(1.- r1*r1);
  TMatrix tmp(F,TMatrix::kMult,C); 
  // Double_t dy2=C(0,0);
  C.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  // C(0,0)+=dy2*sa*sa*r1*r1/(1.- r1*r1);
  // C(1,1)+=dy2*sa*sa*x(4)*x(4)/(1.- r1*r1);
  
  return 1;
}

//_____________________________________________________________________________
void AliTPCtrack::UseClusters() const 
{
  //
  //
  //
  int num_of_clusters=clusters.GetEntriesFast();
  for (int i=0; i<num_of_clusters; i++) {
    //if (i<=14) continue;
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(i);
    c->Use();   
  }
}

//_____________________________________________________________________________
Double_t AliTPCtrack::GetPredictedChi2(const AliTPCcluster *c) const 
{
  //
  // Calculate chi2 for a reconstructed TPC track
  //
  TMatrix H(2,5); H.UnitMatrix();
  TVector m(2);   m(0)=c->fY; m(1)=c->fZ; 
  TMatrix V(2,2); V(0,0)=c->fSigmaY2; V(0,1)=0.; V(1,0)=0.; V(1,1)=c->fSigmaZ2;
  TVector res=x;  res*=H; res-=m; //res*=-1; 
  TMatrix tmp(H,TMatrix::kMult,C);
  TMatrix R(tmp,TMatrix::kMult,TMatrix(TMatrix::kTransposed,H)); R+=V;
  
  Double_t det=(Double_t)R(0,0)*R(1,1) - (Double_t)R(0,1)*R(1,0);
  if (TMath::Abs(det) < 1.e-10) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Singular matrix !\n";
    return 1e10;
  }
  R(0,1)=R(0,0); R(0,0)=R(1,1); R(1,1)=R(0,1); 
  R(1,0)*=-1; R(0,1)=R(1,0);
  R*=1./det;
  
  //R.Invert();
  
  TVector r=res;
  res*=R;
  return r*res;
}

//_____________________________________________________________________________
int AliTPCtrack::GetLab() const 
{
  //
  //
  //
  int lab = 0;
  struct {
    int lab;
    int max;
  } s[MAX_CLUSTER]={{0,0}};
  
  int i;
  int num_of_clusters=clusters.GetEntriesFast();
  for (i=0; i<num_of_clusters; i++) {
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(i);
    lab=c->fTracks[0]; if (lab<0) lab=-lab;
    int j;
    for (j=0; j<MAX_CLUSTER; j++)
      if (s[j].lab==lab || s[j].max==0) break;
    s[j].lab=lab;
    s[j].max++;
  }
  
  int max=0;
  for (i=0; i<num_of_clusters; i++) 
    if (s[i].max>max) {max=s[i].max; lab=s[i].lab;}
  if (lab>0) lab--;
  
  if (1.-float(max)/num_of_clusters > 0.10) return -lab;
  
  if (num_of_clusters < 6) return lab;
  
  max=0;
  for (i=1; i<=6; i++) {
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(num_of_clusters-i);
    if (lab != TMath::Abs(c->fTracks[0])-1
	&& lab != TMath::Abs(c->fTracks[1])-1
	&& lab != TMath::Abs(c->fTracks[2])-1
	) max++;
  }
  if (max>3) return -lab;
  
  return lab;
}

//_____________________________________________________________________________
void AliTPCtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const 
{
  //
  // Get reconstructed TPC track momentum
  //
  Double_t pt=0.3*FIELD/TMath::Abs(x(2))/100; // GeV/c
  Double_t r=x(2)*ref-x(3);
  Double_t y0=x(0) + sqrt(1.- r*r)/x(2);
  px=-pt*(x(0)-y0)*x(2);    //cos(phi);
  py=-pt*(x(3)-ref*x(2));   //sin(phi);
  pz=pt*x(4);
  Double_t tmp=px*TMath::Cos(fAlpha) - py*TMath::Sin(fAlpha);
  py=px*TMath::Sin(fAlpha) + py*TMath::Cos(fAlpha);
  px=tmp;  
}

//_____________________________________________________________________________
//
//     Classes for internal tracking use
//

//_____________________________________________________________________________
void AliTPCRow::InsertCluster(const AliTPCcluster* c) 
{
  //
  // Insert a cluster in the list
  //
  if (num_of_clusters==MAX_CLUSTER_PER_ROW) {
    cerr<<"AliTPCRow::InsertCluster(): Too many clusters !\n"; return;
  }
  if (num_of_clusters==0) {clusters[num_of_clusters++]=c; return;}
  int i=Find(c->fY);
  memmove(clusters+i+1 ,clusters+i,(num_of_clusters-i)*sizeof(AliTPCcluster*));
  clusters[i]=c; num_of_clusters++;
}

//_____________________________________________________________________________
int AliTPCRow::Find(Double_t y) const 
{
  //
  //
  //
  if (y <= clusters[0]->fY) return 0;
  if (y > clusters[num_of_clusters-1]->fY) return num_of_clusters;
  int b=0, e=num_of_clusters-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (y > clusters[m]->fY) b=m+1;
    else e=m; 
  }
  return m;
}
