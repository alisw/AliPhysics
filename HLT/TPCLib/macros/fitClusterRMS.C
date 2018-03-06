#include "TFile.h"
#include "TNtuple.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TH1.h"
#include "TCanvas.h"

//
//clX:clY:clZ:angle:trkX:trkY:trkZ:trkSinPhi:trkDzDs:trkQPt:trkSigmaY2:trkSigmaZ2:trkSigmaSinPhi2:trkSigmaDzDs2:trkSigmaQPt2
//

class fitEntry
{
 public:
  fitEntry() :fN(0), fSum(0), fRMS(0) {}
  int fN;
  double fSum;
  double fRMS;  
};

class fitMap
{  
public :

  fitMap( int nZ,  double zMin, double zMax, int nA, double aMin, double aMax )
    : fNZ(nZ), fNA(nA), fN(nZ*nA), fZMin(zMin), fZMax(zMax), fAMin(aMin), fAMax(aMax)
  {
    fZScale = (fZMax-fZMin)/fNZ;
    fZScaleInv = fNZ/(fZMax-fZMin);
    fAScale = (fAMax-fAMin)/fNA;
    fAScaleInv = fNA/(fAMax-fAMin);
  }


  int N() const { return fN; }
  int NZ() const { return fNZ; }
  int NA() const { return fNA; }

  int ind( double z, double a ) const
  {
    int indA = (int) ( (a-fAMin)*fAScaleInv );
    int indZ = (int) ( (z-fZMin)*fZScaleInv );
    if( indA<0 || indA>=fNA ) return -1;
    if( indZ<0 || indZ>=fNZ ) return -1;
    return indA*fNZ + indZ;
  }
  
  double A(int ind ) const
  {
    int indA = (int) ind/fNZ;
    return fAMin + (0.5+indA)*fAScale;
  }

  double Z(int ind ) const
  {
    int indZ =  ind % fNZ;
    return fZMin + (0.5+indZ)*fZScale;
  }

private:
  int fNZ;
  int fNA;
  int fN;
  double fZMin;
  double fZMax;
  double fAMin;
  double fAMax;
  double fZScale;
  double fZScaleInv;
  double fAScale;
  double fAScaleInv;
};

float GetRMS( float z, float a, float c[6], float *f=0 )
{
  float a2 = a*a;
  float z2 = z*z;  
  float ff[6] = {1,z,a2,z*z,a2*a2,z*a2}; 
  float ret=0;
  for( int i=0; i<6; i++ ) ret+=ff[i]*c[i];
  if( f ) for( int i=0; i<6; i++ ) f[i] = ff[i];
  if( ret<.01 ) ret = .01;
  return ret;
}

int fitClusterRMS()
{
  gSystem->Load("libAliHLTTPC.so");

  TCanvas *canv = new TCanvas("qa","qa",1000,500);
  canv->Divide(2,2);
  TH1D *qaY = new TH1D("qaY","qaY",1000,-15.,15.);
  TH1D *qaZ = new TH1D("qaZ","qaZ",1000,-15.,15.);
  TH1D *qaYcut = new TH1D("qaYcut","qaYcut",1000,-15.,15.);
  TH1D *qaZcut = new TH1D("qaZcut","qaZcut",1000,-15.,15.);

  double xborder[2] = { 133.9, 198.8 };
  double zMax = 250-0.275;


  const char *fname = "clusterresMC-large.root";
  //const char *fname = "clusterresMC.root";

  TFile *file = new TFile(fname,"READ");

  if ( !file || !file->IsOpen() ){
    printf("Can't open [%s] !\n", fname);
    return 1;
  }

  TNtuple *nt = (TNtuple*) file->FindObjectAny("clusterres");
  if ( !nt ) { 
    printf("Error getting TList clusterres \n");
    return 1;
  }

  // Initial errors

  float Coeff[2][3][6] =    
 { 
   { 
 { 3.45964059234e-02, 5.48224896193e-04, -9.15433093905e-03, -4.29471299412e-07, 1.29013538361e+00, -7.61212781072e-04 } 
,  { 3.79444323480e-02, 6.82618469000e-04, 2.89491593838e-01, -1.38612142564e-06, -1.47421330214e-01, -3.49121750332e-04 } 
,  { 6.79186061025e-02, 1.26425642520e-04, 4.24163460732e-01, -5.37607895978e-08, -9.83583211899e-01, 6.61917729303e-05 } 
   } 
,    { 
 { 2.20002025366e-01, -1.63084571250e-03, -6.21695369482e-02, 5.10420568389e-06, -8.36216099560e-03, 1.31321803201e-03 } 
,  { 4.83078956604e-02, 4.03122794523e-05, 1.65162652731e-01, 7.54337804665e-07, -3.89994867146e-02, 4.65555611299e-04 } 
,  { 1.98864527047e-02, 3.98120639147e-04, 2.95489847660e-01, -5.17226112606e-07, -3.12503017485e-02, -1.46960781422e-04 } 
   } 
 }; 








 



  /*
 float Coeff[2][3][6] =    
    { 
      { 
	{ 6.91527724266e-02, 1.93114537979e-04, 7.59759685025e-03, 3.84449322155e-07, 3.53577643633e-01, 4.26960672485e-04,  }, 
	{ 7.93934836984e-02, 7.68070385675e-05, 6.94227740169e-02, 2.42764173208e-07, 3.35806518793e-01, 2.95838603051e-04,  }, 
	{ 9.71398428082e-02, -1.92082690774e-04, 6.25126063824e-02, 8.83549716946e-07, 4.99859154224e-01, 2.41756439209e-04,  }, 
      }, 
      { 
	{ 5.80004863441e-02, 3.41877675964e-05, 1.76747351885e-01, 9.43496047512e-07, -5.34947626293e-02, 2.04607480555e-04,  }, 
	{ 5.73962479830e-02, -6.51067384752e-05, 2.45050489902e-01, 9.39387859944e-07, -6.58090710640e-02, 1.77358349902e-04,  }, 
	{ 4.98481430113e-02, 1.18782438221e-04, 3.84890645742e-01, 2.01719785764e-07, -1.27915903926e-01, -3.80322249839e-04,  }, 
      }
    }; 
    */  
 

  fitMap map(250,0.,zMax, 100,0.,0.9);
  fitEntry *arr[6];
  for( int i=0; i<6; i++ ){
    arr[i] = new fitEntry[map.N()];  
    if( !arr[i] ){
      cout<<"Not enough memory!!"<<endl;
      return -1;
    }
  }

  float clX, clY, clZ, angle;
  float trkX, trkY, trkZ, trkSinPhi, trkDzDs, trkQPt;
  float trkSigmaY2, trkSigmaZ2, trkSigmaSinPhi2, trkSigmaDzDs2, trkSigmaQPt2;

  nt->SetBranchAddress("clX",&clX);
  nt->SetBranchAddress("clY",&clY);
  nt->SetBranchAddress("clZ",&clZ);
  nt->SetBranchAddress("angle",&angle);
  nt->SetBranchAddress("trkX",&trkX);
  nt->SetBranchAddress("trkY",&trkY);
  nt->SetBranchAddress("trkZ",&trkZ);
  nt->SetBranchAddress("trkSinPhi",&trkSinPhi);
  nt->SetBranchAddress("trkDzDs",&trkDzDs);
  nt->SetBranchAddress("trkQPt",&trkQPt);
  nt->SetBranchAddress("trkSigmaY2",&trkSigmaY2);
  nt->SetBranchAddress("trkSigmaZ2",&trkSigmaZ2);

  for( int i=0; i<map.N(); i++ ){   
    for( int k=0; k<6; k++ ) arr[k][i].fRMS = -1;
  }


  for( int iter=0; iter<1; iter++ ){

    for( int k=0; k<6; k++ ){
      for( int i=0; i<map.N(); i++ ){   
	arr[k][i].fN=0;
	arr[k][i].fSum=0;
      }      
    }
 
    double cy = 0.778/1.007;
    double cz = 0.895/0.987;

    for( int i=0; i<nt->GetEntriesFast(); i++){
      if( i%1000000==0 ){
	cout<<"pass " <<iter<<", processing "<<i<<" out of "<<nt->GetEntriesFast()
	    <<"  ("<<((long)i)*100/nt->GetEntriesFast()<<" %)"<<endl;
	cout<<" cy = "<<cy<<" cz = "<<cz<<endl;
	canv->cd(1);
	qaY->Draw();
	canv->cd(2);
	qaZ->Draw();
	canv->cd(3);
	qaYcut->Draw();
	canv->cd(4);
	qaZcut->Draw();	
	canv->Update();
      }
      int ret = nt->GetEntry(i);
      if( ret<=0 ) {
	cout<<"Wrong entry, ret == "<<ret<<endl;
	continue;
      }
      if( trkY == 0 ) continue; // bug
      if( fabs( clX - trkX)>1.e-4){
	cout<<"Wrong X: track "<<trkX<<" cluster "<<clX<<endl;
	continue;
      }
      int type = 0;
      if( clX < xborder[0] ) type = 0;
      else if ( clX >= xborder[1] ) type = 2;
      else type = 1;

      //if( type!=0 ) continue; 
      if( fabs(trkSigmaY2)>.05*.05 ) continue;
      if( fabs(trkSigmaZ2)>.05*.05 ) continue;
      if( fabs(trkQPt)>1./.2 ) continue;
      double z = 250. - 0.275 - fabs(trkZ);
      double aY = fabs( trkSinPhi );
      double aZ = fabs( trkDzDs );      
      double dY = ( clY - trkY );
      double dZ = ( clZ - trkZ );
      int indY = map.ind( z, aY );
      if( indY<0 ) continue;
      int indZ = map.ind( z, aZ );
      if( indZ<0 ) continue;
      float rmsY = GetRMS( z, aY, Coeff[0][type] );
      float rmsZ = GetRMS( z, aZ, Coeff[1][type] );

      // correct dY, dZ for track errors
      dY = dY / sqrt( rmsY*rmsY + trkSigmaY2 ) * rmsY ;
      dZ = dZ / sqrt( rmsZ*rmsZ + trkSigmaZ2 ) * rmsZ ;

      // fill normalised errors to QA histo
      float rmsY3 = cy*rmsY;
      float rmsZ3 = cz*rmsZ;
      
      if( fabs(dY)<=10*rmsY ) qaY->Fill(dY/rmsY);
      if( fabs(dZ)<=10*rmsZ ) qaZ->Fill(dZ/rmsZ);
      if( fabs(dY)<=3*rmsY ) qaYcut->Fill(dY/rmsY3);
      if( fabs(dZ)<=3*rmsZ ) qaZcut->Fill(dZ/rmsZ3);
      
      double cutY = 10*rmsY;//3*arr[type][indY].fRMS;
      if( cutY<0 ) cutY = 2;
      if( fabs(dY)<=cutY ){
	arr[type][indY].fSum+=dY*dY;
	arr[type][indY].fN++;
      }

      double cutZ = 10*rmsZ;//3*arr[3+type][indZ].fRMS;
      if( cutZ<0 ) cutZ = 2;
      if( fabs(dZ)<=cutZ ){      
	arr[3+type][indZ].fSum+=dZ*dZ;
	arr[3+type][indZ].fN++;
      }
    }

    for( int k=0; k<6; k++ ){
      for( int i=0; i<map.N(); i++ ){
	if( arr[k][i].fN>=1000 ){
	  arr[k][i].fRMS = sqrt(arr[k][i].fSum/arr[k][i].fN);
	} else {
	  arr[k][i].fRMS = -1;
	}
      }
    }
  }
  
  canv->cd(1);
  qaY->Draw();
  canv->cd(2);
  qaZ->Draw();
  canv->cd(3);
  qaYcut->Draw();
  canv->cd(4);
  qaZcut->Draw();

  canv->Update();

  TFile *fileOut = new TFile("rms.root","RECREATE");
  fileOut->cd();
  for( int k=0; k<6; k++ ){
    TString name = "rms";
    name+=k;
    TNtuple *ntRMS = new TNtuple(name.Data(),name.Data(),"z:a:rms"); 
    for( int i=0; i<map.N(); i++ ){
      if( arr[k][i].fRMS > 0 ){	
	ntRMS->Fill(map.Z(i), map.A(i),arr[k][i].fRMS ); 
      }
    }  
    ntRMS->Write();
  }

  //n->Draw("clY-trkY","abs(trkY)!=0&&abs(clY-trkY)<2 && clX<134","h");
  
  //TH3F *hY = new TH3F("hY","residuals Y",100,0,.9,(int)zMax*3,1,zMax-1,100,-2.,2.);
  //clusterres->Draw("abs(trkSinPhi):250-0.275-abs(trkZ):clY-trkY>>hY","abs(trkY)!=0&&clX<134","h");
  // clusterres->Draw("clY-trkY","abs(trkY)!=0&&abs(clY-trkY)<20 && clX<134 && abs(clZ)<20 && abs(abs(trkSinPhi)-.4)<.1 && abs(abs(trkDzDs)-.1)<.1","h");
  return 0;
}
