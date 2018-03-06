#include "TFile.h"
#include "TNtuple.h"
#include "Riostream.h"
#include "TSystem.h"
#include <iomanip>
#include <limits>
//#include "AliHLTTPCPolynomFit.h"


int fitClusterErrors()
{
  gSystem->Load("libAliHLTTPC.so");
  
  const char *fname = "rms.root";

  TFile *file = new TFile(fname,"READ");

  if ( !file || !file->IsOpen() ){
    printf("Can't open [%s] !\n", fname);
    return 1;
  }

  float Coefficients[6][6];

  for( int k=0; k<6; k++ ){
    double corr = 1;
    //if( k<3 ) corr = 2.48 / 1.35; // correction to 3*rms cut
    //else corr = 2.149/ 1.362;    

    TString name = "rms";
    name+=k;

    TNtuple *ntY = (TNtuple*) file->FindObjectAny(name.Data());
    if ( !ntY ) { 
      printf("Error getting TList cluster rms \n");
      return 1;
    }

    float z, a, rms;

    ntY->SetBranchAddress("z",&z);
    ntY->SetBranchAddress("a",&a);
    ntY->SetBranchAddress("rms",&rms);

    // v = c[0] + c[1]*z + c[2]*angle2 + c[3]*z*z +c[4]*angle2*angle2 + c[5]*z*angle2;    
    AliHLTTPCPolynomFit fit(6);
    
    for( int i=0; i<ntY->GetEntriesFast(); i++ ){
      if( i%100==0 ) cout<<"ntuple "<<k<<" of 6: processed "<<i<<" out of "<<ntY->GetEntriesFast()
			 <<"  ("<<((long)i)*100/ntY->GetEntriesFast()<<" %)"<<endl;    
      ntY->GetEntry(i);
      double a2=a*a;
      float f[6] = {1,z,a2,z*z,a2*a2,z*a2};
      fit.AddMeasurement( f, rms*corr);
    }

    fit.Fit( Coefficients[k] ); 
  }

  cout<<std::scientific;
  cout<<std::setprecision( 11 );
  cout<<"fParamS0Par[2][3][6]="<<endl;
  cout<<" { "<<endl;
  int yztype=0;
  for( int i=0; i<2; i++ ){
    if(i>0 ) cout <<", "; 
    cout<<"   { "<<endl;   
    for( int j=0; j<3; j++, yztype++){
      if(j>0 ) cout <<", "; 
      cout<<" { ";   
      for( int k=0; k<6; k++){
	if(k>0 ) cout <<", "; 
	//cout<<fParamS0Par[i][j][k]<<", "; 
	cout<<Coefficients[yztype][k];
      }
      cout<<" } "<<endl;   
    }
    cout<<"   } "<<endl;
  }
  cout<<" }; "<<endl;
  
  return 0;
}

/*
fParamS0Par[2][3][7]=
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
   }, 
 }; 

 */
