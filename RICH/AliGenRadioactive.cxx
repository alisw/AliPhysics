#include "AliGenRadioactive.h"
#include <TPDGCode.h>
#include <AliRun.h>
#include <TDatabasePDG.h>
#include <AliLog.h>

//__________________________________________________________________________________________________
AliGenRadioactive::AliGenRadioactive(Int_t iSrcNucleus,Int_t iNsecondaries):AliGenerator()
{
//Main ctor. Used to define radioactive source.
  Double_t e[100],a[100];//arrays to store experimental points
  Int_t nPoints; 
  switch(iSrcNucleus){
    case kSr90: fPartId=kElectron; nPoints=46;  //experimental part
  a[ 0]=0.08605; a[ 1]=0.0878;  a[ 2]=0.08705; a[ 3]=0.07855; a[ 4]=0.0709; a[ 5]=0.0647;  a[ 6]=0.05015; a[ 7]=0.0372; a[ 8]=0.0268;  a[ 9]=0.0215;
  a[10]=0.0157;  a[11]=0.01685; a[12]=0.01745; a[13]=0.01645; a[14]=0.0175; a[15]=0.01635; a[16]=0.01825; a[17]=0.0177; a[18]=0.01735; a[19]=0.0161;
  a[20]=0.0159;  a[21]=0.0176;  a[22]=0.01605; a[23]=0.0161;  a[24]=0.01495;a[25]=0.01595; a[26]=0.01525; a[27]=0.0138; a[28]=0.0121;  a[29]=0.0101;
  a[30]=0.01175; a[31]=0.01095; a[32]=0.0089;  a[33]=0.0091;  a[34]=0.0625; a[35]=0.0505;  a[36]=0.0475;  a[37]=0.0039; a[38]=0.0031;  a[39]=0.0028;
  a[40]=0.0025;  a[41]=0.0017;  a[42]=4.5e-4;  a[43]=4.5e-4;  a[44]=1.5e-4; a[45]=0;     
  break;
    default:    AliError("Wrong source nucleus specified");  return;   
  }
  for(Int_t i=0;i<nPoints;i++) e[i]=0.001*(i*0.05+0.025); //kinetic energy GeV
  fGenH1=new TH1F("Sr90","Sr90 generator hist",nPoints-1,0,e[nPoints-1]);  
  for(Int_t i=0;i<nPoints;i++) fGenH1->Fill(e[i],a[i]);
  fNpart=iNsecondaries;
}
//__________________________________________________________________________________________________
void AliGenRadioactive::Generate()
{
// Generate one trigger
  Int_t nt=0;
  Double_t ekin=0,p=0,theta=0,phi=0,x=0,y=0,z=0,px=0,py=0,pz=0,polx=0,poly=0,polz=0;
  Double_t m=gAlice->PDGDB()->GetParticle(fPartId)->Mass();
  for(Int_t i=0;i<fNpart;i++){
    x=fOrigin.At(0)+fOsigma.At(0)*(Rndm()-0.5); y=fOrigin.At(1)+fOsigma.At(1)*(Rndm()-0.5); z=fOrigin.At(2)+fOsigma.At(2)*(Rndm()-0.5);
    ekin=fGenH1->GetRandom();  p=TMath::Sqrt(ekin*(2*m+ekin)); 
    theta=Rndm()*fThetaMax; phi=Rndm()*fPhiMax;
    px=p*TMath::Cos(theta)*TMath::Cos(phi);  py=p*TMath::Cos(theta)*TMath::Sin(phi);  pz=p*TMath::Sin(theta);
//   AliDebug(1,Form("Origin=(%5.2f,%5.2f,%5.2f) Ekin=%5.3fMeV,P=%5.3fMeV (%5.2f,%5.2f,%5.2f)",x,y,z,1000*ekin,1000*p,1000*px,1000*py,1000*pz));
    PushTrack(fTrackIt,-1,fPartId,px,py,pz,ekin+m,
                                   x, y, z,0,
                                    polx=0,poly=0,polz=0,kPPrimary,nt);//cm, GeV
  }
}
