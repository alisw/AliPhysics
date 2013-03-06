// $Id$

//_________________________________________________________________________
// Class for easier handling of track-cluster electron PID.
//
// Author: Tomas Aronsson (Yale)
/////////////////////////////////////////////////////////////////////////// 
  
#include <Riostream.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <AliEMCALGeometry.h>
#include <AliESD.h>
#include <AliESD.h>
#include <AliESDCaloCells.h>
#include <AliESDCaloCluster.h>
#include <AliESDtrack.h>
#include <AliESDtrack.h>
#include "AliEMCALClusterParams.h" 

ClassImp(AliEMCALClusterParams)
  
//____________________________________________________________________________
AliEMCALClusterParams::AliEMCALClusterParams(AliESDtrack *trackin,
                                             AliESDCaloCluster *clusin, 
                                             AliEMCALGeometry *geometryin, 
                                             AliESDCaloCells *cellsin) : 
  TObject(),
  fTrack(trackin), 
  fCluster(clusin), 
  fGeom(geometryin),
  fCells(cellsin)
{
  // Constructor.
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetPe() const
{
  Double_t pe=fTrack->Pt()/fCluster->E();
  return pe;
}

//____________________________________________________________________________
Int_t AliEMCALClusterParams::IsElectron() const
{
  Double_t ep=fCluster->E()/fTrack->Pt();
//  if (fCluster->GetNCells()       < 2    ) return 0;
//  if (fCluster->GetNCells()       > 35   ) return 0;
//  if (fCluster->E()               < 0    ) return 0;
//  if (fCluster->GetDispersion()   > 1.08 ) return 0;
//  if (fCluster->GetM20()          > 0.42 ) return 0;
//  if (fCluster->GetM02()          > 0.4  ) return 0;
//  if (fCluster->GetM20()          < 0    ) return 0;
//  if (fCluster->GetM02()          < 0.06 ) return 0;
  Int_t istpce=0;
  Double_t dEdx=0;
  dEdx=fTrack->GetTPCsignal();
  //if(dEdx0>75.&&dEdx0<95&&fRunnumber>141794&&fRunnumber<146861) istpce0=1;
  //if(dEdx0>56.&&dEdx0<73&&fRunnumber>151564&&fRunnumber<155385) istpce0=1;
  //if(dEdx0 > 70.&&fRunnumber<141795) istpce0=1;
  if(dEdx>70.0&&dEdx<95.0) istpce=1;

  if (ep>0.8&&ep<1.2&&istpce) return 1;
  else return 0;  
}

//____________________________________________________________________________
void AliEMCALClusterParams::LoopThroughCells() const
{
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    
    printf("Cell %d id: %d. ",i,celllist[i]);
    
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    printf("iphi,eta: %d %d, ",iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    printf("smod %d tower %d ieta %d iphi %d \n",iSupMod,iTower,iphi,ieta);
  }
}

//____________________________________________________________________________
void AliEMCALClusterParams::PrintClusterParameters() const
{
  printf("N fCells: %d Energy: %f Dispersion %f M02: %f M20: %f p/E: %f \n",
         fCluster->GetNCells(),fCluster->E(),fCluster->GetDispersion(),
         fCluster->GetM02(),fCluster->GetM20(),fTrack->Pt()/fCluster->E());
}

//===================================== UNWEIGHTED PARAMETERS==================================

//____________________________________________________________________________
void AliEMCALClusterParams::GetCentroid(Double_t &xback, Double_t &yback, Double_t &rback) const
{
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();

  //Calculates mean x,y and r of the fCluster
  Double_t etai=0, phii=0; 
  Double_t xsum=0;
  Double_t ysum=0;
  Double_t rsum=0;
  Double_t esum=0;
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;

    //printf("Cell %d id: %d. ",i,celllist[i]);
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    // printf("iphi,eta: %d %d, ",iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    //printf("smod %d tower %d ieta %d iphi %d \n",iSupMod,iTower,iphi,ieta);

    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;
    xsum+=etai*fCells->GetCellAmplitude(celllist[i]);
    ysum+=phii*fCells->GetCellAmplitude(celllist[i]);
    esum+=fCells->GetCellAmplitude(celllist[i]);
    rsum+=sqrt(etai*etai+phii*phii)*fCells->GetCellAmplitude(celllist[i]);
  }
  yback=ysum/esum;
  xback=xsum/esum;
  rback=rsum/esum;
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetR(Double_t x, Double_t y) const
{
  // Takes fCluster, and cell position (x,y) and returns distance of cell from centroid

  Double_t rmean=-99;
  Double_t xmean=-99;
  Double_t ymean=-99;

  GetCentroid(xmean,ymean,rmean);
  return sqrt(pow((x-xmean),2)+pow((y-ymean),2));
}

//_____________________________________________________________________________
Double_t AliEMCALClusterParams::GetRfactor() const
{
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();

  //Calculates mean x,y and r of the fCluster
  Double_t etai=0, phii=0; 
  Double_t rsum=0;
  Double_t esum=0;
  
  for (Int_t i=0;i<nclusfCells;i++) { 
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;
    esum+=fCells->GetCellAmplitude(celllist[i]);
    rsum+=GetR(etai,phii)*fCells->GetCellAmplitude(celllist[i]);//GetR!
  }
  
  rsum=rsum/esum;
  return rsum;
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::ElectronFraction(Double_t r, Double_t tce) const
{
  // In determination of the 'K-factor', we measure the 'dispersion' of the [energy fraction of the cell
  // as a function of cell distance from the centroid] from the the distribution for electrons
  // the functional form for electrons is saved here.
  // tce is TOTAL FCLUSTER ENERGY, r is distance of cell to centroid.

  Double_t parm[20][11] = {{0.968421,-0.952347,0.196982,-495.028,495.028,0.512768,3.21555,-134.368,950.132,0.258585,-0.0819375},
			 {0.964309,-0.955084,0.212321,-401.273,401.273,0.600581,-2.44701,-43.9126,1009.65,0.19388,-0.0527486},
			 {0.949579,-0.938399,0.169353,-398.208,398.208,0.801709,-16.1048,399.842,-2319.64,0.16054,-0.044954},
			 {0.937618,-0.926665,0.173489,-291.715,291.715,0.723635,-5.65824,222.248,-1525.94,0.125206,-0.0108058},
			 {0.922977,-0.91741,0.174686,-269.155,269.155,0.777305,-0.656328,86.8073,-733.079,0.123508,-0.0159208},
			 {0.88858,-0.86354,0.102552,-280.217,280.217,0.746299,7.06471,-93.7562,323.688,0.128045,-0.0306688},
			 {0.859556,-0.83609,0.108273,-275.632,275.632,0.752012,4.78987,-66.0216,262.093,0.102209,-0.0201003},
			 {0.789254,-0.749649,0.0960812,-331.001,331.001,0.703363,6.70229,-129.452,725.871,0.0861534,-0.0127965},
			 {0.813629,-0.784791,0.0847636,-287.469,287.469,0.7138,7.15787,-142.766,820.634,0.0625165,-0.00910334},
			 {0.828972,-0.805683,0.0812484,-522.356,522.356,0.813623,1.77465,-34.6065,147.756,0.0467835,-0.00425002},
			 {0.900904,-0.890754,0.048338,-1180.84,1180.84,0.870175,0.0385033,-6.99334,19.6102,0.0232503,-0.00446416},
			 {0.881545,-0.870425,0.0382378,-750,750,0.857701,0.181933,-7.72458,12.6329,0.0198669,-0.00436687},
			 {0.893767,-0.878626,0.0465724,-750,750,0.867552,0.0868854,-9.17702,27.6996,0.013101,-0.00261109},
			 {0.88893,-0.882229,0.0896708,-1166.94,828.101,0.864913,-0.245158,-1.5693,-10.5061,0.00658265,-0.000176185},
			 {0.875546,-0.867819,0.0433809,-750,750,0.853546,0.103699,-4.53405,0.67608,0.00667051,-0.00110355},
			 {0.879879,-0.874221,0.0510928,-750,750,0.846377,0.246414,-4.28576,-22.6443,0.00537529,-0.000852894},
			 {0.889361,-0.884282,0.0475369,-750,750,0.845565,0.419865,-7.83763,-2.71567,0.00348465,-0.000407448},
			 {0.893711,-0.890189,0.0649599,-847.351,847.351,0.844507,0.379177,-8.07506,7.69156,0.0005,-2.89084e-13},
			 {0.893706,-0.891587,0.0409521,-750,750,0.841632,0.269466,-4.04288,-14.8299,0.00211674,-0.000173327},
			 {0.942361,-0.941082,0.0703767,-750,750,0.83202,0.829046,-11.3403,29.0905,0.00025,5.67537e-05}};
  
  Double_t par[11]={0};
  if (tce>0&&tce<0.6) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[0][l];
  } else if (tce>0.6&&tce<0.8) {
     for (Int_t l=0;l<11;l++)
      par[l]=parm[1][l];
  } else if (tce>0.8&&tce<1.12) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[2][l];
  } else if (tce>1.12&&tce<1.37) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[3][l];
  } else if (tce>1.37&&tce<1.75) {
     for (Int_t l=0;l<11;l++)
      par[l]=parm[4][l];
  } else if (tce>1.75&&tce<2.5) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[5][l];
  } else if (tce>2.5&&tce<3.5) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[6][l];
  } else if (tce>3.5&&tce<4.5) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[7][l];
  } else if (tce>4.5&&tce<5.5) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[8][l];
  } else if (tce>5.5&&tce<8) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[9][l];
  } else if (tce<8&&tce>12) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[10][l];
  } else if (tce>12&&tce<20) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[11][l];
  } else if (tce>20&&tce<40) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[12][l];
  } else if (tce>40&&tce<63) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[13][l];
  } else if (tce>63&&tce<88) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[14][l];
  } else if (tce>88&&tce<113) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[15][l];
  } else if (tce>113&&tce<138) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[16][l];
  } else if (tce>138&&tce<163) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[17][l];
  } else if (tce>163&&tce<188) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[18][l];
  } else if (tce>188) {
    for (Int_t l=0;l<11;l++)
      par[l]=parm[19][l];
  }

  Double_t fr=0.;
  if (r<0.1) {
    fr =par[5]+par[6]*r+par[7]*r*r+par[8]*r*r*r;
  } else if (r>0.1&&r<0.92) {
    fr = par[0]+par[1]*r+par[2]*exp(par[3]-par[4]*r);
  } else if (r>0.92&&r<6) {
    fr = par[9] + par[10]*r;
  } 
  if (r>6) {
    cout<<"Crappy f(r_i) value, more than 6!"<<endl;
  }

  return fr; //returns fraction of energy at r,E: f(r,E), correct?
}

//____________________________________________________________________________
Double_t  AliEMCALClusterParams::GetKfactor() const
{
  // Determines the k-factor, uses the funtion ElectronFraction(r,totalenergyoffCluster)

  Double_t totalFClusterEnergy = fCluster->E();
  //Calculates mean x,y and r of the fCluster
  Double_t etai=0, phii=0; 

  Double_t kfactor=0;
 
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;

    kfactor+=GetR(etai,phii)*pow(fCells->GetCellAmplitude(celllist[i])/(totalFClusterEnergy) - 
                                 ElectronFraction(GetR(etai,phii),totalFClusterEnergy),2);
  }

  return kfactor;
}

//____________________________________________________________________________
Double_t  AliEMCALClusterParams::GetDispersionX() const
{
  Double_t rmean=-99;
  Double_t xmean=-99;
  Double_t ymean=-99;
  Double_t esum=0;
  Double_t dsum=0;

  GetCentroid(xmean,ymean,rmean);

  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    Double_t etai=(Double_t)ieta+1.;
     
    dsum+=(pow(etai-xmean,2))*fCells->GetCellAmplitude(celllist[i]);
    esum+=fCells->GetCellAmplitude(celllist[i]);
  }

  return sqrt(dsum/esum);
}

//____________________________________________________________________________
Double_t  AliEMCALClusterParams::GetDispersionY() const
{
  Double_t rmean=-99;
  Double_t xmean=-99;
  Double_t ymean=-99;
  Double_t esum=0;
  Double_t dsum=0;

  GetCentroid(xmean,ymean,rmean);

  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
 
    Double_t phii=(Double_t)iphi+1.;
     
    dsum+=(pow(phii-ymean,2))*fCells->GetCellAmplitude(celllist[i]);
    esum+=fCells->GetCellAmplitude(celllist[i]);
  }

  return sqrt(dsum/esum);
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetDispersionMax() const
{
  Double_t dispX = GetDispersionX();
  Double_t dispY = GetDispersionY();
  if (dispY > dispX) {
    return dispY;
  }
  else {
    return dispX;
  }
}

//____________________________________________________________________________
void  AliEMCALClusterParams::GetEllipseParameters(Double_t &param1, Double_t &param2) const
{
  Double_t sumxx=0;
  Double_t sumyy=0;
  Double_t sumx=0;
  Double_t sumy=0;
  Double_t sumxy=0;
  Double_t esum=0;

  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    Double_t phii=(Double_t)iphi+1.;
    Double_t etai=(Double_t)ieta+1.;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);

    sumxx += amp*pow(etai,2);
    sumyy += amp*pow(phii,2);
    sumx += amp*etai;
    sumy += amp*phii;
    sumxy += amp*etai*phii;
    esum+=amp;
  }
  sumxx=sumxx/esum;
  sumyy = sumyy/esum;
  sumx = sumx/esum;
  sumy = sumy/esum;
  sumxx=sumxx-sumx*sumx;
  sumyy = sumyy-sumy*sumy;
  sumxy = sumxy/esum - sumx*sumy;

  param1 = 0.5*(sumxx+sumyy) + sqrt(0.25*pow((sumxx-sumyy),2) + sumxy*sumxy);
  param2 = 0.5*(sumxx+sumyy) - sqrt(0.25*pow((sumxx-sumyy),2) + sumxy*sumxy);
}

//____________________________________________________________________________
Double_t  AliEMCALClusterParams::GetDispersion() const
{
  Double_t rmean=-99;
  Double_t xmean=-99;
  Double_t ymean=-99;
  Double_t esum=0;
  Double_t dsum=0;

  GetCentroid(xmean,ymean,rmean);

  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();

  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    Double_t phii=(Double_t)iphi+1.;
    Double_t etai=(Double_t)ieta+1.;
     
    dsum+=(pow(phii-ymean,2)+pow(etai-xmean,2))*amp;
    esum+=amp;
  }

  return sqrt(dsum/esum);
}

//============================LOG WEIGHTED PARAMETERS========================================================

//____________________________________________________________________________
void  AliEMCALClusterParams::GetWeightedCentroid(Double_t &xback, Double_t &yback, Double_t &rback) const
{
  Float_t logWeight = 4.5;
  Double_t totalFClusterEnergy=fCluster->E();

  //Calculates mean x,y and r of the fCluster
  Double_t etai=0, phii=0; 
  Double_t xsum=0;
  Double_t ysum=0;
  Double_t rsum=0;
  Double_t esum=0;

  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;
    Double_t w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy ));
    xsum+=etai*w;
    ysum+=phii*w;
    esum+=w;
    rsum+=sqrt(etai*etai+phii*phii)*w;

  }
  yback=ysum/esum;
  xback=xsum/esum;
  rback=rsum/esum;
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedR(Double_t x, Double_t y) const
{
  //takes fCluster, and cell position (x,y) and returns distance of cell from
  //centroid
  Double_t rmean=-99;
  Double_t xmean=-99;
  Double_t ymean=-99;

  GetWeightedCentroid(xmean,ymean,rmean);
  return sqrt(pow((x-xmean),2)+pow((y-ymean),2));
}

//_____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedRfactor() const
{
  Double_t rsum=0;
  Double_t esum=0;
  Double_t etai=0, phii=0; 

  Float_t logWeight = 4.5;

  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;
    Double_t w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy ));
    esum+=w;
    rsum+=GetWeightedR(etai,phii)*w;
  }
  rsum=rsum/esum;

  return rsum;
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::ElectronfractionWeighted(Double_t r, Double_t tce) const
{
  // In determination of the 'K-factor', we measure the 'dispersion' of the [energy fraction of the cell
  // as a function of cell distance from the centroid] from the the distribution for electrons
  // the functional form for electrons is saved here.

  Double_t parm[20][10] = {{9.65749,-25.2283,146.575,-207.63,-3.53453,86.296,-176.854,110.51,8.21828,0.097827},
			{9.69895,-19.8042,112.42,-156.762,16.8874,-27.4961,28.4271,-9.31945,8.18897,-0.213499},
			{9.58816,-15.0126,92.8682,-137.142,19.3093,-42.9118,59.9102,-29.9932,8.10929,-0.118502},
			{9.66433,-17.3556,105.064,-154.824,19.3254,-43.845,63.0948,-32.7479,7.3532,0.437176},
			{9.63564,-12.3028,76.0095,-113.329,16.337,-29.9111,42.2855,-22.8323,7.56083,0.14301},
			{9.57174,-8.92634,58.7273,-91.3475,9.81528,2.85217,-10.6939,4.76399,7.93151,-0.320945},
			{9.49951,-3.10053,18.9051,-26.8457,11.5665,-6.0263,3.44864,-2.37582,7.36371,-0.057224},
			{9.43673,-2.63179,13.7039,-19.2528,9.83544,0.91246,-6.30083,2.46559,7.36203,-0.146645},
			{9.5378,-2.35688,12.6874,-19.3792,7.64049,12.0476,-24.2138,11.5242,7.2368,-0.337072},
			{9.52255,-0.989157,7.61241,-14.1598,6.28615,18.4146,-33.4616,15.6595,6.66628,-0.0979268},
			{9.92079,-1.74431,9.58206,-16.8904,20.9936,-35.2392,29.6651,-9.03178,5.74954,-0.082485},
			{9.88709,-1.11991,6.64989,-13.221,15.0524,-15.8677,9.13398,-1.97541,5.5779,-0.146283},
			{9.91925,-1.71077,9.34349,-16.919,15.913,-16.3563,8.1726,-1.38982,5.04338,-0.120024},
			{9.90981,-2.53254,13.2516,-21.0921,15.1033,-14.4442,6.87168,-1.18726,4.25787,0.0881321},
			{9.7964,-0.0431799,2.45754,-9.00884,15.2986,-15.7366,8.37136,-1.66767,3.83555,0.0409475},
			{9.82242,-0.718889,5.19571,-11.6806,15.4565,-16.0309,8.56511,-1.70815,3.46559,0.0840781},
			{9.88062,-1.21927,7.07515,-13.7956,14.2168,-13.5604,6.90633,-1.3678,3.06548,0.162603},
			{9.98456,-2.58626,12.6175,-20.0417,12.9291,-10.2085,4.02983,-0.606252,3.34551,0.040163},
			{9.89892,-1.45526,8.17534,-14.9205,12.5899,-9.74255,4.08869,-0.718114,2.51632,0.245529},
			{11.2705,-5.29437,0.810794,-0.00331736,10.2653,-0.725442,-3.84838,1.23406,3.85982,-0.188577}};

  Double_t par[10]={0};
  if (tce>0&&tce<0.6) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[0][l];
  } else if (tce>0.6&&tce<0.8) {
     for (Int_t l=0;l<10;l++)
      par[l]=parm[1][l];
  } else if (tce>0.8&&tce<1.12) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[2][l];
  } else if (tce>1.12&&tce<1.37) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[3][l];
  } else if (tce>1.37&&tce<1.75) {
     for (Int_t l=0;l<10;l++)
      par[l]=parm[4][l];
  } else if (tce>1.75&&tce<2.5) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[5][l];
  } else if (tce>2.5&&tce<3.5) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[6][l];
  } else if (tce>3.5&&tce<4.5) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[7][l];
  } else if (tce>4.5&&tce<5.5) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[8][l];
  } else if (tce>5.5&&tce<8) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[9][l];
  } else if (tce<8&&tce>12) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[10][l];
  } else if (tce>12&&tce<20) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[8][l];
  } else if (tce>20&&tce<40) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[12][l];
  } else if (tce>40&&tce<63) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[13][l];
  } else if (tce>63&&tce<88) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[14][l];
  } else if (tce>88&&tce<113) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[15][l];
  } else if (tce>113&&tce<138) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[16][l];
  } else if (tce>138&&tce<163) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[17][l];
  } else if (tce>163&&tce<188) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[18][l];
  } else if (tce>188) {
    for (Int_t l=0;l<10;l++)
      par[l]=parm[19][l];
  }

  Double_t fr=0.;
  if (r<0.5) {
    fr = par[4]+par[5]*r+par[6]*r*r+par[7]*r*r*r;
  } else if (r>0.5&&r<1.0) {
    fr = par[0]+par[1]*r+par[2]*r*r+par[3]*r*r*r;
  } else if (r>1.0&&r<6) {
    fr = par[8] + par[9]*r;
  }
  if (r>6) {
    cout<<"Crappy f(r_i) value, more than 6!"<<endl;
  }

  //cout<<"Weighter fr is: "<<fr<<" and tce,r: "<<tce<<" "<<r<<endl;
  return fr;
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedKfactor() const
{
  Double_t logWeight = 4.5;
  //determines the k-factor
  Double_t kfactor=0;

  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    Double_t etai=(Double_t)ieta+1.;
    Double_t phii=(Double_t)iphi+1.;

    Double_t w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy ));
    kfactor+=GetWeightedR(etai,phii)*pow(w - ElectronfractionWeighted(GetWeightedR(etai,phii),totalFClusterEnergy),2);
  }

  return kfactor;
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedDispersionX() const
{
  Float_t logWeight = 4.5;
  // Calculates the dispersion of the shower at the origin of the RecPoint
  // in cell units - Nov 16,2006

  Double_t d = 0., wtot = 0., w = 0.;
  Int_t nstat=0;
	
  // Calculates the dispersion in cell units 
  Double_t etai=0, etaMean=0.0; 

  // Calculate mean values
  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();

  Int_t ncell=0;//cell counter
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    ncell++;
    
    etai=(Double_t)ieta+1.;
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy )); //Calc weight
    if (w>0.0) {
      etaMean += etai*w;
      wtot    += w;
    }
  }

  if (wtot>0)
  etaMean = etaMean/wtot;

  // Calculate dispersion
  // Loop over fCells in the newly created fCluster
  Int_t ncell1=0;//cell counter
  nstat=0;

  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    ncell1++;
    etai=(Double_t)ieta+1.;
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy));
    
    if (w>0.0) {
      nstat++;
      d += w*((etai-etaMean)*(etai-etaMean)); //Add squares
    }
  }
  
  if ( wtot > 0 && nstat>1) d /= wtot;
  else                      d = 0.; 

  return TMath::Sqrt(d);
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedDispersionY() const
{
  Float_t logWeight = 4.5;
  // Calculates the dispersion of the shower at the origin of the RecPoint
  // in cell units - Nov 16,2006

  Double_t d = 0., wtot = 0., w = 0.;
  Int_t nstat=0;
	
  // Calculates the dispersion in cell units 
  Double_t phii=0, phiMean=0.0; 

  // Calculate mean values
  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();

  Int_t ncell=0;//cell counter
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    ncell++;
    
    phii=(Double_t)iphi+1.;
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy )); //Calc weight
    if (w>0.0) {
      phiMean += phii*w;
      wtot    += w;
    }
  }

  if (wtot>0)
  phiMean = phiMean/wtot;

  // Calculate dispersion
  // Loop over fCells in the newly created fCluster
  Int_t ncell1=0;//cell counter
  nstat=0;
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
     
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    ncell1++;

    phii=(Double_t)iphi+1.;
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy));
    
    if (w>0.0) {
      nstat++;
      d += w*((phii-phiMean)*(phii-phiMean)); //Add squares
    }
  }
  
  if ( wtot > 0 && nstat>1) d /= wtot;
  else                      d = 0.; 

  return TMath::Sqrt(d);
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedDispersionMax() const
{
  Double_t dispX = GetWeightedDispersionX();
  Double_t dispY = GetWeightedDispersionY();
  if (dispY > dispX) {
    return dispY;
  }
  else {
    return dispX;
  }
}

//____________________________________________________________________________
void AliEMCALClusterParams::GetWeightedEllipseParameters(Double_t &param1, Double_t &param2) const
{
  Float_t logWeight=4.5;

  Double_t wtot = 0.;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
	
  Double_t etai =0, phii=0, w=0; 

  Int_t ncell=0;//cell counter

  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);

    ncell++;
    etai = phii = 0.; 
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    etai=(Double_t)ieta;
    phii=(Double_t)iphi;

    //Weight!	  
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy ) ); //Energy of tower/total clus E, in GeV
    dxx  += w * etai * etai;
    x    += w * etai;
    dzz  += w * phii * phii;
    z    += w * phii; 
    dxz  += w * etai * phii; 
    wtot += w;
  }

  if ( wtot > 0 ) { 
    dxx /= wtot;
    x   /= wtot;
    dxx -= x * x;
    dzz /= wtot;
    z   /= wtot;
    dzz -= z * z;
    dxz /= wtot;
    dxz -= x * z;
    param1 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ) ;
    param2 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ) ;
  } else { 
    param1= 0.;
    param2= 0.;
  }
}

//____________________________________________________________________________
Double_t AliEMCALClusterParams::GetWeightedDispersion(Double_t &dispersionback) const
{
  Float_t logWeight = 4.5;
  // Calculates the dispersion of the shower at the origin of the RecPoint
  // in cell units - Nov 16,2006

  Double_t d = 0., wtot = 0., w = 0.;
  Int_t nstat=0;

	
  // Calculates the dispersion in cell units 
  Double_t etai, phii, etaMean=0.0, phiMean=0.0; 

  // Calculate mean values
  Int_t ncell=0;//cell counter

  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
  
  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);

    ncell++;
    etai = phii = 0.; 
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
 
    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy )); //Calc weight
    if (w>0.0) {
      phiMean += phii*w;
      etaMean += etai*w;
      wtot    += w;
    }
    //digit->Delete();
  }

  if (wtot>0.0) {
    phiMean = phiMean/wtot;
    etaMean = etaMean/wtot;
  }

  // Calculate dispersion
  Int_t ncell1=0;//cell counter
  nstat=0;

  for (Int_t i=0;i<nclusfCells;i++) {
    Int_t iSupMod = -1;
    Int_t iTower  = -1;
    Int_t iIphi   = -1;
    Int_t iIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);

    ncell1++;
    etai = phii = 0.; 
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    etai=(Double_t)ieta+1.;
    phii=(Double_t)iphi+1.;
    w = TMath::Max(0.,logWeight+TMath::Log(amp/totalFClusterEnergy));
    
    if (w>0.0) {
      nstat++;
      d += w*((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean)); //Add squares
    }
    //delete digit; 
  }
  
  //delete digit;
  
  if ( wtot > 0 && nstat>1) d /= wtot;
  else                      d = 0.; 

  //printf("AliEMCALRecPoint::EvalDispersion() : Dispersion %f \n",d);
  dispersionback=d;
  return TMath::Sqrt(d);

}

//____________________________________________________________________________
void AliEMCALClusterParams::RecalculateClusterShowerShapeParameters(Double_t &m02, Double_t &m20, Double_t &dispersion) const
{
  // Calculates new center of gravity in the local EMCAL-module coordinates 
  // and tranfers into global ALICE coordinates
  // Calculates Dispersion and main axis

  Double_t fW0=4.5;
  
  Int_t nstat  = 0;
  Float_t wtot = 0.;
  Double_t eCell = 0.;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Double_t etai = -1.;
  Double_t phii = -1.;
  
  Double_t w     = 0.;
  Double_t d     = 0.;
  Double_t dxx   = 0.;
  Double_t dzz   = 0.;
  Double_t dxz   = 0.;  
  Double_t xmean = 0.;
  Double_t zmean = 0.;

  Double_t totalFClusterEnergy=fCluster->E();
  Int_t nclusfCells=fCluster->GetNCells();
  UShort_t *celllist;
  celllist=fCluster->GetCellsAbsId();
    
  //Loop on fCells
  for (Int_t i=0;i<nclusfCells;i++) {
    //Get from the absid the supermodule, tower and eta/phi numbers

    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    Double_t amp=fCells->GetCellAmplitude(celllist[i]);
    eCell =amp;

    if (totalFClusterEnergy > 0 && eCell > 0) {
      w  = TMath::Max( 0., fW0 + TMath::Log( eCell / totalFClusterEnergy ));
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;		
      if (w > 0.0) {
        wtot += w;
        nstat++;		        
        //Shower shape
        dxx  += w * etai * etai;
        xmean+= w * etai;
        dzz  += w * phii * phii;
        zmean+= w * phii; 
        dxz  += w * etai * phii; 
      }
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", eCell, fCluster->E()));
  }//cell loop
  
  //Normalize to the weight	
  if (wtot > 0) {
    xmean /= wtot;
    zmean /= wtot;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));
  
  //Calculate dispersion	
  for (Int_t i=0;i<nclusfCells;i++) {
    //Get from the absid the supermodule, tower and eta/phi numbers
    fGeom->GetCellIndex(celllist[i],iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);

    eCell  = fCells->GetCellAmplitude(celllist[i]);
    
    if (totalFClusterEnergy > 0 && eCell > 0) {
      w  = TMath::Max( 0., fW0 + TMath::Log( eCell / totalFClusterEnergy ));
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;		
      if (w > 0.0)  d +=  w*((etai-xmean)*(etai-xmean)+(phii-zmean)*(phii-zmean)); 
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", eCell, fCluster->E()));
  }// cell loop
  
  //Normalize to the weigth and set shower shape parameters
  if (wtot > 0 && nstat > 1) {
    d /= wtot;
    dxx /= wtot;
    dzz /= wtot;
    dxz /= wtot;
    dxx -= xmean * xmean;
    dzz -= zmean * zmean;
    dxz -= xmean * zmean;
    m02=(0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ));
    m20=(0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz ));
  }
  else{
    d=0.;
    m02=0;
    m20=0;
  }	
  
  if (d>=0)
    dispersion=TMath::Sqrt(d);
  else    
    dispersion=0;
}
