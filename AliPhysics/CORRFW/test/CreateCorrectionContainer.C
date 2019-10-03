AliCFContainer* CreateCorrectionContainer(const char* string) {

  const UInt_t nVar   = 5 ; //number of variables on the grid 
                            //(pt=0, eta=1, phi=2, Z_Vertex=3, event multiplicity=4)
  const UInt_t nStep  = 2 ; //number of selection steps :2 (ESD=0,AOD=1)
  const UInt_t nBin[nVar] = //number of bins for each variable
    {
      10 , //bins in pt
      20 , //bins in eta
      20 , //bins in phi
      20 , //bins in Z_Vertex
      20   //bins in multiplicity
    };
  
  //values for bin lower bounds
  Double_t valuesMinMax[nVar][2] ;
  valuesMinMax[0][0] =   0.0 ; // pt min (GeV/c)
  valuesMinMax[0][1] =  10.0 ; // pt max (GeV/c)
  valuesMinMax[1][0] = - 2.5 ; // eta min
  valuesMinMax[1][1] =   2.5 ; // eta max
  valuesMinMax[2][0] =   0.0 ; // phi min (rad)
  valuesMinMax[2][1] = 2*TMath::Pi() ; // phi max (rad)
  valuesMinMax[3][0] = -20.0 ; // zVtx min (cm)
  valuesMinMax[3][1] =  20.0 ; // zVtx max (cm)
  valuesMinMax[4][0] =   0.0 ; // zVtx min (cm)
  valuesMinMax[4][1] = 200.0 ; // zVtx max (cm)

  //arrays for lower bounds :
  Double_t *binLim[nVar] ;
  for (UInt_t iVar=0; iVar<nVar; iVar++) binLim[iVar] = new Double_t[nBin[iVar]+1];

  // fill bin limits using uniform binning between the defined min and max values.
  for (UInt_t iVar=0; iVar<nVar; iVar++) 
    for (UInt_t iBin=0; iBin<=nBin[iVar]; iBin++) 
      binLim[iVar][iBin] = valuesMinMax[iVar][0] + (valuesMinMax[iVar][1]-valuesMinMax[iVar][0])/nBin[iVar]*iBin ;

  //for non-uniform binning :

  //pt
  binLim[0][0]  =  0.0 ;
  binLim[0][1]  =  0.2 ;
  binLim[0][2]  =  0.4 ;
  binLim[0][3]  =  0.7 ;
  binLim[0][4]  =  1.0 ;
  binLim[0][5]  =  1.5 ;
  binLim[0][6]  =  2.0 ;
  binLim[0][7]  =  3.0 ;
  binLim[0][8]  =  4.0 ;
  binLim[0][9]  =  6.0 ;
  binLim[0][10] = 10.0 ;
  
  //multiplicity
  binLim[4][0]  =   0.0 ;
  binLim[4][1]  =   2.0 ;
  binLim[4][2]  =   4.0 ;
  binLim[4][3]  =   6.0 ;
  binLim[4][4]  =   8.0 ;
  binLim[4][5]  =  10.0 ;
  binLim[4][6]  =  15.0 ;
  binLim[4][7]  =  20.0 ;
  binLim[4][8]  =  25.0 ;
  binLim[4][9]  =  30.0 ;
  binLim[4][10] =  35.0 ;
  binLim[4][11] =  40.0 ;
  binLim[4][12] =  50.0 ;
  binLim[4][13] =  60.0 ;
  binLim[4][14] =  70.0 ;
  binLim[4][15] =  80.0 ;
  binLim[4][16] = 100.0 ;
  binLim[4][17] = 120.0 ;
  binLim[4][18] = 140.0 ;
  binLim[4][19] = 170.0 ;
  binLim[4][20] = 200.0 ;
  
  // create a container
  char contName[50];
  char contDesc[50];
  sprintf(contName,"container_%s",string);
  sprintf(contDesc,"correction container for %s",string);
  AliCFContainer* container = new AliCFContainer (contName,contDesc,nStep,nVar,(Int_t*)nBin) ;
  // set the bin limits
  for (UInt_t iVar=0; iVar<nVar; iVar++) container->SetBinLimits(iVar,binLim[iVar]);
  return container;
}
