//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************
#include "AliRICHParam.h"
#include <iostream.h>

ClassImp(AliRICHParam)
Bool_t   AliRICHParam::fgIsWireSag            =kTRUE;
Bool_t   AliRICHParam::fgIsResolveClusters    =kTRUE;
Double_t AliRICHParam::fgAngleRot             =-60;
Int_t    AliRICHParam::fgHV[kNsectors]        ={2150,2100,2050,2000,2150,2150};
Int_t    AliRICHParam::fgNsigmaTh             =4;
Float_t  AliRICHParam::fgSigmaThMean          =1.5;
Float_t  AliRICHParam::fgSigmaThSpread        =0.5;      
Float_t  AliRICHParam::fSigmaThMap[kNCH][kNpadsX][kNpadsY];

void AliRICHParam::GenSigmaThMap()
{
// Generate the map of thresholds sigmas for all pads of all chambers 
  for(Int_t iChamber=0;iChamber<kNCH;iChamber++)
    for(Int_t ipadX=0;ipadX<NpadsX();ipadX++)
      for(Int_t ipadY=0;ipadY<NpadsY();ipadY++) 
        fSigmaThMap[iChamber][ipadX][ipadY] = SigmaThMean()+(1.-2*gRandom->Rndm())*SigmaThSpread();
  Info("GenSigmaThMap"," Threshold map generated for all RICH chambers");
}
//__________________________________________________________________________________________________
void AliRICHParam::Print()
{
  cout<<"\nPads in chamber ("<<NpadsX()<<','<<NpadsY()<<") in sector ("<<NpadsXsec()<<','<<NpadsYsec()<<')'<<endl;
  cout<<"PC size ("<<PcSizeX()<<','<<PcSizeY()<<") sector size ("<<SectorSizeX()<<','<<SectorSizeY()<<") pad size ("<<PadSizeX()<<','
      <<PadSizeY()<<") Dead zone "<<DeadZone()<<endl;
  cout<<"Anode wire pitch "<<WirePitch()<<" Anode-Cathode gap "<<AnodeCathodeGap()<<" Protection wires-cathode gap "<<ProximityGap()<<endl<<endl;
}
//__________________________________________________________________________________________________
