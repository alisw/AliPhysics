/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  Wrapper for the set of mag.field parameterizations by Chebyshev polinomials  //
//  To obtain the field in cartesian coordinates/components use                  //
//    Field(float* xyz, float* bxyz);                                            //
//  For cylindrical coordinates/components:                                      //
//    FieldCyl(float* rphiz, float* brphiz)                                      //
//                                                                               //
//  For the moment only the solenoid part is parameterized in the volume defined //
//  by R<500, -550<Z<550 cm                                                      //
//                                                                               //
//  The region R<423 cm,  -343.3<Z<481.3 for 30kA and -343.3<Z<481.3 for 12kA    //
//  is parameterized using measured data while outside the Tosca calculation     //
//  is used (matched to data on the boundary of the measurements)                //
//                                                                               //
//  If the querried point is outside the validity region no the return values    //
//  for the field components are set to 0.                                       //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

#include "AliMagFCheb.h"

ClassImp(AliMagFCheb)



//__________________________________________________________________________________________
AliMagFCheb::AliMagFCheb(const char* inputFile) :
    TNamed("Field Map", inputFile),
    fNParamsSol(0),
    fNSegZSol(0),
    fNParamsDip(0),
    fSegZSol(0),
    fSegRSol(0),
    fNSegRSol(0),
    fSegZIdSol(0),
    fMinZSol(0.),
    fMaxZSol(0.),
    fMaxRSol(0.)
{
  Init0();
  LoadData(inputFile);
}

//__________________________________________________________________________________________
AliMagFCheb::~AliMagFCheb()
{
  if (fNParamsSol) {
    delete fParamsSol;
    delete[] fSegZSol;
    delete[] fSegRSol;
    delete[] fNSegRSol;
    delete[] fSegZIdSol;
  }
  //
  // Dipole part ...
  if (fNParamsDip) {
    delete fParamsDip;
  }
}

//__________________________________________________________________________________________
void AliMagFCheb::Init0()
{
  // Solenoid part
  fNParamsSol = 0;
  fNSegZSol   = 0;
  //
  fSegZSol    = 0;
  fSegRSol    = 0;
  //
  fNSegRSol   = 0;
  fSegZIdSol  = 0;
  //
  fMinZSol = fMaxZSol = fMaxRSol = fMaxRSol = 0;
  fParamsSol = 0;
  //
  // Dipole part ...
  fNParamsDip = 0;
  fParamsDip  = 0;
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::AddParamSol(AliCheb3D* param)
{
  // adds new parameterization piece for Sol
  // NOTE: pieces must be added strictly in increasing R then increasing Z order
  //
  if (!fParamsSol) fParamsSol = new TObjArray();
  fParamsSol->Add(param);
  fNParamsSol++;
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::AddParamDip(AliCheb3D* param)
{
  // adds new parameterization piece for Dipole
  //
  if (!fParamsDip) fParamsDip = new TObjArray();
  fParamsDip->Add(param);
  fNParamsDip++;
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::BuildTableSol()
{
  // build the indexes for each parameterization of Solenoid
  //
  const float kSafety=0.001;
  //
  fSegRSol   = new Float_t[fNParamsSol];
  float *tmpbufF  = new float[fNParamsSol+1];
  int   *tmpbufI  = new int[fNParamsSol+1];
  int   *tmpbufI1 = new int[fNParamsSol+1];
  //
  // count number of Z slices and number of R slices in each Z slice
  for (int ip=0;ip<fNParamsSol;ip++) {
    if (ip==0 || (GetParamSol(ip)->GetBoundMax(2)-GetParamSol(ip-1)->GetBoundMax(2))>kSafety) { // new Z slice
      tmpbufF[fNSegZSol] = GetParamSol(ip)->GetBoundMax(2);                                     // 
      tmpbufI[fNSegZSol] = 0;
      tmpbufI1[fNSegZSol++] = ip;
    }
    fSegRSol[ip] = GetParamSol(ip)->GetBoundMax(0);  // upper R
    tmpbufI[fNSegZSol-1]++;
  }
  //
  fSegZSol   = new Float_t[fNSegZSol];
  fSegZIdSol = new Int_t[fNSegZSol];
  fNSegRSol  = new Int_t[fNSegZSol];
  for (int iz=0;iz<fNSegZSol;iz++) {
    fSegZSol[iz]   = tmpbufF[iz];
    fNSegRSol[iz]  = tmpbufI[iz];
    fSegZIdSol[iz] = tmpbufI1[iz];
  }
  //
  fMinZSol = GetParamSol(0)->GetBoundMin(2);
  fMaxZSol = GetParamSol(fNParamsSol-1)->GetBoundMax(2);
  fMaxRSol = GetParamSol(fNParamsSol-1)->GetBoundMax(0);
  //
  delete[] tmpbufF;
  delete[] tmpbufI;
  delete[] tmpbufI1;
  //
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::Field(Float_t *xyz, Float_t *b) const
{
  // compute field in cartesian coordinates
  float rphiz[3];
  if (xyz[2]>GetMaxZSol() || xyz[2]<GetMinZSol()) {for (int i=3;i--;) b[i]=0; return;}
  //
  // Sol region
  // convert coordinates to cyl system
  rphiz[0] = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  rphiz[1] = TMath::ATan2(xyz[1],xyz[0]);
  rphiz[2] = xyz[2];
  if (rphiz[0]>GetMaxRSol()) {for (int i=3;i--;) b[i]=0; return;}
  //
  FieldCylSol(rphiz,b);
  //
  // convert field to cartesian system
  float btr = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]);
  float psiPLUSphi = rphiz[1] + TMath::ATan2(b[1],b[0]);
  b[0] = btr*TMath::Cos(psiPLUSphi);
  b[1] = btr*TMath::Sin(psiPLUSphi);
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::FieldCylSol(Float_t *rphiz, Float_t *b) const
{
  // compute Solenoid field in Cylindircal coordinates
  // note: the check for the point being inside the parameterized region is done outside
  float &r = rphiz[0];
  float &z = rphiz[2];
  int SolZId = 0;
  while (z>fSegZSol[SolZId]) ++SolZId;    // find Z segment
  int SolRId = fSegZIdSol[SolZId];        // first R segment for this Z
  while (r>fSegRSol[SolRId]) ++SolRId;    // find R segment
  GetParamSol( SolRId )->Eval(rphiz,b);
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::Print(Option_t *) const
{
  printf("Alice magnetic field parameterized by Chebyshev polynomials\n");
  printf("Segmentation for Solenoid (%+.2f<Z<%+.2f cm | R<%.2f cm)\n",fMinZSol,fMaxZSol,fMaxRSol);
  //
  for (int iz=0;iz<fNSegZSol;iz++) {
    AliCheb3D* param = GetParamSol( fSegZIdSol[iz] );
    printf("*** Z Segment %2d (%+7.2f<Z<%+7.2f)\t***\n",iz,param->GetBoundMin(2),param->GetBoundMax(2));
    for (int ir=0;ir<fNSegRSol[iz];ir++) {
      param = GetParamSol( fSegZIdSol[iz]+ir );
      printf("    R Segment %2d (%+7.2f<R<%+7.2f, Precision: %.1e) (ID=%2d)\n",ir, param->GetBoundMin(0),param->GetBoundMax(0),
	     param->GetPrecision(),fSegZIdSol[iz]+ir);
    }
  }
}

//_______________________________________________
#ifdef  _INC_CREATION_ALICHEB3D_
void AliMagFCheb::SaveData(const char* outfile) const
{
  // writes coefficients data to output text file
  TString strf = outfile;
  gSystem->ExpandPathName(strf);
  FILE* stream = fopen(strf,"w+");
  //
  // Sol part
  fprintf(stream,"# Set of Chebyshev parameterizations for ALICE magnetic field\nSTART %s\n",GetName());
  fprintf(stream,"START SOLENOID\n#Number of pieces\n%d\n",fNParamsSol);
  for (int ip=0;ip<fNParamsSol;ip++) GetParamSol(ip)->SaveData(stream);
  fprintf(stream,"#\nEND SOLENOID\n");
  //
  // Dip part
  fprintf(stream,"START DIPOLE\n#Number of pieces\n%d\n",fNParamsDip);
  for (int ip=0;ip<fNParamsDip;ip++) GetParamDip(ip)->SaveData(stream);
  fprintf(stream,"#\nEND DIPOLE\n");
  //
  fprintf(stream,"#\nEND %s\n",GetName());
  //
  fclose(stream);
  //
}
#endif

//_______________________________________________
void AliMagFCheb::LoadData(const char* inpfile)
{
  // read coefficients data from the text file
  //
  TString strf = inpfile;
  gSystem->ExpandPathName(strf);
  FILE* stream = fopen(strf,"r");
  if (!stream) {
    printf("Did not find input file %s\n",strf.Data());
    return;
  }
  //
  TString buffs;
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START")) {Error("LoadData","Expected: \"START <name>\", found \"%s\"\nStop\n",buffs.Data());exit(1);}
  if (buffs.First(' ')>0) SetName(buffs.Data()+buffs.First(' ')+1);
  //
  // Solenoid part
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START SOLENOID")) {Error("LoadData","Expected: \"START SOLENOID\", found \"%s\"\nStop\n",buffs.Data());exit(1);}
  AliCheb3DCalc::ReadLine(buffs,stream); // nparam
  int nparSol = buffs.Atoi(); 
  //
  for (int ip=0;ip<nparSol;ip++) {
    AliCheb3D* cheb = new AliCheb3D();
    cheb->LoadData(stream);
    AddParamSol(cheb);
  }
  //
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("END SOLENOID")) {Error("LoadData","Expected \"END SOLENOID\", found \"%s\"\nStop\n",buffs.Data());exit(1);}
  //
  // Dipole part
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START DIPOLE")) {Error("LoadData","Expected: \"START DIPOLE\", found \"%s\"\nStop\n",buffs.Data());exit(1);}
  AliCheb3DCalc::ReadLine(buffs,stream); // nparam
  int nparDip = buffs.Atoi();  
  //
  for (int ip=0;ip<nparDip;ip++) {
    AliCheb3D* cheb = new AliCheb3D();
    cheb->LoadData(stream);
    AddParamDip(cheb);
  }
  //
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("END DIPOLE")) {Error("LoadData","Expected \"END DIPOLE\", found \"%s\"\nStop\n",GetName(),buffs.Data());exit(1);}
  //
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("END") || !buffs.Contains(GetName())) {Error("LoadData","Expected: \"END %s\", found \"%s\"\nStop\n",GetName(),buffs.Data());exit(1);}
  //
  fclose(stream);
  BuildTableSol();
  //  BuildDipTable();
  printf("Loaded magnetic field \"%s\" from %s\n",GetName(),strf.Data());
  //
}
