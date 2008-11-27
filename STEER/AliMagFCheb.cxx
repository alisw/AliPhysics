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

#include "AliMagFCheb.h"
#include <TSystem.h>

ClassImp(AliMagFCheb)

//__________________________________________________________________________________________
AliMagFCheb::AliMagFCheb() : 
  fNParamsSol(0),
  fNSegZSol(0),
  fNParamsTPCInt(0),
  fNSegZTPCInt(0),
  fNParamsDip(0),
//
  fNZSegDip(0),
  fNYSegDip(0),
  fNXSegDip(0),
//
  fSegZSol(0),
  fSegRSol(0),
//
  fSegZTPCInt(0),
  fSegRTPCInt(0),
//
  fSegZDip(0),
  fSegYDip(0),
  fSegXDip(0),
//
  fNSegRSol(0),
  fSegZIdSol(0),
//
  fNSegRTPCInt(0),
  fSegZIdTPCInt(0),
//
  fBegSegYDip(0),
  fNSegYDip(0),
  fBegSegXDip(0),
  fNSegXDip(0),
  fSegIDDip(0),
//
  fMinZSol(1e6),
  fMaxZSol(-1e6),
  fMaxRSol(-1e6), 
//
  fMinZDip(1e6),
  fMaxZDip(-1e6),
//
  fMinZTPCInt(1e6),
  fMaxZTPCInt(-1e6),
  fMaxRTPCInt(-1e6), 
//
  fParamsSol(0),
  fParamsDip(0),
  fParamsTPCInt(0)
//
{
  // default constructor
}

//__________________________________________________________________________________________
AliMagFCheb::AliMagFCheb(const AliMagFCheb& src) : 
  TNamed(src),
  fNParamsSol(0),
  fNSegZSol(0),
  fNParamsTPCInt(0),
  fNSegZTPCInt(0),
  fNParamsDip(0),
//
  fNZSegDip(0),
  fNYSegDip(0),
  fNXSegDip(0),
//
  fSegZSol(0),
  fSegRSol(0),
//
  fSegZTPCInt(0),
  fSegRTPCInt(0),
//
  fSegZDip(0),
  fSegYDip(0),
  fSegXDip(0),
//
  fNSegRSol(0),
  fSegZIdSol(0),
//
  fNSegRTPCInt(0),
  fSegZIdTPCInt(0),
//
  fBegSegYDip(0),
  fNSegYDip(0),
  fBegSegXDip(0),
  fNSegXDip(0),
  fSegIDDip(0),
//
  fMinZSol(1e6),
  fMaxZSol(-1e6),
  fMaxRSol(-1e6), 
//
  fMinZDip(1e6),
  fMaxZDip(-1e6),
//
  fMinZTPCInt(1e6),
  fMaxZTPCInt(-1e6),
  fMaxRTPCInt(-1e6), 
//
  fParamsSol(0),
  fParamsDip(0),
  fParamsTPCInt(0)
{
  // copy constructor
  CopyFrom(src);
}

//__________________________________________________________________________________________
void AliMagFCheb::CopyFrom(const AliMagFCheb& src) 
{ 
  // copy method
  Clear();
  SetName(src.GetName());
  SetTitle(src.GetTitle());
  fNParamsSol    = src.fNParamsSol;
  fNSegZSol      = src.fNSegZSol;
  fNParamsTPCInt = src.fNParamsTPCInt;
  fNSegZTPCInt   = src.fNSegZTPCInt; 
  fNParamsDip    = src.fNParamsDip;
  //
  fNZSegDip      = src.fNZSegDip;
  fNYSegDip      = src.fNYSegDip;
  fNXSegDip      = src.fNXSegDip;  
  //
  fMinZSol       = src.fMinZSol; 
  fMaxZSol       = src.fMaxZSol;
  fMaxRSol       = src.fMaxRSol; 
  //
  fMinZDip       = src.fMinZDip;
  fMaxZDip       = src.fMaxZDip;
  //
  fMinZTPCInt    = src.fMinZTPCInt;
  fMaxZTPCInt    = src.fMaxZTPCInt;
  fMaxRTPCInt    = src.fMaxRTPCInt; 
  // 
  if (src.fNParamsSol) {
    memcpy(fSegZSol  = new Float_t[fNSegZSol], src.fSegZSol, sizeof(Float_t)*fNSegZSol);
    memcpy(fSegRSol  = new Float_t[fNParamsSol], src.fSegRSol, sizeof(Float_t)*fNParamsSol);
    memcpy(fNSegRSol = new Int_t[fNSegZSol], src.fNSegRSol, sizeof(Int_t)*fNSegZSol);
    memcpy(fSegZIdSol= new Int_t[fNSegZSol], src.fSegZIdSol, sizeof(Int_t)*fNSegZSol);
    fParamsSol       = new TObjArray(fNParamsSol);
    for (int i=0;i<fNParamsSol;i++) fParamsSol->AddAtAndExpand(new AliCheb3D(*src.GetParamSol(i)),i);
  }
  //
  if (src.fNParamsDip) {
    memcpy(fSegZDip   = new Float_t[fNZSegDip], src.fSegZDip, sizeof(Float_t)*fNZSegDip);
    memcpy(fSegYDip   = new Float_t[fNYSegDip], src.fSegYDip, sizeof(Float_t)*fNYSegDip);
    memcpy(fSegXDip   = new Float_t[fNXSegDip], src.fSegZDip, sizeof(Float_t)*fNXSegDip);
    memcpy(fBegSegYDip= new Int_t[fNZSegDip], src.fBegSegYDip, sizeof(Int_t)*fNZSegDip);
    memcpy(fNSegYDip  = new Int_t[fNZSegDip], src.fNSegYDip, sizeof(Int_t)*fNZSegDip);
    memcpy(fBegSegXDip= new Int_t[fNYSegDip], src.fBegSegXDip, sizeof(Int_t)*fNYSegDip);
    memcpy(fNSegXDip  = new Int_t[fNYSegDip], src.fNSegXDip, sizeof(Int_t)*fNYSegDip);
    memcpy(fSegIDDip  = new Int_t[fNXSegDip], src.fSegIDDip, sizeof(Int_t)*fNXSegDip);
    fParamsDip        = new TObjArray(fNParamsDip);
    for (int i=0;i<fNParamsDip;i++) fParamsDip->AddAtAndExpand(new AliCheb3D(*src.GetParamDip(i)),i);
  }
  //
  if (src.fNParamsTPCInt) {
    memcpy(fSegZTPCInt  = new Float_t[fNSegZTPCInt], src.fSegZTPCInt, sizeof(Float_t)*fNSegZTPCInt);
    memcpy(fSegRTPCInt  = new Float_t[fNParamsTPCInt], src.fSegRTPCInt, sizeof(Float_t)*fNParamsTPCInt);
    memcpy(fNSegRTPCInt = new Int_t[fNSegZTPCInt], src.fNSegRTPCInt, sizeof(Int_t)*fNSegZTPCInt);
    memcpy(fSegZIdTPCInt= new Int_t[fNSegZTPCInt], src.fSegZIdTPCInt, sizeof(Int_t)*fNSegZTPCInt);
    fParamsTPCInt       = new TObjArray(fNParamsTPCInt);
    for (int i=0;i<fNParamsTPCInt;i++) fParamsTPCInt->AddAtAndExpand(new AliCheb3D(*src.GetParamTPCInt(i)),i);
  }
  //
}

//__________________________________________________________________________________________
AliMagFCheb& AliMagFCheb::operator=(const AliMagFCheb& rhs)
{
  // assignment
  if (this != &rhs) {  
    Clear();
    CopyFrom(rhs);
  }
  return *this;  
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::Clear(const Option_t *)
{
  // clear all dynamic parts
  if (fNParamsSol) {
    delete   fParamsSol;
    delete[] fSegZSol;
    delete[] fSegRSol;
    delete[] fNSegRSol;
    delete[] fSegZIdSol;
  }
  //
  if (fNParamsTPCInt) {
    delete   fParamsTPCInt;
    delete[] fSegZTPCInt;
    delete[] fSegRTPCInt;
    delete[] fNSegRTPCInt;
    delete[] fSegZIdTPCInt;
  }
  // 
  if (fNParamsDip) {
    delete   fParamsDip;
    delete[] fSegZDip;
    delete[] fSegYDip;
    delete[] fSegXDip;
    delete[] fBegSegYDip;
    delete[] fNSegYDip;
    delete[] fBegSegXDip;
    delete[] fNSegXDip;
    delete[] fSegIDDip;
  }
  fNParamsSol = fNParamsTPCInt = fNParamsDip = fNZSegDip = fNYSegDip = fNXSegDip = 0;
  fNSegZSol = fNSegZTPCInt = 0;
  fMinZSol = fMinZDip = fMinZTPCInt = 1e6;
  fMaxZSol = fMaxZDip = fMaxZTPCInt = fMaxRSol = fMaxRTPCInt = -1e6;
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::Field(Float_t *xyz, Float_t *b) const
{
  // compute field in cartesian coordinates. If point is outside of the parameterized region
  // get it at closest valid point
  static float rphiz[3];
  //
#ifndef _BRING_TO_BOUNDARY_  // exact matching to fitted volume is requested
  if ( !(xyz[2]>=GetMinZSol()&&xyz[2]<=GetMaxZSol()) && 
       !(xyz[2]>=GetMinZDip()&&xyz[2]<=GetMaxZDip())  ) {for (int i=3;i--;) b[i]=0; return;}
#endif
  //
  if (xyz[2]<fMaxZDip) {    // dipole part?
#ifndef _BRING_TO_BOUNDARY_
    AliCheb3D* par = GetParamDip(FindDipSegment(xyz));
    if (par->IsInside(xyz)) {par->Eval(xyz,b); return;}
    for (int i=3;i--;) b[i]=0; return;
#else
    GetParamDip(FindDipSegment(xyz))->Eval(xyz,b); return;  
#endif
  }
  //
  // Sol region: convert coordinates to cyl system
  CartToCyl(xyz,rphiz);
#ifndef _BRING_TO_BOUNDARY_
  if (rphiz[0]>GetMaxRSol()) {for (int i=3;i--;) b[i]=0; return;}
#endif
  //
  FieldCylSol(rphiz,b);
  //
  // convert field to cartesian system
  CylToCartCylB(rphiz, b,b);
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::GetTPCInt(Float_t *xyz, Float_t *b) const
{
  // compute TPC region field integral in cartesian coordinates.
  // If point is outside of the parameterized region get it at closeset valid point
  static float rphiz[3];
  //
  // TPCInt region
  // convert coordinates to cyl system
  CartToCyl(xyz,rphiz);
#ifndef _BRING_TO_BOUNDARY_
  if ( (rphiz[2]>GetMaxZTPCInt()||rphiz[2]<GetMinZTPCInt()) ||
       rphiz[0]>GetMaxRTPCInt()) {for (int i=3;i--;) b[i]=0; return;}
#endif
  //
  GetTPCIntCyl(rphiz,b);
  //
  // convert field to cartesian system
  CylToCartCylB(rphiz, b,b);
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::FieldCylSol(const Float_t *rphiz, Float_t *b) const
{
  // compute Solenoid field in Cylindircal coordinates
  // note: if the point is outside the volume get the field in closest parameterized point
  const float &r = rphiz[0];
  const float &z = rphiz[2];
  int SolZId = 0;
  while (z>fSegZSol[SolZId] && SolZId<fNSegZSol-1) ++SolZId;    // find Z segment
  int SolRId = fSegZIdSol[SolZId];        // first R segment for this Z
  int SolRMax = SolRId + fNSegRSol[SolZId];
  while (r>fSegRSol[SolRId] && SolRId<SolRMax-1) ++SolRId;    // find R segment
  GetParamSol( SolRId )->Eval(rphiz,b);
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::GetTPCIntCyl(Float_t *rphiz, Float_t *b) const
{
  // compute field integral in TPC region in Cylindircal coordinates
  // note: the check for the point being inside the parameterized region is done outside
  const float &r = rphiz[0];
  const float &z = rphiz[2];
  int tpcIntZId = 0;
  while (z>fSegZTPCInt[tpcIntZId] && tpcIntZId<fNSegZTPCInt) ++tpcIntZId;    // find Z segment
  int tpcIntRId = fSegZIdTPCInt[tpcIntZId];        // first R segment for this Z
  int tpcIntRIdMax = tpcIntRId + fNSegRTPCInt[tpcIntZId];
  while (r>fSegRTPCInt[tpcIntRId] && tpcIntRId<tpcIntRIdMax) ++tpcIntRId;    // find R segment
  GetParamTPCInt( tpcIntRId )->Eval(rphiz,b);
  //
}


//__________________________________________________________________________________________
void AliMagFCheb::Print(Option_t *) const
{
  // print info
  printf("Alice magnetic field parameterized by Chebyshev polynomials\n");
  printf("Segmentation for Solenoid (%+.2f<Z<%+.2f cm | R<%.2f cm)\n",fMinZSol,fMaxZSol,fMaxRSol);
  //
  if (fParamsSol) fParamsSol->Print();
  /*
  for (int iz=0;iz<fNSegZSol;iz++) {
    AliCheb3D* param = GetParamSol( fSegZIdSol[iz] );
    printf("*** Z Segment %2d (%+7.2f<Z<%+7.2f)\t***\n",iz,param->GetBoundMin(2),param->GetBoundMax(2));
    for (int ir=0;ir<fNSegRSol[iz];ir++) {
      param = GetParamSol( fSegZIdSol[iz]+ir );
      printf("    R Segment %2d (%+7.2f<R<%+7.2f, Precision: %.1e) (ID=%2d)\n",ir, param->GetBoundMin(0),
	     param->GetBoundMax(0),param->GetPrecision(),fSegZIdSol[iz]+ir);
    }
  }
  */
  //
  printf("Segmentation for TPC field integral (%+.2f<Z<%+.2f cm | R<%.2f cm)\n",fMinZTPCInt,fMaxZTPCInt,fMaxRTPCInt);
  //
  if (fParamsTPCInt) fParamsTPCInt->Print();
  /*
  for (int iz=0;iz<fNSegZTPCInt;iz++) {
    AliCheb3D* param = GetParamTPCInt( fSegZIdTPCInt[iz] );
    printf("*** Z Segment %2d (%+7.2f<Z<%+7.2f)\t***\n",iz,param->GetBoundMin(2),param->GetBoundMax(2));
    for (int ir=0;ir<fNSegRTPCInt[iz];ir++) {
      param = GetParamTPCInt( fSegZIdTPCInt[iz]+ir );
      printf("    R Segment %2d (%+7.2f<R<%+7.2f, Precision: %.1e) (ID=%2d)\n",ir, param->GetBoundMin(0),
	     param->GetBoundMax(0),param->GetPrecision(),fSegZIdTPCInt[iz]+ir);
    }
  }
  */
  //
  printf("Segmentation for Dipole (%+.2f<Z<%+.2f cm)\n",fMinZDip,fMaxZDip);
  if (fParamsDip) fParamsDip->Print();
  //
  

}
#ifdef  _INC_CREATION_ALICHEB3D_
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
  if (!buffs.BeginsWith("START")) {
    Error("LoadData","Expected: \"START <name>\", found \"%s\"\nStop\n",buffs.Data());
    exit(1);
  }
  if (buffs.First(' ')>0) SetName(buffs.Data()+buffs.First(' ')+1);
  //
  // Solenoid part    -----------------------------------------------------------
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START SOLENOID")) {
    Error("LoadData","Expected: \"START SOLENOID\", found \"%s\"\nStop\n",buffs.Data());
    exit(1);
  }
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
  if (!buffs.BeginsWith("END SOLENOID")) {
    Error("LoadData","Expected \"END SOLENOID\", found \"%s\"\nStop\n",buffs.Data());
    exit(1);
  }
  //
  // TPCInt part     -----------------------------------------------------------
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START TPCINT")) {
    Error("LoadData","Expected: \"START TPCINT\", found \"%s\"\nStop\n",buffs.Data());
    exit(1);
  }
  AliCheb3DCalc::ReadLine(buffs,stream); // nparam
  int nparTPCInt = buffs.Atoi(); 
  //
  for (int ip=0;ip<nparTPCInt;ip++) {
    AliCheb3D* cheb = new AliCheb3D();
    cheb->LoadData(stream);
    AddParamTPCInt(cheb);
  }
  //
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("END TPCINT")) {
    Error("LoadData","Expected \"END TPCINT\", found \"%s\"\nStop\n",buffs.Data());
    exit(1);
  }
  //
  // Dipole part    -----------------------------------------------------------
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("START DIPOLE")) {
    Error("LoadData","Expected: \"START DIPOLE\", found \"%s\"\nStop\n",buffs.Data());
    exit(1);
  }
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
  if (!buffs.BeginsWith("END DIPOLE")) {
    Error("LoadData","Expected \"END DIPOLE\", found \"%s\"\nStop\n",GetName(),buffs.Data());
    exit(1);
  }
  //
  AliCheb3DCalc::ReadLine(buffs,stream);
  if (!buffs.BeginsWith("END") || !buffs.Contains(GetName())) {
    Error("LoadData","Expected: \"END %s\", found \"%s\"\nStop\n",GetName(),buffs.Data());
    exit(1);
  }
  //
  // ---------------------------------------------------------------------------
  fclose(stream);
  BuildTableSol();
  BuildTableTPCInt();
  BuildTableDip();
  printf("Loaded magnetic field \"%s\" from %s\n",GetName(),strf.Data());
  //
}
#endif

//_________________________________________________________________________
Int_t AliMagFCheb::FindDipSegment(const float *xyz) const
{
  // find the segment containing point xyz. If it is outside find the closest segment 
  int xid,yid,zid = TMath::BinarySearch(fNZSegDip,fSegZDip,xyz[2]); // find zsegment
  int ysegBeg = fBegSegYDip[zid];
  //
  for (yid=0;yid<fNSegYDip[zid];yid++) if (xyz[1]<fSegYDip[ysegBeg+yid]) break;
  if ( --yid < 0 ) yid = 0;
  yid +=  ysegBeg;
  //
  int xsegBeg = fBegSegXDip[yid];
  for (xid=0;xid<fNSegXDip[yid];xid++) if (xyz[0]<fSegXDip[xsegBeg+xid]) break;
  if ( --xid < 0) xid = 0;
  xid +=  xsegBeg;
  //
  return fSegIDDip[xid];
}

//_______________________________________________
#ifdef  _INC_CREATION_ALICHEB3D_


//__________________________________________________________________________________________
AliMagFCheb::AliMagFCheb(const char* inputFile) : 
  fNParamsSol(0),
  fNSegZSol(0),
  fNParamsTPCInt(0),
  fNSegZTPCInt(0),
  fNParamsDip(0),
//
  fNZSegDip(0),
  fNYSegDip(0),
  fNXSegDip(0),
//
  fSegZSol(0),
  fSegRSol(0),
//
  fSegZTPCInt(0),
  fSegRTPCInt(0),
//
  fSegZDip(0),
  fSegYDip(0),
  fSegXDip(0),
//
  fNSegRSol(0),
  fSegZIdSol(0),
//
  fNSegRTPCInt(0),
  fSegZIdTPCInt(0),
//
  fBegSegYDip(0),
  fNSegYDip(0),
  fBegSegXDip(0),
  fNSegXDip(0),
  fSegIDDip(0),
//
  fMinZSol(0),
  fMaxZSol(0),
  fMaxRSol(0), 
//
  fMinZDip(0),
  fMaxZDip(0),
//
  fMinZTPCInt(0),
  fMaxZTPCInt(0),
  fMaxRTPCInt(0), 
//
  fParamsSol(0),
  fParamsDip(0),
  fParamsTPCInt(0)
//
//
{
  // construct from coeffs from the text file
  LoadData(inputFile);
}

//__________________________________________________________________________________________
void AliMagFCheb::AddParamSol(const AliCheb3D* param)
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
void AliMagFCheb::AddParamTPCInt(const AliCheb3D* param)
{
  // adds new parameterization piece for TPCInt
  // NOTE: pieces must be added strictly in increasing R then increasing Z order
  //
  if (!fParamsTPCInt) fParamsTPCInt = new TObjArray();
  fParamsTPCInt->Add(param);
  fNParamsTPCInt++;
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::AddParamDip(const AliCheb3D* param)
{
  // adds new parameterization piece for Dipole
  //
  if (!fParamsDip) fParamsDip = new TObjArray();
  fParamsDip->Add(param);
  fNParamsDip++;
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::ResetTPCInt()
{
  // clean TPC field integral (used for update)
  if (!fNParamsTPCInt) return;
  delete fParamsTPCInt;  
  delete[] fSegZTPCInt;
  delete[] fSegRTPCInt;
  delete[] fNSegRTPCInt;
  delete[] fSegZIdTPCInt;
  //
  fNParamsTPCInt = 0; 
  fNSegZTPCInt   = 0;
  fSegZTPCInt    = 0; 
  fSegRTPCInt    = 0; 
  fNSegRTPCInt   = 0; 
  fSegZIdTPCInt  = 0;
  fMinZTPCInt    = 0; 
  fMaxZTPCInt    = 0; 
  fMaxRTPCInt    = 0; 
  fParamsTPCInt  = 0;
  //
}

//__________________________________________________
void AliMagFCheb::BuildTableDip()
{
  // build lookup table for dipole
  //
  TArrayF segY,segX;
  TArrayI begSegYDip,begSegXDip;
  TArrayI nsegYDip,nsegXDip;
  TArrayI segID;
  float *tmpSegZ,*tmpSegY,*tmpSegX;
  //
  // create segmentation in Z
  fNZSegDip = SegmentDipDimension(&tmpSegZ, fParamsDip, fNParamsDip, 2, 1,-1, 1,-1, 1,-1) - 1;
  fNYSegDip = 0;
  fNXSegDip = 0;
  //
  // for each Z slice create segmentation in Y
  begSegYDip.Set(fNZSegDip);
  nsegYDip.Set(fNZSegDip);
  float xyz[3];
  for (int iz=0;iz<fNZSegDip;iz++) {
    printf("\nZSegment#%d  %+e : %+e\n",iz,tmpSegZ[iz],tmpSegZ[iz+1]);
    int ny = SegmentDipDimension(&tmpSegY, fParamsDip, fNParamsDip, 1, 
				 1,-1, 1,-1, tmpSegZ[iz],tmpSegZ[iz+1]) - 1;
    segY.Set(ny + fNYSegDip);
    for (int iy=0;iy<ny;iy++) segY[fNYSegDip+iy] = tmpSegY[iy];
    begSegYDip[iz] = fNYSegDip;
    nsegYDip[iz] = ny;
    printf(" Found %d YSegments, to start from %d\n",ny, begSegYDip[iz]);
    //
    // for each slice in Z and Y create segmentation in X
    begSegXDip.Set(fNYSegDip+ny);
    nsegXDip.Set(fNYSegDip+ny);
    xyz[2] = (tmpSegZ[iz]+tmpSegZ[iz+1])/2.; // mean Z of this segment
    //
    for (int iy=0;iy<ny;iy++) {
      int isg = fNYSegDip+iy;
      printf("\n   YSegment#%d  %+e : %+e\n",iy, tmpSegY[iy],tmpSegY[iy+1]);
      int nx = SegmentDipDimension(&tmpSegX, fParamsDip, fNParamsDip, 0, 
				   1,-1, tmpSegY[iy],tmpSegY[iy+1], tmpSegZ[iz],tmpSegZ[iz+1]) - 1;
      //
      segX.Set(nx + fNXSegDip);
      for (int ix=0;ix<nx;ix++) segX[fNXSegDip+ix] = tmpSegX[ix];
      begSegXDip[isg] = fNXSegDip;
      nsegXDip[isg] = nx;
      printf("   Found %d XSegments, to start from %d\n",nx, begSegXDip[isg]);
      //
      segID.Set(fNXSegDip+nx);
      //
      // find corresponding params
      xyz[1] = (tmpSegY[iy]+tmpSegY[iy+1])/2.; // mean Y of this segment
      //
      for (int ix=0;ix<nx;ix++) {
	xyz[0] = (tmpSegX[ix]+tmpSegX[ix+1])/2.; // mean X of this segment
	for (int ipar=0;ipar<fNParamsDip;ipar++) {
	  AliCheb3D* cheb = (AliCheb3D*) fParamsDip->At(ipar);
	  if (!cheb->IsInside(xyz)) continue;
	  segID[fNXSegDip+ix] = ipar;
	  break;
	}
      }
      fNXSegDip += nx;
      //
      delete[] tmpSegX;
    }
    delete[] tmpSegY;
    fNYSegDip += ny;
  }
  //
  fMinZDip = tmpSegZ[0];
  fMaxZDip = tmpSegZ[fNZSegDip];
  fSegZDip    = new Float_t[fNZSegDip];
  for (int i=fNZSegDip;i--;) fSegZDip[i] = tmpSegZ[i];
  delete[] tmpSegZ;
  //
  fSegYDip    = new Float_t[fNYSegDip];
  fSegXDip    = new Float_t[fNXSegDip];
  fBegSegYDip = new Int_t[fNZSegDip];
  fNSegYDip   = new Int_t[fNZSegDip];
  fBegSegXDip = new Int_t[fNYSegDip];
  fNSegXDip   = new Int_t[fNYSegDip];
  fSegIDDip   = new Int_t[fNXSegDip];
  //
  for (int i=fNYSegDip;i--;) fSegYDip[i] = segY[i];
  for (int i=fNXSegDip;i--;) fSegXDip[i] = segX[i];
  for (int i=fNZSegDip;i--;) {fBegSegYDip[i] = begSegYDip[i]; fNSegYDip[i] = nsegYDip[i];}
  for (int i=fNYSegDip;i--;) {fBegSegXDip[i] = begSegXDip[i]; fNSegXDip[i] = nsegXDip[i];}
  for (int i=fNXSegDip;i--;) {fSegIDDip[i]   = segID[i];}
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::BuildTableSol()
{
  // build the indexes for each parameterization of Solenoid
  //
  const float kSafety=0.001;
  //
  if (fNParamsSol<1) return;
  fNSegZSol = 0;
  fMaxRSol = 0;
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
    if (fMaxRSol<fSegRSol[ip]) fMaxRSol = fSegRSol[ip];
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
  //
  delete[] tmpbufF;
  delete[] tmpbufI;
  delete[] tmpbufI1;
  //
  //
}

//__________________________________________________________________________________________
void AliMagFCheb::BuildTableTPCInt()
{
  // build the indexes for each parameterization of TPC field integral
  //
  const float kSafety=0.001;
  //
  if (fNParamsTPCInt<1) return;
  fNSegZTPCInt = 0;
  fSegRTPCInt   = new Float_t[fNParamsTPCInt];
  float *tmpbufF  = new float[fNParamsTPCInt+1];
  int   *tmpbufI  = new int[fNParamsTPCInt+1];
  int   *tmpbufI1 = new int[fNParamsTPCInt+1];
  //
  // count number of Z slices and number of R slices in each Z slice
  for (int ip=0;ip<fNParamsTPCInt;ip++) {
    if (ip==0 || (GetParamTPCInt(ip)->GetBoundMax(2)-GetParamTPCInt(ip-1)->GetBoundMax(2))>kSafety) { // new Z slice
      tmpbufF[fNSegZTPCInt] = GetParamTPCInt(ip)->GetBoundMax(2);                                     // 
      tmpbufI[fNSegZTPCInt] = 0;
      tmpbufI1[fNSegZTPCInt++] = ip;
    }
    fSegRTPCInt[ip] = GetParamTPCInt(ip)->GetBoundMax(0);  // upper R
    tmpbufI[fNSegZTPCInt-1]++;
  }
  //
  fSegZTPCInt   = new Float_t[fNSegZTPCInt];
  fSegZIdTPCInt = new Int_t[fNSegZTPCInt];
  fNSegRTPCInt  = new Int_t[fNSegZTPCInt];
  for (int iz=0;iz<fNSegZTPCInt;iz++) {
    fSegZTPCInt[iz]   = tmpbufF[iz];
    fNSegRTPCInt[iz]  = tmpbufI[iz];
    fSegZIdTPCInt[iz] = tmpbufI1[iz];
  }
  //
  fMinZTPCInt = GetParamTPCInt(0)->GetBoundMin(2);
  fMaxZTPCInt = GetParamTPCInt(fNParamsTPCInt-1)->GetBoundMax(2);
  fMaxRTPCInt = GetParamTPCInt(fNParamsTPCInt-1)->GetBoundMax(0);
  //
  delete[] tmpbufF;
  delete[] tmpbufI;
  delete[] tmpbufI1;
  //
  //
}

void AliMagFCheb::SaveData(const char* outfile) const
{
  // writes coefficients data to output text file
  TString strf = outfile;
  gSystem->ExpandPathName(strf);
  FILE* stream = fopen(strf,"w+");
  //
  // Sol part    ---------------------------------------------------------
  fprintf(stream,"# Set of Chebyshev parameterizations for ALICE magnetic field\nSTART %s\n",GetName());
  fprintf(stream,"START SOLENOID\n#Number of pieces\n%d\n",fNParamsSol);
  for (int ip=0;ip<fNParamsSol;ip++) GetParamSol(ip)->SaveData(stream);
  fprintf(stream,"#\nEND SOLENOID\n");
  //
  // TPCInt part ---------------------------------------------------------
  fprintf(stream,"# Set of Chebyshev parameterizations for ALICE magnetic field\nSTART %s\n",GetName());
  fprintf(stream,"START TPCINT\n#Number of pieces\n%d\n",fNParamsTPCInt);
  for (int ip=0;ip<fNParamsTPCInt;ip++) GetParamTPCInt(ip)->SaveData(stream);
  fprintf(stream,"#\nEND TPCINT\n");
  //
  // Dip part   ---------------------------------------------------------
  fprintf(stream,"START DIPOLE\n#Number of pieces\n%d\n",fNParamsDip);
  for (int ip=0;ip<fNParamsDip;ip++) GetParamDip(ip)->SaveData(stream);
  fprintf(stream,"#\nEND DIPOLE\n");
  //
  fprintf(stream,"#\nEND %s\n",GetName());
  //
  fclose(stream);
  //
}

Int_t AliMagFCheb::SegmentDipDimension(float** seg,const TObjArray* par,int npar, int dim, 
				       float xmn,float xmx,float ymn,float ymx,float zmn,float zmx)
{
  // find all boundaries in deimension dim for boxes in given region.
  // if mn>mx for given projection the check is not done for it.
  float *tmpC = new float[2*npar];
  int *tmpInd = new int[2*npar];
  int nseg0 = 0;
  for (int ip=0;ip<npar;ip++) {
    AliCheb3D* cheb = (AliCheb3D*) par->At(ip);
    if (xmn<xmx && (cheb->GetBoundMin(0)>(xmx+xmn)/2 || cheb->GetBoundMax(0)<(xmn+xmx)/2)) continue;
    if (ymn<ymx && (cheb->GetBoundMin(1)>(ymx+ymn)/2 || cheb->GetBoundMax(1)<(ymn+ymx)/2)) continue;
    if (zmn<zmx && (cheb->GetBoundMin(2)>(zmx+zmn)/2 || cheb->GetBoundMax(2)<(zmn+zmx)/2)) continue;
    //
    tmpC[nseg0++] = cheb->GetBoundMin(dim);
    tmpC[nseg0++] = cheb->GetBoundMax(dim);    
  }
  // range Dim's boundaries in increasing order
  TMath::Sort(nseg0,tmpC,tmpInd,kFALSE);
  // count number of really different Z's
  int nseg = 0;
  float cprev = -1e6;
  for (int ip=0;ip<nseg0;ip++) {
    if (TMath::Abs(cprev-tmpC[ tmpInd[ip] ])>1e-4) {
      cprev = tmpC[ tmpInd[ip] ];
      nseg++;
    }
    else tmpInd[ip] = -1; // supress redundant Z
  }
  // 
  *seg  = new float[nseg]; // create final Z segmenations
  nseg = 0;
  for (int ip=0;ip<nseg0;ip++) if (tmpInd[ip]>=0) (*seg)[nseg++] = tmpC[ tmpInd[ip] ];
  //
  delete[] tmpC;
  delete[] tmpInd;
  return nseg;
}

#endif

