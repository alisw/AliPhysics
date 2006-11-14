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


/***********************************************************************
*  This code defines the reconstructed v0 visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/

#include "Track.h"
#include "V0.h"
#include "MCHelixLine.hi"

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

#include <Reve/RGTopFrame.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include <vector>

using namespace Reve;



/***********************************************************************
*
*  V0 class
*
************************************************************************/

const Float_t V0::fgkMassPion2 = 0.13956995*0.13956995;
const Float_t V0::fgkMassProton2 = 0.93827231*0.93827231;

ClassImp(Reve::V0)

V0::V0() :
  RenderElement(),
  TPolyMarker3D(1),
  fV_neg(),
  fP_neg(),
  fV_pos(),
  fP_pos(),
  fV_v0(),
  fV0_birth(),
  fBeta_neg(0),
  fBeta_pos(0),
  fLabel_neg(0),
  fLabel_pos(0),
  fPathMarksNeg(),
  fPathMarksPos(),
  fRnrStyle(0),
  fPolyLineNeg(),
  fPolyLinePos(),
  fPolyLineV0(),
  fESDIndex(-1),
  fDaughterDCA(999),
  fCosPointingAngle(999),
  fDecayLength(0)
{}


V0::V0(Reve::RecTrack* tNeg, Reve::RecTrack* tPos,
       Reve::RecV0* v0, TrackRnrStyle* rs) :
  RenderElement(),
  TPolyMarker3D(1),
  fV_neg(v0->V_neg),
  fP_neg(v0->P_neg ),
  fV_pos(v0->V_pos ),
  fP_pos(v0->P_pos),
  fV_v0(v0->V_ca),
  fV0_birth(v0->V0_birth),
  fBeta_neg(tNeg->beta),
  fBeta_pos(tPos->beta),
  fLabel_neg(v0->d_label[0]),
  fLabel_pos(v0->d_label[1]),
  fPathMarksNeg(),
  fPathMarksPos(),
  fRnrStyle(rs),
  fPolyLineNeg(),
  fPolyLinePos(),
  fPolyLineV0(),
  fESDIndex(-1),
  fDaughterDCA(999),
  fCosPointingAngle(999),
  fDecayLength(0)
{
  fMarkerColor = fRnrStyle->GetColor();
  fPolyLineV0.SetLineColor(fMarkerColor);
  fPolyLinePos.SetLineColor(2);  // red
  fPolyLineNeg.SetLineColor(7);  // light blue

  fMainColorPtr = &fMarkerColor;
  fMarkerStyle = 20;
  fMarkerColor = 5;
  fMarkerSize = 0.3;
}


V0::~V0()
{
  for (vpPathMark_i i=fPathMarksNeg.begin(); i!=fPathMarksNeg.end(); ++i)
    delete *i;
  for (vpPathMark_i i=fPathMarksPos.begin(); i!=fPathMarksPos.end(); ++i)
    delete *i;
}


void V0::Reset(TPolyLine3D* polyLine) {
  //polyLine->SetPolyLine(n_points);
  polyLine->SetPolyLine(0);
}

//______________________________________________________________________
void V0::SetDecayLength(Float_t primx, Float_t primy, Float_t primz) {

  Float_t dx = fV_v0.x-primx;
  Float_t dy = fV_v0.y-primy;
  Float_t dz = fV_v0.z-primz;

  fDecayLength = sqrt(dx*dx+dy*dy+dz*dz);

  // This is probably wrong but I can only do this for now
  Float_t distNorm = fDecayLength/GetMomentum();
  fV0_birth.x = fV_v0.x - distNorm*GetPx();
  fV0_birth.y = fV_v0.y - distNorm*GetPy();
  fV0_birth.z = fV_v0.z - distNorm*GetPz();
}



//______________________________________________________________________
void V0::MakeTrack(vpPathMark_t& pathMark, Reve::Vector& vtx,  Reve::Vector& p,
		   Int_t charge, Float_t beta, TPolyLine3D& polyLine) {

  TrackRnrStyle& RS((fRnrStyle != 0) ? *fRnrStyle : TrackRnrStyle::fgDefStyle);

  Float_t px = p.x, py = p.y, pz = p.z;  

  MCVertex  mc_v0;
  mc_v0.x = vtx.x;
  mc_v0.y = vtx.y; 
  mc_v0.z = vtx.z; 
  mc_v0.t = 0;

  std::vector<MCVertex> track_points;
  Bool_t decay = kFALSE;

  if ((TMath::Abs(vtx.z) > RS.fMaxZ) || (vtx.x*vtx.x + vtx.y*vtx.y > RS.fMaxR*RS.fMaxR)) 
    goto make_polyline;
  
  if (TMath::Abs(RS.fMagField) > 1e-5) {

    // Charged particle in magnetic field

    Float_t a = RS.fgkB2C * RS.fMagField * charge;
   
    MCHelix helix(fRnrStyle, &mc_v0, TMath::C()*beta, &track_points, a); //m->cm
    helix.Init(TMath::Sqrt(px*px+py*py), pz);
   
    if(!pathMark.empty()){
      for(std::vector<Reve::PathMark*>::iterator i=pathMark.begin();
	  i!=pathMark.end(); ++i) {

	Reve::PathMark* pm = *i;
        
	if(RS.fFitDaughters &&  pm->type == Reve::PathMark::Daughter){
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto helix_bounds;

          //printf("%s fit daughter  \n", fName.Data()); 
	  helix.LoopToVertex(p.x, p.y, p.z, pm->V.x, pm->V.y, pm->V.z);
	  p.x -=  pm->P.x;
	  p.y -=  pm->P.y;
	  p.z -=  pm->P.z;
	}
	if(RS.fFitDecay &&  pm->type == Reve::PathMark::Decay){
	  
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto helix_bounds;
	  helix.LoopToVertex(p.x, p.y, p.z, pm->V.x, pm->V.y, pm->V.z);
          decay = true;
          break;
	}
      }
    }
  helix_bounds:
    //go to bounds
    if(!decay || RS.fFitDecay == kFALSE){
      helix.LoopToBounds(px,py,pz);
      // printf("%s loop to bounds  \n",fName.Data() );
    }

  } else {

    // Neutral particle or no field

    MCLine line(fRnrStyle, &mc_v0, TMath::C()*beta, &track_points);
   
    if(!pathMark.empty()) {
      for(std::vector<Reve::PathMark*>::iterator i=pathMark.begin();
	  i!=pathMark.end(); ++i) {
	Reve::PathMark* pm = *i;

	if(RS.fFitDaughters &&  pm->type == Reve::PathMark::Daughter){
          if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto line_bounds;
	  line.GotoVertex(pm->V.x, pm->V.y, pm->V.z);
	  p.x -=  pm->P.x;
	  p.y -=  pm->P.y;
	  p.z -=  pm->P.z;
	}

	if(RS.fFitDecay &&  pm->type == Reve::PathMark::Decay){
	  if(TMath::Abs(pm->V.z) > RS.fMaxZ 
	     || TMath::Sqrt(pm->V.x*pm->V.x + pm->V.y*pm->V.y) > RS.fMaxR )
	    goto line_bounds;
	  line.GotoVertex(pm->V.x, pm->V.y, pm->V.z);
          decay = true;
	  break;
	}
      }
    }

  line_bounds:
    if(!decay || RS.fFitDecay == kFALSE)
      line.GotoBounds(px,py,pz);

  }
make_polyline:
  Reset(&polyLine);
  for(std::vector<MCVertex>::iterator i=track_points.begin();
      i!=track_points.end(); ++i) {
    polyLine.SetNextPoint(i->x,i->y, i->z);
  }

}


//______________________________________________________________________
void V0::MakeV0path() {
  
  MCVertex  mc_v0;
  mc_v0.x = fV_v0.x;
  mc_v0.y = fV_v0.y; 
  mc_v0.z = fV_v0.z; 
  mc_v0.t = 0;

  std::vector<MCVertex> track_points;
  MCLine line(fRnrStyle, &mc_v0, TMath::C()*0.99, &track_points);

 line.GotoVertex(fV0_birth.x,fV0_birth.y,fV0_birth.z);

  Reset(&fPolyLineV0);
  for(std::vector<MCVertex>::iterator i=track_points.begin();
      i!=track_points.end(); ++i) {
    fPolyLineV0.SetNextPoint(i->x,i->y, i->z);
  }
}


//______________________________________________________________________
void V0::MakeV0()
{
  SetNextPoint(fV_v0.x, fV_v0.y, fV_v0.z);
  MakeTrack(fPathMarksNeg, fV_neg, fP_neg, -1, fBeta_neg, fPolyLineNeg);
  MakeTrack(fPathMarksPos, fV_pos, fP_pos,  1, fBeta_pos, fPolyLinePos);
  MakeV0path();
  //fN = fPolyLineNeg.GetN();
}


//______________________________________________________________________
Float_t const V0::GetAlphaArmenteros() {

  Float_t  posXv0 = fP_pos.x*GetPx() + fP_pos.y*GetPy() + fP_pos.z*GetPz();
  Float_t  negXv0 = fP_neg.x*GetPx() + fP_neg.y*GetPy() + fP_neg.z*GetPz();

  if (posXv0+negXv0 > 1.e-39)
    return (posXv0-negXv0)/(posXv0+negXv0);
  else return -999;
}

//______________________________________________________________________
  Float_t const V0::GetPtArmenteros() {

  Float_t  posXv0 = fP_pos.x*GetPx() + fP_pos.y*GetPy() + fP_pos.z*GetPz();
  Float_t  v0mom2  = GetP2();

  if (v0mom2 > 1.e-39)
    return  TMath::Sqrt( GetPosP2() - posXv0*posXv0/v0mom2 ) ;
  else return -999;
}





/***********************************************************************
*
*  V0List class
*
************************************************************************/

ClassImp(Reve::V0List)

//______________________________________________________________________
V0List::V0List() :
  RenderElementListBase(),
  fTitle(),
  fRnrStyle(0),
  fRnrDaughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fNegColor(0),
  fPosColor(0)
{
  for (Int_t i=0; i<fgkNcutVar; i++)
    fHist[i] = 0;
  for (Int_t i=0; i<fgkNcutVar2D; i++)
    fHist2D[i] = 0;
}

//______________________________________________________________________
V0List::V0List(TrackRnrStyle* rs) :
  RenderElementListBase(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fNegColor(0),
  fPosColor(0)
{
  Init();
}

//______________________________________________________________________
V0List::V0List(const Text_t* name, TrackRnrStyle* rs) :
  RenderElementListBase(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fNegColor(0),
  fPosColor(0)
{
  Init();
  SetName(name);
}

//______________________________________________________________________
void V0List::Init()
{

  if (fRnrStyle== 0) fRnrStyle = new TrackRnrStyle;
  SetMainColorPtr(&fRnrStyle->fColor);

  fMin[0]  = 0;   fMax[0] = 10; // pt
  fMin[1]  = 0;   fMax[1] = 10; // K0s mass
  fMin[2]  = 0;   fMax[2] = 10; // lambda mass
  fMin[3]  = 0;   fMax[3] = 10; // anti-lambda mass
  fMin[4]  = 0;   fMax[4] = 10; // daughter DCA
  fMin[5]  = 0.8; fMax[5] = 1;  // cos of pointing angle

  fMin[6]  =  0;   fMax[6] = 200;  // radius
  fMin[7]  = -2;   fMax[7] =   2;  // PseudoRapidity
  fMin[8]  =  0;   fMax[8] =  10;  // NegPt
  fMin[9]  = -2;   fMax[9] =   2;  // NegPseudoRapidity
  fMin[10] =  0;   fMax[10] = 10;  // PosPt
  fMin[11] = -2;   fMax[11] =  2;  // PosPseudoRapidity
  fMin[12] =  0;   fMax[12] =  1e5;  // index, would be good to be able to adjust it

  char *ch = "ptV0";
  fHist[0] = new TH1F(ch,ch, 100, fMin[0], fMax[0]);
  ch = "K0sMass";
  fHist[1] = new TH1F(ch,ch, 100, fMin[1], fMax[1]);
  ch = "LambdaMass";
  fHist[2] = new TH1F(ch,ch, 100, fMin[2], fMax[2]);
  ch = "AntiLambdaMass";
  fHist[3] = new TH1F(ch,ch, 100, fMin[3], fMax[3]);
  ch = "daughterDCA";
  fHist[4] = new TH1F(ch,ch, 100, fMin[4], fMax[4]);
  ch = "cosPointingAngle";
  fHist[5] = new TH1F(ch,ch, 100, fMin[5], fMax[5]);
  ch = "radius";
  fHist[6] = new TH1F(ch,ch, 100, fMin[6], fMax[6]);
  ch = "PseudoRapidity";
  fHist[7] = new TH1F(ch,ch, 100, fMin[7], fMax[7]);

  ch = "NegPt";
  fHist[8] = new TH1F(ch,ch, 100, fMin[8], fMax[8]);
  ch = "NegPseudoRapidity";
  fHist[9] = new TH1F(ch,ch, 100, fMin[9], fMax[9]);
  ch = "PosPt";
  fHist[10] = new TH1F(ch,ch, 100, fMin[10], fMax[10]);
  ch = "PosPseudoRapidity";
  fHist[11] = new TH1F(ch,ch, 100, fMin[11], fMax[11]);
  ch = "v0Index";
  fHist[12] = new TH1F(ch,ch, 100, fMin[12], fMax[12]);


  fMinX[0] = -1.2;
  fMaxX[0] = 1.2;
  fMinY[0] = 0;
  fMaxY[0] = 0.4;
  ch = "ArmenterosPodolansky";
  fHist2D[0] = new TH2F(ch,ch, 70, fMinX[0], fMaxX[0], 70,
			fMinY[0], fMaxY[0]);

  for (Int_t i=0; i<fgkNcutVar; i++) {
    fHist[i]->GetXaxis()->SetLabelSize(0.07);
    fHist[i]->GetYaxis()->SetLabelSize(0.07);
    fHist[i]->SetStats(0);
  }
  for (Int_t i=0; i<fgkNcutVar2D; i++) {
    fHist2D[i]->GetXaxis()->SetLabelSize(0.07);
    fHist2D[i]->GetYaxis()->SetLabelSize(0.07);
    fHist2D[i]->SetStats(0);
  }
}

//______________________________________________________________________
V0List::~V0List() {

  for (Int_t i=0; i<fgkNcutVar; i++)
    if (fHist[i]) delete fHist[i];
  for (Int_t i=0; i<fgkNcutVar2D; i++)
    if (fHist2D[i]) delete fHist2D[i];

}

//______________________________________________________________________
void V0List::Paint(Option_t* option) {
  if(fRnrElement) {

    if(fRnrV0vtx) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((V0*)(*i))->Paint(option);
	}
      }
    }

    if(fRnrDaughters) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((V0*)(*i))->PaintDaughters(option);
	}
      }
    }

    if(fRnrV0path) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((V0*)(*i))->PaintPath(option);
	}
      }
    }
  }
}


//______________________________________________________________________

void V0List::AddElement(RenderElement* el) {

  static const Exc_t eH("V0List::AddElement ");
  if (dynamic_cast<V0*>(el)  == 0)
    throw(eH + "new element not a V0.");
  RenderElementListBase::AddElement(el);
}



void V0List::SetRnrV0vtx(Bool_t rnr) {
  fRnrV0vtx = rnr;
  gReve->Redraw3D();
}

void V0List::SetRnrV0path(Bool_t rnr) {
  fRnrV0path = rnr;
  gReve->Redraw3D();
}

void V0List::SetRnrDaughters(Bool_t rnr) {
  fRnrDaughters = rnr;
  gReve->Redraw3D();
}


//______________________________________________________________________

void V0List::MakeV0s() {
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((V0*)(*i))->MakeV0();
  }
  gReve->Redraw3D();
}


void V0List::MakeMarkers() {
  gReve->Redraw3D();
}


//_________________________________________________________________________
void V0List::AdjustHist(Int_t iHist) {

  if (! fHist[iHist]) return;
  
  TString name = fHist[iHist]->GetName();
  Int_t nBin = fHist[iHist]->GetXaxis()->GetNbins();
  delete fHist[iHist];
  fHist[iHist] = new TH1F(name.Data(), name.Data(), nBin, GetMin(iHist),
			  GetMax(iHist));
  fHist[iHist]->GetXaxis()->SetLabelSize(0.07);
  fHist[iHist]->GetYaxis()->SetLabelSize(0.07);
  fHist[iHist]->SetStats(0);

}


//______________________________________________________________________
void V0List::UnFill(V0* v0) {

    Int_t bin = fHist[0]->GetXaxis()->FindBin(v0->GetPt());
    fHist[0]->SetBinContent( bin, fHist[0]->GetBinContent(bin)-1 );

    bin = fHist[1]->GetXaxis()->FindBin( v0->GetK0mass() );
    fHist[1]->SetBinContent( bin, fHist[1]->GetBinContent(bin)-1 );

    bin = fHist[2]->GetXaxis()->FindBin( v0->GetLamMass() );
    fHist[2]->SetBinContent( bin, fHist[2]->GetBinContent(bin)-1 );

    bin = fHist[3]->GetXaxis()->FindBin( v0->GetAntiLamMass() );
    fHist[3]->SetBinContent( bin, fHist[3]->GetBinContent(bin)-1 );

    bin = fHist[4]->GetXaxis()->FindBin( v0->GetDaughterDCA() );
    fHist[4]->SetBinContent( bin, fHist[4]->GetBinContent(bin)-1 );

    bin = fHist[5]->GetXaxis()->FindBin( v0->GetCosPointingAngle() );
    fHist[5]->SetBinContent( bin, fHist[5]->GetBinContent(bin)-1 );


    bin = fHist[6]->GetXaxis()->FindBin( v0->GetRadius() );// radius
    fHist[6]->SetBinContent( bin, fHist[6]->GetBinContent(bin)-1 );

    bin = fHist[7]->GetXaxis()->FindBin( v0->GetPseudoRapidity() ); // PseudoRapidity
    fHist[7]->SetBinContent( bin, fHist[7]->GetBinContent(bin)-1 );

    bin = fHist[8]->GetXaxis()->FindBin( v0->GetNegPt() );// NegPt
    fHist[8]->SetBinContent( bin, fHist[8]->GetBinContent(bin)-1 );

    bin = fHist[9]->GetXaxis()->FindBin( v0->GetNegPseudoRapidity() );// NegPseudoRapidity
    fHist[9]->SetBinContent( bin, fHist[9]->GetBinContent(bin)-1 );

    bin = fHist[10]->GetXaxis()->FindBin( v0->GetPosPt() ); // PosPt
    fHist[10]->SetBinContent( bin, fHist[10]->GetBinContent(bin)-1 );

    bin = fHist[11]->GetXaxis()->FindBin( v0->GetPosPseudoRapidity() ); // PosPseudoRapidity
    fHist[11]->SetBinContent( bin, fHist[11]->GetBinContent(bin)-1 );

    bin = fHist[12]->GetXaxis()->FindBin( v0->GetESDIndex() ); // ESD index
    fHist[12]->SetBinContent( bin, fHist[12]->GetBinContent(bin)-1 );

    //---
          bin  = fHist2D[0]->GetXaxis()->FindBin( v0->GetAlphaArmenteros() );
    Int_t binY = fHist2D[0]->GetYaxis()->FindBin( v0->GetPtArmenteros() );
    fHist2D[0]->SetBinContent( bin, binY, fHist2D[0]->GetBinContent(bin,binY)-1 );
}


//______________________________________________________________________
void V0List::Filter(V0* v0) {

  Float_t pt = v0->GetPt();
  if ((pt<fMin[0])||(pt>fMax[0])) return;

  Float_t k0sMass = v0->GetK0mass();
  if ( (k0sMass<fMin[1])||(k0sMass>fMax[1]) ) return;

  Float_t lamMass = v0->GetLamMass();
  if ( (lamMass<fMin[2])||(lamMass>fMax[2]) ) return;

  Float_t alamMass = v0->GetAntiLamMass();
  if ( (alamMass<fMin[3])||(alamMass>fMax[3]) ) return;

  Float_t daughtDCA = v0->GetDaughterDCA();
  if ( (daughtDCA<fMin[4])||(daughtDCA>fMax[4]) ) return;

  Float_t cosPointing = v0->GetCosPointingAngle();
  if ( (cosPointing<fMin[5])||(cosPointing>fMax[5]) ) return;


  Float_t radius = v0->GetRadius();
  if ( (radius<fMin[6])||(radius>fMax[6]) ) return;

  Float_t pseudoRapidity = v0->GetPseudoRapidity();
  if ( (pseudoRapidity<fMin[7])||(pseudoRapidity>fMax[7]) ) return;

  Float_t negPt = v0->GetNegPt();
  if ( (negPt<fMin[8])||(negPt>fMax[8]) ) return;

  Float_t negPseudoRapidity = v0->GetNegPseudoRapidity();
  if ( (negPseudoRapidity<fMin[9])||(negPseudoRapidity>fMax[9]) ) return;

  Float_t posPt = v0->GetPosPt();
  if ( (posPt<fMin[10])||(posPt>fMax[10]) ) return;

  Float_t posPseudoRapidity = v0->GetPosPseudoRapidity();
  if ( (posPseudoRapidity<fMin[11])||(posPseudoRapidity>fMax[11]) ) return;

   Float_t esdIndex = v0->GetESDIndex();
   if ( (esdIndex<fMin[12])||(esdIndex>fMax[12]) ) return;

   Float_t alphaArm = v0->GetAlphaArmenteros();
//   if ( (alphaArm<fMinX[0])||(alphaArm>fMaxX[0]) ) return;

   Float_t ptArm = v0->GetPtArmenteros();
//   if ( (ptArm<fMinY[0])||(ptArm>fMaxY[0]) ) return;

  v0->SetRnrElement(kTRUE);
  fHist[0]->Fill(pt);
  fHist[1]->Fill(k0sMass);
  fHist[2]->Fill(lamMass);
  fHist[3]->Fill(alamMass);
  fHist[4]->Fill(daughtDCA);
  fHist[5]->Fill(cosPointing);
  fHist[6]->Fill(radius);
  fHist[7]->Fill(pseudoRapidity);
  fHist[8]->Fill(negPt);
  fHist[9]->Fill(negPseudoRapidity);
  fHist[10]->Fill(posPt);
  fHist[11]->Fill(posPseudoRapidity);
  fHist[12]->Fill(esdIndex);
  fHist2D[0]->Fill(alphaArm, ptArm);
}

//______________________________________________________________________
void V0List::FilterAll() {

  for (Int_t i=0; i<fgkNcutVar; i++)
    fHist[i]->Reset();

  for (Int_t i=0; i<fgkNcutVar2D; i++)
    fHist2D[i]->Reset();
  
  V0* myV0;
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    Filter(myV0);
  }
}


//______________________________________________________________________
void V0List::GetV0IndexRange(Int_t &imin, Int_t &imax) {

  Int_t index;
  V0* myV0;
  lpRE_i i=fChildren.begin();
  myV0 = (V0*)(*i);
  index = myV0->GetESDIndex();
  imin = index;
  imax = index;

  for(; i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    index = myV0->GetESDIndex();
    if (index<imin) imin = index;
    if (index>imax) imax = index;
  }
}


//______________________________________________________________________
void V0List::PtFilter(Float_t min, Float_t max) {

  fMin[0] = min;
  fMax[0] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetPt();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}


//______________________________________________________________________
void V0List::K0sMFilter(Float_t min, Float_t max) {

  fMin[1] = min;
  fMax[1] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetK0mass();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::LamMFilter(Float_t min, Float_t max) {

  fMin[2] = min;
  fMax[2] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetLamMass();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}



//______________________________________________________________________
void V0List::ALamMFilter(Float_t min, Float_t max) {

  fMin[3] = min;
  fMax[3] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetAntiLamMass();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::CosPointingFilter(Float_t min, Float_t max) {

  fMin[5] = min;
  fMax[5] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetCosPointingAngle();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::DaughterDCAFilter(Float_t min, Float_t max) {

  fMin[4] = min;
  fMax[4] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetDaughterDCA();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::RadiusFilter(Float_t min, Float_t max) {

  fMin[6] = min;
  fMax[6] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetRadius();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::EtaFilter(Float_t min, Float_t max) {

  fMin[7] = min;
  fMax[7] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetPseudoRapidity();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::NegPtFilter(Float_t min, Float_t max) {

  fMin[8] = min;
  fMax[8] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetNegPt();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::NegEtaFilter(Float_t min, Float_t max) {

  fMin[9] = min;
  fMax[9] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetNegPseudoRapidity();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::PosPtFilter(Float_t min, Float_t max) {

  fMin[10] = min;
  fMax[10] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetPosPt();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::PosEtaFilter(Float_t min, Float_t max) {

  fMin[11] = min;
  fMax[11] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetPosPseudoRapidity();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

//______________________________________________________________________
void V0List::IndexFilter(Float_t min, Float_t max) {

  fMin[12] = min;
  fMax[12] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  V0* myV0;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myV0 = (V0*)(*i);
    val = myV0->GetESDIndex();
    wasSelected = myV0->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myV0);
	myV0->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myV0);
    }
  }
}

