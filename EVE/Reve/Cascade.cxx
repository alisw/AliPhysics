

/***********************************************************************
*  This code defines the reconstructed cascades visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/

#include "Track.h"
#include "Cascade.h"
#include "MCHelixLine.hi"

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

// Updates
#include <Reve/RGTopFrame.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include <vector>

using namespace Reve;



/***********************************************************************
*
*  Cascade class
*
************************************************************************/

const Float_t Cascade::fgkMassPion2 = 0.13956995*0.13956995;
const Float_t Cascade::fgkMassKaon2 = 0.493677*0.493677;
const Float_t Cascade::fgkMassProton2 = 0.93827231*0.93827231;
const Float_t Cascade::fgkMassLambda2 = 1.115683*1.115683;

ClassImp(Reve::Cascade)


Cascade::Cascade() :
  RenderElement(),
  TPolyMarker3D(1),
  fV_neg(),
  fP_neg(),
  fV_pos(),
  fP_pos(),
  fV_bach(),
  fP_bach(),
  fV_decay(),
  fV_birth(),
  fPathMarksNeg(),
  fPathMarksPos(),
  fPathMarksBach(),
  fRnrStyle(0),
  fPolyLineNeg(),
  fPolyLinePos(),
  fPolyLineBach(),
  fPolyLineV0(),
  fPolyLineCas(),
  fBeta_neg(0),
  fBeta_pos(0),
  fBeta_bach(0),
  fESDIndex(-1),
  fDCA_v0_Bach(999),
  fCasCosPointingAngle(999),
  fCasDecayLength(999)
{
  fPolyLinePos.SetLineColor(2);  // red
  fPolyLineNeg.SetLineColor(7);  // light blue
  fPolyLineBach.SetLineColor(4);  //  blue

  fMarkerStyle = 20;
  fMarkerColor = 5;
  fMarkerSize = 0.3;
}



Cascade::Cascade(TrackRnrStyle* rs) :
  RenderElement(),
  TPolyMarker3D(1),
  fV_neg(),
  fP_neg(),
  fV_pos(),
  fP_pos(),
  fV_bach(),
  fP_bach(),
  fV_decay(),
  fV_birth(),
  fPathMarksNeg(),
  fPathMarksPos(),
  fPathMarksBach(),
  fRnrStyle(rs),
  fPolyLineNeg(),
  fPolyLinePos(),
  fPolyLineBach(),
  fPolyLineV0(),
  fPolyLineCas(),
  fBeta_neg(0),
  fBeta_pos(0),
  fBeta_bach(0),
  fESDIndex(-1),
  fDCA_v0_Bach(999),
  fCasCosPointingAngle(999),
  fCasDecayLength(999)
{
  fMarkerColor = fRnrStyle->GetColor();
  fPolyLineV0.SetLineColor(fMarkerColor);
  fPolyLinePos.SetLineColor(2);  // red
  fPolyLineNeg.SetLineColor(7);  // light blue
  fPolyLineBach.SetLineColor(4);  //  blue

  fMainColorPtr = &fMarkerColor;
  fMarkerStyle = 20;
  fMarkerColor = 5;
  fMarkerSize = 0.3;
}


Cascade::~Cascade()
{
  for (vpPathMark_i i=fPathMarksNeg.begin(); i!=fPathMarksNeg.end(); ++i)
    delete *i;
  for (vpPathMark_i i=fPathMarksPos.begin(); i!=fPathMarksPos.end(); ++i)
    delete *i;
  for (vpPathMark_i i=fPathMarksBach.begin(); i!=fPathMarksBach.end(); ++i)
    delete *i;
}
void Cascade::Reset(TPolyLine3D* polyLine) {
  //polyLine->SetPolyLine(n_points);
  polyLine->SetPolyLine(0);
}


//______________________________________________________________________
void Cascade::SetDecayLength(Float_t primx, Float_t primy, Float_t primz) {


  Float_t dx = fV_decay.x-primx;
  Float_t dy = fV_decay.y-primy;
  Float_t dz = fV_decay.z-primz;

  fCasDecayLength = sqrt(dx*dx+dy*dy+dz*dz);
  // This is probably wrong but I can only do this for now
  Float_t distNorm = fCasDecayLength/GetMomentum();
  fV_birth.x = fV_decay.x - distNorm*GetPx();
  fV_birth.y = fV_decay.y - distNorm*GetPy();
  fV_birth.z = fV_decay.z - distNorm*GetPz();

  fV_bach.x = fV_decay.x;
  fV_bach.y = fV_decay.y;
  fV_bach.z = fV_decay.z;
}


//______________________________________________________________________
void Cascade::MakeTrack(vpPathMark_t& pathMark, Reve::Vector& vtx,  Reve::Vector& p,
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
void Cascade::MakeV0path() {
  
  MCVertex  mc_v0;
  mc_v0.x = (fV_neg.x+fV_pos.x)/2;
  mc_v0.y = (fV_neg.y+fV_pos.y)/2;
  mc_v0.z = (fV_neg.z+fV_pos.z)/2;
  mc_v0.t = 0;

  std::vector<MCVertex> track_points;
  MCLine line(fRnrStyle, &mc_v0, TMath::C()*0.99, &track_points);

 line.GotoVertex(fV_decay.x,fV_decay.y,fV_decay.z);

  Reset(&fPolyLineV0);
  for(std::vector<MCVertex>::iterator i=track_points.begin();
      i!=track_points.end(); ++i) {
    fPolyLineV0.SetNextPoint(i->x,i->y, i->z);
  }

}

//______________________________________________________________________
void Cascade::MakeCasPath() {
  
  MCVertex  mc_v0;
  mc_v0.x = fV_birth.x;
  mc_v0.y = fV_birth.y;
  mc_v0.z = fV_birth.z;
  mc_v0.t = 0;

  std::vector<MCVertex> track_points;
  MCLine line(fRnrStyle, &mc_v0, TMath::C()*0.99, &track_points);

 line.GotoVertex(fV_decay.x,fV_decay.y,fV_decay.z);

  Reset(&fPolyLineCas);
  for(std::vector<MCVertex>::iterator i=track_points.begin();
      i!=track_points.end(); ++i) {
    fPolyLineCas.SetNextPoint(i->x,i->y, i->z);
  }
}


//______________________________________________________________________
void Cascade::MakeCascade()
{
  SetNextPoint(fV_neg.x, fV_neg.y, fV_neg.z);
  SetNextPoint(fV_decay.x, fV_decay.y, fV_decay.z);

  MakeTrack(fPathMarksNeg, fV_neg, fP_neg, -1, fBeta_neg, fPolyLineNeg);
  MakeTrack(fPathMarksPos, fV_pos, fP_pos,  1, fBeta_pos, fPolyLinePos);
  if (fBeta_bach>0)
    MakeTrack(fPathMarksBach, fV_bach, fP_bach,  1, fBeta_bach, fPolyLineBach);
  else 
    MakeTrack(fPathMarksBach, fV_bach, fP_bach, -1, -fBeta_bach, fPolyLineBach);
  MakeV0path();
  MakeCasPath();
}




//______________________________________________________________________
Float_t const Cascade::GetCasAlphaArmenteros() {

  Float_t px = GetPx(), py = GetPy(), pz = GetPz();
  Float_t posXcas, negXcas;

  if (fBeta_bach>0) {
    posXcas = fP_bach.x*px + fP_bach.y*py + fP_bach.z*pz;
    negXcas = (fP_neg.x+fP_pos.x)*px + (fP_neg.y+fP_pos.y)*py + (fP_neg.z+fP_pos.z)*pz;
  } else {
    posXcas = (fP_neg.x+fP_pos.x)*px + (fP_neg.y+fP_pos.y)*py + (fP_neg.z+fP_pos.z)*pz;
    negXcas = fP_bach.x*px + fP_bach.y*py + fP_bach.z*pz;
  }

  if (posXcas + negXcas > 1.e-39)
    return (posXcas - negXcas)/(posXcas + negXcas);
  else return -999;
}


//______________________________________________________________________
Float_t const Cascade::GetCasPtArmenteros() {

  Float_t px = GetPx(), py = GetPy(), pz = GetPz();
  Float_t p2 = px*px + py*py + pz*pz;
  if (p2 < 1.e-39) return  -999;
  
  Float_t posXcas, posP2;
  
  if (fBeta_bach>0) {
    posXcas = fP_bach.x*px + fP_bach.y*py + fP_bach.z*pz;
    posP2 = GetBachP2();
  } else {
    posXcas = (fP_neg.x+fP_pos.x)*px + (fP_neg.y+fP_pos.y)*py + (fP_neg.z+fP_pos.z)*pz;
    posP2 = GetV0P2();
  }
  return sqrt( posP2 - posXcas*posXcas/p2 );
}



/***********************************************************************
*
*  CascadeList class
*
************************************************************************/

ClassImp(Reve::CascadeList)


//______________________________________________________________________
CascadeList::CascadeList(TrackRnrStyle* rs) :
  RenderElementListBase(),
  fTitle(),
  fRnrStyle(rs),
  fRnrBach(kTRUE),
  fRnrV0Daughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fRnrCasVtx(kTRUE),
  fRnrCasPath(kTRUE),
  fNegColor(0),
  fPosColor(0),
  fBachColor(0)
{
  Init();
}


//______________________________________________________________________
CascadeList::CascadeList(const Text_t* name, TrackRnrStyle* rs) :
  RenderElementListBase(),
  fTitle(),
  fRnrStyle(rs),
  fRnrBach(kTRUE),
  fRnrV0Daughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fRnrCasVtx(kTRUE),
  fRnrCasPath(kTRUE),
  fNegColor(0),
  fPosColor(0),
  fBachColor(0)
{
  Init();
  SetName(name);
}


//______________________________________________________________________
void CascadeList::Init()
{

  if (fRnrStyle== 0) fRnrStyle = new TrackRnrStyle;
  SetMainColorPtr(&fRnrStyle->fColor);


  fMin[0]  =  0;     fMax[0]  = 5; // Xi mass
  fMin[1]  =  0;     fMax[1]  = 5; // Omega mass
  fMin[2]  =  0;     fMax[2]  = 1e5; // Index
  fMin[3]  =  0.8;   fMax[3]  = 1; // cosPointingAngle
  fMin[4]  =  0;     fMax[4]  = 5; // bachV0DCA
  fMin[5]  =  0;     fMax[5]  = 100; // radius
  fMin[6]  =  0;     fMax[6]  = 10; // Pt
  fMin[7]  = -2;     fMax[7]  = 2; // PseudoRapidity
  fMin[8]  =  0;     fMax[8]  = 10; // negPt
  fMin[9]  = -2;     fMax[9]  = 2; // negEta
  fMin[10] =  0;    fMax[10]  = 10; // posPt
  fMin[11] = -2;    fMax[11]  = 2; // posEta
  fMin[12] =  0;    fMax[12]  = 10; // bachPt
  fMin[13] = -2;    fMax[13]  = 2; // backEta

  char *ch = "XiMass";
  fHist[0] = new TH1F(ch,ch, 100, fMin[0], fMax[0]);
  ch = "OmegaMass";
  fHist[1] = new TH1F(ch,ch, 100, fMin[1], fMax[1]);
  ch = "Index";
  fHist[2] = new TH1F(ch,ch, 100, fMin[2], fMax[2]);

  ch = "cosPointingAngle";
  fHist[3] = new TH1F(ch,ch, 100, fMin[3], fMax[3]);
  ch = "bachV0DCA";
  fHist[4] = new TH1F(ch,ch, 100, fMin[4], fMax[4]);
  ch = "radius";
  fHist[5] = new TH1F(ch,ch, 100, fMin[5], fMax[5]);
  ch = "Pt";
  fHist[6] = new TH1F(ch,ch, 100, fMin[6], fMax[6]);
  ch = "PseudoRapidity";
  fHist[7] = new TH1F(ch,ch, 100, fMin[7], fMax[7]);

  ch = "negPt";
  fHist[8] = new TH1F(ch,ch, 100, fMin[8], fMax[8]);
  ch = "negEta";
  fHist[9] = new TH1F(ch,ch, 100, fMin[9], fMax[9]);
  ch = "posPt";
  fHist[10] = new TH1F(ch,ch, 100, fMin[10], fMax[10]);
  ch = "posEta";
  fHist[11] = new TH1F(ch,ch, 100, fMin[11], fMax[11]);
  ch = "bachPt";
  fHist[12] = new TH1F(ch,ch, 100, fMin[12], fMax[12]);
  ch = "backEta";
  fHist[13] = new TH1F(ch,ch, 100, fMin[13], fMax[13]);

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
void CascadeList::Paint(Option_t* option) {
  if(fRnrElement) {

    if(fRnrBach) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((Cascade*)(*i))->PaintBachelor(option);
	}
      }
    }

    if(fRnrV0Daughters) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((Cascade*)(*i))->PaintV0Daughters(option);
	}
      }
    }

    if(fRnrV0path) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((Cascade*)(*i))->PaintV0Path(option);
	}
      }
    }

    if(fRnrCasVtx) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((Cascade*)(*i))->Paint(option);
	}
      }
    }

    if(fRnrCasPath) {
      for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrElement()) {
	  ((Cascade*)(*i))->PaintCasPath(option);
	}
      }
    }

  } // end if(fRnrElement)
}


//______________________________________________________________________

void CascadeList::AddElement(RenderElement* el)
{
  static const Exc_t eH("CascadeList::AddElement ");
  if (dynamic_cast<Cascade*>(el)  == 0)
    throw(eH + "new element not a Cascade.");
  RenderElementListBase::AddElement(el);
}



void CascadeList::SetRnrV0vtx(Bool_t rnr)
{
  fRnrV0vtx = rnr;
  gReve->Redraw3D();
}

void CascadeList::SetRnrV0path(Bool_t rnr)
{
  fRnrV0path = rnr;
  gReve->Redraw3D();
}

void CascadeList::SetRnrV0Daughters(Bool_t rnr)
{
  fRnrV0Daughters = rnr;
  gReve->Redraw3D();
}


void CascadeList::SetRnrCasPath(Bool_t rnr)
{
  fRnrCasPath = rnr;
  gReve->Redraw3D();
}

void CascadeList::SetRnrCasVtx(Bool_t rnr)
{
  fRnrCasVtx = rnr;
  gReve->Redraw3D();
}

void CascadeList::SetRnrBachelor(Bool_t rnr)
{
  fRnrBach = rnr;
  gReve->Redraw3D();
}


//______________________________________________________________________

void CascadeList::MakeCascades()
{
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((Cascade*)(*i))->MakeCascade();
  }
  gReve->Redraw3D();
}

//_________________________________________________________________________
void CascadeList::AdjustHist(Int_t iHist) {

  if ((iHist<0)||(iHist>=fgkNcutVar)) return;
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
void CascadeList::UnFill(Cascade* cas) {


    Int_t bin = fHist[0]->GetXaxis()->FindBin(cas->GetXiMass());
    fHist[0]->SetBinContent( bin, fHist[0]->GetBinContent(bin)-1 );

    bin = fHist[1]->GetXaxis()->FindBin( cas->GetOmegaMass() );
    fHist[1]->SetBinContent( bin, fHist[1]->GetBinContent(bin)-1 );


    bin = fHist[2]->GetXaxis()->FindBin( cas->GetESDIndex() );
    fHist[2]->SetBinContent( bin, fHist[2]->GetBinContent(bin)-1 );

    bin = fHist[3]->GetXaxis()->FindBin( cas->GetCasCosPointingAngle() );
    fHist[3]->SetBinContent( bin, fHist[3]->GetBinContent(bin)-1 );

    bin = fHist[4]->GetXaxis()->FindBin( cas->GetDCA_v0_Bach() );
    fHist[4]->SetBinContent( bin, fHist[4]->GetBinContent(bin)-1 );

    bin = fHist[5]->GetXaxis()->FindBin( cas->GetRadius() );
    fHist[5]->SetBinContent( bin, fHist[5]->GetBinContent(bin)-1 );

    bin = fHist[6]->GetXaxis()->FindBin( cas->GetPt() );
    fHist[6]->SetBinContent( bin, fHist[6]->GetBinContent(bin)-1 );

    bin = fHist[7]->GetXaxis()->FindBin( cas->GetPseudoRapidity() );
    fHist[7]->SetBinContent( bin, fHist[7]->GetBinContent(bin)-1 );
    //---

    bin = fHist[8]->GetXaxis()->FindBin( cas->GetNegPt() );
    fHist[8]->SetBinContent( bin, fHist[8]->GetBinContent(bin)-1 );

    bin = fHist[9]->GetXaxis()->FindBin( cas->GetNegPseudoRapidity() );
    fHist[9]->SetBinContent( bin, fHist[9]->GetBinContent(bin)-1 );

    bin = fHist[10]->GetXaxis()->FindBin( cas->GetPosPt() );
    fHist[10]->SetBinContent( bin, fHist[10]->GetBinContent(bin)-1 );

    bin = fHist[11]->GetXaxis()->FindBin( cas->GetPosPseudoRapidity() );
    fHist[11]->SetBinContent( bin, fHist[11]->GetBinContent(bin)-1 );

    bin = fHist[12]->GetXaxis()->FindBin( cas->GetBachPt() );
    fHist[12]->SetBinContent( bin, fHist[12]->GetBinContent(bin)-1 );

    bin = fHist[13]->GetXaxis()->FindBin( cas->GetBachPseudoRapidity() );
    fHist[13]->SetBinContent( bin, fHist[13]->GetBinContent(bin)-1 );

    //---
          bin  = fHist2D[0]->GetXaxis()->FindBin( cas->GetCasAlphaArmenteros() );
    Int_t binY = fHist2D[0]->GetYaxis()->FindBin( cas->GetCasPtArmenteros() );
    fHist2D[0]->SetBinContent( bin, binY, fHist2D[0]->GetBinContent(bin,binY)-1 );
}


//______________________________________________________________________
void CascadeList::Filter(Cascade* cas) {

  Float_t xiMass = cas->GetXiMass();
  if ((xiMass<fMin[0])||(xiMass>fMax[0])) return;

  Float_t omegaMass = cas->GetOmegaMass();
  if ( (omegaMass<fMin[1])||(omegaMass>fMax[1]) ) return;


  Float_t index = cas->GetESDIndex();
  if ( (index<fMin[2])||(index>fMax[2]) ) return;

  Float_t cosPointingAngle = cas->GetCasCosPointingAngle();
  if ( (cosPointingAngle<fMin[3])||(cosPointingAngle>fMax[3]) ) return;

  Float_t bachV0DCA = cas->GetDCA_v0_Bach();
  if ( (bachV0DCA<fMin[4])||(bachV0DCA>fMax[4]) ) return;

  Float_t radius = cas->GetRadius();
  if ( (radius<fMin[5])||(radius>fMax[5]) ) return;

  Float_t pt = cas->GetPt();
  if ( (pt<fMin[6])||(pt>fMax[6]) ) return;

  Float_t pseudoRapidity = cas->GetPseudoRapidity();
  if ( (pseudoRapidity<fMin[7])||(pseudoRapidity>fMax[7]) ) return;

  Float_t negPt = cas->GetNegPt();
  if ( (negPt<fMin[8])||(negPt>fMax[8]) ) return;

  Float_t negEta = cas->GetNegPseudoRapidity();
  if ( (negEta<fMin[9])||(negEta>fMax[9]) ) return;

  Float_t posPt = cas->GetPosPt();
  if ( (posPt<fMin[10])||(posPt>fMax[10]) ) return;

  Float_t posEta = cas->GetPosPseudoRapidity();
  if ( (posEta<fMin[11])||(posEta>fMax[11]) ) return;

  Float_t bachPt = cas->GetBachPt();
  if ( (bachPt<fMin[12])||(bachPt>fMax[12]) ) return;

  Float_t bachEta = cas->GetBachPseudoRapidity();
  if ( (bachEta<fMin[13])||(bachEta>fMax[13]) ) return;

  cas->SetRnrElement(kTRUE);
  fHist[0]->Fill(xiMass);
  fHist[1]->Fill(omegaMass);
  fHist[2]->Fill(index);
  fHist[3]->Fill(cosPointingAngle);
  fHist[4]->Fill(bachV0DCA);
  fHist[5]->Fill(radius);
  fHist[6]->Fill(pt);
  fHist[7]->Fill(pseudoRapidity);
  fHist[8]->Fill(negPt);
  fHist[9]->Fill(negEta);
  fHist[10]->Fill(posPt);
  fHist[11]->Fill(posEta);
  fHist[12]->Fill(bachPt);
  fHist[13]->Fill(bachEta);

  fHist2D[0]->Fill(cas->GetCasAlphaArmenteros(), cas->GetCasPtArmenteros() );
}

//______________________________________________________________________
void CascadeList::FilterAll() {

  for (Int_t i=0; i<fgkNcutVar; i++)
    fHist[i]->Reset();

  for (Int_t i=0; i<fgkNcutVar2D; i++)
    fHist2D[i]->Reset();
  
  Cascade* myCas;
  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    Filter(myCas);
  }
}


//______________________________________________________________________
void CascadeList::GetCasIndexRange(Int_t &imin, Int_t &imax) {

  Int_t index;
  Cascade* myCas;
  lpRE_i i = fChildren.begin();
  myCas = (Cascade*)(*i);
  index = myCas->GetESDIndex();
  imin = index;
  imax = index;

  for(; i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    index = myCas->GetESDIndex();
    if (index<imin) imin = index;
    if (index>imax) imax = index;
  }
}

//______________________________________________________________________
void CascadeList::XiMassFilter(Float_t min, Float_t max) {

  fMin[0] = min;
  fMax[0] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetXiMass();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}

//______________________________________________________________________
void CascadeList::OmegaMassFilter(Float_t min, Float_t max) {

  fMin[1] = min;
  fMax[1] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetOmegaMass();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}

//______________________________________________________________________
void CascadeList::IndexFilter(Float_t min, Float_t max) {

  fMin[2] = min;
  fMax[2] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetESDIndex();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}

//______________________________________________________________________
void CascadeList::CosPointingFilter(Float_t min, Float_t max) {

  fMin[3] = min;
  fMax[3] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetCasCosPointingAngle();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::BachV0DCAFilter(Float_t min, Float_t max) {

  fMin[4] = min;
  fMax[4] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetDCA_v0_Bach();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::RadiusFilter(Float_t min, Float_t max) {

  fMin[5] = min;
  fMax[5] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetRadius();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::PtFilter(Float_t min, Float_t max) {

  fMin[6] = min;
  fMax[6] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetPt();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::PseudoRapFilter(Float_t min, Float_t max) {

  fMin[7] = min;
  fMax[7] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetPseudoRapidity();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::NegPtFilter(Float_t min, Float_t max) {

  fMin[8] = min;
  fMax[8] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetNegPt();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::NegEtaFilter(Float_t min, Float_t max) {

  fMin[9] = min;
  fMax[9] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetNegPseudoRapidity();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::PosPtFilter(Float_t min, Float_t max) {

  fMin[10] = min;
  fMax[10] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetPosPt();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}

//______________________________________________________________________
void CascadeList::PosEtaFilter(Float_t min, Float_t max) {

  fMin[11] = min;
  fMax[11] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetPosPseudoRapidity();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}


//______________________________________________________________________
void CascadeList::BachPtFilter(Float_t min, Float_t max) {

  fMin[12] = min;
  fMax[12] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetBachPt();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}

//______________________________________________________________________
void CascadeList::BachEtaFilter(Float_t min, Float_t max) {

  fMin[13] = min;
  fMax[13] = max;

  Float_t val;
  Bool_t wasSelected;
  Bool_t isSelected;
  Cascade* myCas;

  for(lpRE_i i=fChildren.begin(); i!=fChildren.end(); ++i) {

    myCas = (Cascade*)(*i);
    val = myCas->GetBachPseudoRapidity();
    wasSelected = myCas->GetRnrElement();
    isSelected = ( (val>=min) && (val<=max) );

    if (wasSelected) {
      if (! isSelected) {
	UnFill(myCas);
	myCas->SetRnrElement(isSelected);
      }
    } else {
      if (isSelected) Filter(myCas);
    }
  }
}

