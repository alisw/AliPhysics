#include "NLTProjector.h"
#include "RGTopFrame.h"
#include "NLTPolygonSet.h"
#include "PODs.h"
#include "PointSet.h"
#include "Track.h"

#include "TBuffer3D.h"
#include "TColor.h"
#include "TPointSet3D.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

#include <list>

using namespace Reve;

namespace {
  struct Seg {
    Int_t v1;
    Int_t v2;

    Seg(Int_t i1=-1, Int_t i2=-1):v1(i1), v2(i2){};
  };
   
  typedef std::list<Seg>::iterator It_t;    
}

//______________________________________________________________________________
Vector*  NLTProjection::Project(Vector* origPnts, Int_t Npnts, Bool_t copy)
{
  Vector* pnts = 0; 
  if(copy) 
  {
    pnts = (Vector* )malloc(Npnts*sizeof(Vector));
    memcpy(pnts, origPnts, Npnts*sizeof(Vector));
  }
  else
  { 
    pnts =  origPnts;
  }
  return pnts;
}

//______________________________________________________________________________
Vector*  RhoZ::Project(Vector* origPnts, Int_t Npnts, Bool_t copy)
{
  Vector* pnts = NLTProjection::Project(origPnts, Npnts, copy);
  Float_t R, NR, y, z;
  for(Int_t i = 0; i<Npnts; i++) 
  {
    R = pnts[i].R();
    NR =R/(1+R*fDistortion);
    y = ( pnts[i].y > 0) ? NR : -NR;
    z = pnts[i].z;
    pnts[i].Set(z/(1+TMath::Abs(z*fDistortion)), y, 0);
  }
  return pnts;
}

//______________________________________________________________________________
Bool_t RhoZ::AcceptSegment(Vector& v1, Vector& v2) 
{
  return (v1.y*v2.y < 0) ? kFALSE : kTRUE;
}

//______________________________________________________________________________
Vector*  CircularFishEye::Project(Vector* origPnts, Int_t Npnts, Bool_t copy ) 
{
  Vector* pnts = NLTProjection::Project(origPnts, Npnts, copy);
  
  Float_t R, NR, phi;
  for(Int_t i = 0; i<Npnts; i++) {
    R = pnts[i].R();
    phi = pnts[i].Phi();
    NR = R/(1+R*fDistortion);
    pnts[i].Set(NR*TMath::Cos(phi), NR*TMath::Sin(phi), 0);
  }

  return pnts;
}

/**************************************************************************/
/**************************************************************************/
//______________________________________________________________________
// NLTProjector
//

ClassImp(NLTProjector)

NLTProjector::NLTProjector():
  TNamed("NLTProjector",""),
  fProjection(0),
  fEps(0.05),fIdxMap(0), 
  fNRPnts(0), fRPnts(0)
{
}

//______________________________________________________________________________
NLTProjector::~NLTProjector()
{
  CleanUp();

  if(fProjection) delete fProjection;
}

//______________________________________________________________________________
void NLTProjector::CleanUp()
{
  if(fIdxMap)
  {
    delete [] fIdxMap; 
    fIdxMap = 0;
  }
  if(fRPnts)
  {
    delete [] fRPnts; 
    fNRPnts  = 0;
    fRPnts = 0;
  }
}

/**************************************************************************/
void NLTProjector::SetProjection(NLTProjection::Type_e type, Float_t distort)
{
  static const Exc_t eH("NLTProjector::SetProjection ");
  switch (type)
  {
    case NLTProjection::RhoZ:
    {
      fProjection  = new RhoZ();
      break;
    }
    case NLTProjection::CFishEye:
    {
      fProjection  = new CircularFishEye();
      break;
    }
    default:
      throw(eH + "projection type not valid.");
      break;
  }
  fProjection->fDistortion = distort;
}

//______________________________________________________________________________
void NLTProjector::DumpBuffer(TBuffer3D* buff)
{
  Int_t* bpols = buff->fPols;
  
  for(UInt_t pi = 0; pi< buff->NbPols(); pi++) 
  {
    UInt_t Nseg = bpols[1]; 
    printf("%d polygon of %d has %d segments \n", pi,buff->NbPols(),Nseg);
    
    Int_t* seg =  &bpols[2];
    for(UInt_t a=0; a<Nseg; a++)
    {
      Int_t a1 = buff->fSegs[3*seg[a]+ 1];
      Int_t a2 = buff->fSegs[3*seg[a]+ 2];
      printf("(%d, %d) \n", a1, a2);
      printf("ORIG points :(%f, %f, %f)  (%f, %f, %f)\n", 
      	     buff->fPnts[3*a1],buff->fPnts[3*a1+1], buff->fPnts[3*a1+2],
      	     buff->fPnts[3*a2],buff->fPnts[3*a2+1], buff->fPnts[3*a2+2]);
    }
    printf("\n");
    bpols += (Nseg+2);
  }
}

//______________________________________________________________________________
void NLTProjector::CheckPoint(Int_t idx, Float_t* bbox)
{
  if(fRPnts[idx].x < bbox[0]) bbox[0] = fRPnts[idx].x;   
  if(fRPnts[idx].x > bbox[1]) bbox[1] = fRPnts[idx].x;

  if(fRPnts[idx].y < bbox[2]) bbox[2] = fRPnts[idx].y;   
  if(fRPnts[idx].y > bbox[3]) bbox[3] = fRPnts[idx].y;

  if(fRPnts[idx].z < bbox[4]) bbox[4] = fRPnts[idx].z;   
  if(fRPnts[idx].z > bbox[5]) bbox[5] = fRPnts[idx].z;
}

//______________________________________________________________________________
Bool_t NLTProjector::IsFirstIdxHead(Int_t s0, Int_t s1, TBuffer3D* buff)
{
  Int_t v0 = buff->fSegs[3*s0 + 1];
  Int_t v2 = buff->fSegs[3*s1 + 1];
  Int_t v3 = buff->fSegs[3*s1 + 2];
  if(v0 != v2 && v0 != v3 )
    return kTRUE;
  else 
    return kFALSE;
}

//______________________________________________________________________________
void  NLTProjector::ReducePoints(Vector* pnts, Int_t N)
{
  fIdxMap   = new Int_t[N];  
  Int_t* ra = new Int_t[N];  // list of reduced vertices
  fNRPnts = 0;
  
  for(UInt_t v = 0; v < (UInt_t)N; v++)
  {
    fIdxMap[v] = -1;
    for(Int_t k = 0; k<fNRPnts; k++) 
    {
      if(pnts[v].SquareDistance(pnts[ra[k]]) < fEps*fEps)
      {
	fIdxMap[v] = k; 
	break;
      }
    } 
    // have not find a point inside epsilon, add new point in scaled array
    if(fIdxMap[v] == -1)
    {
      fIdxMap[v] = fNRPnts;
      ra[fNRPnts] = v;
      fNRPnts++;
    }
    // printf("(%f, %f) vertex map %d -> %d \n", pnts[v*2], pnts[v*2 + 1], v, fIdxMap[v]);
  }
  
  // create an array of scaled points
  fRPnts = new Vector[fNRPnts];
  for(Int_t i = 0; i<fNRPnts; i++)
    fRPnts[i].Set(pnts[ra[i]].x,  pnts[ra[i]].y,  pnts[ra[i]].z);
  
  delete [] ra;  
  //  printf("reduced %d points of %d\n", fNRPnts, N);
}

//______________________________________________________________________________
void  NLTProjector::MakePolygonsFromBP(TBuffer3D* buff, std::list<NLTPolygon>& pols)
{
  //  printf("START NLTProjector::MakePolygonsFromBP\n");
  Int_t* bpols = buff->fPols;
  for(UInt_t pi = 0; pi< buff->NbPols(); pi++) 
  {
    std::list<Int_t>  pp; // points in current polygon 
    UInt_t Nseg = bpols[1]; 
    Int_t* seg =  &bpols[2];

    // start idx in the fist segment depends of second segment 
    Int_t  tail, head;
    Bool_t h = IsFirstIdxHead(seg[0], seg[1], buff);
    if(h) {
      head = fIdxMap[buff->fSegs[3*seg[0] + 1]];
      tail = fIdxMap[buff->fSegs[3*seg[0] + 2]];
    }
    else {
      head = fIdxMap[buff->fSegs[3*seg[0] + 2]];
      tail = fIdxMap[buff->fSegs[3*seg[0] + 1]];
    }
    Float_t bbox[] = {0., 0., 0., 0., 0., 0.};
    pp.push_back(head);
    CheckPoint(head, bbox);
    // printf("start idx head %d, tail %d\n", head, tail);
    
    std::list<Seg> segs;  
    for(UInt_t s=1; s< Nseg; s++)
      segs.push_back(Seg(buff->fSegs[3*seg[s] + 1],buff->fSegs[3*seg[s] + 2]));
    Bool_t accepted = kFALSE; 
    for(std::list<Seg>::iterator it = segs.begin(); it != segs.end(); it++ )
    { 
      Int_t mv1 = fIdxMap[(*it).v1];
      Int_t mv2 = fIdxMap[(*it).v2];         
      accepted = fProjection->AcceptSegment(fRPnts[mv1], fRPnts[mv2]);
      if(accepted == kFALSE)
      {
	pp.clear();
	break;
      }	  
      if(tail != pp.back()) 
      {
	pp.push_back(tail);
	CheckPoint(tail, bbox);
      }
      tail = (mv1 == tail) ? mv2 :mv1;
    }

    // render class loops indices last and first should not be equal
    if(accepted && pp.front() == pp.back())
      pp.pop_front();
    

    if( pp.size()>2 && (bbox[1]-bbox[0])>fEps && (bbox[3]-bbox[2])>fEps)
    {
      Bool_t duplicate = kFALSE;
      for (std::list<NLTPolygon>::iterator poi = pols.begin(); poi!= pols.end(); poi++)
      {
        NLTPolygon P = *poi;
        if(pp.size() != (UInt_t)P.fNPnts) 
	  continue;      
	Int_t matched = 0;
	for (std::list<Int_t>::iterator u = pp.begin(); u!= pp.end(); u++) {
          if (*u == P.fPnts[*u]) matched++;
          else continue;
        }
        if (matched == P.fNPnts) {
          duplicate = kTRUE;
          break;
        }
      }

      if(duplicate == kFALSE) 
      {
	Int_t* pv = new Int_t[pp.size()];
	Int_t count=0;
	// printf("%d NLTPolygon points %d \n", pols.size(), pp.size());
        // printf("bbox (%f, %f) (%f, %f)\n", bbox[1], bbox[0], bbox[3], bbox[2]);
	for( std::list<Int_t>::iterator u = pp.begin(); u!= pp.end(); u++){
	  pv[count] = *u;
	  count++;
	}
	pols.push_back(NLTPolygon(pp.size(), pv));
      }
    } // end of NLTPolygon
    bpols += (Nseg+2);
  }
}// MakePolygonsFromBP

//______________________________________________________________________________
void  NLTProjector::MakePolygonsFromBS(TBuffer3D* buff, std::list<NLTPolygon>& pols)
{
  // create your own list of segments according to reduced and projected points
  std::list<Seg> segs;  
  std::list<Seg>::iterator it;
  for(UInt_t s=0; s< buff->NbSegs(); s++)
  {
    Bool_t duplicate = kFALSE;
    Int_t vo1, vo2;   // idx from buff segment
    Int_t vor1, vor2; // mapped idx 
    vo1 =  buff->fSegs[3*s + 1];
    vo2 =  buff->fSegs[3*s + 2]; //... skip color info
    vor1 = fIdxMap[vo1];
    vor2 = fIdxMap[vo2];
    if(vor1 == vor2) continue;
    // check duplicate
    for(it = segs.begin(); it != segs.end(); it++ ){
      Int_t vv1 = (*it).v1;
      Int_t vv2 = (*it).v2;
      if((vv1 == vor1 && vv2 == vor2 )||(vv1 == vor2 && vv2 == vor1 )){
	duplicate = kTRUE;
	continue;
      }
    }
    if(duplicate == kFALSE && fProjection->AcceptSegment(fRPnts[vor1], fRPnts[vor2]))
    {
      segs.push_back(Seg(vor1, vor2));
    }
  }

  while(segs.empty() == kFALSE)
  {
    // printf("Start building polygon %d from %d segments in POOL \n", pols.size(), segs.size());
    std::list<Int_t> pp; // points in current polygon
    pp.push_back(segs.front().v1);
    Int_t tail = segs.front().v2;
    segs.pop_front();
    Bool_t match = kTRUE;
    while(match && segs.empty() == kFALSE)
    {
      // printf("second loop search tail %d \n",tail);
      for(It_t k=segs.begin(); k!=segs.end(); ++k){
	Int_t cv1 = (*k).v1;
	Int_t cv2 = (*k).v2;
	if( cv1 == tail || cv2 == tail){
	  // printf("found point %d  in %d,%d  \n", tail, cv1, cv2);
	  pp.push_back(tail);
	  if(cv1 == tail) 
	    tail = cv2;
	  else 
	    tail = cv1;

	  It_t to_erase = k--;
	  segs.erase(to_erase);
	  match = kTRUE;
	  break;
	}
	else 
	{
	  match = kFALSE;
	}
      } // end for loop in the segment pool
      if(tail ==  pp.front())
	break;
    };

    if(pp.size() > 2){      
      // printf("CLOSING polygon %d with %d points \n", pols.size(),pp.size());
      Int_t* pv = new Int_t[pp.size()];
      Int_t count=0;
      for( std::list<Int_t>::iterator u = pp.begin(); u!= pp.end(); u++){
	pv[count] = *u;
	count++;
      }
      pols.push_back(NLTPolygon(pp.size(), pv));
    }
  } // while segment pool empty
}//MakePolygonsFromBS

 //______________________________________________________________________________
NLTPolygonSet*  NLTProjector::ProjectGeoShape(TBuffer3D* buff, Int_t useBuffPols)
{ 
  // DumpBuffer(buff);  
  // rewrite from Double_t to Vector array
  Vector*  pnts  = new Vector[buff->NbPnts()];
  for(Int_t i = 0; i<(Int_t)buff->NbPnts(); i++) 
    pnts[i].Set(buff->fPnts[3*i],buff->fPnts[3*i+1], buff->fPnts[3*i+2]);

  fProjection->Project(pnts, buff->NbPnts(), kFALSE);
  ReducePoints(pnts, buff->NbPnts());

  // build polygons
  std::list<NLTPolygon>   pols;
  if(useBuffPols == -1) 
  { 
    MakePolygonsFromBP(buff, pols); 
    if(pols.empty())
    {
      // printf("useBuffPols == -1 call MakePolygonsFromBS \n");
      MakePolygonsFromBS(buff, pols);
    }
  }
  if(useBuffPols == 1)
  {
    MakePolygonsFromBP(buff, pols); 
  }
  if(useBuffPols == 0) 
  {
    MakePolygonsFromBS(buff, pols); 
  }
 
  NLTPolygonSet* ps = new NLTPolygonSet();
  if(pols.empty() == kFALSE) 
  {
    // points
    ps->SetPoints(fRPnts,fNRPnts);
    // polygons
    NLTPolygon* nltp = new NLTPolygon[pols.size()];
    Int_t pc = 0;
    for( std::list<NLTPolygon>::iterator pi = pols.begin(); pi!= pols.end(); pi++)
    {
      nltp[pc].fNPnts = (*pi).fNPnts;
      nltp[pc].fPnts = (*pi).fPnts;
      pc++;
    }
    ps->SetPolygons(nltp, pols.size());
    // ps->Dump();
  }
  fRPnts = 0;
  CleanUp();
  return ps;
}

//______________________________________________________________________________
void NLTProjector::ProjectPointSet(PointSet* orig)
{
  // keep original rewrite points into Vector array
  Float_t* op = orig->GetP();
  Vector* vp = new Vector[orig->GetN()];
  for( Int_t k=0; k<orig->GetN(); k++)
    vp[k].Set(op[3*k], op[3*k+1], op[3*k+2]);
  Vector* pnts = fProjection->Project(vp, orig->GetN(), kFALSE);
  ReducePoints(pnts, orig->GetN());
  orig->Reset(fNRPnts);
  for(Int_t i = 0; i< fNRPnts; i++){
    orig->SetPoint(i, fRPnts[i].x, fRPnts[i].y, fRPnts[3].z);
  }

  CleanUp();
}
  
//______________________________________________________________________________
void  NLTProjector::ProjectTrackList(TrackList* orig_tl)
{
  // fill a track list from original track list 
  Reve::RenderElement::List_i i = orig_tl->BeginChildren();
  Int_t Nc =  orig_tl->GetNChildren();
  for(Int_t c=0; c<Nc; c++)
  {
    // printf("project track with %d points \n", orig->GetN());
    Track* orig = dynamic_cast<Track*>(*i);
    Int_t start = 0;
    // keep original copy points to Vector array
    Float_t* op = orig->GetP();
    Vector* pnts = new Vector[orig->GetN()];
    for( Int_t k=0; k<orig->GetN(); k++)
      pnts[k].Set(op[3*k], op[3*k+1], op[3*k+2]);

    for( Int_t k=0; k<orig->GetN()-1; k++)
    {
      if(fProjection->AcceptSegment(pnts[k], pnts[(k+1)]) == kFALSE){
        //printf("split track at %d\n", k);
	Vector* tp = fProjection->Project(&pnts[start], k-start, kFALSE);
        ReducePoints(tp, k-start);
	Track* nt = new Track();  
        nt->Reset(fNRPnts);
	nt->SetMainColor(orig->GetMainColor());
        for(Int_t j = 0; j<fNRPnts; j++)
          nt->SetPoint(j, fRPnts[j].x, fRPnts[j].y, fRPnts[j].z);
        orig_tl->AddElement(nt);
	CleanUp();
        start = k+1;
      }
    }
    if(start != orig->GetN()-1) {
      // printf("finish last track %d %d\n", start,orig->GetN()-1);
      Track* nt = new Track();
      Vector* tp = fProjection->Project(&pnts[start], orig->GetN()-start, kFALSE);
      ReducePoints(tp, orig->GetN()-start);
      nt->Reset(fNRPnts);
      nt->SetMainColor(orig->GetMainColor());
      for(Int_t j = 0; j<fNRPnts; j++)
	nt->SetPoint(j, fRPnts[j].x, fRPnts[j].y, fRPnts[j].z);
      orig_tl->AddElement(nt);
      CleanUp();
    }
    // destroy unprojected
    ++i;
  } // end loop tracks

  // remove original tracks
  for(Int_t c=0; c<Nc; c++)
  {
    i = orig_tl->BeginChildren();
    (*i)->Destroy();
  }
}
