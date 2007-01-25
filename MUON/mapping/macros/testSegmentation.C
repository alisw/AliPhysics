// $Id$
// $MpId: testSegmentation.C,v 1.1 2005/09/19 19:02:53 ivana Exp $

#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONVGeometryDESegmentation.h"
#include "AliMpPlaneType.h"
#include "AliMUONSt345SlatSegmentationV2.h"
//#include "AliMpSt345Reader.h"
#include "AliMpSegmentationManager.h"
#include "AliMpSlat.h"
#include "AliMpPCB.h"

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <sstream>

std::map<int,std::string> detElemIdToSlatTypeMap;

//_____________________________________________________________________________
bool equal(double a, double b, double epsilon=1e-5)
{
  double d = fabs(a-b)/b;
  if ( d < epsilon ) return true;
  return false;
}

//_____________________________________________________________________________
bool error(const char* msg)
{
  std::cout << msg << std::endl;
  return false;
}

//_____________________________________________________________________________
bool readDetElemIdToSlatType(const char* file = "/Users/aphecetc/Workspace/aliroot/MUON/mapping/data/station345/DetElemIdToSlatType.dat")
{
  std::ifstream in(file);
  if (!in.good()) return false;

  char line[80];

  while ( in.getline(line,80) )
    {
      if ( !isdigit(line[0]) ) continue;
      std::string sline(line);
      
      std::string::size_type pos = sline.find_first_of(' ');
      int detelemid = std::atoi(sline.substr(0,pos).c_str());
      detElemIdToSlatTypeMap[detelemid] = sline.substr(pos+1);
    }

  in.close();

  return true;
}

//_____________________________________________________________________________
const char* SlatType(int detelemid)
{
  if ( detElemIdToSlatTypeMap.empty() ) readDetElemIdToSlatType();

  std::map<int,std::string>::const_iterator it = 
    detElemIdToSlatTypeMap.find(detelemid);

  if ( it != detElemIdToSlatTypeMap.end() )
    {
      return it->second.c_str();
    }
  else
    {
      return 0;
    }
}

//_____________________________________________________________________________
bool testAllDetelem()
{
  // test we can read the map file (detelemid <-> slat type)
  bool ok = readDetElemIdToSlatType();
  if (!ok)
    {
      return false;
    }

  // now loops over all detelemid and check that:
  // 1. We can actually read the slat files, both bending and non-bending
  // 2. Bending and non-bending slats do have the same x-size

  std::map<int,std::string>::const_iterator it;
  std::set<std::string> slattypes;
  std::map<std::string,std::string> problems;

  for ( it = detElemIdToSlatTypeMap.begin(); it != detElemIdToSlatTypeMap.end();
	++it )
    {
      slattypes.insert(it->second);
      std::cout << it->first << " : " << it->second <<  " ";
      AliMpSlat* b = static_cast<AliMpSlat*>
	(AliMpSegmentationManager::Segmentation(it->first,AliMp::kBendingPlane));
      AliMpSlat* nb = static_cast<AliMpSlat*>
	(AliMpSegmentationManager::Segmentation(it->first,AliMp::kNonBendingPlane));
      std::cout << " B : " << b << " NB : " << nb;
      Double_t bx = 0.0;
      Double_t nbx = 0.0;
      std::ostringstream pb;

      if (b) 
	{
	  bx = b->DX()*2;
	}
      else
	{
	  std::cout << " Missing BENDING. ";
	}
     if (nb) 
	{
	  nbx = nb->DX()*2;
	}
     else
       {
	 std::cout << " Missing NON-BENDING.";
       }
     std::cout << std::endl;
     if (!b || !nb) 
       {
	 return false;
       }
     if ( bx != nbx )
       {
	 // Find which pcb(s) have a different size
	 if ( b->GetSize() != nb->GetSize() )
	   {
	     std::cout << "Not the same number of PCBs !" << std::endl;
	     return false;
	   }

	 pb << " DIFFERENT SIZES ! Bending= " << bx
	    << " Non-Bending= " << nbx
	    << " delta = " << bx - nbx;

	 for ( size_t i = 0; i < b->GetSize(); ++ i )
	   {
	     Double_t delta = b->GetPCB(i)->DX() - nb->GetPCB(i)->DX();
	     if ( fabs(delta*2) > 1E-3 )
	       {
		 pb << " DELTA(" << i << ")=" << delta*2;
	       }
	   }	 
	 pb << std::ends;
	 problems[it->second] = pb.str();
       }
    }

  //
  std::cout << "Number of detelem (per plane): " 
	    << detElemIdToSlatTypeMap.size()
	    << std::endl;
  std::cout << "Number of slat types : " << slattypes.size()
	    << std::endl;
  std::cout << "(Potential) Problems:" << std::endl;

  std::map<std::string,std::string>::const_iterator sit;

  for ( sit = problems.begin(); sit != problems.end(); ++sit ) 
    {
      std::cout << sit->first << " : " << sit->second << std::endl;
    }

  return ok;
}

//_____________________________________________________________________________
bool getSegmentation(int detElemId, 
		     AliMpPlaneType bendingornot,
		     AliMUONVGeometryDESegmentation*& oldSegm,
		     AliMUONVGeometryDESegmentation*& newSegm)
{
  AliMUON* muon = static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  if (!muon) 
    {
      gAlice->Init("${ALICE_ROOT}/MUON/Config.C");
      muon = static_cast<AliMUON*>(gAlice->GetModule("MUON"));
      if (!muon) return error("Cannot get MUON !");
    }

  int ichamber = detElemId/100 - 1;
  //  std::cout << "muon->Chamber(" << ichamber << ")" << std::endl;
  AliMUONChamber& chamber = muon->Chamber(ichamber);

  int icathode = 1;
  if ( bendingornot == AliMp::kNonBendingPlane ) icathode = 2;

  AliMUONGeometrySegmentation* gs = chamber.SegmentationModel2(icathode);
  if (!gs) return error(Form("Cannot get cathode %d",icathode));

  oldSegm = const_cast<AliMUONVGeometryDESegmentation*>(gs->GetDESegmentation(detElemId));
  if (!oldSegm) 
    {
      return error(Form("Cannot get segmentation for detElemId",detElemId));
    }

  newSegm =
    new AliMUONSt345SlatSegmentationV2(detElemId,bendingornot);

  return true;
}

//_____________________________________________________________________________
Int_t countPads(AliMUONVGeometryDESegmentation* seg)
{
  if (!seg) return 0;

  Int_t npads = 0;

  for ( Int_t ix = 1; ix <= seg->Npx(); ++ix )
    {
      for ( Int_t iy = 1; iy <= seg->Npy(); ++iy )
	{
	  if ( seg->HasPad(ix,iy) ) ++npads;
	}
    }
  return npads;
}

//_____________________________________________________________________________
bool testIC(Int_t d1, Int_t d2, AliMpPlaneType bendingornot)
{
  readDetElemIdToSlatType();
  AliMUONVGeometryDESegmentation* o;
  AliMUONVGeometryDESegmentation* n;

  Int_t ntested = 0;

  for ( Int_t d = d1; d <= d2; ++d )
    {
      if ( detElemIdToSlatTypeMap.find(d) == detElemIdToSlatTypeMap.end() ) continue;
      std::cout << d << std::endl;
      bool ok = getSegmentation(d,bendingornot,o,n);
      if (!ok) return false;

      Int_t nx = std::max(o->Npx(),n->Npx());
      Int_t ny = std::max(o->Npy(),n->Npy());
      
      for ( Int_t ix = 1; ix <= nx; ++ix )
	{
	  for ( Int_t iy = 1; iy <= ny; ++iy )
	    {
	      float xn,yn,zn;
	      float xo,yo,zo;
	      if ( !o->HasPad(ix,iy) || !n->HasPad(ix,iy) ) continue;
	      o->GetPadC(ix,iy,xo,yo,zo);
	      n->GetPadC(ix,iy,xn,yn,zn);
	      ++ntested;
	      if ( !equal(xn,xo) || !equal(yn,yo) )
		{
		  printf("%4d (%4d,%4d) -> OLD (%e,%e) NEW (%e,%e) DELTA (%e,%e)\n",d,ix,iy,xo,yo,xn,yn,xn-xo,yn-yo);
		}
	      // Now tries to go back to ix,iy from positions
	      Int_t nix,niy;
	      Int_t oix,oiy;
	      o->GetPadI(xo,yo,zo,oix,oiy);
	      n->GetPadI(xn,yn,zn,nix,niy);
	      std::string msg;
	      if ( ix != oix || iy != oiy )
		{
		  msg += "OLD";
		}
	      if ( ix != nix || iy != niy )
		{
		  msg += "NEW";
		}
	      if ( !msg.empty() )
		{
		  printf("Circular test failed for %s : (%4d,%4d) -> OLD (%e,%e) NEW (%e,%e) -> OLD (%4d,%4d) NEW (%4d,%4d)\n",msg.c_str(),ix,iy,xo,yo,xn,yn,oix,oiy,nix,niy);		  
		  return false;
		}
	    }
	}      
      delete n;
    }
  std::cout << "Number of tested pads  = " << ntested << std::endl;

  return true;
}

//_____________________________________________________________________________
void countPads(Int_t d1, Int_t d2)
{
  readDetElemIdToSlatType();

  Int_t onpads = 0;
  Int_t nnpads = 0;
  AliMUONVGeometryDESegmentation* o;
  AliMUONVGeometryDESegmentation* n;

  AliMpPlaneType pt[] = { AliMp::kNonBendingPlane, AliMp::kBendingPlane };

  for ( int ipt = 0; ipt < 2; ++ipt )
    {
      std::cout << ( (pt[ipt] == AliMp::kNonBendingPlane )?"NonBending"
		     : "Bending" )
		<< std::endl;
      Int_t p_nnpads = 0;
      Int_t p_onpads = 0;
      for ( Int_t d = d1; d <= d2; ++d )
	{
	  if ( detElemIdToSlatTypeMap.find(d) == detElemIdToSlatTypeMap.end() )
	    continue; // skip non-existing detElemId
	  bool ok = getSegmentation(d,pt[ipt],o,n);
	  if (!ok) return;
	  Int_t nn = countPads(n);
	  Int_t oo = countPads(o);
	  p_nnpads += nn;
	  p_onpads += oo;
	  printf("%4d OLD %5d NEW %5d Delta %6d OLDIXMAX %3d NEWIXMAX %3d\n",d,oo,nn,nn-oo,(o?o->Npx():0),n->Npx());
	  delete n;
	}

      std::cout << "OLD:" << p_onpads << std::endl
		<< "NEW:" << p_nnpads << std::endl;
      nnpads += p_nnpads;
      onpads += p_onpads;
    }

  std::cout << "-----" << std::endl
	    << "OLD:" << onpads << std::endl
	    << "NEW:" << nnpads << std::endl;
}
