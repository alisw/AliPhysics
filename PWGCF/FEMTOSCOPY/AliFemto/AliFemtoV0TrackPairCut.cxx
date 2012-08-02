/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityPairCut - a pair cut which checks for some pair     //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityPairCut.cxx 50722 2011-07-21 15:18:38Z akisiel $
 *
 * Author: Adam Kisiel, Ohio State, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/

#include "AliFemtoV0TrackPairCut.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoV0TrackPairCut)
#endif

//__________________
AliFemtoV0TrackPairCut::AliFemtoV0TrackPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareQualityMax(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(0),
  fTrackTPCOnly(0),
  fDataType(kAOD),
  fDTPCMin(0),
  fDTPCExitMin(0),
  fKstarCut(0),
  fFirstParticleType(kLambda), 
  fSecondParticleType(kProton),
  fMinAvgSepTrackPos(0),
  fMinAvgSepTrackNeg(0)
{
  // Default constructor
  // Nothing to do
}
//__________________
AliFemtoV0TrackPairCut::~AliFemtoV0TrackPairCut(){
  /* no-op */
}
AliFemtoV0TrackPairCut& AliFemtoV0TrackPairCut::operator=(const AliFemtoV0TrackPairCut& cut)
{
  if (this != &cut) {
   
    AliFemtoPairCut::operator=(cut);
    fNPairsPassed = 0;
    fNPairsFailed =0;
    fV0Max = 1.0;
    fShareQualityMax = 1.0;
    fShareFractionMax = 1.0;
    fRemoveSameLabel = 0;
    fTrackTPCOnly = 0;
    fDataType = kAOD;
    fDTPCMin = 0;
    fDTPCExitMin = 0;
    fMinAvgSepTrackPos = 0;
    fMinAvgSepTrackNeg = 0;
  }

  return *this;
}
//__________________
bool AliFemtoV0TrackPairCut::Pass(const AliFemtoPair* pair){
  // Check for pairs that are possibly shared/double reconstruction

  bool temp = true;
  //Track1 - V0
  //Track2 - track

  if(!(pair->Track1()->V0() && pair->Track2()->Track()))
    {
      return false;
    }
  if(fTrackTPCOnly)
    {
      if( (-(pair->Track1()->V0()->IdNeg()+1)) ==pair->Track2()->TrackId() || (-(pair->Track1()->V0()->IdPos()+1)) ==pair->Track2()->TrackId())
	{
	  return false;
	}
    }
  else
    {
      if(pair->Track1()->V0()->IdNeg()==pair->Track2()->TrackId() || pair->Track1()->V0()->IdPos()==pair->Track2()->TrackId())
	{
	  return false;
	}
    }
  
  bool tempTPCEntrancePos = true;
  bool tempTPCEntranceNeg = true;
  bool tempTPCExitPos = true;
  bool tempTPCExitNeg = true;
  if(fDataType==kESD || fDataType==kAOD)
    {
      double distx = pair->Track1()->V0()->NominalTpcEntrancePointPos().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
      double disty = pair->Track1()->V0()->NominalTpcEntrancePointPos().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
      double distz = pair->Track1()->V0()->NominalTpcEntrancePointPos().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
      double distPos = sqrt(distx*distx + disty*disty + distz*distz);

      distx = pair->Track1()->V0()->NominalTpcEntrancePointNeg().x() - pair->Track2()->Track()->NominalTpcEntrancePoint().x();
      disty = pair->Track1()->V0()->NominalTpcEntrancePointNeg().y() - pair->Track2()->Track()->NominalTpcEntrancePoint().y();
      distz = pair->Track1()->V0()->NominalTpcEntrancePointNeg().z() - pair->Track2()->Track()->NominalTpcEntrancePoint().z();
      double distNeg = sqrt(distx*distx + disty*disty + distz*distz);

      double distExitX = pair->Track1()->V0()->NominalTpcExitPointPos().x() - pair->Track2()->Track()->NominalTpcExitPoint().x();
      double distExitY = pair->Track1()->V0()->NominalTpcExitPointPos().y() - pair->Track2()->Track()->NominalTpcExitPoint().y();
      double distExitZ = pair->Track1()->V0()->NominalTpcExitPointPos().z() - pair->Track2()->Track()->NominalTpcExitPoint().z();
      double distExitPos = sqrt(distExitX*distExitX + distExitY*distExitY + distExitZ*distExitZ);

      distExitX = pair->Track1()->V0()->NominalTpcExitPointNeg().x() - pair->Track2()->Track()->NominalTpcExitPoint().x();
      distExitY = pair->Track1()->V0()->NominalTpcExitPointNeg().y() - pair->Track2()->Track()->NominalTpcExitPoint().y();
      distExitZ = pair->Track1()->V0()->NominalTpcExitPointNeg().z() - pair->Track2()->Track()->NominalTpcExitPoint().z();
      double distExitNeg = sqrt(distExitX*distExitX + distExitY*distExitY + distExitZ*distExitZ);

      tempTPCEntrancePos = distPos > fDTPCMin;
      tempTPCEntranceNeg = distNeg > fDTPCMin;
      tempTPCExitPos = distExitPos > fDTPCExitMin;
      tempTPCExitNeg = distExitNeg > fDTPCExitMin;
    }
 
  if(!tempTPCEntrancePos || !tempTPCEntranceNeg || !tempTPCExitPos || !tempTPCExitNeg) return false;

  //reject merged trakcs in TPC
  

  //temp = dist > fDTPCMin;
  //koniec kopii
  //if(!temp)
  //return false;

  //kopia z AliFemtoShareQualityPairCut.cxx
  Int_t nh = 0;
  Int_t an = 0;
  Int_t ns = 0;

 
  
  if ((fShareFractionMax < 1.0) && ( fShareQualityMax < 1.0)) {
    for (unsigned int imap=0; imap<pair->Track1()->V0()->TPCclustersPos().GetNbits(); imap++) {
      // If both have clusters in the same row
      if (pair->Track1()->V0()->TPCclustersPos().TestBitNumber(imap) && 
	  pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// Do they share it ?
	if (pair->Track1()->V0()->TPCsharingPos().TestBitNumber(imap) &&
	    pair->Track2()->Track()->TPCsharing().TestBitNumber(imap))
	  {
	    //	  cout << "A shared cluster !!!" << endl;
	    //	cout << "imap idx1 idx2 " 
	    //	     << imap << " "
	    //	     << tP1idx[imap] << " " << tP2idx[imap] << endl;
	    an++;
	    nh+=2;
	    ns+=2;
	  }
	
	// Different hits on the same padrow
	else {
	  an--;
	  nh+=2;
	}
      }
      else if (pair->Track1()->V0()->TPCclustersPos().TestBitNumber(imap) ||
	       pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// One track has a hit, the other does not
	an++;
	nh++;
      }
    }
    
    Float_t hsmval = 0.0;
    Float_t hsfval = 0.0;
    
    if (nh >0) {
      hsmval = an*1.0/nh;
      hsfval = ns*1.0/nh;
    }
    //  if (hsmval > -0.4) {
    //   cout << "Pair quality: " << hsmval << " " << an << " " << nh << " " 
    //        << (pair->Track1()->Track()) << " " 
    //        << (pair->Track2()->Track()) << endl;
    //   cout << "Bits: " << pair->Track1()->Track()->TPCclusters().GetNbits() << endl;
    //  }
    //   if (hsfval > 0.0) {
    //     cout << "Pair sharity: " << hsfval << " " << ns << " " << nh << "    " << hsmval << " " << an << " " << nh << endl;
    //   }
    
    temp = (hsmval < fShareQualityMax) && (hsfval < fShareFractionMax);
    if(!temp) return false;

    nh = 0;
    an = 0;
    ns = 0;

    for (unsigned int imap=0; imap<pair->Track1()->V0()->TPCclustersNeg().GetNbits(); imap++) {
      // If both have clusters in the same row
      if (pair->Track1()->V0()->TPCclustersNeg().TestBitNumber(imap) && 
	  pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// Do they share it ?
	if (pair->Track1()->V0()->TPCsharingNeg().TestBitNumber(imap) &&
	    pair->Track2()->Track()->TPCsharing().TestBitNumber(imap))
	  {
	    //	  cout << "A shared cluster !!!" << endl;
	    //	cout << "imap idx1 idx2 " 
	    //	     << imap << " "
	    //	     << tP1idx[imap] << " " << tP2idx[imap] << endl;
	    an++;
	    nh+=2;
	    ns+=2;
	  }
	
	// Different hits on the same padrow
	else {
	  an--;
	  nh+=2;
	}
      }
      else if (pair->Track1()->V0()->TPCclustersNeg().TestBitNumber(imap) ||
	       pair->Track2()->Track()->TPCclusters().TestBitNumber(imap)) {
	// One track has a hit, the other does not
	an++;
	nh++;
      }
    }
    
    hsmval = 0.0;
    hsfval = 0.0;
    
    if (nh >0) {
      hsmval = an*1.0/nh;
      hsfval = ns*1.0/nh;
    }
    //  if (hsmval > -0.4) {
    //   cout << "Pair quality: " << hsmval << " " << an << " " << nh << " " 
    //        << (pair->Track1()->Track()) << " " 
    //        << (pair->Track2()->Track()) << endl;
    //   cout << "Bits: " << pair->Track1()->Track()->TPCclusters().GetNbits() << endl;
    //  }
    //   if (hsfval > 0.0) {
    //     cout << "Pair sharity: " << hsfval << " " << ns << " " << nh << "    " << hsmval << " " << an << " " << nh << endl;
    //   }
    
    temp = (hsmval < fShareQualityMax) && (hsfval < fShareFractionMax);


  }
  else
    temp = true;
  //koniec kopii


  //avg sep pair cut
  double avgSep=0;
  AliFemtoThreeVector first, second, tmp;
  for(int i=0; i<8 ;i++)
    {
      tmp = pair->Track1()->V0()->NominalTpcPointPos(i);
      //cout<<"X pos: "<<tmp.x()<<endl;
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());
      
      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z()); 

      avgSep += TMath::Sqrt(((double)first.x()-(double)second.x())*((double)first.x()-(double)second.x())+((double)first.y()-(double)second.y())*((double)first.y()-second.y())+((double)first.z()-(double)second.z())*((double)first.z()-(double)second.z()));
    }
  avgSep /= 8;

  if(avgSep<fMinAvgSepTrackPos) return false;

  avgSep = 0;
  
  for(int i=0; i<8 ;i++)
    {
      tmp = pair->Track1()->V0()->NominalTpcPointNeg(i);
      //cout<<"X pos: "<<tmp.x()<<endl;
      first.SetX((double)(tmp.x()));
      first.SetY((double)tmp.y());
      first.SetZ((double)tmp.z());
      
      tmp = pair->Track2()->Track()->NominalTpcPoint(i);
      second.SetX((double)tmp.x());
      second.SetY((double)tmp.y());
      second.SetZ((double)tmp.z()); 

      avgSep += TMath::Sqrt(((double)first.x()-(double)second.x())*((double)first.x()-(double)second.x())+((double)first.y()-(double)second.y())*((double)first.y()-second.y())+((double)first.z()-(double)second.z())*((double)first.z()-(double)second.z()));
    }
  avgSep /= 8;

  if(avgSep<fMinAvgSepTrackNeg) return false;




  //Qinv cut (we are trying to get rid of antiresidual correlation between primary protons)
  if(fKstarCut > 0)
    {

      //2 daughters of first V0
      double mom1PosX = pair->Track1()->V0()->MomPosX();
      double mom1PosY = pair->Track1()->V0()->MomPosY();
      double mom1PosZ = pair->Track1()->V0()->MomPosZ();
      double mom1NegX = pair->Track1()->V0()->MomNegX();
      double mom1NegY = pair->Track1()->V0()->MomNegY();
      double mom1NegZ = pair->Track1()->V0()->MomNegZ();


      //double PionMass = 0.13956995;
      //double KaonMass = 0.493677;
      double ProtonMass = 0.938272;
      //double LambdaMass = 1.115683;

      AliFemtoLorentzVector fFourMomentum1; // Particle momentum
      AliFemtoLorentzVector fFourMomentum2; // Particle momentum

      AliFemtoThreeVector temp1;
      double ener1=0;
      if(fFirstParticleType == 0) //lambda
	{
	  if(fSecondParticleType == 2) //proton
	    {
	      //temp1 = ::sqrt(mom1PosX*mom1PosX+mom1PosY*mom1PosY+mom1PosZ*mom1PosZ);
	      temp1.SetX(mom1PosX);
	      temp1.SetY(mom1PosY);
	      temp1.SetZ(mom1PosZ);
	      ener1 = ::sqrt(temp1.Mag2()+ProtonMass*ProtonMass);	  
	    }
	}
      else if(fFirstParticleType == 1) //antilambda
	{
	  if(fSecondParticleType == 3) //antiproton
	    {
	      //temp1 = ::sqrt(mom1NegX*mom1NegX+mom1NegY*mom1NegY+mom1NegZ*mom1NegZ);
	      temp1.SetX(mom1NegX);
	      temp1.SetY(mom1NegY);
	      temp1.SetZ(mom1NegZ);
	      ener1 = ::sqrt(temp1.Mag2()+ProtonMass*ProtonMass);
	    }
	}
      fFourMomentum1.SetVect(temp1);
      fFourMomentum1.SetE(ener1);

      //AliFemtoLorentzVector fFourMomentum2; // Particle momentum
      AliFemtoThreeVector temp2;
      double ener2=0;

      //2 daughters of second V0
      temp2 = pair->Track2()->Track()->P();
 
 
      if(fSecondParticleType == 2 || fSecondParticleType == 3) //proton
	{
	    ener2 = ::sqrt(temp2.Mag2()+ProtonMass*ProtonMass);
	}

      fFourMomentum2.SetVect(temp2);
      fFourMomentum2.SetE(ener2);

      //QInv calculation
      AliFemtoLorentzVector tDiff = (fFourMomentum1-fFourMomentum2);

      double tQinv = fabs(-1.*tDiff.m());   // note - qInv() will be negative for identical pairs...
   
      //cout<<"tQinv/2: "<<tQinv/2<<endl;
      if(tQinv/2 < fKstarCut)
	{

	  temp = false;
	}
    }
  return temp;
}
//__________________
AliFemtoString AliFemtoV0TrackPairCut::Report(){
  // Prepare the report from the execution
  string stemp = "AliFemtoV0 Pair Cut - remove shared and split pairs\n";  char ctemp[100];
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

void AliFemtoV0TrackPairCut::SetV0Max(Double_t aV0Max) {
  fV0Max = aV0Max;
}

Double_t AliFemtoV0TrackPairCut::GetAliFemtoV0Max() const {
  return fV0Max;
}


TList *AliFemtoV0TrackPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoV0TrackPairCut.sharequalitymax=%f", fV0Max);
  snprintf(buf, 200, "AliFemtoV0TrackPairCut.sharefractionmax=%f", fShareFractionMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void     AliFemtoV0TrackPairCut::SetRemoveSameLabel(Bool_t aRemove)
{
  fRemoveSameLabel = aRemove;
}

void AliFemtoV0TrackPairCut::SetTPCOnly(Bool_t tpconly)
{
  fTrackTPCOnly = tpconly;
}

void AliFemtoV0TrackPairCut::SetShareQualityMax(Double_t aShareQualityMax) {
  fShareQualityMax = aShareQualityMax;
}

void AliFemtoV0TrackPairCut::SetShareFractionMax(Double_t aShareFractionMax) {
  fShareFractionMax = aShareFractionMax;
}

void AliFemtoV0TrackPairCut::SetDataType(AliFemtoDataType type)
{
  fDataType = type;
}

void AliFemtoV0TrackPairCut::SetTPCEntranceSepMinimum(double dtpc)
{
  fDTPCMin = dtpc;
}

void AliFemtoV0TrackPairCut::SetTPCExitSepMinimum(double dtpc)
{
  fDTPCExitMin = dtpc;
}

void AliFemtoV0TrackPairCut::SetKstarCut(double kstar, AliFemtoParticleType firstParticle, AliFemtoParticleType secondParticle)
{
  fKstarCut = kstar; 
  fFirstParticleType = firstParticle;  //for kstar - first particle type (V0 type) 
  fSecondParticleType = secondParticle;
}

void AliFemtoV0TrackPairCut::SetMinAvgSeparation(int type, double minSep)
{
  if(type == 0) //Track-Pos
    fMinAvgSepTrackPos = minSep;
  else if(type == 1) //Track-Neg
    fMinAvgSepTrackNeg = minSep;
}
