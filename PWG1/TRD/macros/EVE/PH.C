#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TProfile.h>

#include <TRD/AliTRDgeometry.h>
#include <TRD/AliTRDcalibDB.h>
#include <TRD/AliTRDcluster.h>
#include <TRD/AliTRDseedV1.h>
#include <TRD/AliTRDtrackV1.h>
#endif

TH1* PH(const AliTRDtrackV1* track)
{
  if (!track)  return 0x0;
  
  AliTRDcluster* cls = 0;
  AliTRDseedV1 *tracklet = 0x0;
  
  TProfile* ph = 0x0;
  if(!(ph = (TProfile*)gROOT->FindObject("PH"))){
    Int_t ntb = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
    ph = new TProfile("PH", "Average PH", ntb, -.5, ntb-.5);
    ph->GetXaxis()->SetTitle("drift time [1/100ns]");
    ph->GetYaxis()->SetTitle("<PH> [a.u.]");
  } else ph->Reset();


  for (Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++) {
    if(!(tracklet = track->GetTracklet(ily))) continue;
    if(!tracklet->IsOK()) continue;
    
    for (Int_t icl = 0; icl < AliTRDseed::knTimebins; icl++) {
      if(!(cls = tracklet->GetClusters(icl))) continue;
            
	    ph->Fill(cls->GetLocalTimeBin(), cls->GetQ());
	  }
  }

  return ph;
}

