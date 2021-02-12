#include <math.h>

namespace {

	bool trackPassesCut0(AliAODTrack *t)
	{
		return true;
	}

	bool trackPassesCut1(AliAODTrack *t)
	{
		return true;
	}

	bool trackPassesCut2(AliAODTrack *t)
	{
		bool passCut = true;

		// SetMinNClustersTPC(50)
		// SetCutGeoNcrNcl(2.0, 130.0, 1.5, 0.0, 0.0)
		// SetMaxChi2PerClusterTPC(4)
		// SetAcceptKinkDaughters(kFALSE)
		// SetRequireTPCRefit(kTRUE)
		// SetRequireITSRefit(kTRUE)
		// SetClusterRequirementITS(kSPD, kAny)
		// SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1")
		// SetMaxChi2TPCConstrainedGlobal(36)
		// SetMaxdcATovertexZ(2)
		// SetDCAToVertex2D(kFALSE)
		// SetRequireSigmaToVertex(kFALSE)
		// SetMaxChi2PerClusterITS(36)

		return passCut;
	}

	bool trackPassesCut3(AliAODTrack *t)
	{
		return true;
	}

	bool trackPassesCut4(AliAODTrack *t, AliVEvent *event)
	{
		bool passCut = true;

		Double_t dz[2] = { NAN, NAN };
		Double_t cov[3] = { NAN, NAN, NAN };	
		t->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), kVeryBig, dz, cov);

		// SetRequireITSPureStandAlone(kTRUE)
		if (passCut && (!(t->GetStatus() & 0x10000000))) { passCut = false; }
		// SetPtRange(0.15, 1e+15)
		if (passCut && (t->Pt() < 0.15)) { passCut = false; }
		if (passCut && (t->Pt() > 1e+15)) { passCut = false; }
		// SetMinNClustersITS(4)
		if (passCut && (t->GetITSNcls() < 4)) { passCut = false; }
		// SetMaxDCAToVertexXY(2.4)
		// SetMaxDCAToVertexZ(3.2)
		// SetDCAToVertex2D(kTRUE)
		if (passCut && (hypot(dz[0], dz[1]) > hypot(2.4, 3.2))) { passCut = false; }
		// SetMaxChi2PerClusterITS(36)
		if (passCut && (t->GetITSchi2()) > 36) { passCut = false; }

		return passCut;
	}

	UInt_t get_local_track_cut_bits(AliAODTrack *t, AliVEvent *event)
	{
		UInt_t _local_track_cut_bits = 0;

		_local_track_cut_bits |= trackPassesCut0(t) << 0;
		_local_track_cut_bits |= trackPassesCut1(t) << 1;
		_local_track_cut_bits |= trackPassesCut2(t) << 2;
		_local_track_cut_bits |= trackPassesCut3(t) << 3;
		_local_track_cut_bits |= trackPassesCut4(t, event) << 4;

		return _local_track_cut_bits;
	}

	bool trackPassesTPCCuts(AliAODTrack *t) {
		bool passCut = true;
		passCut = passCut && trackPassesCut0(t);
		passCut = passCut && trackPassesCut1(t);
		passCut = passCut && trackPassesCut2(t);
		passCut = passCut && trackPassesCut3(t);

		return passCut;
	}
}
