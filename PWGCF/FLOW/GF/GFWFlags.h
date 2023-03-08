#ifndef GFWFLAGS__H
#define GFWFLAGS__H
namespace GFWFlags {
  enum kLocalEventFlags {
    klEventCuts =      BIT(0), //Event cuts (AliEventCuts)
    klVtxOK =          BIT(1), //Vertex rezolution
    klVtxZ10 =         BIT(2), //Vtx. z < 10 cm
    klVtxZ5 =          BIT(3), //Vtx. z < 5 cm
    klVtxZ7 =          BIT(4), //Vtx. z < 7 cm
    klVtxZ9 =          BIT(5)  //Vtx. z < 9 cm
  };
  enum kLocalTrackFlags {
    klFB32 =           BIT(0), //FB32
    klFB64 =           BIT(1), //FB64
    klFB256 =          BIT(2), //FB256
    klFB512 =          BIT(3), //FB512
    klFB96 =           BIT(4), //FB96 ( 32 || 64)
    klFB96Tuned =      BIT(5), //FB96 with fraction of shared TPC Clusters < 0.4 (for compatibility with FB768Tuned)
    klFB768 =          BIT(6), //FB768 (256 || 512)
    klFB768Tuned =     BIT(7), //FB768 where 512 FB requires a hit on first SDD (for compatibility with FB96Tuned)
    klSharedClusters = BIT(8), //fraction of shared TPC clusters < 0.4
    klHitOnSDD =       BIT(9), //hit on first layer of SDD
    klNoSPD =          BIT(10), //hit on first layer of SDD
    klDCAz20 =         BIT(11),//DCA z < 2 cm
    klDCAz10 =         BIT(12),//DCA z < 1 cm
    klDCAz05 =         BIT(13),//DCA z < 0.5 cm
    klDCAxy2010 =      BIT(14),//DCAxy cut tuned to 2010
    klDCAxy2011 =      BIT(15),//DCAxy cut tuned to 2011
    klDCAxy8Sigma =    BIT(16),//DCAxy cut on 8 sigma, 2011
    klDCAxy4Sigma =    BIT(17),//DCAxy cut on 4 sigma, 2011
    klDCAxy10Sigma =   BIT(18),//DCAxy cut on 10 sigma, 2011
    klTPCchi2PC25 =    BIT(19),//TPC chi2/cluster < 2.5
    klTPCchi2PC20 =    BIT(20),//TPC chi2/cluster < 2.0
    klTPCchi2PC30 =    BIT(21),//TPC chi2/cluster < 3.0
    klTPCchi2PC40 =    BIT(22),//TPC chi2/cluster < 4.0
    klNTPCcls70 =      BIT(23),//Number of TPC clusters 70
    klNTPCcls80 =      BIT(24),//Number of TPC clusters 80
    klNTPCcls90 =      BIT(25),//Number of TPC clusters 90
    klNTPCcls100 =     BIT(26)//Number of TPC clusters 100
  };
  enum EventFlags { kNominal=0, kVtx9, kVtx7, kVtx5, kAllEvFlags}; //Better keep these as uint_t so that we reuse them as array indeces
  enum TrackFlags { kFB96=0, kFB768,
                    kDCAz10, kDCAz05,
                    kDCA4Sigma, kDCA10Sigma,
                    kChiSq2, kChiSq3,
                    kNTPC80, kNTPC90, kNTPC100,
                    kFB96Tuned, kFB768Tuned, //These are for developing purposes and shouldn't be used!
                    kFB768DCAz,
                    kFB768DCAxyLow,
                    kFB768DCAxyHigh,
                    kFB768ChiSq2, kFB768ChiSq3,
                    kFB768nTPC, kFB96MergedDCA,
                    kChiSq25, //For testing purposes on 15_pass1 data
                    kRebinnedNUA, //For testing effects of NUA rebinning
                    kAllTrFlags
                  };
  static inline const Int_t BitIndex(const UInt_t &lFlag) {
    if(lFlag==0) return -1;
    for(Int_t i=0;i<sizeof(lFlag)*8;i++) if(lFlag&(1<<i)) return i;
    return -1;
  };
  static inline const TString GetSystPF(UInt_t lEv, UInt_t lTr) { return TString(Form("_Ev%i_Tr%i",lEv,lTr)); };
  const Int_t gNTrackFlags=22;
  const Int_t gNEventFlags=4;
};
#endif
