#include "AliMuonCompactCluster.h"
#include "TString.h"
#include <cassert>

std::ostream& operator<<(std::ostream& out, const AliMuonCompactCluster& cl)
{
    out << Form("DE %04d (B %6d NB %6d)",
            cl.DetElemId(),cl.BendingManuIndex(),cl.NonBendingManuIndex());
    return out;
}

/// Return the detection element identifier of this cluster
///
/// we use those big constant vectors to avoid depending on the re-creation
/// of the CompactMapping just for that... 
/// Note however that the code for those vectors has been generated
/// by CompactMapping::GenerateStaticOffsets method (and thus, ahem...
/// cut and pasted here). Not pretty, but the DE mapping is not changing
/// often ;-) )
Int_t AliMuonCompactCluster::DetElemId() const
{

    const std::vector<int> detectionElementIds = {100,101,102,103,200,201,202,203,300,301,302,303,400,401,402,403,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025};

    const std::vector<int> detectionElementIdOffsets = {0,451,902,1353,1804,2255,2706,3157,3608,4051,4494,4937,5380,5823,6266,6709,7152,7230,7325,7408,7459,7493,7527,7578,7661,7756,7834,7929,8012,8063,8097,8131,8182,8265,8360,8440,8537,8622,8673,8707,8741,8792,8877,8974,9054,9151,9236,9287,9321,9355,9406,9491,9588,9674,9784,9895,9964,10016,10043,10061,10079,10106,10158,10227,10338,10448,10534,10644,10755,10824,10876,10903,10921,10939,10966,11018,11087,11198,11308,11394,11504,11615,11684,11736,11763,11781,11799,11826,11878,11947,12058,12168,12254,12364,12475,12544,12596,12623,12641,12659,12686,12738,12807,12918,13028,13114,13224,13344,13422,13483,13519,13546,13573,13609,13670,13748,13868,13978,14064,14174,14294,14372,14433,14469,14496,14523,14559,14620,14698,14818,14928,15014,15124,15244,15322,15383,15419,15446,15473,15509,15570,15648,15768,15878,15964,16074,16194,16272,16333,16369,16396,16423,16459,16520,16598,16718};

    Int_t absManuIndex = BendingManuIndex();
    if ( absManuIndex < 0 ) 
    {
        absManuIndex = NonBendingManuIndex();
    }

    auto p = std::lower_bound(detectionElementIdOffsets.begin(),detectionElementIdOffsets.end(),absManuIndex);

    int ix = p-detectionElementIdOffsets.begin()-1;

    assert(ix>=0 && ix < 156);

    return detectionElementIds[ix];
}

