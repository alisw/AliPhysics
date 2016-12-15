#ifndef ALIMUONCOMPACTQUICKACCEFFCHECKER_H
#define ALIMUONCOMPACTQUICKACCEFFCHECKER_H

struct AliMuonCompactQuickAccEffChecker {

    static void CompareWithFullAccEff(const char* fullacceff="EfficiencyJPsiRun_from_astrid.root", const char* quickacceff="lhc15pp.monocathodes-accepted.root", const char* q2="lhc15pp.monocathodes-not-accepted.root");

};

#endif
