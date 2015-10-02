/*
 MCH DA for online BP occupancy evolution

 Contact: Laurent Aphecetche <laurent.aphecetche@subatech.in2p3.fr>
 Link: 
 Run Type: PHYSICS STANDALONE
 DA Type: MON
 Inputs files: mchbpevo.conf
 Number of events needed: all (or at least as much as possible...)
 Reference Run: 233721
 Link: 15000233721019.1104.root
 Output Files: mchbpevo.root, to be exported to the DAQ FXS
 Trigger types used: PHYSICS_EVENT
 */

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

///
/// MUON TRACKER DA to compute the evolution of the bus patches
/// hit count, which can be used to alert in case of high occupancy
///
/// In the end, this DA produces a Root file containing
/// histograms with the hit count all the buspatches (that were seen in the data flow)
/// of MCH as a function of time
/// (and the number of events seen, so we can later on compute the occupancy)
///
/// The time resolutions used are defined in the detDB configuration file mchbpevo.conf
/// if it exists or otherwise defaults to 60s
///
/// Note that, to avoid having to access the mapping (from the OCDB) of MCH, we
/// do rely on some mapping information (list of bus patch ids and their number of
/// pads) to be present in the configuration file for this DA.
/// We also provide defaults in this code in case the configuration is not properly setup.
///
/// $Id$

#include "AliDAQ.h"
#include "AliMergeableCollection.h"
#include "AliMpConstants.h"
#include "AliMUONBusPatchEvolution.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReaderDate.h"
#include "daqDA.h"
#include "event.h"
#include "monitor.h"
#include "Riostream.h"
#include "signal.h"
#include "TDatime.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TObjString.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTimeStamp.h"

const char* OUTPUT_FILE = "mchbpevo.root";
const char* DAVERSION = "MCHPBEVOda v0.2 ($Id$)";
const char* CONFIG_FILE = "mchbpevo.conf";
int DEFAULT_MAX_DURATION = 5*60*60; // 5 hours, in seconds
int DEFAULT_TIME_RESOLUTION = 60; // seconds
int DEFAULT_NOF_EVENTS_REQUIRED_FOR_DECISION = 10000; // ten thousand events needed before any decision can be taken
float DEFAULT_OCCUPANCY_THRESHOLD = 0.30; // by default a bus patch with more than 30% occupancy for the last 1000 is declared suspicious ;-)

//_________________________________________________________________________________________________
Int_t AssertHistogramsRange(AliMergeableCollection& hc, Int_t timeFromRef,
													  const std::vector<int>& timeResolutions,
													  const std::map<int,int>& busPatches)
{
	/// Insure the histograms in hc are big enough to accomodate timeFromRef
	/// timeFromRef can be negative or positive
	/// Return 1 if a re-allocation has been performed, 0 otherwise

	assert(timeResolutions.size()>=1);
	assert(busPatches.size()==888);

	TH1* htest = hc.Histo(Form("/BUSPATCH/HITS/%ds/BP%04d",timeResolutions[0],1));

	assert(htest!=0x0);

	TAxis* axis = htest->GetXaxis();

	Int_t expansionTime(0);

	if ( timeFromRef < axis->GetXmin() )
	{
		expansionTime = timeFromRef - axis->GetXmin();
	}

	if ( timeFromRef > axis->GetXmax() )
	{
		expansionTime = timeFromRef - axis->GetXmax();
	}

	if ( TMath::Abs(expansionTime) < timeResolutions[0] )
	{
		// current binning is OK. Nothing to do.
		return 0;
	}

	if ( TMath::Abs(expansionTime) > 12*3600 )
	{
		std::cout << DAVERSION << " ERROR : Time range expansion bigger than 12 hours. Not doing it." << std::endl;
		return 0;
	}

	Int_t currentDuration = TMath::Nint(axis->GetXmax()-axis->GetXmin());

	std::cout << Form("%s INFO : I will expand my time range by %10d s (current duration %10d s, new one %10d s)",
			DAVERSION,expansionTime,currentDuration,currentDuration+TMath::Abs(expansionTime))<< std::endl;

	// increase the time range by expansionTime, for all histograms

	for ( std::vector<int>::size_type itr = 0; itr < timeResolutions.size(); ++itr )
	{
		int timeReso = timeResolutions[itr];

		std::map<int,int>::const_iterator it;

		for ( it = busPatches.begin(); it != busPatches.end(); ++it )
		{
			int busPatchId = it->first;

			TString path;

			path.Form("/BUSPATCH/HITS/%ds/BP%04d",timeReso,busPatchId);

			TH1* h = hc.Histo(path.Data());

			assert(h!=0x0);

			TH1* hnew = AliMUONBusPatchEvolution::ExpandTimeAxis(*h,  expansionTime,  timeReso);

			hc.Remove(path.Data());

			hc.Adopt(hc.GetIdentifier(path),hnew);
		}
	}


	if  ( expansionTime < 0 )
	{
		TTimeStamp origin;

		htest = hc.Histo(Form("/BUSPATCH/HITS/%ds/BP%04d",timeResolutions[0],1));

		AliMUONBusPatchEvolution::GetTimeOffset(*htest,origin);

		std::cout << DAVERSION << " INFO : New refTime = " << origin.AsString() << std::endl;
	}

	return 1;
}

//_________________________________________________________________________________________________
void FillBusPatchesMap(std::map<int,int>& buspatches)
{
	/// Fill the map of buspatches, giving the relation between the bus patch id and
	/// the number of pads readout by this bus patch.
	///
	/// In principle (for maximum flexibility) that information should come from
	/// the configuration file. But the mapping of the detector is not a thing
	/// that is changing often, if ever, so we provide a very sensible default here.
	///
	/// Ain't pretty, but the code below can be cut and paste from the
	/// execution of the MakeConfigCodeForBPEVOda function in
	/// $ALICE_ROOT/../src/MUON/DA/MCHBPEVOdaUtils.C macro
	///
	buspatches[1236]=992;
	buspatches[1]=1272;
	buspatches[1237]=1184;
	buspatches[2]=1336;
	buspatches[3]=1528;
	buspatches[4]=1656;
	buspatches[5]=1592;
	buspatches[6]=1528;
	buspatches[7]=1528;
	buspatches[8]=1152;
	buspatches[9]=808;
	buspatches[10]=704;
	buspatches[11]=512;
	buspatches[13]=1280;
	buspatches[12]=664;
	buspatches[14]=1344;
	buspatches[15]=1536;
	buspatches[16]=1664;
	buspatches[17]=1664;
	buspatches[18]=1536;
	buspatches[19]=1536;
	buspatches[20]=1152;
	buspatches[21]=800;
	buspatches[22]=624;
	buspatches[23]=512;
	buspatches[24]=744;
	buspatches[25]=1272;
	buspatches[26]=1336;
	buspatches[27]=1528;
	buspatches[28]=1656;
	buspatches[29]=1592;
	buspatches[30]=1528;
	buspatches[31]=1528;
	buspatches[32]=1152;
	buspatches[33]=808;
	buspatches[34]=704;
	buspatches[35]=512;
	buspatches[37]=1280;
	buspatches[36]=664;
	buspatches[38]=1344;
	buspatches[39]=1536;
	buspatches[40]=1664;
	buspatches[41]=1664;
	buspatches[42]=1536;
	buspatches[43]=1536;
	buspatches[44]=1152;
	buspatches[45]=800;
	buspatches[46]=624;
	buspatches[47]=512;
	buspatches[48]=744;
	buspatches[1238]=1024;
	buspatches[1239]=576;
	buspatches[1240]=992;
	buspatches[1241]=1184;
	buspatches[1242]=512;
	buspatches[1243]=720;
	buspatches[1244]=912;
	buspatches[1245]=480;
	buspatches[1246]=608;
	buspatches[1301]=480;
	buspatches[1303]=720;
	buspatches[1302]=608;
	buspatches[1305]=576;
	buspatches[1304]=912;
	buspatches[1306]=992;
	buspatches[1307]=1184;
	buspatches[1309]=1152;
	buspatches[1308]=512;
	buspatches[1310]=992;
	buspatches[1311]=1184;
	buspatches[1313]=2176;
	buspatches[1312]=1024;
	buspatches[1314]=1264;
	buspatches[1315]=1456;
	buspatches[1317]=1664;
	buspatches[1316]=2176;
	buspatches[1318]=1184;
	buspatches[1319]=992;
	buspatches[1320]=2304;
	buspatches[1321]=832;
	buspatches[1322]=1664;
	buspatches[1323]=992;
	buspatches[1324]=1184;
	buspatches[1325]=1600;
	buspatches[1326]=1664;
	buspatches[1327]=1184;
	buspatches[1328]=992;
	buspatches[1329]=2304;
	buspatches[1331]=2176;
	buspatches[1330]=832;
	buspatches[1332]=1264;
	buspatches[1333]=1456;
	buspatches[1335]=1152;
	buspatches[1334]=2176;
	buspatches[1336]=992;
	buspatches[101]=1272;
	buspatches[1339]=576;
	buspatches[102]=1336;
	buspatches[103]=1528;
	buspatches[104]=1656;
	buspatches[105]=1592;
	buspatches[106]=1528;
	buspatches[107]=1528;
	buspatches[108]=1152;
	buspatches[109]=808;
	buspatches[110]=704;
	buspatches[111]=512;
	buspatches[113]=1280;
	buspatches[112]=664;
	buspatches[114]=1344;
	buspatches[115]=1536;
	buspatches[116]=1664;
	buspatches[117]=1664;
	buspatches[118]=1536;
	buspatches[119]=1536;
	buspatches[120]=1152;
	buspatches[121]=800;
	buspatches[122]=624;
	buspatches[123]=512;
	buspatches[124]=744;
	buspatches[125]=1272;
	buspatches[126]=1336;
	buspatches[127]=1528;
	buspatches[128]=1656;
	buspatches[129]=1592;
	buspatches[130]=1528;
	buspatches[131]=1528;
	buspatches[132]=1152;
	buspatches[133]=808;
	buspatches[134]=704;
	buspatches[135]=512;
	buspatches[137]=1280;
	buspatches[136]=664;
	buspatches[138]=1344;
	buspatches[139]=1536;
	buspatches[140]=1664;
	buspatches[141]=1664;
	buspatches[142]=1536;
	buspatches[143]=1536;
	buspatches[144]=1152;
	buspatches[145]=800;
	buspatches[146]=624;
	buspatches[147]=512;
	buspatches[148]=744;
	buspatches[1345]=480;
	buspatches[1346]=608;
	buspatches[1343]=720;
	buspatches[1344]=912;
	buspatches[1340]=992;
	buspatches[1341]=1184;
	buspatches[1342]=512;
	buspatches[1337]=1184;
	buspatches[1338]=1024;
	buspatches[1401]=480;
	buspatches[1402]=608;
	buspatches[1403]=720;
	buspatches[1404]=912;
	buspatches[1405]=576;
	buspatches[1406]=992;
	buspatches[1407]=1184;
	buspatches[1408]=512;
	buspatches[1409]=1152;
	buspatches[1410]=992;
	buspatches[1411]=1184;
	buspatches[1412]=1024;
	buspatches[1413]=2176;
	buspatches[1414]=1264;
	buspatches[1415]=1456;
	buspatches[1416]=2176;
	buspatches[1417]=1664;
	buspatches[1418]=1184;
	buspatches[1419]=992;
	buspatches[1420]=2304;
	buspatches[1421]=832;
	buspatches[1422]=1664;
	buspatches[1423]=992;
	buspatches[1424]=1184;
	buspatches[1425]=1600;
	buspatches[1426]=1664;
	buspatches[1427]=1184;
	buspatches[1428]=992;
	buspatches[1429]=2304;
	buspatches[1430]=832;
	buspatches[1431]=2176;
	buspatches[1432]=1264;
	buspatches[1433]=1456;
	buspatches[1434]=2176;
	buspatches[1435]=1152;
	buspatches[1436]=992;
	buspatches[201]=1272;
	buspatches[1437]=1184;
	buspatches[202]=1336;
	buspatches[203]=1528;
	buspatches[204]=1656;
	buspatches[205]=1592;
	buspatches[206]=1528;
	buspatches[207]=1528;
	buspatches[208]=1152;
	buspatches[209]=808;
	buspatches[210]=704;
	buspatches[211]=512;
	buspatches[213]=1280;
	buspatches[212]=664;
	buspatches[214]=1344;
	buspatches[215]=1536;
	buspatches[216]=1664;
	buspatches[217]=1664;
	buspatches[218]=1536;
	buspatches[219]=1536;
	buspatches[220]=1152;
	buspatches[221]=800;
	buspatches[222]=624;
	buspatches[223]=512;
	buspatches[224]=744;
	buspatches[225]=1272;
	buspatches[226]=1336;
	buspatches[227]=1528;
	buspatches[228]=1656;
	buspatches[229]=1592;
	buspatches[230]=1528;
	buspatches[231]=1528;
	buspatches[232]=1152;
	buspatches[233]=808;
	buspatches[234]=704;
	buspatches[235]=512;
	buspatches[237]=1280;
	buspatches[236]=664;
	buspatches[238]=1344;
	buspatches[239]=1536;
	buspatches[240]=1664;
	buspatches[241]=1664;
	buspatches[242]=1536;
	buspatches[243]=1536;
	buspatches[244]=1152;
	buspatches[245]=800;
	buspatches[246]=624;
	buspatches[247]=512;
	buspatches[248]=744;
	buspatches[1438]=1024;
	buspatches[1439]=576;
	buspatches[1440]=992;
	buspatches[1441]=1184;
	buspatches[1442]=512;
	buspatches[1443]=720;
	buspatches[1444]=912;
	buspatches[1445]=480;
	buspatches[1446]=608;
	buspatches[1501]=480;
	buspatches[1503]=720;
	buspatches[1502]=608;
	buspatches[1505]=576;
	buspatches[1504]=912;
	buspatches[1506]=992;
	buspatches[1507]=1184;
	buspatches[1509]=1152;
	buspatches[1508]=512;
	buspatches[1510]=992;
	buspatches[1511]=1184;
	buspatches[1513]=2176;
	buspatches[1512]=1024;
	buspatches[1514]=1264;
	buspatches[1515]=1456;
	buspatches[1517]=1664;
	buspatches[1516]=2176;
	buspatches[1518]=1184;
	buspatches[1519]=992;
	buspatches[1520]=2304;
	buspatches[1521]=832;
	buspatches[1522]=1664;
	buspatches[1523]=992;
	buspatches[1524]=1184;
	buspatches[1525]=1600;
	buspatches[1526]=1664;
	buspatches[1527]=1184;
	buspatches[1528]=992;
	buspatches[1529]=2304;
	buspatches[1531]=2176;
	buspatches[1530]=832;
	buspatches[1532]=1264;
	buspatches[1533]=1456;
	buspatches[1535]=1152;
	buspatches[1534]=2176;
	buspatches[1536]=992;
	buspatches[301]=1272;
	buspatches[1539]=576;
	buspatches[302]=1336;
	buspatches[303]=1528;
	buspatches[304]=1656;
	buspatches[305]=1592;
	buspatches[306]=1528;
	buspatches[307]=1528;
	buspatches[308]=1152;
	buspatches[309]=808;
	buspatches[310]=704;
	buspatches[311]=512;
	buspatches[313]=1280;
	buspatches[312]=664;
	buspatches[314]=1344;
	buspatches[315]=1536;
	buspatches[316]=1664;
	buspatches[317]=1664;
	buspatches[318]=1536;
	buspatches[319]=1536;
	buspatches[320]=1152;
	buspatches[321]=800;
	buspatches[322]=624;
	buspatches[323]=512;
	buspatches[324]=744;
	buspatches[325]=1272;
	buspatches[326]=1336;
	buspatches[327]=1528;
	buspatches[328]=1656;
	buspatches[329]=1592;
	buspatches[330]=1528;
	buspatches[331]=1528;
	buspatches[332]=1152;
	buspatches[333]=808;
	buspatches[334]=704;
	buspatches[335]=512;
	buspatches[337]=1280;
	buspatches[336]=664;
	buspatches[338]=1344;
	buspatches[339]=1536;
	buspatches[340]=1664;
	buspatches[341]=1664;
	buspatches[342]=1536;
	buspatches[343]=1536;
	buspatches[344]=1152;
	buspatches[345]=800;
	buspatches[346]=624;
	buspatches[347]=512;
	buspatches[348]=744;
	buspatches[1545]=480;
	buspatches[1546]=608;
	buspatches[1543]=720;
	buspatches[1544]=912;
	buspatches[1540]=992;
	buspatches[1541]=1184;
	buspatches[1542]=512;
	buspatches[1537]=1184;
	buspatches[1538]=1024;
	buspatches[1601]=720;
	buspatches[1602]=912;
	buspatches[1603]=304;
	buspatches[1604]=720;
	buspatches[1605]=912;
	buspatches[1606]=240;
	buspatches[1607]=1152;
	buspatches[1608]=720;
	buspatches[1609]=912;
	buspatches[1610]=1024;
	buspatches[1611]=1728;
	buspatches[1612]=720;
	buspatches[1613]=912;
	buspatches[1614]=1536;
	buspatches[1615]=2752;
	buspatches[1616]=992;
	buspatches[1617]=1184;
	buspatches[1618]=2688;
	buspatches[1619]=1664;
	buspatches[1620]=1184;
	buspatches[1621]=992;
	buspatches[1622]=2304;
	buspatches[1623]=832;
	buspatches[1624]=1664;
	buspatches[1625]=992;
	buspatches[1626]=1184;
	buspatches[1627]=1600;
	buspatches[1628]=1664;
	buspatches[1629]=1184;
	buspatches[1630]=992;
	buspatches[1631]=2304;
	buspatches[1632]=832;
	buspatches[1633]=2752;
	buspatches[1634]=992;
	buspatches[1635]=1184;
	buspatches[1636]=2688;
	buspatches[401]=1196;
	buspatches[1637]=1728;
	buspatches[402]=1260;
	buspatches[403]=1446;
	buspatches[404]=1594;
	buspatches[405]=1472;
	buspatches[406]=1404;
	buspatches[407]=1275;
	buspatches[408]=938;
	buspatches[409]=832;
	buspatches[410]=799;
	buspatches[411]=1145;
	buspatches[413]=1200;
	buspatches[412]=625;
	buspatches[414]=1264;
	buspatches[415]=1450;
	buspatches[416]=1594;
	buspatches[417]=1472;
	buspatches[418]=1408;
	buspatches[419]=1265;
	buspatches[420]=935;
	buspatches[421]=832;
	buspatches[422]=736;
	buspatches[423]=1247;
	buspatches[424]=544;
	buspatches[425]=1196;
	buspatches[426]=1260;
	buspatches[427]=1446;
	buspatches[428]=1594;
	buspatches[429]=1472;
	buspatches[430]=1404;
	buspatches[431]=1275;
	buspatches[432]=938;
	buspatches[433]=832;
	buspatches[434]=799;
	buspatches[435]=1145;
	buspatches[437]=1200;
	buspatches[436]=625;
	buspatches[438]=1264;
	buspatches[439]=1450;
	buspatches[440]=1594;
	buspatches[441]=1472;
	buspatches[442]=1408;
	buspatches[443]=1265;
	buspatches[444]=935;
	buspatches[445]=832;
	buspatches[446]=736;
	buspatches[447]=1247;
	buspatches[448]=544;
	buspatches[1638]=720;
	buspatches[1639]=912;
	buspatches[1640]=1536;
	buspatches[1641]=1152;
	buspatches[1642]=720;
	buspatches[1643]=912;
	buspatches[1644]=1024;
	buspatches[1645]=304;
	buspatches[1646]=720;
	buspatches[1647]=912;
	buspatches[1648]=240;
	buspatches[1649]=720;
	buspatches[1650]=912;
	buspatches[1701]=720;
	buspatches[1703]=304;
	buspatches[1702]=912;
	buspatches[1704]=720;
	buspatches[1705]=912;
	buspatches[1707]=1152;
	buspatches[1706]=240;
	buspatches[1708]=720;
	buspatches[1709]=912;
	buspatches[1711]=1728;
	buspatches[1710]=1024;
	buspatches[1712]=720;
	buspatches[1713]=912;
	buspatches[1715]=2752;
	buspatches[1714]=1536;
	buspatches[1716]=992;
	buspatches[1717]=1184;
	buspatches[1719]=1664;
	buspatches[1718]=2688;
	buspatches[1720]=1184;
	buspatches[1721]=992;
	buspatches[1722]=2304;
	buspatches[1723]=832;
	buspatches[1724]=1664;
	buspatches[1725]=992;
	buspatches[1726]=1184;
	buspatches[1727]=1600;
	buspatches[1728]=1664;
	buspatches[1729]=1184;
	buspatches[1730]=992;
	buspatches[1731]=2304;
	buspatches[1733]=2752;
	buspatches[1732]=832;
	buspatches[1734]=992;
	buspatches[1735]=1184;
	buspatches[1737]=1728;
	buspatches[501]=1196;
	buspatches[1738]=720;
	buspatches[502]=1260;
	buspatches[503]=1446;
	buspatches[504]=1594;
	buspatches[505]=1472;
	buspatches[506]=1404;
	buspatches[507]=1275;
	buspatches[508]=938;
	buspatches[509]=832;
	buspatches[510]=799;
	buspatches[511]=1145;
	buspatches[513]=1200;
	buspatches[512]=625;
	buspatches[514]=1264;
	buspatches[515]=1450;
	buspatches[516]=1594;
	buspatches[517]=1472;
	buspatches[518]=1408;
	buspatches[519]=1265;
	buspatches[520]=935;
	buspatches[521]=832;
	buspatches[522]=736;
	buspatches[523]=1247;
	buspatches[524]=544;
	buspatches[525]=1196;
	buspatches[526]=1260;
	buspatches[527]=1446;
	buspatches[528]=1594;
	buspatches[529]=1472;
	buspatches[530]=1404;
	buspatches[531]=1275;
	buspatches[532]=938;
	buspatches[533]=832;
	buspatches[534]=799;
	buspatches[535]=1145;
	buspatches[537]=1200;
	buspatches[536]=625;
	buspatches[538]=1264;
	buspatches[539]=1450;
	buspatches[540]=1594;
	buspatches[541]=1472;
	buspatches[542]=1408;
	buspatches[543]=1265;
	buspatches[544]=935;
	buspatches[545]=832;
	buspatches[546]=736;
	buspatches[547]=1247;
	buspatches[548]=544;
	buspatches[1749]=720;
	buspatches[1750]=912;
	buspatches[1745]=304;
	buspatches[1746]=720;
	buspatches[1747]=912;
	buspatches[1748]=240;
	buspatches[1741]=1152;
	buspatches[1742]=720;
	buspatches[1743]=912;
	buspatches[1744]=1024;
	buspatches[1739]=912;
	buspatches[1740]=1536;
	buspatches[1736]=2688;
	buspatches[1801]=720;
	buspatches[1802]=912;
	buspatches[1803]=304;
	buspatches[1804]=720;
	buspatches[1805]=912;
	buspatches[1806]=240;
	buspatches[1807]=1152;
	buspatches[1808]=720;
	buspatches[1809]=912;
	buspatches[1810]=1024;
	buspatches[1811]=1728;
	buspatches[1812]=720;
	buspatches[1813]=912;
	buspatches[1814]=1536;
	buspatches[1815]=2752;
	buspatches[1816]=992;
	buspatches[1817]=1184;
	buspatches[1818]=2688;
	buspatches[1819]=1664;
	buspatches[1820]=1184;
	buspatches[1821]=992;
	buspatches[1822]=2304;
	buspatches[1823]=832;
	buspatches[1824]=1664;
	buspatches[1825]=992;
	buspatches[1826]=1184;
	buspatches[1827]=1600;
	buspatches[1828]=1664;
	buspatches[1829]=1184;
	buspatches[1830]=992;
	buspatches[1831]=2304;
	buspatches[1832]=832;
	buspatches[1833]=2752;
	buspatches[1834]=992;
	buspatches[1835]=1184;
	buspatches[1836]=2688;
	buspatches[601]=1196;
	buspatches[1837]=1728;
	buspatches[602]=1260;
	buspatches[603]=1446;
	buspatches[604]=1594;
	buspatches[605]=1472;
	buspatches[606]=1404;
	buspatches[607]=1275;
	buspatches[608]=938;
	buspatches[609]=832;
	buspatches[610]=799;
	buspatches[611]=1145;
	buspatches[613]=1200;
	buspatches[612]=625;
	buspatches[614]=1264;
	buspatches[615]=1450;
	buspatches[616]=1594;
	buspatches[617]=1472;
	buspatches[618]=1408;
	buspatches[619]=1265;
	buspatches[620]=935;
	buspatches[621]=832;
	buspatches[622]=736;
	buspatches[623]=1247;
	buspatches[624]=544;
	buspatches[625]=1196;
	buspatches[626]=1260;
	buspatches[627]=1446;
	buspatches[628]=1594;
	buspatches[629]=1472;
	buspatches[630]=1404;
	buspatches[631]=1275;
	buspatches[632]=938;
	buspatches[633]=832;
	buspatches[634]=799;
	buspatches[635]=1145;
	buspatches[637]=1200;
	buspatches[636]=625;
	buspatches[638]=1264;
	buspatches[639]=1450;
	buspatches[640]=1594;
	buspatches[641]=1472;
	buspatches[642]=1408;
	buspatches[643]=1265;
	buspatches[644]=935;
	buspatches[645]=832;
	buspatches[646]=736;
	buspatches[647]=1247;
	buspatches[648]=544;
	buspatches[1838]=720;
	buspatches[1839]=912;
	buspatches[1840]=1536;
	buspatches[1841]=1152;
	buspatches[1842]=720;
	buspatches[1843]=912;
	buspatches[1844]=1024;
	buspatches[1845]=304;
	buspatches[1846]=720;
	buspatches[1847]=912;
	buspatches[1848]=240;
	buspatches[1849]=720;
	buspatches[1850]=912;
	buspatches[1901]=720;
	buspatches[1903]=304;
	buspatches[1902]=912;
	buspatches[1904]=720;
	buspatches[1905]=912;
	buspatches[1907]=1152;
	buspatches[1906]=240;
	buspatches[1908]=720;
	buspatches[1909]=912;
	buspatches[1911]=1728;
	buspatches[1910]=1024;
	buspatches[1912]=720;
	buspatches[1913]=912;
	buspatches[1915]=2752;
	buspatches[1914]=1536;
	buspatches[1916]=992;
	buspatches[1917]=1184;
	buspatches[1919]=1664;
	buspatches[1918]=2688;
	buspatches[1920]=1184;
	buspatches[1921]=992;
	buspatches[1922]=2304;
	buspatches[1923]=832;
	buspatches[1924]=1664;
	buspatches[1925]=992;
	buspatches[1926]=1184;
	buspatches[1927]=1600;
	buspatches[1928]=1664;
	buspatches[1929]=1184;
	buspatches[1930]=992;
	buspatches[1931]=2304;
	buspatches[1933]=2752;
	buspatches[1932]=832;
	buspatches[1934]=992;
	buspatches[1935]=1184;
	buspatches[1937]=1728;
	buspatches[701]=1196;
	buspatches[1938]=720;
	buspatches[702]=1260;
	buspatches[703]=1446;
	buspatches[704]=1594;
	buspatches[705]=1472;
	buspatches[706]=1404;
	buspatches[707]=1275;
	buspatches[708]=938;
	buspatches[709]=832;
	buspatches[710]=799;
	buspatches[711]=1145;
	buspatches[713]=1200;
	buspatches[712]=625;
	buspatches[714]=1264;
	buspatches[715]=1450;
	buspatches[716]=1594;
	buspatches[717]=1472;
	buspatches[718]=1408;
	buspatches[719]=1265;
	buspatches[720]=935;
	buspatches[721]=832;
	buspatches[722]=736;
	buspatches[723]=1247;
	buspatches[724]=544;
	buspatches[725]=1196;
	buspatches[726]=1260;
	buspatches[727]=1446;
	buspatches[728]=1594;
	buspatches[729]=1472;
	buspatches[730]=1404;
	buspatches[731]=1275;
	buspatches[732]=938;
	buspatches[733]=832;
	buspatches[734]=799;
	buspatches[735]=1145;
	buspatches[737]=1200;
	buspatches[736]=625;
	buspatches[738]=1264;
	buspatches[739]=1450;
	buspatches[740]=1594;
	buspatches[741]=1472;
	buspatches[742]=1408;
	buspatches[743]=1265;
	buspatches[744]=935;
	buspatches[745]=832;
	buspatches[746]=736;
	buspatches[747]=1247;
	buspatches[748]=544;
	buspatches[1949]=720;
	buspatches[1950]=912;
	buspatches[1945]=304;
	buspatches[1946]=720;
	buspatches[1947]=912;
	buspatches[1948]=240;
	buspatches[1941]=1152;
	buspatches[1942]=720;
	buspatches[1943]=912;
	buspatches[1944]=1024;
	buspatches[1939]=912;
	buspatches[1940]=1536;
	buspatches[1936]=2688;
	buspatches[801]=1024;
	buspatches[802]=1152;
	buspatches[803]=1536;
	buspatches[804]=1728;
	buspatches[805]=1088;
	buspatches[806]=1472;
	buspatches[807]=1648;
	buspatches[808]=1088;
	buspatches[809]=448;
	buspatches[810]=2160;
	buspatches[811]=2048;
	buspatches[812]=576;
	buspatches[813]=832;
	buspatches[814]=1408;
	buspatches[815]=960;
	buspatches[816]=1072;
	buspatches[817]=1536;
	buspatches[818]=1024;
	buspatches[819]=1152;
	buspatches[820]=1536;
	buspatches[821]=1728;
	buspatches[822]=1088;
	buspatches[823]=1536;
	buspatches[824]=1728;
	buspatches[825]=1088;
	buspatches[826]=448;
	buspatches[827]=2240;
	buspatches[828]=2112;
	buspatches[829]=576;
	buspatches[830]=832;
	buspatches[901]=1024;
	buspatches[903]=1536;
	buspatches[902]=1152;
	buspatches[905]=1088;
	buspatches[904]=1728;
	buspatches[906]=1472;
	buspatches[907]=1648;
	buspatches[909]=448;
	buspatches[908]=1088;
	buspatches[910]=2160;
	buspatches[911]=2048;
	buspatches[912]=576;
	buspatches[913]=832;
	buspatches[914]=1024;
	buspatches[915]=1152;
	buspatches[916]=1536;
	buspatches[917]=1728;
	buspatches[918]=1088;
	buspatches[919]=1536;
	buspatches[920]=1728;
	buspatches[921]=1088;
	buspatches[922]=448;
	buspatches[923]=2240;
	buspatches[924]=2112;
	buspatches[925]=576;
	buspatches[927]=1408;
	buspatches[926]=832;
	buspatches[928]=1024;
	buspatches[929]=1152;
	buspatches[930]=1536;
	buspatches[1001]=448;
	buspatches[1002]=2160;
	buspatches[1003]=2048;
	buspatches[1004]=576;
	buspatches[1005]=832;
	buspatches[1006]=1088;
	buspatches[1007]=1472;
	buspatches[1008]=1648;
	buspatches[1009]=1088;
	buspatches[1010]=1536;
	buspatches[1011]=1728;
	buspatches[1012]=1024;
	buspatches[1013]=1152;
	buspatches[1014]=1408;
	buspatches[1015]=1024;
	buspatches[1016]=1152;
	buspatches[1017]=1536;
	buspatches[1018]=448;
	buspatches[1019]=2240;
	buspatches[1020]=2112;
	buspatches[1021]=576;
	buspatches[1022]=832;
	buspatches[1023]=1088;
	buspatches[1024]=1536;
	buspatches[1025]=1728;
	buspatches[1026]=1088;
	buspatches[1027]=1536;
	buspatches[1028]=1728;
	buspatches[1029]=1024;
	buspatches[1030]=1152;
	buspatches[1101]=1408;
	buspatches[1102]=960;
	buspatches[1103]=1072;
	buspatches[1104]=1536;
	buspatches[1105]=448;
	buspatches[1106]=2160;
	buspatches[1107]=2048;
	buspatches[1108]=576;
	buspatches[1109]=832;
	buspatches[1110]=1088;
	buspatches[1111]=1472;
	buspatches[1112]=1648;
	buspatches[1113]=1088;
	buspatches[1114]=1536;
	buspatches[1115]=1728;
	buspatches[1116]=1024;
	buspatches[1117]=1152;
	buspatches[1118]=448;
	buspatches[1119]=2240;
	buspatches[1120]=2112;
	buspatches[1121]=576;
	buspatches[1123]=1088;
	buspatches[1122]=832;
	buspatches[1124]=1536;
	buspatches[1125]=1728;
	buspatches[1127]=1536;
	buspatches[1126]=1088;
	buspatches[1129]=1024;
	buspatches[1128]=1728;
	buspatches[1130]=1152;
	buspatches[1201]=480;
	buspatches[1202]=608;
	buspatches[1203]=720;
	buspatches[1204]=912;
	buspatches[1205]=576;
	buspatches[1206]=992;
	buspatches[1207]=1184;
	buspatches[1208]=512;
	buspatches[1209]=1152;
	buspatches[1210]=992;
	buspatches[1211]=1184;
	buspatches[1212]=1024;
	buspatches[1213]=2176;
	buspatches[1214]=1264;
	buspatches[1215]=1456;
	buspatches[1216]=2176;
	buspatches[1217]=1664;
	buspatches[1218]=1184;
	buspatches[1219]=992;
	buspatches[1220]=2304;
	buspatches[1221]=832;
	buspatches[1222]=1664;
	buspatches[1223]=992;
	buspatches[1224]=1184;
	buspatches[1225]=1600;
	buspatches[1226]=1664;
	buspatches[1227]=1184;
	buspatches[1228]=992;
	buspatches[1229]=2304;
	buspatches[1230]=832;
	buspatches[1231]=2176;
	buspatches[1232]=1264;
	buspatches[1233]=1456;
	buspatches[1234]=2176;
	buspatches[1235]=1152;
}

//_________________________________________________________________________________________________
void FillCollection(AliMergeableCollection& hc, UInt_t timeOrigin, int maxDuration, int timeResolution, const std::map<int,int>& buspatches)
{
	/// Create the BusPatch and Event histograms
	///

	assert(buspatches.size()==888);

	TTimeStamp origin(static_cast<time_t>(timeOrigin),0);

	std::cout << DAVERSION << " INFO : Using refTime = " << origin.AsString() << std::endl;

	Double_t xmin =  0;
	Double_t xmax = xmin + maxDuration*1.0;

	int nbins = TMath::Nint((xmax-xmin)/timeResolution);

	// basis for all plot = per buspatch

	std::map<int,int>::const_iterator it;

	for ( it = buspatches.begin(); it != buspatches.end(); ++it)
	{
		int busPatchId = it->first;
		int nofPads = it->second;

		TH1* h = new TH1F(Form("BP%04d",busPatchId),Form("Number of hits in %d s bins",timeResolution),nbins,xmin,xmax);

		h->GetXaxis()->SetTimeDisplay(1);
		h->GetXaxis()->SetTimeFormat("%d/%m/%y %H:%M");
		h->GetXaxis()->SetTimeOffset(timeOrigin,"gmt");
		hc.Adopt(Form("/BUSPATCH/HITS/%ds",timeResolution),h);
	}

	// number of events needed for normalization

	TH1* h = new TH1F(Form("Nevents%ds",timeResolution),Form("Number of events %d s bins",timeResolution),nbins,xmin,xmax);

	h->GetXaxis()->SetTimeDisplay(1);
	h->GetXaxis()->SetTimeFormat("%d/%m/%y %H:%M");
	h->GetXaxis()->SetTimeOffset(timeOrigin,"gmt");

	hc.Adopt("",h);
}

//______________________________________________________________________________
void ReportFaultyBusPatches(AliMergeableCollection& hc,
		const std::map<int,int>& buspatches,
		int timeResolution,
		int requiredEvents,
		float occupancyThreshold)
{
	/// Loop over all bus patches and compute their occupancy in the last nevents

	// find how many bins should be considered, by finding in the Nevents histogram
	// how many of the latest bins are required to get an integral >= requiredEvents

	AliMUONBusPatchEvolution bpevo(hc,buspatches);

	std::map<int,double> faultyBP;

	Bool_t ok = bpevo.GetFaultyBusPatches(timeResolution,requiredEvents,occupancyThreshold,faultyBP);

	if (!ok) return;

	std::map<int,double>::const_iterator it;

	if ( faultyBP.empty() ) return;

	std::cout << Form("%s WARNING : %3d bus patches above occupancy threshold (of %7.2f %%) : ",DAVERSION,(int)faultyBP.size(),occupancyThreshold*100.0);

	for ( it = faultyBP.begin(); it != faultyBP.end(); ++it )
	{
		int busPatchId = it->first;
		std::cout << Form("%4d ",busPatchId);
	}

	std::cout << std::endl;
	std::cout << Form("%s INFO : Faulty bus patches occupancies : ",DAVERSION);


	for ( it = faultyBP.begin(); it != faultyBP.end(); ++it )
	{
		int busPatchId = it->first;
		double busPatchOccupancy = it->second;
		std::cout << Form("%4d (%7.2f %%) ",busPatchId,busPatchOccupancy*100.0);
	}

	std::cout << std::endl;
}

//______________________________________________________________________________
int main(int argc, char **argv) 
{
	/// Main method.

	signal(SIGSEGV,SIG_DFL); // to be able to get core dumps...

	TStopwatch timers;
	timers.Start(kTRUE);

	ios::sync_with_stdio();

	std::cout << "Running " << DAVERSION << std::endl;

	if ( argc != 2 )
	{
		std::cout << "Wrong number of arguments" << std::endl;
		std::cout << "Usage : " << argv[0] << " datasource" << std::endl;
		return -1;
	}

	TString dataSource=argv[1];

	UInt_t refTime = 0;
	TString outputFileName = OUTPUT_FILE;

	// needed for streamer application
	gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
			"*",
			"TStreamerInfo",
			"RIO",
			"TStreamerInfo()");

	Int_t numberOfEvents(0);
	Int_t numberOfPhysicsEvent(0);
	Int_t numberOfBadEvents(0);
	Int_t numberOfUsedEvents(0);
	Int_t numberOfIncompleteEvents(0);
	Int_t numberOfFlushedEvents(0);
	Int_t numberOfEventsWithMCH(0);
	Int_t numberOfEventsSinceLastDecision(0);

	UInt_t runNumber(0);

	std::vector<int> timeResolutions;
	std::map<int,int> busPatches;
	int maxDuration;
	int nofEventsRequiredForDecision;
	float occupancyThreshold;

	int status=daqDA_DB_getFile(CONFIG_FILE,CONFIG_FILE);

	if (status)
	{
		std::cout << DAVERSION << " WARNING: Could not get configuration file " << CONFIG_FILE
				<< " from DAQ detector database. Will use (hopefully sensible) defaults" << std::endl;

		timeResolutions.push_back(DEFAULT_TIME_RESOLUTION);
		maxDuration = DEFAULT_MAX_DURATION;
		nofEventsRequiredForDecision = DEFAULT_NOF_EVENTS_REQUIRED_FOR_DECISION;
		occupancyThreshold = DEFAULT_OCCUPANCY_THRESHOLD;
		FillBusPatchesMap(busPatches);
	}
	else
	{
		TEnv env(CONFIG_FILE);

		occupancyThreshold = env.GetValue("occupancyThreshold",DEFAULT_OCCUPANCY_THRESHOLD);
		nofEventsRequiredForDecision = env.GetValue("nofEventsRequiredForDecision",DEFAULT_NOF_EVENTS_REQUIRED_FOR_DECISION);

		TString tr = env.GetValue("timeResolutions","");
		TObjArray* a = tr.Tokenize(" ");
		TObjString* s;
		TIter next(a);
		while ( ( s = static_cast<TObjString*>(next())))
		{
			if (s->String().Atoi())
			{
				timeResolutions.push_back(s->String().Atoi());
			}
		}
		delete a;

		if ( timeResolutions.size() == 0 )
		{
			std::cout << DAVERSION << " WARNING : Configuration file " << CONFIG_FILE << " contains an incorrect timeResolutions setting. "
					<< " Using default of " << DEFAULT_TIME_RESOLUTION << " s instead"
					<< std::endl;
			timeResolutions.push_back(DEFAULT_TIME_RESOLUTION);
		}

		maxDuration = env.GetValue("maxDuration",DEFAULT_MAX_DURATION);

		TString sbp = env.GetValue("bp","");

		TObjArray* b = sbp.Tokenize(" ");
		TIter nextBP(b);

		int total(0);

		while ( ( s = static_cast<TObjString*>(nextBP())))
		{
			TString bpinfo = s->String();
			int ix = bpinfo.Index(':');
			if (ix<=0) continue;
			int bpid = TString(bpinfo(0,ix)).Atoi();
			int npads = TString(bpinfo(ix+1,bpinfo.Length()-ix)).Atoi();
			total += npads;
			busPatches[bpid]=npads;
		}

		delete b;

		if (busPatches.size() != 888)
		{
			std::cout << DAVERSION << " ERROR : Configuration file is incorrect : wrong number of bus patches : "
					<< " got " << busPatches.size() << " but expected 888. Cannot work like that." << std::endl;
			return -3;
		}

		if (busPatches.size() > 0 && total != 1064008)
		{
			std::cout << DAVERSION << " ERROR : Configuration file is incorrect : wrong total number of pads ! Got " << total
					<< " but expected 1064008. Cannot work like that." << std::endl;
			return -2;
		}

		if ( busPatches.size() == 0 )
		{
			std::cout << DAVERSION << " WARNING : Configuration file " << CONFIG_FILE << " does not contain the list of buspath ids (and their "
					" corresponding number of pads). Using default values instead"
					<< std::endl;
			FillBusPatchesMap(busPatches);
		}
	}

	/// We sort the time resolutions so the finer one will be used in AssertHistogramsRange
	std::sort(timeResolutions.begin(),timeResolutions.end());

	AliMergeableCollection* hc(0x0);

	AliRawReaderDate* rawReader(0x0);

	// define data source :
	status=monitorSetDataSource(dataSource.Data());
	if (status!=0)
	{
		std::cout << DAVERSION << " ERROR : monitorSetDataSource() failed: " << monitorDecodeError(status) << endl;
		return -1;
	}

	// Declare monitoring program
	status=monitorDeclareMp("MUON_TRK_BPEVO");
	if (status!=0)
	{
		std::cout << DAVERSION << " ERROR : monitorDeclareMp() failed: " << monitorDecodeError(status) << endl;
		return -1;
	}

	// Define wait event timeout - 1s max
	monitorSetNowait();
	monitorSetNoWaitNetworkTimeout(1000);

	int nofRealloc(0);

	for(;;)
	{
		struct eventHeaderStruct *event(0x0);
		eventTypeType eventT;

		status=monitorGetEventDynamic((void **)&event);
		if (status!=0)
		{
			std::cout << DAVERSION << " ERROR : " << monitorDecodeError(status) << std::endl;
			delete event;
			break;
		}

		/* check shutdown condition */
		if (daqDA_checkShutdown())
		{
			std::cout << DAVERSION << " WARNING : I am requested to stop, so I will stop..." << std::endl;
			delete event;
			break;
		}

		/* retry if got no event */
		if (event==NULL) continue;

		++numberOfEvents;

		eventT=event->eventType;

		if ((eventT == END_OF_RUN)||(eventT == END_OF_RUN_FILES))
		{
			delete event;
			break;
		}

		if (eventT != PHYSICS_EVENT)
		{
			delete event;
			continue;
		}

		if( TEST_SYSTEM_ATTRIBUTE(event->eventTypeAttribute, ATTR_INCOMPLETE_EVENT))
		{
			delete event;
			++numberOfIncompleteEvents;
			continue;
		}

		if( TEST_SYSTEM_ATTRIBUTE(event->eventTypeAttribute, ATTR_FLUSHED_EVENT))
		{
			delete event;
			++numberOfFlushedEvents;
			continue;
		}

		++numberOfPhysicsEvent;

		rawReader = new AliRawReaderDate((void*)event);

		if ( rawReader->GetRunNumber() != runNumber )
		{
			if ( runNumber != 0 )
			{
				std::cout << DAVERSION << " ERROR : Uh oh. That's bad... Changing of run number ???" << std::endl;
				delete event;
				delete rawReader;
				return -9999;
			}

			runNumber = rawReader->GetRunNumber();
		}

		// *without* using the MUONTRK event decoder, we update the number
		// of events where we have information about MUONTRK.
		// this should be what is called subevents in the logbook
		// we do *not* use the MCH decoder on purpose, to not mistakenly
		// believe there's no event if the data is corrupted (and thus the decoder
		// "sees" nothing).

		Bool_t mchThere(kFALSE);

		for ( int iDDL = 0; iDDL < AliDAQ::NumberOfDdls("MUONTRK") && !mchThere; ++iDDL )
		{
			rawReader->Reset();
			rawReader->Select("MUONTRK",iDDL,iDDL);
			if (rawReader->ReadHeader() )
			{
				if (rawReader->GetEquipmentSize() ) mchThere = kTRUE;
			}
		}

		if ( mchThere )
		{
			++numberOfEventsWithMCH;
		}
		else
		{
			delete event;
			delete rawReader;
			continue;
		}

		rawReader->Reset();

		// now do our real work with the MCH decoder

		AliMUONRawStreamTrackerHP stream(rawReader);

		stream.DisableWarnings();
		stream.First();

		Int_t buspatchId;
		UShort_t  manuId;
		UChar_t manuChannel;
		UShort_t adc;

		std::map<int,int> bpValues;

		while ( stream.Next(buspatchId,manuId,manuChannel,adc,kTRUE) )
		{
			bpValues[buspatchId]++;
		}

		if ( bpValues.empty() )
		{
			static int emptyevent(0);

			emptyevent++;

			std::cout << DAVERSION << " WARNING : Empty event" << emptyevent << std::endl;
		}

		Bool_t badEvent = stream.HasPaddingError() || stream.HasGlitchError();

		if ( !badEvent )
		{
			if (!hc)
			{
				hc = new AliMergeableCollection("bpevo");

				assert(refTime==0);

				refTime = rawReader->GetTimestamp();

				for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
				{
					FillCollection(*hc, refTime,maxDuration,timeResolutions[is],busPatches);
				}
			}
			else
			{
				assert(refTime>0);
				nofRealloc += AssertHistogramsRange(*hc,rawReader->GetTimestamp()-refTime,timeResolutions,busPatches);
			}

			++numberOfUsedEvents;

			++numberOfEventsSinceLastDecision;

			if ( numberOfEventsSinceLastDecision > nofEventsRequiredForDecision )
			{
				ReportFaultyBusPatches(*hc,busPatches,timeResolutions[timeResolutions.size()-1],nofEventsRequiredForDecision,occupancyThreshold);
				numberOfEventsSinceLastDecision=0;
			}

			for ( std::map<int,int>::const_iterator it = bpValues.begin(); it != bpValues.end(); ++it )
			{
				const int& buspatchId = it->first;
				const int& bpvalue = it->second;

				TString bpName = Form("BP%04d",buspatchId);

				for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
				{
					TString hname;

					hname.Form("/BUSPATCH/HITS/%ds/%s",timeResolutions[is],bpName.Data());
					TH1* h = hc->Histo(hname.Data());

					if (!h)
					{
						std::cout << DAVERSION << " ERROR : histogram " << hname.Data() << " not found" << std::endl;
						continue;
					}

					h->Fill(rawReader->GetTimestamp()-refTime,bpvalue);
				}
			}

			for ( std::vector<int>::size_type is = 0; is < timeResolutions.size(); ++is )
			{
				TString hname;

				hname.Form("Nevents%ds",timeResolutions[is]);
				TH1* h = hc->Histo(hname.Data());

				if (!h)
				{
					std::cout << DAVERSION << " ERROR : histogram " << hname.Data() << " not found" << std::endl;
					continue;
				}

				h->Fill(rawReader->GetTimestamp()-refTime);
			}

		}
		else
		{
			++numberOfBadEvents;
		}

		delete event;

		delete rawReader;
	}


	std::cout << DAVERSION << Form(" INFO : %12d events processed : %12d physics %d used ones %d bad ones [ %d with MCH information ] [ %d incomplet %d flushed ]",
			numberOfEvents,numberOfPhysicsEvent,numberOfUsedEvents,numberOfBadEvents,numberOfEventsWithMCH,
			numberOfIncompleteEvents,numberOfFlushedEvents) << std::endl;

	std::cout << DAVERSION << " INFO : Number of histogram re-allocation : " << nofRealloc << std::endl;

	if ( hc )
	{
		TFile* f = TFile::Open(outputFileName,"RECREATE");
		hc->Write();
		delete f;
	}

  /* store the result file on FXS */
  if (daqDA_FES_storeFile(outputFileName,"BPEVO")) return -9;

	timers.Stop();

	std::cout << DAVERSION <<  Form(" INFO : Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime()) << std::endl;

	return 0;
}
