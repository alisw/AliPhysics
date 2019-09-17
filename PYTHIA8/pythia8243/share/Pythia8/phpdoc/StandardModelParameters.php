<html>
<head>
<title>Standard-Model Parameters</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='StandardModelParameters.php'>
 
<h2>Standard-Model Parameters</h2> 
<ol id="toc">
  <li><a href="#section0">The strong coupling</a></li>
  <li><a href="#section1">The electromagnetic coupling</a></li>
  <li><a href="#section2">The electroweak couplings</a></li>
  <li><a href="#section3">The quark weak-mixing matrix</a></li>
  <li><a href="#section4">The CoupSM class</a></li>
</ol>

 
<a name="section0"></a> 
<h3>The strong coupling</h3> 
 
The <code>AlphaStrong</code> class is used to provide a first- or 
second-order running <i>alpha_strong</i> (or, trivially, a 
zeroth-order fixed one). Formulae are the standard ones found in 
[<a href="Bibliography.php#refYao06" target="page">Yao06</a>]. The second-order expression used, eq. (9.5), 
may be somewhat different in other approaches (with differences 
formally of higher order), so do not necessarily expect perfect 
agreement, especially not at small <i>Q^2</i> scales. The starting 
<i>alpha_strong</i> value is defined at the <i>M_Z</i> mass scale. 
The <i>Lambda</i> values are matched at the <i>c</i>, <i>b</i> 
and <i>t</i> flavour thresholds, 
such that <i>alpha_strong</i> is continuous. 
For second-order matching an approximate iterative method is used. 
 
<p/> 
For backwards compatibility, 
the following global switch determines whether 5- or 6-flavour running 
will be used above the <i>t</i> threshold: 
<br/><br/><table><tr><td><strong>StandardModel:alphaSnfmax  </td><td>  &nbsp;&nbsp;(<code>default = <strong>6</strong></code>; <code>minimum = 5</code>; <code>maximum = 6</code>)</td></tr></table>
<br/>
<input type="radio" name="1" value="5"><strong>5 </strong>: Use 5-flavour running for all scales above the  <ei>b</ei> flavour threshold (old default).<br/>
<input type="radio" name="1" value="6" checked="checked"><strong>6 </strong>: Use 6-flavour running above the <ei>t</ei> threshold  (new default).<br/>
 
<p/> 
Since we allow <i>alpha_strong</i> to vary separately for 
hard processes, timelike showers, spacelike showers and  multiparton 
interactions, all other relevant values are set in each of these classes. 
The default behaviour is everywhere first-order running. 
 
<p/> 
The <i>alpha_strong</i> calculation is initialized by 
<code>init( value, order, nfmax)</code>, where <code>value</code> 
is the <i>alpha_strong</i> value at <i>M_Z</i>, <code>order</code> 
is the order of the running, 0, 1 or 2, and <code>nfmax</code> 
is the highest number of flavours to include in the running. Thereafter 
the value can be calculated by <code>alphaS(scale2)</code>, where 
<code>scale2</code> is the <i>Q^2</i> scale in GeV^2. 
 
<p/> 
By default the charm, bottom and top threshold-matching mass values 
are chosen to be 1.5, 4.8 and 171 GeV, respectively. The 
<code>setThresholds(double mc, double mb, double mt)</code> 
method can be invoked to select other values. To take effect, this 
must be done before the <code>AlphaStrong::init()</code> method is called, 
since this is where the flavour-dependent <i>Lambda_i</i> values are 
calculated and stored. If in doubt, better call it once again. 
 
<p/> 
For applications inside shower programs, a second-order <code>alpha_s</code> 
value can be obtained as the product of the two functions 
<code>alphaS1Ord(scale2)</code> and <code>alphaS2OrdCorr(scale2)</code>, 
where the first gives a simple first-order running (but with the 
second-order <i>Lambda</i>) and the second the correction factor, 
below unity, for the second-order terms. This allows a compact handling 
of evolution equations. 
 
<p/> 
Resummation arguments [<a href="Bibliography.php#refCat91" target="page">Cat91</a>] show that a set of 
universal QCD corrections can be absorbed in coherent parton showers by 
applying the so-called CMW rescaling of the MSbar value of 
<i>Lambda_QCD</i>. This can be accomplished via a fourth (optional) 
boolean argument to <code>init( value, order, nfmax, useCMW)</code>, 
with default value <code>useCMW = false</code>. When set to 
<code>true</code>, the translation amounts to an <i>N_F</i>-dependent 
rescaling of <i>Lambda_QCD</i>, relative to its MSbar value, by 
a factor 1.661 for NF=3, 1.618 for NF=4, 1.569 for NF=5, 
and 1.513 for NF=6. When using this option, 
be aware that the original CMW arguments were derived using two-loop running 
and that the CMW rescaling may need be taken into account in the context of 
matrix-element matching. Note also that this option has only been made 
available for timelike and spacelike showers, not for hard processes. 
 
<a name="section1"></a> 
<h3>The electromagnetic coupling</h3> 
 
The <code>AlphaEM</code> class is used to generate a running 
<i>alpha_em</i>. The input <code>StandardModel:alphaEMmZ</code> 
value at the <i>M_Z</i> mass is matched to a low-energy behaviour 
with running starting at the electron mass threshold. The matching 
is done by fitting an effective running coefficient in the region 
between the light-quark threshold and the charm/tau threshold. This 
procedure is approximate, but good enough for our purposes. 
 
<p/> 
Since we allow <i>alpha_em</i> to vary separately for 
hard processes, timelike showers, spacelike showers and  multiparton 
interactions, the choice between using a fixed or a running 
<i>alpha_em</i> can be made in each of these classes. 
The default behaviour is everywhere first-order running. 
The actual values assumed at zero momentum transfer and 
at <i>M_Z</i> are only set here, however. 
 
<br/><br/><table><tr><td><strong>StandardModel:alphaEM0 </td><td></td><td> <input type="text" name="2" value="0.00729735" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.00729735</strong></code>; <code>minimum = 0.0072973</code>; <code>maximum = 0.0072974</code>)</td></tr></table>
The <i>alpha_em</i> value at vanishing momentum transfer 
(and also below <i>m_e</i>). 
   
 
<br/><br/><table><tr><td><strong>StandardModel:alphaEMmZ </td><td></td><td> <input type="text" name="3" value="0.00781751" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.00781751</strong></code>; <code>minimum = 0.00780</code>; <code>maximum = 0.00783</code>)</td></tr></table>
The <i>alpha_em</i> value at the <i>M_Z</i> mass scale. 
Default is taken from [<a href="Bibliography.php#refYao06" target="page">Yao06</a>]. 
   
 
<p/> 
The <i>alpha_em</i> calculation is initialized by 
<code>init(order)</code>, where <code>order</code> is the order of 
the running, 0 or 1, with -1 a special option to use the fix value 
provided at <i>M_Z</i>.   Thereafter the value can be 
calculated by <code>alphaEM(scale2)</code>, where 
<code>scale2</code> is the <i>Q^2</i> scale in GeV^2. 
 
<a name="section2"></a> 
<h3>The electroweak couplings</h3> 
 
There are two degrees of freedom that can be set, related to the 
electroweak mixing angle: 
 
<br/><br/><table><tr><td><strong>StandardModel:sin2thetaW </td><td></td><td> <input type="text" name="4" value="0.2312" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2312</strong></code>; <code>minimum = 0.225</code>; <code>maximum = 0.240</code>)</td></tr></table>
The sine-squared of the weak mixing angle, as used in all <i>Z^0</i> 
and <i>W^+-</i> masses and couplings, except for the vector couplings 
of fermions to the <i>Z^0</i>, see below. Default is the MSbar value 
from [<a href="Bibliography.php#refYao06" target="page">Yao06</a>]. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:sin2thetaWbar </td><td></td><td> <input type="text" name="5" value="0.2315" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2315</strong></code>; <code>minimum = 0.225</code>; <code>maximum = 0.240</code>)</td></tr></table>
The sine-squared of the weak mixing angle, as used to derive the vector 
couplings of fermions to the <i>Z^0</i>, in the relation 
<i>v_f = a_f - 4 e_f sin^2(theta_W)bar</i>. Default is the 
effective-angle value from [<a href="Bibliography.php#refYao06" target="page">Yao06</a>]. 
   
 
<p/> 
The Fermi constant is not much used in the currently coded matrix elements, 
since it is redundant, but it is available: 
 
<br/><br/><table><tr><td><strong>StandardModel:GF </td><td></td><td> <input type="text" name="6" value="1.16637e-5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.16637e-5</strong></code>; <code>minimum = 1.0e-5</code>; <code>maximum = 1.3e-5</code>)</td></tr></table>
The Fermi coupling constant, in units of GeV<i>^-2</i>. 
   
 
<a name="section3"></a> 
<h3>The quark weak-mixing matrix</h3> 
 
The absolute values of the Cabibbo-Kobayashi-Maskawa matrix elements are 
set by the following nine real values taken from [<a href="Bibliography.php#refYao06" target="page">Yao06</a>] - 
currently the CP-violating phase is not taken into account in this 
parametrization. It is up to the user to pick a consistent unitary 
set of new values whenever changes are made. 
 
<br/><br/><table><tr><td><strong>StandardModel:Vud </td><td></td><td> <input type="text" name="7" value="0.97383" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.97383</strong></code>; <code>minimum = 0.973</code>; <code>maximum = 0.975</code>)</td></tr></table>
The <i>V_ud</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vus </td><td></td><td> <input type="text" name="8" value="0.2272" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2272</strong></code>; <code>minimum = 0.224</code>; <code>maximum = 0.230</code>)</td></tr></table>
The <i>V_us</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vub </td><td></td><td> <input type="text" name="9" value="0.00396" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.00396</strong></code>; <code>minimum = 0.0037</code>; <code>maximum = 0.0042</code>)</td></tr></table>
The <i>V_ub</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vcd </td><td></td><td> <input type="text" name="10" value="0.2271" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2271</strong></code>; <code>minimum = 0.224</code>; <code>maximum = 0.230</code>)</td></tr></table>
The <i>V_cd</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vcs </td><td></td><td> <input type="text" name="11" value="0.97296" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.97296</strong></code>; <code>minimum = 0.972</code>; <code>maximum = 0.974</code>)</td></tr></table>
The <i>V_cs</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vcb </td><td></td><td> <input type="text" name="12" value="0.04221" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.04221</strong></code>; <code>minimum = 0.0418</code>; <code>maximum = 0.0426</code>)</td></tr></table>
The <i>V_cb</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vtd </td><td></td><td> <input type="text" name="13" value="0.00814" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.00814</strong></code>; <code>minimum = 0.006</code>; <code>maximum = 0.010</code>)</td></tr></table>
The <i>V_td</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vts </td><td></td><td> <input type="text" name="14" value="0.04161" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.04161</strong></code>; <code>minimum = 0.039</code>; <code>maximum = 0.043</code>)</td></tr></table>
The <i>V_ts</i> CKM matrix element. 
   
 
<br/><br/><table><tr><td><strong>StandardModel:Vtb </td><td></td><td> <input type="text" name="15" value="0.9991" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9991</strong></code>; <code>minimum = 0.99907</code>; <code>maximum = 0.9992</code>)</td></tr></table>
The <i>V_tb</i> CKM matrix element. 
   
 
<a name="section4"></a> 
<h3>The CoupSM class</h3> 
 
The <code><?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>Pythia</a></code> class contains a 
public instance <code>coupSM</code> of the <code>CoupSM</code> class. 
This class contains one instance each of the <code>AlphaStrong</code> 
and <code>AlphaEM</code> classes, and additionally stores the weak couplings 
and the quark mixing matrix mentioned above. This class is used especially 
in the calculation of cross sections and resonance widths, but could also 
be used elsewhere. Specifically, as already mentioned, there are separate 
<code>AlphaStrong</code> and <code>AlphaEM</code> instances for timelike 
and spacelike showers and for multiparton interactions, while weak couplings 
and the quark mixing matrix are only stored here. With the exception of the 
first two methods below, which are for internal use, the subsequent ones 
could also be used externally. 
 
<a name="anchor1"></a>
<p/><strong> CoupSM::CoupSM() &nbsp;</strong> <br/>
the constructor does nothing. Internal. 
   
 
<a name="anchor2"></a>
<p/><strong> void CoupSM::init(Settings& settings, Rndm* rndmPtr) &nbsp;</strong> <br/>
this is where the <code>AlphaStrong</code> and <code>AlphaEM</code> 
instances are initialized, and weak couplings and the quark mixing matrix 
are read in and set. This is based on the values stored on this page and 
among the <?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>Couplings and Scales</a>. 
Internal. 
   
 
<a name="anchor3"></a>
<p/><strong> double CoupSM::alphaS(double scale2) &nbsp;</strong> <br/>
the <i>alpha_strong</i> value at the quadratic scale <code>scale2</code>. 
   
 
<a name="anchor4"></a>
<p/><strong> double CoupSM::alphaS1Ord(double scale2) &nbsp;</strong> <br/>
a first-order overestimate of the full second-order <i>alpha_strong</i> 
value at the quadratic scale <code>scale2</code>. 
   
 
<a name="anchor5"></a>
<p/><strong> double CoupSM::alphaS2OrdCorr(double scale2) &nbsp;</strong> <br/>
a multiplicative correction factor, below unity, that brings the 
first-order overestimate above into agreement with the full second-order 
<i>alpha_strong</i> value at the quadratic scale <code>scale2</code>. 
   
 
<a name="anchor6"></a>
<p/><strong> double CoupSM::Lambda3() &nbsp;</strong> <br/>
   
<a name="anchor7"></a>
<strong> double CoupSM::Lambda4() &nbsp;</strong> <br/>
   
<a name="anchor8"></a>
<strong> double CoupSM::Lambda5() &nbsp;</strong> <br/>
the three-, four-, and five-flavour <i>Lambda</i> scale. 
   
 
<a name="anchor9"></a>
<p/><strong> double CoupSM::alphaEM(double scale2) &nbsp;</strong> <br/>
the <i>alpha_em</i> value at the quadratic scale <code>scale2</code>. 
   
 
<a name="anchor10"></a>
<p/><strong> double CoupSM::sin2thetaW() &nbsp;</strong> <br/>
   
<a name="anchor11"></a>
<strong> double CoupSM::cos2thetaW() &nbsp;</strong> <br/>
the sine-squared and cosine-squared of the weak mixing angle, as used in 
the gauge-boson sector. 
   
 
<a name="anchor12"></a>
<p/><strong> double CoupSM::sin2thetaWbar() &nbsp;</strong> <br/>
the sine-squared of the weak mixing angle, as used to derive the vector 
couplings of fermions to the <i>Z^0</i>. 
   
 
<a name="anchor13"></a>
<p/><strong> double CoupSM::GF() &nbsp;</strong> <br/>
the Fermi constant of weak decays, in GeV<i>^-2</i>. 
   
 
<a name="anchor14"></a>
<p/><strong> double CoupSM::ef(int idAbs) &nbsp;</strong> <br/>
the electrical charge of a fermion, by the absolute sign of the PDF code, 
i.e. <code>idAbs</code> must be in the range between 1 and 18. 
   
 
<a name="anchor15"></a>
<p/><strong> double CoupSM::vf(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor16"></a>
<strong> double CoupSM::af(int idAbs) &nbsp;</strong> <br/>
the vector and axial charges of a fermion, by the absolute sign of the PDF 
code (<i>a_f = +-1, v_f = a_f - 4. * sin2thetaWbar * e_f</i>). 
   
 
<a name="anchor17"></a>
<p/><strong> double CoupSM::t3f(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor18"></a>
<strong> double CoupSM::lf(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor19"></a>
<strong> double CoupSM::rf(int idAbs) &nbsp;</strong> <br/>
the weak isospin, left- and righthanded charges of a fermion, by the 
absolute sign of the PDF code (<i>t^3_f = a_f/2, l_f = (v_f + a_f)/2, 
r_f = (v_f - a_f)/2</i>; you may find other conventions in the literature 
that differ by a factor of 2). 
   
 
<a name="anchor20"></a>
<p/><strong> double CoupSM::ef2(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor21"></a>
<strong> double CoupSM::vf2(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor22"></a>
<strong> double CoupSM::af2(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor23"></a>
<strong> double CoupSM::efvf(int idAbs) &nbsp;</strong> <br/>
   
<a name="anchor24"></a>
<strong> double CoupSM::vf2af2(int idAbs) &nbsp;</strong> <br/>
common quadratic combinations of the above couplings: 
<i>e_f^2, v_f^2, a_f^2, e_f * v_f, v_f^2 + a_f^2</i>. 
   
 
<a name="anchor25"></a>
<p/><strong> double CoupSM::VCKMgen(int genU, int genD) &nbsp;</strong> <br/>
   
<a name="anchor26"></a>
<strong> double CoupSM::V2CKMgen(int genU, int genD) &nbsp;</strong> <br/>
the CKM mixing element,or the square of it, for 
up-type generation index <code>genU</code> 
(<i>1 = u, 2 = c, 3 = t, 4 = t'</i>) and 
down-type generation index <code>genD</code> 
(<i>1 = d, 2 = s, 3 = b, 4 = b'</i>). 
   
 
<a name="anchor27"></a>
<p/><strong> double CoupSM::VCKMid(int id1, int id2) &nbsp;</strong> <br/>
   
<a name="anchor28"></a>
<strong> double CoupSM::V2CKMid(int id1, int id2) &nbsp;</strong> <br/>
the CKM mixing element,or the square of it, for 
flavours <code>id1</code> and <code>id2</code>, both in the 
range from <i>-18</i> to <i>+18</i>. The sign is here not 
checked (so it can be used both for <i>u + dbar &rarr; W+</i> 
and <i>u &rarr; d + W+</i>, say), but impossible flavour combinations 
evaluate to zero. The neutrino sector is numbered by flavor 
eigenstates, so there is no mixing in the lepton-neutrino system. 
   
 
<a name="anchor29"></a>
<p/><strong> double CoupSM::V2CKMsum(int id) &nbsp;</strong> <br/>
the sum of squared CKM mixing element that a given flavour can couple to, 
excluding the top quark and fourth generation. Is close to unity 
for the first two generations. Returns unity for the lepton-neutrino 
sector. 
   
 
<a name="anchor30"></a>
<p/><strong> int CoupSM::V2CKMpick(int id) &nbsp;</strong> <br/>
picks a random CKM partner quark or lepton (with the same sign as 
<code>id</code>) according to the respective squared elements, again 
excluding the top quark and fourth generation from the list of 
possibilities. Unambiguous choice for the lepton-neutrino sector. 
   
 
<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "6")
{
$data = "StandardModel:alphaSnfmax = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.00729735")
{
$data = "StandardModel:alphaEM0 = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.00781751")
{
$data = "StandardModel:alphaEMmZ = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.2312")
{
$data = "StandardModel:sin2thetaW = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.2315")
{
$data = "StandardModel:sin2thetaWbar = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1.16637e-5")
{
$data = "StandardModel:GF = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.97383")
{
$data = "StandardModel:Vud = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0.2272")
{
$data = "StandardModel:Vus = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "0.00396")
{
$data = "StandardModel:Vub = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.2271")
{
$data = "StandardModel:Vcd = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "0.97296")
{
$data = "StandardModel:Vcs = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.04221")
{
$data = "StandardModel:Vcb = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.00814")
{
$data = "StandardModel:Vtd = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.04161")
{
$data = "StandardModel:Vts = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.9991")
{
$data = "StandardModel:Vtb = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
