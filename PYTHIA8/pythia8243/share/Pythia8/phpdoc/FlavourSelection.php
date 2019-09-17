<html>
<head>
<title>Flavour Selection</title>
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

<form method='post' action='FlavourSelection.php'>
 
<h2>Flavour Selection</h2> 
<ol id="toc">
  <li><a href="#section0">Flavour Selection for Gaussian <ei>pT</ei> Distribution</a></li>
  <li><a href="#section1">Flavour Selection for Thermal <ei>pT</ei> Distribution</a></li>
</ol>

 
The <code>StringFlav</code> class handles the choice of a new flavour 
in the fragmentation process, and the production of a new hadron 
from a set of input flavours. It is mainly used by the string 
fragmentation machinery (including ministrings), but also e.g. 
in some particle decays and for some beam-remnant cases. The basic 
concepts are in agreement with [<a href="Bibliography.php#refAnd83" target="page">And83</a>]. An alternative 
"thermal model" is described further below. 
 
<br/><br/><hr/> 
<a name="section0"></a> 
<h3>Flavour Selection for Gaussian <i>pT</i> Distribution</h3> 
 
The relative production rates of different particle species is 
influenced by the parameters below. Some have only an impact on 
one specific quantity, but most directly or indirectly have 
consequences for many observables. Therefore the values to use have 
to be viewed in the context of a complete <?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>tune</a>. 
 
<h4>New flavours</h4> 
 
The main parameters of the selection of a new flavour are 
 
<br/><br/><table><tr><td><strong>StringFlav:probStoUD </td><td></td><td> <input type="text" name="1" value="0.217" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.217</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
the suppression of <i>s</i> quark production relative to ordinary 
<i>u</i> or <i>d</i> one. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:probQQtoQ </td><td></td><td> <input type="text" name="2" value="0.081" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.081</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
the suppression of diquark production relative to quark production, 
i.e. of baryon relative to meson production. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:probSQtoQQ </td><td></td><td> <input type="text" name="3" value="0.915" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.915</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
the suppression of strange diquark production relative to light 
diquark production, over and above the one already given by 
<code>probStoU</code>. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:probQQ1toQQ0 </td><td></td><td> <input type="text" name="4" value="0.0275" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0275</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
the suppression of spin 1 diquark production relative to spin 0 one, 
apart from the factor of 3 enhancement of spin 1 from counting the 
number of states. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:probQQ1toQQ0join </td><td></td><td> <input type="text" name="5" value="0.5,0.7,0.9,1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5,0.7,0.9,1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
when two already produced quarks are to be combined to a diquark, 
e.g. in the junction framework, these numbers give the suppression 
of spin 1 diquark production relative to spin 0 one, apart from the 
factor of 3 enhancement of spin 1 from counting the number of states. 
The four components give the suppression when the heaviest quark is 
<i>u/d</i>, <i>s</i>, <i>c</i> or <i>b</i>, respectively. 
These parameters are seldom used and currently not constrained by any 
data, so very much a guesswork. Character-string input of this vector 
should be as a comma-separated list, without any blanks. 
   
 
<h4>Standard-meson production</h4> 
 
The bulk of the particle production corresponds to the lowest-lying 
pseudoscalar and vector multiplets. Their production rates are 
determined by the parameters in this section. 
 
<p/> 
For a given set of flavours, produced according to the probabilities 
outlined above, the ratio of vector-to-pseudocalar meson production 
is described by the parameters below. 
The maximum allowed rate for each case has been set according to 
spin-counting rules, but we expect the real rates to be lower, 
especially for lighter mesons, owing to the vector-pseudoscalar 
mass splitting. 
 
<br/><br/><table><tr><td><strong>StringFlav:mesonUDvector </td><td></td><td> <input type="text" name="6" value="0.50" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.50</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative production ratio vector/pseudoscalar for light 
(<i>u</i>, <i>d</i>) mesons. 
   
<br/><br/><table><tr><td><strong>StringFlav:mesonSvector </td><td></td><td> <input type="text" name="7" value="0.55" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.55</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative production ratio vector/pseudoscalar for strange mesons. 
   
<br/><br/><table><tr><td><strong>StringFlav:mesonCvector </td><td></td><td> <input type="text" name="8" value="0.88" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.88</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative production ratio vector/pseudoscalar for charm mesons. 
   
<br/><br/><table><tr><td><strong>StringFlav:mesonBvector </td><td></td><td> <input type="text" name="9" value="2.20" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.20</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative production ratio vector/pseudoscalar for bottom mesons. 
   
 
<p/> 
Inside each light-quark meson nonet, an octet-singlet mixing angle 
describes the mixing of the two flavour-diagonal isoscalar = 0 states. 
(For terminology and details see [<a href="Bibliography.php#refYao06" target="page">Yao06</a>], chapter 14 on the 
quark model.) 
This angle is needed to specify the probability for such a <i>q qbar</i> 
state to project onto a specific meson. More transparent formulae are 
obtained by introducing the angle <i>alpha = theta + 54.7</i> degrees: 
<br/><i> 
   f  = (uubar + ddbar)/sqrt(2) * sin(alpha) + ssbar * cos(alpha)<br/> 
   f' = (uubar + ddbar)/sqrt(2) * cos(alpha) - ssbar * sin(alpha) 
</i><br/> 
 
<br/><br/><table><tr><td><strong>StringFlav:thetaPS </td><td></td><td> <input type="text" name="10" value="-15." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-15.</strong></code>; <code>minimum = -90.</code>; <code>maximum = 90.</code>)</td></tr></table>
gives the mixing angle <i>theta_PS</i> in the pseudoscalar meson sector 
(which is rather poorly determined), expressed in degrees. 
Here <i>f</i> is associated with <i>eta'</i> and <i>f'</i> with 
<i>eta</i>. (This standard but counterintuitive choice is fixed up 
in the code by replacing <i>alpha &rarr; 90^0 - alpha</i> so that 
<i>eta &harr; eta'</i>; relative signs do not matter since we are 
interested in probabilities only.) 
   
 
<br/><br/><table><tr><td><strong>StringFlav:thetaV </td><td></td><td> <input type="text" name="11" value="36." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>36.</strong></code>; <code>minimum = -90.</code>; <code>maximum = 90.</code>)</td></tr></table>
gives the mixing angle <i>theta_V</i> in the vector meson sector 
(which is somewhat better determined), expressed in degrees. 
Here <i>f</i> is associated with <i>omega</i> and <i>f'</i> 
with <i>phi</i>. 
   
 
<p/> 
Further, the simple model overestimates the production of <i>eta</i> 
and, in particular, <i>eta'</i> mesons, which can be rectified by 
 
<br/><br/><table><tr><td><strong>StringFlav:etaSup </td><td></td><td> <input type="text" name="12" value="0.60" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.60</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the additional suppression of <i>eta</i> production, multiplying the 
normal production probability. Thus 0 means no <i>eta</i> at all 
are produced, while 1 means full rate. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:etaPrimeSup </td><td></td><td> <input type="text" name="13" value="0.12" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.12</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the additional suppression of <i>eta'</i> production, multiplying the 
normal production probability. Thus 0 means no <i>eta'</i> at all 
are produced, while 1 means full rate. 
   
 
<h4>Excited-meson production</h4> 
 
Several excited mesons, ie. with radial or orbital excitations, have been 
observed at non-negligible production rates. Extrapolated to all states 
a fair fraction of all particle production might proceed through such 
states. There are big uncertainties, however, since these excited 
mesons in many cases are extremely poorly known. This also means that 
the modeling of their production and decay is very primitive, and 
even that the inclusion of the production of such states may lead to a 
degraded agreement with data. Currently the default is that all such 
production is switched off. 
 
<p/> 
Parameters are provided to switch them on. By demand, this machinery 
has been made more flexible than in the past. Therefore one parameter is 
provided for each combination of heaviest flavour 
(<i>u/d</i>, <i>s</i>, <i>c</i> or <i>b</i>) and 
multiplet produced. In each case the production rate is normalized to 
that of the lowest-lying pseudoscalar of the same flavour content, as for 
the vector-meson rates introduced above. The multiplets available are the 
four obtained for one unit of orbital angular momentum, in the 
nonrelativistic classification. Using <i>J</i> to denote the sum of 
quark spin <i>S</i> and orbital angular momentum <i>L</i>, i.e. what 
would normally be called the spin of the meson, one has: 
<ul> 
<li>a pseudovector multiplet with <i>L=1, S=0, J=1</i>;</li> 
<li>a scalar multiplet with <i>L=1, S=1, J=0</i>;</li> 
<li>a pseudovector multiplet with <i>L=1, S=1, J=1</i>;</li> 
<li>a tensor multiplet with <i>L=1, S=1, J=2</i>.</li> 
</ul> 
 
The maximum allowed rate for each case has been set according to 
spin-counting rules, but we expect the real rates to be significantly 
lower, owing to mass suppression. 
 
<br/><br/><table><tr><td><strong>StringFlav:mesonUDL1S0J1 </td><td></td><td> <input type="text" name="14" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=0,J=1)</i>/pseudoscalar 
for light (<i>u</i>, <i>d</i>) mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonUDL1S1J0 </td><td></td><td> <input type="text" name="15" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the relative scalar production ratio 
<i>(L=1,S=1,J=0)</i>/pseudoscalar 
for light (<i>u</i>, <i>d</i>) mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonUDL1S1J1 </td><td></td><td> <input type="text" name="16" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=1,J=1)</i>/pseudoscalar 
for light (<i>u</i>, <i>d</i>) mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonUDL1S1J2 </td><td></td><td> <input type="text" name="17" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 5.</code>)</td></tr></table>
the relative tensor production ratio 
<i>(L=1,S=1,J=2)</i>/pseudoscalar 
for light (<i>u</i>, <i>d</i>) mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonSL1S0J1 </td><td></td><td> <input type="text" name="18" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=0,J=1)</i>/pseudoscalar 
for strange mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonSL1S1J0 </td><td></td><td> <input type="text" name="19" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the relative scalar production ratio 
<i>(L=1,S=1,J=0)</i>/pseudoscalar 
for strange mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonSL1S1J1 </td><td></td><td> <input type="text" name="20" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=1,J=1)</i>/pseudoscalar 
for strange mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonSL1S1J2 </td><td></td><td> <input type="text" name="21" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 5.</code>)</td></tr></table>
the relative tensor production ratio 
<i>(L=1,S=1,J=2)</i>/pseudoscalar 
for strange mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonCL1S0J1 </td><td></td><td> <input type="text" name="22" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=0,J=1)</i>/pseudoscalar 
for charm mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonCL1S1J0 </td><td></td><td> <input type="text" name="23" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the relative scalar production ratio 
<i>(L=1,S=1,J=0)</i>/pseudoscalar 
for charm mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonCL1S1J1 </td><td></td><td> <input type="text" name="24" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=1,J=1)</i>/pseudoscalar 
for charm mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonCL1S1J2 </td><td></td><td> <input type="text" name="25" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 5.</code>)</td></tr></table>
the relative tensor production ratio 
<i>(L=1,S=1,J=2)</i>/pseudoscalar 
for charm mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonBL1S0J1 </td><td></td><td> <input type="text" name="26" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=0,J=1)</i>/pseudoscalar 
for bottom mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonBL1S1J0 </td><td></td><td> <input type="text" name="27" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the relative scalar production ratio 
<i>(L=1,S=1,J=0)</i>/pseudoscalar 
for bottom mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonBL1S1J1 </td><td></td><td> <input type="text" name="28" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 3.</code>)</td></tr></table>
the relative pseudovector production ratio 
<i>(L=1,S=1,J=1)</i>/pseudoscalar 
for bottom mesons. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:mesonBL1S1J2 </td><td></td><td> <input type="text" name="29" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 5.</code>)</td></tr></table>
the relative tensor production ratio 
<i>(L=1,S=1,J=2)</i>/pseudoscalar 
for bottom mesons. 
   
 
<p/> 
In addition, an octet-singlet mixing angle is needed for each multiplet, 
as for the pseudoscalar and vector multiplets above. Only for the 
tensor multiplet does any determination exist; for the other multiplets 
default has been chose so that <i>ssbar</i> does not mix with the light 
quarks, and so that the <i>ssbar</i> state is the heavier of the two. 
 
<br/><br/><table><tr><td><strong>StringFlav:thetaL1S0J1 </td><td></td><td> <input type="text" name="30" value="35.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>35.3</strong></code>; <code>minimum = -90.</code>; <code>maximum = 90.</code>)</td></tr></table>
gives the mixing angle <i>theta</i> in the <i>(L=1,S=0,J=1)</i> 
pseudovector meson sector, expressed in degrees. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:thetaL1S1J0 </td><td></td><td> <input type="text" name="31" value="35.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>35.3</strong></code>; <code>minimum = -90.</code>; <code>maximum = 90.</code>)</td></tr></table>
gives the mixing angle <i>theta</i> in the <i>(L=1,S=1,J=0)</i> 
scalar meson sector, expressed in degrees. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:thetaL1S1J1 </td><td></td><td> <input type="text" name="32" value="35.3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>35.3</strong></code>; <code>minimum = -90.</code>; <code>maximum = 90.</code>)</td></tr></table>
gives the mixing angle <i>theta</i> in the <i>(L=1,S=1,J=1)</i> 
pseudovector meson sector, expressed in degrees. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:thetaL1S1J2 </td><td></td><td> <input type="text" name="33" value="28.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>28.0</strong></code>; <code>minimum = -90.</code>; <code>maximum = 90.</code>)</td></tr></table>
gives the mixing angle <i>theta</i> in the <i>(L=1,S=1,J=2)</i> 
tensor meson sector, expressed in degrees. 
   
 
<h4>Baryon production</h4> 
 
The relative rate of baryon production is mainly given by the quark 
and diquark production parameters above, plus SU(6) Clebsch-Gordans. 
The one modifiable parameter related to these coefficients is 
 
<br/><br/><table><tr><td><strong>StringFlav:decupletSup </td><td></td><td> <input type="text" name="34" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
the suppression, relative to default SU(6) factors, of decuplet 
baryon production. Default corresponds to no suppression, while 0 
corresponds to no decuplet production at all. 
   
 
<p/> 
In addition, if popcorn production is allowed, wherein a set of mesons 
(<i>M</i>) may be produced in between the baryon (<i>B</i>) and 
the antibaryon (<i>Bbar</i>), a set of further parameters is introduced. 
Currently only the simplest scenario is implemented, wherein at most 
one intermediate meson may be produced. 
 
<br/><br/><table><tr><td><strong>StringFlav:popcornRate </td><td></td><td> <input type="text" name="35" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 2.0</code>)</td></tr></table>
gives the relative rates of <i>B Bbar</i> and <i>B M Bbar</i> 
production, roughly as 
<br/><i> 
Prob(B M Bbar) / (Prob(B Bbar) + Prob(B M Bbar)) = 
popcornRate / (0.5 + popcornRate) 
</i><br/> 
(the complete expression depends on all the quark and diquark production 
parameters and is therefore not so useful). 
   
 
<br/><br/><table><tr><td><strong>StringFlav:popcornSpair </td><td></td><td> <input type="text" name="36" value="0.9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.0</code>)</td></tr></table>
extra suppression for having an <i>s sbar</i> pair shared between 
the <i>B</i> and <i>Bbar</i> in a <i>B M Bbar</i> configuration. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:popcornSmeson </td><td></td><td> <input type="text" name="37" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.0</code>)</td></tr></table>
extra suppression for having a strange meson <i>M</i> in a 
<i>B M Bbar</i> configuration. 
   
 
<p/> 
Finally, there are some indications that leading-baryon production 
may be further suppressed. A proper description should probably be 
based on a suppression of early production times [<a href="Bibliography.php#refEde97" target="page">Ede97</a>], 
but we here only implement a simpler version where production near 
the end of a string, as defined by rank, is suppressed. The more 
detailed studies suggest that leading <i>c</i> and <i>b</i> baryon 
production will be less suppressed, so we leave it open to set 
light- and heavy-baryon suppression separately. 
 
<br/><br/><strong>StringFlav:suppressLeadingB</strong>  <input type="radio" name="38" value="on"><strong>On</strong>
<input type="radio" name="38" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Suppress leading-baryon production. 
<br/><code>option </code><strong> off</strong> : No suppression.   
<br/><code>option </code><strong> on</strong> : Suppress the production of a diquark in the string 
breaking closest to a quark end of a string, by either of the factors 
below. This suppresses the production of first-rank baryons by the same 
amount. Indirectly also the second-rank and, if popcorn production is 
switched on, third-rank (anti)baryon production is affected.    
   
 
<br/><br/><table><tr><td><strong>StringFlav:lightLeadingBSup </td><td></td><td> <input type="text" name="39" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.0</code>)</td></tr></table>
extra suppression of leading-baryon production for a light-quark 
jet, i.e. <i>d</i>, <i>u</i> or <i>s</i>, when 
<code>suppressLeadingB = on</code>. Thus 0 means no leading-baryon 
production at all, while 1 means full rate. 
   
 
<br/><br/><table><tr><td><strong>StringFlav:heavyLeadingBSup </td><td></td><td> <input type="text" name="40" value="0.9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.0</code>)</td></tr></table>
extra suppression of leading-baryon production for a heavy-quark 
jet, i.e. <i>c</i> or <i>b</i>, when 
<code>suppressLeadingB = on</code>. Thus 0 means no leading-baryon 
production at all, while 1 means full rate. 
   
 
<br/><br/><hr/> 
<a name="section1"></a> 
<h3>Flavour Selection for Thermal <i>pT</i> Distribution</h3> 
 
If the hadronic <i>pT</i> is generated according to the non-default 
thermal distribution, i.e. if <code>StringPT:thermalModel = on</code>, 
the choice of a new flavour in the fragmentation process, and the 
production of a new hadron from a set of input flavours, depends mainly on 
the hadron mass [<a href="Bibliography.php#refFis16" target="page">Fis16</a>]. For a given <i>pT</i> value the new 
flavour is chosen according to 
<br/><i> 
  exp( -mT_had/T) = exp( - sqrt( pT_had^2 + mT_had^2 )/T). 
</i><br/> 
Here <i>T</i> is primarily given by <code>StringPT:temperature</code>, 
but can be further modified in the context of closely packed strings, 
<code>StringPT:closePacking = on</code>. 
Additional factors are included from theory arguments, for instance 
the ratio of vector-to-pseudocalar meson production is set according 
to spin-counting rules. 
Note that the octet-singlet mixing angles in the light-quark meson 
nonets are taken from the parameters above. 
Currently popcorn production has not been implemented, i.e. a baryon 
and an antibaryon are nearest neighbours in the flavour fragmentation 
chain, and share the flavours of one diquark. 
In addition the following two factors are introduced to provide an 
improved description of the flavour composition, although not as good 
as obtained in the default Gaussian scenario, with its bigger selection 
of free parameters. 
 
<p/> 
<br/><br/><table><tr><td><strong>StringFlav:BtoMratio </td><td></td><td> <input type="text" name="41" value="0.357" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.357</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.0</code>)</td></tr></table>
Ratio of the relative rate of baryon to meson production, i.e. every 
baryon Clebsch-Gordan coefficient gets multiplied by this factor. 
   
 
<p/> 
<br/><br/><table><tr><td><strong>StringFlav:StrangeSuppression </td><td></td><td> <input type="text" name="42" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.01</code>; <code>maximum = 1.0</code>)</td></tr></table>
Extra suppression factor for strange quarks. Note that in case of more 
than one strange quark in the hadron the factor gets squared or tripled 
respectively. 
   
 
<p/> 
The following parameters are used to determine which hadrons to choose 
from. By default only the pseudoscalar and vector meson nonet (L=0) 
and baryons with u/d/s quarks are included. For an already-existing 
heavier flavour, say c or b, this corresponds to picking only u/d/s 
for the new quark(s). 
<br/><b>Note:</b> The computer time for selecting the flavour of new 
hadrons goes linearly with the number of hadrons included. Therefore 
we recommend sticking to the default options as heavier hadrons are 
produced less likely anyway. 
 
<a name="anchor1"></a>
<p/><code>mode&nbsp; </code><strong> StringFlav:nQuark &nbsp;</strong> 
 (<code>default = <strong>3</strong></code>; <code>minimum = 3</code>; <code>maximum = 5</code>)<br/>
Selects the newly produced quark flavours that may be included in hadrons. 
The default corresponds to only include u/d/s quarks. 
   
 
<br/><br/><strong>StringFlav:mesonNonetL1</strong>  <input type="radio" name="43" value="on"><strong>On</strong>
<input type="radio" name="43" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Switch on to include the pseudovector, scalar, pseudovector, and tensor 
nonet (L=1). 
   
 
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

if($_POST["1"] != "0.217")
{
$data = "StringFlav:probStoUD = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.081")
{
$data = "StringFlav:probQQtoQ = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.915")
{
$data = "StringFlav:probSQtoQQ = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.0275")
{
$data = "StringFlav:probQQ1toQQ0 = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.5,0.7,0.9,1.0")
{
$data = "StringFlav:probQQ1toQQ0join = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.50")
{
$data = "StringFlav:mesonUDvector = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.55")
{
$data = "StringFlav:mesonSvector = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0.88")
{
$data = "StringFlav:mesonCvector = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "2.20")
{
$data = "StringFlav:mesonBvector = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "-15.")
{
$data = "StringFlav:thetaPS = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "36.")
{
$data = "StringFlav:thetaV = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.60")
{
$data = "StringFlav:etaSup = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.12")
{
$data = "StringFlav:etaPrimeSup = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.0")
{
$data = "StringFlav:mesonUDL1S0J1 = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.0")
{
$data = "StringFlav:mesonUDL1S1J0 = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.0")
{
$data = "StringFlav:mesonUDL1S1J1 = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.0")
{
$data = "StringFlav:mesonUDL1S1J2 = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.0")
{
$data = "StringFlav:mesonSL1S0J1 = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.0")
{
$data = "StringFlav:mesonSL1S1J0 = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.0")
{
$data = "StringFlav:mesonSL1S1J1 = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.0")
{
$data = "StringFlav:mesonSL1S1J2 = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "0.0")
{
$data = "StringFlav:mesonCL1S0J1 = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0.0")
{
$data = "StringFlav:mesonCL1S1J0 = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.0")
{
$data = "StringFlav:mesonCL1S1J1 = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "0.0")
{
$data = "StringFlav:mesonCL1S1J2 = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "0.0")
{
$data = "StringFlav:mesonBL1S0J1 = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "0.0")
{
$data = "StringFlav:mesonBL1S1J0 = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "0.0")
{
$data = "StringFlav:mesonBL1S1J1 = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "0.0")
{
$data = "StringFlav:mesonBL1S1J2 = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "35.3")
{
$data = "StringFlav:thetaL1S0J1 = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "35.3")
{
$data = "StringFlav:thetaL1S1J0 = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "35.3")
{
$data = "StringFlav:thetaL1S1J1 = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "28.0")
{
$data = "StringFlav:thetaL1S1J2 = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "1.0")
{
$data = "StringFlav:decupletSup = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "0.5")
{
$data = "StringFlav:popcornRate = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
if($_POST["36"] != "0.9")
{
$data = "StringFlav:popcornSpair = ".$_POST["36"]."\n";
fwrite($handle,$data);
}
if($_POST["37"] != "0.5")
{
$data = "StringFlav:popcornSmeson = ".$_POST["37"]."\n";
fwrite($handle,$data);
}
if($_POST["38"] != "off")
{
$data = "StringFlav:suppressLeadingB = ".$_POST["38"]."\n";
fwrite($handle,$data);
}
if($_POST["39"] != "0.5")
{
$data = "StringFlav:lightLeadingBSup = ".$_POST["39"]."\n";
fwrite($handle,$data);
}
if($_POST["40"] != "0.9")
{
$data = "StringFlav:heavyLeadingBSup = ".$_POST["40"]."\n";
fwrite($handle,$data);
}
if($_POST["41"] != "0.357")
{
$data = "StringFlav:BtoMratio = ".$_POST["41"]."\n";
fwrite($handle,$data);
}
if($_POST["42"] != "0.5")
{
$data = "StringFlav:StrangeSuppression = ".$_POST["42"]."\n";
fwrite($handle,$data);
}
if($_POST["43"] != "off")
{
$data = "StringFlav:mesonNonetL1 = ".$_POST["43"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
