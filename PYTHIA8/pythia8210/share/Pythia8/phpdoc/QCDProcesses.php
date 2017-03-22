<html>
<head>
<title>QCD Processes</title>
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

<form method='post' action='QCDProcesses.php'>
 
<h2>QCD Processes</h2> 
 
This section is subdivided into soft and hard QCD processes, with 
open charm and bottom production set aside as a special part of the 
latter, and three-jet topologies as a special subset. Kindly note 
that there is a considerable amount of overlap between the soft and 
hard process classes, so that you are likely to double-count 
if you include both in a run. 
 
<h3>Soft QCD processes</h3> 
 
As a rule, the processes in this class should not be mixed with 
the simulation of other processes. All by themselves, they are 
intended to represent the total cross section of hadron collisions, 
with the exception of the "rare processes" that one wishes to study 
separately. In particular, jet physics at all scales occurs as part 
of the minimum-bias description. 
 
<p/> 
We here use the "minimum bias" expression as a shorthand for 
inelastic, nondiffractive events. Strictly speaking, "minimum bias" 
represents an experimental procedure of accepting "everything", with 
some non-universal cuts to exclude elastic and diffractive topologies. 
In practice, the experimental minimum-bias sample may then contain 
some contamination of what is in PYTHIA classified as diffractive, 
especially (high-mass) double diffractive. 
 
<p/> 
Some options to modify these cross sections are found on the 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>Total Cross Sections</a> page. 
 
<br/><br/><strong>SoftQCD:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of all soft QCD processes, 
as listed separately in the following. 
   
 
<br/><br/><strong>SoftQCD:nonDiffractive</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The inelastic nondiffrative part of the total cross section, i.e. 
what would often be called the "minimum-bias component". 
The formalism is based on an <?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?> 
eikonalized description</a> of all the hard QCD processes, so 
includes them in combination with low-<i>pT</i> events. 
Code 101.<br/> 
Since the current description is handled by the multiparton-interactions 
machinery as part of the parton-level processing, no hard process at 
all is defined at the process-level part of the event generation. 
Fortunately, in this case a special 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>codeSub()</a></code> 
method provides information on the first, i.e. hardest, subprocess 
selected by the multiparton-interactions machinery. 
<br/><b>Note</b>: this event class is almost equivalent to the 
minimum-bias component of the total cross section. "Minimum-bias" 
usually refers to the experimental selection procedure, however, 
while "(inelastic) non-diffractive" better relates to the way events 
are generated in the program code. (Although also what separates 
diffractive from nondiffractive physics can be a matter of definition, 
especially once colour reconnection is to be modelled.) 
   
 
<br/><br/><strong>SoftQCD:elastic</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Elastic scattering <i>A B &rarr; A B</i>. 
Code 102. It is possible to include <?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?> 
Coulomb corrections</a>, but by default this is off. 
   
 
<br/><br/><strong>SoftQCD:singleDiffractive</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Single diffractive scattering <i>A B &rarr; X B</i> and 
<i>A B &rarr; A X</i>. See page on <?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?> 
Diffraction</a> for details. Codes 103 and 104. 
   
 
<br/><br/><strong>SoftQCD:doubleDiffractive</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Double diffractive scattering <i>A B &rarr; X_1 X_2</i>. 
See page on <?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?>Diffraction</a> 
for details. Code 105. 
   
 
<br/><br/><strong>SoftQCD:centralDiffractive</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Central diffractive scattering <i>A B &rarr; A X B</i> 
(a.k.a. double-Pomeron exchange, DPE). See pages on 
<?php $filepath = $_GET["filepath"];
echo "<a href='TotalCrossSections.php?filepath=".$filepath."' target='page'>";?>Total Cross Sections</a> 
and on <?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?>Diffraction</a> for details. 
In particular note the <code>SigmaTotal:zeroAXB</code> flag, 
which is on in most tunes, meaning no central diffraction, and 
that therefore would need to be reset to off after the selection 
of a tune (even the default one) to get central diffraction. 
Code 106. 
   
 
<br/><br/><strong>SoftQCD:inelastic</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
All of the above processes, except for elastic. Codes 101, 
103, 104, 105 and 106. 
   
 
<h3>Hard QCD processes</h3> 
 
This group contains the processes for QCD jet production above 
some minimum <i>pT</i> threshold. The <i>pT_min</i> cut cannot be put 
too low, or else unreasonably large jet cross sections will be obtained. 
This is because the divergent perturbative QCD cross section is used 
in this process group, without any regularization modifications. 
An eikonalized description, intended to be valid at all <i>pT</i>, 
is instead included as part of the multiparton-interactions framework, 
specifically in <code>SoftQCD:nonDiffractive</code> above. 
<br/><b>Warning 1</b>: you <b>must</b> remember to set the 
<code>PhaseSpace:pTHatMin</code> value if you use any of these 
processes; there is no sensible default. 
<br/><b>Warning 2</b>: you <b>must not</b> mix processes from the 
<code>SoftQCD</code> and <code>HardQCD</code> process groups, since 
this is likely to lead to double-counting. 
 
<br/><br/><strong>HardQCD:all</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of all hard QCD processes, 
as listed separately in the following. 
   
 
<br/><br/><strong>HardQCD:gg2gg</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; g g</i>. 
Code 111. 
   
 
<br/><br/><strong>HardQCD:gg2qqbar</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; q qbar</i>, where <i>q</i> by default 
is a light quark (<i>u, d, s</i>) (see below). 
Code 112. 
   
 
<br/><br/><strong>HardQCD:qg2qg</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g &rarr; q g</i> and <i>qbar g &rarr; qbar g</i>. 
Code 113. 
   
 
<br/><br/><strong>HardQCD:qq2qq</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q' &rarr; q q'</i>, <i>q qbar' &rarr; q qbar'</i>, 
<i>qbar qbar' &rarr; qbar qbar'</i>, where <i>q'</i> and <i>q</i> 
may agree, but the outgoing flavours equals the incoming ones 
Code 114. 
   
 
<br/><br/><strong>HardQCD:qqbar2gg</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; g g</i>. 
Code 115. 
   
 
<br/><br/><strong>HardQCD:qqbar2qqbarNew</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; q' qbar'</i>, where <i>q'</i> 
by default is a light quark (<i>u, d, s</i>) (see below). 
Code 116. 
   
 
<br/><br/><table><tr><td><strong>HardQCD:nQuarkNew  </td><td></td><td> <input type="text" name="15" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed outgoing new quark flavours in the above 
<i>g g &rarr; q qbar</i> and <i>q qbar &rarr; q' qbar'</i> processes, 
where quarks are treated as massless in the matrix-element expressions 
(but correctly in the phase space). It is thus assumed that <i>c cbar</i> 
and <i>b bbar</i> are added separately with masses taken into account, 
using the processes below. A change to 4 would also include <i>c cbar</i> 
in the massless approximation, etc. In order to avoid double-counting 
the processes below should then not be used simultaneously. 
   
 
<h3>Hard QCD processes: heavy-flavour subset</h3> 
 
These processes form a natural part of the above class, but can 
also be generated separately. Formally the heavy-quark mass makes 
these matrix elements finite in the <i>pT &rarr; 0</i> limit, but at 
high energies one may still question the validity of the expressions 
at low <i>pT</i> values, like for the other hard-QCD processes. 
Also as above, an eikonalized description, intended to be valid at all 
<i>pT</i>, is included as part of the multiparton-interactions framework. 
<br/>Note that the processes below only represent the "tip of the iceberg" 
of charm and bottom production at high energies, where flavour excitation 
and shower branchings provide major additional sources. All these sources 
come together in the descriptions offered by 
<code>SoftQCD:nonDiffractive</code> and <code>HardQCD:all</code>. 
 
<br/><br/><strong>HardQCD:gg2ccbar</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; c cbar</i>. 
Code 121. 
   
 
<br/><br/><strong>HardQCD:qqbar2ccbar</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; c cbar</i>. 
Code 122. 
   
 
<br/><br/><strong>HardQCD:hardccbar</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Sum of the previous two event types. 
Codes 121 and 122. 
   
 
<br/><br/><strong>HardQCD:gg2bbbar</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; b bbar</i>. 
Code 123. 
   
 
<br/><br/><strong>HardQCD:qqbar2bbbar</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; b bbar</i>. 
Code 124. 
   
 
<br/><br/><strong>HardQCD:hardbbbar</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Sum of the previous two event types. 
Codes 123 and 124. 
   
 
<h3>Hard QCD three-parton processes</h3> 
 
Three-parton final states are generated by showers off two-parton 
processes. Topologies then cannot be specified beforehand, beyond 
what is provided by the two-parton hard process. For some checks 
it may be convenient to have access to the dedicated three-parton 
final states, which is what this set of processes allows. 
Cross sections have been taken from [<a href="Bibliography.php" target="page">Ber81</a>]. 
<br/>Note that the processes in this section are  <it>not</it> 
affected by the <code>HardQCD:all</code> switch. In fact, it would 
be double-counting to include both the <code>HardQCD:all</code> and 
the <code>HardQCD:3parton</code> processes in a run or study. 
<br/><b>Warning:</b> this section is still incomplete, e.g. the 
selection of colour flow is very simple, and so it should only 
be used with caution. 
 
<br/><br/><strong>HardQCD:3parton</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of all hard QCD processes with three 
partons in the final state, as listed separately in the following. 
   
 
<br/><br/><strong>HardQCD:gg2ggg</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; g g g</i>. 
Code 131. 
   
 
<br/><br/><strong>HardQCD:qqbar2ggg</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; g g g</i>. 
Code 132. 
   
 
<br/><br/><strong>HardQCD:qg2qgg</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g &rarr; q g g</i> and <i>qbar g &rarr; qbar g g</i>. 
Code 133. 
   
 
<br/><br/><strong>HardQCD:qq2qqgDiff</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q' &rarr; q q' g</i>, <i>q qbar' &rarr; q qbar' g</i>, 
and <i>qbar qbar' &rarr; qbar qbar' g</i>. 
Code 134. 
   
 
<br/><br/><strong>HardQCD:qq2qqgSame</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q q &rarr; q q g</i> and 
<i>qbar qbar &rarr; qbar qbar g</i> 
(<i>q qbar &rarr; q qbar g</i> scatterings are considered separately 
below, see <code>HardQCD:qqbar2qqbargSame</code>). 
Code 135. 
   
 
<br/><br/><strong>HardQCD:qqbar2qqbargDiff</strong>  <input type="radio" name="28" value="on"><strong>On</strong>
<input type="radio" name="28" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; q' qbar' g</i>, where <i>q'</i> 
by default is a light quark (<i>u, d, s</i>) 
(see <code>HardQCD:nQuarkNew</code> above). 
Code 136. 
   
 
<br/><br/><strong>HardQCD:qqbar2qqbargSame</strong>  <input type="radio" name="29" value="on"><strong>On</strong>
<input type="radio" name="29" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar &rarr; q qbar g</i>. 
Code 137. 
   
 
<br/><br/><strong>HardQCD:gg2qqbarg</strong>  <input type="radio" name="30" value="on"><strong>On</strong>
<input type="radio" name="30" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g &rarr; q qbar g</i>, where <i>q</i> by 
default is a light quark (<i>u, d, s</i>) 
(see <code>HardQCD:nQuarkNew</code> above). 
Code 138. 
   
 
<br/><br/><strong>HardQCD:qg2qqqbarDiff</strong>  <input type="radio" name="31" value="on"><strong>On</strong>
<input type="radio" name="31" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g &rarr; q q' qbar'</i> and 
<i>qbar g &rarr; qbar qbar' q'</i>, where <i>q'</i> 
by default is a light quark (<i>u, d, s</i>) 
(see <code>HardQCD:nQuarkNew</code> above). 
Code 139. 
   
 
<br/><br/><strong>HardQCD:qg2qqqbarSame</strong>  <input type="radio" name="32" value="on"><strong>On</strong>
<input type="radio" name="32" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g &rarr; q q qbar</i> and 
<i>qbar g &rarr; qbar qbar q</i>. 
Code 140. 
   
 
 
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

if($_POST["1"] != "off")
{
$data = "SoftQCD:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "SoftQCD:nonDiffractive = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "SoftQCD:elastic = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "SoftQCD:singleDiffractive = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "SoftQCD:doubleDiffractive = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "SoftQCD:centralDiffractive = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "SoftQCD:inelastic = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "HardQCD:all = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "HardQCD:gg2gg = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "HardQCD:gg2qqbar = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "HardQCD:qg2qg = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "HardQCD:qq2qq = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "HardQCD:qqbar2gg = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "HardQCD:qqbar2qqbarNew = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "3")
{
$data = "HardQCD:nQuarkNew = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "HardQCD:gg2ccbar = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "HardQCD:qqbar2ccbar = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "HardQCD:hardccbar = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "HardQCD:gg2bbbar = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "HardQCD:qqbar2bbbar = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "HardQCD:hardbbbar = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "HardQCD:3parton = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "HardQCD:gg2ggg = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "HardQCD:qqbar2ggg = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "HardQCD:qg2qgg = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "HardQCD:qq2qqgDiff = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "HardQCD:qq2qqgSame = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "off")
{
$data = "HardQCD:qqbar2qqbargDiff = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "off")
{
$data = "HardQCD:qqbar2qqbargSame = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "off")
{
$data = "HardQCD:gg2qqbarg = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "off")
{
$data = "HardQCD:qg2qqqbarDiff = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "off")
{
$data = "HardQCD:qg2qqqbarSame = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
