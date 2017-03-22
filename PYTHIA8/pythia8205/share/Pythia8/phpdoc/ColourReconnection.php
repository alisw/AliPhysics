<html>
<head>
<title>Colour reconnection</title>
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

<form method='post' action='ColourReconnection.php'>
 
<h2>Colour Reconnection</h2> 
 
The colour flows in the separate subprocesses defined in the 
multiparton-interactions scenario are tied together via the assignment 
of colour flow in the beam remnant. This is not an unambiguous 
procedure, and currently two different methods are implemented. In the first 
model the colour flow is reconstructed by how a PS could have 
constructed the configuration. In the second model, the full QCD colour 
calculation is taken into account, however the dynamical effects are modeled 
loosely, only an overall saturation is taken into account. The idea is to 
later account for other dynamical effects through colour reconnections. 
 
<p/> 
A simple "minimal" procedure of colour flow only via the beam remnants 
does not result in a scenario in agreement with data, however, 
notably not a sufficiently steep rise of 
<i>&lt;pT&gt;(n_ch)</i>. The true origin of this behaviour and the 
correct mechanism to reproduce it remains one of the big unsolved issues 
at the borderline between perturbative and nonperturbative QCD. Since no final 
answer is known, several models are implemented. The different models also 
rely on the two different colour assignments in the beam remnant. There are 
two, somewhat motivated, models implemented: the original PYTHIA scheme and 
a new scheme that tries to incorporate more of the colour knowledge from QCD. 
 
<p/> 
The original PYTHIA scheme relies on the PS-like colour configuration of the 
beam remnant. This is combined with an additional step, wherein the gluons 
of a lower-<i>pT</i> MPI system are merged with the ones in a higher-pT MPI. 
A more detailed description of the merging can be found below. 
Relative to the other models it tests fewer reconnection possibilities, 
and therefore tends to be reasonably fast. 
 
<p/> 
The new scheme [<a href="Bibliography.php" target="page">Chr14a</a>]relies on the full QCD colour configuration 
in the beam remnant. This is followed up by a colour reconnection, where the 
potential string energy is minimized (ie. the <i>lambda</i> measure is 
minimized). The QCD colour rules are also incorporated in the colour 
reconnection, and determine the probability that a reconnection is allowed. 
The model also allows the creation of junction structures. 
 
<p/> 
In addition to the two models described above, a simple model is implemented, 
wherein gluons can be moved from one location to another so as to reduce the 
total string length. This is one out of a range of simple models developed 
to study potential colour reconnection effects e.g. on top mass 
[<a href="Bibliography.php" target="page">Arg14</a>], not from the point of view of having the most realistic 
description, but in order to probe the potential worst-case spread of 
predictions. All of these models are made available separately in 
<code>include/Pythia8Plugins/ColourReconnectionHooks.h</code>, with the 
setup illustrated in <code>examples/main29.cc</code>, but only the 
gluon-move one is sufficiently general and realistic that it has been 
included among the standard options here. 
 
<br/><br/><strong>ColourReconnection:reconnect</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow or not a system to be merged with another one. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:mode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Determine which model is used for colour reconnection. Beware that 
different <code>BeamRemnants:remnantMode</code> should be used for 
different reconnection schemes. 
<br/>
<input type="radio" name="2" value="0" checked="checked"><strong>0 </strong>:  The MPI-based original Pythia 8 scheme.  <br/>
<input type="radio" name="2" value="1"><strong>1 </strong>:  The new more QCD based scheme.  <br/>
<input type="radio" name="2" value="2"><strong>2 </strong>:  The new gluon-move model.  <br/>
 
<h3>The MPI-based scheme</h3> 
 
In this scheme partons are classified by which MPI system they belong to. 
The colour flow of two such systems can be fused, and if so the partons 
of the lower-<i>pT</i> system are added to the strings defined by the 
higher-<i>pT</i> system in such a way as to give the smallest total 
string length. The bulk of these lower-<i>pT</i> partons are gluons, 
and this is what the scheme is optimized to handle. 
 
<p/> 
In more detail, an MPI system with a scale <i>pT</i> of the hard 
interaction (normally <i>2 &rarr; 2</i>) can be merged with one of 
a harder scale with a probability that is 
<i>pT0_Rec^2 / (pT0_Rec^2 + pT^2)</i>, where <i>pT0_Rec</i> is 
<code>range</code> times <i>pT0</i>, the latter being the same 
energy-dependent dampening parameter as used for MPIs. 
Thus it is easy to merge a low-<i>pT</i> system with any other, 
but difficult to merge two high-<i>pT</i> ones with each other. 
 
<br/><br/><table><tr><td><strong>ColourReconnection:range </td><td></td><td> <input type="text" name="3" value="1.8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.8</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
The <code>range</code> parameter defined above. The higher this number is 
the more reconnections can occur. For values above unity the reconnection 
rate tends to saturate, since then most systems are already connected with 
each other. This is why 10 is set as an effective upper limit, beyond 
which it is not meaningful to let the parameter go. 
   
 
<p/> 
The reconnection procedure is applied iteratively. Thus first the 
reconnection probability <i>P = pT0_Rec^2 / (pT0_Rec^2 + pT^2)</i> of the 
lowest-<i>pT</i> system is found, and gives the probability for merger with 
the second-lowest one. If not merged, it is tested with the third-lowest 
one, and so on. For the <i>m</i>'th higher system the reconnection 
probability thus becomes <i>(1 - P)^(m-1) P</i>. That is, there is 
no explicit dependence on the higher <i>pT</i> scale, but implicitly 
there is via the survival probability of not already having been merged 
with a lower-<i>pT</i> system. Also note that the total reconnection 
probability for the lowest-<i>pT</i> system in an event with <i>n</i> 
systems becomes <i>1 - (1 - P)^(n-1)</i>. Once the fate of the 
lowest-<i>pT</i> system has been decided, the second-lowest is considered 
with respect to the ones above it, then the third-lowest, and so on. 
 
<p/> 
Once it has been decided which systems should be joined, the actual merging 
is carried out in the opposite direction. That is, first the hardest 
system is studied, and all colour dipoles in it are found (including to 
the beam remnants, as defined by the holes of the incoming partons). 
Next each softer system to be merged is studied in turn. Its gluons are, 
in decreasing <i>pT</i> order, inserted on the colour dipole <i>i,j</i> 
that gives the smallest <i>(p_g p_i)(p_g p_j)/(p_i p_j)</i>, i.e. 
minimizes the "disturbance" on the existing dipole, in terms of 
<i>pT^2</i> or <i>Lambda</i> measure (string length). The insertion 
of the gluon means that the old dipole is replaced by two new ones. 
Also the (rather few) quark-antiquark pairs that can be traced back to 
a gluon splitting are treated in close analogy with the gluon case. 
Quark lines that attach directly to the beam remnants cannot be merged 
but are left behind. 
 
<p/> 
The joining procedure can be viewed as a more sophisticated variant of 
the one introduced already in [<a href="Bibliography.php" target="page">Sjo87</a>]. Clearly it is ad hoc. 
It hopefully captures some elements of truth. The lower <i>pT</i> scale 
a system has the larger its spatial extent and therefore the larger its 
overlap with other systems. It could be argued that one should classify 
individual initial-state partons by <i>pT</i> rather than the system 
as a whole. However, for final-state radiation, a soft gluon radiated off 
a hard parton is actually produced at late times and therefore probably 
less likely to reconnect. In the balance, a classification by system 
<i>pT</i> scale appears sensible as a first try. 
 
<p/> 
Note that the reconnection is carried out before resonance decays are 
considered by default. Colour inside a resonance therefore is not 
reconnected. The 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='MasterSwitches.php?filepath=".$filepath."' target='page'>";?>PartonLevel:earlyResDec</a></code> 
can be switched on to perform resonance decays before colour reconnection, 
and then the partons from resonance decays could be affected. 
Ideally the time scales of resonance decays and of colour reconnection 
should be picked dynamically, but this is not yet the case. Notably 
the <i>W</i>, <i>Z</i> and <i>t</i> have intermediate decay time 
scales, somewhat but not much shorter than typical hadronization times, 
whereas the <i>H</i> is much more long-lived. 
 
<h3>The newer scheme</h3> 
 
<br/><br/><table><tr><td><strong>ColourReconnection:m0 </td><td></td><td> <input type="text" name="4" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 5.</code>)</td></tr></table>
This is the variable used in the lambda measure for the string length. 
See the different choices of lambda measure for exact formulaes. This variable 
is in the new model also used as a cut for forming pseudo particles that are 
not colour reconnected. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:nColours  </td><td></td><td> <input type="text" name="5" value="9" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>9</strong></code>; <code>minimum = 1</code>; <code>maximum = 30</code>)</td></tr></table>
The number of reconnection colours, this should not be confused with the 
standard number of QCD colours. Each string is given an integer number between 
0 and <code>nColours - 1</code>. Only strings with the same number are allowed 
to do a normal string reconnection. The default value provides 
the standard QCD probability that a triplet and an anti-triplet is in a 
singlet state. The probability for two strings to form a junction structure is 
in QCD given by the product of two triplets, which gives a probability of 1/3. 
Therefore the number of reconnection colours for junction formation is 
<code>iColours % 3</code>, where iColours refer to the integer of the string. 
The behaviour of junction formation therefore only changes slightly with this 
variable. 
   
 
<br/><br/><strong>ColourReconnection:sameNeighbourColours</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
In the normal colour reconnection two neighbouring strings are not allowed 
to have the same colour. Similar two strings orginating from a gluon split is 
not allowed to reconnect. The physics motivation for this is that it would 
require colour singlet gluons, and therefore for ordinary physics studies this 
should be turned off. But for testing of extreme scenarios (i.e. 1 colour), 
this variable needs to be turned on, since it is not possible to have 
different neighbouring colours. 
   
 
<br/><br/><strong>ColourReconnection:allowJunctions</strong>  <input type="radio" name="7" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="7" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This switch disables the formation of junctions in the colour reconnection. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:lambdaForm  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
This allows to switch between different options for what 
<ei>lambda</ei>-measure to use. 
The formula shown are how much each end of a dipole or junction contribute to 
the total <ei>lambda</ei>-measure. The energies are defined in respectively 
the dipole or junction restframe. 
<br/>
<input type="radio" name="8" value="0" checked="checked"><strong>0 </strong>:  <ei>lambda = ln (1 + sqrt(2) E/m0)</ei>  <br/>
<input type="radio" name="8" value="1"><strong>1 </strong>:  <ei>lambda = ln (1 + 2 E/m0)</ei>  <br/>
<input type="radio" name="8" value="2"><strong>2 </strong>:  <ei>lambda = ln (2 E/m0)</ei>  <br/>
 
<br/><br/><table><tr><td><strong>ColourReconnection:minimumGainJun </td><td></td><td> <input type="text" name="9" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = -100</code>; <code>maximum = 100</code>)</td></tr></table>
The minimum <i>lambda</i> has to decrease in order to create a junction 
antijunction pair. 
   
 
<br/><br/><strong>ColourReconnection:allowDoubleJunRem</strong>  <input type="radio" name="10" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="10" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
This parameter tells whether or not to allow a directly connected 
junction-antijunction pair to split into two strings. The lambda measure of 
the junction system is compared to that of the two possible string 
configurations. If the chosen configuration is the junction system, a q-qbar 
system is inserted between the junctions by removing some energy/momentum from 
the other legs. 
   
 
<h3>The gluon-move scheme</h3> 
 
This approach contains two steps, a first "move" one and an optional 
second "flip" one. Both are intended to reduce the total "string length" 
<i>lambda</i> measure of an event. For multiparton topologies the 
correct <i>lambda</i> measure can become quite cumbersome, so here it 
is approximated by the sum of <i>lambda</i> contributions from each pair 
of partons connected by a colour string piece. For two partons <i>i</i> 
and <i>j</i> with invariant mass <i>m_ij</i> this contribution 
is defined as <i>lambda_ij = ln(1 + m^2_ij / m2Lambda)</i>. 
The 1 is added ad hoc to avoid problems in the <i>m_ij &rarr; 0</i> 
limit, problems which mainly comes from the approximate treatment, 
and <i>m2Lambda</i> is a parameter set below. 
 
<p/> 
In the move step all final gluons are identified, alternatively only a 
fraction <i>fracGluon</i> of them, and also all colour-connected 
parton pairs. For each gluon and each pair it is calculated how the total 
<i>lambda</i> would shift if the gluon would be removed from its current 
location and inserted in between the pair. The gluon move that gives the 
largest negative shift, if any, is then carried out. Alternatively, only 
shifts more negative than <i>dLambdaCut</i> are considered. The procedure 
is iterated so long as allowed moves can be found. 
 
<p/> 
There is some fine print. If a colour singlet subsystem consists of two 
gluons only then it is not allowed to move any of them, since that then 
would result in in a colour singlet gluon. Also, at most as many moves 
are made as there are gluons, which normally should be enough. A specific 
gluon may be moved more than once, however. Finally, a gluon directly 
connected to a junction cannot be moved, and also no gluon can be inserted 
between it and the junction. This is entirely for practical reasons, but 
should not be a problem, since junctions are rare in this model. 
 
<p/> 
The gluon-move steps will not break the connection between string endpoints, 
in the sense that a quark and an antiquark that are colour-connected via 
a number of gluons will remain so, only the number and identity of the 
intermediate gluons may change. Such a scenario may be too restrictive. 
Therefore an optional second flip step is introduced. In it all such 
colour chains are identified, omitting closed gluon loops. The lambda 
change is defined by what happens if the two colour lines are crossed 
somewhere, e.g. such that two systems <i>q1 - g1 - qbar1</i> and 
<i>q2 - g2 - qbar2</i> are flipped to <i>q1 - g1 - g2 - qbar2</i> 
and <i>q2 - qbar1</i>. The flip that gives the largest <i>lambda</i> 
reduction is carried out, again with <i>dLambdaCut</i> offering a 
possibility to restrict the options. As with the move step, the procedure 
is repeated so long as it is allowed. An important restriction is 
imposed, however, that a given system is only allowed to flip once, 
and not with itself. The practical reason is that repeated flips could 
split off closed gluon loops quite easily, which tends to result in 
bad agreement with data. 
 
<p/> 
As an option, singlet subsystems containing a junction may or may not 
be allowed to take part in the flip step. Since the number of junction 
systems is limited in this model the differences are not so important. 
 
<br/><br/><table><tr><td><strong>ColourReconnection:m2Lambda </td><td></td><td> <input type="text" name="11" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 16.</code>)</td></tr></table>
The <i>m2Lambda</i> parameter used in the definition of the approximate 
<i>lambda</i> measure above. It represents an approximate hadronic 
mass-square scale, cf. <i>m0</i> in the previous model. Its value is 
uncertain up to factors of 2, but the <i>lambda</i> change induced by 
a potential move or flip is rather insensitive to the precise value, 
owing to large cancellations. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:fracGluon </td><td></td><td> <input type="text" name="12" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
The probability that a given gluon will be considered for being moved. 
It thus gives the average fraction of gluons being considered. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:dLambdaCut </td><td></td><td> <input type="text" name="13" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
Restrict gluon moves and colour flips to those that reduce <i>lambda</i> 
by more than this amount. The larger this number, the fewer moves and flips 
will be performed, but those that remain are the ones most likely to produce 
large effects. 
   
 
<br/><br/><table><tr><td><strong>ColourReconnection:flipMode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Performing the flip step or not. 
<br/>
<input type="radio" name="14" value="0" checked="checked"><strong>0 </strong>:  No flip handling.  <br/>
<input type="radio" name="14" value="1"><strong>1 </strong>:  Allow flips, but not for strings in junction topologies.  <br/>
<input type="radio" name="14" value="2"><strong>2 </strong>:  Allow flips, including for strings in junction topologies.  <br/>
 
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

if($_POST["1"] != "on")
{
$data = "ColourReconnection:reconnect = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0")
{
$data = "ColourReconnection:mode = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1.8")
{
$data = "ColourReconnection:range = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0.5")
{
$data = "ColourReconnection:m0 = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "9")
{
$data = "ColourReconnection:nColours = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "ColourReconnection:sameNeighbourColours = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "on")
{
$data = "ColourReconnection:allowJunctions = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0")
{
$data = "ColourReconnection:lambdaForm = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1")
{
$data = "ColourReconnection:minimumGainJun = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "on")
{
$data = "ColourReconnection:allowDoubleJunRem = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1.")
{
$data = "ColourReconnection:m2Lambda = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1.")
{
$data = "ColourReconnection:fracGluon = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.")
{
$data = "ColourReconnection:dLambdaCut = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0")
{
$data = "ColourReconnection:flipMode = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
