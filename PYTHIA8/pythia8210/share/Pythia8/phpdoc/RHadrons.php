<html>
<head>
<title>R-hadrons</title>
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

<form method='post' action='RHadrons.php'>
 
<h2>R-hadrons</h2> 
 
When a coloured SUSY particle is longer-lived than typical 
hadronization scales, i.e. around c*tau > 1 fm, or equivalently 
width Gamma < 0.2 GeV, it will have time to hadronize into a colour 
singlet hadronic state, a R-hadron. Currently a set of such 
R-hadrons have been implemented for the case of a long-lived 
gluino, stop or sbottom. Needless to say, the normal case would be 
that only one of them will be long-lived enough to form R-hadrons. 
 
<p/> 
For simplicity all gluino-mesons are assumed to have light-flavour 
spin 1, since those are the lightest and favoured by spin-state 
counting. Further, all gluino-baryons are bookkept as having 
light-flavour spin 3/2, and flavours are listed in descending order. 
This is more for convenience of notation, however, since the normal 
baryon octet e.g. has no uuu = "p++" state. When a diquark is 
extracted, a mixture of spin 0 and spin 1 is allowed. Names and codes 
are essentially in agreement with the PDG conventions, e.g. 
<br/>1000993 <code>R0(~g g)</code> (or gluinoball) 
<br/>1009213 <code>R+(~g u dbar)</code> (or gluino-rho+) 
<br/>1092214 <code>R+(~g uud)</code> (or gluino-Delta+) 
<br/>For internal bookkeeping of momenta, the code 1009002, 
<code>Rtemp(~g q)</code>, is used to denote the intermediate 
state formed when only one of the two string pieces attached to 
the gluino has broken. 
 
<p/> 
For the stop- and sbottom-hadrons the spin counting is simpler, 
since it is entirely given by the constituent quark or diquark spin. 
Again names and codes follow PDG conventions, e.g. 
<br/>1000612 <code>R+(~t dbar)</code> 
<br/>1006211 <code>R+(~t ud0)</code> 
 
<p/> 
The spin and electromagnetic charge of the new particle plays only 
a minor role in the hadronization process, that can be neglected 
to first approximation. Therefore it is possible to use the same 
R-hadrons framework instead for other BSM scenarios with long-lived 
coloured particles, e.g. with massive extra-dimensions copies 
of gluons and quarks, or with leptoquarks. This can be regulated by 
the switches below. Note that the codes and names of the R-hadrons 
is not changed when the heavy particle involved is switched, for 
reasons of administrative simplicity. R-hadron mass spectra and 
other relevant particle data is automatically updated to reflect 
the change, however. 
 
<br/><br/><strong>RHadrons:allow</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allows the gluino, stop and sbottom to hadronize if their respective 
widths are below the limit <code>RHadrons:maxWidth</code>. 
   
 
<br/><br/><table><tr><td><strong>RHadrons:maxWidth </td><td></td><td> <input type="text" name="2" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The maximum width of the gluino for which it is possible to form 
R-hadrons, provided that <code>RHadrons:allow</code> is on. 
   
 
<p/><code>mode&nbsp; </code><strong> RHadrons:idGluino &nbsp;</strong> 
 (<code>default = <strong>1000021</strong></code>)<br/>
The gluino identity code. For other scenarios than SUSY this code 
could be changed to represent another long-lived uncharged colour 
octet particle, that then would be treated in the same spirit. 
Could be set to 0 to forbid any gluino R-hadron formation even when 
the above two criteria, <code>RHadrons:allow</code> 
and <code>RHadrons:maxWidth</code>, are met. 
   
 
<p/><code>mode&nbsp; </code><strong> RHadrons:idStop &nbsp;</strong> 
 (<code>default = <strong>1000006</strong></code>)<br/>
The lightest stop identity code. For other scenarios than SUSY this 
code could be changed to represent another long-lived charge 2/3 
colour triplet particle, that then would be treated in the same 
spirit. As above it could be set to 0 to forbid any stop R-hadron 
formation. 
   
 
<p/><code>mode&nbsp; </code><strong> RHadrons:idSbottom &nbsp;</strong> 
 (<code>default = <strong>1000005</strong></code>)<br/>
The lightest sbottom identity code. For other scenarios than SUSY this 
code could be changed to represent another long-lived charge -1/3 
colour triplet particle, that then would be treated in the same 
spirit. As above it could be set to 0 to forbid any sbottom R-hadron 
formation. 
   
 
<br/><br/><strong>RHadrons:allowDecay</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allows the R-hadrons to decay or not. If the gluino/stop/sbottom is 
stable or too long-lived to decay inside the detector this switch 
has no real function, since then no decays will be performed anyway. 
If the sparticle is so short-lived that it decays before reaching 
the beam pipe then having the decay on is the logical choice. 
So the interesting region is when the decays happens after the 
R-hadron has passed through part of the detector, and changed its 
momentum and quite possibly its flavour content before it is to 
decay. Then normal decays should be switched off, and the R-hadron 
tracked through matter by a program like GEANT 
[<a href="Bibliography.php" target="page">Kra04,Mac07</a>]. After that, the new R-hadron info can be 
overwritten into the event record and the 
<code>Pythia::forceRHadronDecay()</code> method can be called 
to force this modified R-hadron to decay. 
   
 
<br/><br/><strong>RHadrons:setMasses</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use simple mass formulae to construct all available R-hadron masses 
based on the currently initialized gluino/squark masses and the 
constituent masses of the other partons in the hadron. If you switch 
this off, it is your responsibility to set each of the R-hadron masses 
on your own, and set them in an internally consistent way. If you 
mess up on this you may generate accordingly crazy results. 
Specifically, it is to be assumed that none of the R-hadrons has a 
mass below its constituent sparticle, i.e. that the light degrees 
of freedom and the additional confinement gluon field gives a net 
positive contribution to the R-hadron mass. 
   
 
<br/><br/><table><tr><td><strong>RHadrons:probGluinoball </td><td></td><td> <input type="text" name="5" value="0.1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
The fraction of produced gluino R-hadrons that are contain a "valence" 
gluon, with the rest containing a meson or baryon quark flavour content. 
   
 
<br/><br/><table><tr><td><strong>RHadrons:mOffsetCloud </td><td></td><td> <input type="text" name="6" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Extra mass (in GeV) added to each of the one or two extra constituent 
masses in an R-hadron, to calculate the mass of a R-hadron. The same 
offset is also used when the R-hadron momentum and mass is split 
between the squark or gluino and the one or two light (di)quarks, 
one for a squark and two for a gluino. Thus once or twice this amount 
represents a part of the nominal squark or gluino mass that will not 
decay weakly, since it is taken to correspond to the cloud of gluons 
that surround the squark or gluino. 
   
 
<br/><br/><table><tr><td><strong>RHadrons:mCollapse </td><td></td><td> <input type="text" name="7" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
A colour singlet system with an invariant mass less than this amount, 
above the R-hadron mass with the given flavour content, is assumed to 
collapse to this single R-hadron, whereas a full fragmentation handling 
is applied above this mass. 
   
 
<br/><br/><table><tr><td><strong>RHadrons:diquarkSpin1 </td><td></td><td> <input type="text" name="8" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
Probability that a diquark extracted from the flavour code of a gluino 
R-hadron should be assigned spin 1, with the rest being spin 0. Does 
not apply for two identical quarks, where spin 1 is only possibility. 
Note that gluino R-hadron codes for simplicity are assigned as if spin 
is 1 always, and so give no guidance. For stop and sbottom the diquark 
spin is preserved in the particle code, so there is no corresponding 
issue. 
   
 
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
$data = "RHadrons:allow = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0.2")
{
$data = "RHadrons:maxWidth = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "RHadrons:allowDecay = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "RHadrons:setMasses = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0.1")
{
$data = "RHadrons:probGluinoball = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.2")
{
$data = "RHadrons:mOffsetCloud = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1.0")
{
$data = "RHadrons:mCollapse = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0.5")
{
$data = "RHadrons:diquarkSpin1 = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2015 Torbjorn Sjostrand --> 
