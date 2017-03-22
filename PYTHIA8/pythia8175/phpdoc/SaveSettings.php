<html>
<head>
<title>Save Settings</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>

<?php
if($_GET["filename"] != "")
{
$_POST["filename"] = $_GET["filename"];
}
?>
<script language="JavaScript1.2">
<?php echo "setTimeout(\"top.index.location='Index.php?filename=".$_POST["filename"]."'\", 1)"; ?>
</script>
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

<form method='post' action='SaveSettings.php'>

<h2>Save Settings</h2>

The information on this webpage is only valid if you access the PHP
dynamic webpages via a web browser, and does not apply to the static
HTML equivalents. With PHP, all of the settings in the PYTHIA program 
are represented by radio buttons or fill-in boxes, that makes it easy 
for you to  construct a file with your desired changes. This file can 
then be read into PYTHIA by your main program to steer the whole run.

<h3>Basic instructions</h3>

The functionality of the PHP option is described in the following. 

<p/>
<table border="2" cellpadding="5"><td>
<?php

echo "<table border = 0><tr><td align=left valign=center>";

if($_POST["resetval"] == 1)
{
  //DELETE OLD FILE ON RESET
  unlink($_POST["rmfile"]);
}

if (!$_POST["filename"])
{
  echo "Please choose a (temporary) file name:</br>";
  echo "<form method='POST' action='SaveParameters.php'>";
  echo "<input type='text' name='filename'/>";
  echo "<input type='submit' value='Submit'/>";
  echo "</form>";
}
else if ($_POST["filename"] != "")
{
  $filename = $_POST["filename"];
  $filepath = "files/".$_POST["filename"];
  $filecheck = file_exists($filepath);
  
  if(($_POST["finish"] == 1) || ($_GET["returning"] == 1))
    {
      $filecheck = 0;
      $fileopen = 1;
    }
  
  if ($filecheck == 1)
    {
      echo "<font color='red'>File in use.. choose again:</br></font>";
      echo "<form method='POST' action='SaveParameters.php'>";
      echo "<input type='text' name='filename'/>";
      echo "<input type='submit' value='Submit'/>";
      echo "</form>";
    }
  else
    {
      if(!$_POST["finish"])
	{
	  echo "Filename chosen:<b>  "; 
          echo $_POST["filename"];
          echo "</br></b>";
	  echo "<table border = 0 valign = top><tr><td valign=\"top\"><form method='POST' action='SaveParameters.php'>";
	  echo "<input type='hidden' name='filename' value='".$filename."'/>";
	  echo "<input type='hidden' name='finish' value='1'/>";
	  echo "<input type='submit' value='Finish File'/>";
	  echo "</form></td>";
	  echo "<td valign=\"top\"><form method='POST' action='SaveParameters.php'>";
	  echo "<input type='hidden' name='filename' value=''/>";
	  echo "<input type='hidden' name='rmfile' value='".$filepath."'/>";
	  echo "<input type='hidden' name='resetval' value='1'/>";
	  echo "<input type='submit' value='RESET'/>";
	  echo "</form></td></tr></table>";
	}
      if($fileopen != 1)
	{
	  //CREATE & OPEN FILE & CHANGE PERMISSIONS TO 666
	  $handle = fopen($filepath, 'w');
	  chmod($filepath, 0666);
	  fclose($handle);

	}
      if ($_POST["finish"])
	{
	  echo "<b><li><a href='".$filepath."'>Right Click to Save Target: ".$filename."</a></b></br>";
	  echo "<form method='POST' action='SaveParameters.php'>";
	  echo "<input type='hidden' name='filename' value=''/>";
	  echo "<input type='hidden' name='rmfile' value='".$filepath."'/>";
	  echo "<input type='hidden' name='resetval' value='1'/>";
	  echo "<input type='submit' value='RESET'/>";
	  echo "</form>";
	}
    }
}

echo "</td><td></td></tr></table>";

?>
</td></table>

<ul>

<p/><li>
To begin with, you must specify a <b>(temporary) file name</b> in the 
box above. If the filename already exists on the server, you will be 
requested to pick a new name.</li>

<p/><li>
Once you have <b>Submit</b>ted your filename, you can browse through the
pages and make your selections. The values currently selected when you 
load the page are the default values.</li>

<p/><li>
When you have finished making your changes to a particular page, 
you <b>must</b> click on <b>Save Settings</b> at the <b>bottom</b> of 
the page. This will write the changes to your temporary file. If you make 
a mistake, just repeat the procedure for that category again.<br> 

<p/><li>
When you have finished all the changes you need, return to this page 
and click <b>Finish File</b>.</li> 

<p/><li>
You will then get up a link, that you are asked to <b>right-click</b> 
with your mouse (or equivalent).</li> 

<p/><li>
In the menu that appears, pick the option <b>Save Link As</b>
(or equivalent).</li> 

<p/><li>
You will now get up a file browser, for you to pick and <b>Save</b>
the location and file name (the latter by default the same as the 
temporary file name).</li> 

<p/><li>
At any time, if you click the <b>RESET</b> button, your temporary
file will be erased and you can start anew.</li> 

<p/><li>
Before you use a file, be sure to <b>check it visually</b> to confirm 
that you saved what you intended to. Minor corrections are easily made 
in a text editor.
</li> 

</ul>

<p/>
<h3>Supplementary notes</h3>


The documentation files exist in three versions.
<ol>

<p/><li> 
As a set of <code>.xml</code> files, in the <code>xmldoc/</code>
subdirectory. These are the master copies that no user ever should 
touch, but that are used to generate the variants below.</li>  

<p/><li> 
As a set of <code>.html</code> files, in the <code>htmldoc/</code>
subdirectory. You can open your own locally installed copy of the 
<code>Welcome.html</code> file in your web browser and thereafter 
navigate among all the pages. You can learn which parameters are free 
to be changed, but not change anything, except by brute-force 
cut-and-paste to a file of your own.</li> 

<p/><li>
As a set of <code>.php</code> files, in the <code>phpdoc/</code>
subdirectory. For these files to provide the functionality described 
above they have to accessed via a webserver. The one where you have 
your homepage should work fine. Alternatively you can use pages already 
available on another server.</li>

</ol>

<p/>
A few further comments about the operation of the PHP option:
<ul>

<p/><li>
To set up the PHP files on your webserver, you have to install the whole 
<code>phpdoc/</code> subdirectory there. In addition to the  
<code>.php</code> code this includes a few more files, plus a 
subdirectory named <code>files</code> where the temporary files 
are stored. This subdirectory must have public write access to work
(<code>chmod a+w files</code> if not).</li>

<p/><li>
The "temporary" files stored in <code>files</code> actually remain 
unless the RESET button is used. The good news is that this makes
it possible to recover a file that otherwise might be lost. The bad 
news is that the <code>files</code> directory may need to be cleaned 
up from time to time. (But typically the files are pretty small, so
this should not be a major problem.)</li>

<p/><li>
When you click the <b>Save Settings</b> button on the bottom of a page
all changed settings are written on the temporary file in the format
<pre>
name-of-flag/mode/parameter/word = value
</pre>  
with one variable per line. Thereafter all the settings on the page
are restored to their default values.</li>

<p/><li>
You can return to a page to do some further changes and save those.
If you change the same parameter twice, it is the last value that
counts. (Both values are stored in the file, with the more recent
lower down, and then PYTHIA does the changes sequentially.) However
remember that unchanged values are not stored, so if you want to 
restore some default value it may be simpler to edit the file 
afterwards.</li>   

<p/><li>
The changeable flags/modes/parameters/words are mainly in the  
"Setup Run Tasks" section of the index, but a few (less
frequently used ones) can also be found lower down, in the 
"Study Output" and "Link to Other Programs" pages. 

<p/><li>
It is not (yet) possible to modify particle data within the PHP-based
setup approach. This is a more difficult task, since e.g. the
modifications one may want to do in a decay table can be quite
interrelated.

</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->

