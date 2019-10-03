function toggleDisplay(node,name,forceHide)
{
    if (node.style.display == "none" && forceHide == false) {
	node.style.display = "table-row";
	return true;
    }
    else {
	node.style.display = "none";
	toggle(name, true);
    }
    return false;
}
function toggleVis(node,name,forceHide)
{
    if (node.style.visibility == "collapse" && forceHide == false) {
	node.style.visibility = "visible";
	return true;
    }
    else {
	node.style.visibility = "collapse";
	toggle(name, true);
    }
    return false;
}

function hideAll()
{
    var i = 1;
    while (true) { 
	var nam = 'r' + i;
	var row = document.getElementById(nam);
	if (row == null) break;
	// console.log('Closing ' + nam);
	toggle(nam, true);
	i++;
    }
}
/*function expandAll()
  {
  hideAll();
  var i = 1;
  while (true) { 
  var nam = 'r' + i;
  var row = document.getElementById(nam);
  if (row == null) break;
  // console.log('Closing ' + nam);
  toggle(nam, false);
  i++;
  }
  }*/
    
function expandAll()
{
    hideAll();
    var max = document.getElementById("maxExpand").value;
    console.log("Got maximum level to be " + max);
    expandRecursive(null, 1, max);
}

function changeDepth(lvl)
{
    document.getElementById("currentMax").value = lvl;
}
function expandRecursive(parent,lvl,max)
{
    if (lvl > max) return;
    
    console.log('Expanding recursively @ lvl=' + lvl + '/' + max + ' ' + parent);

    var i = 1;
    var p = lvl == 1 ? "r" : parent + ".";
    console.log("Prefix=" + p);
 
    while (true) { 
	var nam = p + i;
	var row = document.getElementById(nam);
	if (row == null) break;
	console.log("Toggle " + nam);
	toggle(nam, false);
	i++;
	
	if (lvl+1 > max) continue;
	console.log("Recursing into " + nam);
	expandRecursive(nam,lvl+1,max);
    }
}
    
function toggle(node,forceHide)
{
    var i = 1;
    var c = 0;
    // console.log('toggle ' + node);
    while (true) {
	var sub = node + "." + i;
	var row = document.getElementById(sub);
	// console.log('got element ' + sub + ': ' + row);
	if (row) {
	    if (toggleDisplay(row, sub, forceHide)) c++;
	}
	else 
	    break;
	i++;
    } // while 
    var arn = node.replace('r','a');
    var arr = document.getElementById(arn);
    // console.log('Got arrow ' + arn + ': ' + arr + ' c=' + c);
    if (arr) {
	// console.log("inner HTML of " + arn + ': ' + arr.innerHTML);
	if (arr.innerHTML != "&nbsp;") {
	    if (c > 0) arr.innerHTML = "&blacktriangledown;";
	    else       arr.innerHTML = "&blacktriangleright;"
	}
    }   
}

function showSub(sub) { 
    console.log("Will show " + sub + "/index.html");
    var frm = document.getElementById("frame");
    if (frm == null) {
	console.log("No frame, will show in parent")
	location.assign(sub + "/index.html");
	return;
    }

    console.log("Got frame, will show in that")
    var ifrm = document.getElementById("iframe");
    ifrm.src = sub + "/index.html";
    
    frm.style.visibility = "visible";
}
function closeDisplay() { 
    var frm = document.getElementById("frame");
    if (frm == null) {
	console.log("No frame to close");
	return;
    }

    console.log("Closing frame");
    frm.style.visibility = "hidden";
}

function browseRootFile(filename) {
    var jsRoot = document.location.origin + "/jsRoot/";
    var path   = document.location.origin + document.location.pathname;
    path.replace(/index.html/,"");
    if (path[path.lenght-1] != '/') path += "/";
    path += filename;

    window.open(jsRoot + "?url=" + path, "_blank", 
		"location=no,menubar=no,status=no,titlebar=no");
}
	
//
// EOF
//
