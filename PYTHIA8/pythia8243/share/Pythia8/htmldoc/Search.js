// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten, May 2017, based on search script by Eadaoin Ilten.

// Function to search the Pythia 8 HTML documentation.
$(function() {
    
    // Search result format.
    var resultFormat = _.template(
	'<a href="javascript:link(\'<%=result.link%>\')" ' +
	    'style="white-space:nowrap"><%=result.name%></a>');
    
    // Search input.
    var searchInput = $('#search-input');
    searchInput.keyup(searchIndex).on('search', searchIndex);
    
    // Search result.
    var searchResult = {
	$field: $('#search-input'), $list: $('#search-result'),
	length: 0, nohide: false, input: ''
    }
    
    // Load the search index.
    var load = document.createElement('script');
    load.src = 'Index.js';
    document.head.appendChild(load);

    // Function to search the index.
    function searchIndex() {

	// Check if the search input has been updated.
	var input = searchResult.$field.val();
	if (input == searchResult.input) return true;
	else searchResult.input = input;
	var results = [];

	// Perform the search if the input is not empty.
	if (input) {
	    regexp = new RegExp(
		_.map((input ? '' + input : '').split(/\s/), function(str) {
		    return str.replace(
			    /[\\^$.*+?|{[()]/g, "\\$&")}).join('.*'), "i");
	    results = _.filter(index, function(result) {
		return regexp.test(result.text)});
	}

	// Format the search results.
	searchResult.$list.html(_.map(results, function(result) {
	    return resultFormat({result:result})}).join('<br/>'));
	searchResult.length = results.length;
	if (!results.length) return false;
	return true;
    }
})
