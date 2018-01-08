/* JS file to run the results page of MIPhy.
-- I modified jsPhyloSVG, encapsulating the entire thing in a load function. This allows it to be reloaded without keeping data from a previous tree.
-- Majority of processing time is in jsPhyloSVG (84%), 8% from searchHighlights, 4% drawing species markers.
  -- Using raphael 2.2 is ~20% slower than raphael 1.4.
-- RANGER-DTL does similar thing minus ILS, but extremely fast. See if there's anything I can use from their algorithm.

   * Major to-do *
-- Add loading text or image, for when page is first loading, and when reclustering.
-- If one cluster has any more than 50% of the sequences, that cluster shape is way too big.
  -- Should add 1 concave node at the tree centre somehow.
-- If a web session has timed out, perhaps give the user an option to set up a new instance, keeping the current settings of the page.
-- Was trying the web version; had an instance of vert cyps open for a while, then the session closed. don't know why. I was searching and zooming mostly.
-- If a user's computer is quite laggy, the program might exit before the html loads. Provide option to set this wait.
-- Some mechanism of setting parameters for each species independently.
-- Method to select a node, and open a new tab of that subtree. Should be using the same data as the original (changing parameters would recluster the entire tree), but only displaying that sub-tree.
-- Internet Explorer and MS Edge hang when displaying results page. If I can't solve it, should at least inform users.
-- MIPhy captures and displays all warnings and errors from the browser. Can I at least title these as browser errors?
  -- Seems to happen frequently if I run it with the WiFi turned off.

   * Technical to-do *
-- The more times a tree is re-drawn (either by changing weights or by changing colours), the more extreme the action of the +/- zoom buttons. No idea why, but fix.
-- If there are many species, the legend is tall, and the tree should hug the top. If their names are long, it is wide, and the tree should hug the left side of the screen.
  -- The program should make this decision based on the measured lengths (or maybe a drawing option).
-- Mouseover on the search highlight shape should be the same as mouseover on the cluster (ideally), or perhaps just that sequence.
-- Change order of numbered dots in cluster tooltop to match order of species tree legend.
-- Change order of event counts in cluster tooltip to: D, L, I.
-- If you have a sequence highlighted from the search, and redraw from the "display options", the search highlight disappears. You have to click the magnifying glass again to show up.
-- Add more colours to the list.
-- If the species names are long, the summary box has issues.
-- If you try to reset params or hit previous when the instance is closed, the values still change before the alert.
-- Add a checkbox by the zoom controls to disable scroll zoom.
  -- DO THIS. Could also completely remove that, instead set to right click+mouseup or down.
  -- Very sensitive scrolling on a mac laptop, but not on mine.
-- If the species names have bad characters (like a space or !), crashes with a bad error. No need to fix, but make the error comprehensible.
-- Large data sets take a few seconds to visualize, so default 'results' page should have a message reading 'Calculating...', or similar.
-- Clean up displaySummaryStats (iterates too many times right now).
-- Change Instability to Imbalance?
-- Make sure any images loaded are alraedy at their correct width/height. saves time resizing them.
-- When saving svg and opening in inkscape, some weirdness.
  -- The legend is tied to the document, can't be moved. Resizing the document moves it around inexplicably.
  -- Part of the tree extends beyond the borders of the image object itself. Probably to do with the shifting; just resize some aspect.
*/

// Option settings (some are overwritten by calculateDefaultSizes()):
var opts = {
  'sizes': {
    'marker_radius':4, 'highlight_padding':8, 'highlight_shift':5, 'tree':585,
    'species_bar':10, 'bar_chart_height':30, 'bar_chart_buffer':3, 'search_buffer':5,
    'inner_label_buffer':3
  },
  'colours': {
    'clusters':'#FFE0B2', 'highlight':'#FF9800', 'search':'#67F410', // 59E800 00DC58
    'bar_chart':'#333333', 'tree_names_background':'#FFFFFF',
    'species':[
      '#6aad27','#990000','#ff9900','#8b7100','#000099','#A1EF00','#00ffff',
      '#eaeaea','#4f4f4f','#5A0370','#F2E800','#CC08FF','#b5b5b5'
    ]
  },
  'fonts': {
    'tree_names':13, 'legend_names':12, 'family':'Helvetica, Arial, sans-serif'
  }
};
// Page settings:
var server_url='http://'+window.location.host, session_id='',
    maintain_wait=2000, instance_closed=false, maintain_interval;
// Globals from Python backend:
var species=[], sequenceIDs=[], sequence_species={}, species_colours={},
    cluster_list=[], init_weights=[], num_sequences, original_tree_data,
    tree_data, species_tree_data, use_coords, web_version;
// Useful globals:
var r_paper, legend_paper, panZoom, seqs = {}, clusters = {},
    cur_params = [], prev_params = [];
var rad = (Math.PI / 180);

function daemonURL(url) {
  return server_url + '/daemon' + url;
}

function setupPage() {
  session_id = location.search.slice(1);
  $.ajax({
    url: daemonURL('/get-data'),
    type: 'POST',
    data: {'session_id': session_id},
    success: function(data_obj) {
      parseData(data_obj);
      calculateDefaultSizes();
      displayTree();
      setupSummaryStatsPane();
      setupParametersPane();
      setupOptionsPane();
      setupExportPane();
      reclusterTree();
      panZoom = svgPanZoom('#figureSvg', {
        fit: false,
        center: false
      });
      $.ajax({
        url: daemonURL('/page-loaded'),
        type: 'POST',
        data: {'session_id': session_id},
        success: function() { maintain_interval = setInterval(maintainServer, maintain_wait); },
        error: function(error) { processError(error, "Error sending that the page had loaded"); }
      });
    },
    error: function(error) { processError(error, "Error loading data from the server"); }
  })
}
function parseData(data_obj) {
  // change init_weights to an object, not an array.
  var data = $.parseJSON(data_obj);
  maintain_wait = data.maintainwait;
  species_tree_data = data.speciestree;
  original_tree_data = data.treedata;
  sequence_species = data.seqspecies;
  init_weights = data.initweights;
  sequenceIDs = data.sequencenames;
  species = data.specieslist;
  use_coords = data.usecoords;
  web_version = data.webversion;
  cur_params = prev_params = init_weights;
  for (var i=0; i<species.length; ++i) {
    species_colours[species[i]] = opts.colours.species[i % opts.colours.species.length];
  }
  num_sequences = sequenceIDs.length;
  seqs = {};
  for (var i=0; i<num_sequences; ++i) {
    seqs[sequenceIDs[i]] = {};
  }
}
function calculateDefaultSizes() {
  var max_default_font=13, max_highlight_padding=50;
  if (num_sequences > 400) {
    opts.fonts.tree_names = 0;
    opts.sizes.species_bar = 15;
  } else {
    opts.fonts.tree_names = Math.min(Math.floor(300.0/num_sequences)+9, max_default_font);
  }
  opts.sizes.marker_radius = Math.min(Math.floor(200.0/num_sequences)+2, 4);
  opts.sizes.highlight_padding = Math.min(Math.floor(400.0/num_sequences)+2, max_highlight_padding);
  opts.sizes.tree = opts.sizes.tree + (max_default_font-opts.fonts.tree_names)*15; // 20
}
function clearTree() {
  clearClusters();
  r_paper.remove();
  legend_paper.remove();
  $("#canvasSvg").empty();
  $("#legendSvg").empty();
  $("#treeGroup").empty();
}
function clearClusters() {
  var toRemove = ['groupBG', 'groupLines', 'clustBG', 'barChart'];
  for (var clustID in clusters) {
    if (clusters.hasOwnProperty(clustID)) {
      for (var i=0; i<toRemove.length; ++i) {
        if (clusters[clustID].hasOwnProperty(toRemove[i])) {
          clusters[clustID][toRemove[i]].remove();
  } } } }
  $("#clustersTable tbody").empty();
}
function redrawTree() {
  if (r_paper) { clearTree(); }
  displayTree();
  displayClusters();
  setupInteractions();
}
function reclusterTree() {
  $.ajax({
    url: daemonURL('/cluster-tree'),
    type: 'POST',
    data: {'session_id': session_id, 'ILS':cur_params[0], 'dups':cur_params[1],
      'loss':cur_params[2], 'spread':cur_params[3]},
    success: function(data) {
      cluster_list = $.parseJSON(data);
      if (clusters) { clearClusters(); }
      displayClusters(); // Finishes 'seqs' and sets 'clusters'.
      setupInteractions(); // Sets up tooltips, and all click/hover functions.
      displaySummaryStats();
    },
    error: function(error){
      processError(error, "Unable to recluster the tree");
    }
  });
}
function displayTree() {
  loadPhyloSVG(); // Reloads jsPhyloSVG.
  addSpeciesBarToTree();

  Smits.PhyloCanvas.Render.Parameters.jsOverride = 1;
  Smits.PhyloCanvas.Render.Style.text["font-size"] = opts.fonts.tree_names;
  Smits.PhyloCanvas.Render.Style.text["font-family"] = opts.fonts.family;

  var maxLabelLength = getMaxLabelLength(Object.keys(seqs));
  var innerCircle = opts.sizes.tree/2.0 + opts.sizes.marker_radius + opts.sizes.inner_label_buffer - 1;
  var outerRad = innerCircle + maxLabelLength + Smits.PhyloCanvas.Render.Parameters.Circular.bufferOuterLabels + opts.sizes.species_bar + opts.sizes.bar_chart_buffer + opts.sizes.bar_chart_height;
  var canvas_size = outerRad * 2;

  Smits.PhyloCanvas.Render.Style.connectedDash['stroke'] = 'none';
  Smits.PhyloCanvas.Render.Parameters.Circular.bufferRadius = (canvas_size-opts.sizes.tree)/canvas_size;
  Smits.PhyloCanvas.Render.Parameters.Circular.bufferInnerLabels = opts.sizes.inner_label_buffer+opts.sizes.marker_radius+1;
  var dataObject = {phyloxml: tree_data};
  var phylocanvas = new Smits.PhyloCanvas(
    dataObject,
    'svgCanvas',
    canvas_size, canvas_size,
    'circular'
  );
  $("#svgCanvas > svg").attr("id", "treeSvg");
  r_paper = phylocanvas.getSvg().svg;
  drawSpeciesMarkers(); // Fills out some of 'seqs' global.
  drawSearchHighlights();
  var shiftTree = displayLegendTree();
  $("#figureSvg").attr({'width':canvas_size, 'height':canvas_size + shiftTree.y});
  $("#treeSvg").attr({'x':shiftTree.x, 'y':shiftTree.y});
  $("#figureSvg").append($("#legendSvg"));
  $("#treeGroup").append($("#treeSvg"));
  $("#treeGroup").parent().prepend($("#treeGroup"));

  //r_paper.circle(treeDrawingParams['cx'],treeDrawingParams['cy'], innerCircle).attr({stroke:'blue'});
  //r_paper.circle(treeDrawingParams['cx'],treeDrawingParams['cy'], outerRad).attr({stroke:'red'});
  //r_paper.circle(treeDrawingParams['cx'],treeDrawingParams['cy'], treeDrawingParams['barChartRadius'] + opts.sizes.bar_chart_buffer + opts.sizes.bar_chart_height).attr({stroke:'green'});

  var eID, leftBox = '#summaryPane', eIDs = [leftBox, '#treeContainer'];
  setTimeout(function() {
    var header_width = parseFloat($('#mainBody').css('padding-left')) + parseFloat($(leftBox).css('width')) +
        parseFloat($(leftBox).css('margin-right')) + canvas_size - parseFloat($('#header').css('padding-left'));
    for (var i=0; i<eIDs.length; ++i) {
      eID = eIDs[i];
      header_width += parseFloat($(eID).css('padding-left')) + parseFloat($(eID).css('padding-right')) +
          parseFloat($(eID).css('border-top-width'))*2 + parseFloat($(eID).css('margin-left'));
    }
    //$('#header').css('width', header_width+'px');

    //$('#header').css('width', '100%');
    //$('#navMenu').css('width', '100%');
  }, 0);
}
function addSpeciesBarToTree() {
  var xmlDoc = $.parseXML(original_tree_data),
      $xml = $(xmlDoc);
  var ns = $xml.children()[0].prefix, prefix = '';
  if (ns != null) { prefix = ns + '\\:'; }
  var $render = $(prefix+"phyloxml "+prefix+"phylogeny render", xmlDoc);
  var $charts = $("charts", $render);
  $charts.append('<species type="binary" bufferInner="0" thickness="'+opts.sizes.species_bar+'"></species>');
  var $styles = $("styles", $render);
  var styleStr = '';
  for (var spc in species_colours) {
    styleStr += '<'+spc+'_style fill="'+species_colours[spc]+'" stroke="#FFF" stroke-width="1"></'+spc+'_style>';
  }
  $styles.append(styleStr);
  var serializer = new XMLSerializer();
  tree_data = serializer.serializeToString(xmlDoc);
}
function drawSpeciesMarkers() {
  var txt, seqNode, seqID, nodeCoords,
      leafX, leafY, textX, textY;
  var svg = $("#svgCanvas > svg");
  svg.find("text").each(function() {
    txt = $(this);
    seqID = txt.text();
    nodeCoords = parseLeafTextCoords(this);
    leafX=nodeCoords[0], leafY=nodeCoords[1], textX=nodeCoords[2], textY=nodeCoords[3];
    seqs[seqID]['leaf_coords'] = [leafX, leafY];
    //seqs[seqID]['text_coords'] = [textX, textY];
    seqs[seqID]['text_coords'] = moveAwayFromCentre([textX, textY], 5);

    seqNode = r_paper.circle(leafX, leafY, opts.sizes.marker_radius);
    seqNode.attr({fill: species_colours[sequence_species[seqID]], 'stroke-width':0.5});
    $(seqNode.node).attr("class","sequenceNode");
    $(seqNode.node).prop('seqID', seqID);
    txt.attr({'class':'sequenceText', 'font':opts.fonts.tree_names+'px '+opts.fonts.family});
    txt.prop('seqID', seqID);
  });
}
function parseLeafTextCoords(a_obj) {
  var coordsStr = $(a_obj).prev().attr("d");
  var L_ind = coordsStr.indexOf("L");
  var leafCoords = coordsStr.slice(1, L_ind).split(",");
  var txtCoords = coordsStr.slice(L_ind+1).split(",");
  return [parseFloat(leafCoords[0]), parseFloat(leafCoords[1]),
      parseFloat(txtCoords[0]), parseFloat(txtCoords[1])];
}
function drawSearchHighlights() {
  var angleOffset = treeDrawingParams['scaleAngle']/2;
  var seqID, angle, pathStr, leafX, leafY, textX, textY,
      nameHighlight, nodeHighlight, lineHighlight, highlightSet;
  for (var i=0; i<num_sequences; ++i) {
    seqID = treeDrawingParams['seqs'][i][0];
    angle = treeDrawingParams['seqs'][i][1];
    highlightSet = r_paper.set();
    // Highlight around the sequence name.
    pathStr = sectorPathString(
      treeDrawingParams['minBGRadius']+opts.sizes.marker_radius+1,
      treeDrawingParams['barChartRadius']+opts.sizes.bar_chart_buffer+opts.sizes.bar_chart_height+opts.sizes.search_buffer,
      angle-angleOffset, angle+angleOffset);
    nameHighlight = r_paper.path(pathStr);
    $(nameHighlight.node).attr("class","sequenceNode");
    $(nameHighlight.node).prop('seqID', seqID);
    // Highlight around the tree node.
    leafX = seqs[seqID]['leaf_coords'][0], leafY=seqs[seqID]['leaf_coords'][1];
    nodeHighlight = r_paper.circle(leafX, leafY, opts.sizes.marker_radius*1.5+1);
    // Highlight connecting the tree node and the sequence name.
    textX = seqs[seqID]['text_coords'][0], textY = seqs[seqID]['text_coords'][1];
    lineHighlight = r_paper.path('M'+leafX+','+leafY+' L'+textX+','+textY);
    // Grouping the highlights, and storing the object.
    highlightSet.push(nameHighlight, nodeHighlight, lineHighlight);
    highlightSet.attr({'stroke-width':0, fill:opts.colours.search}).toBack().hide();
    lineHighlight.attr({'stroke-width':2, stroke:opts.colours.search});
    seqs[seqID]['searchHighlight'] = highlightSet;
  }
}
function displayLegendTree() {
  Smits.PhyloCanvas.Render.Parameters.Rectangular.alignPadding = 10;
  Smits.PhyloCanvas.Render.Parameters.Rectangular.paddingY = 10;
  Smits.PhyloCanvas.Render.Style.textSecantBg['fill'] = 'none';
  Smits.PhyloCanvas.Render.Style.textSecantBg['stroke'] = 'none';
  Smits.PhyloCanvas.Render.Style.text["font-size"] = opts.fonts.legend_names; // Must be before calculating maxLabelLength
  Smits.PhyloCanvas.Render.Style.text["font-family"] = opts.fonts.family;

  var maxLabelLength = getMaxLabelLength(species) + 10;
  var width = 100 + maxLabelLength,
      height = species.length*15 + 20;
  var legXOffset=1, legYOffset=8, legPad=5, legBorder=2, legTitleHeight=20, legTitleSize=20;
  var boxWidth=width+legPad*2, boxHeight=height+legTitleHeight+legPad*2,
      legendWidth=boxWidth+legXOffset+legBorder, legendHeight=boxHeight+legYOffset+legBorder;
  var dataObject = {newick: species_tree_data};
  Smits.PhyloCanvas.Render.Parameters.Rectangular.bufferX = maxLabelLength;
  var legendPhylocanvas = new Smits.PhyloCanvas(
    dataObject,
    'legendSvgCanvas',
    width, height
  );
  legend_paper = legendPhylocanvas.getSvg().svg;
  $("#legendSvgCanvas > svg").attr("id", "legendSvg");
  var legendSvg = $("#legendSvg");
  drawLegendSpeciesMarkers(legendSvg);
  var legendBox = new Raphael('figureSvg', legendWidth, legendHeight);
  legendBox.rect(legXOffset,legYOffset, boxWidth, boxHeight)
    .attr({fill:'#FFF', stroke:'black', 'stroke-width':legBorder});
  legendBox.text(legXOffset+legBorder/2.0+legPad+2, legYOffset+legBorder/2.0+legPad+legTitleHeight/2.0, 'Legend')
    .attr(Smits.PhyloCanvas.Render.Style.text)
    .attr({'font-size':legTitleSize, font:legTitleSize+'px', 'text-anchor':'start'});
  $("#legendSvg").attr({'x':legXOffset+legPad+legBorder/2.0, 'y':legYOffset+legBorder/2.0+legPad+legTitleHeight});
  return calculateLegendOverlap(legendWidth, legendHeight);
}
function drawLegendSpeciesMarkers(legendSvg) {
  legendSvg.find("text").each(function() {
    var txt = $(this);
    var spc = txt.text();
    var nodeCoords = parseLeafTextCoords(this);
    var leafX=nodeCoords[0], leafY=nodeCoords[1];
    legend_paper.circle(leafX, leafY, 4).attr({fill:species_colours[spc]});
    txt.attr({'font':opts.fonts.legend_names+'px'});
  });
}
function calculateLegendOverlap(legendWidth, legendHeight) {
  // Calculates the x- and y-shift needed so that the tree is as high as possible without overlapping the legend, as well as being left aligned (if it fully fits within the canvas).
  var treeRadius = treeDrawingParams['barChartRadius']+opts.sizes.bar_chart_buffer+opts.sizes.bar_chart_height,
      fromEdge = r_paper.height/2.0-treeRadius,
      x_shift = 0;
  if (fromEdge > 0) {
    legendWidth += fromEdge;
    x_shift = -fromEdge;
  }
  var allowedHeight = (treeRadius+fromEdge) - Math.sqrt((legendWidth-fromEdge)*(2*treeRadius-legendWidth+fromEdge));
  return {'x':x_shift, 'y':Math.max(legendHeight-allowedHeight, 0)};
}
function getMaxLabelLength(orig_names) {
  var names = orig_names.slice(), max = 0, toCheck = Math.min(names.length, 10);
  if (toCheck == 10) {
    names.sort(function(a, b) { return b.length - a.length; });
  }
  var paper = new Raphael('footerDiv', 1000,1000);
  for (var i=0; i<toCheck; ++i) {
    var t = paper.text(0,0, names[i]).attr(Smits.PhyloCanvas.Render.Style.text);
    var w = t.getBBox().width;
    t.remove();
    if (w > max) { max = w; }
  }
  paper.remove();
  return max;
}
function displayClustersList() {
  // Fills out the rest of seqs and clusters globals.
  //Formats the description strings for each sequence and each cluster.
  var score, breakdown, clustLength, seqIDs, seqID, clustID;
  var spcStrs=[], spcCounts={}, numSpc=0, spc, spcColour, spcCountArr=[];
  var tbody = $("#clustersTable > tbody");
  var eventsStr, varStr, breakdownStr, spcCountStr, clustInd = 0;
  $("#numClusters").text(cluster_list.length + " MIGs");
  for (var i=0; i<cluster_list.length; ++i) {
    score = roundFloat(cluster_list[i][0], 2).toFixed(2);
    breakdown = cluster_list[i][1];
    clustLength = cluster_list[i][2];
    seqIDs = cluster_list[i][3];
    clustID = 'group_'+clustInd;
    eventsStr = 'D=<b>'+breakdown[1]+'</b>, I=<b>'+breakdown[0]+'</b>, L=<b>'+breakdown[2]+'</b>';
    varStr = breakdown[3]!=null ? roundFloat(breakdown[3], 2).toFixed(2) : 'N/A';
    breakdownStr = '<b>Score:</b> '+score+'<br /><b>Events:</b> '+eventsStr+'<br /><b>Relative spread:</b> '+varStr;

    for (var j=0; j<species.length; ++j) { spcCounts[species[j]] = 0; }
    spcStrs = [];
    for (var j=0; j<seqIDs.length; ++j) {
      spc = sequence_species[seqIDs[j]];
      spcCounts[spc] += 1;
      spcColour = species_colours[spc];
      spcStrs.push(speciesDot(spcColour, seqIDs[j]));
    }
    numSpc = 0, spcCountArr = [];
    for (var j=0; j<species.length; ++j) {
      if (spcCounts[species[j]] != 0) {
        numSpc += 1;
      }
      spcColour = species_colours[species[j]];
      spcCountArr.push(speciesDot(spcColour, spcCounts[species[j]]));
    }
    spcCountStr = '<b>Species ('+numSpc+'):</b> '+spcCountArr.join(' ');
    // Set the strings and summary stats.
    clusters[clustID] = {};
    clusters[clustID]['score'] = parseFloat(score);
    clusters[clustID]['seqIDs'] = seqIDs;
    clusters[clustID]['hover'] = breakdownStr+'<br />'+spcCountStr+'<br /><b>Sequences ('+seqIDs.length+'):</b><br />'+spcStrs.join('<br />');
    for (var k=0; k<clustLength; ++k) {
      seqID = seqIDs[k];
      seqs[seqID]['clustID'] = clustID;
      seqs[seqID]['hover'] = '<b>Sequence:</b> '+seqID+'<br /><b>Species:</b> '+sequence_species[seqID]+'<br />'+breakdownStr;
    }
    tbody.append("<tr class='clusterRow' clusterID='"+clustID+"'><td>"+score+"</td><td>"+clustLength+"</td></tr>");
    clustInd += 1;
  }
}
function displayClusters() {
  clusters = {};
  displayClustersList();
  var seqID, angle, clustID,
      clustLengths = [];
  var angleOffset = treeDrawingParams['scaleAngle']/2,
      curCluster = seqs[treeDrawingParams['seqs'][0][0]]['clustID'],
      maxScore = null, startAngle = null, nextCluster;
  for (var clustID in clusters) {
    var score = clusters[clustID]['score'];
    if (maxScore == null || score > maxScore) {
      maxScore = score;
    }
  }
  for (var i=0; i<num_sequences; ++i) {
    seqID = treeDrawingParams['seqs'][i][0];
    angle = treeDrawingParams['seqs'][i][1];
    clustID = seqs[seqID]['clustID'];
    if (startAngle == null) { startAngle = angle; }
    if (i == num_sequences-1) { nextCluster = null; }
    else { nextCluster = seqs[treeDrawingParams['seqs'][i+1][0]]['clustID']; }
    if (curCluster != nextCluster) {
      drawGroupsAndBars(clustID, startAngle, angle, angleOffset, maxScore);
      curCluster = nextCluster;
      startAngle = null;
      clustLengths.push([clustID, clusters[clustID]['seqIDs'].length]);
    }
  }
  var seqIDs, clustBG;
  clustLengths.sort(function(a,b) {return a[1]-b[1]});
  for (var i=0; i<clustLengths.length; ++i) {
    clustID = clustLengths[i][0];
    seqIDs = clusters[clustID]['seqIDs'];
    clustBG = fullHighlightObject(clustLengths[i][1], seqIDs);
    clustBG.attr({fill:opts.colours.clusters, stroke:'black', "stroke-width":0.5}).toBack();
    triggerClusterListEvents(clustBG, clustID);
    clusters[clustID]['clustBG'] = clustBG;
  }
}
function drawGroupsAndBars(clustID, startAngle, angle, angleOffset, maxScore) {
  var groupBG, groupLines, barChart;
  var pathStr = sectorPathString(
    treeDrawingParams['minBGRadius']+opts.sizes.marker_radius+1,
    treeDrawingParams['barChartRadius']+opts.sizes.bar_chart_buffer,
    startAngle-angleOffset, angle+angleOffset);
  groupBG = r_paper.path(pathStr).attr({fill:opts.colours.tree_names_background, stroke:'none'}).toBack();
  groupLines = r_paper.path(pathStr).attr({fill:'none', stroke:'black', "stroke-width":1});
  triggerClusterListEvents(groupBG, clustID);
  clusters[clustID]['groupBG'] = groupBG;
  clusters[clustID]['groupLines'] = groupLines;
  var height = roundFloat(clusters[clustID]['score']/maxScore * opts.sizes.bar_chart_height, 4),
      min_radius = treeDrawingParams['barChartRadius']+opts.sizes.bar_chart_buffer;
  if (height > 0) {
    pathStr = sectorPathString(
      min_radius, min_radius+height,
      startAngle-angleOffset*0.9, angle+angleOffset*0.9);
    barChart = r_paper.path(pathStr).attr({fill:opts.colours.bar_chart, stroke:'none'});
    clusters[clustID]['barChart'] = barChart;
  }
}
function triggerClusterListEvents(obj, clustID) {
  obj.mouseover(function() {
    $('.clusterRow[clusterID="'+clustID+'"]').mouseenter();
  }).mouseout(function() {
    $('.clusterRow[clusterID="'+clustID+'"]').mouseleave();
  }).click(function() {
    $('.clusterRow[clusterID="'+clustID+'"]').click();
  });
}
function setupInteractions() {
  $('.sequenceText, .sequenceNode').on({
    "click": function() {
      var settings = {items: '.sequenceText, .sequenceNode',
        content: seqs[$(this)[0].seqID]['hover'],
        tooltipClass: 'sequence-tooltip',
        position: {my:'right bottom', collision:'fit'}};
      startTooltip.call(this, settings);
      $(this).tooltip("open");
    },
    "mouseenter": function() {
      var clustID = seqs[$(this)[0].seqID]['clustID'];
      $('.clusterRow[clusterID="'+clustID+'"]').mouseenter();
    },
    "mouseleave": function() {
      var clustID = seqs[$(this)[0].seqID]['clustID'];
      $('.clusterRow[clusterID="'+clustID+'"]').mouseleave();
      if ($(this).tooltip('instance')) { $(this).tooltip("disable"); }
    }
  });
  $('.clusterRow').on({
    "click": function() {
      var settings = {items: '.clusterRow',
        content: clusters[this.getAttribute('clusterID')]['hover'],
        tooltipClass: 'cluster-tooltip',
        position: {my:'left+10 center', at:'right bottom', of:$("#clustersTable"), collision:'fit'}};
      startTooltip.call(this, settings);
      $(this).tooltip("open");
    },
    "mouseenter": function() {
      var clustID = $(this).attr('clusterID');
      clusters[clustID]['clustBG'].attr({fill:opts.colours.highlight});
      if (opts.fonts.tree_names != 0) {
        clusters[clustID]['groupBG'].attr({fill:opts.colours.clusters});
      }
      $(this).addClass('hover');
    },
    "mouseleave": function() {
      var clustID = $(this).attr('clusterID');
      clusters[clustID]['clustBG'].attr({fill:opts.colours.clusters});
      if (opts.fonts.tree_names != 0) {
        clusters[clustID]['groupBG'].attr({fill:opts.colours.tree_names_background});
      }
      if ($(this).tooltip('instance')) { $(this).tooltip("disable"); }
      $(this).removeClass('hover');
    }
  });
  $('#searchButton').click(function() {
    searchFunction();
  });
  $("#clearSearchButton").click(function() {
    $("#searchInput").attr('value', '');
    searchFunction();
  });
  $('#zoomOutButton').click(function() {
    panZoom.zoomOut();
  });
  $('#zoomInButton').click(function() {
    panZoom.zoomIn();
  });
  $('#zoomResetButton').click(function() {
    panZoom.resetZoom();
    panZoom.resetPan();
  });
}
function startTooltip(settings) {
  // Causes tooltip to remain visible if mouse is over the tooltip or starting element.
  // Taken from: http://stackoverflow.com/questions/16660576/only-close-tooltip-if-mouse-is-not-over-target-or-tooltip
  var fadein_time = 100, fadeout_time = 400;
  $(this).tooltip({
    items: settings.items,
    content: settings.content,
    tooltipClass: settings.tooltipClass,
    position: settings.position,
    show: null, // show immediately
    //open: function(event, ui) {
      // close any lingering tooltips:
      // var $id = $(ui.tooltip).attr('id'); $('div.ui-tooltip').not('#' + $id).remove();
      // or: $('.selector').not(this).tooltip('close');
    //},
    close: function(event, ui) {
      ui.tooltip.hover(
        function() { $(this).stop(true).fadeTo(fadein_time, 1); },
        function() { $(this).fadeOut(fadeout_time, function() { $(this).remove(); }); }
      );
    }
  });
}
function searchFunction() {
  var query = $('#searchInput').attr('value').trim().toLowerCase();
  var seqID, ind, numHits=0;
  for (var i=0; i<num_sequences; ++i) {
    seqID = sequenceIDs[i];
    ind = seqID.toLowerCase().indexOf(query);
    if (query == '' || ind == -1) {
      seqs[seqID]['searchHighlight'].hide();
    } else {
      numHits += 1;
      seqs[seqID]['searchHighlight'].show();
    }
  }
  if (query == '') { $("#searchHitsText").text(''); }
  else { $("#searchHitsText").text(numHits+' hits'); }
  return false;
}
function displaySummaryStats() {
  var avg = 0, avgClstr = 0, spcName, seqID, score;
  for (var i=0; i<cluster_list.length; ++i) {
    avgClstr += cluster_list[i][0];
  }
  avgClstr = roundFloat(avgClstr / cluster_list.length, 2).toFixed(2);
  $('#summaryNumClst').text(cluster_list.length);
  $('#summaryAvgClstIns').text(avgClstr);
  for (var i=0; i<species.length; ++i) {
    var count = 0, avgSpc = 0;
    for (var j=0; j<num_sequences; ++j) {
      seqID = sequenceIDs[j];
      if (sequence_species[seqID] == species[i]) {
        count += 1;
        score = clusters[ seqs[seqID]['clustID'] ]['score'];
        avgSpc += score;
        avg += score;
      }
    }
    avgSpc = roundFloat(avgSpc/count, 2).toFixed(2);
    spcName = cleanString(species[i]);
    $('#summary'+spcName+'AvgIns').text(avgSpc);
  }
  avg = roundFloat(avg/num_sequences, 2).toFixed(2);
  $('#summaryAvgIns').text(avg);
}
function setupSummaryStatsPane() {
  $('#summaryStatsTable').empty();
  $('#summaryStatsTable').append('<tr><td>MIGs:</td><td><b><span id="summaryNumClst"></span></b></td><td><b><span id="summaryAvgClstIns"></span></b></td></tr>');
  $('#summaryStatsTable').append('<tr><td>Sequences:</td><td><b>'+num_sequences+'</b></td><td><b><span id="summaryAvgIns"></span></b></td></tr>');
  for (var i=0; i<species.length; ++i) {
    spcCount = 0;
    for (var j=0; j<num_sequences; ++j) {
      if (sequence_species[sequenceIDs[j]] == species[i]) {
        spcCount += 1;
      }
    }
    $('#summaryStatsTable').append('<tr><td>'+species[i]+':</td><td><b>'+spcCount+'</b></td><td><b><span id="summary'+species[i]+'AvgIns"></span></b></td></tr>');
  }
}
function setupParametersPane() {
  var paraMin=0, paraMax=10;
  var ilsSpin = $("#ils-spinner").spinner({
    min: paraMin, max: paraMax,
    numberFormat: 'N1', step: 0.1
  });
  ilsSpin.spinner('value', init_weights[0]);
  var dupsSpin = $("#dups-spinner").spinner({
    min: paraMin, max: paraMax,
    numberFormat: 'N1', step: 0.1
  });
  dupsSpin.spinner('value', init_weights[1]);
  var lossSpin = $("#loss-spinner").spinner({
    min: paraMin, max: paraMax,
    numberFormat: 'N1', step: 0.1
  });
  lossSpin.spinner('value', init_weights[2]);
  var spreadSpin = $("#spread-spinner").spinner({
    min: paraMin, max: paraMax,
    numberFormat: 'N1', step: 0.1
  });
  if (use_coords) {
    spreadSpin.spinner('value', init_weights[3]);
  } else {
    spreadSpin.spinner('value', 0);
    spreadSpin.spinner('disable');
  }

  $("#paramReclusterButton").click(function() {
    if (validate(ilsSpin, 'ILS') && validate(dupsSpin, 'Duplications') &&
      validate(lossSpin, 'Gene loss') && validate(spreadSpin, 'Variance')) {
      var temp = [ilsSpin.spinner('value'), dupsSpin.spinner('value'),
        lossSpin.spinner('value'), spreadSpin.spinner('value')];
      if (!arraysEqual(cur_params, temp)) {
        prev_params = cur_params;
        cur_params = temp;
        reclusterTree();
      }
    }
  });
  $("#paramPreviousButton").click(function() {
    if (validate(ilsSpin, 'ILS') && validate(dupsSpin, 'Duplication') &&
      validate(lossSpin, 'Gene loss') && validate(spreadSpin, 'Spread')) {
      var temp = cur_params;
      ilsSpin.spinner('value', prev_params[0]);
      dupsSpin.spinner('value', prev_params[1]);
      lossSpin.spinner('value', prev_params[2]);
      spreadSpin.spinner('value', prev_params[3]);
      cur_params = prev_params;
      prev_params = temp;
      reclusterTree();
    }
  });
  $("#paramResetButton").click(function() {
    ilsSpin.spinner('value', init_weights[0]);
    dupsSpin.spinner('value', init_weights[1]);
    lossSpin.spinner('value', init_weights[2]);
    spreadSpin.spinner('value', init_weights[3]);
    if (!arraysEqual(cur_params, init_weights)) {
      prev_params = cur_params;
      cur_params = init_weights;
      reclusterTree();
    }
  });
}
function setupOptionsPane() {
  var fontSpin = $("#font-spinner").spinner({
    min: 0, max: null,
    numberFormat: 'n'
  });
  fontSpin.spinner('value', opts.fonts.tree_names);
  var nodeSpin = $("#node-spinner").spinner({
    min: 0, max: null,
    numberFormat: 'n'
  });
  nodeSpin.spinner('value', opts.sizes.marker_radius);
  var highSpin = $("#highlight-spinner").spinner({
    min: 0, max: null,
    numberFormat: 'n'
  });
  highSpin.spinner('value', opts.sizes.highlight_padding);
  var treeSpin = $("#tree-spinner").spinner({
    min: 0, max: null,
    numberFormat: 'N1', step: 5
  });
  treeSpin.spinner('value', opts.sizes.tree);
  var barSpin = $("#bar-spinner").spinner({
    min: 0, max: null,
    numberFormat: 'N1', step: 1
  });
  barSpin.spinner('value', opts.sizes.bar_chart_height);
  var spcBarSpin = $("#spc-bar-spinner").spinner({
    min: 0, max: null,
    numberFormat: 'N1', step: 1
  });
  spcBarSpin.spinner('value', opts.sizes.species_bar);
  $("#background-colour")[0].jscolor.fromString(opts.colours.tree_names_background);
  $("#cluster-colour")[0].jscolor.fromString(opts.colours.clusters);
  $("#highlight-colour")[0].jscolor.fromString(opts.colours.highlight);
  $("#bar-colour")[0].jscolor.fromString(opts.colours.bar_chart);
  $("#search-colour")[0].jscolor.fromString(opts.colours.search);

  var $spcSelect = $("#species-colour-select");
  for (var i=0; i<species.length; ++i) {
    $spcSelect.append('<option value="'+species[i]+'">'+species[i]+'</option>');
  }
  $spcSelect.children().first().attr('selected', true);
  $spcSelect.change();
  $("#redrawButton").click(function() {
    if (
      validate(fontSpin, 'Font size') &&
      validate(nodeSpin, 'Marker size') &&
      validate(highSpin, 'Highlight size') &&
      validate(treeSpin, 'Tree width') &&
      validate(barSpin, 'Bar chart height') &&
      validate(spcBarSpin, 'Species bar height')
    ) {
      $("#figureSvg").attr('opacity', 0.3);
      setTimeout(function() {
        opts.fonts.tree_names = fontSpin.spinner("value");
        opts.sizes.marker_radius = nodeSpin.spinner("value");
        opts.sizes.highlight_padding = highSpin.spinner("value");
        opts.sizes.tree = treeSpin.spinner("value");
        opts.sizes.bar_chart_height = barSpin.spinner("value");
        opts.sizes.species_bar = spcBarSpin.spinner("value");
        opts.colours.tree_names_background = '#' + $("#background-colour")[0].value;
        opts.colours.highlight = '#' + $("#highlight-colour")[0].value;
        opts.colours.clusters = '#' + $("#cluster-colour")[0].value;
        opts.colours.bar_chart = '#' + $("#bar-colour")[0].value;
        opts.colours.search = '#' + $("#search-colour")[0].value;
        redrawTree();
        $("#figureSvg").attr('opacity', 1.0);
      }, 25);
    }
  });
  $("#saveButton").click(function() {
    // Blob method might be outdated, but appears necessary on Chrome with the large trees. See http://stackoverflow.com/questions/24610694/export-html-table-to-csv-in-google-chrome-browser/24611096#24611096
    var svgData = $("#figureSvg")[0].outerHTML;
    if (web_version == true) { // Just download directly from javascript.
      var svgBlob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"});
      var svgUrl = URL.createObjectURL(svgBlob);
      var downloadLink = document.createElement("a");
      downloadLink.href = svgUrl;
      downloadLink.download = "MIPhy_tree.svg";
      document.body.appendChild(downloadLink);
      downloadLink.click();
      document.body.removeChild(downloadLink);
    } else { // Send to flask server, so the user can choose save location.
      $.ajax({
        url: daemonURL('/save-svg-locally'),
        type: 'POST',
        data: {'session_id': session_id, 'svgData':svgData},
        error: function(error) { processError(error, "Error saving .svg file"); }
      });
    }
  });

}
function validate(spinner, description) {
  if (spinner.spinner("isValid")) {
    return true;
  } else {
    var min = spinner.spinner("option", "min"),
        max = spinner.spinner("option", "max"),
        step = spinner.spinner("option", "step"), msg;
    if (max) { msg = description+" must be between "+min+" and "+max; }
    else { msg = description+" must be greater than "+min; }
    alert(msg+", and be a multiple of "+step+".");
    return false;
  }
}
function speciesSelectChange(spc) {
  var colour;
  if (species_colours.hasOwnProperty(spc)) {
    colour = species_colours[spc];
  } else {
    colour = '#FFFFFF';
  }
  $("#species-colour")[0].jscolor.fromString(colour);
}
function speciesColourChange(jscol) {
  var spc = $('#species-colour-select')[0].value,
      colour = '#' + jscol;
  species_colours[spc] = colour;
}
function setupExportPane() {
  var pane = $("#exportCheckboxes");
  for (var i=0; i<species.length; ++i) {
    var spcChk = $('<label><input type="checkbox" class="exportSpeciesCheckbox" checked />'+species[i]+'</label>');
    spcChk.children().prop('species', species[i]);
    pane.append(spcChk);
  }
  $("#exportSelectAll").change(function() {
    if ($(this).is(':checked')) {
      $(".exportSpeciesCheckbox").prop('checked', true);
    } else {
      $(".exportSpeciesCheckbox").prop('checked', false);
    }
  });
  $("#exportButton").click(function() {
    var checkedSpecies = [], data = [], csvStr = '';
    $(".exportSpeciesCheckbox").filter(':checked').each(function() {
      checkedSpecies.push($(this).prop('species'));
    });
    if (checkedSpecies.length > 0) {
      for (var i=0; i<num_sequences; ++i) {
        var seqID = sequenceIDs[i];
        var spc = sequence_species[seqID];
        if ($.inArray(spc, checkedSpecies) > -1) {
          var clustID = seqs[seqID]['clustID'];
          var score = clusters[clustID]['score'];
          data.push([seqID, spc, clustID, score]);
        }
      }
    }
    data.sort(function(a,b) {
      return ( (a[3]!=b[3]) ? (a[3]-b[3]) : ( // Sort based on score (lowest first). If equal:
        (a[2]!=b[2]) ? ( (a[2]<b[2]) ? 1 : -1 ) : ( // Sort based on clusterID (highest first). If equal:
          (a[1]!=b[1]) ? ( (a[1]<b[1]) ? -1 : 1 ) : ( // Sort based on species (lowest first). If equal:
            (a[0]<b[0]) ? -1 : ( (a[0]>b[0]) ? 1 : 0 ) // Sort based on sequence name (lowest first).
        ))));
    });

    for (var i=0; i<data.length; ++i) {
      csvStr += data[i][0]+','+data[i][1]+','+data[i][2]+','+data[i][3]+'\n';
    }
    if (web_version == true) {
      var csvBlob = new Blob([csvStr], {type:'text/csv'});
      var csvUrl = URL.createObjectURL(csvBlob);
      var downloadLink = document.createElement("a");
      downloadLink.href = csvUrl;
      downloadLink.download = "MIPhy_scores.csv";
      document.body.appendChild(downloadLink);
      downloadLink.click();
      document.body.removeChild(downloadLink);
    } else {
      $.ajax({
        url: daemonURL('/save-csv-locally'),
        type: 'POST',
        data: {'session_id':session_id, 'csvStr':csvStr},
        error: function(error) { processError(error, "Error saving .csv file"); }
      });
    }
  });
}

function closeInstance() {
  instance_closed = true;
  clearInterval(maintain_interval);
  $.ajax({
    url: daemonURL('/instance-closed'),
    type: 'POST',
    data: {'session_id': session_id},
    error: function(error) {
      console.log("Error closing your instance:");
      console.log(error);
    }
  });
}
// Called once the document has loaded.
$(document).ready(function(){
  setupPage();
});

// This is continually called to maintain the background server.
function maintainServer() {
  if (!instance_closed) {
    $.ajax({
      url: daemonURL('/maintain-server'),
      type: 'POST',
      data: {'session_id': session_id},
      error: function(error) {
        console.log('connection to MIPhy server lost.');
        instance_closed = true;
        clearInterval(maintain_interval);
      }
    });
  }
}

$(window).bind('beforeunload', function() {
  closeInstance();
});

function processError(error, message) {
  console.log('Error occurred. The error object:');
  console.log(error);
  if (error.status == 559) {
    alert(message+", as the server didn't recognize the given session ID. This generally means your session has timed out.");
  } else if (error.status == 0) {
    if (web_version) {
      alert(message+", as no response was received. This may mean the web server is down.");
    } else {
      alert(message+", as no response was received. This generally means the program has stopped.");
    }
  } else {
    alert(message+"; the server returned code "+error.status);
  }
}
// = = = = =  Private methods  = = = = =
function speciesDot(colour, text) {
  //return '<span class="speciesDot" style="color:'+colour+';">&#9899</span>'+text; Increase text size.
  return '<span class="speciesDot" style="color:'+colour+';">&#9679</span>'+text;
}
function cleanString(oldString) {
  // Sanitizes a string so it can be used in a jQuery selector.
  var newStr = oldString.replace(/([.|,])/g, "\\$1");
  return newStr;
}
function roundFloat(num, num_dec) {
  var x = Math.pow(10, num_dec);
  return Math.round(num * x) / x;
}
function arraysEqual(a1, a2) {
  if (a1.length != a2.length) { return false; }
  for (var i=0; i<a1.length; ++i) {
    if (a1[i] != a2[i]) { return false; }
  }
  return true;
}
function moveAwayFromCentre(point, distance) {
  // Given point=[x,y], coordinates on the tree svg, returns the coordinates of a point
  // on the line from that point to the centre, 'distance' further away. If a negative
  // distance is given, the point will be closer to the centre.
  var v, len, u, centreX = treeDrawingParams['cx'], centreY = treeDrawingParams['cy'];
  v = [centreX-point[0], centreY-point[1]];
  len = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
  u = [v[0]/len, v[1]/len];
  return [point[0]-distance*u[0], point[1]-distance*u[1]];
}
// Functions to calculate shape of cluster highlight:
function fullHighlightObject(length, seqIDs) {
  // Adapted from http://stackoverflow.com/questions/13802203/draw-a-border-around-an-arbitrarily-positioned-set-of-shapes-with-raphaeljs
  var coords = [], hull, path,
      v, u, len, shift = opts.sizes.highlight_shift,
      treeCentre = [treeDrawingParams['cx'], treeDrawingParams['cy']];
  for (var ind=0; ind<length; ++ind) {
    var leafCoords = seqs[seqIDs[ind]]['leaf_coords'];
    if (shift > 0) { // Looks good, especially important when all nodes same distance from centre.
      leafCoords = moveAwayFromCentre(leafCoords, -shift);
    }
    coords.push(leafCoords);
    coords.push(seqs[seqIDs[ind]]['text_coords']);
  }
  if (coords.length == 2) { hull = coords; }
  else { hull = convexHull(coords); }
  if (opts.sizes.highlight_padding > 0) { hull = expandHull(hull); }
  var point, pathStr = 'M';
  for (var i=0; i<hull.length; ++i) {
    point = hull[i];
    pathStr = pathStr + (i == 0 ? 'M' : 'L') + point[0] + ',' + point[1] + ' ';
  }
  pathStr += 'Z';
  path = r_paper.path(pathStr);
  return path;
}
function convexHull(coords) {
  // Gift wrapping algorithm. expandHull() relies on this proceeding in a
  // particular direction, and the order of the points in the returned hull.
  var left, point;
  for (var i = 0; i < coords.length; i++) {
    point = coords[i];
    if (!left || point[0] < left[0]) { left = point; }
  }
  var hull = [left], p, q;
  for (var i = 0; i < hull.length; i++) {
    p = hull[i];
    q = next_hull_pt(coords, p);
    if (q[0] != hull[0][0] || q[1] != hull[0][1]) { hull.push(q); }
  }
  return hull;
}
function expandHull(hull) {
  // Makes highlight encompass the tree leaves, instead of just passing through
  // their centres. In total, triples the size of the hull.
  var p1,p2,p3,
      shift, scale,
      newHull = [];
  for (var i=0; i<hull.length; ++i) {
    p1 = hull[i];
    if (i==hull.length-2) { p2 = hull[i+1]; p3 = hull[0]; }
    else if (i==hull.length-1) { p2 = hull[0]; p3 = hull[1]; }
    else { p2 = hull[i+1]; p3 = hull[i+2]; }
    //Shifts line segment out by a distance.
    shift = [p1[1]-p2[1], p2[0]-p1[0]]; //perpendicular to p1p2
    scale = opts.sizes.highlight_padding / Math.sqrt(dist(shift,[0,0]));
    shift = [shift[0]*scale, shift[1]*scale];
    //[p1[0]+shift[0], p1[1]+shift[1]] dont have to check, because of convex hull directionality.
    newHull.push([p1[0]-shift[0], p1[1]-shift[1]]);
    newHull.push([p2[0]-shift[0], p2[1]-shift[1]]);
    //Adds a new point between the newly created line segments. More rounded, especially on small clusters. Found by normalizing p1p2 to the length of p2p3, then a line from the mean of those segments to p2.
    scale = Math.sqrt(dist(p1, p2)) / Math.sqrt(dist(p2, p3));
    shift = [(p2[0]-p3[0]+(p2[0]-p1[0])/scale)/2.0, (p2[1]-p3[1]+(p2[1]-p1[1])/scale)/2.0];
    scale = opts.sizes.highlight_padding / Math.sqrt(dist(shift,[0,0]));
    shift = [shift[0]*scale, shift[1]*scale];
    newHull.push([p2[0]+shift[0], p2[1]+shift[1]]);
  }
  return newHull;
}
function next_hull_pt(coords, p) {
  var q = p, r, t;
  for (var i = 0; i < coords.length; i++) {
    r = coords[i];
    t = turn(p, q, r);
    if (t == -1 || t == 0 && dist(p, r) > dist(p, q)) { q = r; }
  }
  return q;
}
function turn(p, q, r) {
  var x = (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]);
  if (x > 0) { return 1; }
  else if (x < 0) { return -1; }
  else { return 0; }
}
function dist(p, q) {
  var dx = q[0] - p[0];
  var dy = q[1] - p[1];
  return dx * dx + dy * dy;
}
// Functions to draw background of sequences in a cluster:
function sectorPathString(r1, r2, y1, y2) {
  // Adapted from sector() and secant() from jsphylosvg.js
  var coords1 = secPosition(r1, y1), coords2 = secPosition(r2, y2);
  return [["M", coords1[0], coords1[1]], secant(r1, y1, y2, 0),
            ["L", coords2[0], coords2[1]], secant(r2, y2, y1, 1), ['Z']];
}
function secant(r, startAngle, endAngle, invSecant){
  var endPos = secPosition(r, endAngle);
  var n, inv = 0;
  if(Math.abs(normalizeAngle(endAngle-startAngle)) > 180) {
    n = 1;
  } else {
    n = -1;
  }
  if(invSecant){
    n *= -1;
    inv = 1;
  }
  return ["A", r, r, 0, n < 1 ? 0 : 1, inv, endPos[0], endPos[1]];
}
function secPosition(r, deg){
  deg += treeDrawingParams['initStartAngle'];
  return [roundFloat(treeDrawingParams['cx'] + r * Math.sin(deg * rad), 4),
          roundFloat(treeDrawingParams['cy'] + r * Math.cos(deg * rad), 4)];
}
function normalizeAngle(ang){
  while(ang > 360 || ang < 0) {
    if(ang > 360){ ang -= 360; }
    else if (ang < 0){ ang += 360; }
  }
  return ang;
}
