// window.location.origin
var server_url='http://'+window.location.host, results_url=server_url+'/results',
    session_id='', upload_success=false, results_processed=false,
    maxUploadSize = 10 * 1024*1024; // 10MB;
var init_weights = [0.5, 1.0, 1.0, 1.0];


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
  spreadSpin.spinner('value', init_weights[3]);
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
  $("#paramResetButton").click(function() {
    ilsSpin.spinner('value', init_weights[0]);
    dupsSpin.spinner('value', init_weights[1]);
    lossSpin.spinner('value', init_weights[2]);
    spreadSpin.spinner('value', init_weights[3]);
  });
  $("#spreadRefinementCheck").change(function() {
    if ($(this).is(':checked')) {
      spreadSpin.spinner('enable');
    } else {
      spreadSpin.spinner('disable');
    }
  });
  $("#uploadButton").click(function() {
    if (validateUploadValues(ilsSpin, dupsSpin, lossSpin, spreadSpin) == false) {
      return false;
    }
    var form_data = new FormData($('#upload-files')[0]); // The 2 uploaded files
    form_data.append('usecoords', $("#spreadRefinementCheck")[0].checked); // The state of the 'Spread refinement' check box
    var param_data = {'ILS':ilsSpin.spinner('value'), 'dups':dupsSpin.spinner('value'),
      'loss':lossSpin.spinner('value'), 'spread':spreadSpin.spinner('value')};
    results_processed = false;
    $("#summaryText").html('<p>Uploading and processing files...</p>');
    $.ajax({
      type: 'POST',
      url: daemonURL('/upload-files'),
      data: form_data,
      contentType: false,
      cache: false,
      processData: false,
      success: function(data_obj) {
        var data = $.parseJSON(data_obj);
        session_id = data.idnum;
        var num_sequences = data.numseqs;
        var action_msg = data.actionmsg;
        var action_info = data.actioninfo;
        if (action_msg == 'over seq limit') {
          alert('You uploaded a tree containing '+num_sequences+' sequences, which is above the server limit of '+action_info+'. If you want to process this tree you will have to download MIPhy and run it locally.');
          return false;
        } else if (action_msg == 'over refine limit') {
          alert('You uploaded a tree containing '+num_sequences+' sequences, which is above the server refine limit of '+action_info+', so Spread refinement has been disabled. If you want to include this step, please download MIPhy and run it locally.');
          //spreadSpin.spinner('disable');
          $("#spreadRefinementCheck").trigger('click');
        } else if (action_msg) {
          alert('Processing stopped, as an unknown message was returned from the server: '+action_msg);
          return false;
        }
        processFiles(param_data);
      },
      error: function(error) {
        results_processed = false;
        $("#summaryText").html('<p>No files uploaded.</p>');
        processError(error, "Error uploading files");
      }
    });
  });
}

function processFiles(param_data) {
  console.log('param', param_data);
  var process_data = $.extend({'session_id':session_id}, param_data);
  console.log('proc', process_data);
  //data: {'session_id':session_id},
  results_processed = false;
  $.ajax({
    url: daemonURL('/process-data'),
    type: 'POST',
    data: process_data,
    success: function(data_obj) {
      var data = $.parseJSON(data_obj);
      var num_sequences = data.numseqs;
      var num_species = data.numspc;
      var num_clusters = data.numclstrs;
      results_processed = true;
      $("#summaryText").html(
        '<p>Upload and analysis successful.<br /><b>'+num_sequences+'</b> sequences,<br />from <b>'+
        num_species+'</b> species,<br />in <b>'+num_clusters+'</b> clusters.</p>'
      );
    },
    error: function(error) {
      results_processed=false;
      processError(error, "Error processing data");
    }
  });
}


$(document).ready(function(){
  setupParametersPane();

  $("#viewResultsButton").click(function(){
    if (results_processed == true) {
      window.open(results_url+'?'+session_id, '_blank');
    } else {
      alert("You must upload and process valid files before viewing the MIPhy results.");
    }
  });
});

function processError(error, message) {
  console.log('Error occurred. The error object:');
  console.log(error);
  if (error.status == 550) {
    alert("Error validating the uploaded files: "+error.responseText);
  } else if (error.status == 551) {
    alert("Unknown error while validating the uploaded files: "+error.responseText);
  } else if (error.status == 552) {
    alert("Error sending files or Spread refinement option to the server: "+error.responseText);
  } else {
    alert(message+"; the server returned code "+error.status);
  }
}
function daemonURL(url) {
  return server_url + '/daemon' + url;
}
function validateUploadValues(ilsSpin, dupsSpin, lossSpin, spreadSpin) {
  if ($('#upload-tree-input')[0].files.length != 1 || $('#upload-info-input')[0].files.length != 1) {
      alert("You must select 1 gene tree and 1 info file to proceed.");
      return false;
  }
  var file_size = $('#upload-tree-input')[0].files[0].size +
                  $('#upload-info-input')[0].files[0].size;
  if (file_size > maxUploadSize) {
    alert("You're attempting to upload files that are larger than the allowed total file size of "+maxUploadSize/(1024*1024)+"MB.");
    return false;
  }
  if (!validate(ilsSpin, 'ILS') || !validate(dupsSpin, 'Duplications') ||
  !validate(lossSpin, 'Gene loss') || !validate(spreadSpin, 'Variance')) {
    return false;
  }
  return true;
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
