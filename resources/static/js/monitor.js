var server_url='http://'+window.location.host;

function daemonURL(url) {
  return server_url + '/daemon' + url;
}

$(document).ready(function(){
  $("#currentSessionsButton").click(function(){
      $.ajax({
          url: daemonURL("/get-sessions"),
          type: 'GET',
          success: function(data_obj) {
            var data = $.parseJSON(data_obj);
            console.log('session info returned', data);
          },
          error: function(error){ console.log(error); }
        });
  });
});
