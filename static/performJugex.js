$(function() {
    $('a#calculate').bind('click', function() {
        $.getJSON('/_jugex', {
                  }, function(data) {
                      var elem = document.getElementById("myBar");
                      var width = 10;
                      var id = setInterval(frame, 10);
                      function frame() {
                          if (width >= 100) {
                              clearInterval(id);
                          } else {
                              width++;
                              elem.style.width = width + '%';
                              elem.innerHTML = width * 1  + '%';
                          }
                      }
                      var json = JSON.parse(data.result);
                      for(var i=0; i < json.genes.length; i++)
                          console.log(json.genes[i].name+' '+json.genes[i].pval);
                      var table = document.createElement("TABLE");
                      table.border = "1";
                      var columnCount = 2;
                      var row = table.insertRow(-1);
                      var headerCell = document.createElement("TH");
                      headerCell.innerHTML = 'Name of Gene';
                      row.appendChild(headerCell);
                      headerCell = document.createElement("TH");
                      headerCell.innerHTML = 'P Value';
                      row.appendChild(headerCell);
                      for(var i=0; i<json.genes.length; i++){
                          row = table.insertRow(-1);
                          cell = row.insertCell(-1);
                          cell.innerHTML = json.genes[i].name;
                          cell = row.insertCell(-1);
                          cell.innerHTML = json.genes[i].pval;
                      }
                      var dvTable = document.getElementById("dvTable");
                      dvTable.innerHTML = "";
                      dvTable.appendChild(table);
                  });
        return false;
    });
});

$(function(){
    $.ajax({
        url: '/_autocomplete'
        }).done(function (data){
            $('#gene_autocomplete').autocomplete({
                source: data,
                minLength: 2,
                select: function(event, ui) {
                    console.log("You selected: " + ui.item.label);
                }
            });
        });
    });
