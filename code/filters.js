var photchoices = ["Y", "J", "H", "K"];
var specchoices = ["YJ", "HK"];

function setOptions(chosen,change){
    var selbox = change;
    selbox.options.length = 0;

    if (chosen == "0") {
        for (var i=0; i<photchoices.length; i+=1){
            selbox.options[selbox.options.length] = new Option(photchoices[i]);
        }
    }
    if ((chosen == "30") || (chosen == "4000")) {
        for (var i=0; i<specchoices.length; i+=1){
            selbox.options[selbox.options.length] = new Option(specchoices[i], specchoices[i]);
        }
    }
} 