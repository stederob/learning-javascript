var numrighe = 10;
var numColonne = 10;
var matrice = "", riga, r, c;

for (r = 1; r <= numrighe; r++) {
       riga = "";
       for (c = 1; c <= numColonne; c++) {
            if (c == 1) {
                riga += "<td>" + r + "</td>";
            } else {
                riga += "<td>" + (r * c) + "</td>";
            }
        }
        matrice += "<tr>" + riga + "</tr>";
    }
document.write("<table>"+matrice+"</table>");
