var numrighe = 10;
var numColonne = 10;
var matrice = "", riga, r, c;

for (r = 1; r <= numrighe; r++) {
       riga = "";
       for (c = 1; c <= numColonne; c++) {
            if (r == c) {
                riga += "<td> 1, </td>";
            } else {
                riga += "<td> 0, </td>";
            }
        }
        matrice += "<tr>" + riga + "</tr>";
    }
document.write("<table>"+matrice+"</table>");
