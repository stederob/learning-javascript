var n = 10;
var riga="";
for (i=1; i<=n; i++){
     for(j=1;j<=n; j++){
         riga += i*j +"\t";
     }
     riga +="\n";
}
console.log(riga);
