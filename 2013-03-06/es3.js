var n = 10;
var riga="";
for (i=1; i<=n; i++){
     for(j=1;j<=n; j++){
	 if(i==j)
            riga += 1;
	 else
            riga += 0;
	 if(j==10)
            riga += "\t";
         else
            riga += ",\t";
     }
     riga +="\n";
}
console.log(riga);
