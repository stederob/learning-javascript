var n = 10;
var riga="";
for (i=1; i<=n; i++){
     for(j=1;j<=n; j++){
	 if(j==10)
            riga += i*j +"\t";
	 else
            riga += i*j +",\t";
     }
     riga +="\n";
}
console.log(riga);
