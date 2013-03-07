function identity(n){
var riga="\n";
for (i=1; i<=n; i++){
     for(j=1;j<=n; j++){
	 if(i==j)
            riga += 1 +"\t";
	 else
            riga += 0 +"\t";
     }
     riga +="\n";
}
return riga;
}
