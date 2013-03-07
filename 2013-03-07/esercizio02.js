function fibonacci(n){
var a=1;
var b=1;
var somma;
for(i=0; i<n; i++){
   somma = a+b;
   a= b;
   b= somma;
}
fibonacci[1]=1;
fibonacci[2]=1;
return b;
}
