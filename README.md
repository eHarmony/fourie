Author: Will Dabney (@amarack) during his internship at eHarmony.

# Example feature exppansions

nvars - number of variables/attributes to expand/interact
order - how much to expand

We get (order +1) ^ nvars variables

For expansion of 1 variable:
```
nvars=1 order=0 const
nvars=1 order=1 const, cos(1* pi * x)
nvars=1 order=2 const, cos(2* pi * x)
```
For 2 variables we get this expansion:

nvars=2
```
const
cos pi y
cos pi x
cos pi (x+y)
cos pi (2y)
cos pi (2x)
cos pi (x + 2y)
cos pi (x + 2y)
cos pi (2x + 2y)
```
need to scale x between 0 and 1
or -1,1 and need to include sin() features

```
./fourie order

./fourie 3
-1,-2
1,2
0.1,1
[expanded feats]

./fourie 3 --vw
FOO,BAR


WDABNEYMacBookPro15:fourie wdabney$ ./fourie 3
-1,-1
1,1
0.5232,-0.232
1,0.356412,-0.745941,-0.888136,-0.7324,-0.897197,0.0928573,0.963388,0.0728206,0.957802,0.609924,-0.523034,0.625733,-0.505793,-0.986274,-0.197246
-0.0232,-0.123213
1,0.192336,-0.926013,-0.548549,0.0364344,-0.97367,-0.410979,0.815577,-0.997345,-0.263287,0.896066,0.607979,-0.10911,0.954484,0.476274,-0.771275
^C
WDABNEYMacBookPro15:fourie wdabney$ ./fourie 3 --nonorm
0.53,0.2323,0.71232
1,-0.61865,-0.234545,0.908852,0.745313,-0.984903,0.473307,0.39928,0.110983,-0.849473,0.940069,-0.313674,-0.579879,-0.281343,0.927985,-0.866852,-0.0941084,-0.72396,0.989864,-0.500798,-0.733896,-0.079649,0.832446,-0.950336,-0.999856,0.605233,0.251002,-0.915797,-0.756516,0.981825,-0.458296,-0.414776,-0.982287,0.754911,0.0482359,-0.814593,-0.607182,0.999895,-0.629988,-0.220411,0.0772064,0.735558,-0.987312,0.486043,0.722268,0.0965474,-0.841726,0.944919,0.278991,0.581873,-0.998943,0.654118,0.848178,-0.108548,-0.713871,0.991821,0.985325,-0.743677,-0.0651729,0.824316,0.620573,-0.999997,0.616723,0.236926
WDABNEYMacBookPro15:fourie wdabney$ ./fourie 3 --vw 
FOO,BAR
-1,-1
1,1
0.02343324 |namespace X:0.3423 Y:9343.234324 FOO:0.234 BAR:-0.343 Z:2.3434
0.02343324 |namespace X:0.3423 Y:9343.234324 Z:2.3434 FOURIER0:1 FOURIER1:0.513092 FOURIER2:-0.473473 FOURIER3:-0.998963 FOURIER4:-0.359345 FOURIER5:-0.985378 FOURIER6:-0.651834 FOURIER7:0.316477 FOURIER8:-0.741742 FOURIER9:0.19509 FOURIER10:0.94194 FOURIER11:0.771513 FOURIER12:0.892428 FOURIER13:0.845168 FOURIER14:-0.0251303 FOURIER15:-0.870957 
```
