REPO=PATH/TO/betterr
LINK=-L/usr/lib/R/site-library/RInside/lib/libRInside.so -L-lR
MATRIX=-L/usr/lib/R/library/Matrix/libs/Matrix.so
QUADPROG=-L/usr/local/lib/R/site-library/quadprog/libs/quadprog.so
FUNCPTR=-L/usr/local/lib/R/site-library/funcptr/libs/funcptr.so

vector:
	ldmd2 -i testvector.d -of=vector.exe $(LINK)

data:
	ldmd2 -i testdataframe.d -of=dataframe.exe $(LINK)

list:
	ldmd2 -i testlist.d -of=list.exe $(LINK)

lm:
	ldmd2 -i testlm.d -of=lm.exe $(LINK)

base:
	ldmd2 -i testbase.d -of=base.exe $(LINK)

matrix:
	ldmd2 -i testmatrix.d -of=matrix.exe $(LINK)

ts:
	ldmd2 -i testts.d -of=ts.exe $(LINK)

random:
	ldmd2 -i testrandom.d -of=random.exe $(LINK)

array:
	ldmd2 -i testarray.d -of=array.exe $(LINK)

plot:
	ldmd2 -i testplot.d -of=plot.exe $(LINK)

optim:
	ldmd2 -i testoptim.d -of=optim.exe $(LINK)

crantest:
	ldmd2 -i testcran.d -of=cran.exe $(LINK) $(MATRIX)
	
blastest:
	ldmd2 -i testblas.d -of=blas.exe $(LINK) -L-lopenblas -version=OPENBLAS -version=CBLAS_STD_COMPLEX

quadprog:
	ldmd2 -i testqp.d -of=qp.exe $(LINK) $(QUADPROG)

pair:
	ldmd2 -i testpair.d -of=pair.exe $(LINK) $(MATRIX)

memory:
	ldmd2 -i memory.d -of=memory.exe $(LINK)

testmir:
	ldmd2 -i testmir.d -of=mir.exe

funptr:
	ldmd2 -i testfunptr.d -of=funptr.exe $(LINK) $(FUNCPTR)

dep:
	mkdir -p betterr
	cp $(REPO)/*.d $(shell pwd)/betterr

cblas:
	cp $(REPO)/cblas $(shell pwd)/cblas

cran:
	cp $(REPO)/cran $(shell pwd)/cran

mir:
	cp $(REPO)/mir $(shell pwd)/mir

prng:
	ldmd2 -i testgslprng.d -L-lgsl
