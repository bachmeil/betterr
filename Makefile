REPO=~/repos/betterr
LINK=-L/usr/lib/R/site-library/RInside/lib/libRInside.so -L-lR
MATRIX=-L/usr/lib/R/library/Matrix/libs/Matrix.so
QUADPROG=-L/usr/local/lib/R/site-library/quadprog/libs/quadprog.so
FUNCPTR=-L/usr/local/lib/R/site-library/funcptr/libs/funcptr.so

dep:
	mkdir -p betterr
	cp $(REPO)/*.d $(shell pwd)/betterr

cblas:
	cp $(REPO)/cblas $(shell pwd)/cblas

cran:
	cp $(REPO)/cran $(shell pwd)/cran

mir:
	cp $(REPO)/mir $(shell pwd)/mir
