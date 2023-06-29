import betterr.r;
import betterr.plot;
import betterr.vector;

void main() {
	startR();
	Vector x = [1.7, -0.2, 3.4, -0.2, 6.8, 1.2];
	auto p = Plot(x.name);
	p.main = "Plot of random points";
	p.sub = "Just some test data, actually";
	p.type = "b";
	p.xlab = "x values";
	p.ylab = "y values";
	p.create("myplot");
	closeR();
}
