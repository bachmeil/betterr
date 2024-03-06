/* Use the GNU Scientific Library for random number generation.
 * Includes both sequential and parallel generation.
 * Sequential generation should be faster if you're only running one process,
 * as will be the case if you're interacting with R.
 * Parallel generation is done using L'Ecuyer's code ported from Java.
 * ImportC is used to compile the C headers for GSL.
 */
import gslheaders, prng;
import std.parallelism, std.stdio;

void main() {
	StreamGen sg;
	
	// Create an array of GSLStream structs
	Stream[] streams;
	foreach(_; 0..totalCPUs) {
		streams ~= sg.createStream();
	}
	
	// This function has to be static for it to used with a template
	static double[] simulate(Stream s) {
		auto r = GslRng(s);
		double[] result;
		foreach(_; 0..20) {
			result ~= gsl_ran_ugaussian(r);
		}
		return result;
	}
	
	// Work over the streams in parallel
	auto result = taskPool.map!simulate(streams);
	writeln(result);
  writeln("Total CPUs: ", totalCPUs);
}
