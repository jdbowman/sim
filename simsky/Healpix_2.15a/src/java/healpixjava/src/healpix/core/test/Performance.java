package healpix.core.test;
import healpix.core.HealpixIndex;
import healpix.core.base.set.LongRangeSet;
import healpix.tools.SpatialVector;


/**  
 * measures performance
 */
public class Performance {

	/** helper class to print time */
	static public class StopWatch{
		
		long start = 0;
		public void start(){ 
			start = System.currentTimeMillis();
		}
		
		public void printTime(String label){ 
			long time  = System.currentTimeMillis() - start;
			System.out.println(label+" "+time+" ms");
		}

	}
	
	
	static StopWatch sw = new StopWatch();
	static HealpixIndex t = new HealpixIndex();
	static SpatialVector centre = new SpatialVector(1,1,1);
	static int nside = 0;
	static double radius = 0;
	static LongRangeSet result = null;

	public static void main(String[] args) {
		centre.normalized();
		nside = t.calculateNSide(60);
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
	
		nside = t.calculateNSide(60);
		radius = Math.toRadians(10);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("10 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
		
		nside = t.calculateNSide(1);
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");

		nside = t.calculateNSide(1);
		radius = Math.toRadians(10);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("10 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
		
		nside = 1048576; //highest res available with long ranges
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");

		nside = 1048576; //highest res available with long ranges
		radius = Math.toRadians(10);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("10 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
		
		nside = 536870912; //highest res documented in C++ original code. Order=29
		radius = Math.toRadians(0.5);
		sw.start();
		result = t.query_disc(nside, centre, radius, 0,1);
		sw.printTime("0.5 degrees at NSIDE="+nside+"  have "+result.size()+" pixels and took");
	}
}
