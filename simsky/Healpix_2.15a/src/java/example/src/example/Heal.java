package example;

// get access to the healpix lib core part
import healpix.core.*;


public class  Heal {

    // args are the command line args
    public static void main(String[] args) throws Exception{

	System.out.println("Hello Eric !");

	String raStr = "90.9";
	String decStr = "-20.2";

	// process command line
	int ind =0;
	if (args.length > ind) {
		raStr = args[ind++];
	}
	if (args.length > ind) {
		raStr = args[ind++];
	}

	double ra = Double.parseDouble(raStr);
	double dec = Double.parseDouble(decStr);

	// do healpix

	int pix  = Healpix.ang2pix_nest(1024,Math.toRadians(ra),Math.toRadians(dec));
	System.out.println("ra:"+ra+" dec:"+dec+" pix:"+pix);

    }

}

