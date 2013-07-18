/*
 * HEALPix Java code supported by the Gaia project.
 * Copyright (C) 2006-2011 Gaia Data Processing and Analysis Consortium
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

package healpix.core.base.set.test;

import healpix.core.HealpixIndex;
import healpix.core.base.set.LongList;
import healpix.core.base.set.LongSet;
import healpix.tools.SpatialVector;
import junit.framework.TestCase;


public class QueryDiscTest extends TestCase {
	public void testQueryDisc () {
	int nside = 32;
	
	int inclusive = 0;//false;
	double radius = Math.PI;
	double radius1 = Math.PI/2.;
    HealpixIndex pt = new HealpixIndex();
    long npix = HealpixIndex.nside2Npix(nside);
    double res = HealpixIndex.getPixRes(nside); // pixel size in radians
    System.out.println("res="+res);
    double pixSize = Math.toRadians(res/3600.0); // pixel size in radians
    System.out.println("pixSize="+pixSize+" rad");

    
    LongList fullSky = new LongList(pt.query_disc(nside, new SpatialVector(0., 0., 1.), radius,  0, inclusive));
    LongList firstHalfSky = new LongList(pt.query_disc(nside, new SpatialVector(0., 0., 1.), radius1, 0,  inclusive));
    LongList secondHalfSky = new LongList(pt.query_disc(nside, new SpatialVector(0., 0., -1.), radius1, 0, inclusive));
    firstHalfSky.addAll(secondHalfSky);
    LongSet pixHalfsUnique = new LongSet(firstHalfSky);
    LongList pixHalfsList = new LongList(pixHalfsUnique);
    pixHalfsList = pixHalfsList.sort();
    fullSky = fullSky.sort();

    long listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {

    assertEquals(fullSky.get(i),pixHalfsList.get(i));
    }
    


   firstHalfSky = new LongList(pt.query_disc(nside, new SpatialVector(1., 0., 0.), radius1, 0,inclusive));
   secondHalfSky = new LongList(pt.query_disc(nside, new SpatialVector(-1., 0., 0.),radius1, 0, inclusive));
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = new LongSet(firstHalfSky);
    pixHalfsList = new LongList(pixHalfsUnique);
    
    pixHalfsList = pixHalfsList.sort();
    System.out.println("full size="+fullSky.size()+" half size="+pixHalfsList.size());
    listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {
//        System.out.println( "i="+i+" "+fullSky.get(i)+" "+pixHalfsList.get(i));
        assertEquals(fullSky.get(i),pixHalfsList.get(i));
        }


    firstHalfSky = new LongList(pt.query_disc(nside, new SpatialVector(0., 1., 0.), radius1,  0,inclusive));
    secondHalfSky = new LongList(pt.query_disc(nside, new SpatialVector(0., -1., 0.), radius1, 0, inclusive));
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = new LongSet(firstHalfSky);
    pixHalfsList = new LongList(pixHalfsUnique);
    pixHalfsList = pixHalfsList.sort();
    System.out.println("full size="+fullSky.size()+" half size="+pixHalfsList.size());
    listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());

    for ( int i=0; i< listL; i++) {

        assertEquals(fullSky.get(i),pixHalfsList.get(i));
        }
        
}
}
