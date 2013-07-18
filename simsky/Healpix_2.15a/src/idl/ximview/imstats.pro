; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2008   J. P. Leahy
;
;
;  This file is part of Ximview and of HEALPix
;
;  Ximview and HEALPix are free software; you can redistribute them and/or modify
;  them under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  Ximview and HEALPix are distributed in the hope that they will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;
; -----------------------------------------------------------------------------
FUNCTION imstats, image, xpix, ypix, boxsize, ASTROM = astrom, $
                  LUN = lun, UNIT = unit
;+
; NAME:
;       IMSTATS
;
; PURPOSE:
;       Prints basic statistics on either (1) an image box surrounding
;       specified position or (2) a region defined by a coordinate
;       list. Blanks (NaNs) are ignored. Statitistics are also
;       calculated ignoring exactly zero pixels.
;
;       Results are:
;        * Position and values of minimum and maximum pixels (World
;          Coordinates as well as pixel coordinates if ASTROM is set)  
;        * Mean, median, and standard deviation of values
;        * Integrated flux density (if there is enough astrometric
;          information to derive a pixel size.)
;
; CATEGORY:
;       Image processing, Statistics
;
; CALLING SEQUENCE:
;   
;       Result = IMSTATS(Image, Xpix, Ypix, Boxsize)
;
; INPUTS:
;     Image:  2-D image from which statistics are extracted
;
;     Xpix:   (1) X coordinate of box centre
;             (2) list of X coordinates for all pixels in region
;
;     Ypix:    Y coordinates as above
;
; OPTIONAL INPUTS:
;     Boxsize: 2-element array specifying size of box, or one value if
;              box is square. If specified, Xpix & Ypix must be scalar.
;
; KEYWORD PARAMETERS:
;     UNIT:    intensity unit of image
;
;     LUN:     logical unit number for log file
;
;     ASTROM:  astrolib astrometry structure (if any, otherwise 0)
;
; OUTPUTS:
;     Returns coordinates of statistics box BLC and TRC, after
;     adjustment for edges of image. 
;
; SIDE EFFECTS:
;     Prints results to screen and to log file if LUN is set.
;
; EXAMPLE:
;
;          test = RANDOMN(seed,50,50)
;          boxcoord = IMSTATS(test, 30, 30, 33)
;
;  Prints:  
;
;    IMSTATS: Image statistics for    14 <= x <=   46,    14 <= y <=   46
;             Includes   1089 sky pixels of which  1089 contain valid data
; 
;             Maximum:   3.396E+00 at pixel (   14,   14)
;             Minimum:  -3.212E+00 at pixel (   44,   23)
;
;             Mean:          -1.669E-02  Median: -2.421E-03
;             Standard Dev:   9.850E-01  Unit:   unknown
;     
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, 2008
;-
COMPILE_OPT IDL2
ON_ERROR, 2

r2d = 180d0 / !dpi
line = ''
T = SIZE(image)
IF T[1] LT 3 OR T[2] LT 3 THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'Image too small to analyse!'
    GOTO, QUIT
ENDIF

llbb = KEYWORD_SET(astrom) ; we can work out longs & lats

IF N_ELEMENTS(unit) EQ 0 THEN unit = 'unknown'
XS = SIZE(xpix)  & YS = SIZE(ypix)
IF ARRAY_EQUAL(XS, YS) EQ 0B THEN MESSAGE, 'Mismatched X and Y pixel arrays'
npix = N_ELEMENTS(xpix)
is_box = N_PARAMS() EQ 4

IF is_box THEN BEGIN  ; Convert box to pixel list.
    IF npix GT 1 THEN $
      MESSAGE, 'Coordinate list and box size should not both be supplied'

    B = SIZE(boxsize)
    IF B[1] EQ 1 THEN boxsize = REPLICATE(boxsize,2)

    xypix = [xpix, ypix]
    half = boxsize / 2
    xy0 = ((xypix - half) > 0) < (T[1:2] - 1)
    xy1 = ((xypix + half) > 0) < (T[1:2] - 1)
    dxy = xy1 - xy0 + 1
    nvalid = dxy[0]*dxy[1]

    coords = [[xy0], [xy0[0], xy1[1]], [xy1], [xy1[0], xy0[1]]] ; Returned

    IF nvalid EQ 0 THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Box outside image, no data to analyse'
        GOTO, QUIT
    ENDIF
    box = image[xy0[0]:xy1[0],xy0[1]:xy1[1]]

    x = LINDGEN(dxy[0]) + xy0[0]  &  y = LINDGEN(dxy[1]) + xy0[1]
    x = REBIN(x,dxy[0],dxy[1])  & y = TRANSPOSE(REBIN(y,dxy[1],dxy[0]))

    line = [' ', STRING(xy0[0], xy1[0], xy0[1], xy1[1], FORMAT = $
"('IMSTATS: Image statistics for ',I5,' <= x <=',I5,', ',I5,' <= y <=',I5)")]
ENDIF ELSE BEGIN                ; Pixel list supplied
    impix = T[1]*T[2]
    idx = xpix + ypix*T[1]
    box = image[idx]
    good = WHERE(idx GE 0 AND idx LT impix, nvalid)

                                ; Find bounding box of region
    x01 = MINMAX(xpix)
    y01 = MINMAX(ypix)

    coords = [[x01[0], y01[0]], [x01[0], y01[1]], $
              [x01[1], y01[1]], [x01[1], y01[0]]]

    IF nvalid EQ 0 THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Region outside image, no data to analyse'
        GOTO, QUIT
    ENDIF

    IF nvalid LT npix THEN BEGIN
        box = box[good]
        x = xpix[good]  &  y = ypix[good]
    ENDIF ELSE BEGIN
        x = xpix  &  y = ypix
    ENDELSE
    good = 0

    line = [' ', STRING( npix, FORMAT = $
"('IMSTATS: Image statistics for region with ',I8,' pixels')")]
ENDELSE

IF llbb THEN BEGIN
    XY2AD, x, y, astrom, ll, bb

    valid = WHERE(FINITE(ll) AND FINITE(bb), nvalid)
ENDIF ELSE valid = LINDGEN(nvalid)

IF nvalid GT 0 THEN good = WHERE(FINITE(box[valid]), ngood) ELSE ngood = 0

line = [line, STRING(nvalid, ngood, FORMAT = $
"( '         Includes ',I6,' sky pixels of which',I6,' contain valid data')")]

IF ngood EQ 0 THEN BEGIN
    line = [line, 'IMSTATS: No further analysis possible!']
    GOTO, QUIT
ENDIF

line = [line, ' ']

boxmax = MAX(box, /NAN)
maxpix = WHERE(box EQ boxmax, nmax)
IF nmax GT 1 THEN BEGIN
    line = [line, 'IMSTATS: Multiple maxima! Quoting the first found']
    maxpix = maxpix[0]
ENDIF
boxmin = MIN(box, /NAN)
minpix = WHERE(box EQ boxmin, nmin)
IF nmin GT 1 THEN BEGIN
    line = [line, 'IMSTATS: Multiple minima! Quoting the first found']
    minpix = minpix[0]
ENDIF

maxcoor = [x[maxpix], y[maxpix]]
mincoor = [x[minpix], y[minpix]]
x = 0  &  y = 0

;fmax = numunit(boxmax, unit, PRECISION = 4)
;fmin = numunit(boxmin, unit, PRECISION = 4)

IF llbb THEN BEGIN
    XY2AD, maxcoor[0], maxcoor[1], astrom, lmax, bmax
    XY2AD, mincoor[0], mincoor[1], astrom, lmin, bmin


    hpfm   = "(9(' '),A,':',E12.3,' at pixel (',I5,',',I5,')," + $
    " long:',F8.3,' lat:',F7.3)"

    line = [line, STRING('Maximum', boxmax, maxcoor, lmax, bmax, $
                          FORMAT = hpfm)]

    wcs = STRMID(astrom.CTYPE[0],5)
    hpx = wcs EQ 'HPX'
    IF hpx THEN BEGIN
        nside = T[1] / 5
        theta = (90d0 - bmax)/r2d  & phi = lmax/r2d
        ang2pix_ring, nside, theta, phi, maxpixr
        ang2pix_nest, nside, theta, phi, maxpixn
        theta = (90d0 - bmin)/r2d  & phi = lmin/r2d
        ang2pix_ring, nside, theta, phi, minpixr
        ang2pix_nest, nside, theta, phi, minpixn

        hpfmh  = "(13(' '),'at HEALPix pixels: ',I9,' (nest); '" + $
                                               ",I9,' (ring)')"

        line = [line, STRING(maxpixn, maxpixr, FORMAT = hpfmh)]
    ENDIF

    line = [line, '', STRING('Minimum', boxmin, mincoor, lmin, bmin, $
                         FORMAT = hpfm)]

    IF hpx THEN line = [line, STRING( minpixn, minpixr, FORMAT = hpfmh)]

ENDIF ELSE BEGIN
    normfm = "(9(' '),A,':',E12.3,' at pixel (',I5,',',I5,')')"
    line = [line, STRING( 'Maximum', boxmax, maxcoor, FORMAT=normfm)]
    line = [line, STRING( 'Minimum', boxmin, mincoor, FORMAT=normfm)]
ENDELSE


box = box[valid]
box = box[good]
boxmean = MEAN(box)
boxmed  = MEDIAN(box)
boxstd  = STDDEV(box)

line = [line, ' ', STRING( boxmean, boxmed, FORMAT = $
                "(9(' '),'Mean:        ',E12.3,'  Median:',E11.3)")]
line = [line, STRING( boxstd, unit, FORMAT = $
                      "(9(' '),'Standard Dev:',E12.3,'  Unit:   ',A)")]

; Check whether units are flux per pixel
start = STREGEX(unit,'pixel',/FOLD_CASE)
IF start GE 0 THEN len = 5 ELSE BEGIN
    start = STREGEX(unit,'pix',/FOLD_CASE)
    IF start GE 0 THEN len = 3
ENDELSE
IF start GE 0 THEN BEGIN        ; Make sure units are per pixel
                                ; Harder than it looks!
    solidus = STREGEX(unit,'/')
    perpix = solidus GE 0 AND solidus LT start
    bytes = BYTE(unit)

    brace = STREGEX(unit,'(')
    IF brace GT -1 THEN BEGIN
        ecarb = STREGEX(unit,')')
    ENDIF
ENDIF


IF llbb THEN BEGIN
    IF hpx THEN BEGIN
        pixarea = !dpi / (3L*nside^2)
        pixadif = 0d0
    ENDIF ELSE BEGIN
        area = DBLARR(4)
        FOR ipix = 0,3 DO area[ipix] = get_pix_area(coords[*,ipix], astrom)
        pixarea = MEAN(area)
        pixadif = MAX(area) - MIN(area)
    ENDELSE
    flux = boxmean * ngood * pixarea
    fracerr = pixadif / pixarea

    fluxstr = numunit(flux, unit, PRECISION = 4)
    line = [line, STRING( fluxstr, FORMAT = $
                  "(9(' '),'Integrated flux density: ',A,' sr')")]
    IF fracerr GT 1e-4 THEN line = [line, STRING( 100.*fracerr, FORMAT = $
              "(9(' '),'Pixel area changes by',F6.2,'% across region;')"), $
              '         representative value used.']
ENDIF

zeroes = WHERE(box EQ 0.0, nzero)
IF nzero GT 0 THEN BEGIN
    line = [line, ' ', STRING( nzero, FORMAT = $
                    "(9(' '),'Found',I6,' exact zeroes. Excluding them:')")]
    good = WHERE(box NE 0.0, ngood)
    IF good[0] NE -1 THEN BEGIN
        box = box[good]
        boxmean =  MEAN(box)
        boxmed  = MEDIAN(box)
        boxstd  = STDDEV(box)

        line = [line, STRING( boxmean, boxmed, FORMAT = $
               "(9(' '),'Mean:        ',E12.3,'  Median:',E12.3)" )]
        line = [line, STRING( boxstd, FORMAT = $
                                    "( 9(' '),'Standard Dev:',E12.3)")]
        IF llbb THEN BEGIN
            flux = boxmean * ngood * pixarea
            line = [line, STRING( flux, unit, FORMAT = $
                  "(9(' '),'Integrated flux density:',E12.3,' ',A,' sr')")]
        ENDIF
    ENDIF ELSE line = [line, STRING('No pixels left!',FORMAT = "(9(' '),A)")]
ENDIF

QUIT:

IF line[0] NE '' THEN BEGIN
    line = [line,' ']
    PRINT, line, FORMAT = "(A)"  ; Should put each string on its own line
    IF N_ELEMENTS(lun) NE 0 THEN PRINTF, lun, line, FORMAT="(A)"
ENDIF

RETURN, coords

END

