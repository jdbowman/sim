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
FUNCTION grid_header, header, nside, proj, ndims, dims, blc, trc, llshift, $
                      UNIT = unit
;
; J. P. Leahy 2008
;
; Updates fits header from HParray to HPX/XPH projection. Keywords
; assume that the gridded data is stored in the primary HDU.
;
; Standard keywords that will be deduced from the data by fitswriters
; (e.g. BITPIX) are not updated, except for NAXIS, NAXISn, which are
; used for astrometry.
;
; Inputs:
;   header:  FITS header array
;   nside:   nside of FITS image
;   proj:    "Projection" i.e. ordering of grid: 'GRID' = HPX
;            'NPOLE' or 'SPOLE' for XPH = butterfly.
;   ndims:   Number of extra dimensions
;   dims:    Size of each extra dimension (if any)
;   blc:     Bottom left corner relative to full-sky grid
;   trc:     Top right corner relative to full-sky grid.
;   llshift: number of panels shifted in negative longitude
;   unit:    Brightness unit. If specified, overrides header value, if
;            any.
;
; Returns: updated header
;
COMPILE_OPT IDL2, HIDDEN

coord = [['GLON', 'ELON', 'RA--', 'SLON', 'TLON', 'XLON'],$
         ['GLAT', 'ELAT', 'DEC-', 'SLAT', 'TLAT', 'XLAT']]

sys = SXPAR(header,'COORDSYS', COUNT=count)
IF count EQ 0 THEN sys = SXPAR(header,'SKYCOORD', COUNT=count)

IF count GT 0 THEN BEGIN
    sys1 = STRMID(STRUPCASE(sys),0,1)
    CASE sys1 OF 
        'G': systype = 0        ; Galactic
        'E': systype = 1        ; Ecliptic
        'C': systype = 2        ; Celestial = Equatiorial = RA/Dec
        'Q': systype = 2        ; as above
        'S': systype = 3        ; Supergalactic
        'T': systype = 4        ; Terrestrial
        ELSE: BEGIN
            systype = 5
            PRINT, 'GRID_HEADER: unrecognised coordinate system code: ' + sys
        END
    ENDCASE
ENDIF ELSE systype = 5

IF STRCMP(proj, 'GRID', 4, /FOLD_CASE) THEN BEGIN
    projcode = '-HPX'
    npanel   = 5
    crval2 = 0d0
ENDIF ELSE BEGIN
    projcode = '-XPH'
    npanel   = 4
    crval2 = proj EQ 'NPOLE' ? 90d0 : -90d0
ENDELSE
ngrid = nside*npanel
IF N_ELEMENTS(llshift) EQ 0 THEN llshift = 0
IF N_ELEMENTS(blc)     LT 2 THEN blc = [0,0]
IF N_ELEMENTS(trc)     LT 2 THEN trc = [ngrid, ngrid] - blc
nx = trc[0] - blc[0]
ny = trc[1] - blc[1]

; Remove out of date info from header
junk  = SXPAR(header, 'TFORM*', COUNT = count)
IF count GT 0 THEN BEGIN
    istr = STRTRIM( STRING(INDGEN(count)+1), 2)
    SXDELPAR, header, 'TFORM'+istr
ENDIF
SXDELPAR, header, 'PIXTYPE'     ; data is no longer a pixel array
                                ; Other HEALPix keywords left to
                                ; provide history

; Try to recover brightness unit
IF ~N_ELEMENTS(unit) THEN BEGIN
    unit = SXPAR(header, 'TUNIT*', COUNT=is_unit)
    multi = is_unit GT 1
    IF multi THEN BEGIN
        test = unit[0]
        null = WHERE(test EQ unit,count)
        unit = test
        IF count EQ is_unit THEN multi = 0B ; all same unit
    ENDIF ELSE IF ~is_unit THEN unit = SXPAR(header,'BUNIT', COUNT=is_unit)

    IF ~is_unit THEN unit = 'unknown'
ENDIF

; naxis = SXPAR(header, 'NAXIS') ; Should be meaningless 2 for binary table
naxis = 1
types = SXPAR(header,'CTYPE*', COUNT = count)
recorded = count GT 0
IF recorded THEN BEGIN
    naxis = count
    nax   = SXPAR(header,'NAXIS*', COUNT = c0)
    rpixs = SXPAR(header,'CRPIX*', COUNT = c1)
    delts = SXPAR(header,'CDELT*', COUNT = c2)
    rvals = SXPAR(header,'CRVAL*', COUNT = c3)
    recorded = MAX([c0,c1,c2,c3]) GT 0
ENDIF

IF recorded THEN BEGIN ; Re-write axes, shifted by one.
    fcode = "(A,I1)"
    FOR i=0,ndims-1 DO BEGIN
        i2 = i+3
        SXADDPAR, header, STRING('NAXIS', i2, FORMAT = fcode), nax[i2-1]
        SXADDPAR, header, STRING('CTYPE', i2, FORMAT = fcode), types[i2-1]
        SXADDPAR, header, STRING('CRPIX', i2, FORMAT = fcode), rpixs[i2-1]
        SXADDPAR, header, STRING('CDELT', i2, FORMAT = fcode), delts[i2-1]
        SXADDPAR, header, STRING('CRVAL', i2, FORMAT = fcode), rvals[i2-1]
    ENDFOR
ENDIF 
; Standard header parameters for HPX/XPH projection
SXADDPAR, header, 'NAXIS', 2 + naxis - 1
SXADDPAR, header, 'NAXIS1', nx
SXADDPAR, header, 'NAXIS2', ny
SXADDPAR, header, 'CTYPE1', coord[systype,0] + projcode
SXADDPAR, header, 'CTYPE2', coord[systype,1] + projcode
SXADDPAR, header, 'CRPIX1', npanel*nside/2d0 + 0.5d0 - blc[0]
SXADDPAR, header, 'CRPIX2', npanel*nside/2d0 + 0.5d0 - blc[1]
SXADDPAR, header, 'CDELT1', -90d0/nside
SXADDPAR, header, 'CDELT2',  90d0/nside
SXADDPAR, header, 'CRVAL1', llshift*90d0
SXADDPAR, header, 'CRVAL2', crval2

IF STRCMP(proj, 'GRID', 4, /FOLD_CASE) THEN BEGIN
    SXADDPAR, header, 'PC1_1',  0.5d0 ; Pixel grid at 45 deg relative
    SXADDPAR, header, 'PC1_2',  0.5d0 ; to the natural (x,y) coords
    SXADDPAR, header, 'PC2_1', -0.5d0 ; for the projection.
    SXADDPAR, header, 'PC2_2',  0.5d0 ;
    SXADDPAR, header, 'PV2_1', 4
    SXADDPAR, header, 'PV2_2', 3
ENDIF

IF ~multi THEN SXADDPAR, header, 'BUNIT', unit[0] ; otherwise use TUNITi

ttype = SXPAR(header, 'TTYPE*', COUNT = count)
IF STRCMP(ttype[0], 'PIXEL', 5, /FOLD_CASE) THEN BEGIN ; Remove pixel column
    tunit = SXPAR(header, 'TUNIT*', COUNT = cunit)
    IF cunit LT count THEN tunit = [tunit, STRARR(count-cunit)]
    FOR ii = 0, count-2 DO BEGIN
        istr = STRTRIM( STRING(ii+1), 2)
        SXADDPAR, header, 'TTYPE'+istr, ttype[ii+1]
        SXADDPAR, header, 'TUNIT'+istr, tunit[ii+1]
    ENDFOR
    istr = STRTRIM( STRING(count), 2)
    SXDELPAR, header, 'TTYPE'+istr
    SXDELPAR, header, 'TUNIT'+istr
ENDIF

; Non-standard header for internal use
SXADDPAR, header, 'PROJECTI', proj

SXADDHIST, 'Converted from HEALPix array by IDL hpgrid '+SYSTIME(), header

RETURN, header
END

