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
; Module ximview: widget-based image inspection tool
;
; J. P. Leahy 2008
;
; This file contains an IDL documentation block after the ximview
; procedure is declared.
;
; Where is the data?
;  .TOP uservalue = "state" structure: contains widget IDs for key
;  widgets, and various parameters which are fixed or change rarely
;  during execution: values of input parameters (except the actual image),
;  widget geometry, SIZE output for input map (.IMSIZE),
;  colour indices, etc.
;
;  state.TABARR is a pointer to a structure array, each member of
;  which contains the data for each tabbed screen, including the
;  pointer to the byte-scaled image, widget IDs and scale information,
;  and usually the pointer to the original image.
;
;  state.TABS uservalue is a structure ("mode") containing all the
;  variable geometry parameters, switches for different operation modes,
;  zoom details.
;
;  state.LABEL uservalue is the image itself, unless supplied as
;  pointers to heap variables, in which these are stored via the
;  .IM_PTR tag in the relevant tab structure.
;
;  state.OLD_GRAPH is a pointer to a structure containing the
;  former graphics state (restored when the widget loses keyboard
;  focus or is killed.)
;
;  the .LUT tag in the tab structures is a pointer to a structure
;  containing the R,G,B arrays for the colour table asssociated with
;  the tab, along with the colour indices for absent pixels and line
;  graphics.
;--------------------------------------------------------------------------
;
; Utility routines:
;
FUNCTION colour_schemes, index
; Returns the name of the colour scheme associated with index.
;
COMPILE_OPT IDL2, HIDDEN

schemes = ['Rainbow', 'Heat', 'Blue-yellow-white', 'Greyscale', $
           'Red-black-blue', 'Cyclic']

IF index EQ -1 THEN RETURN, schemes ELSE RETURN, schemes[index]
END

FUNCTION scale_funs, index
; Returns the name of the scaling function with index.
;
COMPILE_OPT IDL2, HIDDEN

trfunc = ['Linear','Asinh','Sqrt', 'Hist eq']

IF index EQ -1 THEN RETURN, trfunc ELSE RETURN, trfunc[index]
END

PRO set_scale, scale_pars, str
;
; Sets RANGE in tab structure str based on scale_pars and the colour table
; in effect.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global, windev, redraw_req, colmap, badcol, syscol

mindat = scale_pars[0]
maxdat = scale_pars[1]
zero   = scale_pars[2]
sdev   = scale_pars[3]

; Find which way up the intensity scale is:
ratio = (maxdat-zero) / ABS(mindat - zero)
top = 0
IF ratio GT 5 THEN top = 1 ELSE IF ratio LT 0.2 THEN top = -1
CASE str.TRFUNC OF
    1: range = [mindat, maxdat] ; Asinh scale
    3: BEGIN                    ; Histogram equalization
                                ; Avoid under-digitization on initial
                                ; analysis: of 3000 bins we want at
                                ; least 15 to cover +/- 1 sigma.
        range = ratio GT 1 ? mindat + [0, 400*sdev] : maxdat - [0, 400*sdev]
    END
    ELSE: BEGIN                 ; Func 0 or 2: linear or sqrt
        offsymm = [25., 12., 25., 12., 9., 25., 50.]
        IF str.TRFUNC EQ 2 THEN offsymm = (2./3.)*offsymm^2
        offasym = offsymm - 3.
        range = top*offasym[str.COLTAB] + [-1, 1]*offsymm[str.COLTAB]
        range = range*sdev + zero
    END
ENDCASE
range = (range > mindat) < maxdat
beta = [2., 2., 2., 2., 1.5, 2, 2]

str.ABSRANGE = [mindat, maxdat]
str.RANGE    = range
str.BETA     = sdev*beta[str.COLTAB]
str.ZERO     = zero
str.MODE     = zero
str.SDEV     = sdev
str.IZERO    = scale_image(zero, range, str.WRAP, badcol, str.BOT, str.TOP, $
                           str.TRFUNC, zero, str.BETA)

END

PRO update_labels, state, tabpref, newname, newtabs, tab_lab, name
;
; Updates top title and tab labels when new data tabs are added
;
; Inputs:
;    state:   structure with main global variables
;    tabpref: User-supplied prefix for each tab label
;    newname: Name string built by parse_input
; Input/output
;    newtabs: On input, array of structures to describe new tabs.
;             on output, array describing all tabs.
;    tab_lab: On input, array of tab labels from parse_input
;             on output, updated
; Output:
;    name:    New top title
;
COMPILE_OPT IDL2, HIDDEN

oldtabs = *state.TABARR
ntab = N_ELEMENTS(oldtabs)
ncol = N_ELEMENTS(newtabs)
tab_arr = [oldtabs, newtabs]
oldtabs = 0  &  newtabs = 0

            ; Get labels
oldlabs = get_tab_uvals(state.TABS)

IF ~STREGEX(tabpref, '^ *$', /BOOLEAN) THEN BEGIN
                                ;   Check to see if we have default labels
    IF STRCMP(tab_lab[0], 'Map 1') THEN tab_lab[*] = tabpref $
    ELSE tab_lab = tabpref + ' ' + tab_lab
    tab_lab = [oldlabs, tab_lab]
    donew = 0B
ENDIF ELSE BEGIN
    tab_lab = [oldlabs, tab_lab]
    donew = 1B
ENDELSE
                                ; Work out new plot title and tab
                                ; prefix from name structure
files  = tab_arr.FILE
files  = files[ UNIQ( files, BSORT(files) ) ]
nfiles = N_ELEMENTS(files)

IF nfiles GT 1 THEN BEGIN
                                ; Eliminate repetition of 'band'
    freqs = tab_arr.FREQ
    freqs = freqs[ UNIQ( freqs, BSORT(freqs) ) ]
    nfreq = N_ELEMENTS(freq)
    IF nfreq GT 1 THEN FOR itab = 0, ntab+ncol-1 DO BEGIN
        freq = STRSPLIT(freqs[itab],'band', /REGEX, /EXTRACT)
        tab_arr[itab].FREQ = freq
    ENDFOR

    newname = STRJOIN(files,', ') + ': '
    len1 = STRLEN(newname)
    start = WHERE( STRCMP( TAG_NAMES(tab_arr), 'file', /FOLD_CASE) )
    new_lab = STRARR(ntab+ncol)
    nntail = ''
    FOR i = 1,6 DO BEGIN
        param = tab_arr.(i + start)
        par   = param[ UNIQ( param, BSORT(param) ) ]
        npar  = N_ELEMENTS(par)
        IF npar EQ 1 THEN nntail += par ELSE IF donew THEN BEGIN
            new_lab = new_lab + ' ' + STRTRIM(param,2)
            tab_arr.(i+start) = ''
        ENDIF
    ENDFOR
    IF donew THEN tab_lab = new_lab + ' ' + tab_lab
    len2 = STRLEN(nntail)

    IF len1+len2 GT 80 THEN newname = files[0]+' and others: '
    newname += nntail
ENDIF             ; Else newname supplied by parse_input should be OK.
name = newname

tab_lab = STRTRIM(STRCOMPRESS(tab_lab), 2)
newtabs = tab_arr

END

FUNCTION match_tab, state, T, proj, is_astrom, astrom, csystem
;
; Check that newly loaded data matches other tabs in size and astrometry
;
; Inputs
;   state:     State structure describing existing tabs
;   T:         SIZE array for new data
;   proj:      HEALPix projection (if any)
;   is_astrom: General astrometry is available
;   astrom:    Astrolib WCS astrometry structure
;   csystem:   String describing coordinate system
;
COMPILE_OPT IDL2, HIDDEN

bad1 = ~ARRAY_EQUAL(T[1:2], state.IMSIZE[1:2])
IF bad1 THEN MESSAGE, /INFORMATIONAL, $
  'New array size does not match data already loaded'

bad2 = proj NE state.PROJ
IF bad2 THEN  MESSAGE, /INFORMATIONAL, $
  'New array HP projection:' + proj + ' does not match data already loaded'

bad3 = is_astrom NE state.IS_ASTROM
IF bad3 THEN  MESSAGE, /INFORMATIONAL, 'New array astrom usage:' $
    + is_astrom + ' does not match data already loaded'

IF is_astrom THEN BEGIN         ; check astrom structure
    good = 1
    ostrom = *state.ASTROM
    FOR i=0,4 DO good *= ARRAY_EQUAL(astrom.(i), ostrom.(i))
    good *= astrom.LONGPOLE EQ ostrom.LONGPOLE
    good *= astrom.LATPOLE  EQ ostrom.LATPOLE
    good *= ARRAY_EQUAL(astrom.PV2, ostrom.PV2)
    IF ~ARRAY_EQUAL(astrom.CTYPE, ostrom.CTYPE) THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Warning: new data coordinate type: ' + $
          STRJOIN(astrom.CTYPE,', ')
        MESSAGE, /INFORMATIONAL, '         does not match previous:  ' + $
          STRJOIN(ostrom.CTYPE,', ')
    ENDIF
    bad4 = ~good
    IF bad4 THEN MESSAGE, /INFORMATIONAL, $
        'New array astrometry cannot match data already loaded'
ENDIF ELSE bad4 = 0B

IF csystem NE state.CSYSTEM THEN  BEGIN
    MESSAGE, /INFORMATIONAL, 'Warning: new array coordinate system: '+csystem
    MESSAGE, /INFORMATIONAL, '         does not match previous one: ' + $
      state.CSYSTEM
ENDIF

RETURN, ~(bad1 || bad2 || bad3 || bad4)

END

FUNCTION form_unit, newunit
; Formats unit into centred 12-letter string, unless "unknown"
;
COMPILE_OPT IDL2, HIDDEN

IF ~STRCMP(newunit, 'unknown', 6, /FOLD_CASE) THEN BEGIN
    ounit = newunit
    lunit = STRLEN(ounit)
    pad = 12 - lunit
    p1 = pad/2
    uformat = "("
    p2 = pad - p1 - 1
    IF p2 GT 1 THEN uformat += STRING(p2, FORMAT = "(I1)")+ "(' '),"
    uformat += "' ',A"
    IF p1 GT 0 THEN uformat += STRING(p1, FORMAT = "(',',I1)")+"(' ')"
    uformat += ")"
    mid = STRING(ounit, FORMAT = uformat)
ENDIF ELSE mid = ' Brightness '

RETURN, STRING(mid, FORMAT = "(A12)")
END

PRO fill_tab, tab_arr, namestr, unit, polcode, $
              trfunc, coltab, bot, top, TEMPORARY = temporary
; Fills in entries in the array of structures which describe individual tabs
;
; Inputs:
;   tab_arr:   array of structures to partially fill in
;   namestr:   structure of "name" elements
;   unit:      string or string array containing flux units
;   polcode:   array of polarization codes
;   temporary: set to 1 if data can be deleted when the program closes
;   trfunc:    Graphics transfer function
;   coltab:    Colour table
;   bot:       Colour index that maps to bottom of range.
;   top:       Colour index that maps to top of range.
;
COMPILE_OPT IDL2, HIDDEN

temporary = KEYWORD_SET(temporary)

tab_arr.TEMPORARY  = temporary

; Labelling details
tab_arr.UNIT       = unit
tab_arr.POLCODE    = polcode
tab_arr.FREQCODE   = STRTRIM(namestr.FREQ,2) ; Retained permanently
tab_arr.FILE       = namestr.FILE
tab_arr.TELESCOPE  = namestr.TELESCOPE
tab_arr.INSTRUMENT = namestr.INSTRUMENT
tab_arr.CREATOR    = namestr.CREATOR
tab_arr.OBJECT     = namestr.OBJECT
tab_arr.FREQ       = namestr.FREQ
tab_arr.STOKES     = namestr.STOKES

; Transfer fn details:
IF N_PARAMS() GT 5 THEN BEGIN
    tab_arr.TRFUNC = trfunc
    tab_arr.COLTAB = coltab
    tab_arr.BOT    = bot
    tab_arr.TOP    = top
ENDIF

END

PRO fill_gores, nside, imsize, astrptr, temp
;
; Set off-sky pixel values to "absent". In future should deal with
; projections other than HPX.
;
; Inputs:
;   nside: HEALpix parameter
;   imsize: SIZE array for image
;   astrptr: Pointer to astrometry structure
;   temp:  Array of pointers to byte image arrays.
;
COMPILE_OPT IDL2, HIDDEN

ntab = N_ELEMENTS(temp)
astrom = *astrptr

proj = STRMID(astrom.CTYPE[0],5,3)

CASE proj OF
    'HPX': BEGIN
        nsmin1 = nside - 1L
        nx = imsize[1]
        ny = imsize[2]
        blc = 2.5*nside + 0.5 - astrom.CRPIX
        trc = blc + [nx, ny] - 1
        cropped = ~ARRAY_EQUAL(blc, [0,0]) || ~ARRAY_EQUAL(trc, 5*nside*[1,1])

        IF ~cropped THEN BEGIN
            ix  = REBIN( INDGEN(nside), nside, nside)
            iy  = TRANSPOSE(ix)
            idx = TEMPORARY(ix) + nx*TEMPORARY(iy)
            idx = REFORM(idx, nside*nside, /OVERWRITE)
            absar = REPLICATE(!P.background, nside*nside)
        ENDIF
                                ; BLC coordinates for empty panels:
        xc = nside*[2L, 3L, 4L, 3L, 4L, 0L, 4L, 0L, 1L, 0L, 1L, 2L]
        yc = nside*[0L, 0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L]
        x0 = (xc - blc[0]) > 0
        y0 = (yc - blc[1]) > 0
        x1 = (xc+nside - 1L - blc[0]) < (nx - 1)
        y1 = (yc+nside - 1L - blc[1]) < (ny - 1)
        offsets = x0 + nx*y0
        dx = x1 - x0 + 1  &  dy = y1 - y0 + 1
                                ; Locate each missing facet and blank
        FOR i=0,11 DO BEGIN
            IF cropped THEN BEGIN
                nsx = dx[i]  &  nsy = dy[i]
                IF nsx LT 0 || nsy LT 0 THEN CONTINUE
                ix  = REBIN( INDGEN(nsx), nsx, nsy)
                iy  = TRANSPOSE( nsx EQ nsy ? ix : REBIN(INDGEN(nsy),nsy,nsx))
                idx = TEMPORARY(ix) + nx*TEMPORARY(iy)
                npix = nsx*nsy
                idx = REFORM(idx, npix, /OVERWRITE)
                absar = REPLICATE(!P.background, npix)
            ENDIF
            indices = idx + offsets[i]
            FOR j=0,ntab-1 DO BEGIN
                tptr = temp[j]
                IF PTR_VALID(tptr) THEN (*tptr)[indices] = absar
            ENDFOR
        ENDFOR
    END
    ELSE: ; Don't know what to do with other projections
ENDCASE
END

FUNCTION get_log_name
; Finds a unique name for the log file.
;
COMPILE_OPT IDL2, HIDDEN

existing = FILE_SEARCH('ximview_*.log')
IF existing EQ '' THEN logfile = 'ximview_1.log' ELSE BEGIN
                                ; Find maximum index among existing files:
    existing = STRMID(existing, 8) ; removes 'ximview_' prefix
                                ; now remove '.log' suffix:
    nums    = STRSPLIT( STRJOIN(existing), '.log', /REGEX, /EXTRACT)
    newnum  = STRTRIM( STRING( MAX(FIX(nums)) + 1 ), 2)
    logfile = 'ximview_'+newnum+'.log'
ENDELSE

RETURN, logfile
END

PRO check_ast, astrom, header, is_astrom
;
; Checks to make sure astrometry data is intelligible by IDL astrolib
; astrometry routines. Converts two old AIPS cases if possible, as
; described in Section 6 of the WCS paper.
;
COMPILE_OPT IDL2, HIDDEN

d2r = !dpi/180d0

wcs = STRMID(astrom.CTYPE[0], 5)
IF wcs EQ 'GLS' THEN BEGIN
    rotas = SXPAR(header,'CROTA*', COUNT = count)
    IF MAX(ABS(rotas)) EQ 0 THEN BEGIN
        wcs = 'SFL'
        astrom.CRPIX = astrom.CRPIX + astrom.CRVAL/astrom.CDELT
        astrom.CTYPE[0] = STRMID(astrom.CTYPE[0],0,4) + wcs
        astrom.CTYPE[1] = STRMID(astrom.CTYPE[1],0,4) + wcs
    ENDIF
ENDIF

hit = WHERE(wcs EQ  ['DEF', 'AZP', 'SIN', 'ZEA', 'CEA', 'SFL', 'PAR', $
                     'MOL', 'AIT', 'COE', 'BON', 'CSC', 'QSC', 'TSC', $
                     'TAN', 'STG', 'ARC', 'ZPN', 'AIR', 'CYP', 'CAR', $
                     'MER', 'COP', 'COD', 'COO', 'PCO', 'SZP', 'HPX'])
is_astrom = hit[0] GE 0

END

FUNCTION ramp, nlevel, DOWN = down
; Calculates a linear ramp from 0B to 255B over nlevel levels
;
COMPILE_OPT IDL2, HIDDEN

ramp = BYTE( (255L * LINDGEN(nlevel)) / (nlevel-1) )
IF KEYWORD_SET(down) THEN ramp = REVERSE(ramp, /OVERWRITE)
RETURN, ramp
END

PRO ximview_lut, coltab, zero, decomp, bot, top
;
; Set up colour table for XIMVIEW, with special values.
;
; Inputs:
;   coltab:     requested colour table
;   zero:       colour level corresponding to zero
;               (for blue-black-red scale only)
;
; Outputs:
;   decomp:     required value of device 'decomposed' state
;   bot, top:   regular colour scale runs [bot:top]
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

; Revise colour table
black_colour  = [  0B,   0B,   0B]
white_colour  = [255B, 255B, 255B]

black = 0B
ngrey = N_ELEMENTS(syscol)
grey = syscol

 ; Set special grey levels
IF ngrey LT 1 THEN grey = [106B]
IF ngrey LT 2 THEN grey = [grey, 192B]
absent = grey[0]  &  badcol = grey[1]

; If DirectColor mode, leave gaps in colour table for syscols to avoid
; flashing, otherwise just absent & badcol to facilitate RGB images.

grey = redraw_req ? grey[0:1] : [0B, grey]
ngrey = N_ELEMENTS(grey)

null = WHERE(grey LT zero, zgcount)
zero = zero - zgcount

bad_colour    = REPLICATE(badcol,3) ;  Neutral grey, as in MOLLVIEW et al.
absent_colour = REPLICATE(absent,3) ; Darker grey for off-sky pixels

ncol = !D.table_size
white = ncol - 1
line  = white - 1
bot = redraw_req ? 0B : 1B
top = (white - 4) < (white - ngrey - 2)

decomp = 0B

ntop = top + 1
tail = BYTARR(ncol-ntop)

CASE coltab OF
    0: BEGIN                    ; Rainbow pseudo colour
        basic = 39  ; Standard "Rainbow-white". Saturates at level 235
        line_colour   = [255B, 175B, 175B] ;  pink.
        LOADCT, basic, /SILENT  ; Expensive because reads from file!
        TVLCT, rnew,gnew,bnew, /GET
        rnew[top] = white_colour
        gnew[top] = white_colour
        bnew[top] = white_colour
    END
    1: BEGIN ; Heat (black-red-white)
        line_colour   = [175B, 255B, 175B] ; pale green
        rnew = [ramp(172), REPLICATE(255B, ntop-172), tail]
        gnew = [BYTARR(116), ramp(ntop-116), tail]
        bnew = [BYTARR(186), ramp(ntop-186), tail]
    END
    2: BEGIN              ; Blue-yellow-white (colourblind equivalent)
        line_colour   = [255B, 175B, 175B] ;  pink.
        parab = BYTE(255*(1. - ((3.0/ntop)*FINDGEN(2*ntop/3) - 1.0)^2))
        n2 = ntop - (2*ntop/3) - 1
        bnew = [parab, parab[0:n2], tail]
        rnew = [BYTARR(ntop/3), ramp(ntop/4), $
                REPLICATE(255B, ntop - (ntop/3) - (ntop/4)), tail]
        gnew = rnew
    END
    3: BEGIN                    ; Greyscale
        absent_colour = [100B, 100B, 200B] ; blue-grey
        bad_colour    = [200B, 255B, 200B]
        line_colour   = [100B, 255B, 100B] ; green
        rnew = [ramp(ntop), tail]
        gnew = rnew  &  bnew = rnew
    END
    4: BEGIN                    ; Blue-black-red
        line_colour   = [175B, 255B, 175B] ; pale green
        pos_range = top - zero
        neg_range = zero
        peak = pos_range > neg_range
        p2  = peak/2
        p2b = peak - p2

        ramp0 = ramp(p2)
        ramp1 = [ramp0, REPLICATE(255B, p2b)]
        ramp2 = [BYTARR(p2b), TEMPORARY(ramp0)]
        IF neg_range GT 0 THEN BEGIN
            rnew =  [REVERSE(ramp1[0:neg_range-1]), 0B]
            gnew =  [REVERSE(ramp2[0:neg_range-1]), 0B]
            bnew =  [REVERSE(ramp2[0:neg_range-1]), 0B]
        ENDIF ELSE BEGIN
            rnew = [0B]  &  gnew = [0B]  &  bnew = [0B]
        ENDELSE
        IF pos_range GT 0 THEN BEGIN
            rnew = [rnew, ramp2[0:pos_range-1], tail]
            gnew = [gnew, ramp2[0:pos_range-1], tail]
            bnew = [bnew, ramp1[0:pos_range-1], tail]
        ENDIF
    END
    5: BEGIN                    ; Cyclic
        line_colour   = [255B, 255B, 255B] ; white
        null = FLTARR(ncol - ntop)
        one = [REPLICATE(1., ntop), null]
        hue = [360.* FINDGEN(ntop) / top, null] ; hue is in degrees

        TVLCT, hue, one, one, /HSV
        TVLCT, rnew, gnew, bnew, /GET
    END
    ELSE: MESSAGE, 'Unknown colour table'
ENDCASE

DEVICE, DECOMPOSED = decomp

r = rnew  &  g = gnew  &  b = bnew

; Leave gaps in colour table for special greys, including absent and bad:
i1 = line-1  &  i2 = line-2
FOR ii = 0,ngrey-1 DO BEGIN
    i0 = grey[ii]
    r[i0+1:i1] = r[i0:i2]
    g[i0+1:i1] = g[i0:i2]
    b[i0+1:i1] = b[i0:i2]
    r[i0] = i0  &  b[i0] = i0  &  g[i0] = i0
ENDFOR

top  += ngrey
zero += zgcount
ntop = top + 1
new_colours = TRANSPOSE([[absent_colour], [bad_colour], [line_colour], $
                         [white_colour]])
cols = [absent, badcol, line, white]
r[cols] = new_colours[*,0]
g[cols] = new_colours[*,1]
b[cols] = new_colours[*,2]

TVLCT, r, g, b

!P.color      = line
!P.background = absent

END

PRO get_centre, zoom_factor, xhalf, yhalf
;
; Finds effective central pixel on view window
;
COMPILE_OPT IDL2, HIDDEN

xhalf = !D.x_vsize / 2
yhalf = !D.y_vsize / 2

zfac = ROUND(zoom_factor)
IF zfac GT 1 THEN BEGIN ; make xhalf, yhalf land on a pixel centre
    nbigpix = xhalf/zfac
    xhalf = zfac*nbigpix + zfac/2
    nbigpix = yhalf/zoom_factor
    yhalf = zfac*nbigpix + zfac/2
ENDIF

END

PRO overview, temp, zoom_factor, xpix, ypix, xhalf, yhalf, $
              x_centre, y_centre, resamp, corner, NOCENTRE = no_centre
;
; Plots overview with FOV of zoomed-in field outlined.
;
; Inputs:
;     temp:         pointer to Byte array with image data
;                   (in RGB mode, array of 3 pointers)
;     zoom_factor:  what it says
;     xpix, ypix:   Image coordinates at which to centre cursor
;     xhalf, yhalf: Display coords for display centre
;     x_centre, y_centre:  Image coords for display centre
;     no_centre:    Don't move the cursor to the centre pixel.
;
; Outputs:
;     resamp:       Resample factors [x,y]
;     corner:       Coordinates of BLC on display
;
COMPILE_OPT IDL2, HIDDEN

do_centre = ~KEYWORD_SET(no_centre)

DEVICE, /CURSOR_CROSSHAIR

nchan = N_ELEMENTS(temp)
FOR i = 0, nchan-1 DO IF PTR_VALID(temp[i]) THEN T = SIZE(*(temp[i]))

IF N_ELEMENTS(T) EQ 0 THEN MESSAGE, $
  'Internal error: no valid pointers received'

                                ; Box showing current FOV
xlo  = x_centre - xhalf/zoom_factor
xhi  = xlo  + !D.x_vsize/zoom_factor
ylo  = y_centre - yhalf/zoom_factor
yhi  = ylo  + !D.y_vsize/zoom_factor
xbox = [xlo,xhi,xhi,xlo,xlo]
ybox = [ylo,ylo,yhi,yhi,ylo]

ERASE
resamp = REPLICATE( divup(T[1],!D.x_vsize) > divup(T[2],!D.y_vsize), 2)
IF resamp[0] EQ 1 THEN BEGIN  ; Try zooming in
    resamp = REPLICATE( (!D.x_vsize/T[1]) < (!D.y_vsize/T[2]) , 2)

    xysiz   = T[1:2] * resamp
    xcorner = (!D.x_vsize - xysiz[0] )/2
    ycorner = (!D.y_vsize - xysiz[1] )/2
    xbox    = xbox*resamp[0] + xcorner
    ybox    = ybox*resamp[1] + ycorner
    xtv     = xpix*resamp[0] + xcorner
    ytv     = ypix*resamp[1] + ycorner

    block = xysiz[0]*xysiz[1]
    image = BYTARR(nchan*block)

    FOR i=0,nchan-1 DO IF PTR_VALID(temp[i]) THEN image[i*block] = $
      REFORM(REBIN(*temp[i], xysiz[0], xysiz[1], /SAMPLE), block)

    resamp = 1.0 / resamp
ENDIF ELSE BEGIN
    xcorner = (!D.x_vsize - T[1] / resamp[0])/2
    ycorner = (!D.y_vsize - T[2] / resamp[1])/2
    xbox    = xbox/resamp[0] + xcorner
    ybox    = ybox/resamp[1] + ycorner
    xtv     = xpix/resamp[0] + xcorner
    ytv     = ypix/resamp[1] + ycorner

    xysiz   = DIVUP(T[1:2], resamp)
    block = xysiz[0]*xysiz[1]
    image = BYTARR(nchan*block)

    FOR i=0,nchan-1 DO IF PTR_VALID(temp[i]) THEN image[i*block] = $
      REFORM((*temp[i])[0:*:resamp[0], 0:*:resamp[1]], block)
ENDELSE

image = REFORM(image, xysiz[0], xysiz[1], nchan, /OVERWRITE)

IF nchan EQ 3 THEN BEGIN
    DEVICE, /DECOMPOSED
    TV, image, xcorner, ycorner, TRUE = 3
    DEVICE, DECOMPOSED = 0
ENDIF ELSE TV, image, xcorner, ycorner

                                ; Plot box
PLOTS, xbox, ybox, /DEVICE
                                ; Put cursor on current point
IF do_centre THEN TVCRS, xtv, ytv

corner = [xcorner, ycorner]

END

FUNCTION invert_scale, ivalue, str
;
; Finds the image value corresponding to given byte value (if possible).
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

value = ivalue
bad = WHERE(ivalue EQ syscol, count)
IF bad[0] GT -1 THEN value[bad] =  !values.F_NAN
IF count EQ N_ELEMENTS(ivalue) THEN RETURN, value

bot = str.BOT
top = str.TOP
ngrey = N_ELEMENTS(syscol)
FOR ig=0,ngrey-1 DO BEGIN
    idx = WHERE(ivalue GT syscol[ig])
    IF idx[0] NE -1 THEN value[idx] -= 1
ENDFOR
IF bot GT 0 THEN value -= bot
top = top - ngrey - bot

r1 = str.RANGE[0]  &  r2 = str.RANGE[1]
CASE str.TRFUNC OF
    0: scale = top / (str.RANGE[1] - str.RANGE[0])
    1: BEGIN
        asr = ASINH((str.RANGE - str.ZERO)/str.BETA)
        scale = top / (asr[1] - asr[0])
        r1 = str.ZERO
    END
    2:  scale = top / SQRT(r2 - r1)
    3: MESSAGE, /INFORMATIONAL, 'Cannot invert Histogram equalization'
ENDCASE

value /= scale
CASE str.TRFUNC OF
    0: ; No action for linear
    1: BEGIN  ; Asinh
        expval = EXP(value + asr[0])
        value = str.BETA*(expval - 1d0/expval)/2d0
    END
    2: value = value^2
    3: value = !values.F_NAN
ENDCASE
value += r1

RETURN, value

END

PRO set_colour_bar, str
;
; Draws an intensity scale bar in the current graphics window.
;
; Input:
;    str: structure describing tab to label
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

CASE str.COLLAB[0] OF
    'mono': GOTO, PSEUDOCOLOUR  ; carry on
    'HSV':  hsv_label, str
    ELSE:   rgb_label, str
ENDCASE
RETURN

PSEUDOCOLOUR:

WSET, str.SCALE_INDEX
ERASE

nxpix = !D.x_vsize - 40  &  nypix = 10

absr  = str.ABSRANGE     &  range  = str.RANGE
wrap  = str.WRAP         &  trfunc = str.TRFUNC
bot   = str.BOT          &  top    = str.TOP
zero  = str.ZERO         &  beta   = str.BETA

absent = !P.background

temp = [[absr], [range]]
r1 = MIN(temp[0,*])
r2 = MAX(temp[1,*])
IF wrap LE 0 THEN r1 = range[0]
IF wrap EQ 0 THEN r2 = range[1]

scale = (r2 - r1) * FINDGEN(nxpix) / (nxpix-1) + r1
scale = scale_image(scale, range, wrap, badcol, bot, top, 0)

x0 = 20               &  x1 = x0 + nxpix - 1
y0 = !D.y_vsize - 13  &  y1 = y0 + nypix - 1
TV, REBIN(scale, nxpix, nypix), x0, y0

                                ; Now the hard part: label the scale
position = CONVERT_COORD([x0,x1], [y0,y1], /DEVICE, /TO_NORMAL)

xold = !X  &  yold = !Y
!Y.window = [position[1,0], position[1,1]]
!X.window = [position[0,0], position[0,1]]

dbyte = top - bot
IF trfunc EQ 1 || trfunc EQ 2 THEN BEGIN
                                ; Choose tick values for non-linear scales
    iz2 = scale_image(0.0, range, wrap, badcol, bot, top, $
                       trfunc, zero, beta)
    rough = iz2[0] + (LINDGEN(13) - 6)*(dbyte / 5)
                                ; Purge out-of-range values:
    good = WHERE(rough GE bot AND rough LE top, ngood)
    rough = rough[good]
    idzero = WHERE(rough EQ iz2[0])
                                ; convert rough values to round numbers
    imrough = invert_scale(rough, str)
    IF idzero GE 0 && idzero LE ngood-1 THEN BEGIN
        imrough[idzero] = 2.0   ; (avoid zero and unities)
        ntoav = idzero < (ngood - 1 - idzero)
        IF ntoav GT 0 THEN BEGIN
                                ; symmetrize range around zero
            negs = -REVERSE(imrough[idzero-ntoav:idzero-1])
            pos  =  imrough[idzero+1:idzero+ntoav]
            avs  =  0.5*(pos + negs)
            imrough[idzero-ntoav:idzero-1] = -REVERSE(avs)
            imrough[idzero+1:idzero+ntoav] = avs
        ENDIF
    ENDIF
    logs    = FLOOR( ALOG10( ABS(imrough) ) )
    imscale = imrough / 10d0^logs
    leads   = 1.0*ROUND(imscale)
    unities = WHERE(ABS(imscale) GT 9.5 OR ABS(imscale) LT 1.8)
    IF unities[0] NE -1 THEN $
      leads[unities] = ROUND(imrough[unities] / 10d0^(logs[unities]-1)) / 10d0
    ticks = leads * 10d0^logs
                                ; Restore zero to actual zero:
    IF idzero GT 0 && idzero LT ngood - 1 THEN ticks[idzero] = 0.0

                                ; Purge again:
    good = WHERE(ticks GE range[0] AND ticks LE range[1], nticks)
    ticks   = ticks[good]
                                ; Get extra precision for duplicates
    doubles = WHERE(ticks[0:nticks-2] EQ ticks[1:*])
                                ; Find position on linear (ie. colour
                                ; index) scale:
    iticks = scale_image(ticks, range, wrap, badcol, bot, top, $
                         trfunc, zero, beta)
    tnames = STRTRIM(STRING(ticks, FORMAT="(G9.2)"),2)
ENDIF

                                ; Restore original unit
absmax = MAX(ABS(absr)) * str.MULT^2
test = numunit(absmax, str.UNIT, OUT_UNIT = ounit, /FORCE)

!X.title = ounit
CASE trfunc OF
    0: AXIS, XRANGE = [r1, r2], XSTYLE = 1, XTICKLEN = 0.4
    1: BEGIN                    ; Asinh
        AXIS, XRANGE = [bot, top], XSTYLE = 1, XTICKLEN = 0.4, $
          XTICKS = nticks-1, XTICKV = iticks, XTICKNAME = tnames
    END
    2: BEGIN                    ; Sqrt
        AXIS, XRANGE = [bot, top], XSTYLE = 1, XTICKLEN = 0.4, $
          XTICKS = nticks-1, XTICKV = iticks, XTICKNAME = tnames
    END
    3: BEGIN                    ; Histogram equalization.
                                ;God knows what the LUT is, just mark
                                ;the beginning and end
        ticks = [r1, r2]
        AXIS, XRANGE = [r1, r2], XSTYLE = 1, XTICKLEN = 0.4, $
          XTICKS = 1, XTICKV = ticks
    END
ENDCASE
AXIS, XAXIS = 1, XTICKS = 1, XSTYLE = 0
!X = xold  &  !Y = yold

END

PRO rgb_label, str
;
; Draws labels for the R, G, B channels in the space usually used for
; the pseudo-colour scale bar.
;
; Inputs:
;     str:        Structure describing tab
;
COMPILE_OPT IDL2, HIDDEN

WSET, str.SCALE_INDEX
xsize = !D.x_vsize

DEVICE, GET_CURRENT_FONT = oldfont
oldfontcode = !P.font
oldxchar = !D.x_ch_size  &  oldychar = !D.y_ch_size
!P.font = 1
absent = !P.background
DEVICE,/DECOMPOSED
!P.background = 0
ERASE

DEVICE, SET_FONT = "Helvetica Bold", /TT_FONT
DEVICE, SET_CHARACTER_SIZE = [10,16]

strspace = (xsize - 40) / 3
y0 = 0.5*(!D.y_vsize - !D.y_ch_size)
FOR i=0,2 DO XYOUTS, 20 +(i+0.5)*strspace, y0, str.COLLAB[i], $
  COLOR = 255L*(256L^i), /DEVICE, ALIGNMENT = 0.5

                                ; Restore graphics/font state
DEVICE, DECOMPOSED = 0
!P.background = absent    &  !P.font = oldfontcode
DEVICE, SET_CHARACTER_SIZE = [oldxchar, oldychar]
CASE !P.font OF
    -1:                         ; Hershey don't need explicit setting
    0: DEVICE, SET_FONT = oldfont ; Device font
    1: DEVICE, SET_FONT = oldfont, /TT_FONT ; True type
ENDCASE

END

FUNCTION im2tv, xpix, ypix, x_centre, y_centre, xhalf, yhalf, zoom, zfac
; Converts from image pixel to screen pixel displayed by XIMVIEW
;
COMPILE_OPT IDL2, HIDDEN

IF zoom LE 0 THEN BEGIN
    xtv = (xpix - x_centre)/zfac + xhalf
    ytv = (ypix - y_centre)/zfac + yhalf
ENDIF ELSE BEGIN
    xtv = (xpix - x_centre)*zfac + xhalf
    ytv = (ypix - y_centre)*zfac + yhalf
ENDELSE

RETURN, [xtv, ytv]
END

FUNCTION tv2im, xtv, ytv, x_centre, y_centre, xhalf, yhalf, zoom, zfac
; Converts from screen pixel displayed by XIMVIEW to image pixel
;
COMPILE_OPT IDL2, HIDDEN

IF zoom LE 0 THEN BEGIN
    xshift = xtv - xhalf  &  yshift = ytv - yhalf
    xpix = xshift*zfac + x_centre
    ypix = yshift*zfac + y_centre
ENDIF ELSE BEGIN
    xpix = (xtv/zfac) - (xhalf/zfac) + x_centre
    ypix = (ytv/zfac) - (yhalf/zfac) + y_centre
ENDELSE

RETURN, [xpix, ypix]

END

PRO marker, xpix, ypix, x_centre, y_centre, xhalf, yhalf, zoom, zfac
; Draws a marker at *image* pixel xpix, ypix
;
COMPILE_OPT IDL2, HIDDEN

coord = im2tv(xpix, ypix, x_centre, y_centre, xhalf, yhalf, zoom, zfac)
xtv = coord[0]  &  ytv = coord[1]

PLOTS, [-12,-4]+xtv, [  0, 0]+ytv, /DEVICE
PLOTS, [ 12, 4]+xtv, [  0, 0]+ytv, /DEVICE
PLOTS, [  0, 0]+xtv, [-12,-4]+ytv, /DEVICE
PLOTS, [  0, 0]+xtv, [ 12, 4]+ytv, /DEVICE
;phase = 2*!dpi*findgen(21)/20d0
; x = 8*COS(phase) & y = 8*SIN(phase)

PLOTS, [ 8., 7.608, 6.472, 4.702, 2.472, 0.,-2.472,-4.702,-6.472,-7.608,-8., $
        -7.608,-6.472,-4.702,-2.472, 0., 2.472, 4.702, 6.472, 7.608, 8.] $
  + xtv, [ 0., 2.472, 4.702, 6.472, 7.608, 8., 7.608, 6.472, 4.702, 2.472, 0.,$
           -2.472,-4.702,-6.472,-7.608,-8.,-7.608,-6.472,-4.702,-2.472, 0.] $
  + ytv, /DEVICE

END

PRO pix_print, state, log, start
;
; Prints pixel details to readout and possibly terminal and logfile
;
; Inputs:
;         state:       structures with global parameters
;         log:         Write out to terminal and logfile
;         start:       start time of this cycle, for diagnostic
;                      printing.
;
COMPILE_OPT IDL2, HIDDEN

r2d = 180d / !dpi
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
xpix        = mode.XPIX
ypix        = mode.YPIX
zoom_factor = mode.ZOOM_FACTOR

current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
ntab    = WIDGET_INFO(state.TABS, /TAB_NUMBER)
blinking = mode.BLINK && current EQ ntab-1
iscreen = blinking ? WHERE((*state.TABARR).SCREEN EQ (*state.BLINK_SEQ)[0]) $
                   : 0

nside = state.NSIDE
                                ;  Monochrome or RGB?
mono = (*state.TABARR)[iscreen].COLLAB[0] EQ 'mono'
IF mono THEN BEGIN
    im_ptr = (*state.TABARR)[iscreen].IM_PTR
    mult   = (*state.TABARR)[iscreen].MULT

    IF ~im_ptr THEN BEGIN
        WIDGET_CONTROL, state.LABEL, GET_UVALUE = image, /NO_COPY
        imval = image[xpix,ypix]
        WIDGET_CONTROL, state.LABEL, SET_UVALUE = image, /NO_COPY
    ENDIF ELSE imval = (*im_ptr)[xpix,ypix]
    imval = imval * mult
ENDIF

IF state.IS_ASTROM THEN BEGIN ; Get long and lat if available from astrom.
    XY2AD, xpix, ypix, *state.ASTROM, ll, bb
    goodpix = FINITE(ll) && FINITE(bb)
ENDIF ELSE goodpix = 1B

IF nside GT 0 THEN BEGIN        ; Get HEALPix pixel numbers
    is_grid = state.PROJ EQ 'GRID'
    np2 = is_grid ? 2.5 : 2.0
    blc = np2*nside + 0.5 - (*state.ASTROM).CRPIX

    nfull = is_grid ? 5*nside : 4*nside

    IF is_grid THEN BEGIN
        llshift = ROUND((*state.ASTROM).CRVAL[0] / 90d0)
        blc = blc - llshift*nside
    ENDIF

    xshift = (xpix + blc[0])
    yshift = (ypix + blc[1])

    IF is_grid THEN BEGIN
        IF xshift GE nfull || yshift GE nfull THEN BEGIN
            xshift -= state.NS4
            yshift -= state.NS4
        ENDIF ELSE IF xshift LT 0 || yshift LT 0 THEN BEGIN
            xshift += state.NS4
            yshift += state.NS4
        ENDIF
    ENDIF
    idx = grid2hp_index(nside, ROUND(xshift), ROUND(yshift), state.PROJ, $
                        /BOTH, R_OFF = *state.RING0, /SILENT)
    rix = idx[0]  &  nix = idx[1]
    idx = 0
    ngood = nix[0] GE 0
    rgood = rix[0] GE 0
ENDIF

IF state.IS_ASTROM && nside GT 0 THEN BEGIN ; check
    IF goodpix NE ngood OR goodpix NE rgood THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Astrometry mismatch'
        MESSAGE, /INFORMATIONAL, STRING(xpix, ypix, ll, bb, nix, rix, $
                             FORMAT = "('Params:',2I6,2F12.6,2I10)")
    ENDIF
ENDIF ELSE IF ~state.IS_ASTROM && nside GT 0 THEN BEGIN ; Must be XPH order
    PIX2ANG_NEST, nside, nix, theta, phi
    ll = phi*r2d  &  bb = 90d - theta*r2d
    goodpix = ngood
ENDIF

; Now construct the output string:

pos = STRING(xpix, ypix, FORMAT = "(I5,1X,I5,1X)")

IF nside GT 0 THEN pos += goodpix ? STRING(nix, rix, FORMAT = "(2I9)") $
                                  : '                  '
IF state.IS_ASTROM || nside GT 0 THEN BEGIN
    pos += goodpix ? STRING(ll, bb, FORMAT = "(1X,2F9.3)") $
                   : STRING(' No sky pixel', FORMAT = "(A,6(' '))")
ENDIF

IF mono THEN pos += goodpix ? STRING(imval, FORMAT = "(3X,F8.3,'  ')") $
                            : '           ' ELSE pos += ' [RGB image] '

pos += STRING(zoom_factor, FORMAT = "(1X,F6.3,1X)")

dt = SYSTIME(1) - start

IF state.VERBOSE && nix ge 0 THEN BEGIN
;  nest2ring, nside, nix, rix2
;  pos +=STRING(rix-rix2, FORMAT = "(I10)") ; F10.5 for dt 
  pos +=STRING(dt, FORMAT = "(F10.5)")
ENDIF

IF log THEN BEGIN
    PRINT, pos
    PRINTF, state.LOGLUN, pos
ENDIF

WIDGET_CONTROL, state.READOUT, SET_VALUE = pos

END

PRO ximview_resize, top, state
;
; Resizes widget. Called when the program notices
; that the window has been re-sized.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

dx0 = state.NEWWIDTH  - state.XSIZE ; NB really these are SCR_sizes for unix
dy0 = state.NEWHEIGHT - state.YSIZE ;

IF dx0 EQ 0 AND dy0 EQ 0 THEN RETURN

tabarr = state.TABARR
ntab = N_ELEMENTS(*tabarr)
str = (*tabarr)[0]
;
; Check that we have graphics focus and if not, grab it for a bit
IF ~state.FOCUS THEN swap_lut, str, old_graph
;
; Increase graphics size by change in widget size (up to a point):
;
geom = WIDGET_INFO(str.DRAW, /GEOMETRY)
xdraw = geom.DRAW_XSIZE  &  ydraw = geom.DRAW_YSIZE
geom = WIDGET_INFO(state.PAD1, /GEOMETRY)
xp1 = geom.XSIZE
geom = WIDGET_INFO(state.PAD2, /GEOMETRY)
xp2 = geom.XSIZE

IF state.NEWWIDTH EQ -1 THEN BEGIN  ; Code to reset to default size
    newx = 512  &  newy = 512
    dx = newx - xdraw
    xtra = 0  &  ytra = 0
ENDIF ELSE BEGIN
    geom = WIDGET_INFO(state.READLAB, /GEOMETRY)
    len  = geom.XSIZE > 256 ; string length or min allowed view region

    dx = MAX([dx0, len-xdraw, -(xp1+xp2)])
    dy = dy0 > (256-ydraw)
    maxdx = state.MAXWIN[0] - xdraw
    maxdy = state.MAXWIN[1] - ydraw
    dx = dx < maxdx
    dy = dy < maxdy

    newx = xdraw+dx  &  newy = ydraw + dy
ENDELSE

dx1 = dx/2 > (-xp1)

is_unix = STRCMP(!version.OS_FAMILY,'UNIX', 4,/ FOLD_CASE)
IF dx NE 0 AND is_unix THEN WIDGET_CONTROL, top, UPDATE = 0

FOR itab = 0,ntab-1 DO WIDGET_CONTROL, (*tabarr)[itab].DRAW,  $
  DRAW_XSIZE = newx, DRAW_YSIZE = newy

IF dx NE 0 THEN BEGIN
    FOR itab = 0,ntab-1 DO $
      WIDGET_CONTROL, (*tabarr)[itab].SCALE, DRAW_XSIZE = xdraw+dx
    WIDGET_CONTROL, state.PAD1, XSIZE = xp1 + dx1
    WIDGET_CONTROL, state.PAD2, XSIZE = xp2 + dx - dx1
ENDIF

IF mode.BLINK THEN BEGIN
    draw = WIDGET_INFO(mode.BBASE, /CHILD)
    scale = WIDGET_INFO(draw, /SIBLING)
    WIDGET_CONTROL, draw, DRAW_XSIZE = newx, DRAW_YSIZE = newy
    IF dx NE 0 THEN WIDGET_CONTROL, scale, DRAW_XSIZE = newx
ENDIF

IF dx NE 0 AND is_unix THEN WIDGET_CONTROL, top, /UPDATE;
; Record revised geometry:
;
geom = WIDGET_INFO(top,/GEOMETRY)
IF is_unix THEN BEGIN
    state.XSIZE     = geom.XSIZE
    state.NEWWIDTH  = geom.XSIZE
    state.YSIZE     = geom.YSIZE + state.MBAR
    state.NEWHEIGHT = geom.YSIZE + state.MBAR
ENDIF ELSE BEGIN
    state.XSIZE    = geom.XSIZE  &  state.YSIZE     = geom.YSIZE
    state.NEWWIDTH = geom.XSIZE  &  state.NEWHEIGHT = geom.YSIZE
ENDELSE

WIDGET_CONTROL, top, SET_UVALUE = state
;
; Re-draw screens
;
WSET, str.WINDOW

get_centre, mode.ZOOM_FACTOR, xhalf, yhalf
xcent = mode.X_CENTRE  &  ycent = mode.Y_CENTRE
zoomf = mode.ZOOM_FACTOR

null3 = REPLICATE(PTR_NEW(),3)
IF dx NE 0 || mode.OVERVIEW THEN FOR itab = 0,ntab-1 DO BEGIN
    stri = (*tabarr)[itab]
    mono = stri.COLLAB[0] EQ 'mono'

    lutptr = stri.LUT
    IF redraw_req THEN TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
    !P.background = (*lutptr).ABSENT
    !P.color      = (*lutptr).LINE

    IF mode.OVERVIEW THEN BEGIN
        WSET, stri.WINDOW
        tptr = mono ? stri.BYTE_PTR : stri.RGB
                                ; Unlike normal, don't centre cursor
                                ; in middle of window... a bad idea if
                                ; you are currently dragging the edge
                                ; of the window!
        overview, tptr, zoomf, mode.XPIX, mode.YPIX, $
          xhalf, yhalf, xcent, ycent, resamp, corner, /NOCENTRE
    ENDIF
    IF dx NE 0 THEN set_colour_bar, stri
ENDFOR
;
; Re-draw zoomed-in image:
;
IF mode.OVERVIEW THEN BEGIN
    mode.RESAMP = resamp
    mode.CORNER = corner
    mode.NEW_VIEW = 2
ENDIF ELSE BEGIN
    ierr = 0
    zoom = mode.ZOOM       &  zfac = mode.ZFAC
    mode.XPIX = xcent      &  mode.YPIX = ycent

    tptr = (*tabarr).BYTE_PTR
    coord = gscroll(tptr, xcent, ycent, xhalf, yhalf, $
                    zoomf, ierr, 2, done, DO_WRAP = mode.ROLL)
    IF ierr NE 0 && ierr NE 6 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)

                                ; Disable panning if image fits in screen
    bltv = tv2im(0, 0, xcent, ycent, xhalf, yhalf, zoom, zfac)
    trtv = tv2im(!D.x_vsize, !D.y_vsize, xcent, ycent, xhalf, yhalf, $
                 zoom, zfac)
    mode.PAN =  MAX(bltv) GT 0 OR MAX(state.IMSIZE[1:2] - trtv) GT 1

    mode.DONE = done

    IF mode.XPT NE -1 THEN marker, mode.XPT, mode.YPT, xcent, ycent, $
                                   xhalf, yhalf, zoom, zfac

                                ; Request another go if loading not finished
    IF done EQ 0 THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5
ENDELSE
mode.XHALF = xhalf  &  mode.YHALF = yhalf

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode
WIDGET_CONTROL, state.TABS, /SENSITIVE

IF ~state.FOCUS THEN restore_lut, old_graph

END

PRO update_screen, tabarr, itab, mode, done
;  Redraws one tab
;
;  Inputs:
;     tabarr: Pointer to structure array describing tabs
;     itab:   element of array required
;     mode:   Usual mode structure
;     done:   usual flag.
;
COMPILE_OPT IDL2, HIDDEN

str = (*tabarr)[itab]

mono = str.COLLAB[0] EQ 'mono'
set_colour_bar, str

WSET, str.WINDOW

xhalf = mode.XHALF  &  yhalf = mode.YHALF
IF mode.OVERVIEW THEN BEGIN
    done = 1
    tptr = mono ? str.BYTE_PTR : str.RGB
    overview, tptr, mode.ZOOM_FACTOR, mode.XPIX, mode.YPIX, $
      xhalf, yhalf, mode.X_CENTRE, mode.Y_CENTRE, resamp, corner
    mode.RESAMP = resamp  &  mode.CORNER = corner
    mode.NEW_VIEW = 2
ENDIF ELSE BEGIN
    ierr = 0
    zoom = mode.ZOOM       &  zfac = mode.ZFAC
    xcent = mode.X_CENTRE  &  ycent = mode.Y_CENTRE
    mode.XPIX = xcent      &  mode.YPIX = ycent

    tptr = (*tabarr).BYTE_PTR
    coord = gscroll(tptr, mode.XPIX, mode.YPIX, xhalf, yhalf, $
                    mode.ZOOM_FACTOR, ierr, 2, done, DO_WRAP = mode.ROLL)
    IF ierr NE 0 && ierr NE 6 THEN MESSAGE, $
      'GSCROLL error '+STRING(ierr)
ENDELSE

END

PRO parse_header, T, header, column, roll, verbose, $
                 astrom, is_astrom, csystem, proj, ounit, title, nside, ns4
;
; Interprets data read in via parse_input
;
; Inputs:
;   T:       Size array for the input or data
;   header:  FITS header
;   column:  List of columns/slices extracted from original file
;   roll:    Try to force interpretation as HPX grid
;   verbose: Pixel printout will include timing info
;
; Outputs:
;   astrom:    WCS astrometry structure from header
;   is_astrom: Use the astrom structure
;   csystem:   String descriping coordinate system
;   proj:      Projection: 'GRID', 'NPOLE', 'SPOLE', or '' (N/A)
;   ounit:     Intensity unit from header, possibly with revised prefix
;   title:     string to use for pixel printout title
;   nside:     HEALPix nside, if any
;   ns4:       4*nside (offset for rolling)
;
COMPILE_OPT IDL2, HIDDEN

unit = SXPAR(header,'BUNIT', COUNT=got_unit)
IF ~got_unit THEN BEGIN ; Look for TUNITi cards:
    unit = SXPAR(header,'TUNIT*', COUNT=got_unit)
    IF got_unit GE MAX(column) THEN unit = unit[column-1] ELSE $
      IF got_unit GE 1 THEN unit = unit[0] ELSE unit = 'unknown'
ENDIF
unit = STRTRIM(unit,2)

EXTAST, header, astrom, noparam
is_astrom = noparam GT 0
nside = 0B  &  ns4 = 0B  & proj = ''

IF is_astrom THEN BEGIN
                                ; Turn off astrom if WCS type not recognised
    check_ast, astrom, header, is_astrom
                                ; Get overall astrometry details
    ii = WHERE( STRCMP(astrom.CTYPE[0], $
                       ['GLON', 'ELON', 'RA--', 'SLON', 'TLON'],4) )
    IF ii EQ -1 THEN ii = 5
    csystem = (['Galactic', 'Ecliptic', 'Equatorial', 'Supergalactic', $
                'Terrestrial', 'unknown'])[ii]

    equinox = SXPAR(header, 'EQUINOX')
    IF equinox EQ 0 THEN equinox = SXPAR(header,'EPOCH')

    radesys = SXPAR(header, 'RADESYS', COUNT=count)
    IF count EQ 0 THEN radesys = SXPAR(header, 'RADECSYS', COUNT=count)
    IF count EQ 0 THEN BEGIN  ; FITS defaults are quite specific here!
        IF equinox LE 0 THEN radesys = 'ICRS (?)' ELSE $
          radesys = equinox GT 1984.0 ? 'FK5 (?)' : 'FK4 (?)'
    ENDIF

    IF ii LE 3 THEN csystem = csystem+' '+radesys
    CASE STRMID(radesys,3) OF
        'FK4': borj = 'B'
        'FK5': borj = 'J'
        ELSE:  borj = ''
    ENDCASE
    IF equinox GT 0 THEN csystem = csystem + STRING(borj, equinox, FORMAT= $
                                                    "('  Equinox ',A,F6.1)")

                                ; Find out about specific WCS systems
                                ; (a) for rolling through +/- 180 deg
                                ; (b) for getting pixel areas

    wcs = STRMID(astrom.CTYPE[0],5)
    cylindrical = WHERE(['CYP', 'CEA', 'CAR', 'MER'] EQ wcs)
    equiareal = WHERE(['ZEA', 'CEA', 'SFL', 'GLS', 'PAR', $
                       'MOL', 'AIT', 'COE', 'BON', 'QSC'] EQ wcs)
    equiareal = equiareal GT -1
    IF equiareal THEN pixarea = get_pix_area(astrom.CRPIX, astrom) ELSE BEGIN
                                ; maybe pixel area is independent of
                                ; native longitude:
        lon_eq = WHERE(['TAN','STG', 'ARC', 'ZPN', 'AIR', 'CYP', 'CAR', $
                        'MER', 'COP', 'COD', 'COO'] EQ wcs)
        IF ~lon_eq THEN BEGIN
            CASE wcs OF
                'AZP': lon_eq = astrom.PV2[1] EQ 0
                'SZP': lon_eq = ABS(astrom.PV2[2]) EQ 90d0
                'SIN': lon_eq = astrom.PV2[0] EQ 0d0 AND astrom.PV2[1] EQ 0d0
                ELSE: ; it's a mess.
            ENDCASE
        ENDIF
    ENDELSE

; Cylindrical, equiareal, pixarea and lon_eq are not used or returned at
; present... for future expansion.

    hpx = wcs EQ 'HPX'
    IF hpx THEN BEGIN ; Find nside. Don't assume that image is not cropped.
        ad2xy, [90d0, 0d0], [0d0, 0d0], astrom, x, y
        nside = ROUND(x[1] - x[0])
        proj = 'GRID'
                      ; Enable roll if we have full sky coverage:
        roll = T[1] EQ 5*nside
    ENDIF ELSE IF wcs EQ 'XPH' THEN BEGIN
                      ; Temporary... until XPH astrom is coded.
        nside = T[1]/4
        npix = NSIDE2NPIX(nside)
        IF npix EQ -1 THEN nside = SXPAR(header, 'NSIDE')
        roll = 0B
        proj = astrom.CRVAL[1] GT 0. ? 'NPOLE' : 'SPOLE'
        is_astrom = 0B
        equiareal = 1B
    ENDIF
    IF nside GT 0 THEN BEGIN
        pixarea = !dpi / (3L*nside^2)
        equiareal = 1B
        ns4 = 4*nside
    ENDIF
ENDIF ELSE BEGIN
    hpx = 0B
    astrom = 0B
    csystem = ''
ENDELSE

IF roll THEN BEGIN  ; check format is right:
    error = 0
    error = error && T[1] NE T[2]
    IF ~hpx THEN nside = T[1]/5
    npix = NSIDE2NPIX(nside)    ; Is this a valid NSIDE?
    error = error || npix EQ -1 ;
    ns4 = 4*nside
    IF error THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Roll requested but input map format is wrong'
        roll = 0B
        IF ~hpx THEN BEGIN
            nside = 0B  &  ns4 = 0B
        ENDIF
    ENDIF ELSE IF ~is_astrom THEN BEGIN ; Make HPX astrom structure
        proj = 'GRID'
        cdelt = [-1,1]*(90d0/nside)
        astrom = {naxis: T[1:2], cd: 0.5*[[1,-1],[1,1]], cdelt: cdelt, $
                  crpix: 0.5d0*(T[1:2] + 1), crval: [0d0, 0d0], $
                  ctype: ['XLON-HPX', 'XLAT-HPX'], $
                  longpole: 0d0, latpole: 90d0, pv2: [4d0, 3d0]}
        is_astrom = 1B
    ENDIF
ENDIF

ounit    = STRARR(N_ELEMENTS(column))
ounit[*] = unit

; Set title for intensity readout: units or "Brightness"
mid = form_unit(ounit[0])
tail = '  Zoom  '
IF verbose THEN tail = tail +'  dt (sec)'

thead = 'X pix Y pix '
IF nside GT 0 THEN thead = thead + ' NEST pix RING pix' ; HEALPix grid.
IF is_astrom OR nside GT 0 THEN thead = thead + ' longitude latitude'
title = {head: thead, unit: mid, tail: tail}

END

PRO prep_screen, state, mode, tabarr, oldgraph
;
; Prepares for creation of new tab. Cancels blinking, stashes old
; graphics state, shifts tab arrays to put screen 0 at front,
;
; Inputs:
;     state:    usual state structure
;     mode:     usual mode structure
; Outputs:
;     tabarr:   array of structures for each tab
;     oldgraph: structure containing old graphics state, if any.
;
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, state.TABS, SENSITIVE = 0

IF mode.BLINK THEN BEGIN ; Switch off blinking
    WIDGET_CONTROL, mode.BBASE, /DESTROY
    mode.BLINK = 0B
    mode.BWIN = -1  &  mode.BSWIN = -1   & mode.BBASE = -1
    WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
    gscroll_setpar, /HIDDEN
ENDIF

tabarr = *state.TABARR

; Check that we have graphics focus and if not, grab it for a bit
IF ~state.FOCUS THEN swap_lut, tabarr[0], oldgraph

; Temporarily make screen zero current (without updating visible screen)
WSET, tabarr[0].WINDOW
gscroll_newscreen, 0, tabarr, mode.ZOOM_FACTOR, mode.X_CENTRE, mode.Y_CENTRE, $
  mode.XHALF, mode.YHALF, done, 1B

*state.TABARR = tabarr

END

;--------------------------------------------------------------------------
;
; Cleanup routines:
;
PRO ximview_cleanup, id
; Ensures widget dies tidily
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, id, GET_UVALUE = state

ximview_tidy, state

PRINT, ''
PRINT, 'Ximview finished'

END
PRO ximview_tidy, state
; Does the real clean up. Called by ximview_cleanup and also by the
; catch block in the main program.
;
COMPILE_OPT IDL2, HIDDEN

IF N_ELEMENTS(state) EQ 0 THEN RETURN

IF state.LOGLUN GT 0 THEN FREE_LUN, state.LOGLUN, /FORCE
;
; Free pointers:
;
tabarr = *state.TABARR
ntab = N_ELEMENTS(tabarr)
FOR itab = 0, ntab-1 DO BEGIN
    str = tabarr[itab]
    IF PTR_VALID(str.BYTE_PTR) THEN PTR_FREE, str.BYTE_PTR
    IF PTR_VALID(str.LUT)      THEN PTR_FREE, str.LUT
    FOR i=0,2 DO IF PTR_VALID(str.RGB[i]) THEN PTR_FREE, str.RGB[i]
    IF str.TEMPORARY AND PTR_VALID(str.IM_PTR) THEN PTR_FREE, str.IM_PTR
ENDFOR

IF PTR_VALID(state.BLINK_SEQ) THEN PTR_FREE, state.BLINK_SEQ
IF PTR_VALID(state.RING0)     THEN PTR_FREE, state.RING0
IF PTR_VALID(state.ASTROM)    THEN PTR_FREE, state.ASTROM
IF PTR_VALID(state.TABARR)    THEN PTR_FREE, state.TABARR
;
; Restores graphics state if possible:
;
gscroll_set = state.MAXWIN[0] GT 0
IF gscroll_set THEN gscroll_tidy

scale = *state.OLD_GRAPH
IF N_ELEMENTS(scale) GT 0 THEN restore_lut, scale

IF PTR_VALID(state.OLD_GRAPH) THEN PTR_FREE, state.OLD_GRAPH

END
;--------------------------------------------------------------------------
;
; Event handlers:
;
PRO ximview_rdimage_event, event
; Handles events from _rdimage dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1
WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE STRMID(tag,0,12) OF
    'FILESEL_EVEN': BEGIN
        CASE event.DONE OF
            0: RETURN
            1: BEGIN
                WIDGET_CONTROL, info.INDID, GET_VALUE = index
                result = {cancel: 0B, index: index, file: event.VALUE}
            END
            2: result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    'WIDGET_TEXT_':             ; do nothing
    ELSE: MESSAGE, /INFORMATIONAL, 'Unexpected event type '+tag
ENDCASE
END

PRO ximview_rdimage, event
;  Reads an IDL-compatible image file and loads it in a new tab
;  (not yet activated)
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)

query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, /COLUMN, $
                    TITLE = 'Enter image file name')
base  = WIDGET_BASE(query, /ROW)
void  = WIDGET_LABEL(base, VALUE = 'Index No. (multi-image files):')
indID = WIDGET_TEXT(base, VALUE = '0', XSIZE = 3)
fileID = CW_FILESEL(query, /IMAGE_FILTER)

info = {retID: event.ID, indID: indID, fileID: fileID}
WIDGET_CONTROL, event.ID, GET_UVALUE = save
WIDGET_CONTROL, query, SET_UVALUE = info
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_rdimage', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result
WIDGET_CONTROL, event.ID, SET_UVALUE = save

IF result.CANCEL THEN RETURN

index = FIX(result.INDEX)
image = READ_IMAGE(result.FILE, red, green, blue, IMAGE_INDEX = index)

imsize = SIZE(image)
trucol = imsize[0] EQ 3
dims = imsize[1:imsize[0]]
IF trucol THEN dims = dims[WHERE(dims NE 3)]

; Make nominal FITS header,
mkhdr, header, imsize[imsize[0]+1], imsize[1:imsize[0]]
first = 0B

; Parameters usually set by parse_input:

; Strip directory info from filename (if any)
sep = PATH_SEP()
IF sep EQ '\' THEN sep = '\\'
firstchar = STRSPLIT(result.FILE, '.*'+sep, /REGEX)
path = STRMID(result.FILE, 0, firstchar-1)
file = STRMID(result.FILE, firstchar)
newname = file

howto = 1
tablab = 'Image'
namestr = {path: path, file: file, telescope: '', instrument: '', $
           creator: '', object: '', freq: '', stokes: ''}

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

IF ~first THEN IF ~ARRAY_EQUAL(dims,state.IMSIZE[1:2]) THEN $
  MESSAGE, 'New image does not match size of those already loaded'

                                ; Turn off blinking, make tab 0
                                ; current:
prep_screen, state, mode, tab_arr, oldgraph

                                ; set graphics focus
IF ~state.FOCUS THEN BEGIN
    state.FOCUS = 1
    WIDGET_CONTROL, event.TOP, SET_UVALUE = state
    WIDGET_CONTROL, str.DRAW, /INPUT_FOCUS
    *state.OLD_GRAPH = oldgraph
ENDIF

xsize = !D.x_vsize  &  ysize = !D.y_vsize

WIDGET_CONTROL, /HOURGLASS

ntab = N_ELEMENTS(tab_arr)
column = 1
ncol = 1
str = tab_arr[0]

lutstr  = *str.LUT
scaled = imsize[imsize[0]+1] EQ 1 ; Bytes already
lut = N_ELEMENTS(red) NE 0
IF lut THEN BEGIN
    lutstr.R = red  & lutstr.G = green  &  lutstr.B = blue
    str.LUT = PTR_NEW(lutstr)
ENDIF

; Scale image if necessary
str.ZERO = 0.0  &  str.BETA = 1.0  &  str.MODE = 0.0  &  str.SDEV = 0.0

IF scaled THEN BEGIN
    IF trucol THEN BEGIN
        rgbptr = PTRARR(3)
        inter = WHERE(imsize[1:3] EQ 3)
        CASE inter[0] OF
            0: FOR i=0,2 DO rgbptr[i] = PTR_NEW( REFORM(image[i,*,*]))
            1: FOR i=0,2 DO rgbptr[i] = PTR_NEW( REFORM(image[*,i,*]))
            2: FOR i=0,2 DO rgbptr[i] = PTR_NEW( REFORM(image[*,*,i]))
        ENDCASE
        image = PTR_NEW()
    ENDIF ELSE *str.BYTE_PTR = TEMPORARY(image)
    howto = 2
ENDIF ELSE BEGIN
    iptr = PTR_NEW(TEMPORARY(image))
    tab_arr[itab].IM_PTR = iptr
    str.BYTE_PTR = PTR_NEW(scale_image(*iptr, range,0, badcol, str.BOT, $
                                       str.TOP, ABS_RANGE = ar) )
    str.RANGE = range
    str.ABSRANGE = ar
    howto = 2
ENDELSE

IF trucol THEN BEGIN
    collab = ['Three-colour image', file, '']
    make_rgb_tab, str, ntab, xsize, ysize, rgbptr, collab, tablab, start, $
      state, mode
                                ; Label plot
    rgb_label, str
ENDIF ELSE BEGIN
    make_tabs, dummy, '', column, '', '', 1B, state, mode, str, first, start, $
      ntab, image, header, newname, howto, tablab, namestr, 0, ncol, line, $
      title, extradims, mismatch

    IF mismatch THEN RETURN

    tabarr = *state.TABARR
                                ; Draw initial screens:
    fill_screens, ncol, ntab, tabarr, mode, first, state, start

                                ; Update readout label
    title_string = title.HEAD + form_unit((*tabarr)[ntab].UNIT) + title.TAIL
WIDGET_CONTROL, state.READLAB, SET_VALUE = title_string
ENDELSE

; Enable blinking if possible
IF ncol + ntab GE 2 THEN BEGIN
    WIDGET_CONTROL, state.BLINK,  /SENSITIVE
    WIDGET_CONTROL, state.FRAMES, /SENSITIVE
    *state.BLINK_SEQ = INDGEN(ncol+ntab)
ENDIF

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode
WIDGET_CONTROL, state.TABS, /SENSITIVE

END

PRO ximview_newlog, event
; Starts a new log file
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

CASE label OF
    'Overwrite old file': logfile = 'ximview.log'
    'New sequence #': logfile = get_log_name()
    'Named...': BEGIN
        logfile = get_user_datum(event, 'Enter File name:', 30)
        IF STRCMP(logfile, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

FREE_LUN, state.LOGLUN
OPENW, loglun, logfile, /GET_LUN

; Write header lines to logfile here..
WIDGET_CONTROL, state.LABEL,   GET_VALUE = name
WIDGET_CONTROL, state.READLAB, GET_VALUE = title
nside = state.NSIDE
PRINTF,loglun,  state.VERSION, SYSTIME(), name, FORMAT = $
  "('XIMVIEW Version ',A,' Restarted log at ',A,//'Dataset: ',A)"
line = ['']
IF state.IS_ASTROM THEN line = [line, 'Coordinate system: '+state.CSYSTEM, '']
IF nside GT 0 THEN line = [line, 'Seems to be HEALPix with N_side' + $
                                 STRING(nside,FORMAT="(I5)"), '']
line = [line,title]
PRINTF, loglun, line, FORMAT="(A)"

state.LOGLUN = loglun
WIDGET_CONTROL, event.TOP, SET_UVALUE = state

END

PRO ximview_deltab, event
; Deletes a tab
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

label = WIDGET_INFO(event.ID, /NAME)
IF label EQ 'BUTTON' THEN WIDGET_CONTROL, event.ID,  GET_VALUE = label

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
tabarr = *state.TABARR

WIDGET_CONTROL, state.TABS, SENSITIVE = 0

current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
ntab = N_ELEMENTS(tabarr)
funny = current GE ntab  ; non-standard tab (blink)

IF ntab LE 1 THEN BEGIN
    ok = DIALOG_MESSAGE("Can't delete last tab", DIALOG_PARENT = event.TOP)
    RETURN
ENDIF

tabs = get_tab_uvals(state.TABS, tabid)
CASE label OF
    'Current tab': IF funny THEN BEGIN
        deadtab = WIDGET_INFO(event.TOP, FIND_BY_UNAME = tabs[current])
        id = -1
        index = current
    ENDIF ELSE BEGIN
        id = 0
        deadtab = tabarr[0].BASE
        index   = tabarr[0].SCREEN
    ENDELSE
    'Specify': BEGIN
        tabs = [tabs,'Cancel']
        index = get_user_item(event, 'Select tab to delete', tabs)
        IF index EQ ntab THEN RETURN ELSE deadtab = tabid[index]
        id = WHERE(tabarr.SCREEN EQ index)
    END
    'BASE' : BEGIN ; Directed via PAD2
        index = event.TAB
        deadtab = tabid[index]
        id = WHERE(tabarr.SCREEN EQ index)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

WIDGET_CONTROL, deadtab, /DESTROY

WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
current = WIDGET_INFO(state.TABS, /TAB_CURRENT)

; Update screen numbers for tabs which fall after the one we just deleted:
idx = WHERE(tabarr.SCREEN GT index)
IF idx[0] NE -1 THEN tabarr[idx].SCREEN -= 1

IF id GE 0 THEN BEGIN
    tabarr[id].SCREEN = -1

; delete old pixmaps, structures and pointers, and/or bulk data;
; update new current screen:
    str = tabarr[id]
    IF PTR_VALID(str.LUT)  THEN PTR_FREE, str.LUT
    IF PTR_VALID(str.BYTE_PTR) THEN BEGIN ; not RGB tab either
        IF str.TEMPORARY THEN IF PTR_VALID(str.IM_PTR) THEN $
          PTR_FREE, str.IM_PTR ELSE BEGIN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE = image, /NO_COPY
            image = 0
        ENDELSE
        PTR_FREE, str.BYTE_PTR
    ENDIF ELSE IF str.TEMPORARY THEN FOR i=0,2 DO IF PTR_VALID(str.RGB[i]) $
      THEN PTR_FREE, str.RGB[i]
ENDIF

gscroll_newscreen, current, tabarr, mode.ZOOM_FACTOR, mode.X_CENTRE, $
  mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, mode.OVERVIEW

*state.TABARR = tabarr

IF redraw_REQ THEN DEVICE, DECOMPOSED = 0B $
              ELSE DEVICE, DECOMPOSED = tabarr[0].DECOMPOSED

IF ~funny THEN ntab = ntab - 1
; disable blinking if now only one tab
IF ntab LT 2 THEN BEGIN
    WIDGET_CONTROL, state.BLINK,  SENSITIVE = 0
    WIDGET_CONTROL, state.FRAMES, SENSITIVE = 0
ENDIF

mode.DONE = done
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
                                ; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, tabarr[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS, /SENSITIVE

END

PRO ximview_2png, event
; Dumps current screen as a PNG file
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

name = GET_USER_DATUM(event, 'Enter name for PNG file:', 30, 'ximview.png')
IF name EQ 'Cancel' THEN RETURN

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
str = (*state.TABARR)[0]

IF ~state.FOCUS THEN BEGIN
    WIDGET_CONTROL, str.DRAW, /INPUT_FOCUS
    swap_lut, str, oldgraph
    *state.OLD_GRAPH = oldgraph
ENDIF

WSET, str.WINDOW
image = TVRD(TRUE = 1)
WSET, str.SCALE_INDEX
scalebar = TVRD(TRUE = 1)
WSET, str.WINDOW
pic = [[[scalebar]],[[image]]]
WRITE_PNG, name, pic

END

PRO ximview_reset, event
; Restores possibly damaged program state: draw panel to sensitive,
; pan mode to true.
;
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
mode.PAN = 1B
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
WIDGET_CONTROL, state.TABS, /SENSITIVE

END

PRO ximview_dump, event
; Prints info on current tab
;
COMPILE_OPT IDL2, HIDDEN
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

str = (*state.TABARR)[0]
PRINT, 'Tab structure:'
HELP, str, /STRUCTURE
PRINT, str
PRINT, 'LUT structure:'
HELP, *str.LUT, /STRUCTURE

END

PRO ximview_exit, event
COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, event.TOP, /DESTROY

END

PRO ximview_clear_mark, event
; Unsets the marked point and re-draws screen to remove it if in view.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
mode.XPT = -1  &  mode.YPT = -1
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
gscroll_refresh, mode.XHALF, mode.YHALF

END

PRO ximview_newlut, event
; Sets new colour table
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.ID,  GET_VALUE = label
WIDGET_CONTROL, event.TOP, GET_UVALUE = state

IF ~colmap THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
      'Colour table manipulation is not possible on this system'
    RETURN
ENDIF

str = (*state.TABARR)[0]
schemes = colour_schemes(-1)
index = (WHERE(STRCMP(label, schemes, /FOLD_CASE)))[0]
IF index NE str.COLTAB THEN BEGIN
    WIDGET_CONTROL, state.TABS, SENSITIVE = 0
    done = 1

    IF ~state.FOCUS THEN swap_lut, dummy, old_graph

    ximview_lut, index, str.IZERO, decomp, bot, top
    str.COLTAB     = index
    str.DECOMPOSED = decomp
    str.BOT        = bot
    str.TOP        = top
    TVLCT, r, g, b, /GET
    *str.LUT = {r:r, g:g, b:b, line: !P.color, absent: !P.background}
    (*state.TABARR)[0] = str
    gscroll_setpar, /BLANK

    WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

    IF state.GLOBAL_COLOUR THEN BEGIN
        ntab = N_ELEMENTS(*state.TABARR)
        (*state.TABARR).COLTAB     = index
        (*state.TABARR).DECOMPOSED = decomp
        (*state.TABARR).BOT        = bot
        (*state.TABARR).TOP        = top
        FOR i=0,ntab-1 DO *(*state.TABARR)[i].LUT = *str.LUT

        IF redraw_REQ THEN FOR i = 0,ntab-1 DO $
          update_screen, state.TABARR, i, mode, done

    ENDIF ELSE IF redraw_req THEN $
      update_screen, state.TABARR, 0, mode, done

    IF done EQ 0 THEN WIDGET_CONTROL, str.DRAW, TIMER=0.5
    mode.DONE = done

    IF ~state.FOCUS THEN restore_lut, old_graph

    WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
    WIDGET_CONTROL, state.TABS, /SENSITIVE
ENDIF

END

PRO ximview_colour, event
; Switches between individual and global colour tables
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP, GET_UVALUE = state

oldglob = state.GLOBAL_COLOUR
CASE label OF
    'Same for all':         newglob = 1B
    'Separate on each tab': newglob = 0B
ENDCASE
state.GLOBAL_COLOUR = newglob
WIDGET_CONTROL, event.TOP, SET_UVALUE = state

IF newglob && newglob NE oldglob THEN BEGIN
    WIDGET_CONTROL, /HOURGLASS
    ntab = N_ELEMENTS(*state.TABARR)
    lutptr = (*state.TABARR)[0].LUT
    lutstr = *lutptr
    WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
    FOR i = 1,ntab-1 DO IF PTR_VALID((*state.TABARR)[i].BYTE_PTR) THEN BEGIN
        lutptr = (*state.TABARR)[i].LUT
        *lutptr = lutstr
        IF redraw_req THEN update_screen, state.TABARR, i, mode, done
    ENDIF

    IF redraw_req THEN BEGIN
        mode.DONE = done
        WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
                                ; Request another go if loading not finished
        IF done EQ 0 THEN WIDGET_CONTROL, (*state.TABARR)[0].DRAW, TIMER = 0.5
    ENDIF
ENDIF

END

PRO ximview_goto, event
; Sets view centre to a specified position
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
WIDGET_CONTROL, state.TABS, SENSITIVE = 0

pixvals = STRCOMPRESS(STRING(mode.X_CENTRE, mode.Y_CENTRE))
bnames = 'image pixel'
coindex = [0]
nside = state.NSIDE
IF nside NE 0 THEN BEGIN
    bnames = bnames+'|HP ring pixel|HP nest pixel'
    coindex = [coindex, 1, 2]
ENDIF
IF state.IS_ASTROM $ ; || state.NSIDE NE 0
  THEN  BEGIN
    bnames = bnames+'|(long\, lat)'
    coindex = [coindex, 3]
ENDIF
desc = ['0, TEXT, ' + pixvals + $
        ', LABEL_LEFT=Enter position:, WIDTH=20, TAG=text', $
        '0, BUTTON,' +bnames+ ', ROW, EXCLUSIVE, SET_VALUE=0, TAG=coord', $
        '1, BASE,, ROW', '0, BUTTON, Done, QUIT',$
        '2, BUTTON, Cancel, QUIT, TAG=cancel']
form = CW_FORM(desc, /COLUMN, TITLE = 'Set centre of view')

IF form.CANCEL THEN RETURN

xhalf = mode.XHALF  &  yhalf = mode.YHALF
TVCRS, xhalf, yhalf

pix = FLOAT(STRSPLIT(form.TEXT,', ',/EXTRACT, COUNT = count))
r2d = 180d0 / !dpi

CASE coindex[form.COORD] OF
    0: BEGIN
        IF count NE 2 THEN BEGIN
            MESSAGE, /INFORMATIONAL, $
              'Enter two numbers for pixel coordinates'
            RETURN
        ENDIF
        xpix = pix[0]  &  ypix = pix[1]
    END
    1: BEGIN
        IF count NE 1 THEN BEGIN
            MESSAGE, /INFORMATIONAL, 'Enter only one number for HP pixel'
            RETURN
        ENDIF
        pix = ROUND(pix[0])
        order = 'RING'
        PIX2ANG_RING, nside, pix, theta, phi
        ll = phi[0]*r2d  &  bb = 90d - theta[0]*r2d
    END
    2: BEGIN
        IF count NE 1 THEN BEGIN
            MESSAGE, /INFORMATIONAL, 'Enter only one number for HP pixel'
            RETURN
        ENDIF
        pix = ROUND(pix[0])
        order = 'NESTED'
        PIX2ANG_NEST, nside, pix, theta, phi
        ll = phi[0]*r2d  &  bb = 90d - theta[0]*r2d
    END
    3: BEGIN
        IF count NE 2 THEN BEGIN
            MESSAGE, /INFORMATIONAL, $
              'Enter two numbers for coordinates'
            RETURN
        ENDIF
        ll = pix[0]  &  bb = pix[1]
    END
ENDCASE

IF form.COORD GT 0 THEN BEGIN
    IF state.IS_ASTROM THEN AD2XY, ll, bb, *state.ASTROM, xpix, ypix $
    ELSE BEGIN ; must be XPH grid
; Not yet available, this branch should not be reachable at present.
;        coord = HP2XPH(nside, pix, order, state.PROJ)
;        xpix = coord[0]  &  ypix = coord[1]
    ENDELSE
ENDIF

xpix = (ROUND(xpix) > 0) < (state.IMSIZE[1]-1)
ypix = (ROUND(ypix) > 0) < (state.IMSIZE[2]-1)

; Now we've found the point, update the screen

IF mode.OVERVIEW THEN BEGIN
    mode.OVERVIEW = 0B
    WIDGET_CONTROL, state.ZOOMCOL, /SENSITIVE
    WIDGET_CONTROL, state.READOUT, /SENSITIVE
ENDIF

ierr = 0
temp = (*state.TABARR).BYTE_PTR

coord = gscroll(temp, xpix, ypix, xhalf, yhalf, $
                mode.ZOOM_FACTOR, ierr, 1, done, DO_WRAP = mode.ROLL)
IF ierr NE 0 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)

mode.OXTV     = xhalf  &  mode.OYTV     = yhalf
mode.X_CENTRE = xpix   &  mode.Y_CENTRE = ypix
mode.XPIX     = xpix   &  mode.YPIX     = ypix
mode.DONE     = done   &  mode.NEW_VIEW = 0
mode.DRAG     = 0B

IF mode.XPT GE 0 THEN $
  marker, mode.XPT, mode.YPT, xpix, ypix, xhalf, yhalf, mode.ZOOM, mode.ZFAC

; If panels remain to be loaded, send timer event to DRAW window
; requesting re-draw:
IF done EQ 0 THEN WIDGET_CONTROL, (*state.TABARR)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS, SET_UVALUE=mode
pix_print, state, 0, start

WIDGET_CONTROL, state.TABS, /SENSITIVE
END

PRO ximview_autoscale_all, event
; Autoscales all tabs
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

ntab = N_ELEMENTS(*state.TABARR)
IF ntab EQ -1 THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'No image to scale'
    RETURN
ENDIF

col = INDGEN(ntab)
tablab = get_tab_uvals(state.TABS)

WIDGET_CONTROL, state.TABS, SENSITIVE = 0

do_plot = 0

FOR itab=0,ntab-1 DO BEGIN
    str = (*state.TABARR)[itab]
    IF ~PTR_VALID(str.BYTE_PTR) THEN CONTINUE ; funny tab

    name = tablab[str.SCREEN]
    MESSAGE, /INFORMATIONAL, STRING(name, FORMAT = "('Tab: ',A)")

    IF str.MODE EQ 0.0 THEN BEGIN
        IF ~str.IM_PTR THEN BEGIN
            howto = 1
            imap = str.SCREEN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
        ENDIF ELSE BEGIN
            howto = 3
            data = str.IM_PTR
            imap = 0
        ENDELSE
                                ;Get mode and rms estimate if not done already
        scaling_params, dummy, data, imap, 0, howto, col, $
          ar, zero, sdev, PLOT = do_plot
        IF howto EQ 1 THEN $
          WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY

        MESSAGE, /INFORMATIONAL, 'Found mode: ' + numunit(zero,str.UNIT) + $
          ' and estimated rms: ' + numunit(sdev,str.UNIT)
    ENDIF ELSE BEGIN
        ar   = str.ABSRANGE
        zero = str.MODE
        sdev = str.SDEV
    ENDELSE
                                ; Set min and max appropriate for pseudo-col:
    set_scale, [ar, zero, sdev], str
    range = str.RANGE
    (*state.TABARR)[itab] = str

                                ; Rescale byte images:
    IF ~str.IM_PTR THEN BEGIN
        WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
        *str.BYTE_PTR = scale_image(data, str.RANGE, str.WRAP, badcol, $
                                    str.BOT, str.TOP, str.TRFUNC, str.ZERO, $
                                    str.BETA, ABS_RANGE = str.ABSRANGE)
        WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
    ENDIF ELSE *str.BYTE_PTR = scale_image(*str.IM_PTR, str.RANGE, str.WRAP, $
                                   badcol, str.BOT, str.TOP, str.TRFUNC, $
                                   str.ZERO, str.BETA, ABS_RANGE= str.ABSRANGE)
ENDFOR

nside = state.NSIDE
IF nside NE 0 && state.IS_ASTROM THEN $
  fill_gores, nside, state.IMSIZE, state.ASTROM, (*state.TABARR).BYTE_PTR

FOR itab=0,ntab-1 DO BEGIN
    IF redraw_req THEN BEGIN
        lutptr = (*state.TABARR)[itab].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    update_screen, state.TABARR, itab, mode, done
ENDFOR

mode.DONE = done
                                ; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, (*state.TABARR)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS, /SENSITIVE
WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

END

PRO ximview_grid, event
; Displays coordinate grids
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0

MESSAGE, /INFORMATIONAL, 'Coordinate grids not yet enabled'
END

PRO ximview_vo, event
; Plots data from the Virtual Observatory
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 0

MESSAGE, /INFORMATIONAL, 'Catalog access not yet enabled'
END

PRO ximview_imstat_opts, event
; Sets options for imstats
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
CASE label OF
    'Region of Interest': state.ROI = 1
    'Box':                state.ROI = 0
    'Set box size': BEGIN
; adjust imstats window
        default = STRJOIN(STRTRIM(STRING(state.STATBOX),2),', ')
        text = get_user_datum(event, $
                              'Box size in pixels (x,y; or x for square):', $
                              12, default)
        IF STRCMP(text, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
        strings = STRSPLIT(text,' ,',/EXTRACT)
        IF N_ELEMENTS(strings) EQ 1 THEN strings = REPLICATE(strings,2)
        state.STATBOX = FIX(strings)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

WIDGET_CONTROL, event.TOP, SET_UVALUE=state
END

PRO ximview_maxfit_opts, event
; Sets options for maxfit
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
CASE label OF
    'Find extremum': state.PEAK = 1
    'Find maximum':  state.PEAK = 0
    'Find minimum':  state.PEAK = 2
    'Set box size': BEGIN
; adjust imstats window
        text = get_user_datum(event, $
                              'Peakfit search box size (pixels)', $
                              3, STRING(state.MAXBOX,FORMAT = "(I3)"))
        IF STRCMP(text, 'CANCEL', 6, /FOLD_CASE) THEN RETURN
        state.MAXBOX = FIX(text)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

WIDGET_CONTROL, event.TOP, SET_UVALUE=state
END

PRO ximview_setprop_event, event
;  Processes events from set property dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE tag OF
    'WIDGET_DROPLIST': BEGIN
        CASE event.ID OF
            info.TABLIST: WIDGET_CONTROL, info.CURRENT, $
              SET_VALUE = info.PROPLIST[event.INDEX]
            info.CHOOSE: WIDGET_CONTROL, event.ID, SET_UVALUE = event.INDEX
        ENDCASE
    END
    'WIDGET_COMBOBOX': ;   WIDGET_CONTROL, event.ID, SET_UVALUE = event.STR
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Accept': BEGIN     ; Extract results from widget
                tab = WIDGET_INFO(info.TABLIST, /DROPLIST_SELECT)
;                WIDGET_CONTROL, info.CHOOSE, GET_UVALUE = property
;                IF property EQ -1 THEN $
                type = WIDGET_INFO(info.CHOOSE, /NAME)
                IF type EQ 'DROPLIST' THEN $
                  property = WIDGET_INFO(info.CHOOSE, /DROPLIST_SELECT) ELSE $
                  property = WIDGET_INFO(info.CHOOSE, /COMBOBOX_GETTEXT)
                IF info.CHOOSE2 NE -1 THEN BEGIN
                    WIDGET_CONTROL, info.CHOOSE2, GET_VALUE = prop2
                    IF ~STREGEX(prop2, '^ *$', /BOOLEAN) THEN property = prop2
                ENDIF
                result = {cancel: 0B, tab: tab, property: property}
            END
            'Cancel': result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Unrecognised event type: ' + tag
ENDCASE

END

PRO ximview_setprop, event, tablab, scr0, property, proplist, options, $
             cancel, tab, new, EDITABLE = editable
; Launches a dialog to choose a property for a specified tab from a
; list (optionally editable, ie you can write your own value if not in
; the list)
;
; Inputs:
;     event:    Used to hold results
;     tablab:   Labels for the tabs available to be set
;     scr0:     Index of tablab corresponding to the current tab
;     property: Name of property to set
;     proplist: Currently-assigned properties for each tab
;     options:  List of options for the property
;     editable: True if the user can enter a new value as well as
;               choosing from a list
; Outputs:
;     cancel:   True if operation abandoned
;     tab:      tab to set (index of tablab)
;     new:      New value of property.
;
COMPILE_OPT IDL2, HIDDEN

editable = KEYWORD_SET(editable)

is_unix = STRCMP(!version.OS_FAMILY, 'UNIX', 4, /FOLD_CASE)

query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, /COLUMN, $
                    TITLE = 'Set tab ' + property)
base    = WIDGET_BASE(query, /ROW)
void    = WIDGET_LABEL(base, VALUE ='Tab to set')
tablist = WIDGET_DROPLIST(base, VALUE = tablab, UVALUE = scr0)
void    = WIDGET_LABEL(base, VALUE = ' Current ' + property + ':')
current = WIDGET_LABEL(base, VALUE = proplist[scr0], /DYNAMIC_RESIZE)
WIDGET_CONTROL, tablist, SET_DROPLIST_SELECT = scr0

now = WHERE(proplist[scr0] EQ options)
propset = now[0] NE -1
base = WIDGET_BASE(query, /ROW)
void = WIDGET_LABEL(base, VALUE = ' New '+ property + ':')
IF editable THEN BEGIN
    IF is_unix && !VERSION.release EQ 6.0 THEN BEGIN
        choose = WIDGET_DROPLIST(base, VALUE = options, UVALUE = now[0])
        IF propset THEN WIDGET_CONTROL, choose, SET_DROPLIST_SELECT = now[0]
        void = WIDGET_LABEL(base, VALUE = 'or:')
        choose2 = WIDGET_TEXT(base, VALUE = '', xsize = 10, /EDITABLE)
    ENDIF ELSE BEGIN
        choose = WIDGET_COMBOBOX(base, VALUE = options, $
                                 UVALUE = proplist[scr0], /EDITABLE)
        IF propset THEN WIDGET_CONTROL, choose, SET_COMBOBOX_SELECT = now[0]
        choose2 = -1
    ENDELSE
ENDIF ELSE BEGIN
    choose = WIDGET_DROPLIST(base, VALUE = options, UVALUE = now[0])
    IF propset THEN WIDGET_CONTROL, choose, SET_DROPLIST_SELECT = now[0]
    choose2 = -1
ENDELSE

bbase   = WIDGET_BASE(query, /ROW)
void   = WIDGET_BUTTON(bbase, VALUE = 'Accept')
void   = WIDGET_BUTTON(bbase, VALUE = 'Cancel')

info = {retID: event.ID, tablist: tablist, proplist: proplist, $
        current: current, choose: choose, choose2: choose2}

WIDGET_CONTROL, event.ID, GET_UVALUE = save
WIDGET_CONTROL, query, SET_UVALUE = info
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_setprop', query
WIDGET_CONTROL, event.ID, GET_UVALUE = result
WIDGET_CONTROL, event.ID, SET_UVALUE = save

cancel = result.CANCEL
IF cancel THEN RETURN

tab = result.TAB  &  new = result.PROPERTY

END

PRO ximview_setpol, event
; Sets the polarization code in a tab to user's preferred value.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
scodes = ['YX', 'XY', 'YY', 'XX', 'LR', 'RL', 'LL', 'RR', 'unknown', $
          'I', 'Q', 'U', 'V', 'Pol Intensity', 'Frac. Pol', 'Pol Angle']
polcodes = (*state.TABARR).POLCODE
polstr   = scodes[polcodes+8]
screens  = (*state.TABARR).SCREEN
good     = WHERE( PTR_VALID( (*state.TABARR).BYTE_PTR ))
scr0     = screens[good[0]]
iscreen  = BSORT(screens)
good     = WHERE( PTR_VALID( (*state.TABARR)[iscreen].BYTE_PTR ))
iscreen  = iscreen[good]
tabcode  = polstr[iscreen]
tablab   = get_tab_uvals(state.TABS)
tablab   = tablab[screens[iscreen]]
scr0     = (WHERE(screens[iscreen] EQ scr0))[0]

ximview_setprop, event, tablab, scr0, 'polarization', tabcode, scodes, $
  cancel, tab, polcode

IF cancel THEN RETURN

polcode -= 8
(*state.TABARR)[iscreen[tab]].POLCODE = polcode
MESSAGE, /INFORMATIONAL, 'Setting tab '+tablab[tab]+' polarization to ' $
  + scodes[polcode+8]

END

PRO ximview_setfreq, event
; Sets the frequency code in a tab to user's preferred value.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
freqcode = (*state.TABARR).FREQCODE
screens  = (*state.TABARR).SCREEN
good     = WHERE( PTR_VALID( (*state.TABARR).BYTE_PTR ))
scr0     = screens[good[0]]
iscreen  = BSORT(screens)
good     = WHERE( PTR_VALID( (*state.TABARR)[iscreen].BYTE_PTR ))
iscreen  = iscreen[good]
freqcode  = freqcode[iscreen]
tablab   = get_tab_uvals(state.TABS)
tablab   = tablab[screens[iscreen]]
scr0     = (WHERE(screens[iscreen] EQ scr0))[0]

nulls = WHERE(freqcode EQ '')
IF nulls[0] NE -1 THEN freqcode[nulls] = '          '
ximview_setprop, event, tablab, scr0, 'frequency', freqcode, freqcode, $
  cancel, tab, freq, /EDITABLE

IF cancel THEN RETURN

IF SIZE(freq, /TYPE) NE 7 THEN freq = freqcode[freq]
freq = STRTRIM(freq, 2)
(*state.TABARR)[iscreen[tab]].FREQCODE = freq
IF FREQ EQ '' THEN freq = '<null>'
MESSAGE, /INFORMATIONAL, 'Setting tab ' + tablab[tab] + ' frequency to ' $
  + freq

END

PRO ximview_help, event
; Launches help windows
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

; Find help directory
info = ROUTINE_INFO('ximview', /SOURCE)
dir = STRSPLIT(info.PATH, 'ximview.pro', /REGEX, /EXTRACT) + 'docs/'
unhelp = ['Abandon hope all ye who enter here.', $
          "Lasciate ogne speranza, voi ch'intrate", $
          'Send requests for help to /dev/null', $
          'ZOMG!!! wtf?', 'Please write the help page you need.']
nun = N_ELEMENTS(unhelp)
WIDGET_CONTROL, event.ID, GET_VALUE = label
WIDGET_CONTROL, event.TOP,  GET_UVALUE = state

assistant = 0B
CASE label OF
    'Help': BEGIN
        file = 'help.txt'  &  title = "Don't Panic!"
        done = 'Done with Ximview HELP'
        height = 50
        assistant = !version.RELEASE GE 6.2
;        ok = DIALOG_MESSAGE(unhelp[ FIX(RANDOMU(seed,1)*nun) ], $
;                            DIALOG_PARENT = event.TOP, TITLE = title)
    END
    'Release Notes': BEGIN
        file = 'release_notes.txt' &  title = 'Release Notes'
        done = 'Done with Release notes'  &  height = 24
;        ok = DIALOG_MESSAGE('E E F G G F E D C C D E Ee di D', $
;                            DIALOG_PARENT = event.TOP, TITLE = title)
    END
    'About': BEGIN
        file = 'about.txt' &  title = 'About Ximview'
        done = 'Done with About Ximview'  &  height = 24
;        ok = DIALOG_MESSAGE('Ximview is about 10 cm across', $
;                            DIALOG_PARENT = event.TOP, $
;                            TITLE = title)
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + label + ' not yet available.'
ENDCASE

IF assistant THEN BEGIN
    ONLINE_HELP, BOOK='ximview.adp'
ENDIF ELSE BEGIN
    XDISPLAYFILE, dir+file, GROUP = event.TOP, DONE_BUTTON = done, $
      HEIGHT = height, TITLE = title
ENDELSE
END

PRO ximview_event, event
;
; Processes events from the top level base: resize, timer, keyboard focus.
; Also swallows events froll the draw window which are actually
; processed by ximview_scroll and possibly CW_DEFROI.
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = state

name = TAG_NAMES(event, /STRUCTURE_NAME)
CASE name OF
    'WIDGET_BUTTON': BEGIN
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        CASE label OF
            'Restore default screen size': BEGIN
                state.NEWWIDTH  = -1
                state.NEWHEIGHT = -1
                                ; Disable events from the draw window
                                ; during resize:
                WIDGET_CONTROL, state.TABS, SENSITIVE = 0
                WIDGET_CONTROL, event.TOP, SET_UVALUE = state
                ximview_resize, event.TOP, state
            END
            ELSE:
        ENDCASE
    END
    '': BEGIN                   ; Probably menu event
        WIDGET_CONTROL, event.ID, GET_VALUE = label
        option = event.VALUE
        MESSAGE,/INFORMATIONAL, 'Option ' + option + ' from menu ' $
          + label +  ' is not yet implemented.'
    END
    'WIDGET_BASE': BEGIN        ; Resize event.

; Record current size of tlb. Don't update immediately as many resize
; events are produced when moving the edge of a window.
        state.NEWWIDTH  = event.X
        state.NEWHEIGHT = event.Y
                                ; Disable events from the draw window
                                ; during resize:
        WIDGET_CONTROL, state.TABS, SENSITIVE = 0
        WIDGET_CONTROL, event.TOP, SET_UVALUE = state
                                ; Send off timer event. Last one to
                                ; arrive is acted on (!).
        WIDGET_CONTROL, event.TOP, TIMER = 0.5
    END
    'WIDGET_TIMER': BEGIN       ;... and act here.
        ximview_resize, event.TOP, state
    END
    'WIDGET_KBRD_FOCUS' : IF event.ENTER THEN BEGIN
                                ; When window becomes active, set
                                ; graphics state to what we want and
                                ; save previous values.

; Return immediately if somehow we already have focus:
        IF state.FOCUS THEN RETURN
        state.FOCUS = 1
        WIDGET_CONTROL, event.TOP, SET_UVALUE = state

        WIDGET_CONTROL, state.TABS, SENSITIVE = 0 ; avoid confusion

        str = (*state.tabarr)[0]

        swap_lut, str, oldgraph
        *state.OLD_GRAPH = oldgraph
        WSET, str.WINDOW

        WIDGET_CONTROL, state.TABS, SENSITIVE = 1
    ENDIF ELSE BEGIN            ; Lost keyboard focus
        IF ~state.FOCUS THEN RETURN
        state.FOCUS = 0
        WIDGET_CONTROL, event.TOP, SET_UVALUE = state
                                ; Restore old colour table etc
        restore_lut, *state.OLD_GRAPH
    ENDELSE
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: '+name
ENDCASE

END

PRO ximview_overview, event
;
; Makes overview plot on request from "overview" button.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

WIDGET_CONTROL, event.TOP,   GET_UVALUE = state

ximview_resize, event.top, state ; Checks the window has not been re-sized

tabarr = *state.tabarr
ntab = N_ELEMENTS(tabarr)
temp = tabarr.BYTE_PTR
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

; Switch to overview mode from zoom mode:
mode.OVERVIEW = 1
WIDGET_CONTROL, state.ZOOMCOL, SENSITIVE = 0
WIDGET_CONTROL, state.READOUT, SENSITIVE = 0

; Set cursor on marked point if available, otherwise centre of FOV.
IF mode.XPT GE 0 THEN BEGIN
    mode.XPIX = mode.XPT  &  mode.YPIX = mode.YPT
ENDIF ELSE BEGIN
    mode.XPIX = mode.X_CENTRE  &  mode.YPIX = mode.Y_CENTRE
ENDELSE

FOR itab = 0,ntab-1 DO BEGIN
    WSET, tabarr[itab].WINDOW

    IF redraw_req THEN BEGIN
        lutptr = tabarr[itab].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    null3 = REPLICATE(PTR_NEW(),3)
    IF ARRAY_EQUAL(tabarr[itab].RGB, null3) THEN tptr = tabarr[itab].BYTE_PTR $
    ELSE tptr = tabarr[itab].RGB
    overview, tptr, mode.ZOOM_FACTOR, mode.XPIX, mode.YPIX, $
      mode.XHALF, mode.YHALF, mode.X_CENTRE, mode.Y_CENTRE, resamp, corner
ENDFOR
mode.RESAMP = resamp  &  mode.CORNER = corner
mode.NEW_VIEW = 1

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

END

FUNCTION ximview_zoom, event
;
; Changes zoom factor & redraws view as appropriate
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP,   GET_UVALUE = state

ximview_resize, event.TOP, state ; Just checking (usually)

WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode

tabarr = state.TABARR
temp = (*tabarr).BYTE_PTR
zoom = mode.ZOOM

WIDGET_CONTROL, event.ID, GET_UVALUE = uval
CASE uval OF
    'IN' :  zoom += 1          ; zoom in
    '1:1' : zoom  = 0          ; zoom factor 1
    'OUT' : zoom -= 1          ; zoom out
ENDCASE
zoom_factor = 2.^zoom

xhalf = mode.XHALF     &  yhalf = mode.YHALF
xcent = mode.X_CENTRE  &  ycent = mode.Y_CENTRE
ierr = 0
coord = gscroll(temp, xcent, ycent, xhalf, yhalf, $
                zoom_factor, ierr, 1, done, DO_WRAP = state.ROLL)

lowerr = ierr MOD 8
IF lowerr EQ 6 || lowerr EQ 7 THEN BEGIN
; Can't zoom any more.
    WIDGET_CONTROL, event.ID, SENSITIVE = 0
    state.MAXZOOM = event.ID
    WIDGET_CONTROL, event.TOP, SET_UVALUE = state
    IF lowerr EQ 7 THEN BEGIN
; Max zoom exceeded: shouldn't get here because we hit ierr 6 first.
        MESSAGE, /INFORMATIONAL, 'Attempt to exceed max zoom'
        RETURN, ierr
    ENDIF
ENDIF ELSE IF state.MAXZOOM NE 0 THEN BEGIN ; restore zoom option
    WIDGET_CONTROL, state.MAXZOOM, /SENSITIVE
    state.MAXZOOM = 0
    WIDGET_CONTROL, event.TOP, SET_UVALUE = state
ENDIF

IF lowerr NE 0 && lowerr LT 6 THEN $
  MESSAGE, 'GSCROLL error '+STRING(ierr)

; rolling disabled if ierr = 8
mode.ROLL = ierr-lowerr NE 8 ? state.ROLL : 0

mode.ZOOM = zoom
zfac = 2^ABS(zoom)
mode.ZFAC = zfac
mode.ZOOM_FACTOR = zoom_factor
mode.XHALF = xhalf  &  mode.YHALF = yhalf

;  Mark current point:
IF mode.XPT GE 0 THEN marker, mode.XPT, mode.YPT, xcent, ycent, $
                              xhalf, yhalf, zoom, zfac

                                ; Disable panning if image fits in screen
bltv = tv2im(0, 0, xcent, ycent, xhalf, yhalf, zoom, zfac)
trtv = tv2im(!D.x_vsize, !D.y_vsize, xcent, ycent, xhalf, yhalf, zoom, zfac)

mode.PAN =  MAX(bltv) GT 0 OR MAX(state.IMSIZE[1:2] - trtv) GT 1

mode.DONE = done
; Request another go if loading not finished
IF done EQ 0 THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5

WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

                                ; Record updated zoom state
pix_print, state, 0, start

RETURN, ierr
END

PRO ximview_tab, event
; Processes tab change events and also swallows events from draw window
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)

WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS,  GET_UVALUE = mode
xpix = mode.XPIX  &  ypix = mode.YPIX
zoomf = mode.ZOOM_FACTOR

name = TAG_NAMES(event, /STRUCTURE_NAME)
SWITCH name OF
    'WIDGET_TIMER':          ; Do nothing, processed by ximview_scroll
    'WIDGET_DRAW': BREAK     ; Do nothing, processed by ximview_scroll
    'WIDGET_TAB':            ; Actual tab event or..
    '' :  BEGIN              ; Anonymous event from ximview_blink
        tabarr = state.TABARR
                                ; Check we don't have a temporary tab:
        IF event.TAB GT MAX((*tabarr).SCREEN) THEN RETURN

                                ;   Check for bizarre counting bug:
        current = event.TAB EQ 1 ? WIDGET_INFO(state.TABS, /TAB_CURRENT) $
                                 : event.TAB

        IF ~state.FOCUS THEN swap_lut, (*tabarr)[0], old_graph

        WSET, (*tabarr)[0].WINDOW

        gscroll_newscreen, current, (*tabarr), zoomf, $
          mode.X_CENTRE, mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, $
          mode.OVERVIEW

        draw = (*tabarr)[0].DRAW
        mode.DONE = done
        IF done EQ 0B THEN WIDGET_CONTROL, draw, TIMER = 0.5

                                ; Disable analysis tools for RGB
                                ; window
        imstats = WIDGET_INFO(state.BLINK, /SIBLING)
        peakfit = WIDGET_INFO(imstats, /SIBLING)
        sens = PTR_VALID((*tabarr)[0].BYTE_PTR)
        WIDGET_CONTROL, imstats, SENSITIVE = sens
        WIDGET_CONTROL, peakfit, SENSITIVE = sens

        IF redraw_req THEN DEVICE, DECOMPOSED = 0 $
                      ELSE DEVICE, DECOMPOSED = (*tabarr)[0].DECOMPOSED

                                ;   Update units on readout label
        mid = form_unit( (*tabarr)[0].UNIT)

        title_string = state.TITLE.HEAD + mid + state.TITLE.TAIL
        WIDGET_CONTROL, state.READLAB, SET_VALUE=title_string

                                ;  Mark current point
        IF ~mode.OVERVIEW THEN BEGIN
            IF mode.XPT GE 0 THEN marker, mode.XPT, mode.YPT, mode.X_CENTRE, $
              mode.Y_CENTRE, mode.XHALF, mode.YHALF, mode.ZOOM, mode.ZFAC

            pix_print, state, 0, start
        ENDIF

        IF ~state.FOCUS THEN restore_lut, old_graph

        BREAK
    END
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: '+name
ENDSWITCH
WIDGET_CONTROL, state.TABS,  SET_UVALUE = mode

END

FUNCTION ximview_scroll, event
; Processes events from main draw widget
; Jobs:
;  (0) Overview or zoom mode?
;      If overview, set xpix, ypix on button press and enter zoom
;      mode.
;  Otherwise:
;  (1) Is button 1 down or not?
;      Yes: drag image using gscroll subroutines
;      No:  Record current cursor position
;  (2) Has button 2 been clicked (act on button down)
;      If so, mark spot and record in log
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)

returnable = event

WIDGET_CONTROL, event.TOP,  GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
tabarr = state.TABARR
str = (*tabarr)[0]

xhalf = mode.XHALF  &  yhalf = mode.YHALF

name = TAG_NAMES(event, /STRUCTURE_NAME)
CASE name OF
    'WIDGET_TIMER': BEGIN     ; gscroll didn't finish loading panel
        IF mode.DONE THEN RETURN, 0 ; seems to have been loaded in mean time
;
; Check that we have graphics focus and if not, grab it for a bit
;
        IF ~state.FOCUS THEN BEGIN
            swap_lut, str, oldgraph
            new_view = 1 ; We might as well take our time...
        ENDIF ELSE new_view = 0

        temp = (*tabarr).BYTE_PTR
        coord = gscroll(temp, mode.X_CENTRE, mode.Y_CENTRE, xhalf, yhalf, $
                        mode.ZOOM_FACTOR, ierr, new_view, done, $
                        DO_WRAP = mode.ROLL)

        IF ierr NE 0 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)
        mode.DONE = done
        returnable = 0 ; don't pass this up the chain.
    END
    'WIDGET_DRAW': BEGIN
;
; Grab graphics focus if necessary
        IF ~state.FOCUS THEN BEGIN
            WIDGET_CONTROL, state.TABS, SENSITIVE = 0 ; avoid confusion

            state.FOCUS = 1
            WIDGET_CONTROL, event.TOP, SET_UVALUE = state

; Explicitly get keyboard focus so we are sure to let it go again
; later. The widget concerned must be the draw window rather than the
; tlb because you can't set INPUT_FOCUS on a base, only on buttons,
; draw and text widgets. However, this should transfer to the tlb
; automatically.
            WIDGET_CONTROL, event.ID, /INPUT_FOCUS
            swap_lut, str, oldgraph
            *state.OLD_GRAPH = oldgraph
            WIDGET_CONTROL, state.TABS, SENSITIVE = 1
        ENDIF
        WSET, str.WINDOW
        IF ~redraw_req THEN DEVICE, DECOMPOSED = str.DECOMPOSED

        xtv = event.X  & ytv = event.Y

        xpmax = state.IMSIZE[1]-1
        ypmax = state.IMSIZE[2]-1
        zfac  = mode.ZFAC  &  zoom = mode.ZOOM


        IF ~mode.OVERVIEW THEN BEGIN ; We are in pan/zoom mode

            x_centre = mode.X_CENTRE  &  y_centre = mode.Y_CENTRE

            drag = mode.PAN AND (event.PRESS OR (mode.DRAG && ~event.RELEASE))
            mode.DRAG = drag

            IF drag THEN BEGIN
                               ; Set cursor to gripping hand
                cursor_grip
;               DEVICE, CURSOR_STANDARD=52

                temp = (*tabarr).BYTE_PTR

; Calculate new centre pixel = shifted by opposite amount from cursor
;
                oxtv = mode.OXTV  &  oytv = mode.OYTV
                IF zoom LE 0 THEN BEGIN
                    xshift = xtv - oxtv
                    yshift = ytv - oytv
                    xpix = ((x_centre - xshift*zfac) > 0) < xpmax
                    ypix = ((y_centre - yshift*zfac) > 0) < ypmax
                ENDIF ELSE BEGIN
                    xpix = ((x_centre - (xtv/zfac) + (oxtv/zfac)) > 0) < xpmax
                    ypix = ((y_centre - (ytv/zfac) + (oytv/zfac)) > 0) < ypmax
                ENDELSE
                IF mode.ROLL THEN BEGIN
                    ns4 = state.NS4
                    IF xpix + ypix GT (state.IMSIZE[1] + ns4) THEN BEGIN
                        xpix -= ns4 & ypix -= ns4
                    ENDIF ELSE IF xpix + ypix LT state.NSIDE THEN BEGIN
                        xpix += ns4 & ypix += ns4
                    ENDIF
                ENDIF
                ierr = 0
                coord = gscroll(temp, xpix, ypix, xhalf, yhalf, $
                        mode.ZOOM_FACTOR, ierr, 0, done, DO_WRAP=mode.ROLL)
                IF ierr NE 0 THEN MESSAGE, 'GSCROLL error '+STRING(ierr)
                mode.DONE = done
                mode.X_CENTRE = xpix
                mode.Y_CENTRE = ypix

;  Get position of pixel under cursor
                coord = tv2im(oxtv, oytv, x_centre, y_centre, xhalf, yhalf, $
                              zoom, zfac)
                xpix = coord[0]  &  ypix = coord[1]
                IF mode.ROLL THEN BEGIN
                    IF xpix+ypix GT (state.IMSIZE[1] + ns4) THEN BEGIN
                        xpix -= ns4 & ypix -= ns4
                    ENDIF ELSE IF xpix + ypix LT state.NSIDE THEN BEGIN
                        xpix += ns4 & ypix += ns4
                    ENDIF
                ENDIF

                xpix = (xpix > 0) < xpmax
                ypix = (ypix > 0) < ypmax
            ENDIF ELSE BEGIN
                DEVICE, /CURSOR_CROSSHAIR

                coord = tv2im(xtv, ytv, x_centre, y_centre, xhalf, yhalf, $
                              zoom, zfac)
                xpix = coord[0]  &  ypix = coord[1]
                IF mode.ROLL THEN BEGIN
                    ns4 = state.NS4
                    IF xpix+ypix GT (state.IMSIZE[1] + ns4) THEN BEGIN
                        xpix -= ns4 & ypix -= ns4
                    ENDIF ELSE IF xpix + ypix LT state.NSIDE THEN BEGIN
                        xpix += ns4 & ypix += ns4
                    ENDIF
                ENDIF
                xpix = (xpix > 0) < xpmax
                ypix = (ypix > 0) < ypmax

                IF event.PRESS EQ 2 THEN BEGIN ; mark spot
                    mode.XPT = xpix
                    mode.YPT = ypix
                ENDIF

                done = 1 ; No gscroll update therefore done by definition!
            ENDELSE

; Print pixel details to readout and possibly terminal and logfile:
            logit = event.PRESS EQ 2

            mode.OXTV = xtv   &  mode.OYTV = ytv
            mode.XPIX = xpix  &  mode.YPIX = ypix

; End of code for zoom mode
        ENDIF ELSE IF event.PRESS THEN BEGIN ; Click in Overview mode...
            mode.OVERVIEW = 0   ; Switch to zoom mode...
            WIDGET_CONTROL, state.ZOOMCOL, /SENSITIVE
            WIDGET_CONTROL, state.READOUT, /SENSITIVE

            IF mode.RESAMP[0] GT 1 THEN BEGIN
                resamp = ROUND(mode.RESAMP)
                xpix = ((xtv - mode.CORNER[0])*resamp[0] > 0) < xpmax
                ypix = ((ytv - mode.CORNER[1])*resamp[1] > 0) < ypmax
            ENDIF ELSE BEGIN
                resamp = ROUND(1.0/mode.RESAMP)
                xpix = ((xtv - mode.CORNER[0])/resamp[0] > 0) < xpmax
                ypix = ((ytv - mode.CORNER[1])/resamp[1] > 0) < ypmax
            ENDELSE

            temp = (*tabarr).BYTE_PTR

            ierr = 0
            coord = gscroll(temp, xpix, ypix, xhalf, yhalf, mode.ZOOM_FACTOR, $
                            ierr, mode.NEW_VIEW, done, DO_WRAP=state.ROLL)

            lowerr = ierr MOD 8
            IF lowerr EQ 6 THEN BEGIN ; Tiny map: max zoom already.
                in = WIDGET_INFO(state.ZOOMCOL, FIND_BY_UNAME = 'IN')
                WIDGET_CONTROL, in, SENSITIVE = 0
                state.MAXZOOM = in
                WIDGET_CONTROL, event.TOP, SET_UVALUE = state
            ENDIF

            IF lowerr NE 0 && lowerr NE 6 THEN $
              MESSAGE, 'GSCROLL error '+STRING(ierr)

; rolling disabled if ierr = 8
            mode.ROLL = ierr NE 8 ? state.ROLL : 0
            mode.DONE = done

            xcurs = xhalf  &  ycurs = yhalf
            TVCRS, xcurs, ycurs

                                ; Disable panning if image fits in screen
            bltv = tv2im(0, 0, xpix, ypix, xhalf, yhalf, zoom, zfac)
            trtv = tv2im(!D.x_vsize, !D.y_vsize, xpix, ypix, xhalf, yhalf, $
                         zoom, zfac)

            mode.PAN =  MAX(bltv) GT 0 OR MAX(state.IMSIZE[1:2] - trtv) GT 1

            mode.OXTV  = xcurs    &  mode.OYTV  = ycurs
            mode.XHALF = xhalf    &  mode.YHALF = yhalf
            mode.X_CENTRE = xpix  &  mode.Y_CENTRE = ypix
            mode.XPIX = xpix      &  mode.YPIX = ypix
            mode.NEW_VIEW = 0

            logit = 0
        ENDIF ELSE RETURN, 0  ; overview mode, no click yet.

        IF logit && str.SCREEN NE state.LASTTAB THEN BEGIN
            WIDGET_CONTROL, str.BASE, GET_UVALUE = tname
            tname = 'Now on tab: ' + tname
            PRINT, tname
            PRINTF, state.LOGLUN, tname
            state.LASTTAB = str.SCREEN
            WIDGET_CONTROL, event.TOP, SET_UVALUE = state
        ENDIF

        pix_print, state, logit, start

    END                         ; End of "WIDGET_DRAW" case
    ELSE: MESSAGE, /INFORMATIONAL, 'Unknown event type received: '+name
ENDCASE

                                ; Mark current point
lutptr = str.LUT
TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
!P.background = (*lutptr).ABSENT
!P.color      = (*lutptr).LINE
IF mode.XPT GE 0 && ~mode.OVERVIEW THEN $
  marker, mode.XPT, mode.YPT, mode.X_CENTRE, $
  mode.Y_CENTRE, mode.XHALF, mode.YHALF, mode.ZOOM, mode.ZFAC

EMPTY

IF ~state.FOCUS THEN restore_lut, oldgraph

; If panels remain to be loaded, send timer event to DRAW window
; requesting re-draw:
IF done EQ 0 THEN WIDGET_CONTROL, event.ID, TIMER = 0.5

WIDGET_CONTROL, state.TABS, SET_UVALUE = mode

RETURN, returnable
END

;---------------------------------------------------------------------------
;
; Widget creation routines
;
;
PRO create_tab, base, tablab, uval, uname, tabarr, itab, xsize, ysize, lut
; Creates the widgets that populate a Ximview tab and fills in the mandatory
; fields of the tab descriptor structure
;
str = tabarr[itab]
str.LUT  = N_ELEMENTS(lut) GT 0 ? PTR_NEW(lut) : PTR_NEW(/ALLOCATE_HEAP)

str.BASE = WIDGET_BASE(base, TITLE = tablab, UNAME = uname, $
                       UVALUE = uval, /COLUMN, XPAD = 0, YPAD = 0, SPACE = 0)
str.DRAW = WIDGET_DRAW(str.BASE, XSIZE = xsize, YSIZE = ysize, RETAIN = 2, $
                       EVENT_FUNC = 'ximview_scroll', $
                       /BUTTON_EVENTS, /MOTION_EVENTS)
str.SCALE = WIDGET_DRAW(str.BASE, XSIZE = xsize, YSIZE = 45, RETAIN = 2)
str.SCREEN = itab
WIDGET_CONTROL, str.DRAW, GET_VALUE = index
str.WINDOW = index
WIDGET_CONTROL, str.SCALE, GET_VALUE = index
str.SCALE_INDEX = index

tabarr[itab] = str

END

PRO make_ximview, name, title, state, tlb, ntab, tablab
;
; Sets up widget geometry and installs defaults
;
COMPILE_OPT IDL2, HIDDEN

is_unix = STRCMP(!version.OS_FAMILY,'UNIX', 4,/ FOLD_CASE)

xsize = 512  & ysize = 512

tlb = WIDGET_BASE(TITLE = 'XIMVIEW', MBAR = menu, /COLUMN, $
                  UNAME = 'XIMVIEW', /TLB_SIZE_EVENTS, /KBRD_FOCUS_EVENTS)
file     = WIDGET_BUTTON(menu, VALUE = 'File', /MENU)
; options  = WIDGET_BUTTON(menu, VALUE = 'Options', /MENU)
display  = WIDGET_BUTTON(menu, VALUE = 'Display', /MENU)
frames   = WIDGET_BUTTON(menu, VALUE = 'Tabs', /MENU)
analysis = WIDGET_BUTTON(menu, VALUE = 'Analysis', /MENU)
help     = WIDGET_BUTTON(menu, VALUE = 'Help', /HELP, /MENU)

desc = ['0\Load FITS\ximview_fitsload', $
        '0\Load image file\ximview_rdimage', $
        '1\New logfile\ximview_newlog', '0\Overwrite old file', $
        '0\New sequence #', '2\Named...', $
        '0\Write PNG image\ximview_2png', $
        '0\Reset\ximview_reset', $
;        '0\Tab info\ximview_dump', $  ; Debug item
        '0\Exit\ximview_exit']
file_pd = CW_PDMENU(file, desc, /MBAR, /RETURN_NAME)
desc = ['1\Left Mouse', '0\Pan', '2\No action', $
        '1\Middle Mouse', '0\Mark Point', '0\Box Imstats', '0\Peakfit', $
        '2\No action', $
        '1\Right Mouse', '0\Overview','0\Zoom in', '0\Zoom out', $
        '0\Last zoom', '0\Last view button', '2\No action']
; opt_pd  = CW_PDMENU(options, desc, /MBAR, /RETURN_FULL_NAME)
cs   = colour_schemes(-1)
ncs  = N_ELEMENTS(cs)
cs   = ['0\'+cs[0:ncs-2], '2\'+cs[ncs-1]]
fun  = scale_funs(-1)
nfun = N_ELEMENTS(fun)
fun  = ['0\'+fun[0:nfun-2], '2\'+fun[nfun-1]]
desc = ['0\Adjust scaling\ximview_scale', $
        '0\Auto scale all tabs\ximview_autoscale_all', $
        '1\Colour table\ximview_newlut', cs, $
        '1\Colour handling\ximview_colour', $
        '0\Same for all', '2\Separate on each tab', $
        '0\Set view centre\ximview_goto', $
        '0\Restore default screen size\ximview_event']
disp_pd = CW_PDMENU(display, desc, /MBAR, /RETURN_FULL_NAME)
desc = ['0\Blink setup\ximview_blink', $
        '0\Red-Green-Blue\ximview_rgb', $
        '1\Polarization\ximview_pol', '2\Colour', $ ; '0\Vector', '2\LIC', $
        '1\Delete tab\ximview_deltab', '0\Current tab', '2\Specify']
fram_pd = CW_PDMENU(frames, desc, /MBAR, /RETURN_FULL_NAME)
desc = ['1\Imstats\ximview_imstat_opts', $
        '0\Box', '0\Region of Interest', '2\Set box size', $
        '1\Peakfit\ximview_maxfit_opts', $
        '0\Find extremum', '0\Find maximum', '0\Find minimum',  $
        '2\Set box size', $
        '0\Clear mark\ximview_clear_mark', $
        '1\Set image properties', $
        '0\polarization\ximview_setpol', '2\frequency\ximview_setfreq'] ;, $
;        '1\Grid\ximview_grid', '0\On/Off', '2\Grid interval',$
;        '0\Catalog\ximview_vo', '0\Profile', $
;        '0\Ruler', '0\Protractor', '0\Circles', '0\Mark direction', $
;        '1\Instrument FOVs','2\Planck']
ana_pd  = CW_PDMENU(analysis, desc, /MBAR, /RETURN_FULL_NAME)
desc = ['0\Help\ximview_help', '0\Release Notes\ximview_help', $
        '0\About\ximview_help']
help_pd = CW_PDMENU(help, desc, /MBAR, /HELP, /RETURN_FULL_NAME)


base1 = WIDGET_BASE(tlb, /ROW)

zoomcol = WIDGET_BASE(base1,/COLUMN, EVENT_FUNC = 'ximview_zoom', $
                      /ALIGN_CENTER)
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'Z')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'o')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'o')
zoomt = WIDGET_LABEL(zoomcol, VALUE = 'm')
zoomin    = WIDGET_BUTTON(zoomcol, VALUE = 'in', UVALUE = 'IN', UNAME = 'IN')
;,$ TOOLTIP = 'Zoom in x 2')
zoomreset = WIDGET_BUTTON(zoomcol, VALUE = '1:1', UVALUE = '1:1', $
                          TOOLTIP= 'Reset to 1 image pixel per screen pixel')
zoomout = WIDGET_BUTTON(zoomcol, VALUE = 'out', UVALUE = 'OUT', UNAME = 'OUT')
;,$ TOOLTIP = 'Zoom out x 2')


rightcol = WIDGET_BASE(base1,/COLUMN)
label    = WIDGET_LABEL(rightcol, value=name, /DYNAMIC_RESIZE)
                                ; Set up a separate tab with draw &
                                ; scale windows for each image to display
tabs     = WIDGET_TAB(rightcol, EVENT_PRO = 'ximview_tab')

null = PTR_NEW()
tab_str = {base: 0L, draw: 0L, scale: 0L, $     ; Widget IDs
           window: 0L, scale_index: 0L, $       ; Draw window indices
           im_ptr: null, $                      ; -> original image
           byte_ptr: null, $                    ; -> scaled images
           lut: null, $                         ; -> colour scale structure
           rgb: REPLICATE(null, 3), $           ; -> RGB screens, if any.
           absrange: [0.0, 0.0], $              ; Full intensity range
           sdev: 0.0, $                         ; estimated "noise" rms
           mode: 0.0, $                         ; modal intensity.
           range: [0.0, 0.0], $                 ; intensity range displayed
           beta: 1.0, $                         ; parameter for Asinh scaling
           zero: 0.0, $                         ; zero for display
           unit: 'unknown', $                   ; intensity unit
           mult: 1.0, $                         ; Multiplier for flux
           screen: 0, $                         ; tab #
           polcode: 0, $                        ; Polarization code
           freqcode: '', $                      ; Frequency descriptor
           wrap: 0, $       ; controls display of intensities outside "range"
           trfunc: 0, coltab: 0, $              ; display options
           bot: 0B, top: 0B, $                  ; byte range for LUT
           izero: 0B, $                         ; LUT level of zero
           decomposed: 0B, $                    ; is device decomposed?
           temporary: 0B, $                     ; Delete data at end
           collab: ['mono','',''], $            ; Labels for colour channels
           file: '', telescope: '', $           ; Components of name string
           instrument: '', object: '', creator: '', freq: '', stokes: ''}

tab_arr = REPLICATE(tab_str, ntab)
FOR itab = 0,ntab-1 DO create_tab, tabs, tablab[itab], tablab[itab], 'mono', $
  tab_arr, itab, xsize, ysize

IF is_unix THEN BEGIN
    readlab = WIDGET_LABEL(rightcol, value = title, $
                           /ALIGN_LEFT, /DYNAMIC_RESIZE)
    readout = WIDGET_LABEL(rightcol, value = 'No pixel assigned', $
                           /ALIGN_LEFT, /DYNAMIC_RESIZE)
ENDIF ELSE BEGIN
    font = 'lucida console*10'
    readlab = WIDGET_LABEL(rightcol, value = title, /ALIGN_LEFT, $
                           FONT = font, /DYNAMIC_RESIZE)
    readout = WIDGET_LABEL(rightcol, value = 'No pixel assigned', $
                           FONT = font, /ALIGN_LEFT, /DYNAMIC_RESIZE)
ENDELSE


button_row = WIDGET_BASE(tlb, /ROW)
overview   = WIDGET_BUTTON(button_row, VALUE = 'Overview', $
                           EVENT_PRO = 'ximview_overview', $
             TOOLTIP = 'Return to initial view of whole image in one panel,')

; Use pad1 as a convenient widget when we want explicitly to request scaling:
pad1       = WIDGET_BASE(button_row, EVENT_PRO = 'ximview_scale')
blink      = WIDGET_BUTTON(button_row, VALUE='Blink on/off', $
                           EVENT_PRO = 'ximview_blink')
imstats       = WIDGET_BUTTON(button_row, VALUE = 'Imstats', $
                              EVENT_FUNC = 'ximview_imstats', $
                    TOOLTIP = 'Prints statistics for box around marked point')
maxfit     = WIDGET_BUTTON(button_row, VALUE = 'Peakfit', $
                           EVENT_PRO  ='ximview_maxfit',  $
                       TOOLTIP = 'Fits peak near marked point')
; pad2 is another convenient peg:
pad2  = WIDGET_BASE(button_row, EVENT_PRO = 'ximview_deltab')
exit  = WIDGET_BUTTON(button_row, VALUE = 'Exit', EVENT_PRO = 'ximview_exit')


; Find out potential geometry and adjust a bit:
ogeom = WIDGET_INFO(overview,/GEOMETRY)
WIDGET_CONTROL, exit, XSIZE=ogeom.XSIZE

bgeom = WIDGET_INFO(base1, /GEOMETRY)
brgeom = WIDGET_INFO(button_row, /GEOMETRY)
xgeom = WIDGET_INFO(exit,/GEOMETRY)
length = xgeom.XOFFSET + xgeom.XSIZE + brgeom.XPAD
dx  = FIX(bgeom.XSIZE - bgeom.XPAD - length)
WIDGET_CONTROL, pad1, XSIZE = dx/2
WIDGET_CONTROL, pad2, XSIZE = dx - (dx/2)

mode = {drag: 0B, $             ; true if currently dragging
        pan: 1B,  $             ; true if panning enabled
        overview: 1B, $         ; true if in overview mode
        roll: state.ROLL, $     ; true if rolling is currently enabled
        new_view: 1B, $         ; true if re-load of pixmaps is needed
        done: 1B, $             ; true if all pixmaps up to date
        blink: 0B, $            ; true if blinking is going on
        tab: 0B, $              ; currently visible tab
        bbase: -1L, $           ; base widget of blink tab
        bwin:  -1L, $           ; Window index of main blink window
        bswin: -1L, $           ; Window index of blink scalebar window
        period: 0.5, $          ; current blink time/tab (seconds)
        do_scale: 1, $          ; blink scale as well as image
        zoom_factor: 1.0, $     ; current zoom factor
        zoom: 0, $              ; log_2 of zoom factor
        zfac: 1, $              ; zoom factor or 1/zoom factor
        resamp: [1.0, 1.0], $   ; resampling factors for overview mode.
        corner: [0, 0], $       ; screen coord of image BLC in overview mode
        x_centre: 0, y_centre: 0, $ ; image coord of pixel at screen centre
        xpix: 0, ypix: 0, $     ; image coords of pixel under cursor
        xpt: -1., ypt: -1., $   ; image coords of marked pixel
        xhalf: 0, yhalf: 0, $   ; screen coord of screen centre
        oxtv: xsize/2, oytv: ysize/2 $ ; screen coords of last cursor pos.
       }

tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
DEVICE, GET_SCREEN_SIZE = screen
xmargin = (screen[0]- tlb_geom.SCR_XSIZE)/2
ymargin = (screen[1]- tlb_geom.SCR_YSIZE)/2
xoff = 100 < xmargin

WIDGET_CONTROL, tlb, TLB_SET_XOFFSET = xoff
WIDGET_CONTROL, tlb, /REALIZE
WIDGET_CONTROL, tlb, SENSITIVE = 0   ; Wait until finished setting up
WIDGET_CONTROL, tabs, SET_UVALUE = mode

FOR itab = 0, ntab-1 DO BEGIN
    WIDGET_CONTROL, tab_arr[itab].draw, GET_VALUE = index
    tab_arr[itab].window = index
    WIDGET_CONTROL, tab_arr[itab].scale, GET_VALUE = index
    tab_arr[itab].scale_index = index
ENDFOR

*state.BLINK_SEQ = LINDGEN(ntab)

tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)

IF is_unix THEN BEGIN
    state.XSIZE     = tlb_geom.SCR_XSIZE ; Use SCR values as they include menu
    state.YSIZE     = tlb_geom.SCR_YSIZE ; bar and also the "size" returned by
    state.NEWWIDTH  = tlb_geom.SCR_XSIZE ; re-size events reflects SCR
    state.NEWHEIGHT = tlb_geom.SCR_YSIZE ; values under Linux.
    state.MBAR      = tlb_geom.SCR_YSIZE - tlb_geom.YSIZE
ENDIF ELSE BEGIN
    state.XSIZE     = tlb_geom.XSIZE ; But under MS WIndows they don't
    state.YSIZE     = tlb_geom.YSIZE
    state.NEWWIDTH  = tlb_geom.XSIZE
    state.NEWHEIGHT = tlb_geom.YSIZE
ENDELSE
state.TABS    = tabs     &  state.tabarr  = PTR_NEW(tab_arr)
state.READOUT = readout  &  state.LABEL   = label
state.READLAB = readlab  &  state.ZOOMCOL = zoomcol
state.PAD1    = pad1     &  state.PAD2    = pad2
state.BLINK   = blink    &  state.frames  = frames

END

PRO make_tabs, input, proj, column, roll, name, temporary, state, mode, str, $
               first, start, noldtab, data, header, newname, howto, tablab, $
               namestr, polcodes, ntab, line, title, extradims, mismatch
; Wrapper routine:
;   Interpret header (parse_header)
;   Create tab structures to describe tabs & fill in from header
;   Create tabs on widget
;   Store data in widget uservalues and/or the heap.
;   Set up or update GSCROLL
;   Set up mode if this is the first load.
;
; Outputs:
;  ntab, line, title, extradims, mismatch
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

verbose = state.VERBOSE

IF namestr.PATH NE '' THEN state.PATH = namestr.PATH

ntab = N_ELEMENTS(column)

IF N_ELEMENTS(name) EQ 0 THEN name = ''
log_header = newname + ' ' + name
PRINTF, state.LOGLUN, log_header, FORMAT = "('Dataset: ',A)"
IF first && name EQ '' THEN name = newname

IF howto LT 2 THEN BEGIN        ; swap to pointer usage if possible
    IF ntab EQ 1 THEN BEGIN
        temp = PTR_NEW(/ALLOCATE_HEAP)
        CASE howto OF
            0: *temp = input
            1: *temp = TEMPORARY(data)
        ENDCASE
        data  = temp
        howto = 2
        temporary = 1B ; we can delete heap variable
    ENDIF ELSE IF ~first THEN BEGIN
        used = PTR_VALID((*state.TABARR).IM_PTR)
        IF MIN(used) EQ 0 THEN $
          MESSAGE, 'Internal error: bulk data slot already in use'
    ENDIF
ENDIF

CASE howto OF
    0: T = SIZE(input) ; input data was 2+D array
    1: T = SIZE(data)  ; input data converted to 2-D array
    2: T = SIZE(*data) ; 2-D array is stored via a pointer
    3: T = SIZE(*data[0]) ; array of 2-D images via pointer array
    ELSE: MESSAGE, 'Unknown howto value!'
ENDCASE

IF howto EQ 3 THEN BEGIN
    extras = N_ELEMENTS(data)
    IF extras EQ 1 THEN howto = 2 ELSE BEGIN
        IF T[0] NE 2 THEN MESSAGE, 'Pointer array should be to 2-D images'
    ENDELSE
ENDIF

extradims = T[0] GT 2
intype    = T[T[0]+1]

IF extradims THEN BEGIN
    IF howto NE 0 THEN $
      MESSAGE, STRING(howto, FORMAT = $
                      "('Internal error: howto =',I2,' and >2 dimensions')")
    T = SIZE(input[*,*,0])
ENDIF

CASE howto < 2 OF
    0: WIDGET_CONTROL, state.LABEL, SET_UVALUE = input
    1: WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
    2: ; Data stored via pointer in state structure.
ENDCASE

; Extract info from header
parse_header, T, header, column, roll, state.VERBOSE, $
  astrom, is_astrom, csystem, proj, unit, title, nside, ns4

IF first THEN BEGIN
    line = ['']
    IF is_astrom THEN line = [line, 'Coordinate system: '+csystem, '']
    IF nside GT 0 THEN line = [line, 'Seems to be HEALPix with N_side' + $
                               STRING(nside,FORMAT="(I5)"), '']

                                ; Store info in state structure
    state.TITLE  = title     &  state.PROJ  = proj   &  state.ROLL = roll
    state.IMSIZE = T         &  state.NSIDE = nside  &  state.NS4  = ns4
    state.CSYSTEM = csystem  &  state.IS_ASTROM = is_astrom
    IF SIZE(astrom, /TYPE) NE 0 THEN state.ASTROM = PTR_NEW(astrom)
;    IF nside NE 0 THEN state.RING0 =
    mismatch = 0B
ENDIF ELSE BEGIN             ; Check for agreement with existing maps:
    mismatch = ~match_tab(state, T, proj, is_astrom, astrom, csystem)
    IF mismatch THEN BEGIN
        IF howto GE 2 && temporary THEN  FOR i=0,ntab-1 DO $
          IF PTR_VALID(data[i]) THEN PTR_FREE, data[i]
        RETURN
    ENDIF
ENDELSE
astrom = 0

; Make tab array for new tabs and fill in standard tab array entries
str.MULT = 1.0  &  str.TRFUNC = 0
str.SDEV = 0.0  &  str.MODE   = 0.0
str.ZERO = 0.0  &  str.IZERO  = 0    &  str.BETA  = 1.0
str.ABSRANGE = [0.0, 0.0]            &  str.RANGE = [0.0, 0.0]
str.COLLAB = ['mono','','']          &  str.RGB   = PTRARR(3)

; Set colour table (same as last used) and remember it:
ximview_lut, str.COLTAB, str.IZERO, decomp, bot, top
TVLCT, r, g, b, /GET
lutstr = {r:r, g:g, b:b, line: !P.color, absent: !P.background}

str.BOT = bot  &  str.TOP = top  &  *str.LUT = lutstr
str.DECOMPOSED = decomp


temp    = PTRARR(ntab, /ALLOCATE_HEAP)
newtabs = REPLICATE(str, ntab)
lutstr = *str.LUT

; Fill in state parameters for each tab:
IF howto LE 1 THEN temporary = 1B ; always safe to delete stored version.
fill_tab, newtabs, namestr, unit, polcodes, TEMPORARY = temporary
newtabs.BYTE_PTR = temp
newtabs.IM_PTR =  howto LT 2 ? PTR_NEW() : data

                        ; Update top title and tab labels, and
                        ; merge old & new arrays
IF first THEN topname = name ELSE $
  update_labels, state, name, newname, newtabs, tablab, topname
tab_arr = TEMPORARY(newtabs)    ; Now contains old and new

WIDGET_CONTROL, state.LABEL, SET_VALUE = topname

                        ; Change labels of old tabs:
FOR itab=0,noldtab-1 DO BEGIN
    WIDGET_CONTROL, tab_arr[itab].BASE, BASE_SET_TITLE = tablab[itab]
    WIDGET_CONTROL, tab_arr[itab].BASE, SET_UVALUE = tablab[itab]
ENDFOR

nnewtab = first ? ntab-1 : ntab

WSET, str.WINDOW
xsize = !D.x_vsize  & ysize = !D.y_vsize
FOR icol = 0,nnewtab-1 DO BEGIN    ; Make the new tabs
    itab = icol + noldtab
    create_tab, state.TABS, tablab[itab], tablab[itab], 'mono', $
      tab_arr, itab, xsize, ysize, lutstr
ENDFOR

IF first THEN BEGIN
    xpanel = 128L            ; default size of panel on pixmap.
    maxzoom = xpanel / 4L

    noldtab = 0              ; we want to start with tab 0 from now on
                                ; Set graphics focus to main panel:
    WSET, tab_arr[0].WINDOW
                                ; Initialize gscroll common and pixmap windows:
    ierr = 0
    gscroll_setup, ntab, maxwin, ierr, WINDOW = tab_arr.WINDOW, $
      IMAGE = tab_arr.BYTE_PTR, LUT = tab_arr.LUT, REDRAW = redraw_req, $
      /HIDDEN
    IF ierr NE 0 THEN MESSAGE, 'GSCROLL_SETUP error ' + STRING(ierr)
    state.MAXWIN = maxwin

    IF verbose THEN MESSAGE,/INFORMATIONAL, $
      STRING(SYSTIME(1)-start, "Scrolling set up", $
             FORMAT = "(F7.3,' seconds: ', A)")
    IF verbose THEN HELP, /MEMORY

; Set up for initial overview:
    xpix = T[1]/2  &  ypix = T[2]/2
    x_centre = xpix  &  y_centre = ypix

; Set initial zoom: largest that gets whole map on screen, or 1,
; whichever is larger.
    WSET, tab_arr[0].WINDOW
    ratio = (!D.x_vsize/T[1]) < (!D.y_vsize/T[2])
    ratio = ratio < maxzoom
    zoom = ratio GT 0 ? FIX(ALOG(ratio) / ALOG(2.0)) : 0
    zoom_factor = 2.^zoom

; Find effective central pixel on view window:
    get_centre, zoom_factor, xhalf, yhalf

    WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

    mode.ZOOM     = zoom      &  mode.ZOOM_FACTOR = zoom_factor
    mode.ZFAC     = 2^ABS(zoom)
    mode.OXTV     = xhalf     &  mode.OYTV     = yhalf
    mode.XPIX     = xpix      &  mode.YPIX     = ypix
    mode.XHALF    = xhalf     &  mode.YHALF    = yhalf
    mode.X_CENTRE = x_centre  &  mode.Y_CENTRE = y_centre

ENDIF ELSE FOR icol = 0,ntab-1 DO BEGIN
                                ; Update new low-level structures:
    itab = icol + noldtab
    WSET, tab_arr[itab].WINDOW
    gscroll_addscreen, tab_arr[itab].BYTE_PTR, tab_arr[itab].LUT
ENDFOR

*state.TABARR = tab_arr

END

PRO scale_tabs, ntab, noldtab, column, state, auto_scale, scale_pars, howto, $
                extradims, input, range, start
; Scales a bunch of new images
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

verbose = state.VERBOSE
fmt = "(F7.3,' seconds: ',A)"

IF N_ELEMENTS(range) EQ 0 && ~auto_scale THEN range = 'Full'
                                ; tells scale_image to use full data range
tabarr = state.TABARR

FOR icol = 0,ntab-1 DO BEGIN
    itab = icol + noldtab
    incol = column[icol] - 1
    str = (*tabarr)[itab]
    bptr = str.BYTE_PTR  &  iptr = str.IM_PTR
    bot = str.BOT          &  top = str.TOP    &  wrap = str.WRAP
    & trfunc = str.TRFUNC  &  zero = str.ZERO  &  beta = str.BETA

    IF auto_scale THEN BEGIN
        set_scale, scale_pars[*, icol], str
        range = str.RANGE
    ENDIF

    CASE howto < 2 OF
        0: *bptr = extradims ? scale_image(input[*,*,incol], range, wrap, $
                                           badcol, bot, top, ABS_RANGE = ar, $
                                           trfunc, zero, beta, /VERBOSE) $
                             : scale_image(input, range, wrap, badcol, $
                                           bot, top, ABS_RANGE = ar, trfunc, $
                                           zero, beta, /VERBOSE)
        1: BEGIN
            WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, /NO_COPY
            *bptr = scale_image(data, range, wrap, badcol, bot, top, trfunc, $
                                ABS_RANGE = ar, zero, beta, /VERBOSE)
            WIDGET_CONTROL, state.LABEL, SET_UVALUE = data, /NO_COPY
        END
        2: *bptr = scale_image(*iptr, range, wrap, badcol, bot, top, trfunc, $
                               ABS_RANGE = ar, zero, beta, /VERBOSE)
    ENDCASE

    str.ABSRANGE = ar
    str.RANGE    = range

    absmax = MAX(ABS(ar))
    test = NUMUNIT(absmax, str.UNIT, MULTIPLIER = mult, OUT_UNIT = ounit, $
                   /FORCE)
    str.UNIT = ounit
    str.MULT = mult

    (*tabarr)[itab] = str
ENDFOR

IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING(SYSTIME(1)-start, 'Images scaled', FORMAT = fmt)
IF verbose THEN HELP, /MEMORY

                               ; Set off-sky pixels to absent:
IF state.NSIDE NE 0 && state.IS_ASTROM THEN BEGIN
    fill_gores, state.NSIDE, state.IMSIZE, state.ASTROM, (*tabarr).BYTE_PTR

    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING(SYSTIME(1)-start, "absent pixels noted", FORMAT = fmt)
ENDIF

END

PRO fill_screens, ntab, noldtab, tabarr, mode, first, state, start
; Fills DRAW windows of newly-created tabs
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

FOR i=0,ntab-1 DO BEGIN
    itab = i + noldtab
    IF redraw_req THEN BEGIN
        lutptr = (*tabarr)[itab].LUT
        TVLCT, (*lutptr).R, (*lutptr).G, (*lutptr).B
        !P.background = (*lutptr).ABSENT
        !P.color      = (*lutptr).LINE
    ENDIF
    update_screen, tabarr, itab, mode, done
ENDFOR

IF ~first THEN BEGIN
                                ; Switch to newly-loaded screen
    WIDGET_CONTROL, state.TABS, SET_TAB_CURRENT = noldtab

    WSET, (*tabarr)[0].WINDOW   ; doesn't matter which, just need size
    gscroll_newscreen, noldtab, *tabarr, mode.ZOOM_FACTOR, $
      mode.X_CENTRE, mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, mode.OVERVIEW

    mode.DONE = done
                                ; Request another go if loading not finished
    IF done EQ 0B THEN WIDGET_CONTROL, (*tabarr)[0].DRAW, TIMER = 0.5

                                ;  Mark current point
    IF ~mode.OVERVIEW THEN BEGIN
        IF mode.XPT GE 0 THEN marker, mode.XPT, mode.YPT, mode.X_CENTRE, $
          mode.Y_CENTRE, mode.XHALF, mode.YHALF, mode.ZOOM, mode.ZFAC

        pix_print, state, 0, start
    ENDIF
ENDIF

END

PRO make_rgb_tab, str, ntab, xsize, ysize, rgbptr, collab, tablab, $
                  start, state, mode
;
; Creates a tab for an RGB image, fills in the tab structure
; Fills in the image and labels the plot.
;
; Inputs
;     str:          Structure to describe tab to be created
;     xsize, ysize: Current size of main draw window
;     rgbptr:       Array of 3 pointers to the rgb byte images
;     collab:       Array of 3 labels for the rgb channels
;     tablab:       Tab title label
;     start:        Start systime of this event
;     state, mode:  Usual system structures
;
COMPILE_OPT IDL2, HIDDEN

str.DECOMPOSED = 1
str.RGB        = rgbptr
str.UNIT       = ''
str.COLLAB     = collab
; Cancel existing pointers
str.IM_PTR = PTR_NEW()  &  str.BYTE_PTR = PTR_NEW()
ramp = LINDGEN(255)
lutstr = {r:ramp, g: ramp, b: ramp, line: !P.color, absent: !P.background}
; For the sake of DirectColor visuals

tabarr = TEMPORARY(*state.TABARR)
ntab   = N_ELEMENTS(tabarr)
tabarr = [tabarr, str]
; Create a tab for the RGB image:
create_tab, state.TABS, tablab, tablab, tablab, tabarr, ntab, xsize, ysize, $
  lutstr

str = tabarr[ntab]
                                ; Update low-level structures:
WSET, str.WINDOW
gscroll_addscreen, str.BYTE_PTR, str.LUT, str.RGB

*state.TABARR = TEMPORARY(tabarr)
                                ; Display tab
WSET, str.WINDOW

zoomf    = mode.ZOOM_FACTOR
xhalf    = mode.XHALF     &  yhalf    = mode.YHALF
x_centre = mode.X_CENTRE  &  y_centre = mode.Y_CENTRE

gscroll_newscreen, ntab, *state.TABARR, zoomf, $
  x_centre, y_centre, xhalf, yhalf, done, mode.OVERVIEW

IF mode.OVERVIEW THEN BEGIN
    overview, rgbptr, zoomf, mode.XPIX, mode.YPIX, xhalf, yhalf, $
      x_centre, y_centre
ENDIF ELSE BEGIN  ;  Mark current point
    IF mode.XPT GE 0 THEN marker, mode.XPT, mode.YPT, x_centre, $
      y_centre, xhalf, yhalf, mode.ZOOM, mode.ZFAC

    pix_print, state, 0, start
ENDELSE
                                ; Enable blinking
*state.BLINK_SEQ = INDGEN(ntab+1)
IF ntab + 1 GE 2 THEN BEGIN
    WIDGET_CONTROL, state.BLINK,  /SENSITIVE
    WIDGET_CONTROL, state.FRAMES, /SENSITIVE
ENDIF
                             ; Switch to newly-loaded screen
WIDGET_CONTROL, state.TABS, SET_TAB_CURRENT = ntab

mode.DONE = done
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode
IF done EQ 0B THEN WIDGET_CONTROL, str.DRAW, TIMER = 0.5

END

PRO ximview, input, range, proj, order, COLUMN = column, EXTENSION = exten, $
             WRAP = wrap, ROLL = roll, NAME = name, $
             NPOLE = npole, SPOLE = spole, RING = ring, NESTED=nest, $
             LOG = log, TEMPORARY = temporary, VERBOSE = verbose
;+
; NAME:
;       XIMVIEW
;
; PURPOSE:
;       Inspection and basic analysis tool for large images including
;       those stored in HEALPix format. Supports FITS WCS coordinate
;       and unit descriptions. Runs transparently over X-connections
;       (but requires a fast link to operate effectively).
;
;       Multi-dimensional images are displayed on a set of tabbed
;       screens, which support blinking etc.
;
; CATEGORY:
;       Widgets.
;
; CALLING SEQUENCE:
;
;       XIMVIEW, Input, [Range, Proj, Order]
;
; INPUTS:
;       Input:  Any of the following:
;               o Image array (2 or more dimensions)
;               o HEALPix array (1 or more dimension with the first
;                 being a valid HEALPix size)
;               o Structure containing a FITS header as tag (0),
;                 plus a 2+D image as tag (1), or a set of 1D
;                 HEALPix arrays as tags (1), (2),..., or a set of
;                 pointers to 1D or 2D arrays as tags (1), (2),...
;                 Arrays on all tags must be the same size.
;                 CUT4 format arrays are accepted.
;               o Name of a FITS file containing an image (in the
;                 primary HDU or in an IMAGE extension) or a set of
;                 HEALPix arrays including CUT4 format (in a binary
;                 table extension). The ".fits" or ".FITS" extensions
;                 may be omitted.
;               o A single pointer to any of the above.
;               o An array of pointers to 1D HEALPix or 2D image
;                 arrays.
;
; OPTIONAL INPUTS:
;       Range:  Sets the range of image intensities to map to the
;               displayed colour table. Options are:
;               o  Scalar: maximum image intensity to plot. The
;                   minimum to plot is the minimum in the data.
;               o  [minimum, maximum]
;               o '*' or 'AUTO': auto-scale image based on mode and
;                 robust estimate of standard deviation (separately
;                 for each channel).
;               Default is full range of image.
;
;       Proj:   'GRID', 'NPOLE', or 'SPOLE' if HEALPix array is
;               supplied. Default: 'GRID' = HPX projection.
;
;       Order:  'NESTED' or 'RING' if HEALPix array is
;               supplied. Default: header value if any, otherwise
;               'RING'
;
; KEYWORD PARAMETERS:
;       NPOLE:  Set for 'NPOLE' projection (alternative to Proj input)
;
;       SPOLE:  Set for 'SPOLE' projection (alternative to Proj input)
;
;       RING:   Set for RING order   (alternative to Order input)
;
;       NESTED: Set for NESTED order (alternative to Order input)
;
;       COLUMN: Single value or array of:
;               o Numbers or names of the binary table column to read,
;               o the plane index (coordinate on 3rd dimension) if
;                 the input is a 3-D array. (Dimensions higher than 3
;                 are "collapsed" into the 3rd dimension).
;               o Coordinate on 2nd dimension if the input is an array
;                 of (1-D) HEALPix arrays
;               If not specified, all the columns/planes/arrays are
;               read in.
;
;       EXTENSION:  FITS extension containing data to read (primary
;               HDU is considered extension 0). Default: extension 0 if
;               it contains data, otherwise extension 1.
;
;       WRAP:   = 0 (unset): intensities outside Range are displayed
;                            as min or max colour, as appropriate.
;               < 0: Intensities greater than max set by Range use
;                    "wrapped" colours, starting again from the bottom
;                    of the colour table. Intensities less than the
;                    min set by Range are displayed as the min colour.
;                >0: colour table is wrapped at both ends;
;
;       ROLL:   Enable rolling of image through +/-180 longitude
;               (enabled automatically if program determines the image
;               is a full-sphere HPX projection).
;
;       NAME:   Title for image.
;               Default: constructs one from filename and/or header.
;
;       LOG:    Set (/LOG or LOG=1) to give logfile a unique name of
;               the form "ximview_n.log" where n is an integer. (n=m+1
;               where "ximview_m.log" is the file in the default
;               directory with the largest index m).
;               OR: set = <string>, the name of the logfile.
;
;       TEMPORARY: Set to overwrite input array to save space.
;
;       VERBOSE:   Set to print diagnostics
;
; OUTPUTS:
;       Logfile: see LOG above.
;
; COMMON BLOCKS:
;       GR_GLOBAL:    Contains global graphics state parameters
;       GRID_GSCROLL: Used by low-level GSCROLL library.
;       XY2PIX:       HEALPix low-level library. Only used when
;                     displaying HEALPix images.
;
; SIDE EFFECTS:
;       If your plot device Visual Class is "DirectColor" a private
;       colour map is (usually) used when the cursor is in Ximview's
;       image display area, causing all elements on your screen to
;       change colour while the cursor remains there. To avoid such
;       "flashing", set DEVICE, TRUE_COLOR=24 at the start of your
;       session.
;
;       A similar effect may occur for Visual Class "PseudoColor"
;       except that colours will change only when the widget gains or
;       loses "keyboard focus" (i.e. it becomes the active window).
;
; RESTRICTIONS:
;       Only one instance of Ximview is allowed at a time.
;
;       XIMVIEW uses the HEALPix IDL library for displaying HEALPix
;       images. For best behaviour the XIMVIEW directories should
;       occur before the HEALPix ones.
;
; PROCEDURE:
;       Ximview prints some basic instructions when it starts up.
;       Detailed instructions are available from the HELP button.
;
;       For a description of HEALPix data, and to download the
;       software see:
;       o  http://healpix.jpl.nasa.gov/
;       o  Gorski et al. 2005, Astrophysical Journal, vol 622, p. 759
;
;       HEALPix datasets are converted to 2D arrays in one of the
;       projections described by Calabretta and Roukema (2007, Monthly
;       Notices of the Royal Astronomical Society, vol 381, p. 865.)
;
; EXAMPLE:
;       Display test-card image, scaled min (=0) to 100:
;
;               XIMVIEW, DIST(200,200), 100
;
;       To display the first three HEALPix arrays stored in the first
;       extension to file 'wmap_band_iqumap_r9_3yr_K_v2.fits'
;       (available from http://lambda.gsfc.nasa.gov/), with auto-scaling:
;
;               XIMVIEW, 'wmap_band_iqumap_r9_3yr_K_v2', '*', COL=[1,2,3]
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, Jan-Feb 2008 (to v.0.3)
;       March 2008      v0.4: added RGB and HSV display, CUT4 files,
;                             better colour handling, numerous minor
;                             improvements (see Release Notes).
;       April 2008      v0.4.2: Bug fixes
;       July  2008      v0.5: Bug fixes
;       November 2008   v0.6: Bug fix (HP2HPX), added scale to PNG output.
;       August 2009     v0.6.2: Bug fix in grid2hp_index & related progs.
;-
COMPILE_OPT IDL2

; Global parameters describing graphics state
COMMON gr_global, windev, redraw_req, colmap, badcol, syscol

version = '0.6.2'

start = SYSTIME(1)

temporary = KEYWORD_SET(temporary)
restore_size = 0
ntab = 0

first = XREGISTERED('ximview', /NOSHOW) EQ 0

;ON_ERROR, 0                     ; For debug:
;error_status = 0
ON_ERROR, 2
CATCH, error_status
IF error_status NE 0 THEN BEGIN
; Should arrive here if error occurs before XMANAGER is called and
; starts to handle errors
    CATCH, /CANCEL
    HELP, /LAST_MESSAGE

    IF temporary && N_ELEMENTS(data) GT 0 THEN FOR i=0,ntab-1 DO $
      IF PTR_VALID(data[i]) THEN PTR_FREE, data[i]

    IF restore_size THEN $ ; Restore any dummy axes in input
      input = REFORM(input, input_size[1:input_size[0]], /OVERWRITE)

    IF first THEN BEGIN
        ximview_tidy, state

        IF N_ELEMENTS(tlb) GT 0 THEN WIDGET_CONTROL, tlb, /DESTROY

        MESSAGE, 'Problem initializing Ximview'
    ENDIF ELSE BEGIN            ; Clear up any pointers and pixmaps
                                ; related to new tabs
        tabarr = *state.TABARR
        ntab = N_ELEMENTS(tabarr)
        FOR i=noldtab,ntab-1 DO BEGIN
            id = i
            str = tabarr[i]
            deadtab = str.BASE
            index   = str.SCREEN
            WIDGET_CONTROL, deadtab, /DESTROY
            str.SCREEN = -1
            IF PTR_VALID(str.LUT) THEN PTR_FREE, str.LUT
            IF PTR_VALID(str.BYTE_PTR) THEN PTR_FREE, str.BYTE_PTR
            IF str.TEMPORARY && PTR_VALID(str.IM_PTR) THEN $
              PTR_FREE, str.IM_PTR

            current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
            gscroll_newscreen, current, tabarr, mode.ZOOM_FACTOR, $
              mode.X_CENTRE, mode.Y_CENTRE, mode.XHALF, mode.YHALF, done, 1B
        ENDFOR
        *state.TABARR = tabarr
        MESSAGE, /INFORMATIONAL, 'Problem loading new data'
        WIDGET_CONTROL, state.TABS, /SENSITIVE
    ENDELSE

    RETURN
ENDIF                           ; Catch block ends here

IF ~first THEN BEGIN
    MESSAGE, /INFORMATIONAL, $
        'Data will be added to existing Ximview window'
        ; Find Ximview

    test = WIDGET_INFO(/MANAGED)
    IF test[0] EQ 0 THEN MESSAGE, 'No widgets currently being managed'

    FOR i = 0,N_ELEMENTS(test)-1 DO BEGIN
        uname = WIDGET_INFO(test[i], /UNAME)
        IF uname EQ 'XIMVIEW' THEN GOTO, GOTIT
    ENDFOR

    MESSAGE, 'Ximview is not currently being managed'

GOTIT:

    tlb = test[i]
    WIDGET_CONTROL, tlb, GET_UVALUE = state
    WIDGET_CONTROL, state.TABS, GET_UVALUE = mode
    current = WIDGET_INFO(state.TABS, /TAB_CURRENT)
ENDIF

fmt = "('XIMVIEW: ',F7.3,' seconds: ',A)"

; Sort out input parameters

IF N_PARAMS() LT 1 THEN BEGIN
    PRINT, 'Syntax:'
    PRINT, 'XIMVIEW, input, [range, [proj, [order, ]]]'
    PRINT, '         [{/NPOLE | /SPOLE}, {/RING | /NESTED}'
    PRINT, '         COLUMN=, EXTENSION =,  WRAP=, ROLL=, NAME=, LOG='
    PRINT, '         /TEMPORARY, /VERBOSE]'
    RETURN
ENDIF

ring = KEYWORD_SET(ring)
nest = KEYWORD_SET(nest)
order_set = N_ELEMENTS(order) GT 0

CASE order_set + ring + nest OF
    0: ; Do nothing, may not be healpix
    1: IF ~order_set THEN order = nest ? 'NESTED' : 'RING'
    ELSE:  MESSAGE, 'Please specify ordering only once.'
ENDCASE

verbose = KEYWORD_SET(verbose)
IF verbose THEN HELP, /MEMORY

; Auto-scale range if requested:
auto_scale = 0B
IF SIZE(range,/TYPE) EQ 7 THEN BEGIN
    rtext = STRTRIM(STRUPCASE(range),2)
    SWITCH rtext OF
        '*':
        'AUTO': BEGIN
            auto_scale = 1B
            BREAK
        END
        ELSE: MESSAGE, /INFORMATIONAL, 'Option ' + rtext + $
          ' not yet available.'
    ENDSWITCH
ENDIF

IF first THEN BEGIN ; Initial setup: logfile, widget, common:
    is_unix = STRCMP(!version.OS_FAMILY, 'UNIX', 4, /FOLD_CASE)
    ; Find help directory
    info = ROUTINE_INFO('ximview', /SOURCE)
    helpdir = STRSPLIT(info.PATH, 'ximview.pro', /REGEX, /EXTRACT) $
        + 'docs'
    pathsep = path_sep(/SEARCH_PATH)
    dirs = STRSPLIT(!help_path, pathsep, /EXTRACT)
    ndir = N_ELEMENTS(dirs)
    IF dirs[ndir-1] NE helpdir THEN !help_path = !help_path + pathsep + helpdir

; Save current state in structure old_graph:
    windev = is_unix ? 'X' : 'WIN'
    swap_lut, dummy, old_graph

    IF ~KEYWORD_SET(log) THEN logfile = 'ximview.log' ELSE BEGIN
        ltype = SIZE(log,/TYPE) ; ltype = 7 for a string
        logfile = ltype EQ 7 ? log : get_log_name()
    ENDELSE
    OPENW, loglun, logfile, /GET_LUN
    PRINTF, loglun,  version, SYSTIME(), FORMAT = $
      "('XIMVIEW Version ',A,' started at ',A)"

; Inputs involving input data:

    roll = KEYWORD_SET(roll)

    proj_set = N_ELEMENTS(proj) GT 0
    npole = KEYWORD_SET(npole) & spole = KEYWORD_SET(spole)
    IF proj_set + npole + spole GT 1 THEN $
        MESSAGE, 'Please specify projection only once!'
    IF npole THEN proj = 'NPOLE'
    IF spole THEN proj = 'SPOLE'
    coltab = 0
    noldtab = 0

    wd = FILE_EXPAND_PATH('')
    title = {head: ' ', unit: ' ', tail: ' '} ; Dummy value
    state = {zoomcol: 0L, tabs: 0L, readout: 0L, $         ; Widget IDs
             label: 0L, readlab: 0L, pad1: 0L, pad2: 0L, $ ; Widget IDs
             blink: 0L, frames: 0L, $                      ; Widget IDs
             maxzoom: 0L, $                                ; Variable widget ID
             tabarr:  PTR_NEW(), $    ; -> structure array for each tab
             lasttab: -1L, $          ; last tab # referred to in log file
             title: title,       $    ; structure with elements of read label
             loglun: loglun,     $    ; Logical unit number for log file
             version: version,   $    ; Ximview version #
             proj: '', roll: roll, verbose: verbose, $ ; Input parameters
             focus: 1B, $             ; True if we have keyboard focus
             path: wd, $              ; Path for data (not always log files)
             maxwin: [0S, 0S], $      ; max size of draw windows
             global_colour: 1B, $     ; True if all tabs should use same LUT
             old_graph: PTR_NEW(), $  ; -> old graphics state
             blink_seq: PTR_NEW(-1), $ ; pointer to blink sequence array
             xsize: 0S, ysize: 0S, mbar: 0S, $    ; tlb geometry
             newwidth: 0S, newheight: 0S, $       ; tlb geometry
             imsize: LONARR(5), nside: 0L, ns4: 0L, $ ; image geometry
             ring0: PTR_NEW(/ALLOCATE_HEAP), $    ; -> scratch array
             roi: 0B, statbox: [33S, 33S], $      ; Parameters for Imstats
             maxbox: 7S, peak: 1S, $              ; Parameters for  Peakfit
             is_astrom: 0B, $                ; generic astrometry available
             astrom: PTR_NEW(), csystem: ''} ; astrometry

; Remember old graphics device and colour table:
    state.OLD_GRAPH = PTR_NEW(old_graph)
    make_ximview, 'Welcome to Ximview', ' ', state, tlb, 1, 'Loading...'
    str = (*state.TABARR)[0]
    noldtab = 1

;set up gr_global common:
    sys_cols = WIDGET_INFO(tlb, /SYSTEM_COLORS)
    nsyscol = N_TAGS(sys_cols)
                                ; Find grey levels amoung system colours
    cols = [0]
    FOR i=0,nsyscol-1 DO BEGIN
        col = sys_cols.(i)      ; NB: integer triplet not byte triplet
        IF col[0] NE -1 THEN BEGIN
            test = WHERE(col[0] EQ cols)
            IF test EQ -1 THEN cols = [cols, col[0]]
        ENDIF
    ENDFOR
    sys_cols = 0
    cols  = cols[ UNIQ( cols, BSORT(cols) ) ]
    IF verbose THEN PRINT, 'System grey-levels:', cols

    IF is_unix THEN greys = WHERE(cols NE 0B AND cols NE 255B, ngrey) $
               ELSE ngrey = 0
    IF ngrey GT 0 THEN greys = cols[greys]

; Set grey levels for non-data pixels. colour index = grey level so
; that the CI does not need to be changed for decomposed (RGB) images:
    IF ngrey LT 1 THEN greys = [106] ; Default off-sky level
    IF ngrey LT 2 THEN greys = [greys, 192] ; Default bad pixel level
    syscol = BYTE(greys)

; Find out about colour maps: Do we have one? Do we have to
; re-draw image if the map is changed?
    DEVICE, DECOMPOSED = 0  ; Affects results under TrueColor
    colmap     = COLORMAP_APPLICABLE(redraw_req)
    redraw_req = colmap && redraw_req

    IF ~redraw_req THEN greys = greys[0:1] ; only keep non-data levels.

; Set default colour table and remember it:
;    ximview_lut, coltab, 0, decomp, bot, top
;    TVLCT, r, g, b, /GET
;    lutstr = {r:r, g:g, b:b, line: !P.color, absent: !P.background}
;
;    str.BOT = bot  &  str.TOP = top  &  *str.LUT = lutstr

    IF verbose THEN PRINT, SYSTIME(1)-start, "Graphics state set", FORMAT=fmt
    IF verbose THEN HELP, /MEMORY

ENDIF ELSE BEGIN                ; Loading new data into existing widget
    proj   = state.PROJ
    roll   = state.ROLL

    str = (*state.TABARR)[0]
    prep_screen, state, mode, otabarr, old_graph
    noldtab = N_ELEMENTS(otabarr)
    otabarr = 0
                                ; grab graphics focus if we don't
                                ; already have it:
    IF ~state.FOCUS THEN BEGIN
        state.FOCUS = 1
        WIDGET_CONTROL, tlb, SET_UVALUE = state
        WIDGET_CONTROL, str.DRAW, /INPUT_FOCUS
        *state.OLD_GRAPH = old_graph
    ENDIF
ENDELSE

IF N_ELEMENTS(wrap) NE 0 THEN str.WRAP = wrap

;WIDGET_CONTROL, /HOURGLASS
; Suppress any dummy axes in input
input_size = SIZE(input)
IF input_size[0] GT 0 THEN input = REFORM(input,/OVERWRITE)
restore_size = ~temporary && input_size[0] GT 0

                                ; Parse the input parameter
parse_input, input, ORDER = order, PROJ = proj, TEMPORARY = temporary, $
  COLUMN = column, data, header, newname, howto, tablab, scale_pars, $
  namestr, polcodes, GET_SCALE = auto_scale, EXTENSION = exten, $
  VERBOSE = verbose

IF verbose THEN PRINT, SYSTIME(1)-start, "Input converted", FORMAT=fmt
IF verbose THEN HELP, /MEMORY

; create widget tab(s), and stash data in widget uservalues/heap
make_tabs, input, proj, column, roll, name, temporary, state, mode, $
  str, first, start, noldtab, data, header, newname, howto, $
  tablab, namestr, polcodes, ntab, line, title, extradims, mismatch

IF mismatch THEN MESSAGE, 'New data does not match data already loaded'


; Scale the image
scale_tabs, ntab, noldtab, column, state, auto_scale, scale_pars, howto, $
  extradims, input, range, start

tabarr = state.TABARR

; Draw initial screens:
fill_screens, ntab, noldtab, tabarr, mode, first, state, start

                                ; Update readout label
title_string = title.HEAD + form_unit((*tabarr)[noldtab].UNIT) + title.TAIL
WIDGET_CONTROL, state.READLAB, SET_VALUE = title_string

IF first THEN BEGIN
    PRINT, ''
    PRINT, 'XIMVIEW:  Click cursor to define centre of zoom/scroll, then:'
    PRINT, '          Mouse button 1 to drag image,'
    PRINT, '                       2 to mark point and record value'
    PRINT, 'Maxfit and Imstats buttons work around last marked point.'
    PRINT, ''

    IF verbose THEN PRINT, SYSTIME(1)-start, "Starting event loop", FORMAT=fmt
    IF verbose THEN HELP, /MEMORY

    line = [line, title_string]
    PRINTF, state.LOGLUN, line, FORMAT="(A)"
    PRINT, line, FORMAT="(A)"
ENDIF

; Store global data
WIDGET_CONTROL, tlb, SET_UVALUE = state
WIDGET_CONTROL, state.TABS, SET_UVALUE = mode


; Enable buttons etc:
*state.BLINK_SEQ = INDGEN(ntab+noldtab)
WIDGET_CONTROL, tlb, /SENSITIVE
IF mode.OVERVIEW THEN WIDGET_CONTROL, state.ZOOMCOL, SENSITIVE = 0
sens =  ntab + noldtab GE 2

WIDGET_CONTROL, state.BLINK,  SENSITIVE = sens
WIDGET_CONTROL, state.FRAMES, SENSITIVE = sens


CATCH, /CANCEL

IF first THEN XMANAGER, 'ximview', tlb, CLEANUP='ximview_cleanup', /NO_BLOCK $
         ELSE WIDGET_CONTROL, state.TABS, /SENSITIVE

END
