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
PRO tvdisp, image, resamp, range, tvpos, WRAP=wrap, VERBOSE=verbose
;+
; NAME:
;       TVDISP
;
; PURPOSE:
;       Displays an image on the TV, optionally resampled and scaled,
;       by default centred in the field of view.
;
; CATEGORY:
;       Direct Graphics, Image Display
;
; CALLING SEQUENCE:
;
;       TVDISP, Image, Resamp, Range, TVpos
;
; INPUTS:
;       Image:   2-D array to display (any numerical data type). If
;                more than 2 dimensions, only the first 2-D "slice" is
;                displayed. 
;
; OPTIONAL INPUTS:
;       Resamp:  positive integer or array of two integers giving
;                resampling factor, e.g. if Resamp = 3 every third row
;                and column is displayed, if = [3,6] every third row
;                and sixth column. Default: [1,1]
;
;       Range:   One or two values specifying range of intensity. If
;                two values are given, the bottom of the colour scale
;                is mapped to Range[0] and the top to Range[1]. If one
;                value is given, the range is Data Min to Range. If
;                not specified, the range is Data Min to Data Max.
;
;       TVpos    One or two integers specifying position on the TV
;                screen. One value is passed to TV as its "position"
;                argument, e.g. TVpos = 0 puts the image in the top
;                left corner of he window (if it fits). If two values
;                are specified they are the screen pixel coordinates
;                of the image bottom left corner. If not specified,
;                the image is centred.
;
; KEYWORD PARAMETERS:
;       WRAP:   = 0 (unset): intensities outside Range are displayed
;                            as min or max colour, as appropriate.
;               < 0: Intensities greater than max set by Range use
;                    "wrapped" colours, starting again from the bottom
;                    of the colour table. Intensities less than the
;                    min set by Range are displayed as the min colour.
;                >0: colour table is wrapped at both ends;
;
;       VERBOSE: Prints max and min of image before and after scaling.
;
; SIDE EFFECTS:
;       If no graphics window is available, a new 512x512 pixel window
;       is created.
;
; EXAMPLE:
;
;       TVDISP, DIST(1024,1024), 2
;
; MODIFICATION HISTORY:
;       Written by:  J. P. Leahy November 2007
;
;-
COMPILE_OPT IDL2

verbose = KEYWORD_SET(verbose)

; Get size of display window

IF !D.window EQ -1 THEN WINDOW, xsize=512, ysize=512
xsize = !D.x_vsize
ysize = !D.y_vsize

IF verbose THEN PRINT, MINMAX(image, /NAN), $
  FORMAT="('Minimum & Maximum on original image:', 2(1X,E10.3))"

rx = 1 & ry = 1
S = SIZE(resamp)
IF S[1] GT 0 THEN BEGIN
    IF S[0] GT 0 AND S[1] GT 1 THEN BEGIN
        rx = FIX(resamp[0])
        ry = FIX(resamp[1])
    ENDIF ELSE BEGIN
        rx = FIX(resamp[0])
        ry = rx
    ENDELSE
ENDIF 

wrap = KEYWORD_SET(wrap)
R = SIZE(image[*,*]) ; Implicitly restricted to first plane if image is 3D

ntvp = N_ELEMENTS(tvpos)

docorner = 1
x0 = 0 & y0 = 0
CASE ntvp OF
    2: BEGIN
        xout = !D.x_vsize - tvpos[0] + 1
        yout = !D.y_vsize - tvpos[1] + 1
        x1 = xout*rx < R[1]-1 & y1 = yout*ry < R[2]-1
        xcorner = tvpos[0] & ycorner = tvpos[1]
    END
    1: BEGIN
        xout = !D.x_vsize
        yout = !D.y_vsize
        x1 = xout*rx < R[1]-1 & y1 = yout*ry < R[2]-1
        docorner = 0
    END
    0: BEGIN                     ; Centre image on screen
        xout = !D.x_vsize & yout = !D.y_vsize
        x0 = (R[1] - xout*rx)/2 > 0
        x1 = (R[1] + xout*rx)/2 < R[1]-1
        y0 = (R[2] - yout*ry)/2 > 0
        y1 = (R[2] + yout*ry)/2 < R[2]-1
        xcorner = (!D.x_vsize - (x1-x0+1)/rx )/2
        ycorner = (!D.y_vsize - (y1-y0+1)/ry )/2
    END
ENDCASE

; Colour table?
temp = S[1] EQ 0 ? image[x0:x1,y0:y1] $
                 : image[x0:x1:rx,y0:y1:ry]

temp = scale_image(temp, range, wrap, VERBOSE=verbose) 

IF  docorner THEN TV, temp, xcorner, ycorner ELSE TV, temp, tvpos

END
