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
PRO gscroll_addscreen, imptr, lutptr, rgbptr
;
; Adds new pixmap and screen structure
;
; Inputs:
;    imptr:  Pointer to byte image for the new screen
;    lutptr: Pointer to structure describing colour Look-up-table for
;            new screen
;    rgbptr: Array of 3 pointers to R, G, B channels (only set for
;            3-colour screens).
;
COMPILE_OPT IDL2, HIDDEN
COMMON GRID_GSCROLL

null = WHERE(screens.viewwin EQ old_window, old)

view_x = !D.x_vsize  &  view_y = !D.y_vsize

str = screens[0]
IF ~old THEN WINDOW, /FREE, TITLE = 'GSCROLL 0', XSIZE = view_x, YSIZE = view_y
str.viewwin = !D.window
WINDOW, /FREE, XSIZE = xwidth, YSIZE = ywidth, /PIXMAP

str.PIXWIN = !D.window
str.PIXMAP_panel.idx = -1
str.STARTED = 1
str.IMPTR   = imptr
str.LUT     = lutptr
str.IS_RGB  = N_ELEMENTS(rgbptr) GT 0
IF str.IS_RGB THEN str.RGB = rgbptr

screens = [screens, str]

END
