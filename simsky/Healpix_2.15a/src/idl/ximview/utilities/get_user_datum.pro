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
FUNCTION get_user_datum, event, question, field_size, default
;+
; NAME:
;       GET_USER_DATUM
;
; PURPOSE:
;       Launches a modal dialog to get a piece of info from the user
;
; CATEGORY:
;       Widget Routines, Compound
;
; CALLING SEQUENCE:
;
;       Result = GET_USER_DATUM(Event, Question, Field_size, Default)
;
; INPUTS:
;       Event:      Unused. Present for compatibility with old version.
;       Question:   Text for label
;       Field_size: Size of text box for answer, in characters.
;
; OPTIONAL INPUTS:
;       Default:    Default value placed in box.
;
; OUTPUTS:
;       Result is either the contents of the text box or 'Cancel' if
;       Cancel button was pressed.
;
; EXAMPLE:
;
;       name = GET_USER_DATUM('Please enter your name',30)
;
; MODIFICATION HISTORY:
;       Written by:      J P Leahy March 2008 
;                        (complete re-write of earlier version).
;-
COMPILE_OPT IDL2, HIDDEN

; Escape backslashes and commas in the input strings:
question = STRJOIN(STRSPLIT(question,'\',/EXTRACT,/PRESERVE_NULL),'\\')
question = STRJOIN(STRSPLIT(question,',',/EXTRACT,/PRESERVE_NULL),'\,')

IF N_ELEMENTS(default) EQ 0 THEN default = ''
default = STRJOIN(STRSPLIT(default,'\',/EXTRACT,/PRESERVE_NULL),'\\')
default = STRJOIN(STRSPLIT(default,',',/EXTRACT,/PRESERVE_NULL),'\,')

desc = ['0, LABEL,'+question+' ,', $
        '0, TEXT, '+default+', TAG=answer, WIDTH='+STRING(field_size), $
        '1, BASE,, ROW', '0, BUTTON, Accept, QUIT, TAG=OK', $
        '2, BUTTON, Cancel, QUIT']

result = CW_FORM(desc, /COLUMN, TITLE = 'Enter data')
RETURN, result.OK ? result.ANSWER : 'Cancel' 

END

