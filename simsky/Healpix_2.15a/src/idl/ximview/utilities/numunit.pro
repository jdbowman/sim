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
;  Module NUMUNIT
;  
;  J. P. Leahy 2008
;
;  Contents:
;  get_power:  utility routine
;  numunit:    main function.
; 
;  A standard IDL documentation block follows the declaration of
;  numunit.
;
PRO get_power, string, power, irest, valid
;
; Gets the exponent of a valid FITS unit "raised to the power of" expression
;
; Input: string: contains unit from the point the "power of" part starts
; Output: irest: first character position in string after the
;                expression
;         power: the power
;         valid: 2 if there is a valid integer power expression
;                1 if there is valid fractional power expression
;                0 if there is not a valid expression
;
COMPILE_OPT IDL2, HIDDEN

irest = 0  & power = 0 & valid = 0
IF STRCMP(string,'**',2) THEN ksta = 2 ELSE $
  IF STRCMP(string,'^',1) THEN ksta = 1 ELSE ksta = 0

brace = STRCMP(STRMID(string,ksta),'(',1)
IF brace THEN BEGIN
    ksta = ksta + 1
    kend = STREGEX(string,')') - 1
    IF kend LT ksta THEN $
        MESSAGE, /INFORMATIONAL, 'Mismatched parentheses in unit string'

    ktext = STRMID(string,ksta,kend-ksta+1)
    dot = STREGEX(ktext,'\.',/BOOLEAN)
    solidus = STREGEX(ktext,'/')
    ratio = solidus GT -1
    irest = kend + 2
ENDIF ELSE BEGIN
    dot = 0  & ratio = 0
; Find first char not 0 to 9 or + or -
    klen = STREGEX(STRMID(string,ksta),'[^0-9+-]') 
    ktext = STRMID(string,ksta,klen)
    irest = ksta + klen
ENDELSE
valid = 2 - dot - ratio

IF irest EQ 0 OR valid EQ 0 THEN BEGIN
    valid = 0
    irest = 0
    RETURN
ENDIF

CATCH, error_code
IF error_code NE 0 THEN BEGIN ; Problem converting string.
    CATCH, /CANCEL
    power = 0
    valid = 0
    irest = 0
ENDIF

IF valid EQ 2 THEN power = FIX(ktext) ELSE IF dot THEN $
  power = FLOAT(ktext) ELSE BEGIN
    p1 = FLOAT(STRMID(ktext,0,solidus))
    p2 = FLOAT(STRMID(ktext,solidus+1))
    power = p1/p2
ENDELSE

RETURN

END

FUNCTION numunit, numin, unit, PRECISION = precision, DECIMALS = decimals, $
    FORCE = force, MULTIPLIER = multiplier, OUT_UNIT = out_unit
;+
; NAME:
;       NUMUNIT
;
;
; PURPOSE: 
;       Parses (number, unit) pairs to a string suitable for
;       printing. A multiplier prefix is chosen for the unit to put
;       the number in the range [1,1000). Assumes unit is encoded in
;       the FITS convention (See the FITS standard V3.0 at
;       fits.gsfc.nasa.gov), generalised in that

;       o  "<unit>s" is acceptable as well as "<unit>" 
;          (unless unit is a single character),
;       o  Spelled-out names of units with short (usually one-letter)
;          symbols are allowed.
;       o  unit names which have been upper-cased are allowed,
;          provided they are not too ambiguous. For upper-case unit
;          descriptors, we also search for spelled-out multiplier
;          prefixes, eg. "MILLI".
;       o  Commas are an acceptable separator between unit elements.
;          (The implication is that the text after the comma is a
;          comment, as in "mK,thermodynamic")
;       If the leading component of the unit (E.g. "count" in
;       "count/s") is one listed in the FITS standard as being
;       inappropriate for multiplier prefixes, or if it is not
;       recognised, the unit is returned unchanged unless /FORCE is
;       specified, in which case a numerical prefix is supplied.
;
;       Optional outputs give the conversion factor between input and
;       output units (a power of ten) and the output unit. These can
;       be used for formatting tabular data etc.
;
; CATEGORY:
;       Input/Output, String Processing
;
;
; CALLING SEQUENCE:
;
;       Result = NUMUNIT(Number, Unit)
;
; INPUTS:
;       Number:  Number to be formatted (assumed float or double)
;       Unit:    String containing units of number, may already
;                contain a prefix
;
; KEYWORD PARAMETERS:
;       PRECISION:  Number of significant figures (integer, < 10).
;                   Default: 3.
;       DECIMALS:   Number of decimal places (integer, < 9)
;                   (alternative to PRECISION).
;       FORCE:      Re-scale even units which should not have
;                   prefixes, by including "10^n" in the unit string.
;
; OUTPUTS:
;      Result is a string containing rescaled number and unit.
;
; OPTIONAL OUTPUTS:
;       MULTIPLIER: Scale factor between input and output number:
;                   output = input * multiplier
;       OUT_UNIT:   String containing just the output unit (including
;                   prefix if any). 
;
; EXAMPLE:
;       Simple use:
;
;       IDL> PRINT, NUMUNIT(0.00001281,'  Ohms')
;       12.8 uOhm
;   
;       The illegitimate terminal "s" has been removed. (ohm-seconds
;       should be coded as "Ohm s", "Ohm.s" or "Ohm*s"). 
;
;       IDL> PRINT, NUMUNIT(1.7377e7,'JY/PIXEL', DEC=2)
;         17.38 MJy/pixel
;
;       (Such upper-case units may occur in old FITS headers)
;       Note that the leading blanks are not stripped when DECIMALS is
;       set, to facilitate lining results up in columns.
;
;       To choose a format and unit for a table, run NUMUNIT on the
;       maximum value to ensure there are no overflows, and specify
;       plenty of precision to allow for smaller values:
;
;       IDL> data = 100.0*EXP(RANDOMN(seed, 100))
;       IDL> PRINT, NUMUNIT(MAX(data),'Hz', PREC=5, MULT=mult, OUT=unit)
;       1.9922 kHz
;       IDL> PRINT, unit, mult
;       kHz   0.00100000
;       IDL> PRINT, 'Frequency', '('+unit+')', FORMAT = "(A)"
;       Frequency
;       (kHz)
;       IDL> PRINT, mult*data, FORMAT = "(F7.4)"
;        0.0763
;        0.0816
;       ...etc
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, Jan-Feb 2008
;
;-
COMPILE_OPT IDL2
breaks = ' *.,/^0123456789+-()'
known = ['m','g','s','rad','sr','K','A','mol','cd', $
         'Hz','J','W','V','N','Pa','C','Ohm','S','F','Wb','T','H','lm','lx', $
         'a','yr','eV','pc','Jy','mag','R','G','barn','bit','byte', $
         'metre', 'gram', 'second','sec', 'kelvin', 'ampere', 'mole', $
         'joule', 'watt', 'volt', 'newton', 'coulomb', 'siemens', 'farad', $
         'tesla', 'henry', 'jansky', 'rayleigh', 'gauss']
nopref = ['deg','arcmin','arcsec','mas','min','h','d','erg','Ry','solMass', $
          'u', 'solLum', 'Angstrom', 'solRad','AU','lyr','count','ct', $
          'photon','ph','pixel','pix','D','Sun','chan','bin','voxel','adu', $
          'beam', 'unknown', 'hour', 'day', 'rydberg', 'debye']
unambig = ['m','rad','sr','K','mol','cd', $
           'Hz','J','W','N','Pa','C','Ohm','F','Wb','T','lm','lx', $
           'yr','pc','Jy','mag','R','barn','bit','byte', $
           'metre', 'gram', 'second','sec', 'kelvin', 'ampere', 'mole', $
           'joule', 'watt', 'newton', 'coulomb', 'seivert', 'farad', 'tesla', $
           'henry', 'jansky', 'rayleigh', 'gauss', $
           'deg','arcmin','arcsec','mas','min','erg','Ry','solMass', $
           'u', 'solLum', 'Angstrom', 'solRad','AU','lyr','count','ct', $
           'photon','ph','pixel','pix','Sun','chan','bin','voxel','adu', $
           'beam','log','ln','exp','sqrt', 'hour', 'day', 'rydberg', 'debye']
ambig = ['g','G','s','S','A','a','H','h','d','D','eV']
prefix = ['y','z','a','f','p','n','u','m','','k','M','G','T','P','E','Z','Y']
altpref = ['c','d','','da','h']
prefuc = ['YOCTO','ZEPTO','ATTO','FEMTO','PICO','NANO','MICRO','MILLI','', $
          'KILO','MEGA','GIGA','TERA','PETA','EXA','ZETTA','YOTTA']
apfuc = ['CENTI','DECI','','DECA','HECTO']
functions = ['log','ln','exp','sqrt']

mk = WHERE( STRLEN(known) GT 1)
knowns = known[mk]+'s'
mk = WHERE( STRLEN(nopref) GT 1)
noprefs = nopref[mk]+'s'

prefout = ''
untail = ''
power = 0 ; Default prefix power of ten
number = numin
numout = number
premult = 0B

decset  = N_ELEMENTS(decimals) GT 0
IF decset  THEN IF decimals GT 9 THEN MESSAGE, 'Max 9 decimal places'
precset = N_ELEMENTS(precision) GT 0
IF precset THEN IF precision GT 10 THEN MESSAGE, 'Max 10 significant figures'
IF decset AND precset THEN MESSAGE, 'Set only one of DECIMALS or PRECISION'

IF decset + precset EQ 0 THEN BEGIN
    precset = 1B
    precision = 3
ENDIF

force = KEYWORD_SET(force)
;un = parse_unit(unit)
;FUNCTION parse_unit, unit
;
; Tries to interpret a unit string
;
un = STRTRIM(unit, 2)

; Check for power of ten in unit string
IF STRCMP(un,'10',2) THEN BEGIN
    get_power, STRMID(un,2),kpow, irest, valid
    IF valid EQ 2 THEN BEGIN
        number = number * (10d0)^kpow
        un = STRTRIM(STRMID(un,irest+2),1)
    ENDIF
ENDIF

inpower = number EQ 0.0 ? 0 : FLOOR(ALOG10(ABS(number)))

uc = un EQ STRUPCASE(un) ; If true, unit is in CAPS (very confusing!)

IF uc THEN BEGIN ; Try to convert to regular case where possible:
    ucp  = STRUPCASE([prefix, altpref])
    un0 = STRSPLIT(un, breaks, /EXTRACT, COUNT = count)
    IF count EQ 0 THEN BEGIN ; null unit
        untail = ''
        gotit = -1
        GOTO, NOSCALE
    ENDIF
    FOR j = 0,count-1 DO BEGIN
        gotit = -1
        unlen = STRLEN(un0[j])
        FOR i=0,N_ELEMENTS(unambig)-1 DO BEGIN
            start = STREGEX(un0[j],unambig[i]+'$',LENGTH = len, /FOLD_CASE)
            IF start+len NE unlen THEN CONTINUE
            IF start GT 0 THEN BEGIN
                pref0 = STRMID(un0[j],0,start)
                IF start GT 1 THEN BEGIN
                    gotit = WHERE(pref0 EQ prefuc)
                    IF gotit NE -1 THEN pref0 = prefix[gotit] ELSE BEGIN
                        gotit = WHERE(pref0 EQ apfuc)
                        IF gotit NE -1 THEN pref0 = altpref[gotit] ELSE BEGIN
                            gotit = WHERE(pref0 EQ 'DA')
                            IF gotit NE -1 THEN pref0 = 'da'
                        ENDELSE
                    ENDELSE
                ENDIF ELSE BEGIN
                                ; one letter prefix
                    gotit = WHERE(pref0 EQ ucp, ngot)
                    IF ngot EQ 1 THEN pref0 = ([prefix,altpref])[gotit]
                ENDELSE
                IF gotit[0] NE -1 THEN BEGIN
                    un0[j] = pref0+unambig[i]
                    BREAK
                ENDIF
            ENDIF ELSE IF start NE -1 THEN BEGIN
                un0[j] = unambig[i] ; no prefix
                gotit = 0
                BREAK
            ENDIF
        ENDFOR
        IF gotit[0] EQ -1 THEN GOTO, NOCHANGE ; never found clear match
    ENDFOR
                                ; Now put string back together
    start = STRSPLIT(un, breaks, LENGTH = len)
    IF start[0] GT 0 THEN new = STRMID(un,0,start[0]) ELSE new = ''
    FOR i =0,count-2 DO BEGIN
        gap0 = start[i] + len[i]
        lengap = start[i+1] - gap0
        new = new + un0[i] + STRMID(un,gap0,lengap)
    ENDFOR
    new = new + un0[count-1]
    gap0 = start[count-1]+len[count-1]
    lengap = STRLEN(un) - gap0
    IF lengap GT 0 THEN new = new + STRMID(un,gap0,lengap)
    un = new
    uc = 0B
ENDIF

NOCHANGE:

; Check for functions at start of unit string
brace = STRPOS(un,'(')
IF brace NE -1 THEN BEGIN
    un0 = STRSPLIT(un,'()',/EXTRACT)
    fun = WHERE(un0[0] EQ STRUPCASE(functions))
    IF fun NE -1 THEN BEGIN
        un0 = un
        GOTO, NOSCALE
    ENDIF
ENDIF

un0 = STRSPLIT(un, breaks, /EXTRACT)
un0 = un0[0]
untail = STRMID(un,STRLEN(un0))

; Check to see if leading unit is raised to a funny power
get_power, untail, unpower, irest, valid
IF valid AND unpower NE 1d0 THEN GOTO, NOSCALE ; Can't use prefix multipliers

gotit = WHERE(un0 EQ [known, nopref])
bare = gotit NE -1 OR STRLEN(un0) EQ 1

IF bare EQ 0B AND uc THEN BEGIN
    gotit = WHERE(un0 EQ STRUPCASE([known, nopref]))
    IF gotit[0] NE -1 THEN un0 = ([known, nopref])[gotit[0]]
    bare = gotit[0] NE -1
ENDIF

IF bare EQ 0B THEN BEGIN ; try unauthorised plurals
    gotit = WHERE(un0 EQ [knowns, noprefs])
    bare = gotit NE -1
    IF bare THEN un0 = STRMID(un0,0,STRLEN(un0)-1) ; Drop final 's'
ENDIF

IF bare EQ 0B THEN BEGIN ; unit not recognised, look for scaling prefix
    IF uc THEN BEGIN
        FOR i=0,16 DO BEGIN
            IF STRMID(un0,0,STRLEN(prefuc[i])) EQ prefuc[i] THEN BEGIN
                power = 3*(i - 8)
                un0 = STRMID(un0,STRLEN(prefuc[i]))
                GOTO, PREFSET
            ENDIF
        ENDFOR
        FOR i=0,4 DO BEGIN
            IF STRMID(un0,0,STRLEN(apfuc[i])) EQ apfuc[i] THEN BEGIN
                power = i-2
                un0 = STRMID(un0,STRLEN(apfuc[i]))
                GOTO, PREFSET
            ENDIF
        ENDFOR
    ENDIF

    pref = (WHERE(STRCMP(un0,prefix,1)))[0]
    IF pref NE -1 THEN power = 3*(pref - 8) ELSE BEGIN
        pref = WHERE(STRCMP(un0,altpref,1))
        IF pref[0] EQ 1 THEN pref = STRCMP(un0,'da',2) ? 3 : 1
        IF pref NE -1 THEN power = pref - 2
    ENDELSE
                                ; remove prefix from unit string
    plen = power EQ 1 ? 2 : 1
    oldpref = STRMID(un0,0,plen)
    un0 = STRMID(un0,plen)
ENDIF

gotit = (WHERE(un0 EQ [known,knowns]))[0]
;gotit = WHERE(un0 EQ nopref)
IF gotit GT -1 THEN GOTO, PREFSET

recog = (WHERE(un0 EQ [nopref, noprefs]))[0]

IF recog EQ -1 THEN BEGIN
    un0 = oldpref+un0           ; safer!
    power = 0
ENDIF

NOSCALE:
    IF force THEN BEGIN
        premult = 1B
        GOTO, PREFSET
    ENDIF

    prefout = ''
    numout = number
    IF inpower LT 0 OR inpower GE 3 THEN BEGIN
        IF precset THEN decimals = precision - 1
        fmt = STRING(8+decimals, decimals, FORMAT = "('(E',I2,'.',I1)")
        GOTO, FINISH
    ENDIF

PREFSET:
    number = number * 10d0^(power)
    pon3 = FLOOR((inpower+power)/3d0)
    IF premult THEN BEGIN
        pref_index = 0
        grams = gotit[0] EQ 1
        IF grams THEN BEGIN
            number = number * 1d-3
            inpower = inpower - 3
            pon3 = pon3 - 1
            un0 = 'kg'
        ENDIF
        pow =pon3 * 3
        IF pow EQ 0 THEN prefout = '' ELSE BEGIN
            prefout = '10^' + STRTRIM(STRING(pow),1) + ' '
        ENDELSE
    ENDIF ELSE BEGIN
        pref_index = pon3 + 8
        IF pref_index LT 0 OR pref_index GT 16 THEN BEGIN
            number = number * 10d0^(-power)
            GOTO, NOSCALE                           ; out of prefixes
        ENDIF
        prefout = prefix[pref_index]
    ENDELSE

    outpower = 3*pon3
    numout = number * 10d0^(-outpower)

    outlog = numout EQ 0.0 ? 0 : FLOOR(ALOG10(ABS(numout)))
    IF precset THEN decimals = precision - outlog - 1

                                ; Watch out for rounding up:
    frac = numout * 10d0^(decimals)
    frac = frac - FIX(frac)

    IF ABS(frac) GT 0.5 && FIX(ABS(numout))+1 GE 1000 THEN BEGIN
        IF pref_index + 1 GT 16 THEN BEGIN
            number = number * 10d0^(-power)
            GOTO, NOSCALE
        ENDIF
        numout  /= 1d3
        prefout = prefix[pref_index + 1]
        IF precset THEN decimals = precision - FLOOR(ALOG10(ABS(numout)+1)) - 1
    ENDIF

    fmt = decimals LE 0 ? "(I4" : $
      STRING(5+decimals, decimals, FORMAT = "('(F',I2,'.',I1)")

    IF decimals LE 0 THEN numout = ROUND(numout) 
FINISH:

out_unit = prefout + un0 + untail

fmt = fmt +",' ',A)"

numstring = STRING(numout, out_unit, FORMAT = fmt)
IF precset THEN numstring = STRTRIM(numstring, 2)

multexp = ROUND(ALOG10(numout / numin))
multiplier = 10.0^multexp

RETURN, numstring

END
