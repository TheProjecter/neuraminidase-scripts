# BEGIN FUNCTIONS ###########
# Notation function for 3D structure
def notation_3D(j):
    #------------
    # notacion 3D
    if j<=169:
        # Normal
        new_notation = str(j+1)
    elif j==170:
        # Se resta 1 y se agrega 'a'
        new_notation = str(j)+'a'
    elif j<=307:
        # Se resta 1
        new_notation = str(j)
    elif j<=333:
        # Normal
        new_notation = str(j+1)
    elif j<=335:
        # Se suma 1
        new_notation = str(j+2)
    elif j<=389:
        # Se suma 3
        new_notation = str(j+4)
    elif j<=408:
        # Se suma 4
        new_notation = str(j+5)
    elif j<=412:
        # Se enumera como subtipos 412 (a,b,c,d)
        new_notation = str(412)+chr(j-312)
    elif j<=433:
        # Normal
        new_notation = str(j+1)
    else: 
        # Se suma 1
        new_notation = str(j+2)
    return new_notation

""" Determinate the region of the value using a regions_limits tuple and the value to locate
Args:
    value: value to determine the region
    regions_limits: tuple with the boundaries of the regions
Returns:
    An integer number with the ubication in the tuple or array??
"""
def whatRegion(value,regions_limits):
    for i in range(len(regions_limits)):
        if(value<=regions_limits[i]):
            return i    
    return len(regions_limits)



### Deprecated ####
# Notation function for active site a partir de la notacion_3D
# TODO: definir notacion de sitio activo usar sets.
def whatRegion_AS(x):
    
    if  x < 83 :
	    region = 0
    elif  x < 118 :
	    region = 1
    elif  x == 118 :
	    region = 2
    elif  x == 119 :
	    region = 3
    elif  x < 151 :
	    region = 4
    elif  x == 151 :
	    region = 5
    elif  x == 152 :
	    region = 6
    elif  x < 156 :
	    region = 7
    elif  x == 156 :
	    region = 8
    elif  x < 177 :
	    region = 9
    elif  x == 177 :
	    region = 10
    elif  x < 223 :
	    region = 11
    elif  x == 223 :
	    region = 12
    elif  x < 225 :
	    region = 13
    elif  x == 225 :
	    region = 14
    elif  x < 228 :
	    region = 15
    elif  x == 228 :
	    region = 16
    elif  x < 277 :
	    region = 17
    elif  x == 277 :
	    region = 18
    elif  x == 278 :
	    region = 19
    elif  x < 293 :
	    region = 20
    elif  x == 293 :
	    region = 21
    elif  x < 295 :
	    region = 22
    elif  x == 295 :
	    region = 23
    elif  x < 368 :
	    region = 24
    elif  x == 368 :
	    region = 25
    elif  x < 402 :
	    region = 26
    elif  x == 402 :
	    region = 27
    else     :
	    region = 28

    return region

### Deprecated ####    
# Secondary Structure region. i.e. alpha, beta, loop
def whatRegion_2D(x):
    if x <= 81 :
 	    region = 0
    elif x <= 95 :
 	    region = 1
    elif x <= 103 :
 	    region = 2
    elif x == 104 :
 	    region = 3
    elif x <= 109 :
 	    region = 4
    elif x <= 113 :
 	    region = 5
    elif x <= 118 :
 	    region = 6
    elif x <= 120 :
 	    region = 7
    elif x <= 125 :
 	    region = 8
    elif x <= 127 :
 	    region = 9
    elif x <= 136 :
 	    region = 10
    elif x == 137 :
 	    region = 11
    elif x <= 140 :
 	    region = 12
    elif x <= 155 :
 	    region = 13
    elif x <= 161 :
 	    region = 14
    elif x == 172 :
 	    region = 15
    elif x <= 178 :
 	    region = 16
    elif x == 179 :
 	    region = 17
    elif x <= 186 :
 	    region = 18
    elif x <= 188 :
 	    region = 19
    elif x <= 198 :
 	    region = 20
    elif x <= 200 :
 	    region = 21
    elif x <= 209 :
 	    region = 22
    elif x == 210 :
 	    region = 23
    elif x <= 218 :
 	    region = 24
    elif x <= 231 :
 	    region = 25
    elif x <= 234 :
 	    region = 26
    elif x <= 236 :
 	    region = 27
    elif x <= 244 :
 	    region = 28
    elif x <= 251 :
 	    region = 29
    elif x <= 259 :
 	    region = 30
    elif x <= 261 :
 	    region = 31
    elif x <= 269 :
 	    region = 32
    elif x <= 279 :
 	    region = 33
    elif x <= 284 :
 	    region = 34
    elif x <= 286 :
 	    region = 35
    elif x <= 292 :
 	    region = 36
    elif x <= 300 :
 	    region = 37
    elif x <= 307 :
 	    region = 38
    elif x <= 310 :
 	    region = 39
    elif x <= 317 :
 	    region = 40
    elif x <= 350 :
 	    region = 41
    elif x <= 353 :
 	    region = 42
    elif x <= 355 :
 	    region = 43
    elif x <= 361 :
 	    region = 44
    elif x <= 367 :
 	    region = 45
    elif x <= 377 :
 	    region = 46
    elif x <= 387 :
 	    region = 47
    elif x <= 393 :
 	    region = 48
    elif x <= 397 :
 	    region = 49
    elif x <= 400 :
 	    region = 50
    elif x == 401 :
 	    region = 51
    elif x <= 408 :
 	    region = 52
    elif x == 409 :
 	    region = 53
    elif x <= 413 :
 	    region = 54
    elif x <= 418 :
 	    region = 55
    elif x <= 430 :
 	    region = 56
    elif x <= 436 :
 	    region = 57
    elif x <= 448 :
 	    region = 58
    elif x <= 469 :
 	    region = 59
    else   :
         region = 60

    return region
# END FUNCTIONS ###########
