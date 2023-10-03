import numpy as np

# write an entry of a Fortran binary file
def writeentry(f,e,dtype,endian):
    if dtype == 'i4':
        e = np.array( e, dtype = 'i4' )
    if endian == 'big':
        #b = e.astype( '>'+dtype ).tobytes( order = 'F' )
        b = e.astype( '>'+dtype ).tostring( order = 'F' )
        l = np.array( len(b), dtype = '>i4' ) # entry length
    else:
        #b = e.astype( dtype ).tobytes( order = 'F' )
        b = e.astype( '<'+dtype ).tostring( order = 'F' )
        l = np.array( len(b), dtype = 'i4' ) # entry length
    f.write( l )
    f.write( b )
    f.write( l )


