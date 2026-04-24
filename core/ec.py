from core.fieldArith import modinv
from gmpy2 import f_mod, mul
from base64 import urlsafe_b64encode as enc64
from base64 import urlsafe_b64decode as dec64

# P is a very large prime number that defines the field boundaries.
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
b = 7  # b=7 comes from the curve equation: y^2 = x^3 + 7 
# G is the generator point - the starting point on the curve everyone agrees on.(industry standard)
x = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
y = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
G = (x, y)
# n is the total number of points available on the curve .
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

def pointAdd(p: tuple, q: tuple):
    #addition -   In ECC, adding two points isn't simple math. You draw a line through 
    #P and Q, find where it hits the curve again, and reflect it on x - axis .
    
    if p == None: return q # Adding nothing to Q
    if q == None: return p # Adding nothing to P

    lam = 0 # Lambda represents the slope of the line
    if p[0] == q[0]:
        if p[1] != q[1]: # Points are vertical mirror images so when added they cancel out
            return None
        # case : if two points conincide then we find tangent at that point.
        lam = f_mod((3 * mul(p[0], p[0])) * modinv(2 * p[1], P), P)
    else:
        # case: 2 differnet points -  Find the slope between two different points.
        lam = f_mod((q[1] - p[1]) * modinv(q[0] - p[0], P), P)

    # Calculate the new point (x, y) based on the slope
    x_new = f_mod(mul(lam, lam) - p[0] - q[0], P)
    y_new = f_mod(mul(lam, p[0] - x_new) - p[1], P)
    return (x_new, y_new)

def scalMul(k: int, p: tuple):
    #Multiplies a point P by a number k by double And add algorithm which 
    #is immposible to undo
    
    result = None
    binSum = p
    while(k > 0):
        if((k & 1) == 1): # If the bit is 1, add the current power-of-two point
            result = pointAdd(binSum, result)
        binSum = pointAdd(binSum, binSum) # Double the point for the next bit
        k >>= 1
    return result

def onCurv(p: tuple):
    #Checks if a point actually on curve y^2 = x^3 + 7.
    left = mul(p[1], p[1])
    right = mul(mul(p[0], p[0]), p[0]) + b
    return f_mod(left - right, P) == 0

def xyTob64(p: tuple):
    #Converts a (x, y) math coordinate into a Base64 string.
    raw_bytes = p[0].to_bytes(32, 'big') + p[1].to_bytes(32, 'big')
    return str(enc64(raw_bytes).decode('utf-8'))

def xyFromb64(raw: str):
    #Turns a Base64 string back into (x, y) math coordinates.
    raw_bytes = dec64(raw)
    return (int.from_bytes(raw_bytes[:32], 'big'), int.from_bytes(raw_bytes[32:], 'big'))
