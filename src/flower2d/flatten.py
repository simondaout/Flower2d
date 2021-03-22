import numpy as np
import math,sys

# Returns a flat list containing all elements in the sequence in
# depth-first order
#
# By default, it flattens all the sequence recursively.  Set
# recursive to false, in order to only flatten one level of
# elements

def flatten(sequence,recursive=True):
    result=[]
    if recursive:
        stack=[]
        i = iter(sequence)
        while True:
            try: 
                e=next(i)
                if hasattr(e,"__iter__") and not isinstance(e,str):
                    stack.append(i)
                    i=iter(e)
                else:
                    result.append(e)
            except StopIteration:
                try:
                    i = stack.pop()
                except IndexError:
                    return result
    
    else:
        for e in sequence:
            if hasattr(e,"__iter__") and not isinstance(e,str):
                result.extend(e)
            else:
                result.append(e)
        return result        





