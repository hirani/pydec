__all__ = ['flatten']      
      
      
def flatten(l, ltypes=(list, tuple)):
    """
    Recursively flatten a list
    """
    i = 0
    while i < len(l):
        if isinstance(l[i],ltypes) and not l[i]: # skip empty lists/tuples
            l.pop(i)
            continue
        while isinstance(l[i], ltypes):
            l[i:i+1] = list(l[i])
        i += 1
    return l
