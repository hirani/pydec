__all__ = ['file_extension']

def file_extension(filename):
    """
    Return the extension of a file (if any)
    """
    parts = filename.split(".")
    if len(parts) >= 2:
        return parts[-1]
    else:
        return ""
        
        
