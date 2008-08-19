from pydec.testing import *

from pydec.io.misc import *

def test_file_extension():
    cases = []
    cases.append(("",""))
    cases.append(("a",""))
    cases.append(("a.txt","txt"))
    cases.append(("helloworld",""))
    cases.append(("somefile.exe","exe"))
    cases.append(("one.two.three","three"))
    cases.append(("one.",""))
    cases.append(("one..two..three","three"))
    cases.append(("a-b-.2.34,5.3","3"))
    
    for f,e in cases:
        assert_equal(file_extension(f),e)
    
