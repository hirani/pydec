__all__ = ['read_mesh','write_mesh']

from pydec.mesh import simplicial_mesh
import pydec.io.arrayio
#from xml.dom.ext import PrettyPrint
from xml.dom import minidom
from io import IOBase

import os

class PyMeshException(Exception):  pass
class PyMeshIOException(PyMeshException): pass



mesh_str_type_pairs = [('simplicial_mesh',simplicial_mesh)]
mesh_str_to_type    = dict(mesh_str_type_pairs)
mesh_type_to_str    = dict([(t,s) for (s,t) in mesh_str_type_pairs])



def read_arrays(node_list,filepath):
    array_dict = dict()
    
    for node in node_list:
        file_node = node.getElementsByTagName('file')[0]
        file_name = str(file_node.attributes['name'].value)
        file_name = os.path.join(filepath,file_name)
        
        data = pydec.io.arrayio.read_array(file_name)

        array_dict[str(node.nodeName)] = data
        
    return array_dict



def read_mesh(fid):
    """
    Read a mesh from a given open file or filename.

    Examples:
      my_mesh = read_mesh('torus.xml')
    or
      fid = open('torus.xml')
      my_mesh = read_mesh(fid)
    
    """
    if not hasattr(fid, "read"): fid = open(fid)
        
    xmldoc = minidom.parse(fid)

    mesh_node = xmldoc.firstChild

    if mesh_node.tagName != 'mesh':
        raise PyMeshIOException('Invalid XML root node')


    (filepath, filename) = os.path.split(fid.name)   

    children = [child for child in xmldoc.firstChild.childNodes if child.nodeType == child.ELEMENT_NODE]
    
    array_dict = read_arrays(children, filepath)

    if mesh_node.hasAttribute('type'):
        mesh_str = str(mesh_node.attributes['type'].value)
    else:
        mesh_str = 'simplicial_mesh'

    mesh_type = mesh_str_to_type[mesh_str]

    return mesh_type(array_dict)


def write_mesh(fid, mesh, format='binary'):
    """
    Write a mesh to a given file or filename.


    Examples:
       write_mesh('torus.xml',my_mesh)
    or
       write_mesh('torus.xml',my_mesh,format='ascii')
    or
       fid = open('torus.xml')
       write_mesh(fid,my_mesh,format='basic')
       
    """
    
    if not hasattr(fid, "read"): fid = open(fid, 'w')
    
    (filepath, filename) = os.path.split(fid.name)
    basename = filename.split('.')[0]
    
    
    xmldoc = minidom.Document()
    
    mesh_node = xmldoc.appendChild(xmldoc.createElement('mesh'))
    mesh_node.setAttribute('type', mesh_type_to_str[type(mesh)])
    for key,value in mesh.items():
        data_filename = basename + '.' + key
        
        data_node      = mesh_node.appendChild(xmldoc.createElement(key))
        data_file_node = data_node.appendChild(xmldoc.createElement('file'))
        data_file_node.setAttribute('name',data_filename)
        
        pydec.io.arrayio.write_array(os.path.join(filepath,data_filename),value,format)


    xmldoc.writexml(fid,indent='',addindent='\t',newl='\n')



