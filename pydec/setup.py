#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pydec', parent_package, top_path)

    config.add_subpackage('dec')
    config.add_subpackage('fem')
    config.add_subpackage('io')    
    config.add_subpackage('math')
    config.add_subpackage('mesh')
    config.add_subpackage('util')
    config.add_subpackage('testing')
    config.add_subpackage('vis')

    #config.make_svn_version_py()  # installs __svn_version__.py
    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
