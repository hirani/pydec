#!/usr/bin/env python3
"""
A setup.py script to use setuptools, which gives egg goodness, etc.
"""

from setuptools import setup
with open('setup.py') as f:
    code = compile(f.read(), 'setup.py', 'exec')
    exec(code)
    
