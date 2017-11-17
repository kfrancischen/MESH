#!/bin/bash

OBJDIR="$1"
LIBFILE="$2"

cat <<SETUPPY > setup.py
from distutils.core import setup, Extension
MESHmodule = Extension('MESH',
	sources = [
		'src/main_python.cpp'
	],
	libraries = [
		'mesh',
		'stdc++',
    'openblas'
	],
	library_dirs = [
    '$OBJDIR'
  ],
	extra_link_args = [
		'$LIBFILE',
    '-lgomp'
	],
	include_dirs = [
		'src/arma'
	],
  extra_compile_args=[
    '-std=c++11', 
    '-fopenmp',
		'-Wno-strict-prototypes'
  ]
)
setup(name = 'MESH',
  license = "GNU GPL v3",
	version = '1.2.0',
  author='Kaifeng Chen',
  author_email='kfchen@stanford.edu',
  url='https://kfrancischen.github.io/MESH/',
	description = 'Multilayer Electromagnetic Solver for Heat transfer',
	ext_modules = [MESHmodule]
)
SETUPPY