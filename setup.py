import setuptools  # this is the "magic" import
from setuptools import find_packages
from numpy.distutils.core import setup, Extension

from MapMaker import __version__


ext = Extension(name='MapMaker.Tools.fBinning', sources=['MapMaker/Tools/fBinning.f90'])

setup(
    name = 'MapMaker',
    version = __version__,

    packages = find_packages(),
    ext_modules = [ ext ],

    # metadata for upload to PyPI
    author = "Stuart Harper, Shifan Zuo",
    author_email = "sfzuo@bao.ac.cn",
    description = "Map-making package for single dish astronomical data.",
    # license = "GPL v3.0",
    url = "https://github.com/zuoshifan/MapMaker.git",
)