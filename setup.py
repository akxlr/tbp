#!/usr/bin/env python

import os
import stat
import subprocess
from sys import platform
from setuptools import setup
from distutils.command.build import build
from shutil import copyfile

VERSION = (0, 1, 0)
BASEPATH = os.path.dirname(os.path.abspath(__file__))
LIBDAI_PATH = os.path.join(BASEPATH, 'libdai')
README = os.path.join(BASEPATH, 'README.md')


class TBPBuild(build):
    def run(self):
        build.run(self)

        # Build libDAI to produce libdai/utils/dfgmarg binary (dfgmarg is a C++ program that takes a .dfg file and
        # returns marginals)

        # First, copy platform-specific Makefile.<platform> to Makefile.conf
        if platform == "linux" or platform == "linux2":
            makefile = 'Makefile.LINUX'
        elif platform == "darwin":
            makefile = 'Makefile.MACOSX64'
        else:
            # libDAI has makefiles for most platforms, but these are untested so we raise an exception instead.
            raise Exception("Unexpected platform {} - are you using Linux or OSX?".format(platform))
        copyfile(os.path.join(LIBDAI_PATH, makefile), os.path.join(LIBDAI_PATH, 'Makefile.conf'))

        # Build libDAI with make
        print('Building libDAI')
        subprocess.run('make', cwd=LIBDAI_PATH)

        # Due to the way setuptools works, the binary needs to be contained in the actual tbp package directory in
        # order to be copied to the user's system by the package_data argument below - so copy it there
        dest = os.path.join(self.build_lib, 'tbp/dfgmarg')
        copyfile(os.path.join(LIBDAI_PATH, 'utils/dfgmarg'), dest)
        # Set executable permissions
        st = os.stat(dest)
        os.chmod(dest, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

def setup_tbp():

    # Convert README.md to the required rst format (only used when uploading to PyPI)
    try:
        from pypandoc import convert
        long_description = convert(README, 'rst')
    except ImportError:
        print("Warning: pypandoc module not found, could not convert Markdown to RST")
        long_description = open(README, 'r').read()

    setup(
        name='tbp',
        packages=['tbp'],
        entry_points={
            'console_scripts': ['tbp=tbp.cli:main'],
        },
        install_requires=[
            'numpy',
            'tensorly',
        ],
        # The only part of libDAI actually installed to the user's system is the dfgmarg binary - ensure the
        # binary is included in the installed package dir
        package_data={
            'tbp': ['dfgmarg'],
        },
        version='{}.{}.{}'.format(*VERSION),
        description='Tensor Belief Propagation - algorithm for approximate inference in discrete graphical models',
        long_description=long_description,
        license='MIT',
        author='Andrew Wrigley',
        author_email='andrew@wrigley.io',
        maintainer='Andrew Wrigley',
        maintainer_email='andrew@wrigley.io',
        url='https://github.com/akxlr/tbp',
        download_url='https://github.com/akxlr/tbp/archive/{}.{}.{}.tar.gz'.format(*VERSION),
        keywords=[
            'tensor',
            'belief',
            'propagation',
            'tensor belief propagation',
            'graphical models',
            'machine learning',
            'inference',
        ],
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Operating System :: Unix',
            'Programming Language :: C++',
            'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering :: Artificial Intelligence',
        ],
        cmdclass={
            'build': TBPBuild,
        }
    )


if __name__ == '__main__':
    setup_tbp()


