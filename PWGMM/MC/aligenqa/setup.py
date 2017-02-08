from setuptools import setup
from glob import glob

install_requires = ['rootpy==0.8.2', 'argparse']
# tests_require = ['nose']

setup(
    name='aligenqa',
    version='0.0.1',
    description="Provides tools for post-processing of MC multiplicity estimator studies",
    author='Christian Bourjau',
    author_email='christian.bourjau@cern.ch',
    packages=['aligenqa'],
    long_description=open('README.md').read(),
    scripts=glob('scripts/*'),
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Development Status :: 3 - Alpha",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    dependency_links=[
        'https://github.com/rootpy/rootpy/archive/84325d6ceae7858f284e4d92b8db2acda6993fa6.zip#egg=rootpy-0.8.2.dev0'
    ],
    install_requires=install_requires,
    # extras_require={'test': tests_require},
    test_suite='nose.collector',
)
