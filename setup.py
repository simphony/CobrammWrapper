from setuptools import setup, find_packages

from packageinfo import VERSION, NAME, OSP_CORE_MIN, OSP_CORE_MAX

# Read description
with open('README.md', 'r') as readme:
    README_TEXT = readme.read()


# main setup configuration class
setup(
    name=NAME,
    version=VERSION,
    author='Computational Chemistry CHIMIND, UNIBO',
    description='Simulation wrapper for COBRAMM',
    keywords='simphony, cuds, cobramm',
    long_description=README_TEXT,
    install_requires=[
        'osp-core>=' + OSP_CORE_MIN + ', <' + OSP_CORE_MAX,
        'numpy',
        'scipy',
        'matplotlib',
    ],
    packages=find_packages(),
    test_suite='tests',
    entry_points={
        'wrappers':
            'wrapper = osp.wrappers.cobrammwrapper:CobrammSession'},
)
