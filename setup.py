# setup.py

from setuptools import setup, find_packages

setup(
    name='elastic_cp2k',
    version='0.1.0',
    description='A package for CP2K elastic constants calculation.',
    author='Keming Zhu',
    author_email='zhukeming0127@gmail.com',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'ase',
        'pandas',
        'pyyaml'
    ],
    entry_points={
        'console_scripts': [
            'elastic_cp2k=elastic_cp2k.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
    ],
)
