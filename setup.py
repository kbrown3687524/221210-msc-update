from setuptools import setup
setup(
    name='amia',
    version='0.1.0',    
    description='Automated Mutation Introduction and Analysis',
    author='Keaghan Brown',
    author_email='3687524@myuwc,.ac.za',
    packages=['amia'],
    install_requires=['biopython>=1.79',
                      'tabulate>=0.8.10',
                      'pymol>=2.3.0',
                      'MDAnalysis[analysis]>=2.2.0',
                      'MDAnalysisTests',
                      'pandas>=1.5.1'
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.10.2',
    ],
)
