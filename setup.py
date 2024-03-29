from setuptools import setup

setup(
    name='NESSI',
    version='0.0.1',    
    description='Numerical-Empirical-Sun-as-a-Star-Integrator',
    url='https://github.com/apy-github/NESSI',
    author='Adur Pastor Yabar',
    author_email='adur.pastor@astro.su.se',
    license='GPLv3',
    packages=['nessi',],
    package_data={'': ['data/*.npz', 'data/*.fits']},
    install_requires=['scipy',
                      'numpy',
                      'matplotlib',
                      'astropy',
                      ],
    )

#    package_dir={'': ''},
#    classifiers=[
#        'Development Status :: 1 - Planning',
#        'Intended Audience :: Science/Research',
#        'License :: OSI Approved :: BSD License',  
#        'Operating System :: POSIX :: Linux',        
#        'Programming Language :: Python :: 2',
#        'Programming Language :: Python :: 2.7',
#        'Programming Language :: Python :: 3',
#        'Programming Language :: Python :: 3.4',
#        'Programming Language :: Python :: 3.5',
#    ],

