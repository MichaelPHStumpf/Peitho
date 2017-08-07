from setuptools import setup

setup(name='peitho',
	version='0.1',
	description='Perfecting Experiments with Information Theory',
	url='https://github.com/MichaelPHStumpf/Peitho',
	download_url='https://github.com/MichaelPHStumpf/Peitho/archive/0.1.tar.gz',
	author='Leander Dony, Scott Ward, Jonas Mackerodt',
	author_email='m.stumpf@imperial.ac.uk',
	license='MIT',
	packages=['peitho'],
	install_requires = [
		'pycuda',
		'numpy',
		'matplotlib',
		'python-libsbml'
	],
	entry_points = {
		'console_scripts': ['peitho=peitho.main.main:main']
		},
	keywords = ['information theory','entropy','experimental design'],
	classifiers = ['Development Status :: 3 - Alpha',
	'Environment :: Console',
	'Intended Audience :: Science/Research',
	'Topic :: Scientific/Engineering :: Bio-Informatics',
	'License :: OSI Approved :: MIT License',
	'Programming Language :: Python :: 2.7'],
	zip_safe=False)
