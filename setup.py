from setuptools import setup

setup(name='pyclamps',
      version='0.0.1',
      description='Package for processing CLAMPS data',
      url='',
      author='Tyler Bell',
      author_email='tyler.bell@ou.edu',
      license='MIT',
      packages=['pyclamps'],
      install_requires=[
          'cmocean',
          'matplotlib',
          'netCDF4',
          'numpy',
          'xarray',
      ],
      zip_safe=False)