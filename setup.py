from setuptools import setup
from irap.__init__ import __version__

setup(name='irap',
    version=__version__,
    description='an intelligent protein analysis toolkit based on raac and pssm ',
    url='https://github.com/KingoftheNight/irap',
    author='liangyc',
    author_email='1694822092@qq.com',
    license='BSD 2-Clause',
    packages=['irap'],
    install_requires=[
        'numpy',
        'matplotlib',
        'scikit-learn',
        'seaborn==0.10.1',
        'pyecharts'
        ],
    entry_points={
        'console_scripts': [
        'irap=irap.__main__:rpct_main',
            ]
        },
    python_requires=">=3.6",
    include_package_data=True,
    zip_safe=True)
