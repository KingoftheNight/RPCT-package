from setuptools import setup

setup(
    name="pyrpct",
    version="1.0",
    author="Liang YC",
    author_email="1694822092@qq.com",
    description="A RAAC-PSSM-based Protein Prediction Toolkit",
    long_description="none",
    long_description_content_type="text/markdown",
    url="https://github.com/KingoftheNight/RPCT-package",
    packages=["pyrpct"],
    install_requires=[
        'numpy',
        'matplotlib',
        'scikit-learn',
        'seaborn',
        'pyecharts'
        ],
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=True)
