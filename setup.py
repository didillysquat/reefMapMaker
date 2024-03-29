try:
    from setuptools_conda import dist_conda
    cmdclass = {'dist_conda': dist_conda}
except ImportError:
    cmdclass = {}
import setuptools
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="reefmapmaker",
    version='0.1.10',
    author="Benjamin C C Hume",
    author_email="didillysquat@gmail.com",
    description="Script to plot maps with reference coral reefs annotated.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/didillysquat/reefMapMaker",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux"
    ],
    license='GPL-3.0-only',
    python_requires='>=3.6',
    install_requires=['cartopy', 'pandas', 'matplotlib', 'numpy', 'xlrd', 'scipy'],
    scripts=['scripts/reefmapmaker'],
    project_urls={
        'Bug Reports': 'https://github.com/didillysquat/reefMapMaker/issues',
        'Source': 'https://github.com/didillysquat/reefMapMaker',
    }
)
