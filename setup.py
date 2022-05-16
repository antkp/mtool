
from setuptools import setup


setup(
    name="mtool-1.46_ALPHA",
    version="1.46",
    description="Read me",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/antkp/mtool.git",
    author="p. h.",
    author_email="saburow@gmx.net",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
    ],
    packages=["reader"],
    include_package_data=True,
    install_requires=[
        "feedparser", "html2text", "importlib_resources", "typing"
    ],
    entry_points={"console_scripts": ["realpython=reader.__main__:main"]},
)